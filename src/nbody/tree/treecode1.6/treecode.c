/*
 * treecode.c: hierarchical N-body simulation code.
 * Copyright (c) 2015 by Joshua E. Barnes, Kyoto, Japan.
 */

#define global				// don't default global to extern
#
#include "treecode.h"

#include <sys/time.h>
#include <sys/times.h>
#include <sys/resource.h>

 #include <immintrin.h>

//  Default values for input parameters.
//  ____________________________________

string defv[] = {
#if defined(QUICKSCAN)
				";Hierarchical N-body code (quick scan)",
#elif defined(EXTGRAV)
				";Hierarchical N-body code (+ ext grav)",
#else
				";Hierarchical N-body code"
#  if defined(SOFTCORR)
                                " (soft corr)",
#  else
                                " (no soft corr)",
#  endif
#endif
    "in=",			";Input file with initial conditions",
    "out=",			";Output file pattern for N-body frames.",
				";Directive (eg, %d) formats step number.",
    "restore=",			";Continue run from existing state file",
    "save=",			";State file pattern, written each step.",
				";Directive (eg, %d) alternates 0 and 1.",
    "dtime=1/32",		";Timestep for leapfrog integration",
    "eps=0.0125",		";Smoothing length for force calculation",
#if !defined(QUICKSCAN)
    "theta=1.0",		";Hierarchical force accuracy parameter",
#endif
    "usequad=true",		";If false, don't use quadrapole moments",
    "options=",			";Comma-separated list of program options.",
				";new-tout: reschedule output times.",
				";reset-time: set time to zero.",
				";bh86,sw94,theta-eff: opening criteria.",
    "outputs=" PosTag "," VelTag,
				";Body data arrays written to output.",
				";Other choices: "
				  MassTag "," PhiTag "," AccTag ".",
    "tstop=2.0",		";Time to stop integration",
    "dtout=1/4",		";Time interval between data outputs",
    "nstatic=0",		";Number of static bodies in array.",
				";If pos (neg), count from start (end).",
    "nbody=65536",		";Number of bodies generated for test run.",
				";If no input file given, make Plummer model.",
    "timesteps=100",    ";Number of timesteps to run for.",
    "seed=123",			";Random number seed for test run",
    "stream=",			";Output file pattern for frame stream",
    "log=",			";Output file name for calculation log.",
				";Defaults to stdout unless stream=-.",
#if defined(EXTGRAV)
    "gravgsp=",			";Input GSP for external gravity field",
#endif
    "VERSION=1.6",		";Joshua Barnes  7 May 2015",
    NULL,
};

//  Prototypes for local procedures.
//  ________________________________

local void startrun(void);			// initialize system state
local void newrun(void);			// start new simulation
local void oldrun(void);			// resume old simulation
local void testdata(void);			// generate test data
local void stepsystem(void);			// advance by one time-step
local void treeforce(void);			// do force calculation

#if defined(EXTGRAV)
#include "gsp.h"
gsprof *gravgsp = NULL;				// GSP for external field
#endif

/**
 * Return time in seconds in Epoch.
 */
double wtime() {
  struct rusage r;
  struct timeval t;
  gettimeofday( &t, (struct timezone *)0 );
  return t.tv_sec + t.tv_usec*1.0e-6;
}

//  main: toplevel routine for hierarchical N-body code.
//  ____________________________________________________

int main(int argc, string argv[]) {
  initparam(argv, defv);			// initialize param access
  headline = defv[0] + 1;			// use default headline
  startrun();					// get params & input data
  startoutput();				// activate output code
  if (nstep == 0) {				// if data just initialized
    treeforce_initial_0 = wtime();
    treeforce();				// calculate initial forces
    treeforce_initial_1 = wtime();
    output();					// generate initial output
  }
  if (dtime != 0.0) {		// if time steps requested
    // TODO: make this work in timesteps?
    treeforce_0 = wtime();
    while (nstep <= timesteps) { // while not past tstop
      stepsystem();       // advance step by step
      output();         // output results each time
    }
    // while (tstop - tnow > 0.01 * dtime) {	// while not past tstop
    //   stepsystem();				// advance step by step
    //   output();					// output results each time
    // }
    treeforce_1 = wtime();
  }
  finaloutput();
  bodyptr p;
  bodyptr q;
  float phi = 0.0f;
  for (p = bodytab; p < bodytab+nbody; p++) {// loop over all bodies
    for (q = bodytab; q < bodytab+nbody; q++) {// loop over all bodies
    //  printf("Pos(p) = (%.8f,%.8f,%.8f)\n",  Pos(p)[0], Pos(p)[1], Pos(p)[2]);
    //  printf("Pos(q) = (%.8f,%.8f,%.8f)\n",  Pos(q)[0], Pos(q)[1], Pos(q)[2]);
      float rx = Pos(q)[0] - Pos(p)[0];
      float ry = Pos(q)[1] - Pos(p)[1];
      float rz = Pos(q)[2] - Pos(p)[2];
      float r2 = rx*rx + ry*ry + rz*rz + eps;
      float r2inv = 1.0 / sqrt(r2);
      float r6inv = r2inv * r2inv * r2inv;
      float mass = Mass(q);
      phi += mass * r6inv;
    }
  }
  printf(" Answer = %f\n", phi);
  return (0);					// end with proper status
}

//  startrun: startup hierarchical N-body code.
//  ___________________________________________
// This runs once. 
local void startrun(void) {
  startrun_time_0 = wtime();
  bodyptr p1, p2, p;
  stream gravstr;

  define_body(sizeof(body), Precision, NDIM);	// setup phat body struct
  define_body_offset(PosTag,  BodyOffset(Pos));
  define_body_offset(VelTag,  BodyOffset(Vel));
  define_body_offset(MassTag, BodyOffset(Mass));
  define_body_offset(PhiTag,  BodyOffset(Phi));
  define_body_offset(AccTag,  BodyOffset(Acc));
  infile = getparam("in");			// set I/O file names
  outfile = getparam("out");
  savefile = getparam("save");
  if (strnull(getparam("restore")))	{	// starting a new run?
    newrun();
  } else {						// else resume old run
    oldrun();
  }
  if (ABS(nstatic) > nbody) {			// check nstatic is OK
    error("%s: absurd value for nstatic\n", getargv0());
  }
  p1 = bodytab + MAX(nstatic, 0);		// set dynamic body range
  p2 = bodytab + nbody + MIN(nstatic, 0);
  testcalc = TRUE;				// determine type of calc:
  for (p = p1; p < p2; p++) {
    testcalc = testcalc && (Mass(p) == 0);	// look for dynamic masses
  }
  strfile = getparam("stream");
  logfile = getparam("log");
#if defined(EXTGRAV)
  if (! strnull(getparam("gravgsp"))) {		// was GSP file given?
    gravstr = stropen(getparam("gravgsp"), "r");
    get_history(gravstr);
    gravgsp = get_gsprof(gravstr);		// read external field GSP
    strclose(gravstr);
  }
#endif

  startrun_time_1 = wtime();
}

//  newrun, oldrun: initialize parameters and bodies.
//  _________________________________________________

local void newrun(void) {
  eps = getdparam("eps");			// get input parameters
  dtime = getdparam("dtime");
  nstatic = getiparam("nstatic");
#if !defined(QUICKSCAN)
  theta = getdparam("theta");
#endif
  usequad = getbparam("usequad");
  tstop = getdparam("tstop");
  dtout = getdparam("dtout");
  options = getparam("options");
  outputs = getparam("outputs");
  if (! strnull(infile)) {			// if data file was given
    inputdata();				// then read inital data
  } else {					// else make initial data
    nbody = getiparam("nbody");			// get number of bodies
    // I added this
    timesteps = getiparam("timesteps");     // get number of timesteps
    init_random(getiparam("seed"));		// set random number gen.
    testdata();					// and make plummer model
  }
  rsize = 1.0;					// start root w/ unit cube
  nstep = 0;					// begin counting steps
  tout = tnow;					// schedule first output
}

local void oldrun(void) {
  restorestate(getparam("restore"));		// read in old state file
  if (getparamstat("eps") & ARGPARAM)	{	// was eps given new value?
    eps = getdparam("eps");			// use command line value
  }
  if (getparamstat("nstatic") & ARGPARAM)	{ // likewise for others...
    nstatic = getiparam("nstatic");
  }
#if !defined(QUICKSCAN) 
  if (getparamstat("theta") & ARGPARAM) {
    theta = getdparam("theta");
  }
#endif
  if (getparamstat("usequad") & ARGPARAM) {
    usequad = getbparam("usequad");
  }
  if (getparamstat("options") & ARGPARAM) {
    options = getparam("options");
  }
  if (getparamstat("outputs") & ARGPARAM) {
    outputs = getparam("outputs");
  }
  if (getparamstat("tstop") & ARGPARAM) {
    tstop = getdparam("tstop");
  }
  if (getparamstat("dtout") & ARGPARAM) {
    dtout = getdparam("dtout");
  }
  if (scanopt(options, "new-tout")) {		// if output time reset
    tout = tnow + dtout;			// then offset from now
  }
}

//  testdata: generate Plummer model initial conditions for test runs,
//  scaled to units such that M = -4E = G = 1 (Henon, Hegge, etc).
//  See Aarseth, SJ, Henon, M, & Wielen, R (1974) Astr & Ap, 37, 183.
//  __________________________________________________________________

#define MFRAC  0.999				// cut off 1-MFRAC of mass

local void testdata(void) {
  real rsc, vsc, r, v, x, y;
  bodyptr p;

  float scale = 1.0f;
  float vscale = 1.0f;
  float mscale = 1.0f;

  if (nbody < 1) {			// check for silly values
    error("%s: absurd value for nbody\n", getargv0());
  }
  bodytab = (bodyptr) allocate(nbody * sizeof(body));
						// alloc space for bodies
  rsc = (3 * PI) / 16;				// set length scale factor
  vsc = rsqrt(1.0 / rsc);			// find speed scale factor
  int i = 0;
  masses = malloc(sizeof(int)*nbody);
  for (p = bodytab; p < bodytab+nbody; p++) {	// loop over particles
    Type(p) = BODY;				// tag as a body
    //Mass(p) = 1.0 / nbody;			// set masses equal
    // Set mass randomly, like in brute
    float mass = (rand() / (float) RAND_MAX) * mscale;
    //float mass = 0;
    Mass(p) = mass;
    masses[i] = mass;
    x = xrandom(0.0, MFRAC);			// pick enclosed mass
    r = 1.0 / rsqrt(rpow(x, -2.0/3.0) - 1);	// find enclosing radius
    pickshell(Pos(p), NDIM, rsc * r);		// pick position vector
    do {					// select from fn g(x)
      x = xrandom(0.0, 1.0);			// for x in range 0:1
      y = xrandom(0.0, 0.1);			// max of g(x) is 0.092
    } while (y > x*x * rpow(1 - x*x, 3.5));	// using von Neumann tech
    v = x * rsqrt(2.0 / rsqrt(1 + r*r));	// find resulting speed
    pickshell(Vel(p), NDIM, vsc * v);		// pick velocity vector
    i++;
  }
  tnow = 0.0;					// set elapsed model time
}

//  stepsystem: advance N-body system using simple leap-frog.
//  _________________________________________________________

local void stepsystem(void) {
  bodyptr p1, p2, p;

  p1 = bodytab + MAX(nstatic, 0);		// set dynamic body range
  p2 = bodytab + nbody + MIN(nstatic, 0);
  for (p = p1; p < p2; p++) {			// loop over body range	
    ADDMULVS(Vel(p), Acc(p), 0.5 * dtime);	// advance v by 1/2 step
    ADDMULVS(Pos(p), Vel(p), dtime);		// advance r by 1 step
  }
  treeforce();
  for (p = p1; p < p2; p++) {			// loop over body range	
    ADDMULVS(Vel(p), Acc(p), 0.5 * dtime);	// advance v by 1/2 step
  }
  nstep++;					// count another time step
  tnow = tnow + dtime;				// finally, advance time
}

//  treeforce: supervise force calculation.
//  _______________________________________
local void treeforce(void) {
  bodyptr p1, p2, p;
  real r, mr3i;

  p1 = bodytab + MAX(nstatic, 0);		// set dynamic body range
  p2 = bodytab + nbody + MIN(nstatic, 0);
  for (p = bodytab; p < bodytab+nbody; p++)	// loop over all bodies
    Update(p) = (testcalc ? p1 <= p && p < p2 : TRUE);
						// flag bodies to update
  maketree(bodytab, nbody);			// construct tree structure
  gravcalc();					// compute current forces
  forcereport();				// print force statistics
#if defined(EXTGRAV)
  for (p = bodytab; p < bodytab+nbody; p++)	// loop over all bodies
    if (Update(p) && gravgsp != NULL) {		// update in extern field?
      r = absv(Pos(p));				// get distance from origin
      mr3i = - mass_gsp(gravgsp, r) / (r*r*r);
      ADDMULVS(Acc(p), Pos(p), mr3i);		// add extern acc and phi
      Phi(p) += phi_gsp(gravgsp, r);
    }
#endif
}
