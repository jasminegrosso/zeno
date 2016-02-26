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

global float *masses;
global float *x_pos;
global float *y_pos;
global float *z_pos;
global float *vx;
global float *vy;
global float *vz;
global float *ax;
global float *ay;
global float *az;

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
  global int x[nbody];
  global int y[nbody];
  global int z[nbody];
  startoutput();				// activate output code
  if (nstep == 0) {				// if data just initialized
    treeforce_initial_0 = wtime();
    treeforce();				// calculate initial forces
    treeforce_initial_1 = wtime();
    output();					// generate initial output
  }
  if (dtime != 0.0)				// if time steps requested
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
  printf("startrun\n");
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
  if (strnull(getparam("restore")))		// starting a new run?
    newrun();
  else						// else resume old run
    oldrun();
  if (ABS(nstatic) > nbody)			// check nstatic is OK
    error("%s: absurd value for nstatic\n", getargv0());
  p1 = bodytab + MAX(nstatic, 0);		// set dynamic body range
  p2 = bodytab + nbody + MIN(nstatic, 0);

  testcalc = TRUE;				// determine type of calc:
  //START OF FIRST OPTIMISATION
 /* int vector_size = (nbody/4)*4;
  p = p1;

  __m128 tc = _mm_setzero_ps();

  float zero = 0.0;
  __m128 zero_v = _mm_load_ps1(&zero);
  int new = 0;

  float *a = _mm_malloc(sizeof(float) * nbody,64);

  for (int i = 0; i < vector_size; i+=4) {
    // Load masses into vector.
    __m128 mass_v = _mm_load_ps(masses+i);

    // Want to compare the masses with 0, if true store true in testcalc.
    // 0xffffffff : 0 (true:false)
    __m128 result = _mm_cmpeq_ps(mass_v, zero_v);
    // Add up all the values in the result, if it's equal to -4, 
    // then all of the values of Mass(p) evaluated to true. If not, 
    // one of them was false so testcalc should be false. 
    int f[4];
    _mm_storeu_ps((float*) f, result);

    int total = (int)f[0] + (int)f[1] + (int)f[2] + (int)f[3];
    // If masses don't add up to -4, then one of them wasn't 0, so
    // break out of loop. 
    if (total != -4) {
      testcalc = FALSE;
      break;
    }

    p += 4;
  }
  _mm_free(a);
  if (testcalc) {
    for (; p < p2; p++) {
      testcalc = testcalc && (Mass(p) == 0);
    }
  }*/
  //END OF FIRST OPTIMISAION
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
  if (! strnull(infile))			// if data file was given
    inputdata();				// then read inital data
  else {					// else make initial data
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
  if (getparamstat("eps") & ARGPARAM)		// was eps given new value?
    eps = getdparam("eps");			// use command line value
  if (getparamstat("nstatic") & ARGPARAM)	// likewise for others...
    nstatic = getiparam("nstatic");
#if !defined(QUICKSCAN)
  if (getparamstat("theta") & ARGPARAM)
    theta = getdparam("theta");
#endif
  if (getparamstat("usequad") & ARGPARAM)
    usequad = getbparam("usequad");
  if (getparamstat("options") & ARGPARAM)
    options = getparam("options");
  if (getparamstat("outputs") & ARGPARAM)
    outputs = getparam("outputs");
  if (getparamstat("tstop") & ARGPARAM)
    tstop = getdparam("tstop");
  if (getparamstat("dtout") & ARGPARAM)
    dtout = getdparam("dtout");
  if (scanopt(options, "new-tout"))		// if output time reset
    tout = tnow + dtout;			// then offset from now
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

  if (nbody < 1)				// check for silly values
    error("%s: absurd value for nbody\n", getargv0());
  bodytab = (bodyptr) allocate(nbody * sizeof(body));
						// alloc space for bodies
  rsc = (3 * PI) / 16;				// set length scale factor
  vsc = rsqrt(1.0 / rsc);			// find speed scale factor
  int i = 0;
  masses = malloc(sizeof(int)*nbody);
  x_pos = malloc(sizeof(float)*nbody);
  y_pos = malloc(sizeof(float)*nbody);
  z_pos = malloc(sizeof(float)*nbody);
  vx = malloc(sizeof(float)*nbody);
  vy = malloc(sizeof(float)*nbody);
  vz = malloc(sizeof(float)*nbody);
  ax = malloc(sizeof(float)*nbody);
  ay = malloc(sizeof(float)*nbody);
  az = malloc(sizeof(float)*nbody);

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
    //pickshell(Pos(p), NDIM, rsc * r);		// pick position vector
   
    real rad = rsc * r;
    real rsq, rscale;
    int j;

    do {
      rsq = 0.0;
      for (j = 0; j < NDIM; j++) {
        Pos(p)[j] = xrandom(-1.0, 1.0);
        rsq = rsq + Pos(p)[j] * Pos(p)[j];
      }
    } while (rsq > 1.0);
    rscale = rad / rsqrt(rsq);
    for (j = 0; j < NDIM; j++) {
      Pos(p)[j] = Pos(p)[j] * rscale;
    }

    x_pos[i] = Pos(p)[0];
    y_pos[i] = Pos(p)[1];
    z_pos[i] = Pos(p)[2];

    do {					// select from fn g(x)
      x = xrandom(0.0, 1.0);			// for x in range 0:1
      y = xrandom(0.0, 0.1);			// max of g(x) is 0.092
    } while (y > x*x * rpow(1 - x*x, 3.5));	// using von Neumann tech
    v = x * rsqrt(2.0 / rsqrt(1 + r*r));	// find resulting speed

    //pickshell(Vel(p), NDIM, vsc * v);		// pick velocity vector

    rad = vsc * v;

    do {
      rsq = 0.0;
      for (j = 0; j < NDIM; j++) {
        Vel(p)[j] = xrandom(-1.0, 1.0);
        rsq = rsq + Vel(p)[j] * Vel(p)[j];
      }
    } while (rsq > 1.0);
    rscale = rad / rsqrt(rsq);
    for (j = 0; j < NDIM; j++) {
      Vel(p)[j] = Vel(p)[j] * rscale;
    }

    vx[i] = Vel(p)[0];
    vy[i] = Vel(p)[1];
    vz[i] = Vel(p)[2];

    i++;
  }
  // printf("Positions initialised to: \n");
  // for (int k = 0; k < nbody; k++) {
  //   printf("%f %f %f\n", x_pos[k], y_pos[k], z_pos[k]);
  // }
  // printf("Velocities initialised to: \n");
  // for (int k = 0; k < nbody; k++) {
  //   printf("%f %f %f\n", vx[k], vy[k], vz[k]);
  // }
  // printf("Masses initialised to: \n");
  // for (int j = 0; j < nbody; j++) {
  //   printf("%f\n", masses[j]);
  // }
  tnow = 0.0;					// set elapsed model time
}

//  stepsystem: advance N-body system using simple leap-frog.
//  _________________________________________________________

local void stepsystem(void) {
  bodyptr p1, p2, p;

  p1 = bodytab + MAX(nstatic, 0);		// set dynamic body range
  p2 = bodytab + nbody + MIN(nstatic, 0);
  p = p1;

  int vector_size = (nbody/4)*4;

  float scale = (0.5 * dtime);
  __m128 scale_v = _mm_load_ps1(&scale);

  __m128 dtime_v = _mm_load_ps1(&dtime);

  for (int i = 0; i < vector_size; i+=4) {
    __m128 x_pos_v = _mm_load_ps(x_pos+i); 
    __m128 y_pos_v = _mm_load_ps(y_pos+i);
    __m128 z_pos_v = _mm_load_ps(z_pos+i);

    __m128 vx_v = _mm_load_ps(vx+i); 
    __m128 vy_v = _mm_load_ps(vy+i);
    __m128 vz_v = _mm_load_ps(vz+i);

    __m128 ax_v = _mm_load_ps(ax+i); 
    __m128 ay_v = _mm_load_ps(ay+i);
    __m128 az_v = _mm_load_ps(az+i);

    __m128 calculated_vel_x = _mm_mul_ps(ax_v, scale_v);
    __m128 calculated_vel_y = _mm_mul_ps(ay_v, scale_v);
    __m128 calculated_vel_z = _mm_mul_ps(az_v, scale_v);

    __m128 x_vel_sum_v = _mm_add_ps(vx_v, calculated_vel_x);
    __m128 y_vel_sum_v = _mm_add_ps(vy_v, calculated_vel_y);
    __m128 z_vel_sum_v = _mm_add_ps(vz_v, calculated_vel_z);

    _mm_store_ps(vx+i, x_vel_sum_v); 
    _mm_store_ps(vy+i, y_vel_sum_v); 
    _mm_store_ps(vz+i, z_vel_sum_v); 

    Vel(p)[0] = vx[i];
    Vel(p)[1] = vy[i];
    Vel(p)[2] = vz[i];

    Vel(p+1)[0] = vx[i+1];
    Vel(p+1)[1] = vy[i+1];
    Vel(p+1)[2] = vz[i+1];

    Vel(p+2)[0] = vx[i+2];
    Vel(p+2)[1] = vy[i+2];
    Vel(p+2)[2] = vz[i+2];

    Vel(p+3)[0] = vx[i+3];
    Vel(p+3)[1] = vy[i+3];
    Vel(p+3)[2] = vz[i+3];

    // (Vel(p))[0] += (Acc(p))[0] * scale;           
    // (Vel(p))[1] += (Acc(p))[1] * scale;               
    // (Vel(p))[2] += (Acc(p))[2] * scale;  

    __m128 calculated_x = _mm_mul_ps(vx_v, dtime_v);
    __m128 calculated_y = _mm_mul_ps(vx_v, dtime_v);
    __m128 calculated_z = _mm_mul_ps(vx_v, dtime_v);

    __m128 x_pos_sum_v = _mm_add_ps(x_pos_v, calculated_x);
    __m128 y_pos_sum_v = _mm_add_ps(y_pos_v, calculated_y);
    __m128 z_pos_sum_v = _mm_add_ps(z_pos_v, calculated_z);

    _mm_store_ps(x_pos+i, x_pos_sum_v); 
    _mm_store_ps(y_pos+i, y_pos_sum_v); 
    _mm_store_ps(z_pos+i, z_pos_sum_v); 

    Pos(p)[0] = x_pos[i];
    Pos(p)[1] = y_pos[i];
    Pos(p)[2] = z_pos[i];

    Pos(p+1)[0] = x_pos[i+1];
    Pos(p+1)[1] = y_pos[i+1];
    Pos(p+1)[2] = z_pos[i+1];

    Pos(p+2)[0] = x_pos[i+2];
    Pos(p+2)[1] = y_pos[i+2];
    Pos(p+2)[2] = z_pos[i+2];

    Pos(p+3)[0] = x_pos[i+3];
    Pos(p+3)[1] = y_pos[i+3];
    Pos(p+3)[2] = z_pos[i+3];

    // (Pos(p))[0] += (Vel(p))[0] * (dtime);           
    // (Pos(p))[1] += (Vel(p))[1] * (dtime);               
    // (Pos(p))[2] += (Vel(p))[2] * (dtime);        

    p += 4;

    //ADDMULVS(Vel(p), Acc(p), 0.5 * dtime);	// advance v by 1/2 step
    //ADDMULVS(Pos(p), Vel(p), dtime);		// advance r by 1 step
  }
  for (; p < p2; p++) {
    (Vel(p))[0] += (Acc(p))[0] * (0.5 * dtime);           
    (Vel(p))[1] += (Acc(p))[1] * (0.5 * dtime);               
    (Vel(p))[2] += (Acc(p))[2] * (0.5 * dtime);    

    (Pos(p))[0] += (Vel(p))[0] * (dtime);           
    (Pos(p))[1] += (Vel(p))[1] * (dtime);               
    (Pos(p))[2] += (Vel(p))[2] * (dtime);      
  }

  treeforce();

  p = p1;

  for (int i = 0; i < vector_size; i+=4) {
    //ADDMULVS(Vel(p), Acc(p), 0.5 * dtime);
    __m128 vx_v = _mm_load_ps(vx+i); 
    __m128 vy_v = _mm_load_ps(vy+i);
    __m128 vz_v = _mm_load_ps(vz+i);

    __m128 ax_v = _mm_load_ps(ax+i); 
    __m128 ay_v = _mm_load_ps(ay+i);
    __m128 az_v = _mm_load_ps(az+i);

    __m128 calculated_v_x = _mm_mul_ps(ax_v, scale_v);
    __m128 calculated_v_y = _mm_mul_ps(ay_v, scale_v);
    __m128 calculated_v_z = _mm_mul_ps(az_v, scale_v);

    __m128 x_v_sum_v = _mm_add_ps(vx_v, calculated_v_x);
    __m128 y_v_sum_v = _mm_add_ps(vy_v, calculated_v_y);
    __m128 z_v_sum_v = _mm_add_ps(vz_v, calculated_v_z);

    _mm_store_ps(vx+i, x_v_sum_v); 
    _mm_store_ps(vy+i, y_v_sum_v); 
    _mm_store_ps(vz+i, z_v_sum_v);

    Vel(p)[0] = vx[i];
    Vel(p)[1] = vy[i];
    Vel(p)[2] = vz[i];

    Vel(p+1)[0] = vx[i+1];
    Vel(p+1)[1] = vy[i+1];
    Vel(p+1)[2] = vz[i+1];

    Vel(p+2)[0] = vx[i+2];
    Vel(p+2)[1] = vy[i+2];
    Vel(p+2)[2] = vz[i+2];

    Vel(p+3)[0] = vx[i+3];
    Vel(p+3)[1] = vy[i+3];
    Vel(p+3)[2] = vz[i+3];

    p += 4;
  }  
  for (; p < p2; p++) {			// loop over body range	
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
  for (p = bodytab; p < bodytab+nbody; p++)	{// loop over all bodies
    Update(p) = (testcalc ? p1 <= p && p < p2 : TRUE);
  }
						// flag bodies to update
  maketree(bodytab, nbody);			// construct tree structure
  gravcalc();					// compute current forces
  // Add accelerations to array.
  int i = 0;
  for (p = bodytab; p < bodytab+nbody; p++) {// loop over all bodies
    ax[i] = Acc(p)[0];
    ay[i] = Acc(p)[1];
    az[i] = Acc(p)[2];
    i++;
  }
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
