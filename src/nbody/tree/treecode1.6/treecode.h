/*
 * treecode.h: define various things for treecode.c and treeio.c.
 * Copyright (c) 2015 by Joshua E. Barnes, Kyoto, Japan.
 */

#ifndef _treecode_h
#define _treecode_h

#include "stdinc.h"
#include "mathfns.h"
#include "vectmath.h"
#include "getparam.h"
#include "datatypes.h"
#include "treedefs.h"
#include "fixbody.h"

//  Parameters, state variables, and diagnostics for N-body integration.
//  ____________________________________________________________________

global string infile;			// file name for snapshot input
global string outfile;			// file pattern for snapshot output
global string logfile;			// file name for log file
global string savefile;			// file pattern for state output
global string strfile;			// file pattern for frame stream
global real dtime;			// basic integration timestep
global real dtout;			// data output timestep
global real tstop;			// time to stop calculation
global string outputs;			// list of data to output
global string headline;			// message describing calculation
global real tnow;			// current value of time
global real tout;			// time of next output
global int nstep;			// number of time-steps
global int nbody;			// number of bodies in system
global int nstatic;			// number of static bodies
global bodyptr bodytab;			// points to array of bodies
global bool testcalc;			// indicate test-particle calc.


// my stuff
global double wtime();
global double startrun_time_0;
global double startrun_time_1;
global double treeforce_initial_0;
global double treeforce_initial_1;
global double treeforce_0;
global double treeforce_1;

global int timesteps;     // number of timesteps

//  Prototypes for I/O routines.
//  ____________________________

void inputdata(void);			// read initial data file
void startoutput(void);			// open files for output
void forcereport(void);			// report on force calculation
void output(void);			// perform output operation
void savestate(string);			// save system state		
void restorestate(string);		// restore system state
void finaloutput();

#endif // ! _treecode_h
