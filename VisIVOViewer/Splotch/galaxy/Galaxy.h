# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <string>
# include <iostream>
# include "cxxsupport/paramfile.h"
#include "cxxsupport/bstream.h"
#ifdef HDF5
# include "hdf5.h"
#endif

using namespace std;

long RDiscFunc (paramfile &params, string ComponentName, long number_of_points, long tot, float * coordx,
                 float * coordy, float * coordz, float * II, long nnx, long nny);

long RDiscFuncTirific (paramfile &params, string ComponentName, long number_of_points, long tot, float * coordx,
                 float * coordy, float * coordz, float * II, long nnx, long nny);

long GaussRFunc (paramfile &params, string ComponentName, long number_of_points, float * coordx,
                 float * coordy, float * coordz);

long GaussRDiscFunc (paramfile &params, string ComponentName, long number_of_points, long tot, float * coordx, float * coordy, float * coordz, long nnx, long nny);

long GaussRGlobFunc (paramfile &params, string ComponentName, long number_of_points, long ntot, float * coordx, float * coordy, float * coordz, float * III, long nx, long ny);

void CalculateDensity (float * hsml, float * rho, float * xcoord, float * ycoord,
                       float * zcoord, long numofpart, float smooth);

void CalculateColours (paramfile &params, string ComponentName, long npart, 
                       float * cred, float * cgreen, float * cblue, float * ciii, 
                       float * Red, float * Green, float * Blue, float * III, float * xcoord, 
                       float * ycoord, long nxxx, long nyyy);

long GlobularCluster (paramfile &params, string ComponentName, long number_of_points, long ntot,
                      float * coordx, float * coordy, float * coordz);

long ReadImages (paramfile &params, string infile_rgb, string infile_mask, long numx, long numy, float * RRR,
                 float * GGG, float * BBB, float * III, float * xx, float * yy, long nwant);

const int NUM_OF_FIELDS = 11;

