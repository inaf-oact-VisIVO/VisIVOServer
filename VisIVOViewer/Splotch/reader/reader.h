#ifndef READER_H
#define READER_H

#include "splotch/splotchutils.h"

void gadget_reader(paramfile &params, int interpol_mode,
  std::vector<particle_sim> &p, std::vector<MyIDType> &id,
  std::vector<vec3f> &vel, int snr, double &time, double &boxsize);
void gadget_hdf5_reader(paramfile &params, int interpol_mode,
  std::vector<particle_sim> &p, std::vector<MyIDType> &id,
  std::vector<vec3f> &vel, int snr, double &time, double &redshift, double &boxsize);
void gadget_millenium_reader(paramfile &params, std::vector<particle_sim> &p, int snr, double *time);
void bin_reader_tab (paramfile &params, std::vector<particle_sim> &points);
void bin_reader_block (paramfile &params, std::vector<particle_sim> &points);

#ifdef SPLVISIVO
bool visivo_reader(paramfile &params, std::vector<particle_sim> &points,VisIVOServerOptions &opt);
#else
void visivo_reader();
#endif

long bin_reader_block_mpi (paramfile &params, std::vector<particle_sim> &points, float *maxr, float *minr, int mype, int npes);
void mesh_reader(paramfile &params, std::vector<particle_sim> &points);
void hdf5_reader(paramfile &params, std::vector<particle_sim> &points);
void tipsy_reader(paramfile &params, std::vector<particle_sim> &points);
void galaxy_reader(paramfile &params, std::vector<particle_sim> &points);
void h5part_reader(paramfile &params, std::vector<particle_sim> &points);
void ramses_reader(paramfile &params, std::vector<particle_sim> &points);
#endif
