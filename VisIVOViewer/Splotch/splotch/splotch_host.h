#ifndef SPLOTCH_HOST_H
#define SPLOTCH_HOST_H

#include "splotch/splotchutils.h"

#ifdef SPLVISIVO
void host_rendering (paramfile &params, std::vector<particle_sim> &particle_data,
  arr2<COLOUR> &pic, const vec3 &campos, const vec3 &lookat, const vec3 &sky,
  std::vector<COLOURMAP> &amap, float b_brightness, VisIVOServerOptions &opt);
#else
void host_rendering (paramfile &params, std::vector<particle_sim> &particle_data,
  arr2<COLOUR> &pic, const vec3 &campos, const vec3 &lookat, const vec3 &sky,
  std::vector<COLOURMAP> &amap, float b_brightness);
#endif

#endif
