#ifndef CURENDER_H
#define CURENDER_H

#include "cuda/splotch_cuda2.h"
#include "cuda/splotch_cuda.h"


using namespace std;

int cu_draw_chunk(int mydevID, cu_particle_sim *d_particle_data, int nParticle, arr2<COLOUR> &Pic_host, cu_gpu_vars* gv, bool a_eq_e, float64 grayabsorb, int xres, int yres, bool doLogs);
int add_device_image(arr2<COLOUR> &Pic_host, cu_gpu_vars* gv, int xres, int yres);

#endif
