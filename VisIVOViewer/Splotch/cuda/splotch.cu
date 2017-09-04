/*
Try accelerating splotch with CUDA. July 2009.
Copyright things go here.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <cuda.h>
#include <cutil_inline.h>

#include "cxxsupport/lsconstants.h"
#include "cxxsupport/string_utils.h"
#include "splotch/splotchutils.h"
#include "kernel/transform.h"

#include "cuda/splotch_kernel.cu"
#include "cuda/splotch_cuda.h"
#include "cuda/CuPolicy.h"


using namespace std;

template<typename T> T findParamWithoutChange
  (paramfile *param, std::string &key, T &deflt)
  {
  return param->param_present(key) ? param->find<T>(key) : deflt;
  }

#define CLEAR_MEM(p) if(p) {cutilSafeCall(cudaFree(p)); p=0;}


void getCuTransformParams(cu_param &para_trans,
paramfile &params, vec3 &campos, vec3 &lookat, vec3 &sky)
  {
  int xres = params.find<int>("xres",800),
      yres = params.find<int>("yres",xres);
  double fov = params.find<double>("fov",45); //in degrees
  double fovfct = tan(fov*0.5*degr2rad);
  float64 xfac=0.0, dist=0.0;

  sky.Normalize();
  vec3 zaxis = (lookat-campos).Norm();
  vec3 xaxis = crossprod (sky,zaxis).Norm();
  vec3 yaxis = crossprod (zaxis,xaxis);
  TRANSFORM trans;
  trans.Make_General_Transform
        (TRANSMAT(xaxis.x,xaxis.y,xaxis.z,
                  yaxis.x,yaxis.y,yaxis.z,
                  zaxis.x,zaxis.y,zaxis.z,
                  0,0,0));
  trans.Invert();
  TRANSFORM trans2;
  trans2.Make_Translation_Transform(-campos);
  trans2.Add_Transform(trans);
  trans=trans2;
  bool projection = params.find<bool>("projection",true);

  if (!projection)
    {
    float64 dist= (campos-lookat).Length();
    float64 xfac=1./(fovfct*dist);
    cout << " Field of fiew: " << 1./xfac*2. << endl;
    }

  float minrad_pix = params.find<float>("minrad_pix",1.);

  //retrieve the parameters for transformation
  for (int i=0; i<12; i++)
    para_trans.p[i] =trans.Matrix().p[i];
  para_trans.projection=projection;
  para_trans.xres=xres;
  para_trans.yres=yres;
  para_trans.fovfct=fovfct;
  para_trans.dist=dist;
  para_trans.xfac=xfac;
  para_trans.minrad_pix=minrad_pix;
  }


void cu_init(int devID, int nP, cu_gpu_vars* pgv, paramfile &fparams, vec3 &campos, vec3 &lookat, vec3 &sky)
  {
  cudaSetDevice (devID); // initialize cuda runtime
  
  //allocate device memory for particle data
  size_t s = pgv->policy->GetSizeDPD(nP);
  //one more space allocated for the dumb
  cutilSafeCall(cudaMalloc((void**) &pgv->d_pd, s +sizeof(cu_particle_sim)));
  
  //now prepare memory for d_particle_splotch.
  //one more for dums
  s = nP* sizeof(cu_particle_splotch);
  cutilSafeCall( cudaMalloc((void**) &pgv->d_ps_render, s+sizeof(cu_particle_splotch)));

  size_t size = pgv->policy->GetFBufSize() <<20;
  cutilSafeCall( cudaMalloc((void**) &pgv->d_fbuf, size)); 

  //retrieve parameters
  cu_param tparams;
  getCuTransformParams(tparams,fparams,campos,lookat,sky);

  tparams.zmaxval   = fparams.find<float>("zmax",1.e23);
  tparams.zminval   = fparams.find<float>("zmin",0.0);
  tparams.ptypes    = fparams.find<int>("ptypes",1);

  for(int itype=0; itype<tparams.ptypes; itype++)
    {
    tparams.brightness[itype] = fparams.find<double>("brightness"+dataToString(itype),1.);
    tparams.col_vector[itype] = fparams.find<bool>("color_is_vector"+dataToString(itype),false);
    }
  tparams.rfac=1.5;

  //dump parameters to device
  cutilSafeCall( cudaMemcpyToSymbol(dparams, &tparams, sizeof(cu_param) ));
  }

void cu_allocate_particles(unsigned int nP, cu_gpu_vars* pgv)
  {
  //now resize d_particle_splotch.
  //one more for dums
  size_t s = (nP+1)* sizeof(cu_particle_splotch);
  cutilSafeCall( cudaMalloc((void**) &pgv->d_ps_render, s)); 
  }

void cu_copy_particles_to_device(cu_particle_sim* h_pd, unsigned int n, cu_gpu_vars* pgv)
  {
  //copy particle data to device
  size_t s = pgv->policy->GetSizeDPD(n);
  cutilSafeCall(cudaMemcpy(pgv->d_pd, h_pd, s, cudaMemcpyHostToDevice) );
  }


void cu_transform (unsigned int n, cu_particle_splotch *h_ps, cu_gpu_vars* pgv)
  {

  //Get block dim and grid dim from pgv->policy object
  dim3 dimGrid, dimBlock;
  pgv->policy->GetDimsBlockGrid(n, &dimGrid, &dimBlock);

  //call device transformation
  k_transform<<<dimGrid,dimBlock>>>(pgv->d_pd, pgv->d_ps_render, n);

  //copy the result out
  size_t size = n* sizeof(cu_particle_splotch);
  cutilSafeCall(cudaMemcpy(h_ps, pgv->d_ps_render, size, cudaMemcpyDeviceToHost) );

  }

void cu_init_colormap(cu_colormap_info h_info, cu_gpu_vars* pgv)
  {
  //allocate memories for colormap and ptype_points and dump host data into it
  size_t size =sizeof(cu_color_map_entry)*h_info.mapSize;
  cutilSafeCall(cudaMemcpyToSymbol(dmap, h_info.map, size) );
  //type
  size =sizeof(int)*h_info.ptypes;
  cutilSafeCall(cudaMemcpyToSymbol(ptype_points, h_info.ptype_points, size) );

  //set fields of global variable pgv->d_colormap_info
  pgv->colormap_size   = h_info.mapSize;
  pgv->colormap_ptypes = h_info.ptypes;
  }

void cu_colorize(int n, cu_gpu_vars* pgv)
  {

  //fetch grid dim and block dim and call device
  dim3 dimGrid, dimBlock;
  pgv->policy->GetDimsBlockGrid(n, &dimGrid, &dimBlock);

  cudaEvent_t start,stop;
  cutilSafeCall( cudaEventCreate(&start));
  cutilSafeCall( cudaEventCreate(&stop));
  cutilSafeCall( cudaEventRecord( start, 0));

  k_colorize<<<dimGrid,dimBlock>>>(n, pgv->d_ps_render, pgv->colormap_size, pgv->colormap_ptypes);

  cutilSafeCall( cudaEventRecord( stop, 0));
  cutilSafeCall( cudaEventSynchronize(stop));
 // float elapsedTime;
 // cutilSafeCall( cudaEventElapsedTime(&elapsedTime,start,stop));
  cutilSafeCall( cudaEventDestroy(start));
  cutilSafeCall( cudaEventDestroy(stop));

  //particle_splotch memory on device will be freed in cu_end
  }


void cu_copy_particles_to_render(cu_particle_splotch *p,
  int n, cu_gpu_vars* pgv)
  {
  //copy filtered particles into device
  size_t size = n *sizeof(cu_particle_splotch);
  cutilSafeCall(cudaMemcpy(pgv->d_ps_render, p,size,
    cudaMemcpyHostToDevice) );
  }

void cu_render1
  (int nP, bool a_eq_e, float grayabsorb, cu_gpu_vars* pgv)
  {
  cudaEvent_t start,stop;
  cutilSafeCall( cudaEventCreate(&start));
  cutilSafeCall( cudaEventCreate(&stop));
  cutilSafeCall( cudaEventRecord( start, 0));
 
  //get dims from pgv->policy object first
  dim3 dimGrid, dimBlock;
  pgv->policy->GetDimsBlockGrid(nP, &dimGrid, &dimBlock);

  //call device
  k_render1<<<dimGrid, dimBlock>>>(pgv->d_ps_render, nP,
    pgv->d_fbuf, a_eq_e, grayabsorb, pgv->colormap_size, pgv->colormap_ptypes);

  cutilSafeCall( cudaEventRecord( stop, 0));
  cutilSafeCall( cudaEventSynchronize(stop));
//  float elapsedTime;
//  cutilSafeCall( cudaEventElapsedTime(&elapsedTime,start,stop));
  cutilSafeCall( cudaEventDestroy(start));
  cutilSafeCall( cudaEventDestroy(stop));
  }


void cu_get_fbuf
  (void *h_fbuf, bool a_eq_e, unsigned long n, cu_gpu_vars* pgv)
  {
  size_t size;
  if (a_eq_e)
    size =n* sizeof(cu_fragment_AeqE);
  else
    size =n* sizeof(cu_fragment_AneqE);

  cutilSafeCall( cudaMemcpy(h_fbuf, pgv->d_fbuf,size,
    cudaMemcpyDeviceToHost)) ;
  }

void cu_end(cu_gpu_vars* pgv)
  {
  CLEAR_MEM((pgv->d_pd));
  CLEAR_MEM((pgv->d_ps_render));
  CLEAR_MEM((pgv->d_fbuf));

  cudaThreadExit();

  delete pgv->policy;
  }

int cu_get_chunk_particle_count(paramfile &params, CuPolicy* policy, size_t psize, float pfactor)
  {
   int gMemSize = policy->GetGMemSize();
   int fBufSize = policy->GetFBufSize();
   if (gMemSize <= fBufSize) return 0;

  // float factor =params.find<float>("particle_mem_factor", 3);
   int spareMem = 10;
   int arrayParticleSize = gMemSize - fBufSize - spareMem;

   return (int) (arrayParticleSize/psize/pfactor)*(1<<20);
  }
