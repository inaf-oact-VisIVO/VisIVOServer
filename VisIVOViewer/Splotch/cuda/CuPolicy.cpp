

#include "CuPolicy.h"

CuPolicy::CuPolicy(paramfile &Param)
  {
    cudaDeviceProp deviceProp;
    CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, 0));
    m_blockSize = deviceProp.maxThreadsPerBlock;
    gmsize = deviceProp.totalGlobalMem;
    maxregion = Param.find<int>("max_region", 1024);
//    fbsize = Param.find<int>("fragment_buffer_size", 128);
    fbsize = GetGMemSize()/4;
  }

//Get size of device particle data
size_t CuPolicy::GetSizeDPD(int n)
  { return n* sizeof(cu_particle_sim); }

int CuPolicy::GetMaxRegion()
  { return maxregion; }

int CuPolicy::GetFBufSize() // return dimension in terms of Megabytes
  {
     return fbsize; 
  }

int CuPolicy::GetGMemSize() // return dimension in terms of Megabytes
  { 
    int M = 1<<20;
    int size = gmsize/M;
    return size; 
  }

void CuPolicy::GetDimsBlockGrid(int n, dim3 *dimGrid, dim3 *dimBlock)
  {
    *dimBlock = dim3(m_blockSize);
    int nBlock = (n+m_blockSize-1)/m_blockSize;
    *dimGrid =dim3(nBlock); 
  }
