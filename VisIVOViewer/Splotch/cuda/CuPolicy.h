/*
Copyright things go here.
*/
/*
CuPolicy is the class that knows the overall state of cuda application.
All 'magic numbers' are out of this class.
*/
#ifndef CUPOLICY_H
#define CUPOLICY_H

#include "cxxsupport/paramfile.h"
#include "cuda/splotch_cuda.h"
#include <cuda.h>
#include <cutil_inline.h>

class CuPolicy
  {
  private:
    int m_blockSize, maxregion, fbsize;
    size_t gmsize;
  public:
    CuPolicy(paramfile &Param);

    int GetMaxRegion();
    int GetFBufSize();
    int GetGMemSize();
    size_t GetSizeDPD(int n);
    void GetDimsBlockGrid(int n, dim3 *dimGrid, dim3 *dimBlock);
  };

#endif //CUPOLICY_H
