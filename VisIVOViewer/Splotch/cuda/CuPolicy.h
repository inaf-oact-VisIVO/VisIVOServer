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

#ifdef __CUDACC__
#include <cuda.h>
#else
struct dim3;
#endif

using namespace std;

class CuPolicy
  {
  private:
    int m_gridSize, p_blockSize;
    pair <int,int> res, tile_size;
    int boundary_width, x_num_tiles, y_num_tiles;
    size_t gmsize;
  public:
    CuPolicy(int xres, int yres, paramfile &params);

    void GetTileInfo(int *tile_sidex, int *tiley, int *width, int *nxtiles, int *nytiles);
    int GetNumTiles();
    size_t GetGMemSize();
    size_t GetImageSize();
    int GetBlockSize();
    int GetMaxGridSize();
    void GetDimsBlockGrid(int n, dim3 *dimGrid, dim3 *dimBlock);
  };

#endif //CUPOLICY_H
