find_package(CUDA REQUIRED)
SET(PRJ_SRC
         splotch.cu
         splotch_cuda2.cc
         CuPolicy.cpp
         deviceQuery.cu
        
)

# splotch_kernel.cu
# splotch_cuda.h
# splotch_cuda2.h
# CuPolicy.h