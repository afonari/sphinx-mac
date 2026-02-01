// ---------------------------------------------------------------------------
//
//           The general purpose cross platform C/C++ framework
//
//                       S x A c c e l e r a t e
//
//           Home:       https://www.sxlib.de
//           License:    Apache 2
//           Authors:    see src/AUTHORS
//
// ---------------------------------------------------------------------------

#include <SxVector.h>
#include <string.h>
#ifdef USE_GPU
#   include <cuda_runtime_api.h>
#endif

SxAllocMem* SxAllocMem::create (ssize_t n, int memTypeIn)
{
   SX_CHECK (n > 0, n);
   // allocate handle (separate cacheline)
   void *handle = SxAllocation::getAligned (roundCL(sizeof(SxAllocMem)),
                                            roundCL(1));
   //SX_CHECK (handle);
   // construct SxAllocMem on this memory
   return new (handle) SxAllocMem((size_t)n, memTypeIn);
}

inline
SxAllocMem::SxAllocMem (size_t n, int memTypeIn)
   : allocSize (n), memType(memTypeIn), refCounter (1)
{
   if (memType == NormalRAM)  {
      // --- memtype = 0  => normal RAM
      mem = SxAllocation::get (n);
#     ifndef NDEBUG
         // fill with nan
         ::memset (mem, 0xff, allocSize);
#     endif
#ifdef USE_GPU
   } else if (memType & GPUMemory)  {
      // --- switch to desired GPU device
      bool currentGPU = memType & CurrentGPUMemory;
      int currentDevice = -1;
      int cudaStatus = cudaGetDevice (&currentDevice);
      SX_CHECK(cudaStatus == cudaSuccess, cudaStatus);
      int gpuId = currentGPU ? currentDevice : getLocationId ();
      if (currentGPU)  {
         const_cast<int&>(memType) += currentDevice << 8;
      } else if (currentDevice != gpuId)  {
#ifndef NDEBUG
         int nGPU = -1;
         cudaStatus = cudaGetDeviceCount (&nGPU);
         SX_CHECK(cudaStatus == cudaSuccess, cudaStatus);
         SX_CHECK (gpuId >= 0 && gpuId < nGPU, gpuId, nGPU);
#endif
         cudaStatus = cudaSetDevice (gpuId);
         SX_CHECK(cudaStatus == cudaSuccess, cudaStatus);
      }
      // allocate
      cudaStatus = cudaMalloc (&mem, n);
      SX_CHECK(cudaStatus == cudaSuccess, cudaStatus);
      // --- back to original GPU device
      if (currentDevice != gpuId)  {
         cudaStatus = cudaSetDevice (currentDevice);
         SX_CHECK(cudaStatus == cudaSuccess, cudaStatus);
      }
#endif
   } else {
      // not supported
      SX_EXIT;
   }
}

SxAllocMem::~SxAllocMem ()
{
   if (memType == NormalRAM)  {
      if (mem) free (mem);
#ifdef USE_GPU
   } else if (memType & GPUMemory) {
      bool currentGPU = memType & CurrentGPUMemory;
      int currentDevice = -1;
      int cudaStatus = cudaGetDevice (&currentDevice);
      SX_CHECK(cudaStatus == cudaSuccess, cudaStatus);
      int gpuId = getLocationId ();
      if (currentGPU)  {
         SX_CHECK (gpuId == currentDevice);
      } else if (currentDevice != gpuId)  {
         cudaStatus = cudaSetDevice (gpuId);
         SX_CHECK(cudaStatus == cudaSuccess, cudaStatus);
      }
      // deallocate
      cudaStatus = cudaFree (mem);
      SX_CHECK(cudaStatus == cudaSuccess, cudaStatus);
      // --- back to original GPU device
      if (currentDevice != gpuId)  {
         cudaStatus = cudaSetDevice (currentDevice);
         SX_CHECK(cudaStatus == cudaSuccess, cudaStatus);
      }
#endif
   } else {
      // not supported
      SX_EXIT;
   }
}
