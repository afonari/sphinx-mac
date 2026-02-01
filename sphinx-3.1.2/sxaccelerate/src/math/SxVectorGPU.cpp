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

#include <SxError.h>
#include <SxVectorGPU.h>
#include <SxVector.h>

/* --- This file defines functions needed to work with vectors living
       on the GPU.
       It needs cuda.
*/
#ifdef USE_GPU
#include <cuda_runtime_api.h>

#define CU_EXEC( CU_CALL )                    \
do {                                          \
   cudaError_t status = CU_CALL;              \
   if (status != cudaSuccess)                 \
      sxCudaCrash (__FILE__, __LINE__, status); \
} while (0)

void sxCudaCrash (const char* file, int line, cudaError_t status)
{
   printf("CUDA error %s at %s:%d\n'%s'\n", cudaGetErrorName(status),
          file, line, cudaGetErrorString (status));
   SX_EXIT;
}

#ifndef NDEBUG
inline
#endif
cudaMemcpyKind getCopyKind (int srcMemType, int dstMemType)
{
#ifndef NDEBUG
   int currentDevice = -1;
   int cudaStatus = cudaGetDevice (&currentDevice);
#endif
   SX_CHECK (cudaStatus == cudaSuccess, cudaStatus);
   SX_CHECK (SxAllocMem::isGPU (dstMemType) || SxAllocMem::isGPU (srcMemType));
   SX_CHECK (   SxAllocMem::isRAM (srcMemType)
             || SxAllocMem::getLocationId (srcMemType) == currentDevice,
             SxAllocMem::getLocationId (srcMemType), currentDevice);
   SX_CHECK (   SxAllocMem::isRAM (dstMemType)
             || SxAllocMem::getLocationId (dstMemType) == currentDevice,
             SxAllocMem::getLocationId (dstMemType), currentDevice);
   return    SxAllocMem::isRAM (dstMemType) ? cudaMemcpyDeviceToHost
          : (SxAllocMem::isRAM (srcMemType) ? cudaMemcpyHostToDevice
                                            : cudaMemcpyDeviceToDevice  );
}

void memCopyGPU (      void *dst, ssize_t dstColStride, int dstMemType,
                 const void *src, ssize_t srcColStride, int srcMemType,
                 ssize_t bytesPerCol, ssize_t nCol)
{

   CU_EXEC (cudaMemcpy2D (dst, dstColStride, src, srcColStride,
                          bytesPerCol, nCol,
                          getCopyKind (srcMemType, dstMemType)   ) );
}

void memCopyGPU (      void *dst, int dstMemType,
                 const void *src, int srcMemType, ssize_t nBytes)
{
   CU_EXEC (cudaMemcpy (dst, src, nBytes, getCopyKind (srcMemType, dstMemType)));
}

#endif
