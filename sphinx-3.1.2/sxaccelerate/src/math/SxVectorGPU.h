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

#ifndef _SX_VECTOR_GPU_H_
#define _SX_VECTOR_GPU_H_

#include <SxMath.h>

/// \brief Copy matrix to/from/within GPU
void memCopyGPU (      void *dst, ssize_t dstColStride, int dstMemType,
                 const void *src, ssize_t srcColStride, int srcMemType,
                 ssize_t bytesPerCol, ssize_t nCol);

/// \brief Copy vector to/from/within GPU
void memCopyGPU (      void *dst, int dstMemType,
                 const void *src, int srcMemType, ssize_t nBytes);

#endif /* _SX_VECTOR_GPU_H_ */
