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

#include <SxMath.h>
#include <SxComplex.h>
#ifdef USE_GPU
#include <SxComputeGPU.h>

// --- The next 980 lines were generated from snippets/SxComputeGPU.cpp snippet BINOP
/** \brief Compute res = x + y (contiguous vector)
  @param resPtr target vector
  @param n      vector size
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param gpuId  GPU id, if -1: use the current one
*/
void addGPU (float* resPtr, ssize_t n,
      const float* xPtr, const float* yPtr, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x + y (contiguous columns)
  @param resPtr target vector
  @param nRows        number of rows (contiguous in memory)
  @param nCols        number of columns (with stride >= nRows)
  @param xPtr         1st operand
  @param xColStride   column stride for x
  @param yPtr         2nd operand
  @param yColStride   column stride for y
  @param resPtr       result
  @param resColStride column stride for res
  @param gpuId  GPU id, if -1: use the current one
*/
void addGPU (ssize_t nRows, ssize_t nCols,
                  const float* xPtr, ssize_t xColStride,
                  const float* yPtr, ssize_t yColStride,
                  float* resPtr, ssize_t resColStride,
                  int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute x += y (contiguous vector)
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param n      vector size
  @param gpuId  GPU id, if -1: use the current one
*/
void addInPlaceGPU (float* xPtr, const float* yPtr,
      ssize_t n, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x - y (contiguous vector)
  @param resPtr target vector
  @param n      vector size
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param gpuId  GPU id, if -1: use the current one
*/
void subtractGPU (float* resPtr, ssize_t n,
      const float* xPtr, const float* yPtr, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x - y (contiguous columns)
  @param resPtr target vector
  @param nRows        number of rows (contiguous in memory)
  @param nCols        number of columns (with stride >= nRows)
  @param xPtr         1st operand
  @param xColStride   column stride for x
  @param yPtr         2nd operand
  @param yColStride   column stride for y
  @param resPtr       result
  @param resColStride column stride for res
  @param gpuId  GPU id, if -1: use the current one
*/
void subtractGPU (ssize_t nRows, ssize_t nCols,
                  const float* xPtr, ssize_t xColStride,
                  const float* yPtr, ssize_t yColStride,
                  float* resPtr, ssize_t resColStride,
                  int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute x -= y (contiguous vector)
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param n      vector size
  @param gpuId  GPU id, if -1: use the current one
*/
void subtractInPlaceGPU (float* xPtr, const float* yPtr,
      ssize_t n, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x * y (contiguous vector)
  @param resPtr target vector
  @param n      vector size
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param gpuId  GPU id, if -1: use the current one
*/
void multiplyGPU (float* resPtr, ssize_t n,
      const float* xPtr, const float* yPtr, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x * y (contiguous columns)
  @param resPtr target vector
  @param nRows        number of rows (contiguous in memory)
  @param nCols        number of columns (with stride >= nRows)
  @param xPtr         1st operand
  @param xColStride   column stride for x
  @param yPtr         2nd operand
  @param yColStride   column stride for y
  @param resPtr       result
  @param resColStride column stride for res
  @param gpuId  GPU id, if -1: use the current one
*/
void multiplyGPU (ssize_t nRows, ssize_t nCols,
                  const float* xPtr, ssize_t xColStride,
                  const float* yPtr, ssize_t yColStride,
                  float* resPtr, ssize_t resColStride,
                  int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute x *= y (contiguous vector)
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param n      vector size
  @param gpuId  GPU id, if -1: use the current one
*/
void multiplyInPlaceGPU (float* xPtr, const float* yPtr,
      ssize_t n, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x / y (contiguous vector)
  @param resPtr target vector
  @param n      vector size
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param gpuId  GPU id, if -1: use the current one
*/
void divideGPU (float* resPtr, ssize_t n,
      const float* xPtr, const float* yPtr, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x / y (contiguous columns)
  @param resPtr target vector
  @param nRows        number of rows (contiguous in memory)
  @param nCols        number of columns (with stride >= nRows)
  @param xPtr         1st operand
  @param xColStride   column stride for x
  @param yPtr         2nd operand
  @param yColStride   column stride for y
  @param resPtr       result
  @param resColStride column stride for res
  @param gpuId  GPU id, if -1: use the current one
*/
void divideGPU (ssize_t nRows, ssize_t nCols,
                  const float* xPtr, ssize_t xColStride,
                  const float* yPtr, ssize_t yColStride,
                  float* resPtr, ssize_t resColStride,
                  int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute x /= y (contiguous vector)
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param n      vector size
  @param gpuId  GPU id, if -1: use the current one
*/
void divideInPlaceGPU (float* xPtr, const float* yPtr,
      ssize_t n, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x + y (contiguous vector)
  @param resPtr target vector
  @param n      vector size
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param gpuId  GPU id, if -1: use the current one
*/
void addGPU (double* resPtr, ssize_t n,
      const double* xPtr, const double* yPtr, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x + y (contiguous columns)
  @param resPtr target vector
  @param nRows        number of rows (contiguous in memory)
  @param nCols        number of columns (with stride >= nRows)
  @param xPtr         1st operand
  @param xColStride   column stride for x
  @param yPtr         2nd operand
  @param yColStride   column stride for y
  @param resPtr       result
  @param resColStride column stride for res
  @param gpuId  GPU id, if -1: use the current one
*/
void addGPU (ssize_t nRows, ssize_t nCols,
                  const double* xPtr, ssize_t xColStride,
                  const double* yPtr, ssize_t yColStride,
                  double* resPtr, ssize_t resColStride,
                  int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute x += y (contiguous vector)
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param n      vector size
  @param gpuId  GPU id, if -1: use the current one
*/
void addInPlaceGPU (double* xPtr, const double* yPtr,
      ssize_t n, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x - y (contiguous vector)
  @param resPtr target vector
  @param n      vector size
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param gpuId  GPU id, if -1: use the current one
*/
void subtractGPU (double* resPtr, ssize_t n,
      const double* xPtr, const double* yPtr, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x - y (contiguous columns)
  @param resPtr target vector
  @param nRows        number of rows (contiguous in memory)
  @param nCols        number of columns (with stride >= nRows)
  @param xPtr         1st operand
  @param xColStride   column stride for x
  @param yPtr         2nd operand
  @param yColStride   column stride for y
  @param resPtr       result
  @param resColStride column stride for res
  @param gpuId  GPU id, if -1: use the current one
*/
void subtractGPU (ssize_t nRows, ssize_t nCols,
                  const double* xPtr, ssize_t xColStride,
                  const double* yPtr, ssize_t yColStride,
                  double* resPtr, ssize_t resColStride,
                  int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute x -= y (contiguous vector)
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param n      vector size
  @param gpuId  GPU id, if -1: use the current one
*/
void subtractInPlaceGPU (double* xPtr, const double* yPtr,
      ssize_t n, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x * y (contiguous vector)
  @param resPtr target vector
  @param n      vector size
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param gpuId  GPU id, if -1: use the current one
*/
void multiplyGPU (double* resPtr, ssize_t n,
      const double* xPtr, const double* yPtr, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x * y (contiguous columns)
  @param resPtr target vector
  @param nRows        number of rows (contiguous in memory)
  @param nCols        number of columns (with stride >= nRows)
  @param xPtr         1st operand
  @param xColStride   column stride for x
  @param yPtr         2nd operand
  @param yColStride   column stride for y
  @param resPtr       result
  @param resColStride column stride for res
  @param gpuId  GPU id, if -1: use the current one
*/
void multiplyGPU (ssize_t nRows, ssize_t nCols,
                  const double* xPtr, ssize_t xColStride,
                  const double* yPtr, ssize_t yColStride,
                  double* resPtr, ssize_t resColStride,
                  int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute x *= y (contiguous vector)
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param n      vector size
  @param gpuId  GPU id, if -1: use the current one
*/
void multiplyInPlaceGPU (double* xPtr, const double* yPtr,
      ssize_t n, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x / y (contiguous vector)
  @param resPtr target vector
  @param n      vector size
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param gpuId  GPU id, if -1: use the current one
*/
void divideGPU (double* resPtr, ssize_t n,
      const double* xPtr, const double* yPtr, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x / y (contiguous columns)
  @param resPtr target vector
  @param nRows        number of rows (contiguous in memory)
  @param nCols        number of columns (with stride >= nRows)
  @param xPtr         1st operand
  @param xColStride   column stride for x
  @param yPtr         2nd operand
  @param yColStride   column stride for y
  @param resPtr       result
  @param resColStride column stride for res
  @param gpuId  GPU id, if -1: use the current one
*/
void divideGPU (ssize_t nRows, ssize_t nCols,
                  const double* xPtr, ssize_t xColStride,
                  const double* yPtr, ssize_t yColStride,
                  double* resPtr, ssize_t resColStride,
                  int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute x /= y (contiguous vector)
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param n      vector size
  @param gpuId  GPU id, if -1: use the current one
*/
void divideInPlaceGPU (double* xPtr, const double* yPtr,
      ssize_t n, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x + y (contiguous vector)
  @param resPtr target vector
  @param n      vector size
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param gpuId  GPU id, if -1: use the current one
*/
void addGPU (SxComplex8* resPtr, ssize_t n,
      const SxComplex8* xPtr, const SxComplex8* yPtr, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x + y (contiguous columns)
  @param resPtr target vector
  @param nRows        number of rows (contiguous in memory)
  @param nCols        number of columns (with stride >= nRows)
  @param xPtr         1st operand
  @param xColStride   column stride for x
  @param yPtr         2nd operand
  @param yColStride   column stride for y
  @param resPtr       result
  @param resColStride column stride for res
  @param gpuId  GPU id, if -1: use the current one
*/
void addGPU (ssize_t nRows, ssize_t nCols,
                  const SxComplex8* xPtr, ssize_t xColStride,
                  const SxComplex8* yPtr, ssize_t yColStride,
                  SxComplex8* resPtr, ssize_t resColStride,
                  int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute x += y (contiguous vector)
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param n      vector size
  @param gpuId  GPU id, if -1: use the current one
*/
void addInPlaceGPU (SxComplex8* xPtr, const SxComplex8* yPtr,
      ssize_t n, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x - y (contiguous vector)
  @param resPtr target vector
  @param n      vector size
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param gpuId  GPU id, if -1: use the current one
*/
void subtractGPU (SxComplex8* resPtr, ssize_t n,
      const SxComplex8* xPtr, const SxComplex8* yPtr, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x - y (contiguous columns)
  @param resPtr target vector
  @param nRows        number of rows (contiguous in memory)
  @param nCols        number of columns (with stride >= nRows)
  @param xPtr         1st operand
  @param xColStride   column stride for x
  @param yPtr         2nd operand
  @param yColStride   column stride for y
  @param resPtr       result
  @param resColStride column stride for res
  @param gpuId  GPU id, if -1: use the current one
*/
void subtractGPU (ssize_t nRows, ssize_t nCols,
                  const SxComplex8* xPtr, ssize_t xColStride,
                  const SxComplex8* yPtr, ssize_t yColStride,
                  SxComplex8* resPtr, ssize_t resColStride,
                  int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute x -= y (contiguous vector)
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param n      vector size
  @param gpuId  GPU id, if -1: use the current one
*/
void subtractInPlaceGPU (SxComplex8* xPtr, const SxComplex8* yPtr,
      ssize_t n, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x * y (contiguous vector)
  @param resPtr target vector
  @param n      vector size
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param gpuId  GPU id, if -1: use the current one
*/
void multiplyGPU (SxComplex8* resPtr, ssize_t n,
      const SxComplex8* xPtr, const SxComplex8* yPtr, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x * y (contiguous columns)
  @param resPtr target vector
  @param nRows        number of rows (contiguous in memory)
  @param nCols        number of columns (with stride >= nRows)
  @param xPtr         1st operand
  @param xColStride   column stride for x
  @param yPtr         2nd operand
  @param yColStride   column stride for y
  @param resPtr       result
  @param resColStride column stride for res
  @param gpuId  GPU id, if -1: use the current one
*/
void multiplyGPU (ssize_t nRows, ssize_t nCols,
                  const SxComplex8* xPtr, ssize_t xColStride,
                  const SxComplex8* yPtr, ssize_t yColStride,
                  SxComplex8* resPtr, ssize_t resColStride,
                  int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute x *= y (contiguous vector)
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param n      vector size
  @param gpuId  GPU id, if -1: use the current one
*/
void multiplyInPlaceGPU (SxComplex8* xPtr, const SxComplex8* yPtr,
      ssize_t n, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x / y (contiguous vector)
  @param resPtr target vector
  @param n      vector size
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param gpuId  GPU id, if -1: use the current one
*/
void divideGPU (SxComplex8* resPtr, ssize_t n,
      const SxComplex8* xPtr, const SxComplex8* yPtr, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x / y (contiguous columns)
  @param resPtr target vector
  @param nRows        number of rows (contiguous in memory)
  @param nCols        number of columns (with stride >= nRows)
  @param xPtr         1st operand
  @param xColStride   column stride for x
  @param yPtr         2nd operand
  @param yColStride   column stride for y
  @param resPtr       result
  @param resColStride column stride for res
  @param gpuId  GPU id, if -1: use the current one
*/
void divideGPU (ssize_t nRows, ssize_t nCols,
                  const SxComplex8* xPtr, ssize_t xColStride,
                  const SxComplex8* yPtr, ssize_t yColStride,
                  SxComplex8* resPtr, ssize_t resColStride,
                  int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute x /= y (contiguous vector)
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param n      vector size
  @param gpuId  GPU id, if -1: use the current one
*/
void divideInPlaceGPU (SxComplex8* xPtr, const SxComplex8* yPtr,
      ssize_t n, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x + y (contiguous vector)
  @param resPtr target vector
  @param n      vector size
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param gpuId  GPU id, if -1: use the current one
*/
void addGPU (SxComplex16* resPtr, ssize_t n,
      const SxComplex16* xPtr, const SxComplex16* yPtr, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x + y (contiguous columns)
  @param resPtr target vector
  @param nRows        number of rows (contiguous in memory)
  @param nCols        number of columns (with stride >= nRows)
  @param xPtr         1st operand
  @param xColStride   column stride for x
  @param yPtr         2nd operand
  @param yColStride   column stride for y
  @param resPtr       result
  @param resColStride column stride for res
  @param gpuId  GPU id, if -1: use the current one
*/
void addGPU (ssize_t nRows, ssize_t nCols,
                  const SxComplex16* xPtr, ssize_t xColStride,
                  const SxComplex16* yPtr, ssize_t yColStride,
                  SxComplex16* resPtr, ssize_t resColStride,
                  int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute x += y (contiguous vector)
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param n      vector size
  @param gpuId  GPU id, if -1: use the current one
*/
void addInPlaceGPU (SxComplex16* xPtr, const SxComplex16* yPtr,
      ssize_t n, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x - y (contiguous vector)
  @param resPtr target vector
  @param n      vector size
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param gpuId  GPU id, if -1: use the current one
*/
void subtractGPU (SxComplex16* resPtr, ssize_t n,
      const SxComplex16* xPtr, const SxComplex16* yPtr, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x - y (contiguous columns)
  @param resPtr target vector
  @param nRows        number of rows (contiguous in memory)
  @param nCols        number of columns (with stride >= nRows)
  @param xPtr         1st operand
  @param xColStride   column stride for x
  @param yPtr         2nd operand
  @param yColStride   column stride for y
  @param resPtr       result
  @param resColStride column stride for res
  @param gpuId  GPU id, if -1: use the current one
*/
void subtractGPU (ssize_t nRows, ssize_t nCols,
                  const SxComplex16* xPtr, ssize_t xColStride,
                  const SxComplex16* yPtr, ssize_t yColStride,
                  SxComplex16* resPtr, ssize_t resColStride,
                  int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute x -= y (contiguous vector)
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param n      vector size
  @param gpuId  GPU id, if -1: use the current one
*/
void subtractInPlaceGPU (SxComplex16* xPtr, const SxComplex16* yPtr,
      ssize_t n, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x * y (contiguous vector)
  @param resPtr target vector
  @param n      vector size
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param gpuId  GPU id, if -1: use the current one
*/
void multiplyGPU (SxComplex16* resPtr, ssize_t n,
      const SxComplex16* xPtr, const SxComplex16* yPtr, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x * y (contiguous columns)
  @param resPtr target vector
  @param nRows        number of rows (contiguous in memory)
  @param nCols        number of columns (with stride >= nRows)
  @param xPtr         1st operand
  @param xColStride   column stride for x
  @param yPtr         2nd operand
  @param yColStride   column stride for y
  @param resPtr       result
  @param resColStride column stride for res
  @param gpuId  GPU id, if -1: use the current one
*/
void multiplyGPU (ssize_t nRows, ssize_t nCols,
                  const SxComplex16* xPtr, ssize_t xColStride,
                  const SxComplex16* yPtr, ssize_t yColStride,
                  SxComplex16* resPtr, ssize_t resColStride,
                  int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute x *= y (contiguous vector)
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param n      vector size
  @param gpuId  GPU id, if -1: use the current one
*/
void multiplyInPlaceGPU (SxComplex16* xPtr, const SxComplex16* yPtr,
      ssize_t n, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x / y (contiguous vector)
  @param resPtr target vector
  @param n      vector size
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param gpuId  GPU id, if -1: use the current one
*/
void divideGPU (SxComplex16* resPtr, ssize_t n,
      const SxComplex16* xPtr, const SxComplex16* yPtr, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x / y (contiguous columns)
  @param resPtr target vector
  @param nRows        number of rows (contiguous in memory)
  @param nCols        number of columns (with stride >= nRows)
  @param xPtr         1st operand
  @param xColStride   column stride for x
  @param yPtr         2nd operand
  @param yColStride   column stride for y
  @param resPtr       result
  @param resColStride column stride for res
  @param gpuId  GPU id, if -1: use the current one
*/
void divideGPU (ssize_t nRows, ssize_t nCols,
                  const SxComplex16* xPtr, ssize_t xColStride,
                  const SxComplex16* yPtr, ssize_t yColStride,
                  SxComplex16* resPtr, ssize_t resColStride,
                  int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute x /= y (contiguous vector)
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param n      vector size
  @param gpuId  GPU id, if -1: use the current one
*/
void divideInPlaceGPU (SxComplex16* xPtr, const SxComplex16* yPtr,
      ssize_t n, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x + y (contiguous vector)
  @param resPtr target vector
  @param n      vector size
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param gpuId  GPU id, if -1: use the current one
*/
void addGPU (int* resPtr, ssize_t n,
      const int* xPtr, const int* yPtr, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x + y (contiguous columns)
  @param resPtr target vector
  @param nRows        number of rows (contiguous in memory)
  @param nCols        number of columns (with stride >= nRows)
  @param xPtr         1st operand
  @param xColStride   column stride for x
  @param yPtr         2nd operand
  @param yColStride   column stride for y
  @param resPtr       result
  @param resColStride column stride for res
  @param gpuId  GPU id, if -1: use the current one
*/
void addGPU (ssize_t nRows, ssize_t nCols,
                  const int* xPtr, ssize_t xColStride,
                  const int* yPtr, ssize_t yColStride,
                  int* resPtr, ssize_t resColStride,
                  int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute x += y (contiguous vector)
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param n      vector size
  @param gpuId  GPU id, if -1: use the current one
*/
void addInPlaceGPU (int* xPtr, const int* yPtr,
      ssize_t n, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x - y (contiguous vector)
  @param resPtr target vector
  @param n      vector size
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param gpuId  GPU id, if -1: use the current one
*/
void subtractGPU (int* resPtr, ssize_t n,
      const int* xPtr, const int* yPtr, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x - y (contiguous columns)
  @param resPtr target vector
  @param nRows        number of rows (contiguous in memory)
  @param nCols        number of columns (with stride >= nRows)
  @param xPtr         1st operand
  @param xColStride   column stride for x
  @param yPtr         2nd operand
  @param yColStride   column stride for y
  @param resPtr       result
  @param resColStride column stride for res
  @param gpuId  GPU id, if -1: use the current one
*/
void subtractGPU (ssize_t nRows, ssize_t nCols,
                  const int* xPtr, ssize_t xColStride,
                  const int* yPtr, ssize_t yColStride,
                  int* resPtr, ssize_t resColStride,
                  int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute x -= y (contiguous vector)
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param n      vector size
  @param gpuId  GPU id, if -1: use the current one
*/
void subtractInPlaceGPU (int* xPtr, const int* yPtr,
      ssize_t n, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x * y (contiguous vector)
  @param resPtr target vector
  @param n      vector size
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param gpuId  GPU id, if -1: use the current one
*/
void multiplyGPU (int* resPtr, ssize_t n,
      const int* xPtr, const int* yPtr, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x * y (contiguous columns)
  @param resPtr target vector
  @param nRows        number of rows (contiguous in memory)
  @param nCols        number of columns (with stride >= nRows)
  @param xPtr         1st operand
  @param xColStride   column stride for x
  @param yPtr         2nd operand
  @param yColStride   column stride for y
  @param resPtr       result
  @param resColStride column stride for res
  @param gpuId  GPU id, if -1: use the current one
*/
void multiplyGPU (ssize_t nRows, ssize_t nCols,
                  const int* xPtr, ssize_t xColStride,
                  const int* yPtr, ssize_t yColStride,
                  int* resPtr, ssize_t resColStride,
                  int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute x *= y (contiguous vector)
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param n      vector size
  @param gpuId  GPU id, if -1: use the current one
*/
void multiplyInPlaceGPU (int* xPtr, const int* yPtr,
      ssize_t n, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x / y (contiguous vector)
  @param resPtr target vector
  @param n      vector size
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param gpuId  GPU id, if -1: use the current one
*/
void divideGPU (int* resPtr, ssize_t n,
      const int* xPtr, const int* yPtr, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute res = x / y (contiguous columns)
  @param resPtr target vector
  @param nRows        number of rows (contiguous in memory)
  @param nCols        number of columns (with stride >= nRows)
  @param xPtr         1st operand
  @param xColStride   column stride for x
  @param yPtr         2nd operand
  @param yColStride   column stride for y
  @param resPtr       result
  @param resColStride column stride for res
  @param gpuId  GPU id, if -1: use the current one
*/
void divideGPU (ssize_t nRows, ssize_t nCols,
                  const int* xPtr, ssize_t xColStride,
                  const int* yPtr, ssize_t yColStride,
                  int* resPtr, ssize_t resColStride,
                  int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief Compute x /= y (contiguous vector)
  @param xPtr   1st operand
  @param yPtr   2nd operand
  @param n      vector size
  @param gpuId  GPU id, if -1: use the current one
*/
void divideInPlaceGPU (int* xPtr, const int* yPtr,
      ssize_t n, int gpuId)
{
   // placeholder, needs to be implemented
   SX_EXIT;
}
// --- BINOP

// --- The next 108 lines were generated from snippets/SxComputeGPU.cpp snippet matmult
/** \brief GPU BLAS sgemm wrapper

  C <- (alpha) A B + beta C       (alpha is always 1, beta is a scalar)

  @param resMat   => C   (M x N)
  @param beta     => beta
  @param aMat     => A   (M x K)
  @param bMat     => B   (K x N)
  @param M        => M (number of rows of A & C)
  @param K        => K (number of cols of A, number of rows of B)
  @param N        => N (number of cols of B and C)
  @param lda      => column stride for A
  @param ldb      => column stride for B
  @param ldc      => column stride for C
  @param gpuId    GPU device id where this is to be done; -1 means
                  currently active device
 */
SX_EXPORT_MATH void
matmultGPU (float *resMat, const float beta, const float *aMat,
            const float *bMat, int M, int K,
            int N, int lda, int ldb, int ldc, int gpuId)

{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief GPU BLAS dgemm wrapper

  C <- (alpha) A B + beta C       (alpha is always 1, beta is a scalar)

  @param resMat   => C   (M x N)
  @param beta     => beta
  @param aMat     => A   (M x K)
  @param bMat     => B   (K x N)
  @param M        => M (number of rows of A & C)
  @param K        => K (number of cols of A, number of rows of B)
  @param N        => N (number of cols of B and C)
  @param lda      => column stride for A
  @param ldb      => column stride for B
  @param ldc      => column stride for C
  @param gpuId    GPU device id where this is to be done; -1 means
                  currently active device
 */
SX_EXPORT_MATH void
matmultGPU (double *resMat, const double beta, const double *aMat,
            const double *bMat, int M, int K,
            int N, int lda, int ldb, int ldc, int gpuId)

{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief GPU BLAS cgemm wrapper

  C <- (alpha) A B + beta C       (alpha is always 1, beta is a scalar)

  @param resMat   => C   (M x N)
  @param beta     => beta
  @param aMat     => A   (M x K)
  @param bMat     => B   (K x N)
  @param M        => M (number of rows of A & C)
  @param K        => K (number of cols of A, number of rows of B)
  @param N        => N (number of cols of B and C)
  @param lda      => column stride for A
  @param ldb      => column stride for B
  @param ldc      => column stride for C
  @param gpuId    GPU device id where this is to be done; -1 means
                  currently active device
 */
SX_EXPORT_MATH void
matmultGPU (SxComplex8 *resMat, const SxComplex8 &beta, const SxComplex8 *aMat,
            const SxComplex8 *bMat, int M, int K,
            int N, int lda, int ldb, int ldc, int gpuId)

{
   // placeholder, needs to be implemented
   SX_EXIT;
}

/** \brief GPU BLAS zgemm wrapper

  C <- (alpha) A B + beta C       (alpha is always 1, beta is a scalar)

  @param resMat   => C   (M x N)
  @param beta     => beta
  @param aMat     => A   (M x K)
  @param bMat     => B   (K x N)
  @param M        => M (number of rows of A & C)
  @param K        => K (number of cols of A, number of rows of B)
  @param N        => N (number of cols of B and C)
  @param lda      => column stride for A
  @param ldb      => column stride for B
  @param ldc      => column stride for C
  @param gpuId    GPU device id where this is to be done; -1 means
                  currently active device
 */
SX_EXPORT_MATH void
matmultGPU (SxComplex16 *resMat, const SxComplex16 &beta, const SxComplex16 *aMat,
            const SxComplex16 *bMat, int M, int K,
            int N, int lda, int ldb, int ldc, int gpuId)

{
   // placeholder, needs to be implemented
   SX_EXIT;
}
// --- matmult
#endif
