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

#ifndef _SX_SINGULAR_VAL_H_
#define _SX_SINGULAR_VAL_H_

#include <SxMath.h>
#include <SxVector.h>

/** \brief Container for singular value decomposition of a matrix

    \author Christoph Freysoldt, freysoldt@mpie.de */
template<class T>
class SX_EXPORT_MATH SxSingularValDecomp
{
   public:
      /// The eigenvalues
      SxVector<typename SxTypeMapper<T>::TReal> vals;
      /** \brief The left and right transformations

          If the singular value decomposition reads
          \f$ A = U \Sigma V^\dagger \f$
          "left" is U and "right" is V.

          This is analogous to eigenvectors: multiplying A
          from the right with column <i> of V yields the
          corresponding column <i> of U multiplied with the
          singular value <i>.
        */
      SxVector<T> left, right;

      /// compute routine
      void compute (SxVector<T> &&mat, EIGCMD cmd);

      /// Constructor
      SxSingularValDecomp (const SxVector<T> &mat, EIGCMD cmd = All)
      {
         compute (SxVector<T> (mat), cmd);
      }

      /// Constructor from an about-to-die matrix
      SxSingularValDecomp (SxVector<T> &&mat, EIGCMD cmd = All)
      {
         compute (std::move (mat), cmd);
      }

};

template<class T>
void SxSingularValDecomp<T>::compute (SxVector<T> &&mat, EIGCMD cmd)
{
   SX_CHECK (cmd == All || cmd == ValuesOnly, cmd);
   SX_CHECK (mat.getSize () > 0);
   SX_CHECK (mat.isRAM ()); // no GPU support yet
   int m = int(mat.getNRows ());
   int n = int(mat.getNCols ());
   SX_CHECK(ssize_t(n) * ssize_t(m) == mat.getSize ()); // otherwise: int overflow
   int minNM = min(n,m);

   vals.resize (minNM);
   if (cmd == All)  {
      left.reformat (m, minNM);
      right.reformat (minNM, n);
   } else {
      left.resize (0);
      right.resize (0);
   }

   singularValueDecomp (mat.elements, m, n, vals.elements, left.elements, right.elements, false);
   // we want V rather than V^H
   right = right.adjoint ();

   int nNonZero = minNM;
   for ( ; nNonZero>0 ; --nNonZero)
      if (vals(nNonZero-1) > 0.) break;
   if (nNonZero < minNM)  {
      // replace val/left/right by reference to nonZero part
      // this does not copy the data.
      SxVecRef<typename SxTypeMapper<T>::TReal> valsNonZero = vals.getRef(nNonZero);
      vals = valsNonZero;
      if (cmd == All)  {
         SxVecRef<T> leftNonZero = left.template getRef<Compact> (0, m, nNonZero);
         left = leftNonZero;
         SxVecRef<T> rightNonZero = right.template getRef<Compact> (0, n, nNonZero);
         right = rightNonZero;
      }
   }
}

#endif /* _SX_SINGULAR_VAL_H_ */
