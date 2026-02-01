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

#include <SxUtil.h>
#include <SxVector.h>
#include <SxMathLib.h>
#include <SxConstants.h>  /* for complex I */
#include <SxEigensystem.h>
#include <SxSingularValDecomp.h>
#include <stdio.h>
#include <iostream>

#define TYPE_M1 SxComplex8
#define TYPE_M2 SxComplex16
#define TYPE_M3 float
#define TYPE_M4 double

template<class T>
void print (const SxSingularValDecomp<T> &svd, const SxVector<T> &mat)
{
   typedef SxVecRef<typename SxTypeMapper<T>::TReal> RealVec;
   cout << "values: " << svd.vals << endl;
   cout << "left: " << svd.left << endl;
   SxVector<T> ovlpL = svd.left.overlap (svd.left);
   cout << RealVec(ovlpL) << endl;
   cout << (ovlpL - RealVec(ovlpL)).normSqr () << endl;
   cout << "right: " << svd.right << endl;
   SxVector<T> ovlpR = svd.right.overlap (svd.right);
   cout << RealVec(ovlpR) << endl;
   cout << (ovlpR - RealVec(ovlpR)).normSqr () << endl;

   cout << "Test vals" << endl;
   for (int i = 0; i < svd.vals.getSize (); ++i)  {
      cout << ((mat ^ svd.right.colRef (i))
               - svd.vals(i) * svd.left.colRef (i)).normSqr ()
           << endl;
   }
}

/**
  \example sxsingval.cpp
  \author  C. Freysoldt, based on sxeigen.cpp by Sixten Boeck
  */
int main ()
{
   initSPHInXMath ();

   int n = 4;
   SxVector<TYPE_M1> m1(n,n);
   SxVector<TYPE_M2> m2(n,n);
   SxVector<TYPE_M3> m3(n,n);
   SxVector<TYPE_M4> m4(n,n);

   // --- set up matrix m1 and m3
   int r, c;
   for (r=0; r < n; r++)  {
      for (c=r; c < n; c++)  {
         m1(r,c) = TYPE_M2(sqrt((double)(r*c+1)),
                           cbrt((double)(r*c+1)));
         m1(c,r) = m1(r,c).conj();
         m3(r,c) = ((double)(r+c+1));
         m3(c,r) = m3(r,c);
         if (r==c)  m1(r,c).im = 0.;
      }
   }
   // --- initialize matrix m2/4. m2/4 = m1/3, but trigonally packed
   //     (use only upper-right triangle)
   for (r=0; r < n; r++)
      for (c=0; c < n; c++)  {
         m2(r,c) = TYPE_M2(sqrt((double)(r*c+1)),
                           cbrt((double)(r*c+1)));
         m4(r,c) = sqrt((double)(r*c+1));
         if (r==c)  m2(r,c).im = 0.;
      }

   // --- print matrices
   cout << "m1\n" << m1 << endl;
   cout << "m2\n" << m2 << endl;
   cout << "m3\n" << m3 << endl;
   cout << "m4\n" << m4 << endl;

   // --- compute eigensystems
   SxEigensystem<TYPE_M1> eig1(m1);
   SxEigensystem<TYPE_M3> eig3(m3);
   SxEigensystem<TYPE_M2> eig2(m2);
   SxEigensystem<TYPE_M4> eig4(m4);

   cout << "e1: "; eig1.print(true);
   cout << "e2: "; eig2.print(true);
   cout << "e3: "; eig3.print(true);
   cout << "e4: "; eig4.print(true);

   SxSingularValDecomp<TYPE_M1> svd1(m1);
   SxSingularValDecomp<TYPE_M3> svd3(m3);
   SxSingularValDecomp<TYPE_M2> svd2(m2);
   SxSingularValDecomp<TYPE_M4> svd4(m4);

   cout << "Singular value decomp 1" << endl;
   print (svd1, m1);
   cout << "Singular value decomp 2" << endl;
   cout << "Norm of eigenvals" << sqrt(eig2.vals.absSqr ()) << endl;
   print (svd2, m2);
   cout << "Singular value decomp 3" << endl;
   print (svd3, m3);
   cout << "Singular value decomp 4" << endl;
   print (svd4, m4);

   return 0;
}

