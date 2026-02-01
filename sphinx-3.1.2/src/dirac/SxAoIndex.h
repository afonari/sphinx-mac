
// ---------------------------------------------------------------------------
//
//      The ab-initio based multiscale library
//
//                  S / P H I / n X
//
//      Copyright:  Max-Planck-Institut fuer Eisenforschung GmbH
//                  40237 Duesseldorf, Germany
//
//      Contact:    https://sxrepo.mpie.de
//      Authors:    see sphinx/AUTHORS
//      License:    see sphinx/LICENSE
//
// ---------------------------------------------------------------------------

#ifndef _SX_AO_INDEX_H_
#define _SX_AO_INDEX_H_

#include <SxDiracLib.h>

/**
  \brief This is a container class for atomic orbital indices n,l,m

  \author C.Freysoldt, freysoldt@mpie.de
  */
class SX_EXPORT_DIRAC SxAoIndex  {
   public:
      int n, l, m;
      SxAoIndex () : n(-1), l(-1), m(-1) {}
      SxAoIndex (int n_, int l_, int m_) : n(n_), l(l_), m(m_) {}
      // The standard copy and assignment operators work fine
      // because there's no memory management here.
      // i.e. these two are automatically defined:
      // SxAoIndex (const SxAoIndex);
      // operator= (const &SxAoIndex);

      inline void set (int nIn, int lIn, int mIn)
      {
         n=nIn; l=lIn; m=mIn;
      }

      inline bool operator== (const SxAoIndex &in) const
      {
         return (in.n == n && in.l == l && in.m == m);
      }
      inline bool operator!= (const SxAoIndex &in) const
      {
         return (in.n != n || in.l != l || in.m != m);
      }
      /** \brief Ordering comparison (hierarchy 1) n 2) l 3) m).

        \note
        (n1, l1, m1) < (n2, l2, m2) means
        -# n1 != n2 => n1 < n2
        -# n1 == n2 && l1 != l2 => l1 < l2
        -# n1 == n2 && l1 == l2 => m1 < m2
        */
      inline bool operator< (const SxAoIndex &in) const
      {
         return (n == in.n) ? (l == in.l ? (m < in.m)
                                         : (l < in.l) )
                            :              (n < in.n) ;
      }
      /** \brief Ordering comparison (hierarchy 1) n 2) l 3) m).

        \note
        (n1, l1, m1) > (n2, l2, m2) means
        -# n1 != n2 => n1 > n2
        -# n1 == n2 && l1 != l2 => l1 > l2
        -# n1 == n2 && l1 == l2 => m1 > m2
        */
      inline bool operator> (const SxAoIndex &in) const
      {
         return (n == in.n) ? (l == in.l ? (m > in.m)
                                         : (l > in.l) )
                            :              (n > in.n) ;
      }
};

inline
void operator<< (SxAuxData &aux, const SxAoIndex &nlm)
{
   aux.setOrb (nlm.n, nlm.l, nlm.m);
}

/**
  \brief This is a container class for atomic orbital indices is,ia,io
  */
class SxOrbitalIndex  {
   public:
      int is, ia, io;
      SxOrbitalIndex () : is(-1), ia(-1), io(-1) {}
      SxOrbitalIndex (int is_, int ia_, int io_)
        : is(is_), ia(ia_), io(io_) {}
      // The standard copy and assignment operators work fine
      // because there's no memory management here.
      // i.e. these two are automatically defined:
      // SxOrbitalIndex (const SxOrbitalIndex);
      // operator= (const &SxOrbitalIndex);

      inline void set (int iSpecies, int iAtom, int iOrbital)
      {
         is = iSpecies;
         ia = iAtom;
         io = iOrbital;
      }

      inline bool operator== (const SxOrbitalIndex &in) const
      {
         return (in.is == is && in.ia == ia && in.io == io);
      }
      inline bool operator!= (const SxOrbitalIndex &in) const
      {
         return (in.is != is || in.ia != ia || in.io != io);
      }
      /** \brief Ordering comparison (hierarchy 1) is 2) ia 3) io).

        \note
        (is, ia, io) < (is2, ja, jo) means
        -# is != js => is < js
        -# is == js && ia != ja => ia < ja
        -# is == js && ia == ja => io < jo
        */
      inline bool operator< (const SxOrbitalIndex &in) const
      {
         return (is == in.is) ? (ia == in.ia ? (io < in.io)
                                             : (ia < in.ia) )
                            :                  (ia < in.is) ;
      }
      /** \brief Ordering comparison (hierarchy 1) is 2) ia 3) io).

        \note
        (is, ia, io) > (js, ja, jo) means
        -# is != js => is > js
        -# is == js && ia != ja => ia > ja
        -# is == js && ia == ja => io > jo
        */
      inline bool operator> (const SxOrbitalIndex &in) const
      {
         return (is == in.is) ? (ia == in.ia ? (io > in.io)
                                             : (ia > in.ia) )
                            :                  (ia > in.is) ;
      }
};

#endif /* _SX_AO_INDEX_H_ */
