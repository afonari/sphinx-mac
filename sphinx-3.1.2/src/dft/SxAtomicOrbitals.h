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

#ifndef _SX_ATOMIC_ORBITALS_H_
#define _SX_ATOMIC_ORBITALS_H_

#include <SxDFT.h>
#include <SxVector.h>
#include <SxArray.h>
#include <SxAtomicStructure.h>
#include <SxRadBasis.h>
#include <SxRadialBasis.h>
#include <SxAoIndex.h>

/** \brief Container of atomic orbitals \f$
              \langle r | \mu_{i_s,i_a,n,l,m} \rangle
           \f$

    \b SxAtomicOrbitals = S/PHI/nX Atomic Orbitals

    This class is a container for atomic orbitals sampled on a radial
    grid. The atomic orbitals can be specified by the index of the
    species \f$i_s\f$ and the orbital type (different n, l).

    \sa      \ref page_dirac
    \ingroup group_dft
    \ingroup group_dirac
    \author  Sixten Boeck, Christoph Freysoldt, Björn Lange
  */
class SX_EXPORT_DFT SxAtomicOrbitals
{
   public:

      SxAtomicOrbitals ();
      SxAtomicOrbitals (const SxArray<SxArray<SxVector<double> > > &in,
                        const SxConstPtr<SxBasis> &basisPtrIn,
                        bool splineRepIn = false);

      /** \brief Constructor from list of l-values, waves to be set later
        @param lPhi         list of l-values :iSpecies,:iOrbType
        @param basisPtrIn   pointer to radial basis (optional)
        @param splineRepIn  set to true if wave functions will be provided in
                            spline representation
      */
      SxAtomicOrbitals (const SxArray<SxArray<int> > &lPhi,
                        const SxPtr<SxBasis> &basisPtrIn = SxPtr<SxBasis> (),
                        bool splineRepIn = false);

      /** \brief Constructor from map, waves to be set later
        @param map          reference orbital map
        @param basisPtrIn   pointer to radial basis (optional)
        @param splineRepIn  set to true if wave functions will be provided in
                            spline representation
      */
      SxAtomicOrbitals (const SxArray<SxArray<SxAoIndex> > &map,
                        const SxPtr<SxBasis> &basisPtrIn = SxPtr<SxBasis> (),
                        bool splineRepIn = false);
      /// \brief Copy Constructor
      SxAtomicOrbitals (const SxAtomicOrbitals &);

      /// \brief Read In Constructor
      SxAtomicOrbitals (SxBinIO &io);

      /// \brief Destructor
      virtual ~SxAtomicOrbitals ();

      /// \brief Assignment operator
      SxAtomicOrbitals& operator= (const SxAtomicOrbitals &) = default;

      void operator+= (const SxAtomicOrbitals &);
      void operator-= (const SxAtomicOrbitals &);

      SxAtomicOrbitals operator+ (const SxAtomicOrbitals &) const;
      SxAtomicOrbitals operator- (const SxAtomicOrbitals &) const;

      SxAtomicOrbitals operator* (double skalar) const;
      SxAtomicOrbitals operator/ (double skalar) const;

      void sumMPI() const;

      /// \brief Get number of Species
      inline int getNSpecies () const
      {
         return (int)muSet.getSize ();
      }

      /// \brief Get number of orbital types for one species
      inline int getNOrbTypes (int iSpecies) const
      {
         return (int)muSet(iSpecies).getSize ();
      }

      /// \brief Get number of orbital types for all species
      int getNOrbTypes () const;

      /// \brief Get number of orbitals for one species
      inline int getNOrbs (ssize_t iSpecies) const
      {
         return (int)refOrbMap(iSpecies).getSize ();
      }

      /** \brief Get an orbital shape
          @param is    species id
          @param iot   orbital type (m is not set!)
      */
      template <class Idx1, class Idx2> // could be int/ssize_t/SxAutoLoop&
      inline SxVector<double> &operator() (Idx1 &is, Idx2 &iot)
      {
         return muSet(is)(iot);
      }
      /** \brief Get an orbital shape (const).
          @param is    species id
          @param iot   orbital type (m is not set!)
      */
      template <class Idx1, class Idx2> // could be int/ssize_t/SxAutoLoop&
      inline const SxVector<double> &operator() (Idx1 &is, Idx2 &iot) const
      {
         return muSet(is)(iot);
      }

      /** \brief Get an orbital with all of is,n,l,m set in auxData (ia=-1)
          @param is    species id
          @param io    local orbital id for (n,l,m)
          @note This is in contrast to operator(is, iot) which yields the m-independent shape
          See also @refOrbMap
      */
      template <class Idx1, class Idx2> // could be int/ssize_t/SxAutoLoop&
      inline const SxVecRef<double> getOrb (Idx1 &is, Idx2 &io) const
      {
         const SxAoIndex &nlm = refOrbMap(is)(io);
         SxVecRef<double> res = const_cast<SxVector<double>&>(muSet(is)(nlm.n));
         res.auxData.setAtom (is, -1);
         res.auxData << nlm;
         return res;
      }

      int getIOT (int is, int n, int l) const;

      /** \brief Extract a single orbital

          In principle only \em iSpecies and \em l are needed to extract
          a single orbital. However, for some projections (e.g. onto
          plane-waves) a Dirac vector also need to know it's atomic index and
          the quantum numbers. This function extracts an orbital and
          initializes the necessary auxillary Dirac data. */
      const SxVecRef<double>
      operator() (int is, int ia, int iot, int l, int m) const;

      /// \brief Set all functions to a value (usually zero)
      void set (double val)
      {
         SX_CHECK (getNSpecies () > 0);
         SX_LOOP2 (iSpecies,iot)
            muSet(iSpecies)(iot).set(val);
      }
      SxArray<SxQuantumNumbers> getReducedOrbitalMap () const;

      SxArray<SxQuantumNumbers> getOrbitalMap (const SxAtomicStructure &structure) const;

      int getNOrbitals (const SxAtomicStructure &structure) const;

      void write (SxBinIO &io) const;

      void write (const SxString &filename) const;

      void read (SxBinIO &io);

      void read (const SxString &file);

      void readSiesta (const SxString &file,
                       const SxVecRef<double> &radFunc);

      void writeSiesta (const SxArray<SxString> &file,
                        const SxArray<SxVector<double> > &basis);


      void print (const SxString &file) const;

      void setup (const SxSymbolTable *table);

      void addOrbital (const SxVecRef<double> &orbitalIn);

      /** \brief Read an orbital from a file
        @param fp an open file, containing lines of "r psi(r)"
        @param is  species
        @param n   main index
        @param l   l channel
        @param iot  if set, read to orbital with index iot.
                    Append orbital otherwise.
        @param ignoreComments if # line comments should be ignored
        */
      void readOrbital (FILE *fp, int is, int n, int l, int iot=-1,
                        bool ignoreComments = true);
      /** \brief Read orbitals from a file
        @param fp   an open file, containing lines of "r psi(r)" preceded by
                    comment lines that give l= (and optionally n=)
        @param is   species
       */
      void readOrbitals (FILE *fp, int is);

      double getNormSqr (int is, int iot) const;
      double getNormSqrSum () const;

      double dot (const SxAtomicOrbitals &orbitalsIn, int is, int iot, int js, int jot) const;
      double dot (const SxAtomicOrbitals &orbitalsIn) const;

      double sum (const SxAtomicOrbitals &orbitalsIn, int is, int iot, int js, int jot) const;
      double sum (const SxAtomicOrbitals &orbitalsIn) const;

      void normalize ();

      /// Set the radial basis pointer
      void setBasis (const SxConstPtr<SxBasis> &basisPtrIn);

      const SxRadBasis& getRadBasis () const { return *getRadBasisPtr (); }
      SxConstPtr<SxRadBasis> getRadBasisPtr () const {
         return SxConstPtr<SxRadBasis>(basisPtr);
      }

      SxConstPtr<SxRadialBasis> getRadGBasisPtr () const { return basisPtr; };

      /** \brief Get the radial mesh for species iSpecies */
      const SxVecRef<double> &getRadFunc (int iSpecies) const
      {
         SX_CHECK (basisPtr);
         if (const SxRadBasis *rad
               = dynamic_cast<const SxRadBasis*>(basisPtr.getPtr ()))
            return rad->radFunc(iSpecies);
         else if (const SxRadialBasis *radial
               = dynamic_cast<const SxRadialBasis*>(basisPtr.getPtr ()))
            return radial->getRadFunc ();
         SX_EXIT; // wrong basis?
      }

      void toSpline () const;

      void toVec () const;

      SxVector<double> toVec (int is, int iot) const;

      bool isSpline () const {return splineRep;};

      void createFuncLMap ();

      SxVector<double>& getFuncL (int is, int l, int ifl);

      const SxVector<double>& getFuncL (int is, int l, int ifl) const;

      SxArray<SxArray<int> > getFuncPerL () const;

      int getFuncPerL (int is, int l) const;

      SxArray<SxArray<SxVector<double> > > orthogonalize ();

      SxArray<SxArray<SxVector<double> > > getOverlap () const;

      SxVector<double> getOverlap (int is, int l) const;

      SxArray<SxArray<SxArray<int> > > funcLMap;

      int getOrbitalIdx (int is, int ia, int iot, int l, int m,
                         const SxArray<SxQuantumNumbers> &map) const;

      int getLMax () const;

      inline int getLMax (const int iSpecies) const
      {
         return int(funcLMap(iSpecies).getSize ()) - 1;
      }

      void orthogonalizeOn(SxAtomicOrbitals &basis);

      void rotate(const SxArray<SxArray<SxVector<double> > > &rotMat);

      void refine (SxArray<int> &factor);

      SxVector<double>
      compressWave(SxString &file, int iState, int iSpecies, int l);

      #ifdef USE_HDF5
      void writeHDF5(const SxString &name);
      #endif

   protected:

      mutable bool splineRep;

      SxConstPtr<SxBasis> basisPtr;

      /** \brief The sampling points of each orbitals.

        <b>Storage order</b>:
        -# iSpecies
        -# l
        -# r */
      mutable SxArray<SxArray<SxVector<double> > >  muSet;  // :is,:iot,:r

   public:
      /// Reference orbital map iot -> n,l,m
      SxArray<SxArray<SxAoIndex> > refOrbMap;

      /// Get the muSet
      const SxArray<SxArray<SxVector<double> > > & getMuSet () const
      {
         return muSet;
      }

      /// Inquire l-number of some orbital type
      int getL (int is, int iot) const
      {
         SX_CHECK (is >= 0 && is < muSet.getSize (), is, muSet.getSize ());
         SX_CHECK (iot >= 0 && iot < muSet(is).getSize (),
                   iot, muSet(is).getSize ());
         return muSet(is)(iot).auxData.l;
      }

};

inline SX_EXPORT_DFT
SxAtomicOrbitals operator* (double skalar, const SxAtomicOrbitals &in)
{
   return in * skalar;
}

#endif /* _SX_ATOMIC_ORBITALS_H_ */
