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

#ifndef _SX_AO_BASIS_H_
#define _SX_AO_BASIS_H_

#include <SxPrecision.h>
#include <SxBasis.h>
#include <SxGBasis.h>
#include <SxVector.h>
#include <SxArray.h>
#include <SxGkBasis.h>
#include <SxPtr.h>
#include <SxDiracLib.h>
#include <SxOverlap.h>
#include <SxPWOverlap.h>
#include <SxAoIndex.h>


/** \brief Abstract base class for AO reference orbital containers
  */
class SxRefAOContainerBase {
   public:
      /// Destructor
      virtual ~SxRefAOContainerBase () {}
      /// Get reference orbital
      virtual SxVecRef<double>
      getOrb (ssize_t iSpecies, ssize_t iRef) const = 0;
      /// Get l of reference orbital
      virtual int getL (ssize_t iSpecies, ssize_t iRef) const = 0;
      /// Get number of species
      virtual int getNSpecies () const = 0;
      /// Get number of reference orbitals
      virtual int getNOrb (ssize_t iSpecies) const = 0;
};

/** Template container class for array (is) of array (iRef)
  */
template<class ArrayT>
class SxRefAOContainer : public SxRefAOContainerBase
{
   private:
      const ArrayT &psi;
   public:
      SxRefAOContainer (const ArrayT &psiIn)
         : psi(psiIn)
      { /* empty */ }
      /// Destructor
      virtual ~SxRefAOContainer () {}

      /// Get reference orbital
      SxVecRef<double> getOrb (ssize_t iSpecies, ssize_t iRef) const override
      {
         return const_cast<SxVecRef<double>&>
                ( static_cast<const SxVecRef<double>&>(psi(iSpecies)(iRef)) );
      }
      /// Get l of reference orbital
      int getL (ssize_t iSpecies, ssize_t iRef) const override
      {
         return getOrb(iSpecies, iRef).auxData.l;
      }
      /// Get number of species
      int getNSpecies () const override { return int(psi.getSize ()); }
      /// Get number of reference orbitals
      int getNOrb (ssize_t iSpecies) const override
      {
         return int(psi(iSpecies).getSize ());
      }
};


/** This class describes an abstract atomic orbital basis. The form
    of the atomic orbitals is not specified.

    \b SxAOBasis = S/PHI/nX Atomic Orbital Basis

    \ingroup group_dft
    \author  Christoph Freysoldt
    */
class SX_EXPORT_DIRAC SxAOBasis : public SxBasis
{
   public:

      typedef SxComplex16         CoeffType;
      typedef SxVector<CoeffType> TPsi;
      
      /// Return basis type
      virtual SxString getType () const { return "|mu>"; }
      
      /** \brief Caching behaviour for overlap matrices and their inverse
        */
      enum Caching {
         /// undefined behaviour
         Unknown, 
         /// recompute on demand
         Recompute, 
         /// cache one matrix
         CacheCurrentK, 
         /// cache all matrices
         CacheAll };
      
      /// Empty constructor
      SxAOBasis ();

      /** \brief Lightweight wrapper to get a AO ref orbital container
                 for array of array of psi */
      template<class ArrayT>
      SxRefAOContainer<ArrayT> wrapContainer (const ArrayT& psi)
      {
         return SxRefAOContainer<ArrayT>(psi);
      }

      /** \brief Constructor from pseudowaves
          \sa set
        */
      SxAOBasis (const SxGkBasis &gk, 
                 const SxArray<SxArray<SxVector<double> > > &psi,
                 const SxConstPtr<SxOverlapBase> SPtrIn = SxPtr<SxPWOverlap>::create ());

      /** \brief Constructor from pseudowaves
         @param gk       G+k basis
         @param psi      radial orbitals (:iSpecies)
         @param l-values l-values for these
        */
      SxAOBasis (const SxGkBasis &gk, 
                 const SxArray<SxVector<double> > &psi,
                 const SxArray<SxArray<int> >     &lPsi);

      /** \brief Set to pseudowaves
          \param gk  G+k basis
          \param psi radial wavefunctions is:l:ir

          @note See also #wrapContainer
        */
      void set (const SxGkBasis &gk, 
                const SxRefAOContainerBase &psi);
      /** \brief Single-species single-k "light" constructor */
      SxAOBasis (const SxVecRef<SxComplex16> &YlmGl, int iSpecies);

   protected:
      /// Initialize timer
      void init ();
      /// Compute reference orbitals
      void computeRefOrbitals(const SxGkBasis &gk,
                              const SxRadBasis &rad,
                              const SxArray<SxArray<SxVector<double> > > &psiRad);
   public:
      /// Caching behaviour for reference orbitals
      enum Caching cacheRefOrb;
   protected:
      /// k-point that is cached for reference orbitals
      mutable int refOrbCachedK;
   public:
      /// Pointer to Overlap Operator
      SxConstPtr<SxOverlapBase> SPtr;


      /// Destructor
      virtual ~SxAOBasis ();

      REGISTER_PROJECTOR (SxAOBasis, identity);
      REGISTER_PROJECTOR (SxGBasis, toPWBasis);

      /** \brief Identity projection, does not do anything.
        \note This allows to write
        \code
psiAO = (ao | psi);
        \endcode
        even if psi is in an AO basis already.
        */
      SxVecRef<CoeffType> identity  (const SxAOBasis *,
                                     const SxVecRef<CoeffType> &) const;
      /**
        \brief Project to plane-wave basis

       */
      SxVector<PrecCoeffG> toPWBasis (const SxGBasis *,
                                      const SxVecRef<CoeffType> &) const;

   protected:
#     ifdef USE_SXGEMMM
      /// sxpgemm3m implementation
      SxVector<PrecCoeffG>
      toPWBasisGemm3m (const SxGBasis *gPtr, int ik,
                       const SxVecRef<SxComplex16> &psiAO) const;
#     else
      /// BLAS implementation
      SxVector<PrecCoeffG>
      toPWBasisBLAS   (const SxGBasis *gPtr, int ik,
                       const SxVecRef<SxComplex16> &psiAO) const;
#     endif

   public:
      /**
        \brief Project from plane-wave basis

        \note The projection routine is located here, SxGBasis::toAO just 
              calls this function.
       */
      SxAOBasis::TPsi fromPWBasis (const SxVecRef<CoeffType> &psiG) const;

      /**
        \brief Project from plane-wave basis

       Projection from SxGBasis or SxPAWBasis
       */
      SxAOBasis::TPsi fromPWBasis (const PsiRef &psiS, int ik) const;

      /** \brief Derivative of orbital projections with respect to atomic
                 position
        */
      SxArray<SxVector<SxComplex16> > gradProject (const PsiRef &psi) const;

      /** \brief Derivative of orbital projections with respect to atomic
                 position
          @param psiS  wavefunction with this basis' overlap operator applied
                      (if applicable)
          @param ik   Gk index of reference orbitals to use

          @note This routine is the interface for projections from a SxPAWBasis
                in the absence of an overlap operator in the AO basis. In this
                case, we simply ignore the PAW projection coefficients behind
                the G+k coefficients.

                The routine does NOT check that the PAW basis matches the
                G+k basis of the reference orbitals.

        */
      SxArray<SxVector<SxComplex16> >
      gradProject (const PsiRef &psiS, int ik) const;

      /** \brief The reference orbitals in a |G+k>-Basis ik:is:(ig,iOrb)
          
          For each species, the orbitals are stored once and shifted to
          the appropiate atomic position whenever needed. This reduces
          the necessary amount of memory.
          The iOrb index is a condensed index (n,l,m) and can be decrypted
          with refOrbMap
        */
      SxArray<SxArray<PsiG> > refOrbitals;

      /** Additional phase per orbital */
      SxArray<SxVector<SxComplex16> > extraPhase;

/** \brief Maps iRefOrb to (n,l,m) for each species (is:iRefOrb)

          \example
          \code
SxAoIndex nlm = aoBasis.refOrbMap(is)(iRefOrb);
int n = nlm.n;
int l = nlm.l;
int m = nlm.m;
\endcode
        */
      SxArray<SxArray<SxAoIndex> > refOrbMap;

/** \brief Maps iOrb to (is,ia,io) for each orbital (mu)

          \example
          \code
int mu;
SxOrbitalIndex sao = aoBasis.orbitalMap(mu);
SxAoIndex nlm = aoBasis.refOrbMap(sao.is)(sao.io);
int is = sao.is;
int ia = sao.ia;
int n = nlm.n;
int l = nlm.l;
int m = nlm.m;
\endcode
*/
      SxArray<SxOrbitalIndex> orbitalMap;

      int getNSpecies () const { return int(refOrbMap.getSize ()); }
      int getNOrb () const { return int(orbitalMap.getSize ()); }
      /// Return number of atomic orbitals
      virtual ssize_t getNElements () const  { return getNOrb (); }
      int getLMax () const;

      /// Blocksize for blocked projections
      int blockSize;
/** \brief Return all orbitals in |G+k> basis
        */
      PsiG getAOinG (int ik) const;
      /** \brief Return a single orbital in |G+k> basis
        */
      PsiG getAOinG(int ik, int iOrb) const;

   protected:
      //\name Overlap 
      //\{
      /** \brief The overlap matrices ik:(mu,nu)

          This stores the matrices
          \f[ S^k_{\mu,\nu} = \int dr^3 (\chi^k_{\mu})^*(r) \chi^k_nu(r) \f]
        
        */
      mutable SxArray<SxVecRef<PrecCoeffG> > overlap;
      /// Caching behaviour for overlap matrices
      enum Caching cacheOverlap;
      /// k-point that is cached for overlap
      mutable int overlapCachedK;

      /** \brief The overlap matrices ik:(mu,nu)

          This stores the matrices
          \f[ S^k_{\mu,\nu} = \int dr^3 (\chi^k_{\mu})^*(r) \chi^k_nu(r) \f]
        
        */
      mutable SxArray<SxVecRef<PrecCoeffG> > invOverlap;
      /// Caching behaviour for inverse overlap matrices
      enum Caching cacheInverse;
      /// k-point that is cached for inverse overlap
      mutable int invOverlapCachedK;
      
   public:
      /// Get overlap matrix for k-point ik
      SxVecRef<SxComplex16> getOverlap (int ik) const;
      /// Get inverse overlap matrix for k-point ik
      SxVecRef<SxComplex16> getInverseOverlap (int ik) const;
      /**  \brief Compute the overlap matrix
        \note You may consider caching the overlap matrices instead
              of computing them whenever needed.
        \sa setOverlapCaching, getOverlap 
        */
       SxVector<SxComplex16> calculateOverlap (int ik) const;


      /** \brief Set caching behaviour for overlap matrices
          @param mode the mode
          @param nk number of k-points to cache if mode is CacheAll
                 and k-dependent reference orbitals are not cached.
        */
      void setOverlapCaching (enum Caching mode, int nk = -1);
      /** \brief Set caching behaviour for inverse overlap matrices
          @param mode the mode
          @param nk number of k-points to cache if mode is CacheAll
                 and k-dependent reference orbitals are not cached.
       */
      void setInvOverlapCaching (enum Caching mode, int nk = -1);
      //\}

};

#endif /* _SX_AO_BASIS_H_ */

