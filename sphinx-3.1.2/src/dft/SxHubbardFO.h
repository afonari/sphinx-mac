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

#ifndef _SX_HUBBARD_FO_H_
#define _SX_HUBBARD_FO_H_

#include <SxDFT.h>
#include <SxHubbardXO.h>
#include <SxAtomicOrbitals.h>

class SxHubbardU;

/** \brief FO projectors for DFT + U

    \author Christoph Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DFT SxHubbardFO : public SxHubbardXO
{

   public:
      /// Mapping class for symmetrizing block density matrix
      class SymMap;

   public:
      /// Constructor
      SxHubbardFO (int siteOffsetIn = 0);

      /// Destructor
      virtual ~SxHubbardFO () = default;

      // --- abstract base class interface
      /// Read and setup for homonuclear diatomics
      virtual void read (const SxSymbolTable *table,
                         const SxAtomicStructure &structure,
                         const SxPtr<SxPAWPot> &potPtr) override;
      /// Setup AO projectors in G+k space
      virtual void
      setupProjGk (const SxGkBasis &gk,
                   const SxPtr<SxPartialWaveBasis> &pBasis) override;

      /// Remove intermediate data used for setup
      virtual void finalize () { /* empty */ }

      /// Return number of atomic orbitals for site iSite
      virtual int getNOrbSite (ssize_t iSite) const
      {
         return nFragOrb(iSite);
         //return (int)siteU(iSite).getNCols ();
      }

      /// Return number of correlated orbitals for site iSite
      virtual int getNOrbHam (ssize_t iSite) const
      {
         return (int)siteU(iSite).getNCols ();
      }

      /// Number of atoms for a given site
      virtual int getNAtoms (ssize_t iSite) const override
      {
         return int(atomIdx(iSite).getSize ());
      }

      /** \brief Get global atom id for a given atom
          @param iSite   which site
          @param iAtom   atom id within that site
          \note This must be consistent with atom ordering in trafoForces
         */
      virtual int getIAtom (ssize_t iSite, ssize_t iAtom) const override
      {
         return atomIdx(iSite)(iAtom);
      }

      /** \brief Number of orbitals associated with a specific atom
          @param iSite    which site
          @param iAtom    atom id within that site
        */
      virtual int getNOrbAtom (ssize_t iSite, ssize_t iAtom) const override;

      /** \brief Apply the +U Hamiltonian
          \note The return vector must be projectable to the G+k basis.
        */
      virtual PsiG apply (const SxVecRef<SxComplex16> &psi) const override;

      /** \brief Compute the Hubbard energy and Hamiltonian
        @param hubbardU the Hubbard parent object containing the U parameters,
                        energy, etc.
        @param Pij the AO block density matrix
        @param structure the atomic structure
        */
      virtual void compute (SxHubbardU *hubbardU,
                            const SxBlockDensityMatrix& Pij,
                            const SxAtomicStructure &structure) override;

      /** \brief Add contribution of one k-point to the AO block density matrix
        @param Pij The block density matrix to be computed
        @param waves wave functions for current k-point
        @param weight k-point weight
        @param focc occupation number for each state
        */
      virtual void addToRho (SxBlockDensityMatrix *Pij,
                             const SxVecRef<PrecCoeffG> &waves,
                             double weight,
                             const SxVecRef<double> &focc) const override;
      /** \brief Compute contribution to gradient of the AO block density matrix

          @param fermi Fermi occupations
          @param P     <AO|Psi> projections
          @param gradP <d/dtau AO | Psi> projection gradients
          @return gradient matrices for all sites. The first ("spin") index is
                  used for the direction.

          @note k-point and spin must be set in the projection's (P) auxData.
        */
      virtual SxBlockDensityMatrix
      computeGradPij (const SxFermi &fermi,
                      const SxVecRef<SxComplex16> &P,
                      const SxArray<SxVector<SxComplex16> > &gradP) const override;

      /// Get symmetry map
      virtual SxPtr<SxHubbardXO::SymMap>
      getSymMap (const SxAtomicStructure &structure) const override;

      // --- end of abstract base class interface

   public:
      /// List of central atom indices
      SxArray<int> centralAtoms;
      /// Maximum distance between central atom and ligands
      double maxDist;
      /// Minimum number of atoms in fragment
      ssize_t minNumAtoms;

      /// Orbital energies per ao type
      SxArray<double> orbEnergy;

      /// Wolfsberg-Helmholz coupling factor for interatomic matrix elements
      double kappa;

      /// map fragment orbital to AO basis function type (:iFragOrbTot)
      SxArray<int> aoType;
      /// map fragment orbital to global projector index (:iFragOrbTot)
      SxArray<int> mapFragOrbProj;

      /// size of fragment orbital basis (:iSite)
      SxArray<int> nFragOrb;

      /// site-specific offsets for global fragment basis index (:iSite)
      /// -> iFragOrbTot
      SxArray<int> siteStart;

      /// Atom indices of each site (:iSite)(:iAtom)
      SxArray<SxArray<int> > atomIdx;

      /// Orbital indices (:iSite)(:iFragOrb)
      SxArray<SxArray<SxOrbitalIndex> > orbMap;

      /// map fragment eigenstates to site orbital (:iSite,:iSiteOrb)
      SxArray<int> siteOrbIdx;

      /// Orbital character threshold when selecting by character
      double selectThreshold;
      /// AO type when selecting by character
      int selectType;

      /// rotation from fragment orbital basis to site orbitals
      SxArray<SxVector<SxComplex16> > siteU;

      SxVector<SxComplex16>
      extractFromMatGlobal (ssize_t iSite, const SxVecRef<SxComplex16> &M,
                            const SxVecRef<SxComplex16> &phase) const;

      // Whether to include 3-center PAW corrections
      bool fullPAWNorm;

      // Whether to do global projector normalization
      bool globalProj;

      /// orbitals in radial G basis (:is)(:iol)
      SxAtomicOrbitals refOrbsG,
      /// orbitals projectors radial G basis (:is)(:iol)
                       refOrbsProjG,
      /// PAW projectors in radial G basis (:is)(:ipl)
                       projG;

      /// offset of species in global orbital/projector list ipg=(is,ia,iot)
      SxArray<int> offsetProjGlobal;
      /// map local orbital number iol=(n,l,m) to orbital type (iot)
      SxArray<SxArray<int> > mapType;

      /// Find fragments in the structure
      void findFragments (const SxAtomicStructure &structure);

      /// Get number of sites (fragments)
      int getNSite () const
      {
         SX_CHECK (nFragOrb.getSize () == atomIdx.getSize (),
                   nFragOrb.getSize (), atomIdx.getSize ());
         return int(atomIdx.getSize ());
      }

      SxVector<double> getPAWProjections (ssize_t iSite, const SxPAWPot &pawPot,
                                          const SxAtomicStructure &structure,
                                          const SxAtomicOrbitals &shapeG) const;

      /** \brief Compute the AO basis overlap matrix for one site
          @param iSite       site to compute
          @param pawPot      PAW potentials
          @param structure   current atomic structure
          @param orbitalMap  maps to (is,ia,io) for each orbital (mu)
                             in the global orbital/projector list
        */
      SxVector<double>
      computeSiteOverlap (ssize_t iSite, const SxPAWPot &pawPot,
                          const SxAtomicStructure &structure) const;

      /// Get fragment non-orthogonal Hamiltonian from overlap matrix
      SxVector<double> getHam (ssize_t iSite, const SxVector<double> &siteS) const;

      /// Set up the site rotation matrices
      void setupSiteU (const SxPtr<SxPAWPot> potPtr,
                       const SxAtomicStructure &structure);

      // Get overlap between fragment AO projectors and AO orbitals
      SxVector<double> getSiteT (ssize_t iSite, const SxPAWPot &pawPot,
                                 const SxAtomicStructure &structure);
      /** \brief Projection testing function: project to single site within
                 a global projector basis
          @param psiAO global projections
          @param iSite which site
          @param phase exp(ikR) phases for atoms in this site relative to
                 the representative atom in the structure
          @return projected wave, given in the global basis
        */
      SxVector<SxComplex16>
      projectToSite (const SxVecRef<SxComplex16> &psiAO, ssize_t iSite,
                     const SxVecRef<SxComplex16> &phase) const;

   protected:
      /// Caching behavior for global <chi|proj>^-1
      SxAOBasis::Caching cacheGlobalOvlp;
      /// Cache for global chi|proj
      mutable SxArray<SxVecRef<SxComplex16> > cachedOvlp;
   public:
      /// Get global inverse for turning (co-variant) projections into
      /// (contravariant) expansion coefficients
      const SxVecRef<SxComplex16> getInvOvlp (int ik) const;
      /// Set the caching mode for getInvOvlp
      void setCaching (SxAOBasis::Caching mode);

      /// Get k-vector
      Coord getK (int ik) const;

      /** \brief Get k-dependent contribution to density matrix
        */
      SxVector<SxComplex16>
      getProjMatrix (const PsiG &waves, const SxVecRef<double>& focc) const;

};

#endif /* _SX_HUBBARD_FO_H_ */
