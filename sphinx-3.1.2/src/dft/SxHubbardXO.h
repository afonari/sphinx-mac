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

#ifndef _SX_HUBBARD_XO_H_
#define _SX_HUBBARD_XO_H_

#include <SxDFT.h>
#include <SxNaturalCubicSpline.h>
#include <SxPAWPot.h>
#include <SxPartialWaveBasis.h>
#include <SxRadialBasis.h>
#include <SxFermi.h>
#include <SxBlockDensityMatrix.h>

class SxHubbardU;

/** \brief XO projectors for DFT + U

    \author Christoph Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DFT SxHubbardXO {

   public:
      /// Mapping class for symmetrizing block density matrix
      class SymMap
      {
         public:
            /// Get rotated site id upon rotation iSym
            virtual int rotateSite (int iSite, int iSym) = 0;
            /// Get rotated orbital type of site iSite upon rotation iSym
            virtual int rotatedType (int iSite, int iOrbType, int iSym) = 0;
            /// Get number of orbital types (across all atoms) for one site
            virtual int getNOrbType (int iSite) = 0;
            /// Get the reference site (smallest equivalent site id)
            virtual int getRefSite (int iSite) = 0;
            /// Get orbital offset for particular type at particular site
            virtual int getOffset (int iSite, int iOrbType) = 0;
            /// Get l-channel for an orbital type
            virtual int getL (int iSite, int iOrbType) = 0;

            /// Virtual destructor
            virtual ~SymMap () = default;
      };

      /// Radial G basis for atomic-orbital setup
      SxPtr<SxRadialBasis> radGPtr;

   //protected:
      /// Lattice shifts for molecules across cell boundary
      SxArray<Coord> latShift;

   public:
      /// Atomic orbital projectors
      SxPtr<SxAOBasis> aoProj;

      /** \brief Hubbard Hamiltonian in the AO projector space :iSpin, :iSite
        */
      SxBlockDensityMatrix hamProj;

      /// Force contribution from AO->XO transformation
      SxAtomicStructure trafoForce;

      /// Number of radial points
      int nRad;
      /// Produce debug output
      bool verbose;
      /// Prefix for verbose output files
      SxString prefix;

   protected:
      /// Number of sites
      int nSite;

   public:
      int getNSite () const { return nSite; }

      /// Offset in the complete list of Hubbard U sites
      int siteOffset;

      /// Constructor
      SxHubbardXO (int siteOffsetIn = 0)
         : nRad(200), verbose(false), nSite(0), siteOffset(siteOffsetIn)
      {
         // empty
      }

      /// Destructor
      virtual ~SxHubbardXO () = default;

      // --- abstract base class interface
      /// Read and setup for homonuclear diatomics
      virtual void read (const SxSymbolTable *table,
                         const SxAtomicStructure &structure,
                         const SxPtr<SxPAWPot> &potPtr) = 0;
      /// Setup AO projectors in G+k space
      virtual void setupProjGk (const SxGkBasis &gk,
                                const SxPtr<SxPartialWaveBasis> &pBasis) = 0;

      /// Remove intermediate data used for setup
      virtual void finalize () { /* empty */ }

      /// Return number of atomic orbitals for site iSite
      virtual int getNOrbSite (ssize_t iSite) const = 0;

      /// Return number of correlated orbitals for site iSite
      virtual int getNOrbHam (ssize_t iSite) const
      {
         return getNOrbSite (iSite);
      }

      /// Number of atoms for a given site
      virtual int getNAtoms (ssize_t iSite) const = 0;

      /** \brief Get global atom id for a given atom
          @param iSite   which site
          @param iAtom   atom id within that site
          \note This must be consistent with atom ordering in trafoForces
         */
      virtual int getIAtom (ssize_t iSite, ssize_t iAtom) const = 0;

      /** \brief Number of orbitals associated with a specific atom
          @param iSite    which site
          @param iAtom    atom id within that site
        */
      virtual int getNOrbAtom (ssize_t iSite, ssize_t iAtom) const = 0;

      /** \brief Apply the +U Hamiltonian
          \note The return vector must be projectable to the G+k basis.
        */
      virtual PsiG apply (const SxVecRef<SxComplex16> &psi) const
      {
         return hamProj ^ ( *aoProj | psi );
      }

      /** \brief Compute the Hubbard energy and Hamiltonian
        @param hubbardU the Hubbard parent object containing the U parameters,
                        energy, etc.
        @param Pij the AO block density matrix
        @param structure the atomic structure
        */
      virtual void compute (SxHubbardU *hubbardU,
                            const SxBlockDensityMatrix& Pij,
                            const SxAtomicStructure &structure) = 0;

      /** \brief Add contribution of one k-point to the AO block density matrix
        @param Pij The block density matrix to be computed
        @param waves wave functions for current k-point
        @param weight k-point weight
        @param focc occupation number for each state
        */
      virtual void addToRho (SxBlockDensityMatrix *Pij,
                             const SxVecRef<PrecCoeffG> &waves,
                             double weight,
                             const SxVecRef<double> &focc) const = 0;

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
                      const SxArray<SxVector<SxComplex16> > &gradP) const = 0;

      /// Get symmetry map
      virtual SxPtr<SymMap>
      getSymMap (const SxAtomicStructure &structure) const = 0;

      // --- end of abstract base class interface

      /// Symmetrize the AO block density matrix
      void symmetrize (SxBlockDensityMatrix *Pij,
                       const SxAtomicStructure &structure,
                       const SxYlmRotGroup &ylmRot) const;
      /** \brief symmetrize the gradient of the AO block density matrix

          @param Pij       pointer to the AO block density matrix gradient
          @param structure atomic structure
          @param ylmRot    symmetries in spherical harmonics basis
        */
      void symmetrizeGradPij (SxBlockDensityMatrix *Pij,
                              const SxAtomicStructure &structure,
                              const SxYlmRotGroup &ylmRot) const;

      /** \brief Get forces from the displacement of the AO projectors

          @param gradPij gradient of block density matrix with respect to
                         displacement of the left-hand atom
          @param iSpin   spin-channel this block density matrix refers to
          @return The corresponding forces.
      */
      SxAtomicStructure
      getForce (const SxBlockDensityMatrix &gradPij, ssize_t iSpin) const;

      // --- auxiliary function for setup
      /// Real-space truncate a shape provided in radial G space
      SxVector<double> truncateShapeG (const SxVecRef<double> &shape,
                                       double rCut,
                                       double width) const;
};

#endif /* _SX_HUBBARD_XO_H_ */
