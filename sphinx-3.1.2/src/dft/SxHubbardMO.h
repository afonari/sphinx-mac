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

#ifndef _SX_HUBBARD_MO_H_
#define _SX_HUBBARD_MO_H_

#include <SxDFT.h>
#include <SxNaturalCubicSpline.h>
#include <SxPAWPot.h>
#include <SxPartialWaveBasis.h>
#include <SxRadialBasis.h>
#include <SxFermi.h>
#include <SxBlockDensityMatrix.h>
#include <SxHubbardXO.h>

class SxHubbardU;

/** \brief MO projectors for DFT + U

    \author Christoph Freysoldt, freysoldt@mpie.de */
class SX_EXPORT_DFT SxHubbardMO : public SxHubbardXO {
   protected:
      /// dummy atomic structure for setupG (contains setup cell)
      SxAtomicStructure setupStructure;
      /// G-basis for MO setup
      SxGBasis setupG;

      /// PAW partial wave basis (with projection to setupG)
      SxPtr<SxPartialWaveBasis> pBasis;

      /// Full AO for MO normalization
      SxVector<double> shapeAO;

   protected:
      /// AO projector for site projection
      SxVector<double> shapeProj;

      /// Atom indices, molecule by molecule
      SxArray<int> atomIdx;

   public:
      /// Constructor
      SxHubbardMO (int siteOffsetIn = 0);

      /// Destructor
      virtual ~SxHubbardMO () = default;

   protected:
      /** \brief Setup up atomic orbital and Hubbard projector
          @param ao  Atomic orbital in radial real space
          @param rCut cutoff radius for the truncated projector
          @param truncWidth transition for cutoff
          @param pawPotPtr PAW potential
          @param eCut plane-wave energy cutoff

          The species identity (is) and l quantum number is taken from ao.
          The m value (>=0) of ao determines the rotational quantum number
          of the MO constructed from the ao.
        */
      void setupAO (const SxVecRef<double> &ao,
                    double rCut, double truncWidth,
                    const SxPtr<SxPAWPot> &pawPotPtr,
                    double eCut);
      /** \brief Set up the (simple cubic) setup box and its G basis

          @param rSize cubic lattice constant
          @param eCut  plane-wave cutoff
        */
      void setupBox (double rSize, double eCut);

      /** \brief Find the molecules according to MO Group in input file
          @param moGroup HubbardU.MO group
          @param structure structure
          @param speciesInfo needed for chemical symbols, use PAW potential
          @return species id
          */
      int findMolecules (const SxSymbolTable* moGroup,
                         const SxAtomicStructure &structure,
                         const SxSpeciesData &speciesInfo);
      /** \brief Find the atoms according to AO Group in input file
          @param aoGroup HubbardU.AO group
          @param structure structure
          @param speciesInfo needed for chemical symbols, use PAW potential
          @return species id
          */
      int findAtoms(const SxSymbolTable* aoGroup,
                    const SxAtomicStructure &structure,
                    const SxSpeciesData &speciesInfo);
   public:
      /// Read and setup for homonuclear diatomics
      virtual void read (const SxSymbolTable *table,
                 const SxAtomicStructure &structure,
                 const SxPtr<SxPAWPot> &potPtr) override;
      /// Read and setup for single atoms
      void readAO (const SxSymbolTable *table,
                   const SxAtomicStructure &structure,
                   const SxPtr<SxPAWPot> &potPtr);
      /// Remove intermediate data used for setup
      virtual void finalize () override;
   protected:
      /// Minimum distance
      double distFrom;
      /// dr for pNorm interpolation
      double dDist;
      /// Constituent AO's l value
      int l;
      /// MO rotational momentum
      int mMO;
      /// Phase between the atoms
      double sign;

      /// Number of atomic orbitals per site
      int nAoPerSite;

   public:
      /// Return number of atomic orbitals for site iSite
      virtual int getNOrbSite (ssize_t /* iSite */) const override
      {
         return nAoPerSite;
      }
      /// Number of interpolation points
      int nInterpolate;
      /// Get l-value of orbitals
      int getL () const { return l; }
   protected:
      /// MO projector normalization as function of distance
      SxNaturalCubicSpline pNorm;

      /// Get projector norm
      double getPNorm (double dist)  {
         //return 0.7;
         return pNorm.getValYExtra ((dist - distFrom) / dDist);
      }

      /// Get derivative of projector norm
      double getPNormDeriv (double dist)  {
         //return 0.;
         return pNorm.getDerivYExtra ((dist - distFrom) / dDist) / dDist;
      }

      /** \brief Orbital rotation from the standard z-axis orientation
          to a given molecular axis.

          @param axis the molecular axis (interatomic vector)
          @return a 2x(2L+1) matrix, containing the rotation of the
                  atomic orbital projectors (2L+1 many) to the -M and
                  +M site. The value of M (rotational quantum number)
                  is taken from shapeProj, which gets its value from
                  the #setup routine.
          */
      SxVector<double> getRot (const Coord &axis) const;

      /// Validate distance
      bool validateDist (double dist) const;

      /** Set up normalization table
          The Hubbard MO projector must be normalized such that
          <MO | P_mm | MO > = 1
          if MO is the normalized orbital for the site index m.


          Since P_mm' = |p_m> N <p_m'|, where |p_m> is a
          modified orbital (shapeProj), the normalization constant
          is
          N = <MO|\hat S|MO>/<MO|p_m> <p_m|MO>

          Note that our shapeProj is set up from two truncated atomic
          orbitals that include the PAW overlap operator for the
          atom they are associated with. In contrast, the MO is set up
          from the full atomic orbital (shapeAO) and uses the
          overlap operator for both atoms.

          For homonuclear diatomic molecules, N depends only on the
          interatomic distance, since the shape perpendicular to the axis
          is dictated by symmetry (the rotational quantum number).

          The number of points is nInterpolate, which defaults to 100.

        */
      void setupNormalization (double distFromIn,
                               double distTo);

   public:
      /// Setup AO projectors in G+k space
      virtual void setupProjGk (const SxGkBasis &gk,
                                const SxPtr<SxPartialWaveBasis>&) override;

      /** \brief Compute the Hubbard energy and Hamiltonian
        @param hubbardU the Hubbard parent object containing the U parameters,
                        energy, etc.
        @param Pij the AO block density matrix
        @param structure the atomic structure
        */
      virtual void compute (SxHubbardU *hubbardU,
                    const SxBlockDensityMatrix& Pij,
                    const SxAtomicStructure &structure) override;
   protected:

      /** \brief Compute the Hubbard energy and Hamiltonian for homonuclear
                 diatomics
        @param hubbardU the Hubbard parent object containing the U parameters,
                        energy, etc.
        @param Pij the AO block density matrix
        @param structure the atomic structure
        */
      void computeMO (SxHubbardU *hubbardU,
                      const SxBlockDensityMatrix& Pij,
                      const SxAtomicStructure &structure);
      /** \brief Compute the Hubbard energy and Hamiltonian for single atoms
        @param hubbardU the Hubbard parent object containing the U parameters,
                        energy, etc.
        @param Pij the AO block density matrix
        @param structure the atomic structure
        */
      void computeAO (SxHubbardU *hubbardU,
                      const SxBlockDensityMatrix& Pij,
                      const SxAtomicStructure &structure);
   public:
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

      class SymMap;

      /// Number of atoms for a given site
      virtual int getNAtoms (ssize_t iSite) const override
      {
         return (latShift.getSize () > 0) ? 2 : 1;
      }

      /** \brief Get global atom id for a given atom
          @param iSite  which site
          @param ia     atom id within that site
          \note This must be consistent with atom ordering in trafoForces
         */
      virtual int getIAtom (ssize_t iSite, ssize_t ia) const override
      {
         return int (latShift.getSize () > 0 ? (2 * iSite + ia) : ia);
      }

      /** \brief Number of orbitals associated with a specific atom
          @param iSite    which site
          @param iAtom    atom id within that site
        */
      virtual int getNOrbAtom (ssize_t iSite, ssize_t ia) const override
      {
         return 2 * l + 1;
      }

};

#endif /* _SX_HUBBARD_MO_H_ */
