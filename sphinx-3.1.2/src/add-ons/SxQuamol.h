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

#ifndef _SX_QUAMOL_H_
#define _SX_QUAMOL_H_

#include <SxUtil.h>
#include <SxVector.h>
#include <SxParser.h>
#include <SxAtomicStructure.h>
#include <SxPseudoPot.h>
#include <SxGkBasis.h>
#include <SxRadBasis.h>
#include <SxRadialBasis.h>
#include <SxAtomicOrbitals.h>
#include <SxQuantumNumbers.h>
#include <SxConstants.h>
#include <SxTimer.h>
#include <SxError.h>
#include <SxCLI.h>
#include <SxFermi.h>
#include <SxPW.h>
#include <SxYlm.h>
#include <SxPWOverlap.h>
#include <SxPAWOverlap.h>
#include <SxExt.h>

typedef SxVector<SxComplex16> SxOrbitals; //(iOrbital:ig)
typedef SxVector<SxComplex16> SxDMatC16; //(iOrbital:ig)
typedef SxVector<SxComplex16> SxWaveState; //(ig)

class SX_EXPORT_EXT SxQuamol   {
   public:
      /// \brief Standart Constructor
      SxQuamol ();
      /// \brief Standart Destructor
      ~SxQuamol ();
      /// \brief set Function
      void set ( SxAtomicOrbitals &radialsIn,
                 SxConstPtr<SxRadialBasis> radGBasisPtr,
                 SxConstPtr<SxPW> wavePtrIn,
                 SxConstPtr<SxFermi> fermiPtrIn,
                 SxConstPtr<SxOverlapBase> SPtrIn);

      // Variables
      /// \brief QUasi AtoMic OrbitaLS (QUAMOLS) in radGBasis representation
      SxAtomicOrbitals radials;

      double functionalValue;
      SxAtomicOrbitals grad;

      double spillage;
      SxAtomicOrbitals spillageGrad;
      double sigma;

      double eKin;
      SxAtomicOrbitals eKinGrad;
      double zeta;
      bool adaptiveZeta;

      double locVal;
      SxAtomicOrbitals locGrad;
      double kappa;
      double rStart;
      bool adaptiveKappa;

      /// \brief maximal radial extension from Radialbasis of initial guess
      double rMax;
      void computeLocRGGprime (SxConstPtr<SxRadialBasis> radRBasisPtr);
      SxList<SxVector<double> > locRGGPrime;

      /// \brief Waves
      SxConstPtr<SxPW> wavesPtr;
      /// \brief Fermi
      SxConstPtr<SxFermi> fermiPtr;
      /// \brief Overlap Operator
      SxConstPtr<SxOverlapBase> SPtr;
      /// \brief Convergence criterium for functional change
      double dF;
      /// \brief Convergence criterium for gradient length
      double dRes;
      /// \brief Maximum number of steps for spillage optimization
      int maxSteps;
      double relLineMin;
      /** \brief Flag to print radials gradients and directions 
           in every optimization step
      **/
      bool print;

      void computeFunctional (const SxAtomicOrbitals &functions,
                              bool calcGrad = true, bool calcKinLoc = true);
      void computeSpillage (const SxAtomicOrbitals &functions,
                            bool calcGrad = true);
      void computeKineticEnergy (const SxAtomicOrbitals &functions,
                            bool calcGrad = true);
      void computeLocalization (const SxAtomicOrbitals &functions,
                            bool calcGrad = true);
      SxAtomicOrbitals calcNumericGrad (const SxAtomicOrbitals &functions,
            double h = 1e-6);
      double lineMin (const SxAtomicOrbitals &dir,
                      const SxVector3<double> &x,
                      SxVector3<double> &y);
      SxComplex<double> parabelFit (double x0, double x1, double x2,
                                        double N0, double N1, double N2);
      /// \brief Expand radialsG into Gk Basis
      SxOrbitals expandRadialsG (const SxAtomicOrbitals &functions, int ik);
      SxAtomicOrbitals getOrbitals (const SxAtomicOrbitals &functions,
            SxConstPtr<SxRadBasis> radBasisPtr);
      SxAtomicOrbitals getOrbitals (const SxAtomicOrbitals &functions,
            SxConstPtr<SxRadialBasis> radRBasisPtr);

      void checkTrafo (const SxAtomicOrbitals &functionsIn);

      void refine (SxArray<int> &factor);

      SxAtomicOrbitals setOrthogonal (
            const SxAtomicOrbitals &functions,
            const SxAtomicOrbitals &referenz);

      SxAtomicOrbitals setOrthogonal (
            const SxAtomicOrbitals &functions);

      void completenessProfile(const SxAtomicOrbitals &functions);
      void getNormContribution(const SxAtomicOrbitals &functions,
            SxConstPtr<SxRadBasis> radBasisPtr);

      /// \brief optimization routine
      void compute ();

      int printStep;

      bool printLine;

      bool checkGrad;

      SxArray<SxArray<SxVector<double> > > dRgdRGk;

      void setDRDR (int ik);

      SxArray<SxArray<bool> > fixedOrbitals;

      void setFixedList(const SxSymbolTable *table);

      void printSpillagePerState (const SxAtomicOrbitals &functions);

      void updateGuess (SxVector3<double> &x, SxVector3<double> &y, double val);
      
      void calcResidues (const SxAtomicOrbitals &functions);

};

namespace Timer {
   enum QuamolTimer {
      setup,
      computeFunctional,
      computeSpillage,
      computeKineticEnergy,
      computeLocalization,
      calcNumericGrad,
      stepClock,
      lineMin,
      radG2mu,
      checkTrafo,
      getROrbitals,
      DRDataDRGrid
   };
}

SX_REGISTER_TIMERS(Timer::QuamolTimer)
{
   using namespace Timer;
   regTimer (setup,                "Setup");
   regTimer (computeFunctional,    "compute F");
   regTimer (computeSpillage,      "compute Spillage");
   regTimer (computeKineticEnergy, "compute Ekin");
   regTimer (computeLocalization,  "compute Loc");
   regTimer (calcNumericGrad,      "numeric gradient");
   regTimer (stepClock,            "Step Clock");
   regTimer (lineMin,              "Line Minimization");
   regTimer (radG2mu,              "RadG -> orbitalsG");
   regTimer (checkTrafo,           "trafo check");
   regTimer (getROrbitals,         "get orbitals in R");
   regTimer (DRDataDRGrid,         "DRData-Grid");
}

#endif /* _SX_QUAMOL_H_ */
