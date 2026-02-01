#include <SxDFTConfig.h>
#include <SxAOBasis.h>
#include <SxProjMatrix.h>
#include <SxPWOverlap.h>
#include <SxRadialBasis.h>
#include <SxRadBasis.h>
#ifdef USE_SXGEMMM
#include <SxGemmm.h>
#endif

class SxAOBasisProj
   : public SxProjMatrix<PrecCoeffG>::SaveProjections
{
   public:
      const SxArray<SxOrbitalIndex> &orbitalMap;
      const SxArray<PsiG> &refOrbitals;
      mutable SxVecRef<PrecCoeffG> phase;
      SxAOBasisProj(const SxAOBasis &ao, int ik)
         : SxProjMatrix<PrecCoeffG>::SaveProjections (
               ao.getNOrb (),
               ao.blockSize,
               (int)ao.refOrbitals(ik)(ao.orbitalMap(0).is).getNRows ()),
           orbitalMap(ao.orbitalMap),
           refOrbitals(ao.refOrbitals(ik))
      { /* empty */ }

      virtual ~SxAOBasisProj () {}

      virtual void getProjector(int iOrb, SxVecRef<PrecCoeffG> *target) const
      {
         SX_CHECK (iOrb >= 0 && iOrb < orbitalMap.getSize (),
                   iOrb, orbitalMap.getSize ());
         const SxOrbitalIndex &idx = orbitalMap(iOrb);
         const SxVecRef<PrecCoeffG> &refOrb
            = refOrbitals(idx.is).colRef (idx.io);

         if (phase.getSize () <= 0
             || phase.auxData.is != idx.is
             || phase.auxData.ia != idx.ia)
         {
            phase.unref ();
            phase = refOrb.getBasis<SxGBasis> ().getPhaseFactors (idx.is, idx.ia);
         }
         SX_CHECK (phase.auxData.is == idx.is,
                   phase.auxData.is, idx.is);
         SX_CHECK (phase.auxData.ia == idx.ia,
                   phase.auxData.ia, idx.ia);

#ifdef USE_OPENMP
#pragma omp parallel for if (nElements > sxChunkSize)
#endif
         for (int ig = 0; ig < nElements; ++ig)
            (*target)(ig) = phase(ig) * refOrb(ig);
      }

      HAS_TARGET_GETPROJECTOR;
      NO_GETFACTOR;
};

SxAOBasis::SxAOBasis ()
{
   cacheRefOrb = Unknown;
   init ();
}

namespace Timer {
   // timer ids
   enum AOTimers{ AoRefOrbitalSetup, AoOverlap, AoSInversion, AoProjection,
                  AoGradient, AoGradProject, AoTotalTime} ;
}

SX_REGISTER_TIMERS(Timer::AOTimers)
{
   using namespace Timer;
   regTimer (AoRefOrbitalSetup,"Ref. orbital setup");
   regTimer (AoOverlap,"Overlap setup");
   regTimer (AoSInversion,"Overlap inversion");
   regTimer (AoProjection,"AO projection");
   regTimer (AoGradient,"AO gradient");
   regTimer (AoGradProject,"AO projection d/dR");
   regTimer (AoTotalTime, "AOBasis total");
}

void SxAOBasis::init ()  {
   cacheOverlap = cacheInverse = Unknown;
   refOrbCachedK = overlapCachedK = invOverlapCachedK = -1;
   blockSize = 64;
}

SxAOBasis::~SxAOBasis ()
{
   deregisterAll ();
}

SxAOBasis::SxAOBasis (const SxGkBasis &gk,
                      const SxArray<SxArray<SxVector<double> > > &psiRad,
                      const SxConstPtr<SxOverlapBase> SPtrIn)
   : SPtr (SPtrIn)
{
   init ();
   set (gk, wrapContainer (psiRad) );
}

SxAOBasis::SxAOBasis (const SxVecRef<SxComplex16> &YlmGl, int iSpecies)
{
   init ();
   const SxAtomicStructure &structure = YlmGl.getBasis<SxGBasis> ().getTau ();
   int nRef = (int)YlmGl.getNCols ();
   int nOrb = nRef * structure.getNAtoms (iSpecies);

   orbitalMap.resize (nOrb);
   for (int ia = 0; ia < structure.getNAtoms (iSpecies); ++ia)  {
      for (int io = 0; io < nRef; ++io)  {
         orbitalMap(io + nRef * ia).set (iSpecies, ia, io);
      }
   }

   refOrbitals.resize (1);
   refOrbitals(0).resize (structure.getNSpecies ());
   refOrbitals(0)(iSpecies) = YlmGl;
   cacheRefOrb = CacheAll;

}

void SxAOBasis::set (const SxGkBasis &gk,
                     const SxRefAOContainerBase &psi)
{
   SX_CLOCK (Timer::AoTotalTime);
   if (!SPtr) SPtr = SxPtr<SxPWOverlap>::create ();
   int nSpecies = psi.getNSpecies ();
   SX_CHECK (gk.getNk () > 0, gk.getNk ());
   SX_CHECK (nSpecies > 0, nSpecies);
   cacheRefOrb = CacheAll;
   const SxVector<int> &nAtoms = gk.getTau ().atomInfo->nAtoms;

   // --- setup orbital map io -> (n,l,m) for each species
   refOrbMap.resize (nSpecies);

   int nOrb = 0;
   for (int is = 0; is < nSpecies; ++is)  {
      int nOrbType = psi.getNOrb (is);
      int nOrbSpecies = 0;
      for (int iOrbType = 0; iOrbType < nOrbType; ++iOrbType)
         nOrbSpecies += 2 * psi.getL(is,iOrbType) + 1;
      refOrbMap(is).resize (nOrbSpecies);
      int io = 0;
      for (int iOrbType = 0; iOrbType < nOrbType; ++iOrbType) {
         int l = psi.getL(is,iOrbType);
         for (int m = -l; m <= l; ++m, ++io)  {
            refOrbMap(is)(io).set (iOrbType, l, m);
         }
      }
      nOrb += io * nAtoms(is);
   }

   // -- setup orbital map iOrb -> (is, ia, io)
   orbitalMap.resize (nOrb);
   for (int is = 0, iOrb = 0; is < nSpecies; ++is)
      for (int ia = 0; ia < nAtoms(is); ++ia)
         for (int io = 0; io < refOrbMap(is).getSize (); ++io, ++iOrb)
            orbitalMap(iOrb).set (is, ia, io);

   // --- setup reference orbitals
   SX_CLOCK (Timer::AoRefOrbitalSetup);
   ssize_t nk = gk.getNk ();
   refOrbitals.resize (nk);

   double gMax = 0.;
   for (int ik = 0; ik < nk; ++ik)  {
      SX_MPI_LEVEL("waves-k");
      if (SxLoopMPI::myWork(ik))  {
         refOrbitals(ik).resize(nSpecies);
         for (int is = 0; is < nSpecies; ++is)  {
            ssize_t nOrbital = refOrbMap(is).getSize ();
            if (nOrbital > 0)  {
               refOrbitals(ik)(is).reformat (gk(ik).ng, nOrbital);
               refOrbitals(ik)(is).setBasis (&gk(ik));
            }
         }
         gMax = max(gMax, sqrt(gk(ik).g2(gk(ik).ng - 1)));
      }
   }
   if (gMax == 0.) return; // no k-points

   gMax *= 1.1;
   double cutVol = FOUR_PI / 3. * gMax * gMax * gMax;
   double dg = 0.003;
   int ng1 = int(cutVol / gk.getTau ().cell.getReciprocalCell ().volume);
   int ng2 = int(gMax/dg) + 1;
   SxRadialBasis *radG = NULL;
   if (ng1 * nk > ng2) {
      radG = new SxRadialBasis (0., gMax, ng2, SxRadialBasis::Linear);
      radG->realSpace = false;
   }

   SX_LOOP(is)  {
      if (refOrbMap(is).getSize () == 0) continue;
      // reference orbital with basis used for the |G+k> projection below
      SxVecRef<double> proj;
      proj.auxData.n = -1;
      SX_LOOP(io) {
         const SxAoIndex &nlm = refOrbMap(is)(io);
         // handle to reference orbital with basis as provided by container
         SxVecRef<double> psiRef = psi.getOrb(is, nlm.n);
         psiRef.auxData.is = int (is);
         psiRef.auxData.setOrb (nlm.n, nlm.l, nlm.m);
         // try to get a radial G basis
         const SxRadialBasis *myRadG
            = dynamic_cast<const SxRadialBasis*>(psiRef.getBasisPtr ());
         if (!myRadG) myRadG = radG;

         // now set proj
         if (myRadG) {
            if (proj.auxData.n != nlm.n)  {
               // set up spline for each new ref
               proj.unref ();
               proj = myRadG->toSpline (*myRadG | psiRef);
            }
         } else {
            // use psiRef as is
            proj.unref ();
            proj = psiRef;
         }
         proj.auxData.is = (int)is;
         proj.auxData.setOrb (nlm.n, nlm.l, nlm.m);

         // --- now compute for all |G+k> as needed
         for (int ik = 0; ik < nk; ++ik)  {
            if (refOrbitals(ik).getSize () > 0)
               refOrbitals(ik)(is).colRef(io) <<= gk(ik) | proj;
         }
      }

   }
   // clean overlap cache
   overlap.resize (0);
   invOverlap.resize (0);
   overlapCachedK = invOverlapCachedK = -1;
}

namespace {
   class SxMatrixAO : public SxRefAOContainerBase
   {
      private:
         const SxArray<SxVector<double> > &psi;
         const SxArray<SxArray<int> > &lPhi;
      public:
         SxMatrixAO (const SxArray<SxVector<double> >&psiIn,
                     const SxArray<SxArray<int> > &lIn)
            : psi(psiIn), lPhi(lIn)
         { /* empty */ }
         /// Destructor
         virtual ~SxMatrixAO () {}

         /// Get reference orbital
         SxVecRef<double>
         getOrb (ssize_t iSpecies, ssize_t iRef) const override
         {
            return const_cast<SxVector<double>&>(psi(iSpecies)).colRef(iRef);
         }
         /// Get l of reference orbital
         int getL (ssize_t iSpecies, ssize_t iRef) const override
         {
            return lPhi(iSpecies)(iRef);
         }
         /// Get number of species
         int getNSpecies () const override
         {
            return int(psi.getSize ());
         }
         /// Get number of reference orbitals
         int getNOrb (ssize_t iSpecies) const override
         {
            return int(psi(iSpecies).getNCols ());
         }
   };
}

SxAOBasis::SxAOBasis (const SxGkBasis &gk,
                      const SxArray<SxVector<double> > &psi,
                      const SxArray<SxArray<int> > &lPsi)
{
   init ();
   set (gk, SxMatrixAO (psi, lPsi));
}

SxVector<SxComplex16> SxAOBasis::calculateOverlap (int ik) const
{
#ifndef NDEBUG
   if (cacheRefOrb == CacheAll)  {
      int nk = int(refOrbitals.getSize ());
      SX_CHECK (nk > 0);
      SX_CHECK (ik >= 0 && ik < nk, ik, nk);
   } else if (cacheRefOrb == CacheCurrentK)  {
      SX_CHECK (refOrbCachedK == ik, refOrbCachedK, ik);
   }
#endif
   int nOrb = int(orbitalMap.getSize ());

   SX_CLOCK (Timer::AoTotalTime);
   SX_CLOCK (Timer::AoOverlap);
   SxVector<SxComplex16> ovlp;
   SxOverlap S(SPtr);
   const SxGBasis *gBasis = dynamic_cast<const SxGBasis*>
                           (refOrbitals(ik)(0).getBasisPtr ());
   SX_CHECK (gBasis);
   int ng = gBasis->ng;
   // --- high-mem
   /*
   //SxAOBasisProj aop(*this, ik);
   TPsi aoMu.reformat(ng, nOrb);
   for (iOrb = 0; iOrb < nOrb; ++iOrb)
      aoMu.colRef(iOrb) <<= S | getAOinG(ik, iOrb);
   ovlp = aop.getProjection(aoMu);
   */

   // ---medium-mem
   ovlp.reformat(nOrb, nOrb);
   int nblock = 64;
   TPsi muBlock;
   for (int iOrb = 0, ib = 0; iOrb < nOrb; /* inside */)  {
      if (ib == 0)  {
         if (iOrb + nblock > nOrb) nblock = nOrb - iOrb;
         muBlock.reformat (ng, nblock);
         muBlock.setBasis (gBasis);
         muBlock.auxData.ik = ik;
      }
      muBlock.colRef (ib) <<= getAOinG(ik, iOrb + ib);
      if (++ib == nblock)  {
         SxIdx idx(iOrb * nOrb, (iOrb + nblock) * nOrb - 1);
         ovlp(idx) <<= fromPWBasis (muBlock);
         iOrb += nblock;
         ib = 0;
      }
   }

   return ovlp;
}


SxVecRef<SxAOBasis::CoeffType>
SxAOBasis::identity  (const SxAOBasis *basis,
                      const SxVecRef<CoeffType> &psiAO) const
{
   SX_CHECK (basis == this);
   return const_cast<SxVecRef<CoeffType>&>(psiAO);
}

/** \brief Multiply a row-block of psiAO with an orbital-dependent extra phase

  res = psiAO(ipl,i) * phase(ipl)

  for ipl = iOrbStart .. iOrbStart+nOrbBlock-1

  This is an auxiliary function for toPWBasis.
  */
SxVector<SxComplex16>
getAOExtraPhaseBlock(const SxVecRef<SxComplex16> &psiAO,
                     const SxVecRef<SxComplex16> &phase,
                     ssize_t iOrbStart,
                     ssize_t nOrbBlock)
{
   // --- multiply psiAO by orbital-dependent phase
   ssize_t nPsi = psiAO.getNCols (), nOrb = psiAO.getNRows ();
   SxVector<SxComplex16> psiAOphase (nOrbBlock, nPsi);
   for (ssize_t iPsi = 0, iPh = 0; iPsi < nPsi; ++iPsi)  {
      ssize_t iAO = iOrbStart + iPsi * nOrb;
      for (ssize_t jOrb = 0; jOrb < nOrbBlock; ++jOrb)  {
         psiAOphase(iPh++) = psiAO(iAO++) * phase(iOrbStart + jOrb);
      }
   }
   return psiAOphase;
}


SxVector<PrecCoeffG>
SxAOBasis::toPWBasis (const SxGBasis *gPtr,
                      const SxVecRef<SxComplex16> &psiAO) const
{
   SX_CLOCK (Timer::AoTotalTime);
   SX_CLOCK (Timer::AoGradient);
   // no projection of empty vectors
   SX_CHECK (psiAO.getSize () > 0);

   // --- find k-point
   SX_CHECK (cacheRefOrb == CacheAll);
   int ik = psiAO.auxData.ik;
   int nk = int(refOrbitals.getSize ());
   if (ik >= 0 && ik < nk)  {
#ifndef NDEBUG
      int is = orbitalMap(0).is; // we might start at is > 0
      SX_CHECK (gPtr == refOrbitals(ik)(is).getBasisPtr ());
#endif
   } else if (cacheRefOrb == CacheAll)  {
      int is = orbitalMap(0).is; // we might start at is > 0
      SX_CHECK (is>=0, is);
      for (ik = 0; ik < nk; ++ik)  {
         if (refOrbitals(ik).getSize () == 0) continue; // different MPI task
         if (refOrbitals(ik)(is).getBasisPtr () == gPtr)
            break;
      }
      if (ik == nk)  {
         cout << "Projection from unregistered G-basis in SxAOBasis\n";
         SX_EXIT;
      }
   }
#ifndef NDEBUG
      if (extraPhase.getSize () > 0)  {
         SX_CHECK (extraPhase(ik).getSize () == getNElements (),
                   extraPhase(ik).getSize (), getNElements ());
      }
#endif

   // now execute the actual implementation
#ifdef USE_SXGEMMM
   return toPWBasisGemm3m(gPtr, ik, psiAO);
#else
   return toPWBasisBLAS(gPtr, ik, psiAO);
#endif

}

#ifdef USE_SXGEMMM
SxVector<PrecCoeffG>
SxAOBasis::toPWBasisGemm3m(const SxGBasis *gPtr,
                           int ik,
                           const SxVecRef<SxComplex16> &psiAO) const
{
   /* --- Direct algorithm (sxgemm3m algorithm)
      sxgemm3m is a direct blocked algorithm that calculates
      p(ig,ipl) * T(ig,ia) * psiAO(ipl, ia, iPsi)
      exploiting CPU vectorization and level caches.
      Its performance appears superior on modern CPUs.
   */
   int nPsi = (int)psiAO.getNCols ();
   int nOrb = getNOrb ();
   const SxArray<PsiG> &proj = refOrbitals(ik);

   SxVector<PrecCoeffG> res(gPtr->ng, nPsi);
   SxVector<PrecCoeffG> T;
   res.set (0.);
   // --- loop over species
   for (int iOrb = 0; iOrb < nOrb; /* empty */) {
      int is = orbitalMap(iOrb).is;

      // number of projectors for this species
      int npl = (int)proj(is).getNCols ();

      // --- get number of atoms
      int nAtoms, iOrbStart = iOrb;
      for (nAtoms = 0; iOrb < nOrb; iOrb += npl, ++nAtoms)  {
         if (orbitalMap(iOrb).is != is) break;
         SX_CHECK (orbitalMap(iOrb).io == orbitalMap(iOrbStart).io);
      }

      SxVector<double> projReal = proj(is).real ();
      int nAtomBlock = -1;
      // --- loop over atom blocks (to limit memory for T)
      for (int ia = 0; ia < nAtoms; ia += nAtomBlock)  {
         // adjust size of atom block for final block
         nAtomBlock = min(nAtoms-ia, blockSize);

         // --- collect atomic phase factors exp[-i(G+k)tau(ia)]
         T.reformat (gPtr->ng, nAtomBlock);
         for (int iab = 0; iab < nAtomBlock; ++iab)  {
            int iAtom = orbitalMap(iOrbStart + npl * iab).ia;
            SxVecRef<SxComplex16> Tcol = T.colRef (iab);
            gPtr->getPhaseFactors(is, iAtom, Tcol);
         }

         // --- now do summation over ipl and ia
         if (extraPhase.getSize () > 0)  {
            // multiply psiAO by orbital-dependent phase
            SxVector<SxComplex16> psiAOphase
               = getAOExtraPhaseBlock(psiAO, extraPhase(ik),
                                      iOrbStart, nAtomBlock * npl);
            // do sxgemm summation
            sxpgemm3m(nPsi, nAtomBlock, npl, gPtr->ng,
                      T.elements, projReal.elements,
                      psiAOphase.elements, npl, nAtomBlock * npl,
                      res.elements);
         } else {
            sxpgemm3m(nPsi, nAtomBlock, npl, gPtr->ng,
                      T.elements, projReal.elements,
                      psiAO.elements + iOrbStart, npl, nOrb,
                      res.elements);
         }
         iOrbStart += nAtomBlock * npl;
      }
   }

   SX_VALIDATE_VECTOR(res);

   res.auxData = psiAO.auxData;
   res.auxData.ik = ik;
   res.setBasis (gPtr);
   return res;
}

#else

SxVector<PrecCoeffG>
SxAOBasis::toPWBasisBLAS (const SxGBasis *gPtr,
                           int ik,
                           const SxVecRef<SxComplex16> &psiAO) const
{
   /* This is the BLAS algorithm for

      res(ig,i) = sum_ia,ipl proj(ig,ipl) * T(ig, ia) * psiAO(ipl,ia,i)

      ig  is the plane-wave index
      ipl runs over local orbitals located at atom ia
      ia  runs over all atoms
      i   runs over all states

      It collects the projectors and uses zgemm calls to do
      the data reduction.
      sxgemm3m is an alternative direct algorithm, see below.
    */
   const int SX_PWPROJ_BLOCK_MINPSI = 1;
   int nPsi = (int)psiAO.getNCols ();
   if (nPsi > SX_PWPROJ_BLOCK_MINPSI)  {
      // --- use projector blocking
      SxAOBasisProj aop(*this, ik);
      SxVecRef<SxComplex16> psiAOphase;
      if (extraPhase.getSize () > 0)  {
         psiAOphase = getAOExtraPhaseBlock (psiAO, extraPhase(ik),
                                            0, getNElements ());
      } else {
         psiAOphase = const_cast<SxVecRef<SxComplex16>&>(psiAO);
      }
      SxVector<PrecCoeffG> res = aop.gradient (psiAOphase);
      res.auxData = psiAO.auxData;
      res.setBasis (gPtr);
      res.auxData.ik = ik;
      return res;
   } else if (nPsi > 1)  {
      // use algorithm below for each psi
      SxVector<SxComplex16> res (gPtr->ng, nPsi);
      res.auxData = psiAO.auxData;
      res.setBasis (gPtr);
      res.auxData.ik = ik;

      for (int i = 0; i < nPsi; ++i)
         res.colRef(i) <<= toPWBasis (gPtr, psiAO.colRef(i));
      return res;
   }

   // --- fall-back algorithm for a single psi (rarely used)
   /* For a single psi, we follow this strategy:

      For
      \sum_{is,ia,ipl} T(is,ia) |p(is,ipl)>  A(is, ipl, ia)
      we first compute

      Variant 1: (more local orbitals than atoms)
      |pA(is,ia)> = \sum_{ipl} |p(is,ipl)> A(is, ipl, ia)
      as a matrix operation.

      Variant 2: (more atoms than local orbitals)
      |phaseNl(is,ipl)> = \sum_{ia} |T(is,ia)> A(is, ipl, ia)
      as a matrix operation.

      If the number of atoms becomes large (>64), we use blocking
      to limit the size of pA.
   */

   int nOrb = getNOrb ();
   const SxArray<PsiG> &proj = refOrbitals(ik);

   SxVector<PrecCoeffG> res(gPtr->ng);
   res.set (0.);
   // --- loop over species
   for (int iOrb = 0; iOrb < nOrb; /* empty */) {
      int is = orbitalMap(iOrb).is;

      // number of projectors for this species
      int npl = (int)proj(is).getNCols ();

      // --- get number of atoms
      int nAtoms, iOrbStart = iOrb;
      for (nAtoms = 0; iOrb < nOrb; iOrb += npl, ++nAtoms)  {
         if (orbitalMap(iOrb).is != is) break;
         SX_CHECK (orbitalMap(iOrb).io == orbitalMap(iOrbStart).io);
      }

      int nAtomBlock = -1;
      // --- compute gradient (with blocking on atoms)
      if (npl >= nAtoms)  {
         // --- variant 1: sum first over ipl, then over atoms
         // psiAOphase(ipl,a) = psiAO(ipl,a) * phase(ipl,a)
         // [if extraPhase, otherwise phase=1]
         //
         // (1) pA(ig,a) = sum_ipl proj(ig,ipl) psiAOphase(ipl,a)  (BLAS 3)
         // (2) res(ig) = sum_a pA(ig,a) * T(ig,a)                 (vector loop)
         for (int ia = 0; ia < nAtoms; ia += nAtomBlock)  {
            // adjust size of atom block for final block
            nAtomBlock = min(nAtoms-ia, blockSize);
            // --- compute new pA block
            SxVecRef<SxComplex16> psiAOphase;
            if (extraPhase.getSize () > 0)  {
               psiAOphase = getAOExtraPhaseBlock (psiAO, extraPhase(ik),
                                                  iOrbStart, npl * nAtomBlock);
               psiAOphase.reshape (npl, nAtomBlock);
            } else {
               psiAOphase = psiAO.getRef<Compact> (iOrbStart, npl, nAtomBlock);
            }
            SxVector<PrecCoeffG> pA = proj(is) ^ psiAOphase;

            // apply translation for each atom
            for (int iab = 0; iab < nAtomBlock; ++iab)  {
               int iAtom = orbitalMap(iOrbStart + npl * iab).ia;
               SxVecRef<PrecCoeffG> T     = gPtr->getPhaseFactors (is, iAtom),
                                    pAcol = pA.colRef(iab);
#ifdef USE_OPENMP
#pragma omp parallel for simd if (res.getSize () > sxChunkSize)
#endif
               for (int ig = 0; ig < res.getSize (); ++ig)
                  res(ig) += T(ig) * pAcol(ig);
            }
            iOrbStart += nAtomBlock * npl;
         }
         SX_CHECK (iOrbStart  == iOrb, iOrb - iOrbStart, nAtoms, npl);
      } else {
         // --- variant 2: sum first over atoms, then over ipl
         // psiAOphase(ipl,a) = psiAO(ipl,a) * phase(ipl,a)
         // [if extraPhase, otherwise phase=1]
         //
         // (1) phaseNl(ig,ipl) = sum_ia T(ig,a) psiAOphase(ipl,a)  (BLAS 3)
         // (2) res(ig) = sum_ipl proj(ig,ipl) * phaseNl(ig,ipl)  (vector loop)
         SxVector<PrecCoeffG> phaseNl(gPtr->ng, npl);
         phaseNl.set (0.);
         SxVector<PrecCoeffG> T;
         for (int ia = 0; ia < nAtoms; ia += nAtomBlock)  {
            // adjust size of atom block for final block
            nAtomBlock = min(nAtoms-ia, blockSize);

            // --- collect atomic phase factors exp[-i(G+k)tau(ia)]
            T.reformat (gPtr->ng, nAtomBlock);
            for (int iab = 0; iab < nAtomBlock; ++iab)  {
               int iAtom = orbitalMap(iOrbStart + npl * iab).ia;
               SxVecRef<SxComplex16> Tcol = T.colRef (iab);
               gPtr->getPhaseFactors(is, iAtom, Tcol);
            }

            if (extraPhase.getSize () > 0)  {
               const SxVector<SxComplex16> &phase = extraPhase(ik);
               // --- multiply psiAO by orbital-dependent phase
               //     and reorder to (iAtom, ipl)
               // same as getAOExtraPhaseBlock(...).transpose () in one go
               SxVector<SxComplex16> psiAOphase (nAtomBlock, npl);
               for (ssize_t iab = 0; iab < nAtomBlock; ++iab)
                  for (ssize_t ipl = 0; ipl < npl; ++ipl)
                     psiAOphase(iab, ipl) = psiAO(iOrbStart + ipl + iab * npl)
                                          * phase(iOrbStart + ipl + iab * npl);
               phaseNl += T ^ psiAOphase;
            } else {
               // sum over atoms
               phaseNl += T ^ psiAO.getRef<Compact>(iOrbStart,npl,nAtomBlock,1)
                                   .transpose ();
            }
            iOrbStart += nAtomBlock * npl;
         }
         SX_CHECK (iOrbStart  == iOrb,
                   iOrb - iOrbStart, nAtoms, npl);

         // --- sum over local projectors
         int ng = gPtr->ng;
         for (int ipl = 0; ipl < npl; ++ipl)  {
            const SxVecRef<PrecCoeffG> &phiRef = proj(is).colRef(ipl);
#ifdef USE_OPENMP
#pragma omp parallel for if (ng > sxChunkSize)
#endif
            for (int ig = 0; ig < ng; ig++)
               res(ig) += phiRef(ig) * phaseNl(ig + ng * ipl);
         }
      }
   }

   SX_VALIDATE_VECTOR(res);

   res.auxData = psiAO.auxData;
   res.auxData.ik = ik;
   res.setBasis (gPtr);
   return res;
}
#endif


SxAOBasis::TPsi
SxAOBasis::fromPWBasis (const SxVecRef<SxComplex16> &psiG) const
{
   // no projection of empty vectors
   SX_CHECK (psiG.getSize () > 0);

   // --- find k-point
   SX_CHECK (cacheRefOrb == CacheAll);
   int ik = psiG.auxData.ik;
   int nk = int(refOrbitals.getSize ());
   int is = orbitalMap(0).is; // first species, may be > 0
   if (ik >= 0 && ik < nk)  {
      SX_CHECK (psiG.getBasisPtr () == refOrbitals(ik)(is).getBasisPtr ());
   } else if (cacheRefOrb == CacheAll)  {
      SX_CHECK (is>=0, is);
      for (ik = 0; ik < nk; ++ik)  {
         if (refOrbitals(ik)(is).getBasisPtr () == psiG.getBasisPtr ())
            break;
      }
      if (ik == nk)  {
         cout << "Projection from unregistered G-basis in SxAOBasis\n";
         SX_EXIT;
      }
   }
   return SPtr ? fromPWBasis (SPtr->apply(psiG), ik)
               : fromPWBasis (psiG             , ik);
}

void applyPhaseConj(SxVector<SxComplex16> *resPtr,
                    const SxVecRef<SxComplex16> &phase)
{
   SX_CHECK (resPtr->getNRows () == phase.getSize (),
             resPtr->getNRows (), phase.getSize ());
   // index-loop implementation
   SX_LOOP2(iState, iOrb)
      (*resPtr)(iOrb,iState) *= phase(iOrb).conj ();
   /*
   // iterator implementation
   SxVector<SxComplex16>::Iterator resIt = resPtr->begin ();
   for (ssize_t iState = 0; iState < res->getNCols (); ++iState)
      for (auto phaseIt : phase)
         *resIt++ *= (*phaseIt).conj ();
   */
}

#ifdef USE_SXGEMMM
// --- sxgemmm version
SxVector<SxComplex16>
SxAOBasis::fromPWBasis (const SxVecRef<PrecCoeffG> &psi, int ik) const
{
   SX_CLOCK (Timer::AoTotalTime);
   SX_CLOCK (Timer::AoProjection);

   /* --- Direct algorithm (sxgemmm algorithm)
      We also have a direct blocked algorithm that calculates
      p(ig,ipl) * T(ig,ia).conj () * psi(ig, iState)
      exploiting CPU vectorization and level caches.
      Its performance appears superior on modern CPUs.
   */
   const SxArray<PsiG> &proj = refOrbitals(ik);

   int nOrb = getNOrb ();
   ssize_t nStates = psi.getNCols ();
   SxVector<SxComplex16> res(nOrb, nStates);
   // --- loop over species
   for (int iOrb = 0; iOrb < nOrb; /* empty */) {
      int is = orbitalMap(iOrb).is;

      // number of projectors for this species
      int npl = (int)proj(is).getNCols ();

      // --- get number of atoms
      int nAtoms, iOrbStart = iOrb;
      for (nAtoms = 0; iOrb < nOrb; iOrb += npl, ++nAtoms)  {
         if (orbitalMap(iOrb).is != is) break;
         SX_CHECK (orbitalMap(iOrb).io == orbitalMap(iOrbStart).io);
      }

      const SxGBasis &gk = proj(is).getBasis<SxGBasis> ();
      SX_CHECK((psi.getBasisPtr () == &gk)
               || !dynamic_cast<const SxGBasis*>(psi.getBasisPtr ()));

      SX_CHECK (proj(is).imag ().normSqr () < 1e-20);
      SxVector<double> projReal = proj(is).real ();

      // --- perform projections (with blocking on atoms)
      SxVector<PrecCoeffG> conjT;
      int ng = gk.ng;
      int nAtomBlock = -1;
      for (int ia = 0; ia < nAtoms; ia += nAtomBlock)  {
         // adjust size of atom block for final block
         nAtomBlock = min(nAtoms-ia, blockSize);
         conjT.reformat (ng, nAtomBlock);

         // --- collect conjT = T.conj () for all atoms in block
         for (int iab = 0; iab < nAtomBlock; ++iab) {
            // get the real atom index
            int iAtom = orbitalMap(iOrbStart + npl * iab).ia;
            SxVecRef<PrecCoeffG> bcol = conjT.colRef(iab),
                                 T = gk.getPhaseFactors(is, iAtom);
#           ifdef USE_OPENMP
#           pragma omp parallel for if (ng > sxChunkSize)
#           endif
            for (int ig = 0; ig < ng; ig++)
               bcol(ig) = T(ig).conj ();
         }

         // project current block
         sxgemmm(nStates, nAtomBlock, npl, ng,
                 psi.elements, psi.getNRows (), conjT.elements,
                 projReal.elements, res.elements + iOrbStart, npl, nOrb);

         iOrbStart += nAtomBlock * npl;
      }
      SX_CHECK (iOrbStart == iOrb, iOrb - iOrbStart, nAtoms, npl);
   }
   SX_VALIDATE_VECTOR(res);
   res.auxData = psi.auxData;
   res.setBasis (this);
   if (extraPhase.getSize () > 0)
      applyPhaseConj (&res, extraPhase(ik));
   return res;
}

#else

// --- BLAS version
SxVector<SxComplex16>
SxAOBasis::fromPWBasis (const SxVecRef<PrecCoeffG> &psi, int ik) const
{
   SX_CLOCK (Timer::AoTotalTime);
   SX_CLOCK (Timer::AoProjection);

   /* This is the BLAS algorithm
      It collects the projectors and uses zgemm calls to do
      the data reduction.
    */
   const int SX_PWPROJ_BLOCK_MINPSI = 1;
   if (psi.getNCols () > SX_PWPROJ_BLOCK_MINPSI)  {
      // --- use projector blocking
      SxAOBasisProj aop(*this, ik);

      SxVector<SxComplex16> result = aop.getProjectionFromExtended(psi);

      SX_VALIDATE_VECTOR (result);
      result.auxData = psi.auxData;
      result.setBasis (this);
      if (extraPhase.getSize () > 0)
         applyPhaseConj (&result, extraPhase(ik));
      return result;
   } else if (psi.getNCols () > 1)  {
      int nPsi = (int)psi.getNCols ();
      // use algorithm below for each psi
      SxVector<SxComplex16> res (getNOrb (), nPsi);

      for (int i = 0; i < nPsi; ++i)
         res.colRef(i) <<= fromPWBasis (psi.colRef(i));
      res.auxData = psi.auxData;
      res.setBasis (this);
      return res;
   }

   /* --- Use atom-blocking (BLAS algorithm)
      For a single psi, we follow this strategy:

      First note that
      <p(is,ipl) * T(is,ia)|psi> = <p(is,ipl) | T*(is,ia) * psi>

      We now perform the right-hand formula as matrix operation for
      each species:
      <T * p|psi> (npl * nAtoms) = p(npl x ng) ^ psiT(ng x nAtoms)

      If the number of atoms becomes large (>64), we use blocking
      to limit the size of psiT.

   */
   const SxArray<PsiG> &proj = refOrbitals(ik);
   int nOrb = getNOrb ();
   ssize_t nStates = psi.getNCols ();
   SxVector<SxComplex16> res(nOrb * nStates);
   res.reshape (nOrb, nStates);
   // --- loop over species
   for (int iOrb = 0; iOrb < nOrb; /* empty */) {
      int is = orbitalMap(iOrb).is;

      const SxGBasis &gk = proj(is).getBasis<SxGBasis> ();
      SX_CHECK((psi.getBasisPtr () == &gk)
               || !dynamic_cast<const SxGBasis*>(psi.getBasisPtr ()));

      // number of projectors for this species
      int npl = (int)proj(is).getNCols ();

      // --- get number of atoms
      int nAtoms, iOrbStart = iOrb;
      for (nAtoms = 0; iOrb < nOrb; iOrb += npl, ++nAtoms)  {
         if (orbitalMap(iOrb).is != is) break;
         SX_CHECK (orbitalMap(iOrb).io == orbitalMap(iOrbStart).io);
      }

      // --- shortcut for nAtoms = 1
      if (nAtoms == 1)  {
         PsiG psiTrans = gk.getPhaseFactors(is, orbitalMap(iOrbStart).ia).conj()
                       * psi;
         res (SxIdx(iOrbStart,iOrb-1)) <<= proj(is).overlap (psiTrans);
         continue;
      }

      // --- perform projections (with blocking on atoms)
      SxVector<PrecCoeffG> atomBlock;
      int ng = gk.ng;
      int nAtomBlock = -1;
      for (int ia = 0; ia < nAtoms; ia += nAtomBlock)  {
         // adjust size of atom block for final block
         nAtomBlock = min(nAtoms-ia, blockSize);
         atomBlock.reformat (ng, nAtomBlock);
         // collect atomBlock = psi * T.conj () for all atoms in block
         for (int iab = 0; iab < nAtomBlock; ++iab) {
            // get the real atom index
            int iAtom = orbitalMap(iOrbStart + npl * iab).ia;
            SxVecRef<PrecCoeffG> bcol = atomBlock.colRef(iab),
                                 T = gk.getPhaseFactors(is, iAtom);
#           ifdef USE_OPENMP
#           pragma omp parallel for if (ng > sxChunkSize)
#           endif
            for (int ig = 0; ig < ng; ig++)
               bcol(ig) = psi(ig) * T(ig).conj ();
         }

         // project current block
         SxIdx currentOrbs(iOrbStart, iOrbStart + nAtomBlock * npl - 1);
         res(currentOrbs) <<= proj(is).overlap (atomBlock);

         iOrbStart += nAtomBlock * npl;
      }
      SX_CHECK (iOrbStart == iOrb, iOrb - iOrbStart, nAtoms, npl);
   }
   SX_VALIDATE_VECTOR(res);
   res.auxData = psi.auxData;
   res.setBasis (this);
   if (extraPhase.getSize () > 0)
      applyPhaseConj (&res, extraPhase(ik));
   return res;
}
#endif

int SxAOBasis::getLMax () const
{
   int result = 0;
   int nOrbs = getNOrb ();
   for (int iOrb = 0; iOrb < nOrbs; iOrb++)  {
      int is = orbitalMap(iOrb).is;
      int io = orbitalMap(iOrb).io;
      int l = refOrbMap(is)(io).l;
      if (result < l) result = l;
   }

   return result;
}

PsiG SxAOBasis::getAOinG (int ik) const
{
   SX_CHECK (cacheRefOrb == CacheAll);
   if (cacheRefOrb == CacheCurrentK)  {
      SX_CHECK (refOrbCachedK == ik, refOrbCachedK, ik);
      ik = 0;
   }
   SX_CHECK (ik >= 0 && ik < refOrbitals.getSize (),
             ik, refOrbitals.getSize ());
   const SxGBasis *gPtr
      = dynamic_cast<const SxGBasis *> (refOrbitals(ik)(0).getBasisPtr());
   ssize_t ng = gPtr->g2.getSize ();
   ssize_t nOrbitals = orbitalMap.getSize ();
   PsiG result(ng,nOrbitals);
   for (int iOrb = 0; iOrb < nOrbitals; iOrb++)  {
     result.colRef(iOrb) <<= getAOinG (ik,iOrb);
   }
   result.setBasis(gPtr);
   result.auxData.ik = ik;

   return result;
}

PsiG SxAOBasis::getAOinG (int ik, int iOrb) const
{
   SX_CHECK (cacheRefOrb == CacheAll);
   int trueK = ik;
   if (cacheRefOrb == CacheCurrentK)  {
      SX_CHECK (refOrbCachedK == ik, refOrbCachedK, ik);
      ik = 0;
   }
   SX_CHECK (ik >= 0 && ik < refOrbitals.getSize (),
             ik, refOrbitals.getSize ());
   SX_CHECK (iOrb >= 0 && iOrb < orbitalMap.getSize (),
             iOrb, orbitalMap.getSize ());
   const SxOrbitalIndex &idx = orbitalMap(iOrb);
   const PsiG &refOrb = refOrbitals(ik)(idx.is);

   const SxGBasis *gPtr
      = dynamic_cast<const SxGBasis *> (refOrb.getBasisPtr ());
   SX_CHECK (gPtr);
   PsiG res;
   res = gPtr->getPhaseFactors(idx.is,idx.ia) * refOrb.colRef(idx.io);
   if (extraPhase.getSize () > 0) res *= extraPhase(ik)(iOrb);
   const SxAoIndex &nlm = refOrbMap(idx.is)(idx.io);
   res.auxData.ik = trueK;
   res.auxData.setAtom (idx.is, idx.ia);
   res.auxData.setOrb (nlm.n, nlm.l, nlm.m);

   SX_VALIDATE_VECTOR (res);
   return res;
}

SxVecRef<SxComplex16> SxAOBasis::getOverlap (int ik) const
{
   // --- single k caching
   if (cacheOverlap == CacheCurrentK)  {
      if (ik != overlapCachedK)  {
         overlapCachedK = ik;
         overlap(0).unref ();
         overlap(0) = calculateOverlap(ik);
      }
      return overlap(0);
   }
   // --- all k caching
   if (cacheOverlap == CacheAll)  {
      SX_CHECK (ik >= 0 && ik < overlap.getSize (),
                ik, overlap.getSize ());
      if (overlap(ik).getSize () == 0)
         overlap(ik) = calculateOverlap(ik);
      return overlap(ik);
   }
   // no caching
   SX_CHECK (cacheOverlap == Unknown || cacheOverlap == Recompute);

   if (cacheOverlap == Unknown)  {
      cout << "Warning: undefined caching behaviour for overlap matrices\n"
              "         in SxAOBasis::getOverlap. Recomputing overlap...\n";
   }

   return calculateOverlap(ik);
}

SxVecRef<SxComplex16> SxAOBasis::getInverseOverlap (int ik) const
{
   SX_CLOCK (Timer::AoTotalTime);
   // --- single k caching
   if (cacheInverse == CacheCurrentK)  {
      if (ik != invOverlapCachedK)  {
         invOverlapCachedK = ik;
         invOverlap(0).unref ();
         // get the overlap matrix
         const SxVecRef<SxComplex16> &S = getOverlap(ik);
         SX_CLOCK (Timer::AoSInversion);
         invOverlap(0) = SxVector<SxComplex16> (S).inverse ();
      }
      return invOverlap(0);
   }
   // --- all k caching
   if (cacheInverse == CacheAll)  {
      SX_CHECK (ik >= 0 && ik < invOverlap.getSize (),
                ik, overlap.getSize ());
      if (invOverlap(ik).getSize () == 0)  {
         SX_CLOCK (Timer::AoSInversion);
         invOverlap(ik) = SxVector<SxComplex16> (getOverlap(ik)).inverse ();
      }
      return invOverlap(ik);
   }
   // no caching
   SX_CHECK (cacheInverse == Unknown || cacheInverse == Recompute);

   if (cacheInverse == Unknown)  {
      cout << "Warning: undefined caching behaviour for overlap matrices\n"
              "         in SxAOBasis::getOverlap. Recomputing overlap...\n";
   }

   SxVector<SxComplex16> res = getOverlap(ik);
   SX_CLOCK (Timer::AoSInversion);
   return std::move(res).inverse ();
}

void SxAOBasis::setOverlapCaching (enum Caching mode, int nk)
{
   SX_CHECK (mode == Recompute ||
             mode == CacheCurrentK ||
             mode == CacheAll);
   if (mode == Recompute)  {
      overlap.resize (0);
      overlapCachedK = -1;
   }
   if (mode == CacheCurrentK)  {
      overlap.resize (1);
      overlapCachedK = -1;
   }
   if (mode == CacheAll)  {
      if (nk <= 0) nk = int(refOrbitals.getSize ());
      if (cacheOverlap == CacheCurrentK && overlapCachedK != -1)  {
         SxVecRef<SxComplex16> ovlp = overlap(0);
         overlap.resize (nk);
         overlap(overlapCachedK).unref (); // needed for nk=1
         overlap(overlapCachedK) = ovlp;
      }
      overlapCachedK = -1;
   }
   cacheOverlap = mode;
}

void SxAOBasis::setInvOverlapCaching (enum Caching mode, int nk)
{
   SX_CHECK (mode == Recompute ||
             mode == CacheCurrentK ||
             mode == CacheAll);
   if (mode == Recompute)  {
      invOverlap.resize (0);
      invOverlapCachedK = -1;
   }
   if (mode == CacheCurrentK)  {
      invOverlap.resize (1);
      invOverlapCachedK = -1;
   }
   if (mode == CacheAll)  {
      if (nk <= 0) nk = int(refOrbitals.getSize ());
      if (cacheInverse == CacheCurrentK && invOverlapCachedK != -1)  {
         SxVecRef<SxComplex16> ovlp = overlap(0);
         invOverlap.resize (nk);
         invOverlap(invOverlapCachedK).unref (); // needed for nk=1
         invOverlap(invOverlapCachedK) = ovlp;
      } else {
         overlap.resize (nk);
      }
      invOverlapCachedK = -1;
   }
   cacheInverse = mode;
}

class SxAOBasisProjGrad
: public SxProjMatrix<PrecCoeffG>::SaveProjections
{
   public:
      const SxArray<SxOrbitalIndex> &orbitalMap;
      const SxArray<PsiG> &refOrbitals;
      mutable SxVector<PrecCoeffG> T; // cached phase factors ("atom translation")

      const SxGBasis &gk;

      SxAOBasisProjGrad (const SxAOBasis &ao, int ik)
         : SxProjMatrix<PrecCoeffG>::SaveProjections (
              3 * ao.getNOrb (),
              ao.blockSize,
              (int)ao.refOrbitals(ik)(ao.orbitalMap(0).is).getNRows ()),
         orbitalMap(ao.orbitalMap),
         refOrbitals(ao.refOrbitals(ik)),
         gk(refOrbitals(ao.orbitalMap(0).is).getBasis<SxGBasis> ())
      {
         // empty
      }
      virtual ~SxAOBasisProjGrad () {}

      HAS_TARGET_GETPROJECTOR;
      NO_GETFACTOR;
      virtual void getProjector(int i, SxVecRef<PrecCoeffG> *target) const
      {
         int iOrb = i / 3;
         int iDir = i - 3 * iOrb;
         SX_CHECK (iOrb >= 0 && iOrb < orbitalMap.getSize (),
                   iOrb, orbitalMap.getSize ());
         const SxOrbitalIndex &idx = orbitalMap(iOrb);
         const PsiRef &refOrb = refOrbitals(idx.is).colRef (idx.io);

         if (T.getSize () <= 0
             || T.auxData.is != idx.is
             || T.auxData.ia != idx.ia)
         {
            SxVecRef<PrecCoeffG> phi = refOrb.getBasis<SxGBasis> ().getPhaseFactors (idx.is, idx.ia);
            T.reformat (nElements, 3);
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
            for (int jDir = 0; jDir < 3; ++jDir)  {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
               for (int ig = 0; ig < nElements; ++ig)  {
                  T(ig + jDir * nElements) = I * phi(ig) 
                                           * gk.gVec(ig + jDir * nElements);
               }
            }
            T.auxData.setAtom (idx.is, idx.ia);
         }

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
         for (int ig = 0; ig < nElements; ++ig)
            (*target)(ig) =  T(ig + iDir * nElements) * refOrb(ig);
      }
};

SxArray<SxVector<SxComplex16> >
SxAOBasis::gradProject (const PsiRef &psi) const
{
   //SX_CLOCK (Timer::GradProject);
   SX_CHECK(psi.getSize () >0);

   // --- find k-point
   SX_CHECK (cacheRefOrb == CacheAll);
   int ik = psi.auxData.ik;
   int nk = int(refOrbitals.getSize ());
   PsiRef psiS = SPtr->apply (psi);
   if (ik >= 0 && ik < nk)  {
      SX_CHECK (psiS.getBasisPtr () == refOrbitals(ik)(0).getBasisPtr ());
   } else if (cacheRefOrb == CacheAll)  {
      for (ik = 0; ik < nk; ++ik)  {
         if (refOrbitals(ik)(0).getBasisPtr () == psiS.getBasisPtr ())
            break;
      }
      if (ik == nk)  {
         cout << "Projection from unregistered G-basis in SxAOBasis\n";
         SX_EXIT;
      }
   }
   return gradProject(psiS, ik);
}

SxArray<SxVector<SxComplex16> >
SxAOBasis::gradProject (const PsiRef &psiS, int ik) const
{
   SX_CLOCK (Timer::AoTotalTime);
   SX_CLOCK (Timer::AoGradProject);
   int nPsi = (int)psiS.getNCols ();
   if (nPsi==0) nPsi = 1;
   int nOrb = getNOrb ();

   SxArray<SxVector<SxComplex16> > result(3);
   SxVecRef<SxComplex16> allProj;
   allProj = SxAOBasisProjGrad(*this, ik).getProjectionFromExtended (psiS);
   allProj.reshape (3, nOrb * nPsi);
   for (int iDir = 0; iDir < 3; ++iDir)  {
      result(iDir) = allProj.rowRef (iDir); // performs a copy!
      result(iDir).reshape (nOrb, nPsi);
      if (extraPhase.getSize () > 0)
         applyPhaseConj (&result(iDir), extraPhase(ik));
   }
   return result;
}



