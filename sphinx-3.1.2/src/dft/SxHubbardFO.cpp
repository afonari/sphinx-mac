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

#include <SxHubbardFO.h>
#include <SxNeighbors.h>
#include <SxEigensystem.h>
#include <SxDiagMat.h>
#include <SxSimpleParser.h>
#include <SxTextIO.h>
#include <SxPAWOverlap.h>
#include <SxHubbardU.h>

/** \brief Get overlap <chi1|chi2> between two AO-like functions

    @param f1    radial part of function chi1 in G-space
    @param f2    radial part of function chi2 in G-space
    @param jsb   precomputed Bessel functions for |l1-l2|..(l1 + l2), only
                 every second L is needed, i.e. |l1-l2|+0,2,4,6...
    @param Ylm   precomputed Y_lm(r1 - r2), orientation matters (not r2-r1),
                 must include the Y_lm normalization factor
    @param clebsch Clebsch-Gordan coefficients until l1+l2
    @return the overlap integral

    Assuming that
    \f[ \chi_i(\mathbf g) = f_i(|\mathbf g|) Y_{l_im_i}(\hat \mathbf g) \f]
    This is computed as
    \f[
    \langle \chi_1 | \chi_2 \rangle = \sum_L \int g^2 dg f_1(g)f_2(g)j_L(g|r_{21}|)
                                             \langle l_1m_1\, LM|l_2m_2\rangle Y_{LM}(\hat \mathbf r_{21})
    \f]
    with \f$\mathbf r_{21} = \mathbf r_1 - \mathbf r_2
  */
double getOverlap (const SxVecRef<double> &f1,
                   const SxVecRef<double> &f2,
                   const SxArray<SxVector<double> > &jsb,
                   const SxVector<double> &Ylm,
                   const SxYlm::SxClebschTable &clebsch)
{
   if (Ylm.getSize () == 0)  {
      if (f1.auxData.l != f2.auxData.l || f1.auxData.m != f2.auxData.m)
         return 0.;
      return tr(f1 * f2);
   }
   return SxRadialBasis::getOverlap (f1.auxData.l, f1.auxData.m,
                                     f2.auxData.l, f2.auxData.m,
                                     SxRadialBasis::computeMatElem (f1,f2,jsb),
                                     Ylm, clebsch);
}

SxHubbardFO::SxHubbardFO (int siteOffsetIn)
   : SxHubbardXO (siteOffsetIn),
     selectThreshold(0.7),
     selectType(-1),
     fullPAWNorm(false), globalProj(true),
     cacheGlobalOvlp(SxAOBasis::CacheAll)
{
   // empty
}

void SxHubbardFO::read (const SxSymbolTable *table,
                        const SxAtomicStructure &structure,
                        const SxPtr<SxPAWPot> &potPtr)
{
   // --- parse the symbol table
   SxArray<double> rCut, truncWidth;
   SxList<double> eList;
   SxList<int> centralAtomsList;
   int selectL = -1;
   SYMBOLPARSE(table) SYMBOLGROUP("fragment") {
      SYMBOLGROUP ("AOBasis")  {
         // read orbitals in real space (changed to radial G basis below)
         refOrbsG.setup (SYMBOLGROUP_TABLE);

         // --- read further species-dependent stuff
         int nSpecies = SYMBOL_COUNT("species");
         rCut.resize (nSpecies);
         truncWidth.resize (nSpecies);
         int iSpecies = 0;
         FOREACH_SYMBOLGROUP("species")  {
            FOREACH_SYMBOLGROUP("orbital")  {
               eList << double(SYMBOLGET("energy"));
            }
            rCut(iSpecies)       = SYMBOLGET ("truncate") || -1.;
            truncWidth(iSpecies) = SYMBOLGET ("truncWidth") || 0.5;
            iSpecies++;
         }
      }
      orbEnergy = eList;
      kappa = SYMBOLGET ("kappa");
      // maximum distance between central atom and ligands
      maxDist = SYMBOLGET ("maxDist");
      // minimum number of atoms in fragment
      minNumAtoms = (SYMBOLGET ("minLigands") || 0) + 1;

      // fragments consists of a central atom and nearby ligands
      SYMBOLGROUP ("center")  {
         if (HAVE_SYMBOL("label"))  {
            SxString label = SYMBOLGET("label");
            if (structure.hasLabels ())  {
               const SxArray<SxString> &labels = structure.getLabels ();
               SX_LOOP(iAtom)  {
                  if (labels(iAtom) == label)
                     centralAtomsList << int(iAtom);
               }
            } else {
               cout << "Central atoms for fragments to be selected by label '"
                    << label << "', but structure has no labels" << endl;
               SX_QUIT;
            }
         } else {
            SxString chemSym = SYMBOLGET("element");
            int iSpecies = int(structure.getElements ().findPos (chemSym));
            if (iSpecies != -1)  {
               for (int ia = 0; ia < structure.getNAtoms (iSpecies); ++ia)
                     centralAtomsList << structure.getIAtom (iSpecies, ia);
            } else {
               cout << "Central atoms for fragments to be selected by element '"
                    << chemSym << "', but there is no such element." << endl;
               cout << "Elements are " << structure.getElements () << endl;
               SX_QUIT;
            }
         }
      }
      // whether to use 3-body PAW corrections for <pi(@atom 1)|S(@atom 3)|chi(@atom 2)>
      fullPAWNorm = SYMBOLGET("fullPAWNorm").toBool ();
      // whether to normalize with respect to global <pi|chi> vs. intra-fragment <pi|chi>
      globalProj = ! SYMBOLGET("fragmentOnlyNorm").toBool ();
      // select orbitals
      if (HAVE_SYMBOL("selectByIndex"))  {
         siteOrbIdx = SYMBOLGET("selectByIndex");
         for (auto &orbIdx : siteOrbIdx) orbIdx--;
         selectType = -1;
      } else {
         selectL = SYMBOLGET("selectByL");
         selectThreshold = SYMBOLGET("threshold") || 0.7;
         cout << "Selecting fragment orbitals if l=" << selectL
              << " character of central atom is >= " << selectThreshold << endl;
      }
      nRad = SYMBOLGET("nRadGrid") || 200;
      verbose = SYMBOLGET("verbose").toBool ();
   }
   if (centralAtomsList.getSize () == 0)  {
      cout << "No single central atom found for fragment orbitals" << endl;
      SX_QUIT;
   }
   centralAtoms = centralAtomsList;
   cout << "Using " << (fullPAWNorm ? 3 : 2)
        << "-center PAW normalization (fullPAWNorm flag)" << endl;

   if (!radGPtr)  {
      double eCut = SxGBasis::getECut (table->topLevel ());
      double gMax = sqrt(eCut), dg = 0.003;
      int ng = int(gMax/dg) + 1;
      radGPtr = SxPtr<SxRadialBasis>::create (0., gMax, ng, false);
   }
   // --- transform orbitals to radial G space
   SxConstPtr<SxRadBasis> radBasisPtr = refOrbsG.getRadBasisPtr ();
   SX_LOOP2(is,iot) {
      // --- set up reference orbitals in G space
      refOrbsG(is,iot) = *radGPtr | refOrbsG(is,iot);
   }
   refOrbsG.setBasis (radGPtr);

   // --- set up PAW projectors in radial G space
   // ugly, but easy: destructor call + placement new.
   projG.~SxAtomicOrbitals ();
   new (&projG) SxAtomicOrbitals(potPtr->lPhi, radGPtr);
   const SxPAWPot &pawPot = *potPtr;
   SX_LOOP(is)  {
      const SxVector<double> &projRad = pawPot.pPS(is);
      for (int n = 0; n < pawPot.lPhi(is).getSize (); ++n)  {
         int l = pawPot.lPhi(is)(n);
         SxVecRef<double> proj
            = const_cast<SxVector<double>&>(projRad).colRef (n);
         proj.auxData.setAtom ((int)is, -1);
         proj.auxData.setOrb (n, l, NONE_M);
         projG(is, n) = *radGPtr | proj;
      }
   }

   // --- set up refOrbsProjG
   // resize refOrbsProjG
   // ugly, but easy: destructor call + placement new.
   refOrbsProjG.~SxAtomicOrbitals ();
   new (&refOrbsProjG) SxAtomicOrbitals (refOrbsG.refOrbMap, radGPtr);
   {
      // for plotting...
      SX_LOOP (is) {
         double width = truncWidth(is);
         double rMax = rCut(is) > 0. ? (rCut(is) + 6 * width) : radBasisPtr->getRMax ((int)is);
         SxRadialBasis R(0., rMax, nRad, true);
         SX_LOOP(iType)  {
            const SxVecRef<double> &refOrb = refOrbsG (is,iType);
            SxVecRef<double> &proj = refOrbsProjG (is,iType);
            int l = refOrb.auxData.l;
            prefix = structure.getElements ()(is)
                   + SxString::sprintf("%d-l%d-", int(iType), l);
            if (rCut(is) > 0.)  {
               proj = truncateShapeG (refOrb, rCut(is), width);
            } else {
               proj = SxVector<double>(refOrb); // copy full shape
            }
            proj.auxData.m = NONE_M;

            if (!fullPAWNorm)  {
               // --- PAW norm corrections go into projectors
               // Compute <p_i | chiTrunc>
               SxVector<double> pTrunc(pawPot.getNProjType (is));
               for (int ipt = 0; ipt < pTrunc.getSize (); ++ipt)  {
                  pTrunc(ipt) = (pawPot.lPhi(is)(ipt) == refOrb.auxData.l)
                              ? tr(projG(is,ipt) * proj)
                              : 0.;
               }
               // S_{ij} <p_i | chiTrunc>
               pTrunc = pawPot.deltaS(is) ^ pTrunc;
               // |p_j> S_{ij} <p_i | chiTrunc>
               SX_LOOP(ipt)
                  if (projG(is,ipt).auxData.l == refOrb.auxData.l)
                     proj.plus_assign_ax (pTrunc(ipt), projG(is,ipt));

               if (verbose) {
                  const SxVector<double> &r = R.getRadFunc ();

                  SxVector<double> shapeR = R | proj;
                  SxTextIO(prefix + "projS.dat").writeXYPlot (r, shapeR);
                  SxTextIO(prefix + "shapeTruncSG.dat")
                     .writeXYPlot (radGPtr->getRadFunc (), proj);
               }
            }
         }
      }
   }

   // --- setup
   // * offsetProjGlobal ... offset of different species in projector list
   // * mapType          ... orbital "type" for each fragment orbital
   offsetProjGlobal.resize (structure.getNSpecies ());
   mapType.resize (structure.getNSpecies ());

   int nAOType = 0;
   int offProj = 0;
   int isCentral = structure.getISpecies (centralAtoms(0));
   SX_LOOP(ia)  {
      if (structure.getISpecies (centralAtoms(ia)) != isCentral)  {
         cout << "Not all central atoms have same species. This won't work." << endl;
         SX_QUIT;
      }
   }
   SX_LOOP(is)  {

      int nRef = refOrbsProjG.getNOrbs (is);
      if (nRef == 0) continue;
      // offset of species in global orbital/projector list ipg=(is,ia,iot)
      offsetProjGlobal(is) = offProj;
      offProj += nRef * structure.getNAtoms (int(is));

      mapType(is).resize (nRef);
      // --- loop over reference orbitals
      SxVecRef<double> psiGRef;
      SX_LOOP(io) {
         SxAoIndex nlm = refOrbsProjG.refOrbMap(is)(io);

         if (is == isCentral && nlm.l == selectL)
            selectType = nAOType;

         // --- extract orbital information from refOrbMap
         // map local orbital number iol=(n,l,m) to orbital type (iot)
         mapType(is)(io) = nAOType;
         if (nlm.l == nlm.m) nAOType++;
      }

   }
   // now find the fragments
   findFragments (structure);

   // --- and compute the site rotation matrices
   setupSiteU (potPtr, structure);
}

void SxHubbardFO::setupProjGk (const SxGkBasis &gkBasis,
                               const SxPtr<SxPartialWaveBasis> &pBasis)
{
   // --- create projector basis
   SxPtr<SxOverlapBase> SPtr;
   if (fullPAWNorm)
     SPtr =  SxPtr<SxPAWOverlap>::create (pBasis);
   else
     SPtr = SxPtr<SxPWOverlap>::create ();
   aoProj = SxPtr<SxAOBasis>::create (gkBasis, refOrbsProjG.getMuSet (), SPtr);
   if (globalProj)
      //setCaching (SxAOBasis::CacheCurrentK);
      setCaching (SxAOBasis::CacheAll);
   cout << "Using " << (globalProj ? "global" : "per-fragment")
        << " proj->chi transformation (fragmentOnlyNorm flag)" << endl;

}

void SxHubbardFO::findFragments (const SxAtomicStructure &structure)
{
   SxGrid grid(structure, 10);
   SxNeighbors nn;
   SxStack<Coord> latShiftList;
   SxStack<int> nFragOrbList, aoTypeList;
   SxStack<int> mapFragOrbProjList;
   SxStack<SxArray<int> > atomIdxList;

   SX_LOOP(iCenter)  {
      int iTl = centralAtoms(iCenter);
      nn.compute (grid, structure, structure.getAtom (iTl), maxDist,
                  SxNeighbors::StoreAbs | SxNeighbors::IncludeZeroDistance);
      if (nn.getSize () < minNumAtoms) continue;

      int nFragOrbitals = 0;

      atomIdxList << SxArray<int> ();
      atomIdxList.top ().resize (nn.getSize ());
      for (int in = 0; in < nn.getSize (); ++in)  {
         int jTl = nn.idx(in);
         int ja, js = structure.getISpecies (jTl, &ja);
         // only consider central atom and ligands of right species
         if (jTl != iTl && refOrbsG.getNOrbTypes(js) == 0) continue;

         // --- now set up orbital info
         int nAO = refOrbsG.getNOrbs (js);
         nFragOrbitals += nAO;
         int offsetJOrb = offsetProjGlobal(js) + ja * nAO;
         for (int iAO = 0; iAO < nAO; ++iAO)  {
            // lattice shift
            latShiftList << nn.absPositions(in) - structure.getAtom (jTl);
            mapFragOrbProjList << offsetJOrb + iAO;
            aoTypeList << mapType(js)(iAO);
         }
         cout << structure.getElements ()(js) << " @ " << nn.absPositions(in)
              << " with "  << nAO << " orbitals" << endl;
         atomIdxList.top ()(in) = jTl;
      }
      // update fragment orbital list
      nFragOrbList << nFragOrbitals;
   }
   // tell SxHubbardXO the number of sites
   SxHubbardXO::nSite = (int) atomIdxList.getSize ();
   if (nSite == 0)  {
      cout << "No sites found for fragment orbitals" << endl;
      cout << "Minimum ligands = " << (minNumAtoms-1) << endl;
      SX_QUIT;
   }

   // --- now set up final arrays
   // fragment orbital basis stuff
   aoType = aoTypeList;
   latShift = latShiftList;
   mapFragOrbProj = mapFragOrbProjList;

   // site specific stuff
   nFragOrb = std::move (nFragOrbList);
   atomIdx = std::move (atomIdxList);
   // set up offsets into global fragment orbitals
   siteStart.resize (getNSite ());
   for (int i = 0, off = 0; i < getNSite (); off += nFragOrb(i++))
      siteStart(i) = off;

   // set up orbital map
   orbMap.resize (nSite);
   SX_LOOP(iSite)  {
      orbMap(iSite).resize (nFragOrb(iSite));
      int ifo = 0;
      SX_LOOP(ifa)  {
         int ia, is = structure.getISpecies (atomIdx(iSite)(ifa), &ia);
         int nAO = refOrbsProjG.getNOrbs (is);
         for (int io = 0; io < nAO; ++io)
            orbMap(iSite)(ifo++).set (is, ia, io);
      }
   }
}


void computeYlmNormalized(int lMax, const Coord &dR, SxVector<double>* Ylm)
{
   SX_CHECK (Ylm);
   Ylm->resize (sqr(lMax + 1));
   SxYlm::getYlmArray(lMax, dR, Ylm);
   for (int l = 0, lm = 0; l <= lMax; ++l)
      for (int m = -l; m <= l; ++m, ++lm)
         (*Ylm)(lm) *= SxYlm::getYlmNormFactor(l,m);
}

SxVector<double>
SxHubbardFO::getPAWProjections (ssize_t iSite, const SxPAWPot &pawPot,
                                     const SxAtomicStructure &structure,
                                     const SxAtomicOrbitals &shapeG) const
{
   // --- setup <pPAW|AOps>, i.e., overlap between PAW projectors and
   //     pseudo-AO basis within the current fragment
   int nFO = nFragOrb(iSite);
   int nProjSite = 0;
   for (int jTl : atomIdx(iSite))
      nProjSite += pawPot.getNProj (structure.getISpecies (jTl));

   SxVector<double> pFO(nProjSite, nFO);
   const SxRadialBasis &gRad = *projG.getRadGBasisPtr ();
   // loop over AO basis
   SX_LOOP(ifo)  {
      int ipSite = 0;
      int latOff = siteStart(iSite);
      const SxOrbitalIndex &orb = orbMap(iSite)(ifo);
      Coord R1 = latShift(latOff + ifo);
      // loop over atoms within fragments
      SX_LOOP(jAtom)  {
         int jTl = atomIdx(iSite)(jAtom);
         int js = structure.getISpecies (jTl);
         Coord R2 = latShift(latOff);
         Coord dR = structure.getAtom(jTl) + R2
                  - structure.getAtom(orb.is, orb.ia) - R1;
         int lMax = pawPot.lMax(js) + refOrbsG.getLMax (orb.is);
         double dist = dR.norm ();
         SxArray<SxVector<double> > jsb(lMax + 1);
         SxVector<double> Ylm;
         if (dist > 1e-7)  {
            SX_LOOP (l)
               jsb(l) = SxRadialBasis::jsb ((int)l, dist * gRad.getRadFunc ());
            computeYlmNormalized (lMax, dR, &Ylm);
         }

         // loop over PAW projectors of this atom
         SX_LOOP(ipl)  {
            pFO(ipSite++, ifo) = getOverlap (projG.getOrb(js, ipl),
                                             shapeG.getOrb(orb.is, orb.io),
                                             jsb, Ylm, pawPot.clebschGordan);
         }
         latOff += shapeG.getNOrbs (js); // nAO for atom J
      }
   }
   return pFO;
}


SxVector<double>
SxHubbardFO::computeSiteOverlap (ssize_t iSite, const SxPAWPot &pawPot,
                                 const SxAtomicStructure &structure) const
{
   int nFO = nFragOrb(iSite);
   SxVector<double> pFO = getPAWProjections (iSite, pawPot, structure,
                                             refOrbsG);

   // --- now set off the overlap matrix of the fragment's AO basis
   //     including PAW corrections
   SxVector<double> siteS(nFO, nFO);
   int latOff = siteStart(iSite);
   const SxRadialBasis &gRad = *projG.getRadGBasisPtr ();
   SX_LOOP2(ifo, jfo)  {

      // --- PAW overlap correction <mu_i|p_k><p_k|mu_j>
      int ip = 0;
      double overlapCorr = 0.;
      int pOff = 0;
      SX_LOOP(kAtom)  {
         int kTl = atomIdx(iSite)(kAtom);
         int ks = structure.getISpecies (kTl);
         int pOffNext = pOff + refOrbsG.getNOrbs (ks);
         if (    (pOff > ifo || pOffNext <= ifo)
              && (pOff > jfo || pOffNext <= jfo))
         {
            // avoid 3-body PAW norm corrections
            ip += pawPot.getNProj (ks);
            pOff = pOffNext;
            continue;
         }
         pOff = pOffNext;
         const SxArray<int> &lPhi = pawPot.lPhi(ks);
         const SxVector<double> deltaS = pawPot.deltaS(ks);
         int jp0 = ip;
         SX_LOOP(ipt)  {
            int nm = 2 * lPhi(ipt) + 1;
            int jp = jp0;
            SX_LOOP(jpt)  {
               if (lPhi(ipt) == lPhi(jpt))  {
                  double sum = 0.;
                  for (int mm = 0; mm < nm; ++mm)
                     sum += pFO(ip + mm, ifo) * pFO(jp + mm, jfo);
                  overlapCorr += sum * deltaS(ipt, jpt);
               }
               jp += 2 * lPhi(jpt) + 1;
            }
            ip += nm;
         }
      }

      // --- compute plane-wave overlap + PAW corrections
      const SxOrbitalIndex &orb1 = orbMap(iSite)(ifo);
      const SxOrbitalIndex &orb2 = orbMap(iSite)(jfo);
      Coord dR = structure.getAtom(orb1.is, orb1.ia)
               - structure.getAtom(orb2.is, orb2.ia)
               + latShift (latOff + ifo)
               - latShift (latOff + jfo);
      double dist = dR.norm ();
      int lMax = refOrbsG.getLMax (orb1.is) + refOrbsG.getLMax (orb2.is);
      SxArray<SxVector<double> > jsb(lMax + 1);
      SxVector<double> Ylm;
      if (dist > 1e-7)  {
         SX_LOOP (l)
            jsb(l) = SxRadialBasis::jsb ((int)l, dist * gRad.getRadFunc ());
         computeYlmNormalized (lMax, dR, &Ylm);
      }

      //                raw overlap
      siteS(ifo, jfo) = getOverlap(refOrbsG.getOrb(orb1.is,orb1.io),
                                   refOrbsG.getOrb(orb2.is,orb2.io),
                                   jsb, Ylm, pawPot.clebschGordan)
      //                PAW correction
                      + overlapCorr ;
   }
   return siteS;
}


SxVector<double>
SxHubbardFO::getHam (ssize_t iSite, const SxVector<double> &siteS) const
{
   int nFO = nFragOrb(iSite);
   int offset = siteStart(iSite);
   // set up site Hamiltonian according to extended Hueckel theory
   SxVector<SxComplex16> H(nFO,nFO);
   for (int io = 0; io < nFO; ++io)  {
      double e1 = orbEnergy(aoType(offset + io));
      for (int jo = 0; jo < nFO; ++jo)  {
        double e2 = orbEnergy(aoType(offset + jo));
        // Wolfsberg-Helmholz approximation
        H(io,jo) = 0.5 * (e1 + e2) * siteS(io,jo) * (io == jo ? 1 : kappa);
      }
   }
   return H;
}

void SxHubbardFO::setupSiteU (const SxPtr<SxPAWPot> potPtr,
                              const SxAtomicStructure &structure)
{
   siteU.resize (getNSite ());
   SX_LOOP (iSite)  {

      // --- setup overlap matrix
      SxVector<double> siteS = computeSiteOverlap (iSite, *potPtr, structure);
      cout << "siteS diag: " << siteS.diag () << endl;
      cout << (siteS - siteS.transpose ()).normSqr () << endl;

      // set up site Hamiltonian according to extended Hueckel theory
      SxVector<double> H = getHam (iSite, siteS);
      cout << (H - H.transpose ()).normSqr () << endl;

      // --- solve the generalized eigenvalue Problem H psi = epsilon S psi
      // (1) diagonalize S
      SxSymEigensystem<double> eig (siteS);
      cout << "S eigenvalues = " << eig.vals << endl;

      // (2) set up U = S^{-1/2} from eigenbasis (Loewdin orthogonalization)
      SxDiagMat<double> Ieps(1./sqrt(eig.vals));
      SxVector<double> U = eig.vecs ^ Ieps ^ eig.vecs.adjoint();

      H = U ^ H ^ U;
      // (3) solve normal eigenvalue problem
      eig.compute (std::move(H), All);

      int nFO = nFragOrb(iSite);
      int offset = siteStart(iSite);

      // --- inspect orbital character of eigenstates
      SxVector<double> vecNO = siteS ^ U ^ eig.vecs;
      SxArray<double> normType(orbEnergy.getSize ());
      SxList<int> selectOrbIdx;
      SX_LOOP(i)  {
         cout << (i + 1) << " " << eig.vals(i) << ": ";
         //cout << eig.vecs.colRef(i).absSqr () << endl;
         normType.set (0.);
         SX_LOOP(j)  {
            double c = eig.vecs(j,i);
            normType(aoType(offset + j)) += c*c;
         }
         cout << "norm by type: " << normType << endl;
         // if we are selecting by l-character, select now
         if (selectType >=0 && normType(selectType) >= selectThreshold)  {
            selectOrbIdx.append (int(i));
         }
      }
      if (selectType >= 0)
         siteOrbIdx = selectOrbIdx;

      // --- select the relevant eigenstates
      // TODO: use orbital character to decide
      SxVector<double> siteOrbs(nFO, siteOrbIdx.getSize ());
      SX_LOOP (i)  {
         siteOrbs.colRef(i) <<= eig.vecs.colRef(siteOrbIdx(i));
         cout << "Selecting state " << (siteOrbIdx(i)+1) << endl;
      }
      SxVector<double> finalU = siteS ^ U ^ siteOrbs;
      if (!globalProj)  {
         // include fragment normalization (proj|chi)^-1 in siteU
         finalU = getSiteT (iSite, *potPtr, structure).inverse () ^ finalU;
      }
      siteU(iSite) = finalU; // type-casts to SxComplex16
   }
   if (selectType >= 0)
      siteOrbIdx.resize (0);
}

SxVector<SxComplex16>
SxHubbardFO::extractFromMatGlobal (ssize_t iSite, const SxVecRef<SxComplex16> &M,
                                        const SxVecRef<SxComplex16> &phase) const
{
   int nOrb = nFragOrb(iSite);
   SxVector<SxComplex16> siteM(nOrb, nOrb);
   SX_LOOP2(jfo,ifo)  {
      siteM(ifo,jfo) = M(mapFragOrbProj(siteStart(iSite) + int(ifo)),
                         mapFragOrbProj(siteStart(iSite) + int(jfo)))
                     * phase(ifo).conj () * phase(jfo);
   }
   return siteM;
}

SxVector<double> SxHubbardFO::getSiteT (ssize_t iSite, const SxPAWPot &pawPot,
                                             const SxAtomicStructure &structure)
{
   int nFO = nFragOrb(iSite);
   SxVector<double> pFO = getPAWProjections (iSite, pawPot, structure,
                                             refOrbsG);
   SxVector<double> pFP = getPAWProjections (iSite, pawPot, structure,
                                             refOrbsProjG);

   // --- now set off the overlap matrix of the fragment's AO basis
   SxVector<double> siteT(nFO, nFO); // ifo/jfp
   int latOff = siteStart(iSite);
   const SxRadialBasis &gRad = *projG.getRadGBasisPtr ();
   SX_LOOP2(ifo, jfp)  {

      // --- PAW overlap correction <mu_i|p_k><p_k|mu_j>
      double overlapCorr = 0.;
      if (fullPAWNorm)  {
         // --- overlap is not included in fragment orbital projectors
         int ip = 0;
         SX_LOOP(kAtom)  {
            int kTl = atomIdx(iSite)(kAtom);
            int ks = structure.getISpecies (kTl);
            const SxArray<int> &lPhi = pawPot.lPhi(ks);
            const SxVector<double> deltaS = pawPot.deltaS(ks);
            int jp0 = ip;
            SX_LOOP(ipt)  {
               int nm = 2 * lPhi(ipt) + 1;
               int jp = jp0;
               SX_LOOP(jpt)  {
                  if (lPhi(ipt) == lPhi(jpt))  {
                     double sum = 0.;
                     for (int mm = 0; mm < nm; ++mm)
                        sum += pFO(ip + mm, ifo) * pFP(jp + mm, jfp);
                     overlapCorr += sum * deltaS(ipt, jpt);
                  }
                  jp += 2 * lPhi(jpt) + 1;
               }
               ip += nm;
            }
         }
      }
      // --- compute plane-wave overlap
      const SxOrbitalIndex &orb1 = orbMap(iSite)(ifo);
      const SxOrbitalIndex &orb2 = orbMap(iSite)(jfp);
      Coord dR = structure.getAtom(orb1.is, orb1.ia)
               - structure.getAtom(orb2.is, orb2.ia)
               + latShift (latOff + ifo)
               - latShift (latOff + jfp);
      double dist = dR.norm ();
      int lMax = refOrbsG.getLMax (orb1.is) + refOrbsProjG.getLMax (orb2.is);
      SxArray<SxVector<double> > jsb(lMax + 1);
      SxVector<double> Ylm;
      if (dist > 1e-7)  {
         SX_LOOP (l)
            jsb(l) = SxRadialBasis::jsb ((int)l, dist * gRad.getRadFunc ());
         computeYlmNormalized (lMax, dR, &Ylm);
      }

      //                raw overlap
      siteT(ifo, jfp) = getOverlap(refOrbsG.getOrb(orb1.is, orb1.io),
                                   refOrbsProjG.getOrb(orb2.is, orb2.io),
                                   jsb, Ylm, pawPot.clebschGordan)
      //                PAW correction
                      + overlapCorr;
   }
   return siteT;
}

SxVector<SxComplex16>
SxHubbardFO::projectToSite (const SxVecRef<SxComplex16> &psiAO,
                            ssize_t iSite,
                            const SxVecRef<SxComplex16> &phase) const
{
   int nFO = nFragOrb(iSite);
   ssize_t nStates = psiAO.getNCols ();
   ssize_t offset = siteStart(iSite);
   SxVector<SxComplex16> psiFO(nFO, nStates);
   SX_LOOP2(iState,ifo)  {
      psiFO(ifo,iState) = psiAO(mapFragOrbProj(offset + ifo),iState)
                         * phase(ifo).conj ();
   }
   const SxVecRef<SxComplex16> &U = siteU(iSite);
   const SxVecRef<SxComplex16> x = U.overlap (psiFO);
   //cout << x.overlap (x).diag () << endl;
   psiFO = U ^ (U.overlap (psiFO));
   SxVector<SxComplex16> res(psiAO.getNRows (), nStates);
   res.auxData = psiAO.auxData;
   res.set (0.);
   SX_LOOP2(iState,ifo)  {
      res(mapFragOrbProj(offset + ifo),iState) += psiFO(ifo,iState)
                                                * phase(ifo);
   }
   return res;
}

const SxVecRef<SxComplex16> SxHubbardFO::getInvOvlp (int ik) const
{
   if (cachedOvlp.getSize () > 0)  {
      if (cacheGlobalOvlp == SxAOBasis::CacheAll && cachedOvlp(ik).getSize () > 0)
         return cachedOvlp(ik);
      else if (cacheGlobalOvlp == SxAOBasis::CacheCurrentK
            && cachedOvlp(0).auxData.ik == ik)
         return cachedOvlp(0);
   }
   const SxGBasis &gk = aoProj->refOrbitals(ik)(0).getBasis<SxGBasis> ();
   PsiG aoBasis(gk.getNElements (), aoProj->getNOrb ());
   SX_LOOP(iog)  {
      const SxOrbitalIndex &orb = aoProj->orbitalMap(iog);
      const SxAoIndex &nlm = refOrbsG.refOrbMap (orb.is)(orb.io);
      aoBasis.colRef (iog) <<= gk | refOrbsG(orb.is, orb.ia, nlm.n, nlm.l, nlm.m);
   }
   aoBasis.setBasis (gk);
   // do the actual computation
   SxVector<SxComplex16> overlap
      = (*aoProj | aoBasis).adjoint ().inverse ();
   overlap.auxData.ik = ik;
   overlap.setBasis (NULL); // no well-defined basis
   ssize_t nkCache;
   switch (cacheGlobalOvlp)  {
      case SxAOBasis::CacheAll: nkCache = aoProj->refOrbitals.getSize (); break;
      case SxAOBasis::CacheCurrentK: nkCache = 1; ik = 0; break;
      default: nkCache = 0;
   }


   if (nkCache > 0) {
      if (cachedOvlp.getSize () != nkCache)
         cachedOvlp.resize (nkCache);
      cachedOvlp(ik).unref ();
      cachedOvlp(ik) = overlap;
   }

   return overlap;
}

void SxHubbardFO::setCaching (SxAOBasis::Caching mode)
{
   cacheGlobalOvlp = mode;
   ssize_t nkCache;
   switch (mode)  {
      case SxAOBasis::CacheAll      :
         nkCache = aoProj->refOrbitals.getSize ();
         break;
      case SxAOBasis::CacheCurrentK:
         nkCache = 1;
         break;
      default:
         nkCache = 0;
   }
   cachedOvlp.resize (nkCache);
}

SxVector<SxComplex16>
SxHubbardFO::getProjMatrix (const PsiG &waves,
                            const SxVecRef<double>& focc) const
{
   int ik = waves.auxData.ik;

   int nStates = int(waves.getNCols ());

   // get number of occupied states by ignoring empty states from the end
   int nOccStates = nStates;
   for ( ; nOccStates > 1; nOccStates--)
      if (fabs(focc(nOccStates-1)) > 1e-12) break;

   // --- get waves/focc for occupied states
   SxVecRef<PrecCoeffG> occWaves;
   int ng = int(waves.getNRows ());

   occWaves = waves(SxIdx(0, ng * nOccStates - 1));
   occWaves.reshape (ng, nOccStates);
   occWaves.auxData = waves.auxData;

   // get (S^-1) <mu|psi>
   SxAOBasis::TPsi aoPsi;
   //       < p | psi >
   aoPsi =  (*aoProj | occWaves);
   if (globalProj)
      aoPsi = getInvOvlp (ik) ^ aoPsi;

   // multiply aoPsi with sqrt(focc(iState))
   for (int i = 0; i < nOccStates; ++i)
      aoPsi.colRef (i) *= (focc(i) > 0. ? SxComplex16(1.) : I )
                        * sqrt(fabs(focc(i)));

   // return <mu|psi(i)>focc(i)<psi(i)|mu>
   return (aoPsi ^ aoPsi.adjoint ());
}

/// Number of orbitals associated with a specific atom
int SxHubbardFO::getNOrbAtom (ssize_t iSite, ssize_t ia) const
{
   SX_CHECK (trafoForce.getNAtoms () > 0);
   int is = trafoForce.getISpecies (atomIdx(iSite)(ia));
   return refOrbsProjG.getNOrbs (is);
}

PsiG SxHubbardFO::apply (const SxVecRef<SxComplex16> &psi) const
{
   // TODO: possibly pre-compute Hamiltonian in full projector space
   // |pi>H<pj| = <pi|chi>^{-1} [sum_s SELECT_s^{-1} U_s^adj H_s SELECT_s] <chi|pi>^{-1}
   int ik = psi.auxData.ik;
   int iSpin = psi.auxData.iSpin;
   ssize_t nStates = psi.getNCols ();
   SxAOBasis::TPsi aoPsi =  *aoProj | psi;
   if (globalProj)
      aoPsi = getInvOvlp (ik) ^ aoPsi;

   PsiG res(aoPsi.getNRows (), nStates);
   res.set (0.);

   Coord kVec = getK (ik);
   SX_LOOP (iSite) {
      ssize_t nFO = nFragOrb(iSite);
      int offset = siteStart(iSite);

      // --- setup phases for current site
      SxVector<SxComplex16> phase(nFO);
      SX_LOOP(ifo)  {
         double kR = kVec ^ latShift(offset + ifo);
         phase(ifo) = SxComplex16::phase(-kR);
      }

      // --- extract orbitals for current fragment and apply latShift phase
      SxVector<SxComplex16> psiFrag(nFO, nStates);
      for (ssize_t iState = 0; iState < nStates; ++iState)  {
         for (ssize_t ifo = 0; ifo < nFO; ++ifo)  {
            psiFrag(ifo, iState) = aoPsi(mapFragOrbProj(offset + ifo), iState)
                                 * phase(ifo).conj ();
         }
      }
      // apply site Hamiltonian
      psiFrag = siteU(iSite) ^ SxVector<SxComplex16>(hamProj(iSpin,iSite))
                             ^ siteU(iSite).overlap (psiFrag);
      // --- put back orbitals for current fragment and apply latShift phase
      for (ssize_t iState = 0; iState < nStates; ++iState)  {
         for (ssize_t ifo = 0; ifo < nFO; ++ifo)  {
            res(mapFragOrbProj(offset + ifo), iState) += psiFrag(ifo, iState)
                                                       * phase(ifo);
         }
      }
   }
   if (globalProj)
      res = getInvOvlp (ik).adjoint () ^ res;
   res.auxData = aoPsi.auxData;
   res.setBasis (*aoProj);

   return res;
}

/** \brief Compute the Hubbard energy and Hamiltonian
  @param hubbardU the Hubbard parent object containing the U parameters,
                  energy, etc.
  @param Pij the AO block density matrix
  @param structure the atomic structure
  */
void SxHubbardFO::compute (SxHubbardU *hubbardU,
                      const SxBlockDensityMatrix& Pij,
                      const SxAtomicStructure &structure)
{
   double fFull = Pij.getNSpin () == 1 ? 2. : 1.;
   double fInv = 1. / fFull;
   // rescale energy to single spin
   hubbardU->energy *= fInv;
   hubbardU->eDoubleCounting *= fInv;

   SX_LOOP2(iSpin, iSite)  {
      SxVector<double> U = siteU(iSite).real ();
      ssize_t nm = U.getNCols ();
      SxVector<double> D = U.overlap(Pij(iSpin, iSite) ^ U );
      D.auxData.iSpin = char(iSpin);
      SX_CHECK (D.getNCols () == nm, D.getNCols (), nm);
      SX_CHECK (D.getNRows () == nm, D.getNRows (), nm);
      if (verbose) {
         SxSymEigensystem<double> eig(D);
         cout << "occ eigenspace (site " << (iSite + siteOffset + 1) << "):"
              << endl;
         for (int mm = 0; mm < nm; ++mm)
            cout << eig.vals(mm) << ": " << eig.vecs.colRef(mm) << endl;
      }
      // compute HubbardU stuff from renormalized D (for nSpin=1 <=> f=2)
      hamProj(iSpin, iSite)
         = hubbardU->computeIncr((int)iSite + siteOffset, fInv * D);
   }
   hubbardU->energy *= fFull;
   hubbardU->eDoubleCounting *= fFull;
}

Coord SxHubbardFO::getK (int ik) const
{
   int is = 0;
   while (aoProj->refOrbMap(is).getSize () == 0) ++is;
   return aoProj->refOrbitals(ik)(is).getBasis<SxGBasis> ().getK ();
}

/** \brief Add contribution of one k-point to the AO block density matrix
  @param Pij The block density matrix to be computed
  @param waves wave functions for current k-point
  @param weight k-point weight
  @param focc occupation number for each state
  */
void SxHubbardFO::addToRho (SxBlockDensityMatrix *Pij,
                       const SxVecRef<PrecCoeffG> &waves,
                       double weight,
                       const SxVecRef<double> &focc) const
{
   SX_CHECK (Pij);
   int iSpin = waves.auxData.iSpin;
   Coord kVec = getK (waves.auxData.ik);

   SX_CHECK (iSpin >=0 && iSpin < Pij->getNSpin (), iSpin, Pij->getNSpin ());
   SxVector<SxComplex16> P = getProjMatrix (waves, focc);
   SX_LOOP (iSite) {
      int nFO = nFragOrb(iSite);
      int offset = siteStart(iSite);

      // --- setup phases for current site
      SxVector<SxComplex16> phase(nFO);
      SX_LOOP(ifo)  {
         double kR = kVec ^ latShift(offset + ifo);
         phase(ifo) = SxComplex16::phase(-kR);
      }

      SxVector<SxComplex16> siteP = extractFromMatGlobal(iSite, P, phase);
      // P = U.adjoint () ^ P ^ U
      //siteP = siteU(iSite).overlap(siteP ^ siteU(iSite) );
      (*Pij)(iSpin,iSite).plus_assign_ax (weight, siteP.real ());
   }
}
/** \brief Compute contribution to gradient of the AO block density matrix

    @param fermi Fermi occupations
    @param P     <AO|Psi> projections
    @param gradP <d/dtau AO | Psi> projection gradients
    @return gradient matrices for all sites. The first ("spin") index is
            used for the direction.

    @note k-point and spin must be set in the projection's (P) auxData.
  */
SxBlockDensityMatrix
SxHubbardFO::computeGradPij (const SxFermi &fermi,
                             const SxVecRef<SxComplex16> &P,
                             const SxArray<SxVector<SxComplex16> > &gradP) const
{
   SX_EXIT;
}

/** Rotation mapping class.

  This is used for symmetrizing the block-density matrix.

  */
class SxHubbardFO::SymMap : public SxHubbardXO::SymMap
{
   private:
      // mapped site :iSite,iSym
      SxArray<SxVector<int> > mapSiteSym;
      // mapped orbital type :iOrbType, iSym
      SxArray<SxArray2<int> > mapOrbSym;
      // offset of orbital type iOrbType in total orbital list
      SxArray<SxArray<int> > orbOffset;
      // l of orbital type iOrbType
      SxArray<SxArray<int> > lType;
   public:
   /// Constructor
      SymMap (const SxHubbardFO &parent, const SxAtomicStructure &structure);

      /// Get rotated site id upon rotation iSym
      virtual int rotateSite (int iSite, int iSym)
      {
         return mapSiteSym(iSite)(iSym);
      }
      /// Get rotated orbital type of site iSite upon rotation iSym
      virtual int rotatedType (int iSite, int iOrbType, int iSym)
      {
         return mapOrbSym(iSite)(iOrbType, iSym);
      }
      /// Get number of orbital types (across all atoms) for one site
      virtual int getNOrbType (int iSite)
      {
         return (int)orbOffset(iSite).getSize ();
      }
      /// Get the reference site (smallest equivalent site id)
      virtual int getRefSite (int iSite)
      {
         return mapSiteSym(iSite).minval ();
      }
      /// Get orbital offset for particular type at particular site
      virtual int getOffset (int iSite, int iOrbType)
      {
         return orbOffset(iSite)(iOrbType);
      }
      /// Get l-channel for an orbital type
      virtual int getL (int iSite, int iOrbType)
      {
         return lType(iSite)(iOrbType);
      }

      /// Virtual destructor
      virtual ~SymMap () = default;
};

SxHubbardFO::SymMap::SymMap (const SxHubbardFO &parent,
                             const SxAtomicStructure &structure)
{
   // --- identify orbital types
   const SxArray<SxOrbitalIndex> &orbs = parent.aoProj->orbitalMap;
   const SxArray<SxArray<SxAoIndex> > &refOrbMap = parent.aoProj->refOrbMap;

   int nSite = parent.getNSite ();
   SxArray<int> nOrbType (nSite);
   orbOffset.resize (nSite);
   lType.resize (nSite);
   SX_LOOP(iSite)  {
      nOrbType(iSite) = 0;
      int offset = parent.siteStart (iSite);
      // count orbital types
      for (int io = 0; io < parent.nFragOrb(iSite); ++io)  {
         const SxOrbitalIndex &orb = orbs(parent.mapFragOrbProj(offset + io));
         if (refOrbMap(orb.is)(orb.io).m == 0)  {
            nOrbType(iSite)++;
         }
      }
      // find orbital offsets and l-value
      orbOffset(iSite).resize (nOrbType(iSite));
      lType(iSite).resize (nOrbType(iSite));
      for (int io = 0, iot = 0; io < parent.nFragOrb(iSite); ++io)  {
         const SxOrbitalIndex &orb = orbs(parent.mapFragOrbProj(offset + io));
         const SxAoIndex &nlm = refOrbMap(orb.is)(orb.io);
         if (nlm.m == -nlm.l)  {
            lType(iSite)(iot) = nlm.l;
            orbOffset(iSite)(iot++) = io;
         }
      }
   }

   // --- set up symmetry mapping for atoms
   const SxSymGroup &syms = *structure.cell.symGroupPtr;
   SxGrid grid (structure, 10);
   int nSym = syms.getNSymmorphic ();

   mapSiteSym.resize (nSite);
   mapOrbSym.resize (nSite);
   SX_LOOP(iSite)  {
      mapSiteSym(iSite).resize (nSym);
      mapOrbSym(iSite).reformat (nOrbType(iSite), nSym);
      int offset = parent.siteStart (iSite);
      Coord center = structure.getAtom (parent.centralAtoms(iSite));
      // --- setup symmetry mappings for all symmetries
      for (int iSym = 0; iSym < nSym; iSym++)  {
         // --- rotate center and identify rotated site
         const SymMat &S = syms.getSymmorphic (iSym);
         Coord rotCenter = S ^ center;
         int iaRot = structure.find (rotCenter, grid);
         SX_CHECK (iaRot >=0);
         int iSiteRot = (int)parent.centralAtoms.findPos (iaRot);
         mapSiteSym(iSite)(iSym) = iSiteRot;
         int offRot = parent.siteStart(iSiteRot);

         // --- loop over atoms within site
         int iot = 0;
         for (int ia = 0; ia < parent.getNAtoms(iSite); ++ia) {
            int iTl = parent.atomIdx(iSite)(ia);
            int is = structure.getISpecies (iTl);
            int nOrbTypeAtom = parent.refOrbsProjG.getNOrbTypes (is);
            // --- find rotated atom
            Coord posRel = structure.getAtom (iTl)
                         + parent.latShift(offset + orbOffset(iSite)(iot))
                         - center;
            Coord rotPos = (S ^ posRel)
                         + structure.getAtom (parent.centralAtoms(iSiteRot));
            // brute force search, since we must respect latShift
            int jot = 0;
            for (int ja = 0; ja < parent.getNAtoms(iSite); ++ja)  {
               int jTl = parent.atomIdx(iSiteRot)(ja);
               int js = structure.getISpecies (jTl);
               if (is == js)  {
                  Coord pos = structure.getAtom(jTl)
                            + parent.latShift(offRot + orbOffset(iSiteRot)(jot));
                  if ((pos - rotPos).normSqr () < 1e-6) break;
               }
               jot += parent.refOrbsProjG.getNOrbTypes (js);
            }
            SX_CHECK (jot < nOrbType(iSite), jot, nOrbType(iSite));
            // --- set up symmetry map for orbital types for that atom
            for (int i = 0; i < nOrbTypeAtom; i++)
               mapOrbSym(iSite)(iot++, iSym) = jot++;
         }
         SX_CHECK (iot == nOrbType(iSite), iot, nOrbType(iSite));
      }
   }
}

/// Get symmetry map
SxPtr<SxHubbardXO::SymMap>
SxHubbardFO::getSymMap (const SxAtomicStructure &structure) const
{
   return SxPtr<SxHubbardFO::SymMap>::create (*this, structure);
}

