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

#include <SxHubbardXO.h>
#include <SxTextIO.h>

SxVector<double>
SxHubbardXO::truncateShapeG (const SxVecRef<double> &shape,
                             double rCut,
                             double width) const
{
   const SxRadialBasis *radG
      = dynamic_cast<const SxRadialBasis*>(shape.getBasisPtr ());
   SX_CHECK (radG);
   SxRadialBasis R(0., rCut + 6. * width, nRad, true);
   const SxVecRef<double> &r = R.getRadFunc ();

   SxVector<double> shapeR = R | shape;
   if (verbose) {
     SxTextIO(prefix + "shape.dat").writeXYPlot (r, shapeR);
   }

   // ref. 1, Eq. (19)
   SX_LOOP(i)
      shapeR(i) *=  0.5 * erfc((r(i) - rCut)/width);

   if (verbose) {
      SxTextIO(prefix + "proj.dat").writeXYPlot (r, shapeR);
   }

   return (*radG) | shapeR;
}

void SxHubbardXO::symmetrize (SxBlockDensityMatrix *Pij,
                              const SxAtomicStructure &structure,
                              const SxYlmRotGroup &ylmRot) const
{
   SX_CHECK (Pij->getNSite () == getNSite (),
             Pij->getNSite (), getNSite ());
   // --- TODO
   // symmetrization function in SxBlockDensityMatrix
   // sym (iSite, L1, offset1, L2, offset2, ylmRot, mapSite, mapSpin)
   // => symmetrize subblock 2L1+1 x 2L2+1 (at offset offset1,offset2)
   //    across all equivalent sites, rotating the sites and potentially
   //    the spin according to the maps

   SxPtr<SymMap> symMap = getSymMap (structure);

   SX_CHECK (structure.cell.symGroupPtr);
   int nSym = structure.cell.symGroupPtr->getNSymmorphic ();

   for (int iSite = 0; iSite < getNSite (); ++iSite) {
      // check if symmetrization has been performed already
      if (symMap->getRefSite(iSite) < iSite) continue;

      for (int iSpin = 0; iSpin < Pij->getNSpin (); ++iSpin)  {
         for (int iot = 0; iot < symMap->getNOrbType (iSite); ++iot)  {
            int l1 = symMap->getL (iSite, iot), nm1 = 2 * l1 + 1;
            for (int jot = 0; jot < symMap->getNOrbType (iSite); ++jot)  {
               int l2 = symMap->getL (iSite, jot), nm2 = 2 * l2 + 1;
               SxVector<double> sym (nm1, nm2), rotRho(nm1,nm2);
               sym.set (0.);
               // --- construct symmetrized molDensMat for this site
               for (int iSym = 0; iSym < nSym; ++iSym)  {
                  int rotSite = symMap->rotateSite (iSite, iSym);
                  int iotR = symMap->rotatedType(iSite, iot, iSym);
                  int jotR = symMap->rotatedType(iSite, jot, iSym);
                  int offset  = symMap->getOffset (rotSite, iotR),
                      offset2 = symMap->getOffset (rotSite, jotR);
                  // --- collect rotRho
                  const SxVector<double> &rotD = (*Pij)(iSpin,rotSite);
                  const SxVector<double> &Dl1  = ylmRot(iSym)(l1),
                                         &Dl2  = ylmRot(iSym)(l2);
                  for (int m = 0; m < nm1; ++m)
                     for (int m2 = 0; m2 < nm2; ++m2)
                        rotRho(m,m2) = rotD(offset + m, offset2 + m2);
                  // --- symmetrize
                  // ref. 1, Eq. (45)
                  sym += Dl1.transpose () ^ rotRho ^ Dl2;
               }
               sym /= double(nSym);
               // --- distribute symmetrized matrix
               for (int iSym = 0; iSym < nSym; ++iSym)  {
                  // rotate with symmetry
                  const SxVector<double> &Dl1  = ylmRot(iSym)(l1),
                                         &Dl2  = ylmRot(iSym)(l2);
                  // ref. 1, Eq. (46)
                  rotRho = Dl1 ^ sym ^ Dl2.transpose ();

                  int rotSite = symMap->rotateSite (iSite, iSym);
                  int iotR = symMap->rotatedType(iSite, iot, iSym);
                  int jotR = symMap->rotatedType(iSite, jot, iSym);
                  int offset  = symMap->getOffset (rotSite, iotR),
                      offset2 = symMap->getOffset (rotSite, jotR);
                  // write out
                  SxVector<double> &rotD = (*Pij)(iSpin,rotSite);
                  for (int m = 0; m < nm1; ++m)
                     for (int m2 = 0; m2 < nm2; ++m2)
                        rotD(offset + m, offset2 + m2) = rotRho(m,m2);
               }
            }
         }
      }
   }
}

void SxHubbardXO::symmetrizeGradPij (SxBlockDensityMatrix *Pij,
                                     const SxAtomicStructure &structure,
                                     const SxYlmRotGroup &ylmRot) const
{
   SX_CHECK (Pij->getNSite () == getNSite (),
             Pij->getNSite (), getNSite ());
   // note: for gradients, we use the "spin" index of SxBlockDensityMatrix
   // for the direction
   SX_CHECK(Pij->getNSpin () == 3, Pij->getNSpin ());

   SxPtr<SymMap> symMap = getSymMap (structure);

   SX_CHECK (structure.cell.symGroupPtr);
   const SxSymGroup &syms = *structure.cell.symGroupPtr;
   int nSym = syms.getNSymmorphic ();

   for (int iSite = 0; iSite < getNSite (); ++iSite) {
      // check if symmetrization has been performed already
      if (symMap->getRefSite(iSite) < iSite) continue;

      SxArray<SxVector<double> > sym(3);
      for (int iot = 0; iot < symMap->getNOrbType (iSite); ++iot)  {
         int l1 = symMap->getL (iSite, iot), nm1 = 2 * l1 + 1;
         for (int jot = 0; jot < symMap->getNOrbType (iSite); ++jot)  {
            int l2 = symMap->getL (iSite, jot), nm2 = 2 * l2 + 1;
            SxVector<double> rotRho(nm1,nm2);
            SX_LOOP(iDir)  {
               sym(iDir).reformat (nm1, nm2);
               sym(iDir).set (0.);
            }
            for (int iSym = 0; iSym < nSym; ++iSym)  {
               int rotSite = symMap->rotateSite (iSite, iSym);
               int iotR = symMap->rotatedType(iSite, iot, iSym);
               int jotR = symMap->rotatedType(iSite, jot, iSym);
               int offset  = symMap->getOffset (rotSite, iotR),
                   offset2 = symMap->getOffset (rotSite, jotR);
               const SymMat &S = syms.getSymmorphic (iSym);
               // --- collect rotRho
               const SxVector<double> &Dl1  = ylmRot(iSym)(l1),
                                      &Dl2  = ylmRot(iSym)(l2);
               for (int iDir = 0; iDir < 3; ++iDir)  {
                  const SxVector<double> &rotD = (*Pij)(iDir,rotSite);
                  for (int m = 0; m < nm1; ++m)
                     for (int m2 = 0; m2 < nm2; ++m2)
                        rotRho(m,m2) = rotD(offset + m, offset2 + m2);
                  // ref. 1, Eq. (47), sum over kk'
                  rotRho = Dl1.transpose () ^ rotRho ^ Dl2;
                  // ref. 1, Eq. (47), sum over beta
                  for (int jDir = 0; jDir < 3; ++jDir)
                     sym(jDir).plus_assign_ax (S(iDir,jDir), rotRho);
               }
            }
            SX_LOOP(iDir) sym(iDir) /= double(nSym);
            // --- distribute symmetrized matrix
            for (int iSym = 0; iSym < nSym; ++iSym)  {
               // rotate with symmetry
               const SymMat &S = syms.getSymmorphic (iSym);
               const SxVector<double> &Dl1  = ylmRot(iSym)(l1),
                                      &Dl2  = ylmRot(iSym)(l2);
               SX_LOOP(iDir)  {
                  rotRho.set (0.);
                  // ref. 1, Eq. (48), sum over beta
                  SX_LOOP(jDir)
                     rotRho.plus_assign_ax (S(iDir,jDir), sym(jDir));

                  // ref. 1, Eq. (48), sum over kk'
                  rotRho = Dl1 ^ rotRho ^ Dl2.transpose ();

                  int rotSite = symMap->rotateSite (iSite, iSym);
                  int iotR = symMap->rotatedType(iSite, iot, iSym);
                  int jotR = symMap->rotatedType(iSite, jot, iSym);
                  int offset  = symMap->getOffset (rotSite, iotR),
                      offset2 = symMap->getOffset (rotSite, jotR);
                  // write out
                  SxVector<double> &rotD = (*Pij)(iDir,rotSite);
                  for (int m = 0; m < nm1; ++m)
                     for (int m2 = 0; m2 < nm2; ++m2)
                        rotD(offset + m, offset2 + m2) = rotRho(m,m2);
               }
            }
         }
      }
   }
}

SxAtomicStructure
SxHubbardXO::getForce (const SxBlockDensityMatrix &gradPij,
                       ssize_t iSpin) const
{
   SxAtomicStructure f2 = trafoForce.getNewStr ();
   f2.set (0., 0., 0.);
   SX_LOOP2(iDir, iSite)  {
      for (int ia = 0, io = 0; ia < getNAtoms(iSite); ++ia)  {
         int iTl = getIAtom(iSite, ia);
         for (int iAO = 0; iAO < getNOrbAtom(iSite, ia); iAO++, ++io)  {
            // ref. 1, Eq. (40)
            f2.ref(iTl)(iDir) += dot(gradPij(iDir,iSite).rowRef(io),
                                       hamProj(iSpin,iSite).rowRef(io)
                                     + hamProj(iSpin,iSite).colRef(io));
         }
      }
   }
   return f2;
}
