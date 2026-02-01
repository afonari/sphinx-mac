//
//      Contact:    Sixten Boeck, boeck@mpie.de
//                  Algorithm Design and Modeling Group
//                  Computational Materials Design
//                  Max-Planck-Institute for Iron Research
//                  40237 Duesseldorf, Germany
//
//      Authors:    see src/AUTHORS
//
// ---------------------------------------------------------------------------

#ifndef _SX_RADIAL_BASIS_H_
#define _SX_RADIAL_BASIS_H_

#include <SxPrecision.h>
#include <SxQuantumNumbers.h>
#include <SxBasis.h>
#include <SxVector.h>
#include <SxArray.h>
#include <SxTimer.h>
#include <SxDiracLib.h>
#include <SxYlm.h>


class SX_EXPORT_DIRAC SxRadialBasis : public SxBasis
{
   public:
      /// Standard coefficient type
      typedef double              CoeffType;
      /// Standard vector type
      typedef SxVector<CoeffType> TPsi;

      /// If true, this is radial R. If false, this is radial G.
      bool realSpace;

      /// \brief Type of the grid
      enum gridType {
         /// \f$ r_i = r_0 + i \cdot \frac{r_{\rm max} - r_0}{N-1} \f$
         Linear,
         /// \f$ r_i = r_0 + \sqrt i \cdot \frac{r_{\rm max} - r_0}{\sqrt{N-1}} \f$
         Quadratic,
         /// quadratic grid for the first 1/3 of the total range; other 2/3 are linear
         Mixed,
         /// \f$ r_i = r_0 + (r_{\rm max} - r_0) * \frac{2i}{i + N-1} \f$
         Hyperbel,
         /// \f$ r_i = r_0 * e^{\lambda i}\f$ with \f$ \lambda = \frac{1}{N-1}\ln \frac{r_{\rm max}}{r_0}\f$
         Logarithmic};

      /// Empty constructor
      SxRadialBasis ();
      /// Destructor
      ~SxRadialBasis ();
      /** \brief Constructor
        @param radMin     value for first grid point \f$r_0\f$
        @param radMax     value for last grid point \f$r_{\rm max}\f$
        @param nPoints    number of points in grid
        @param realSpace_ specify where this basis lives: real-space (true)
                          or reciprocal space (false)
        @param mode       type of grid
      */
      SxRadialBasis (double radMin, double radMax, int nPoints, bool realSpace_, enum gridType mode = Linear);
      /** \brief Constructor
        @param basis      explicit grid points
        @param realSpace_ specify where this basis lives: real-space (true)
                          or reciprocal space (false)
      */
      SxRadialBasis (const SxVector<CoeffType> &basis, bool realSpace_);
      SxRadialBasis (const SxVecRef<CoeffType> &dRadIn, double factor, double rad0, bool realSpace_);

      /** \brief Set basis
        @param radMin     value for first grid point \f$r_0\f$
        @param radMax     value for last grid point \f$r_{\rm max}\f$
        @param nPoints    number of points in grid
        @param mode       type of grid

        \note realSpace must be set separately!
      */
      void set (double radMin, double radMax, int nPoints, enum gridType mode);
      /** \brief Set basis
        @param basis      explicit grid points

        \note realSpace must be set separately!
      */
      void set (const SxVector<CoeffType> &basis);
      void set (const SxVecRef<CoeffType> &dRadIn, double factor, double rad0);

      /// Get all the grid points
      const SxVector<CoeffType> & getRadFunc () const
      {
         return rad;
      }
      /// Get minimum radial grid value
      inline double getRadMin () const { return rad(0); }
      /// Get maximum radial grid value
      inline double getRadMax () const { return rad(rad.getSize () - 1); }
      /** \brief Get the grid differential \f$ r_{i+1}-r_i \f$

        \note the last differential is always set to zero
      */
      const SxVector<CoeffType> & getRadDR () const
      {
         return dRad;
      }

      /// Get type of basis
      virtual SxString getType () const {
         return realSpace ? "|r>" : "|g>";
      }

      /// Get the number of grid points
      virtual ssize_t getNElements () const { return rad.getSize (); }

      REGISTER_PROJECTOR (SxRadialBasis, radToRad);
      REGISTER_PROJECTOR (SxRadBasis, toRadBasis);
      REGISTER_PROJECTOR (SxGBasis, toGBasis);

      /** \brief Transformation routine

        \note At present, only identity (input = output) and Fourier transforms are supported
        */
      SxVecRef<CoeffType> radToRad (const SxRadialBasis *basisPtr,
                                    const SxVecRef<CoeffType> &vec) const;
      /// Transform to real-space radial basis
      SxVector<double> toRadBasis (const SxRadBasis *radBasisPtr,
                                   const SxVecRef<double> &vec) const;
      /// Transform between real and reciprocal space via spherical-Bessel transform
      SxVector<double> fourier (const SxRadialBasis *radBasisPtr,
                                const SxVecRef<double> &vec) const;
      /// Project from radial G to full G via interpolation
      SxVector<PrecCoeffG> toGBasis (const SxGBasis *gBasisPtr,
                                     const SxVecRef<CoeffType> &vec) const;

      /** \brief Projector from \b r space to some anonymous basis

          This is the versatile interface: we try to dynamically determine
          the basis to project to.

          \sa \ref page_dirac */
      virtual SxVecRef<double>
      projectTo (const SxBasis *, const SxVecRef<double> &) const override;
      // int(... R^2dR)
      double integrate (const SxVecRef<CoeffType> &integrand, bool useSpline = true) const;
      /// 4-point Simpson integration
      inline double simpsonTerm (double x1, double x2, double x3, double x4, double f1, double f2, double f3, double f4) const;
      /// 3-point Simpson integration
      inline double simpsonTerm (double x1, double x2, double x3, double f1, double f2, double f3) const;
      /// 2-point integration (trapezoidal rule)
      inline double simpsonTerm (double x1, double x2, double f1, double f2) const;

      /// Trace implementation: integrate over radial space
      virtual Real8 tr(const SxVecRef<double> &x) const
      {
         return integrate (x, false);
      }


      /// Compute spline coefficients
      SxVector<double> toSpline (const SxVecRef<double> &vec) const;
      /// Interpolate spline coefficients to the basis
      SxVector<double> toVec (const SxVecRef<double> &vec) const;

      /// Spherical Bessel function
      static SxVector<double> jsb (int l, const SxVecRef<double> &z);

      bool isLogarithmic () const
      {
          SX_CHECK (rad.getSize () > 0);
          if (rad(0) < 1e-14) return false;
          double logdr1 = log(rad(1)/rad(0));
          double logdrN = log(rad(rad.getSize () - 1)/rad(0))
                        / double(rad.getSize () - 1);
          return fabs(logdr1 - logdrN) < 1e-6;
      }

   protected:

      /// The actual grid points
      SxVector<double> rad;
      /** \brief The differences between grid points \f$r_{i+1} - r_i\f$. 
        \note The last difference for $i=N-1$ is set to zero
      */
      SxVector<double> dRad;

   public:
      /** \brief Get overlap <chi1|chi2> between two AO-like functions

          @param l1,m1 angular quantum numbers of chi1
          @param l2,m2 angular quantum numbers of chi2
          @param matElem matrix elements of \f$ M_L = \int g^2 dg f_1(g) j_L(g r) f_2(g) \f$,
                         only L=|l1-l2| .. |l1+l2| in steps of 2 are needed
          @param Ylm   precomputed Y_lm(r1 - r2), orientation matters (not r2-r1),
                       must include the Y_lm normalization factor
          @param clebsch Clebsch-Gordan coefficients until l1+l2
          @return the overlap integral

          Assuming that
          \f[ \chi_i(\mathbf g) = f_i(|\mathbf g|) Y_{l_im_i}(\hat \mathbf g) \f]
          This is computed as
          \f[
          \langle \chi_1 | \chi_2 \rangle
             = \sum_L M_L \langle l_1m_1\, LM|l_2m_2\rangle Y_{LM}(\hat \mathbf r_{21})
          \f]
          with \f$\mathbf r_{21} = \mathbf r_1 - \mathbf r_2\f$
          and \f$ M_L = \int g^2 dg f_1(g) j_L(g r) f_2(g) \f$.
        */
      static double getOverlap(int l1, int m1, int l2, int m2,
                               const SxVecRef<double> &matElem,
                               const SxVecRef<double> &Ylm,
                               const SxYlm::SxClebschTable &clebsch);

      /** \brief Get reduced matrix element for overlap <chi1|chi2>
                 between two AO-like functions, see #getOverlap

          @param f1    radial part in g-space of chi1
          @param f2    radial part in g-space of chi2
          @param jsb   precomputed Bessel functions \f$j_L(gr_{21})\f$ for
                       \f$|l_1-l_2|\dots(l_1 + l_2)\$f, only
                       every second L is needed, i.e. |l1-l2|+0,2,4,6...
          \f$r_21\f$ is the distance between the atoms.

          @return \f$ M_L = \int g^2 dg f_1(g) j_L(g r_{21}) f_2(g) \f$
          */
      static SxVector<double> computeMatElem(const SxVecRef<double> &f1,
                                             const SxVecRef<double> &f2,
                                             const SxArray<SxVector<double> > &jsb);
      /** \brief Compute array of spherical Bessel functions \f$j_l (g r)\f$.

          @param lmax maximum l value
          @param dist distance, i.e., scaling factor r if basis is g.
          @return Bessel functions \f$j_L(gr_{21})\f$
      */
      SxArray<SxVector<double> > computeJsb (int lmax, double dist);
};

namespace Timer {
   enum RadialBasisTimer {
      RadFourier, RadR2Rad, RadG2Rad, RadG2Gk, IntegrateRad,
      VecToSpline, SplineToVec, PhaseFactors, SplineGetY
   };
}

SX_REGISTER_TIMERS (Timer::RadialBasisTimer)
{
   using namespace Timer;
   regTimer (RadFourier,   "<g|r> or <r|g>");
   regTimer (RadR2Rad,     "RadR to Rad Basis");
   regTimer (RadG2Rad,     "RadG to Rad Basis");
   regTimer (RadG2Gk,      "RadG to Gk Basis");
   regTimer (IntegrateRad, "radial Integration");
   regTimer (VecToSpline,  "Vector to Spline");
   regTimer (SplineToVec,  "Spline to Vector");
   regTimer (PhaseFactors, "calc phasefactors");
   regTimer (SplineGetY,   "Spline to |G+k|");
}

#endif /* _SX_RADIAL_BASIS_H_ */
