/**
***************************************************************************
* @file brick/geometry/bullseye2D_impl.hh
*
* Source file defining the Bullseye2D class.
*
* Copyright (C) 2008-2013 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_GEOMETRY_BULLSEYE2D_IMPL_HH
#define BRICK_GEOMETRY_BULLSEYE2D_IMPL_HH

// This file is included by bullseye2D.hh, and should not be directly included
// by user code, so no need to include bullseye2D.hh here.
// 
// #include <brick/geometry/bullseye2D.hh>

#include <brick/linearAlgebra/linearAlgebra.hh>
#include <brick/common/constants.hh>
#include <brick/common/mathFunctions.hh>
#include <brick/numeric/array2D.hh>
#include <brick/numeric/maxRecorder.hh>
#include <brick/numeric/utilities.hh>
#include <brick/numeric/vector2D.hh>

namespace brick {

  namespace geometry {
    
    // The default constructor initializes to the unit circle.
    template <class Type>
    Bullseye2D<Type>::
    Bullseye2D()
      : m_origin(0.0, 0.0),
        m_semimajorAxis(1.0, 0.0),
        m_semiminorAxis(0.0, 1.0),
        m_scales(3)
    {
      m_scales[0] = 0.5;
      m_scales[1] = 1.0;
      m_scales[1] = 1.5;
    }

    
    // Construct a bullseye, explicitly setting parameters.
    template <class Type>
    template <class Iter>
    Bullseye2D<Type>::
    Bullseye2D(Ellipse2D<Type> const& ellipse,
               Iter scalesBegin, Iter scalesEnd)
      : m_origin(ellipse.getOrigin()),
        m_semimajorAxis(ellipse.getSemimajorAxis()),
        m_semiminorAxis(ellipse.getSemiminorAxis()),
        m_scales(scalesEnd - scalesBegin)
    {
      std::copy(scalesBegin, scalesEnd, m_scales.begin());
    }

    
    // The copy constructor deep copies its argument.
    template <class Type>
    Bullseye2D<Type>::
    Bullseye2D(Bullseye2D<Type> const& source)
      : m_origin(source.m_origin),
        m_semimajorAxis(source.m_semimajorAxis),
        m_semiminorAxis(source.m_semiminorAxis),
        m_scales(source.m_scales)
    {
      // Empty.
    }


    // The assignment operator deep copies its argument.
    template <class Type>
    Bullseye2D<Type>&
    Bullseye2D<Type>::
    operator=(Bullseye2D<Type> const& source)
    {
      if(&source != this) {
        m_origin = source.m_origin;
        m_semimajorAxis = source.m_semimajorAxis;
        m_semiminorAxis = source.m_semiminorAxis;
        m_scales = source.m_scales;
      }
      return *this;
    }


    // Estimate bullseye parameters from a series of points on the
    // bullseye.
    template <class Type>
    template<class PointsIterType, class CountsIterType>
    Type
    Bullseye2D<Type>::
    estimate(PointsIterType pointsBeginIter, PointsIterType pointsEndIter,
             CountsIterType countsBeginIter, CountsIterType countsEndIter,
             bool computeResidual)
    {
      // This line may be O(N), depending on what type of iterator.
      unsigned int numberOfPoints = pointsEndIter - pointsBeginIter;

      // Do the estimation.
      brick::numeric::Array1D<Type> residualVector(numberOfPoints);
      this->estimate(pointsBeginIter, pointsEndIter,
                     countsBeginIter, countsEndIter,
                     residualVector.begin(), computeResidual);

      // Compute RMS residual.
      Type residual = static_cast<Type>(0);
      residual = brick::numeric::dot<Type>(residualVector, residualVector);
      residual = brick::common::squareRoot(residual);
      return residual;
    }


    // This member function duplicates the functionality of the other
    // estimate() member function, but returns algebraic residuals for
    // each input point.
    template <class Type>
    template<class PointsIterType, class CountsIterType, class ResidualIter>
    void
    Bullseye2D<Type>::
    estimate(PointsIterType pointsBeginIter, PointsIterType pointsEndIter,
             CountsIterType countsBeginIter, CountsIterType countsEndIter,
             ResidualIter residualIter, bool computeResidual)
    {
      unsigned int const numberOfRings = countsEndIter - countsBeginIter;
      if(numberOfRings == 0) {
        BRICK_THROW(brick::common::ValueException,
                    "Bullseye2D::estimate()",
                    "Sequence of counts must have at least one element.");
      }

      // We assume each ring of the bullseye is elliptical, with the
      // implicit parameterization
      // 
      //   F(x, y) = a * x^2 + b * x * y + c * y^2 + d * x + e * y + f = 0,
      //
      // with the ellipse-specific constraint
      //
      //   b^2 - 4 * a * c < 0.
      //
      // All rings share the same a, b, c, d, and e parameters.  Each
      // ring has its own f, which controls its scaling.  So, the
      // entire bullseye can be described by the parameter vector p =
      // [a, b, c, d, e, f_0, f_1, f_2, ...]^T.
      // 
      // Note that the ellipse parameterization can be scaled
      // arbitrarily, allowing us to restate the constraint
      //
      //   4 * a * c - b^2 = 1.
      //
      // These equations can be written in matrix form
      //
      //   D * p = 0; p^T * C * p = 1,
      //
      // D is the N by (5 + number-of-rings) "design matrix,"
      //
      // @code
      //   N = [[x_0^2, x_0 * y_0, y_0^2, x_0, y_0, 1, 0, 0, ...],
      //        [x_1^2, x_1 * y_1, y_1^2, x_1, y_1, 1, 0, 0, ...],
      //        ...
      //        [x_k^2, x_k * y_k, y_k^2, x_k, y_k, 0, 1, 0, ...],
      //        ...
      //        [x_m^2, x_m * y_m, y_m^2, x_m, y_m, 0, 0, 1, ...],
      //        ...
      //        ],
      // @endcode
      //
      // where k is the index of the first point in the second ring,
      // and m is the index of the first point in the third ring.
      // 
      // C is the (5 + number-of-rings) x (5 + number-of-rings) 
      // constraint matrix, identically zero except for a -1 in the
      // second element of the second row, and 2s in the third element
      // of the first row and the first element of the third row.
      // 
      // Minimizing F(x, y) in the least squares sense, we have
      //
      //   pHat = min_over_p(p^T * S * p)
      //
      // with the constraint that
      //
      //   p^T * C * p = 1,
      //
      // where the "scatter matrix" S = D^T * D.
      //
      // This is readily solved using Lagrange multipliers.  Here, we
      // decompose the solution following Halir & Flasser (see
      // reference [1] in the documentation of this function).

      // This line may be O(N), depending on what type of iterator.
      unsigned int numberOfPoints = pointsEndIter - pointsBeginIter;

      // Left half of the design matrix.
      brick::numeric::Array2D<Type> D1(numberOfPoints, 3);

      // Right half of the design matrix.
      brick::numeric::Array2D<Type> D2(numberOfPoints, 2 + numberOfRings);

      // Fill in the design matrix.
      unsigned int rowNumber = 0;
      unsigned int ringNumber = 0;
      unsigned int ringEndIndex = *countsBeginIter;
      while(pointsBeginIter != pointsEndIter) {
        // Keep track of which ring this input point belongs to, so we
        // can fill in the design matrix appropritely.
        while(rowNumber >= ringEndIndex) {
          ++ringNumber;
          ++countsBeginIter;
          if(!(countsBeginIter != countsEndIter)) {
            BRICK_THROW(brick::common::IndexException,
                        "Bullseye2D::estimate()",
                        "There are more input points than specified by the "
                        "countsBeginIter - countsEndIter input sequence.");
          }
          ringEndIndex += *countsBeginIter;
        }

        // The first five columns of the design matrix are easy.
        brick::numeric::Vector2D<Type> samplePoint = *pointsBeginIter;
        D1(rowNumber, 0) = samplePoint.x() * samplePoint.x();
        D1(rowNumber, 1) = samplePoint.x() * samplePoint.y();
        D1(rowNumber, 2) = samplePoint.y() * samplePoint.y();
        D2(rowNumber, 0) = samplePoint.x();
        D2(rowNumber, 1) = samplePoint.y();

        // The last few depend on which ring the input pointn belongs to.
        for(unsigned int jj = 0; jj < numberOfRings; ++jj) {
          D2(rowNumber, 2 + jj) = Type(0);
        }
        D2(rowNumber, 2 + ringNumber) = Type(1);

        ++pointsBeginIter;
        ++rowNumber;
      }

      // Sanity check the point counts to make sure everything looks
      // right.
      if(rowNumber != ringEndIndex) {
        BRICK_THROW(brick::common::IndexException,
                    "Bullseye2D::estimate()",
                    "There are fewer input points than specified by the "
                    "countsBeginIter - countsEndIter input sequence.");
      }
      
      // Compute the four quadrants of the scatter matrix.  Actually,
      // the upper right and lower left quadrants are identical, so only
      // three matrices to compute here.
      brick::numeric::Array2D<Type> D1Transpose = D1.transpose();
      brick::numeric::Array2D<Type> S1 =
        brick::numeric::matrixMultiply<Type>(D1Transpose, D1);
      brick::numeric::Array2D<Type> S2 =
        brick::numeric::matrixMultiply<Type>(D1Transpose, D2);
      brick::numeric::Array2D<Type> S3 =
        brick::numeric::matrixMultiply<Type>(D2.transpose(), D2);

      // Rather than solving for the whole parameter vector at once,
      // solve first for just the first three elements.  The remaining
      // elements will be recovered using the matrix T.  Justification
      // for this, and for the next few steps, is in the paper, but
      // isn't reproduced here (don't worry, though, it's just
      // algebra, nothing fancy).
      brick::numeric::Array2D<Type> TT;
      try {
        TT = (Type(-1)
              * brick::numeric::matrixMultiply<Type>(
                brick::linearAlgebra::inverse(S3), S2.transpose()));
      } catch(brick::common::ValueException) {
        BRICK_THROW(brick::common::ValueException,
                    "Bullseye2D::estimate()",
                    "Input points are not sufficient to estimate bullseye.");
      }
      
      // After some algebra, Halir & Flasser reduce this the solution
      // for the first three elements to an eigenproblem with 3x3
      // matrix M.
      brick::numeric::Array2D<Type> M0 =
        S1 + brick::numeric::matrixMultiply<Type>(S2, TT);
      brick::numeric::Array2D<Type> MM(M0.rows(), M0.columns());
      MM.getRow(0).copy(M0.getRow(2) * Type(0.5));
      MM.getRow(1).copy(M0.getRow(1) * Type(-1));
      MM.getRow(2).copy(M0.getRow(0) * Type(0.5));

      // The Lagrange solution for the first three parameters.
      brick::numeric::Array1D< std::complex<Type> > eigenvalues;
      brick::numeric::Array2D< std::complex<Type> > eigenvectors;
      brick::linearAlgebra::eigenvectors(
        MM, eigenvalues, eigenvectors);

      // The eigenvector we want is the one with the minimum positive
      // eigenvalue, but we're not going to look at the eigenvalues
      // directly because if we have a very good fit, that eigenvalue
      // might be slightly negative due to numerical issues.  Instead,
      // we evaluate the constraint a^T * C * a, looking eigenvector
      // that gives a positive value (positive because our original
      // constraint set this quantity to 1, and also because for the
      // correct solution our eigensystem guarantees that the
      // constraint value is equal to the eigenvalue.  The claim is
      // that there will only ever be one eigenvector with a positive
      // constraint value, although I have not verified this.
      // Assuming it to be true, we simply take the eigenvector with
      // the largest constraint value.  Another concern, I currently
      // assume all eigenvalues/eigenvectors are real.
      brick::numeric::MaxRecorder<Type, unsigned int> maxRecorder;
      for(unsigned int ii = 0; ii < eigenvalues.size(); ++ii) {
        Type constraintValue =
          4.0 * eigenvectors(0, ii).real() * eigenvectors(2, ii).real()
          - eigenvectors(1, ii).real() * eigenvectors(1, ii).real();
        maxRecorder.test(constraintValue, ii);
      }

      unsigned int selectedColumn = maxRecorder.getPayload();
      brick::numeric::Array1D<Type> firstThreeParameters(eigenvectors.rows());
      for(unsigned int ii = 0; ii < eigenvectors.rows(); ++ii) {
        firstThreeParameters[ii] = eigenvectors(ii, selectedColumn).real();
      }

      // Solve for the remaining parameters using matrix T, as
      // mentioned above.
      brick::numeric::Array1D<Type> remainingParameters =
        brick::numeric::matrixMultiply<Type>(TT, firstThreeParameters);

      // Combine all into one vector.
      brick::numeric::Array1D<Type> allParameters(
        firstThreeParameters.size() + remainingParameters.size());
      std::copy(firstThreeParameters.begin(), firstThreeParameters.end(),
                allParameters.begin());
      std::copy(remainingParameters.begin(), remainingParameters.end(),
                allParameters.begin() + firstThreeParameters.size());

      // Convert into the parameterization we prefer.
      this->convertAlgebraicToTrigonometric(
        allParameters, m_origin, m_semimajorAxis, m_semiminorAxis,
        m_scales);

      // Compute the algebraic residuals.
      if(computeResidual) {
        brick::numeric::Array1D<Type> residualTerm1 =
          brick::numeric::matrixMultiply<Type>(D1, firstThreeParameters);
        brick::numeric::Array1D<Type> residualTerm2 =
          brick::numeric::matrixMultiply<Type>(D2, remainingParameters);
        brick::numeric::Array1D<Type> residualVector =
          residualTerm1 + residualTerm2;
        std::copy(residualVector.begin(), residualVector.end(), residualIter);
      }
    }


    // Returns an ellipse that describes each of the rings of the
    // bullseye.
    template <class Type>
    Ellipse2D<Type>
    Bullseye2D<Type>::
    getEllipse()
    {
      return Ellipse2D<Type>(m_origin, m_semimajorAxis, m_semiminorAxis);
    }

    
    // Convert from implicit bullseye representation to trigonometric
    // parameters.
    template <class Type>
    void
    Bullseye2D<Type>::
    convertAlgebraicToTrigonometric(
      brick::numeric::Array1D<Type> const& algebraicParameters,
      brick::numeric::Vector2D<Type>& origin,
      brick::numeric::Vector2D<Type>& semimajorAxis,
      brick::numeric::Vector2D<Type>& semiminorAxis,
      std::vector<Type>& scales)
    {
      // This conversion is largely copied from an article on ellipse
      // parameterization written by Robert J. Lopez, and published on
      // the maplesoft tips & techniques page.  We extend it slightly
      // to handle bullseyes, and use it without proof here, although
      // R. J. Lopez's article includes a derivation, which follows
      // these steps: substitute u = (x - h) and v = (y - k) into the
      // algebraic distance equation.  When h and k are chosen so that
      // the ellipse is centered at the origin, then the algebraic
      // distance equation will have no coefficients for the linear
      // terms in u and v.  Setting these coefficients to zeros and
      // solving for h and k gives the origin.
      //
      // Next, substituting for h and k to get a quadratic with no
      // linear terms, we have the equation for a non-axis-aligned
      // ellipse at the origin.  The constant element of this
      // polynomial happens to look like a determinant over a simple
      // quadratic, but that's not really important.  Lopez only has
      // to worry about one constant element, because he's dealing
      // with a single ellipse.  Our bullseye looks like several
      // concentric ellipses, so we'll have several of these constant
      // terms, which control the scale of the ellipses.
      //
      // An axis-aligned ellipse has no u*v cross-term.  Again, Lopez
      // makes a variable substitution: u = C*X - S*Y; v = S*X + C*Y,
      // where S and C are sin and cos of a rotation angle.  Isolating
      // the XY term, setting it to zero, and applying double-angle
      // equalities gives sin(2*Theta) / cos(2*Theta) = b / (a - c).
      //
      // This final equation is manupulated to get expressions for
      // sin(Theta) and cos(Theta), which are substituted back into
      // the rotated ellipse equation and use to solve for the axes.

      // Sanity check.
      if(algebraicParameters.size() < 6) {
        BRICK_THROW(brick::common::ValueException,
                    "Bullseye2D::convertAlgebraicToTrigonometric()",
                    "Parameter vector must have at least six elements.");
      }

      unsigned int const numberOfRings = algebraicParameters.size() - 5;

      // First extract parameterization,
      Type aa = algebraicParameters[0];
      Type bb = algebraicParameters[1];
      Type cc = algebraicParameters[2];
      Type dd = algebraicParameters[3];
      Type ee = algebraicParameters[4];

      std::vector<Type> ffVector(numberOfRings);
      for(unsigned int ii = 0; ii < numberOfRings; ++ii) {
        ffVector[ii] = algebraicParameters[5 + ii]; 
      }

      // Compute center of ellipse.
      Type phi = bb * bb - 4 * aa * cc;
      origin.setValue((2 * cc * dd - bb * ee) / phi,
                      (2 * aa * ee - bb * dd) / phi);

      // Compute sine and cosine of ellipse rotation angle.  These
      // define directions of the major and minor axes of the ellipse.
      Type CC = brick::common::constants::rootOverTwo;
      Type SS = brick::common::constants::rootOverTwo;
      if(aa != cc) {
        Type aMinusC = aa - cc;
        Type lambda = (brick::common::absoluteValue(aMinusC)
                       / brick::common::squareRoot(
                         bb * bb + aMinusC * aMinusC));
        CC = brick::common::squareRoot((1.0 + lambda) / 2);
        SS = brick::common::squareRoot((1.0 - lambda) / 2);
      }

      // We'll compute scaling for each ring, then commit to semimajor
      // and semiminor axes later.
      std::vector<Type> alphaVector(numberOfRings);
      std::vector<Type> betaVector(numberOfRings);
      Type PP = aa * CC * CC + cc * SS * SS + bb * CC * SS;
      Type QQ = aa * SS * SS + cc * CC * CC - bb * CC * SS;
      for(unsigned int ringNumber = 0; ringNumber < numberOfRings;
          ++ringNumber) {
        Type gamma = (-2 * aa * ee*ee
                      + 8 * aa * cc * ffVector[ringNumber]
                      + 2 * bb * ee * dd
                      - 2 * ffVector[ringNumber] * bb * bb
                      - 2 * cc * dd * dd);
        Type rho = gamma / (2 * phi);
        alphaVector[ringNumber] = brick::common::squareRoot(rho / PP);
        betaVector[ringNumber] = brick::common::squareRoot(rho / QQ);
      }

      // Now pinck a major and minor axis that are approximately
      // scaled to match the input data.  We'll just take the average
      // of the scales of the rings.
      Type alpha = Type(0);
      Type beta = Type(0);
      for(unsigned int ringNumber = 0; ringNumber < numberOfRings;
          ++ringNumber) {
        alpha += alphaVector[ringNumber] / numberOfRings;
        beta += betaVector[ringNumber] / numberOfRings;
      }
      if(alpha >= beta) {
        semimajorAxis.setValue(alpha * CC, alpha * SS);
        semiminorAxis.setValue(-beta * SS, beta * CC);
      } else {
        semimajorAxis.setValue(-beta * SS, beta * CC);
        semiminorAxis.setValue(alpha * CC, alpha * SS);
      }

      // And record the scales of the the different rings.
      for(unsigned int ringNumber = 0; ringNumber < numberOfRings;
          ++ringNumber) {
        scales[ringNumber] = ((alphaVector[ringNumber] + betaVector[ringNumber])
                              / (alpha + beta));
      }
    }


    /* ======= Non-member functions. ======= */

    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream, const Bullseye2D<Type>& bullseye)
    {
      stream << "Bullseye2D{ "
             << bullseye.getOrigin() << ", "
             << bullseye.getSemimajorAxis() << ", "
             << bullseye.getSemiminorAxis() << " }";
      return stream;
    }
    
  } // namespace geometry
    
} // namespace brick

#endif /* #ifndef BRICK_GEOMETRY_BULLSEYE2D_IMPL_HH */
