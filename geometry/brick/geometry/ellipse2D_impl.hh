/**
***************************************************************************
* @file brick/geometry/ellipse2D_impl.hh
*
* Source file defining the Ellipse2D class.
*
* Copyright (C) 2008 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_GEOMETRY_ELLIPSE2D_IMPL_HH
#define BRICK_GEOMETRY_ELLIPSE2D_IMPL_HH

// This file is included by ellipse2D.hh, and should not be directly included
// by user code, so no need to include ellipse2D.hh here.
// 
// #include <brick/geometry/ellipse2D.hh>

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
    Ellipse2D<Type>::
    Ellipse2D()
      : m_origin(0.0, 0.0),
        m_semimajorAxis(1.0, 0.0),
        m_semiminorAxis(0.0, 1.0)
    {
      // Empty.
    }

    
    // This constructor initializes the ellipse using explicitly
    // specified values.
    template <class Type>
    Ellipse2D<Type>::
    Ellipse2D(brick::numeric::Vector2D<Type> const& origin,
              brick::numeric::Vector2D<Type> const& semimajorAxis,
              Type ratio)
      : m_origin(origin),
        m_semimajorAxis(semimajorAxis),
        m_semiminorAxis(semimajorAxis.y() * ratio, semimajorAxis.x() * ratio)
    {
      if(ratio > 1.0) {
        std::swap(m_semimajorAxis, m_semiminorAxis);
      }
    }

    
    // The copy constructor deep copies its argument.
    template <class Type>
    Ellipse2D<Type>::
    Ellipse2D(Ellipse2D<Type> const& source)
      : m_origin(source.m_origin),
        m_semimajorAxis(source.m_semimajorAxis),
        m_semiminorAxis(source.m_semiminorAxis)
    {
      // Empty.
    }


    // The assignment operator deep copies its argument.
    template <class Type>
    Ellipse2D<Type>&
    Ellipse2D<Type>::
    operator=(Ellipse2D<Type> const& source)
    {
      if(&source != this) {
        m_origin = source.m_origin;
        m_semimajorAxis = source.m_semimajorAxis;
        m_semiminorAxis = source.m_semiminorAxis;
      }
      return *this;
    }


    // Estimate ellipse parameters from a series of points on the
    // ellipse.
    template <class Type>
    template<class IterType>
    void
    Ellipse2D<Type>::
    estimate(IterType beginIter, IterType endIter)
    {
      // This algorithm is based on the ellipse parameterization
      // 
      //   F(x, y) = a * x^2 + b * x * y + c * y^2 + d * x + e * y + f = 0,
      //
      // with the ellipse-specific constraint
      //
      //   b^2 - 4 * a * c < 0.
      //
      // Note that the ellipse parameterization can be scaled
      // arbitrarily, allowing us to restate the constraint
      //
      //   4 * a * c - b^2 = 1.
      //
      // These equations can be written in matrix form
      //
      //   D * a = 0; a^T * C * a = 1,
      //
      // where we reuse the variable a to represent the entire
      // 6-element parameterization, [a, b, c, d, e, f]^T.  D is the
      // Nx6 "design matrix," [[x_0^2, x_0 * y_0, y_0^2, x_0, y_0, 1],
      // [x_1^2, x_1 * y_1, y_1^2, x_1, y_1, 1], ...].  C is the 6x6
      // constraint matrix, identically zero except for a -1 in the
      // second element of the second row, and 2s in the third element
      // of the first row and the first element of the third row.
      // 
      // Minimizing F(x, y) in the least squares sense, we have
      //
      //   aHat = min_over_a(a^T * S * a)
      //
      // with the constraint that
      //
      //   a^T * C * a = 1,
      //
      // where the "scatter matrix" S = D^T * D.
      //
      // This is readily solved using Lagrange multipliers.  Here, we
      // decompose the solution following Halir & Flasser (see
      // reference [1] in the documentation of this function).

      // This line may be O(N), depending on what type of iterator.
      unsigned int numberOfPoints = endIter - beginIter;

      // Left half of the design matrix.
      brick::numeric::Array2D<Type> D1(numberOfPoints, 3);

      // Right half of the design matrix.
      brick::numeric::Array2D<Type> D2(numberOfPoints, 3);

      // Fill in the design matrix.
      unsigned int rowNumber = 0;
      while(beginIter != endIter) {
        brick::numeric::Vector2D<Type> samplePoint = *beginIter;
        D1(rowNumber, 0) = samplePoint.x() * samplePoint.x();
        D1(rowNumber, 1) = samplePoint.x() * samplePoint.y();
        D1(rowNumber, 2) = samplePoint.y() * samplePoint.y();
        D2(rowNumber, 0) = samplePoint.x();
        D2(rowNumber, 1) = samplePoint.y();
        D2(rowNumber, 2) = Type(1);

        ++beginIter;
        ++rowNumber;
      }

      // Compute the four quadrants of the scatter matrix.  Actually,
      // the upper right and lower left quadrants are identical, so on
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
      brick::numeric::Array2D<Type> TT = (
        Type(-1)
        * brick::numeric::matrixMultiply<Type>(
          brick::linearAlgebra::inverse(S3), S2.transpose()));
      
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
        allParameters, m_origin, m_semimajorAxis, m_semiminorAxis);
    }


    // Convert from implicit ellipse representation to trigonometric
    // parameters.
    template <class Type>
    void
    Ellipse2D<Type>::
    convertAlgebraicToTrigonometric(
      brick::numeric::Array1D<Type> const& algebraicParameters,
      brick::numeric::Vector2D<Type>& origin,
      brick::numeric::Vector2D<Type>& semimajorAxis,
      brick::numeric::Vector2D<Type>& semiminorAxis)
    {
      // This conversion is largely copied from an article written by
      // Robert J. Lopez, and published on the maplesoft tips &
      // techniques page.  We use it without proof here, although
      // R. J. Lopez's article includes a derivation, which follows
      // these steps: substitute u = (x - h) and v = (y - k) into the
      // algebraic distance equation.  When h and k are chosen so that
      // the ellipse is centered at the origin, thn the algebraic
      // distance equation will have no coefficients for the linear
      // terms in u and v.  Setting these coefficients to zeros and
      // solving for h and k gives the origin.
      //
      // Next, substituding for h and k to get a quadratic with no
      // linear terms, we have the equation for a non-axis-aligned
      // ellipse at the origin.  The constant element of this
      // polynomial happens to look like a determinant over a simple
      // quadratic, but that's not really important.  This explains
      // the ifdef'd out code, below.
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
      if(algebraicParameters.size() != 6) {
        BRICK_THROW(brick::common::ValueException,
                    "Ellipse2D::convertAlgebraicToTrigonometric()",
                    "Parameter vector must have six elements.");
      }

      // First extract parameterization.
      Type aa = algebraicParameters[0];
      Type bb = algebraicParameters[1];
      Type cc = algebraicParameters[2];
      Type dd = algebraicParameters[3];
      Type ee = algebraicParameters[4];
      Type ff = algebraicParameters[5];

      // Compute center of ellipse.
#if 0
      brick::numeric::Array2D<Type> AA(3, 3);
      AA(0,0) = 2 * aa;
      AA(0,1) = bb;
      AA(0,2) = dd;
      AA(1,0) = bb;
      AA(1,1) = 2 * cc;
      AA(1,2) = ee;
      AA(2,0) = dd;
      AA(2,1) = ee;
      AA(2,2) = 2 * ff;
      Type gamma = brick::linearAlgebra::determinant(AA);
#else
      Type gamma = (-2 * aa * ee*ee
                    + 8 * aa * cc * ff
                    + 2 * bb * ee * dd
                    - 2 * ff * bb * bb
                    - 2 * cc * dd * dd);
#endif
      Type phi = bb * bb - 4 * aa * cc;
      Type rho = gamma / (2 * phi);

      origin.setValue((2 * cc * dd - bb * ee) / phi,
                      (2 * aa * ee - bb * dd) / phi);


      // Compute axes.
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
      Type PP = aa * CC * CC + cc * SS * SS + bb * CC * SS;
      Type QQ = aa * SS * SS + cc * CC * CC - bb * CC * SS;
      
      Type alpha = brick::common::squareRoot(rho / PP);
      Type beta = brick::common::squareRoot(rho / QQ);
      
      if(alpha >= beta) {
        semimajorAxis.setValue(alpha * CC, alpha * SS);
        semiminorAxis.setValue(-beta * SS, beta * CC);
      } else {
        semimajorAxis.setValue(-beta * SS, beta * CC);
        semiminorAxis.setValue(alpha * CC, alpha * SS);
      }
    }


    /* ======= Non-member functions. ======= */

    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream, const Ellipse2D<Type>& ellipse)
    {
      stream << "Ellipse2D{ "
             << ellipse.getOrigin() << ", "
             << ellipse.getSemimajorAxis() << ", "
             << ellipse.getSemiminorAxis() << " }";
      return stream;
    }
    
  } // namespace geometry
    
} // namespace brick

#endif /* #ifndef BRICK_GEOMETRY_CIRCLE2D_IMPL_HH */
