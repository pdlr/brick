/**
***************************************************************************
* @file brick/computerVision/threePointAlgorithm_impl.hh
*
* Header file defining inline and template functions from
* threePointAlgorithm.hh
*
* Copyright (C) 2009,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_THREEPOINTALGORITHM_IMPL_HH
#define BRICK_COMPUTERVISION_THREEPOINTALGORITHM_IMPL_HH

// This file is included by threePointAlgorithm.hh, and should not be
// directly included by user code, so no need to include
// threePointAlgorithm.hh here.
//
// #include <brick/computerVision/threePointAlgorithm.hh>

#include <cmath>
#include <limits>
#include <brick/common/complexNumber.hh>
#include <brick/common/mathFunctions.hh>
#include <brick/computerVision/registerPoints3D.hh>
#include <brick/numeric/solveQuartic.hh>
#include <brick/numeric/utilities.hh>

namespace brick {

  namespace computerVision {

    // This function implements the "three point perspective pose
    // estimation algorithm" of Grunert[1][2] for recovering the
    // camera-frame coordinates of the corners of a triangle of known
    // size, given the projections of those corners in the camera
    // image.
    template <class FloatType, class IterType>
    unsigned int
    threePointAlgorithm(brick::numeric::Vector3D<FloatType> const& w0,
                        brick::numeric::Vector3D<FloatType> const& w1,
                        brick::numeric::Vector3D<FloatType> const& w2,
                        brick::numeric::Vector2D<FloatType> const& u0,
                        brick::numeric::Vector2D<FloatType> const& u1,
                        brick::numeric::Vector2D<FloatType> const& u2,
                        CameraIntrinsicsPinhole<FloatType> const& intrinsics,
                        IterType p0OutputIter,
                        IterType p1OutputIter,
                        IterType p2OutputIter,
                        FloatType epsilon)
    {
      unsigned int numberOfSolutions = 0;

      // Following the conventions of the paper, we define j0, j1, j2
      // to be unit vectors in the camera coordinate frame that point
      // toward the three world points w0, w1, and w2, respectively.
      brick::numeric::Vector3D<FloatType> j0 =
        intrinsics.reverseProject(u0).getDirectionVector();
      brick::numeric::Vector3D<FloatType> j1 =
        intrinsics.reverseProject(u1).getDirectionVector();
      brick::numeric::Vector3D<FloatType> j2 =
        intrinsics.reverseProject(u2).getDirectionVector();

      // Define alpha to be the angle between j1 and j2, beta to be
      // the angle between j0 and j2, and gamma to be the angle
      // between j0 and j1.
      FloatType cosAlpha = brick::numeric::dot<FloatType>(j1, j2);
      FloatType cosBeta = brick::numeric::dot<FloatType>(j0, j2);
      FloatType cosGamma = brick::numeric::dot<FloatType>(j0, j1);

      // Similarly, define a, b, and c to be the distance between the
      // 3D points that define alpha, beta, and gamma, respectively.
      // We only need the squares of these distances, so we save a
      // sqrt() call for each by computing only the square.
      FloatType a2 = brick::numeric::magnitudeSquared<FloatType>(w1 - w2);
      FloatType b2 = brick::numeric::magnitudeSquared<FloatType>(w0 - w2);
      FloatType c2 = brick::numeric::magnitudeSquared<FloatType>(w0 - w1);

      // If we define s0, s1, and s2 to be the distances from the
      // camera focus to each of the three points, then we have:
      //
      //   p0 = s0 * j0
      //   p1 = s1 * j1
      //   p2 = s2 * j2
      //
      // Where p0, p1, and p2 are the positions of the three 3D points
      // in camera coordinates.  The law of cosines gives us:
      //
      //   s1^2 + s2^2 - 2*s1*s2*cos(alpha) = a^2
      //   s0^2 + s2^2 - 2*s0*s2*cos(beta) = b^2
      //   s0^2 + s1^2 - 2*s0*s1*cos(gamma) = c^2

      // If we choose k1, k2 so that s1 = k1*s0, and s2 = k2*s0 and
      // substitute into (and rearrange) each of the law of cosines
      // equations, we get:
      //
      //   s0^2 = a^2 / (k1^2 + k2^2 + 2*k1*k2*cos(alpha))
      //   s0^2 = b^2 / (1 + k2^2 + 2*k2*cos(beta))
      //   s0^2 = c^2 / (1 + k1^2 + 2*k1*cos(gamma))
      //
      // Combining these to eliminate k2, we obtain an expression for
      // k1 in terms of k2, and a quartic equation in k2.  These are
      // equations 8 and 9 in [1], and are not reproduced in this
      // comment (although the code below implements first the quartic
      // equation, and then the expression for k1).
      FloatType s0Array0[4];
      FloatType s1Array0[4];
      FloatType s2Array0[4];
      FloatType condition0;
      unsigned int newSolutions0 = solveThreePointAlgorithmQuarticSystem(
        cosAlpha, cosBeta, cosGamma, a2, b2, c2, epsilon,
        &(s0Array0[0]), &(s1Array0[0]), &(s2Array0[0]),
        condition0);
      FloatType s0Array1[4];
      FloatType s1Array1[4];
      FloatType s2Array1[4];
      FloatType condition1;
      unsigned int newSolutions1 = solveThreePointAlgorithmQuarticSystem(
        cosBeta, cosAlpha, cosGamma, b2, a2, c2, epsilon,
        &(s1Array1[0]), &(s0Array1[0]), &(s2Array1[0]),
        condition1);
      FloatType s0Array2[4];
      FloatType s1Array2[4];
      FloatType s2Array2[4];
      FloatType condition2;
      unsigned int newSolutions2 = solveThreePointAlgorithmQuarticSystem(
        cosBeta, cosGamma, cosAlpha, b2, c2, a2, epsilon,
        &(s1Array2[0]), &(s2Array2[0]), &(s0Array2[0]),
        condition2);

      FloatType* s0Array;
      FloatType* s1Array;
      FloatType* s2Array;
      unsigned int newSolutions;
      if((condition0 >= condition1) && (condition0 >= condition2) ) {
        s0Array = s0Array0;
        s1Array = s1Array0;
        s2Array = s2Array0;
        newSolutions = newSolutions0;
      } else if((condition1 >= condition0) && (condition1 >= condition2) ) {
        s0Array = s0Array1;
        s1Array = s1Array1;
        s2Array = s2Array1;
        newSolutions = newSolutions1;
      } else {
        s0Array = s0Array2;
        s1Array = s1Array2;
        s2Array = s2Array2;
        newSolutions = newSolutions2;
      }

      for(unsigned int ii = 0; ii < newSolutions; ++ii) {
        *p0OutputIter = s0Array[ii] * j0;
        *p1OutputIter = s1Array[ii] * j1;
        *p2OutputIter = s2Array[ii] * j2;

        ++p0OutputIter;
        ++p1OutputIter;
        ++p2OutputIter;
        ++numberOfSolutions;
      }
      return numberOfSolutions;
    }


    template<class FloatType, class InIter3D, class InIter2D>
    brick::numeric::Transform3D<FloatType>
    threePointAlgorithmRobust(
      InIter3D worldPointsBegin,
      InIter3D worldPointsEnd,
      InIter2D imagePointsBegin,
      CameraIntrinsicsPinhole<FloatType> const& intrinsics,
      size_t iterations,
      FloatType inlierProportion,
      FloatType& score,
      brick::random::PseudoRandom& pRandom)
    {
      // State variables so we'll remember the correct essential
      // matrix once we find it.
      FloatType bestErrorSoFar = std::numeric_limits<FloatType>::max();
      brick::numeric::Transform3D<FloatType> selectedCandidate;

      // Sanity check arguments.
      size_t numberOfPoints = worldPointsEnd - worldPointsBegin;
      if(numberOfPoints < 3) {
        BRICK_THROW(brick::common::ValueException,
                    "threePointAlgorithmRobust()",
                    "Input sequence must have at least three elements.");
      }

      // Copy input points into local buffers.
      std::vector< brick::numeric::Vector3D<FloatType> > worldPoints(numberOfPoints);
      std::vector< brick::numeric::Vector2D<FloatType> > imagePoints(numberOfPoints);
      std::copy(worldPointsBegin, worldPointsEnd, worldPoints.begin());
      std::copy(imagePointsBegin, imagePointsBegin + numberOfPoints,
                imagePoints.begin());

      // Make a buffer to hold points in camera space (and from which
      // to compute residual errors.
      std::vector< brick::numeric::Vector3D<FloatType> > cameraPoints(numberOfPoints);

      // Start the algorithm!
      for(size_t ii = 0; ii < iterations; ++ii) {

        // Select three points.
        for(size_t jj = 0; jj < 3; ++jj) {
          int selectedIndex = pRandom.uniformInt(jj, numberOfPoints);
          if(selectedIndex != static_cast<int>(jj)) {
            std::swap(worldPoints[jj], worldPoints[selectedIndex]);
            std::swap(imagePoints[jj], imagePoints[selectedIndex]);
          }
        }

        // Get candidate cameraPoints.
        brick::numeric::Vector3D<FloatType> testPoints0_cam[4];
        brick::numeric::Vector3D<FloatType> testPoints1_cam[4];
        brick::numeric::Vector3D<FloatType> testPoints2_cam[4];
        unsigned int numberOfSolutions = threePointAlgorithm(
          worldPoints[0], worldPoints[1], worldPoints[2],
          imagePoints[0], imagePoints[1], imagePoints[2], intrinsics,
          testPoints0_cam, testPoints1_cam, testPoints2_cam);

        // Test each candidate solution.
        for(size_t jj = 0; jj < numberOfSolutions; ++jj) {

          // Recover the camTworld transform corresponding to this
          // solution.
          cameraPoints[0] = testPoints0_cam[jj];
          cameraPoints[1] = testPoints1_cam[jj];
          cameraPoints[2] = testPoints2_cam[jj];
          brick::numeric::Transform3D<FloatType> camTworld =
            registerPoints3D<FloatType>(
            worldPoints.begin(), worldPoints.begin() + 3, cameraPoints.begin());

          // Transform all world points into camera coordinates.
          std::transform(worldPoints.begin(), worldPoints.end(),
                         cameraPoints.begin(), camTworld.getFunctor());

          // Project all camera points into image coordinates and
          // compute residuals..
          std::vector<FloatType> residualVector(numberOfPoints);
          for(size_t kk = 0; kk < cameraPoints.size(); ++kk) {
            brick::numeric::Vector2D<FloatType> testPoint_image =
              intrinsics.project(cameraPoints[kk]);
            residualVector[kk] = brick::numeric::magnitudeSquared<FloatType>(
              testPoint_image - imagePoints[kk]);
          }

          // Compute robust error statistic.
          //
          // Note(xxx): Better not to sort here, since it changes the
          // algorithm to O(NlogN).
          std::sort(residualVector.begin(), residualVector.end());
          int testIndex = static_cast<int>(
            inlierProportion * (residualVector.size() - 1) + 0.5);
          if(testIndex >= static_cast<int>(residualVector.size())) {
            testIndex = residualVector.size() - 1;
          }
          if(testIndex < 0) {
            testIndex = 0;
          }
          FloatType errorValue = residualVector[testIndex];

          // Remember candidate if it's the best so far.
          if(errorValue < bestErrorSoFar) {
            selectedCandidate = camTworld;
            bestErrorSoFar = errorValue;
          }
        }
      }
      score = bestErrorSoFar;
      return selectedCandidate;
    }


    template <class FloatType, class OutIter>
    unsigned int
    solveThreePointAlgorithmQuarticSystem(
      FloatType cosAlpha, FloatType cosBeta, FloatType cosGamma,
      FloatType a2, FloatType b2, FloatType c2,
      FloatType epsilon,
      OutIter s0Iter, OutIter s1Iter, OutIter s2Iter,
      FloatType& condition)
    {
      unsigned int numberOfSolutions = 0;
      condition = std::numeric_limits<FloatType>::max();

      FloatType cos2Alpha = cosAlpha * cosAlpha;
      FloatType cos2Beta = cosBeta * cosBeta;
      FloatType cos2Gamma = cosGamma * cosGamma;
      FloatType a2OverB2 = a2 / b2;
      FloatType a2MinusC2OverB2 = (a2 - c2) / b2;
      FloatType a2MinusC2OverB2Sq = a2MinusC2OverB2 * a2MinusC2OverB2;
      FloatType a2MinusC2OverB2Minus1 = a2MinusC2OverB2 - 1.0;
      FloatType a2PlusC2OverB2 = (a2 + c2) / b2;
      FloatType b2MinusA2OverB2 = (b2 - a2) / b2;
      FloatType b2MinusC2OverB2 = (b2 - c2) / b2;
      FloatType c2OverB2 = c2 / b2;
      FloatType oneMinusA2PlusC2OverB2 = (1.0 - a2PlusC2OverB2);
      FloatType oneMinusA2MinusC2OverB2 = (1.0 - a2MinusC2OverB2);
      FloatType onePlusA2MinusC2OverB2 = (1.0 + a2MinusC2OverB2);

      FloatType A0 = (onePlusA2MinusC2OverB2 * onePlusA2MinusC2OverB2
                   - 4.0 * a2OverB2 * cos2Gamma);

      FloatType A1 = 4.0 * ((-a2MinusC2OverB2 * onePlusA2MinusC2OverB2 * cosBeta)
                         + (2.0 * a2OverB2 * cos2Gamma * cosBeta)
                         - (oneMinusA2PlusC2OverB2 * cosAlpha * cosGamma));

      FloatType A2 = 2.0 * (a2MinusC2OverB2Sq - 1.0
                         + 2.0 * a2MinusC2OverB2Sq * cos2Beta
                         + 2.0 * b2MinusC2OverB2 * cos2Alpha
                         - 4.0 * a2PlusC2OverB2 * cosAlpha * cosBeta * cosGamma
                         + 2.0 * b2MinusA2OverB2 * cos2Gamma);

      FloatType A3 = 4.0 * (a2MinusC2OverB2 * oneMinusA2MinusC2OverB2 * cosBeta
                         - oneMinusA2PlusC2OverB2 * cosAlpha * cosGamma
                         + 2.0 * c2OverB2 * cos2Alpha * cosBeta);

      FloatType A4 = (a2MinusC2OverB2Minus1 * a2MinusC2OverB2Minus1
                   - 4.0 * c2OverB2 * cos2Alpha);

      // Now we solve for the roots of the quartic, which tell us
      // valid values of k2.  Each root corresponds to a scale factor
      // that's consistent with the observed data.
      brick::common::ComplexNumber<FloatType> k2Roots[4];
      brick::numeric::solveQuartic(
        A3 / A4, A2 / A4, A1 / A4, A0 / A4,
        k2Roots[0], k2Roots[1], k2Roots[2], k2Roots[3]);

      // For real value of k1, there's a corresponding value of k2 (see
      // Eq. 8 in [1]).
      for(unsigned int ii = 0; ii < 4; ++ii) {
        bool isReal = (brick::common::absoluteValue(
                         k2Roots[ii].getImaginaryPart())
                       <= epsilon);
        if(isReal) {
          FloatType k2 = k2Roots[ii].getRealPart();
          FloatType numerator = (a2MinusC2OverB2Minus1 * k2 * k2
                              - 2.0 * a2MinusC2OverB2 * cosBeta * k2
                              + onePlusA2MinusC2OverB2);
          FloatType denominator = 2.0 * (cosGamma - k2 * cosAlpha);
          FloatType magDenominator = brick::common::absoluteValue(denominator);
          if(magDenominator < condition) {
            condition = magDenominator;
          }
          if(magDenominator < epsilon) {
            continue;
          }
          FloatType k1 = numerator / denominator;

          // Now that we have k1 and k2, recover the distance from the
          // focus to each of the observed points.
          FloatType s0Sq = c2 / (1.0 + k1 * k1 - 2.0 * k1 * cosGamma);
          if(s0Sq < 0.0) {
            continue;
          }
          *s0Iter = brick::common::squareRoot(s0Sq);
          *s1Iter = k1 * (*s0Iter);
          *s2Iter = k2 * (*s0Iter);

          ++s0Iter;
          ++s1Iter;
          ++s2Iter;
          ++numberOfSolutions;
        }
      }
      if(numberOfSolutions == 0) {
        condition = 0.0;
      }
      return numberOfSolutions;
    }

  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_THREEPOINTALGORITHM_IMPL_HH */
