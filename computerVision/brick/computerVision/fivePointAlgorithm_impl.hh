/**
***************************************************************************
* @file brick/computerVision/fivePointAlgorithm_impl.hh
*
* Header file defining inline and function templates declared in
* fivePointAlgorithm.hh.
*
* Copyright (C) 2009,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_FIVEPOINTALGORITHM_IMPL_HH
#define BRICK_COMPUTERVISION_FIVEPOINTALGORITHM_IMPL_HH

// This file is included by fivePointAlgorithm.hh, and should not be
// directly included by user code, so no need to include
// fivePointAlgorithm.hh here.
//
// #include <brick/computerVision/fivePointAlgorithm.hh>

#include <cmath>
#include <complex>
#include <limits>
#include <brick/computerVision/cameraIntrinsicsPinhole.hh>
#include <brick/computerVision/threePointAlgorithm.hh>
#include <brick/geometry/ray2D.hh>
#include <brick/geometry/utilities2D.hh>
#include <brick/geometry/utilities3D.hh>
#include <brick/linearAlgebra/linearAlgebra.hh>
#include <brick/numeric/subArray2D.hh>
#include <brick/numeric/utilities.hh>

namespace brick {

  namespace computerVision {


    template<class FloatType, class Iterator>
    std::vector< brick::numeric::Array2D<FloatType> >
    fivePointAlgorithm(Iterator sequence0Begin, Iterator sequence0End,
                       Iterator sequence1Begin)
    {
      // The five (or more) pairs of input points give at least five
      // linear constraints on the essential matrix, E.  Since E has
      // nine elements, this means there's a four-dimentional null
      // space of these constraints.  We solve for an orthogonal basis
      // of this null space here.

      // For each pair of points q, q', the essential matrix, E,
      // satisfies the equation:
      //
      //   transpose(q') * E * q = 0
      //
      // where the points q are drawn from sequence0, and the points
      // q' are drawn from sequence1.
      //
      // We rearrange this equation to get
      //
      //   ||e_00*q_x*q'_x + e_01*q_y*q'_x + e_02*q'_x
      //     + e_10*q_x*q'_y + e_11*q_y*q'_y + e_12*q'_y
      //     + e_20*q_x + e_21*q_y + e_22|| = 0
      //
      // where e_ij is the element of E in the i^th row and the j^th
      // column, q_x & q_y are the x & y components of q,
      // respectively, and q'_x & q'_y are the x & y components of q',
      // respectively.
      //
      // or,
      //
      //   ||A * vec(E)|| = 0
      //
      // With the matrix A as specified in the code below.
      size_t numberOfCorrespondences = sequence0End - sequence0Begin;
      brick::numeric::Array2D<FloatType> AMatrix(numberOfCorrespondences, 9);
      for(size_t rowIndex = 0; rowIndex < numberOfCorrespondences; ++rowIndex) {
        brick::numeric::Array1D<FloatType> currentRow = AMatrix.getRow(rowIndex);
        const brick::numeric::Vector2D<FloatType>& qq = *sequence0Begin;
        const brick::numeric::Vector2D<FloatType>& qPrime = *sequence1Begin;
        currentRow[0] = qq.x() * qPrime.x();
        currentRow[1] = qq.y() * qPrime.x();
        currentRow[2] = qPrime.x();
        currentRow[3] = qq.x() * qPrime.y();
        currentRow[4] = qq.y() * qPrime.y();
        currentRow[5] = qPrime.y();
        currentRow[6] = qq.x();
        currentRow[7] = qq.y();
        currentRow[8] = 1.0;
        ++sequence0Begin;
        ++sequence1Begin;
      }

      // Following [1], we solve for the four dimensional null space
      // using SVD.
      brick::numeric::Array2D<FloatType> uArray;
      brick::numeric::Array1D<FloatType> sigmaArray;
      brick::numeric::Array2D<FloatType> vTransposeArray;
      brick::linearAlgebra::singularValueDecomposition(
        AMatrix, uArray, sigmaArray, vTransposeArray, true);
      brick::numeric::Array2D<FloatType> E0Array(3, 3);
      brick::numeric::Array2D<FloatType> E1Array(3, 3);
      brick::numeric::Array2D<FloatType> E2Array(3, 3);
      brick::numeric::Array2D<FloatType> E3Array(3, 3);
      std::copy(vTransposeArray.getRow(5).begin(),
                vTransposeArray.getRow(5).end(),
                E0Array.begin());
      std::copy(vTransposeArray.getRow(6).begin(),
                vTransposeArray.getRow(6).end(),
                E1Array.begin());
      std::copy(vTransposeArray.getRow(7).begin(),
                vTransposeArray.getRow(7).end(),
                E2Array.begin());
      std::copy(vTransposeArray.getRow(8).begin(),
                vTransposeArray.getRow(8).end(),
                E3Array.begin());

      // Let E = x*E0 + y*E1 + z*E2 + w*E3, and assume w == 1 (Yuck.
      // We will fix this soon).  We have (equations (2) and (3) of [1])
      //
      //   det(E) = 0
      //
      // and
      //
      //   2 * E * transpose(E) * E - trace(E * transpose(E)) * E = 0
      //
      // These two equations give us ten cubic constraints on x, y,
      // and z.  As described in the paper, we write these constraints
      // as the matrix product of a 10x20 coefficient matrix, M, and a
      // monomials vector X.
      //
      //   M * X = [0]
      //
      //   X = [x^3, x^2*y, x^2*z, x*y^2, x*y*z, x*z^2, y^3, y^2*z, y*z^2, z^3,
      //        x^2, x*y, x*z, y^2, y*z, z^2, x, y, z, 1]
      //
      // Although computing M just involves doing a bunch of matrix
      // multiplications, tranposes, determinants, etc., the
      // computation is quite tedious.  In fact, even typing in the
      // result is quite tedious.  For this reason, we use
      // automatically generated C code exported from the computer
      // algebra system Maxima.
      brick::numeric::Array2D<FloatType> M = generateFivePointConstraintMatrix(
        E0Array, E1Array, E2Array, E3Array);

#if 0
      brick::numeric::Array2D<FloatType> diff = M - MnL;
      brick::numeric::Array2D<bool> flags(diff.rows(), diff.columns());
      for(size_t nn = 0; nn < diff.size(); ++nn) {
        flags[nn] = std::fabs(diff[nn]) < 1.0E-9;
      }

      std::cout << "===================" << std::endl;
      FloatType target = MnL(0, 0);
      for(size_t mm = 0; mm < M.rows(); ++mm) {
        for(size_t nn = 0; nn < M.columns(); ++nn) {
          if(std::fabs(M(mm, nn) - target) < 1.0E-9) {
            std::cout << " [" << mm << ", " << nn << "] matches." << std::endl;
          }
        }
      }
      std::cout << "===================" << std::endl;
#endif

      // The paper calls for Gauss-Jordan elimination to reduce the
      // first 10 columns of M to the identity matrix, so that the
      // rows of the eliminated matrix form a Groebner basis of the
      // ideal of the 10 original cubic constraints (see [3] for more
      // information on Groebner bases, ideals, and algebraic geometry
      // in general).  We don't have a Gauss-Jordan elimination
      // routine coded up, so we use what we have available.

      // Extract the first 10 columns of M.
      brick::numeric::Array2D<FloatType> M0 = brick::numeric::subArray(
        M, brick::numeric::Slice(), brick::numeric::Slice(0, 10));

      // Invert them (presumably using Gauss-Jordan elimination.
      brick::numeric::Array2D<FloatType> M0Inv =
        brick::linearAlgebra::inverse(M0);

      // Extract the 10th through 20th columns of M.
      brick::numeric::Array2D<FloatType> M1 = brick::numeric::subArray(
        M, brick::numeric::Slice(), brick::numeric::Slice(10, 0));

      // And manipulate them so that they look like they went through
      // the same elimination process.  Following the paper, we call
      // the result of this manipulation B, although the real Brobner
      // basis is the set of polynomials whose coefficients are drawn
      // from the 10x20 matrix resulting from the elimination (first
      // 10 columns are 10x10 identity matrix, remaining columns equal
      // to our 10x10 matrix B.
      brick::numeric::Array2D<FloatType> B =
        brick::numeric::matrixMultiply<FloatType>(M0Inv, M1);

      // Since we eliminated the left half of M, we're safe to assume
      // that the leading terms of the Groebner basis we just computed
      // do not include the monomials x, y, or z.  This means that x,
      // y, and z will be in the quotient ring of complex third-order
      // polynomials over the ideal of the polynomials expressed in M.
      //
      // Theorem xxx in [3] gives us that, at the solutions, x, of the
      // polynomials in M,
      //
      //   f(x) * transpose(u(x)) = transpose(u(x)) * A_f
      //
      // Where f(x) is an arbitrary polynomial in the quotient ring,
      // u(x) is the vector of basis monomials of the quotient ring,
      // and A_f is the action matrix associated with f(x).
      //
      // This is useful because it means that, at the solutions to our
      // system of third order polynomials, the values of the
      // monomials in u(x) correspond to the left eigenvectors of A_f.
      // But we've already agreed that x, y, and z are three of the
      // monomials in u(x), so if we can just find A_f for an
      // arbitrary polynomial in the quotient ring, then we can take
      // its left eigenvectors and grab the elements corresponding to
      // x, y, and z to find our solutions!
      //
      // We'll use the action matrix for the polynomial x, and
      // transpose it so that we can use our normal right-eigenvector
      // solver to find the eigenvectors.

      // Warning(xxx): We'll use the action matrix from Stewenius &
      // Nister's paper, which doesn't appear to correspond to
      // suggested degree-then-lexicographic order of monomials.  The
      // right solution to this is to generate the correct action
      // matrix below, but instead we temporarily shuffle the columns
      // of our constraint matrix to match the action matrix.  This
      // shuffling happens under the hood in the call to
      // generateFivePointConstraintMatrix(), above.
      brick::numeric::Array2D<FloatType> At(10, 10);
      At = 0.0;
      for(size_t ii = 0; ii < 10; ++ii) {
        At(0, ii) = -B(0, ii);
        At(1, ii) = -B(1, ii);
        At(2, ii) = -B(2, ii);
        At(3, ii) = -B(4, ii);
        At(4, ii) = -B(5, ii);
        At(5, ii) = -B(7, ii);
      }
      At(6, 0) = 1.0;
      At(7, 1) = 1.0;
      At(8, 3) = 1.0;
      At(9, 6) = 1.0;

      // Solutions for x, y, z can be found from the 10 eigenvectors
      // of this (tranposed) action matrix.  We normalize the final
      // element of each eigenvector to 1.0.
      brick::numeric::Array1D< std::complex<FloatType> > eigenvalues;
      brick::numeric::Array2D< std::complex<FloatType> > eigenvectors;
      brick::linearAlgebra::eigenvectors(At, eigenvalues, eigenvectors);

      // Now that we have solutions for x, y, and z, we can go back
      // and use them to reconstruct potential solutions for the
      // essential matrix.
      std::vector< brick::numeric::Array2D<FloatType> > result;
      for(size_t ii = 0; ii < 10; ++ii) {
        if(eigenvalues[ii].imag() == 0) {
          FloatType xx = eigenvectors(6, ii).real() / eigenvectors(9, ii).real();
          FloatType yy = eigenvectors(7, ii).real() / eigenvectors(9, ii).real();
          FloatType zz = eigenvectors(8, ii).real() / eigenvectors(9, ii).real();
          Array2D<FloatType> solution =
            (xx * E0Array) + (yy * E1Array) + (zz * E2Array) + E3Array;
          result.push_back(solution);
        }
      }

      return result;
    }


    template<class FloatType, class Iterator>
    brick::numeric::Array2D<FloatType>
    fivePointAlgorithmRobust(Iterator sequence0Begin, Iterator sequence0End,
                             Iterator sequence1Begin,
                             size_t iterations,
                             FloatType inlierProportion,
                             FloatType& score,
                             brick::random::PseudoRandom pRandom)
    {
      // State variables so we'll remember the correct essential
      // matrix once we find it.
      FloatType bestErrorSoFar = std::numeric_limits<FloatType>::max();
      brick::numeric::Array2D<FloatType> selectedCandidate(3, 3);
      selectedCandidate = 0.0;

      // Copy input points into local buffers.
      size_t numberOfPoints = sequence0End - sequence0Begin;
      std::vector< brick::numeric::Vector2D<FloatType> >
        qVector(numberOfPoints);
      std::vector< brick::numeric::Vector2D<FloatType> > qPrimeVector(numberOfPoints);
      std::copy(sequence0Begin, sequence0End, qVector.begin());
      std::copy(sequence1Begin, sequence1Begin + numberOfPoints,
                qPrimeVector.begin());

      for(size_t ii = 0; ii < iterations; ++ii) {
        // Select five points.
        for(size_t jj = 0; jj < 5; ++jj) {
          int selectedIndex = pRandom.uniformInt(jj, numberOfPoints);
          if(selectedIndex != static_cast<int>(jj)) {
            std::swap(qVector[jj], qVector[selectedIndex]);
            std::swap(qPrimeVector[jj], qPrimeVector[selectedIndex]);
          }
        }

        // Get candidate essential matrices.
        std::vector< brick::numeric::Array2D<FloatType> > EVector =
          fivePointAlgorithm<FloatType>(qVector.begin(), qVector.begin() + 5,
                                        qPrimeVector.begin());

        // Test each candidate.
        for(size_t jj = 0; jj < EVector.size(); ++jj) {
          Array2D<FloatType> EE = EVector[jj];

          // Compute residuals for all input points.
          std::vector<FloatType> residualVector(numberOfPoints);
          for(size_t kk = 0; kk < numberOfPoints; ++kk) {
            residualVector[kk] = checkEpipolarConstraint(
              EE, qVector[kk], qPrimeVector[kk]);
          }

          // Compute robust error statistic.
          //
          // Note(xxx): Better not to sort here, since it changes the
          // algorithm to O(NlogN).
          std::sort(residualVector.begin(), residualVector.end());
          int testIndex = static_cast<int>(
            inlierProportion * residualVector.size() + 0.5);
          FloatType errorValue = residualVector[testIndex];

          // Remember candidate if it's the best so far.
          if(errorValue < bestErrorSoFar) {
            selectedCandidate = EE;
            bestErrorSoFar = errorValue;
          }
        }
      }
      score = bestErrorSoFar;
      return selectedCandidate;
    }


    template<class FloatType, class Iterator>
    void
    fivePointAlgorithmRobust(Iterator sequence0Begin, Iterator sequence0End,
                             Iterator sequence1Begin,
                             Iterator sequence2Begin,
                             size_t iterations,
                             FloatType inlierProportion,
                             brick::numeric::Array2D<FloatType>& cam2Ecam0,
                             brick::numeric::Transform3D<FloatType>& cam0Tcam2,
                             brick::numeric::Transform3D<FloatType>& cam1Tcam2,
                             FloatType& score,
                             brick::random::PseudoRandom pRandom)
    {
      // Since we're using calibrated image points, all cameras have
      // the same intrinsics.
      CameraIntrinsicsPinhole<FloatType> intrinsics(
        1, 1, 1.0, 1.0, 1.0, 0.0, 0.0);

      // State variables so we'll remember the correct essential
      // matrices once we find them.
      FloatType bestErrorSoFar = std::numeric_limits<FloatType>::max();
      brick::numeric::Array2D<FloatType> selectedCam2Ecam0(3, 3);
      brick::numeric::Transform3D<FloatType> selectedCam0Tcam2;
      brick::numeric::Transform3D<FloatType> selectedCam1Tcam2;
      selectedCam2Ecam0 = 0.0;

      // Copy input points into local buffers.
      size_t numberOfPoints = sequence0End - sequence0Begin;
      std::vector< brick::numeric::Vector2D<FloatType> > points2D_cam0(numberOfPoints);
      std::vector< brick::numeric::Vector2D<FloatType> > points2D_cam1(numberOfPoints);
      std::vector< brick::numeric::Vector2D<FloatType> > points2D_cam2(numberOfPoints);
      std::copy(sequence0Begin, sequence0End, points2D_cam0.begin());
      std::copy(sequence1Begin, sequence1Begin + numberOfPoints,
                points2D_cam1.begin());
      std::copy(sequence2Begin, sequence2Begin + numberOfPoints,
                points2D_cam2.begin());

      // Allocate storage for temporary values prior to starting loop.
      std::vector< brick::numeric::Vector3D<FloatType> > points3D_cam0(numberOfPoints);
      std::vector< brick::numeric::Vector3D<FloatType> > points3D_cam1(numberOfPoints);
      std::vector< brick::numeric::Vector3D<FloatType> > points3D_cam2(numberOfPoints);
      std::vector<FloatType> residualVector(numberOfPoints);

      // Loop over many sets of random samples.
      for(size_t ii = 0; ii < iterations; ++ii) {

        // Select five points.
        for(size_t jj = 0; jj < 5; ++jj) {
          int selectedIndex = pRandom.uniformInt(jj, numberOfPoints);
          if(selectedIndex != static_cast<int>(jj)) {
            std::swap(points2D_cam0[jj], points2D_cam0[selectedIndex]);
            std::swap(points2D_cam1[jj], points2D_cam1[selectedIndex]);
            std::swap(points2D_cam2[jj], points2D_cam2[selectedIndex]);
          }
        }

        // Get candidate essential matrices.
        std::vector< brick::numeric::Array2D<FloatType> > EVector =
          fivePointAlgorithm<FloatType>(
            points2D_cam0.begin(), points2D_cam0.begin() + 5,
            points2D_cam2.begin());

        // Test each candidate.
        for(size_t jj = 0; jj < EVector.size(); ++jj) {
          Array2D<FloatType> EE = EVector[jj];

          // Recover relative motion between cameras, assuming EE is
          // correct.
          brick::numeric::Transform3D<FloatType> c2Tc0;
          try {
            c2Tc0 = getCameraMotionFromEssentialMatrix(
              EE, points2D_cam0[0], points2D_cam2[0]);
          } catch(brick::common::ValueException&) {
            // Input points were on parallel rays!  No point in evaluating
            // this candidate.
            continue;
          }

          // Given relative motion, recover 3D position of each input
          // point in camera 2 coordinates.
          for(size_t kk = 0; kk < numberOfPoints; ++kk) {
            points3D_cam2[kk] = triangulateCalibratedImagePoint(
              c2Tc0, points2D_cam2[kk], points2D_cam0[kk]);
          }

          // Recover 3D position and orientation of camera 1.
          FloatType internalScore;
          brick::numeric::Transform3D<FloatType> c1Tc2 = threePointAlgorithmRobust(
            points3D_cam2.begin(), points3D_cam2.begin() + 5,
            points2D_cam1.begin(), intrinsics, 1, 1.0, internalScore, pRandom);

          // We expect the model to fit these five points better than
          // it fits other points in the set.  If internalScore
          // doesn't beat the average residual over all points for our
          // best so far, then there's no point in continuing with
          // this candidate.
          if(internalScore >= bestErrorSoFar) {
            continue;
          }

          // Compute the 3D position of each input point in camera 1
          // coordinates.
          for(size_t kk = 0; kk < numberOfPoints; ++kk) {
            points3D_cam1[kk] = c1Tc2 * points3D_cam2[kk];
          }

          // Recover 3D position of each of the input points in camera
          // 0 coordinates.
          brick::numeric::Transform3D<FloatType> c0Tc2 = c2Tc0.invert();
          for(size_t kk = 0; kk < numberOfPoints; ++kk) {
            points3D_cam0[kk] = c0Tc2 * points3D_cam2[kk];
          }

          // Project 3D points into each image, and compute residual.
          for(size_t kk = 0; kk < numberOfPoints; ++kk) {
            brick::numeric::Vector2D<FloatType> residualVec0(
              points3D_cam0[kk].x() / points3D_cam0[kk].z(),
              points3D_cam0[kk].y() / points3D_cam0[kk].z());
            residualVec0 -= points2D_cam0[kk];
            brick::numeric::Vector2D<FloatType> residualVec1(
              points3D_cam1[kk].x() / points3D_cam1[kk].z(),
              points3D_cam1[kk].y() / points3D_cam1[kk].z());
            residualVec1 -= points2D_cam1[kk];
            brick::numeric::Vector2D<FloatType> residualVec2(
              points3D_cam2[kk].x() / points3D_cam2[kk].z(),
              points3D_cam2[kk].y() / points3D_cam2[kk].z());
            residualVec2 -= points2D_cam2[kk];

            residualVector[kk] =
              (brick::numeric::magnitudeSquared<FloatType>(residualVec0)
               + brick::numeric::magnitudeSquared<FloatType>(residualVec1)
               + brick::numeric::magnitudeSquared<FloatType>(residualVec2)) / 3.0;
          }

          // Compute robust error statistic.
          //
          // Note(xxx): Better not to sort here, since it changes the
          // algorithm to O(NlogN).
          std::sort(residualVector.begin(), residualVector.end());
          int testIndex = static_cast<int>(
            inlierProportion * residualVector.size() + 0.5);
          FloatType errorValue = residualVector[testIndex];

          // Remember candidate if it's the best so far.
          if(errorValue < bestErrorSoFar) {
            selectedCam2Ecam0 = EE;
            selectedCam0Tcam2 = c0Tc2;
            selectedCam1Tcam2 = c1Tc2;
            bestErrorSoFar = errorValue;
          }
        }
      }
      score = bestErrorSoFar;
      cam2Ecam0 = selectedCam2Ecam0;
      cam0Tcam2 = selectedCam0Tcam2;
      cam1Tcam2 = selectedCam1Tcam2;
    }


    template <class FloatType>
    FloatType
    checkEpipolarConstraint(
      brick::numeric::Array2D<FloatType> const& fundamentalMx,
      brick::numeric::Vector2D<FloatType>& point0,
      brick::numeric::Vector2D<FloatType>& point1)
    {
      // Compute matrix-vector product F * q.
      FloatType lineCoefficient0 =
        fundamentalMx[0] * point0.x() + fundamentalMx[1] * point0.y()
        + fundamentalMx[2];
      FloatType lineCoefficient1 =
        fundamentalMx[3] * point0.x() + fundamentalMx[4] * point0.y()
        + fundamentalMx[5];
      FloatType lineCoefficient2 =
        fundamentalMx[6] * point0.x() + fundamentalMx[7] * point0.y()
        + fundamentalMx[8];

      // Get the epipolar line corresponding to point0.
      brick::geometry::Ray2D<FloatType> ray(lineCoefficient0, lineCoefficient1, lineCoefficient2);

      // Find the squared distance between q' and this ray.
      brick::numeric::Vector2D<FloatType> closestPoint =
        brick::geometry::findClosestPoint(point1, ray);
      return brick::numeric::magnitudeSquared<FloatType>(closestPoint - point1);
    }


    template <class FloatType>
    brick::numeric::Transform3D<FloatType>
    getCameraMotionFromEssentialMatrix(
      brick::numeric::Array2D<FloatType> const& EE,
      brick::numeric::Vector2D<FloatType> const& testPointCamera0,
      brick::numeric::Vector2D<FloatType> const& testPointCamera1)
    {
      // Note(xxx): Need to put a comment here explaining why this
      // works.  Ref R. Hartley and A. Zisserman, Multiple View
      // Geometry in Computer Vision, Cambridge University Press, ISBN
      // 0-521-62304-9, 2000.  Or ref R. Tsai and T. Huang, Uniqueness
      // and Estimation of Three-Dimensional Motion Parameters of
      // Rigid Objects with Curved Surfaces, IEEE Transactions on
      // Pattern Analysis and Machine Intelligence, 6(1):13-27, 1984.
      brick::numeric::Array2D<FloatType> uArray;
      brick::numeric::Array1D<FloatType> sigmaArray;
      brick::numeric::Array2D<FloatType> vTransposeArray;
      brick::linearAlgebra::singularValueDecomposition(
        EE, uArray, sigmaArray, vTransposeArray);

#if 0
      // Nister's paper outlines the following steps, which are faster
      // than what we're actually using below, and exactly equivalent,
      // but I haven't had time to test and make sure I've implemented
      // them correctly.  Note that the right-handedness check is
      // missing from this code.

      // Multiply by D == [[0, 1, 0], [-1, 0, 0], [0, 0, 1]].
      vTransposeArray.getRow(0) *= -1.0;
      std::swap(vTransposeArray[0], vTransposeArray[3]);
      std::swap(vTransposeArray[1], vTransposeArray[4]);
      std::swap(vTransposeArray[2], vTransposeArray[5]);
      brick::numeric::Array2D<FloatType> rotation0 =
        matrixMultiply<FloatType>(uArray, vTransposeArray);

      brick::numeric::Transform3D<FloatType> c1Tc0(
        rotation0[0], rotation0[1], rotation0[2], uArray(0, 2),
        rotation0[3], rotation0[4], rotation0[5], uArray(1, 2),
        rotation0[6], rotation0[7], rotation0[8], uArray(2, 2),
        0.0,          0.0,          0.0,          1.0);

      // Note: this call might throw...
      brick::numeric::Vector3D<FloatType> testPoint3D0 = triangulateCalibratedImagePoint(
        c1Tc0, testPointCamera1, testPointCamera0);
      brick::numeric::Vector3D<FloatType> testPoint3D1 = c1Tc0 * testPoint3D0;
      if((testPoint3D0.z() < 0.0) && (testPoint3D1.z() < 0.0)) {
        // 3D position is behind both cameras.
        c1Tc0.setTransform(
          -rotation0[0], -rotation0[1], -rotation0[2], -uArray(0, 2),
          -rotation0[3], -rotation0[4], -rotation0[5], -uArray(1, 2),
          -rotation0[6], -rotation0[7], -rotation0[8], -uArray(2, 2),
          0.0,          0.0,          0.0,          1.0);
      } else if((testPoint3D0.z() * testPoint3D1.z()) < 0.0) {
        // 3D position is behind one camera, but not both.
        Transform3D<FloatType> Ht(-1.0, 0.0, 0.0, 0.0,
                       0.0, -1.0, 0.0, 0.0,
                       0.0, 0.0, -1.0, 0.0,
                       2*vTransposeArray(2, 0), 2*vTransposeArray(2, 1),
                       2*vTransposeArray(2, 1), 1.0);
        c1Tc0 = Ht * c1Tc0;
        testPoint3D0 = Ht * testPoint3D0;
        if((testPoint3D0.z() < 0.0) || (testPoint3D1.z() < 0.0)) {
          c1Tc0.setTransform(
            -c1Tc0.value<0, 0>(), -c1Tc0.value<0, 1>(),
            -c1Tc0.value<0, 2>(), -c1Tc0.value<0, 3>(),
            -c1Tc0.value<1, 0>(), -c1Tc0.value<1, 1>(),
            -c1Tc0.value<1, 2>(), -c1Tc0.value<1, 3>(),
            -c1Tc0.value<2, 0>(), -c1Tc0.value<2, 1>(),
            -c1Tc0.value<2, 2>(), -c1Tc0.value<2, 3>(),
            0.0, 0.0, 0.0, 1.0);
        }
      }
      return c1Tc0;
#else
      // Here's a bonehead version of the above.

      // Multiply by D == [[0, 1, 0], [-1, 0, 0], [0, 0, 1]].
      vTransposeArray.getRow(0) *= -1.0;
      std::swap(vTransposeArray[0], vTransposeArray[3]);
      std::swap(vTransposeArray[1], vTransposeArray[4]);
      std::swap(vTransposeArray[2], vTransposeArray[5]);
      brick::numeric::Array2D<FloatType> rotation0 =
        matrixMultiply<FloatType>(uArray, vTransposeArray);

      // Undo previous contortion, and multiply by D == [[0, -1, 0],
      // [1, 0, 0], [0, 0, 1]].
      vTransposeArray.getRow(0) *= -1.0;
      vTransposeArray.getRow(1) *= -1.0;
      brick::numeric::Array2D<FloatType> rotation1 =
        matrixMultiply<FloatType>(uArray, vTransposeArray);

      // SVD could easily return -1 * the rotatin we want.  Of course, this
      // would make the resulting coordinate system left handed.  Check for
      // this here.
      brick::numeric::Vector3D<FloatType> e0(rotation0(0, 0), rotation0(1, 0), rotation0(2, 0));
      brick::numeric::Vector3D<FloatType> e1(rotation0(0, 1), rotation0(1, 1), rotation0(2, 1));
      brick::numeric::Vector3D<FloatType> e2(rotation0(0, 2), rotation0(1, 2), rotation0(2, 2));
      bool isRightHanded = brick::numeric::dot<FloatType>(brick::numeric::cross(e0, e1), e2) > 0.0;
      if(!isRightHanded) {
        rotation0 *= -1.0;
        rotation1 *= -1.0;
      }

      // Now we have some candidates for rotation.  We'll test them
      // using the test point that was passed in.
      brick::numeric::Transform3D<FloatType> c1Tc0_0(
        rotation0[0], rotation0[1], rotation0[2], uArray(0, 2),
        rotation0[3], rotation0[4], rotation0[5], uArray(1, 2),
        rotation0[6], rotation0[7], rotation0[8], uArray(2, 2),
        0.0,          0.0,          0.0,          1.0);
      brick::numeric::Transform3D<FloatType> c1Tc0_1(
        rotation0[0], rotation0[1], rotation0[2], -uArray(0, 2),
        rotation0[3], rotation0[4], rotation0[5], -uArray(1, 2),
        rotation0[6], rotation0[7], rotation0[8], -uArray(2, 2),
        0.0,          0.0,          0.0,          1.0);

      brick::numeric::Vector3D<FloatType> testPoint3D1 = triangulateCalibratedImagePoint(
        c1Tc0_0, testPointCamera1, testPointCamera0);
      brick::numeric::Vector3D<FloatType> testPoint3D0 = c1Tc0_0.invert() * testPoint3D1;
      if((testPoint3D0.z() > 0.0) && (testPoint3D1.z() > 0.0)) {
        return c1Tc0_0;
      }
      if((testPoint3D0.z() < 0.0) && (testPoint3D1.z() < 0.0)) {
        return c1Tc0_1;
      }

      brick::numeric::Transform3D<FloatType> c1Tc0_2(
        rotation1[0], rotation1[1], rotation1[2], uArray(0, 2),
        rotation1[3], rotation1[4], rotation1[5], uArray(1, 2),
        rotation1[6], rotation1[7], rotation1[8], uArray(2, 2),
        0.0,          0.0,          0.0,          1.0);
      brick::numeric::Transform3D<FloatType> c1Tc0_3(
        rotation1[0], rotation1[1], rotation1[2], -uArray(0, 2),
        rotation1[3], rotation1[4], rotation1[5], -uArray(1, 2),
        rotation1[6], rotation1[7], rotation1[8], -uArray(2, 2),
        0.0,          0.0,          0.0,          1.0);
      // Note: this call might throw...
      testPoint3D1 = triangulateCalibratedImagePoint(
        c1Tc0_2, testPointCamera1, testPointCamera0);
      testPoint3D0 = c1Tc0_2.invert() * testPoint3D1;
      if((testPoint3D0.z() > 0.0) && (testPoint3D1.z() > 0.0)) {
        return c1Tc0_2;
      }
      return c1Tc0_3;
#endif
    }


    template <class FloatType>
    brick::numeric::Vector3D<FloatType>
    triangulateCalibratedImagePoint(
      brick::numeric::Transform3D<FloatType> const& c0Tc1,
      brick::numeric::Vector2D<FloatType> const& testPointCamera0,
      brick::numeric::Vector2D<FloatType> const& testPointCamera1)
    {
      // This triangulation is much less efficient that the approach
      // described in [2], but is conceptually clearer, and has nicer
      // error characteristics.

      // Each test point defines a ray in the coordinate system of the
      // corresponding camera.
      brick::geometry::Ray3D<FloatType> ray0_0(
        brick::numeric::Vector3D<FloatType>(0.0, 0.0, 0.0),
        brick::numeric::Vector3D<FloatType>(testPointCamera0.x(), testPointCamera0.y(), 1.0));
      brick::geometry::Ray3D<FloatType> ray1_1(
        brick::numeric::Vector3D<FloatType>(0.0, 0.0, 0.0),
        brick::numeric::Vector3D<FloatType>(testPointCamera1.x(), testPointCamera1.y(), 1.0));

      // Transform one of the two rays from camera1 coordinates to
      // camera0 coordinates.
      brick::geometry::Ray3D<FloatType> ray1_0 = c0Tc1 * ray1_1;

      // And find the closest point of approach between the two rays.
      FloatType distance0;
      FloatType distance1;
      FloatType residual;
      return brick::geometry::findIntersect(
        ray0_0, ray1_0, distance0, distance1, residual);
    }


    // This function is used internally by fivePointAlgorithm() to
    // generate a 10x20 matrix of coefficients of polynomial
    // constraints.
    template <class FloatType>
    brick::numeric::Array2D<FloatType>
    generateFivePointConstraintMatrix(
      brick::numeric::Array2D<FloatType> const& E0Array,
      brick::numeric::Array2D<FloatType> const& E1Array,
      brick::numeric::Array2D<FloatType> const& E2Array,
      brick::numeric::Array2D<FloatType> const& E3Array)
    {

      /*
        =====================================================
        ============== Start of long comment ================
        =====================================================

        # Here's maxima code that automatically generates our
        # constraint equations.  Save this as maximaCode.ma, and
        # execute the command
        #
        #   maxima -b maximaCode.ma
        #
        # to generate an intermedate text file called
        # maximaConstraints.txt, then procede to the python code
        # below.

        # ---------------- Start of Maxima code ------------------

        E0:matrix([e0rc00, e0rc01, e0rc02],
                  [e0rc10, e0rc11, e0rc12],
                  [e0rc20, e0rc21, e0rc22]);
        E1:matrix([e1rc00, e1rc01, e1rc02],
                  [e1rc10, e1rc11, e1rc12],
                  [e1rc20, e1rc21, e1rc22]);
        E2:matrix([e2rc00, e2rc01, e2rc02],
                  [e2rc10, e2rc11, e2rc12],
                  [e2rc20, e2rc21, e2rc22]);
        E3:matrix([e3rc00, e3rc01, e3rc02],
                  [e3rc10, e3rc11, e3rc12],
                  [e3rc20, e3rc21, e3rc22]);
        EE:x0*E0 + x1*E1 + x2*E2 + E3;

        C0:determinant(EE);
        C0E:expand(C0);

        EExEET:EE . transpose(EE);
        ETrace:EExEET[1][1] + EExEET[2][2] + EExEET[3][3];
        CC:2 * (EExEET . EE) - (ETrace * EE);
        CCE:expand(CC);

        Constraints:[CCE[1][1], CCE[1][2], CCE[1][3],
                     CCE[2][1], CCE[2][2], CCE[2][3],
                     CCE[3][1], CCE[3][2], CCE[3][3],
                     C0E];

        appendfile("maximaConstraints.txt");
        grind(Constraints);

        # ----------------- End of Maxima code -------------------

        # Here's python code that parses maximaConstraints.txt
        # (generated above) to create code for computing the
        # constraints matrix.  Save this code as pythonCode.py, and
        # execute the command
        #
        #   python pythonCode.py > cppSnippet.cpp
        #
        # Then copy the contents of cppSnippet.cpp into the cut &
        # paste section below.

        # ---------------- Start of Python code ------------------
        import re
        constraintsRe = re.compile("\[.*\]", re.DOTALL)

        # This regex works as follows:
        #   (                  - Starts a group representing the coefficient.
        #     [+-]?              - Optional prefix of '+' or '-', as in
        #                          '-2*e0rc00'.
        #     (?:\d\*)?          - Optional leading multiplier, as in
        #                          '2*e0rc00'.
        #     e[e0123rc^*]*[0123] - A product of e?rc?? terms, as in
        #                           'e0rc00^2*e1rc21'.
        #   )                  - Ends the group representing the coefficient.
        #   (?:                - Starts a group.
        #     \*                 - Group has to start with a "*", as in
        #                          '*x0^3'.
        #     (x[x0123*^]*)      - A group representing the monomial, as in
        #                          'x0^2*x1'.
        #   )?                 - Closes the group, and makes it optional.
        termRe = re.compile(
          "([+-]?(?:\d\*)?e[e0123rc^*]*[0123])(?:\*(x[x0123*^]*))?")
        expectedMonomials = ['x0^3', 'x0^2*x1', 'x0^2*x2',
                             'x0*x1^2', 'x0*x1*x2', 'x0*x2^2',
                             'x1^3', 'x1^2*x2', 'x1*x2^2',
                             'x2^3',
                             'x0^2', 'x0*x1', 'x0*x2',
                             'x1^2', 'x1*x2', 'x2^2',
                             'x0', 'x1', 'x2',
                             None]

        inputFile = open('maximaConstraints.txt')
        inputText = inputFile.read()

        # This array represents the 10x20 constraint matrix.
        coefficientArray = map(lambda x: [''] * 20, range(10))

        # We're expecting to find 10 equations: 9 polynomials for the nine
        # terms in the "2 * E * transpose(E) * E ..." constraint, and one
        # polynomial for the determinant constraint.
        constraintsText = constraintsRe.search(inputText).group(0)[1:-1]
        constraintsList = constraintsText.split(',')

        # Iterate over all constraints.
        constraintNumber = 0
        for constraint in constraintsList:
          coefficientMap = {}

          # Reset this row of the constraint matrix.
          for monomial in expectedMonomials:
            coefficientMap[monomial] = ''
          # end for

          # Examine all terms of this polynomial, split each term into
          # coefficient and monomial parts (where coefficient might be
          # something like "e0rc12^2*e2rc21", and monomial might be
          # something like "x1*x2^2"), and add the coefficient to the
          # total coefficient for the monomial (since there may be
          # more than one term associated with each monomial).
          startPos = 0
          while 1:
            termMatch = termRe.search(constraint, startPos)
            if not termMatch:
              break
            # end if
            coefficient = termMatch.group(1)
            monomial = termMatch.group(2)
            if not coefficientMap.has_key(monomial):
              raise IOError, 'Unexpected monomial'
            # end if
            coefficientMap[monomial] = (
              coefficientMap[monomial] + coefficient.replace('^', '_'))
            startPos = termMatch.end()
          # end while

          # Now put the accumulated coefficients into the constraint matrix.
          for ii in range(len(expectedMonomials)):
            coefficient = coefficientMap[expectedMonomials[ii]]
            if coefficient != '':
              coefficientArray[constraintNumber][ii] = coefficient
            else:
              coefficientArray[constraintNumber][ii] = '0.0'
            # end if
          # end for
          constraintNumber = constraintNumber + 1
        # end for

        # Send C++ code to standard out.
        for ii in range(len(coefficientArray)):
          for jj in range(len(expectedMonomials)):
            print ('    AMatrix(%d, %d) = %s;'
                   % (ii, jj, coefficientArray[ii][jj]))
          # end for
        # end for

        # ----------------- End of Python code -------------------

        =====================================================
        ==============  End of long comment  ================
        =====================================================
       */

      // First set up the variables for the auto-generated code...

      // Individual matrix elements.  The name eXrcYZ means the (Y, Z)
      // element of matrix E'X'Array.
      FloatType e0rc00 = E0Array(0, 0);
      FloatType e0rc01 = E0Array(0, 1);
      FloatType e0rc02 = E0Array(0, 2);
      FloatType e0rc10 = E0Array(1, 0);
      FloatType e0rc11 = E0Array(1, 1);
      FloatType e0rc12 = E0Array(1, 2);
      FloatType e0rc20 = E0Array(2, 0);
      FloatType e0rc21 = E0Array(2, 1);
      FloatType e0rc22 = E0Array(2, 2);
      FloatType e1rc00 = E1Array(0, 0);
      FloatType e1rc01 = E1Array(0, 1);
      FloatType e1rc02 = E1Array(0, 2);
      FloatType e1rc10 = E1Array(1, 0);
      FloatType e1rc11 = E1Array(1, 1);
      FloatType e1rc12 = E1Array(1, 2);
      FloatType e1rc20 = E1Array(2, 0);
      FloatType e1rc21 = E1Array(2, 1);
      FloatType e1rc22 = E1Array(2, 2);
      FloatType e2rc00 = E2Array(0, 0);
      FloatType e2rc01 = E2Array(0, 1);
      FloatType e2rc02 = E2Array(0, 2);
      FloatType e2rc10 = E2Array(1, 0);
      FloatType e2rc11 = E2Array(1, 1);
      FloatType e2rc12 = E2Array(1, 2);
      FloatType e2rc20 = E2Array(2, 0);
      FloatType e2rc21 = E2Array(2, 1);
      FloatType e2rc22 = E2Array(2, 2);
      FloatType e3rc00 = E3Array(0, 0);
      FloatType e3rc01 = E3Array(0, 1);
      FloatType e3rc02 = E3Array(0, 2);
      FloatType e3rc10 = E3Array(1, 0);
      FloatType e3rc11 = E3Array(1, 1);
      FloatType e3rc12 = E3Array(1, 2);
      FloatType e3rc20 = E3Array(2, 0);
      FloatType e3rc21 = E3Array(2, 1);
      FloatType e3rc22 = E3Array(2, 2);

      // Squares of individual matrix elements.
      FloatType e0rc00_2 = e0rc00 * e0rc00;
      FloatType e0rc01_2 = e0rc01 * e0rc01;
      FloatType e0rc02_2 = e0rc02 * e0rc02;
      FloatType e0rc10_2 = e0rc10 * e0rc10;
      FloatType e0rc11_2 = e0rc11 * e0rc11;
      FloatType e0rc12_2 = e0rc12 * e0rc12;
      FloatType e0rc20_2 = e0rc20 * e0rc20;
      FloatType e0rc21_2 = e0rc21 * e0rc21;
      FloatType e0rc22_2 = e0rc22 * e0rc22;
      FloatType e1rc00_2 = e1rc00 * e1rc00;
      FloatType e1rc01_2 = e1rc01 * e1rc01;
      FloatType e1rc02_2 = e1rc02 * e1rc02;
      FloatType e1rc10_2 = e1rc10 * e1rc10;
      FloatType e1rc11_2 = e1rc11 * e1rc11;
      FloatType e1rc12_2 = e1rc12 * e1rc12;
      FloatType e1rc20_2 = e1rc20 * e1rc20;
      FloatType e1rc21_2 = e1rc21 * e1rc21;
      FloatType e1rc22_2 = e1rc22 * e1rc22;
      FloatType e2rc00_2 = e2rc00 * e2rc00;
      FloatType e2rc01_2 = e2rc01 * e2rc01;
      FloatType e2rc02_2 = e2rc02 * e2rc02;
      FloatType e2rc10_2 = e2rc10 * e2rc10;
      FloatType e2rc11_2 = e2rc11 * e2rc11;
      FloatType e2rc12_2 = e2rc12 * e2rc12;
      FloatType e2rc20_2 = e2rc20 * e2rc20;
      FloatType e2rc21_2 = e2rc21 * e2rc21;
      FloatType e2rc22_2 = e2rc22 * e2rc22;
      FloatType e3rc00_2 = e3rc00 * e3rc00;
      FloatType e3rc01_2 = e3rc01 * e3rc01;
      FloatType e3rc02_2 = e3rc02 * e3rc02;
      FloatType e3rc10_2 = e3rc10 * e3rc10;
      FloatType e3rc11_2 = e3rc11 * e3rc11;
      FloatType e3rc12_2 = e3rc12 * e3rc12;
      FloatType e3rc20_2 = e3rc20 * e3rc20;
      FloatType e3rc21_2 = e3rc21 * e3rc21;
      FloatType e3rc22_2 = e3rc22 * e3rc22;

      // Cubes of individual matrix elements.
      FloatType e0rc00_3 = e0rc00_2 * e0rc00;
      FloatType e0rc01_3 = e0rc01_2 * e0rc01;
      FloatType e0rc02_3 = e0rc02_2 * e0rc02;
      FloatType e0rc10_3 = e0rc10_2 * e0rc10;
      FloatType e0rc11_3 = e0rc11_2 * e0rc11;
      FloatType e0rc12_3 = e0rc12_2 * e0rc12;
      FloatType e0rc20_3 = e0rc20_2 * e0rc20;
      FloatType e0rc21_3 = e0rc21_2 * e0rc21;
      FloatType e0rc22_3 = e0rc22_2 * e0rc22;
      FloatType e1rc00_3 = e1rc00_2 * e1rc00;
      FloatType e1rc01_3 = e1rc01_2 * e1rc01;
      FloatType e1rc02_3 = e1rc02_2 * e1rc02;
      FloatType e1rc10_3 = e1rc10_2 * e1rc10;
      FloatType e1rc11_3 = e1rc11_2 * e1rc11;
      FloatType e1rc12_3 = e1rc12_2 * e1rc12;
      FloatType e1rc20_3 = e1rc20_2 * e1rc20;
      FloatType e1rc21_3 = e1rc21_2 * e1rc21;
      FloatType e1rc22_3 = e1rc22_2 * e1rc22;
      FloatType e2rc00_3 = e2rc00_2 * e2rc00;
      FloatType e2rc01_3 = e2rc01_2 * e2rc01;
      FloatType e2rc02_3 = e2rc02_2 * e2rc02;
      FloatType e2rc10_3 = e2rc10_2 * e2rc10;
      FloatType e2rc11_3 = e2rc11_2 * e2rc11;
      FloatType e2rc12_3 = e2rc12_2 * e2rc12;
      FloatType e2rc20_3 = e2rc20_2 * e2rc20;
      FloatType e2rc21_3 = e2rc21_2 * e2rc21;
      FloatType e2rc22_3 = e2rc22_2 * e2rc22;
      FloatType e3rc00_3 = e3rc00_2 * e3rc00;
      FloatType e3rc01_3 = e3rc01_2 * e3rc01;
      FloatType e3rc02_3 = e3rc02_2 * e3rc02;
      FloatType e3rc10_3 = e3rc10_2 * e3rc10;
      FloatType e3rc11_3 = e3rc11_2 * e3rc11;
      FloatType e3rc12_3 = e3rc12_2 * e3rc12;
      FloatType e3rc20_3 = e3rc20_2 * e3rc20;
      FloatType e3rc21_3 = e3rc21_2 * e3rc21;
      FloatType e3rc22_3 = e3rc22_2 * e3rc22;

      // Big matrix to accept all of the constraint coefficients.
      brick::numeric::Array2D<FloatType> AMatrix(10, 20);


      // ============= Begin cut & paste section =============

    AMatrix(0, 0) = -e0rc00*e0rc22_2+2*e0rc02*e0rc20*e0rc22-e0rc00*e0rc21_2+2*e0rc01*e0rc20*e0rc21+e0rc00*e0rc20_2-e0rc00*e0rc12_2+2*e0rc02*e0rc10*e0rc12-e0rc00*e0rc11_2+2*e0rc01*e0rc10*e0rc11+e0rc00*e0rc10_2+e0rc00*e0rc02_2+e0rc00*e0rc01_2+e0rc00_3;
    AMatrix(0, 1) = -2*e0rc00*e0rc22*e1rc22+2*e0rc02*e0rc20*e1rc22-2*e0rc00*e0rc21*e1rc21+2*e0rc01*e0rc20*e1rc21+2*e0rc02*e0rc22*e1rc20+2*e0rc01*e0rc21*e1rc20+2*e0rc00*e0rc20*e1rc20-2*e0rc00*e0rc12*e1rc12+2*e0rc02*e0rc10*e1rc12-2*e0rc00*e0rc11*e1rc11+2*e0rc01*e0rc10*e1rc11+2*e0rc02*e0rc12*e1rc10+2*e0rc01*e0rc11*e1rc10+2*e0rc00*e0rc10*e1rc10+2*e0rc20*e0rc22*e1rc02+2*e0rc10*e0rc12*e1rc02+2*e0rc00*e0rc02*e1rc02+2*e0rc20*e0rc21*e1rc01+2*e0rc10*e0rc11*e1rc01+2*e0rc00*e0rc01*e1rc01-e0rc22_2*e1rc00-e0rc21_2*e1rc00+e0rc20_2*e1rc00-e0rc12_2*e1rc00-e0rc11_2*e1rc00+e0rc10_2*e1rc00+e0rc02_2*e1rc00+e0rc01_2*e1rc00+3*e0rc00_2*e1rc00;
    AMatrix(0, 2) = -2*e0rc00*e0rc22*e2rc22+2*e0rc02*e0rc20*e2rc22-2*e0rc00*e0rc21*e2rc21+2*e0rc01*e0rc20*e2rc21+2*e0rc02*e0rc22*e2rc20+2*e0rc01*e0rc21*e2rc20+2*e0rc00*e0rc20*e2rc20-2*e0rc00*e0rc12*e2rc12+2*e0rc02*e0rc10*e2rc12-2*e0rc00*e0rc11*e2rc11+2*e0rc01*e0rc10*e2rc11+2*e0rc02*e0rc12*e2rc10+2*e0rc01*e0rc11*e2rc10+2*e0rc00*e0rc10*e2rc10+2*e0rc20*e0rc22*e2rc02+2*e0rc10*e0rc12*e2rc02+2*e0rc00*e0rc02*e2rc02+2*e0rc20*e0rc21*e2rc01+2*e0rc10*e0rc11*e2rc01+2*e0rc00*e0rc01*e2rc01-e0rc22_2*e2rc00-e0rc21_2*e2rc00+e0rc20_2*e2rc00-e0rc12_2*e2rc00-e0rc11_2*e2rc00+e0rc10_2*e2rc00+e0rc02_2*e2rc00+e0rc01_2*e2rc00+3*e0rc00_2*e2rc00;
    AMatrix(0, 3) = -e0rc00*e1rc22_2+2*e0rc02*e1rc20*e1rc22+2*e0rc20*e1rc02*e1rc22-2*e0rc22*e1rc00*e1rc22-e0rc00*e1rc21_2+2*e0rc01*e1rc20*e1rc21+2*e0rc20*e1rc01*e1rc21-2*e0rc21*e1rc00*e1rc21+e0rc00*e1rc20_2+2*e0rc22*e1rc02*e1rc20+2*e0rc21*e1rc01*e1rc20+2*e0rc20*e1rc00*e1rc20-e0rc00*e1rc12_2+2*e0rc02*e1rc10*e1rc12+2*e0rc10*e1rc02*e1rc12-2*e0rc12*e1rc00*e1rc12-e0rc00*e1rc11_2+2*e0rc01*e1rc10*e1rc11+2*e0rc10*e1rc01*e1rc11-2*e0rc11*e1rc00*e1rc11+e0rc00*e1rc10_2+2*e0rc12*e1rc02*e1rc10+2*e0rc11*e1rc01*e1rc10+2*e0rc10*e1rc00*e1rc10+e0rc00*e1rc02_2+2*e0rc02*e1rc00*e1rc02+e0rc00*e1rc01_2+2*e0rc01*e1rc00*e1rc01+3*e0rc00*e1rc00_2;
    AMatrix(0, 4) = -2*e0rc00*e1rc22*e2rc22+2*e0rc02*e1rc20*e2rc22+2*e0rc20*e1rc02*e2rc22-2*e0rc22*e1rc00*e2rc22-2*e0rc00*e1rc21*e2rc21+2*e0rc01*e1rc20*e2rc21+2*e0rc20*e1rc01*e2rc21-2*e0rc21*e1rc00*e2rc21+2*e0rc02*e1rc22*e2rc20+2*e0rc01*e1rc21*e2rc20+2*e0rc00*e1rc20*e2rc20+2*e0rc22*e1rc02*e2rc20+2*e0rc21*e1rc01*e2rc20+2*e0rc20*e1rc00*e2rc20-2*e0rc00*e1rc12*e2rc12+2*e0rc02*e1rc10*e2rc12+2*e0rc10*e1rc02*e2rc12-2*e0rc12*e1rc00*e2rc12-2*e0rc00*e1rc11*e2rc11+2*e0rc01*e1rc10*e2rc11+2*e0rc10*e1rc01*e2rc11-2*e0rc11*e1rc00*e2rc11+2*e0rc02*e1rc12*e2rc10+2*e0rc01*e1rc11*e2rc10+2*e0rc00*e1rc10*e2rc10+2*e0rc12*e1rc02*e2rc10+2*e0rc11*e1rc01*e2rc10+2*e0rc10*e1rc00*e2rc10+2*e0rc20*e1rc22*e2rc02+2*e0rc22*e1rc20*e2rc02+2*e0rc10*e1rc12*e2rc02+2*e0rc12*e1rc10*e2rc02+2*e0rc00*e1rc02*e2rc02+2*e0rc02*e1rc00*e2rc02+2*e0rc20*e1rc21*e2rc01+2*e0rc21*e1rc20*e2rc01+2*e0rc10*e1rc11*e2rc01+2*e0rc11*e1rc10*e2rc01+2*e0rc00*e1rc01*e2rc01+2*e0rc01*e1rc00*e2rc01-2*e0rc22*e1rc22*e2rc00-2*e0rc21*e1rc21*e2rc00+2*e0rc20*e1rc20*e2rc00-2*e0rc12*e1rc12*e2rc00-2*e0rc11*e1rc11*e2rc00+2*e0rc10*e1rc10*e2rc00+2*e0rc02*e1rc02*e2rc00+2*e0rc01*e1rc01*e2rc00+6*e0rc00*e1rc00*e2rc00;
    AMatrix(0, 5) = -e0rc00*e2rc22_2+2*e0rc02*e2rc20*e2rc22+2*e0rc20*e2rc02*e2rc22-2*e0rc22*e2rc00*e2rc22-e0rc00*e2rc21_2+2*e0rc01*e2rc20*e2rc21+2*e0rc20*e2rc01*e2rc21-2*e0rc21*e2rc00*e2rc21+e0rc00*e2rc20_2+2*e0rc22*e2rc02*e2rc20+2*e0rc21*e2rc01*e2rc20+2*e0rc20*e2rc00*e2rc20-e0rc00*e2rc12_2+2*e0rc02*e2rc10*e2rc12+2*e0rc10*e2rc02*e2rc12-2*e0rc12*e2rc00*e2rc12-e0rc00*e2rc11_2+2*e0rc01*e2rc10*e2rc11+2*e0rc10*e2rc01*e2rc11-2*e0rc11*e2rc00*e2rc11+e0rc00*e2rc10_2+2*e0rc12*e2rc02*e2rc10+2*e0rc11*e2rc01*e2rc10+2*e0rc10*e2rc00*e2rc10+e0rc00*e2rc02_2+2*e0rc02*e2rc00*e2rc02+e0rc00*e2rc01_2+2*e0rc01*e2rc00*e2rc01+3*e0rc00*e2rc00_2;
    AMatrix(0, 6) = -e1rc00*e1rc22_2+2*e1rc02*e1rc20*e1rc22-e1rc00*e1rc21_2+2*e1rc01*e1rc20*e1rc21+e1rc00*e1rc20_2-e1rc00*e1rc12_2+2*e1rc02*e1rc10*e1rc12-e1rc00*e1rc11_2+2*e1rc01*e1rc10*e1rc11+e1rc00*e1rc10_2+e1rc00*e1rc02_2+e1rc00*e1rc01_2+e1rc00_3;
    AMatrix(0, 7) = -2*e1rc00*e1rc22*e2rc22+2*e1rc02*e1rc20*e2rc22-2*e1rc00*e1rc21*e2rc21+2*e1rc01*e1rc20*e2rc21+2*e1rc02*e1rc22*e2rc20+2*e1rc01*e1rc21*e2rc20+2*e1rc00*e1rc20*e2rc20-2*e1rc00*e1rc12*e2rc12+2*e1rc02*e1rc10*e2rc12-2*e1rc00*e1rc11*e2rc11+2*e1rc01*e1rc10*e2rc11+2*e1rc02*e1rc12*e2rc10+2*e1rc01*e1rc11*e2rc10+2*e1rc00*e1rc10*e2rc10+2*e1rc20*e1rc22*e2rc02+2*e1rc10*e1rc12*e2rc02+2*e1rc00*e1rc02*e2rc02+2*e1rc20*e1rc21*e2rc01+2*e1rc10*e1rc11*e2rc01+2*e1rc00*e1rc01*e2rc01-e1rc22_2*e2rc00-e1rc21_2*e2rc00+e1rc20_2*e2rc00-e1rc12_2*e2rc00-e1rc11_2*e2rc00+e1rc10_2*e2rc00+e1rc02_2*e2rc00+e1rc01_2*e2rc00+3*e1rc00_2*e2rc00;
    AMatrix(0, 8) = -e1rc00*e2rc22_2+2*e1rc02*e2rc20*e2rc22+2*e1rc20*e2rc02*e2rc22-2*e1rc22*e2rc00*e2rc22-e1rc00*e2rc21_2+2*e1rc01*e2rc20*e2rc21+2*e1rc20*e2rc01*e2rc21-2*e1rc21*e2rc00*e2rc21+e1rc00*e2rc20_2+2*e1rc22*e2rc02*e2rc20+2*e1rc21*e2rc01*e2rc20+2*e1rc20*e2rc00*e2rc20-e1rc00*e2rc12_2+2*e1rc02*e2rc10*e2rc12+2*e1rc10*e2rc02*e2rc12-2*e1rc12*e2rc00*e2rc12-e1rc00*e2rc11_2+2*e1rc01*e2rc10*e2rc11+2*e1rc10*e2rc01*e2rc11-2*e1rc11*e2rc00*e2rc11+e1rc00*e2rc10_2+2*e1rc12*e2rc02*e2rc10+2*e1rc11*e2rc01*e2rc10+2*e1rc10*e2rc00*e2rc10+e1rc00*e2rc02_2+2*e1rc02*e2rc00*e2rc02+e1rc00*e2rc01_2+2*e1rc01*e2rc00*e2rc01+3*e1rc00*e2rc00_2;
    AMatrix(0, 9) = -e2rc00*e2rc22_2+2*e2rc02*e2rc20*e2rc22-e2rc00*e2rc21_2+2*e2rc01*e2rc20*e2rc21+e2rc00*e2rc20_2-e2rc00*e2rc12_2+2*e2rc02*e2rc10*e2rc12-e2rc00*e2rc11_2+2*e2rc01*e2rc10*e2rc11+e2rc00*e2rc10_2+e2rc00*e2rc02_2+e2rc00*e2rc01_2+e2rc00_3;
    AMatrix(0, 10) = -2*e0rc00*e0rc22*e3rc22+2*e0rc02*e0rc20*e3rc22-2*e0rc00*e0rc21*e3rc21+2*e0rc01*e0rc20*e3rc21+2*e0rc02*e0rc22*e3rc20+2*e0rc01*e0rc21*e3rc20+2*e0rc00*e0rc20*e3rc20-2*e0rc00*e0rc12*e3rc12+2*e0rc02*e0rc10*e3rc12-2*e0rc00*e0rc11*e3rc11+2*e0rc01*e0rc10*e3rc11+2*e0rc02*e0rc12*e3rc10+2*e0rc01*e0rc11*e3rc10+2*e0rc00*e0rc10*e3rc10+2*e0rc20*e0rc22*e3rc02+2*e0rc10*e0rc12*e3rc02+2*e0rc00*e0rc02*e3rc02+2*e0rc20*e0rc21*e3rc01+2*e0rc10*e0rc11*e3rc01+2*e0rc00*e0rc01*e3rc01-e0rc22_2*e3rc00-e0rc21_2*e3rc00+e0rc20_2*e3rc00-e0rc12_2*e3rc00-e0rc11_2*e3rc00+e0rc10_2*e3rc00+e0rc02_2*e3rc00+e0rc01_2*e3rc00+3*e0rc00_2*e3rc00;
    AMatrix(0, 11) = -2*e0rc00*e1rc22*e3rc22+2*e0rc02*e1rc20*e3rc22+2*e0rc20*e1rc02*e3rc22-2*e0rc22*e1rc00*e3rc22-2*e0rc00*e1rc21*e3rc21+2*e0rc01*e1rc20*e3rc21+2*e0rc20*e1rc01*e3rc21-2*e0rc21*e1rc00*e3rc21+2*e0rc02*e1rc22*e3rc20+2*e0rc01*e1rc21*e3rc20+2*e0rc00*e1rc20*e3rc20+2*e0rc22*e1rc02*e3rc20+2*e0rc21*e1rc01*e3rc20+2*e0rc20*e1rc00*e3rc20-2*e0rc00*e1rc12*e3rc12+2*e0rc02*e1rc10*e3rc12+2*e0rc10*e1rc02*e3rc12-2*e0rc12*e1rc00*e3rc12-2*e0rc00*e1rc11*e3rc11+2*e0rc01*e1rc10*e3rc11+2*e0rc10*e1rc01*e3rc11-2*e0rc11*e1rc00*e3rc11+2*e0rc02*e1rc12*e3rc10+2*e0rc01*e1rc11*e3rc10+2*e0rc00*e1rc10*e3rc10+2*e0rc12*e1rc02*e3rc10+2*e0rc11*e1rc01*e3rc10+2*e0rc10*e1rc00*e3rc10+2*e0rc20*e1rc22*e3rc02+2*e0rc22*e1rc20*e3rc02+2*e0rc10*e1rc12*e3rc02+2*e0rc12*e1rc10*e3rc02+2*e0rc00*e1rc02*e3rc02+2*e0rc02*e1rc00*e3rc02+2*e0rc20*e1rc21*e3rc01+2*e0rc21*e1rc20*e3rc01+2*e0rc10*e1rc11*e3rc01+2*e0rc11*e1rc10*e3rc01+2*e0rc00*e1rc01*e3rc01+2*e0rc01*e1rc00*e3rc01-2*e0rc22*e1rc22*e3rc00-2*e0rc21*e1rc21*e3rc00+2*e0rc20*e1rc20*e3rc00-2*e0rc12*e1rc12*e3rc00-2*e0rc11*e1rc11*e3rc00+2*e0rc10*e1rc10*e3rc00+2*e0rc02*e1rc02*e3rc00+2*e0rc01*e1rc01*e3rc00+6*e0rc00*e1rc00*e3rc00;
    AMatrix(0, 12) = -2*e0rc00*e2rc22*e3rc22+2*e0rc02*e2rc20*e3rc22+2*e0rc20*e2rc02*e3rc22-2*e0rc22*e2rc00*e3rc22-2*e0rc00*e2rc21*e3rc21+2*e0rc01*e2rc20*e3rc21+2*e0rc20*e2rc01*e3rc21-2*e0rc21*e2rc00*e3rc21+2*e0rc02*e2rc22*e3rc20+2*e0rc01*e2rc21*e3rc20+2*e0rc00*e2rc20*e3rc20+2*e0rc22*e2rc02*e3rc20+2*e0rc21*e2rc01*e3rc20+2*e0rc20*e2rc00*e3rc20-2*e0rc00*e2rc12*e3rc12+2*e0rc02*e2rc10*e3rc12+2*e0rc10*e2rc02*e3rc12-2*e0rc12*e2rc00*e3rc12-2*e0rc00*e2rc11*e3rc11+2*e0rc01*e2rc10*e3rc11+2*e0rc10*e2rc01*e3rc11-2*e0rc11*e2rc00*e3rc11+2*e0rc02*e2rc12*e3rc10+2*e0rc01*e2rc11*e3rc10+2*e0rc00*e2rc10*e3rc10+2*e0rc12*e2rc02*e3rc10+2*e0rc11*e2rc01*e3rc10+2*e0rc10*e2rc00*e3rc10+2*e0rc20*e2rc22*e3rc02+2*e0rc22*e2rc20*e3rc02+2*e0rc10*e2rc12*e3rc02+2*e0rc12*e2rc10*e3rc02+2*e0rc00*e2rc02*e3rc02+2*e0rc02*e2rc00*e3rc02+2*e0rc20*e2rc21*e3rc01+2*e0rc21*e2rc20*e3rc01+2*e0rc10*e2rc11*e3rc01+2*e0rc11*e2rc10*e3rc01+2*e0rc00*e2rc01*e3rc01+2*e0rc01*e2rc00*e3rc01-2*e0rc22*e2rc22*e3rc00-2*e0rc21*e2rc21*e3rc00+2*e0rc20*e2rc20*e3rc00-2*e0rc12*e2rc12*e3rc00-2*e0rc11*e2rc11*e3rc00+2*e0rc10*e2rc10*e3rc00+2*e0rc02*e2rc02*e3rc00+2*e0rc01*e2rc01*e3rc00+6*e0rc00*e2rc00*e3rc00;
    AMatrix(0, 13) = -2*e1rc00*e1rc22*e3rc22+2*e1rc02*e1rc20*e3rc22-2*e1rc00*e1rc21*e3rc21+2*e1rc01*e1rc20*e3rc21+2*e1rc02*e1rc22*e3rc20+2*e1rc01*e1rc21*e3rc20+2*e1rc00*e1rc20*e3rc20-2*e1rc00*e1rc12*e3rc12+2*e1rc02*e1rc10*e3rc12-2*e1rc00*e1rc11*e3rc11+2*e1rc01*e1rc10*e3rc11+2*e1rc02*e1rc12*e3rc10+2*e1rc01*e1rc11*e3rc10+2*e1rc00*e1rc10*e3rc10+2*e1rc20*e1rc22*e3rc02+2*e1rc10*e1rc12*e3rc02+2*e1rc00*e1rc02*e3rc02+2*e1rc20*e1rc21*e3rc01+2*e1rc10*e1rc11*e3rc01+2*e1rc00*e1rc01*e3rc01-e1rc22_2*e3rc00-e1rc21_2*e3rc00+e1rc20_2*e3rc00-e1rc12_2*e3rc00-e1rc11_2*e3rc00+e1rc10_2*e3rc00+e1rc02_2*e3rc00+e1rc01_2*e3rc00+3*e1rc00_2*e3rc00;
    AMatrix(0, 14) = -2*e1rc00*e2rc22*e3rc22+2*e1rc02*e2rc20*e3rc22+2*e1rc20*e2rc02*e3rc22-2*e1rc22*e2rc00*e3rc22-2*e1rc00*e2rc21*e3rc21+2*e1rc01*e2rc20*e3rc21+2*e1rc20*e2rc01*e3rc21-2*e1rc21*e2rc00*e3rc21+2*e1rc02*e2rc22*e3rc20+2*e1rc01*e2rc21*e3rc20+2*e1rc00*e2rc20*e3rc20+2*e1rc22*e2rc02*e3rc20+2*e1rc21*e2rc01*e3rc20+2*e1rc20*e2rc00*e3rc20-2*e1rc00*e2rc12*e3rc12+2*e1rc02*e2rc10*e3rc12+2*e1rc10*e2rc02*e3rc12-2*e1rc12*e2rc00*e3rc12-2*e1rc00*e2rc11*e3rc11+2*e1rc01*e2rc10*e3rc11+2*e1rc10*e2rc01*e3rc11-2*e1rc11*e2rc00*e3rc11+2*e1rc02*e2rc12*e3rc10+2*e1rc01*e2rc11*e3rc10+2*e1rc00*e2rc10*e3rc10+2*e1rc12*e2rc02*e3rc10+2*e1rc11*e2rc01*e3rc10+2*e1rc10*e2rc00*e3rc10+2*e1rc20*e2rc22*e3rc02+2*e1rc22*e2rc20*e3rc02+2*e1rc10*e2rc12*e3rc02+2*e1rc12*e2rc10*e3rc02+2*e1rc00*e2rc02*e3rc02+2*e1rc02*e2rc00*e3rc02+2*e1rc20*e2rc21*e3rc01+2*e1rc21*e2rc20*e3rc01+2*e1rc10*e2rc11*e3rc01+2*e1rc11*e2rc10*e3rc01+2*e1rc00*e2rc01*e3rc01+2*e1rc01*e2rc00*e3rc01-2*e1rc22*e2rc22*e3rc00-2*e1rc21*e2rc21*e3rc00+2*e1rc20*e2rc20*e3rc00-2*e1rc12*e2rc12*e3rc00-2*e1rc11*e2rc11*e3rc00+2*e1rc10*e2rc10*e3rc00+2*e1rc02*e2rc02*e3rc00+2*e1rc01*e2rc01*e3rc00+6*e1rc00*e2rc00*e3rc00;
    AMatrix(0, 15) = -2*e2rc00*e2rc22*e3rc22+2*e2rc02*e2rc20*e3rc22-2*e2rc00*e2rc21*e3rc21+2*e2rc01*e2rc20*e3rc21+2*e2rc02*e2rc22*e3rc20+2*e2rc01*e2rc21*e3rc20+2*e2rc00*e2rc20*e3rc20-2*e2rc00*e2rc12*e3rc12+2*e2rc02*e2rc10*e3rc12-2*e2rc00*e2rc11*e3rc11+2*e2rc01*e2rc10*e3rc11+2*e2rc02*e2rc12*e3rc10+2*e2rc01*e2rc11*e3rc10+2*e2rc00*e2rc10*e3rc10+2*e2rc20*e2rc22*e3rc02+2*e2rc10*e2rc12*e3rc02+2*e2rc00*e2rc02*e3rc02+2*e2rc20*e2rc21*e3rc01+2*e2rc10*e2rc11*e3rc01+2*e2rc00*e2rc01*e3rc01-e2rc22_2*e3rc00-e2rc21_2*e3rc00+e2rc20_2*e3rc00-e2rc12_2*e3rc00-e2rc11_2*e3rc00+e2rc10_2*e3rc00+e2rc02_2*e3rc00+e2rc01_2*e3rc00+3*e2rc00_2*e3rc00;
    AMatrix(0, 16) = -e0rc00*e3rc22_2+2*e0rc02*e3rc20*e3rc22+2*e0rc20*e3rc02*e3rc22-2*e0rc22*e3rc00*e3rc22-e0rc00*e3rc21_2+2*e0rc01*e3rc20*e3rc21+2*e0rc20*e3rc01*e3rc21-2*e0rc21*e3rc00*e3rc21+e0rc00*e3rc20_2+2*e0rc22*e3rc02*e3rc20+2*e0rc21*e3rc01*e3rc20+2*e0rc20*e3rc00*e3rc20-e0rc00*e3rc12_2+2*e0rc02*e3rc10*e3rc12+2*e0rc10*e3rc02*e3rc12-2*e0rc12*e3rc00*e3rc12-e0rc00*e3rc11_2+2*e0rc01*e3rc10*e3rc11+2*e0rc10*e3rc01*e3rc11-2*e0rc11*e3rc00*e3rc11+e0rc00*e3rc10_2+2*e0rc12*e3rc02*e3rc10+2*e0rc11*e3rc01*e3rc10+2*e0rc10*e3rc00*e3rc10+e0rc00*e3rc02_2+2*e0rc02*e3rc00*e3rc02+e0rc00*e3rc01_2+2*e0rc01*e3rc00*e3rc01+3*e0rc00*e3rc00_2;
    AMatrix(0, 17) = -e1rc00*e3rc22_2+2*e1rc02*e3rc20*e3rc22+2*e1rc20*e3rc02*e3rc22-2*e1rc22*e3rc00*e3rc22-e1rc00*e3rc21_2+2*e1rc01*e3rc20*e3rc21+2*e1rc20*e3rc01*e3rc21-2*e1rc21*e3rc00*e3rc21+e1rc00*e3rc20_2+2*e1rc22*e3rc02*e3rc20+2*e1rc21*e3rc01*e3rc20+2*e1rc20*e3rc00*e3rc20-e1rc00*e3rc12_2+2*e1rc02*e3rc10*e3rc12+2*e1rc10*e3rc02*e3rc12-2*e1rc12*e3rc00*e3rc12-e1rc00*e3rc11_2+2*e1rc01*e3rc10*e3rc11+2*e1rc10*e3rc01*e3rc11-2*e1rc11*e3rc00*e3rc11+e1rc00*e3rc10_2+2*e1rc12*e3rc02*e3rc10+2*e1rc11*e3rc01*e3rc10+2*e1rc10*e3rc00*e3rc10+e1rc00*e3rc02_2+2*e1rc02*e3rc00*e3rc02+e1rc00*e3rc01_2+2*e1rc01*e3rc00*e3rc01+3*e1rc00*e3rc00_2;
    AMatrix(0, 18) = -e2rc00*e3rc22_2+2*e2rc02*e3rc20*e3rc22+2*e2rc20*e3rc02*e3rc22-2*e2rc22*e3rc00*e3rc22-e2rc00*e3rc21_2+2*e2rc01*e3rc20*e3rc21+2*e2rc20*e3rc01*e3rc21-2*e2rc21*e3rc00*e3rc21+e2rc00*e3rc20_2+2*e2rc22*e3rc02*e3rc20+2*e2rc21*e3rc01*e3rc20+2*e2rc20*e3rc00*e3rc20-e2rc00*e3rc12_2+2*e2rc02*e3rc10*e3rc12+2*e2rc10*e3rc02*e3rc12-2*e2rc12*e3rc00*e3rc12-e2rc00*e3rc11_2+2*e2rc01*e3rc10*e3rc11+2*e2rc10*e3rc01*e3rc11-2*e2rc11*e3rc00*e3rc11+e2rc00*e3rc10_2+2*e2rc12*e3rc02*e3rc10+2*e2rc11*e3rc01*e3rc10+2*e2rc10*e3rc00*e3rc10+e2rc00*e3rc02_2+2*e2rc02*e3rc00*e3rc02+e2rc00*e3rc01_2+2*e2rc01*e3rc00*e3rc01+3*e2rc00*e3rc00_2;
    AMatrix(0, 19) = -e3rc00*e3rc22_2+2*e3rc02*e3rc20*e3rc22-e3rc00*e3rc21_2+2*e3rc01*e3rc20*e3rc21+e3rc00*e3rc20_2-e3rc00*e3rc12_2+2*e3rc02*e3rc10*e3rc12-e3rc00*e3rc11_2+2*e3rc01*e3rc10*e3rc11+e3rc00*e3rc10_2+e3rc00*e3rc02_2+e3rc00*e3rc01_2+e3rc00_3;
    AMatrix(1, 0) = -e0rc01*e0rc22_2+2*e0rc02*e0rc21*e0rc22+e0rc01*e0rc21_2+2*e0rc00*e0rc20*e0rc21-e0rc01*e0rc20_2-e0rc01*e0rc12_2+2*e0rc02*e0rc11*e0rc12+e0rc01*e0rc11_2+2*e0rc00*e0rc10*e0rc11-e0rc01*e0rc10_2+e0rc01*e0rc02_2+e0rc01_3+e0rc00_2*e0rc01;
    AMatrix(1, 1) = -2*e0rc01*e0rc22*e1rc22+2*e0rc02*e0rc21*e1rc22+2*e0rc02*e0rc22*e1rc21+2*e0rc01*e0rc21*e1rc21+2*e0rc00*e0rc20*e1rc21+2*e0rc00*e0rc21*e1rc20-2*e0rc01*e0rc20*e1rc20-2*e0rc01*e0rc12*e1rc12+2*e0rc02*e0rc11*e1rc12+2*e0rc02*e0rc12*e1rc11+2*e0rc01*e0rc11*e1rc11+2*e0rc00*e0rc10*e1rc11+2*e0rc00*e0rc11*e1rc10-2*e0rc01*e0rc10*e1rc10+2*e0rc21*e0rc22*e1rc02+2*e0rc11*e0rc12*e1rc02+2*e0rc01*e0rc02*e1rc02-e0rc22_2*e1rc01+e0rc21_2*e1rc01-e0rc20_2*e1rc01-e0rc12_2*e1rc01+e0rc11_2*e1rc01-e0rc10_2*e1rc01+e0rc02_2*e1rc01+3*e0rc01_2*e1rc01+e0rc00_2*e1rc01+2*e0rc20*e0rc21*e1rc00+2*e0rc10*e0rc11*e1rc00+2*e0rc00*e0rc01*e1rc00;
    AMatrix(1, 2) = -2*e0rc01*e0rc22*e2rc22+2*e0rc02*e0rc21*e2rc22+2*e0rc02*e0rc22*e2rc21+2*e0rc01*e0rc21*e2rc21+2*e0rc00*e0rc20*e2rc21+2*e0rc00*e0rc21*e2rc20-2*e0rc01*e0rc20*e2rc20-2*e0rc01*e0rc12*e2rc12+2*e0rc02*e0rc11*e2rc12+2*e0rc02*e0rc12*e2rc11+2*e0rc01*e0rc11*e2rc11+2*e0rc00*e0rc10*e2rc11+2*e0rc00*e0rc11*e2rc10-2*e0rc01*e0rc10*e2rc10+2*e0rc21*e0rc22*e2rc02+2*e0rc11*e0rc12*e2rc02+2*e0rc01*e0rc02*e2rc02-e0rc22_2*e2rc01+e0rc21_2*e2rc01-e0rc20_2*e2rc01-e0rc12_2*e2rc01+e0rc11_2*e2rc01-e0rc10_2*e2rc01+e0rc02_2*e2rc01+3*e0rc01_2*e2rc01+e0rc00_2*e2rc01+2*e0rc20*e0rc21*e2rc00+2*e0rc10*e0rc11*e2rc00+2*e0rc00*e0rc01*e2rc00;
    AMatrix(1, 3) = -e0rc01*e1rc22_2+2*e0rc02*e1rc21*e1rc22+2*e0rc21*e1rc02*e1rc22-2*e0rc22*e1rc01*e1rc22+e0rc01*e1rc21_2+2*e0rc00*e1rc20*e1rc21+2*e0rc22*e1rc02*e1rc21+2*e0rc21*e1rc01*e1rc21+2*e0rc20*e1rc00*e1rc21-e0rc01*e1rc20_2-2*e0rc20*e1rc01*e1rc20+2*e0rc21*e1rc00*e1rc20-e0rc01*e1rc12_2+2*e0rc02*e1rc11*e1rc12+2*e0rc11*e1rc02*e1rc12-2*e0rc12*e1rc01*e1rc12+e0rc01*e1rc11_2+2*e0rc00*e1rc10*e1rc11+2*e0rc12*e1rc02*e1rc11+2*e0rc11*e1rc01*e1rc11+2*e0rc10*e1rc00*e1rc11-e0rc01*e1rc10_2-2*e0rc10*e1rc01*e1rc10+2*e0rc11*e1rc00*e1rc10+e0rc01*e1rc02_2+2*e0rc02*e1rc01*e1rc02+3*e0rc01*e1rc01_2+2*e0rc00*e1rc00*e1rc01+e0rc01*e1rc00_2;
    AMatrix(1, 4) = -2*e0rc01*e1rc22*e2rc22+2*e0rc02*e1rc21*e2rc22+2*e0rc21*e1rc02*e2rc22-2*e0rc22*e1rc01*e2rc22+2*e0rc02*e1rc22*e2rc21+2*e0rc01*e1rc21*e2rc21+2*e0rc00*e1rc20*e2rc21+2*e0rc22*e1rc02*e2rc21+2*e0rc21*e1rc01*e2rc21+2*e0rc20*e1rc00*e2rc21+2*e0rc00*e1rc21*e2rc20-2*e0rc01*e1rc20*e2rc20-2*e0rc20*e1rc01*e2rc20+2*e0rc21*e1rc00*e2rc20-2*e0rc01*e1rc12*e2rc12+2*e0rc02*e1rc11*e2rc12+2*e0rc11*e1rc02*e2rc12-2*e0rc12*e1rc01*e2rc12+2*e0rc02*e1rc12*e2rc11+2*e0rc01*e1rc11*e2rc11+2*e0rc00*e1rc10*e2rc11+2*e0rc12*e1rc02*e2rc11+2*e0rc11*e1rc01*e2rc11+2*e0rc10*e1rc00*e2rc11+2*e0rc00*e1rc11*e2rc10-2*e0rc01*e1rc10*e2rc10-2*e0rc10*e1rc01*e2rc10+2*e0rc11*e1rc00*e2rc10+2*e0rc21*e1rc22*e2rc02+2*e0rc22*e1rc21*e2rc02+2*e0rc11*e1rc12*e2rc02+2*e0rc12*e1rc11*e2rc02+2*e0rc01*e1rc02*e2rc02+2*e0rc02*e1rc01*e2rc02-2*e0rc22*e1rc22*e2rc01+2*e0rc21*e1rc21*e2rc01-2*e0rc20*e1rc20*e2rc01-2*e0rc12*e1rc12*e2rc01+2*e0rc11*e1rc11*e2rc01-2*e0rc10*e1rc10*e2rc01+2*e0rc02*e1rc02*e2rc01+6*e0rc01*e1rc01*e2rc01+2*e0rc00*e1rc00*e2rc01+2*e0rc20*e1rc21*e2rc00+2*e0rc21*e1rc20*e2rc00+2*e0rc10*e1rc11*e2rc00+2*e0rc11*e1rc10*e2rc00+2*e0rc00*e1rc01*e2rc00+2*e0rc01*e1rc00*e2rc00;
    AMatrix(1, 5) = -e0rc01*e2rc22_2+2*e0rc02*e2rc21*e2rc22+2*e0rc21*e2rc02*e2rc22-2*e0rc22*e2rc01*e2rc22+e0rc01*e2rc21_2+2*e0rc00*e2rc20*e2rc21+2*e0rc22*e2rc02*e2rc21+2*e0rc21*e2rc01*e2rc21+2*e0rc20*e2rc00*e2rc21-e0rc01*e2rc20_2-2*e0rc20*e2rc01*e2rc20+2*e0rc21*e2rc00*e2rc20-e0rc01*e2rc12_2+2*e0rc02*e2rc11*e2rc12+2*e0rc11*e2rc02*e2rc12-2*e0rc12*e2rc01*e2rc12+e0rc01*e2rc11_2+2*e0rc00*e2rc10*e2rc11+2*e0rc12*e2rc02*e2rc11+2*e0rc11*e2rc01*e2rc11+2*e0rc10*e2rc00*e2rc11-e0rc01*e2rc10_2-2*e0rc10*e2rc01*e2rc10+2*e0rc11*e2rc00*e2rc10+e0rc01*e2rc02_2+2*e0rc02*e2rc01*e2rc02+3*e0rc01*e2rc01_2+2*e0rc00*e2rc00*e2rc01+e0rc01*e2rc00_2;
    AMatrix(1, 6) = -e1rc01*e1rc22_2+2*e1rc02*e1rc21*e1rc22+e1rc01*e1rc21_2+2*e1rc00*e1rc20*e1rc21-e1rc01*e1rc20_2-e1rc01*e1rc12_2+2*e1rc02*e1rc11*e1rc12+e1rc01*e1rc11_2+2*e1rc00*e1rc10*e1rc11-e1rc01*e1rc10_2+e1rc01*e1rc02_2+e1rc01_3+e1rc00_2*e1rc01;
    AMatrix(1, 7) = -2*e1rc01*e1rc22*e2rc22+2*e1rc02*e1rc21*e2rc22+2*e1rc02*e1rc22*e2rc21+2*e1rc01*e1rc21*e2rc21+2*e1rc00*e1rc20*e2rc21+2*e1rc00*e1rc21*e2rc20-2*e1rc01*e1rc20*e2rc20-2*e1rc01*e1rc12*e2rc12+2*e1rc02*e1rc11*e2rc12+2*e1rc02*e1rc12*e2rc11+2*e1rc01*e1rc11*e2rc11+2*e1rc00*e1rc10*e2rc11+2*e1rc00*e1rc11*e2rc10-2*e1rc01*e1rc10*e2rc10+2*e1rc21*e1rc22*e2rc02+2*e1rc11*e1rc12*e2rc02+2*e1rc01*e1rc02*e2rc02-e1rc22_2*e2rc01+e1rc21_2*e2rc01-e1rc20_2*e2rc01-e1rc12_2*e2rc01+e1rc11_2*e2rc01-e1rc10_2*e2rc01+e1rc02_2*e2rc01+3*e1rc01_2*e2rc01+e1rc00_2*e2rc01+2*e1rc20*e1rc21*e2rc00+2*e1rc10*e1rc11*e2rc00+2*e1rc00*e1rc01*e2rc00;
    AMatrix(1, 8) = -e1rc01*e2rc22_2+2*e1rc02*e2rc21*e2rc22+2*e1rc21*e2rc02*e2rc22-2*e1rc22*e2rc01*e2rc22+e1rc01*e2rc21_2+2*e1rc00*e2rc20*e2rc21+2*e1rc22*e2rc02*e2rc21+2*e1rc21*e2rc01*e2rc21+2*e1rc20*e2rc00*e2rc21-e1rc01*e2rc20_2-2*e1rc20*e2rc01*e2rc20+2*e1rc21*e2rc00*e2rc20-e1rc01*e2rc12_2+2*e1rc02*e2rc11*e2rc12+2*e1rc11*e2rc02*e2rc12-2*e1rc12*e2rc01*e2rc12+e1rc01*e2rc11_2+2*e1rc00*e2rc10*e2rc11+2*e1rc12*e2rc02*e2rc11+2*e1rc11*e2rc01*e2rc11+2*e1rc10*e2rc00*e2rc11-e1rc01*e2rc10_2-2*e1rc10*e2rc01*e2rc10+2*e1rc11*e2rc00*e2rc10+e1rc01*e2rc02_2+2*e1rc02*e2rc01*e2rc02+3*e1rc01*e2rc01_2+2*e1rc00*e2rc00*e2rc01+e1rc01*e2rc00_2;
    AMatrix(1, 9) = -e2rc01*e2rc22_2+2*e2rc02*e2rc21*e2rc22+e2rc01*e2rc21_2+2*e2rc00*e2rc20*e2rc21-e2rc01*e2rc20_2-e2rc01*e2rc12_2+2*e2rc02*e2rc11*e2rc12+e2rc01*e2rc11_2+2*e2rc00*e2rc10*e2rc11-e2rc01*e2rc10_2+e2rc01*e2rc02_2+e2rc01_3+e2rc00_2*e2rc01;
    AMatrix(1, 10) = -2*e0rc01*e0rc22*e3rc22+2*e0rc02*e0rc21*e3rc22+2*e0rc02*e0rc22*e3rc21+2*e0rc01*e0rc21*e3rc21+2*e0rc00*e0rc20*e3rc21+2*e0rc00*e0rc21*e3rc20-2*e0rc01*e0rc20*e3rc20-2*e0rc01*e0rc12*e3rc12+2*e0rc02*e0rc11*e3rc12+2*e0rc02*e0rc12*e3rc11+2*e0rc01*e0rc11*e3rc11+2*e0rc00*e0rc10*e3rc11+2*e0rc00*e0rc11*e3rc10-2*e0rc01*e0rc10*e3rc10+2*e0rc21*e0rc22*e3rc02+2*e0rc11*e0rc12*e3rc02+2*e0rc01*e0rc02*e3rc02-e0rc22_2*e3rc01+e0rc21_2*e3rc01-e0rc20_2*e3rc01-e0rc12_2*e3rc01+e0rc11_2*e3rc01-e0rc10_2*e3rc01+e0rc02_2*e3rc01+3*e0rc01_2*e3rc01+e0rc00_2*e3rc01+2*e0rc20*e0rc21*e3rc00+2*e0rc10*e0rc11*e3rc00+2*e0rc00*e0rc01*e3rc00;
    AMatrix(1, 11) = -2*e0rc01*e1rc22*e3rc22+2*e0rc02*e1rc21*e3rc22+2*e0rc21*e1rc02*e3rc22-2*e0rc22*e1rc01*e3rc22+2*e0rc02*e1rc22*e3rc21+2*e0rc01*e1rc21*e3rc21+2*e0rc00*e1rc20*e3rc21+2*e0rc22*e1rc02*e3rc21+2*e0rc21*e1rc01*e3rc21+2*e0rc20*e1rc00*e3rc21+2*e0rc00*e1rc21*e3rc20-2*e0rc01*e1rc20*e3rc20-2*e0rc20*e1rc01*e3rc20+2*e0rc21*e1rc00*e3rc20-2*e0rc01*e1rc12*e3rc12+2*e0rc02*e1rc11*e3rc12+2*e0rc11*e1rc02*e3rc12-2*e0rc12*e1rc01*e3rc12+2*e0rc02*e1rc12*e3rc11+2*e0rc01*e1rc11*e3rc11+2*e0rc00*e1rc10*e3rc11+2*e0rc12*e1rc02*e3rc11+2*e0rc11*e1rc01*e3rc11+2*e0rc10*e1rc00*e3rc11+2*e0rc00*e1rc11*e3rc10-2*e0rc01*e1rc10*e3rc10-2*e0rc10*e1rc01*e3rc10+2*e0rc11*e1rc00*e3rc10+2*e0rc21*e1rc22*e3rc02+2*e0rc22*e1rc21*e3rc02+2*e0rc11*e1rc12*e3rc02+2*e0rc12*e1rc11*e3rc02+2*e0rc01*e1rc02*e3rc02+2*e0rc02*e1rc01*e3rc02-2*e0rc22*e1rc22*e3rc01+2*e0rc21*e1rc21*e3rc01-2*e0rc20*e1rc20*e3rc01-2*e0rc12*e1rc12*e3rc01+2*e0rc11*e1rc11*e3rc01-2*e0rc10*e1rc10*e3rc01+2*e0rc02*e1rc02*e3rc01+6*e0rc01*e1rc01*e3rc01+2*e0rc00*e1rc00*e3rc01+2*e0rc20*e1rc21*e3rc00+2*e0rc21*e1rc20*e3rc00+2*e0rc10*e1rc11*e3rc00+2*e0rc11*e1rc10*e3rc00+2*e0rc00*e1rc01*e3rc00+2*e0rc01*e1rc00*e3rc00;
    AMatrix(1, 12) = -2*e0rc01*e2rc22*e3rc22+2*e0rc02*e2rc21*e3rc22+2*e0rc21*e2rc02*e3rc22-2*e0rc22*e2rc01*e3rc22+2*e0rc02*e2rc22*e3rc21+2*e0rc01*e2rc21*e3rc21+2*e0rc00*e2rc20*e3rc21+2*e0rc22*e2rc02*e3rc21+2*e0rc21*e2rc01*e3rc21+2*e0rc20*e2rc00*e3rc21+2*e0rc00*e2rc21*e3rc20-2*e0rc01*e2rc20*e3rc20-2*e0rc20*e2rc01*e3rc20+2*e0rc21*e2rc00*e3rc20-2*e0rc01*e2rc12*e3rc12+2*e0rc02*e2rc11*e3rc12+2*e0rc11*e2rc02*e3rc12-2*e0rc12*e2rc01*e3rc12+2*e0rc02*e2rc12*e3rc11+2*e0rc01*e2rc11*e3rc11+2*e0rc00*e2rc10*e3rc11+2*e0rc12*e2rc02*e3rc11+2*e0rc11*e2rc01*e3rc11+2*e0rc10*e2rc00*e3rc11+2*e0rc00*e2rc11*e3rc10-2*e0rc01*e2rc10*e3rc10-2*e0rc10*e2rc01*e3rc10+2*e0rc11*e2rc00*e3rc10+2*e0rc21*e2rc22*e3rc02+2*e0rc22*e2rc21*e3rc02+2*e0rc11*e2rc12*e3rc02+2*e0rc12*e2rc11*e3rc02+2*e0rc01*e2rc02*e3rc02+2*e0rc02*e2rc01*e3rc02-2*e0rc22*e2rc22*e3rc01+2*e0rc21*e2rc21*e3rc01-2*e0rc20*e2rc20*e3rc01-2*e0rc12*e2rc12*e3rc01+2*e0rc11*e2rc11*e3rc01-2*e0rc10*e2rc10*e3rc01+2*e0rc02*e2rc02*e3rc01+6*e0rc01*e2rc01*e3rc01+2*e0rc00*e2rc00*e3rc01+2*e0rc20*e2rc21*e3rc00+2*e0rc21*e2rc20*e3rc00+2*e0rc10*e2rc11*e3rc00+2*e0rc11*e2rc10*e3rc00+2*e0rc00*e2rc01*e3rc00+2*e0rc01*e2rc00*e3rc00;
    AMatrix(1, 13) = -2*e1rc01*e1rc22*e3rc22+2*e1rc02*e1rc21*e3rc22+2*e1rc02*e1rc22*e3rc21+2*e1rc01*e1rc21*e3rc21+2*e1rc00*e1rc20*e3rc21+2*e1rc00*e1rc21*e3rc20-2*e1rc01*e1rc20*e3rc20-2*e1rc01*e1rc12*e3rc12+2*e1rc02*e1rc11*e3rc12+2*e1rc02*e1rc12*e3rc11+2*e1rc01*e1rc11*e3rc11+2*e1rc00*e1rc10*e3rc11+2*e1rc00*e1rc11*e3rc10-2*e1rc01*e1rc10*e3rc10+2*e1rc21*e1rc22*e3rc02+2*e1rc11*e1rc12*e3rc02+2*e1rc01*e1rc02*e3rc02-e1rc22_2*e3rc01+e1rc21_2*e3rc01-e1rc20_2*e3rc01-e1rc12_2*e3rc01+e1rc11_2*e3rc01-e1rc10_2*e3rc01+e1rc02_2*e3rc01+3*e1rc01_2*e3rc01+e1rc00_2*e3rc01+2*e1rc20*e1rc21*e3rc00+2*e1rc10*e1rc11*e3rc00+2*e1rc00*e1rc01*e3rc00;
    AMatrix(1, 14) = -2*e1rc01*e2rc22*e3rc22+2*e1rc02*e2rc21*e3rc22+2*e1rc21*e2rc02*e3rc22-2*e1rc22*e2rc01*e3rc22+2*e1rc02*e2rc22*e3rc21+2*e1rc01*e2rc21*e3rc21+2*e1rc00*e2rc20*e3rc21+2*e1rc22*e2rc02*e3rc21+2*e1rc21*e2rc01*e3rc21+2*e1rc20*e2rc00*e3rc21+2*e1rc00*e2rc21*e3rc20-2*e1rc01*e2rc20*e3rc20-2*e1rc20*e2rc01*e3rc20+2*e1rc21*e2rc00*e3rc20-2*e1rc01*e2rc12*e3rc12+2*e1rc02*e2rc11*e3rc12+2*e1rc11*e2rc02*e3rc12-2*e1rc12*e2rc01*e3rc12+2*e1rc02*e2rc12*e3rc11+2*e1rc01*e2rc11*e3rc11+2*e1rc00*e2rc10*e3rc11+2*e1rc12*e2rc02*e3rc11+2*e1rc11*e2rc01*e3rc11+2*e1rc10*e2rc00*e3rc11+2*e1rc00*e2rc11*e3rc10-2*e1rc01*e2rc10*e3rc10-2*e1rc10*e2rc01*e3rc10+2*e1rc11*e2rc00*e3rc10+2*e1rc21*e2rc22*e3rc02+2*e1rc22*e2rc21*e3rc02+2*e1rc11*e2rc12*e3rc02+2*e1rc12*e2rc11*e3rc02+2*e1rc01*e2rc02*e3rc02+2*e1rc02*e2rc01*e3rc02-2*e1rc22*e2rc22*e3rc01+2*e1rc21*e2rc21*e3rc01-2*e1rc20*e2rc20*e3rc01-2*e1rc12*e2rc12*e3rc01+2*e1rc11*e2rc11*e3rc01-2*e1rc10*e2rc10*e3rc01+2*e1rc02*e2rc02*e3rc01+6*e1rc01*e2rc01*e3rc01+2*e1rc00*e2rc00*e3rc01+2*e1rc20*e2rc21*e3rc00+2*e1rc21*e2rc20*e3rc00+2*e1rc10*e2rc11*e3rc00+2*e1rc11*e2rc10*e3rc00+2*e1rc00*e2rc01*e3rc00+2*e1rc01*e2rc00*e3rc00;
    AMatrix(1, 15) = -2*e2rc01*e2rc22*e3rc22+2*e2rc02*e2rc21*e3rc22+2*e2rc02*e2rc22*e3rc21+2*e2rc01*e2rc21*e3rc21+2*e2rc00*e2rc20*e3rc21+2*e2rc00*e2rc21*e3rc20-2*e2rc01*e2rc20*e3rc20-2*e2rc01*e2rc12*e3rc12+2*e2rc02*e2rc11*e3rc12+2*e2rc02*e2rc12*e3rc11+2*e2rc01*e2rc11*e3rc11+2*e2rc00*e2rc10*e3rc11+2*e2rc00*e2rc11*e3rc10-2*e2rc01*e2rc10*e3rc10+2*e2rc21*e2rc22*e3rc02+2*e2rc11*e2rc12*e3rc02+2*e2rc01*e2rc02*e3rc02-e2rc22_2*e3rc01+e2rc21_2*e3rc01-e2rc20_2*e3rc01-e2rc12_2*e3rc01+e2rc11_2*e3rc01-e2rc10_2*e3rc01+e2rc02_2*e3rc01+3*e2rc01_2*e3rc01+e2rc00_2*e3rc01+2*e2rc20*e2rc21*e3rc00+2*e2rc10*e2rc11*e3rc00+2*e2rc00*e2rc01*e3rc00;
    AMatrix(1, 16) = -e0rc01*e3rc22_2+2*e0rc02*e3rc21*e3rc22+2*e0rc21*e3rc02*e3rc22-2*e0rc22*e3rc01*e3rc22+e0rc01*e3rc21_2+2*e0rc00*e3rc20*e3rc21+2*e0rc22*e3rc02*e3rc21+2*e0rc21*e3rc01*e3rc21+2*e0rc20*e3rc00*e3rc21-e0rc01*e3rc20_2-2*e0rc20*e3rc01*e3rc20+2*e0rc21*e3rc00*e3rc20-e0rc01*e3rc12_2+2*e0rc02*e3rc11*e3rc12+2*e0rc11*e3rc02*e3rc12-2*e0rc12*e3rc01*e3rc12+e0rc01*e3rc11_2+2*e0rc00*e3rc10*e3rc11+2*e0rc12*e3rc02*e3rc11+2*e0rc11*e3rc01*e3rc11+2*e0rc10*e3rc00*e3rc11-e0rc01*e3rc10_2-2*e0rc10*e3rc01*e3rc10+2*e0rc11*e3rc00*e3rc10+e0rc01*e3rc02_2+2*e0rc02*e3rc01*e3rc02+3*e0rc01*e3rc01_2+2*e0rc00*e3rc00*e3rc01+e0rc01*e3rc00_2;
    AMatrix(1, 17) = -e1rc01*e3rc22_2+2*e1rc02*e3rc21*e3rc22+2*e1rc21*e3rc02*e3rc22-2*e1rc22*e3rc01*e3rc22+e1rc01*e3rc21_2+2*e1rc00*e3rc20*e3rc21+2*e1rc22*e3rc02*e3rc21+2*e1rc21*e3rc01*e3rc21+2*e1rc20*e3rc00*e3rc21-e1rc01*e3rc20_2-2*e1rc20*e3rc01*e3rc20+2*e1rc21*e3rc00*e3rc20-e1rc01*e3rc12_2+2*e1rc02*e3rc11*e3rc12+2*e1rc11*e3rc02*e3rc12-2*e1rc12*e3rc01*e3rc12+e1rc01*e3rc11_2+2*e1rc00*e3rc10*e3rc11+2*e1rc12*e3rc02*e3rc11+2*e1rc11*e3rc01*e3rc11+2*e1rc10*e3rc00*e3rc11-e1rc01*e3rc10_2-2*e1rc10*e3rc01*e3rc10+2*e1rc11*e3rc00*e3rc10+e1rc01*e3rc02_2+2*e1rc02*e3rc01*e3rc02+3*e1rc01*e3rc01_2+2*e1rc00*e3rc00*e3rc01+e1rc01*e3rc00_2;
    AMatrix(1, 18) = -e2rc01*e3rc22_2+2*e2rc02*e3rc21*e3rc22+2*e2rc21*e3rc02*e3rc22-2*e2rc22*e3rc01*e3rc22+e2rc01*e3rc21_2+2*e2rc00*e3rc20*e3rc21+2*e2rc22*e3rc02*e3rc21+2*e2rc21*e3rc01*e3rc21+2*e2rc20*e3rc00*e3rc21-e2rc01*e3rc20_2-2*e2rc20*e3rc01*e3rc20+2*e2rc21*e3rc00*e3rc20-e2rc01*e3rc12_2+2*e2rc02*e3rc11*e3rc12+2*e2rc11*e3rc02*e3rc12-2*e2rc12*e3rc01*e3rc12+e2rc01*e3rc11_2+2*e2rc00*e3rc10*e3rc11+2*e2rc12*e3rc02*e3rc11+2*e2rc11*e3rc01*e3rc11+2*e2rc10*e3rc00*e3rc11-e2rc01*e3rc10_2-2*e2rc10*e3rc01*e3rc10+2*e2rc11*e3rc00*e3rc10+e2rc01*e3rc02_2+2*e2rc02*e3rc01*e3rc02+3*e2rc01*e3rc01_2+2*e2rc00*e3rc00*e3rc01+e2rc01*e3rc00_2;
    AMatrix(1, 19) = -e3rc01*e3rc22_2+2*e3rc02*e3rc21*e3rc22+e3rc01*e3rc21_2+2*e3rc00*e3rc20*e3rc21-e3rc01*e3rc20_2-e3rc01*e3rc12_2+2*e3rc02*e3rc11*e3rc12+e3rc01*e3rc11_2+2*e3rc00*e3rc10*e3rc11-e3rc01*e3rc10_2+e3rc01*e3rc02_2+e3rc01_3+e3rc00_2*e3rc01;
    AMatrix(2, 0) = +e0rc02*e0rc22_2+2*e0rc01*e0rc21*e0rc22+2*e0rc00*e0rc20*e0rc22-e0rc02*e0rc21_2-e0rc02*e0rc20_2+e0rc02*e0rc12_2+2*e0rc01*e0rc11*e0rc12+2*e0rc00*e0rc10*e0rc12-e0rc02*e0rc11_2-e0rc02*e0rc10_2+e0rc02_3+e0rc01_2*e0rc02+e0rc00_2*e0rc02;
    AMatrix(2, 1) = +2*e0rc02*e0rc22*e1rc22+2*e0rc01*e0rc21*e1rc22+2*e0rc00*e0rc20*e1rc22+2*e0rc01*e0rc22*e1rc21-2*e0rc02*e0rc21*e1rc21+2*e0rc00*e0rc22*e1rc20-2*e0rc02*e0rc20*e1rc20+2*e0rc02*e0rc12*e1rc12+2*e0rc01*e0rc11*e1rc12+2*e0rc00*e0rc10*e1rc12+2*e0rc01*e0rc12*e1rc11-2*e0rc02*e0rc11*e1rc11+2*e0rc00*e0rc12*e1rc10-2*e0rc02*e0rc10*e1rc10+e0rc22_2*e1rc02-e0rc21_2*e1rc02-e0rc20_2*e1rc02+e0rc12_2*e1rc02-e0rc11_2*e1rc02-e0rc10_2*e1rc02+3*e0rc02_2*e1rc02+e0rc01_2*e1rc02+e0rc00_2*e1rc02+2*e0rc21*e0rc22*e1rc01+2*e0rc11*e0rc12*e1rc01+2*e0rc01*e0rc02*e1rc01+2*e0rc20*e0rc22*e1rc00+2*e0rc10*e0rc12*e1rc00+2*e0rc00*e0rc02*e1rc00;
    AMatrix(2, 2) = +2*e0rc02*e0rc22*e2rc22+2*e0rc01*e0rc21*e2rc22+2*e0rc00*e0rc20*e2rc22+2*e0rc01*e0rc22*e2rc21-2*e0rc02*e0rc21*e2rc21+2*e0rc00*e0rc22*e2rc20-2*e0rc02*e0rc20*e2rc20+2*e0rc02*e0rc12*e2rc12+2*e0rc01*e0rc11*e2rc12+2*e0rc00*e0rc10*e2rc12+2*e0rc01*e0rc12*e2rc11-2*e0rc02*e0rc11*e2rc11+2*e0rc00*e0rc12*e2rc10-2*e0rc02*e0rc10*e2rc10+e0rc22_2*e2rc02-e0rc21_2*e2rc02-e0rc20_2*e2rc02+e0rc12_2*e2rc02-e0rc11_2*e2rc02-e0rc10_2*e2rc02+3*e0rc02_2*e2rc02+e0rc01_2*e2rc02+e0rc00_2*e2rc02+2*e0rc21*e0rc22*e2rc01+2*e0rc11*e0rc12*e2rc01+2*e0rc01*e0rc02*e2rc01+2*e0rc20*e0rc22*e2rc00+2*e0rc10*e0rc12*e2rc00+2*e0rc00*e0rc02*e2rc00;
    AMatrix(2, 3) = +e0rc02*e1rc22_2+2*e0rc01*e1rc21*e1rc22+2*e0rc00*e1rc20*e1rc22+2*e0rc22*e1rc02*e1rc22+2*e0rc21*e1rc01*e1rc22+2*e0rc20*e1rc00*e1rc22-e0rc02*e1rc21_2-2*e0rc21*e1rc02*e1rc21+2*e0rc22*e1rc01*e1rc21-e0rc02*e1rc20_2-2*e0rc20*e1rc02*e1rc20+2*e0rc22*e1rc00*e1rc20+e0rc02*e1rc12_2+2*e0rc01*e1rc11*e1rc12+2*e0rc00*e1rc10*e1rc12+2*e0rc12*e1rc02*e1rc12+2*e0rc11*e1rc01*e1rc12+2*e0rc10*e1rc00*e1rc12-e0rc02*e1rc11_2-2*e0rc11*e1rc02*e1rc11+2*e0rc12*e1rc01*e1rc11-e0rc02*e1rc10_2-2*e0rc10*e1rc02*e1rc10+2*e0rc12*e1rc00*e1rc10+3*e0rc02*e1rc02_2+2*e0rc01*e1rc01*e1rc02+2*e0rc00*e1rc00*e1rc02+e0rc02*e1rc01_2+e0rc02*e1rc00_2;
    AMatrix(2, 4) = +2*e0rc02*e1rc22*e2rc22+2*e0rc01*e1rc21*e2rc22+2*e0rc00*e1rc20*e2rc22+2*e0rc22*e1rc02*e2rc22+2*e0rc21*e1rc01*e2rc22+2*e0rc20*e1rc00*e2rc22+2*e0rc01*e1rc22*e2rc21-2*e0rc02*e1rc21*e2rc21-2*e0rc21*e1rc02*e2rc21+2*e0rc22*e1rc01*e2rc21+2*e0rc00*e1rc22*e2rc20-2*e0rc02*e1rc20*e2rc20-2*e0rc20*e1rc02*e2rc20+2*e0rc22*e1rc00*e2rc20+2*e0rc02*e1rc12*e2rc12+2*e0rc01*e1rc11*e2rc12+2*e0rc00*e1rc10*e2rc12+2*e0rc12*e1rc02*e2rc12+2*e0rc11*e1rc01*e2rc12+2*e0rc10*e1rc00*e2rc12+2*e0rc01*e1rc12*e2rc11-2*e0rc02*e1rc11*e2rc11-2*e0rc11*e1rc02*e2rc11+2*e0rc12*e1rc01*e2rc11+2*e0rc00*e1rc12*e2rc10-2*e0rc02*e1rc10*e2rc10-2*e0rc10*e1rc02*e2rc10+2*e0rc12*e1rc00*e2rc10+2*e0rc22*e1rc22*e2rc02-2*e0rc21*e1rc21*e2rc02-2*e0rc20*e1rc20*e2rc02+2*e0rc12*e1rc12*e2rc02-2*e0rc11*e1rc11*e2rc02-2*e0rc10*e1rc10*e2rc02+6*e0rc02*e1rc02*e2rc02+2*e0rc01*e1rc01*e2rc02+2*e0rc00*e1rc00*e2rc02+2*e0rc21*e1rc22*e2rc01+2*e0rc22*e1rc21*e2rc01+2*e0rc11*e1rc12*e2rc01+2*e0rc12*e1rc11*e2rc01+2*e0rc01*e1rc02*e2rc01+2*e0rc02*e1rc01*e2rc01+2*e0rc20*e1rc22*e2rc00+2*e0rc22*e1rc20*e2rc00+2*e0rc10*e1rc12*e2rc00+2*e0rc12*e1rc10*e2rc00+2*e0rc00*e1rc02*e2rc00+2*e0rc02*e1rc00*e2rc00;
    AMatrix(2, 5) = +e0rc02*e2rc22_2+2*e0rc01*e2rc21*e2rc22+2*e0rc00*e2rc20*e2rc22+2*e0rc22*e2rc02*e2rc22+2*e0rc21*e2rc01*e2rc22+2*e0rc20*e2rc00*e2rc22-e0rc02*e2rc21_2-2*e0rc21*e2rc02*e2rc21+2*e0rc22*e2rc01*e2rc21-e0rc02*e2rc20_2-2*e0rc20*e2rc02*e2rc20+2*e0rc22*e2rc00*e2rc20+e0rc02*e2rc12_2+2*e0rc01*e2rc11*e2rc12+2*e0rc00*e2rc10*e2rc12+2*e0rc12*e2rc02*e2rc12+2*e0rc11*e2rc01*e2rc12+2*e0rc10*e2rc00*e2rc12-e0rc02*e2rc11_2-2*e0rc11*e2rc02*e2rc11+2*e0rc12*e2rc01*e2rc11-e0rc02*e2rc10_2-2*e0rc10*e2rc02*e2rc10+2*e0rc12*e2rc00*e2rc10+3*e0rc02*e2rc02_2+2*e0rc01*e2rc01*e2rc02+2*e0rc00*e2rc00*e2rc02+e0rc02*e2rc01_2+e0rc02*e2rc00_2;
    AMatrix(2, 6) = +e1rc02*e1rc22_2+2*e1rc01*e1rc21*e1rc22+2*e1rc00*e1rc20*e1rc22-e1rc02*e1rc21_2-e1rc02*e1rc20_2+e1rc02*e1rc12_2+2*e1rc01*e1rc11*e1rc12+2*e1rc00*e1rc10*e1rc12-e1rc02*e1rc11_2-e1rc02*e1rc10_2+e1rc02_3+e1rc01_2*e1rc02+e1rc00_2*e1rc02;
    AMatrix(2, 7) = +2*e1rc02*e1rc22*e2rc22+2*e1rc01*e1rc21*e2rc22+2*e1rc00*e1rc20*e2rc22+2*e1rc01*e1rc22*e2rc21-2*e1rc02*e1rc21*e2rc21+2*e1rc00*e1rc22*e2rc20-2*e1rc02*e1rc20*e2rc20+2*e1rc02*e1rc12*e2rc12+2*e1rc01*e1rc11*e2rc12+2*e1rc00*e1rc10*e2rc12+2*e1rc01*e1rc12*e2rc11-2*e1rc02*e1rc11*e2rc11+2*e1rc00*e1rc12*e2rc10-2*e1rc02*e1rc10*e2rc10+e1rc22_2*e2rc02-e1rc21_2*e2rc02-e1rc20_2*e2rc02+e1rc12_2*e2rc02-e1rc11_2*e2rc02-e1rc10_2*e2rc02+3*e1rc02_2*e2rc02+e1rc01_2*e2rc02+e1rc00_2*e2rc02+2*e1rc21*e1rc22*e2rc01+2*e1rc11*e1rc12*e2rc01+2*e1rc01*e1rc02*e2rc01+2*e1rc20*e1rc22*e2rc00+2*e1rc10*e1rc12*e2rc00+2*e1rc00*e1rc02*e2rc00;
    AMatrix(2, 8) = +e1rc02*e2rc22_2+2*e1rc01*e2rc21*e2rc22+2*e1rc00*e2rc20*e2rc22+2*e1rc22*e2rc02*e2rc22+2*e1rc21*e2rc01*e2rc22+2*e1rc20*e2rc00*e2rc22-e1rc02*e2rc21_2-2*e1rc21*e2rc02*e2rc21+2*e1rc22*e2rc01*e2rc21-e1rc02*e2rc20_2-2*e1rc20*e2rc02*e2rc20+2*e1rc22*e2rc00*e2rc20+e1rc02*e2rc12_2+2*e1rc01*e2rc11*e2rc12+2*e1rc00*e2rc10*e2rc12+2*e1rc12*e2rc02*e2rc12+2*e1rc11*e2rc01*e2rc12+2*e1rc10*e2rc00*e2rc12-e1rc02*e2rc11_2-2*e1rc11*e2rc02*e2rc11+2*e1rc12*e2rc01*e2rc11-e1rc02*e2rc10_2-2*e1rc10*e2rc02*e2rc10+2*e1rc12*e2rc00*e2rc10+3*e1rc02*e2rc02_2+2*e1rc01*e2rc01*e2rc02+2*e1rc00*e2rc00*e2rc02+e1rc02*e2rc01_2+e1rc02*e2rc00_2;
    AMatrix(2, 9) = e2rc02*e2rc22_2+2*e2rc01*e2rc21*e2rc22+2*e2rc00*e2rc20*e2rc22-e2rc02*e2rc21_2-e2rc02*e2rc20_2+e2rc02*e2rc12_2+2*e2rc01*e2rc11*e2rc12+2*e2rc00*e2rc10*e2rc12-e2rc02*e2rc11_2-e2rc02*e2rc10_2+e2rc02_3+e2rc01_2*e2rc02+e2rc00_2*e2rc02;
    AMatrix(2, 10) = +2*e0rc02*e0rc22*e3rc22+2*e0rc01*e0rc21*e3rc22+2*e0rc00*e0rc20*e3rc22+2*e0rc01*e0rc22*e3rc21-2*e0rc02*e0rc21*e3rc21+2*e0rc00*e0rc22*e3rc20-2*e0rc02*e0rc20*e3rc20+2*e0rc02*e0rc12*e3rc12+2*e0rc01*e0rc11*e3rc12+2*e0rc00*e0rc10*e3rc12+2*e0rc01*e0rc12*e3rc11-2*e0rc02*e0rc11*e3rc11+2*e0rc00*e0rc12*e3rc10-2*e0rc02*e0rc10*e3rc10+e0rc22_2*e3rc02-e0rc21_2*e3rc02-e0rc20_2*e3rc02+e0rc12_2*e3rc02-e0rc11_2*e3rc02-e0rc10_2*e3rc02+3*e0rc02_2*e3rc02+e0rc01_2*e3rc02+e0rc00_2*e3rc02+2*e0rc21*e0rc22*e3rc01+2*e0rc11*e0rc12*e3rc01+2*e0rc01*e0rc02*e3rc01+2*e0rc20*e0rc22*e3rc00+2*e0rc10*e0rc12*e3rc00+2*e0rc00*e0rc02*e3rc00;
    AMatrix(2, 11) = +2*e0rc02*e1rc22*e3rc22+2*e0rc01*e1rc21*e3rc22+2*e0rc00*e1rc20*e3rc22+2*e0rc22*e1rc02*e3rc22+2*e0rc21*e1rc01*e3rc22+2*e0rc20*e1rc00*e3rc22+2*e0rc01*e1rc22*e3rc21-2*e0rc02*e1rc21*e3rc21-2*e0rc21*e1rc02*e3rc21+2*e0rc22*e1rc01*e3rc21+2*e0rc00*e1rc22*e3rc20-2*e0rc02*e1rc20*e3rc20-2*e0rc20*e1rc02*e3rc20+2*e0rc22*e1rc00*e3rc20+2*e0rc02*e1rc12*e3rc12+2*e0rc01*e1rc11*e3rc12+2*e0rc00*e1rc10*e3rc12+2*e0rc12*e1rc02*e3rc12+2*e0rc11*e1rc01*e3rc12+2*e0rc10*e1rc00*e3rc12+2*e0rc01*e1rc12*e3rc11-2*e0rc02*e1rc11*e3rc11-2*e0rc11*e1rc02*e3rc11+2*e0rc12*e1rc01*e3rc11+2*e0rc00*e1rc12*e3rc10-2*e0rc02*e1rc10*e3rc10-2*e0rc10*e1rc02*e3rc10+2*e0rc12*e1rc00*e3rc10+2*e0rc22*e1rc22*e3rc02-2*e0rc21*e1rc21*e3rc02-2*e0rc20*e1rc20*e3rc02+2*e0rc12*e1rc12*e3rc02-2*e0rc11*e1rc11*e3rc02-2*e0rc10*e1rc10*e3rc02+6*e0rc02*e1rc02*e3rc02+2*e0rc01*e1rc01*e3rc02+2*e0rc00*e1rc00*e3rc02+2*e0rc21*e1rc22*e3rc01+2*e0rc22*e1rc21*e3rc01+2*e0rc11*e1rc12*e3rc01+2*e0rc12*e1rc11*e3rc01+2*e0rc01*e1rc02*e3rc01+2*e0rc02*e1rc01*e3rc01+2*e0rc20*e1rc22*e3rc00+2*e0rc22*e1rc20*e3rc00+2*e0rc10*e1rc12*e3rc00+2*e0rc12*e1rc10*e3rc00+2*e0rc00*e1rc02*e3rc00+2*e0rc02*e1rc00*e3rc00;
    AMatrix(2, 12) = +2*e0rc02*e2rc22*e3rc22+2*e0rc01*e2rc21*e3rc22+2*e0rc00*e2rc20*e3rc22+2*e0rc22*e2rc02*e3rc22+2*e0rc21*e2rc01*e3rc22+2*e0rc20*e2rc00*e3rc22+2*e0rc01*e2rc22*e3rc21-2*e0rc02*e2rc21*e3rc21-2*e0rc21*e2rc02*e3rc21+2*e0rc22*e2rc01*e3rc21+2*e0rc00*e2rc22*e3rc20-2*e0rc02*e2rc20*e3rc20-2*e0rc20*e2rc02*e3rc20+2*e0rc22*e2rc00*e3rc20+2*e0rc02*e2rc12*e3rc12+2*e0rc01*e2rc11*e3rc12+2*e0rc00*e2rc10*e3rc12+2*e0rc12*e2rc02*e3rc12+2*e0rc11*e2rc01*e3rc12+2*e0rc10*e2rc00*e3rc12+2*e0rc01*e2rc12*e3rc11-2*e0rc02*e2rc11*e3rc11-2*e0rc11*e2rc02*e3rc11+2*e0rc12*e2rc01*e3rc11+2*e0rc00*e2rc12*e3rc10-2*e0rc02*e2rc10*e3rc10-2*e0rc10*e2rc02*e3rc10+2*e0rc12*e2rc00*e3rc10+2*e0rc22*e2rc22*e3rc02-2*e0rc21*e2rc21*e3rc02-2*e0rc20*e2rc20*e3rc02+2*e0rc12*e2rc12*e3rc02-2*e0rc11*e2rc11*e3rc02-2*e0rc10*e2rc10*e3rc02+6*e0rc02*e2rc02*e3rc02+2*e0rc01*e2rc01*e3rc02+2*e0rc00*e2rc00*e3rc02+2*e0rc21*e2rc22*e3rc01+2*e0rc22*e2rc21*e3rc01+2*e0rc11*e2rc12*e3rc01+2*e0rc12*e2rc11*e3rc01+2*e0rc01*e2rc02*e3rc01+2*e0rc02*e2rc01*e3rc01+2*e0rc20*e2rc22*e3rc00+2*e0rc22*e2rc20*e3rc00+2*e0rc10*e2rc12*e3rc00+2*e0rc12*e2rc10*e3rc00+2*e0rc00*e2rc02*e3rc00+2*e0rc02*e2rc00*e3rc00;
    AMatrix(2, 13) = +2*e1rc02*e1rc22*e3rc22+2*e1rc01*e1rc21*e3rc22+2*e1rc00*e1rc20*e3rc22+2*e1rc01*e1rc22*e3rc21-2*e1rc02*e1rc21*e3rc21+2*e1rc00*e1rc22*e3rc20-2*e1rc02*e1rc20*e3rc20+2*e1rc02*e1rc12*e3rc12+2*e1rc01*e1rc11*e3rc12+2*e1rc00*e1rc10*e3rc12+2*e1rc01*e1rc12*e3rc11-2*e1rc02*e1rc11*e3rc11+2*e1rc00*e1rc12*e3rc10-2*e1rc02*e1rc10*e3rc10+e1rc22_2*e3rc02-e1rc21_2*e3rc02-e1rc20_2*e3rc02+e1rc12_2*e3rc02-e1rc11_2*e3rc02-e1rc10_2*e3rc02+3*e1rc02_2*e3rc02+e1rc01_2*e3rc02+e1rc00_2*e3rc02+2*e1rc21*e1rc22*e3rc01+2*e1rc11*e1rc12*e3rc01+2*e1rc01*e1rc02*e3rc01+2*e1rc20*e1rc22*e3rc00+2*e1rc10*e1rc12*e3rc00+2*e1rc00*e1rc02*e3rc00;
    AMatrix(2, 14) = +2*e1rc02*e2rc22*e3rc22+2*e1rc01*e2rc21*e3rc22+2*e1rc00*e2rc20*e3rc22+2*e1rc22*e2rc02*e3rc22+2*e1rc21*e2rc01*e3rc22+2*e1rc20*e2rc00*e3rc22+2*e1rc01*e2rc22*e3rc21-2*e1rc02*e2rc21*e3rc21-2*e1rc21*e2rc02*e3rc21+2*e1rc22*e2rc01*e3rc21+2*e1rc00*e2rc22*e3rc20-2*e1rc02*e2rc20*e3rc20-2*e1rc20*e2rc02*e3rc20+2*e1rc22*e2rc00*e3rc20+2*e1rc02*e2rc12*e3rc12+2*e1rc01*e2rc11*e3rc12+2*e1rc00*e2rc10*e3rc12+2*e1rc12*e2rc02*e3rc12+2*e1rc11*e2rc01*e3rc12+2*e1rc10*e2rc00*e3rc12+2*e1rc01*e2rc12*e3rc11-2*e1rc02*e2rc11*e3rc11-2*e1rc11*e2rc02*e3rc11+2*e1rc12*e2rc01*e3rc11+2*e1rc00*e2rc12*e3rc10-2*e1rc02*e2rc10*e3rc10-2*e1rc10*e2rc02*e3rc10+2*e1rc12*e2rc00*e3rc10+2*e1rc22*e2rc22*e3rc02-2*e1rc21*e2rc21*e3rc02-2*e1rc20*e2rc20*e3rc02+2*e1rc12*e2rc12*e3rc02-2*e1rc11*e2rc11*e3rc02-2*e1rc10*e2rc10*e3rc02+6*e1rc02*e2rc02*e3rc02+2*e1rc01*e2rc01*e3rc02+2*e1rc00*e2rc00*e3rc02+2*e1rc21*e2rc22*e3rc01+2*e1rc22*e2rc21*e3rc01+2*e1rc11*e2rc12*e3rc01+2*e1rc12*e2rc11*e3rc01+2*e1rc01*e2rc02*e3rc01+2*e1rc02*e2rc01*e3rc01+2*e1rc20*e2rc22*e3rc00+2*e1rc22*e2rc20*e3rc00+2*e1rc10*e2rc12*e3rc00+2*e1rc12*e2rc10*e3rc00+2*e1rc00*e2rc02*e3rc00+2*e1rc02*e2rc00*e3rc00;
    AMatrix(2, 15) = +2*e2rc02*e2rc22*e3rc22+2*e2rc01*e2rc21*e3rc22+2*e2rc00*e2rc20*e3rc22+2*e2rc01*e2rc22*e3rc21-2*e2rc02*e2rc21*e3rc21+2*e2rc00*e2rc22*e3rc20-2*e2rc02*e2rc20*e3rc20+2*e2rc02*e2rc12*e3rc12+2*e2rc01*e2rc11*e3rc12+2*e2rc00*e2rc10*e3rc12+2*e2rc01*e2rc12*e3rc11-2*e2rc02*e2rc11*e3rc11+2*e2rc00*e2rc12*e3rc10-2*e2rc02*e2rc10*e3rc10+e2rc22_2*e3rc02-e2rc21_2*e3rc02-e2rc20_2*e3rc02+e2rc12_2*e3rc02-e2rc11_2*e3rc02-e2rc10_2*e3rc02+3*e2rc02_2*e3rc02+e2rc01_2*e3rc02+e2rc00_2*e3rc02+2*e2rc21*e2rc22*e3rc01+2*e2rc11*e2rc12*e3rc01+2*e2rc01*e2rc02*e3rc01+2*e2rc20*e2rc22*e3rc00+2*e2rc10*e2rc12*e3rc00+2*e2rc00*e2rc02*e3rc00;
    AMatrix(2, 16) = +e0rc02*e3rc22_2+2*e0rc01*e3rc21*e3rc22+2*e0rc00*e3rc20*e3rc22+2*e0rc22*e3rc02*e3rc22+2*e0rc21*e3rc01*e3rc22+2*e0rc20*e3rc00*e3rc22-e0rc02*e3rc21_2-2*e0rc21*e3rc02*e3rc21+2*e0rc22*e3rc01*e3rc21-e0rc02*e3rc20_2-2*e0rc20*e3rc02*e3rc20+2*e0rc22*e3rc00*e3rc20+e0rc02*e3rc12_2+2*e0rc01*e3rc11*e3rc12+2*e0rc00*e3rc10*e3rc12+2*e0rc12*e3rc02*e3rc12+2*e0rc11*e3rc01*e3rc12+2*e0rc10*e3rc00*e3rc12-e0rc02*e3rc11_2-2*e0rc11*e3rc02*e3rc11+2*e0rc12*e3rc01*e3rc11-e0rc02*e3rc10_2-2*e0rc10*e3rc02*e3rc10+2*e0rc12*e3rc00*e3rc10+3*e0rc02*e3rc02_2+2*e0rc01*e3rc01*e3rc02+2*e0rc00*e3rc00*e3rc02+e0rc02*e3rc01_2+e0rc02*e3rc00_2;
    AMatrix(2, 17) = +e1rc02*e3rc22_2+2*e1rc01*e3rc21*e3rc22+2*e1rc00*e3rc20*e3rc22+2*e1rc22*e3rc02*e3rc22+2*e1rc21*e3rc01*e3rc22+2*e1rc20*e3rc00*e3rc22-e1rc02*e3rc21_2-2*e1rc21*e3rc02*e3rc21+2*e1rc22*e3rc01*e3rc21-e1rc02*e3rc20_2-2*e1rc20*e3rc02*e3rc20+2*e1rc22*e3rc00*e3rc20+e1rc02*e3rc12_2+2*e1rc01*e3rc11*e3rc12+2*e1rc00*e3rc10*e3rc12+2*e1rc12*e3rc02*e3rc12+2*e1rc11*e3rc01*e3rc12+2*e1rc10*e3rc00*e3rc12-e1rc02*e3rc11_2-2*e1rc11*e3rc02*e3rc11+2*e1rc12*e3rc01*e3rc11-e1rc02*e3rc10_2-2*e1rc10*e3rc02*e3rc10+2*e1rc12*e3rc00*e3rc10+3*e1rc02*e3rc02_2+2*e1rc01*e3rc01*e3rc02+2*e1rc00*e3rc00*e3rc02+e1rc02*e3rc01_2+e1rc02*e3rc00_2;
    AMatrix(2, 18) = +e2rc02*e3rc22_2+2*e2rc01*e3rc21*e3rc22+2*e2rc00*e3rc20*e3rc22+2*e2rc22*e3rc02*e3rc22+2*e2rc21*e3rc01*e3rc22+2*e2rc20*e3rc00*e3rc22-e2rc02*e3rc21_2-2*e2rc21*e3rc02*e3rc21+2*e2rc22*e3rc01*e3rc21-e2rc02*e3rc20_2-2*e2rc20*e3rc02*e3rc20+2*e2rc22*e3rc00*e3rc20+e2rc02*e3rc12_2+2*e2rc01*e3rc11*e3rc12+2*e2rc00*e3rc10*e3rc12+2*e2rc12*e3rc02*e3rc12+2*e2rc11*e3rc01*e3rc12+2*e2rc10*e3rc00*e3rc12-e2rc02*e3rc11_2-2*e2rc11*e3rc02*e3rc11+2*e2rc12*e3rc01*e3rc11-e2rc02*e3rc10_2-2*e2rc10*e3rc02*e3rc10+2*e2rc12*e3rc00*e3rc10+3*e2rc02*e3rc02_2+2*e2rc01*e3rc01*e3rc02+2*e2rc00*e3rc00*e3rc02+e2rc02*e3rc01_2+e2rc02*e3rc00_2;
    AMatrix(2, 19) = +e3rc02*e3rc22_2+2*e3rc01*e3rc21*e3rc22+2*e3rc00*e3rc20*e3rc22-e3rc02*e3rc21_2-e3rc02*e3rc20_2+e3rc02*e3rc12_2+2*e3rc01*e3rc11*e3rc12+2*e3rc00*e3rc10*e3rc12-e3rc02*e3rc11_2-e3rc02*e3rc10_2+e3rc02_3+e3rc01_2*e3rc02+e3rc00_2*e3rc02;
    AMatrix(3, 0) = -e0rc10*e0rc22_2+2*e0rc12*e0rc20*e0rc22-e0rc10*e0rc21_2+2*e0rc11*e0rc20*e0rc21+e0rc10*e0rc20_2+e0rc10*e0rc12_2+2*e0rc00*e0rc02*e0rc12+e0rc10*e0rc11_2+2*e0rc00*e0rc01*e0rc11+e0rc10_3-e0rc02_2*e0rc10-e0rc01_2*e0rc10+e0rc00_2*e0rc10;
    AMatrix(3, 1) = -2*e0rc10*e0rc22*e1rc22+2*e0rc12*e0rc20*e1rc22-2*e0rc10*e0rc21*e1rc21+2*e0rc11*e0rc20*e1rc21+2*e0rc12*e0rc22*e1rc20+2*e0rc11*e0rc21*e1rc20+2*e0rc10*e0rc20*e1rc20+2*e0rc20*e0rc22*e1rc12+2*e0rc10*e0rc12*e1rc12+2*e0rc00*e0rc02*e1rc12+2*e0rc20*e0rc21*e1rc11+2*e0rc10*e0rc11*e1rc11+2*e0rc00*e0rc01*e1rc11-e0rc22_2*e1rc10-e0rc21_2*e1rc10+e0rc20_2*e1rc10+e0rc12_2*e1rc10+e0rc11_2*e1rc10+3*e0rc10_2*e1rc10-e0rc02_2*e1rc10-e0rc01_2*e1rc10+e0rc00_2*e1rc10+2*e0rc00*e0rc12*e1rc02-2*e0rc02*e0rc10*e1rc02+2*e0rc00*e0rc11*e1rc01-2*e0rc01*e0rc10*e1rc01+2*e0rc02*e0rc12*e1rc00+2*e0rc01*e0rc11*e1rc00+2*e0rc00*e0rc10*e1rc00;
    AMatrix(3, 2) = -2*e0rc10*e0rc22*e2rc22+2*e0rc12*e0rc20*e2rc22-2*e0rc10*e0rc21*e2rc21+2*e0rc11*e0rc20*e2rc21+2*e0rc12*e0rc22*e2rc20+2*e0rc11*e0rc21*e2rc20+2*e0rc10*e0rc20*e2rc20+2*e0rc20*e0rc22*e2rc12+2*e0rc10*e0rc12*e2rc12+2*e0rc00*e0rc02*e2rc12+2*e0rc20*e0rc21*e2rc11+2*e0rc10*e0rc11*e2rc11+2*e0rc00*e0rc01*e2rc11-e0rc22_2*e2rc10-e0rc21_2*e2rc10+e0rc20_2*e2rc10+e0rc12_2*e2rc10+e0rc11_2*e2rc10+3*e0rc10_2*e2rc10-e0rc02_2*e2rc10-e0rc01_2*e2rc10+e0rc00_2*e2rc10+2*e0rc00*e0rc12*e2rc02-2*e0rc02*e0rc10*e2rc02+2*e0rc00*e0rc11*e2rc01-2*e0rc01*e0rc10*e2rc01+2*e0rc02*e0rc12*e2rc00+2*e0rc01*e0rc11*e2rc00+2*e0rc00*e0rc10*e2rc00;
    AMatrix(3, 3) = -e0rc10*e1rc22_2+2*e0rc12*e1rc20*e1rc22+2*e0rc20*e1rc12*e1rc22-2*e0rc22*e1rc10*e1rc22-e0rc10*e1rc21_2+2*e0rc11*e1rc20*e1rc21+2*e0rc20*e1rc11*e1rc21-2*e0rc21*e1rc10*e1rc21+e0rc10*e1rc20_2+2*e0rc22*e1rc12*e1rc20+2*e0rc21*e1rc11*e1rc20+2*e0rc20*e1rc10*e1rc20+e0rc10*e1rc12_2+2*e0rc12*e1rc10*e1rc12+2*e0rc00*e1rc02*e1rc12+2*e0rc02*e1rc00*e1rc12+e0rc10*e1rc11_2+2*e0rc11*e1rc10*e1rc11+2*e0rc00*e1rc01*e1rc11+2*e0rc01*e1rc00*e1rc11+3*e0rc10*e1rc10_2-2*e0rc02*e1rc02*e1rc10-2*e0rc01*e1rc01*e1rc10+2*e0rc00*e1rc00*e1rc10-e0rc10*e1rc02_2+2*e0rc12*e1rc00*e1rc02-e0rc10*e1rc01_2+2*e0rc11*e1rc00*e1rc01+e0rc10*e1rc00_2;
    AMatrix(3, 4) = -2*e0rc10*e1rc22*e2rc22+2*e0rc12*e1rc20*e2rc22+2*e0rc20*e1rc12*e2rc22-2*e0rc22*e1rc10*e2rc22-2*e0rc10*e1rc21*e2rc21+2*e0rc11*e1rc20*e2rc21+2*e0rc20*e1rc11*e2rc21-2*e0rc21*e1rc10*e2rc21+2*e0rc12*e1rc22*e2rc20+2*e0rc11*e1rc21*e2rc20+2*e0rc10*e1rc20*e2rc20+2*e0rc22*e1rc12*e2rc20+2*e0rc21*e1rc11*e2rc20+2*e0rc20*e1rc10*e2rc20+2*e0rc20*e1rc22*e2rc12+2*e0rc22*e1rc20*e2rc12+2*e0rc10*e1rc12*e2rc12+2*e0rc12*e1rc10*e2rc12+2*e0rc00*e1rc02*e2rc12+2*e0rc02*e1rc00*e2rc12+2*e0rc20*e1rc21*e2rc11+2*e0rc21*e1rc20*e2rc11+2*e0rc10*e1rc11*e2rc11+2*e0rc11*e1rc10*e2rc11+2*e0rc00*e1rc01*e2rc11+2*e0rc01*e1rc00*e2rc11-2*e0rc22*e1rc22*e2rc10-2*e0rc21*e1rc21*e2rc10+2*e0rc20*e1rc20*e2rc10+2*e0rc12*e1rc12*e2rc10+2*e0rc11*e1rc11*e2rc10+6*e0rc10*e1rc10*e2rc10-2*e0rc02*e1rc02*e2rc10-2*e0rc01*e1rc01*e2rc10+2*e0rc00*e1rc00*e2rc10+2*e0rc00*e1rc12*e2rc02-2*e0rc02*e1rc10*e2rc02-2*e0rc10*e1rc02*e2rc02+2*e0rc12*e1rc00*e2rc02+2*e0rc00*e1rc11*e2rc01-2*e0rc01*e1rc10*e2rc01-2*e0rc10*e1rc01*e2rc01+2*e0rc11*e1rc00*e2rc01+2*e0rc02*e1rc12*e2rc00+2*e0rc01*e1rc11*e2rc00+2*e0rc00*e1rc10*e2rc00+2*e0rc12*e1rc02*e2rc00+2*e0rc11*e1rc01*e2rc00+2*e0rc10*e1rc00*e2rc00;
    AMatrix(3, 5) = -e0rc10*e2rc22_2+2*e0rc12*e2rc20*e2rc22+2*e0rc20*e2rc12*e2rc22-2*e0rc22*e2rc10*e2rc22-e0rc10*e2rc21_2+2*e0rc11*e2rc20*e2rc21+2*e0rc20*e2rc11*e2rc21-2*e0rc21*e2rc10*e2rc21+e0rc10*e2rc20_2+2*e0rc22*e2rc12*e2rc20+2*e0rc21*e2rc11*e2rc20+2*e0rc20*e2rc10*e2rc20+e0rc10*e2rc12_2+2*e0rc12*e2rc10*e2rc12+2*e0rc00*e2rc02*e2rc12+2*e0rc02*e2rc00*e2rc12+e0rc10*e2rc11_2+2*e0rc11*e2rc10*e2rc11+2*e0rc00*e2rc01*e2rc11+2*e0rc01*e2rc00*e2rc11+3*e0rc10*e2rc10_2-2*e0rc02*e2rc02*e2rc10-2*e0rc01*e2rc01*e2rc10+2*e0rc00*e2rc00*e2rc10-e0rc10*e2rc02_2+2*e0rc12*e2rc00*e2rc02-e0rc10*e2rc01_2+2*e0rc11*e2rc00*e2rc01+e0rc10*e2rc00_2;
    AMatrix(3, 6) = -e1rc10*e1rc22_2+2*e1rc12*e1rc20*e1rc22-e1rc10*e1rc21_2+2*e1rc11*e1rc20*e1rc21+e1rc10*e1rc20_2+e1rc10*e1rc12_2+2*e1rc00*e1rc02*e1rc12+e1rc10*e1rc11_2+2*e1rc00*e1rc01*e1rc11+e1rc10_3-e1rc02_2*e1rc10-e1rc01_2*e1rc10+e1rc00_2*e1rc10;
    AMatrix(3, 7) = -2*e1rc10*e1rc22*e2rc22+2*e1rc12*e1rc20*e2rc22-2*e1rc10*e1rc21*e2rc21+2*e1rc11*e1rc20*e2rc21+2*e1rc12*e1rc22*e2rc20+2*e1rc11*e1rc21*e2rc20+2*e1rc10*e1rc20*e2rc20+2*e1rc20*e1rc22*e2rc12+2*e1rc10*e1rc12*e2rc12+2*e1rc00*e1rc02*e2rc12+2*e1rc20*e1rc21*e2rc11+2*e1rc10*e1rc11*e2rc11+2*e1rc00*e1rc01*e2rc11-e1rc22_2*e2rc10-e1rc21_2*e2rc10+e1rc20_2*e2rc10+e1rc12_2*e2rc10+e1rc11_2*e2rc10+3*e1rc10_2*e2rc10-e1rc02_2*e2rc10-e1rc01_2*e2rc10+e1rc00_2*e2rc10+2*e1rc00*e1rc12*e2rc02-2*e1rc02*e1rc10*e2rc02+2*e1rc00*e1rc11*e2rc01-2*e1rc01*e1rc10*e2rc01+2*e1rc02*e1rc12*e2rc00+2*e1rc01*e1rc11*e2rc00+2*e1rc00*e1rc10*e2rc00;
    AMatrix(3, 8) = -e1rc10*e2rc22_2+2*e1rc12*e2rc20*e2rc22+2*e1rc20*e2rc12*e2rc22-2*e1rc22*e2rc10*e2rc22-e1rc10*e2rc21_2+2*e1rc11*e2rc20*e2rc21+2*e1rc20*e2rc11*e2rc21-2*e1rc21*e2rc10*e2rc21+e1rc10*e2rc20_2+2*e1rc22*e2rc12*e2rc20+2*e1rc21*e2rc11*e2rc20+2*e1rc20*e2rc10*e2rc20+e1rc10*e2rc12_2+2*e1rc12*e2rc10*e2rc12+2*e1rc00*e2rc02*e2rc12+2*e1rc02*e2rc00*e2rc12+e1rc10*e2rc11_2+2*e1rc11*e2rc10*e2rc11+2*e1rc00*e2rc01*e2rc11+2*e1rc01*e2rc00*e2rc11+3*e1rc10*e2rc10_2-2*e1rc02*e2rc02*e2rc10-2*e1rc01*e2rc01*e2rc10+2*e1rc00*e2rc00*e2rc10-e1rc10*e2rc02_2+2*e1rc12*e2rc00*e2rc02-e1rc10*e2rc01_2+2*e1rc11*e2rc00*e2rc01+e1rc10*e2rc00_2;
    AMatrix(3, 9) = -e2rc10*e2rc22_2+2*e2rc12*e2rc20*e2rc22-e2rc10*e2rc21_2+2*e2rc11*e2rc20*e2rc21+e2rc10*e2rc20_2+e2rc10*e2rc12_2+2*e2rc00*e2rc02*e2rc12+e2rc10*e2rc11_2+2*e2rc00*e2rc01*e2rc11+e2rc10_3-e2rc02_2*e2rc10-e2rc01_2*e2rc10+e2rc00_2*e2rc10;
    AMatrix(3, 10) = -2*e0rc10*e0rc22*e3rc22+2*e0rc12*e0rc20*e3rc22-2*e0rc10*e0rc21*e3rc21+2*e0rc11*e0rc20*e3rc21+2*e0rc12*e0rc22*e3rc20+2*e0rc11*e0rc21*e3rc20+2*e0rc10*e0rc20*e3rc20+2*e0rc20*e0rc22*e3rc12+2*e0rc10*e0rc12*e3rc12+2*e0rc00*e0rc02*e3rc12+2*e0rc20*e0rc21*e3rc11+2*e0rc10*e0rc11*e3rc11+2*e0rc00*e0rc01*e3rc11-e0rc22_2*e3rc10-e0rc21_2*e3rc10+e0rc20_2*e3rc10+e0rc12_2*e3rc10+e0rc11_2*e3rc10+3*e0rc10_2*e3rc10-e0rc02_2*e3rc10-e0rc01_2*e3rc10+e0rc00_2*e3rc10+2*e0rc00*e0rc12*e3rc02-2*e0rc02*e0rc10*e3rc02+2*e0rc00*e0rc11*e3rc01-2*e0rc01*e0rc10*e3rc01+2*e0rc02*e0rc12*e3rc00+2*e0rc01*e0rc11*e3rc00+2*e0rc00*e0rc10*e3rc00;
    AMatrix(3, 11) = -2*e0rc10*e1rc22*e3rc22+2*e0rc12*e1rc20*e3rc22+2*e0rc20*e1rc12*e3rc22-2*e0rc22*e1rc10*e3rc22-2*e0rc10*e1rc21*e3rc21+2*e0rc11*e1rc20*e3rc21+2*e0rc20*e1rc11*e3rc21-2*e0rc21*e1rc10*e3rc21+2*e0rc12*e1rc22*e3rc20+2*e0rc11*e1rc21*e3rc20+2*e0rc10*e1rc20*e3rc20+2*e0rc22*e1rc12*e3rc20+2*e0rc21*e1rc11*e3rc20+2*e0rc20*e1rc10*e3rc20+2*e0rc20*e1rc22*e3rc12+2*e0rc22*e1rc20*e3rc12+2*e0rc10*e1rc12*e3rc12+2*e0rc12*e1rc10*e3rc12+2*e0rc00*e1rc02*e3rc12+2*e0rc02*e1rc00*e3rc12+2*e0rc20*e1rc21*e3rc11+2*e0rc21*e1rc20*e3rc11+2*e0rc10*e1rc11*e3rc11+2*e0rc11*e1rc10*e3rc11+2*e0rc00*e1rc01*e3rc11+2*e0rc01*e1rc00*e3rc11-2*e0rc22*e1rc22*e3rc10-2*e0rc21*e1rc21*e3rc10+2*e0rc20*e1rc20*e3rc10+2*e0rc12*e1rc12*e3rc10+2*e0rc11*e1rc11*e3rc10+6*e0rc10*e1rc10*e3rc10-2*e0rc02*e1rc02*e3rc10-2*e0rc01*e1rc01*e3rc10+2*e0rc00*e1rc00*e3rc10+2*e0rc00*e1rc12*e3rc02-2*e0rc02*e1rc10*e3rc02-2*e0rc10*e1rc02*e3rc02+2*e0rc12*e1rc00*e3rc02+2*e0rc00*e1rc11*e3rc01-2*e0rc01*e1rc10*e3rc01-2*e0rc10*e1rc01*e3rc01+2*e0rc11*e1rc00*e3rc01+2*e0rc02*e1rc12*e3rc00+2*e0rc01*e1rc11*e3rc00+2*e0rc00*e1rc10*e3rc00+2*e0rc12*e1rc02*e3rc00+2*e0rc11*e1rc01*e3rc00+2*e0rc10*e1rc00*e3rc00;
    AMatrix(3, 12) = -2*e0rc10*e2rc22*e3rc22+2*e0rc12*e2rc20*e3rc22+2*e0rc20*e2rc12*e3rc22-2*e0rc22*e2rc10*e3rc22-2*e0rc10*e2rc21*e3rc21+2*e0rc11*e2rc20*e3rc21+2*e0rc20*e2rc11*e3rc21-2*e0rc21*e2rc10*e3rc21+2*e0rc12*e2rc22*e3rc20+2*e0rc11*e2rc21*e3rc20+2*e0rc10*e2rc20*e3rc20+2*e0rc22*e2rc12*e3rc20+2*e0rc21*e2rc11*e3rc20+2*e0rc20*e2rc10*e3rc20+2*e0rc20*e2rc22*e3rc12+2*e0rc22*e2rc20*e3rc12+2*e0rc10*e2rc12*e3rc12+2*e0rc12*e2rc10*e3rc12+2*e0rc00*e2rc02*e3rc12+2*e0rc02*e2rc00*e3rc12+2*e0rc20*e2rc21*e3rc11+2*e0rc21*e2rc20*e3rc11+2*e0rc10*e2rc11*e3rc11+2*e0rc11*e2rc10*e3rc11+2*e0rc00*e2rc01*e3rc11+2*e0rc01*e2rc00*e3rc11-2*e0rc22*e2rc22*e3rc10-2*e0rc21*e2rc21*e3rc10+2*e0rc20*e2rc20*e3rc10+2*e0rc12*e2rc12*e3rc10+2*e0rc11*e2rc11*e3rc10+6*e0rc10*e2rc10*e3rc10-2*e0rc02*e2rc02*e3rc10-2*e0rc01*e2rc01*e3rc10+2*e0rc00*e2rc00*e3rc10+2*e0rc00*e2rc12*e3rc02-2*e0rc02*e2rc10*e3rc02-2*e0rc10*e2rc02*e3rc02+2*e0rc12*e2rc00*e3rc02+2*e0rc00*e2rc11*e3rc01-2*e0rc01*e2rc10*e3rc01-2*e0rc10*e2rc01*e3rc01+2*e0rc11*e2rc00*e3rc01+2*e0rc02*e2rc12*e3rc00+2*e0rc01*e2rc11*e3rc00+2*e0rc00*e2rc10*e3rc00+2*e0rc12*e2rc02*e3rc00+2*e0rc11*e2rc01*e3rc00+2*e0rc10*e2rc00*e3rc00;
    AMatrix(3, 13) = -2*e1rc10*e1rc22*e3rc22+2*e1rc12*e1rc20*e3rc22-2*e1rc10*e1rc21*e3rc21+2*e1rc11*e1rc20*e3rc21+2*e1rc12*e1rc22*e3rc20+2*e1rc11*e1rc21*e3rc20+2*e1rc10*e1rc20*e3rc20+2*e1rc20*e1rc22*e3rc12+2*e1rc10*e1rc12*e3rc12+2*e1rc00*e1rc02*e3rc12+2*e1rc20*e1rc21*e3rc11+2*e1rc10*e1rc11*e3rc11+2*e1rc00*e1rc01*e3rc11-e1rc22_2*e3rc10-e1rc21_2*e3rc10+e1rc20_2*e3rc10+e1rc12_2*e3rc10+e1rc11_2*e3rc10+3*e1rc10_2*e3rc10-e1rc02_2*e3rc10-e1rc01_2*e3rc10+e1rc00_2*e3rc10+2*e1rc00*e1rc12*e3rc02-2*e1rc02*e1rc10*e3rc02+2*e1rc00*e1rc11*e3rc01-2*e1rc01*e1rc10*e3rc01+2*e1rc02*e1rc12*e3rc00+2*e1rc01*e1rc11*e3rc00+2*e1rc00*e1rc10*e3rc00;
    AMatrix(3, 14) = -2*e1rc10*e2rc22*e3rc22+2*e1rc12*e2rc20*e3rc22+2*e1rc20*e2rc12*e3rc22-2*e1rc22*e2rc10*e3rc22-2*e1rc10*e2rc21*e3rc21+2*e1rc11*e2rc20*e3rc21+2*e1rc20*e2rc11*e3rc21-2*e1rc21*e2rc10*e3rc21+2*e1rc12*e2rc22*e3rc20+2*e1rc11*e2rc21*e3rc20+2*e1rc10*e2rc20*e3rc20+2*e1rc22*e2rc12*e3rc20+2*e1rc21*e2rc11*e3rc20+2*e1rc20*e2rc10*e3rc20+2*e1rc20*e2rc22*e3rc12+2*e1rc22*e2rc20*e3rc12+2*e1rc10*e2rc12*e3rc12+2*e1rc12*e2rc10*e3rc12+2*e1rc00*e2rc02*e3rc12+2*e1rc02*e2rc00*e3rc12+2*e1rc20*e2rc21*e3rc11+2*e1rc21*e2rc20*e3rc11+2*e1rc10*e2rc11*e3rc11+2*e1rc11*e2rc10*e3rc11+2*e1rc00*e2rc01*e3rc11+2*e1rc01*e2rc00*e3rc11-2*e1rc22*e2rc22*e3rc10-2*e1rc21*e2rc21*e3rc10+2*e1rc20*e2rc20*e3rc10+2*e1rc12*e2rc12*e3rc10+2*e1rc11*e2rc11*e3rc10+6*e1rc10*e2rc10*e3rc10-2*e1rc02*e2rc02*e3rc10-2*e1rc01*e2rc01*e3rc10+2*e1rc00*e2rc00*e3rc10+2*e1rc00*e2rc12*e3rc02-2*e1rc02*e2rc10*e3rc02-2*e1rc10*e2rc02*e3rc02+2*e1rc12*e2rc00*e3rc02+2*e1rc00*e2rc11*e3rc01-2*e1rc01*e2rc10*e3rc01-2*e1rc10*e2rc01*e3rc01+2*e1rc11*e2rc00*e3rc01+2*e1rc02*e2rc12*e3rc00+2*e1rc01*e2rc11*e3rc00+2*e1rc00*e2rc10*e3rc00+2*e1rc12*e2rc02*e3rc00+2*e1rc11*e2rc01*e3rc00+2*e1rc10*e2rc00*e3rc00;
    AMatrix(3, 15) = -2*e2rc10*e2rc22*e3rc22+2*e2rc12*e2rc20*e3rc22-2*e2rc10*e2rc21*e3rc21+2*e2rc11*e2rc20*e3rc21+2*e2rc12*e2rc22*e3rc20+2*e2rc11*e2rc21*e3rc20+2*e2rc10*e2rc20*e3rc20+2*e2rc20*e2rc22*e3rc12+2*e2rc10*e2rc12*e3rc12+2*e2rc00*e2rc02*e3rc12+2*e2rc20*e2rc21*e3rc11+2*e2rc10*e2rc11*e3rc11+2*e2rc00*e2rc01*e3rc11-e2rc22_2*e3rc10-e2rc21_2*e3rc10+e2rc20_2*e3rc10+e2rc12_2*e3rc10+e2rc11_2*e3rc10+3*e2rc10_2*e3rc10-e2rc02_2*e3rc10-e2rc01_2*e3rc10+e2rc00_2*e3rc10+2*e2rc00*e2rc12*e3rc02-2*e2rc02*e2rc10*e3rc02+2*e2rc00*e2rc11*e3rc01-2*e2rc01*e2rc10*e3rc01+2*e2rc02*e2rc12*e3rc00+2*e2rc01*e2rc11*e3rc00+2*e2rc00*e2rc10*e3rc00;
    AMatrix(3, 16) = -e0rc10*e3rc22_2+2*e0rc12*e3rc20*e3rc22+2*e0rc20*e3rc12*e3rc22-2*e0rc22*e3rc10*e3rc22-e0rc10*e3rc21_2+2*e0rc11*e3rc20*e3rc21+2*e0rc20*e3rc11*e3rc21-2*e0rc21*e3rc10*e3rc21+e0rc10*e3rc20_2+2*e0rc22*e3rc12*e3rc20+2*e0rc21*e3rc11*e3rc20+2*e0rc20*e3rc10*e3rc20+e0rc10*e3rc12_2+2*e0rc12*e3rc10*e3rc12+2*e0rc00*e3rc02*e3rc12+2*e0rc02*e3rc00*e3rc12+e0rc10*e3rc11_2+2*e0rc11*e3rc10*e3rc11+2*e0rc00*e3rc01*e3rc11+2*e0rc01*e3rc00*e3rc11+3*e0rc10*e3rc10_2-2*e0rc02*e3rc02*e3rc10-2*e0rc01*e3rc01*e3rc10+2*e0rc00*e3rc00*e3rc10-e0rc10*e3rc02_2+2*e0rc12*e3rc00*e3rc02-e0rc10*e3rc01_2+2*e0rc11*e3rc00*e3rc01+e0rc10*e3rc00_2;
    AMatrix(3, 17) = -e1rc10*e3rc22_2+2*e1rc12*e3rc20*e3rc22+2*e1rc20*e3rc12*e3rc22-2*e1rc22*e3rc10*e3rc22-e1rc10*e3rc21_2+2*e1rc11*e3rc20*e3rc21+2*e1rc20*e3rc11*e3rc21-2*e1rc21*e3rc10*e3rc21+e1rc10*e3rc20_2+2*e1rc22*e3rc12*e3rc20+2*e1rc21*e3rc11*e3rc20+2*e1rc20*e3rc10*e3rc20+e1rc10*e3rc12_2+2*e1rc12*e3rc10*e3rc12+2*e1rc00*e3rc02*e3rc12+2*e1rc02*e3rc00*e3rc12+e1rc10*e3rc11_2+2*e1rc11*e3rc10*e3rc11+2*e1rc00*e3rc01*e3rc11+2*e1rc01*e3rc00*e3rc11+3*e1rc10*e3rc10_2-2*e1rc02*e3rc02*e3rc10-2*e1rc01*e3rc01*e3rc10+2*e1rc00*e3rc00*e3rc10-e1rc10*e3rc02_2+2*e1rc12*e3rc00*e3rc02-e1rc10*e3rc01_2+2*e1rc11*e3rc00*e3rc01+e1rc10*e3rc00_2;
    AMatrix(3, 18) = -e2rc10*e3rc22_2+2*e2rc12*e3rc20*e3rc22+2*e2rc20*e3rc12*e3rc22-2*e2rc22*e3rc10*e3rc22-e2rc10*e3rc21_2+2*e2rc11*e3rc20*e3rc21+2*e2rc20*e3rc11*e3rc21-2*e2rc21*e3rc10*e3rc21+e2rc10*e3rc20_2+2*e2rc22*e3rc12*e3rc20+2*e2rc21*e3rc11*e3rc20+2*e2rc20*e3rc10*e3rc20+e2rc10*e3rc12_2+2*e2rc12*e3rc10*e3rc12+2*e2rc00*e3rc02*e3rc12+2*e2rc02*e3rc00*e3rc12+e2rc10*e3rc11_2+2*e2rc11*e3rc10*e3rc11+2*e2rc00*e3rc01*e3rc11+2*e2rc01*e3rc00*e3rc11+3*e2rc10*e3rc10_2-2*e2rc02*e3rc02*e3rc10-2*e2rc01*e3rc01*e3rc10+2*e2rc00*e3rc00*e3rc10-e2rc10*e3rc02_2+2*e2rc12*e3rc00*e3rc02-e2rc10*e3rc01_2+2*e2rc11*e3rc00*e3rc01+e2rc10*e3rc00_2;
    AMatrix(3, 19) = -e3rc10*e3rc22_2+2*e3rc12*e3rc20*e3rc22-e3rc10*e3rc21_2+2*e3rc11*e3rc20*e3rc21+e3rc10*e3rc20_2+e3rc10*e3rc12_2+2*e3rc00*e3rc02*e3rc12+e3rc10*e3rc11_2+2*e3rc00*e3rc01*e3rc11+e3rc10_3-e3rc02_2*e3rc10-e3rc01_2*e3rc10+e3rc00_2*e3rc10;
    AMatrix(4, 0) = -e0rc11*e0rc22_2+2*e0rc12*e0rc21*e0rc22+e0rc11*e0rc21_2+2*e0rc10*e0rc20*e0rc21-e0rc11*e0rc20_2+e0rc11*e0rc12_2+2*e0rc01*e0rc02*e0rc12+e0rc11_3+e0rc10_2*e0rc11-e0rc02_2*e0rc11+e0rc01_2*e0rc11-e0rc00_2*e0rc11+2*e0rc00*e0rc01*e0rc10;
    AMatrix(4, 1) = -2*e0rc11*e0rc22*e1rc22+2*e0rc12*e0rc21*e1rc22+2*e0rc12*e0rc22*e1rc21+2*e0rc11*e0rc21*e1rc21+2*e0rc10*e0rc20*e1rc21+2*e0rc10*e0rc21*e1rc20-2*e0rc11*e0rc20*e1rc20+2*e0rc21*e0rc22*e1rc12+2*e0rc11*e0rc12*e1rc12+2*e0rc01*e0rc02*e1rc12-e0rc22_2*e1rc11+e0rc21_2*e1rc11-e0rc20_2*e1rc11+e0rc12_2*e1rc11+3*e0rc11_2*e1rc11+e0rc10_2*e1rc11-e0rc02_2*e1rc11+e0rc01_2*e1rc11-e0rc00_2*e1rc11+2*e0rc20*e0rc21*e1rc10+2*e0rc10*e0rc11*e1rc10+2*e0rc00*e0rc01*e1rc10+2*e0rc01*e0rc12*e1rc02-2*e0rc02*e0rc11*e1rc02+2*e0rc02*e0rc12*e1rc01+2*e0rc01*e0rc11*e1rc01+2*e0rc00*e0rc10*e1rc01-2*e0rc00*e0rc11*e1rc00+2*e0rc01*e0rc10*e1rc00;
    AMatrix(4, 2) = -2*e0rc11*e0rc22*e2rc22+2*e0rc12*e0rc21*e2rc22+2*e0rc12*e0rc22*e2rc21+2*e0rc11*e0rc21*e2rc21+2*e0rc10*e0rc20*e2rc21+2*e0rc10*e0rc21*e2rc20-2*e0rc11*e0rc20*e2rc20+2*e0rc21*e0rc22*e2rc12+2*e0rc11*e0rc12*e2rc12+2*e0rc01*e0rc02*e2rc12-e0rc22_2*e2rc11+e0rc21_2*e2rc11-e0rc20_2*e2rc11+e0rc12_2*e2rc11+3*e0rc11_2*e2rc11+e0rc10_2*e2rc11-e0rc02_2*e2rc11+e0rc01_2*e2rc11-e0rc00_2*e2rc11+2*e0rc20*e0rc21*e2rc10+2*e0rc10*e0rc11*e2rc10+2*e0rc00*e0rc01*e2rc10+2*e0rc01*e0rc12*e2rc02-2*e0rc02*e0rc11*e2rc02+2*e0rc02*e0rc12*e2rc01+2*e0rc01*e0rc11*e2rc01+2*e0rc00*e0rc10*e2rc01-2*e0rc00*e0rc11*e2rc00+2*e0rc01*e0rc10*e2rc00;
    AMatrix(4, 3) = -e0rc11*e1rc22_2+2*e0rc12*e1rc21*e1rc22+2*e0rc21*e1rc12*e1rc22-2*e0rc22*e1rc11*e1rc22+e0rc11*e1rc21_2+2*e0rc10*e1rc20*e1rc21+2*e0rc22*e1rc12*e1rc21+2*e0rc21*e1rc11*e1rc21+2*e0rc20*e1rc10*e1rc21-e0rc11*e1rc20_2-2*e0rc20*e1rc11*e1rc20+2*e0rc21*e1rc10*e1rc20+e0rc11*e1rc12_2+2*e0rc12*e1rc11*e1rc12+2*e0rc01*e1rc02*e1rc12+2*e0rc02*e1rc01*e1rc12+3*e0rc11*e1rc11_2+2*e0rc10*e1rc10*e1rc11-2*e0rc02*e1rc02*e1rc11+2*e0rc01*e1rc01*e1rc11-2*e0rc00*e1rc00*e1rc11+e0rc11*e1rc10_2+2*e0rc00*e1rc01*e1rc10+2*e0rc01*e1rc00*e1rc10-e0rc11*e1rc02_2+2*e0rc12*e1rc01*e1rc02+e0rc11*e1rc01_2+2*e0rc10*e1rc00*e1rc01-e0rc11*e1rc00_2;
    AMatrix(4, 4) = -2*e0rc11*e1rc22*e2rc22+2*e0rc12*e1rc21*e2rc22+2*e0rc21*e1rc12*e2rc22-2*e0rc22*e1rc11*e2rc22+2*e0rc12*e1rc22*e2rc21+2*e0rc11*e1rc21*e2rc21+2*e0rc10*e1rc20*e2rc21+2*e0rc22*e1rc12*e2rc21+2*e0rc21*e1rc11*e2rc21+2*e0rc20*e1rc10*e2rc21+2*e0rc10*e1rc21*e2rc20-2*e0rc11*e1rc20*e2rc20-2*e0rc20*e1rc11*e2rc20+2*e0rc21*e1rc10*e2rc20+2*e0rc21*e1rc22*e2rc12+2*e0rc22*e1rc21*e2rc12+2*e0rc11*e1rc12*e2rc12+2*e0rc12*e1rc11*e2rc12+2*e0rc01*e1rc02*e2rc12+2*e0rc02*e1rc01*e2rc12-2*e0rc22*e1rc22*e2rc11+2*e0rc21*e1rc21*e2rc11-2*e0rc20*e1rc20*e2rc11+2*e0rc12*e1rc12*e2rc11+6*e0rc11*e1rc11*e2rc11+2*e0rc10*e1rc10*e2rc11-2*e0rc02*e1rc02*e2rc11+2*e0rc01*e1rc01*e2rc11-2*e0rc00*e1rc00*e2rc11+2*e0rc20*e1rc21*e2rc10+2*e0rc21*e1rc20*e2rc10+2*e0rc10*e1rc11*e2rc10+2*e0rc11*e1rc10*e2rc10+2*e0rc00*e1rc01*e2rc10+2*e0rc01*e1rc00*e2rc10+2*e0rc01*e1rc12*e2rc02-2*e0rc02*e1rc11*e2rc02-2*e0rc11*e1rc02*e2rc02+2*e0rc12*e1rc01*e2rc02+2*e0rc02*e1rc12*e2rc01+2*e0rc01*e1rc11*e2rc01+2*e0rc00*e1rc10*e2rc01+2*e0rc12*e1rc02*e2rc01+2*e0rc11*e1rc01*e2rc01+2*e0rc10*e1rc00*e2rc01-2*e0rc00*e1rc11*e2rc00+2*e0rc01*e1rc10*e2rc00+2*e0rc10*e1rc01*e2rc00-2*e0rc11*e1rc00*e2rc00;
    AMatrix(4, 5) = -e0rc11*e2rc22_2+2*e0rc12*e2rc21*e2rc22+2*e0rc21*e2rc12*e2rc22-2*e0rc22*e2rc11*e2rc22+e0rc11*e2rc21_2+2*e0rc10*e2rc20*e2rc21+2*e0rc22*e2rc12*e2rc21+2*e0rc21*e2rc11*e2rc21+2*e0rc20*e2rc10*e2rc21-e0rc11*e2rc20_2-2*e0rc20*e2rc11*e2rc20+2*e0rc21*e2rc10*e2rc20+e0rc11*e2rc12_2+2*e0rc12*e2rc11*e2rc12+2*e0rc01*e2rc02*e2rc12+2*e0rc02*e2rc01*e2rc12+3*e0rc11*e2rc11_2+2*e0rc10*e2rc10*e2rc11-2*e0rc02*e2rc02*e2rc11+2*e0rc01*e2rc01*e2rc11-2*e0rc00*e2rc00*e2rc11+e0rc11*e2rc10_2+2*e0rc00*e2rc01*e2rc10+2*e0rc01*e2rc00*e2rc10-e0rc11*e2rc02_2+2*e0rc12*e2rc01*e2rc02+e0rc11*e2rc01_2+2*e0rc10*e2rc00*e2rc01-e0rc11*e2rc00_2;
    AMatrix(4, 6) = -e1rc11*e1rc22_2+2*e1rc12*e1rc21*e1rc22+e1rc11*e1rc21_2+2*e1rc10*e1rc20*e1rc21-e1rc11*e1rc20_2+e1rc11*e1rc12_2+2*e1rc01*e1rc02*e1rc12+e1rc11_3+e1rc10_2*e1rc11-e1rc02_2*e1rc11+e1rc01_2*e1rc11-e1rc00_2*e1rc11+2*e1rc00*e1rc01*e1rc10;
    AMatrix(4, 7) = -2*e1rc11*e1rc22*e2rc22+2*e1rc12*e1rc21*e2rc22+2*e1rc12*e1rc22*e2rc21+2*e1rc11*e1rc21*e2rc21+2*e1rc10*e1rc20*e2rc21+2*e1rc10*e1rc21*e2rc20-2*e1rc11*e1rc20*e2rc20+2*e1rc21*e1rc22*e2rc12+2*e1rc11*e1rc12*e2rc12+2*e1rc01*e1rc02*e2rc12-e1rc22_2*e2rc11+e1rc21_2*e2rc11-e1rc20_2*e2rc11+e1rc12_2*e2rc11+3*e1rc11_2*e2rc11+e1rc10_2*e2rc11-e1rc02_2*e2rc11+e1rc01_2*e2rc11-e1rc00_2*e2rc11+2*e1rc20*e1rc21*e2rc10+2*e1rc10*e1rc11*e2rc10+2*e1rc00*e1rc01*e2rc10+2*e1rc01*e1rc12*e2rc02-2*e1rc02*e1rc11*e2rc02+2*e1rc02*e1rc12*e2rc01+2*e1rc01*e1rc11*e2rc01+2*e1rc00*e1rc10*e2rc01-2*e1rc00*e1rc11*e2rc00+2*e1rc01*e1rc10*e2rc00;
    AMatrix(4, 8) = -e1rc11*e2rc22_2+2*e1rc12*e2rc21*e2rc22+2*e1rc21*e2rc12*e2rc22-2*e1rc22*e2rc11*e2rc22+e1rc11*e2rc21_2+2*e1rc10*e2rc20*e2rc21+2*e1rc22*e2rc12*e2rc21+2*e1rc21*e2rc11*e2rc21+2*e1rc20*e2rc10*e2rc21-e1rc11*e2rc20_2-2*e1rc20*e2rc11*e2rc20+2*e1rc21*e2rc10*e2rc20+e1rc11*e2rc12_2+2*e1rc12*e2rc11*e2rc12+2*e1rc01*e2rc02*e2rc12+2*e1rc02*e2rc01*e2rc12+3*e1rc11*e2rc11_2+2*e1rc10*e2rc10*e2rc11-2*e1rc02*e2rc02*e2rc11+2*e1rc01*e2rc01*e2rc11-2*e1rc00*e2rc00*e2rc11+e1rc11*e2rc10_2+2*e1rc00*e2rc01*e2rc10+2*e1rc01*e2rc00*e2rc10-e1rc11*e2rc02_2+2*e1rc12*e2rc01*e2rc02+e1rc11*e2rc01_2+2*e1rc10*e2rc00*e2rc01-e1rc11*e2rc00_2;
    AMatrix(4, 9) = -e2rc11*e2rc22_2+2*e2rc12*e2rc21*e2rc22+e2rc11*e2rc21_2+2*e2rc10*e2rc20*e2rc21-e2rc11*e2rc20_2+e2rc11*e2rc12_2+2*e2rc01*e2rc02*e2rc12+e2rc11_3+e2rc10_2*e2rc11-e2rc02_2*e2rc11+e2rc01_2*e2rc11-e2rc00_2*e2rc11+2*e2rc00*e2rc01*e2rc10;
    AMatrix(4, 10) = -2*e0rc11*e0rc22*e3rc22+2*e0rc12*e0rc21*e3rc22+2*e0rc12*e0rc22*e3rc21+2*e0rc11*e0rc21*e3rc21+2*e0rc10*e0rc20*e3rc21+2*e0rc10*e0rc21*e3rc20-2*e0rc11*e0rc20*e3rc20+2*e0rc21*e0rc22*e3rc12+2*e0rc11*e0rc12*e3rc12+2*e0rc01*e0rc02*e3rc12-e0rc22_2*e3rc11+e0rc21_2*e3rc11-e0rc20_2*e3rc11+e0rc12_2*e3rc11+3*e0rc11_2*e3rc11+e0rc10_2*e3rc11-e0rc02_2*e3rc11+e0rc01_2*e3rc11-e0rc00_2*e3rc11+2*e0rc20*e0rc21*e3rc10+2*e0rc10*e0rc11*e3rc10+2*e0rc00*e0rc01*e3rc10+2*e0rc01*e0rc12*e3rc02-2*e0rc02*e0rc11*e3rc02+2*e0rc02*e0rc12*e3rc01+2*e0rc01*e0rc11*e3rc01+2*e0rc00*e0rc10*e3rc01-2*e0rc00*e0rc11*e3rc00+2*e0rc01*e0rc10*e3rc00;
    AMatrix(4, 11) = -2*e0rc11*e1rc22*e3rc22+2*e0rc12*e1rc21*e3rc22+2*e0rc21*e1rc12*e3rc22-2*e0rc22*e1rc11*e3rc22+2*e0rc12*e1rc22*e3rc21+2*e0rc11*e1rc21*e3rc21+2*e0rc10*e1rc20*e3rc21+2*e0rc22*e1rc12*e3rc21+2*e0rc21*e1rc11*e3rc21+2*e0rc20*e1rc10*e3rc21+2*e0rc10*e1rc21*e3rc20-2*e0rc11*e1rc20*e3rc20-2*e0rc20*e1rc11*e3rc20+2*e0rc21*e1rc10*e3rc20+2*e0rc21*e1rc22*e3rc12+2*e0rc22*e1rc21*e3rc12+2*e0rc11*e1rc12*e3rc12+2*e0rc12*e1rc11*e3rc12+2*e0rc01*e1rc02*e3rc12+2*e0rc02*e1rc01*e3rc12-2*e0rc22*e1rc22*e3rc11+2*e0rc21*e1rc21*e3rc11-2*e0rc20*e1rc20*e3rc11+2*e0rc12*e1rc12*e3rc11+6*e0rc11*e1rc11*e3rc11+2*e0rc10*e1rc10*e3rc11-2*e0rc02*e1rc02*e3rc11+2*e0rc01*e1rc01*e3rc11-2*e0rc00*e1rc00*e3rc11+2*e0rc20*e1rc21*e3rc10+2*e0rc21*e1rc20*e3rc10+2*e0rc10*e1rc11*e3rc10+2*e0rc11*e1rc10*e3rc10+2*e0rc00*e1rc01*e3rc10+2*e0rc01*e1rc00*e3rc10+2*e0rc01*e1rc12*e3rc02-2*e0rc02*e1rc11*e3rc02-2*e0rc11*e1rc02*e3rc02+2*e0rc12*e1rc01*e3rc02+2*e0rc02*e1rc12*e3rc01+2*e0rc01*e1rc11*e3rc01+2*e0rc00*e1rc10*e3rc01+2*e0rc12*e1rc02*e3rc01+2*e0rc11*e1rc01*e3rc01+2*e0rc10*e1rc00*e3rc01-2*e0rc00*e1rc11*e3rc00+2*e0rc01*e1rc10*e3rc00+2*e0rc10*e1rc01*e3rc00-2*e0rc11*e1rc00*e3rc00;
    AMatrix(4, 12) = -2*e0rc11*e2rc22*e3rc22+2*e0rc12*e2rc21*e3rc22+2*e0rc21*e2rc12*e3rc22-2*e0rc22*e2rc11*e3rc22+2*e0rc12*e2rc22*e3rc21+2*e0rc11*e2rc21*e3rc21+2*e0rc10*e2rc20*e3rc21+2*e0rc22*e2rc12*e3rc21+2*e0rc21*e2rc11*e3rc21+2*e0rc20*e2rc10*e3rc21+2*e0rc10*e2rc21*e3rc20-2*e0rc11*e2rc20*e3rc20-2*e0rc20*e2rc11*e3rc20+2*e0rc21*e2rc10*e3rc20+2*e0rc21*e2rc22*e3rc12+2*e0rc22*e2rc21*e3rc12+2*e0rc11*e2rc12*e3rc12+2*e0rc12*e2rc11*e3rc12+2*e0rc01*e2rc02*e3rc12+2*e0rc02*e2rc01*e3rc12-2*e0rc22*e2rc22*e3rc11+2*e0rc21*e2rc21*e3rc11-2*e0rc20*e2rc20*e3rc11+2*e0rc12*e2rc12*e3rc11+6*e0rc11*e2rc11*e3rc11+2*e0rc10*e2rc10*e3rc11-2*e0rc02*e2rc02*e3rc11+2*e0rc01*e2rc01*e3rc11-2*e0rc00*e2rc00*e3rc11+2*e0rc20*e2rc21*e3rc10+2*e0rc21*e2rc20*e3rc10+2*e0rc10*e2rc11*e3rc10+2*e0rc11*e2rc10*e3rc10+2*e0rc00*e2rc01*e3rc10+2*e0rc01*e2rc00*e3rc10+2*e0rc01*e2rc12*e3rc02-2*e0rc02*e2rc11*e3rc02-2*e0rc11*e2rc02*e3rc02+2*e0rc12*e2rc01*e3rc02+2*e0rc02*e2rc12*e3rc01+2*e0rc01*e2rc11*e3rc01+2*e0rc00*e2rc10*e3rc01+2*e0rc12*e2rc02*e3rc01+2*e0rc11*e2rc01*e3rc01+2*e0rc10*e2rc00*e3rc01-2*e0rc00*e2rc11*e3rc00+2*e0rc01*e2rc10*e3rc00+2*e0rc10*e2rc01*e3rc00-2*e0rc11*e2rc00*e3rc00;
    AMatrix(4, 13) = -2*e1rc11*e1rc22*e3rc22+2*e1rc12*e1rc21*e3rc22+2*e1rc12*e1rc22*e3rc21+2*e1rc11*e1rc21*e3rc21+2*e1rc10*e1rc20*e3rc21+2*e1rc10*e1rc21*e3rc20-2*e1rc11*e1rc20*e3rc20+2*e1rc21*e1rc22*e3rc12+2*e1rc11*e1rc12*e3rc12+2*e1rc01*e1rc02*e3rc12-e1rc22_2*e3rc11+e1rc21_2*e3rc11-e1rc20_2*e3rc11+e1rc12_2*e3rc11+3*e1rc11_2*e3rc11+e1rc10_2*e3rc11-e1rc02_2*e3rc11+e1rc01_2*e3rc11-e1rc00_2*e3rc11+2*e1rc20*e1rc21*e3rc10+2*e1rc10*e1rc11*e3rc10+2*e1rc00*e1rc01*e3rc10+2*e1rc01*e1rc12*e3rc02-2*e1rc02*e1rc11*e3rc02+2*e1rc02*e1rc12*e3rc01+2*e1rc01*e1rc11*e3rc01+2*e1rc00*e1rc10*e3rc01-2*e1rc00*e1rc11*e3rc00+2*e1rc01*e1rc10*e3rc00;
    AMatrix(4, 14) = -2*e1rc11*e2rc22*e3rc22+2*e1rc12*e2rc21*e3rc22+2*e1rc21*e2rc12*e3rc22-2*e1rc22*e2rc11*e3rc22+2*e1rc12*e2rc22*e3rc21+2*e1rc11*e2rc21*e3rc21+2*e1rc10*e2rc20*e3rc21+2*e1rc22*e2rc12*e3rc21+2*e1rc21*e2rc11*e3rc21+2*e1rc20*e2rc10*e3rc21+2*e1rc10*e2rc21*e3rc20-2*e1rc11*e2rc20*e3rc20-2*e1rc20*e2rc11*e3rc20+2*e1rc21*e2rc10*e3rc20+2*e1rc21*e2rc22*e3rc12+2*e1rc22*e2rc21*e3rc12+2*e1rc11*e2rc12*e3rc12+2*e1rc12*e2rc11*e3rc12+2*e1rc01*e2rc02*e3rc12+2*e1rc02*e2rc01*e3rc12-2*e1rc22*e2rc22*e3rc11+2*e1rc21*e2rc21*e3rc11-2*e1rc20*e2rc20*e3rc11+2*e1rc12*e2rc12*e3rc11+6*e1rc11*e2rc11*e3rc11+2*e1rc10*e2rc10*e3rc11-2*e1rc02*e2rc02*e3rc11+2*e1rc01*e2rc01*e3rc11-2*e1rc00*e2rc00*e3rc11+2*e1rc20*e2rc21*e3rc10+2*e1rc21*e2rc20*e3rc10+2*e1rc10*e2rc11*e3rc10+2*e1rc11*e2rc10*e3rc10+2*e1rc00*e2rc01*e3rc10+2*e1rc01*e2rc00*e3rc10+2*e1rc01*e2rc12*e3rc02-2*e1rc02*e2rc11*e3rc02-2*e1rc11*e2rc02*e3rc02+2*e1rc12*e2rc01*e3rc02+2*e1rc02*e2rc12*e3rc01+2*e1rc01*e2rc11*e3rc01+2*e1rc00*e2rc10*e3rc01+2*e1rc12*e2rc02*e3rc01+2*e1rc11*e2rc01*e3rc01+2*e1rc10*e2rc00*e3rc01-2*e1rc00*e2rc11*e3rc00+2*e1rc01*e2rc10*e3rc00+2*e1rc10*e2rc01*e3rc00-2*e1rc11*e2rc00*e3rc00;
    AMatrix(4, 15) = -2*e2rc11*e2rc22*e3rc22+2*e2rc12*e2rc21*e3rc22+2*e2rc12*e2rc22*e3rc21+2*e2rc11*e2rc21*e3rc21+2*e2rc10*e2rc20*e3rc21+2*e2rc10*e2rc21*e3rc20-2*e2rc11*e2rc20*e3rc20+2*e2rc21*e2rc22*e3rc12+2*e2rc11*e2rc12*e3rc12+2*e2rc01*e2rc02*e3rc12-e2rc22_2*e3rc11+e2rc21_2*e3rc11-e2rc20_2*e3rc11+e2rc12_2*e3rc11+3*e2rc11_2*e3rc11+e2rc10_2*e3rc11-e2rc02_2*e3rc11+e2rc01_2*e3rc11-e2rc00_2*e3rc11+2*e2rc20*e2rc21*e3rc10+2*e2rc10*e2rc11*e3rc10+2*e2rc00*e2rc01*e3rc10+2*e2rc01*e2rc12*e3rc02-2*e2rc02*e2rc11*e3rc02+2*e2rc02*e2rc12*e3rc01+2*e2rc01*e2rc11*e3rc01+2*e2rc00*e2rc10*e3rc01-2*e2rc00*e2rc11*e3rc00+2*e2rc01*e2rc10*e3rc00;
    AMatrix(4, 16) = -e0rc11*e3rc22_2+2*e0rc12*e3rc21*e3rc22+2*e0rc21*e3rc12*e3rc22-2*e0rc22*e3rc11*e3rc22+e0rc11*e3rc21_2+2*e0rc10*e3rc20*e3rc21+2*e0rc22*e3rc12*e3rc21+2*e0rc21*e3rc11*e3rc21+2*e0rc20*e3rc10*e3rc21-e0rc11*e3rc20_2-2*e0rc20*e3rc11*e3rc20+2*e0rc21*e3rc10*e3rc20+e0rc11*e3rc12_2+2*e0rc12*e3rc11*e3rc12+2*e0rc01*e3rc02*e3rc12+2*e0rc02*e3rc01*e3rc12+3*e0rc11*e3rc11_2+2*e0rc10*e3rc10*e3rc11-2*e0rc02*e3rc02*e3rc11+2*e0rc01*e3rc01*e3rc11-2*e0rc00*e3rc00*e3rc11+e0rc11*e3rc10_2+2*e0rc00*e3rc01*e3rc10+2*e0rc01*e3rc00*e3rc10-e0rc11*e3rc02_2+2*e0rc12*e3rc01*e3rc02+e0rc11*e3rc01_2+2*e0rc10*e3rc00*e3rc01-e0rc11*e3rc00_2;
    AMatrix(4, 17) = -e1rc11*e3rc22_2+2*e1rc12*e3rc21*e3rc22+2*e1rc21*e3rc12*e3rc22-2*e1rc22*e3rc11*e3rc22+e1rc11*e3rc21_2+2*e1rc10*e3rc20*e3rc21+2*e1rc22*e3rc12*e3rc21+2*e1rc21*e3rc11*e3rc21+2*e1rc20*e3rc10*e3rc21-e1rc11*e3rc20_2-2*e1rc20*e3rc11*e3rc20+2*e1rc21*e3rc10*e3rc20+e1rc11*e3rc12_2+2*e1rc12*e3rc11*e3rc12+2*e1rc01*e3rc02*e3rc12+2*e1rc02*e3rc01*e3rc12+3*e1rc11*e3rc11_2+2*e1rc10*e3rc10*e3rc11-2*e1rc02*e3rc02*e3rc11+2*e1rc01*e3rc01*e3rc11-2*e1rc00*e3rc00*e3rc11+e1rc11*e3rc10_2+2*e1rc00*e3rc01*e3rc10+2*e1rc01*e3rc00*e3rc10-e1rc11*e3rc02_2+2*e1rc12*e3rc01*e3rc02+e1rc11*e3rc01_2+2*e1rc10*e3rc00*e3rc01-e1rc11*e3rc00_2;
    AMatrix(4, 18) = -e2rc11*e3rc22_2+2*e2rc12*e3rc21*e3rc22+2*e2rc21*e3rc12*e3rc22-2*e2rc22*e3rc11*e3rc22+e2rc11*e3rc21_2+2*e2rc10*e3rc20*e3rc21+2*e2rc22*e3rc12*e3rc21+2*e2rc21*e3rc11*e3rc21+2*e2rc20*e3rc10*e3rc21-e2rc11*e3rc20_2-2*e2rc20*e3rc11*e3rc20+2*e2rc21*e3rc10*e3rc20+e2rc11*e3rc12_2+2*e2rc12*e3rc11*e3rc12+2*e2rc01*e3rc02*e3rc12+2*e2rc02*e3rc01*e3rc12+3*e2rc11*e3rc11_2+2*e2rc10*e3rc10*e3rc11-2*e2rc02*e3rc02*e3rc11+2*e2rc01*e3rc01*e3rc11-2*e2rc00*e3rc00*e3rc11+e2rc11*e3rc10_2+2*e2rc00*e3rc01*e3rc10+2*e2rc01*e3rc00*e3rc10-e2rc11*e3rc02_2+2*e2rc12*e3rc01*e3rc02+e2rc11*e3rc01_2+2*e2rc10*e3rc00*e3rc01-e2rc11*e3rc00_2;
    AMatrix(4, 19) = -e3rc11*e3rc22_2+2*e3rc12*e3rc21*e3rc22+e3rc11*e3rc21_2+2*e3rc10*e3rc20*e3rc21-e3rc11*e3rc20_2+e3rc11*e3rc12_2+2*e3rc01*e3rc02*e3rc12+e3rc11_3+e3rc10_2*e3rc11-e3rc02_2*e3rc11+e3rc01_2*e3rc11-e3rc00_2*e3rc11+2*e3rc00*e3rc01*e3rc10;
    AMatrix(5, 0) = +e0rc12*e0rc22_2+2*e0rc11*e0rc21*e0rc22+2*e0rc10*e0rc20*e0rc22-e0rc12*e0rc21_2-e0rc12*e0rc20_2+e0rc12_3+e0rc11_2*e0rc12+e0rc10_2*e0rc12+e0rc02_2*e0rc12-e0rc01_2*e0rc12-e0rc00_2*e0rc12+2*e0rc01*e0rc02*e0rc11+2*e0rc00*e0rc02*e0rc10;
    AMatrix(5, 1) = +2*e0rc12*e0rc22*e1rc22+2*e0rc11*e0rc21*e1rc22+2*e0rc10*e0rc20*e1rc22+2*e0rc11*e0rc22*e1rc21-2*e0rc12*e0rc21*e1rc21+2*e0rc10*e0rc22*e1rc20-2*e0rc12*e0rc20*e1rc20+e0rc22_2*e1rc12-e0rc21_2*e1rc12-e0rc20_2*e1rc12+3*e0rc12_2*e1rc12+e0rc11_2*e1rc12+e0rc10_2*e1rc12+e0rc02_2*e1rc12-e0rc01_2*e1rc12-e0rc00_2*e1rc12+2*e0rc21*e0rc22*e1rc11+2*e0rc11*e0rc12*e1rc11+2*e0rc01*e0rc02*e1rc11+2*e0rc20*e0rc22*e1rc10+2*e0rc10*e0rc12*e1rc10+2*e0rc00*e0rc02*e1rc10+2*e0rc02*e0rc12*e1rc02+2*e0rc01*e0rc11*e1rc02+2*e0rc00*e0rc10*e1rc02-2*e0rc01*e0rc12*e1rc01+2*e0rc02*e0rc11*e1rc01-2*e0rc00*e0rc12*e1rc00+2*e0rc02*e0rc10*e1rc00;
    AMatrix(5, 2) = +2*e0rc12*e0rc22*e2rc22+2*e0rc11*e0rc21*e2rc22+2*e0rc10*e0rc20*e2rc22+2*e0rc11*e0rc22*e2rc21-2*e0rc12*e0rc21*e2rc21+2*e0rc10*e0rc22*e2rc20-2*e0rc12*e0rc20*e2rc20+e0rc22_2*e2rc12-e0rc21_2*e2rc12-e0rc20_2*e2rc12+3*e0rc12_2*e2rc12+e0rc11_2*e2rc12+e0rc10_2*e2rc12+e0rc02_2*e2rc12-e0rc01_2*e2rc12-e0rc00_2*e2rc12+2*e0rc21*e0rc22*e2rc11+2*e0rc11*e0rc12*e2rc11+2*e0rc01*e0rc02*e2rc11+2*e0rc20*e0rc22*e2rc10+2*e0rc10*e0rc12*e2rc10+2*e0rc00*e0rc02*e2rc10+2*e0rc02*e0rc12*e2rc02+2*e0rc01*e0rc11*e2rc02+2*e0rc00*e0rc10*e2rc02-2*e0rc01*e0rc12*e2rc01+2*e0rc02*e0rc11*e2rc01-2*e0rc00*e0rc12*e2rc00+2*e0rc02*e0rc10*e2rc00;
    AMatrix(5, 3) = +e0rc12*e1rc22_2+2*e0rc11*e1rc21*e1rc22+2*e0rc10*e1rc20*e1rc22+2*e0rc22*e1rc12*e1rc22+2*e0rc21*e1rc11*e1rc22+2*e0rc20*e1rc10*e1rc22-e0rc12*e1rc21_2-2*e0rc21*e1rc12*e1rc21+2*e0rc22*e1rc11*e1rc21-e0rc12*e1rc20_2-2*e0rc20*e1rc12*e1rc20+2*e0rc22*e1rc10*e1rc20+3*e0rc12*e1rc12_2+2*e0rc11*e1rc11*e1rc12+2*e0rc10*e1rc10*e1rc12+2*e0rc02*e1rc02*e1rc12-2*e0rc01*e1rc01*e1rc12-2*e0rc00*e1rc00*e1rc12+e0rc12*e1rc11_2+2*e0rc01*e1rc02*e1rc11+2*e0rc02*e1rc01*e1rc11+e0rc12*e1rc10_2+2*e0rc00*e1rc02*e1rc10+2*e0rc02*e1rc00*e1rc10+e0rc12*e1rc02_2+2*e0rc11*e1rc01*e1rc02+2*e0rc10*e1rc00*e1rc02-e0rc12*e1rc01_2-e0rc12*e1rc00_2;
    AMatrix(5, 4) = +2*e0rc12*e1rc22*e2rc22+2*e0rc11*e1rc21*e2rc22+2*e0rc10*e1rc20*e2rc22+2*e0rc22*e1rc12*e2rc22+2*e0rc21*e1rc11*e2rc22+2*e0rc20*e1rc10*e2rc22+2*e0rc11*e1rc22*e2rc21-2*e0rc12*e1rc21*e2rc21-2*e0rc21*e1rc12*e2rc21+2*e0rc22*e1rc11*e2rc21+2*e0rc10*e1rc22*e2rc20-2*e0rc12*e1rc20*e2rc20-2*e0rc20*e1rc12*e2rc20+2*e0rc22*e1rc10*e2rc20+2*e0rc22*e1rc22*e2rc12-2*e0rc21*e1rc21*e2rc12-2*e0rc20*e1rc20*e2rc12+6*e0rc12*e1rc12*e2rc12+2*e0rc11*e1rc11*e2rc12+2*e0rc10*e1rc10*e2rc12+2*e0rc02*e1rc02*e2rc12-2*e0rc01*e1rc01*e2rc12-2*e0rc00*e1rc00*e2rc12+2*e0rc21*e1rc22*e2rc11+2*e0rc22*e1rc21*e2rc11+2*e0rc11*e1rc12*e2rc11+2*e0rc12*e1rc11*e2rc11+2*e0rc01*e1rc02*e2rc11+2*e0rc02*e1rc01*e2rc11+2*e0rc20*e1rc22*e2rc10+2*e0rc22*e1rc20*e2rc10+2*e0rc10*e1rc12*e2rc10+2*e0rc12*e1rc10*e2rc10+2*e0rc00*e1rc02*e2rc10+2*e0rc02*e1rc00*e2rc10+2*e0rc02*e1rc12*e2rc02+2*e0rc01*e1rc11*e2rc02+2*e0rc00*e1rc10*e2rc02+2*e0rc12*e1rc02*e2rc02+2*e0rc11*e1rc01*e2rc02+2*e0rc10*e1rc00*e2rc02-2*e0rc01*e1rc12*e2rc01+2*e0rc02*e1rc11*e2rc01+2*e0rc11*e1rc02*e2rc01-2*e0rc12*e1rc01*e2rc01-2*e0rc00*e1rc12*e2rc00+2*e0rc02*e1rc10*e2rc00+2*e0rc10*e1rc02*e2rc00-2*e0rc12*e1rc00*e2rc00;
    AMatrix(5, 5) = +e0rc12*e2rc22_2+2*e0rc11*e2rc21*e2rc22+2*e0rc10*e2rc20*e2rc22+2*e0rc22*e2rc12*e2rc22+2*e0rc21*e2rc11*e2rc22+2*e0rc20*e2rc10*e2rc22-e0rc12*e2rc21_2-2*e0rc21*e2rc12*e2rc21+2*e0rc22*e2rc11*e2rc21-e0rc12*e2rc20_2-2*e0rc20*e2rc12*e2rc20+2*e0rc22*e2rc10*e2rc20+3*e0rc12*e2rc12_2+2*e0rc11*e2rc11*e2rc12+2*e0rc10*e2rc10*e2rc12+2*e0rc02*e2rc02*e2rc12-2*e0rc01*e2rc01*e2rc12-2*e0rc00*e2rc00*e2rc12+e0rc12*e2rc11_2+2*e0rc01*e2rc02*e2rc11+2*e0rc02*e2rc01*e2rc11+e0rc12*e2rc10_2+2*e0rc00*e2rc02*e2rc10+2*e0rc02*e2rc00*e2rc10+e0rc12*e2rc02_2+2*e0rc11*e2rc01*e2rc02+2*e0rc10*e2rc00*e2rc02-e0rc12*e2rc01_2-e0rc12*e2rc00_2;
    AMatrix(5, 6) = +e1rc12*e1rc22_2+2*e1rc11*e1rc21*e1rc22+2*e1rc10*e1rc20*e1rc22-e1rc12*e1rc21_2-e1rc12*e1rc20_2+e1rc12_3+e1rc11_2*e1rc12+e1rc10_2*e1rc12+e1rc02_2*e1rc12-e1rc01_2*e1rc12-e1rc00_2*e1rc12+2*e1rc01*e1rc02*e1rc11+2*e1rc00*e1rc02*e1rc10;
    AMatrix(5, 7) = +2*e1rc12*e1rc22*e2rc22+2*e1rc11*e1rc21*e2rc22+2*e1rc10*e1rc20*e2rc22+2*e1rc11*e1rc22*e2rc21-2*e1rc12*e1rc21*e2rc21+2*e1rc10*e1rc22*e2rc20-2*e1rc12*e1rc20*e2rc20+e1rc22_2*e2rc12-e1rc21_2*e2rc12-e1rc20_2*e2rc12+3*e1rc12_2*e2rc12+e1rc11_2*e2rc12+e1rc10_2*e2rc12+e1rc02_2*e2rc12-e1rc01_2*e2rc12-e1rc00_2*e2rc12+2*e1rc21*e1rc22*e2rc11+2*e1rc11*e1rc12*e2rc11+2*e1rc01*e1rc02*e2rc11+2*e1rc20*e1rc22*e2rc10+2*e1rc10*e1rc12*e2rc10+2*e1rc00*e1rc02*e2rc10+2*e1rc02*e1rc12*e2rc02+2*e1rc01*e1rc11*e2rc02+2*e1rc00*e1rc10*e2rc02-2*e1rc01*e1rc12*e2rc01+2*e1rc02*e1rc11*e2rc01-2*e1rc00*e1rc12*e2rc00+2*e1rc02*e1rc10*e2rc00;
    AMatrix(5, 8) = +e1rc12*e2rc22_2+2*e1rc11*e2rc21*e2rc22+2*e1rc10*e2rc20*e2rc22+2*e1rc22*e2rc12*e2rc22+2*e1rc21*e2rc11*e2rc22+2*e1rc20*e2rc10*e2rc22-e1rc12*e2rc21_2-2*e1rc21*e2rc12*e2rc21+2*e1rc22*e2rc11*e2rc21-e1rc12*e2rc20_2-2*e1rc20*e2rc12*e2rc20+2*e1rc22*e2rc10*e2rc20+3*e1rc12*e2rc12_2+2*e1rc11*e2rc11*e2rc12+2*e1rc10*e2rc10*e2rc12+2*e1rc02*e2rc02*e2rc12-2*e1rc01*e2rc01*e2rc12-2*e1rc00*e2rc00*e2rc12+e1rc12*e2rc11_2+2*e1rc01*e2rc02*e2rc11+2*e1rc02*e2rc01*e2rc11+e1rc12*e2rc10_2+2*e1rc00*e2rc02*e2rc10+2*e1rc02*e2rc00*e2rc10+e1rc12*e2rc02_2+2*e1rc11*e2rc01*e2rc02+2*e1rc10*e2rc00*e2rc02-e1rc12*e2rc01_2-e1rc12*e2rc00_2;
    AMatrix(5, 9) = e2rc12*e2rc22_2+2*e2rc11*e2rc21*e2rc22+2*e2rc10*e2rc20*e2rc22-e2rc12*e2rc21_2-e2rc12*e2rc20_2+e2rc12_3+e2rc11_2*e2rc12+e2rc10_2*e2rc12+e2rc02_2*e2rc12-e2rc01_2*e2rc12-e2rc00_2*e2rc12+2*e2rc01*e2rc02*e2rc11+2*e2rc00*e2rc02*e2rc10;
    AMatrix(5, 10) = +2*e0rc12*e0rc22*e3rc22+2*e0rc11*e0rc21*e3rc22+2*e0rc10*e0rc20*e3rc22+2*e0rc11*e0rc22*e3rc21-2*e0rc12*e0rc21*e3rc21+2*e0rc10*e0rc22*e3rc20-2*e0rc12*e0rc20*e3rc20+e0rc22_2*e3rc12-e0rc21_2*e3rc12-e0rc20_2*e3rc12+3*e0rc12_2*e3rc12+e0rc11_2*e3rc12+e0rc10_2*e3rc12+e0rc02_2*e3rc12-e0rc01_2*e3rc12-e0rc00_2*e3rc12+2*e0rc21*e0rc22*e3rc11+2*e0rc11*e0rc12*e3rc11+2*e0rc01*e0rc02*e3rc11+2*e0rc20*e0rc22*e3rc10+2*e0rc10*e0rc12*e3rc10+2*e0rc00*e0rc02*e3rc10+2*e0rc02*e0rc12*e3rc02+2*e0rc01*e0rc11*e3rc02+2*e0rc00*e0rc10*e3rc02-2*e0rc01*e0rc12*e3rc01+2*e0rc02*e0rc11*e3rc01-2*e0rc00*e0rc12*e3rc00+2*e0rc02*e0rc10*e3rc00;
    AMatrix(5, 11) = +2*e0rc12*e1rc22*e3rc22+2*e0rc11*e1rc21*e3rc22+2*e0rc10*e1rc20*e3rc22+2*e0rc22*e1rc12*e3rc22+2*e0rc21*e1rc11*e3rc22+2*e0rc20*e1rc10*e3rc22+2*e0rc11*e1rc22*e3rc21-2*e0rc12*e1rc21*e3rc21-2*e0rc21*e1rc12*e3rc21+2*e0rc22*e1rc11*e3rc21+2*e0rc10*e1rc22*e3rc20-2*e0rc12*e1rc20*e3rc20-2*e0rc20*e1rc12*e3rc20+2*e0rc22*e1rc10*e3rc20+2*e0rc22*e1rc22*e3rc12-2*e0rc21*e1rc21*e3rc12-2*e0rc20*e1rc20*e3rc12+6*e0rc12*e1rc12*e3rc12+2*e0rc11*e1rc11*e3rc12+2*e0rc10*e1rc10*e3rc12+2*e0rc02*e1rc02*e3rc12-2*e0rc01*e1rc01*e3rc12-2*e0rc00*e1rc00*e3rc12+2*e0rc21*e1rc22*e3rc11+2*e0rc22*e1rc21*e3rc11+2*e0rc11*e1rc12*e3rc11+2*e0rc12*e1rc11*e3rc11+2*e0rc01*e1rc02*e3rc11+2*e0rc02*e1rc01*e3rc11+2*e0rc20*e1rc22*e3rc10+2*e0rc22*e1rc20*e3rc10+2*e0rc10*e1rc12*e3rc10+2*e0rc12*e1rc10*e3rc10+2*e0rc00*e1rc02*e3rc10+2*e0rc02*e1rc00*e3rc10+2*e0rc02*e1rc12*e3rc02+2*e0rc01*e1rc11*e3rc02+2*e0rc00*e1rc10*e3rc02+2*e0rc12*e1rc02*e3rc02+2*e0rc11*e1rc01*e3rc02+2*e0rc10*e1rc00*e3rc02-2*e0rc01*e1rc12*e3rc01+2*e0rc02*e1rc11*e3rc01+2*e0rc11*e1rc02*e3rc01-2*e0rc12*e1rc01*e3rc01-2*e0rc00*e1rc12*e3rc00+2*e0rc02*e1rc10*e3rc00+2*e0rc10*e1rc02*e3rc00-2*e0rc12*e1rc00*e3rc00;
    AMatrix(5, 12) = +2*e0rc12*e2rc22*e3rc22+2*e0rc11*e2rc21*e3rc22+2*e0rc10*e2rc20*e3rc22+2*e0rc22*e2rc12*e3rc22+2*e0rc21*e2rc11*e3rc22+2*e0rc20*e2rc10*e3rc22+2*e0rc11*e2rc22*e3rc21-2*e0rc12*e2rc21*e3rc21-2*e0rc21*e2rc12*e3rc21+2*e0rc22*e2rc11*e3rc21+2*e0rc10*e2rc22*e3rc20-2*e0rc12*e2rc20*e3rc20-2*e0rc20*e2rc12*e3rc20+2*e0rc22*e2rc10*e3rc20+2*e0rc22*e2rc22*e3rc12-2*e0rc21*e2rc21*e3rc12-2*e0rc20*e2rc20*e3rc12+6*e0rc12*e2rc12*e3rc12+2*e0rc11*e2rc11*e3rc12+2*e0rc10*e2rc10*e3rc12+2*e0rc02*e2rc02*e3rc12-2*e0rc01*e2rc01*e3rc12-2*e0rc00*e2rc00*e3rc12+2*e0rc21*e2rc22*e3rc11+2*e0rc22*e2rc21*e3rc11+2*e0rc11*e2rc12*e3rc11+2*e0rc12*e2rc11*e3rc11+2*e0rc01*e2rc02*e3rc11+2*e0rc02*e2rc01*e3rc11+2*e0rc20*e2rc22*e3rc10+2*e0rc22*e2rc20*e3rc10+2*e0rc10*e2rc12*e3rc10+2*e0rc12*e2rc10*e3rc10+2*e0rc00*e2rc02*e3rc10+2*e0rc02*e2rc00*e3rc10+2*e0rc02*e2rc12*e3rc02+2*e0rc01*e2rc11*e3rc02+2*e0rc00*e2rc10*e3rc02+2*e0rc12*e2rc02*e3rc02+2*e0rc11*e2rc01*e3rc02+2*e0rc10*e2rc00*e3rc02-2*e0rc01*e2rc12*e3rc01+2*e0rc02*e2rc11*e3rc01+2*e0rc11*e2rc02*e3rc01-2*e0rc12*e2rc01*e3rc01-2*e0rc00*e2rc12*e3rc00+2*e0rc02*e2rc10*e3rc00+2*e0rc10*e2rc02*e3rc00-2*e0rc12*e2rc00*e3rc00;
    AMatrix(5, 13) = +2*e1rc12*e1rc22*e3rc22+2*e1rc11*e1rc21*e3rc22+2*e1rc10*e1rc20*e3rc22+2*e1rc11*e1rc22*e3rc21-2*e1rc12*e1rc21*e3rc21+2*e1rc10*e1rc22*e3rc20-2*e1rc12*e1rc20*e3rc20+e1rc22_2*e3rc12-e1rc21_2*e3rc12-e1rc20_2*e3rc12+3*e1rc12_2*e3rc12+e1rc11_2*e3rc12+e1rc10_2*e3rc12+e1rc02_2*e3rc12-e1rc01_2*e3rc12-e1rc00_2*e3rc12+2*e1rc21*e1rc22*e3rc11+2*e1rc11*e1rc12*e3rc11+2*e1rc01*e1rc02*e3rc11+2*e1rc20*e1rc22*e3rc10+2*e1rc10*e1rc12*e3rc10+2*e1rc00*e1rc02*e3rc10+2*e1rc02*e1rc12*e3rc02+2*e1rc01*e1rc11*e3rc02+2*e1rc00*e1rc10*e3rc02-2*e1rc01*e1rc12*e3rc01+2*e1rc02*e1rc11*e3rc01-2*e1rc00*e1rc12*e3rc00+2*e1rc02*e1rc10*e3rc00;
    AMatrix(5, 14) = +2*e1rc12*e2rc22*e3rc22+2*e1rc11*e2rc21*e3rc22+2*e1rc10*e2rc20*e3rc22+2*e1rc22*e2rc12*e3rc22+2*e1rc21*e2rc11*e3rc22+2*e1rc20*e2rc10*e3rc22+2*e1rc11*e2rc22*e3rc21-2*e1rc12*e2rc21*e3rc21-2*e1rc21*e2rc12*e3rc21+2*e1rc22*e2rc11*e3rc21+2*e1rc10*e2rc22*e3rc20-2*e1rc12*e2rc20*e3rc20-2*e1rc20*e2rc12*e3rc20+2*e1rc22*e2rc10*e3rc20+2*e1rc22*e2rc22*e3rc12-2*e1rc21*e2rc21*e3rc12-2*e1rc20*e2rc20*e3rc12+6*e1rc12*e2rc12*e3rc12+2*e1rc11*e2rc11*e3rc12+2*e1rc10*e2rc10*e3rc12+2*e1rc02*e2rc02*e3rc12-2*e1rc01*e2rc01*e3rc12-2*e1rc00*e2rc00*e3rc12+2*e1rc21*e2rc22*e3rc11+2*e1rc22*e2rc21*e3rc11+2*e1rc11*e2rc12*e3rc11+2*e1rc12*e2rc11*e3rc11+2*e1rc01*e2rc02*e3rc11+2*e1rc02*e2rc01*e3rc11+2*e1rc20*e2rc22*e3rc10+2*e1rc22*e2rc20*e3rc10+2*e1rc10*e2rc12*e3rc10+2*e1rc12*e2rc10*e3rc10+2*e1rc00*e2rc02*e3rc10+2*e1rc02*e2rc00*e3rc10+2*e1rc02*e2rc12*e3rc02+2*e1rc01*e2rc11*e3rc02+2*e1rc00*e2rc10*e3rc02+2*e1rc12*e2rc02*e3rc02+2*e1rc11*e2rc01*e3rc02+2*e1rc10*e2rc00*e3rc02-2*e1rc01*e2rc12*e3rc01+2*e1rc02*e2rc11*e3rc01+2*e1rc11*e2rc02*e3rc01-2*e1rc12*e2rc01*e3rc01-2*e1rc00*e2rc12*e3rc00+2*e1rc02*e2rc10*e3rc00+2*e1rc10*e2rc02*e3rc00-2*e1rc12*e2rc00*e3rc00;
    AMatrix(5, 15) = +2*e2rc12*e2rc22*e3rc22+2*e2rc11*e2rc21*e3rc22+2*e2rc10*e2rc20*e3rc22+2*e2rc11*e2rc22*e3rc21-2*e2rc12*e2rc21*e3rc21+2*e2rc10*e2rc22*e3rc20-2*e2rc12*e2rc20*e3rc20+e2rc22_2*e3rc12-e2rc21_2*e3rc12-e2rc20_2*e3rc12+3*e2rc12_2*e3rc12+e2rc11_2*e3rc12+e2rc10_2*e3rc12+e2rc02_2*e3rc12-e2rc01_2*e3rc12-e2rc00_2*e3rc12+2*e2rc21*e2rc22*e3rc11+2*e2rc11*e2rc12*e3rc11+2*e2rc01*e2rc02*e3rc11+2*e2rc20*e2rc22*e3rc10+2*e2rc10*e2rc12*e3rc10+2*e2rc00*e2rc02*e3rc10+2*e2rc02*e2rc12*e3rc02+2*e2rc01*e2rc11*e3rc02+2*e2rc00*e2rc10*e3rc02-2*e2rc01*e2rc12*e3rc01+2*e2rc02*e2rc11*e3rc01-2*e2rc00*e2rc12*e3rc00+2*e2rc02*e2rc10*e3rc00;
    AMatrix(5, 16) = +e0rc12*e3rc22_2+2*e0rc11*e3rc21*e3rc22+2*e0rc10*e3rc20*e3rc22+2*e0rc22*e3rc12*e3rc22+2*e0rc21*e3rc11*e3rc22+2*e0rc20*e3rc10*e3rc22-e0rc12*e3rc21_2-2*e0rc21*e3rc12*e3rc21+2*e0rc22*e3rc11*e3rc21-e0rc12*e3rc20_2-2*e0rc20*e3rc12*e3rc20+2*e0rc22*e3rc10*e3rc20+3*e0rc12*e3rc12_2+2*e0rc11*e3rc11*e3rc12+2*e0rc10*e3rc10*e3rc12+2*e0rc02*e3rc02*e3rc12-2*e0rc01*e3rc01*e3rc12-2*e0rc00*e3rc00*e3rc12+e0rc12*e3rc11_2+2*e0rc01*e3rc02*e3rc11+2*e0rc02*e3rc01*e3rc11+e0rc12*e3rc10_2+2*e0rc00*e3rc02*e3rc10+2*e0rc02*e3rc00*e3rc10+e0rc12*e3rc02_2+2*e0rc11*e3rc01*e3rc02+2*e0rc10*e3rc00*e3rc02-e0rc12*e3rc01_2-e0rc12*e3rc00_2;
    AMatrix(5, 17) = +e1rc12*e3rc22_2+2*e1rc11*e3rc21*e3rc22+2*e1rc10*e3rc20*e3rc22+2*e1rc22*e3rc12*e3rc22+2*e1rc21*e3rc11*e3rc22+2*e1rc20*e3rc10*e3rc22-e1rc12*e3rc21_2-2*e1rc21*e3rc12*e3rc21+2*e1rc22*e3rc11*e3rc21-e1rc12*e3rc20_2-2*e1rc20*e3rc12*e3rc20+2*e1rc22*e3rc10*e3rc20+3*e1rc12*e3rc12_2+2*e1rc11*e3rc11*e3rc12+2*e1rc10*e3rc10*e3rc12+2*e1rc02*e3rc02*e3rc12-2*e1rc01*e3rc01*e3rc12-2*e1rc00*e3rc00*e3rc12+e1rc12*e3rc11_2+2*e1rc01*e3rc02*e3rc11+2*e1rc02*e3rc01*e3rc11+e1rc12*e3rc10_2+2*e1rc00*e3rc02*e3rc10+2*e1rc02*e3rc00*e3rc10+e1rc12*e3rc02_2+2*e1rc11*e3rc01*e3rc02+2*e1rc10*e3rc00*e3rc02-e1rc12*e3rc01_2-e1rc12*e3rc00_2;
    AMatrix(5, 18) = +e2rc12*e3rc22_2+2*e2rc11*e3rc21*e3rc22+2*e2rc10*e3rc20*e3rc22+2*e2rc22*e3rc12*e3rc22+2*e2rc21*e3rc11*e3rc22+2*e2rc20*e3rc10*e3rc22-e2rc12*e3rc21_2-2*e2rc21*e3rc12*e3rc21+2*e2rc22*e3rc11*e3rc21-e2rc12*e3rc20_2-2*e2rc20*e3rc12*e3rc20+2*e2rc22*e3rc10*e3rc20+3*e2rc12*e3rc12_2+2*e2rc11*e3rc11*e3rc12+2*e2rc10*e3rc10*e3rc12+2*e2rc02*e3rc02*e3rc12-2*e2rc01*e3rc01*e3rc12-2*e2rc00*e3rc00*e3rc12+e2rc12*e3rc11_2+2*e2rc01*e3rc02*e3rc11+2*e2rc02*e3rc01*e3rc11+e2rc12*e3rc10_2+2*e2rc00*e3rc02*e3rc10+2*e2rc02*e3rc00*e3rc10+e2rc12*e3rc02_2+2*e2rc11*e3rc01*e3rc02+2*e2rc10*e3rc00*e3rc02-e2rc12*e3rc01_2-e2rc12*e3rc00_2;
    AMatrix(5, 19) = +e3rc12*e3rc22_2+2*e3rc11*e3rc21*e3rc22+2*e3rc10*e3rc20*e3rc22-e3rc12*e3rc21_2-e3rc12*e3rc20_2+e3rc12_3+e3rc11_2*e3rc12+e3rc10_2*e3rc12+e3rc02_2*e3rc12-e3rc01_2*e3rc12-e3rc00_2*e3rc12+2*e3rc01*e3rc02*e3rc11+2*e3rc00*e3rc02*e3rc10;
    AMatrix(6, 0) = +e0rc20*e0rc22_2+2*e0rc10*e0rc12*e0rc22+2*e0rc00*e0rc02*e0rc22+e0rc20*e0rc21_2+2*e0rc10*e0rc11*e0rc21+2*e0rc00*e0rc01*e0rc21+e0rc20_3-e0rc12_2*e0rc20-e0rc11_2*e0rc20+e0rc10_2*e0rc20-e0rc02_2*e0rc20-e0rc01_2*e0rc20+e0rc00_2*e0rc20;
    AMatrix(6, 1) = +2*e0rc20*e0rc22*e1rc22+2*e0rc10*e0rc12*e1rc22+2*e0rc00*e0rc02*e1rc22+2*e0rc20*e0rc21*e1rc21+2*e0rc10*e0rc11*e1rc21+2*e0rc00*e0rc01*e1rc21+e0rc22_2*e1rc20+e0rc21_2*e1rc20+3*e0rc20_2*e1rc20-e0rc12_2*e1rc20-e0rc11_2*e1rc20+e0rc10_2*e1rc20-e0rc02_2*e1rc20-e0rc01_2*e1rc20+e0rc00_2*e1rc20+2*e0rc10*e0rc22*e1rc12-2*e0rc12*e0rc20*e1rc12+2*e0rc10*e0rc21*e1rc11-2*e0rc11*e0rc20*e1rc11+2*e0rc12*e0rc22*e1rc10+2*e0rc11*e0rc21*e1rc10+2*e0rc10*e0rc20*e1rc10+2*e0rc00*e0rc22*e1rc02-2*e0rc02*e0rc20*e1rc02+2*e0rc00*e0rc21*e1rc01-2*e0rc01*e0rc20*e1rc01+2*e0rc02*e0rc22*e1rc00+2*e0rc01*e0rc21*e1rc00+2*e0rc00*e0rc20*e1rc00;
    AMatrix(6, 2) = +2*e0rc20*e0rc22*e2rc22+2*e0rc10*e0rc12*e2rc22+2*e0rc00*e0rc02*e2rc22+2*e0rc20*e0rc21*e2rc21+2*e0rc10*e0rc11*e2rc21+2*e0rc00*e0rc01*e2rc21+e0rc22_2*e2rc20+e0rc21_2*e2rc20+3*e0rc20_2*e2rc20-e0rc12_2*e2rc20-e0rc11_2*e2rc20+e0rc10_2*e2rc20-e0rc02_2*e2rc20-e0rc01_2*e2rc20+e0rc00_2*e2rc20+2*e0rc10*e0rc22*e2rc12-2*e0rc12*e0rc20*e2rc12+2*e0rc10*e0rc21*e2rc11-2*e0rc11*e0rc20*e2rc11+2*e0rc12*e0rc22*e2rc10+2*e0rc11*e0rc21*e2rc10+2*e0rc10*e0rc20*e2rc10+2*e0rc00*e0rc22*e2rc02-2*e0rc02*e0rc20*e2rc02+2*e0rc00*e0rc21*e2rc01-2*e0rc01*e0rc20*e2rc01+2*e0rc02*e0rc22*e2rc00+2*e0rc01*e0rc21*e2rc00+2*e0rc00*e0rc20*e2rc00;
    AMatrix(6, 3) = +e0rc20*e1rc22_2+2*e0rc22*e1rc20*e1rc22+2*e0rc10*e1rc12*e1rc22+2*e0rc12*e1rc10*e1rc22+2*e0rc00*e1rc02*e1rc22+2*e0rc02*e1rc00*e1rc22+e0rc20*e1rc21_2+2*e0rc21*e1rc20*e1rc21+2*e0rc10*e1rc11*e1rc21+2*e0rc11*e1rc10*e1rc21+2*e0rc00*e1rc01*e1rc21+2*e0rc01*e1rc00*e1rc21+3*e0rc20*e1rc20_2-2*e0rc12*e1rc12*e1rc20-2*e0rc11*e1rc11*e1rc20+2*e0rc10*e1rc10*e1rc20-2*e0rc02*e1rc02*e1rc20-2*e0rc01*e1rc01*e1rc20+2*e0rc00*e1rc00*e1rc20-e0rc20*e1rc12_2+2*e0rc22*e1rc10*e1rc12-e0rc20*e1rc11_2+2*e0rc21*e1rc10*e1rc11+e0rc20*e1rc10_2-e0rc20*e1rc02_2+2*e0rc22*e1rc00*e1rc02-e0rc20*e1rc01_2+2*e0rc21*e1rc00*e1rc01+e0rc20*e1rc00_2;
    AMatrix(6, 4) = +2*e0rc20*e1rc22*e2rc22+2*e0rc22*e1rc20*e2rc22+2*e0rc10*e1rc12*e2rc22+2*e0rc12*e1rc10*e2rc22+2*e0rc00*e1rc02*e2rc22+2*e0rc02*e1rc00*e2rc22+2*e0rc20*e1rc21*e2rc21+2*e0rc21*e1rc20*e2rc21+2*e0rc10*e1rc11*e2rc21+2*e0rc11*e1rc10*e2rc21+2*e0rc00*e1rc01*e2rc21+2*e0rc01*e1rc00*e2rc21+2*e0rc22*e1rc22*e2rc20+2*e0rc21*e1rc21*e2rc20+6*e0rc20*e1rc20*e2rc20-2*e0rc12*e1rc12*e2rc20-2*e0rc11*e1rc11*e2rc20+2*e0rc10*e1rc10*e2rc20-2*e0rc02*e1rc02*e2rc20-2*e0rc01*e1rc01*e2rc20+2*e0rc00*e1rc00*e2rc20+2*e0rc10*e1rc22*e2rc12-2*e0rc12*e1rc20*e2rc12-2*e0rc20*e1rc12*e2rc12+2*e0rc22*e1rc10*e2rc12+2*e0rc10*e1rc21*e2rc11-2*e0rc11*e1rc20*e2rc11-2*e0rc20*e1rc11*e2rc11+2*e0rc21*e1rc10*e2rc11+2*e0rc12*e1rc22*e2rc10+2*e0rc11*e1rc21*e2rc10+2*e0rc10*e1rc20*e2rc10+2*e0rc22*e1rc12*e2rc10+2*e0rc21*e1rc11*e2rc10+2*e0rc20*e1rc10*e2rc10+2*e0rc00*e1rc22*e2rc02-2*e0rc02*e1rc20*e2rc02-2*e0rc20*e1rc02*e2rc02+2*e0rc22*e1rc00*e2rc02+2*e0rc00*e1rc21*e2rc01-2*e0rc01*e1rc20*e2rc01-2*e0rc20*e1rc01*e2rc01+2*e0rc21*e1rc00*e2rc01+2*e0rc02*e1rc22*e2rc00+2*e0rc01*e1rc21*e2rc00+2*e0rc00*e1rc20*e2rc00+2*e0rc22*e1rc02*e2rc00+2*e0rc21*e1rc01*e2rc00+2*e0rc20*e1rc00*e2rc00;
    AMatrix(6, 5) = +e0rc20*e2rc22_2+2*e0rc22*e2rc20*e2rc22+2*e0rc10*e2rc12*e2rc22+2*e0rc12*e2rc10*e2rc22+2*e0rc00*e2rc02*e2rc22+2*e0rc02*e2rc00*e2rc22+e0rc20*e2rc21_2+2*e0rc21*e2rc20*e2rc21+2*e0rc10*e2rc11*e2rc21+2*e0rc11*e2rc10*e2rc21+2*e0rc00*e2rc01*e2rc21+2*e0rc01*e2rc00*e2rc21+3*e0rc20*e2rc20_2-2*e0rc12*e2rc12*e2rc20-2*e0rc11*e2rc11*e2rc20+2*e0rc10*e2rc10*e2rc20-2*e0rc02*e2rc02*e2rc20-2*e0rc01*e2rc01*e2rc20+2*e0rc00*e2rc00*e2rc20-e0rc20*e2rc12_2+2*e0rc22*e2rc10*e2rc12-e0rc20*e2rc11_2+2*e0rc21*e2rc10*e2rc11+e0rc20*e2rc10_2-e0rc20*e2rc02_2+2*e0rc22*e2rc00*e2rc02-e0rc20*e2rc01_2+2*e0rc21*e2rc00*e2rc01+e0rc20*e2rc00_2;
    AMatrix(6, 6) = +e1rc20*e1rc22_2+2*e1rc10*e1rc12*e1rc22+2*e1rc00*e1rc02*e1rc22+e1rc20*e1rc21_2+2*e1rc10*e1rc11*e1rc21+2*e1rc00*e1rc01*e1rc21+e1rc20_3-e1rc12_2*e1rc20-e1rc11_2*e1rc20+e1rc10_2*e1rc20-e1rc02_2*e1rc20-e1rc01_2*e1rc20+e1rc00_2*e1rc20;
    AMatrix(6, 7) = +2*e1rc20*e1rc22*e2rc22+2*e1rc10*e1rc12*e2rc22+2*e1rc00*e1rc02*e2rc22+2*e1rc20*e1rc21*e2rc21+2*e1rc10*e1rc11*e2rc21+2*e1rc00*e1rc01*e2rc21+e1rc22_2*e2rc20+e1rc21_2*e2rc20+3*e1rc20_2*e2rc20-e1rc12_2*e2rc20-e1rc11_2*e2rc20+e1rc10_2*e2rc20-e1rc02_2*e2rc20-e1rc01_2*e2rc20+e1rc00_2*e2rc20+2*e1rc10*e1rc22*e2rc12-2*e1rc12*e1rc20*e2rc12+2*e1rc10*e1rc21*e2rc11-2*e1rc11*e1rc20*e2rc11+2*e1rc12*e1rc22*e2rc10+2*e1rc11*e1rc21*e2rc10+2*e1rc10*e1rc20*e2rc10+2*e1rc00*e1rc22*e2rc02-2*e1rc02*e1rc20*e2rc02+2*e1rc00*e1rc21*e2rc01-2*e1rc01*e1rc20*e2rc01+2*e1rc02*e1rc22*e2rc00+2*e1rc01*e1rc21*e2rc00+2*e1rc00*e1rc20*e2rc00;
    AMatrix(6, 8) = +e1rc20*e2rc22_2+2*e1rc22*e2rc20*e2rc22+2*e1rc10*e2rc12*e2rc22+2*e1rc12*e2rc10*e2rc22+2*e1rc00*e2rc02*e2rc22+2*e1rc02*e2rc00*e2rc22+e1rc20*e2rc21_2+2*e1rc21*e2rc20*e2rc21+2*e1rc10*e2rc11*e2rc21+2*e1rc11*e2rc10*e2rc21+2*e1rc00*e2rc01*e2rc21+2*e1rc01*e2rc00*e2rc21+3*e1rc20*e2rc20_2-2*e1rc12*e2rc12*e2rc20-2*e1rc11*e2rc11*e2rc20+2*e1rc10*e2rc10*e2rc20-2*e1rc02*e2rc02*e2rc20-2*e1rc01*e2rc01*e2rc20+2*e1rc00*e2rc00*e2rc20-e1rc20*e2rc12_2+2*e1rc22*e2rc10*e2rc12-e1rc20*e2rc11_2+2*e1rc21*e2rc10*e2rc11+e1rc20*e2rc10_2-e1rc20*e2rc02_2+2*e1rc22*e2rc00*e2rc02-e1rc20*e2rc01_2+2*e1rc21*e2rc00*e2rc01+e1rc20*e2rc00_2;
    AMatrix(6, 9) = e2rc20*e2rc22_2+2*e2rc10*e2rc12*e2rc22+2*e2rc00*e2rc02*e2rc22+e2rc20*e2rc21_2+2*e2rc10*e2rc11*e2rc21+2*e2rc00*e2rc01*e2rc21+e2rc20_3-e2rc12_2*e2rc20-e2rc11_2*e2rc20+e2rc10_2*e2rc20-e2rc02_2*e2rc20-e2rc01_2*e2rc20+e2rc00_2*e2rc20;
    AMatrix(6, 10) = +2*e0rc20*e0rc22*e3rc22+2*e0rc10*e0rc12*e3rc22+2*e0rc00*e0rc02*e3rc22+2*e0rc20*e0rc21*e3rc21+2*e0rc10*e0rc11*e3rc21+2*e0rc00*e0rc01*e3rc21+e0rc22_2*e3rc20+e0rc21_2*e3rc20+3*e0rc20_2*e3rc20-e0rc12_2*e3rc20-e0rc11_2*e3rc20+e0rc10_2*e3rc20-e0rc02_2*e3rc20-e0rc01_2*e3rc20+e0rc00_2*e3rc20+2*e0rc10*e0rc22*e3rc12-2*e0rc12*e0rc20*e3rc12+2*e0rc10*e0rc21*e3rc11-2*e0rc11*e0rc20*e3rc11+2*e0rc12*e0rc22*e3rc10+2*e0rc11*e0rc21*e3rc10+2*e0rc10*e0rc20*e3rc10+2*e0rc00*e0rc22*e3rc02-2*e0rc02*e0rc20*e3rc02+2*e0rc00*e0rc21*e3rc01-2*e0rc01*e0rc20*e3rc01+2*e0rc02*e0rc22*e3rc00+2*e0rc01*e0rc21*e3rc00+2*e0rc00*e0rc20*e3rc00;
    AMatrix(6, 11) = +2*e0rc20*e1rc22*e3rc22+2*e0rc22*e1rc20*e3rc22+2*e0rc10*e1rc12*e3rc22+2*e0rc12*e1rc10*e3rc22+2*e0rc00*e1rc02*e3rc22+2*e0rc02*e1rc00*e3rc22+2*e0rc20*e1rc21*e3rc21+2*e0rc21*e1rc20*e3rc21+2*e0rc10*e1rc11*e3rc21+2*e0rc11*e1rc10*e3rc21+2*e0rc00*e1rc01*e3rc21+2*e0rc01*e1rc00*e3rc21+2*e0rc22*e1rc22*e3rc20+2*e0rc21*e1rc21*e3rc20+6*e0rc20*e1rc20*e3rc20-2*e0rc12*e1rc12*e3rc20-2*e0rc11*e1rc11*e3rc20+2*e0rc10*e1rc10*e3rc20-2*e0rc02*e1rc02*e3rc20-2*e0rc01*e1rc01*e3rc20+2*e0rc00*e1rc00*e3rc20+2*e0rc10*e1rc22*e3rc12-2*e0rc12*e1rc20*e3rc12-2*e0rc20*e1rc12*e3rc12+2*e0rc22*e1rc10*e3rc12+2*e0rc10*e1rc21*e3rc11-2*e0rc11*e1rc20*e3rc11-2*e0rc20*e1rc11*e3rc11+2*e0rc21*e1rc10*e3rc11+2*e0rc12*e1rc22*e3rc10+2*e0rc11*e1rc21*e3rc10+2*e0rc10*e1rc20*e3rc10+2*e0rc22*e1rc12*e3rc10+2*e0rc21*e1rc11*e3rc10+2*e0rc20*e1rc10*e3rc10+2*e0rc00*e1rc22*e3rc02-2*e0rc02*e1rc20*e3rc02-2*e0rc20*e1rc02*e3rc02+2*e0rc22*e1rc00*e3rc02+2*e0rc00*e1rc21*e3rc01-2*e0rc01*e1rc20*e3rc01-2*e0rc20*e1rc01*e3rc01+2*e0rc21*e1rc00*e3rc01+2*e0rc02*e1rc22*e3rc00+2*e0rc01*e1rc21*e3rc00+2*e0rc00*e1rc20*e3rc00+2*e0rc22*e1rc02*e3rc00+2*e0rc21*e1rc01*e3rc00+2*e0rc20*e1rc00*e3rc00;
    AMatrix(6, 12) = +2*e0rc20*e2rc22*e3rc22+2*e0rc22*e2rc20*e3rc22+2*e0rc10*e2rc12*e3rc22+2*e0rc12*e2rc10*e3rc22+2*e0rc00*e2rc02*e3rc22+2*e0rc02*e2rc00*e3rc22+2*e0rc20*e2rc21*e3rc21+2*e0rc21*e2rc20*e3rc21+2*e0rc10*e2rc11*e3rc21+2*e0rc11*e2rc10*e3rc21+2*e0rc00*e2rc01*e3rc21+2*e0rc01*e2rc00*e3rc21+2*e0rc22*e2rc22*e3rc20+2*e0rc21*e2rc21*e3rc20+6*e0rc20*e2rc20*e3rc20-2*e0rc12*e2rc12*e3rc20-2*e0rc11*e2rc11*e3rc20+2*e0rc10*e2rc10*e3rc20-2*e0rc02*e2rc02*e3rc20-2*e0rc01*e2rc01*e3rc20+2*e0rc00*e2rc00*e3rc20+2*e0rc10*e2rc22*e3rc12-2*e0rc12*e2rc20*e3rc12-2*e0rc20*e2rc12*e3rc12+2*e0rc22*e2rc10*e3rc12+2*e0rc10*e2rc21*e3rc11-2*e0rc11*e2rc20*e3rc11-2*e0rc20*e2rc11*e3rc11+2*e0rc21*e2rc10*e3rc11+2*e0rc12*e2rc22*e3rc10+2*e0rc11*e2rc21*e3rc10+2*e0rc10*e2rc20*e3rc10+2*e0rc22*e2rc12*e3rc10+2*e0rc21*e2rc11*e3rc10+2*e0rc20*e2rc10*e3rc10+2*e0rc00*e2rc22*e3rc02-2*e0rc02*e2rc20*e3rc02-2*e0rc20*e2rc02*e3rc02+2*e0rc22*e2rc00*e3rc02+2*e0rc00*e2rc21*e3rc01-2*e0rc01*e2rc20*e3rc01-2*e0rc20*e2rc01*e3rc01+2*e0rc21*e2rc00*e3rc01+2*e0rc02*e2rc22*e3rc00+2*e0rc01*e2rc21*e3rc00+2*e0rc00*e2rc20*e3rc00+2*e0rc22*e2rc02*e3rc00+2*e0rc21*e2rc01*e3rc00+2*e0rc20*e2rc00*e3rc00;
    AMatrix(6, 13) = +2*e1rc20*e1rc22*e3rc22+2*e1rc10*e1rc12*e3rc22+2*e1rc00*e1rc02*e3rc22+2*e1rc20*e1rc21*e3rc21+2*e1rc10*e1rc11*e3rc21+2*e1rc00*e1rc01*e3rc21+e1rc22_2*e3rc20+e1rc21_2*e3rc20+3*e1rc20_2*e3rc20-e1rc12_2*e3rc20-e1rc11_2*e3rc20+e1rc10_2*e3rc20-e1rc02_2*e3rc20-e1rc01_2*e3rc20+e1rc00_2*e3rc20+2*e1rc10*e1rc22*e3rc12-2*e1rc12*e1rc20*e3rc12+2*e1rc10*e1rc21*e3rc11-2*e1rc11*e1rc20*e3rc11+2*e1rc12*e1rc22*e3rc10+2*e1rc11*e1rc21*e3rc10+2*e1rc10*e1rc20*e3rc10+2*e1rc00*e1rc22*e3rc02-2*e1rc02*e1rc20*e3rc02+2*e1rc00*e1rc21*e3rc01-2*e1rc01*e1rc20*e3rc01+2*e1rc02*e1rc22*e3rc00+2*e1rc01*e1rc21*e3rc00+2*e1rc00*e1rc20*e3rc00;
    AMatrix(6, 14) = +2*e1rc20*e2rc22*e3rc22+2*e1rc22*e2rc20*e3rc22+2*e1rc10*e2rc12*e3rc22+2*e1rc12*e2rc10*e3rc22+2*e1rc00*e2rc02*e3rc22+2*e1rc02*e2rc00*e3rc22+2*e1rc20*e2rc21*e3rc21+2*e1rc21*e2rc20*e3rc21+2*e1rc10*e2rc11*e3rc21+2*e1rc11*e2rc10*e3rc21+2*e1rc00*e2rc01*e3rc21+2*e1rc01*e2rc00*e3rc21+2*e1rc22*e2rc22*e3rc20+2*e1rc21*e2rc21*e3rc20+6*e1rc20*e2rc20*e3rc20-2*e1rc12*e2rc12*e3rc20-2*e1rc11*e2rc11*e3rc20+2*e1rc10*e2rc10*e3rc20-2*e1rc02*e2rc02*e3rc20-2*e1rc01*e2rc01*e3rc20+2*e1rc00*e2rc00*e3rc20+2*e1rc10*e2rc22*e3rc12-2*e1rc12*e2rc20*e3rc12-2*e1rc20*e2rc12*e3rc12+2*e1rc22*e2rc10*e3rc12+2*e1rc10*e2rc21*e3rc11-2*e1rc11*e2rc20*e3rc11-2*e1rc20*e2rc11*e3rc11+2*e1rc21*e2rc10*e3rc11+2*e1rc12*e2rc22*e3rc10+2*e1rc11*e2rc21*e3rc10+2*e1rc10*e2rc20*e3rc10+2*e1rc22*e2rc12*e3rc10+2*e1rc21*e2rc11*e3rc10+2*e1rc20*e2rc10*e3rc10+2*e1rc00*e2rc22*e3rc02-2*e1rc02*e2rc20*e3rc02-2*e1rc20*e2rc02*e3rc02+2*e1rc22*e2rc00*e3rc02+2*e1rc00*e2rc21*e3rc01-2*e1rc01*e2rc20*e3rc01-2*e1rc20*e2rc01*e3rc01+2*e1rc21*e2rc00*e3rc01+2*e1rc02*e2rc22*e3rc00+2*e1rc01*e2rc21*e3rc00+2*e1rc00*e2rc20*e3rc00+2*e1rc22*e2rc02*e3rc00+2*e1rc21*e2rc01*e3rc00+2*e1rc20*e2rc00*e3rc00;
    AMatrix(6, 15) = +2*e2rc20*e2rc22*e3rc22+2*e2rc10*e2rc12*e3rc22+2*e2rc00*e2rc02*e3rc22+2*e2rc20*e2rc21*e3rc21+2*e2rc10*e2rc11*e3rc21+2*e2rc00*e2rc01*e3rc21+e2rc22_2*e3rc20+e2rc21_2*e3rc20+3*e2rc20_2*e3rc20-e2rc12_2*e3rc20-e2rc11_2*e3rc20+e2rc10_2*e3rc20-e2rc02_2*e3rc20-e2rc01_2*e3rc20+e2rc00_2*e3rc20+2*e2rc10*e2rc22*e3rc12-2*e2rc12*e2rc20*e3rc12+2*e2rc10*e2rc21*e3rc11-2*e2rc11*e2rc20*e3rc11+2*e2rc12*e2rc22*e3rc10+2*e2rc11*e2rc21*e3rc10+2*e2rc10*e2rc20*e3rc10+2*e2rc00*e2rc22*e3rc02-2*e2rc02*e2rc20*e3rc02+2*e2rc00*e2rc21*e3rc01-2*e2rc01*e2rc20*e3rc01+2*e2rc02*e2rc22*e3rc00+2*e2rc01*e2rc21*e3rc00+2*e2rc00*e2rc20*e3rc00;
    AMatrix(6, 16) = +e0rc20*e3rc22_2+2*e0rc22*e3rc20*e3rc22+2*e0rc10*e3rc12*e3rc22+2*e0rc12*e3rc10*e3rc22+2*e0rc00*e3rc02*e3rc22+2*e0rc02*e3rc00*e3rc22+e0rc20*e3rc21_2+2*e0rc21*e3rc20*e3rc21+2*e0rc10*e3rc11*e3rc21+2*e0rc11*e3rc10*e3rc21+2*e0rc00*e3rc01*e3rc21+2*e0rc01*e3rc00*e3rc21+3*e0rc20*e3rc20_2-2*e0rc12*e3rc12*e3rc20-2*e0rc11*e3rc11*e3rc20+2*e0rc10*e3rc10*e3rc20-2*e0rc02*e3rc02*e3rc20-2*e0rc01*e3rc01*e3rc20+2*e0rc00*e3rc00*e3rc20-e0rc20*e3rc12_2+2*e0rc22*e3rc10*e3rc12-e0rc20*e3rc11_2+2*e0rc21*e3rc10*e3rc11+e0rc20*e3rc10_2-e0rc20*e3rc02_2+2*e0rc22*e3rc00*e3rc02-e0rc20*e3rc01_2+2*e0rc21*e3rc00*e3rc01+e0rc20*e3rc00_2;
    AMatrix(6, 17) = +e1rc20*e3rc22_2+2*e1rc22*e3rc20*e3rc22+2*e1rc10*e3rc12*e3rc22+2*e1rc12*e3rc10*e3rc22+2*e1rc00*e3rc02*e3rc22+2*e1rc02*e3rc00*e3rc22+e1rc20*e3rc21_2+2*e1rc21*e3rc20*e3rc21+2*e1rc10*e3rc11*e3rc21+2*e1rc11*e3rc10*e3rc21+2*e1rc00*e3rc01*e3rc21+2*e1rc01*e3rc00*e3rc21+3*e1rc20*e3rc20_2-2*e1rc12*e3rc12*e3rc20-2*e1rc11*e3rc11*e3rc20+2*e1rc10*e3rc10*e3rc20-2*e1rc02*e3rc02*e3rc20-2*e1rc01*e3rc01*e3rc20+2*e1rc00*e3rc00*e3rc20-e1rc20*e3rc12_2+2*e1rc22*e3rc10*e3rc12-e1rc20*e3rc11_2+2*e1rc21*e3rc10*e3rc11+e1rc20*e3rc10_2-e1rc20*e3rc02_2+2*e1rc22*e3rc00*e3rc02-e1rc20*e3rc01_2+2*e1rc21*e3rc00*e3rc01+e1rc20*e3rc00_2;
    AMatrix(6, 18) = +e2rc20*e3rc22_2+2*e2rc22*e3rc20*e3rc22+2*e2rc10*e3rc12*e3rc22+2*e2rc12*e3rc10*e3rc22+2*e2rc00*e3rc02*e3rc22+2*e2rc02*e3rc00*e3rc22+e2rc20*e3rc21_2+2*e2rc21*e3rc20*e3rc21+2*e2rc10*e3rc11*e3rc21+2*e2rc11*e3rc10*e3rc21+2*e2rc00*e3rc01*e3rc21+2*e2rc01*e3rc00*e3rc21+3*e2rc20*e3rc20_2-2*e2rc12*e3rc12*e3rc20-2*e2rc11*e3rc11*e3rc20+2*e2rc10*e3rc10*e3rc20-2*e2rc02*e3rc02*e3rc20-2*e2rc01*e3rc01*e3rc20+2*e2rc00*e3rc00*e3rc20-e2rc20*e3rc12_2+2*e2rc22*e3rc10*e3rc12-e2rc20*e3rc11_2+2*e2rc21*e3rc10*e3rc11+e2rc20*e3rc10_2-e2rc20*e3rc02_2+2*e2rc22*e3rc00*e3rc02-e2rc20*e3rc01_2+2*e2rc21*e3rc00*e3rc01+e2rc20*e3rc00_2;
    AMatrix(6, 19) = +e3rc20*e3rc22_2+2*e3rc10*e3rc12*e3rc22+2*e3rc00*e3rc02*e3rc22+e3rc20*e3rc21_2+2*e3rc10*e3rc11*e3rc21+2*e3rc00*e3rc01*e3rc21+e3rc20_3-e3rc12_2*e3rc20-e3rc11_2*e3rc20+e3rc10_2*e3rc20-e3rc02_2*e3rc20-e3rc01_2*e3rc20+e3rc00_2*e3rc20;
    AMatrix(7, 0) = +e0rc21*e0rc22_2+2*e0rc11*e0rc12*e0rc22+2*e0rc01*e0rc02*e0rc22+e0rc21_3+e0rc20_2*e0rc21-e0rc12_2*e0rc21+e0rc11_2*e0rc21-e0rc10_2*e0rc21-e0rc02_2*e0rc21+e0rc01_2*e0rc21-e0rc00_2*e0rc21+2*e0rc10*e0rc11*e0rc20+2*e0rc00*e0rc01*e0rc20;
    AMatrix(7, 1) = +2*e0rc21*e0rc22*e1rc22+2*e0rc11*e0rc12*e1rc22+2*e0rc01*e0rc02*e1rc22+e0rc22_2*e1rc21+3*e0rc21_2*e1rc21+e0rc20_2*e1rc21-e0rc12_2*e1rc21+e0rc11_2*e1rc21-e0rc10_2*e1rc21-e0rc02_2*e1rc21+e0rc01_2*e1rc21-e0rc00_2*e1rc21+2*e0rc20*e0rc21*e1rc20+2*e0rc10*e0rc11*e1rc20+2*e0rc00*e0rc01*e1rc20+2*e0rc11*e0rc22*e1rc12-2*e0rc12*e0rc21*e1rc12+2*e0rc12*e0rc22*e1rc11+2*e0rc11*e0rc21*e1rc11+2*e0rc10*e0rc20*e1rc11-2*e0rc10*e0rc21*e1rc10+2*e0rc11*e0rc20*e1rc10+2*e0rc01*e0rc22*e1rc02-2*e0rc02*e0rc21*e1rc02+2*e0rc02*e0rc22*e1rc01+2*e0rc01*e0rc21*e1rc01+2*e0rc00*e0rc20*e1rc01-2*e0rc00*e0rc21*e1rc00+2*e0rc01*e0rc20*e1rc00;
    AMatrix(7, 2) = +2*e0rc21*e0rc22*e2rc22+2*e0rc11*e0rc12*e2rc22+2*e0rc01*e0rc02*e2rc22+e0rc22_2*e2rc21+3*e0rc21_2*e2rc21+e0rc20_2*e2rc21-e0rc12_2*e2rc21+e0rc11_2*e2rc21-e0rc10_2*e2rc21-e0rc02_2*e2rc21+e0rc01_2*e2rc21-e0rc00_2*e2rc21+2*e0rc20*e0rc21*e2rc20+2*e0rc10*e0rc11*e2rc20+2*e0rc00*e0rc01*e2rc20+2*e0rc11*e0rc22*e2rc12-2*e0rc12*e0rc21*e2rc12+2*e0rc12*e0rc22*e2rc11+2*e0rc11*e0rc21*e2rc11+2*e0rc10*e0rc20*e2rc11-2*e0rc10*e0rc21*e2rc10+2*e0rc11*e0rc20*e2rc10+2*e0rc01*e0rc22*e2rc02-2*e0rc02*e0rc21*e2rc02+2*e0rc02*e0rc22*e2rc01+2*e0rc01*e0rc21*e2rc01+2*e0rc00*e0rc20*e2rc01-2*e0rc00*e0rc21*e2rc00+2*e0rc01*e0rc20*e2rc00;
    AMatrix(7, 3) = +e0rc21*e1rc22_2+2*e0rc22*e1rc21*e1rc22+2*e0rc11*e1rc12*e1rc22+2*e0rc12*e1rc11*e1rc22+2*e0rc01*e1rc02*e1rc22+2*e0rc02*e1rc01*e1rc22+3*e0rc21*e1rc21_2+2*e0rc20*e1rc20*e1rc21-2*e0rc12*e1rc12*e1rc21+2*e0rc11*e1rc11*e1rc21-2*e0rc10*e1rc10*e1rc21-2*e0rc02*e1rc02*e1rc21+2*e0rc01*e1rc01*e1rc21-2*e0rc00*e1rc00*e1rc21+e0rc21*e1rc20_2+2*e0rc10*e1rc11*e1rc20+2*e0rc11*e1rc10*e1rc20+2*e0rc00*e1rc01*e1rc20+2*e0rc01*e1rc00*e1rc20-e0rc21*e1rc12_2+2*e0rc22*e1rc11*e1rc12+e0rc21*e1rc11_2+2*e0rc20*e1rc10*e1rc11-e0rc21*e1rc10_2-e0rc21*e1rc02_2+2*e0rc22*e1rc01*e1rc02+e0rc21*e1rc01_2+2*e0rc20*e1rc00*e1rc01-e0rc21*e1rc00_2;
    AMatrix(7, 4) = +2*e0rc21*e1rc22*e2rc22+2*e0rc22*e1rc21*e2rc22+2*e0rc11*e1rc12*e2rc22+2*e0rc12*e1rc11*e2rc22+2*e0rc01*e1rc02*e2rc22+2*e0rc02*e1rc01*e2rc22+2*e0rc22*e1rc22*e2rc21+6*e0rc21*e1rc21*e2rc21+2*e0rc20*e1rc20*e2rc21-2*e0rc12*e1rc12*e2rc21+2*e0rc11*e1rc11*e2rc21-2*e0rc10*e1rc10*e2rc21-2*e0rc02*e1rc02*e2rc21+2*e0rc01*e1rc01*e2rc21-2*e0rc00*e1rc00*e2rc21+2*e0rc20*e1rc21*e2rc20+2*e0rc21*e1rc20*e2rc20+2*e0rc10*e1rc11*e2rc20+2*e0rc11*e1rc10*e2rc20+2*e0rc00*e1rc01*e2rc20+2*e0rc01*e1rc00*e2rc20+2*e0rc11*e1rc22*e2rc12-2*e0rc12*e1rc21*e2rc12-2*e0rc21*e1rc12*e2rc12+2*e0rc22*e1rc11*e2rc12+2*e0rc12*e1rc22*e2rc11+2*e0rc11*e1rc21*e2rc11+2*e0rc10*e1rc20*e2rc11+2*e0rc22*e1rc12*e2rc11+2*e0rc21*e1rc11*e2rc11+2*e0rc20*e1rc10*e2rc11-2*e0rc10*e1rc21*e2rc10+2*e0rc11*e1rc20*e2rc10+2*e0rc20*e1rc11*e2rc10-2*e0rc21*e1rc10*e2rc10+2*e0rc01*e1rc22*e2rc02-2*e0rc02*e1rc21*e2rc02-2*e0rc21*e1rc02*e2rc02+2*e0rc22*e1rc01*e2rc02+2*e0rc02*e1rc22*e2rc01+2*e0rc01*e1rc21*e2rc01+2*e0rc00*e1rc20*e2rc01+2*e0rc22*e1rc02*e2rc01+2*e0rc21*e1rc01*e2rc01+2*e0rc20*e1rc00*e2rc01-2*e0rc00*e1rc21*e2rc00+2*e0rc01*e1rc20*e2rc00+2*e0rc20*e1rc01*e2rc00-2*e0rc21*e1rc00*e2rc00;
    AMatrix(7, 5) = +e0rc21*e2rc22_2+2*e0rc22*e2rc21*e2rc22+2*e0rc11*e2rc12*e2rc22+2*e0rc12*e2rc11*e2rc22+2*e0rc01*e2rc02*e2rc22+2*e0rc02*e2rc01*e2rc22+3*e0rc21*e2rc21_2+2*e0rc20*e2rc20*e2rc21-2*e0rc12*e2rc12*e2rc21+2*e0rc11*e2rc11*e2rc21-2*e0rc10*e2rc10*e2rc21-2*e0rc02*e2rc02*e2rc21+2*e0rc01*e2rc01*e2rc21-2*e0rc00*e2rc00*e2rc21+e0rc21*e2rc20_2+2*e0rc10*e2rc11*e2rc20+2*e0rc11*e2rc10*e2rc20+2*e0rc00*e2rc01*e2rc20+2*e0rc01*e2rc00*e2rc20-e0rc21*e2rc12_2+2*e0rc22*e2rc11*e2rc12+e0rc21*e2rc11_2+2*e0rc20*e2rc10*e2rc11-e0rc21*e2rc10_2-e0rc21*e2rc02_2+2*e0rc22*e2rc01*e2rc02+e0rc21*e2rc01_2+2*e0rc20*e2rc00*e2rc01-e0rc21*e2rc00_2;
    AMatrix(7, 6) = +e1rc21*e1rc22_2+2*e1rc11*e1rc12*e1rc22+2*e1rc01*e1rc02*e1rc22+e1rc21_3+e1rc20_2*e1rc21-e1rc12_2*e1rc21+e1rc11_2*e1rc21-e1rc10_2*e1rc21-e1rc02_2*e1rc21+e1rc01_2*e1rc21-e1rc00_2*e1rc21+2*e1rc10*e1rc11*e1rc20+2*e1rc00*e1rc01*e1rc20;
    AMatrix(7, 7) = +2*e1rc21*e1rc22*e2rc22+2*e1rc11*e1rc12*e2rc22+2*e1rc01*e1rc02*e2rc22+e1rc22_2*e2rc21+3*e1rc21_2*e2rc21+e1rc20_2*e2rc21-e1rc12_2*e2rc21+e1rc11_2*e2rc21-e1rc10_2*e2rc21-e1rc02_2*e2rc21+e1rc01_2*e2rc21-e1rc00_2*e2rc21+2*e1rc20*e1rc21*e2rc20+2*e1rc10*e1rc11*e2rc20+2*e1rc00*e1rc01*e2rc20+2*e1rc11*e1rc22*e2rc12-2*e1rc12*e1rc21*e2rc12+2*e1rc12*e1rc22*e2rc11+2*e1rc11*e1rc21*e2rc11+2*e1rc10*e1rc20*e2rc11-2*e1rc10*e1rc21*e2rc10+2*e1rc11*e1rc20*e2rc10+2*e1rc01*e1rc22*e2rc02-2*e1rc02*e1rc21*e2rc02+2*e1rc02*e1rc22*e2rc01+2*e1rc01*e1rc21*e2rc01+2*e1rc00*e1rc20*e2rc01-2*e1rc00*e1rc21*e2rc00+2*e1rc01*e1rc20*e2rc00;
    AMatrix(7, 8) = +e1rc21*e2rc22_2+2*e1rc22*e2rc21*e2rc22+2*e1rc11*e2rc12*e2rc22+2*e1rc12*e2rc11*e2rc22+2*e1rc01*e2rc02*e2rc22+2*e1rc02*e2rc01*e2rc22+3*e1rc21*e2rc21_2+2*e1rc20*e2rc20*e2rc21-2*e1rc12*e2rc12*e2rc21+2*e1rc11*e2rc11*e2rc21-2*e1rc10*e2rc10*e2rc21-2*e1rc02*e2rc02*e2rc21+2*e1rc01*e2rc01*e2rc21-2*e1rc00*e2rc00*e2rc21+e1rc21*e2rc20_2+2*e1rc10*e2rc11*e2rc20+2*e1rc11*e2rc10*e2rc20+2*e1rc00*e2rc01*e2rc20+2*e1rc01*e2rc00*e2rc20-e1rc21*e2rc12_2+2*e1rc22*e2rc11*e2rc12+e1rc21*e2rc11_2+2*e1rc20*e2rc10*e2rc11-e1rc21*e2rc10_2-e1rc21*e2rc02_2+2*e1rc22*e2rc01*e2rc02+e1rc21*e2rc01_2+2*e1rc20*e2rc00*e2rc01-e1rc21*e2rc00_2;
    AMatrix(7, 9) = e2rc21*e2rc22_2+2*e2rc11*e2rc12*e2rc22+2*e2rc01*e2rc02*e2rc22+e2rc21_3+e2rc20_2*e2rc21-e2rc12_2*e2rc21+e2rc11_2*e2rc21-e2rc10_2*e2rc21-e2rc02_2*e2rc21+e2rc01_2*e2rc21-e2rc00_2*e2rc21+2*e2rc10*e2rc11*e2rc20+2*e2rc00*e2rc01*e2rc20;
    AMatrix(7, 10) = +2*e0rc21*e0rc22*e3rc22+2*e0rc11*e0rc12*e3rc22+2*e0rc01*e0rc02*e3rc22+e0rc22_2*e3rc21+3*e0rc21_2*e3rc21+e0rc20_2*e3rc21-e0rc12_2*e3rc21+e0rc11_2*e3rc21-e0rc10_2*e3rc21-e0rc02_2*e3rc21+e0rc01_2*e3rc21-e0rc00_2*e3rc21+2*e0rc20*e0rc21*e3rc20+2*e0rc10*e0rc11*e3rc20+2*e0rc00*e0rc01*e3rc20+2*e0rc11*e0rc22*e3rc12-2*e0rc12*e0rc21*e3rc12+2*e0rc12*e0rc22*e3rc11+2*e0rc11*e0rc21*e3rc11+2*e0rc10*e0rc20*e3rc11-2*e0rc10*e0rc21*e3rc10+2*e0rc11*e0rc20*e3rc10+2*e0rc01*e0rc22*e3rc02-2*e0rc02*e0rc21*e3rc02+2*e0rc02*e0rc22*e3rc01+2*e0rc01*e0rc21*e3rc01+2*e0rc00*e0rc20*e3rc01-2*e0rc00*e0rc21*e3rc00+2*e0rc01*e0rc20*e3rc00;
    AMatrix(7, 11) = +2*e0rc21*e1rc22*e3rc22+2*e0rc22*e1rc21*e3rc22+2*e0rc11*e1rc12*e3rc22+2*e0rc12*e1rc11*e3rc22+2*e0rc01*e1rc02*e3rc22+2*e0rc02*e1rc01*e3rc22+2*e0rc22*e1rc22*e3rc21+6*e0rc21*e1rc21*e3rc21+2*e0rc20*e1rc20*e3rc21-2*e0rc12*e1rc12*e3rc21+2*e0rc11*e1rc11*e3rc21-2*e0rc10*e1rc10*e3rc21-2*e0rc02*e1rc02*e3rc21+2*e0rc01*e1rc01*e3rc21-2*e0rc00*e1rc00*e3rc21+2*e0rc20*e1rc21*e3rc20+2*e0rc21*e1rc20*e3rc20+2*e0rc10*e1rc11*e3rc20+2*e0rc11*e1rc10*e3rc20+2*e0rc00*e1rc01*e3rc20+2*e0rc01*e1rc00*e3rc20+2*e0rc11*e1rc22*e3rc12-2*e0rc12*e1rc21*e3rc12-2*e0rc21*e1rc12*e3rc12+2*e0rc22*e1rc11*e3rc12+2*e0rc12*e1rc22*e3rc11+2*e0rc11*e1rc21*e3rc11+2*e0rc10*e1rc20*e3rc11+2*e0rc22*e1rc12*e3rc11+2*e0rc21*e1rc11*e3rc11+2*e0rc20*e1rc10*e3rc11-2*e0rc10*e1rc21*e3rc10+2*e0rc11*e1rc20*e3rc10+2*e0rc20*e1rc11*e3rc10-2*e0rc21*e1rc10*e3rc10+2*e0rc01*e1rc22*e3rc02-2*e0rc02*e1rc21*e3rc02-2*e0rc21*e1rc02*e3rc02+2*e0rc22*e1rc01*e3rc02+2*e0rc02*e1rc22*e3rc01+2*e0rc01*e1rc21*e3rc01+2*e0rc00*e1rc20*e3rc01+2*e0rc22*e1rc02*e3rc01+2*e0rc21*e1rc01*e3rc01+2*e0rc20*e1rc00*e3rc01-2*e0rc00*e1rc21*e3rc00+2*e0rc01*e1rc20*e3rc00+2*e0rc20*e1rc01*e3rc00-2*e0rc21*e1rc00*e3rc00;
    AMatrix(7, 12) = +2*e0rc21*e2rc22*e3rc22+2*e0rc22*e2rc21*e3rc22+2*e0rc11*e2rc12*e3rc22+2*e0rc12*e2rc11*e3rc22+2*e0rc01*e2rc02*e3rc22+2*e0rc02*e2rc01*e3rc22+2*e0rc22*e2rc22*e3rc21+6*e0rc21*e2rc21*e3rc21+2*e0rc20*e2rc20*e3rc21-2*e0rc12*e2rc12*e3rc21+2*e0rc11*e2rc11*e3rc21-2*e0rc10*e2rc10*e3rc21-2*e0rc02*e2rc02*e3rc21+2*e0rc01*e2rc01*e3rc21-2*e0rc00*e2rc00*e3rc21+2*e0rc20*e2rc21*e3rc20+2*e0rc21*e2rc20*e3rc20+2*e0rc10*e2rc11*e3rc20+2*e0rc11*e2rc10*e3rc20+2*e0rc00*e2rc01*e3rc20+2*e0rc01*e2rc00*e3rc20+2*e0rc11*e2rc22*e3rc12-2*e0rc12*e2rc21*e3rc12-2*e0rc21*e2rc12*e3rc12+2*e0rc22*e2rc11*e3rc12+2*e0rc12*e2rc22*e3rc11+2*e0rc11*e2rc21*e3rc11+2*e0rc10*e2rc20*e3rc11+2*e0rc22*e2rc12*e3rc11+2*e0rc21*e2rc11*e3rc11+2*e0rc20*e2rc10*e3rc11-2*e0rc10*e2rc21*e3rc10+2*e0rc11*e2rc20*e3rc10+2*e0rc20*e2rc11*e3rc10-2*e0rc21*e2rc10*e3rc10+2*e0rc01*e2rc22*e3rc02-2*e0rc02*e2rc21*e3rc02-2*e0rc21*e2rc02*e3rc02+2*e0rc22*e2rc01*e3rc02+2*e0rc02*e2rc22*e3rc01+2*e0rc01*e2rc21*e3rc01+2*e0rc00*e2rc20*e3rc01+2*e0rc22*e2rc02*e3rc01+2*e0rc21*e2rc01*e3rc01+2*e0rc20*e2rc00*e3rc01-2*e0rc00*e2rc21*e3rc00+2*e0rc01*e2rc20*e3rc00+2*e0rc20*e2rc01*e3rc00-2*e0rc21*e2rc00*e3rc00;
    AMatrix(7, 13) = +2*e1rc21*e1rc22*e3rc22+2*e1rc11*e1rc12*e3rc22+2*e1rc01*e1rc02*e3rc22+e1rc22_2*e3rc21+3*e1rc21_2*e3rc21+e1rc20_2*e3rc21-e1rc12_2*e3rc21+e1rc11_2*e3rc21-e1rc10_2*e3rc21-e1rc02_2*e3rc21+e1rc01_2*e3rc21-e1rc00_2*e3rc21+2*e1rc20*e1rc21*e3rc20+2*e1rc10*e1rc11*e3rc20+2*e1rc00*e1rc01*e3rc20+2*e1rc11*e1rc22*e3rc12-2*e1rc12*e1rc21*e3rc12+2*e1rc12*e1rc22*e3rc11+2*e1rc11*e1rc21*e3rc11+2*e1rc10*e1rc20*e3rc11-2*e1rc10*e1rc21*e3rc10+2*e1rc11*e1rc20*e3rc10+2*e1rc01*e1rc22*e3rc02-2*e1rc02*e1rc21*e3rc02+2*e1rc02*e1rc22*e3rc01+2*e1rc01*e1rc21*e3rc01+2*e1rc00*e1rc20*e3rc01-2*e1rc00*e1rc21*e3rc00+2*e1rc01*e1rc20*e3rc00;
    AMatrix(7, 14) = +2*e1rc21*e2rc22*e3rc22+2*e1rc22*e2rc21*e3rc22+2*e1rc11*e2rc12*e3rc22+2*e1rc12*e2rc11*e3rc22+2*e1rc01*e2rc02*e3rc22+2*e1rc02*e2rc01*e3rc22+2*e1rc22*e2rc22*e3rc21+6*e1rc21*e2rc21*e3rc21+2*e1rc20*e2rc20*e3rc21-2*e1rc12*e2rc12*e3rc21+2*e1rc11*e2rc11*e3rc21-2*e1rc10*e2rc10*e3rc21-2*e1rc02*e2rc02*e3rc21+2*e1rc01*e2rc01*e3rc21-2*e1rc00*e2rc00*e3rc21+2*e1rc20*e2rc21*e3rc20+2*e1rc21*e2rc20*e3rc20+2*e1rc10*e2rc11*e3rc20+2*e1rc11*e2rc10*e3rc20+2*e1rc00*e2rc01*e3rc20+2*e1rc01*e2rc00*e3rc20+2*e1rc11*e2rc22*e3rc12-2*e1rc12*e2rc21*e3rc12-2*e1rc21*e2rc12*e3rc12+2*e1rc22*e2rc11*e3rc12+2*e1rc12*e2rc22*e3rc11+2*e1rc11*e2rc21*e3rc11+2*e1rc10*e2rc20*e3rc11+2*e1rc22*e2rc12*e3rc11+2*e1rc21*e2rc11*e3rc11+2*e1rc20*e2rc10*e3rc11-2*e1rc10*e2rc21*e3rc10+2*e1rc11*e2rc20*e3rc10+2*e1rc20*e2rc11*e3rc10-2*e1rc21*e2rc10*e3rc10+2*e1rc01*e2rc22*e3rc02-2*e1rc02*e2rc21*e3rc02-2*e1rc21*e2rc02*e3rc02+2*e1rc22*e2rc01*e3rc02+2*e1rc02*e2rc22*e3rc01+2*e1rc01*e2rc21*e3rc01+2*e1rc00*e2rc20*e3rc01+2*e1rc22*e2rc02*e3rc01+2*e1rc21*e2rc01*e3rc01+2*e1rc20*e2rc00*e3rc01-2*e1rc00*e2rc21*e3rc00+2*e1rc01*e2rc20*e3rc00+2*e1rc20*e2rc01*e3rc00-2*e1rc21*e2rc00*e3rc00;
    AMatrix(7, 15) = +2*e2rc21*e2rc22*e3rc22+2*e2rc11*e2rc12*e3rc22+2*e2rc01*e2rc02*e3rc22+e2rc22_2*e3rc21+3*e2rc21_2*e3rc21+e2rc20_2*e3rc21-e2rc12_2*e3rc21+e2rc11_2*e3rc21-e2rc10_2*e3rc21-e2rc02_2*e3rc21+e2rc01_2*e3rc21-e2rc00_2*e3rc21+2*e2rc20*e2rc21*e3rc20+2*e2rc10*e2rc11*e3rc20+2*e2rc00*e2rc01*e3rc20+2*e2rc11*e2rc22*e3rc12-2*e2rc12*e2rc21*e3rc12+2*e2rc12*e2rc22*e3rc11+2*e2rc11*e2rc21*e3rc11+2*e2rc10*e2rc20*e3rc11-2*e2rc10*e2rc21*e3rc10+2*e2rc11*e2rc20*e3rc10+2*e2rc01*e2rc22*e3rc02-2*e2rc02*e2rc21*e3rc02+2*e2rc02*e2rc22*e3rc01+2*e2rc01*e2rc21*e3rc01+2*e2rc00*e2rc20*e3rc01-2*e2rc00*e2rc21*e3rc00+2*e2rc01*e2rc20*e3rc00;
    AMatrix(7, 16) = +e0rc21*e3rc22_2+2*e0rc22*e3rc21*e3rc22+2*e0rc11*e3rc12*e3rc22+2*e0rc12*e3rc11*e3rc22+2*e0rc01*e3rc02*e3rc22+2*e0rc02*e3rc01*e3rc22+3*e0rc21*e3rc21_2+2*e0rc20*e3rc20*e3rc21-2*e0rc12*e3rc12*e3rc21+2*e0rc11*e3rc11*e3rc21-2*e0rc10*e3rc10*e3rc21-2*e0rc02*e3rc02*e3rc21+2*e0rc01*e3rc01*e3rc21-2*e0rc00*e3rc00*e3rc21+e0rc21*e3rc20_2+2*e0rc10*e3rc11*e3rc20+2*e0rc11*e3rc10*e3rc20+2*e0rc00*e3rc01*e3rc20+2*e0rc01*e3rc00*e3rc20-e0rc21*e3rc12_2+2*e0rc22*e3rc11*e3rc12+e0rc21*e3rc11_2+2*e0rc20*e3rc10*e3rc11-e0rc21*e3rc10_2-e0rc21*e3rc02_2+2*e0rc22*e3rc01*e3rc02+e0rc21*e3rc01_2+2*e0rc20*e3rc00*e3rc01-e0rc21*e3rc00_2;
    AMatrix(7, 17) = +e1rc21*e3rc22_2+2*e1rc22*e3rc21*e3rc22+2*e1rc11*e3rc12*e3rc22+2*e1rc12*e3rc11*e3rc22+2*e1rc01*e3rc02*e3rc22+2*e1rc02*e3rc01*e3rc22+3*e1rc21*e3rc21_2+2*e1rc20*e3rc20*e3rc21-2*e1rc12*e3rc12*e3rc21+2*e1rc11*e3rc11*e3rc21-2*e1rc10*e3rc10*e3rc21-2*e1rc02*e3rc02*e3rc21+2*e1rc01*e3rc01*e3rc21-2*e1rc00*e3rc00*e3rc21+e1rc21*e3rc20_2+2*e1rc10*e3rc11*e3rc20+2*e1rc11*e3rc10*e3rc20+2*e1rc00*e3rc01*e3rc20+2*e1rc01*e3rc00*e3rc20-e1rc21*e3rc12_2+2*e1rc22*e3rc11*e3rc12+e1rc21*e3rc11_2+2*e1rc20*e3rc10*e3rc11-e1rc21*e3rc10_2-e1rc21*e3rc02_2+2*e1rc22*e3rc01*e3rc02+e1rc21*e3rc01_2+2*e1rc20*e3rc00*e3rc01-e1rc21*e3rc00_2;
    AMatrix(7, 18) = +e2rc21*e3rc22_2+2*e2rc22*e3rc21*e3rc22+2*e2rc11*e3rc12*e3rc22+2*e2rc12*e3rc11*e3rc22+2*e2rc01*e3rc02*e3rc22+2*e2rc02*e3rc01*e3rc22+3*e2rc21*e3rc21_2+2*e2rc20*e3rc20*e3rc21-2*e2rc12*e3rc12*e3rc21+2*e2rc11*e3rc11*e3rc21-2*e2rc10*e3rc10*e3rc21-2*e2rc02*e3rc02*e3rc21+2*e2rc01*e3rc01*e3rc21-2*e2rc00*e3rc00*e3rc21+e2rc21*e3rc20_2+2*e2rc10*e3rc11*e3rc20+2*e2rc11*e3rc10*e3rc20+2*e2rc00*e3rc01*e3rc20+2*e2rc01*e3rc00*e3rc20-e2rc21*e3rc12_2+2*e2rc22*e3rc11*e3rc12+e2rc21*e3rc11_2+2*e2rc20*e3rc10*e3rc11-e2rc21*e3rc10_2-e2rc21*e3rc02_2+2*e2rc22*e3rc01*e3rc02+e2rc21*e3rc01_2+2*e2rc20*e3rc00*e3rc01-e2rc21*e3rc00_2;
    AMatrix(7, 19) = +e3rc21*e3rc22_2+2*e3rc11*e3rc12*e3rc22+2*e3rc01*e3rc02*e3rc22+e3rc21_3+e3rc20_2*e3rc21-e3rc12_2*e3rc21+e3rc11_2*e3rc21-e3rc10_2*e3rc21-e3rc02_2*e3rc21+e3rc01_2*e3rc21-e3rc00_2*e3rc21+2*e3rc10*e3rc11*e3rc20+2*e3rc00*e3rc01*e3rc20;
    AMatrix(8, 0) = +e0rc22_3+e0rc21_2*e0rc22+e0rc20_2*e0rc22+e0rc12_2*e0rc22-e0rc11_2*e0rc22-e0rc10_2*e0rc22+e0rc02_2*e0rc22-e0rc01_2*e0rc22-e0rc00_2*e0rc22+2*e0rc11*e0rc12*e0rc21+2*e0rc01*e0rc02*e0rc21+2*e0rc10*e0rc12*e0rc20+2*e0rc00*e0rc02*e0rc20;
    AMatrix(8, 1) = +3*e0rc22_2*e1rc22+e0rc21_2*e1rc22+e0rc20_2*e1rc22+e0rc12_2*e1rc22-e0rc11_2*e1rc22-e0rc10_2*e1rc22+e0rc02_2*e1rc22-e0rc01_2*e1rc22-e0rc00_2*e1rc22+2*e0rc21*e0rc22*e1rc21+2*e0rc11*e0rc12*e1rc21+2*e0rc01*e0rc02*e1rc21+2*e0rc20*e0rc22*e1rc20+2*e0rc10*e0rc12*e1rc20+2*e0rc00*e0rc02*e1rc20+2*e0rc12*e0rc22*e1rc12+2*e0rc11*e0rc21*e1rc12+2*e0rc10*e0rc20*e1rc12-2*e0rc11*e0rc22*e1rc11+2*e0rc12*e0rc21*e1rc11-2*e0rc10*e0rc22*e1rc10+2*e0rc12*e0rc20*e1rc10+2*e0rc02*e0rc22*e1rc02+2*e0rc01*e0rc21*e1rc02+2*e0rc00*e0rc20*e1rc02-2*e0rc01*e0rc22*e1rc01+2*e0rc02*e0rc21*e1rc01-2*e0rc00*e0rc22*e1rc00+2*e0rc02*e0rc20*e1rc00;
    AMatrix(8, 2) = +3*e0rc22_2*e2rc22+e0rc21_2*e2rc22+e0rc20_2*e2rc22+e0rc12_2*e2rc22-e0rc11_2*e2rc22-e0rc10_2*e2rc22+e0rc02_2*e2rc22-e0rc01_2*e2rc22-e0rc00_2*e2rc22+2*e0rc21*e0rc22*e2rc21+2*e0rc11*e0rc12*e2rc21+2*e0rc01*e0rc02*e2rc21+2*e0rc20*e0rc22*e2rc20+2*e0rc10*e0rc12*e2rc20+2*e0rc00*e0rc02*e2rc20+2*e0rc12*e0rc22*e2rc12+2*e0rc11*e0rc21*e2rc12+2*e0rc10*e0rc20*e2rc12-2*e0rc11*e0rc22*e2rc11+2*e0rc12*e0rc21*e2rc11-2*e0rc10*e0rc22*e2rc10+2*e0rc12*e0rc20*e2rc10+2*e0rc02*e0rc22*e2rc02+2*e0rc01*e0rc21*e2rc02+2*e0rc00*e0rc20*e2rc02-2*e0rc01*e0rc22*e2rc01+2*e0rc02*e0rc21*e2rc01-2*e0rc00*e0rc22*e2rc00+2*e0rc02*e0rc20*e2rc00;
    AMatrix(8, 3) = +3*e0rc22*e1rc22_2+2*e0rc21*e1rc21*e1rc22+2*e0rc20*e1rc20*e1rc22+2*e0rc12*e1rc12*e1rc22-2*e0rc11*e1rc11*e1rc22-2*e0rc10*e1rc10*e1rc22+2*e0rc02*e1rc02*e1rc22-2*e0rc01*e1rc01*e1rc22-2*e0rc00*e1rc00*e1rc22+e0rc22*e1rc21_2+2*e0rc11*e1rc12*e1rc21+2*e0rc12*e1rc11*e1rc21+2*e0rc01*e1rc02*e1rc21+2*e0rc02*e1rc01*e1rc21+e0rc22*e1rc20_2+2*e0rc10*e1rc12*e1rc20+2*e0rc12*e1rc10*e1rc20+2*e0rc00*e1rc02*e1rc20+2*e0rc02*e1rc00*e1rc20+e0rc22*e1rc12_2+2*e0rc21*e1rc11*e1rc12+2*e0rc20*e1rc10*e1rc12-e0rc22*e1rc11_2-e0rc22*e1rc10_2+e0rc22*e1rc02_2+2*e0rc21*e1rc01*e1rc02+2*e0rc20*e1rc00*e1rc02-e0rc22*e1rc01_2-e0rc22*e1rc00_2;
    AMatrix(8, 4) = +6*e0rc22*e1rc22*e2rc22+2*e0rc21*e1rc21*e2rc22+2*e0rc20*e1rc20*e2rc22+2*e0rc12*e1rc12*e2rc22-2*e0rc11*e1rc11*e2rc22-2*e0rc10*e1rc10*e2rc22+2*e0rc02*e1rc02*e2rc22-2*e0rc01*e1rc01*e2rc22-2*e0rc00*e1rc00*e2rc22+2*e0rc21*e1rc22*e2rc21+2*e0rc22*e1rc21*e2rc21+2*e0rc11*e1rc12*e2rc21+2*e0rc12*e1rc11*e2rc21+2*e0rc01*e1rc02*e2rc21+2*e0rc02*e1rc01*e2rc21+2*e0rc20*e1rc22*e2rc20+2*e0rc22*e1rc20*e2rc20+2*e0rc10*e1rc12*e2rc20+2*e0rc12*e1rc10*e2rc20+2*e0rc00*e1rc02*e2rc20+2*e0rc02*e1rc00*e2rc20+2*e0rc12*e1rc22*e2rc12+2*e0rc11*e1rc21*e2rc12+2*e0rc10*e1rc20*e2rc12+2*e0rc22*e1rc12*e2rc12+2*e0rc21*e1rc11*e2rc12+2*e0rc20*e1rc10*e2rc12-2*e0rc11*e1rc22*e2rc11+2*e0rc12*e1rc21*e2rc11+2*e0rc21*e1rc12*e2rc11-2*e0rc22*e1rc11*e2rc11-2*e0rc10*e1rc22*e2rc10+2*e0rc12*e1rc20*e2rc10+2*e0rc20*e1rc12*e2rc10-2*e0rc22*e1rc10*e2rc10+2*e0rc02*e1rc22*e2rc02+2*e0rc01*e1rc21*e2rc02+2*e0rc00*e1rc20*e2rc02+2*e0rc22*e1rc02*e2rc02+2*e0rc21*e1rc01*e2rc02+2*e0rc20*e1rc00*e2rc02-2*e0rc01*e1rc22*e2rc01+2*e0rc02*e1rc21*e2rc01+2*e0rc21*e1rc02*e2rc01-2*e0rc22*e1rc01*e2rc01-2*e0rc00*e1rc22*e2rc00+2*e0rc02*e1rc20*e2rc00+2*e0rc20*e1rc02*e2rc00-2*e0rc22*e1rc00*e2rc00;
    AMatrix(8, 5) = +3*e0rc22*e2rc22_2+2*e0rc21*e2rc21*e2rc22+2*e0rc20*e2rc20*e2rc22+2*e0rc12*e2rc12*e2rc22-2*e0rc11*e2rc11*e2rc22-2*e0rc10*e2rc10*e2rc22+2*e0rc02*e2rc02*e2rc22-2*e0rc01*e2rc01*e2rc22-2*e0rc00*e2rc00*e2rc22+e0rc22*e2rc21_2+2*e0rc11*e2rc12*e2rc21+2*e0rc12*e2rc11*e2rc21+2*e0rc01*e2rc02*e2rc21+2*e0rc02*e2rc01*e2rc21+e0rc22*e2rc20_2+2*e0rc10*e2rc12*e2rc20+2*e0rc12*e2rc10*e2rc20+2*e0rc00*e2rc02*e2rc20+2*e0rc02*e2rc00*e2rc20+e0rc22*e2rc12_2+2*e0rc21*e2rc11*e2rc12+2*e0rc20*e2rc10*e2rc12-e0rc22*e2rc11_2-e0rc22*e2rc10_2+e0rc22*e2rc02_2+2*e0rc21*e2rc01*e2rc02+2*e0rc20*e2rc00*e2rc02-e0rc22*e2rc01_2-e0rc22*e2rc00_2;
    AMatrix(8, 6) = +e1rc22_3+e1rc21_2*e1rc22+e1rc20_2*e1rc22+e1rc12_2*e1rc22-e1rc11_2*e1rc22-e1rc10_2*e1rc22+e1rc02_2*e1rc22-e1rc01_2*e1rc22-e1rc00_2*e1rc22+2*e1rc11*e1rc12*e1rc21+2*e1rc01*e1rc02*e1rc21+2*e1rc10*e1rc12*e1rc20+2*e1rc00*e1rc02*e1rc20;
    AMatrix(8, 7) = +3*e1rc22_2*e2rc22+e1rc21_2*e2rc22+e1rc20_2*e2rc22+e1rc12_2*e2rc22-e1rc11_2*e2rc22-e1rc10_2*e2rc22+e1rc02_2*e2rc22-e1rc01_2*e2rc22-e1rc00_2*e2rc22+2*e1rc21*e1rc22*e2rc21+2*e1rc11*e1rc12*e2rc21+2*e1rc01*e1rc02*e2rc21+2*e1rc20*e1rc22*e2rc20+2*e1rc10*e1rc12*e2rc20+2*e1rc00*e1rc02*e2rc20+2*e1rc12*e1rc22*e2rc12+2*e1rc11*e1rc21*e2rc12+2*e1rc10*e1rc20*e2rc12-2*e1rc11*e1rc22*e2rc11+2*e1rc12*e1rc21*e2rc11-2*e1rc10*e1rc22*e2rc10+2*e1rc12*e1rc20*e2rc10+2*e1rc02*e1rc22*e2rc02+2*e1rc01*e1rc21*e2rc02+2*e1rc00*e1rc20*e2rc02-2*e1rc01*e1rc22*e2rc01+2*e1rc02*e1rc21*e2rc01-2*e1rc00*e1rc22*e2rc00+2*e1rc02*e1rc20*e2rc00;
    AMatrix(8, 8) = +3*e1rc22*e2rc22_2+2*e1rc21*e2rc21*e2rc22+2*e1rc20*e2rc20*e2rc22+2*e1rc12*e2rc12*e2rc22-2*e1rc11*e2rc11*e2rc22-2*e1rc10*e2rc10*e2rc22+2*e1rc02*e2rc02*e2rc22-2*e1rc01*e2rc01*e2rc22-2*e1rc00*e2rc00*e2rc22+e1rc22*e2rc21_2+2*e1rc11*e2rc12*e2rc21+2*e1rc12*e2rc11*e2rc21+2*e1rc01*e2rc02*e2rc21+2*e1rc02*e2rc01*e2rc21+e1rc22*e2rc20_2+2*e1rc10*e2rc12*e2rc20+2*e1rc12*e2rc10*e2rc20+2*e1rc00*e2rc02*e2rc20+2*e1rc02*e2rc00*e2rc20+e1rc22*e2rc12_2+2*e1rc21*e2rc11*e2rc12+2*e1rc20*e2rc10*e2rc12-e1rc22*e2rc11_2-e1rc22*e2rc10_2+e1rc22*e2rc02_2+2*e1rc21*e2rc01*e2rc02+2*e1rc20*e2rc00*e2rc02-e1rc22*e2rc01_2-e1rc22*e2rc00_2;
    AMatrix(8, 9) = e2rc22_3+e2rc21_2*e2rc22+e2rc20_2*e2rc22+e2rc12_2*e2rc22-e2rc11_2*e2rc22-e2rc10_2*e2rc22+e2rc02_2*e2rc22-e2rc01_2*e2rc22-e2rc00_2*e2rc22+2*e2rc11*e2rc12*e2rc21+2*e2rc01*e2rc02*e2rc21+2*e2rc10*e2rc12*e2rc20+2*e2rc00*e2rc02*e2rc20;
    AMatrix(8, 10) = +3*e0rc22_2*e3rc22+e0rc21_2*e3rc22+e0rc20_2*e3rc22+e0rc12_2*e3rc22-e0rc11_2*e3rc22-e0rc10_2*e3rc22+e0rc02_2*e3rc22-e0rc01_2*e3rc22-e0rc00_2*e3rc22+2*e0rc21*e0rc22*e3rc21+2*e0rc11*e0rc12*e3rc21+2*e0rc01*e0rc02*e3rc21+2*e0rc20*e0rc22*e3rc20+2*e0rc10*e0rc12*e3rc20+2*e0rc00*e0rc02*e3rc20+2*e0rc12*e0rc22*e3rc12+2*e0rc11*e0rc21*e3rc12+2*e0rc10*e0rc20*e3rc12-2*e0rc11*e0rc22*e3rc11+2*e0rc12*e0rc21*e3rc11-2*e0rc10*e0rc22*e3rc10+2*e0rc12*e0rc20*e3rc10+2*e0rc02*e0rc22*e3rc02+2*e0rc01*e0rc21*e3rc02+2*e0rc00*e0rc20*e3rc02-2*e0rc01*e0rc22*e3rc01+2*e0rc02*e0rc21*e3rc01-2*e0rc00*e0rc22*e3rc00+2*e0rc02*e0rc20*e3rc00;
    AMatrix(8, 11) = +6*e0rc22*e1rc22*e3rc22+2*e0rc21*e1rc21*e3rc22+2*e0rc20*e1rc20*e3rc22+2*e0rc12*e1rc12*e3rc22-2*e0rc11*e1rc11*e3rc22-2*e0rc10*e1rc10*e3rc22+2*e0rc02*e1rc02*e3rc22-2*e0rc01*e1rc01*e3rc22-2*e0rc00*e1rc00*e3rc22+2*e0rc21*e1rc22*e3rc21+2*e0rc22*e1rc21*e3rc21+2*e0rc11*e1rc12*e3rc21+2*e0rc12*e1rc11*e3rc21+2*e0rc01*e1rc02*e3rc21+2*e0rc02*e1rc01*e3rc21+2*e0rc20*e1rc22*e3rc20+2*e0rc22*e1rc20*e3rc20+2*e0rc10*e1rc12*e3rc20+2*e0rc12*e1rc10*e3rc20+2*e0rc00*e1rc02*e3rc20+2*e0rc02*e1rc00*e3rc20+2*e0rc12*e1rc22*e3rc12+2*e0rc11*e1rc21*e3rc12+2*e0rc10*e1rc20*e3rc12+2*e0rc22*e1rc12*e3rc12+2*e0rc21*e1rc11*e3rc12+2*e0rc20*e1rc10*e3rc12-2*e0rc11*e1rc22*e3rc11+2*e0rc12*e1rc21*e3rc11+2*e0rc21*e1rc12*e3rc11-2*e0rc22*e1rc11*e3rc11-2*e0rc10*e1rc22*e3rc10+2*e0rc12*e1rc20*e3rc10+2*e0rc20*e1rc12*e3rc10-2*e0rc22*e1rc10*e3rc10+2*e0rc02*e1rc22*e3rc02+2*e0rc01*e1rc21*e3rc02+2*e0rc00*e1rc20*e3rc02+2*e0rc22*e1rc02*e3rc02+2*e0rc21*e1rc01*e3rc02+2*e0rc20*e1rc00*e3rc02-2*e0rc01*e1rc22*e3rc01+2*e0rc02*e1rc21*e3rc01+2*e0rc21*e1rc02*e3rc01-2*e0rc22*e1rc01*e3rc01-2*e0rc00*e1rc22*e3rc00+2*e0rc02*e1rc20*e3rc00+2*e0rc20*e1rc02*e3rc00-2*e0rc22*e1rc00*e3rc00;
    AMatrix(8, 12) = +6*e0rc22*e2rc22*e3rc22+2*e0rc21*e2rc21*e3rc22+2*e0rc20*e2rc20*e3rc22+2*e0rc12*e2rc12*e3rc22-2*e0rc11*e2rc11*e3rc22-2*e0rc10*e2rc10*e3rc22+2*e0rc02*e2rc02*e3rc22-2*e0rc01*e2rc01*e3rc22-2*e0rc00*e2rc00*e3rc22+2*e0rc21*e2rc22*e3rc21+2*e0rc22*e2rc21*e3rc21+2*e0rc11*e2rc12*e3rc21+2*e0rc12*e2rc11*e3rc21+2*e0rc01*e2rc02*e3rc21+2*e0rc02*e2rc01*e3rc21+2*e0rc20*e2rc22*e3rc20+2*e0rc22*e2rc20*e3rc20+2*e0rc10*e2rc12*e3rc20+2*e0rc12*e2rc10*e3rc20+2*e0rc00*e2rc02*e3rc20+2*e0rc02*e2rc00*e3rc20+2*e0rc12*e2rc22*e3rc12+2*e0rc11*e2rc21*e3rc12+2*e0rc10*e2rc20*e3rc12+2*e0rc22*e2rc12*e3rc12+2*e0rc21*e2rc11*e3rc12+2*e0rc20*e2rc10*e3rc12-2*e0rc11*e2rc22*e3rc11+2*e0rc12*e2rc21*e3rc11+2*e0rc21*e2rc12*e3rc11-2*e0rc22*e2rc11*e3rc11-2*e0rc10*e2rc22*e3rc10+2*e0rc12*e2rc20*e3rc10+2*e0rc20*e2rc12*e3rc10-2*e0rc22*e2rc10*e3rc10+2*e0rc02*e2rc22*e3rc02+2*e0rc01*e2rc21*e3rc02+2*e0rc00*e2rc20*e3rc02+2*e0rc22*e2rc02*e3rc02+2*e0rc21*e2rc01*e3rc02+2*e0rc20*e2rc00*e3rc02-2*e0rc01*e2rc22*e3rc01+2*e0rc02*e2rc21*e3rc01+2*e0rc21*e2rc02*e3rc01-2*e0rc22*e2rc01*e3rc01-2*e0rc00*e2rc22*e3rc00+2*e0rc02*e2rc20*e3rc00+2*e0rc20*e2rc02*e3rc00-2*e0rc22*e2rc00*e3rc00;
    AMatrix(8, 13) = +3*e1rc22_2*e3rc22+e1rc21_2*e3rc22+e1rc20_2*e3rc22+e1rc12_2*e3rc22-e1rc11_2*e3rc22-e1rc10_2*e3rc22+e1rc02_2*e3rc22-e1rc01_2*e3rc22-e1rc00_2*e3rc22+2*e1rc21*e1rc22*e3rc21+2*e1rc11*e1rc12*e3rc21+2*e1rc01*e1rc02*e3rc21+2*e1rc20*e1rc22*e3rc20+2*e1rc10*e1rc12*e3rc20+2*e1rc00*e1rc02*e3rc20+2*e1rc12*e1rc22*e3rc12+2*e1rc11*e1rc21*e3rc12+2*e1rc10*e1rc20*e3rc12-2*e1rc11*e1rc22*e3rc11+2*e1rc12*e1rc21*e3rc11-2*e1rc10*e1rc22*e3rc10+2*e1rc12*e1rc20*e3rc10+2*e1rc02*e1rc22*e3rc02+2*e1rc01*e1rc21*e3rc02+2*e1rc00*e1rc20*e3rc02-2*e1rc01*e1rc22*e3rc01+2*e1rc02*e1rc21*e3rc01-2*e1rc00*e1rc22*e3rc00+2*e1rc02*e1rc20*e3rc00;
    AMatrix(8, 14) = +6*e1rc22*e2rc22*e3rc22+2*e1rc21*e2rc21*e3rc22+2*e1rc20*e2rc20*e3rc22+2*e1rc12*e2rc12*e3rc22-2*e1rc11*e2rc11*e3rc22-2*e1rc10*e2rc10*e3rc22+2*e1rc02*e2rc02*e3rc22-2*e1rc01*e2rc01*e3rc22-2*e1rc00*e2rc00*e3rc22+2*e1rc21*e2rc22*e3rc21+2*e1rc22*e2rc21*e3rc21+2*e1rc11*e2rc12*e3rc21+2*e1rc12*e2rc11*e3rc21+2*e1rc01*e2rc02*e3rc21+2*e1rc02*e2rc01*e3rc21+2*e1rc20*e2rc22*e3rc20+2*e1rc22*e2rc20*e3rc20+2*e1rc10*e2rc12*e3rc20+2*e1rc12*e2rc10*e3rc20+2*e1rc00*e2rc02*e3rc20+2*e1rc02*e2rc00*e3rc20+2*e1rc12*e2rc22*e3rc12+2*e1rc11*e2rc21*e3rc12+2*e1rc10*e2rc20*e3rc12+2*e1rc22*e2rc12*e3rc12+2*e1rc21*e2rc11*e3rc12+2*e1rc20*e2rc10*e3rc12-2*e1rc11*e2rc22*e3rc11+2*e1rc12*e2rc21*e3rc11+2*e1rc21*e2rc12*e3rc11-2*e1rc22*e2rc11*e3rc11-2*e1rc10*e2rc22*e3rc10+2*e1rc12*e2rc20*e3rc10+2*e1rc20*e2rc12*e3rc10-2*e1rc22*e2rc10*e3rc10+2*e1rc02*e2rc22*e3rc02+2*e1rc01*e2rc21*e3rc02+2*e1rc00*e2rc20*e3rc02+2*e1rc22*e2rc02*e3rc02+2*e1rc21*e2rc01*e3rc02+2*e1rc20*e2rc00*e3rc02-2*e1rc01*e2rc22*e3rc01+2*e1rc02*e2rc21*e3rc01+2*e1rc21*e2rc02*e3rc01-2*e1rc22*e2rc01*e3rc01-2*e1rc00*e2rc22*e3rc00+2*e1rc02*e2rc20*e3rc00+2*e1rc20*e2rc02*e3rc00-2*e1rc22*e2rc00*e3rc00;
    AMatrix(8, 15) = +3*e2rc22_2*e3rc22+e2rc21_2*e3rc22+e2rc20_2*e3rc22+e2rc12_2*e3rc22-e2rc11_2*e3rc22-e2rc10_2*e3rc22+e2rc02_2*e3rc22-e2rc01_2*e3rc22-e2rc00_2*e3rc22+2*e2rc21*e2rc22*e3rc21+2*e2rc11*e2rc12*e3rc21+2*e2rc01*e2rc02*e3rc21+2*e2rc20*e2rc22*e3rc20+2*e2rc10*e2rc12*e3rc20+2*e2rc00*e2rc02*e3rc20+2*e2rc12*e2rc22*e3rc12+2*e2rc11*e2rc21*e3rc12+2*e2rc10*e2rc20*e3rc12-2*e2rc11*e2rc22*e3rc11+2*e2rc12*e2rc21*e3rc11-2*e2rc10*e2rc22*e3rc10+2*e2rc12*e2rc20*e3rc10+2*e2rc02*e2rc22*e3rc02+2*e2rc01*e2rc21*e3rc02+2*e2rc00*e2rc20*e3rc02-2*e2rc01*e2rc22*e3rc01+2*e2rc02*e2rc21*e3rc01-2*e2rc00*e2rc22*e3rc00+2*e2rc02*e2rc20*e3rc00;
    AMatrix(8, 16) = +3*e0rc22*e3rc22_2+2*e0rc21*e3rc21*e3rc22+2*e0rc20*e3rc20*e3rc22+2*e0rc12*e3rc12*e3rc22-2*e0rc11*e3rc11*e3rc22-2*e0rc10*e3rc10*e3rc22+2*e0rc02*e3rc02*e3rc22-2*e0rc01*e3rc01*e3rc22-2*e0rc00*e3rc00*e3rc22+e0rc22*e3rc21_2+2*e0rc11*e3rc12*e3rc21+2*e0rc12*e3rc11*e3rc21+2*e0rc01*e3rc02*e3rc21+2*e0rc02*e3rc01*e3rc21+e0rc22*e3rc20_2+2*e0rc10*e3rc12*e3rc20+2*e0rc12*e3rc10*e3rc20+2*e0rc00*e3rc02*e3rc20+2*e0rc02*e3rc00*e3rc20+e0rc22*e3rc12_2+2*e0rc21*e3rc11*e3rc12+2*e0rc20*e3rc10*e3rc12-e0rc22*e3rc11_2-e0rc22*e3rc10_2+e0rc22*e3rc02_2+2*e0rc21*e3rc01*e3rc02+2*e0rc20*e3rc00*e3rc02-e0rc22*e3rc01_2-e0rc22*e3rc00_2;
    AMatrix(8, 17) = +3*e1rc22*e3rc22_2+2*e1rc21*e3rc21*e3rc22+2*e1rc20*e3rc20*e3rc22+2*e1rc12*e3rc12*e3rc22-2*e1rc11*e3rc11*e3rc22-2*e1rc10*e3rc10*e3rc22+2*e1rc02*e3rc02*e3rc22-2*e1rc01*e3rc01*e3rc22-2*e1rc00*e3rc00*e3rc22+e1rc22*e3rc21_2+2*e1rc11*e3rc12*e3rc21+2*e1rc12*e3rc11*e3rc21+2*e1rc01*e3rc02*e3rc21+2*e1rc02*e3rc01*e3rc21+e1rc22*e3rc20_2+2*e1rc10*e3rc12*e3rc20+2*e1rc12*e3rc10*e3rc20+2*e1rc00*e3rc02*e3rc20+2*e1rc02*e3rc00*e3rc20+e1rc22*e3rc12_2+2*e1rc21*e3rc11*e3rc12+2*e1rc20*e3rc10*e3rc12-e1rc22*e3rc11_2-e1rc22*e3rc10_2+e1rc22*e3rc02_2+2*e1rc21*e3rc01*e3rc02+2*e1rc20*e3rc00*e3rc02-e1rc22*e3rc01_2-e1rc22*e3rc00_2;
    AMatrix(8, 18) = +3*e2rc22*e3rc22_2+2*e2rc21*e3rc21*e3rc22+2*e2rc20*e3rc20*e3rc22+2*e2rc12*e3rc12*e3rc22-2*e2rc11*e3rc11*e3rc22-2*e2rc10*e3rc10*e3rc22+2*e2rc02*e3rc02*e3rc22-2*e2rc01*e3rc01*e3rc22-2*e2rc00*e3rc00*e3rc22+e2rc22*e3rc21_2+2*e2rc11*e3rc12*e3rc21+2*e2rc12*e3rc11*e3rc21+2*e2rc01*e3rc02*e3rc21+2*e2rc02*e3rc01*e3rc21+e2rc22*e3rc20_2+2*e2rc10*e3rc12*e3rc20+2*e2rc12*e3rc10*e3rc20+2*e2rc00*e3rc02*e3rc20+2*e2rc02*e3rc00*e3rc20+e2rc22*e3rc12_2+2*e2rc21*e3rc11*e3rc12+2*e2rc20*e3rc10*e3rc12-e2rc22*e3rc11_2-e2rc22*e3rc10_2+e2rc22*e3rc02_2+2*e2rc21*e3rc01*e3rc02+2*e2rc20*e3rc00*e3rc02-e2rc22*e3rc01_2-e2rc22*e3rc00_2;
    AMatrix(8, 19) = +e3rc22_3+e3rc21_2*e3rc22+e3rc20_2*e3rc22+e3rc12_2*e3rc22-e3rc11_2*e3rc22-e3rc10_2*e3rc22+e3rc02_2*e3rc22-e3rc01_2*e3rc22-e3rc00_2*e3rc22+2*e3rc11*e3rc12*e3rc21+2*e3rc01*e3rc02*e3rc21+2*e3rc10*e3rc12*e3rc20+2*e3rc00*e3rc02*e3rc20;
    AMatrix(9, 0) = +e0rc00*e0rc11*e0rc22-e0rc01*e0rc10*e0rc22-e0rc00*e0rc12*e0rc21+e0rc02*e0rc10*e0rc21+e0rc01*e0rc12*e0rc20-e0rc02*e0rc11*e0rc20;
    AMatrix(9, 1) = +e0rc00*e0rc11*e1rc22-e0rc01*e0rc10*e1rc22-e0rc00*e0rc12*e1rc21+e0rc02*e0rc10*e1rc21+e0rc01*e0rc12*e1rc20-e0rc02*e0rc11*e1rc20-e0rc00*e0rc21*e1rc12+e0rc01*e0rc20*e1rc12+e0rc00*e0rc22*e1rc11-e0rc02*e0rc20*e1rc11-e0rc01*e0rc22*e1rc10+e0rc02*e0rc21*e1rc10+e0rc10*e0rc21*e1rc02-e0rc11*e0rc20*e1rc02-e0rc10*e0rc22*e1rc01+e0rc12*e0rc20*e1rc01+e0rc11*e0rc22*e1rc00-e0rc12*e0rc21*e1rc00;
    AMatrix(9, 2) = +e0rc00*e0rc11*e2rc22-e0rc01*e0rc10*e2rc22-e0rc00*e0rc12*e2rc21+e0rc02*e0rc10*e2rc21+e0rc01*e0rc12*e2rc20-e0rc02*e0rc11*e2rc20-e0rc00*e0rc21*e2rc12+e0rc01*e0rc20*e2rc12+e0rc00*e0rc22*e2rc11-e0rc02*e0rc20*e2rc11-e0rc01*e0rc22*e2rc10+e0rc02*e0rc21*e2rc10+e0rc10*e0rc21*e2rc02-e0rc11*e0rc20*e2rc02-e0rc10*e0rc22*e2rc01+e0rc12*e0rc20*e2rc01+e0rc11*e0rc22*e2rc00-e0rc12*e0rc21*e2rc00;
    AMatrix(9, 3) = +e0rc00*e1rc11*e1rc22-e0rc01*e1rc10*e1rc22-e0rc10*e1rc01*e1rc22+e0rc11*e1rc00*e1rc22-e0rc00*e1rc12*e1rc21+e0rc02*e1rc10*e1rc21+e0rc10*e1rc02*e1rc21-e0rc12*e1rc00*e1rc21+e0rc01*e1rc12*e1rc20-e0rc02*e1rc11*e1rc20-e0rc11*e1rc02*e1rc20+e0rc12*e1rc01*e1rc20+e0rc20*e1rc01*e1rc12-e0rc21*e1rc00*e1rc12-e0rc20*e1rc02*e1rc11+e0rc22*e1rc00*e1rc11+e0rc21*e1rc02*e1rc10-e0rc22*e1rc01*e1rc10;
    AMatrix(9, 4) = +e0rc00*e1rc11*e2rc22-e0rc01*e1rc10*e2rc22-e0rc10*e1rc01*e2rc22+e0rc11*e1rc00*e2rc22-e0rc00*e1rc12*e2rc21+e0rc02*e1rc10*e2rc21+e0rc10*e1rc02*e2rc21-e0rc12*e1rc00*e2rc21+e0rc01*e1rc12*e2rc20-e0rc02*e1rc11*e2rc20-e0rc11*e1rc02*e2rc20+e0rc12*e1rc01*e2rc20-e0rc00*e1rc21*e2rc12+e0rc01*e1rc20*e2rc12+e0rc20*e1rc01*e2rc12-e0rc21*e1rc00*e2rc12+e0rc00*e1rc22*e2rc11-e0rc02*e1rc20*e2rc11-e0rc20*e1rc02*e2rc11+e0rc22*e1rc00*e2rc11-e0rc01*e1rc22*e2rc10+e0rc02*e1rc21*e2rc10+e0rc21*e1rc02*e2rc10-e0rc22*e1rc01*e2rc10+e0rc10*e1rc21*e2rc02-e0rc11*e1rc20*e2rc02-e0rc20*e1rc11*e2rc02+e0rc21*e1rc10*e2rc02-e0rc10*e1rc22*e2rc01+e0rc12*e1rc20*e2rc01+e0rc20*e1rc12*e2rc01-e0rc22*e1rc10*e2rc01+e0rc11*e1rc22*e2rc00-e0rc12*e1rc21*e2rc00-e0rc21*e1rc12*e2rc00+e0rc22*e1rc11*e2rc00;
    AMatrix(9, 5) = +e0rc00*e2rc11*e2rc22-e0rc01*e2rc10*e2rc22-e0rc10*e2rc01*e2rc22+e0rc11*e2rc00*e2rc22-e0rc00*e2rc12*e2rc21+e0rc02*e2rc10*e2rc21+e0rc10*e2rc02*e2rc21-e0rc12*e2rc00*e2rc21+e0rc01*e2rc12*e2rc20-e0rc02*e2rc11*e2rc20-e0rc11*e2rc02*e2rc20+e0rc12*e2rc01*e2rc20+e0rc20*e2rc01*e2rc12-e0rc21*e2rc00*e2rc12-e0rc20*e2rc02*e2rc11+e0rc22*e2rc00*e2rc11+e0rc21*e2rc02*e2rc10-e0rc22*e2rc01*e2rc10;
    AMatrix(9, 6) = +e1rc00*e1rc11*e1rc22-e1rc01*e1rc10*e1rc22-e1rc00*e1rc12*e1rc21+e1rc02*e1rc10*e1rc21+e1rc01*e1rc12*e1rc20-e1rc02*e1rc11*e1rc20;
    AMatrix(9, 7) = +e1rc00*e1rc11*e2rc22-e1rc01*e1rc10*e2rc22-e1rc00*e1rc12*e2rc21+e1rc02*e1rc10*e2rc21+e1rc01*e1rc12*e2rc20-e1rc02*e1rc11*e2rc20-e1rc00*e1rc21*e2rc12+e1rc01*e1rc20*e2rc12+e1rc00*e1rc22*e2rc11-e1rc02*e1rc20*e2rc11-e1rc01*e1rc22*e2rc10+e1rc02*e1rc21*e2rc10+e1rc10*e1rc21*e2rc02-e1rc11*e1rc20*e2rc02-e1rc10*e1rc22*e2rc01+e1rc12*e1rc20*e2rc01+e1rc11*e1rc22*e2rc00-e1rc12*e1rc21*e2rc00;
    AMatrix(9, 8) = +e1rc00*e2rc11*e2rc22-e1rc01*e2rc10*e2rc22-e1rc10*e2rc01*e2rc22+e1rc11*e2rc00*e2rc22-e1rc00*e2rc12*e2rc21+e1rc02*e2rc10*e2rc21+e1rc10*e2rc02*e2rc21-e1rc12*e2rc00*e2rc21+e1rc01*e2rc12*e2rc20-e1rc02*e2rc11*e2rc20-e1rc11*e2rc02*e2rc20+e1rc12*e2rc01*e2rc20+e1rc20*e2rc01*e2rc12-e1rc21*e2rc00*e2rc12-e1rc20*e2rc02*e2rc11+e1rc22*e2rc00*e2rc11+e1rc21*e2rc02*e2rc10-e1rc22*e2rc01*e2rc10;
    AMatrix(9, 9) = e2rc00*e2rc11*e2rc22-e2rc01*e2rc10*e2rc22-e2rc00*e2rc12*e2rc21+e2rc02*e2rc10*e2rc21+e2rc01*e2rc12*e2rc20-e2rc02*e2rc11*e2rc20;
    AMatrix(9, 10) = +e0rc00*e0rc11*e3rc22-e0rc01*e0rc10*e3rc22-e0rc00*e0rc12*e3rc21+e0rc02*e0rc10*e3rc21+e0rc01*e0rc12*e3rc20-e0rc02*e0rc11*e3rc20-e0rc00*e0rc21*e3rc12+e0rc01*e0rc20*e3rc12+e0rc00*e0rc22*e3rc11-e0rc02*e0rc20*e3rc11-e0rc01*e0rc22*e3rc10+e0rc02*e0rc21*e3rc10+e0rc10*e0rc21*e3rc02-e0rc11*e0rc20*e3rc02-e0rc10*e0rc22*e3rc01+e0rc12*e0rc20*e3rc01+e0rc11*e0rc22*e3rc00-e0rc12*e0rc21*e3rc00;
    AMatrix(9, 11) = +e0rc00*e1rc11*e3rc22-e0rc01*e1rc10*e3rc22-e0rc10*e1rc01*e3rc22+e0rc11*e1rc00*e3rc22-e0rc00*e1rc12*e3rc21+e0rc02*e1rc10*e3rc21+e0rc10*e1rc02*e3rc21-e0rc12*e1rc00*e3rc21+e0rc01*e1rc12*e3rc20-e0rc02*e1rc11*e3rc20-e0rc11*e1rc02*e3rc20+e0rc12*e1rc01*e3rc20-e0rc00*e1rc21*e3rc12+e0rc01*e1rc20*e3rc12+e0rc20*e1rc01*e3rc12-e0rc21*e1rc00*e3rc12+e0rc00*e1rc22*e3rc11-e0rc02*e1rc20*e3rc11-e0rc20*e1rc02*e3rc11+e0rc22*e1rc00*e3rc11-e0rc01*e1rc22*e3rc10+e0rc02*e1rc21*e3rc10+e0rc21*e1rc02*e3rc10-e0rc22*e1rc01*e3rc10+e0rc10*e1rc21*e3rc02-e0rc11*e1rc20*e3rc02-e0rc20*e1rc11*e3rc02+e0rc21*e1rc10*e3rc02-e0rc10*e1rc22*e3rc01+e0rc12*e1rc20*e3rc01+e0rc20*e1rc12*e3rc01-e0rc22*e1rc10*e3rc01+e0rc11*e1rc22*e3rc00-e0rc12*e1rc21*e3rc00-e0rc21*e1rc12*e3rc00+e0rc22*e1rc11*e3rc00;
    AMatrix(9, 12) = +e0rc00*e2rc11*e3rc22-e0rc01*e2rc10*e3rc22-e0rc10*e2rc01*e3rc22+e0rc11*e2rc00*e3rc22-e0rc00*e2rc12*e3rc21+e0rc02*e2rc10*e3rc21+e0rc10*e2rc02*e3rc21-e0rc12*e2rc00*e3rc21+e0rc01*e2rc12*e3rc20-e0rc02*e2rc11*e3rc20-e0rc11*e2rc02*e3rc20+e0rc12*e2rc01*e3rc20-e0rc00*e2rc21*e3rc12+e0rc01*e2rc20*e3rc12+e0rc20*e2rc01*e3rc12-e0rc21*e2rc00*e3rc12+e0rc00*e2rc22*e3rc11-e0rc02*e2rc20*e3rc11-e0rc20*e2rc02*e3rc11+e0rc22*e2rc00*e3rc11-e0rc01*e2rc22*e3rc10+e0rc02*e2rc21*e3rc10+e0rc21*e2rc02*e3rc10-e0rc22*e2rc01*e3rc10+e0rc10*e2rc21*e3rc02-e0rc11*e2rc20*e3rc02-e0rc20*e2rc11*e3rc02+e0rc21*e2rc10*e3rc02-e0rc10*e2rc22*e3rc01+e0rc12*e2rc20*e3rc01+e0rc20*e2rc12*e3rc01-e0rc22*e2rc10*e3rc01+e0rc11*e2rc22*e3rc00-e0rc12*e2rc21*e3rc00-e0rc21*e2rc12*e3rc00+e0rc22*e2rc11*e3rc00;
    AMatrix(9, 13) = +e1rc00*e1rc11*e3rc22-e1rc01*e1rc10*e3rc22-e1rc00*e1rc12*e3rc21+e1rc02*e1rc10*e3rc21+e1rc01*e1rc12*e3rc20-e1rc02*e1rc11*e3rc20-e1rc00*e1rc21*e3rc12+e1rc01*e1rc20*e3rc12+e1rc00*e1rc22*e3rc11-e1rc02*e1rc20*e3rc11-e1rc01*e1rc22*e3rc10+e1rc02*e1rc21*e3rc10+e1rc10*e1rc21*e3rc02-e1rc11*e1rc20*e3rc02-e1rc10*e1rc22*e3rc01+e1rc12*e1rc20*e3rc01+e1rc11*e1rc22*e3rc00-e1rc12*e1rc21*e3rc00;
    AMatrix(9, 14) = +e1rc00*e2rc11*e3rc22-e1rc01*e2rc10*e3rc22-e1rc10*e2rc01*e3rc22+e1rc11*e2rc00*e3rc22-e1rc00*e2rc12*e3rc21+e1rc02*e2rc10*e3rc21+e1rc10*e2rc02*e3rc21-e1rc12*e2rc00*e3rc21+e1rc01*e2rc12*e3rc20-e1rc02*e2rc11*e3rc20-e1rc11*e2rc02*e3rc20+e1rc12*e2rc01*e3rc20-e1rc00*e2rc21*e3rc12+e1rc01*e2rc20*e3rc12+e1rc20*e2rc01*e3rc12-e1rc21*e2rc00*e3rc12+e1rc00*e2rc22*e3rc11-e1rc02*e2rc20*e3rc11-e1rc20*e2rc02*e3rc11+e1rc22*e2rc00*e3rc11-e1rc01*e2rc22*e3rc10+e1rc02*e2rc21*e3rc10+e1rc21*e2rc02*e3rc10-e1rc22*e2rc01*e3rc10+e1rc10*e2rc21*e3rc02-e1rc11*e2rc20*e3rc02-e1rc20*e2rc11*e3rc02+e1rc21*e2rc10*e3rc02-e1rc10*e2rc22*e3rc01+e1rc12*e2rc20*e3rc01+e1rc20*e2rc12*e3rc01-e1rc22*e2rc10*e3rc01+e1rc11*e2rc22*e3rc00-e1rc12*e2rc21*e3rc00-e1rc21*e2rc12*e3rc00+e1rc22*e2rc11*e3rc00;
    AMatrix(9, 15) = +e2rc00*e2rc11*e3rc22-e2rc01*e2rc10*e3rc22-e2rc00*e2rc12*e3rc21+e2rc02*e2rc10*e3rc21+e2rc01*e2rc12*e3rc20-e2rc02*e2rc11*e3rc20-e2rc00*e2rc21*e3rc12+e2rc01*e2rc20*e3rc12+e2rc00*e2rc22*e3rc11-e2rc02*e2rc20*e3rc11-e2rc01*e2rc22*e3rc10+e2rc02*e2rc21*e3rc10+e2rc10*e2rc21*e3rc02-e2rc11*e2rc20*e3rc02-e2rc10*e2rc22*e3rc01+e2rc12*e2rc20*e3rc01+e2rc11*e2rc22*e3rc00-e2rc12*e2rc21*e3rc00;
    AMatrix(9, 16) = +e0rc00*e3rc11*e3rc22-e0rc01*e3rc10*e3rc22-e0rc10*e3rc01*e3rc22+e0rc11*e3rc00*e3rc22-e0rc00*e3rc12*e3rc21+e0rc02*e3rc10*e3rc21+e0rc10*e3rc02*e3rc21-e0rc12*e3rc00*e3rc21+e0rc01*e3rc12*e3rc20-e0rc02*e3rc11*e3rc20-e0rc11*e3rc02*e3rc20+e0rc12*e3rc01*e3rc20+e0rc20*e3rc01*e3rc12-e0rc21*e3rc00*e3rc12-e0rc20*e3rc02*e3rc11+e0rc22*e3rc00*e3rc11+e0rc21*e3rc02*e3rc10-e0rc22*e3rc01*e3rc10;
    AMatrix(9, 17) = +e1rc00*e3rc11*e3rc22-e1rc01*e3rc10*e3rc22-e1rc10*e3rc01*e3rc22+e1rc11*e3rc00*e3rc22-e1rc00*e3rc12*e3rc21+e1rc02*e3rc10*e3rc21+e1rc10*e3rc02*e3rc21-e1rc12*e3rc00*e3rc21+e1rc01*e3rc12*e3rc20-e1rc02*e3rc11*e3rc20-e1rc11*e3rc02*e3rc20+e1rc12*e3rc01*e3rc20+e1rc20*e3rc01*e3rc12-e1rc21*e3rc00*e3rc12-e1rc20*e3rc02*e3rc11+e1rc22*e3rc00*e3rc11+e1rc21*e3rc02*e3rc10-e1rc22*e3rc01*e3rc10;
    AMatrix(9, 18) = +e2rc00*e3rc11*e3rc22-e2rc01*e3rc10*e3rc22-e2rc10*e3rc01*e3rc22+e2rc11*e3rc00*e3rc22-e2rc00*e3rc12*e3rc21+e2rc02*e3rc10*e3rc21+e2rc10*e3rc02*e3rc21-e2rc12*e3rc00*e3rc21+e2rc01*e3rc12*e3rc20-e2rc02*e3rc11*e3rc20-e2rc11*e3rc02*e3rc20+e2rc12*e3rc01*e3rc20+e2rc20*e3rc01*e3rc12-e2rc21*e3rc00*e3rc12-e2rc20*e3rc02*e3rc11+e2rc22*e3rc00*e3rc11+e2rc21*e3rc02*e3rc10-e2rc22*e3rc01*e3rc10;
    AMatrix(9, 19) = +e3rc00*e3rc11*e3rc22-e3rc01*e3rc10*e3rc22-e3rc00*e3rc12*e3rc21+e3rc02*e3rc10*e3rc21+e3rc01*e3rc12*e3rc20-e3rc02*e3rc11*e3rc20;


      // ============= End cut & paste section =============

      // Warning(xxx): The action matrix in Stewenius & Nister's paper
      // doesn't appear to correspond to degree-then-lexicographic
      // order of monomials.  The right solution to this is to
      // generate the correct action matrix in fivePointAlgorithm.h,
      // but instead we temporarily shuffle the columns of our
      // constraint matrix to match the action matrix.
      brick::numeric::Array2D<FloatType> A2Matrix(AMatrix.rows(), AMatrix.columns());
      brick::numeric::Array1D<int> shuffle(
        "[0, 1, 3, 6, 2, 4, 7, 5, 8, 9, 10, 11, 13, 12, 14, 15, 16, 17, 18, 19]"
        );
      for(size_t rowIndex = 0; rowIndex < AMatrix.rows(); ++rowIndex) {
        for(size_t columnIndex = 0; columnIndex < AMatrix.columns();
            ++columnIndex) {
          A2Matrix(rowIndex, columnIndex) =
            AMatrix(rowIndex, shuffle[columnIndex]);
        }
      }
      return A2Matrix;
    }

  } // namespace computerVision

} // namespace brick


#endif /* #ifndef BRICK_COMPUTERVISION_FIVEPOINTALGORITHM_IMPL_HH */
