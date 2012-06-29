/**
***************************************************************************
* @file brick/computerVision/registerPoints3D_impl.hh
*
* Header file defining inline and template functions from
* registerPoints3D.hh
*
* Copyright (C) 1998-2007,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_COMPUTERVISION_REGISTERPOINTS3D_IMPL_HH
#define BRICK_COMPUTERVISION_REGISTERPOINTS3D_IMPL_HH

// This file is included by registerPoints3D.hh, and should not be
// directly included by user code, so no need to include
// registerPoints3D.hh here.
// 
// #include <brick/computerVision/registerPoints3D.hh>

#include <limits>
#include <brick/common/functional.hh>
#include <brick/common/types.hh>
#include <brick/numeric/array1D.hh>
#include <brick/numeric/array2D.hh>
#include <brick/numeric/quaternion.hh>
#include <brick/numeric/rotations.hh>
#include <brick/numeric/utilities.hh>
#include <brick/numeric/vector3D.hh>
#include <brick/linearAlgebra/linearAlgebra.hh>

namespace brick {

  namespace computerVision {
    
    /// @cond privateCode
    namespace privateCode {

      // Private routine to generate the element selection array
      // proposed in Horn's paper.
      template <class FloatType>
      inline const Array2D<FloatType>&
      getHornNCoefficientMatrix()
      {
        static Array2D<FloatType> NCoefficients(
          "[[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], "
          " [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0], "
          " [0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0], "
          " [0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0], "
          " [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0], "
          " [1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0], "
          " [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0], "
          " [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0], "
          " [0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0], "
          " [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0], "
          " [-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0], "
          " [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0], "
          " [0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0], "
          " [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0], "
          " [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0], "
          " [-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0]]");
        return NCoefficients;
      }


      // Private routine containing code that is common to all flavors
      // of estimateTransform3D().
      template <class FloatType>
      brick::numeric::Transform3D<FloatType>
      estimateTransformFromZeroMeanPointArrays(
        Array2D<FloatType> const& translatedFromPoints,
        Array2D<FloatType> const& translatedToPoints,
        Vector3D<FloatType> const& fromMean,
        Vector3D<FloatType> const& toMean)
      {
        // Compute the (not quite) covariance matrix between the two point
        // clouds.  This is a  3x3 matrix.
        Array2D<FloatType> matrixM = matrixMultiply<FloatType>(
          translatedFromPoints.transpose(), translatedToPoints);

        // Select elements as describe by Horn.  The member function
        // ravel() simply returns a flattened (1D) array referencing the
        // elements of the Array2D instance in row-major order.
        Array1D<FloatType> vectorN = matrixMultiply<FloatType>(
          privateCode::getHornNCoefficientMatrix<FloatType>(), matrixM.ravel());
        Array2D<FloatType> matrixN(
          4, 4, vectorN.data(), vectorN.getReferenceCount());

        // Find the largest eigenvector of matrixN.  This is a unit
        // quaternion describing the best fit rotation.
        Array1D<brick::common::Float64> eValues;
        Array2D<brick::common::Float64> eVectors;

        // // Note(xxx): Need to implement more versions of
        // // eigenvectorsSymmetric()!
        // 
        // Array2D<Float64> tempArray(matrixN.rows(), matrixN.columns());
        // tempArray.copy(matrixN);
        // eigenvectorsSymmetric(matrixN, eValues, eVectors);
        brick::linearAlgebra::eigenvectorsSymmetric(matrixN, eValues, eVectors);

        // The function argmax returns the index of the largest element of
        // eValues.
        size_t index0 = argmax(eValues);
        Quaternion<FloatType> q0(static_cast<FloatType>(eVectors(0, index0)),
                                 static_cast<FloatType>(eVectors(1, index0)),
                                 static_cast<FloatType>(eVectors(2, index0)),
                                 static_cast<FloatType>(eVectors(3, index0)));

        // Convert the unit quaternion to a rotation matrix.
        brick::numeric::Transform3D<FloatType> xf = quaternionToTransform3D(q0);

        // And add in best fit translation.
        Vector3D<FloatType> bestFitTranslation = toMean - xf * fromMean;
        xf.setValue(0, 3, bestFitTranslation.x());
        xf.setValue(1, 3, bestFitTranslation.y());
        xf.setValue(2, 3, bestFitTranslation.z());
        return xf;
      }        
      
    } // namespace privateCode
    /// @endcond
  

    template <class FloatType, class InIter0, class InIter1>
    brick::numeric::Transform3D<FloatType>
    registerPoints3D(InIter0 fromPointsBegin, InIter0 fromPointsEnd,
                     InIter1 toPointsBegin)
    {
      std::vector<bool> flags(fromPointsEnd - fromPointsBegin, true);
      return registerPoints3D<FloatType>(
        fromPointsBegin, fromPointsEnd, toPointsBegin, flags.begin());
    }

  
    template <class FloatType, class InIter0, class InIter1, class InIter2>
    brick::numeric::Transform3D<FloatType>
    registerPoints3D(InIter0 fromPointsBegin, InIter0 fromPointsEnd,
                     InIter1 toPointsBegin, InIter2 flagsBegin)
    {
      // Compute the mean of each point cloud, so we can translate its
      // center to the origin.
      InIter0 fromIter = fromPointsBegin;
      InIter1 toIter = toPointsBegin;
      InIter2 flagsIter = flagsBegin;
      Vector3D<FloatType> fromMean(0.0, 0.0, 0.0);
      Vector3D<FloatType> toMean(0.0, 0.0, 0.0);
      size_t count = 0;
      while(fromIter != fromPointsEnd) {
        if(*flagsIter) {
          fromMean += *fromIter;
          toMean += *toIter;
          ++count;
        }
        ++fromIter;
        ++toIter;
        ++flagsIter;
      }
      if(count == 0) {
        BRICK_THROW(brick::common::ValueException, "registerPoints3D()",
                  "No points to register!");
      }
      fromMean /= static_cast<FloatType>(count);
      toMean /= static_cast<FloatType>(count);

      // Now translate each point cloud so that its center of mass is at
      // the origin.  We'll use these translated point clouds to compute
      // the best fit rotation.  We use Nx3 arrays to represent these
      // translated ponits so that we'll be able to conveniently do
      // linear algebra later.
      fromIter = fromPointsBegin;
      toIter = toPointsBegin;
      flagsIter = flagsBegin;
      Array2D<FloatType> translatedFromPoints(count, 3);
      Array2D<FloatType> translatedToPoints(count, 3);
      size_t arrayIndex = 0;
      while(fromIter != fromPointsEnd) {
        if(*flagsIter) {
          Vector3D<FloatType>& fromPoint = *fromIter;
          Vector3D<FloatType>& toPoint = *toIter;
        
          // Copy the first element in the current row.
          translatedFromPoints[arrayIndex] = fromPoint.x() - fromMean.x();
          translatedToPoints[arrayIndex] = toPoint.x() - toMean.x();
        
          // Advance to next element in the current row.
          ++arrayIndex;
          translatedFromPoints[arrayIndex] = fromPoint.y() - fromMean.y();
          translatedToPoints[arrayIndex] = toPoint.y() - toMean.y();
        
          // Advance to next element in the current row.
          ++arrayIndex;
          translatedFromPoints[arrayIndex] = fromPoint.z() - fromMean.z();
          translatedToPoints[arrayIndex] = toPoint.z() - toMean.z();
        
          // Wrap around to next row.
          ++arrayIndex;
        }
        // Advance to next point.
        ++fromIter;
        ++toIter;
        ++flagsIter;
      }

      // Now that the input matrices are constructed, dispatch to a
      // subroutine for the actual transform estimation.
      return privateCode::estimateTransformFromZeroMeanPointArrays(
        translatedFromPoints, translatedToPoints, fromMean, toMean);
    }
  

    template <class FloatType, class InIter0, class InIter1, class InIter2>
    brick::numeric::Transform3D<FloatType>
    registerPoints3D(InIter0 fromPointsBegin, InIter0 fromPointsEnd,
                     InIter1 toPointsBegin, InIter2 weightsBegin,
                     bool /* dummy */)
    {
      // Compute the mean of each point cloud, so we can translate its
      // center to the origin.
      InIter0 fromIter = fromPointsBegin;
      InIter1 toIter = toPointsBegin;
      InIter2 weightsIter = weightsBegin;
      Vector3D<FloatType> fromMean(0.0, 0.0, 0.0);
      Vector3D<FloatType> toMean(0.0, 0.0, 0.0);
      FloatType totalWeight = 0.0;
      unsigned int count = 0;
      while(fromIter != fromPointsEnd) {
        FloatType weight = *weightsIter;
        fromMean += weight * (*fromIter);
        toMean += weight * (*toIter);
        totalWeight += weight;
        ++count;
        ++fromIter;
        ++toIter;
        ++weightsIter;
      }
      if((count == 0) || (totalWeight == 0.0)) {
        BRICK_THROW(brick::common::ValueException, "registerPoints3D()",
                  "No weighted points to register!");
      }
      fromMean /= totalWeight;
      toMean /= totalWeight;

      // Now translate each point cloud so that its center of mass is at
      // the origin.  We'll use these translated point clouds to compute
      // the best fit rotation.  We use Nx3 arrays to represent these
      // translated ponits so that we'll be able to conveniently do
      // linear algebra later.
      fromIter = fromPointsBegin;
      toIter = toPointsBegin;
      weightsIter = weightsBegin;
      Array2D<FloatType> translatedFromPoints(count, 3);
      Array2D<FloatType> translatedToPoints(count, 3);
      size_t arrayIndex = 0;
      while(fromIter != fromPointsEnd) {
        FloatType weight = *weightsIter;
        Vector3D<FloatType> const& fromPoint = *fromIter;
        Vector3D<FloatType> const& toPoint = *toIter;
        
        // Copy the first element in the current row.  We know that,
        // later on, the thing we care about is the matrix product of
        // translatedFromPoints.transpose() and translatedToPoints.
        // This means that, after subtracting out the mean, we can
        // multiply fromPoint and toPoint by the square root of the
        // weight.  Equivalently, we can multiply only one of
        // fromPoint and toPoint by weight, and achieve the same
        // effect without doing a sqrt operation.
        translatedFromPoints[arrayIndex] =
          weight * (fromPoint.x() - fromMean.x());
        translatedToPoints[arrayIndex] = toPoint.x() - toMean.x();
        
        // Advance to next element in the current row.
        ++arrayIndex;
        translatedFromPoints[arrayIndex] = 
          weight * (fromPoint.y() - fromMean.y());
        translatedToPoints[arrayIndex] = toPoint.y() - toMean.y();
        
        // Advance to next element in the current row.
        ++arrayIndex;
        translatedFromPoints[arrayIndex] =
          weight * (fromPoint.z() - fromMean.z());
        translatedToPoints[arrayIndex] = toPoint.z() - toMean.z();
        
        // Wrap around to next row.
        ++arrayIndex;

        // Advance to next point.
        ++fromIter;
        ++toIter;
        ++weightsIter;
      }

      // Now that the input matrices are constructed, dispatch to a
      // subroutine for the actual transform estimation.
      return privateCode::estimateTransformFromZeroMeanPointArrays(
        translatedFromPoints, translatedToPoints, fromMean, toMean);
    }
  

    template <class FloatType, class InIter0, class InIter1, class OutIter0>
    brick::numeric::Transform3D<FloatType>
    registerPoints3D(InIter0 fromPointsBegin, InIter0 fromPointsEnd,
                     InIter1 toPointsBegin, OutIter0 selectedFlagsBegin,
                     FloatType inclusion, FloatType maximumResidual,
                     size_t maximumIterations)
    {
      // Sort out arguments, array sizes, etc.
      size_t numberOfPoints = fromPointsEnd - fromPointsBegin;
      size_t numberToExclude = 0;
      if(inclusion > 0.0 && inclusion < 1.0) {
        size_t numberToInclude = static_cast<size_t>(
          std::floor(inclusion * numberOfPoints));
        if(numberToInclude >= numberOfPoints) {
          BRICK_THROW(brick::common::ValueException, "registerPoints3D()",
                    "Inclusion percentage is too low... all points excluded.");
        }
        numberToExclude = numberOfPoints - numberToInclude;
      }

      // Register repeatedly until we converge.
      brick::numeric::Transform3D<FloatType> resultTransform;
      std::vector<bool> flags(numberOfPoints, true);
      for(size_t iterationIndex = 0;
          iterationIndex < maximumIterations;
          ++iterationIndex) {
        // Do a sub-registration.
        resultTransform = registerPoints3D<FloatType>(
          fromPointsBegin, fromPointsEnd, toPointsBegin, flags.begin());

        // Compute residuals.
        std::vector< Vector3D<FloatType> > transformedPoints(numberOfPoints);
        std::transform(fromPointsBegin, fromPointsEnd,
                       transformedPoints.begin(), resultTransform.getFunctor());

        std::vector<FloatType> squaredResiduals(numberOfPoints);
        typedef FloatType (*magVec3D)(const Vector3D<FloatType>&);
        std::transform(
          transformedPoints.begin(), transformedPoints.end(),
          toPointsBegin, squaredResiduals.begin(),
          composeFunctor_1_2(
            std::ptr_fun(static_cast<magVec3D>(magnitudeSquared)),
            std::minus< Vector3D<FloatType> >()));

        // Figure out how big a residual has to be before we ignore the
        // associated point.
        FloatType residualThreshold =
          (maximumResidual > 0.0)
          ? maximumResidual : std::numeric_limits<FloatType>::max();
        if(numberToExclude != 0) {
          // Find the numberToExclude largest residuals.
          std::vector<FloatType> sortedSquaredResiduals(numberToExclude);
          std::partial_sort_copy(
            squaredResiduals.begin(), squaredResiduals.end(),
            sortedSquaredResiduals.begin(), sortedSquaredResiduals.end(),
            std::greater<FloatType>());
          FloatType threshold = sortedSquaredResiduals[numberToExclude - 1];
          if(threshold < residualThreshold) {
            residualThreshold = threshold;
          }
        }

        // Check which points exceed the residual threshold and should
        // be counted as outliers.
        std::vector<bool> newFlags(numberOfPoints);
        std::transform(
          squaredResiduals.begin(), squaredResiduals.end(), newFlags.begin(),
          std::bind2nd(std::less<FloatType>(), residualThreshold));

        // Any change from last iteration?
        if(std::equal(flags.begin(), flags.end(), newFlags.begin())) {
          // No change.  We're done.
          break;
        }
      
        // Reset flags so the next iteration will exclude points which
        // exceed the threshold.
        std::copy(newFlags.begin(), newFlags.end(), flags.begin());
      } // for

      std::copy(flags.begin(), flags.end(), selectedFlagsBegin);
      return resultTransform;
    }
  
  } // namespace computerVision    

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_REGISTERPOINTS3D_IMPL_HH */
