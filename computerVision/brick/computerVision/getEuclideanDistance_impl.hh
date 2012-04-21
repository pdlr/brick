/**
***************************************************************************
* @file brick/computerVision/getEuclideanDistance_impl.hh
*
* Header file defining inline and template functions declared in
* getEuclideanDistance.hh.
*
* Copyright (C) 2006,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_EUCLIDEANDISTANCE_IMPL_HH
#define BRICK_COMPUTERVISION_EUCLIDEANDISTANCE_IMPL_HH

// This file is included by getEuclideanDistance.hh, and should not be
// directly included by user code, so no need to include
// getEuclideanDistance.hh here.
// 
// #include <brick/computerVision/getEuclideanDistance.hh>

#include <list>
#include <brick/computerVision/imageFormat.hh>
#include <brick/computerVision/imageIO.hh>
#include <brick/computerVision/image.hh>
#include <brick/numeric/index2D.hh>

namespace brick {

  namespace computerVision {
    
    template<class FloatType, ImageFormat FORMAT>
    brick::numeric::Array2D<FloatType>
    getEuclideanDistance(const Image<FORMAT>& inputImage,
                         size_t maxNumberOfPasses=10);


    template<class FloatType, ImageFormat FORMAT>
    brick::numeric::Array2D<FloatType>
    getEuclideanDistance(const Image<FORMAT>& inputImage,
                         size_t maxNumberOfPasses,
                         size_t& numberOfPassesUsed);
  
  } // namespace computerVision

} // namespace brick


/* ============ Definitions of inline & template functions ============ */


#include <cmath>
#include <brick/numeric/index2D.hh>

namespace brick {

  namespace computerVision {
    
    /// @cond privateCode
    namespace privateCode {

      template<class FloatType>
      inline bool
      eucDistPropagate(size_t toIndex,
                       size_t fromIndex,
                       brick::numeric::Array2D<FloatType>& distanceMap,
                       brick::numeric::Array2D<brick::numeric::Index2D>& referentMap,
                       int row,
                       int column)
      {
        // Inexpensive check to see if we _might_ need to update.
        if(distanceMap(fromIndex) < distanceMap(toIndex)) {
          brick::numeric::Index2D& referent = referentMap(fromIndex);
          FloatType deltaU = column - referent.getColumn();
          FloatType deltaV = row - referent.getRow();
          FloatType newDistance = deltaU * deltaU + deltaV * deltaV;
          if(newDistance < distanceMap(toIndex)) {
            distanceMap(toIndex) = newDistance;
            referentMap(toIndex) = referent;
            return true;
          }
        }
        return false;
      }

    } // namespace privateCode
    /// @endcond

  
    template<class FloatType, ImageFormat FORMAT>
    brick::numeric::Array2D<FloatType>
    getEuclideanDistance(const Image<FORMAT>& inputImage,
                         size_t maxNumberOfPasses)
    {
      size_t numberOfPassesUsed;
      return getEuclideanDistance<FloatType>(inputImage, maxNumberOfPasses,
                                             numberOfPassesUsed);
    }


    template<class FloatType, ImageFormat FORMAT>
    brick::numeric::Array2D<FloatType>
    getEuclideanDistance(const Image<FORMAT>& inputImage,
                         size_t maxNumberOfPasses,
                         size_t& numberOfPassesUsed)
    {
      brick::numeric::Array2D<FloatType> distanceMap(inputImage.rows(), inputImage.columns());
      brick::numeric::Array2D<brick::numeric::Index2D> referentMap(inputImage.rows(), inputImage.columns());

      // Initialize the distance matrix.
      FloatType maxSqDistance = (inputImage.rows() * inputImage.rows()
                              + inputImage.columns() * inputImage.columns());
      size_t index0 = 0;
      for(size_t row = 0; row < inputImage.rows(); ++row) {
        for(size_t column = 0; column < inputImage.columns(); ++column) {
          if(inputImage(index0)) {
            distanceMap(index0) = 0.0;
            referentMap(index0) = brick::numeric::Index2D(static_cast<int>(row), static_cast<int>(column));
          } else {
            distanceMap(index0) = maxSqDistance;
          }
          ++index0;
        }
      }

      // This is hard to do efficiently closed-form, so we'll iterate.
      // If numberOfPasses is set to 1, we'll still get a reasonable
      // map.
      size_t columns = distanceMap.columns();
      size_t rows = distanceMap.rows();
      size_t passNumber;
      for(passNumber = 0; passNumber < maxNumberOfPasses; ++passNumber) {
        // This variable will tell us if we can quit early because the
        // distances are all correct.
        bool isChanged = false;
      
        // === Propagate distances East. === 

        // Propagate distances East along top row.
        index0 = 1;
        for(size_t column = 1; column < inputImage.columns(); ++column) {
          isChanged |= privateCode::eucDistPropagate(
            index0, index0 - 1, distanceMap, referentMap,
            0, static_cast<int>(column));
          isChanged |= privateCode::eucDistPropagate(
            index0, index0 + columns - 1, distanceMap, referentMap,
            0, static_cast<int>(column));
          ++index0;
        }

        // Propagate distances East through the bulk of the image.
        for(size_t row = 1; row < inputImage.rows() - 1; ++row) {
          index0 = row * inputImage.columns() + 1;
          for(size_t column = 1; column < inputImage.columns(); ++column) {
            isChanged |= privateCode::eucDistPropagate(
              index0, index0 - columns - 1, distanceMap, referentMap,
              static_cast<int>(row), static_cast<int>(column));
            isChanged |= privateCode::eucDistPropagate(
              index0, index0 - 1, distanceMap, referentMap,
              static_cast<int>(row), static_cast<int>(column));
            isChanged |= privateCode::eucDistPropagate(
              index0, index0 + columns - 1, distanceMap, referentMap,
              static_cast<int>(row), static_cast<int>(column));
            ++index0;
          }
        }

        // Propagate distances East along bottom row.
        index0 = (inputImage.rows() - 1) * inputImage.columns() + 1;
        for(size_t column = 1; column < inputImage.columns(); ++column) {
          isChanged |= privateCode::eucDistPropagate(
            index0, index0 - columns - 1, distanceMap, referentMap,
            static_cast<int>(rows) - 1, static_cast<int>(column));
          isChanged |= privateCode::eucDistPropagate(
            index0, index0 - 1, distanceMap, referentMap,
            static_cast<int>(rows) - 1, static_cast<int>(column));
          ++index0;
        }

        // === Propagate distances South. ===
        index0 = columns + 1;
        for(size_t row = 1; row < inputImage.rows(); ++row) {
          index0 = row * inputImage.columns();

          // Update first pixel of the row.
          isChanged |= privateCode::eucDistPropagate(
            index0, index0 - columns, distanceMap, referentMap,
            static_cast<int>(row), 0);
          isChanged |= privateCode::eucDistPropagate(
            index0, index0 - columns + 1, distanceMap, referentMap,
            static_cast<int>(row), 0);
          ++index0;

          // Update the bulk of the pixels in the current row.
          for(size_t column = 1; column < inputImage.columns() - 1; ++column) {
            isChanged |= privateCode::eucDistPropagate(
              index0, index0 - columns - 1, distanceMap, referentMap,
              static_cast<int>(row), static_cast<int>(column));
            isChanged |= privateCode::eucDistPropagate(
              index0, index0 - columns, distanceMap, referentMap,
              static_cast<int>(row), static_cast<int>(column));
            isChanged |= privateCode::eucDistPropagate(
              index0, index0 - columns + 1, distanceMap, referentMap,
              static_cast<int>(row), static_cast<int>(column));
            ++index0;
          }

          // Update last pixel of the row.
          isChanged |= privateCode::eucDistPropagate(
            index0, index0 - columns - 1, distanceMap, referentMap,
            static_cast<int>(row), static_cast<int>(columns) - 1);
          isChanged |= privateCode::eucDistPropagate(
            index0, index0 - columns, distanceMap, referentMap,
            static_cast<int>(row), static_cast<int>(columns) - 1);
          ++index0;
        }

        // === Propagate distances West. ===

        // Propagate distances West along top row.
        index0 = columns - 2;

        // Trick here... note that column wraps around to a very large
        // number once instead of going negative.
        for(size_t column = inputImage.columns() - 2; column < columns;
            --column) {
          isChanged |= privateCode::eucDistPropagate(
            index0, index0 + 1, distanceMap, referentMap,
            0, static_cast<int>(column));
          isChanged |= privateCode::eucDistPropagate(
            index0, index0 + columns + 1, distanceMap, referentMap,
            0, static_cast<int>(column));
          --index0;
        }
      
        // Propagate distances West through the bulk of the image.
        for(size_t row = 1; row < inputImage.rows() - 1; ++row) {
          index0 = (row + 1) * inputImage.columns() - 2;
          for(size_t column = inputImage.columns() - 2; column < columns;
              --column) {
            isChanged |= privateCode::eucDistPropagate(
              index0, index0 - columns + 1, distanceMap, referentMap,
              static_cast<int>(row), static_cast<int>(column));
            isChanged |= privateCode::eucDistPropagate(
              index0, index0 + 1, distanceMap, referentMap,
              static_cast<int>(row), static_cast<int>(column));
            isChanged |= privateCode::eucDistPropagate(
              index0, index0 + columns + 1, distanceMap, referentMap,
              static_cast<int>(row), static_cast<int>(column));
            --index0;
          }
        }
      
        // Propagate distances West along bottom row.
        index0 = rows * columns - 2;
        for(size_t column = inputImage.columns() - 2; column < columns; --column) {
          isChanged |= privateCode::eucDistPropagate(
            index0, index0 - columns + 1, distanceMap, referentMap,
            static_cast<int>(rows - 1), static_cast<int>(column));
          isChanged |= privateCode::eucDistPropagate(
            index0, index0 + 1, distanceMap, referentMap,
            static_cast<int>(rows - 1), static_cast<int>(column));
          --index0;
        }


        // === Propagate distances North. ===

        // Remember that row wraps around to a very large number instead
        // of going negative.
        for(size_t row = inputImage.rows() - 2; row < rows; --row) {
          index0 = row * columns;

          // Update first pixel of the row.
          isChanged |= privateCode::eucDistPropagate(
            index0, index0 + columns, distanceMap, referentMap,
            static_cast<int>(row), 0);
          isChanged |= privateCode::eucDistPropagate(
            index0, index0 + columns + 1, distanceMap, referentMap,
            static_cast<int>(row), 0);
          ++index0;

          // Update the bulk of the pixels in the current row.
          for(size_t column = 1; column < inputImage.columns() - 1; ++column) {
            isChanged |= privateCode::eucDistPropagate(
              index0, index0 + columns - 1, distanceMap, referentMap,
              static_cast<int>(row), static_cast<int>(column));
            isChanged |= privateCode::eucDistPropagate(
              index0, index0 + columns, distanceMap, referentMap,
              static_cast<int>(row), static_cast<int>(column));
            isChanged |= privateCode::eucDistPropagate(
              index0, index0 + columns + 1, distanceMap, referentMap,
              static_cast<int>(row), static_cast<int>(column));
            ++index0;
          }

          // Update last pixel of the row.
          isChanged |= privateCode::eucDistPropagate(
            index0, index0 + columns - 1, distanceMap, referentMap,
            static_cast<int>(row), static_cast<int>(columns) - 1);
          isChanged |= privateCode::eucDistPropagate(
            index0, index0 + columns, distanceMap, referentMap,
            static_cast<int>(row), static_cast<int>(columns) - 1);
          ++index0;
        }

        // If no pixels needed changing, we're done.
        if(!isChanged) {
          break;
        }
      }

      // Convert squared distance to distance.
      for(index0 = 0; index0 < distanceMap.size(); ++index0) {
        distanceMap(index0) = std::sqrt(distanceMap(index0));
      }

      numberOfPassesUsed = passNumber;
      return distanceMap;
    }
  
  } // namespace computerVision
  
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_EUCLIDEANDISTANCE_IMPL_HH */






