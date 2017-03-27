/**
***************************************************************************
* @file brick/computerVision/connectedComponents_impl.hh
*
* Header file declaring inline and template functions declared in
* connectedComponents.hh.
*
* Copyright (C) 2006,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_CONNECTEDCOMPONENTS_IMPL_HH
#define BRICK_COMPUTERVISION_CONNECTEDCOMPONENTS_IMPL_HH

// This file is included by connectedComponents.hh, and should not be
// directly included by user code, so no need to include
// connectedComponents.hh here.
// 
// #include <brick/computerVision/connectedComponents.hh>

#include <cmath>
#include <memory>
#include <brick/computerVision/disjointSet.hh>

namespace brick {

  namespace computerVision {

    /// @cond privateCode
    namespace privateCode {

      template<ImageFormat FORMAT_OUT>
      void
      populateOutputImage(Image<FORMAT_OUT>& outputImage,
                          brick::numeric::Array2D<size_t> const& labelImage,
                          std::vector<size_t> const& labelArray);
      

    } // namespace privateCode
    /// @endcond

  
    // This function does connected components analysis on a previously
    // segmented image.
    template<ImageFormat FORMAT_OUT, ImageFormat FORMAT_IN>
    Image<FORMAT_OUT>
    connectedComponents(const Image<FORMAT_IN>& inputImage)
    {
      unsigned int numberOfComponents;
      return connectedComponents<FORMAT_OUT>(inputImage, numberOfComponents);
    }


    // This function does connected components analysis on a previously
    // segmented image.
    template<ImageFormat FORMAT_OUT, ImageFormat FORMAT_IN>
    Image<FORMAT_OUT>
    connectedComponents(const Image<FORMAT_IN>& inputImage,
                        unsigned int& numberOfComponents)
    {
      typedef typename Image<FORMAT_IN>::const_iterator InIterator;
      typedef brick::numeric::Array2D<size_t>::iterator LabelIterator;

      // Allocate storage for intermediate and final results.
      Image<FORMAT_OUT> outputImage(inputImage.rows(), inputImage.columns());
      brick::numeric::Array2D<size_t> labelImage(inputImage.rows(), 
                                                 inputImage.columns());

      // This vector will do the accounting of which components abut
      // one another.  We use unique_ptr here because DisjointSet is
      // not copyable.
      std::vector< std::unique_ptr< DisjointSet<size_t> > >
        correspondenceVector;

      // We'll label the very first component as 1. Labels of zero
      // mean background.
      size_t currentLabel = 0;
      correspondenceVector.emplace_back(new DisjointSet<size_t>(0));

      // This variable will be used to keep track of whether or not the
      // previous pixel was part of a blob.
      bool isActive = false;
    
      // Get iterators pointing to the first pixel of the input image, and
      // the first pixel of the label image.
      InIterator inIter = inputImage.begin();
      LabelIterator labelIter = labelImage.begin();

      // Label the first row.
      for(size_t columnIndex = 0; columnIndex < inputImage.columns();
          ++columnIndex) {
        if(!(*inIter)) {
          // The current pixel is background.  Mark it as such.
          *labelIter = 0;
          isActive = false;
        } else if(isActive) {
          // The current pixel is in a blob, and the previous pixel
          // was in a blob.  Adopt the same label as the previous
          // pixel.  Note there's no need to set isActive... it's
          // already true.
          *labelIter = currentLabel;
          // isActive = true;
        } else {
          // The current pixel is in a blob, but the previous pixel
          // was not.  This is a new blob!  Get a new label.
          ++currentLabel;
          correspondenceVector.emplace_back(
            new DisjointSet<size_t>(currentLabel));
          *labelIter = currentLabel;
          isActive = true;              
        }
        // Move to the next pixel.
        ++inIter;
        ++labelIter;
      }

      // Label the remaining rows.
      size_t const numberOfColumns = inputImage.columns();
      size_t previousLabel = 0;
      
      for(size_t rowIndex = 1; rowIndex < inputImage.rows(); ++rowIndex) {
        // At the beginning of each row, pretend there was a zero
        // pixel that we just came from, and update isActive
        // accordingly.  This lets us avoid adding special case code
        // for the first column.
        isActive = false;
        previousLabel = 0;

        // Iterate over the current row.
        for(size_t columnIndex = 0; columnIndex < inputImage.columns();
            ++columnIndex) {

          // Get the label of the pixel one row above the current
          // pixel.
          size_t parentLabel = *(labelIter - numberOfColumns);

          if(!(*inIter)) {
            // The current pixel is background.  This is one of our
            // two most likely cases, so handle it as simply as
            // possible.
            *labelIter = 0;
            isActive = false;
          } else if(isActive) {
            // The current pixel is in the interior of a blob, which
            // is the other most likely case.  Handle it as simply as
            // possible.
            *labelIter = previousLabel;

            // No need to set isActive... it's already true.
            // isActive = true;

            // We may have just joined two blobs.  Record the
            // correspondence, if appropriate.
            if(parentLabel && (parentLabel != previousLabel)) {
              correspondenceVector[previousLabel]->merge(
                  *(correspondenceVector[parentLabel]));
            }
          } else {
            // The current pixel is in a blob, but the previous pixel
            // was not.  This puts on the left edge of a blob.  Set
            // isActive so that we'll remember next iteration that we
            // just labeled a pixel.
            isActive = true;
            
            // Perhaps the pixel in the row above was also part of
            // this blob.
            if(parentLabel) {
              // The pixel in the previous row was inside a blob.
              // Adopt its label.
              *labelIter = parentLabel;
              previousLabel = parentLabel;
            } else {
              // The pixel in the previous row was not in a blob.  The
              // blob the current pixel is in might be new!  Get a new
              // label.
              ++currentLabel;
              correspondenceVector.emplace_back(
                new DisjointSet<size_t>(currentLabel));
              *labelIter = currentLabel;
              previousLabel = currentLabel;
            }
          }
          // Move to the next pixel.
          ++inIter;
          ++labelIter;
        }
      }

      // === Resolve label equivalences. ===

      // Create a look up table which will take tentative labels
      // (assigned above) and map them to finalized labels.  We arrange
      // that finalized labels are 1, 2, 3, etc., without any gaps,
      // but make no guarantees about which blob gets which label.
      std::vector<size_t> labelArray(correspondenceVector.size(), 0);
      size_t outputLabel = 0;
      for(size_t ii = 0; ii < labelArray.size(); ++ii) {
        if(labelArray[ii] == 0) {
          size_t headOfFamily = correspondenceVector[ii]->find().getPayload();
          if(labelArray[headOfFamily] != 0) {
            labelArray[ii] = labelArray[headOfFamily];
          } else {
            labelArray[ii] = outputLabel;
            labelArray[headOfFamily] = outputLabel;
            ++outputLabel;
          }
        }
      }

      // Relabel the blobs.
      privateCode::populateOutputImage(outputImage, labelImage, labelArray);

      if(currentLabel != 0) {
        numberOfComponents = static_cast<unsigned int>(currentLabel - 1);
      } else {
        // This case will never happen, but we leave the code here so
        // the compiler (and collaborators) won't worry.
        numberOfComponents = 0;
      }
    
      return outputImage;
    }    


    /// @cond privateCode
    namespace privateCode {

      template<ImageFormat FORMAT_OUT>
      void
      populateOutputImage(Image<FORMAT_OUT>& outputImage,
                          brick::numeric::Array2D<size_t> const& labelImage,
                          std::vector<size_t> const& labelArray)
      {
        auto labelIter = labelImage.begin();
        auto outIter = outputImage.begin();
        
        while(labelIter != labelImage.end()) {
          *outIter = static_cast<typename ImageFormatTraits<FORMAT_OUT>::PixelType>(labelArray[*labelIter]);
          ++outIter;
          ++labelIter;
        }
        
      }
      
    } // namespace privateCode
    /// @endcond
    
  } // namespace computerVision
    
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_CONNECTEDCOMPONENTS_IMPL_HH */
