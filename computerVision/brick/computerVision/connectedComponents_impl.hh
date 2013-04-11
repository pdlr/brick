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

namespace brick {

  namespace computerVision {

    /// @cond privateCode
    namespace privateCode {

      // This function traverses the blob label graph, starting at a
      // particular node, propagating that node's label so that any blobs
      // which are connected to the starting node inherit the the starting
      // nodes label iff the starting nodes label is less than the current
      // label of the connected blob.
      bool
      propagateLabel(size_t label, size_t node, std::vector<size_t>& labelArray,
                     const std::vector< std::list<size_t> >& neighborsVector);

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
      typedef typename Image<FORMAT_OUT>::iterator OutIterator;
      typedef brick::numeric::Array2D<size_t>::iterator LabelIterator;
    
      Image<FORMAT_OUT> outputImage(inputImage.rows(), inputImage.columns());
      brick::numeric::Array2D<size_t> labelImage(inputImage.rows(), 
                                                 inputImage.columns());
      size_t currentLabel = 0;

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
        if(*inIter) {
          // The pixel value is nonzero, we're in a blob.
          if(!isActive) {
            // If the previous pixel was zero, we're in what might be a
            // new blob, so increment the blob label.
            ++currentLabel;
            isActive = true;
          }
          // Label the pixel with its tentative label.
          *labelIter = currentLabel;
        } else {
          // The pixel value is zero, not in a blob.
          *labelIter = 0;
          isActive = false;
        }
        // Move to the next pixel.
        ++inIter;
        ++labelIter;
      }

      // Label the remaining rows.
      std::vector< std::pair<size_t, size_t> > correspondences;
      LabelIterator chaseIter = labelImage.begin();
      size_t workingLabel = 0;
      size_t previousParentLabel = 0;
      for(size_t rowIndex = 1; rowIndex < inputImage.rows(); ++rowIndex) {
        // At the beginning of each row, pretend there was a zero pixel
        // that we just came from, and update isActive accordingly.
        isActive = false;
        for(size_t columnIndex = 0; columnIndex < inputImage.columns();
            ++columnIndex) {
          if(*inIter) {
            // The pixel value is nonzero, we're in a blob.  Get the
            // blob label of the pixel in the row above (0 means "not in
            // a blob").
            size_t parentLabel = *chaseIter;
            if(!isActive) {
              // The pixel in the previous column was zero.  This blob
              // might be a new one.
              if((parentLabel != 0)) {
                // The pixel in the previous row was part of a blob.
                // Adopt its label.
                workingLabel = parentLabel;
                previousParentLabel = parentLabel;
              } else {
                // The pixel in the previous row was not part of a blob.
                // Get a new label.
                ++currentLabel;
                workingLabel = currentLabel;
                previousParentLabel = 0;
              }
              isActive = true;
            } else {
              // The pixel in the previous column was part of a blob.
              if((parentLabel != 0)
                 && (parentLabel != previousParentLabel)) {
                // We just connected a blob in the previous row with the
                // blob in the previous column.  Record the
                // correspondence.
                correspondences.push_back(
                  std::make_pair(parentLabel, workingLabel));
                previousParentLabel = parentLabel;
              }
            }
            // Assign the appropriate label to this pixel.
            *labelIter = workingLabel;
          } else {
            // We're at a zero pixel.  Update its label accordingly.
            *labelIter = 0;
            isActive = false;
          }
          // Move to the next pixel.
          ++inIter;
          ++labelIter;
          ++chaseIter;
        }
      }

      // === Resolve label equivalences. ===

      // Create an look up table which will take tentative labels
      // (assigned above) and map them to finalized labels.  We'll start
      // by assigning each element in the lookup table an impossible
      // value, and then we'll go back and correct each element in turn.
      std::vector<size_t> labelArray(currentLabel + 1);
      for(size_t label = 0; label < labelArray.size(); ++label) {
        labelArray[label] = currentLabel + 1;
      }

      // Imagine a graph in which each node is a label, and each edge is
      // an equivalence between labels.  Create a list of edges for each
      // node in the graph.
      typedef std::vector< std::pair<size_t, size_t> >::const_iterator
        CorrespondenceIterator;
      std::vector< std::list<size_t> > neighborsVector(labelArray.size());
      CorrespondenceIterator cIter = correspondences.begin();
      while(cIter != correspondences.end()) {
        neighborsVector[cIter->first].push_back(cIter->second);
        neighborsVector[cIter->second].push_back(cIter->first);
        ++cIter;
      }

      // Propagate labels over all edges.
      currentLabel = 0;
      for(size_t node = 0; node < labelArray.size(); ++node) {
        if(privateCode::propagateLabel(
             currentLabel, node, labelArray, neighborsVector)) {
          ++currentLabel;
        }
      }

      // Relabel the blobs.
      labelIter = labelImage.begin();
      OutIterator outIter = outputImage.begin();
      while(labelIter != labelImage.end()) {
		*outIter = static_cast<typename ImageFormatTraits<FORMAT_OUT>::PixelType>(labelArray[*labelIter]);
        ++outIter;
        ++labelIter;
      }

      if(currentLabel != 0) {
        numberOfComponents = static_cast<unsigned int>(currentLabel - 1);
      } else {
        numberOfComponents = 0;
      }
    
      return outputImage;
    }

  } // namespace computerVision
    
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_CONNECTEDCOMPONENTS_IMPL_HH */
