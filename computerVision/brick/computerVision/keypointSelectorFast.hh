/**
***************************************************************************
* @file brick/computerVision/keypointSelectorFast.hh
*
* Header file declaring a class template for selecting stable
* keypoints from an image.
*
* Copyright (C) 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_KEYPOINTSELECTORFAST_HH
#define BRICK_COMPUTERVISION_KEYPOINTSELECTORFAST_HH

#include <limits>
#include <vector>
#include <brick/computerVision/image.hh>
#include <brick/numeric/index2D.hh>

namespace brick {

  namespace computerVision {

    struct KeypointFast {
      int row;
      int column;
      brick::common::UnsignedInt8 featureVector[16];
      bool isPositive;
    };
    

    /**
     ** This class template selects keypoints from an input image
     ** using Rosten's FAST keypoint detector [1].  Unlike other
     ** keypoint selectors in this library, it does not use a scale
     ** space, and operates directly on the input image.
     **/
    class KeypointSelectorFast {
    public:

      // ========= Public member functions. =========

      KeypointSelectorFast();

      
      void
      estimateThreshold(
        Image<GRAY8> const& inImage, 
        unsigned int startRow = 0,
        unsigned int startColumn = 0,
        unsigned int stopRow = std::numeric_limits<unsigned int>::max(),
        unsigned int stopColumn = std::numeric_limits<unsigned int>::max(),
        unsigned int expectedKeypointsPerImage = 500);

      
      std::vector<KeypointFast>
      getKeypoints();


      brick::common::Int16
      getThreshold();

      
      void
      setImage(
        Image<GRAY8> const& inImage,
        unsigned int startRow = 0,
        unsigned int startColumn = 0,
        unsigned int stopRow = std::numeric_limits<unsigned int>::max(),
        unsigned int stopColumn = std::numeric_limits<unsigned int>::max());

    private:

      void
      checkAndRepairRegionOfInterest(Image<GRAY8> const& inImage,
                                     unsigned int pixelMeasurementRadius,
                                     unsigned int& startRow,
                                     unsigned int& startColumn,
                                     unsigned int& stopRow,
                                     unsigned int& stopColumn);

      brick::common::Int16
      measurePixelThreshold(Image<GRAY8> const& image,
                            unsigned int row, unsigned int column);
        
      void
      measurePixelThreshold(Image<GRAY8> const& inImage,
                            unsigned int row, unsigned int column,
                            brick::common::Int16& threshold);

      bool
      testPixel(Image<GRAY8> const& image,
                unsigned int row, unsigned int column,
                const common::Int16 threshold,
                KeypointFast& keypoint);


      bool
      testPixelDetails(Image<GRAY8> const& image,
                       unsigned int row, unsigned int column,
                       const brick::common::Int16 testValue,
                       const brick::common::Int16 threshold,
                       KeypointFast& keypoint,
                       bool isPositive);
      
      std::vector<KeypointFast> m_keypointVector;
      brick::common::Int16 m_threshold;
    };

  } // namespace computerVision
  
} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/keypointSelectorFast_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_KEYPOINTSELECTORFAST_HH */
