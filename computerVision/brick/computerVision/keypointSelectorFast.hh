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

      static const unsigned int numberOfFeatures = 16;
    };
    

    /**
     ** This class template selects keypoints from an input image
     ** using Rosten's FAST keypoint detector [1,2].  Unlike other
     ** keypoint selectors in this library, it does not use a scale
     ** space, and operates directly on the input image.
     **
     ** [1] E. Rosten and T. Drummond, "Fusing points and lines for
     ** high performance tracking.", IEEE International Conference on
     ** Computer Vision, Vol 2, pp 1508--1511", 2005.
     **
     ** [2] E. Rosten and T. Drummond, "Machine learning for
     ** high-speed corner detection", European Conference on Computer
     ** Vision, Vol 1, pp 430-443, 2006.
     **/
    class KeypointSelectorFast {
    public:

      // ========= Public member functions. =========


      /** 
       * Default constructor.
       */
      KeypointSelectorFast();


      /** 
       * This function tries to automatically set the internal
       * threshold of the algorithm based on an input image.  Lower
       * thresholds lead to more keypoints being detected.  If you
       * like, you can also set the threshold explicitly using member
       * function setThreshold().  Note that you may need to adjust
       * the threshold over time as lighting in the image changes.
       * 
       * @param inImage This argument is the image from which to
       * estimate the threshold.
       * 
       * @param startRow This argument specifies the bounding box of
       * the image region to use during estimation.  Omit it to
       * process the whole image.
       * 
       * @param startColumn This argument specifies the bounding box of
       * the image region to use during estimation.  Omit it to
       * process the whole image.
       * 
       * @param stopRow This argument specifies the bounding box of
       * the image region to use during estimation.  Omit it to
       * process the whole image.
       * 
       * @param stopColumn This argument specifies the bounding box of
       * the image region to use during estimation.  Omit it to
       * process the whole image.
       * 
       * @param expectedKeypointsPerImage This argument specifies
       * (roughly) how many keypoints you'd like to the automatically
       * tuned threshold to detect.
       */
      void
      estimateThreshold(
        Image<GRAY8> const& inImage, 
        unsigned int startRow = 0,
        unsigned int startColumn = 0,
        unsigned int stopRow = std::numeric_limits<unsigned int>::max(),
        unsigned int stopColumn = std::numeric_limits<unsigned int>::max(),
        unsigned int expectedKeypointsPerImage = 500);


      /** 
       * Return the keypoints detected during the most recent call to
       * member function setImage().
       * 
       * @return The return value is vector of KeypointFast instances.
       */
      std::vector<KeypointFast>
      getKeypoints() const;


      /** 
       * Return the value of the threshold used in keypoint detection.
       * See member function estimateThreshold().
       * 
       * @return The threshold value.
       */
      brick::common::Int16
      getThreshold() const;


      /** 
       * Process an image to find keypoints.
       * 
       * @param inImage This argument is the image in which to look
       * for keypoints.
       * 
       * @param startRow This argument specifies the bounding box of
       * the image region to use during estimation.  Omit it to
       * process the whole image.
       * 
       * @param startColumn This argument specifies the bounding box of
       * the image region to use during estimation.  Omit it to
       * process the whole image.
       * 
       * @param startRow This argument specifies the bounding box of
       * the image region to use during estimation.  Omit it to
       * process the whole image.
       * 
       * @param startColumn This argument specifies the bounding box of
       * the image region to use during estimation.  Omit it to
       * process the whole image.
       */
      void
      setImage(
        Image<GRAY8> const& inImage,
        unsigned int startRow = 0,
        unsigned int startColumn = 0,
        unsigned int stopRow = std::numeric_limits<unsigned int>::max(),
        unsigned int stopColumn = std::numeric_limits<unsigned int>::max());


      /** 
       * Manually set the internal threshold of the algorithm.  See
       * also member function estimateThreshold().
       * 
       * @param threshold This argument specifies the desired threshold.
       */
      void
      setThreshold(brick::common::Int16 threshold);

      
    private:

      // Make sure bounding box of processing region is sane.
      void
      checkAndRepairRegionOfInterest(Image<GRAY8> const& inImage,
                                     unsigned int pixelMeasurementRadius,
                                     unsigned int& startRow,
                                     unsigned int& startColumn,
                                     unsigned int& stopRow,
                                     unsigned int& stopColumn) const;

      // Find the highest threshold value that would still allow this
      // particular pixel to pass and be selected as a keypoint.
      brick::common::Int16
      measurePixelThreshold(Image<GRAY8> const& image,
                            unsigned int row, unsigned int column) const;
      

      // Check to see if a specific pixel should be selected as a keypoint.
      bool
      testPixel(Image<GRAY8> const& image,
                unsigned int row, unsigned int column,
                const common::Int16 threshold,
                KeypointFast& keypoint) const;


      // This is called by testPixel to do the heavy lifting.
      bool
      testPixelDetails(Image<GRAY8> const& image,
                       unsigned int row, unsigned int column,
                       const brick::common::Int16 testValue,
                       const brick::common::Int16 threshold,
                       KeypointFast& keypoint,
                       bool isPositive) const;

      /* ======== Data members ========= */
      std::vector<KeypointFast> m_keypointVector;
      brick::common::Int16 m_threshold;
    };

  } // namespace computerVision
  
} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/keypointSelectorFast_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_KEYPOINTSELECTORFAST_HH */
