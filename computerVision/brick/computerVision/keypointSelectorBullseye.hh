/**
***************************************************************************
* @file brick/computerVision/keypointSelectorBullseye.hh
*
* Header file declaring a class template for selecting stable
* keypoints from an image.
*
* Copyright (C) 2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_KEYPOINTSELECTORBULLSEYE_HH
#define BRICK_COMPUTERVISION_KEYPOINTSELECTORBULLSEYE_HH

#include <limits>
#include <vector>
#include <brick/computerVision/image.hh>
#include <brick/numeric/index2D.hh>

namespace brick {

  namespace computerVision {

    template <class CoordinateType>
    struct KeypointBullseye {
      CoordinateType row;
      CoordinateType column;
      CoordinateType value;

    private:

    public:
      KeypointBullseye(CoordinateType rowArg,
                     CoordinateType columnArg,
                     CoordinateType valueArg)
        : row(rowArg),
          column(columnArg),
          value(valueArg) {}

    };
    

    /**
     ** This class template looks for bullseye targets in an input
     ** image.  It does not use a scale space, and operates directly
     ** on the input image.
     **/
    template <class FloatType>
    class KeypointSelectorBullseye {
    public:

      // ========= Public member functions. =========


      /** 
       * Default constructor.
       * 
       * @param kappa This argument is the free parameter in the
       * Bullseye corner metric computation.
       * 
       * @param maxNumberOfBullseyes Use this argument to indicate how
       * many bullseye targets you expect to find in the image.
       */
      KeypointSelectorBullseye(unsigned int maxNumberOfBullseyes);


      /** 
       * Return the keypoints detected during the most recent call to
       * member function setImage().
       * 
       * @return The return value is vector of KeypointBullseye instances.
       */
      std::vector< KeypointBullseye<brick::common::Int32> >
      getKeypoints() const;


      /** 
       * Return the keypoints detected during the most recent call to
       * member function setImage().
       * 
       * @param iterator This argument must be a writable iterator
       * pointing to KeypoindBullseye<brick::common::Int32>.
       * 
       * @param threshold Increasing this threshold eliminates weaker
       * corners from the output.
       */
      template <class Iter>
      void
      getKeypoints(Iter iterator, FloatType threshold = 0.0) const;


      /** 
       * Return the keypoints detected during the most recent call to
       * member function setImage(), and use subpixel interpolation to
       * refine their positions.
       * 
       * @return The return value is vector of KeypointBullseye instances.
       */
      std::vector< KeypointBullseye<FloatType> >
      getKeypointsGeneralPosition() const;

      
      /** 
       * Return the keypoints detected during the most recent call to
       * member function setImage(), and use subpixel interpolation to
       * refine their positions.
       * 
       * @param iterator This argument must be a writable iterator
       * pointing to KeypoindBullseye<FloatType>.
       * 
       * @param threshold Increasing this threshold eliminates weaker
       * corners from the output.
       */
      template <class Iter>
      void
      getKeypointsGeneralPosition(Iter iterator) const;
      
      
      /** 
       * Process an image to find keypoints.
       * 
       * @param inImage This argument is the image in which to look
       * for keypoints.
       */
      void
      setImage(Image<GRAY8> const& inImage);


    private:

      typedef brick::common::Int32 AccumulatedType;


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

    };

  } // namespace computerVision
  
} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/keypointSelectorBullseye_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_KEYPOINTSELECTORBULLSEYE_HH */
