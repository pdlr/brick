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

    unsigned int const keypointBullseyeMaxRadius = 64;
    
    template <class CoordinateType>
    struct KeypointBullseye {
      CoordinateType row;
      CoordinateType column;
      CoordinateType value;
      unsigned int horizontalScale;
      unsigned int verticalScale;

      brick::common::UnsignedInt8 leftSpoke[keypointBullseyeMaxRadius];
      brick::common::UnsignedInt8 rightSpoke[keypointBullseyeMaxRadius];
      brick::common::UnsignedInt8 topSpoke[keypointBullseyeMaxRadius];
      brick::common::UnsignedInt8 bottomSpoke[keypointBullseyeMaxRadius];
      
    private:

    public:
      KeypointBullseye(CoordinateType rowArg,
                     CoordinateType columnArg,
                     CoordinateType valueArg)
        : row(rowArg),
          column(columnArg),
          value(valueArg),
          horizontalScale(0),
          verticalScale(0),
          leftSpoke(),
          rightSpoke(),
          topSpoke(),
          bottomSpoke() {}

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
       *
       * @param maxRadius Use this argument to indicate the biggest
       * bullseye you expect to see.  For example, if the biggest
       * bullseye you expect to see is 20 pixels across, set this to
       * 10 (because radius is half of diameter).
       */
      KeypointSelectorBullseye(unsigned int maxNumberOfBullseyes,
                               unsigned int maxRadius);


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
      setImage(Image<GRAY8> const& inImage,
               unsigned int startRow,
               unsigned int startColumn,
               unsigned int stopRow,
               unsigned int stopColumn);


    private:

      typedef brick::common::Int32 AccumulatedType;


      // Accumulate statistics related to the difference between two
      // pixel values.
      void
      accumulateAsymmetrySums(brick::common::Int32 pixel0,
                              brick::common::Int32 pixel1,
                              brick::common::UnsignedInt32& pixelSum,
                              brick::common::UnsignedInt32& pixelSquaredSum,
                              brick::common::UnsignedInt32& asymmetrySum);

      
      // Make sure bounding box of processing region is sane.
      void
      checkAndRepairRegionOfInterest(Image<GRAY8> const& inImage,
                                     unsigned int radius,
                                     unsigned int& startRow,
                                     unsigned int& startColumn,
                                     unsigned int& stopRow,
                                     unsigned int& stopColumn) const;

      // Estimate how much the target is squished along each axis.  Is
      // it circular?  Elliptical?
      void
      estimateScale(Image<GRAY8> const& image,
                    unsigned int radius,
                    unsigned int row, unsigned int column,
                    KeypointBullseye<brick::common::Int32>& keypoint) const;
      

      // See if an image location is plausibly the center of a
      // bullseye by looking for symetry around it.  Note that the
      // symmetry computation currently just compares up with down and
      // left with right.
      FloatType
      evaluateSymmetry(Image<GRAY8> const& image,
                       unsigned int radius,
                       unsigned int row, unsigned int column,
                       KeypointBullseye<brick::common::UnsignedInt32>& keypoint)
        const;

      
      // Find the highest threshold value that would still allow this
      // particular pixel to pass and be selected as a keypoint.
      brick::common::Int16
      measurePixelThreshold(Image<GRAY8> const& image,
                            unsigned int row, unsigned int column) const;


      // Private data members.

      std::vector< KeypointBullseye<FloatType> > m_keypointVector;
      brick::common::UnsignedInt32 m_maxNumberOfBullseyes;
      brick::common::UnsignedInt32 m_maxRadius;
      

    };

  } // namespace computerVision
  
} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/keypointSelectorBullseye_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_KEYPOINTSELECTORBULLSEYE_HH */
