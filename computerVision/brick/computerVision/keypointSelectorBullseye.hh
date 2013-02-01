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
#include <brick/geometry/bullseye2D.hh>
#include <brick/numeric/index2D.hh>
#include <brick/numeric/vector2D.hh>

namespace brick {

  namespace computerVision {

    unsigned int const keypointBullseyeMaxRadius = 64;
    
    template <class CoordinateType, class FloatType = double>
    struct KeypointBullseye {
      CoordinateType row;
      CoordinateType column;
      FloatType asymmetry;
      FloatType bullseyeMetric;
      brick::geometry::Bullseye2D<FloatType> bullseye;

      KeypointBullseye(CoordinateType rowArg,
                       CoordinateType columnArg,
                       FloatType asymmetryArg)
        : row(rowArg),
          column(columnArg),
          asymmetry(asymmetryArg),
          bullseyeMetric(0.0),
          bullseye() {}
    };
    

    /**
     ** This class template looks for bullseye targets in an input
     ** image.  It does not use a scale space, and operates directly
     ** on the input image.  Use it like this:
     **
     ** @code
     ** KeypointSelectorBullseye<double> myKeypointSelector(2, 40, 10);
     ** myKeypointSelector.setImage(myGrayscaleImage);
     ** std::vector< KeypointBullseye<double> > keypoints
     **    = myKeypointSelector.getKeypointsGeneralPosition() const;
     ** @endcode
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
       * @param maxRadius Use this argument to indicate the size of
       * the biggest bullseye you expect to see.  For example, if the
       * biggest bullseye you expect to see is 20 pixels across, set
       * this to 10 (because radius is half of diameter).
       * 
       * @param minRadius Use this argument to indicate the size of
       * the smallest bullseye you expect to see.  For example, if the
       * smallest bullseye you expect to see is 8 pixels across, set
       * this to 4 (because radius is half of diameter).
       */
      KeypointSelectorBullseye(unsigned int maxNumberOfBullseyes,
                               unsigned int maxRadius,
                               unsigned int minRadius);


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
      

      /** 
       * Process a region of an image to find keypoints.
       * 
       * @param inImage This argument is the image in which to look
       * for keypoints.
       * 
       * @param startRow This argument helps to specify the upper-left
       * corner of the region to be searched.
       * 
       * @param startColumn This argument helps to specify the upper-left
       * corner of the region to be searched.
       * 
       * @param stopRow This argument helps to specify the lower-right
       * corner of the region to be searched.
       * 
       * @param stopColumn This argument  helps to specify the lower-right
       * corner of the region to be searched.
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
      accumulateAsymmetrySums(
        brick::common::Int32 pixel0,
        brick::common::Int32 pixel1,
        brick::common::UnsignedInt32& pixelSum,
        brick::common::UnsignedInt32& pixelSquaredSum,
        brick::common::UnsignedInt32& asymmetrySum) const;

      
      // Make sure bounding box of processing region is sane.
      void
      checkAndRepairRegionOfInterest(Image<GRAY8> const& inImage,
                                     unsigned int radius,
                                     unsigned int& startRow,
                                     unsigned int& startColumn,
                                     unsigned int& stopRow,
                                     unsigned int& stopColumn) const;


      bool
      estimateBullseye(
        brick::geometry::Bullseye2D<FloatType>& bullseye,
        std::vector< std::vector< brick::numeric::Vector2D<FloatType> > > const&
          edgePositions,
        unsigned int numberOfTransitions);

      
      // Estimate how much the target is squished along each axis.  Is
      // it circular?  Elliptical?
      void
      estimateScale(Image<GRAY8> const& image,
                    unsigned int radius,
                    unsigned int row, unsigned int column,
                    KeypointBullseye<brick::common::Int32>& keypoint) const;


      
      // Figure out what "normal" is for the asymmetry measure, and
      // pick a threshold that's low enough.  Low enough means that
      // only interesting pixels have a lower score from
      // evaluateAsymmetry() (where low means "symmetrical").
      FloatType estimateAsymmetryThreshold(
        Image<GRAY8> const& inImage,
        unsigned int radius,
        unsigned int startRow,
        unsigned int startColumn,
        unsigned int stopRow,
        unsigned int stopColumn,
        unsigned int numberOfSamples) const;


      // See if an image location is plausibly the center of a
      // bullseye by looking for symetry around it.  Note that the
      // symmetry computation currently just compares up with down,
      // left with right, opposing diagonals.
      bool
      evaluateAsymmetry(Image<GRAY8> const& image,
                        unsigned int radius,
                        unsigned int row, unsigned int column,
                        FloatType& asymmetry) const;

      
      // Compute a measure of bullseye-ness that's more expensive --
      // and more accurate -- than asymmetry, and Fill out a KeyPoint
      // instance with the corresponding information.
      void
      evaluateBullseyeMetric(
        KeypointBullseye<brick::common::Int32>& keypoint,
        Image<GRAY1> const& edgeImage,
        brick::numeric::Array2D<FloatType> const& gradientX,
        brick::numeric::Array2D<FloatType> const& gradientY,
        // unsigned int minRadius,
        unsigned int maxRadius);

      
      // Find the highest threshold value that would still allow this
      // particular pixel to pass and be selected as a keypoint.
      brick::common::Int16
      measurePixelThreshold(Image<GRAY8> const& image,
                            unsigned int row, unsigned int column) const;


      // Insert the new keypoint into a sorted vector, discarding the
      // worst point if the addition would make the vector longer than
      // maxNumberOfBullseyes.
      void
      sortedInsert(
        KeypointBullseye<brick::common::Int32> const& keypoint,
        std::vector< KeypointBullseye<brick::common::Int32> >& keypointVector,
        unsigned int maxNumberOfBullseyes);


      bool
      validateBullseye(brick::geometry::Bullseye2D<FloatType> const& bullseye,
                       // Image<GRAY8> const& inImage,
                       Image<GRAY1> const& edgeImage,
                       brick::numeric::Array2D<FloatType> const& gradientX,
                       brick::numeric::Array2D<FloatType> const& gradientY,
                       unsigned int row,
                       unsigned int column,
                       unsigned int maxRadius,
                       FloatType& goodness);
      

      // Private data members.

      std::vector< brick::numeric::Vector2D<FloatType> >
        m_bullseyePoints;
      std::vector<unsigned int>
        m_bullseyeEdgeCounts;
      std::vector< std::vector< brick::numeric::Vector2D<FloatType> > >
        m_edgePositions;
      
      std::vector< KeypointBullseye<brick::common::Int32> > m_keypointVector;
      brick::common::UnsignedInt32 m_maxNumberOfBullseyes;
      brick::common::UnsignedInt32 m_numberOfTransitions;
      brick::common::UnsignedInt32 m_maxRadius;
      brick::common::UnsignedInt32 m_minRadius;

    };

  } // namespace computerVision
  
} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/keypointSelectorBullseye_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_KEYPOINTSELECTORBULLSEYE_HH */
