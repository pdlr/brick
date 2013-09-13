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
      brick::common::UInt8 darkColor;
      brick::common::UInt8 lightColor;
      brick::geometry::Bullseye2D<FloatType> bullseye;

      KeypointBullseye()
        : row(0),
          column(0),
          asymmetry(0.0),
          bullseyeMetric(0.0),
          darkColor(0),
          lightColor(0),
          bullseye() {}

      KeypointBullseye(CoordinateType rowArg,
                       CoordinateType columnArg)
        : row(rowArg),
          column(columnArg),
          asymmetry(0.0),
          bullseyeMetric(0.0),
          darkColor(0),
          lightColor(0),
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
     **    = myKeypointSelector.getKeypointsGeneralPosition();
     ** @endcode
     **/
    template <class FloatType>
    class KeypointSelectorBullseye {
    public:

      // ========= Public member functions. =========


      /** 
       * Default constructor.
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
       *
       * @param numberOfTransitions This argument tells the selector
       * how many light-to-dark and dark-to-light transitions to
       * expect when moving from the center of the bulllseye toward
       * the edge.  For example, a bullseye with a colored center and
       * one colored ring would have 3 transitions (center to lighter
       * background; background to first ring; and first ring to
       * background).
       *
       * @param isGeneralPositionRequired Use this argument to turn
       * off computation of subpixel bullseye positions if you don't
       * want them.  This saves a little bit of computation.
       */
      KeypointSelectorBullseye(brick::common::UInt32 maxNumberOfBullseyes,
                               brick::common::UInt32 maxRadius,
                               brick::common::UInt32 minRadius,
                               brick::common::UInt32 numberOfTransitions = 3,
                               bool isGeneralPositionRequired = true);


      /** 
       * Given a keypoint, find its subpixel (general position)
       * location.  This function is exposed publically to give the
       * user flexibility in finding locations in modified version of
       * the original image.
       * 
       * @param inputKeypoint This argument is the keypoint to be
       * fine-tuned.
       * 
       * @param inImage This argument is the input image against which
       * to compute the refined position.
       * 
       * @return The return value is a copy of the input keypoint with
       * its row and column updated to subpixel values.
       */
      KeypointBullseye<FloatType, FloatType>
      fineTuneKeypoint(
        KeypointBullseye<brick::common::Int32, FloatType> const& inputKeypoint,
        Image<GRAY8> const& inImage);

      
      /** 
       * Return the keypoints detected during the most recent call to
       * member function setImage().
       * 
       * @return The return value is vector of KeypointBullseye instances.
       */
      std::vector< KeypointBullseye<brick::common::Int32, FloatType> >
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
      std::vector< KeypointBullseye<FloatType, FloatType> >
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
               brick::common::UInt32 startRow,
               brick::common::UInt32 startColumn,
               brick::common::UInt32 stopRow,
               brick::common::UInt32 stopColumn);


    private:

      typedef brick::common::Int32 AccumulatedType;


      // Accumulate statistics related to the difference between two
      // pixel values.
      void
      accumulateAsymmetrySums(
        brick::common::Int32 pixel0,
        brick::common::Int32 pixel1,
        brick::common::UInt32& pixelSum,
        brick::common::UInt32& pixelSquaredSum,
        brick::common::UInt32& asymmetrySum) const;


      bool
      countTransitions(std::vector<brick::common::UInt8> const& spoke,
                       brick::common::UInt32 numberOfTransitions,
                       brick::common::UInt8 minDynamicRange,
                       brick::common::UInt8& darkColor,
                       brick::common::UInt8& lightColor,
                       brick::common::UInt32 minRadius,
                       brick::common::UInt32& actualRadius) const;
      
      
      // Make sure bounding box of processing region is sane.
      void
      checkAndRepairRegionOfInterest(brick::common::UInt32 rows,
                                     brick::common::UInt32 columns,
                                     brick::common::UInt32 radius,
                                     brick::common::UInt32& startRow,
                                     brick::common::UInt32& startColumn,
                                     brick::common::UInt32& stopRow,
                                     brick::common::UInt32& stopColumn) const;


      bool
      estimateBullseye(
        brick::geometry::Bullseye2D<FloatType>& bullseye,
        std::vector< std::vector< brick::numeric::Vector2D<FloatType> > > const&
          edgePositions,
        brick::common::UInt32 numberOfTransitions) const;

      
      // Estimate how much the target is squished along each axis.  Is
      // it circular?  Elliptical?
      void
      estimateScale(
        Image<GRAY8> const& image,
        brick::common::UInt32 radius,
        brick::common::UInt32 row, brick::common::UInt32 column,
        KeypointBullseye<brick::common::Int32, FloatType>& keypoint) const;


      
      // Figure out what "normal" is for the asymmetry measure, and
      // pick a threshold that's low enough.  Low enough means that
      // only interesting pixels have a lower score from
      // evaluateAsymmetry() (where low means "symmetrical").
      FloatType
      estimateAsymmetryThreshold(Image<GRAY8> const& inImage,
                                 brick::common::UInt32 minRadius,
                                 brick::common::UInt32 maxRadius,
                                 brick::common::UInt32 startRow,
                                 brick::common::UInt32 startColumn,
                                 brick::common::UInt32 stopRow,
                                 brick::common::UInt32 stopColumn,
                                 brick::common::UInt32 numberOfSamples) const;


      // See if an image location is plausibly the center of a
      // bullseye by looking for symetry around it.  Note that the
      // symmetry computation currently just compares up with down,
      // left with right, opposing diagonals.
      bool
      evaluateAsymmetry(Image<GRAY8> const& image,
                        brick::common::UInt32 radius,
                        brick::common::UInt32 row, brick::common::UInt32 column,
                        FloatType& asymmetry) const;

      
      inline bool
      testAndRecordEdges(
        Image<GRAY1> const& edgeImage,
        int const& row,
        int const& column,
        std::vector< std::vector< brick::numeric::Vector2D<FloatType> > >&
        edgePositions,
        brick::common::UInt32& edgeCount,
        brick::common::UInt32 const& numberOfTransitions) const
      {
        if(edgeImage(row, column)) {
          edgePositions[edgeCount].push_back(
            brick::numeric::Vector2D<FloatType>(column, row));
          return (++edgeCount >= numberOfTransitions);
        }
        return false;
      }
        

      inline bool
      testAndRecordEdgesDiagonal(
        Image<GRAY1> const& edgeImage,
        int const& row,
        int const& column,
        int const& d0,
        int const& d1,
        std::vector< std::vector< brick::numeric::Vector2D<FloatType> > >&
        edgePositions,
        brick::common::UInt32& edgeCount,
        brick::common::UInt32 const& numberOfTransitions) const
      {
        if(edgeImage(row, column)
           || (edgeImage(row + d0, column)
               && edgeImage(row, column + d1)
               && (!edgeImage(row + d0, column + d1)))) {
          edgePositions[edgeCount].push_back(
            brick::numeric::Vector2D<FloatType>(column, row));
          return (++edgeCount >= numberOfTransitions);
        }
        return false;
      }

      
      // Compute a measure of bullseye-ness that's more expensive --
      // and more accurate -- than asymmetry, and Fill out a KeyPoint
      // instance with the corresponding information.
      void
      evaluateBullseyeMetric(
        KeypointBullseye<brick::common::Int32, FloatType>& keypoint,
        Image<GRAY1> const& edgeImage,
        brick::numeric::Array2D<FloatType> const& gradientX,
        brick::numeric::Array2D<FloatType> const& gradientY,
        brick::common::UInt32 minRadius,
        brick::common::UInt32 maxRadius) const;


      bool
      isPlausibleBullseye(
        KeypointBullseye<brick::common::Int32, FloatType>& keypoint,
        Image<GRAY8> const& inImage,
        brick::common::UInt32 minRadius,
        brick::common::UInt32 maxRadius,
        FloatType asymmetryThreshold,
        bool forceAsymmetry = false) const;

      
      // Insert the new keypoint into a sorted vector, discarding the
      // worst point if the addition would make the vector longer than
      // maxNumberOfBullseyes.
      void
      sortedInsert(
        KeypointBullseye<brick::common::Int32, FloatType> const& keypoint,
        std::vector< KeypointBullseye<brick::common::Int32, FloatType> >& keypointVector,
        brick::common::UInt32 maxNumberOfBullseyes);


      bool
      validateBullseye(brick::geometry::Bullseye2D<FloatType> const& bullseye,
                       // Image<GRAY8> const& inImage,
                       Image<GRAY1> const& edgeImage,
                       brick::numeric::Array2D<FloatType> const& gradientX,
                       brick::numeric::Array2D<FloatType> const& gradientY,
                       brick::common::UInt32 row,
                       brick::common::UInt32 column,
                       brick::common::UInt32 minRadius,
                       brick::common::UInt32 maxRadius,
                       FloatType& goodness) const;
      

      // Private data members.

      std::vector< brick::numeric::Vector2D<FloatType> >
        m_bullseyePoints;
      std::vector<brick::common::UInt32>
        m_bullseyeEdgeCounts;
      std::vector< std::vector< brick::numeric::Vector2D<FloatType> > >
        m_edgePositions;

      bool m_isGeneralPositionRequired;      
      std::vector< KeypointBullseye<brick::common::Int32, FloatType> > m_keypointVector;
      std::vector< KeypointBullseye<FloatType, FloatType> > m_keypointGPVector;
      brick::common::UInt32 m_maxNumberOfBullseyes;
      brick::common::UInt32 m_numberOfTransitions;
      brick::common::UInt32 m_maxRadius;
      brick::common::UInt32 m_minRadius;
      brick::common::UInt8 m_minDynamicRange;

    };

  } // namespace computerVision
  
} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/keypointSelectorBullseye_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_KEYPOINTSELECTORBULLSEYE_HH */
