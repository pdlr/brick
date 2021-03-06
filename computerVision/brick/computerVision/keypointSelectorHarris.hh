/**
***************************************************************************
* @file brick/computerVision/keypointSelectorHarris.hh
*
* Header file declaring a class template for selecting stable
* keypoints from an image.
*
* Copyright (C) 2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_KEYPOINTSELECTORHARRIS_HH
#define BRICK_COMPUTERVISION_KEYPOINTSELECTORHARRIS_HH

#include <limits>
#include <vector>
#include <brick/computerVision/image.hh>
#include <brick/numeric/index2D.hh>

namespace brick {

  namespace computerVision {

    template <class CoordinateType>
    struct KeypointHarris {
      CoordinateType row;
      CoordinateType column;
      CoordinateType value;

    private:
      CoordinateType m_xx;
      CoordinateType m_xy;
      CoordinateType m_yy;

    public:
      KeypointHarris(CoordinateType rowArg,
                     CoordinateType columnArg,
                     CoordinateType valueArg,
                     CoordinateType xxArg,
                     CoordinateType xyArg,
                     CoordinateType yyArg)
        : row(rowArg),
          column(columnArg),
          value(valueArg),
          m_xx(xxArg),
          m_xy(xyArg),
          m_yy(yyArg) {}


      template <class FloatType>
      void
      getCovariance(FloatType& c00, FloatType& c01, FloatType& c11);
    };


    /**
     ** This class template selects keypoints from an input image
     ** using Harris's and Stephen's keypoint detector [1].  Unlike other
     ** keypoint selectors in this library, it does not use a scale
     ** space, and operates directly on the input image.
     **
     ** [1] C. Harris and M. Stephens, "A combined corner and edge
     ** detector", Proceedings of the 4th Alvey Vision Conference,
     ** pp. 147–151, 1988.
     **/
    template <class FloatType>
    class KeypointSelectorHarris {
    public:

      // ========= Public member functions. =========


      /**
       * Default constructor.
       *
       * @param kappa This argument is the free parameter in the
       * Harris corner metric computation.
       *
       * @param sigma This argument is currently ignored.
       */
      KeypointSelectorHarris(FloatType kappa = 0.1,
                             FloatType sigma = 1.6);


      /**
       * Return the keypoints detected during the most recent call to
       * member function setImage().
       *
       * @return The return value is vector of KeypointHarris instances.
       */
      std::vector< KeypointHarris<brick::common::Int32> >
      getKeypoints() const;


      /**
       * Return the keypoints detected during the most recent call to
       * member function setImage().
       *
       * @param iterator This argument must be a writable iterator
       * pointing to KeypoindHarris<brick::common::Int32>.
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
       * @return The return value is vector of KeypointHarris instances.
       */
      std::vector< KeypointHarris<FloatType> >
      getKeypointsGeneralPosition() const;


      /**
       * Return the keypoints detected during the most recent call to
       * member function setImage(), and use subpixel interpolation to
       * refine their positions.
       *
       * @param iterator This argument must be a writable iterator
       * pointing to KeypoindHarris<FloatType>.
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


      // Slight misnomer.  This actually computes products of
      // gradients within the region defined by corner0 and corner1.
      // Also updates corner0 and corner1 to reflect the region over
      // which the gradient calculation is valid.
      void
      computeGradients(
        Image<GRAY_SIGNED32> const& inImage,
        brick::numeric::Array2D<AccumulatedType>& gradientXX,
        brick::numeric::Array2D<AccumulatedType>& gradientXY,
        brick::numeric::Array2D<AccumulatedType>& gradientYY,
        brick::numeric::Index2D& corner0,
        brick::numeric::Index2D& corner1);


      // Give products of gradients, compute the Harris corner
      // indicator at each pixel within the region defined by corner0
      // and corner1.
      void
      computeHarrisIndicators(
        brick::numeric::Array2D<AccumulatedType> const& gradientXX,
        brick::numeric::Array2D<AccumulatedType> const& gradientXY,
        brick::numeric::Array2D<AccumulatedType> const& gradientYY,
        brick::numeric::Array2D<FloatType>& harrisIndicators,
        brick::numeric::Index2D const& corner0,
        brick::numeric::Index2D const& corner1);


      /* ======== Data members ========= */
      brick::numeric::Array2D<FloatType>       m_harrisIndicators;
      brick::numeric::Array2D<AccumulatedType> m_gradientXX;
      brick::numeric::Array2D<AccumulatedType> m_gradientXY;
      brick::numeric::Array2D<AccumulatedType> m_gradientYY;

      // These two members keep track of what portion of the image has
      // valid data from which to compute keypoints.  This search
      // region will dwindle as as an input image is processed because
      // each successive filter leaves a little more deadspace at the
      // edges of the image.
      brick::numeric::Index2D m_searchRegionCorner0;
      brick::numeric::Index2D m_searchRegionCorner1;

      // Parameters of the algorithm itself.
      FloatType m_kappa;
      FloatType m_sigma;
    };

  } // namespace computerVision

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/keypointSelectorHarris_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_KEYPOINTSELECTORHARRIS_HH */
