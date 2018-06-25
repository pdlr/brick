/**
***************************************************************************
* @file brick/computerVision/thresholderSauvola.hh
*
* Header file declaring a class that implements the adaptive text
* thresholding algorithm of Sauvola and Pietikainen.
*
* Copyright (C) 2017 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_THRESHOLDERSAUVOLA_HH
#define BRICK_COMPUTERVISION_THRESHOLDERSAUVOLA_HH

#include <brick/computerVision/image.hh>

#include <brick/numeric/boxIntegrator2D.hh>

namespace brick {

  namespace computerVision {

    /**
     ** This struct template controls the behavior of the
     ** ThresholderSauvola class.  If you need to customize the
     ** thresholderSauvola for your specific application, you can
     ** declare a similar struct, and use it as the second template
     ** argument of ThresholderSauvola.
     **/
    template<ImageFormat Format>
    struct ThresholderSauvolaConfig {
      typedef brick::common::Float64 FloatType;

      // Must be big enough to accumulate one window's worth of
      // squared pixel values.
      typedef uint32_t SumType;

      static uint8_t getBlackValue() {return 0;}
      static FloatType getMaxStdDev() {return 128.0;}
      static uint8_t getWhiteValue() {return 255;}
    };


    /**
     ** This class implements the adaptive text thresholding algorithm
     ** of Sauvola and Pietikainen [1].  It does not implement the
     ** corresponding non-textual algorithm, and it does not implement
     ** the region analysis and switching logic that the authors
     ** intended to select between the text thresholding algorithm and
     ** the non-textual thresholding algorithm.
     **
     ** Our implementation uses integral images to speed up the
     ** computation of local mean and variance, at the expense of
     ** increased memory footprint.
     **
     ** Use this class as follows:
     **
     ** @code
     **   Image<GRAY8> inputImage0 = readPGM8(getTestImageFileNamePGM0());
     **   uint32_t const windowRadius = 32;
     **   Float64 const kappa = 0.5;
     **   ThresholderSauvola<GRAY8> thresholder(windowRadius, kappa);
     **   Image<GRAY8> binaryImage = thresholder(inputImage0);
     ** @endcode
     **
     ** [1] J. Sauvola and M. Pietikainen, “Adaptive document image
     ** binarization,” Pattern Recognition 33(2), pp. 225–236, 2000.
     **/
    template < ImageFormat Format,
               class Config = ThresholderSauvolaConfig<Format> >
    class ThresholderSauvola {
    public:

      typedef typename ImageFormatTraits<Format>::PixelType PixelType;

      typedef typename Config::FloatType FloatType;
      typedef typename Config::SumType SumType;

      /**
       * Constructor.
       *
       * @param windowRadius This argument specifies the size of the
       * image adaptation region.  Thresholds will be computed based
       * on a region of support that is 2*windowRadius + 1 pixels
       * square.  Larger values make the threshold adapt more slowly
       * over the image.
       *
       * @param kappa This argument specifies the threshold
       * sensitivity.  Larger values of kappa make there be less text
       * (black) in the output image.
       */
      ThresholderSauvola(uint32_t windowRadius, FloatType kappa = 0.5);


      /**
       * Destructor.
       */
      virtual
      ~ThresholderSauvola() {}


      /**
       * Computes a thresholded image based on the input image
       * previously set by member function setImage().  Calling this
       * member function without have previously called setImage() is
       * an error.
       *
       * @return The return value is an Image<GRAY8> in which
       * forground (text) pixels are black (i.e., have pixel value 0)
       * and background pixels are white (i.e., have pixel value 255).
       */
      Image<GRAY8>
      computeBinaryImage();


      /**
       * Sets the image to be thresholded, and does preprocessing so
       * that subsequent calls to member function computeBinaryImage()
       * can execute quickly.  Note that inputImage is only shallow
       * copied.  Any changes to the original images that occur after
       * the call to setImage(), but before a call to
       * computeBinaryImage(), will affect the result of
       * computeBinaryImage().
       *
       * @param inputImage This argument is the image to be thresholded.
       */
      void
      setImage(Image<Format> const& inputImage);


      /**
       * Set the sensitivity of the thresholding algorithm, overriding
       * constructor argument and any previous calls and
       * setSensitivity().
       *
       * @param kappa This argument specifies the threshold
       * sensitivity.  Larger values of kappa make there be more text
       * (black) in the output image.
       */
      void
      setSensitivity(FloatType kappa) {
        this->m_kappa = kappa;
      }


      /**
       * Set the size of the adaptation region, overriding constructor
       * argument and any previous calls and setWindowRadius().
       *
       * @param windowRadius This argument specifies the size of the
       * image adaptation region.  Thresholds will be computed based
       * on a region of support that is 2*windowRadius + 1 pixels
       * square.  Larger values make the threshold adapt more slowly
       * over the image.
       */
      void
      setWindowRadius(uint32_t windowRadius) {
        this->m_windowRadius = windowRadius;
      }


      /**
       * This operator thresholds an input image.  It is equivalent to
       * calling member function setImage(), followed by member
       * function computeBinaryImage().
       *
       * @param inputImage This argument is the image to be thresholded.
       *
       * @return The return value is an Image<GRAY8> in which
       * forground (text) pixels are black (i.e., have pixel value 0)
       * and background pixels are white (i.e., have pixel value 255).
       */
      Image<GRAY8>
      operator()(Image<Format> const& inputImage) {
        this->setImage(inputImage);
        return this->computeBinaryImage();
      }

    private:

      // Quick and dirty routine to make sure an image ROI is
      // completely within the image (i.e., doesn't extend past one of
      // the image borders).  Assumes the image is bigger than the
      // ROI, beginCoord is smaller than endCoord, and lowerBound is
      // smaller than upperBound.
      void
      adjustWindowCoordinates(
        int32_t& beginCoord, int32_t& endCoord,
        int32_t lowerBound, int32_t upperBound, int32_t windowSize)
        {
          if(beginCoord < lowerBound) {
            beginCoord = lowerBound;
            endCoord = lowerBound + windowSize;
          } else if(endCoord > upperBound) {
            endCoord = upperBound;
            beginCoord = upperBound - windowSize;
          }
        }

      // This quantity gets computed a couple of times in the
      // implementation, so we abstract it out into a function.
      uint32_t
      getWindowSize() {return 2 * this->m_windowRadius + 1;}

      // These member variables affect the thresholding algorithm.
      FloatType m_kappa;
      uint32_t m_windowRadius;

      // This member variables holds a shallow copy of the input image.
      Image<Format> m_inputImage;

      brick::numeric::BoxIntegrator2D<PixelType, SumType>
        m_sumIntegrator;
      brick::numeric::BoxIntegrator2D<PixelType, SumType>
        m_squaredSumIntegrator;
    };

  } // namespace computerVision

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/thresholderSauvola_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_THRESHOLDERSAUVOLA_HH */
