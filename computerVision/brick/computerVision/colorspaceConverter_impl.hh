/**
***************************************************************************
* @file brick/computerVision/colorspaceConverter_impl.hh
*
* Header file defining ColorspaceConverter class template.
*
* Copyright (C) 2006-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_COLORSPACECONVERTER_IMPL_HH
#define BRICK_COMPUTERVISION_COLORSPACECONVERTER_IMPL_HH

// This file is included by colorspaceConverter.hh, and should not be
// directly included by user code, so no need to include
// colorspaceConverter.hh here.
//
// #include <brick/computerVision/colorspaceConverter.hh>

#include <cmath>

#include <brick/common/mathFunctions.hh>

namespace brick {

  namespace computerVision {

    namespace privateCode {

      // Here we define file-scoped versions of common conversion
      // because dispatching between specializations of
      // ColorspaceConverter (to avoid code duplication) is too hard
      // for the compiler.  Instead, we dispatch to this function,
      // which is not templated.  Life gets much easiser.


      // This conversion follows the wikipedia article "HSL color
      // space" as of 2009-02-19, with the exceptions that all
      // resulting components are scaled [0.0, 1.0], and that the
      // resulting hue value is rotated 60 degrees to match the page's
      // contents as of 2014-08-21.
      inline void
      doColorspaceConversion(const Image<RGB_FLOAT64>::PixelType& inputPixel,
                             Image<HSV_FLOAT64>::PixelType& outputPixel)
      {
        // "value" is the simplest to compute... it's just the max RGB value.
        brick::common::Float64 maxVal = static_cast<brick::common::Float64>(
          std::max(inputPixel.red,
                   std::max(inputPixel.green, inputPixel.blue)));
        outputPixel.value = maxVal;

        // Handle first special case up-front.
        if(maxVal == 0.0) {
          outputPixel.hue = brick::common::Float64(0.0);
          outputPixel.saturation = brick::common::Float64(0.0);
          return;
        }

        // Now we can compute "saturation".
        brick::common::Float64 minVal = static_cast<brick::common::Float64>(
          std::min(inputPixel.red,
                   std::min(inputPixel.green, inputPixel.blue)));
        brick::common::Float64 delta = maxVal - minVal;
        outputPixel.saturation = delta / maxVal;

        // Handle second special case up-front.
        if(delta == 0.0) {
          outputPixel.hue = brick::common::Float64(0.0);
          return;
        }

        // Warning(xxx): our definition of hue goes from 0.0 to 1.0, not
        // 0.0 to 360.0, as is traditional.
        if(inputPixel.red == maxVal) {
          outputPixel.hue =
            brick::common::Float64(1.0 / 6.0)
            + (static_cast<brick::common::Float64>(inputPixel.green)
               - static_cast<brick::common::Float64>(inputPixel.blue)) / (6.0 * delta);
        } else if(inputPixel.green == maxVal) {
          outputPixel.hue =
            brick::common::Float64(0.5)
            + (static_cast<brick::common::Float64>(inputPixel.blue)
               - static_cast<brick::common::Float64>(inputPixel.red)) / (6.0 * delta);
        } else {
          outputPixel.hue =
            brick::common::Float64(5.0 / 6.0)
            + (static_cast<brick::common::Float64>(inputPixel.red)
               - static_cast<brick::common::Float64>(inputPixel.green)) / (6.0 * delta);
        }

#if 1
        // Previous code doesn't correspond to wikipedia anymore.
        // Possibly the 2009-02-19 page has been updated?  New
        // definition of hue is rotated 60 degrees from previous.
        outputPixel.hue -= 1.0 / 6.0;
        outputPixel.hue = ((outputPixel.hue < 0.0)
                           ? (outputPixel.hue + 1.0) : outputPixel.hue);
#endif
      }


      // This conversion inverts the RGB->HSV conversion defined above.
      inline void
      doColorspaceConversion(const Image<HSV_FLOAT64>::PixelType& inputPixel,
                             Image<RGB_FLOAT64>::PixelType& outputPixel)
      {
        brick::common::Float64 chroma =
          inputPixel.value * inputPixel.saturation;
        brick::common::Float64 huePrime =
          inputPixel.hue * 6.0;

        // Make sure huePrime is positive.
        while(huePrime < 0.0) {
          huePrime += 6.0;
        }

        // Compute huePrime % 2
        double remainder = huePrime;
        while(remainder >= 2.0) {
          remainder -= 2.0;
        }

        brick::common::Float64 xValue =
          chroma * (1.0 - brick::common::absoluteValue(remainder - 1));

        brick::common::Float64 red = 0.0;
        brick::common::Float64 green = 0.0;
        brick::common::Float64 blue = 0.0;
        if(huePrime < 1.0) {
          red = chroma;
          green = xValue;
          blue = 0.0;
        } else if (huePrime < 2.0) {
          red = xValue;
          green = chroma;
          blue = 0.0;
        } else if (huePrime < 3.0) {
          red = 0.0;
          green = chroma;
          blue = xValue;
        } else if (huePrime < 4.0) {
          red = 0.0;
          green = xValue;
          blue = chroma;
        } else if (huePrime < 5.0) {
          red = xValue;
          green = 0.0;
          blue = chroma;
        } else {
          red = chroma;
          green = 0.0;
          blue = xValue;
        }

        // Now update to reflect pixel brightness.
        double increment = inputPixel.value - chroma;
        outputPixel.red = red + increment;
        outputPixel.green = green + increment;
        outputPixel.blue = blue + increment;
      }

    } // namespace privateCode


    template<>
    inline
    void
    ColorspaceConverter<GRAY1, GRAY8>::
    operator()(const Image<GRAY1>::PixelType& inputPixel,
               Image<GRAY8>::PixelType& outputPixel)
    {
      if(inputPixel) {
        outputPixel = Image<GRAY8>::PixelType(255);
      } else {
        outputPixel = Image<GRAY8>::PixelType(0);
      }
    }


    template<>
    inline
    void
    ColorspaceConverter<GRAY1, RGB8>::
    operator()(const Image<GRAY1>::PixelType& inputPixel,
               Image<RGB8>::PixelType& outputPixel)
    {
      if(inputPixel) {
        outputPixel.red = ImageFormatTraits<RGB8>::ComponentType(255);
        outputPixel.green = ImageFormatTraits<RGB8>::ComponentType(255);
        outputPixel.blue = ImageFormatTraits<RGB8>::ComponentType(255);
      } else {
        outputPixel.red = ImageFormatTraits<RGB8>::ComponentType(0);
        outputPixel.green = ImageFormatTraits<RGB8>::ComponentType(0);
        outputPixel.blue = ImageFormatTraits<RGB8>::ComponentType(0);
      }
    }


    template<>
    inline
    void
    ColorspaceConverter<GRAY8, GRAY16>::
    operator()(const Image<GRAY8>::PixelType& inputPixel,
               Image<GRAY16>::PixelType& outputPixel)
    {
      outputPixel = Image<GRAY16>::PixelType(inputPixel) << 8;
    }


    template<>
    inline
    void
    ColorspaceConverter<GRAY8, RGB8>::
    operator()(const Image<GRAY8>::PixelType& inputPixel,
               Image<RGB8>::PixelType& outputPixel)
    {
      outputPixel.red = inputPixel;
      outputPixel.green = inputPixel;
      outputPixel.blue = inputPixel;
    }


    template<>
    inline
    void
    ColorspaceConverter<GRAY16, GRAY8>::
    operator()(const Image<GRAY16>::PixelType& inputPixel,
               Image<GRAY8>::PixelType& outputPixel)
    {
      outputPixel = Image<GRAY8>::PixelType(inputPixel >> 8);
    }


    template<>
    inline
    void
    ColorspaceConverter<RGB8, GRAY8>::
    operator()(const Image<RGB8>::PixelType& inputPixel,
               Image<GRAY8>::PixelType& outputPixel)
    {
      double accumulator = (0.3 * inputPixel.red
                            + 0.59 * inputPixel.green
                            + 0.11 * inputPixel.blue);
      outputPixel = static_cast<Image<GRAY8>::PixelType>(accumulator + 0.5);
    }


    template<>
    inline
    void
    ColorspaceConverter<RGB8, GRAY16>::
    operator()(const Image<RGB8>::PixelType& inputPixel,
               Image<GRAY16>::PixelType& outputPixel)
    {
      double accumulator = (0.3 * inputPixel.red
                            + 0.59 * inputPixel.green
                            + 0.11 * inputPixel.blue) * 256;
      outputPixel = static_cast<Image<GRAY16>::PixelType>(accumulator + 0.5);
    }


    template<>
    inline
    void
    ColorspaceConverter<RGB8, GRAY_FLOAT64>::
    operator()(const Image<RGB8>::PixelType& inputPixel,
               Image<GRAY_FLOAT64>::PixelType& outputPixel)
    {
      double accumulator = (0.3 * inputPixel.red
                            + 0.59 * inputPixel.green
                            + 0.11 * inputPixel.blue);
      outputPixel =
        static_cast<Image<GRAY_FLOAT64>::PixelType>(accumulator + 0.5);
    }


    template<>
    inline
    void
    ColorspaceConverter<RGB8, RGB16>::
    operator()(const Image<RGB8>::PixelType& inputPixel,
               Image<RGB16>::PixelType& outputPixel)
    {
      outputPixel.red = static_cast<brick::common::UnsignedInt16>(inputPixel.red) << 8;
      outputPixel.green = static_cast<brick::common::UnsignedInt16>(inputPixel.green) << 8;
      outputPixel.blue = static_cast<brick::common::UnsignedInt16>(inputPixel.blue) << 8;
    }


    template<>
    inline
    void
    ColorspaceConverter<RGB8, RGB_FLOAT32>::
    operator()(const Image<RGB8>::PixelType& inputPixel,
               Image<RGB_FLOAT32>::PixelType& outputPixel)
    {
      outputPixel.red = static_cast<brick::common::Float32>(inputPixel.red);
      outputPixel.green = static_cast<brick::common::Float32>(inputPixel.green);
      outputPixel.blue = static_cast<brick::common::Float32>(inputPixel.blue);
    }


    template<>
    inline
    void
    ColorspaceConverter<RGB8, RGB_FLOAT64>::
    operator()(const Image<RGB8>::PixelType& inputPixel,
               Image<RGB_FLOAT64>::PixelType& outputPixel)
    {
      outputPixel.red = static_cast<brick::common::Float64>(inputPixel.red);
      outputPixel.green = static_cast<brick::common::Float64>(inputPixel.green);
      outputPixel.blue = static_cast<brick::common::Float64>(inputPixel.blue);
    }


    template<>
    inline
    void
    ColorspaceConverter<RGB8, HSV_FLOAT64>::
    operator()(const Image<RGB8>::PixelType& inputPixel,
               Image<HSV_FLOAT64>::PixelType& outputPixel)
    {
      // Dispatch to the static version of this converter.  Note
      // that we're not rescaling our RGB values [0.0, 1.0].  This
      // should have no effect on the S & V values, which are
      // essentually ratios of R, G, and B anyway.  We'll rescale
      // value later, thus saving two floating point operations.
      privateCode::doColorspaceConversion(
        Image<RGB_FLOAT64>::PixelType(inputPixel.red, inputPixel.green,
                                      inputPixel.blue),
        outputPixel);

      // Scale value 0.0 - 1.0, consistent with H and S.
      outputPixel.value /= 255.0;
    }


    template<>
    inline
    void
    ColorspaceConverter<RGB8, YIQ_FLOAT64>::
    operator()(const Image<RGB8>::PixelType& inputPixel,
               Image<YIQ_FLOAT64>::PixelType& outputPixel)
    {
      outputPixel.luma = ((0.299 / 255.0) * inputPixel.red
                          + (0.587 / 255.0) * inputPixel.green
                          + (0.114 / 255.0) * inputPixel.blue);
      outputPixel.inPhase = ((0.595716 / 255.0)* inputPixel.red
                             - (0.274453 / 255.0) * inputPixel.green
                             - (0.321263 / 255.0) * inputPixel.blue);
      outputPixel.quadrature = ((0.211456 / 255.0) * inputPixel.red
                                - (0.522591 / 255.0) * inputPixel.green
                                + (0.311135 / 255.0) * inputPixel.blue);
    }


    template<>
    inline
    void
    ColorspaceConverter<RGB8, BGRA8>::
    operator()(const Image<RGB8>::PixelType& inputPixel,
               Image<BGRA8>::PixelType& outputPixel)
    {
      outputPixel.red = inputPixel.red;
      outputPixel.green = inputPixel.green;
      outputPixel.blue = inputPixel.blue;
      outputPixel.alpha = brick::common::UnsignedInt8(255);
    }


    template<>
    inline
    void
    ColorspaceConverter<RGB8, RGBA8>::
    operator()(const Image<RGB8>::PixelType& inputPixel,
               Image<RGBA8>::PixelType& outputPixel)
    {
      outputPixel.red = inputPixel.red;
      outputPixel.green = inputPixel.green;
      outputPixel.blue = inputPixel.blue;
      outputPixel.alpha = brick::common::UnsignedInt8(255);
    }


    template<>
    inline
    void
    ColorspaceConverter<RGB16, GRAY8>::
    operator()(const Image<RGB16>::PixelType& inputPixel,
               Image<GRAY8>::PixelType& outputPixel)
    {
      // Scale by 256 to avoid overflow.
      double accumulator = (0.3 * inputPixel.red
                            + 0.59 * inputPixel.green
                            + 0.11 * inputPixel.blue) / 256.0;
      outputPixel = static_cast<Image<GRAY8>::PixelType>(accumulator + 0.5);
    }


    template<>
    inline
    void
    ColorspaceConverter<RGB16, GRAY16>::
    operator()(const Image<RGB16>::PixelType& inputPixel,
               Image<GRAY16>::PixelType& outputPixel)
    {
      // Scale by 256 to avoid overflow.
      double accumulator = (0.3 * inputPixel.red
                            + 0.59 * inputPixel.green
                            + 0.11 * inputPixel.blue);
      outputPixel = static_cast<Image<GRAY16>::PixelType>(accumulator + 0.5);
    }


    template<>
    inline
    void
    ColorspaceConverter<RGB16, RGB8>::
    operator()(const Image<RGB16>::PixelType& inputPixel,
               Image<RGB8>::PixelType& outputPixel)
    {
      outputPixel.red = static_cast<brick::common::UnsignedInt8>(inputPixel.red >> 8);
      outputPixel.green = static_cast<brick::common::UnsignedInt8>(inputPixel.green >> 8);
      outputPixel.blue = static_cast<brick::common::UnsignedInt8>(inputPixel.blue >> 8);
    }


    template<>
    inline
    void
    ColorspaceConverter<RGB_FLOAT32, GRAY8>::
    operator()(const Image<RGB_FLOAT32>::PixelType& inputPixel,
               Image<GRAY8>::PixelType& outputPixel)
    {
      brick::common::Float32 accumulator = (0.3 * inputPixel.red
                                            + 0.59 * inputPixel.green
                                            + 0.11 * inputPixel.blue);
      outputPixel =
        static_cast<Image<GRAY8>::PixelType>(accumulator + 0.5);
    }


    template<>
    inline
    void
    ColorspaceConverter<RGB_FLOAT32, RGB8>::
    operator()(const Image<RGB_FLOAT32>::PixelType& inputPixel,
               Image<RGB8>::PixelType& outputPixel)
    {
      outputPixel.red = static_cast<brick::common::UnsignedInt8>(inputPixel.red + 0.5);
      outputPixel.green = static_cast<brick::common::UnsignedInt8>(inputPixel.green + 0.5);
      outputPixel.blue = static_cast<brick::common::UnsignedInt8>(inputPixel.blue + 0.5);
    }


    template<>
    inline
    void
    ColorspaceConverter<RGB_FLOAT64, RGB8>::
    operator()(const Image<RGB_FLOAT64>::PixelType& inputPixel,
               Image<RGB8>::PixelType& outputPixel)
    {
      outputPixel.red = static_cast<brick::common::UnsignedInt8>(inputPixel.red + 0.5);
      outputPixel.green = static_cast<brick::common::UnsignedInt8>(inputPixel.green + 0.5);
      outputPixel.blue = static_cast<brick::common::UnsignedInt8>(inputPixel.blue + 0.5);
    }


    template<>
    inline
    void
    ColorspaceConverter<RGB_FLOAT64, HSV_FLOAT64>::
    operator()(const Image<RGB_FLOAT64>::PixelType& inputPixel,
               Image<HSV_FLOAT64>::PixelType& outputPixel)
    {
      privateCode::doColorspaceConversion(inputPixel, outputPixel);
    }


    template<>
    inline
    void
    ColorspaceConverter<RGB_FLOAT64, YIQ_FLOAT64>::
    operator()(const Image<RGB_FLOAT64>::PixelType& inputPixel,
               Image<YIQ_FLOAT64>::PixelType& outputPixel)
    {
      outputPixel.luma = (0.299 * inputPixel.red
                          + 0.587 * inputPixel.green
                          + 0.114 * inputPixel.blue);
      outputPixel.inPhase = (0.595716 * inputPixel.red
                             - 0.274453 * inputPixel.green
                             - 0.321263 * inputPixel.blue);
      outputPixel.quadrature = (0.211456 * inputPixel.red
                                - 0.522591 * inputPixel.green
                                + 0.311135 * inputPixel.blue);
    }


    template<>
    inline
    void
    ColorspaceConverter<BGRA8, RGB8>::
    operator()(const Image<BGRA8>::PixelType& inputPixel,
               Image<RGB8>::PixelType& outputPixel)
    {
      outputPixel.red = inputPixel.red;
      outputPixel.green = inputPixel.green;
      outputPixel.blue = inputPixel.blue;
    }


    template<>
    inline
    void
    ColorspaceConverter<RGBA8, RGB8>::
    operator()(const Image<RGBA8>::PixelType& inputPixel,
               Image<RGB8>::PixelType& outputPixel)
    {
      outputPixel.red = inputPixel.red;
      outputPixel.green = inputPixel.green;
      outputPixel.blue = inputPixel.blue;
    }


    template<>
    inline
    void
    ColorspaceConverter<HSV_FLOAT64, RGB8>::
    operator()(const Image<HSV_FLOAT64>::PixelType& inputPixel,
               Image<RGB8>::PixelType& outputPixel)
    {
      PixelRGBFloat64 tempResult;
      privateCode::doColorspaceConversion(inputPixel, tempResult);
      outputPixel.red =
        static_cast<common::UInt8>(tempResult.red * 255.0 + 0.5);
      outputPixel.green =
        static_cast<common::UInt8>(tempResult.green * 255.0 + 0.5);
      outputPixel.blue =
        static_cast<common::UInt8>(tempResult.blue * 255.0 + 0.5);
    }



    template<>
    inline
    void
    ColorspaceConverter<HSV_FLOAT64, RGB_FLOAT64>::
    operator()(const Image<HSV_FLOAT64>::PixelType& inputPixel,
               Image<RGB_FLOAT64>::PixelType& outputPixel)
    {
      privateCode::doColorspaceConversion(inputPixel, outputPixel);
    }



  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_COLORSPACECONVERTER_IMPL_HH */
