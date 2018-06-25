/**
***************************************************************************
* @file brick/computerVision/pixelOperations.hh
*
* Header file declaring pixel-related functions that operate on built
* in types.  This lets us have, for example, multiplyPixel() function
* templates for grayscale images, in which the pixels are built in
* types.
*
* Copyright (C) 2015 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_PIXELOPERATIONS_HH
#define BRICK_COMPUTERVISION_PIXELOPERATIONS_HH

#include <stdint.h>

#include <brick/numeric/numericTraits.hh>
#include <brick/common/types.hh>

namespace brick {

  namespace computerVision {

    /**
     * This function is equivalent to operator*, but can be
     * specialized for specific values (e.g., when Multiplier is a
     * power of two) to improve efficiency.
     *
     * @param pixel0 This argument is the pixel to be multiplied.
     *
     * @return The return value is a copy of pixel0 in which each
     * element has been multiplied by Multiplier.
     */
    template<int Multiplier>
    inline brick::common::UInt8
    multiplyPixel(brick::common::UInt8 const& pixel0);

    template<int Multiplier>
    inline brick::common::UInt16
    multiplyPixel(brick::common::UInt16 const& pixel0);

    template<int Multiplier>
    inline brick::common::Float32
    multiplyPixel(brick::common::Float32 const& pixel0);

    template<int Multiplier>
    inline brick::common::Float64
    multiplyPixel(brick::common::Float64 const& pixel0);

    template<int Multiplier>
    inline PixelRGB8
    multiplyPixel(PixelRGB8 const& pixel0);

    template<int Multiplier>
    inline PixelRGB16
    multiplyPixel(PixelRGB16 const& pixel0);

    template<int Multiplier>
    inline PixelRGBFloat32
    multiplyPixel(PixelRGBFloat32 const& pixel0);


    /**
     * This function is equivalent to operator/, but can be
     * specialized for specific values (e.g., when Divisor is a
     * power of two) to improve efficiency.
     *
     * @param pixel0 This argument is the pixel to be divided.
     *
     * @return The return value is a copy of pixel0 in which each
     * element has been divided by Divisor.
     */
    template<int Divisor>
    inline brick::common::UInt8
    dividePixel(brick::common::UInt8 const& pixel0);

    template<int Divisor>
    inline brick::common::UInt16
    dividePixel(brick::common::UInt16 const& pixel0);

    template<int Divisor>
    inline brick::common::Float32
    dividePixel(brick::common::Float32 const& pixel0);

    template<int Divisor>
    inline brick::common::Float64
    dividePixel(brick::common::Float64 const& pixel0);

    template<int Divisor>
    inline PixelRGB8
    dividePixel(PixelRGB8 const& pixel0);

    template<int Divisor>
    inline PixelRGB16
    dividePixel(PixelRGB16 const& pixel0);

    template<int Divisor>
    inline PixelRGBFloat32
    dividePixel(PixelRGBFloat32 const& pixel0);

  } // namespace computerVision

} // namespace brick

/* ============ Definitions of inline & template functions ============ */


namespace brick {

  namespace computerVision {


    // This function is equivalent to operator*, but can be
    // specialized for specific values (e.g., when Multiplier is a
    // power of two) to improve efficiency.
    template<int Multiplier>
    inline brick::common::UInt8
    multiplyPixel(brick::common::UInt8 const& pixel0)
    {
      return pixel0 * static_cast<brick::common::UInt8>(Multiplier);
    }

    template<int Multiplier>
    inline brick::common::UInt16
    multiplyPixel(brick::common::UInt16 const& pixel0)
    {
      return pixel0 * static_cast<brick::common::UInt16>(Multiplier);
    }

    template<>
    inline brick::common::UInt16
    multiplyPixel<2>(brick::common::UInt16 const& pixel0)
    {
      return pixel0 << 1;
    }

    template<int Multiplier>
    inline brick::common::Float32
    multiplyPixel(brick::common::Float32 const& pixel0)
    {
      return pixel0 * static_cast<brick::common::Float32>(Multiplier);
    }

    template<int Multiplier>
    inline brick::common::Float64
    multiplyPixel(brick::common::Float64 const& pixel0)
    {
      return pixel0 * static_cast<brick::common::Float64>(Multiplier);
    }

    template<int Multiplier>
    inline PixelRGB8
    multiplyPixel(PixelRGB8 const& pixel0)
    {
      return PixelRGB8(
        pixel0.red * static_cast<brick::common::UInt8>(Multiplier),
        pixel0.green * static_cast<brick::common::UInt8>(Multiplier),
        pixel0.blue * static_cast<brick::common::UInt8>(Multiplier));
    }

    template<int Multiplier>
    inline PixelRGB16
    multiplyPixel(PixelRGB16 const& pixel0)
    {
      return PixelRGB16(
        pixel0.red * static_cast<brick::common::UInt16>(Multiplier),
        pixel0.green * static_cast<brick::common::UInt16>(Multiplier),
        pixel0.blue * static_cast<brick::common::UInt16>(Multiplier));
    }

    template<>
    inline PixelRGB16
    multiplyPixel<2>(PixelRGB16 const& pixel0)
    {
      return PixelRGB16(pixel0.red << 1, pixel0.green << 1, pixel0.blue << 1);
    }

    template<int Multiplier>
    inline PixelRGBFloat32
    multiplyPixel(PixelRGBFloat32 const& pixel0)
    {
      return PixelRGBFloat32(
        pixel0.red * static_cast<brick::common::Float32>(Multiplier),
        pixel0.green * static_cast<brick::common::Float32>(Multiplier),
        pixel0.blue * static_cast<brick::common::Float32>(Multiplier));
    }


    // This function is equivalent to operator/, but can be
    // specialized for specific values (e.g., when Divisor is a
    // power of two) to improve efficiency.
    template<int Divisor>
    inline brick::common::UInt8
    dividePixel(brick::common::UInt8 const& pixel0)
    {
      return pixel0 / static_cast<brick::common::UInt8>(Divisor);
    }

    template<int Divisor>
    inline brick::common::UInt16
    dividePixel(brick::common::UInt16 const& pixel0)
    {
      return pixel0 / static_cast<brick::common::UInt16>(Divisor);
    }

    template<>
    inline brick::common::UInt16
    dividePixel<4>(brick::common::UInt16 const& pixel0)
    {
      return pixel0 >> 2;
    }

    template<int Divisor>
    inline brick::common::Float32
    dividePixel(brick::common::Float32 const& pixel0)
    {
      return pixel0 / static_cast<brick::common::Float32>(Divisor);
    }

    template<int Divisor>
    inline brick::common::Float64
    dividePixel(brick::common::Float64 const& pixel0)
    {
      return pixel0 / static_cast<brick::common::Float64>(Divisor);
    }

    template<int Divisor>
    inline PixelRGB8
    dividePixel(PixelRGB8 const& pixel0)
    {
      return PixelRGB8(
        pixel0.red / static_cast<brick::common::UInt8>(Divisor),
        pixel0.green / static_cast<brick::common::UInt8>(Divisor),
        pixel0.blue / static_cast<brick::common::UInt8>(Divisor));
    }

    template<int Divisor>
    inline PixelRGB16
    dividePixel(PixelRGB16 const& pixel0)
    {
      return PixelRGB16(
        pixel0.red / static_cast<brick::common::UInt16>(Divisor),
        pixel0.green / static_cast<brick::common::UInt16>(Divisor),
        pixel0.blue / static_cast<brick::common::UInt16>(Divisor));
    }

    template<>
    inline PixelRGB16
    dividePixel<4>(PixelRGB16 const& pixel0)
    {
      return PixelRGB16(pixel0.red >> 2, pixel0.green >> 2, pixel0.blue >> 2);
    }

    template<int Divisor>
    inline PixelRGBFloat32
    dividePixel(PixelRGBFloat32 const& pixel0)
    {
      return PixelRGBFloat32(
        pixel0.red / static_cast<brick::common::Float32>(Divisor),
        pixel0.green / static_cast<brick::common::Float32>(Divisor),
        pixel0.blue / static_cast<brick::common::Float32>(Divisor));
    }

  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_PIXELOPERATIONS_HH */
