/**
***************************************************************************
* @file brick/computerVision/pixelRGBA.hh
*
* Header file declaring the PixelRGBA class template.
*
* Copyright (C) 2005-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_PIXELRGBA_HH
#define BRICK_COMPUTERVISION_PIXELRGBA_HH

#include <brick/common/types.hh>

namespace brick {

  namespace computerVision {

    template<class TYPE>
    struct PixelRGBA
    {
    public:

      /**
       * This constructor makes no guarantees about the color of the
       * pixel.
       */
      PixelRGBA()
        : red(), green(), blue(), alpha() {}


      /**
       * This constructor explicitly sets the pixel value.
       *
       * @param redValue This argument specifies the red value for the
       * new pixel.
       *
       * @param greenValue This argument specifies the green value for
       * the new pixel.
       *
       * @param blueValue This argument specifies the blue value for
       * the new pixel.
       *
       * @param alphaValue This argument specifies the alpha value for
       * the new pixel.
       */
      PixelRGBA(const TYPE& redValue, const TYPE& greenValue,
                const TYPE& blueValue,
                const TYPE& alphaValue)
        : red(redValue), green(greenValue), blue(blueValue),
          alpha(alphaValue) {}


      /**
       * The destructor deletes cleans up for deletion.
       */
      ~PixelRGBA() {}


      /**
       * This member function copies the pixel component values, in
       * order, from consecutive iterator targets, incrementing the
       * iterator after each copy.
       *
       * @param iter This argument is the iterator from which to copy the
       * pixel components.
       *
       * @return The return value is a reference to iter after the
       * copying and incrementing is done.
       */
      template<class Iter>
      inline Iter&
      copyFromIterator(Iter& iter);


      /**
       * This member function assigns the pixel component values, in
       * order, to consecutive iterator targets, incrementing the
       * iterator after each assignment.
       *
       * @param iter This argument is the iterator to which to copy the
       * pixel components.
       *
       * @return The return value is a reference to iter after the
       * copying and incrementing is done.
       */
      template<class Iter>
      inline Iter&
      copyToIterator(Iter& iter);


      /* ====== Public data members ====== */

      TYPE red;
      TYPE green;
      TYPE blue;
      TYPE alpha;


      /* ====== Static member functions ====== */

      /**
       * This static member function indicates whether or not a pixel
       * instance is memory identical to a contiguous array of Component
       * type.
       *
       * @return The return value is true if the pixel structure is not
       * padded by the compiler.
       */
      static inline bool
      isContiguous();

    };


    typedef PixelRGBA<brick::common::UnsignedInt8> PixelRGBA8;
    typedef PixelRGBA<brick::common::UnsignedInt16> PixelRGBA16;
    typedef PixelRGBA<brick::common::Int16> PixelRGBASigned16;
    typedef PixelRGBA<brick::common::Int32> PixelRGBASigned32;
    typedef PixelRGBA<brick::common::Float32> PixelRGBAFloat32;
    typedef PixelRGBA<brick::common::Float64> PixelRGBAFloat64;


    /**
     * This operator subtracts the values of the individual color
     * components of its arguments.
     *
     * @param pixel0 The color component values of pixel1 will be
     * subtracted from the color component values of this pixel.
     *
     * @param pixel1 The color component values of this pixel will be
     * subtracted from the color component values of pixel1.
     *
     * @return The return value is a pixel in which each color component
     * value is the difference of the corresponding values in the two
     * input pixels.
     */
    template<class TYPE>
    inline PixelRGBA<TYPE>
    operator-(const PixelRGBA<TYPE>& pixel0, const PixelRGBA<TYPE>& pixel1);


    /**
     * This operator returns true if the contents of the two argments
     * are identical, false otherwise.
     *
     * @param pixel0 This argument is the first pixel value to be compared.
     *
     * @param pixel1 This argument is the second pixel value to be compared.
     *
     * @return The return value indicates whether the two pixels have
     * identical values.
     */
    template<class TYPE>
    inline bool
    operator==(const PixelRGBA<TYPE>& pixel0, const PixelRGBA<TYPE>& pixel1);

  } // namespace computerVision

} // namespace brick

/* ============ Definitions of inline & template functions ============ */


namespace brick {

  namespace computerVision {

    // This member function copies the pixel component values, in
    // order, from consecutive iterator targets, incrementing the
    // iterator after each copy.
    template<class TYPE>
    template<class Iter>
    inline Iter&
    PixelRGBA<TYPE>::
    copyFromIterator(Iter& iter)
    {
      red = *(iter++);
      green = *(iter++);
      blue = *(iter++);
      alpha = *(iter++);
      return iter;
    }


    // This member function assigns the pixel component values, in
    // order, to consecutive iterator targets, incrementing the
    // iterator after each assignment.
    template<class TYPE>
    template<class Iter>
    inline Iter&
    PixelRGBA<TYPE>::
    copyToIterator(Iter& iter)
    {
      *(iter++) = red;
      *(iter++) = green;
      *(iter++) = blue;
      *(iter++) = alpha;
      return iter;
    }


    // This static member function indicates whether or not a pixel
    // instance is memory identical to a contiguous array of Component
    // type.
    template<class TYPE>
    inline bool
    PixelRGBA<TYPE>::
    isContiguous()
    {
      TYPE dummy;
      TYPE* typePtr = &dummy;
      PixelRGBA<TYPE>* pixelPtr = reinterpret_cast<PixelRGBA<TYPE>*>(typePtr);
      PixelRGBA<TYPE>* pixelPtr2 = pixelPtr + 1;
      return
        ((reinterpret_cast<TYPE*>(&(pixelPtr2->red)) == &(typePtr[4]))
         && (reinterpret_cast<TYPE*>(&(pixelPtr2->green)) == &(typePtr[5]))
         && (reinterpret_cast<TYPE*>(&(pixelPtr2->blue)) == &(typePtr[6]))
         && (reinterpret_cast<TYPE*>(&(pixelPtr2->alpha)) == &(typePtr[7])));
    }


    // This operator subtracts the values of the individual color
    // components of its arguments.
    template<class TYPE>
    inline PixelRGBA<TYPE>
    operator-(const PixelRGBA<TYPE>& pixel0, const PixelRGBA<TYPE>& pixel1)
    {
      return PixelRGBA<TYPE>(
        pixel0.red - pixel1.red,
        pixel0.green - pixel1.green,
        pixel0.blue - pixel1.blue,
        pixel0.alpha - pixel1.alpha);
    }


    // This operator returns true if the contents of the two argments
    // are identical, false otherwise.
    template<class TYPE>
    inline bool
    operator==(const PixelRGBA<TYPE>& pixel0, const PixelRGBA<TYPE>& pixel1)
    {
      return (pixel0.red == pixel1.red
              && pixel0.green == pixel1.green
              && pixel0.blue == pixel1.blue
              && pixel0.alpha == pixel1.alpha);
    }

  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_PIXELRGBA_HH */
