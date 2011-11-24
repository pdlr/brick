/**
***************************************************************************
* @file brick/computerVision/pixelBGRA.hh
*
* Header file declaring the PixelBGRA class template.
*
* Copyright (C) 2005-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_PIXELBGRA_HH
#define BRICK_COMPUTERVISION_PIXELBGRA_HH

#include <brick/common/types.hh>

namespace brick {

  namespace computerVision {
    
    template<class TYPE>
    struct PixelBGRA
    {
    public:

      /** 
       * This constructor makes no guarantees about the color of the
       * pixel.
       */
      PixelBGRA()
        : blue(), green(), red(), alpha() {}

    
      /** 
       * This constructor explicitly sets the pixel value.
       * 
       * @param blueValue This argument specifies the blue value for
       * the new pixel.
       * 
       * @param greenValue This argument specifies the green value for
       * the new pixel.
       * 
       * @param redValue This argument specifies the red value for the
       * new pixel.
       *
       * @param alphaValue This argument specifies the alpha value for
       * the new pixel.
       */
      PixelBGRA(const TYPE& blueValue, const TYPE& greenValue,
                const TYPE& redValue,  const TYPE& alphaValue)
        : blue(blueValue), green(greenValue), red(redValue),
          alpha(alphaValue) {}


      /** 
       * The destructor deletes cleans up for deletion.
       */
      ~PixelBGRA() {}


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

      TYPE blue;
      TYPE green;
      TYPE red;
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


    typedef PixelBGRA<brick::common::UnsignedInt8> PixelBGRA8;
    typedef PixelBGRA<brick::common::UnsignedInt16> PixelBGRA16;
    typedef PixelBGRA<brick::common::Int16> PixelBGRASigned16;
    typedef PixelBGRA<brick::common::Int32> PixelBGRASigned32;
    typedef PixelBGRA<brick::common::Float32> PixelBGRAFloat32;
    typedef PixelBGRA<brick::common::Float64> PixelBGRAFloat64;


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
    inline PixelBGRA<TYPE>
    operator-(const PixelBGRA<TYPE>& pixel0, const PixelBGRA<TYPE>& pixel1);

  
    /** 
     * This operator returns true if the contents of the two argments
     * are identical, false otherwise.
     * 
     * @param pixel0 This argument is the first pixel value to be compablue.
     * 
     * @param pixel1 This argument is the second pixel value to be compablue.
     * 
     * @return The return value indicates whether the two pixels have
     * identical values.
     */
    template<class TYPE>
    inline bool
    operator==(const PixelBGRA<TYPE>& pixel0, const PixelBGRA<TYPE>& pixel1);

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
    PixelBGRA<TYPE>::
    copyFromIterator(Iter& iter)
    {
      blue = *(iter++);
      green = *(iter++);
      red = *(iter++);
      alpha = *(iter++);
      return iter;
    }


    // This member function assigns the pixel component values, in
    // order, to consecutive iterator targets, incrementing the
    // iterator after each assignment.
    template<class TYPE>
    template<class Iter>
    inline Iter&
    PixelBGRA<TYPE>::
    copyToIterator(Iter& iter)
    {
      *(iter++) = blue;
      *(iter++) = green;
      *(iter++) = red;
      *(iter++) = alpha;
      return iter;
    }


    // This static member function indicates whether or not a pixel
    // instance is memory identical to a contiguous array of Component
    // type.
    template<class TYPE>
    inline bool
    PixelBGRA<TYPE>::
    isContiguous()
    {
      TYPE dummy;
      TYPE* typePtr = &dummy;
      PixelBGRA<TYPE>* pixelPtr = reinterpret_cast<PixelBGRA<TYPE>*>(typePtr);
      PixelBGRA<TYPE>* pixelPtr2 = pixelPtr + 1;
      return
        ((reinterpret_cast<TYPE*>(&(pixelPtr2->blue)) == &(typePtr[4]))
         && (reinterpret_cast<TYPE*>(&(pixelPtr2->green)) == &(typePtr[5]))
         && (reinterpret_cast<TYPE*>(&(pixelPtr2->red)) == &(typePtr[6]))
         && (reinterpret_cast<TYPE*>(&(pixelPtr2->alpha)) == &(typePtr[7])));
    }


    // This operator subtracts the values of the individual color
    // components of its arguments.
    template<class TYPE>
    inline PixelBGRA<TYPE>
    operator-(const PixelBGRA<TYPE>& pixel0, const PixelBGRA<TYPE>& pixel1)
    {
      return PixelBGRA<TYPE>(
        pixel0.blue - pixel1.blue,
        pixel0.green - pixel1.green,
        pixel0.red - pixel1.red,
        pixel0.alpha - pixel1.alpha);
    }


    // This operator returns true if the contents of the two argments
    // are identical, false otherwise.
    template<class TYPE>
    inline bool
    operator==(const PixelBGRA<TYPE>& pixel0, const PixelBGRA<TYPE>& pixel1)
    {
      return (pixel0.blue == pixel1.blue
              && pixel0.green == pixel1.green
              && pixel0.red == pixel1.red
              && pixel0.alpha == pixel1.alpha);
    }

  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_PIXELBGRA_HH */
