/**
***************************************************************************
* @file brick/computerVision/pixelHSV.hh
*
* Header file declaring the PixelHSV class template.
*
* Copyright (C) 2005-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_PIXELHSV_HH
#define BRICK_COMPUTERVISION_PIXELHSV_HH

#include <brick/common/types.hh>

namespace brick {

  namespace computerVision {
    
    template<class TYPE>
    struct PixelHSV
    {
    public:

      /** 
       * This constructor makes no guarantees about the color of the
       * pixel.
       */
      PixelHSV()
        : hue(), saturation(), value() {}

    
      /** 
       * This constructor explicitly sets the pixel value.
       * 
       * @param hue This argument specifies the hue value for the new pixel.
       * 
       * @param saturation This argument specifies the saturation
       * value for the new pixel.
       * 
       * @param value This argument specifies the value value for the new pixel.
       */
      PixelHSV(const TYPE& hueValue, const TYPE& saturationValue,
               const TYPE& valueValue)
        : hue(hueValue), saturation(saturationValue), value(valueValue) {}


      /** 
       * The destructor deletes cleans up for deletion.
       */
      ~PixelHSV() {}

    
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

      TYPE hue;
      TYPE saturation;
      TYPE value;


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


    typedef PixelHSV<brick::common::UnsignedInt8> PixelHSV8;
    typedef PixelHSV<brick::common::UnsignedInt16> PixelHSV16;
    typedef PixelHSV<brick::common::Int16> PixelHSVSigned16;
    typedef PixelHSV<brick::common::Int32> PixelHSVSigned32;
    typedef PixelHSV<brick::common::Float32> PixelHSVFloat32;
    typedef PixelHSV<brick::common::Float64> PixelHSVFloat64;


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
    inline PixelHSV<TYPE>
    operator-(const PixelHSV<TYPE>& pixel0, const PixelHSV<TYPE>& pixel1);

  
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
    operator==(const PixelHSV<TYPE>& pixel0, const PixelHSV<TYPE>& pixel1);

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
    PixelHSV<TYPE>::
    copyFromIterator(Iter& iter)
    {
      hue = *(iter++);
      saturation = *(iter++);
      value = *(iter++);
      return iter;
    }


    // This member function assigns the pixel component values, in
    // order, to consecutive iterator targets, incrementing the
    // iterator after each assignment.
    template<class TYPE>
    template<class Iter>
    inline Iter&
    PixelHSV<TYPE>::
    copyToIterator(Iter& iter)
    {
      *(iter++) = hue;
      *(iter++) = saturation;
      *(iter++) = value;
      return iter;
    }


    // This static member function indicates whether or not a pixel
    // instance is memory identical to a contiguous array of Component
    // type.
    template<class TYPE>
    inline bool
    PixelHSV<TYPE>::
    isContiguous()
    {
      TYPE dummy;
      TYPE* typePtr = &dummy;
      PixelHSV<TYPE>* pixelPtr = reinterpret_cast<PixelHSV<TYPE>*>(typePtr);
      PixelHSV<TYPE>* pixelPtr2 = pixelPtr + 1;
      return
        ((reinterpret_cast<TYPE*>(&(pixelPtr2->hue)) == &(typePtr[3]))
         && (reinterpret_cast<TYPE*>(&(pixelPtr2->saturation)) == &(typePtr[4]))
         && (reinterpret_cast<TYPE*>(&(pixelPtr2->value)) == &(typePtr[5])));
    }


    // This operator subtracts the values of the individual color
    // components of its arguments.
    template<class TYPE>
    inline PixelHSV<TYPE>
    operator-(const PixelHSV<TYPE>& pixel0, const PixelHSV<TYPE>& pixel1)
    {
      return PixelHSV<TYPE>(
        pixel0.hue - pixel1.hue,
        pixel0.saturation - pixel1.saturation,
        pixel0.value - pixel1.value);
    }


    // This operator returns true if the contents of the two argments
    // are identical, false otherwise.
    template<class TYPE>
    inline bool
    operator==(const PixelHSV<TYPE>& pixel0, const PixelHSV<TYPE>& pixel1)
    {
      return (pixel0.hue == pixel1.hue
              && pixel0.saturation == pixel1.saturation
              && pixel0.value == pixel1.value);
    }

  } // namespace computerVision
  
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_PIXELHSV_HH */
