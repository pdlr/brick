/**
***************************************************************************
* @file brick/computerVision/pixelRGB.hh
*
* Header file declaring the PixelRGB class template.
*
* Copyright (C) 2005-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_PIXELRGB_HH
#define BRICK_COMPUTERVISION_PIXELRGB_HH

#include <brick/numeric/numericTraits.hh>
#include <brick/common/types.hh>

namespace brick {

  namespace computerVision {
    
    template<class Type>
    struct PixelRGB
    {
      /** 
       * This constructor makes no guarantees about the color of the
       * pixel.
       */
      PixelRGB()
        : red(), green(), blue() {}

    
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
       */
      PixelRGB(Type const& redValue, Type const& greenValue,
               Type const& blueValue)
        : red(redValue), green(greenValue), blue(blueValue) {}


      template <class Scalar>
      explicit PixelRGB(Scalar scalar)
        : red(static_cast<Type>(scalar)),
          green(static_cast<Type>(scalar)),
          blue(static_cast<Type>(scalar)) {}

      
      template <class OtherType>
      PixelRGB(PixelRGB<OtherType> const& other)
        : red(static_cast<Type>(other.red)),
          green(static_cast<Type>(other.green)),
          blue(static_cast<Type>(other.blue)) {}

      
      /** 
       * The destructor deletes cleans up for deletion.
       */
      ~PixelRGB() {}


      template <class Scalar>
      PixelRGB<Type>&
      operator=(Scalar scalar) {
        red = static_cast<Type>(scalar);
        green = static_cast<Type>(scalar);
        blue = static_cast<Type>(scalar);
        return *this;
      }

      
      template <class OtherType>
      PixelRGB<Type>&
      operator=(PixelRGB<OtherType> const& other) {
        red = static_cast<Type>(other.red);
        green = static_cast<Type>(other.green);
        blue = static_cast<Type>(other.blue);
        return *this;
      }
      

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

      Type red;
      Type green;
      Type blue;


      /* ====== Public arithmetic operators ====== */

      /** 
       * This operator adds the color component values of its argument to the 
       * corresponding values of its *this.
       * 
       * @param other The color component values of other will be
       * added to *this.
       * 
       * @return The return value is a reference to *this.
       */
      template <class OtherType>
      inline PixelRGB<Type>&
      operator+=(PixelRGB<OtherType> const& other);


      /** 
       * This operator adds the color component values of its argument to the 
       * corresponding values of its *this.
       * 
       * @param other The color component values of other will be
       * added to *this.
       * 
       * @return The return value is a reference to *this.
       */
      template <class OtherType>
      inline PixelRGB<Type>&
      operator-=(PixelRGB<OtherType> const& other);


      /** 
       * This operator divides each color component of *this by the
       * the value of its argument.
       * 
       * @param scalar This is the value by which to divide each color
       * component of *this.
       * 
       * @return The return value is a reference to *this.
       */
      template <class Scalar>
      inline PixelRGB<Type>&
      operator*=(Scalar const& scalar);


      /** 
       * This operator divides each color component of *this by the
       * the value of its argument.
       * 
       * @param scalar This is the value by which to divide each color
       * component of *this.
       * 
       * @return The return value is a reference to *this.
       */
      template <class Scalar>
      inline PixelRGB<Type>&
      operator/=(Scalar const& scalar);


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


    typedef PixelRGB<brick::common::UnsignedInt8> PixelRGB8;
    typedef PixelRGB<brick::common::UnsignedInt16> PixelRGB16;
    typedef PixelRGB<brick::common::Int16> PixelRGBSigned16;
    typedef PixelRGB<brick::common::Int32> PixelRGBSigned32;
    typedef PixelRGB<brick::common::Float32> PixelRGBFloat32;
    typedef PixelRGB<brick::common::Float64> PixelRGBFloat64;


    /** 
     * This operator multiplies the color component values of a pixel
     * by a scalar.
     * 
     * @param scalar Each color component will be multiplied by this
     * value.
     * 
     * @param pixel1 The color component values of this pixel will be
     * multiplied by the scalar.
     * 
     * @return The return value is the result of the multiplication.
     */
    template<class Scalar, class Type>
    inline PixelRGB<Type>
    operator*(Scalar scalar, PixelRGB<Type> const& pixel1);


    /** 
     * This operator adds the values of the individual color
     * components of its arguments.
     * 
     * @param pixel0 The color component values of pixel1 will be
     * added to the color component values of this pixel.
     * 
     * @param pixel1 The color component values of this pixel will be
     * added to the color component values of pixel0.
     * 
     * @return The return value is a pixel in which each color
     * component value is the sum of the corresponding values in the
     * two input pixels.
     */
    template<class Type>
    inline PixelRGB<Type>
    operator+(PixelRGB<Type> const& pixel0, PixelRGB<Type> const& pixel1);


    /** 
     * This operator subtracts the values of the individual color
     * components of its arguments.
     * 
     * @param pixel0 The color component values of pixel1 will be
     * subtracted from the color component values of this pixel.
     * 
     * @param pixel1 The color component values of this pixel will be
     * subtracted from the color component values of pixel0.
     * 
     * @return The return value is a pixel in which each color component
     * value is the difference of the corresponding values in the two
     * input pixels.
     */
    template<class Type>
    inline PixelRGB<Type>
    operator-(PixelRGB<Type> const& pixel0, PixelRGB<Type> const& pixel1);
    

    /** 
     * This operator multiplies the values of the individual color
     * components of its arguments.
     * 
     * @param pixel0 The color component values of this pixel will be
     * multiplied by the color component values of pixel1.
     * 
     * @param pixel1 The color component values of pixel0 will be
     * multiplied by the color component values of this pixel.
     * 
     * @return The return value is a pixel in which each color
     * component value is the product of the corresponding values in
     * the two input pixels.
     */
    template<class Type>
    inline PixelRGB<Type>
    operator*(PixelRGB<Type> const& pixel0, PixelRGB<Type> const& pixel1);


    /** 
     * This operator divides the values of the individual color
     * components of its arguments.
     * 
     * @param pixel0 The color component values this pixel will be
     * divided by the color component values of pixel1.
     * 
     * @param pixel1 The color component values of pixel0 will be
     * divided by the color component values of this pixel.
     * 
     * @return The return value is a pixel in which each color
     * component value is the dividend of the corresponding values in
     * the two input pixels.
     */
    template<class Type>
    inline PixelRGB<Type>
    operator/(PixelRGB<Type> const& pixel0, PixelRGB<Type> const& pixel1);


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
    template<class Type>
    inline bool
    operator==(PixelRGB<Type> const& pixel0, PixelRGB<Type> const& pixel1);

  } // namespace computerVision    

} // namespace brick
  
/* ============ Definitions of inline & template functions ============ */


namespace brick {

  namespace computerVision {
    
    // This member function copies the pixel component values, in
    // order, from consecutive iterator targets, incrementing the
    // iterator after each copy.
    template<class Type>
    template<class Iter>
    inline Iter&
    PixelRGB<Type>::
    copyFromIterator(Iter& iter)
    {
      red = *(iter++);
      green = *(iter++);
      blue = *(iter++);
      return iter;
    }


    // This member function assigns the pixel component values, in
    // order, to consecutive iterator targets, incrementing the
    // iterator after each assignment.
    template<class Type>
    template<class Iter>
    inline Iter&
    PixelRGB<Type>::
    copyToIterator(Iter& iter)
    {
      *(iter++) = red;
      *(iter++) = green;
      *(iter++) = blue;
      return iter;
    }



    // This operator adds the color component values of its argument to the 
    // corresponding values of its *this.
    template<class Type>
    template <class OtherType>
    inline PixelRGB<Type>&
    PixelRGB<Type>::
    operator+=(PixelRGB<OtherType> const& other)
    {
      red += other.red;
      green += other.green;
      blue += other.blue;
      return *this;
    }

    
    // This operator adds the color component values of its argument to the 
    // corresponding values of its *this.
    template<class Type>
    template <class OtherType>
    inline PixelRGB<Type>&
    PixelRGB<Type>::
    operator-=(PixelRGB<OtherType> const& other)
    {
      red -= other.red;
      green -= other.green;
      blue -= other.blue;
      return *this;
    }

    
    // This operator divides each color component of *this by the
    // the value of its argument.
    template<class Type>
    template <class Scalar>
    inline PixelRGB<Type>&
    PixelRGB<Type>::
    operator*=(Scalar const& scalar)
    {
      red *= scalar;
      green *= scalar;
      blue *= scalar;
      return *this;
    }

    
    // This operator divides each color component of *this by the
    // the value of its argument.
    template<class Type>
    template <class Scalar>
    inline PixelRGB<Type>&
    PixelRGB<Type>::
    operator/=(Scalar const& scalar)
    {
      red /= scalar;
      green /= scalar;
      blue /= scalar;
      return *this;
    }

    
    // This static member function indicates whether or not a pixel
    // instance is memory identical to a contiguous array of Component
    // type.
    template<class Type>
    inline bool
    PixelRGB<Type>::
    isContiguous()
    {
      Type dummy;
      Type* typePtr = &dummy;
      PixelRGB<Type>* pixelPtr = reinterpret_cast<PixelRGB<Type>*>(typePtr);
      PixelRGB<Type>* pixelPtr2 = pixelPtr + 1;
      return
        ((reinterpret_cast<Type*>(&(pixelPtr2->red)) == &(typePtr[3]))
         && (reinterpret_cast<Type*>(&(pixelPtr2->green)) == &(typePtr[4]))
         && (reinterpret_cast<Type*>(&(pixelPtr2->blue)) == &(typePtr[5])));
    }


    // This operator subtracts the values of the individual color
    // components of its arguments.
    template<class Scalar, class Type>
    inline PixelRGB<Type>
    operator*(Scalar scalar, PixelRGB<Type> const& pixel1)
    {
      return PixelRGB<Type>(
          scalar * pixel1.red,
          scalar * pixel1.green,
          scalar * pixel1.blue);
    }


    // This operator subtracts the values of the individual color
    // components of its arguments.
    template<class Type>
    inline PixelRGB<Type>
    operator+(PixelRGB<Type> const& pixel0, PixelRGB<Type> const& pixel1)
    {
      return PixelRGB<Type>(
        pixel0.red + pixel1.red,
        pixel0.green + pixel1.green,
        pixel0.blue + pixel1.blue);
    }


    // This operator subtracts the values of the individual color
    // components of its arguments.
    template<class Type>
    inline PixelRGB<Type>
    operator-(PixelRGB<Type> const& pixel0, PixelRGB<Type> const& pixel1)
    {
      return PixelRGB<Type>(
            pixel0.red - pixel1.red,
            pixel0.green - pixel1.green,
            pixel0.blue - pixel1.blue);
    }


    // This operator subtracts the values of the individual color
    // components of its arguments.
    template<class Type>
    inline PixelRGB<Type>
    operator*(PixelRGB<Type> const& pixel0, PixelRGB<Type> const& pixel1)
    {
      return PixelRGB<Type>(
        pixel0.red * pixel1.red,
        pixel0.green * pixel1.green,
        pixel0.blue * pixel1.blue);
    }


    // This operator subtracts the values of the individual color
    // components of its arguments.
    template<class Type>
    inline PixelRGB<Type>
    operator/(PixelRGB<Type> const& pixel0, PixelRGB<Type> const& pixel1)
    {
      return PixelRGB<Type>(
        pixel0.red / pixel1.red,
        pixel0.green / pixel1.green,
        pixel0.blue / pixel1.blue);
    }


    // This operator returns true if the contents of the two argments
    // are identical, false otherwise.
    template<class Type>
    inline bool
    operator==(PixelRGB<Type> const& pixel0, PixelRGB<Type> const& pixel1)
    {
      return (pixel0.red == pixel1.red
              && pixel0.green == pixel1.green
              && pixel0.blue == pixel1.blue);
    }

  } // namespace computerVision
  
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_PIXELRGB_HH */
