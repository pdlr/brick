/**
***************************************************************************
* @file brick/numeric/iEEEFloat32.hh
*
* Header file declaring the IEEEFloat32 class.
*
* Copyright (C) 2004-2007,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_IEEEFLOAT32_HH
#define BRICK_NUMERIC_IEEEFLOAT32_HH

#include <brick/common/exception.hh>

namespace brick {

  namespace numeric {

    /**
     ** The IEEEFloat32 class is for manipulating 32-bit IEEE floating
     ** point numbers.  This class is useful if, for example, you need
     ** to figure out what specific pattern of bits represents the
     ** number 1.347 in 32 bit IEEE floating point format, and you're
     ** not sure that your machine uses IEEE floating point
     ** representation.
     **
     ** The 32-bit IEEE floating point format is one sign bit (s),
     ** followed by an 8-bit biased exponent, followed by a fraction
     ** from the normalized mantissa (f).  This is broken into bytes as
     ** follows:
     **
     **   | seeeeeee | efffffff | ffffffff | ffffffff |
     **
     ** This arrangements of bits represents a real number, F, according
     ** to the following formula:
     **
     **   F = (-1)^s * 2E(e - 127) * (1.f)_2,
     **
     ** where the notation (1.f)_2 means "the 24 bit binary number
     ** consisting of 1, followed by the 23 bits of f, where the 23 bits
     ** of f represent the fractional part.  That is the first bit of f
     ** has a weight of 2E-1, the second has a weight of 2E-2, and so
     ** on.
     **
     ** WARNING: This class currently does not handle NaN or Inf values.
     ** Furthermore, the value -0.0 is silently converted to 0.0.
     **/
    class IEEEFloat32 {
    public:

      /**
       * This typedef specifies the native type that will be used to do
       * computations internally.  It must have at least as many
       * exponent and mantissa bits as IEEE 32-bit floating point
       * representation.
       */
      typedef float FloatType;


      /**
       * Default constructor initializes to 0.0;
       */
      IEEEFloat32();


      /**
       * This constructor initializes the IEEEFloat32 instance to the value
       * specified by its argument.
       *
       * @param value This argument specifies the value of the float in
       * question.
       */
      IEEEFloat32(FloatType value);


      /**
       * This constructor initializes the IEEEFloat32 instance using its
       * 32 bit binary representation.
       *
       * @param byte0 This argument represents the first 8 bits of the
       * binary representation (the sign bit and the first 7 exponent
       * bits).
       *
       * @param byte1 This argument represents the second 8 bits of the
       * binary representation (the final exponent bit and the first 7
       * mantissa bits).
       *
       * @param byte2 This argument represents the third 8 bits of the
       * binary representation (the 8th - 15th mantissa bits).
       *
       * @param byte3 This argument represents the final 8 bits of the
       * binary representation (the 16th - 23rd mantissa bits).
       */
      IEEEFloat32(unsigned char byte0,
                  unsigned char byte1,
                  unsigned char byte2,
                  unsigned char byte3);


      /**
       * This is the copy constructor.  It deep copies its argument.
       *
       * @param source This argument is the IEEEFloat32 instance to be
       * copied.
       */
      IEEEFloat32(const IEEEFloat32& source);


      /**
       * The destructor destroys the class instance and cleans up any
       * storage.
       */
      ~IEEEFloat32() {}


      /**
       * This conversion operator returns the float as a built-in type.
       *
       * @return The return value is a FloatType instance corresponding
       * to the number described by this class.
       */
      operator
      FloatType() const {return static_cast<FloatType>(m_value);}


      /**
       * This member function returns the requested 8 bits byte from the
       * IEEE floating point representation.
       *
       * @param index0 This argument specifies which byte to return. If
       * its value is zero, the sign bit and 1st 7 exponent bits will be
       * returned.  If its value is one, the final exponent bit and the
       * first 7 mantissa bits will be returned.  If its value is two,
       * the subsequent 8 mantissa bits will be returned.  If its value
       * is two, the final 8 mantissa bits will be returned.
       *
       * @return The return is an unsigned char containing the requested
       * 8 bits of the IEEE floating point representation.
       */
      unsigned char
      getByte(size_t index0);


      /**
       * This member function returns an instance of FloatType having
       * the same value as *this.  It is provided for those times when
       * an implicit type conversion isn't possible, and a static_cast
       * is too clunky.
       *
       * @return The return value is a FloatType instance having the
       * value represented by *this.
       */
      FloatType
      getFloat() {return static_cast<FloatType>(*this);}


      /**
       * This member function sets the IEEEFloat32 instance to the value
       * specified by its argument.
       *
       * @param value This argument specifies the value of the float in
       * question.
       */
      void
      setValue(FloatType value);


      /**
       * This member function sets the IEEEFloat32 instance using the
       * 32 bit binary representation.
       *
       * @param byte0 This argument represents the first 8 bits of the
       * binary representation (the sign bit and the first 7 exponent
       * bits).
       *
       * @param byte1 This argument represents the second 8 bits of the
       * binary representation (the final exponent bit and the first 7
       * mantissa bits).
       *
       * @param byte2 This argument represents the third 8 bits of the
       * binary representation (the 8th - 15th mantissa bits).
       *
       * @param byte3 This argument represents the final 8 bits of the
       * binary representation (the 16th - 23rd mantissa bits).
       */
      void
      setValue(unsigned char byte0,
               unsigned char byte1,
               unsigned char byte2,
               unsigned char byte3);


    private:
      /**
       * This private member function implements the actual conversion
       * from binary representation to FloatType.  On machines with
       * underlying big-endian IEEE floating point representation, you
       * could implement this as follows:
       *
       *   float floatValue;
       *   *reinterpret_cast<unsigned char*>(&floatValue) = byte0;
       *   *(reinterpret_cast<unsigned char*>(&floatValue) + 1) = byte1;
       *   *(reinterpret_cast<unsigned char*>(&floatValue) + 2) = byte2;
       *   *(reinterpret_cast<unsigned char*>(&floatValue) + 3) = byte3;
       *   value = floatValue;
       *
       * @param byte0 This argument is the first byte of the IEEE 32 bit
       * representation.
       *
       * @param byte1 This argument is the second byte of the IEEE 32 bit
       * representation.
       *
       * @param byte2 This argument is the third byte of the IEEE 32 bit
       * representation.
       *
       * @param byte3 This argument is the fourth byte of the IEEE 32 bit
       * representation.
       *
       * @param value This argument returns the recovered floating point
       * value.
       */
      void
      binaryToFloat(unsigned char byte0,
                    unsigned char byte1,
                    unsigned char byte2,
                    unsigned char byte3,
                    FloatType& value);


      /**
       * This private member function verifies that the compiler
       * built-in types have sufficient precision to implement the math
       * in this class.
       */
      void
      checkTypes();


      /**
       * This private member function implements the actual conversion
       * from float to binary representation.  On machines with
       * underlying big-endian IEEE floating point representation, you
       * could implement this as follows:
       *
       *   float floatValue = value;
       *   byte0 = *reinterpret_cast<unsigned char*>(&floatValue);
       *   byte1 = *(reinterpret_cast<unsigned char*>(&floatValue) + 1);
       *   byte2 = *(reinterpret_cast<unsigned char*>(&floatValue) + 2);
       *   byte3 = *(reinterpret_cast<unsigned char*>(&floatValue) + 3);
       *
       * @param value This argument is the floating point value.
       *
       * @param byte0 This argument returns the first byte of the IEEE 32 bit
       * representation.
       *
       * @param byte1 This argument returns the second byte of the IEEE 32 bit
       * representation.
       *
       * @param byte2 This argument returns the third byte of the IEEE 32 bit
       * representation.
       *
       * @param byte3 This argument returns the fourth byte of the IEEE 32 bit
       * representation.
       */
      void
      floatToBinary(FloatType value,
                    unsigned char& byte0,
                    unsigned char& byte1,
                    unsigned char& byte2,
                    unsigned char& byte3);


      /**
       * This member variable stores the floating point value of the
       * class instance.
       */
      FloatType m_value;

      /**
       * This member variable stores the binary representation of the
       * value of the class instance.
       */
      unsigned char m_bytes[4];
    }; // class IEEEFloat32

  } // namespace numeric

} // namespace brick

#endif // #ifndef BRICK_NUMERIC_IEEEFLOAT32_HH
