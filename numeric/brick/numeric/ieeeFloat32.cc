/**
***************************************************************************
* @file brick/numeric/IEEEFloat32.cpp
*
* Source file defining the IEEEFloat32 class.
*
* Copyright (C) 2004-2007,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <algorithm>
#include <cmath>
#include <sstream>

#include <brick/numeric/ieeeFloat32.hh>

namespace brick {

  namespace numeric {

    // Default constructor initializes to 0.0;
    IEEEFloat32::
    IEEEFloat32()
      : m_value(),
        m_bytes()
    {
      this->checkTypes();
      this->setValue(0.0);
    }


    // This constructor initializes the IEEEFloat32 instance to the value
    // specified by its argument.
    IEEEFloat32::
    IEEEFloat32(FloatType value)
      : m_value(),
        m_bytes()
    {
      this->checkTypes();
      this->setValue(value);
    }


    // This constructor initializes the IEEEFloat32 instance using its
    // 32 bit binary representation.
    IEEEFloat32::
    IEEEFloat32(unsigned char byte0,
                unsigned char byte1,
                unsigned char byte2,
                unsigned char byte3)
      : m_value(),
        m_bytes()
    {
      this->checkTypes();
      this->setValue(byte0, byte1, byte2, byte3);
    }


    // This is the copy constructor.  It deep copies its argument.
    IEEEFloat32::
    IEEEFloat32(const IEEEFloat32& source)
      : m_value(source.m_value),
        m_bytes()
    {
      std::copy(&(source.m_bytes[0]), &(source.m_bytes[0]) + 4, &(m_bytes[0]));
    }


    // This member function returns the requested 8 bits byte from the
    // IEEE floating point representation.
    unsigned char
    IEEEFloat32::
    getByte(size_t index0)
    {
      // Check argument.
      if(index0 >= 4) {
        std::ostringstream message;
        message << "Index value, " << index0 << ", is out of bounds.";
        BRICK_THROW(common::IndexException, "IEEEFloat32::getByte(size_t)",
                    message.str().c_str());
      }
      return m_bytes[index0];
    }


    // This member function sets the IEEEFloat32 instance to the value
    // specified by its argument.
    void
    IEEEFloat32::
    setValue(FloatType value)
    {
      // Copy the input argument.
      m_value = value;

      // Handle special cases.
      // Note that this || is redundant.  x == 0.0 should imply x == -0.0.
      if((value == 0.0) || (value == -0.0)) {
        m_bytes[0] = 0x00;
        m_bytes[1] = 0x00;
        m_bytes[2] = 0x00;
        m_bytes[3] = 0x00;
      } else {
        // Do the conversion.
        this->floatToBinary(
          value, m_bytes[0], m_bytes[1], m_bytes[2], m_bytes[3]);
      }
    }


    // This member function sets the IEEEFloat32 instance using the
    // 32 bit binary representation.
    void
    IEEEFloat32::
    setValue(unsigned char byte0,
             unsigned char byte1,
             unsigned char byte2,
             unsigned char byte3)
    {
      // Remember the input arguments.
      m_bytes[0] = byte0;
      m_bytes[1] = byte1;
      m_bytes[2] = byte2;
      m_bytes[3] = byte3;

      // Handle special cases.
      if((byte0 == 0x00)
         && (byte1 == 0x00)
         && (byte2 == 0x00)
         && (byte3 == 0x00)) {
        m_value = 0.0;
      } else if((byte0 == 0x80)
                && (byte1 == 0x00)
                && (byte2 == 0x00)
                && (byte3 == 0x00)) {
        // Treat -0.0 just like 0.0.
        m_bytes[0] = 0x00;
        m_value = 0.0;
      } else {
        // Do the conversion.
        this->binaryToFloat(byte0, byte1, byte2, byte3, m_value);
      }
    }



    // This private member function implements the actual conversion
    // from binary representation to float.
    void
    IEEEFloat32::
    binaryToFloat(unsigned char byte0,
                  unsigned char byte1,
                  unsigned char byte2,
                  unsigned char byte3,
                  IEEEFloat32::FloatType& value)
    {
      // Start by calculating the mantissa using an int.  We'll divide
      // by 2^23 later.  Remember that the 2^0 bit is set by the IEEE
      // floating point definition.
      int mantissaByte0 = static_cast<int>(byte1 | 0x80);
      int mantissaByte1 = static_cast<int>(byte2);
      int mantissaByte2 = static_cast<int>(byte3);
      int mantissaAsInt =
        (mantissaByte0 << 16) | (mantissaByte1 << 8) | mantissaByte2;

      // Now include the sign.
      if(byte0 & 0x80) {
        mantissaAsInt *= -1;
      }

      // Recover the exponent.
      int exponentByte0 = static_cast<int>(byte0 & 0x7f);
      int exponentByte1 = static_cast<int>(byte1 & 0x80);
      int exponentAsInt = (exponentByte0 << 1) | (exponentByte1 >> 7);
      // Remember that IEEE format requires us to offset the exponent.
      exponentAsInt -= 127;

      // Now convert to floating point.  Hope that FloatType has enough
      // precision to do this without roundoff errors (true as long as
      // FloatType has at least as many bits as IEEE 32 bit float for both
      // exponent and mantissa).  We explicitly divide by 2^23 (rather
      // than simply subtracting 23 from the exponent) in case the local
      // representation of FloatType has only 8 exponent bits.
      FloatType mantissa =
        static_cast<FloatType>(mantissaAsInt)
        / static_cast<FloatType>(std::pow(2.0, 23.0));
      value = mantissa * static_cast<FloatType>(
        std::pow(2.0, static_cast<double>(exponentAsInt)));
    }



    // This private member function verifies that the compiler
    // built-in types have sufficient precision to implement the math
    // in this class.
    void
    IEEEFloat32::
    checkTypes()
    {
      if(sizeof(FloatType) < 4) {
        BRICK_THROW(common::RunTimeException, "IEEEFloat32::checkTypes()",
                    "FloatType has insufficient precision.");
      }

      // Probably need some more checking here.
    }


    // This private member function implements the actual conversion
    // from float to binary representation.
    void
    IEEEFloat32::
    floatToBinary(FloatType value,
                  unsigned char& byte0,
                  unsigned char& byte1,
                  unsigned char& byte2,
                  unsigned char& byte3)
    {
      // Extract the sign bit.
      if(value < 0.0) {
        byte0 = 0x80;
        value *= -1;
      } else {
        byte0 = 0x00;
      }

      // Extract the exponent.
      int exponentAsInt = 0;
      while(value >= 2.0) {
        value /= 2.0;
        exponentAsInt += 1;
      }
      while(value < 1.0) {
        value *= 2.0;
        exponentAsInt -= 1;
      }

      // Remember the offset required by IEEE format.
      exponentAsInt += 127;

      // Now we know the mantissa (it's what's left in value).
      int mantissaAsInt = static_cast<int>(value * pow(2.0, 23));

      // Go ahead and assign the bytes.
      byte0 |= ((exponentAsInt & 0xfe) >> 1);
      byte1 = ((exponentAsInt & 0x01) << 7) | ((mantissaAsInt & 0x7f0000) >> 16);
      byte2 = ((mantissaAsInt & 0x00ff00) >> 8);
      byte3 = (mantissaAsInt & 0x0000ff);
    }

  } // namespace numeric

} // namespace brick
