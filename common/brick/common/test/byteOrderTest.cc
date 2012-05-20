/**
***************************************************************************
* @file byteOrderTest.cc
* 
* Source file defining tests for exception trace code.
*
* Copyright (C) 2005 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <iostream>
#include <limits>
#include <vector>
#include <brick/common/byteOrder.hh>
#include <brick/common/types.hh>

namespace brick {

  namespace common {

    // We don't want to introduce a dependency on non-brick code for
    // unit testing, and the brick::test library is not available in
    // this context, so we just hack up some test functions.

    template <class Type>
    std::vector<Type>
    getTestSequence(unsigned int numberOfElements, Type increment)
    {
      std::vector<Type> testSequence;
      Type accumulator(0);
      for(unsigned int ii = 0; ii < numberOfElements; ++ii) {
        testSequence.push_back(accumulator);
        accumulator += increment;
      }
      return testSequence;
    }
    
      
    bool
    testGetByteOrder()
    {
      std::cout << "Testing getByteOrder..." << std::endl;

      std::vector<UnsignedInt16> testSequence =
        getTestSequence<UnsignedInt16>(2000, 27);
      
      for(unsigned int ii = 0; ii < testSequence.size(); ++ii) {
        UnsignedInt16 const testVal = testSequence[ii];
        UnsignedInt8 const highByte(testVal >> 8);
        UnsignedInt8 const lowByte(testVal & 0x0ff);

        UnsignedInt8 const* asBytes = (UnsignedInt8*)(&testVal);
        UnsignedInt8 const firstByte = asBytes[0];
        UnsignedInt8 const secondByte = asBytes[1];

        if(getByteOrder() == BRICK_BIG_ENDIAN) {
          if((firstByte != highByte) || (secondByte != lowByte)) {
            return false;
          }
        } else {
          if((firstByte != lowByte)  || (secondByte != highByte)) {
            return false;
          }
        }
      }

      // If we get this far, then all is well.
      return true;
    }


    bool
    testSwitchByteOrder()
    {
      std::cout << "Testing switchByteOrder..." << std::endl;

      // OK, the code below is duplicated and tedious, but we're going
      // for understandable, not elegant.
      
      // Get test sequences.
      unsigned int const numElements = 2000;
      UnsignedInt16 increment16 =
        (std::numeric_limits<UnsignedInt16>::max() - 1) / numElements;
      std::vector<UnsignedInt16> testSequence16 =
        getTestSequence<UnsignedInt16>(numElements, increment16);

      UnsignedInt32 increment32 =
        (std::numeric_limits<UnsignedInt32>::max() - 1) / numElements;
      std::vector<UnsignedInt32> testSequence32 =
        getTestSequence<UnsignedInt32>(numElements, increment32);

      UnsignedInt64 increment64 =
        (std::numeric_limits<UnsignedInt64>::max() - 1) / numElements;
      std::vector<UnsignedInt64> testSequence64 =
        getTestSequence<UnsignedInt64>(numElements, increment64);

      
      // Byte-swap the test sequences manually.
      std::vector<UnsignedInt16> swappedSequence16(testSequence16.size());
      for(unsigned int ii = 0; ii < testSequence16.size(); ++ii) {
        UnsignedInt16 const testVal = testSequence16[ii];
        UnsignedInt16 const byte0(testVal & 0x00ff);
        UnsignedInt16 const byte1(testVal >> 8);
        UnsignedInt16 const swappedVal = (byte0 << 8) + byte1;
        swappedSequence16[ii] = swappedVal;
      }
      
      std::vector<UnsignedInt32> swappedSequence32(testSequence32.size());
      for(unsigned int ii = 0; ii < testSequence32.size(); ++ii) {
        UnsignedInt32 const testVal = testSequence32[ii];
        UnsignedInt32 const byte0(testVal & 0x000000ff);
        UnsignedInt32 const byte1((testVal & 0x0000ff00) >> 8);
        UnsignedInt32 const byte2((testVal & 0x00ff0000) >> 16);
        UnsignedInt32 const byte3((testVal & 0xff000000) >> 24);
        UnsignedInt32 const swappedVal =
          (byte0 << 24) + (byte1 << 16) + (byte2 << 8) + byte3;
        swappedSequence32[ii] = swappedVal;
      }
      
      std::vector<UnsignedInt64> swappedSequence64(testSequence64.size());
      for(unsigned int ii = 0; ii < testSequence64.size(); ++ii) {
        UnsignedInt64 const testVal = testSequence64[ii];
        UnsignedInt64 const byte0(testVal & 0x00000000000000ffLL);
        UnsignedInt64 const byte1((testVal & 0x000000000000ff00LL) >> 8);
        UnsignedInt64 const byte2((testVal & 0x0000000000ff0000LL) >> 16);
        UnsignedInt64 const byte3((testVal & 0x00000000ff000000LL) >> 24);
        UnsignedInt64 const byte4((testVal & 0x000000ff00000000LL) >> 32);
        UnsignedInt64 const byte5((testVal & 0x0000ff0000000000LL) >> 40);
        UnsignedInt64 const byte6((testVal & 0x00ff000000000000LL) >> 48);
        UnsignedInt64 const byte7((testVal & 0xff00000000000000LL) >> 56);
        UnsignedInt64 const swappedVal =
          (byte0 << 56) + (byte1 << 48) + (byte2 << 40) + (byte3 << 32)
          + (byte4 << 24) + (byte5 << 16) + (byte6 << 8) + byte7;
        swappedSequence64[ii] = swappedVal;
      }

      // Figure out what we want the byteOrder.hh routines to do.
      ByteOrder fromOrder = BRICK_BIG_ENDIAN;
      ByteOrder toOrder = BRICK_LITTLE_ENDIAN;
      if(fromOrder != getByteOrder()) {
        std::swap(fromOrder, toOrder);
      }

      // Use the routines from byteOrder.hh to swap the byte order,
      // and make sure the result is correct.
      std::vector<UnsignedInt16> workspace16 = testSequence16;
      std::vector<UnsignedInt16> resultSequence16(testSequence16.size());
      switchByteOrder(&(workspace16[0]), workspace16.size(),
                      &(resultSequence16[0]), fromOrder, toOrder);
      if(!std::equal(resultSequence16.begin(), resultSequence16.end(),
                     swappedSequence16.begin())) {
        return false;
      }
      
      std::vector<UnsignedInt32> workspace32 = testSequence32;
      std::vector<UnsignedInt32> resultSequence32(testSequence32.size());
      switchByteOrder(&(workspace32[0]), workspace32.size(),
                      &(resultSequence32[0]), fromOrder, toOrder);
      if(!std::equal(resultSequence32.begin(), resultSequence32.end(),
                     swappedSequence32.begin())) {
        return false;
      }
      
      std::vector<UnsignedInt64> workspace64 = testSequence64;
      std::vector<UnsignedInt64> resultSequence64(testSequence64.size());
      switchByteOrder(&(workspace64[0]), workspace64.size(),
                      &(resultSequence64[0]), fromOrder, toOrder);
      if(!std::equal(resultSequence64.begin(), resultSequence64.end(),
                     swappedSequence64.begin())) {
        return false;
      }

      // Make sure the swapping didn't change the original data.
      if(!std::equal(workspace16.begin(), workspace16.end(),
                     testSequence16.begin())) {
        return false;
      }
      if(!std::equal(workspace32.begin(), workspace32.end(),
                     testSequence32.begin())) {
        return false;
      }
      if(!std::equal(workspace64.begin(), workspace64.end(),
                     testSequence64.begin())) {
        return false;
      }

      // Use the routines from byteOrder.hh to _not_ swap the byte order,
      // and make sure the result is correct.
      switchByteOrder(&(workspace16[0]), workspace16.size(),
                      &(resultSequence16[0]), fromOrder, fromOrder);
      if(!std::equal(resultSequence16.begin(), resultSequence16.end(),
                     testSequence16.begin())) {
        return false;
      }
      
      switchByteOrder(&(workspace32[0]), workspace32.size(),
                      &(resultSequence32[0]), fromOrder, fromOrder);
      if(!std::equal(resultSequence32.begin(), resultSequence32.end(),
                     testSequence32.begin())) {
        return false;
      }
      
      switchByteOrder(&(workspace64[0]), workspace64.size(),
                      &(resultSequence64[0]), fromOrder, fromOrder);
      if(!std::equal(resultSequence64.begin(), resultSequence64.end(),
                     testSequence64.begin())) {
        return false;
      }

      // Make sure the swapping didn't change the original data.
      if(!std::equal(workspace16.begin(), workspace16.end(),
                     testSequence16.begin())) {
        return false;
      }
      if(!std::equal(workspace32.begin(), workspace32.end(),
                     testSequence32.begin())) {
        return false;
      }
      if(!std::equal(workspace64.begin(), workspace64.end(),
                     testSequence64.begin())) {
        return false;
      }

      // Now swap the data in place, and make sure it comes out right.
      switchByteOrder(&(workspace16[0]), workspace16.size(),
                      fromOrder, toOrder);
      if(!std::equal(workspace16.begin(), workspace16.end(),
                     swappedSequence16.begin())) {
        return false;
      }
      switchByteOrder(&(workspace32[0]), workspace32.size(),
                      fromOrder, toOrder);
      if(!std::equal(workspace32.begin(), workspace32.end(),
                     swappedSequence32.begin())) {
        return false;
      }
      switchByteOrder(&(workspace64[0]), workspace64.size(),
                      fromOrder, toOrder);
      if(!std::equal(workspace64.begin(), workspace64.end(),
                     swappedSequence64.begin())) {
        return false;
      }

      // If we get this far, then all is well.
      return true;
    }

  } // namespace common

} // namespace brick


// int main(int argc, char** argv)
int main(int, char**)
{
  bool result = true;
  result &= brick::common::testGetByteOrder();
  result &= brick::common::testSwitchByteOrder();
  return (result ? 0 : 1);
}
