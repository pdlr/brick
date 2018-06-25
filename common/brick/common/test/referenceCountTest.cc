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
#include <brick/common/referenceCount.hh>
#include <brick/common/types.hh>

namespace brick {

  namespace common {

    // We don't want to introduce a dependency on non-brick code for
    // unit testing, and the brick::test library is not available in
    // this context, so we just hack up some test functions.


    bool
    checkState(ReferenceCount const& referenceCount,
               bool isCounted, bool isShared, int count)
    {
      if(referenceCount.isCounted() != isCounted) {
        return false;
      }
      if(referenceCount.isShared() != isShared) {
        return false;
      }
      if(referenceCount.getCount() != count) {
        return false;
      }
      return true;
    }


    bool
    testConstructor()
    {
      std::cout << "Testing ReferenceCount::ReferenceCount()..." << std::endl;

      ReferenceCount count0(0);
      if(!checkState(count0, false, false, 0)) {
        return false;
      }

      ReferenceCount count1(1);
      if(!checkState(count1, true, false, 1)) {
        return false;
      }

      ReferenceCount count2(2);
      if(!checkState(count2, true, true, 2)) {
        return false;
      }

      // If we get this far, then all is well.
      return true;
    }


    bool
    testCopyConstructor()
    {
      std::cout << "Testing ReferenceCount::ReferenceCount(other)..."
                << std::endl;

      ReferenceCount count0(0);
      ReferenceCount count1(count0);
      if(!checkState(count0, false, false, 0)) {
        return false;
      }
      if(!checkState(count1, false, false, 0)) {
        return false;
      }

      ReferenceCount count2(1);
      ReferenceCount count3(count2);
      if(!checkState(count2, true, true, 2)) {
        return false;
      }
      if(!checkState(count3, true, true, 2)) {
        return false;
      }
      ReferenceCount count4(count2);
      if(!checkState(count2, true, true, 3)) {
        return false;
      }
      if(!checkState(count3, true, true, 3)) {
        return false;
      }
      if(!checkState(count4, true, true, 3)) {
        return false;
      }

      // If we get this far, then all is well.
      return true;
    }


    bool
    testDestructor()
    {
      std::cout << "Testing ReferenceCount::~ReferenceCount()..."
                << std::endl;

      ReferenceCount count0(0);
      {
        ReferenceCount count1(count0);
      }
      if(!checkState(count0, false, false, 0)) {
        return false;
      }

      ReferenceCount count2(1);
      {
        ReferenceCount count3(count2);
        {
          ReferenceCount count4(count2);
          if(!checkState(count2, true, true, 3)) {
            return false;
          }
          if(!checkState(count3, true, true, 3)) {
            return false;
          }
          if(!checkState(count4, true, true, 3)) {
            return false;
          }
        }
        if(!checkState(count2, true, true, 2)) {
          return false;
        }
        if(!checkState(count3, true, true, 2)) {
          return false;
        }
      }
      if(!checkState(count2, true, false, 1)) {
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
  result &= brick::common::testConstructor();
  result &= brick::common::testCopyConstructor();
  result &= brick::common::testDestructor();
  return (result ? 0 : 1);
}
