/**
***************************************************************************
* @file brick/numeric/test/index3DTest.cpp
*
* Source file defining tests for the Index3D data structure.
*
* Copyright (C) 2009,2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <algorithm>
#include <iostream>
#include <brick/numeric/index3D.hh>
#include <brick/test/testFixture.hh>

#include <iostream>
#include <brick/common/exception.hh>

namespace num = brick::numeric;

// Anonymous namespace for local functions and classes.
namespace {

  class cmp {
  public:
    bool operator()(num::Index3D const& arg0, num::Index3D const& arg1) {
      return arg0.getSlice() < arg1.getSlice();
    }
  };

} // namespace


namespace brick {

  namespace numeric {

    class Index3DTest
      : public test::TestFixture<Index3DTest> {

    public:

      Index3DTest();
      ~Index3DTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testIndex3D();

    private:

      double m_defaultTolerance;

    }; // class Index3DTest


    /* ============== Member Function Definititions ============== */

    Index3DTest::
    Index3DTest()
      : brick::test::TestFixture<Index3DTest>("Index3DTest"),
        m_defaultTolerance(1.0E-10)
    {
      BRICK_TEST_REGISTER_MEMBER(testIndex3D);
    }


    void
    Index3DTest::
    testIndex3D()
    {
      // Currently we have only one test for this class, and it's
      // a legacy test that we used to track down a bug in the
      // assignment operator.  It stays around because the cost of
      // an extra test is very small.
      std::vector<num::Index3D> testPoints;
      testPoints.push_back(num::Index3D(2, 3, 1));
      testPoints.push_back(num::Index3D(5, 3, 3));
      testPoints.push_back(num::Index3D(3, 1, 2));
      testPoints.push_back(num::Index3D(7, 4, 2));
      testPoints.push_back(num::Index3D(-2, 6, 5));

      std::vector<num::Index3D> sortedPoints;
      sortedPoints.push_back(num::Index3D(-2, 6, 5));
      sortedPoints.push_back(num::Index3D(2, 3, 1));
      sortedPoints.push_back(num::Index3D(3, 1, 2));
      sortedPoints.push_back(num::Index3D(5, 3, 3));
      sortedPoints.push_back(num::Index3D(7, 4, 2));

      std::sort(testPoints.begin(), testPoints.end(), cmp());

      for(size_t ii = 0; ii < testPoints.size(); ++ii) {
        BRICK_TEST_ASSERT(testPoints[ii] == sortedPoints[ii]);
      }
    }

  } // namespace numeric

} // namespace brick

#if 0

int main(int argc, char** argv)
{
  brick::numeric::Index3DTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::Index3DTest currentTest;

}

#endif
