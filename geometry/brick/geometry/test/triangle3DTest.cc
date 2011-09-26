/**
***************************************************************************
* @file triangle3DTest.cpp
*
* Source file defining tests for the Triangle3D class.
*
* Copyright (C) 2007 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/functional.hh>
#include <brick/geometry/triangle3D.hh>
#include <brick/numeric/utilities.hh>
#include <brick/test/testFixture.hh>


namespace brick {

  namespace geometry {
    
    class Triangle3DTest : public brick::test::TestFixture<Triangle3DTest> {

    public:

      Triangle3DTest();
      ~Triangle3DTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.

    private:

      const double m_defaultTolerance;
      
    }; // class Triangle3DTest


    /* ============== Member Function Definititions ============== */

    Triangle3DTest::
    Triangle3DTest()
      : brick::test::TestFixture<Triangle3DTest>("Triangle3DTest"),
        m_defaultTolerance(1.0E-12)
    {
      // Empty.
    }

  } // namespace geometry

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::geometry::Triangle3DTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::geometry::Triangle3DTest currentTest;

}

#endif
