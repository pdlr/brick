/**
***************************************************************************
* @file ellipse2DTest.cpp
*
* Source file defining tests for the Ellipse2D class.
*
* Copyright (C) 2007 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/functional.hh>
#include <brick/geometry/ellipse2D.hh>
#include <brick/numeric/utilities.hh>
#include <brick/test/testFixture.hh>


namespace brick {

  namespace geometry {
    
    class Ellipse2DTest : public brick::test::TestFixture<Ellipse2DTest> {

    public:

      Ellipse2DTest();
      ~Ellipse2DTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testEstimate();

    private:

      const double m_defaultTolerance;
      
    }; // class Ellipse2DTest


    /* ============== Member Function Definititions ============== */

    Ellipse2DTest::
    Ellipse2DTest()
      : brick::test::TestFixture<Ellipse2DTest>("Ellipse2DTest"),
        m_defaultTolerance(1.0E-12)
    {
      BRICK_TEST_REGISTER_MEMBER(testEstimate);
    }


    void
    Ellipse2DTest::
    testEstimate()
    {
      brick::numeric::Array1D<double> algebraicParameters(
        "[3.0, -5.0, 7.0, -52.0, 23.0, 10.0]");
      Ellipse2D<double> ellipse2D;
      // ellipse2D.estimate();
    }

  } // namespace geometry

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::geometry::Ellipse2DTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::geometry::Ellipse2DTest currentTest;

}

#endif
