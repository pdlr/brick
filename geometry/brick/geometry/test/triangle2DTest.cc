/**
***************************************************************************
* @file triangle2DTest.cpp
*
* Source file defining tests for the Triangle2D class.
*
* Copyright (C) 2014 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/functional.hh>
#include <brick/common/mathFunctions.hh>
#include <brick/geometry/triangle2D.hh>
#include <brick/numeric/utilities.hh>
#include <brick/test/testFixture.hh>


namespace brick {

  namespace geometry {

    class Triangle2DTest : public brick::test::TestFixture<Triangle2DTest> {

    public:

      Triangle2DTest();
      ~Triangle2DTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testConstructor();
      void testConstructor_Vector2D_Vector2D_Vector2D();
      void testGetArea();

    private:

      const double m_defaultTolerance;

    }; // class Triangle2DTest


    /* ============== Member Function Definititions ============== */

    Triangle2DTest::
    Triangle2DTest()
      : brick::test::TestFixture<Triangle2DTest>("Triangle2DTest"),
        m_defaultTolerance(1.0E-12)
    {
      BRICK_TEST_REGISTER_MEMBER(testConstructor);
      BRICK_TEST_REGISTER_MEMBER(testConstructor_Vector2D_Vector2D_Vector2D);
      BRICK_TEST_REGISTER_MEMBER(testGetArea);
    }


    void
    Triangle2DTest::
    testConstructor()
    {
      Triangle2D<double> triangle;
      numeric::Vector2D<double> vertex0 = triangle.getVertex0();
      numeric::Vector2D<double> vertex1 = triangle.getVertex1();
      numeric::Vector2D<double> vertex2 = triangle.getVertex2();

      BRICK_TEST_ASSERT(
        approximatelyEqual(vertex0.getX(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(vertex0.getY(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(vertex1.getX(), 1.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(vertex1.getY(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(vertex2.getX(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(vertex2.getY(), 1.0, m_defaultTolerance));

    }

    void
    Triangle2DTest::
    testConstructor_Vector2D_Vector2D_Vector2D()
    {
      // Indirectly tested by testGetArea().
    }


    void
    Triangle2DTest::
    testGetArea()
    {
      // For a bunch of different triangles...
      for(brick::common::UInt32 ii = 0; ii < 10; ++ii) {
        for(brick::common::UInt32 jj = 0; jj < 10; ++jj) {
          for(brick::common::UInt32 kk = 0; kk < 10; ++kk) {

            // Three arbitrary vertices.
            numeric::Vector2D<double> vertex0(10.0 * ii + 5.0, 2.0 * jj);
            numeric::Vector2D<double> vertex1(-2.0 * jj - 3.0, 1.0 * kk + 7.0);
            numeric::Vector2D<double> vertex2(7.0 * kk, -2.0 * ii + 1.0);

            // Independent calculation of the area of the triangle
            // defined by the three vertices.  We'll use the "0.5 *
            // base * height" method, so first we need the base and
            // height.
            //
            // ... The base is easy.
            numeric::Vector2D<double> baseVector = vertex1 - vertex0;
            double base = numeric::magnitude<double>(baseVector);

            // ... The height requires a little work.
            numeric::Vector2D<double> heightVector = vertex2 - vertex0;
            numeric::Vector2D<double> baseDirection = baseVector / base;
            double projectionOntoBase =
              numeric::dot<double>(heightVector, baseDirection);
            heightVector -= projectionOntoBase * baseDirection;
            double height = numeric::magnitude<double>(heightVector);

            // ... All done.
            double referenceArea = 0.5 * base * height;

            // Now compute area using Triangle2D.
            Triangle2D<double> triangle(vertex0, vertex1, vertex2);
            double area = triangle.getArea();

            BRICK_TEST_ASSERT(
              approximatelyEqual(area, referenceArea, m_defaultTolerance));
          }
        }
      }
    }

  } // namespace geometry

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::geometry::Triangle2DTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::geometry::Triangle2DTest currentTest;

}

#endif
