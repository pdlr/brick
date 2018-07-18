/**
***************************************************************************
* @file triangle2DTest.cpp
*
* Source file defining tests for the Circle3D class template.
*
* Copyright (C) 2014 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/constants.hh>
#include <brick/common/functional.hh>
#include <brick/common/mathFunctions.hh>
#include <brick/geometry/circle3D.hh>
#include <brick/numeric/utilities.hh>
#include <brick/test/testFixture.hh>


namespace brick {

  namespace geometry {

    class Circle3DTest : public brick::test::TestFixture<Circle3DTest> {

    public:

      Circle3DTest();
      ~Circle3DTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testConstructor();
      void testConstructor_Vector3D_Vector3D_Vector3D();

    private:

      const double m_defaultTolerance;

    }; // class Circle3DTest


    /* ============== Member Function Definititions ============== */

    Circle3DTest::
    Circle3DTest()
      : brick::test::TestFixture<Circle3DTest>("Circle3DTest"),
        m_defaultTolerance(1.0E-12)
    {
      BRICK_TEST_REGISTER_MEMBER(testConstructor);
      BRICK_TEST_REGISTER_MEMBER(testConstructor_Vector3D_Vector3D_Vector3D);
    }


    void
    Circle3DTest::
    testConstructor()
    {
      Circle3D<double> circle;
      numeric::Vector3D<double> origin = circle.getOrigin();
      numeric::Vector3D<double> point0 = circle.getPerimeterPoint(0.0);
      numeric::Vector3D<double> point1 = circle.getPerimeterPoint(
        common::constants::piOverTwo);

      BRICK_TEST_ASSERT(
        approximatelyEqual(origin.getX(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(origin.getY(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(origin.getZ(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(point0.getX(), 1.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(point0.getY(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(point0.getZ(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(point1.getX(), 0.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(point1.getY(), 1.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(point1.getZ(), 0.0, m_defaultTolerance));
    }

    void
    Circle3DTest::
    testConstructor_Vector3D_Vector3D_Vector3D()
    {
      for(brick::common::UInt32 ii = 0; ii < 10; ++ii) {
        for(brick::common::UInt32 jj = 0; jj < 10; ++jj) {
          for(brick::common::UInt32 kk = 0; kk < 10; ++kk) {

            // Three arbitrary vertices.
            numeric::Vector3D<double> origin(
              10.0 * ii + 5.0, 2.0 * jj, -3.0 * kk + 1.0);
            numeric::Vector3D<double> basisVector0(
              -2.0 * jj - 3.0, 1.0 * kk + 7.0, 3.0 * kk + 1.5);
            numeric::Vector3D<double> basisVector1(
              7.0 * kk, -2.0 * ii + 1.0, -1.0 * kk);

            // Make basisVector1 orthogonal to basisVector0, and make
            // it have the same length.
            double magnitude0 = numeric::magnitude<double>(basisVector0);
            double projection =
              numeric::dot<double>(basisVector1, basisVector0 / magnitude0);
            basisVector1 -= (projection / magnitude0) * basisVector0;
            double magnitude1 = numeric::magnitude<double>(basisVector1);
            basisVector1 /= magnitude1;
            basisVector1 *= magnitude0;

            // Create and test a circle.
            Circle3D<double> circle(origin, basisVector0, basisVector1);

            numeric::Vector3D<double> origin1 = circle.getOrigin();
            numeric::Vector3D<double> point0 = circle.getPerimeterPoint(0.0);
            numeric::Vector3D<double> point1 = circle.getPerimeterPoint(
              common::constants::piOverTwo);
            double radius = circle.getRadius();

            BRICK_TEST_ASSERT(
              approximatelyEqual(origin1.getX(), origin.getX(),
                                 m_defaultTolerance));
            BRICK_TEST_ASSERT(
              approximatelyEqual(origin1.getY(), origin.getY(),
                                 m_defaultTolerance));
            BRICK_TEST_ASSERT(
              approximatelyEqual(origin1.getZ(), origin.getZ(),
                                 m_defaultTolerance));
            BRICK_TEST_ASSERT(
              approximatelyEqual(point0.getX(),
                                 origin.getX() + basisVector0.getX(),
                                 m_defaultTolerance));
            BRICK_TEST_ASSERT(
              approximatelyEqual(point0.getY(),
                                 origin.getY() + basisVector0.getY(),
                                 m_defaultTolerance));
            BRICK_TEST_ASSERT(
              approximatelyEqual(point0.getZ(),
                                 origin.getZ() + basisVector0.getZ(),
                                 m_defaultTolerance));
            BRICK_TEST_ASSERT(
              approximatelyEqual(point1.getX(),
                                 origin.getX() + basisVector1.getX(),
                                 m_defaultTolerance));
            BRICK_TEST_ASSERT(
              approximatelyEqual(point1.getY(),
                                 origin.getY() + basisVector1.getY(),
                                 m_defaultTolerance));
            BRICK_TEST_ASSERT(
              approximatelyEqual(point1.getZ(),
                                 origin.getZ() + basisVector1.getZ(),
                                 m_defaultTolerance));

            BRICK_TEST_ASSERT(
              approximatelyEqual(radius, magnitude0, m_defaultTolerance));

          }
        }
      }
    }

  } // namespace geometry

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::geometry::Circle3DTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::geometry::Circle3DTest currentTest;

}

#endif
