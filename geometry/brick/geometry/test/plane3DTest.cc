/**
***************************************************************************
* @file plane3DTest.cpp
*
* Source file defining tests for the Plane3D class.
*
* Copyright (C) 2007 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/functional.hh>
#include <brick/geometry/plane3D.hh>
#include <brick/numeric/utilities.hh>
#include <brick/test/testFixture.hh>


namespace brick {

  namespace geometry {
    
    class Plane3DTest : public TestFixture<Plane3DTest> {

    public:

      Plane3DTest();
      ~Plane3DTest() {}

      void setUp(const std::string& testName) {}
      void tearDown(const std::string& testName) {}

      // Tests.
      void testConstructor__iterator__iterator();

    private:

      const double m_defaultTolerance;
      
    }; // class Plane3DTest


    /* ============== Member Function Definititions ============== */

    Plane3DTest::
    Plane3DTest()
      : TestFixture<Plane3DTest>("Plane3DTest"),
        m_defaultTolerance(1.0E-12)
    {
      BRICK_TEST_REGISTER_MEMBER(testConstructor__iterator__iterator);
    }


    void
    Plane3DTest::
    testConstructor__iterator__iterator()
    {
      Array1D<Vector3D> pointArray("[Vector3D(-1.0, -1.0, 1.0),"
                                   " Vector3D(-1.0,  0.0, 1.0),"
                                   " Vector3D(-1.0,  1.0, 1.0),"
                                   " Vector3D( 0.0, -1.0, 1.0),"
                                   " Vector3D( 0.0,  0.0, 1.0),"
                                   " Vector3D( 0.0,  1.0, 1.0),"
                                   " Vector3D( 0.0,  1.0, 2.0)," // Outlier.
                                   " Vector3D( 1.0, -1.0, 1.0),"
                                   " Vector3D( 1.0,  0.0, 1.0),"
                                   " Vector3D( 1.0,  1.0, 1.0)]");

      Plane3D plane0(pointArray.begin(), pointArray.end(), 0.9);
      BRICK_TEST_ASSERT(
        approximatelyEqual(plane0.getOrigin().z(), 1.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(plane0.getDirectionVector0().z(), 0.0,
                           m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(plane0.getDirectionVector1().z(), 0.0,
                           m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(
          dot(plane0.getDirectionVector0(), plane0.getDirectionVector1()),
          0.0, m_defaultTolerance));
    }
  
  } // namespace geometry

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::geometry::Plane3DTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::geometry::Plane3DTest currentTest;

}

#endif
