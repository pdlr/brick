/**
***************************************************************************
* @file utilities3DTest.cpp
*
* Source file defining tests for dlrGeometry library utilities.
*
* Copyright (C) 2007 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/functional.hh>
#include <brick/geometry/utilities3D.hh>
#include <brick/numeric/rotations.hh>
#include <brick/numeric/utilities.hh>
#include <brick/test/testFixture.hh>


namespace brick {

  namespace geometry {
    
    class Utilities3DTest : public TestFixture<Utilities3DTest> {

    public:

      Utilities3DTest();
      ~Utilities3DTest() {}

      void setUp(const std::string& testName) {}
      void tearDown(const std::string& testName) {}

      // Tests.
      void testCheckIntersect__ray__triangle();
      void testFindIntersect__ray__plane__doubleRef();
      void testFindIntersect__ray3D__ray3D__doubleRef();
      void testOperatorTimes__Transform3D__Plane3D();
      void testOperatorTimes__Transform3D__Ray3D();

    private:

      const double m_defaultTolerance;
      
    }; // class Utilities3DTest


    /* ============== Member Function Definititions ============== */

    Utilities3DTest::
    Utilities3DTest()
      : TestFixture<Utilities3DTest>("Utilities3DTest"),
        m_defaultTolerance(1.0E-12)
    {
      BRICK_TEST_REGISTER_MEMBER(testCheckIntersect__ray__triangle);
      BRICK_TEST_REGISTER_MEMBER(testFindIntersect__ray__plane__doubleRef);
      BRICK_TEST_REGISTER_MEMBER(testFindIntersect__ray3D__ray3D__doubleRef);
      BRICK_TEST_REGISTER_MEMBER(testOperatorTimes__Transform3D__Plane3D);
      BRICK_TEST_REGISTER_MEMBER(testOperatorTimes__Transform3D__Ray3D);
    }


    void
    Utilities3DTest::
    testCheckIntersect__ray__triangle()
    {
      std::vector<Triangle3D> triangles;
      triangles.push_back(Triangle3D(Vector3D(-1.0, 3.0, 0.0),
                                     Vector3D(20.0, 3.0, 0.0),
                                     Vector3D(2.0, 0.0, 0.0)));
      triangles.push_back(Triangle3D(Vector3D(20.0, 3.0, 0.0),
                                     Vector3D(-1.0, 3.0, 0.0),
                                     Vector3D(2.0, 0.0, 0.0)));
      triangles.push_back(Triangle3D(Vector3D(2.0, 0.0, 0.0),
                                     Vector3D(20.0, 3.0, 0.0),
                                     Vector3D(-1.0, 3.0, 0.0)));

      std::vector<Vector3D> interiorPoints;
      interiorPoints.push_back(Vector3D(1.1, 1.0, 0.0));
      interiorPoints.push_back(Vector3D(0.0, 2.9, 0.0));
      interiorPoints.push_back(Vector3D(13.9, 2.0, 0.0));

      std::vector<Vector3D> exteriorPoints;
      exteriorPoints.push_back(Vector3D(0.9, 1.0, 0.0));
      exteriorPoints.push_back(Vector3D(0.0, 3.1, 0.0));
      exteriorPoints.push_back(Vector3D(14.1, 2.0, 0.0));

      std::vector<Vector3D> rayOrigins;
      rayOrigins.push_back(Vector3D(1.0, 2.0, 10.0));
      rayOrigins.push_back(Vector3D(1.0, 2.0, -10.0));
      rayOrigins.push_back(Vector3D(-14.0, 2.0, 2.0));
      
      Vector3D intersect;
      for(size_t triangleIndex = 0; triangleIndex < triangles.size();
          ++triangleIndex) {
        for(size_t originIndex = 0; originIndex < rayOrigins.size();
            ++originIndex) {
          for(size_t pointIndex = 0; pointIndex < interiorPoints.size();
              ++pointIndex) {
            Vector3D direction =
              interiorPoints[pointIndex] - rayOrigins[originIndex];
            Ray3D ray(rayOrigins[originIndex], direction);
            BRICK_TEST_ASSERT(
              checkIntersect(ray, triangles[triangleIndex], intersect));
            double residual = magnitude(intersect - interiorPoints[pointIndex]);
            BRICK_TEST_ASSERT(residual < m_defaultTolerance);
          }
          for(size_t pointIndex = 0; pointIndex < exteriorPoints.size();
              ++pointIndex) {
            Vector3D direction =
              exteriorPoints[pointIndex] - rayOrigins[originIndex];
            Ray3D ray(rayOrigins[originIndex], direction);
            BRICK_TEST_ASSERT(
              !checkIntersect(ray, triangles[triangleIndex], intersect));
          }
        }
      }
    }

    
    void
    Utilities3DTest::
    testFindIntersect__ray__plane__doubleRef()
    {
      Plane3D plane0(Vector3D(1.0, 0.0, 1.0),
                     Vector3D(1.0, 1.0, 1.0),
                     Vector3D(0.0, 0.0, 1.0));
      Ray3D ray0(Vector3D(4.0, 2.0, 10.0), Vector3D(0.0, 0.0, -1.0));
      Ray3D ray1(Vector3D(8.0, 33.0, 13.0), Vector3D(-1.0, -5.0, -2.0));

      double distance0;
      double distance1;
      Vector3D point0 = findIntersect(ray0, plane0, distance0);
      Vector3D point1 = findIntersect(ray1, plane0, distance1);

      Vector3D referencePoint0(4.0, 2.0, 1.0);
      Vector3D referencePoint1(2.0, 3.0, 1.0);
      double referenceDistance0 =
        magnitude(ray0.getOrigin() - referencePoint0);
      double referenceDistance1 =
        magnitude(ray1.getOrigin() - referencePoint1);
      
      BRICK_TEST_ASSERT(approximatelyEqual(point0.x(), referencePoint0.x(),
                                         m_defaultTolerance));
      BRICK_TEST_ASSERT(approximatelyEqual(point0.x(), referencePoint0.x(),
                                         m_defaultTolerance));
      BRICK_TEST_ASSERT(approximatelyEqual(point0.x(), referencePoint0.x(),
                                         m_defaultTolerance));
      BRICK_TEST_ASSERT(approximatelyEqual(distance0, referenceDistance0,
                                         m_defaultTolerance));

      BRICK_TEST_ASSERT(approximatelyEqual(point1.x(), referencePoint1.x(),
                                         m_defaultTolerance));
      BRICK_TEST_ASSERT(approximatelyEqual(point1.x(), referencePoint1.x(),
                                         m_defaultTolerance));
      BRICK_TEST_ASSERT(approximatelyEqual(point1.x(), referencePoint1.x(),
                                         m_defaultTolerance));
      BRICK_TEST_ASSERT(approximatelyEqual(distance1, referenceDistance1,
                                         m_defaultTolerance));
    }


    void
    Utilities3DTest::
    testFindIntersect__ray3D__ray3D__doubleRef()
    {
      Ray3D ray0(Vector3D(1.0, 2.0, 10.0), Vector3D(0.0, 0.0, -10.0));
      Ray3D ray1(Vector3D(2.0, 4.0, -5.0), Vector3D(0.0, -2.0, 0.0));

      double residual;
      double distance0;
      double distance1;
      Vector3D point0 = findIntersect(
        ray0, ray1, distance0, distance1, residual);

      Vector3D referencePoint0(1.5, 2.0, -5.0);
      double referenceDistance0 = 15.0;
      double referenceDistance1 = 2.0;
      double referenceResidual(0.5);
      
      BRICK_TEST_ASSERT(approximatelyEqual(point0.x(), referencePoint0.x(),
                                         m_defaultTolerance));
      BRICK_TEST_ASSERT(approximatelyEqual(point0.y(), referencePoint0.y(),
                                         m_defaultTolerance));
      BRICK_TEST_ASSERT(approximatelyEqual(point0.z(), referencePoint0.z(),
                                         m_defaultTolerance));
      BRICK_TEST_ASSERT(approximatelyEqual(distance0, referenceDistance0,
                                         m_defaultTolerance));
      BRICK_TEST_ASSERT(approximatelyEqual(distance1, referenceDistance1,
                                         m_defaultTolerance));
      BRICK_TEST_ASSERT(approximatelyEqual(residual, referenceResidual,
                                         m_defaultTolerance));
    }

    
    void
    Utilities3DTest::
    testOperatorTimes__Transform3D__Plane3D()
    {
      // Arbitrary points and transform.
      Vector3D origin(1.0, 2.0, 3.0);
      Vector3D endPoint0(2.0, -1.0, 3.0);;
      Vector3D endPoint1(-2.0, -3.0, 4.0);;
      Transform3D xf = rollPitchYawToTransform3D(Vector3D(0.2, -0.1, 0.5));
      xf.setValue<0, 3>(2.0);
      xf.setValue<1, 3>(-1.0);
      xf.setValue<2, 3>(-3.0);

      Vector3D newOrigin = xf * origin;
      Vector3D newEndPoint0 = xf * endPoint0;
      Vector3D newEndPoint1 = xf * endPoint1;

      Plane3D plane0(origin, endPoint0, endPoint1);
      Plane3D plane1(newOrigin, newEndPoint0, newEndPoint1);
      Plane3D plane2 = xf * plane0;
      BRICK_TEST_ASSERT(
        magnitude(plane2.getOrigin() - plane1.getOrigin())
        < m_defaultTolerance);
      BRICK_TEST_ASSERT(
        magnitude(plane2.getDirectionVector0() - plane1.getDirectionVector0())
        < m_defaultTolerance);
      BRICK_TEST_ASSERT(
        magnitude(plane2.getDirectionVector1() - plane1.getDirectionVector1())
        < m_defaultTolerance);
    }


    void
    Utilities3DTest::
    testOperatorTimes__Transform3D__Ray3D()
    {
      // Arbitrary points and transform.
      Vector3D startPoint(1.0, 2.0, 3.0);
      Vector3D endPoint(2.0, -1.0, 3.0);;
      Transform3D xf( 1.0, -0.5,  0.2,  3.0,
                     -0.6,  0.7, -0.7,  4.0,
                      0.2,  5.0,  2.0, -6.0,
                      0.1, -0.3,  1.1,  4.0);

      Vector3D newStartPoint = xf * startPoint;
      Vector3D newEndPoint = xf * endPoint;

      Ray3D ray0(startPoint, endPoint - startPoint, false);
      Ray3D ray1(newStartPoint, newEndPoint - newStartPoint, false);
      Ray3D ray2 = xf * ray0;
      BRICK_TEST_ASSERT(
        magnitude(ray2.getOrigin() - ray1.getOrigin())
        < m_defaultTolerance);
      BRICK_TEST_ASSERT(
        magnitude(ray2.getDirectionVector() - ray1.getDirectionVector())
        < m_defaultTolerance);
    }

    
  } // namespace geometry

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::geometry::Utilities3DTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::geometry::Utilities3DTest currentTest;

}

#endif
