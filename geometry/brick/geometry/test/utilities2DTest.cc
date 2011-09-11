/**
***************************************************************************
* @file utilities2DTest.cpp
*
* Source file defining tests for dlrGeometry library utilities.
*
* Copyright (C) 2009 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <cmath>
#include <brick/numeric/array2D.hh>
#include <brick/numeric/utilities.hh>
#include <brick/geometry/utilities2D.hh>
#include <brick/test/testFixture.hh>

namespace num = brick::numeric;

namespace brick {

  namespace geometry {
    
    class Utilities2DTest : public TestFixture<Utilities2DTest> {

    public:

      Utilities2DTest();
      ~Utilities2DTest() {}

      void setUp(const std::string& testName) {}
      void tearDown(const std::string& testName) {}

      // Tests.
      void testCheckIntersect__lineSegment2D__lineSegment2D__vector2DRef();
      void testCheckIntersect__ray2D__lineSegment2D__vector2DRef__doubleRef();
      void testFindClosestPoint__vector2D__ray2D();
      void testOperatorTimes__Transform2D__LineSegment2D();
      void testOperatorTimes__Transform2D__Ray2D();

    private:

      bool isApproximatelyEqual(double val0, double val1,
                                double tolerance = 1.0E-12);
      bool isApproximatelyEqual(numeric::Vector2D const& point0,
                                numeric::Vector2D const& point1,
                                double tolerance = 1.0E-12);
        
      const double m_defaultTolerance;
      
    }; // class Utilities2DTest


    /* ============== Member Function Definititions ============== */

    Utilities2DTest::
    Utilities2DTest()
      : TestFixture<Utilities2DTest>("Utilities2DTest"),
        m_defaultTolerance(1.0E-12)
    {
      BRICK_TEST_REGISTER_MEMBER(
        testCheckIntersect__lineSegment2D__lineSegment2D__vector2DRef);
      BRICK_TEST_REGISTER_MEMBER(
        testCheckIntersect__ray2D__lineSegment2D__vector2DRef__doubleRef);
      BRICK_TEST_REGISTER_MEMBER(testFindClosestPoint__vector2D__ray2D);
      BRICK_TEST_REGISTER_MEMBER(testOperatorTimes__Transform2D__LineSegment2D);
      BRICK_TEST_REGISTER_MEMBER(testOperatorTimes__Transform2D__Ray2D);
    }


    void
    Utilities2DTest::    
    testCheckIntersect__lineSegment2D__lineSegment2D__vector2DRef()
    {
      // Define a few test segments.
      std::vector<LineSegment2D> testSegments0;
      testSegments0.push_back(
        LineSegment2D(Vector2D(1.0, 2.0), Vector2D(7.0, 2.0)));
      testSegments0.push_back(
        LineSegment2D(Vector2D(1.0, 2.0), Vector2D(1.0, 19.0)));
      testSegments0.push_back(
        LineSegment2D(Vector2D(1.0, 2.0), Vector2D(-10.0, 2.0)));

      // Define a few more test segments.
      std::vector<LineSegment2D> testSegments1;
      testSegments1.push_back(
        LineSegment2D(Vector2D(3.0, 5.0), Vector2D(9.0, -1.0)));
      testSegments1.push_back(
        LineSegment2D(Vector2D(-3.0, 5.0), Vector2D(3.0, -1.0)));
      testSegments1.push_back(
        LineSegment2D(Vector2D(-10.0, 15.0), Vector2D(10.0, 15.0)));

      // Define intersection points (by inspection).
      num::Array2D<bool> validityArray(
        testSegments0.size(), testSegments1.size());
      num::Array2D<Vector2D> intersectArray(
        testSegments0.size(), testSegments1.size());

      validityArray(0, 0) = true;
      intersectArray(0, 0) = Vector2D(6.0, 2.0);
      validityArray(0, 1) = false;
      validityArray(0, 2) = false;
      validityArray(1, 0) = false;
      validityArray(1, 1) = false;
      validityArray(1, 2) = true;
      intersectArray(1, 2) = Vector2D(1.0, 15.0);
      validityArray(2, 0) = false;
      validityArray(2, 1) = true;
      intersectArray(2, 1) = Vector2D(0.0, 2.0);
      validityArray(2, 2) = false;

      // Now check that the function under test agrees.
      for(unsigned int ii = 0; ii < testSegments0.size(); ++ii) {
        for(unsigned int jj = 0; jj < testSegments1.size(); ++jj) {
          Vector2D intersect;
          bool isValid = checkIntersect(
            testSegments0[ii], testSegments1[jj], intersect);

          BRICK_TEST_ASSERT(isValid == validityArray(ii, jj));
          if(isValid) {
            BRICK_TEST_ASSERT(
              this->isApproximatelyEqual(intersect, intersectArray(ii, jj)));
          }
        }
      }
    }


    void
    Utilities2DTest::
    testCheckIntersect__ray2D__lineSegment2D__vector2DRef__doubleRef()
    {
      // Define a few test rays.
      std::vector<Ray2D> testRays;
      testRays.push_back(Ray2D(Vector2D(1.0, 2.0), Vector2D(1.0, 0.0)));
      testRays.push_back(Ray2D(Vector2D(1.0, 2.0), Vector2D(0.0, 1.0)));
      testRays.push_back(Ray2D(Vector2D(1.0, 2.0), Vector2D(-1.0, 0.0)));

      // Define a few test segments.
      std::vector<LineSegment2D> testSegments;
      testSegments.push_back(
        LineSegment2D(Vector2D(3.0, 5.0), Vector2D(9.0, -1.0)));
      testSegments.push_back(
        LineSegment2D(Vector2D(-3.0, 5.0), Vector2D(3.0, -1.0)));
      testSegments.push_back(
        LineSegment2D(Vector2D(-10.0, 15.0), Vector2D(10.0, 15.0)));

      // Define intersection points (by inspection).
      num::Array2D<bool> validityArray(
        testRays.size(), testSegments.size());
      num::Array2D<Vector2D> intersectArray(
        testRays.size(), testSegments.size());
      num::Array2D<double> lambdaArray(
        testRays.size(), testSegments.size());

      validityArray(0, 0) = true;
      intersectArray(0, 0) = Vector2D(6.0, 2.0);
      lambdaArray(0, 0) = 5.0;
      validityArray(0, 1) = false;
      validityArray(0, 2) = false;
      validityArray(1, 0) = false;
      validityArray(1, 1) = false;
      validityArray(1, 2) = true;
      intersectArray(1, 2) = Vector2D(1.0, 15.0);
      lambdaArray(1, 2) = 13.0;
      validityArray(2, 0) = false;
      validityArray(2, 1) = true;
      intersectArray(2, 1) = Vector2D(0.0, 2.0);
      lambdaArray(2, 1) = 1.0;
      validityArray(2, 2) = false;

      // Now check that the function under test agrees.
      for(unsigned int ii = 0; ii < testRays.size(); ++ii) {
        for(unsigned int jj = 0; jj < testSegments.size(); ++jj) {
          Vector2D intersect;
          double lambda;
          bool isValid = checkIntersect(
            testRays[ii], testSegments[jj], intersect, lambda);

          BRICK_TEST_ASSERT(isValid == validityArray(ii, jj));
          if(isValid) {
            BRICK_TEST_ASSERT(
              this->isApproximatelyEqual(lambda, lambdaArray(ii, jj)));
            BRICK_TEST_ASSERT(
              this->isApproximatelyEqual(intersect, intersectArray(ii, jj)));
          }
        }
      }
    }


    void
    Utilities2DTest::
    testFindClosestPoint__vector2D__ray2D()
    {
      // Make a 2x2 box centered at 10, 20.  Define it's four corners
      // and the center points of each side.
      double boxSize = 2.0;
      num::Vector2D centerPoint(10.0, 20.0);
      num::Vector2D lowerLeft =
        centerPoint + num::Vector2D(-boxSize / 2.0, -boxSize / 2.0);
      num::Vector2D upperLeft =
        centerPoint + num::Vector2D(-boxSize / 2.0, boxSize / 2.0);
      num::Vector2D upperRight =
        centerPoint + num::Vector2D(boxSize / 2.0, boxSize / 2.0);
      num::Vector2D lowerRight =
        centerPoint + num::Vector2D(boxSize / 2.0, -boxSize / 2.0);
      num::Vector2D leftPoint =
        centerPoint + num::Vector2D(-boxSize / 2.0, 0.0);
      num::Vector2D topPoint =
        centerPoint + num::Vector2D(0.0, boxSize / 2.0);
      num::Vector2D rightPoint =
        centerPoint + num::Vector2D(boxSize / 2.0, 0.0);
      num::Vector2D bottomPoint =
        centerPoint + num::Vector2D(0.0, -boxSize / 2.0);

      // Rays along each edge of the box, pointing clockwise.
      Ray2D leftEdge(lowerLeft, num::Vector2D(0.0, 1.0));
      Ray2D topEdge(upperLeft, num::Vector2D(1.0, 0.0));
      Ray2D rightEdge(upperRight, num::Vector2D(0.0, -1.0));
      Ray2D bottomEdge(lowerRight, num::Vector2D(-1.0, 0.0));
      
      // Now find the closest point on each edge to the center.
      // Should be the previously defined center points.
      num::Vector2D testPoint;
      testPoint = findClosestPoint(centerPoint, leftEdge);
      BRICK_TEST_ASSERT(this->isApproximatelyEqual(testPoint, leftPoint));
      testPoint = findClosestPoint(centerPoint, topEdge);
      BRICK_TEST_ASSERT(this->isApproximatelyEqual(testPoint, topPoint));
      testPoint = findClosestPoint(centerPoint, rightEdge);
      BRICK_TEST_ASSERT(this->isApproximatelyEqual(testPoint, rightPoint));
      testPoint = findClosestPoint(centerPoint, bottomEdge);
      BRICK_TEST_ASSERT(this->isApproximatelyEqual(testPoint, bottomPoint));
    }


    void
    Utilities2DTest::
    testOperatorTimes__Transform2D__LineSegment2D()
    {
      // Arbitrary points and transform.
      Vector2D startPoint(1.0, 2.0);
      Vector2D endPoint(2.0, -1.0);;
      Transform2D xf( 1.0, -0.5, 3.0,
                     -0.6,  0.7, 4.0,
                      0.1, -0.3, 4.0);

      Vector2D newStartPoint = xf * startPoint;
      Vector2D newEndPoint = xf * endPoint;

      LineSegment2D lineSegment0(startPoint, endPoint);
      LineSegment2D lineSegment1(newStartPoint, newEndPoint);
      LineSegment2D lineSegment2 = xf * lineSegment0;
      BRICK_TEST_ASSERT(
        magnitude(lineSegment2.getVertex0() - lineSegment1.getVertex0())
        < m_defaultTolerance);
      BRICK_TEST_ASSERT(
        magnitude(lineSegment2.getVertex1() - lineSegment1.getVertex1())
        < m_defaultTolerance);
    }

        
    void
    Utilities2DTest::
    testOperatorTimes__Transform2D__Ray2D()
    {
      // Arbitrary points and transform.
      Vector2D startPoint(1.0, 2.0);
      Vector2D endPoint(2.0, -1.0);;
      Transform2D xf( 1.0, -0.5, 3.0,
                     -0.6,  0.7, 4.0,
                      0.1, -0.3, 4.0);

      Vector2D newStartPoint = xf * startPoint;
      Vector2D newEndPoint = xf * endPoint;

      Ray2D ray0(startPoint, endPoint - startPoint, false);
      Ray2D ray1(newStartPoint, newEndPoint - newStartPoint, false);
      Ray2D ray2 = xf * ray0;
      BRICK_TEST_ASSERT(
        magnitude(ray2.getOrigin() - ray1.getOrigin())
        < m_defaultTolerance);
      BRICK_TEST_ASSERT(
        magnitude(ray2.getDirectionVector() - ray1.getDirectionVector())
        < m_defaultTolerance);
    }

        
    bool
    Utilities2DTest::
    isApproximatelyEqual(double val0,
                         double val1,
                         double tolerance)
    {
      return (std::fabs(val0 - val1) < tolerance);
    }

    
    bool
    Utilities2DTest::
    isApproximatelyEqual(numeric::Vector2D const& point0,
                         numeric::Vector2D const& point1,
                         double tolerance)
    {
      return ((std::fabs(point0.x() - point1.x()) < tolerance)
              && (std::fabs(point0.y() - point1.y()) < tolerance));
    }
    
  } // namespace geometry

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::geometry::Utilities2DTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::geometry::Utilities2DTest currentTest;

}

#endif
