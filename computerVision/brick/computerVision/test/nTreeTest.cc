/**
***************************************************************************
* @file brick/computerVision/test/nTreeTest.cc
*
* Source file defining tests for the NTree data structure.
*
* Copyright (C) 2009,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/mathFunctions.hh>
#include <brick/computerVision/nTree.hh>
#include <brick/numeric/index3D.hh>
#include <brick/numeric/vector3D.hh>
#include <brick/numeric/utilities.hh>
#include <brick/test/testFixture.hh>

namespace com = brick::common;
namespace num = brick::numeric;


namespace brick {

  namespace computerVision {

    class NTreeTest
      : public brick::test::TestFixture<NTreeTest> {

    public:

      NTreeTest();
      ~NTreeTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testConstructor();
      void testFind();
      void testFindNearest();

    private:

      num::Vector3D<double> const&
      findNearest(num::Vector3D<double> const& point,
                  std::vector< num::Vector3D<double> > const& candidateVector,
                  double& distance);


      double m_defaultTolerance;

    }; // class NTreeTest


    /* =========== Specialization of KDComparator to work with =========== */
    /* ===========         brick::numeric::Array1D<int>          =========== */

    template <>
    bool
    NComparator< 3, brick::numeric::Array1D<int> >::
    isEqual(brick::numeric::Array1D<int> const& arg0,
            brick::numeric::Array1D<int> const& arg1) const
    {
      if(arg0.size() != arg1.size()) {
        return false;
      }
      for(size_t ii = 0; ii < arg0.size(); ++ii) {
        if(arg0[ii] != arg1[ii]) {
          return false;
        }
      }
      return true;
    }


    /* ============ NTreeTest Member Function Definititions ============ */

    NTreeTest::
    NTreeTest()
      : brick::test::TestFixture<NTreeTest>("NTreeTest"),
        m_defaultTolerance(1.0E-10)
    {
      BRICK_TEST_REGISTER_MEMBER(testConstructor);
      BRICK_TEST_REGISTER_MEMBER(testFind);
      BRICK_TEST_REGISTER_MEMBER(testFindNearest);
    }


    void
    NTreeTest::
    testConstructor()
    {
      std::vector<num::Index3D> inPoints;
      std::vector<num::Index3D> outPoints;

      inPoints.push_back(num::Index3D(2, 3, 1));
      inPoints.push_back(num::Index3D(5, 3, 3));
      inPoints.push_back(num::Index3D(2, 1, 2));
      inPoints.push_back(num::Index3D(7, 4, 2));
      inPoints.push_back(num::Index3D(-2, 6, 5));
      inPoints.push_back(num::Index3D(7, -6, -2));
      inPoints.push_back(num::Index3D(0, 0, 0));
      inPoints.push_back(num::Index3D(3, 4, 5));
      inPoints.push_back(num::Index3D(2, 1, 5));
      inPoints.push_back(num::Index3D(3, 6, 4));
      inPoints.push_back(num::Index3D(-2, 6, 4));

      outPoints.push_back(num::Index3D(1, 3, 2));
      outPoints.push_back(num::Index3D(5, 4, 4));
      outPoints.push_back(num::Index3D(7, 2, 3));
      outPoints.push_back(num::Index3D(2, 2, 3));
      outPoints.push_back(num::Index3D(2, 7, 6));

      NTree<3, num::Index3D> nTree(num::Index3D(-7, -7, -7),
                                   num::Index3D(10, 10, 10));
      for(auto point : inPoints) {
        nTree.add(point);
      }

      for(size_t ii = 0; ii < inPoints.size(); ++ii) {
        BRICK_TEST_ASSERT(nTree.find(inPoints[ii]));
      }
      for(size_t ii = 0; ii < outPoints.size(); ++ii) {
        BRICK_TEST_ASSERT(!(nTree.find(outPoints[ii])));
      }
    }


    void
    NTreeTest::
    testFind()
    {
      std::vector< num::Array1D<int> > inPoints;
      std::vector< num::Array1D<int> > outPoints;

      inPoints.push_back(num::Array1D<int>("[2, 3, 1]"));
      inPoints.push_back(num::Array1D<int>("[5, 3, 3]"));
      inPoints.push_back(num::Array1D<int>("[2, 1, 2]"));
      inPoints.push_back(num::Array1D<int>("[7, 4, 2]"));
      inPoints.push_back(num::Array1D<int>("[-2, 6, 5]"));
      inPoints.push_back(num::Array1D<int>("[7, -6, -2]"));
      inPoints.push_back(num::Array1D<int>("[0, 0, 0]"));
      inPoints.push_back(num::Array1D<int>("[3, 4, 5]"));
      inPoints.push_back(num::Array1D<int>("[2, 1, 5]"));
      inPoints.push_back(num::Array1D<int>("[3, 6, 4]"));
      inPoints.push_back(num::Array1D<int>("[-2, 6, 4]"));

      outPoints.push_back(num::Array1D<int>("[1, 3, 2]"));
      outPoints.push_back(num::Array1D<int>("[5, 4, 4]"));
      outPoints.push_back(num::Array1D<int>("[7, 2, 3]"));
      outPoints.push_back(num::Array1D<int>("[2, 2, 3]"));
      outPoints.push_back(num::Array1D<int>("[2, 7, 6]"));

      NTree< 3, num::Array1D<int> > nTree(num::Array1D<int>("[10, -7, 10]"),
                                          num::Array1D<int>("[-3, 9, 0]"));
      for(auto point : inPoints) {
        nTree.add(point);
      }

      for(size_t ii = 0; ii < inPoints.size(); ++ii) {
        BRICK_TEST_ASSERT(nTree.find(inPoints[ii]));
      }
      for(size_t ii = 0; ii < outPoints.size(); ++ii) {
        BRICK_TEST_ASSERT(!(nTree.find(outPoints[ii])));
      }
    }


    void
    NTreeTest::
    testFindNearest()
    {
      std::vector< num::Vector3D<double> > inPoints;
      std::vector< num::Vector3D<double> > outPoints;

      inPoints.push_back(num::Vector3D<double>(2.1, 3.5, 2.0));
      inPoints.push_back(num::Vector3D<double>(5.0, 3.2, 2.5));
      inPoints.push_back(num::Vector3D<double>(2.4, 1.6, 1.3));
      inPoints.push_back(num::Vector3D<double>(7.7, 4.7, 1.1));
      inPoints.push_back(num::Vector3D<double>(-2.0, 6.3, 5.0));
      inPoints.push_back(num::Vector3D<double>(0.0, 0.0, 0.0));
      inPoints.push_back(num::Vector3D<double>(3.1, 4.7, 5.4));
      inPoints.push_back(num::Vector3D<double>(2.2, 1.6, 5.4));
      inPoints.push_back(num::Vector3D<double>(3.5, 6.9, 4.4));
      inPoints.push_back(num::Vector3D<double>(-2.2, 6.8, 4.3));

      outPoints.push_back(num::Vector3D<double>(1.0, 3.4, 2.0));
      outPoints.push_back(num::Vector3D<double>(5.2, 4.5, 4.5));
      outPoints.push_back(num::Vector3D<double>(7.3, 2.5, 1.1));
      outPoints.push_back(num::Vector3D<double>(2.6, 2.8, 3.0));
      outPoints.push_back(num::Vector3D<double>(2.0, 7.9, 2.5));

      NTree< 3, num::Vector3D<double> > nTree(
        num::Vector3D<double>(-5.0, 8.0, 6.0),
        num::Vector3D<double>(9.0, -1.0, -2.0));
      for(auto point : inPoints) {
        nTree.add(point);
      }

      for(size_t ii = 0; ii < outPoints.size(); ++ii) {
        double distance;
        num::Vector3D<double> nearest = this->findNearest(
          outPoints[ii], inPoints, distance);
        double maybeDistance;
        num::Vector3D<double> maybeNearest = nTree.findNearest(
          outPoints[ii], maybeDistance);
        BRICK_TEST_ASSERT(
          num::magnitude<double>(nearest - maybeNearest) < m_defaultTolerance);
        BRICK_TEST_ASSERT(com::absoluteValue(distance - maybeDistance) < m_defaultTolerance);
      }
    }


    num::Vector3D<double> const&
    NTreeTest::
    findNearest(num::Vector3D<double> const& point,
                std::vector< num::Vector3D<double> > const& candidateVector,
                double& distance)
    {
      num::Vector3D<double> const* nearestPointPtr = &(candidateVector[0]);
      distance = num::magnitudeSquared<double>(point - candidateVector[0]);
      for(size_t ii = 1; ii < candidateVector.size(); ++ii) {
        double newDistance = num::magnitudeSquared<double>(point - candidateVector[ii]);
        if(newDistance < distance) {
          distance = newDistance;
          nearestPointPtr = &candidateVector[ii];
        }
      }
      return *nearestPointPtr;
    }

  } // namespace computerVision

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::computerVision::NTreeTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::NTreeTest currentTest;

}

#endif
