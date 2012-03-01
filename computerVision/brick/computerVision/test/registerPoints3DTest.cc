/**
***************************************************************************
* @file brick/computerVision/test/registerPoints3DTest.cc
*
* Source file defining tests for Horn's method 3D point registration.
*
* Copyright (C) 2006,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/functional.hh>
#include <brick/computerVision/registerPoints3D.hh>
#include <brick/numeric/rotations.hh>
#include <brick/numeric/transform3D.hh>
#include <brick/numeric/vector3D.hh>
#include <brick/test/testFixture.hh>


namespace brick {

  namespace computerVision {
    
    class RegisterPoints3DTest :
      public brick::test::TestFixture<RegisterPoints3DTest> {

    public:

      RegisterPoints3DTest();
      ~RegisterPoints3DTest() {}

      void setUp(const std::string& /* testName */);
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testRegisterPoints3D__In__In__In();
      void testRegisterPoints3D__In__In__In__In();
      void testRegisterPoints3D__In__In__In__Out__double__double__size_t();
    
    private:

      bool
      isApproximatelySameTransform(const Transform3D<double>& transform0,
                                   const Transform3D<double>& transform1);

    
      double m_defaultTolerance;

      size_t m_firstOutlier;
      std::vector<bool> m_flagsVector;
      std::vector<bool> m_flagsVectorReference;
      std::vector< Vector3D<double> > m_fromPoints;
      std::vector< Vector3D<double> > m_fromPointsReference;
      std::vector< Vector3D<double> > m_toPoints;
      std::vector< Vector3D<double> > m_toPointsReference;
      std::vector< Vector3D<double> > m_toPoints2;
      std::vector< Vector3D<double> > m_toPoints2Reference;
    
      Transform3D<double> m_referenceXf;
    
    }; // class RegisterPoints3DTest


    /* ============== Member Function Definititions ============== */

    RegisterPoints3DTest::
    RegisterPoints3DTest()
      : TestFixture<RegisterPoints3DTest>("RegisterPoints3DTest"),
        m_defaultTolerance(1.0E-9),
        m_firstOutlier(2),
        m_flagsVector(),
        m_flagsVectorReference(),
        m_fromPoints(),
        m_fromPointsReference(),
        m_toPoints(),
        m_toPointsReference(),
        m_toPoints2(),
        m_toPoints2Reference(),
        m_referenceXf()    
    {
      BRICK_TEST_REGISTER_MEMBER(testRegisterPoints3D__In__In__In);
      BRICK_TEST_REGISTER_MEMBER(testRegisterPoints3D__In__In__In__In);
      BRICK_TEST_REGISTER_MEMBER(
        testRegisterPoints3D__In__In__In__Out__double__double__size_t);

      // Set up one set of points.
      m_fromPointsReference.push_back(Vector3D<double>(0.0, 0.0, 0.0));
      m_fromPointsReference.push_back(Vector3D<double>(1.0, 0.0, 0.0));
      m_fromPointsReference.push_back(Vector3D<double>(0.5, 0.2, 0.0));
      m_fromPointsReference.push_back(Vector3D<double>(0.0, 3.0, 0.0));
      m_fromPointsReference.push_back(Vector3D<double>(0.0, 0.2, 0.3));
      m_fromPointsReference.push_back(Vector3D<double>(0.0, 0.0, 0.2));
      m_fromPointsReference.push_back(Vector3D<double>(0.1, 0.0, 2.0));
      m_fromPointsReference.push_back(Vector3D<double>(2.0, 0.0, 2.0));

      // Set up a transform.
      Vector3D<double> translation(3.0, -2.5, 1.1);
      Vector3D<double> rotation(0.2, 1.03, -0.15);
      m_referenceXf = rollPitchYawToTransform3D(rotation);
      m_referenceXf.setValue<0, 3>(translation.x());
      m_referenceXf.setValue<1, 3>(translation.y());
      m_referenceXf.setValue<2, 3>(translation.z());

      // Set up 2nd set of points.
      std::transform(m_fromPointsReference.begin(), m_fromPointsReference.end(),
                     std::back_inserter(m_toPointsReference),
                     m_referenceXf.getFunctor());

      // Make some outliers.
      std::copy(m_toPointsReference.begin(), m_toPointsReference.end(),
                std::back_inserter(m_toPoints2Reference));
      m_flagsVectorReference.resize(m_fromPointsReference.size());
      std::fill(m_flagsVectorReference.begin(), m_flagsVectorReference.end(),
                true);
      for(size_t index0 = m_firstOutlier;
          index0 < m_flagsVectorReference.size(); index0 += 2) {
        m_flagsVectorReference[index0] = false;
        m_toPoints2Reference[index0].setValue(1.0, 0.1, 0.1);
      }
    }


    void
    RegisterPoints3DTest::
    setUp(const std::string& /* testName */)
    {
      // Make copies of data so we don't stomp on the original.
      m_flagsVector = m_flagsVectorReference;
      m_fromPoints = m_fromPointsReference;
      m_toPoints = m_toPointsReference;
      m_toPoints2 = m_toPoints2Reference;
    }
  
    
    void
    RegisterPoints3DTest::
    testRegisterPoints3D__In__In__In()
    {
      Transform3D<double> recoveredXf = registerPoints3D<double>(
        m_fromPoints.begin(), m_fromPoints.end(), m_toPoints.begin());
      BRICK_TEST_ASSERT(
        this->isApproximatelySameTransform(recoveredXf, m_referenceXf));
    }


    void
    RegisterPoints3DTest::
    testRegisterPoints3D__In__In__In__In()
    {
      // Make sure we still get the right xf.
      Transform3D<double> recoveredXf = registerPoints3D<double>(
        m_fromPoints.begin(), m_fromPoints.end(), m_toPoints2.begin(),
        m_flagsVector.begin());
      BRICK_TEST_ASSERT(
        this->isApproximatelySameTransform(recoveredXf, m_referenceXf));

      // Make sure an outlier screws us up.
      m_flagsVector[m_firstOutlier] = true;
      recoveredXf = registerPoints3D<double>(
        m_fromPoints.begin(), m_fromPoints.end(), m_toPoints2.begin(),
        m_flagsVector.begin());
      BRICK_TEST_ASSERT(
        !this->isApproximatelySameTransform(recoveredXf, m_referenceXf));
    }

  
    void
    RegisterPoints3DTest::
    testRegisterPoints3D__In__In__In__Out__double__double__size_t()
    {
      // Figure out what percentage of inliers we have, and cheat it up
      // a little to avoid numerical issues.
      size_t numberOfInliers =
        std::count(m_flagsVector.begin(), m_flagsVector.end(), true);
      double inclusion = (numberOfInliers + 0.5) / m_flagsVector.size();

      // Make sure the registration works, and correctly identifies
      // outliers.
      std::vector<bool> flagsVector(m_fromPoints.size());
      Transform3D<double> recoveredXf = registerPoints3D<double>(
        m_fromPoints.begin(), m_fromPoints.end(), m_toPoints2.begin(),
        flagsVector.begin(), inclusion);

      BRICK_TEST_ASSERT(std::equal(flagsVector.begin(), flagsVector.end(),
                                 m_flagsVector.begin()));
      BRICK_TEST_ASSERT(
        this->isApproximatelySameTransform(recoveredXf, m_referenceXf));

      // Now we'll do the test again, only this time we'll specify a
      // threshold for residual values, rather than a percentage of
      // outliers.
      std::vector< Vector3D<double> > transformedPoints(m_fromPoints.size());
      std::transform(m_fromPoints.begin(), m_fromPoints.end(),
                     transformedPoints.begin(), recoveredXf.getFunctor());
      std::vector<double> residuals(m_fromPoints.size());
      for(size_t index0 = 0; index0 < m_fromPoints.size(); ++index0) {
        residuals[index0] =
          magnitude<double>(m_fromPoints[index0] - transformedPoints[index0]);
      }
      std::sort(residuals.begin(), residuals.end());
      double threshold =
        (residuals[numberOfInliers] + residuals[numberOfInliers + 1]) / 2.0;
    
      std::fill(flagsVector.begin(), flagsVector.end(), true);
      recoveredXf = registerPoints3D<double>(
        m_fromPoints.begin(), m_fromPoints.end(), m_toPoints2.begin(),
        flagsVector.begin(), 0.0, threshold);

      BRICK_TEST_ASSERT(std::equal(flagsVector.begin(), flagsVector.end(),
                                 m_flagsVector.begin()));
      BRICK_TEST_ASSERT(
        this->isApproximatelySameTransform(recoveredXf, m_referenceXf));
    }


    bool
    RegisterPoints3DTest::
    isApproximatelySameTransform(const Transform3D<double>& transform0,
                                 const Transform3D<double>& transform1)
    {
      for(size_t rowIndex = 0; rowIndex < 4; ++rowIndex) {
        for(size_t columnIndex = 0; columnIndex < 4; ++columnIndex) {
          if(!approximatelyEqual(transform0(rowIndex, columnIndex),
                                 transform1(rowIndex, columnIndex),
                                 m_defaultTolerance)) {
            return false;
          }
        }
      }
      return true;
    }

  } // namespace computerVision
  
} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::computerVision::RegisterPoints3DTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::RegisterPoints3DTest currentTest;

}

#endif
