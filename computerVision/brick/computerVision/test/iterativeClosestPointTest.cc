/**
***************************************************************************
* @file brick/computerVision/test/iterativeClosestPointTest.cc
*
* Source file defining tests for the IterativeClosestPoint class template.
*
* Copyright (C) 2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <vector>

#include <brick/computerVision/iterativeClosestPoint.hh>
#include <brick/numeric/rotations.hh>
#include <brick/numeric/vector3D.hh>
#include <brick/test/testFixture.hh>

namespace com = brick::common;
namespace num = brick::numeric;


namespace brick {

  namespace computerVision {

    class IterativeClosestPointTest
      : public brick::test::TestFixture<IterativeClosestPointTest> {

    public:

      IterativeClosestPointTest();
      ~IterativeClosestPointTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testGetTransform();
    
    private:

      double m_defaultTolerance;
    
    }; // class IterativeClosestPointTest


    
    /* ====== IterativeClosestPointTest Member Function Definititions ====== */

    IterativeClosestPointTest::
    IterativeClosestPointTest()
      : brick::test::TestFixture<IterativeClosestPointTest>("IterativeClosestPointTest"),
        m_defaultTolerance(1.0E-5)
    {
      BRICK_TEST_REGISTER_MEMBER(testGetTransform);
    }


    void
    IterativeClosestPointTest::
    testGetTransform()
    {
      // Test parameters.
      unsigned int const extent       = 10;
      // double       const maxRotation  = 0.5;  // Radians.
      double       const maxRotation  = 0.1;  // Radians.
      unsigned int const numRotations = 3;
      // double       const maxTranslation  = 0.5;
      double       const maxTranslation  = 0.2;
      unsigned int const numTranslations = 2;

      std::vector< num::Vector3D<double> > modelPoints;

      // Set up a well behaved shape to register.
      for(unsigned int ii = 0; ii < extent; ++ii) {
        for(unsigned int jj = 0; jj < extent; ++jj) {
          double xValue(ii);
          double yValue(jj);
          double zValue = (std::sin(4.0 * xValue / extent)
                           * std::cos(5.0 * yValue / extent));
          modelPoints.push_back(num::Vector3D<double>(xValue, yValue, zValue));
        }
      }

      // Iterate through a bunch of different ground truth rotations
      // and translations.
      for(double roll = -maxRotation; roll <= maxRotation; 
          roll += (2 * maxRotation / numRotations)) {
        for(double pitch = -maxRotation; pitch <= maxRotation; 
            pitch += (2 * maxRotation / numRotations)) {
          for(double yaw = -maxRotation; yaw <= maxRotation; 
              yaw += (2 * maxRotation / numRotations)) {
            for(double tx = -maxTranslation; tx <= maxTranslation; 
                tx += (2 * maxTranslation / numTranslations)) {
              for(double ty = -maxTranslation; ty <= maxTranslation; 
                  ty += (2 * maxTranslation / numTranslations)) {
                for(double tz = -maxTranslation; tz <= maxTranslation; 
                    tz += (2 * maxTranslation / numTranslations)) {

                  std::cout << "." << std::flush;
                  
                  // Build a ground truth coordinate transform.
                  num::Vector3D<double> rollPitchYaw(roll, pitch, yaw);
                  num::Transform3D<double> observedFromModel =
                    num::rollPitchYawToTransform3D<double>(rollPitchYaw);
                  observedFromModel.setValue(0, 3, tx);
                  observedFromModel.setValue(1, 3, ty);
                  observedFromModel.setValue(2, 3, tz);
                  num::Transform3D<double> modelFromObserved = 
                    observedFromModel.invert();

                  // Use the ground truth transform to compute some
                  // hypothetical observations.
                  std::vector< num::Vector3D<double> > observedPoints(
                    modelPoints.size());
                  std::transform(modelPoints.begin(), modelPoints.end(),
                                 observedPoints.begin(), 
                                 observedFromModel.getFunctor());

                  // Recover the coordinate transformation by ICP.
                  num::Transform3D<double> modelFromObservedEstimate;
                  num::Transform3D<double> observedFromModelEstimate;
                  IterativeClosestPoint<3, num::Vector3D<double>, double> icp;
                  icp.setModelPoints(modelPoints.begin(), modelPoints.end());
                  modelFromObservedEstimate = icp.registerPoints(
                    observedPoints.begin(), observedPoints.end(),
                    modelFromObservedEstimate);
                  observedFromModelEstimate = 
                    modelFromObservedEstimate.invert();

                  // Verify that we've recovered the correct transform.
                  double txEstimate = observedFromModelEstimate(0, 3);
                  double tyEstimate = observedFromModelEstimate(1, 3);
                  double tzEstimate = observedFromModelEstimate(2, 3);
                  num::Vector3D<double> rollPitchYawEstimate =
                    num::transform3DToRollPitchYaw<double>(
                      observedFromModelEstimate);

                  BRICK_TEST_ASSERT(
                    com::approximatelyEqual(
                      txEstimate, tx, this->m_defaultTolerance));
                  BRICK_TEST_ASSERT(
                    com::approximatelyEqual(
                      tyEstimate, ty, this->m_defaultTolerance));
                  BRICK_TEST_ASSERT(
                    com::approximatelyEqual(
                      tzEstimate, tz, this->m_defaultTolerance));
                  BRICK_TEST_ASSERT(
                    com::approximatelyEqual(rollPitchYawEstimate.x(), roll, 
                                            this->m_defaultTolerance));
                  BRICK_TEST_ASSERT(
                    com::approximatelyEqual(rollPitchYawEstimate.y(), pitch, 
                                            this->m_defaultTolerance));
                  BRICK_TEST_ASSERT(
                    com::approximatelyEqual(rollPitchYawEstimate.z(), yaw,
                                            this->m_defaultTolerance));
                }
              }
            }
          }
        }
      }
    }
  
  } // namespace computerVision
  
} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::computerVision::IterativeClosestPointTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::IterativeClosestPointTest currentTest;
  
}

#endif
