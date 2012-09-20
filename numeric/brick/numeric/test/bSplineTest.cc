/**
***************************************************************************
* @file bSplineTest.cpp
* 
* Source file defining BSplineTest class.
*
* Copyright (C) 2006,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/functional.hh>
#include <brick/numeric/bSpline.hh>
#include <brick/test/testFixture.hh>


namespace brick {

  namespace numeric {

    class BSplineTest : public brick::test::TestFixture<BSplineTest> {

    public:

      BSplineTest();
      ~BSplineTest() {};

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests of member functions.
      void testConstructor__size_t__bool();
      void testConstructor__size_t__bool__2();

      // Regression tests.
      void testNonperiodicSpline();
      
    private:

      double
      computePeriodicCubicSpline(double sValue,
                                 std::vector<double>& controlPoints);


      double
      computePeriodicQuadraticSpline(double sValue,
                                     std::vector<double>& controlPoints);
    
    
      double
      cubicBasisFunction(double x);
    
    
      double
      quadraticBasisFunction(double x);
    
    }; // class BSplineTest


    /* ============== Member Function Definititions ============== */

    BSplineTest::
    BSplineTest()
      : brick::test::TestFixture<BSplineTest>("BSplineTest")
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testConstructor__size_t__bool);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__size_t__bool__2);
      // BRICK_TEST_REGISTER_MEMBER(testNonperiodicSpline);
    }


    void
    BSplineTest::
    testConstructor__size_t__bool()
    {
      const double rangeStart = -10.0;
      const double rangeStop = 10.0;
      const double rangeStep = 4.5;
      const size_t numberOfNodes = 4;

      BSpline<double> bSpline;
      bSpline.setNumberOfNodes(numberOfNodes);
      std::vector<double> controlPoints(3, 0.0);
      for(double point0 = rangeStart; point0 < rangeStop;
          point0 += rangeStep) {
        controlPoints[0] = point0;
        for(double point1 = rangeStart; point1 < rangeStop;
            point1 += rangeStep) {
          controlPoints[1] = point1;
          for(double point2 = rangeStart; point2 < rangeStop;
              point2 += rangeStep) {
            controlPoints[2] = point2;
            bSpline.setControlPoints(controlPoints);

            for(double sValue = 0.0; sValue < double(numberOfNodes - 1);
                sValue += 0.2) {
              double result = bSpline(sValue);
              double referenceValue = this->computePeriodicQuadraticSpline(
                sValue, controlPoints);
              BRICK_TEST_ASSERT(
                approximatelyEqual(result, referenceValue, 1.0E-11));
            }
          }
        }
      }
    }


    void
    BSplineTest::
    testConstructor__size_t__bool__2()
    {
      const double rangeStart = -10.0;
      const double rangeStop = 10.0;
      const double rangeStep = 4.5;
      const size_t numberOfNodes = 5;

      BSpline<double> bSpline(3);
      bSpline.setNumberOfNodes(numberOfNodes);
      std::vector<double> controlPoints(numberOfNodes - 1, 0.0);
      for(double point0 = rangeStart; point0 < rangeStop;
          point0 += rangeStep) {
        controlPoints[0] = point0;
        for(double point1 = rangeStart; point1 < rangeStop;
            point1 += rangeStep) {
          controlPoints[1] = point1;
          for(double point2 = rangeStart; point2 < rangeStop;
              point2 += rangeStep) {
            controlPoints[2] = point2;
            for(double point3 = rangeStart; point2 < rangeStop;
                point2 += rangeStep) {
              controlPoints[3] = point3;
              bSpline.setControlPoints(controlPoints);

              for(double sValue = 0.0; sValue < double(numberOfNodes - 1);
                  sValue += 0.2) {
                double result = bSpline(sValue);
                double referenceValue = this->computePeriodicCubicSpline(
                  sValue, controlPoints);
                BRICK_TEST_ASSERT(
                  approximatelyEqual(result, referenceValue, 1.0E-11));
              }
            }
          }
        }
      }
    }


    void
    BSplineTest::
    testNonperiodicSpline()
    {
      BSpline<double> zoomInterp(3, false);

      // Build spline lookup table.
      {
        BSpline<double>& spline = zoomInterp;

        std::vector<double> controlPoints(36, 0.0);
        std::vector<size_t> knotMultiplicities(36, 0);

        // spline.setNumberOfNodes(36, false);
        spline.setNumberOfNodes(36);

        for(int32_t ii=0; ii < 36; ++ii) {
          if(0 == ii || 35 == ii) {
            knotMultiplicities[ii] = 2;
          } else {
            knotMultiplicities[ii] = 1;
          }
        }              

        // spline.setKnotMultiplicities(knotMultiplicities);

        // Sony block zoom position from FCB-EX1020_TM.pdf, page 55
        controlPoints[0]  = 0.0;
        controlPoints[1]  = 5743.0;
        controlPoints[2]  = 8176.0;
        controlPoints[3]  = 9597.0;
        controlPoints[4]  = 10560.0;
        controlPoints[5]  = 11266.0;
        controlPoints[6]  = 11819.0;
        controlPoints[7]  = 12270.0;
        controlPoints[8]  = 12650.0;
        controlPoints[9]  = 12978.0;
        controlPoints[10] = 13268.0;
        controlPoints[11] = 13529.0;
        controlPoints[12] = 13768.0;
        controlPoints[13] = 13988.0;
        controlPoints[14] = 14195.0;
        controlPoints[15] = 14195.0;
        controlPoints[16] = 14576.0;
        controlPoints[17] = 14752.0;
        controlPoints[18] = 14921.0;
        controlPoints[19] = 15080.0;
        controlPoints[20] = 15231.0;
        controlPoints[21] = 15372.0;
        controlPoints[22] = 15502.0;
        controlPoints[23] = 15622.0;
        controlPoints[24] = 15731.0;
        controlPoints[25] = 15828.0;
        controlPoints[26] = 15916.0;
        controlPoints[27] = 15996.0;
        controlPoints[28] = 16066.0;
        controlPoints[29] = 16128.0;
        controlPoints[30] = 16184.0;
        controlPoints[31] = 16232.0;
        controlPoints[32] = 16276.0;
        controlPoints[33] = 16317.0;
        controlPoints[34] = 16351.0;
        controlPoints[35] = 16384.0;

        spline.setControlPoints(controlPoints);

        //
        // Test

        std::cout << "Zoom 0 : " << spline(0.0) << std::endl;
        std::cout << "Zoom 16: " << spline(16.0) << std::endl;
        std::cout << "Zoom 35: " << spline(34.0) << std::endl;
        std::cout << "Zoom 16.5: " << spline(16.5) << std::endl;
      }
    }
  

    double
    BSplineTest::
    computePeriodicQuadraticSpline(double sValue,
                                   std::vector<double>& controlPoints)
    {
      int spanIndex2 = static_cast<int>(std::floor(sValue));
      int spanIndex1 = spanIndex2 - 1;
      int spanIndex0 = spanIndex1 - 1;

      int maxSpanIndex = static_cast<int>(controlPoints.size()) - 1;
      while(spanIndex2 < 0) {spanIndex2 += static_cast<int>(controlPoints.size());}
      while(spanIndex1 < 0) {spanIndex1 += static_cast<int>(controlPoints.size());}
      while(spanIndex0 < 0) {spanIndex0 += static_cast<int>(controlPoints.size());}
      while(spanIndex2 > maxSpanIndex) {spanIndex2 -= static_cast<int>(controlPoints.size());}
      while(spanIndex1 > maxSpanIndex) {spanIndex1 -= static_cast<int>(controlPoints.size());}
      while(spanIndex0 > maxSpanIndex) {spanIndex0 -= static_cast<int>(controlPoints.size());}

      double result = 0.0;
      double offset = sValue - floor(sValue);
      result += (controlPoints[spanIndex0]
                 * this->quadraticBasisFunction(offset + 2.0));
      result += (controlPoints[spanIndex1]
                 * this->quadraticBasisFunction(offset + 1.0));
      result += (controlPoints[spanIndex2]
                 * this->quadraticBasisFunction(offset));
      return result;
    }

  
    double
    BSplineTest::
    computePeriodicCubicSpline(double sValue,
                               std::vector<double>& controlPoints)
    {
      int spanIndex3 = static_cast<int>(std::floor(sValue));
      int spanIndex2 = spanIndex3 - 1;
      int spanIndex1 = spanIndex2 - 1;
      int spanIndex0 = spanIndex1 - 1;

      int maxSpanIndex = static_cast<int>(controlPoints.size()) - 1;
      while(spanIndex2 < 0) {spanIndex2 += static_cast<int>(controlPoints.size());}
      while(spanIndex1 < 0) {spanIndex1 += static_cast<int>(controlPoints.size());}
      while(spanIndex0 < 0) {spanIndex0 += static_cast<int>(controlPoints.size());}
      while(spanIndex2 > maxSpanIndex) {spanIndex2 -= static_cast<int>(controlPoints.size());}
      while(spanIndex1 > maxSpanIndex) {spanIndex1 -= static_cast<int>(controlPoints.size());}
      while(spanIndex0 > maxSpanIndex) {spanIndex0 -= static_cast<int>(controlPoints.size());}

      double result = 0.0;
      double offset = sValue - floor(sValue);
      result += (controlPoints[spanIndex0]
                 * this->cubicBasisFunction(offset + 3.0));
      result += (controlPoints[spanIndex1]
                 * this->cubicBasisFunction(offset + 2.0));
      result += (controlPoints[spanIndex2]
                 * this->cubicBasisFunction(offset + 1.0));
      result += (controlPoints[spanIndex3]
                 * this->cubicBasisFunction(offset));
      return result;
    }

  
    double
    BSplineTest::
    cubicBasisFunction(double x)
    {
      double returnValue;
      if(x < 0.0) {
        returnValue = 0;
      } else if(x < 1.0) {
        double x3 = x * x * x;
        returnValue =  (1.0 / 6.0) * x3;
      } else if(x < 2.0) {
        double x2 = x * x;
        double x3 = x2 * x;
        returnValue = (-0.5 * x3) + (2.0 * x2) - (2.0 * x) + (2.0 / 3.0);
      } else if(x < 3.0) {
        double x2 = x * x;
        double x3 = x2 * x;
        returnValue = (0.5 * x3) - (4.0 * x2) + (10.0 * x) - (22.0 / 3.0);
      } else if(x < 4.0) {
        double x2 = x * x;
        double x3 = x2 * x;
        returnValue = (((-1.0 / 6.0) * x3) + (2.0 * x2) - (8.0 * x)
                       + (32.0 / 3.0));
      } else {
        returnValue = 0;
      }
      return returnValue;
    }

  
    double
    BSplineTest::
    quadraticBasisFunction(double x)
    {
      double returnValue;
      if(x < 0.0) {
        returnValue = 0;
      } else if(x < 1.0) {
        returnValue =  0.5 * x * x;
      } else if(x < 2.0) {
        returnValue = (-x * x) + (3.0 * x) - (3.0 / 2.0);
      } else if(x < 3.0) {
        returnValue = (0.5 * x * x) - (3.0 * x) + (9.0 / 2.0);
      } else {
        returnValue = 0;
      }
      return returnValue;
    }

  } // namespace numeric

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::numeric::BSplineTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::BSplineTest currentTest;

}

#endif

