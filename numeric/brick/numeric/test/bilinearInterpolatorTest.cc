/**
***************************************************************************
* @file bilinearInterpolatorTest.cc
* 
* Source file defining BilinearInterpolatorTest class.
*
* Copyright (C) 2004-2005 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/functional.hh>
#include <brick/numeric/bilinearInterpolator.hh>
#include <brick/test/testFixture.hh>


namespace brick {

  namespace numeric {
    
    class BilinearInterpolatorTest :
      public brick::test::TestFixture<BilinearInterpolatorTest> {

    public:

      BilinearInterpolatorTest();
      ~BilinearInterpolatorTest() {};

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests of member functions.
      void testIndexOperator();
    
    private:

     double m_defaultTolerance;
      
    }; // class BilinearInterpolatorTest


    /* ============== Member Function Definititions ============== */

    BilinearInterpolatorTest::
    BilinearInterpolatorTest()
      : brick::test::TestFixture<BilinearInterpolatorTest>("BilinearInterpolatorTest"),
      m_defaultTolerance(1.0E-6)
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testIndexOperator);
    }


    void
    BilinearInterpolatorTest::
    testIndexOperator()
    {
      Array2D<double> myArray(3, 5);
      for(unsigned int rr = 0; rr < myArray.rows(); ++rr) {
        for(unsigned int cc = 0; cc < myArray.columns(); ++cc) {
          myArray(rr, cc) = 2.0 * rr + 3.0 * cc;
        }
      }

      BilinearInterpolator<double> interpolator(myArray);
      for(double yy = 0.0;
          yy < static_cast<double>(myArray.rows() - 1.0);
          yy += 0.1) {

        for(double xx = 0.0;
            xx < static_cast<double>(myArray.columns() - 1.0);
            xx += 0.1) {

          double nominalValue = 2.0 * yy + 3.0 * xx;
          double observedValue = interpolator(yy, xx);

          BRICK_TEST_ASSERT(approximatelyEqual(nominalValue, observedValue,
                                               this->m_defaultTolerance));
        }
      }
    }
    
  } //  namespace numeric

} // namespace brick


#if 1

int main(int /* argc */, char** /* argv */)
{
  brick::numeric::BilinearInterpolatorTest currentTest;
  // bool result = currentTest.run();
  bool result = 0;
  currentTest.testIndexOperator();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::BilinearInterpolatorTest currentTest;

}

#endif
