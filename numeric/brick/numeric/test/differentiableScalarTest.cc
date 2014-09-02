/**
***************************************************************************
* @file differentiableScalarTest.cc
* 
* Source file defining DifferentiableScalarTest class.
*
* Copyright (C) 2014 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/functional.hh>
#include <brick/numeric/array1D.hh>
#include <brick/numeric/differentiableScalar.hh>
#include <brick/test/testFixture.hh>


namespace brick {

  namespace numeric {
    
    class DifferentiableScalarTest
      : public brick::test::TestFixture<DifferentiableScalarTest> {

    public:

      DifferentiableScalarTest();
      ~DifferentiableScalarTest() {};

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests of member functions.
      void testConstructor__void();
      void testConstructor__Type();
      void testConstructor__Type__Iter();
      void testConstructor__DifferentiableScalar();
      void testAssignmentOperator();
      void testTimesEqualsOperator();
      void testDivideEqualsOperator();
      void testPlusEqualsOperator();
      void testMinusEqualsOperator();

      // Tests of non-member functions.
      void testOperatorTimes();
      void testOperatorDivide();
      void testOperatorPlus();
      void testOperatorMinus();
      void testCosine();
      void testSine();
    
    private:

      double m_defaultTolerance;

    }; // class DifferentiableScalarTest


    /* ============== Member Function Definititions ============== */

    DifferentiableScalarTest::
    DifferentiableScalarTest()
      : brick::test::TestFixture<DifferentiableScalarTest>("DifferentiableScalarTest"),
        m_defaultTolerance(1.0E-11)
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testConstructor__void);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__Type);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__Type__Iter);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__DifferentiableScalar);
      BRICK_TEST_REGISTER_MEMBER(testAssignmentOperator);
      BRICK_TEST_REGISTER_MEMBER(testTimesEqualsOperator);
      BRICK_TEST_REGISTER_MEMBER(testDivideEqualsOperator);
      BRICK_TEST_REGISTER_MEMBER(testPlusEqualsOperator);
      BRICK_TEST_REGISTER_MEMBER(testMinusEqualsOperator);
      BRICK_TEST_REGISTER_MEMBER(testOperatorTimes);
      BRICK_TEST_REGISTER_MEMBER(testOperatorDivide);
      BRICK_TEST_REGISTER_MEMBER(testOperatorPlus);
      BRICK_TEST_REGISTER_MEMBER(testOperatorMinus);
      BRICK_TEST_REGISTER_MEMBER(testCosine);
      BRICK_TEST_REGISTER_MEMBER(testSine);
    }


    void
    DifferentiableScalarTest::
    testConstructor__void()
    {
      {
        DifferentiableScalar<double> ds;
        BRICK_TEST_ASSERT(
          approximatelyEqual(ds.getValue(), 0.0, this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(ds.getDerivative(), 1.0, this->m_defaultTolerance));
      }

      {
        DifferentiableScalar<double, 4> ds;
        BRICK_TEST_ASSERT(
          approximatelyEqual(ds.getValue(), 0.0, this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(ds.getPartialDerivative(0), 1.0,
                             this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(ds.getPartialDerivative(1), 0.0,
                             this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(ds.getPartialDerivative(2), 0.0,
                             this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(ds.getPartialDerivative(3), 0.0,
                             this->m_defaultTolerance));
      }
    }

  
    void
    DifferentiableScalarTest::
    testConstructor__Type()
    {
      for(double xx = -3.0; xx < 3.0; xx += 0.124) {

        {
          DifferentiableScalar<double> ds(xx);
          BRICK_TEST_ASSERT(
            approximatelyEqual(ds.getValue(), xx, this->m_defaultTolerance));
          BRICK_TEST_ASSERT(
            approximatelyEqual(ds.getDerivative(), 1.0,
                               this->m_defaultTolerance));
        }

        {
          DifferentiableScalar<double, 4> ds;
          BRICK_TEST_ASSERT(
            approximatelyEqual(ds.getValue(), 0.0, this->m_defaultTolerance));
          BRICK_TEST_ASSERT(
            approximatelyEqual(ds.getPartialDerivative(0), 1.0,
                               this->m_defaultTolerance));
          BRICK_TEST_ASSERT(
            approximatelyEqual(ds.getPartialDerivative(1), 0.0,
                               this->m_defaultTolerance));
          BRICK_TEST_ASSERT(
            approximatelyEqual(ds.getPartialDerivative(2), 0.0,
                               this->m_defaultTolerance));
          BRICK_TEST_ASSERT(
            approximatelyEqual(ds.getPartialDerivative(3), 0.0,
                               this->m_defaultTolerance));
        }
        
      }
    }

  
    void
    DifferentiableScalarTest::
    testConstructor__Type__Iter()
    {
      Array1D<double> partials(2);
      for(double xx = -3.0; xx < 3.0; xx += 0.4) {
        for(double partial0 = -3.0; partial0 < 3.0; partial0 += 0.2) {
          for(double partial1 = -3.0; partial1 < 3.0; partial1 += 0.2) {
            partials[0] = partial0;
            partials[1] = partial1;
            DifferentiableScalar<double, 2> ds(xx, partials.begin());
            BRICK_TEST_ASSERT(
              approximatelyEqual(ds.getValue(), xx, this->m_defaultTolerance));
            BRICK_TEST_ASSERT(
              approximatelyEqual(ds.getPartialDerivative(0),
                                 partials[0], this->m_defaultTolerance));
            BRICK_TEST_ASSERT(
              approximatelyEqual(ds.getPartialDerivative(1),
                                 partials[1], this->m_defaultTolerance));
          }
        }
      }
    }

  
    void
    DifferentiableScalarTest::
    testConstructor__DifferentiableScalar()
    {
      Array1D<double> partials(2);
      for(double xx = -3.0; xx < 3.0; xx += 0.4) {
        for(double partial0 = -3.0; partial0 < 3.0; partial0 += 0.2) {
          for(double partial1 = -3.0; partial1 < 3.0; partial1 += 0.2) {
            partials[0] = partial0;
            partials[1] = partial1;
            DifferentiableScalar<double, 2> other(xx, partials.begin());
            DifferentiableScalar<double, 2> ds(other);
            BRICK_TEST_ASSERT(
              approximatelyEqual(ds.getValue(), xx, this->m_defaultTolerance));
            BRICK_TEST_ASSERT(
              approximatelyEqual(ds.getPartialDerivative(0),
                                 partials[0], this->m_defaultTolerance));
            BRICK_TEST_ASSERT(
              approximatelyEqual(ds.getPartialDerivative(1),
                                 partials[1], this->m_defaultTolerance));
          }
        }
      }
    }
  

    void
    DifferentiableScalarTest::
    testAssignmentOperator()
    {
      Array1D<double> partials(2);
      DifferentiableScalar<double, 2> ds;
      for(double xx = -3.0; xx < 3.0; xx += 0.4) {
        for(double partial0 = -3.0; partial0 < 3.0; partial0 += 0.2) {
          for(double partial1 = -3.0; partial1 < 3.0; partial1 += 0.2) {
            partials[0] = partial0;
            partials[1] = partial1;
            DifferentiableScalar<double, 2> other(xx, partials.begin());
            ds = other;
            BRICK_TEST_ASSERT(
              approximatelyEqual(ds.getValue(), xx, this->m_defaultTolerance));
            BRICK_TEST_ASSERT(
              approximatelyEqual(ds.getPartialDerivative(0),
                                 partials[0], this->m_defaultTolerance));
            BRICK_TEST_ASSERT(
              approximatelyEqual(ds.getPartialDerivative(1),
                                 partials[1], this->m_defaultTolerance));
          }
        }
      }
    }

  
    void
    DifferentiableScalarTest::
    testTimesEqualsOperator()
    {
      double xx = 3.0;
      double yy = 2.0;
      
      // Arg0 represents 2x^2 + y,
      double value0 = 2.0 * xx * xx + yy;
      double d0_dx = 4.0 * xx;
      double d0_dy = 1.0;
      DifferentiableScalar<double, 2> arg0;
      arg0.setValue(value0);
      arg0.setPartialDerivative(0, d0_dx);
      arg0.setPartialDerivative(1, d0_dy);

      // Arg0 represents 2x + 3y + 5.
      double value1 = 2.0 * xx + 3.0 * yy + 5.0;
      double d1_dx = 2.0;
      double d1_dy = 3.0;
      DifferentiableScalar<double, 2> arg1;
      arg1.setValue(value1);
      arg1.setPartialDerivative(0, d1_dx);
      arg1.setPartialDerivative(1, d1_dy);

      // Their product should match 4x^3 + 6(x^2)y + 10x^2 + 2xy + 3y^2 + 5y.
      // d/dx should be 12x^2 + 12xy + 20x + 2y.
      // d/dy should be 6x^2 + 2x + 6y + 5.
      double referenceValue = (4.0 * xx * xx * xx + 6.0 * xx * xx * yy
                               + 10.0 * xx * xx + 2.0 * xx * yy
                               + 3.0 * yy * yy + 5.0 * yy);
      double partial0 = 12.0 * xx * xx + 12.0 * xx * yy + 20.0 * xx + 2.0 * yy;
      double partial1 = 6.0 * xx * xx + 2.0 * xx + 6.0 * yy + 5.0;

      arg0 *= arg1;

      BRICK_TEST_ASSERT(
        approximatelyEqual(arg0.getValue(), referenceValue,
                           this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(arg0.getPartialDerivative(0),
                           partial0, this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(arg0.getPartialDerivative(1),
                           partial1, this->m_defaultTolerance));
    }

  
    void
    DifferentiableScalarTest::
    testDivideEqualsOperator()
    {
      double xx = 3.0;
      double yy = 2.0;

      // Arg0 represents 4x^3 + 6(x^2)y + 10x^2 + 2xy + 3y^2 + 5y.
      // d/dx should be 12x^2 + 12xy + 20x + 2y.
      // d/dy should be 6x^2 + 2x + 6y + 5.
      double value0 = (4.0 * xx * xx * xx + 6.0 * xx * xx * yy
                       + 10.0 * xx * xx + 2.0 * xx * yy
                       + 3.0 * yy * yy + 5.0 * yy);
      double d0_dx = 12.0 * xx * xx + 12.0 * xx * yy + 20.0 * xx + 2.0 * yy;
      double d0_dy = 6.0 * xx * xx + 2.0 * xx + 6.0 * yy + 5.0;
      DifferentiableScalar<double, 2> arg0;
      arg0.setValue(value0);
      arg0.setPartialDerivative(0, d0_dx);
      arg0.setPartialDerivative(1, d0_dy);

      // Arg0 represents 2x + 3y + 5.
      double value1 = 2.0 * xx + 3.0 * yy + 5.0;
      double d1_dx = 2.0;
      double d1_dy = 3.0;
      DifferentiableScalar<double, 2> arg1;
      arg1.setValue(value1);
      arg1.setPartialDerivative(0, d1_dx);
      arg1.setPartialDerivative(1, d1_dy);

      // Their quotient should match 2x^2 + y,
      double referenceValue = 2.0 * xx * xx + yy;
      double partial0 = 4.0 * xx;
      double partial1 = 1.0;

      arg0 /= arg1;
      
      BRICK_TEST_ASSERT(
        approximatelyEqual(arg0.getValue(), referenceValue,
                           this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(arg0.getPartialDerivative(0),
                           partial0, this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(arg0.getPartialDerivative(1),
                           partial1, this->m_defaultTolerance));
    }

  
    void
    DifferentiableScalarTest::
    testPlusEqualsOperator()
    {
      // Not implemented yet.
    }

  
    void
    DifferentiableScalarTest::
    testMinusEqualsOperator()
    {
      // Not implemented yet.
    }

  
    void
    DifferentiableScalarTest::
    testOperatorTimes()
    {
      // Not implemented yet.
    }


    void
    DifferentiableScalarTest::
    testOperatorDivide()
    {
      // Not implemented yet.
    }


    void
    DifferentiableScalarTest::
    testOperatorPlus()
    {
      // Not implemented yet.
    }
  

    void
    DifferentiableScalarTest::
    testOperatorMinus()
    {
      // Not implemented yet.
    }

    
    void
    DifferentiableScalarTest::
    testCosine()
    {
      // Not implemented yet.
    }


    void
    DifferentiableScalarTest::
    testSine()
    {
      // Not implemented yet.
    }
    
  } //  namespace numeric

} // namespace brick


#if 0

int main(int /* argc */, char** /* argv */)
{
  brick::numeric::DifferentiableScalarTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::DifferentiableScalarTest currentTest;

}

#endif
