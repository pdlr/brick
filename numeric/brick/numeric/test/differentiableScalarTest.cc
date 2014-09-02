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
      ~DifferentiableScalarTest() {}

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
      void testGetDerivative();
      void testGetPartialDerivative();
      void testGetValue();

      // Tests of non-member functions.
      void testOperatorTimes();
      void testOperatorDivide();
      void testOperatorPlus();
      void testOperatorMinus();
      void testCosine();
      void testSine();
    
    private:

      void
      getExampleScalars(DifferentiableScalar<double, 2>& twoXSquaredPlusY,
                        DifferentiableScalar<double, 2>& twoXPlusThreeYPlusFive,
                        DifferentiableScalar<double, 2>& fourXCubedEtc);
      
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
      BRICK_TEST_REGISTER_MEMBER(testGetDerivative);
      BRICK_TEST_REGISTER_MEMBER(testGetPartialDerivative);
      BRICK_TEST_REGISTER_MEMBER(testGetValue);
      
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
      DifferentiableScalar<double, 2> twoXSquaredPlusY;
      DifferentiableScalar<double, 2> twoXPlusThreeYPlusFive;
      DifferentiableScalar<double, 2> fourXCubedEtc;
      this->getExampleScalars(twoXSquaredPlusY, twoXPlusThreeYPlusFive,
                              fourXCubedEtc);
      
      twoXSquaredPlusY *= twoXPlusThreeYPlusFive;

      BRICK_TEST_ASSERT(
        approximatelyEqual(twoXSquaredPlusY.getValue(),
                           fourXCubedEtc.getValue(),
                           this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(twoXSquaredPlusY.getPartialDerivative(0),
                           fourXCubedEtc.getPartialDerivative(0),
                           this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(twoXSquaredPlusY.getPartialDerivative(1),
                           fourXCubedEtc.getPartialDerivative(1),
                           this->m_defaultTolerance));
    }

  
    void
    DifferentiableScalarTest::
    testDivideEqualsOperator()
    {
      DifferentiableScalar<double, 2> twoXSquaredPlusY;
      DifferentiableScalar<double, 2> twoXPlusThreeYPlusFive;
      DifferentiableScalar<double, 2> fourXCubedEtc;
      this->getExampleScalars(twoXSquaredPlusY, twoXPlusThreeYPlusFive,
                              fourXCubedEtc);
      
      fourXCubedEtc /= twoXPlusThreeYPlusFive;

      BRICK_TEST_ASSERT(
        approximatelyEqual(twoXSquaredPlusY.getValue(),
                           fourXCubedEtc.getValue(),
                           this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(twoXSquaredPlusY.getPartialDerivative(0),
                           fourXCubedEtc.getPartialDerivative(0),
                           this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(twoXSquaredPlusY.getPartialDerivative(1),
                           fourXCubedEtc.getPartialDerivative(1),
                           this->m_defaultTolerance));
    }

  
    void
    DifferentiableScalarTest::
    testPlusEqualsOperator()
    {
      DifferentiableScalar<double, 2> twoXSquaredPlusY;
      DifferentiableScalar<double, 2> twoXPlusThreeYPlusFive;
      DifferentiableScalar<double, 2> fourXCubedEtc;
      this->getExampleScalars(twoXSquaredPlusY, twoXPlusThreeYPlusFive,
                              fourXCubedEtc);

      DifferentiableScalar<double, 2> arg0;
      DifferentiableScalar<double, 2> arg1;
      DifferentiableScalar<double, 2> arg2;
      this->getExampleScalars(arg0, arg1, arg2);

      arg0 += arg1;

      BRICK_TEST_ASSERT(
        approximatelyEqual(arg0.getValue(),
                           (twoXSquaredPlusY.getValue()
                            + twoXPlusThreeYPlusFive.getValue()),
                           this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(arg0.getPartialDerivative(0),
                           (twoXSquaredPlusY.getPartialDerivative(0)
                            + twoXPlusThreeYPlusFive.getPartialDerivative(0)),
                           this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(arg0.getPartialDerivative(1),
                           (twoXSquaredPlusY.getPartialDerivative(1)
                            + twoXPlusThreeYPlusFive.getPartialDerivative(1)),
                           this->m_defaultTolerance));
    }

  
    void
    DifferentiableScalarTest::
    testMinusEqualsOperator()
    {
      DifferentiableScalar<double, 2> twoXSquaredPlusY;
      DifferentiableScalar<double, 2> twoXPlusThreeYPlusFive;
      DifferentiableScalar<double, 2> fourXCubedEtc;
      this->getExampleScalars(twoXSquaredPlusY, twoXPlusThreeYPlusFive,
                              fourXCubedEtc);

      DifferentiableScalar<double, 2> arg0;
      DifferentiableScalar<double, 2> arg1;
      DifferentiableScalar<double, 2> arg2;
      this->getExampleScalars(arg0, arg1, arg2);

      arg0 -= arg1;

      BRICK_TEST_ASSERT(
        approximatelyEqual(arg0.getValue(),
                           (twoXSquaredPlusY.getValue()
                            - twoXPlusThreeYPlusFive.getValue()),
                           this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(arg0.getPartialDerivative(0),
                           (twoXSquaredPlusY.getPartialDerivative(0)
                            - twoXPlusThreeYPlusFive.getPartialDerivative(0)),
                           this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(arg0.getPartialDerivative(1),
                           (twoXSquaredPlusY.getPartialDerivative(1)
                            - twoXPlusThreeYPlusFive.getPartialDerivative(1)),
                           this->m_defaultTolerance));
    }

  
    void
    DifferentiableScalarTest::
    testGetDerivative()
    {
      double xx = 3.0;
      double yy = 2.0;
      
      DifferentiableScalar<double, 2> twoXSquaredPlusY;
      DifferentiableScalar<double, 2> twoXPlusThreeYPlusFive;
      DifferentiableScalar<double, 2> fourXCubedEtc;
      this->getExampleScalars(twoXSquaredPlusY, twoXPlusThreeYPlusFive,
                              fourXCubedEtc);

      double referenceValue = 2.0 * xx * xx + yy;
      double referencePartial0 = 4.0 * xx;
      double referencePartial1 = 1.0;
      
      BRICK_TEST_ASSERT(
        approximatelyEqual(twoXSquaredPlusY.getValue(),
                           referenceValue,
                           this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(twoXSquaredPlusY.getDerivative(),
                           referencePartial0,
                           this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(twoXSquaredPlusY.getPartialDerivative(0),
                           referencePartial0,
                           this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(twoXSquaredPlusY.getPartialDerivative(1),
                           referencePartial1,
                           this->m_defaultTolerance));
    }

    
    void
    DifferentiableScalarTest::
    testGetPartialDerivative()
    {
      // Tested above.
    }

    
    void
    DifferentiableScalarTest::
    testGetValue()
    {
      // Tested above.
    }
    
    
    void
    DifferentiableScalarTest::
    testOperatorTimes()
    {
      DifferentiableScalar<double, 2> twoXSquaredPlusY;
      DifferentiableScalar<double, 2> twoXPlusThreeYPlusFive;
      DifferentiableScalar<double, 2> fourXCubedEtc;
      this->getExampleScalars(twoXSquaredPlusY, twoXPlusThreeYPlusFive,
                              fourXCubedEtc);

      DifferentiableScalar<double, 2> product0 =
        twoXSquaredPlusY * twoXPlusThreeYPlusFive;

      // Operator*=() has already been tested.
      DifferentiableScalar<double, 2> product1 = twoXSquaredPlusY;
      product1 *= twoXPlusThreeYPlusFive;

      
      BRICK_TEST_ASSERT(
        approximatelyEqual(product0.getValue(),
                           product1.getValue(),
                           this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(product0.getPartialDerivative(0),
                           product1.getPartialDerivative(0),
                           this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(product0.getPartialDerivative(1),
                           product1.getPartialDerivative(1),
                           this->m_defaultTolerance));
    }


    void
    DifferentiableScalarTest::
    testOperatorDivide()
    {
      DifferentiableScalar<double, 2> twoXSquaredPlusY;
      DifferentiableScalar<double, 2> twoXPlusThreeYPlusFive;
      DifferentiableScalar<double, 2> fourXCubedEtc;
      this->getExampleScalars(twoXSquaredPlusY, twoXPlusThreeYPlusFive,
                              fourXCubedEtc);

      DifferentiableScalar<double, 2> quotient0 =
        twoXSquaredPlusY / twoXPlusThreeYPlusFive;

      // Operator/=() has already been tested.
      DifferentiableScalar<double, 2> quotient1 = twoXSquaredPlusY;
      quotient1 /= twoXPlusThreeYPlusFive;

      
      BRICK_TEST_ASSERT(
        approximatelyEqual(quotient0.getValue(),
                           quotient1.getValue(),
                           this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(quotient0.getPartialDerivative(0),
                           quotient1.getPartialDerivative(0),
                           this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(quotient0.getPartialDerivative(1),
                           quotient1.getPartialDerivative(1),
                           this->m_defaultTolerance));
    }


    void
    DifferentiableScalarTest::
    testOperatorPlus()
    {
      DifferentiableScalar<double, 2> twoXSquaredPlusY;
      DifferentiableScalar<double, 2> twoXPlusThreeYPlusFive;
      DifferentiableScalar<double, 2> fourXCubedEtc;
      this->getExampleScalars(twoXSquaredPlusY, twoXPlusThreeYPlusFive,
                              fourXCubedEtc);

      DifferentiableScalar<double, 2> sum0 =
        twoXSquaredPlusY + twoXPlusThreeYPlusFive;

      // Operator+=() has already been tested.
      DifferentiableScalar<double, 2> sum1 = twoXSquaredPlusY;
      sum1 += twoXPlusThreeYPlusFive;

      
      BRICK_TEST_ASSERT(
        approximatelyEqual(sum0.getValue(),
                           sum1.getValue(),
                           this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(sum0.getPartialDerivative(0),
                           sum1.getPartialDerivative(0),
                           this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(sum0.getPartialDerivative(1),
                           sum1.getPartialDerivative(1),
                           this->m_defaultTolerance));
    }
  

    void
    DifferentiableScalarTest::
    testOperatorMinus()
    {
      DifferentiableScalar<double, 2> twoXSquaredPlusY;
      DifferentiableScalar<double, 2> twoXPlusThreeYPlusFive;
      DifferentiableScalar<double, 2> fourXCubedEtc;
      this->getExampleScalars(twoXSquaredPlusY, twoXPlusThreeYPlusFive,
                              fourXCubedEtc);

      DifferentiableScalar<double, 2> difference0 =
        twoXSquaredPlusY - twoXPlusThreeYPlusFive;

      // Operator-=() has already been tested.
      DifferentiableScalar<double, 2> difference1 = twoXSquaredPlusY;
      difference1 -= twoXPlusThreeYPlusFive;

      
      BRICK_TEST_ASSERT(
        approximatelyEqual(difference0.getValue(),
                           difference1.getValue(),
                           this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(difference0.getPartialDerivative(0),
                           difference1.getPartialDerivative(0),
                           this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(difference0.getPartialDerivative(1),
                           difference1.getPartialDerivative(1),
                           this->m_defaultTolerance));
    }

    
    void
    DifferentiableScalarTest::
    testCosine()
    {
      double theta0 = 0.312;
      double theta1 = -0.125;

      // This variable represents 5.0*theta0 + 2.0*theta1.
      DifferentiableScalar<double, 2> fiveTheta0PlusTwoTheta1;
      fiveTheta0PlusTwoTheta1.setValue(5.0 * theta0 + 2.0 * theta1);
      fiveTheta0PlusTwoTheta1.setPartialDerivative(0, 5.0);
      fiveTheta0PlusTwoTheta1.setPartialDerivative(1, 2.0);

      // Do the operation under test. 
      DifferentiableScalar<double, 2> cosStuff = cosine(
        fiveTheta0PlusTwoTheta1);

      // Check the result.
      double omega = 5.0 * theta0 + 2.0 * theta1;
      double referenceValue = brick::common::cosine(omega);
      double referencePartial0 = -5.0 * brick::common::sine(omega);
      double referencePartial1 = -2.0 * brick::common::sine(omega);
      
      BRICK_TEST_ASSERT(
        approximatelyEqual(cosStuff.getValue(),
                           referenceValue,
                           this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(cosStuff.getPartialDerivative(0),
                           referencePartial0,
                           this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(cosStuff.getPartialDerivative(1),
                           referencePartial1,
                           this->m_defaultTolerance));
    }


    void
    DifferentiableScalarTest::
    testSine()
    {
      double theta0 = 0.312;
      double theta1 = -0.125;

      // This variable represents 5.0*theta0 + 2.0*theta1.
      DifferentiableScalar<double, 2> fiveTheta0PlusTwoTheta1;
      fiveTheta0PlusTwoTheta1.setValue(5.0 * theta0 + 2.0 * theta1);
      fiveTheta0PlusTwoTheta1.setPartialDerivative(0, 5.0);
      fiveTheta0PlusTwoTheta1.setPartialDerivative(1, 2.0);

      // Do the operation under test. 
      DifferentiableScalar<double, 2> sinStuff = sine(
        fiveTheta0PlusTwoTheta1);

      // Check the result.
      double omega = 5.0 * theta0 + 2.0 * theta1;
      double referenceValue = brick::common::sine(omega);
      double referencePartial0 = 5.0 * brick::common::cosine(omega);
      double referencePartial1 = 2.0 * brick::common::cosine(omega);
      
      BRICK_TEST_ASSERT(
        approximatelyEqual(sinStuff.getValue(),
                           referenceValue,
                           this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(sinStuff.getPartialDerivative(0),
                           referencePartial0,
                           this->m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(sinStuff.getPartialDerivative(1),
                           referencePartial1,
                           this->m_defaultTolerance));
    }


    void
    DifferentiableScalarTest::
    getExampleScalars(DifferentiableScalar<double, 2>& twoXSquaredPlusY,
                      DifferentiableScalar<double, 2>& twoXPlusThreeYPlusFive,
                      DifferentiableScalar<double, 2>& fourXCubedEtc)
    {
      double const xx = 3.0;
      double const yy = 2.0;

      // Arg0 represents 2x^2 + y,
      double value0 = 2.0 * xx * xx + yy;
      double d0_dx = 4.0 * xx;
      double d0_dy = 1.0;
      twoXSquaredPlusY.setValue(value0);
      twoXSquaredPlusY.setPartialDerivative(0, d0_dx);
      twoXSquaredPlusY.setPartialDerivative(1, d0_dy);

      // Arg1 represents 2x + 3y + 5.
      double value1 = 2.0 * xx + 3.0 * yy + 5.0;
      double d1_dx = 2.0;
      double d1_dy = 3.0;
      twoXPlusThreeYPlusFive.setValue(value1);
      twoXPlusThreeYPlusFive.setPartialDerivative(0, d1_dx);
      twoXPlusThreeYPlusFive.setPartialDerivative(1, d1_dy);

      // Their product should match 4x^3 + 6(x^2)y + 10x^2 + 2xy + 3y^2 + 5y.
      // d/dx should be 12x^2 + 12xy + 20x + 2y.
      // d/dy should be 6x^2 + 2x + 6y + 5.
      double value2 = (4.0 * xx * xx * xx + 6.0 * xx * xx * yy
                       + 10.0 * xx * xx + 2.0 * xx * yy
                       + 3.0 * yy * yy + 5.0 * yy);
      double d2_dx = 12.0 * xx * xx + 12.0 * xx * yy + 20.0 * xx + 2.0 * yy;
      double d2_dy = 6.0 * xx * xx + 2.0 * xx + 6.0 * yy + 5.0;
      fourXCubedEtc.setValue(value2);
      fourXCubedEtc.setPartialDerivative(0, d2_dx);
      fourXCubedEtc.setPartialDerivative(1, d2_dy);
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
