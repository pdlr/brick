/**
***************************************************************************
* @file brick/optimization/test/lossFunctionsTest.cc
*
* Source file defining a test class for symbols declared in
* brick/optimization/lossFunctions.hh.
*
* Copyright (C) 2018 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/optimization/lossFunctions.hh>

#include <brick/common/functional.hh>
#include <brick/numeric/differentiableScalar.hh>
#include <brick/test/testFixture.hh>

namespace brick {

  namespace optimization {

    class LossFunctionsTest : public brick::test::TestFixture<LossFunctionsTest> {

    public:

      LossFunctionsTest();
      ~LossFunctionsTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      void testLossFunctionCauchy();
      void testLossFunctionHuber();
      void testLossFunctionPseudoHuber();
      void testLossFunctionTukeyBiweight();

    private:


      double m_defaultTolerance;

    }; // class LossFunctionsTest


    /* ============== Member Function Definititions ============== */

    LossFunctionsTest::
    LossFunctionsTest()
      : brick::test::TestFixture<LossFunctionsTest>(
        "LossFunctionsTest"),
        m_defaultTolerance(1.0E-7)
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testLossFunctionCauchy);
      BRICK_TEST_REGISTER_MEMBER(testLossFunctionHuber);
      BRICK_TEST_REGISTER_MEMBER(testLossFunctionPseudoHuber);
      BRICK_TEST_REGISTER_MEMBER(testLossFunctionTukeyBiweight);
    }


    void
    LossFunctionsTest::
    testLossFunctionCauchy()
    {
      double constexpr rootTwo = brick::common::constants::rootTwo;
      double constexpr rootTwoOverTwo =
        brick::common::constants::rootTwoOverTwo;
      double constexpr rangeStart = -10.0;
      double constexpr rangeStop = 10.0;
      double constexpr rangeIncrement = 0.01;


      // Spot check to make sure derivative values are sane.
      {
        LossFunctionCauchy<double> lossFunction;

        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getWeight(0.0), 0.0,
                             this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getWeight(rootTwo), rootTwoOverTwo,
                             this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getWeight(-rootTwo), -rootTwoOverTwo,
                             this->m_defaultTolerance));
      }

      // Check that function value and derivative relate correctly.
      {
        LossFunctionCauchy<brick::numeric::DifferentiableScalar<double>>
          lossFunction;

        for(double zz = rangeStart; zz < rangeStop; zz += rangeIncrement) {
          brick::numeric::DifferentiableScalar<double> zzDiff(zz);
          zzDiff.setDerivative(1.0);

          brick::numeric::DifferentiableScalar<double> loss =
            lossFunction.getValue(zzDiff);
          brick::numeric::DifferentiableScalar<double> weight =
            lossFunction.getWeight(zzDiff);

          BRICK_TEST_ASSERT(
            approximatelyEqual(loss.getDerivative(), weight.getValue(),
                               this->m_defaultTolerance));
        }
      }

      // Check that function value and square root relate correctly.
      {
        LossFunctionCauchy<double> lossFunction;

        for(double zz = rangeStart; zz < rangeStop; zz += rangeIncrement) {
          double loss = lossFunction.getValue(zz);
          double residual = lossFunction.getL2Equivalent(zz);

          BRICK_TEST_ASSERT(
            approximatelyEqual(loss, residual * residual,
                               this->m_defaultTolerance));
        }
      }
    }


    void
    LossFunctionsTest::
    testLossFunctionHuber()
    {
      double constexpr rangeStart = -10.0;
      double constexpr rangeStop = 10.0;
      double constexpr rangeIncrement = 0.01;
      double constexpr delta = 2.0;


      // Spot check to make sure function values are sane.
      {
        LossFunctionHuber<double> lossFunction(delta);

        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getValue(0.0), 0.0,
                             this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getValue(1.0), 0.5,
                             this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getValue(2.0), 2.0,
                             this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getValue(4.0), 6.0,
                             this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getValue(-1.0), 0.5,
                             this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getValue(-2.0), 2.0,
                             this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getValue(-4.0), 6.0,
                             this->m_defaultTolerance));
      }

      // Check that function value and derivative relate correctly.
      {
        LossFunctionHuber<brick::numeric::DifferentiableScalar<double>>
          lossFunction;

        for(double zz = rangeStart; zz < rangeStop; zz += rangeIncrement) {
          brick::numeric::DifferentiableScalar<double> zzDiff(zz);
          zzDiff.setDerivative(1.0);

          brick::numeric::DifferentiableScalar<double> loss =
            lossFunction.getValue(zzDiff);
          brick::numeric::DifferentiableScalar<double> weight =
            lossFunction.getWeight(zzDiff);

          BRICK_TEST_ASSERT(
            approximatelyEqual(loss.getDerivative(), weight.getValue(),
                               this->m_defaultTolerance));
        }
      }

      // Check that function value and square root relate correctly.
      {
        LossFunctionHuber<double> lossFunction;

        for(double zz = rangeStart; zz < rangeStop; zz += rangeIncrement) {
          double loss = lossFunction.getValue(zz);
          double residual = lossFunction.getL2Equivalent(zz);

          BRICK_TEST_ASSERT(
            approximatelyEqual(loss, residual * residual,
                               this->m_defaultTolerance));
        }
      }
    }


    void
    LossFunctionsTest::
    testLossFunctionPseudoHuber()
    {
      double constexpr rangeStart = -10.0;
      double constexpr rangeStop = 10.0;
      double constexpr rangeIncrement = 0.01;
      double constexpr delta = 2.0;

      // Spot check to make sure function values are sane.
      {
        LossFunctionPseudoHuber<double> lossFunction(delta);

        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getValue(0.0), 0.0,
                             this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getValue(1.0), 0.472135955,
                             this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getValue(2.0), 1.656854249,
                             this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getValue(4.0), 4.94427191,
                             this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getValue(-1.0), 0.472135955,
                             this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getValue(-2.0), 1.656854249,
                             this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getValue(-4.0), 4.94427191,
                             this->m_defaultTolerance));
      }

      // Check that function value and derivative relate correctly.
      {
        LossFunctionPseudoHuber<brick::numeric::DifferentiableScalar<double>>
          lossFunction;

        for(double zz = rangeStart; zz < rangeStop; zz += rangeIncrement) {
          brick::numeric::DifferentiableScalar<double> zzDiff(zz);
          zzDiff.setDerivative(1.0);

          brick::numeric::DifferentiableScalar<double> loss =
            lossFunction.getValue(zzDiff);
          brick::numeric::DifferentiableScalar<double> weight =
            lossFunction.getWeight(zzDiff);

          BRICK_TEST_ASSERT(
            approximatelyEqual(loss.getDerivative(), weight.getValue(),
                               this->m_defaultTolerance));
        }
      }

      // Check that function value and square root relate correctly.
      {
        LossFunctionPseudoHuber<double> lossFunction;

        for(double zz = rangeStart; zz < rangeStop; zz += rangeIncrement) {
          double loss = lossFunction.getValue(zz);
          double residual = lossFunction.getL2Equivalent(zz);

          BRICK_TEST_ASSERT(
            approximatelyEqual(loss, residual * residual,
                               this->m_defaultTolerance));
        }
      }
    }

    void
    LossFunctionsTest::
    testLossFunctionTukeyBiweight()
    {
      double constexpr rangeStart = -10.0;
      double constexpr rangeStop = 10.0;
      double constexpr rangeIncrement = 0.01;


      // Spot check to make sure derivative values are sane.
      {
        LossFunctionTukeyBiweight<double> lossFunction;

        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getWeight(0.0), 0.0,
                             this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getWeight(1.0), 0.945216049,
                             this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getWeight(-1.0), -0.945216049,
                             this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getWeight(4.0), 1.234567901,
                             this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getWeight(-4.0), -1.234567901,
                             this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getWeight(6.0), 0.0,
                             this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getWeight(-6.0), -0.0,
                             this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getWeight(8.0), 0.0,
                             this->m_defaultTolerance));
        BRICK_TEST_ASSERT(
          approximatelyEqual(lossFunction.getWeight(-8.0), 0.0,
                             this->m_defaultTolerance));
      }

      // Check that function value and derivative relate correctly.
      {
        LossFunctionTukeyBiweight<brick::numeric::DifferentiableScalar<double>>
          lossFunction;

        for(double zz = rangeStart; zz < rangeStop; zz += rangeIncrement) {
          brick::numeric::DifferentiableScalar<double> zzDiff(zz);
          zzDiff.setDerivative(1.0);

          brick::numeric::DifferentiableScalar<double> loss =
            lossFunction.getValue(zzDiff);
          brick::numeric::DifferentiableScalar<double> weight =
            lossFunction.getWeight(zzDiff);

          BRICK_TEST_ASSERT(
            approximatelyEqual(loss.getDerivative(), weight.getValue(),
                               this->m_defaultTolerance));
        }
      }

      // Check that function value and square root relate correctly.
      {
        LossFunctionTukeyBiweight<double> lossFunction;

        for(double zz = rangeStart; zz < rangeStop; zz += rangeIncrement) {
          double loss = lossFunction.getValue(zz);
          double residual = lossFunction.getL2Equivalent(zz);

          BRICK_TEST_ASSERT(
            approximatelyEqual(loss, residual * residual,
                               this->m_defaultTolerance));
        }
      }
    }

  } // namespace optimization

} // namespace brick


#if 0

int main(int /* argc */, char** /* argv */)
{
  brick::optimization::LossFunctionsTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::optimization::LossFunctionsTest currentTest;

}

#endif
