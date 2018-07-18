/**
***************************************************************************
* @file brick/optimization/test/autoGradientFunctionLMTest.cc
*
* Source file defining autoGradientFunctionLMTest class.
*
* Copyright (C) 2018 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/optimization/autoGradientFunctionLM.hh>

#include <brick/common/functional.hh>
#include <brick/test/testFixture.hh>

namespace brick {

  namespace optimization {

    class AutoGradientFunctionLMTest : public brick::test::TestFixture<AutoGradientFunctionLMTest> {

    public:

      AutoGradientFunctionLMTest();
      ~AutoGradientFunctionLMTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      void testConstructor();
      void testApplicationOperator();
      void testComputeGradientAndHessian();

    private:

      // Implements [(5*x^2 + 2*y^3), 8*y, -2*x*y].
      struct MySSDFunction {
        template <class InputIter, class OutputIter>
        void apply(InputIter argsBegin, OutputIter resultBegin) {
          typedef typename std::remove_reference<decltype(*resultBegin)>::type
            Scalar;

          Scalar result0 = 5.0 * (*argsBegin) * (*argsBegin);
          Scalar result1 = 0.0;
          Scalar result2 = -2.0 * (*argsBegin);

          ++argsBegin;
          result0 += 2.0 * (*argsBegin) * (*argsBegin) * (*argsBegin);
          result1 += 8.0 * (*argsBegin);
          result2 *= (*argsBegin);

          *resultBegin = result0;
          ++resultBegin;
          *resultBegin = result1;
          ++resultBegin;
          *resultBegin = result2;
        }

        std::size_t getNumberOfArguments(){return 2;}
        std::size_t getNumberOfErrorTerms(){return 3;}
      };

      double m_defaultTolerance;

    }; // class AutoGradientFunctionLMTest


    /* ============== Member Function Definititions ============== */

    AutoGradientFunctionLMTest::
    AutoGradientFunctionLMTest()
      : brick::test::TestFixture<AutoGradientFunctionLMTest>(
        "AutoGradientFunctionLMTest"),
        m_defaultTolerance(1.0E-7)
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testConstructor);
      BRICK_TEST_REGISTER_MEMBER(testApplicationOperator);
      BRICK_TEST_REGISTER_MEMBER(testComputeGradientAndHessian);
    }


    void
    AutoGradientFunctionLMTest::
    testConstructor()
    {
      // Tested implicitly by testApplicationOperator().
    }


    void
    AutoGradientFunctionLMTest::
    testApplicationOperator()
    {
      MySSDFunction ssdFunction;
      AutoGradientFunctionLM<MySSDFunction, 2> gradFunctor(ssdFunction);
      brick::numeric::Array1D<double> args(2);
      args[0] = 0.5;
      args[1] = 2.0;

      double result = gradFunctor(args);

      double term0 = (5.0 * args[0] * args[0]
                      + 2.0 * args[1] * args[1] * args[1]);
      double term1 = 8.0 * args[1];
      double term2 = -2.0 * args[0] * args[1];
      double referenceResult = term0 * term0 + term1 * term1 + term2 * term2;
      BRICK_TEST_ASSERT(
        approximatelyEqual(result, referenceResult, this->m_defaultTolerance));
    }


    void
    AutoGradientFunctionLMTest::
    testComputeGradientAndHessian()
    {
      MySSDFunction ssdFunction;
      AutoGradientFunctionLM<MySSDFunction, 2> gradFunctor(ssdFunction);
      brick::numeric::Array1D<double> args(2);
      args[0] = 0.5;
      args[1] = 2.0;

      brick::numeric::Array1D<double> dEdX;
      brick::numeric::Array2D<double> d2EdX2;

      gradFunctor.computeGradientAndHessian(args, dEdX, d2EdX2);

      double term0 = (5.0 * args[0] * args[0]
                      + 2.0 * args[1] * args[1] * args[1]);
      double term1 = 8.0 * args[1];
      double term2 = -2.0 * args[0] * args[1];

      brick::numeric::Array1D<double> dEdXReference(2);
      dEdXReference[0] = (20.0 * term0 * args[0] - 4.0 * term2 * args[1]);
      dEdXReference[1] = (12.0 * term0 * args[1] * args[1]
                          + 16.0 * term1
                          - 4.0 * term2 * args[0]);

      for(std::size_t ii = 0; ii < 2; ++ii) {
        BRICK_TEST_ASSERT(
          approximatelyEqual(dEdX[ii], dEdXReference[ii],
                             this->m_defaultTolerance));
      }
    }

  } // namespace optimization

} // namespace brick


#if 0

int main(int /* argc */, char** /* argv */)
{
  brick::optimization::AutoGradientFunctionLMTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::optimization::AutoGradientFunctionLMTest currentTest;

}

#endif
