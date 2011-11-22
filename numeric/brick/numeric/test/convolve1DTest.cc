/**
***************************************************************************
* @file convolve1DTest.cpp
*
* Source file defining Convolve1DTest class.
*
* Copyright (C) 2007 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/functional.hh>
#include <brick/numeric/convolve1D.hh>
#include <brick/test/testFixture.hh>

namespace brick {

  namespace numeric {
  
    template <class Type>
    class Convolve1DTest
      : public brick::test::TestFixture< Convolve1DTest<Type> >
    {

    public:

      // Typedef required why?
      typedef Convolve1DTest<Type> TestFixtureType;
    
    
      Convolve1DTest(const std::string& typeName);
      ~Convolve1DTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      void testConvolve1D_truncateResult();
      void testConvolve1D_padResult();
      void testConvolve1D_padSignal();
      void testConvolve1D_zeroPadSignal();
      void testConvolve1D_reflectSignal();
      void testConvolve1D_wrapSignal();
      void testCorrelate1D_truncateResult();
      void testCorrelate1D_padResult();
      void testCorrelate1D_padSignal();
      void testCorrelate1D_zeroPadSignal();
      void testCorrelate1D_reflectSignal();
      void testCorrelate1D_wrapSignal();

    private:

      template <class Type2>
      bool equivalent(const Array1D<Type2>& arg0,
		      const Array1D<Type2>& arg1,
		      Type2 tolerance);


      Array1D<Type> m_convolve1DKernel;
      Array1D<Type> m_correlate1DKernel;
      Type m_defaultTolerance;
      Type m_fillValue;
      Array1D<Type> m_signal;
      Array1D<Type> m_result_truncateResult;
      Array1D<Type> m_result_padResult;
      Array1D<Type> m_result_padSignal;
      Array1D<Type> m_result_zeroPadSignal;
      Array1D<Type> m_result_reflectSignal;
      Array1D<Type> m_result_wrapSignal;
    
    }; // class Convolve1DTest


    /* ============== Member Function Definititions ============== */

    template <class Type>
    Convolve1DTest<Type>::
    Convolve1DTest(const std::string& typeName)
      : brick::test::TestFixture<Convolve1DTest>(
	std::string("Convolve1DTest<" + typeName + ">")),
	m_convolve1DKernel("[6, 4, 2]"),
	m_correlate1DKernel("[2, 4, 6]"),
	m_defaultTolerance(static_cast<Type>(1.0E-6)),
	m_fillValue(static_cast<Type>(3)),
	m_signal("[1, 2, 3, 4, 5, 6]"),
	m_result_truncateResult("[28, 40, 52, 64]"),
	m_result_padResult("[3, 28, 40, 52, 64, 3]"),
	m_result_padSignal("[24, 22, 28, 40, 52, 64, 52, 42]"),
	m_result_zeroPadSignal("[6, 16, 28, 40, 52, 64, 34, 12]"),
	m_result_reflectSignal("[14, 18, 28, 40, 52, 64, 70, 66]"),
	m_result_wrapSignal("[40, 28, 28, 40, 52, 64, 40, 28]")
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testConvolve1D_truncateResult);
      BRICK_TEST_REGISTER_MEMBER(testConvolve1D_padResult);
      BRICK_TEST_REGISTER_MEMBER(testConvolve1D_padSignal);
      BRICK_TEST_REGISTER_MEMBER(testConvolve1D_zeroPadSignal);
      BRICK_TEST_REGISTER_MEMBER(testConvolve1D_reflectSignal);
      BRICK_TEST_REGISTER_MEMBER(testConvolve1D_wrapSignal);
      BRICK_TEST_REGISTER_MEMBER(testCorrelate1D_truncateResult);
      BRICK_TEST_REGISTER_MEMBER(testCorrelate1D_padResult);
      BRICK_TEST_REGISTER_MEMBER(testCorrelate1D_padSignal);
      BRICK_TEST_REGISTER_MEMBER(testCorrelate1D_zeroPadSignal);
      BRICK_TEST_REGISTER_MEMBER(testCorrelate1D_reflectSignal);
      BRICK_TEST_REGISTER_MEMBER(testCorrelate1D_wrapSignal);
    }


    template <class Type>
    void
    Convolve1DTest<Type>::
    testConvolve1D_truncateResult()
    {
      Array1D<Type> result = convolve1D<Type>(m_convolve1DKernel, m_signal,
					      BRICK_CONVOLVE_TRUNCATE_RESULT,
					      BRICK_CONVOLVE_ROI_VALID);
      BRICK_TEST_ASSERT(
	this->equivalent(result, m_result_truncateResult, m_defaultTolerance));
    }

  
    template <class Type>
    void
    Convolve1DTest<Type>::
    testConvolve1D_padResult()
    {
      Array1D<Type> result = convolve1D<Type>(m_convolve1DKernel, m_signal,
					      BRICK_CONVOLVE_PAD_RESULT,
					      BRICK_CONVOLVE_ROI_SAME,
					      m_fillValue);
      BRICK_TEST_ASSERT(
	this->equivalent(result, m_result_padResult, m_defaultTolerance));
    }


    template <class Type>
    void
    Convolve1DTest<Type>::
    testConvolve1D_padSignal()
    {
      Array1D<Type> result = convolve1D<Type>(m_convolve1DKernel, m_signal,
					      BRICK_CONVOLVE_PAD_SIGNAL,
					      BRICK_CONVOLVE_ROI_FULL,
					      m_fillValue);
      BRICK_TEST_ASSERT(
	this->equivalent(result, m_result_padSignal, m_defaultTolerance));
    }


    template <class Type>
    void
    Convolve1DTest<Type>::
    testConvolve1D_zeroPadSignal()
    {
      Array1D<Type> result = convolve1D<Type>(m_convolve1DKernel, m_signal,
					      BRICK_CONVOLVE_ZERO_PAD_SIGNAL,
					      BRICK_CONVOLVE_ROI_FULL);
      BRICK_TEST_ASSERT(
	this->equivalent(result, m_result_zeroPadSignal, m_defaultTolerance));
    }


    template <class Type>
    void
    Convolve1DTest<Type>::
    testConvolve1D_reflectSignal()
    {
      Array1D<Type> result = convolve1D<Type>(m_convolve1DKernel, m_signal,
					      BRICK_CONVOLVE_REFLECT_SIGNAL,
					      BRICK_CONVOLVE_ROI_FULL);
      BRICK_TEST_ASSERT(
	this->equivalent(result, m_result_reflectSignal, m_defaultTolerance));
    }


    template <class Type>
    void
    Convolve1DTest<Type>::
    testConvolve1D_wrapSignal()
    {
      Array1D<Type> result = convolve1D<Type>(m_convolve1DKernel, m_signal,
					      BRICK_CONVOLVE_WRAP_SIGNAL,
					      BRICK_CONVOLVE_ROI_FULL);
      BRICK_TEST_ASSERT(
	this->equivalent(result, m_result_wrapSignal, m_defaultTolerance));
    }


    template <class Type>
    void
    Convolve1DTest<Type>::
    testCorrelate1D_truncateResult()
    {
      Array1D<Type> result = correlate1D<Type>(m_correlate1DKernel, m_signal,
					       BRICK_CONVOLVE_TRUNCATE_RESULT,
					       BRICK_CONVOLVE_ROI_VALID);
      BRICK_TEST_ASSERT(
	this->equivalent(result, m_result_truncateResult, m_defaultTolerance));
    }

  
    template <class Type>
    void
    Convolve1DTest<Type>::
    testCorrelate1D_padResult()
    {
      Array1D<Type> result = correlate1D<Type>(m_correlate1DKernel, m_signal,
					       BRICK_CONVOLVE_PAD_RESULT,
					       BRICK_CONVOLVE_ROI_SAME,
					       m_fillValue);
      BRICK_TEST_ASSERT(
	this->equivalent(result, m_result_padResult, m_defaultTolerance));
    }


    template <class Type>
    void
    Convolve1DTest<Type>::
    testCorrelate1D_padSignal()
    {
      Array1D<Type> result = correlate1D<Type>(m_correlate1DKernel, m_signal,
					       BRICK_CONVOLVE_PAD_SIGNAL,
					       BRICK_CONVOLVE_ROI_FULL,
					       m_fillValue);
      BRICK_TEST_ASSERT(
	this->equivalent(result, m_result_padSignal, m_defaultTolerance));
    }


    template <class Type>
    void
    Convolve1DTest<Type>::
    testCorrelate1D_zeroPadSignal()
    {
      Array1D<Type> result = correlate1D<Type>(m_correlate1DKernel, m_signal,
					       BRICK_CONVOLVE_ZERO_PAD_SIGNAL,
					       BRICK_CONVOLVE_ROI_FULL);
      BRICK_TEST_ASSERT(
	this->equivalent(result, m_result_zeroPadSignal, m_defaultTolerance));
    }


    template <class Type>
    void
    Convolve1DTest<Type>::
    testCorrelate1D_reflectSignal()
    {
      Array1D<Type> result = correlate1D<Type>(m_correlate1DKernel, m_signal,
					       BRICK_CONVOLVE_REFLECT_SIGNAL,
					       BRICK_CONVOLVE_ROI_FULL);
      BRICK_TEST_ASSERT(
	this->equivalent(result, m_result_reflectSignal, m_defaultTolerance));
    }


    template <class Type>
    void
    Convolve1DTest<Type>::
    testCorrelate1D_wrapSignal()
    {
      Array1D<Type> result = correlate1D<Type>(m_correlate1DKernel, m_signal,
					       BRICK_CONVOLVE_WRAP_SIGNAL,
					       BRICK_CONVOLVE_ROI_FULL);
      BRICK_TEST_ASSERT(
	this->equivalent(result, m_result_wrapSignal, m_defaultTolerance));
    }


    template <class Type>
    template <class Type2>
    bool
    Convolve1DTest<Type>::
    equivalent(const Array1D<Type2>& array0,
	       const Array1D<Type2>& array1,
	       Type2 tolerance)
    {
      if(array0.size() != array1.size()) {
	return false;
      }
      return std::equal(array0.begin(), array0.end(), array1.begin(),
			ApproximatelyEqualFunctor<Type2>(tolerance));
    }

  } //  namespace numeric
  
} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::numeric::Convolve1DTest<double> currentTest0("double");
  brick::numeric::Convolve1DTest<float> currentTest1("float");
  brick::numeric::Convolve1DTest<int> currentTest2("int");
  bool result = (currentTest0.run()
                 && currentTest1.run()
                 && currentTest2.run());
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::Convolve1DTest<double> currentTest0("double");
  brick::numeric::Convolve1DTest<float> currentTest1("float");
  brick::numeric::Convolve1DTest<int> currentTest2("int");

}

#endif

