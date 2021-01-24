/**
***************************************************************************
* @file brick/numeric/test/convolveNDTest.cpp
*
* Source file defining ConvolveNDTest class.
*
* Copyright (C) 2007-2008 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/functional.hh>
#include <brick/numeric/array2D.hh>
#include <brick/numeric/convolve2D.hh>
#include <brick/numeric/convolveND.hh>
#include <brick/test/testFixture.hh>

namespace brick {

  namespace numeric {

    template <class Type>
    class ConvolveNDTest
      : public brick::test::TestFixture< ConvolveNDTest<Type> >
    {

    public:

      // Typedef required why?
      typedef ConvolveNDTest<Type> TestFixtureType;


      ConvolveNDTest(const std::string& typeName);
      ~ConvolveNDTest() {}

      void setUp(const std::string&) {}
      void tearDown(const std::string&) {}

      void testConvolveND_padResult();

    private:

      template <class Type2, size_t Dimension>
      bool equivalent(const ArrayND<Dimension, Type2>& arg0,
		      const ArrayND<Dimension, Type2>& arg1,
		      Type2 tolerance);


      Array1D<Type> m_convolveNDKernel;
      Type m_defaultTolerance;
      Type m_fillValue;
      ArrayND<2, Type> m_signal;
      Array2D<Type> m_signal2D;

    }; // class ConvolveNDTest


    /* ============== Member Function Definititions ============== */

    template <class Type>
    ConvolveNDTest<Type>::
    ConvolveNDTest(const std::string& typeName)
      : brick::test::TestFixture<ConvolveNDTest>(
	std::string("ConvolveNDTest<" + typeName + ">")),
	m_convolveNDKernel("[6, 4, 2, 0, 5]"),
	m_defaultTolerance(static_cast<Type>(1.0E-6)),
	m_fillValue(static_cast<Type>(3)),
	m_signal(),
        m_signal2D("[[1, 2, 3, 4, 5, 6, -1, 11],"
                   " [2, 3, 4, 5, 6, 7,  2,  2],"
                   " [0, 1, 2, 3, 4, 5, -5,  6],"
                   " [3, 2, 1, 1, 2, 3, 10, -3],"
                   " [4, 4, 3, 1, 2, 0,  4,  4],"
                   " [6, 5, 4, 3, 2, 1, -2, -8],"
                   " [1, 1, 5, 7, 0, 1,  0,  9],"
                   " [1, 2, 3, 4, 3, 2,  1,  2],"
                   " [1, 2, 4, 5, 7, 8,  2,  5]]")
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testConvolveND_padResult);

      // Initialize members.
      m_signal.reinit(m_signal2D.shape());
      m_signal.copy(m_signal2D.data());
    }


    template <class Type>
    void
    ConvolveNDTest<Type>::
    testConvolveND_padResult()
    {
      // Run the routine under test.
      ArrayND<2, Type> result0 = convolve<Type, Type>(
        m_convolveNDKernel, m_signal, 0, BRICK_CONVOLVE_PAD_RESULT,
        BRICK_CONVOLVE_ROI_SAME);

      // Calculate the correct answer.
      Array2D<Type> kernel(m_convolveNDKernel.size(), 1);
      kernel.copy(m_convolveNDKernel.data());
      Array2D<Type> reference2D = convolve2D<Type, Type>(
        kernel, m_signal2D, BRICK_CONVOLVE_PAD_RESULT,
        BRICK_CONVOLVE_ROI_SAME, static_cast<Type>(0));
      ArrayND<2, Type> reference(reference2D.shape());
      reference.copy(reference2D.data());

      BRICK_TEST_ASSERT(
	this->equivalent(result0, reference, m_defaultTolerance));

      // Repeat for next axis.
      ArrayND<2, Type> result1 = convolve<Type, Type>(
        m_convolveNDKernel, m_signal, 1, BRICK_CONVOLVE_PAD_RESULT,
        BRICK_CONVOLVE_ROI_SAME);

      // Calculate the correct answer.
      kernel.reshape(1, m_convolveNDKernel.size());
      reference2D = convolve2D<Type, Type>(
        kernel, m_signal2D, BRICK_CONVOLVE_PAD_RESULT,
        BRICK_CONVOLVE_ROI_SAME, static_cast<Type>(0));
      reference.reinit(reference2D.shape());
      reference.copy(reference2D.data());

      BRICK_TEST_ASSERT(
	this->equivalent(result1, reference, m_defaultTolerance));
    }


    template <class Type>
    template <class Type2, size_t Dimension>
    bool
    ConvolveNDTest<Type>::
    equivalent(const ArrayND<Dimension, Type2>& array0,
	       const ArrayND<Dimension, Type2>& array1,
	       Type2 tolerance)
    {
      Array1D<size_t> shape0 = array0.getShape();
      Array1D<size_t> shape1 = array1.getShape();
      if(shape0.size() != shape1.size()) {
        return false;
      }
      if(!std::equal(shape0.begin(), shape0.end(), shape1.begin())) {
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
  brick::numeric::ConvolveNDTest<double> currentTest0("double");
  brick::numeric::ConvolveNDTest<float> currentTest1("float");
  brick::numeric::ConvolveNDTest<int> currentTest2("int");
  bool result = (currentTest0.run()
                 && currentTest1.run()
                 && currentTest2.run());
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::ConvolveNDTest<double> currentTest0("double");
  brick::numeric::ConvolveNDTest<float> currentTest1("float");
  brick::numeric::ConvolveNDTest<int> currentTest2("int");

}

#endif
