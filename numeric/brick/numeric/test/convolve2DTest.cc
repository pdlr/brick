/**
***************************************************************************
* @file convolve2DTest.cpp
*
* Source file defining Convolve2DTest class.
*
* Copyright (C) 2007 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/functional.hh>
#include <brick/numeric/convolve2D.hh>
#include <brick/test/testFixture.hh>

namespace brick {

  namespace numeric {
  
    template <class Type>
    class Convolve2DTest
      : public brick::test::TestFixture< Convolve2DTest<Type> >
    {

    public:

      // Typedef required why?
      typedef Convolve2DTest<Type> TestFixtureType;
    
    
      Convolve2DTest(const std::string& typeName);
      ~Convolve2DTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      void testConvolve2D_truncateResult();
      void testConvolve2D_padResult();
      void testConvolve2D_padSignal();
      void testConvolve2D_zeroPadSignal();
      void testConvolve2D_reflectSignal();
      void testConvolve2D_wrapSignal();
      void testCorrelate2D_truncateResult();
      void testCorrelate2D_padResult();
      void testCorrelate2D_padSignal();
      void testCorrelate2D_zeroPadSignal();
      void testCorrelate2D_reflectSignal();
      void testCorrelate2D_wrapSignal();

    private:

      template <class Type2>
      bool equivalent(const Array2D<Type2>& arg0,
		      const Array2D<Type2>& arg1,
		      Type2 tolerance);


      Array2D<Type> m_convolve2DKernel;
      Array2D<Type> m_correlate2DKernel;
      Type m_defaultTolerance;
      Type m_fillValue;
      Array2D<Type> m_signal;
      Array2D<Type> m_result_truncateResult;
      Array2D<Type> m_result_padResult;
      Array2D<Type> m_result_padSignal;
      Array2D<Type> m_result_zeroPadSignal;
      Array2D<Type> m_result_reflectSignal;
      Array2D<Type> m_result_wrapSignal;
    
    }; // class Convolve2DTest


    /* ============== Member Function Definititions ============== */

    template <class Type>
    Convolve2DTest<Type>::
    Convolve2DTest(const std::string& typeName)
      : brick::test::TestFixture<Convolve2DTest>(
	std::string("Convolve2DTest<" + typeName + ">")),
	m_convolve2DKernel("[[6, 4, 2],"
			   " [1, 5, 0],"
			   " [0, 3, 3],"
			   " [4, 3, 4],"
			   " [2, 3, 4]]"),
	m_correlate2DKernel("[[4, 3, 2],"
			    " [4, 3, 4],"
			    " [3, 3, 0],"
			    " [0, 5, 1],"
			    " [2, 4, 6]]"),
	m_defaultTolerance(static_cast<Type>(1.0E-6)),
	m_fillValue(static_cast<Type>(3)),
	m_signal("[[1, 2, 3, 4, 5, 6],"
		 " [2, 3, 4, 5, 6, 7],"
		 " [0, 1, 2, 3, 4, 5],"
		 " [3, 2, 1, 1, 2, 3],"
		 " [4, 4, 3, 1, 2, 0],"
		 " [6, 5, 4, 3, 2, 1],"
		 " [1, 1, 5, 7, 0, 1],"
		 " [1, 2, 3, 4, 3, 2],"
		 " [1, 2, 4, 5, 7, 8]]"),
	m_result_truncateResult("[[105, 110, 133, 153],"
				" [130, 125, 121, 135],"
				" [118, 139, 107,  96],"
				" [131, 141, 130,  74],"
				" [142, 158, 181, 156]]"),
	m_result_padResult("[[3, 3, 3, 3, 3, 3],"
			   " [3, 3, 3, 3, 3, 3],"
			   " [3, 105, 110, 133, 153, 3],"
			   " [3, 130, 125, 121, 135, 3],"
			   " [3, 118, 139, 107,  96, 3],"
			   " [3, 131, 141, 130,  74, 3],"
			   " [3, 142, 158, 181, 156, 3],"
			   " [3, 3, 3, 3, 3, 3],"
			   " [3, 3, 3, 3, 3, 3]]"),
	m_result_padSignal("[[120, 118, 124, 136, 148, 160, 148, 138],"
			   " [124, 117, 131, 149, 167, 185, 169, 140],"
			   " [113,  97, 104, 128, 152, 176, 177, 145],"
			   " [121,  96,  91, 108, 135, 168, 178, 156],"
			   " [130, 121, 105, 110, 133, 153, 168, 160],"
			   " [137, 142, 130, 125, 121, 135, 138, 152],"
			   " [117, 115, 118, 139, 107,  96, 102, 127],"
			   " [122, 120, 131, 141, 130,  74,  92, 112],"
			   " [132, 123, 142, 158, 181, 156, 117, 116],"
			   " [128, 114, 133, 161, 159, 152, 126, 113],"
			   " [120, 106, 102, 138, 158, 153, 138, 135],"
			   " [120, 114, 114, 137, 161, 174, 160, 148],"
			   " [128, 124, 123, 135, 150, 162, 163, 152]]"),
	m_result_zeroPadSignal("[[  6,  16,  28,  40,  52,  64,  34,  12],"
			       " [ 13,  33,  53,  71,  89, 107,  70,  14],"
			       " [  2,  22,  44,  68,  92, 116,  96,  28],"
			       " [ 22,  42,  64,  81, 108, 141, 118,  51],"
			       " [ 37,  82, 105, 110, 133, 153, 129,  67],"
			       " [ 44, 103, 130, 125, 121, 135,  99,  59],"
			       " [ 24,  76, 118, 139, 107,  96,  63,  34],"
			       " [ 29,  81, 131, 141, 130,  74,  53,  19],"
			       " [ 39,  84, 142, 158, 181, 156,  78,  23],"
			       " [ 17,  45,  97, 125, 123, 116,  69,  14],"
			       " [  6,  19,  48,  84, 104,  99,  66,  36],"
			       " [  6,  18,  42,  65,  89, 102,  70,  40],"
			       " [  2,   7,  18,  30,  45,  57,  52,  32]]"),
	m_result_reflectSignal("[[ 62,  66,  87, 115, 148, 188, 219, 219],"
			       " [ 71,  71, 102, 146, 190, 234, 265, 265],"
			       " [ 60,  60,  91, 135, 179, 223, 254, 254],"
			       " [ 77,  69,  80, 106, 142, 184, 215, 215],"
			       " [108, 102, 105, 110, 133, 153, 172, 184],"
			       " [135, 132, 130, 125, 121, 135, 139, 150],"
			       " [105, 102, 118, 139, 107,  96,  92,  85],"
			       " [120, 129, 131, 141, 130,  74,  72,  74],"
			       " [124, 129, 142, 158, 181, 156, 132, 128],"
			       " [ 90,  94, 131, 175, 193, 202, 177, 167],"
			       " [ 53,  55,  90, 149, 176, 174, 150, 158],"
			       " [ 55,  51, 100, 166, 177, 175, 173, 171],"
			       " [113, 104, 119, 161, 192, 172, 187, 189]]"),
	m_result_wrapSignal("[[126,  87, 125, 165, 175, 180, 126,  87],"
			    " [155, 102, 101, 155, 193, 206, 155, 102],"
			    " [174, 108,  86, 133, 181, 218, 174, 108],"
			    " [194, 132,  82, 111, 153, 198, 194, 132],"
			    " [166, 149, 105, 110, 133, 153, 166, 149],"
			    " [143, 162, 130, 125, 121, 135, 143, 162],"
			    " [ 87, 110, 118, 139, 107,  96,  87, 110],"
			    " [ 82, 100, 131, 141, 130,  74,  82, 100],"
			    " [117, 107, 142, 158, 181, 156, 117, 107],"
			    " [126,  87, 125, 165, 175, 180, 126,  87],"
			    " [155, 102, 101, 155, 193, 206, 155, 102],"
			    " [174, 108,  86, 133, 181, 218, 174, 108],"
			    " [194, 132,  82, 111, 153, 198, 194, 132]]")
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testConvolve2D_truncateResult);
      BRICK_TEST_REGISTER_MEMBER(testConvolve2D_padResult);
      BRICK_TEST_REGISTER_MEMBER(testConvolve2D_padSignal);
      BRICK_TEST_REGISTER_MEMBER(testConvolve2D_zeroPadSignal);
      BRICK_TEST_REGISTER_MEMBER(testConvolve2D_reflectSignal);
      BRICK_TEST_REGISTER_MEMBER(testConvolve2D_wrapSignal);
      BRICK_TEST_REGISTER_MEMBER(testCorrelate2D_truncateResult);
      BRICK_TEST_REGISTER_MEMBER(testCorrelate2D_padResult);
      BRICK_TEST_REGISTER_MEMBER(testCorrelate2D_padSignal);
      BRICK_TEST_REGISTER_MEMBER(testCorrelate2D_zeroPadSignal);
      BRICK_TEST_REGISTER_MEMBER(testCorrelate2D_reflectSignal);
      BRICK_TEST_REGISTER_MEMBER(testCorrelate2D_wrapSignal);
    }


    template <class Type>
    void
    Convolve2DTest<Type>::
    testConvolve2D_truncateResult()
    {
      Array2D<Type> result = convolve2D<Type, Type>(
	m_convolve2DKernel, m_signal, BRICK_CONVOLVE_TRUNCATE_RESULT,
	BRICK_CONVOLVE_ROI_VALID);
      BRICK_TEST_ASSERT(
	this->equivalent(result, m_result_truncateResult, m_defaultTolerance));
    }

  
    template <class Type>
    void
    Convolve2DTest<Type>::
    testConvolve2D_padResult()
    {
      Array2D<Type> result = convolve2D<Type, Type>(m_convolve2DKernel, m_signal,
					      BRICK_CONVOLVE_PAD_RESULT,
					      BRICK_CONVOLVE_ROI_SAME,
					      m_fillValue);
      BRICK_TEST_ASSERT(
	this->equivalent(result, m_result_padResult, m_defaultTolerance));
    }


    template <class Type>
    void
    Convolve2DTest<Type>::
    testConvolve2D_padSignal()
    {
      Array2D<Type> result = convolve2D<Type, Type>(m_convolve2DKernel, m_signal,
					      BRICK_CONVOLVE_PAD_SIGNAL,
					      BRICK_CONVOLVE_ROI_FULL,
					      m_fillValue);
      BRICK_TEST_ASSERT(
	this->equivalent(result, m_result_padSignal, m_defaultTolerance));
    }


    template <class Type>
    void
    Convolve2DTest<Type>::
    testConvolve2D_zeroPadSignal()
    {
      Array2D<Type> result = convolve2D<Type, Type>(m_convolve2DKernel, m_signal,
					      BRICK_CONVOLVE_ZERO_PAD_SIGNAL,
					      BRICK_CONVOLVE_ROI_FULL);
      BRICK_TEST_ASSERT(
	this->equivalent(result, m_result_zeroPadSignal, m_defaultTolerance));
    }


    template <class Type>
    void
    Convolve2DTest<Type>::
    testConvolve2D_reflectSignal()
    {
      Array2D<Type> result = convolve2D<Type, Type>(m_convolve2DKernel, m_signal,
					      BRICK_CONVOLVE_REFLECT_SIGNAL,
					      BRICK_CONVOLVE_ROI_FULL);
      BRICK_TEST_ASSERT(
	this->equivalent(result, m_result_reflectSignal, m_defaultTolerance));
    }


    template <class Type>
    void
    Convolve2DTest<Type>::
    testConvolve2D_wrapSignal()
    {
      Array2D<Type> result = convolve2D<Type, Type>(m_convolve2DKernel, m_signal,
					      BRICK_CONVOLVE_WRAP_SIGNAL,
					      BRICK_CONVOLVE_ROI_FULL);
      BRICK_TEST_ASSERT(
	this->equivalent(result, m_result_wrapSignal, m_defaultTolerance));
    }


    template <class Type>
    void
    Convolve2DTest<Type>::
    testCorrelate2D_truncateResult()
    {
      Array2D<Type> result = correlate2D<Type, Type>(m_correlate2DKernel, m_signal,
					       BRICK_CONVOLVE_TRUNCATE_RESULT,
					       BRICK_CONVOLVE_ROI_VALID);
      BRICK_TEST_ASSERT(
	this->equivalent(result, m_result_truncateResult, m_defaultTolerance));
    }

  
    template <class Type>
    void
    Convolve2DTest<Type>::
    testCorrelate2D_padResult()
    {
      Array2D<Type> result = correlate2D<Type, Type>(m_correlate2DKernel, m_signal,
					       BRICK_CONVOLVE_PAD_RESULT,
					       BRICK_CONVOLVE_ROI_SAME,
					       m_fillValue);
      BRICK_TEST_ASSERT(
	this->equivalent(result, m_result_padResult, m_defaultTolerance));
    }


    template <class Type>
    void
    Convolve2DTest<Type>::
    testCorrelate2D_padSignal()
    {
      Array2D<Type> result = correlate2D<Type, Type>(m_correlate2DKernel, m_signal,
					       BRICK_CONVOLVE_PAD_SIGNAL,
					       BRICK_CONVOLVE_ROI_FULL,
					       m_fillValue);
      BRICK_TEST_ASSERT(
	this->equivalent(result, m_result_padSignal, m_defaultTolerance));
    }


    template <class Type>
    void
    Convolve2DTest<Type>::
    testCorrelate2D_zeroPadSignal()
    {
      Array2D<Type> result = correlate2D<Type, Type>(m_correlate2DKernel, m_signal,
					       BRICK_CONVOLVE_ZERO_PAD_SIGNAL,
					       BRICK_CONVOLVE_ROI_FULL);
      BRICK_TEST_ASSERT(
	this->equivalent(result, m_result_zeroPadSignal, m_defaultTolerance));
    }


    template <class Type>
    void
    Convolve2DTest<Type>::
    testCorrelate2D_reflectSignal()
    {
      Array2D<Type> result = correlate2D<Type, Type>(m_correlate2DKernel, m_signal,
					       BRICK_CONVOLVE_REFLECT_SIGNAL,
					       BRICK_CONVOLVE_ROI_FULL);
      BRICK_TEST_ASSERT(
	this->equivalent(result, m_result_reflectSignal, m_defaultTolerance));
    }


    template <class Type>
    void
    Convolve2DTest<Type>::
    testCorrelate2D_wrapSignal()
    {
      Array2D<Type> result = correlate2D<Type, Type>(
	m_correlate2DKernel, m_signal, BRICK_CONVOLVE_WRAP_SIGNAL,
	BRICK_CONVOLVE_ROI_FULL);
      BRICK_TEST_ASSERT(
	this->equivalent(result, m_result_wrapSignal, m_defaultTolerance));
    }


    template <class Type>
    template <class Type2>
    bool
    Convolve2DTest<Type>::
    equivalent(const Array2D<Type2>& array0,
	       const Array2D<Type2>& array1,
	       Type2 tolerance)
    {
      if(array0.rows() != array1.rows()) {
	return false;
      }
      if(array0.columns() != array1.columns()) {
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
  brick::numeric::Convolve2DTest<double> currentTest0("double");
  brick::numeric::Convolve2DTest<float> currentTest1("float");
  brick::numeric::Convolve2DTest<int> currentTest2("int");
  bool result = (currentTest0.run()
                 && currentTest1.run()
                 && currentTest2.run());
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::Convolve2DTest<double> currentTest0("double");
  brick::numeric::Convolve2DTest<float> currentTest1("float");
  brick::numeric::Convolve2DTest<int> currentTest2("int");

}

#endif
