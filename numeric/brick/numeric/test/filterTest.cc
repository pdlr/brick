/**
***************************************************************************
* @file brick/numeric/test/filterTest.cc
*
* Source file defining tests for functions defined in
* brick/numeric/filter.hh.
*
* Copyright (C) 2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_DEVELOPER
#define BRICK_NUMERIC_DEVELOPER 0
#endif /* #ifndef BRICK_NUMERIC_DEVELOPER */

#include <brick/numeric/filter.hh>
#include <brick/test/testFixture.hh>

#if BRICK_NUMERIC_DEVELOPER
#include <brick/utilities/timeUtilities.hh>
#endif /* #if BRICK_NUMERIC_DEVELOPER */

// Anonymous namespace for locally defined functions.
namespace {

} // namespace

namespace brick {

  namespace numeric {
    
    class FilterTest
      : public brick::test::TestFixture<FilterTest> {

    public:

      FilterTest();
      ~FilterTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testFilterColumns();
      void testFilterRows();

#if BRICK_NUMERIC_DEVELOPER
      void timeFilterColumns();
      void timeFilterRows();
#endif /* #if BRICK_NUMERIC_DEVELOPER */
      
    private:

      Array2D<common::Float64>
      getSignal();

      Array1D<common::Float64>
      getKernel(unsigned int size);

      bool
      isApproximatelyEqual(const Array2D<common::Float64>& array0,
                           const Array2D<common::Float64>& array1);

      void
      referenceFilterColumns(Array2D<common::Float64>& result,
                             Array2D<common::Float64>& signal,
                             Array1D<common::Float64>& kernel);

      void
      referenceFilterRows(Array2D<common::Float64>& result,
                          Array2D<common::Float64>& signal,
                          Array1D<common::Float64>& kernel);


      common::Float64 m_defaultTolerance;
      
    }; // class FilterTest


    /* ============== Member Function Definititions ============== */

    FilterTest::
    FilterTest()
      : TestFixture<FilterTest>("FilterTest"),
        m_defaultTolerance(1.0E-5)
    {
      BRICK_TEST_REGISTER_MEMBER(testFilterColumns);
      BRICK_TEST_REGISTER_MEMBER(testFilterRows);

#if BRICK_NUMERIC_DEVELOPER
      BRICK_TEST_REGISTER_MEMBER(timeFilterColumns);
      BRICK_TEST_REGISTER_MEMBER(timeFilterRows);
#endif /* #if BRICK_NUMERIC_DEVELOPER */
    }


    void
    FilterTest::
    testFilterColumns()
    {
      for(unsigned int ii = 1; ii < 15; ii += 2) {
        Array2D<common::Float64> signal = this->getSignal();
        Array1D<common::Float64> kernel = this->getKernel(ii);

        Array2D<common::Float64> reference(signal.rows(), signal.columns());
        reference = 0.0;
        this->referenceFilterColumns(reference, signal, kernel);
      
        Array2D<common::Float64> result(signal.rows(), signal.columns());
        result = 0.0;
        filterColumns(result, signal, kernel);

        BRICK_TEST_ASSERT(this->isApproximatelyEqual(result, reference));
      }
    }
    

    void
    FilterTest::
    testFilterRows()
    {
      for(unsigned int ii = 1; ii < 15; ii += 2) {
        Array2D<common::Float64> signal = this->getSignal();
        Array1D<common::Float64> kernel = this->getKernel(ii);

        Array2D<common::Float64> reference(signal.rows(), signal.columns());
        reference = 0.0;
        this->referenceFilterRows(reference, signal, kernel);
      
        Array2D<common::Float64> result(signal.rows(), signal.columns());
        result = 0.0;
        filterRows(result, signal, kernel);

        BRICK_TEST_ASSERT(this->isApproximatelyEqual(result, reference));
      }
    }


#if BRICK_NUMERIC_DEVELOPER
    void
    FilterTest::
    timeFilterColumns()
    {
      unsigned int const iterations = 100;

      for(unsigned int ii = 1; ii < 15; ii += 2) {
        Array2D<common::Float64> signal = this->getSignal();
        Array1D<common::Float64> kernel = this->getKernel(ii);

        Array2D<common::Float64> result(signal.rows(), signal.columns());
      
        double t0 = utilities::getCurrentTime();
        for(unsigned int jj = 0; jj < iterations; ++jj) {
          filterColumns(result, signal, kernel);
        }
        double t1 = utilities::getCurrentTime();

        std::cout << "\nAverage ET for filterColumns() with kernel size "
                  << ii << ": " << (t1 - t0) / iterations << std::endl;
      }
    }

    
    void
    FilterTest::
    timeFilterRows()
    {
      unsigned int const iterations = 100;
      
      for(unsigned int ii = 1; ii < 15; ii += 2) {      
        Array2D<common::Float64> signal = this->getSignal();
        Array1D<common::Float64> kernel = this->getKernel(ii);

        Array2D<common::Float64> result(signal.rows(), signal.columns());
      
        double t0 = utilities::getCurrentTime();
        for(unsigned int jj = 0; jj < iterations; ++jj) {
          filterRows(result, signal, kernel);
        }
        double t1 = utilities::getCurrentTime();

        std::cout << "\nAverage ET for filterRows() with kernel size "
                  << ii << ": " << (t1 - t0) / iterations << std::endl;
      }
    }
#endif /* #if BRICK_NUMERIC_DEVELOPER */


    Array2D<common::Float64>
    FilterTest::
    getSignal()
    {
      Array2D<common::Float64> signal(240, 320);
      for(unsigned int rr = 0; rr < signal.rows(); ++rr) {
        for(unsigned int cc = 0; cc < signal.columns(); ++cc) {
          signal(rr, cc) = std::sin((rr + cc) / 50.0) + std::cos(cc / 29.0);
        }
      }
      return signal;
    }

    
    Array1D<common::Float64>
    FilterTest::
    getKernel(unsigned int size)
    {
      Array1D<common::Float64> source(
        "[2.4, 1.67, 8.9, 5.9, 3.3, 0.7, 2.6, 10.0, 1.1, 1.0, 4.9, 5.5, 1.5]");
      Array1D<common::Float64> kernel(size);
      std::copy(source.begin(), source.begin() + size, kernel.begin());
      return kernel;
    }

    
    bool
    FilterTest::
    isApproximatelyEqual(const Array2D<common::Float64>& array0,
                         const Array2D<common::Float64>& array1)
    {
      if(array0.rows() != array1.rows()) {
        return false;
      }
      if(array0.columns() != array1.columns()) {
        return false;
      }
      return std::equal(
        array0.begin(), array0.end(), array1.begin(),
        ApproximatelyEqualFunctor<common::Float64>(m_defaultTolerance));
    }


    void
    FilterTest::
    referenceFilterColumns(Array2D<common::Float64>& result,
                           Array2D<common::Float64>& signal,
                           Array1D<common::Float64>& kernel)
    {
      unsigned int offset = kernel.size() / 2;
      for(unsigned int rr = offset; rr < (signal.rows() - offset); ++rr) {
        for(unsigned int cc = 0; cc < signal.columns(); ++cc) {
          common::Float64 element = 0.0;
          for(unsigned int ii = 0; ii < kernel.size(); ++ii) {
            element += kernel[ii] * signal(rr + ii - offset, cc);
          }
          result(rr, cc) = element;
        }
      }
    }
    

    void
    FilterTest::
    referenceFilterRows(Array2D<common::Float64>& result,
                        Array2D<common::Float64>& signal,
                        Array1D<common::Float64>& kernel)
    {
      unsigned int offset = kernel.size() / 2;
      for(unsigned int rr = 0; rr < signal.rows(); ++rr) {
        for(unsigned int cc = offset; cc < (signal.columns() - offset); ++cc) {
          common::Float64 element = 0.0;
          for(unsigned int ii = 0; ii < kernel.size(); ++ii) {
            element += kernel[ii] * signal(rr, cc + ii - offset);
          }
          result(rr, cc) = element;
        }
      }
    }
    
  } // namespace numeric

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::numeric::FilterTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::FilterTest currentTest;

}

#endif
