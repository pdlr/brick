/**
***************************************************************************
* @file brick/numeric/test/boxIntegrator2DTest.cc
*
* Source file defining BoxIntegrator2DTest class.
*
* Copyright (C) 2006,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/functional.hh>
#include <brick/numeric/boxIntegrator2D.hh>
#include <brick/numeric/subArray2D.hh>
#include <brick/numeric/utilities.hh>
#include <brick/test/testFixture.hh>


namespace brick {

  namespace numeric {

    class BoxIntegrator2DTest
      : public brick::test::TestFixture<BoxIntegrator2DTest> {

    public:

      BoxIntegrator2DTest();
      ~BoxIntegrator2DTest() {};

      void setUp(const std::string& /* testName */);
      void tearDown(const std::string& /* testName */) {}

      // Tests of member functions.
      void testConstructor__void();
      void testConstructor__Array2D();
      void testConstructor__Array2D__Index2D__Index2D();
      void testConstructor__BoxIntegrator2D();
      void testGetIntegral__Index2D__Index2D();
      void testGetIntegral__Index2D__Index2D__bool();
      void testSetArray__Array2D();
      void testSetArray__Array2D__Index2D__Index2D();

    private:

      const double m_defaultTolerance;
      Array2D<double> m_testArray0;
      Array2D<double> m_testArray1;

    }; // class BoxIntegrator2DTest


    /* ============== Member Function Definititions ============== */

    BoxIntegrator2DTest::
    BoxIntegrator2DTest()
      : brick::test::TestFixture<BoxIntegrator2DTest>("BoxIntegrator2DTest"),
        m_defaultTolerance(1.0E-11),
        m_testArray0(),
        m_testArray1()
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testConstructor__void);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__Array2D);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__Array2D__Index2D__Index2D);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__BoxIntegrator2D);
      BRICK_TEST_REGISTER_MEMBER(testGetIntegral__Index2D__Index2D);
      BRICK_TEST_REGISTER_MEMBER(testGetIntegral__Index2D__Index2D__bool);
      BRICK_TEST_REGISTER_MEMBER(testSetArray__Array2D);
      BRICK_TEST_REGISTER_MEMBER(testSetArray__Array2D__Index2D__Index2D);
    }


    void
    BoxIntegrator2DTest::
    setUp(const std::string& /* testName */)
    {
      m_testArray0.reinit(100, 120);
      m_testArray1.reinit(m_testArray0.rows(), m_testArray0.columns());
      for(size_t row = 0; row < m_testArray0.rows(); ++row) {
        for(size_t column = 0; column < m_testArray0.columns(); ++column) {
          m_testArray0(row, column) = row;
            // std::sin(row / 20.0 + column / 30.0) + 5.0;
          m_testArray1(row, column) = column;
            // std::sin(row / 30.0 + column / 20.0) - 5.0;
        }
      }
    }


    void
    BoxIntegrator2DTest::
    testConstructor__void()
    {
      // Nothing really to test.
      BoxIntegrator2D<double, double> boxIntegrator2D;
      BRICK_TEST_ASSERT(true);
    }


    void
    BoxIntegrator2DTest::
    testConstructor__Array2D()
    {
      BoxIntegrator2D<double, double> boxIntegrator2D(m_testArray0);
      for(int regionRows = 6; regionRows < 12; regionRows += 2) {
        for(int regionColumns = 6; regionColumns < 12; regionColumns += 2) {

          for(size_t row0 = 0;
              row0 < m_testArray0.rows() - regionRows;
              row0 += m_testArray0.rows() / 10) {
            for(size_t column0 = 0;
                column0 < m_testArray0.columns() - regionColumns;
                column0 += m_testArray0.columns() / 10) {

              Array2D<double> roi =
                subArray(m_testArray0, Slice(row0, row0 + regionRows),
                         Slice(column0, column0 + regionColumns));
              double referenceValue = sum<double>(roi.ravel());

              double testValue = boxIntegrator2D.getIntegral(
                Index2D(row0, column0),
                Index2D(row0 + regionRows, column0 + regionColumns));

              BRICK_TEST_ASSERT(approximatelyEqual(testValue, referenceValue,
                                                 m_defaultTolerance));
            }
          }
        }
      }
    }


    void
    BoxIntegrator2DTest::
    testConstructor__Array2D__Index2D__Index2D()
    {
      for(int regionRows = 6; regionRows < 12; regionRows += 2) {
        for(int regionColumns = 6; regionColumns < 12; regionColumns += 2) {

          for(int rowOffset = 0;
              rowOffset < static_cast<int>(m_testArray0.rows() / 3);
              rowOffset += (m_testArray0.rows() / 10)) {
            for(int columnOffset = 0;
                columnOffset < static_cast<int>(m_testArray0.columns() / 3);
                columnOffset += (m_testArray0.columns() / 10)) {

              BoxIntegrator2D<double, double> boxIntegrator2D(
                m_testArray0, Index2D(rowOffset, columnOffset),
                Index2D(m_testArray0.rows(), m_testArray0.columns()));

              for(size_t row0 = 0;
                  row0 < (m_testArray0.rows() - rowOffset
                          - regionRows);
                  row0 += m_testArray0.rows() / 10) {
                for(size_t column0 = 0;
                    column0 < (m_testArray0.columns() - columnOffset
                               - regionColumns);
                    column0 += m_testArray0.columns() / 10) {

                  Array2D<double> roi = subArray(
                    m_testArray0,
                    Slice(row0 + rowOffset, row0 + rowOffset + regionRows),
                    Slice(column0 + columnOffset, column0 + columnOffset
                          + regionColumns));
                  double referenceValue = sum<double>(roi.ravel());

                  double testValue = boxIntegrator2D.getIntegral(
                    Index2D(row0, column0),
                    Index2D(row0 + regionRows, column0 + regionColumns));

                  BRICK_TEST_ASSERT(approximatelyEqual(testValue, referenceValue,
                                                     m_defaultTolerance));
                }
              }
            }
          }
        }
      }
    }


    void
    BoxIntegrator2DTest::
    testConstructor__BoxIntegrator2D()
    {
      BoxIntegrator2D<double, double> boxIntegrator2D0(m_testArray0);
      BoxIntegrator2D<double, double> boxIntegrator2D1(boxIntegrator2D0);
      for(int regionRows = 6; regionRows < 12; regionRows += 2) {
        for(int regionColumns = 6; regionColumns < 12; regionColumns += 2) {

          for(size_t row0 = 0;
              row0 < m_testArray0.rows() - regionRows;
              row0 += m_testArray0.rows() / 10) {
            for(size_t column0 = 0;
                column0 < m_testArray0.columns() - regionColumns;
                column0 += m_testArray0.columns() / 10) {

              double testValue0 = boxIntegrator2D0.getIntegral(
                Index2D(row0, column0),
                Index2D(row0 + regionRows, column0 + regionColumns));
              double testValue1 = boxIntegrator2D1.getIntegral(
                Index2D(row0, column0),
                Index2D(row0 + regionRows, column0 + regionColumns));

              BRICK_TEST_ASSERT(approximatelyEqual(testValue0, testValue1,
                                                 m_defaultTolerance));
            }
          }
        }
      }
    }


    void
    BoxIntegrator2DTest::
    testGetIntegral__Index2D__Index2D()
    {
      // Already tested by constructor tests.
    }


    void
    BoxIntegrator2DTest::
    testGetIntegral__Index2D__Index2D__bool()
    {
      for(int regionRows = 6; regionRows < 12; regionRows += 2) {
        for(int regionColumns = 6; regionColumns < 12; regionColumns += 2) {

          for(int rowOffset = 0;
              rowOffset < static_cast<int>(m_testArray0.rows() / 3);
              rowOffset += static_cast<int>(m_testArray0.rows() / 10)) {
            for(int columnOffset = 0;
                columnOffset < static_cast<int>(m_testArray0.columns() / 3);
                columnOffset += static_cast<int>(m_testArray0.columns() / 10)) {

              BoxIntegrator2D<double, double> boxIntegrator2D(
                m_testArray0, Index2D(rowOffset, columnOffset),
                Index2D(m_testArray0.rows(), m_testArray0.columns()));

              for(size_t row0 = rowOffset;
                  row0 < (m_testArray0.rows() - regionRows);
                  row0 += m_testArray0.rows() / 10) {
                for(size_t column0 = columnOffset;
                    column0 < (m_testArray0.columns() - regionColumns);
                    column0 += m_testArray0.columns() / 10) {

                  Array2D<double> roi = subArray(
                    m_testArray0, Slice(row0, row0 + regionRows),
                    Slice(column0, column0 + regionColumns));
                  double referenceValue = sum<double>(roi.ravel());

                  double testValue = boxIntegrator2D.getIntegral(
                    Index2D(row0, column0),
                    Index2D(row0 + regionRows, column0 + regionColumns),
                    true);

                  BRICK_TEST_ASSERT(approximatelyEqual(testValue, referenceValue,
                                                     m_defaultTolerance));
                }
              }
            }
          }
        }
      }
    }


    void
    BoxIntegrator2DTest::
    testSetArray__Array2D()
    {
      BoxIntegrator2D<double, double> boxIntegrator2D0(m_testArray0);
      BoxIntegrator2D<double, double> boxIntegrator2D1(boxIntegrator2D0);
      boxIntegrator2D1.setArray(m_testArray1);
      for(int regionRows = 6; regionRows < 12; regionRows += 2) {
        for(int regionColumns = 6; regionColumns < 12; regionColumns += 2) {

          for(size_t row0 = 0;
              row0 < m_testArray0.rows() - regionRows;
              row0 += m_testArray0.rows() / 10) {
            for(size_t column0 = 0;
                column0 < m_testArray0.columns() - regionColumns;
                column0 += m_testArray0.columns() / 10) {

              Array2D<double> roi0 =
                subArray(m_testArray0, Slice(row0, row0 + regionRows),
                         Slice(column0, column0 + regionColumns));
              double referenceValue0 = sum<double>(roi0.ravel());
              Array2D<double> roi1 =
                subArray(m_testArray1, Slice(row0, row0 + regionRows),
                         Slice(column0, column0 + regionColumns));
              double referenceValue1 = sum<double>(roi1.ravel());

              double testValue0 = boxIntegrator2D0.getIntegral(
                Index2D(row0, column0),
                Index2D(row0 + regionRows, column0 + regionColumns));
              double testValue1 = boxIntegrator2D1.getIntegral(
                Index2D(row0, column0),
                Index2D(row0 + regionRows, column0 + regionColumns));

              BRICK_TEST_ASSERT(approximatelyEqual(testValue0, referenceValue0,
                                                 m_defaultTolerance));
              BRICK_TEST_ASSERT(approximatelyEqual(testValue1, referenceValue1,
                                                 m_defaultTolerance));
            }
          }
        }
      }
    }


    void
    BoxIntegrator2DTest::
    testSetArray__Array2D__Index2D__Index2D()
    {
      for(int regionRows = 6; regionRows < 12; regionRows += 2) {
        for(int regionColumns = 6; regionColumns < 12; regionColumns += 2) {

          for(int rowOffset = 0;
              rowOffset < static_cast<int>(m_testArray0.rows() / 3);
              rowOffset += (m_testArray0.rows() / 10)) {
            for(int columnOffset = 0;
                columnOffset < static_cast<int>(m_testArray0.columns() / 3);
                columnOffset += (m_testArray0.columns() / 10)) {

              BoxIntegrator2D<double, double> boxIntegrator2D0(m_testArray0);
              BoxIntegrator2D<double, double> boxIntegrator2D1(
                boxIntegrator2D0);
              boxIntegrator2D1.setArray(
                m_testArray0, Index2D(rowOffset, columnOffset),
                Index2D(m_testArray0.rows(), m_testArray0.columns()));

              for(size_t row0 = 0;
                  row0 < (m_testArray0.rows() - rowOffset
                          - regionRows);
                  row0 += m_testArray0.rows() / 10) {
                for(size_t column0 = 0;
                    column0 < (m_testArray0.columns() - columnOffset
                               - regionColumns);
                    column0 += m_testArray0.columns() / 10) {

                  Array2D<double> roi = subArray(
                    m_testArray0,
                    Slice(row0 + rowOffset, row0 + rowOffset + regionRows),
                    Slice(column0 + columnOffset, column0 + columnOffset
                          + regionColumns));
                  double referenceValue = sum<double>(roi.ravel());

                  double testValue = boxIntegrator2D1.getIntegral(
                    Index2D(row0, column0),
                    Index2D(row0 + regionRows, column0 + regionColumns));

                  BRICK_TEST_ASSERT(approximatelyEqual(testValue, referenceValue,
                                                     m_defaultTolerance));
                }
              }
            }
          }
        }
      }
    }


  } // namespace numeric

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::numeric::BoxIntegrator2DTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
  return 0;
}

#else

namespace {

  brick::numeric::BoxIntegrator2DTest currentTest;

}

#endif
