/**
***************************************************************************
* @file brick/numeric/test/utilitiesTest.cpp
*
* Source file defining UtilitiesTest class.
*
* Copyright (C) 2005 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/common/functional.hh>
#include <brick/numeric/utilities.hh>
#include <brick/test/testFixture.hh>

namespace brick {

  namespace numeric {
    template <class Type>
    class UtilitiesTest
      : public test::TestFixture< UtilitiesTest<Type> >
    {

    public:

      typedef UtilitiesTest<Type> TestFixtureType;


      UtilitiesTest(const std::string& typeName);
      ~UtilitiesTest() {}

      void setUp(const std::string&) {}
      void tearDown(const std::string&) {}

      void testAbs0();
      void testAbs1();
      void testAllTrueEtc();
      void testArgmax0();
      void testArgmax1();
      void testArgmax2D0();
      void testArgmax2D1();
      void testArgmin0();
      void testArgmin1();
      void testArgmin2D0();
      void testArgmin2D1();
      void testArgsort0();
      void testAxisMax();
      void testAxisMaximum0();
      void testAxisMaximum1();
      void testAxisMin();
      void testAxisMinimum0();
      void testAxisMinimum1();
      void testAxisSum0();
      void testAxisSum1();
      void testAxisSum2();
      void testColumnIndices();
      void testCompress0();
      void testCompress1();
      void testCount();
      void testGetCentroid();
      void testGetMeanAndCovariance();
      void testMaximum();
      void testNormalizedCorrelation();
      void testSum0();
      void testTake_Array1D_Array1D();
      void testTake_Array2D_Array1D();
      void testTake_Array2D_Array1D_unsignedInt();

    private:

      template <class Type2>
      bool equivalent(const Array1D<Type2>& arg0,
                      const Array1D<Type2>& arg1,
                      Type2 tolerance);

      template <class Type2>
      bool equivalent(const Array2D<Type2>& arg0, const Array2D<Type2>& arg1,
                      Type2 tolerance);

      Array1D<double>
      getRandomSample(size_t dimensionality);

      double m_defaultTolerance;

    }; // class UtilitiesTest


    /* ============== Member Function Definititions ============== */

    template <class Type>
    UtilitiesTest<Type>::
    UtilitiesTest(const std::string& typeName)
      : test::TestFixture<UtilitiesTest>(
        std::string("UtilitiesTest<" + typeName + ">")),
        m_defaultTolerance(1.0E-6)
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testAbs0);
      BRICK_TEST_REGISTER_MEMBER(testAbs1);
      BRICK_TEST_REGISTER_MEMBER(testAllTrueEtc);
      BRICK_TEST_REGISTER_MEMBER(testArgmax0);
      BRICK_TEST_REGISTER_MEMBER(testArgmax1);
      BRICK_TEST_REGISTER_MEMBER(testArgmax2D0);
      BRICK_TEST_REGISTER_MEMBER(testArgmax2D1);
      BRICK_TEST_REGISTER_MEMBER(testArgmin0);
      BRICK_TEST_REGISTER_MEMBER(testArgmin1);
      BRICK_TEST_REGISTER_MEMBER(testArgmin2D0);
      BRICK_TEST_REGISTER_MEMBER(testArgmin2D1);
      BRICK_TEST_REGISTER_MEMBER(testArgsort0);
      BRICK_TEST_REGISTER_MEMBER(testAxisMax);
      BRICK_TEST_REGISTER_MEMBER(testAxisMaximum0);
      BRICK_TEST_REGISTER_MEMBER(testAxisMaximum1);
      BRICK_TEST_REGISTER_MEMBER(testAxisMin);
      BRICK_TEST_REGISTER_MEMBER(testAxisMinimum0);
      BRICK_TEST_REGISTER_MEMBER(testAxisMinimum1);
      BRICK_TEST_REGISTER_MEMBER(testAxisSum0);
      BRICK_TEST_REGISTER_MEMBER(testAxisSum1);
      BRICK_TEST_REGISTER_MEMBER(testAxisSum2);
      BRICK_TEST_REGISTER_MEMBER(testColumnIndices);
      BRICK_TEST_REGISTER_MEMBER(testCompress0);
      BRICK_TEST_REGISTER_MEMBER(testCompress1);
      BRICK_TEST_REGISTER_MEMBER(testCount);
      BRICK_TEST_REGISTER_MEMBER(testGetCentroid);
      BRICK_TEST_REGISTER_MEMBER(testGetMeanAndCovariance);
      BRICK_TEST_REGISTER_MEMBER(testMaximum);
      BRICK_TEST_REGISTER_MEMBER(testSum0);
      BRICK_TEST_REGISTER_MEMBER(testTake_Array1D_Array1D);
      BRICK_TEST_REGISTER_MEMBER(testTake_Array2D_Array1D);
      BRICK_TEST_REGISTER_MEMBER(testTake_Array2D_Array1D_unsignedInt);
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testAbs0()
    {
      Array1D<Type> testArray0("[-4, -3, -2, -1, 0, 1, 2, 3]");
      Array1D<Type> testArray1("[4, 3, 2, 1, 0, 1, 2, 3]");
      Array1D<Type> testArray2 = abs(testArray0);
      BRICK_TEST_ASSERT(
        this->equivalent(
          testArray1, testArray2, static_cast<Type>(m_defaultTolerance)));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testAbs1()
    {
      Array2D<Type> testArray0("[[-4, -3, -2],"
                               " [-1, 0, 1]]");
      Array2D<Type> testArray1("[[4, 3, 2],"
                               " [1, 0, 1]]");
      Array2D<Type> testArray2 = abs(testArray0);
      BRICK_TEST_ASSERT(
        this->equivalent(
          testArray1, testArray2, static_cast<Type>(m_defaultTolerance)));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testAllTrueEtc()
    {
      Array1D<Type> testArray0("[-2, -1, -10, 5, 0, 3, 1]");
      Array1D<Type> testArray1("[-2, -1, -10, 5, -1, 3, 1]");
      Array1D<Type> testArray2("[0, 0, 1, 0, 0, 0, 0, 0]");
      Array1D<Type> testArray3("[0, 0, 0, 0, 0, 0, 0, 0]");

      BRICK_TEST_ASSERT(!allTrue(testArray0));
      BRICK_TEST_ASSERT(!allFalse(testArray0));
      BRICK_TEST_ASSERT(anyTrue(testArray0));
      BRICK_TEST_ASSERT(anyFalse(testArray0));

      BRICK_TEST_ASSERT(allTrue(testArray1));
      BRICK_TEST_ASSERT(!allFalse(testArray1));
      BRICK_TEST_ASSERT(anyTrue(testArray1));
      BRICK_TEST_ASSERT(!anyFalse(testArray1));

      BRICK_TEST_ASSERT(!allTrue(testArray2));
      BRICK_TEST_ASSERT(!allFalse(testArray2));
      BRICK_TEST_ASSERT(anyTrue(testArray2));
      BRICK_TEST_ASSERT(anyFalse(testArray2));

      BRICK_TEST_ASSERT(!allTrue(testArray3));
      BRICK_TEST_ASSERT(allFalse(testArray3));
      BRICK_TEST_ASSERT(!anyTrue(testArray3));
      BRICK_TEST_ASSERT(anyFalse(testArray3));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testArgmax0()
    {
      Array1D<Type> testArray0("[-2, -1, -10, 5, 2, 3, 0]");
      BRICK_TEST_ASSERT(argmax(testArray0) == 3);
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testArgmax1()
    {
      Array1D<Type> testArray0("[-2, -1, -10, 5, 2, 3, 0]");
      BRICK_TEST_ASSERT(argmax(testArray0, std::greater<Type>()) == 2);
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testArgmax2D0()
    {
      Array2D<Type> testArray0("[[-2, -1,   -10,   5,  2,   3,  0],"
                               " [ 7, -11,  -10,   5, -1,  12,  6],"
                               " [-4, -100, -10,   9,  4,  -3,  5],"
                               " [ 2,  1,    4,   -14,  6,   3,  7]]");
      BRICK_TEST_ASSERT(argmax2D(testArray0) == Index2D(1, 5));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testArgmax2D1()
    {
      Array2D<Type> testArray0("[[-2, -1,   -10,   5,  2,   3,  0],"
                               " [ 7, -11,  -10,   5, -1,  12,  6],"
                               " [-4, -100, -10,   9,  4,  -3,  5],"
                               " [ 2,  1,    4,   -14,  6,   3,  7]]");
      BRICK_TEST_ASSERT(argmax2D(testArray0, std::greater<Type>())
                        == Index2D(2, 1));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testArgmin0()
    {
      Array1D<Type> testArray0("[-2, -1, -10, 5, 2, 3, 0]");
      BRICK_TEST_ASSERT(argmin(testArray0, std::greater<Type>()) == 3);
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testArgmin1()
    {
      Array1D<Type> testArray0("[-2, -1, -10, 5, 2, 3, 0]");
      BRICK_TEST_ASSERT(argmin(testArray0, std::greater<Type>()) == 3);
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testArgmin2D0()
    {
      Array2D<Type> testArray0("[[-2, -1,   -10,   5,  2,   3,  0],"
                               " [ 7, -11,  -10,   5, -1,  12,  6],"
                               " [-4, -100, -10,   9,  4,  -3,  5],"
                               " [ 2,  1,    4,   -14,  6,   3,  7]]");
      BRICK_TEST_ASSERT(argmin2D(testArray0) == Index2D(2, 1));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testArgmin2D1()
    {
      Array2D<Type> testArray0("[[-2, -1,   -10,   5,  2,   3,  0],"
                               " [ 7, -11,  -10,   5, -1,  12,  6],"
                               " [-4, -100, -10,   9,  4,  -3,  5],"
                               " [ 2,  1,    4,   -14,  6,   3,  7]]");
      BRICK_TEST_ASSERT(argmin2D(testArray0, std::greater<Type>())
                        == Index2D(1, 5));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testArgsort0()
    {
      Array1D<Type> testArray0("[-4, 3, -2, 10, 6, 5, 7, 1]");
      Array1D<size_t> referenceArray0("[0, 2, 7, 1, 5, 4, 6, 3]");
      Array1D<size_t> resultArray0 = argsort(testArray0);
      BRICK_TEST_ASSERT(
        this->equivalent(resultArray0, referenceArray0, size_t(0)));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testAxisMax()
    {
      Array2D<Type> testArray0("[[-2, -1, -10],"
                               " [5, 2, 3],"
                               " [10, 0, -2],"
                               " [1, 2, 2]]");
      Array1D<Type> maxOverRows("[10, 2, 3]");
      Array1D<Type> maxOverColumns("[-1, 5, 10, 2]");
      BRICK_TEST_ASSERT(
        this->equivalent(axisMax(testArray0, 0), maxOverRows,
                         static_cast<Type>(m_defaultTolerance)));
      BRICK_TEST_ASSERT(
        this->equivalent(axisMax(testArray0, 1), maxOverColumns,
                         static_cast<Type>(m_defaultTolerance)));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testAxisMaximum0()
    {
      Array2D<Type> testArray0("[[-2, -1, -10],"
                               " [5, 2, 3],"
                               " [10, 0, -2],"
                               " [1, 2, 2]]");
      Array1D<Type> maxOverRows("[10, 2, 3]");
      Array1D<Type> maxOverColumns("[-1, 5, 10, 2]");
      BRICK_TEST_ASSERT(
        this->equivalent(axisMaximum(testArray0, 0), maxOverRows,
                         static_cast<Type>(m_defaultTolerance)));
      BRICK_TEST_ASSERT(
        this->equivalent(axisMaximum(testArray0, 1), maxOverColumns,
                         static_cast<Type>(m_defaultTolerance)));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testAxisMaximum1()
    {
      Array2D<Type> testArray0("[[-2, -1, -10],"
                               " [5, 2, 3],"
                               " [10, 0, -2],"
                               " [1, 2, 2]]");
      Array1D<Type> minOverRows("[-2, -1, -10]");
      Array1D<Type> minOverColumns("[-10, 2, -2, 1]");
      BRICK_TEST_ASSERT(
        this->equivalent(axisMaximum(testArray0, 0, std::greater<Type>()),
                         minOverRows, static_cast<Type>(m_defaultTolerance)));
      BRICK_TEST_ASSERT(
        this->equivalent(axisMaximum(testArray0, 1, std::greater<Type>()),
                         minOverColumns, static_cast<Type>(m_defaultTolerance)));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testAxisMin()
    {
      Array2D<Type> testArray0("[[-2, -1, -10],"
                               " [5, 2, 3],"
                               " [10, 0, -2],"
                               " [1, 2, 2]]");
      Array1D<Type> minOverRows("[-2, -1, -10]");
      Array1D<Type> minOverColumns("[-10, 2, -2, 1]");
      BRICK_TEST_ASSERT(
        this->equivalent(axisMin(testArray0, 0), minOverRows,
                         static_cast<Type>(m_defaultTolerance)));
      BRICK_TEST_ASSERT(
        this->equivalent(axisMin(testArray0, 1), minOverColumns,
                         static_cast<Type>(m_defaultTolerance)));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testAxisMinimum0()
    {
      Array2D<Type> testArray0("[[-2, -1, -10],"
                               " [5, 2, 3],"
                               " [10, 0, -2],"
                               " [1, 2, 2]]");
      Array1D<Type> minOverRows("[-2, -1, -10]");
      Array1D<Type> minOverColumns("[-10, 2, -2, 1]");
      BRICK_TEST_ASSERT(
        this->equivalent(axisMinimum(testArray0, 0), minOverRows,
                         static_cast<Type>(m_defaultTolerance)));
      BRICK_TEST_ASSERT(
        this->equivalent(axisMinimum(testArray0, 1), minOverColumns,
                         static_cast<Type>(m_defaultTolerance)));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testAxisMinimum1()
    {
      Array2D<Type> testArray0("[[-2, -1, -10],"
                               " [5, 2, 3],"
                               " [10, 0, -2],"
                               " [1, 2, 2]]");
      Array1D<Type> maxOverRows("[10, 2, 3]");
      Array1D<Type> maxOverColumns("[-1, 5, 10, 2]");
      BRICK_TEST_ASSERT(
        this->equivalent(axisMinimum(testArray0, 0, std::greater<Type>()),
                         maxOverRows, static_cast<Type>(m_defaultTolerance)));
      BRICK_TEST_ASSERT(
        this->equivalent(axisMinimum(testArray0, 1, std::greater<Type>()),
                         maxOverColumns, static_cast<Type>(m_defaultTolerance)));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testAxisSum0()
    {
      typedef typename ArithmeticTraits<Type, Type>::SumType SumType;

      Array2D<Type> testArray0("[[-2, -1, -10],"
                               " [5, 2, 3],"
                               " [10, 0, -2],"
                               " [1, 2, 2]]");
      Array1D<SumType> sumOverRows("[14, 3, -7]");
      Array1D<SumType> sumOverColumns("[-13, 10, 8, 5]");
      Array1D<SumType> computedSum0 = axisSum<SumType>(testArray0, 0);
      Array1D<SumType> computedSum1 = axisSum<SumType>(testArray0, 1);

      BRICK_TEST_ASSERT(
        this->equivalent(computedSum0, sumOverRows,
                         static_cast<SumType>(m_defaultTolerance)));
      BRICK_TEST_ASSERT(
        this->equivalent(computedSum1, sumOverColumns,
                         static_cast<SumType>(m_defaultTolerance)));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testAxisSum1()
    {
      Array2D<Type> testArray0("[[-2, -1, -10],"
                               " [5, 2, 3],"
                               " [10, 0, -2],"
                               " [1, 2, 2]]");
      Array1D<double> sumOverRows("[14, 3, -7]");
      Array1D<double> sumOverColumns("[-13, 10, 8, 5]");
      Array1D<double> computedSum0 = axisSum<double>(testArray0, 0);
      Array1D<double> computedSum1 = axisSum<double>(testArray0, 1);

      BRICK_TEST_ASSERT(
        this->equivalent(computedSum0, sumOverRows,
                         static_cast<double>(m_defaultTolerance)));
      BRICK_TEST_ASSERT(
        this->equivalent(computedSum1, sumOverColumns,
                         static_cast<double>(m_defaultTolerance)));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testAxisSum2()
    {
      // Warning(xxx): Can't make this compile.
//     Array2D<Type> testArray0("[[-2, -1, -10],"
//                              " [5, 2, 3],"
//                              " [10, 0, -2],"
//                              " [1, 2, 2]]");
//     Array1D<double> sumOverRows("[14, 3, -7]");
//     Array1D<double> sumOverColumns("[-13, 10, 8, 5]");

//     Array1D<double> computedSum0 =
//       axisSum(testArray0, 0, Double, std::plus<double>());
//     Array1D<double> computedSum1 =
//       axisSum(testArray0, 1, Double, std::plus<double>());

//     BRICK_TEST_ASSERT(
//       this->equivalent(computedSum0, sumOverRows,
//                        static_cast<double>(m_defaultTolerance)));
//     BRICK_TEST_ASSERT(
//       this->equivalent(computedSum1, sumOverColumns,
//                        static_cast<double>(m_defaultTolerance)));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testColumnIndices()
    {
      Array2D<Type> indices("[[0, 1, 2, 3],"
                            " [0, 1, 2, 3],"
                            " [0, 1, 2, 3]]");
      Array2D<Type> computedIndices =
        columnIndices<Type>(indices.rows(), indices.columns());

      BRICK_TEST_ASSERT(
        this->equivalent(computedIndices, indices,
                         static_cast<Type>(m_defaultTolerance)));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testCompress0()
    {
      Array1D<Type> inputArray("[1, 2, 3, 4, 5]");
      Array1D<bool> mask("[1, 1, 0, 1, 0]");
      Array1D<Type> target("[1, 2, 4]");
      Array1D<Type> computedArray = compress(mask, inputArray);

      BRICK_TEST_ASSERT(
        this->equivalent(computedArray, target,
                         static_cast<Type>(m_defaultTolerance)));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testCompress1()
    {
      Array1D<Type> inputArray("[1, 2, 3, 4, 5]");
      Array1D<bool> mask("[1, 1, 0, 1, 0]");
      Array1D<Type> target("[1, 2, 4]");
      size_t numTrue = 3;
      Array1D<Type> computedArray = compress(mask, inputArray, numTrue);

      BRICK_TEST_ASSERT(
        this->equivalent(computedArray, target,
                         static_cast<Type>(m_defaultTolerance)));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testCount()
    {
      Array1D<Type> inputArray("[1, 2, 0, 4, 0, 5, -1, 0]");
      size_t numTrue = 5;
      BRICK_TEST_ASSERT(count(inputArray) == numTrue);
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testGetCentroid()
    {
      Array1D<Type> inputArray0("[0, 2, 0, 2, 0, 0, 0, 0]");
      double centroid = getCentroid<double>(inputArray0);
      BRICK_TEST_ASSERT(approximatelyEqual(centroid, 2.0, m_defaultTolerance));

      Array1D<Type> inputArray1("[-5, -1, 2, 6, 6, 2, -1, -5]");
      centroid = getCentroid<double>(inputArray1);
      BRICK_TEST_ASSERT(approximatelyEqual(centroid, 3.5, m_defaultTolerance));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testGetMeanAndCovariance()
    {
      const size_t numberOfSamples = 10000;
      const size_t dimensionality = 3;

      // Build samples.
      Array2D<double> sampleArray(numberOfSamples, dimensionality);
      for(size_t sampleIndex = 0; sampleIndex < numberOfSamples; ++sampleIndex) {
        Array1D<double> newSample = this->getRandomSample(dimensionality);
        sampleArray.row(sampleIndex).copy(newSample);
      }

      // Compute sample mean.
      Array1D<double> referenceMean(dimensionality);
      referenceMean = 0.0;
      for(size_t sampleIndex = 0; sampleIndex < numberOfSamples; ++sampleIndex) {
        referenceMean += sampleArray.row(sampleIndex);
      }
      referenceMean /= static_cast<double>(numberOfSamples);

      // Compute sample covariance estimate.
      Array2D<double> referenceCovariance(dimensionality, dimensionality);
      referenceCovariance = 0.0;
      for(size_t sampleIndex = 0; sampleIndex < numberOfSamples; ++sampleIndex) {
        Array1D<double> sampleOffset =
          sampleArray.row(sampleIndex) - referenceMean;
        referenceCovariance += outerProduct<double>(
          sampleOffset, sampleOffset);
      }
      referenceCovariance /= static_cast<double>(numberOfSamples - 1);

      // Run routine under test.
      Array1D<double> testMean;
      Array2D<double> testCovariance;
      getMeanAndCovariance(sampleArray, testMean, testCovariance);

      // Verify result.
      BRICK_TEST_ASSERT(
        this->equivalent(
          testMean, referenceMean, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        this->equivalent(
          testCovariance, referenceCovariance, m_defaultTolerance));

      // Run routine under test for column-major samples.
      getMeanAndCovariance(sampleArray.transpose(), testMean, testCovariance, 1);

      // Verify result.
      BRICK_TEST_ASSERT(
        this->equivalent(
          testMean, referenceMean, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        this->equivalent(
          testCovariance, referenceCovariance, m_defaultTolerance));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testMaximum()
    {
      // Build inputArray..
      Array2D<double> inputArray(
        "[[1.0,  4.0,   2.0,  -3.0, 21.0],"
        " [6.0,  4.0,   7.0,   6.0, -2.0],"
        " [2.0,  5.0, -10.0,  -3.0, 10.0],"
        " [20.0, 4.0,   2.0,   3.0, -2.0]]");

      BRICK_TEST_ASSERT(
        approximatelyEqual(
          maximum(inputArray.getRegion(Index2D(0, 0), Index2D(3, 4))),
          7.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(
          maximum(inputArray.getRegion(Index2D(0, 0), Index2D(3, 5))),
          21.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(
          maximum(inputArray.getRegion(Index2D(0, 0), Index2D(4, 4))),
          20.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(
          maximum(inputArray.getRegion(Index2D(0, 1), Index2D(4, 4))),
          7.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(
          maximum(inputArray.getRegion(Index2D(1, 0), Index2D(4, 4))),
          20.0, m_defaultTolerance));
      BRICK_TEST_ASSERT(
        approximatelyEqual(
          maximum(inputArray.getRegion(Index2D(1, 1), Index2D(4, 5))),
          10.0, m_defaultTolerance));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testNormalizedCorrelation()
    {
      Array1D<Type> inputArray0("[1, 2, 0, 4, 0, 5, -1, 0]");
      Array1D<Type> inputArray1("[7, 9, 2, 4, 0, 3,  6, 7]");

      // Reference computation of normalized correlation.
      double referenceValue;
      Array1D<double> x0(inputArray0.size());
      Array1D<double> y0(inputArray1.size());
      x0.copy(inputArray0);
      y0.copy(inputArray1);
      x0 -= mean<double>(x0);
      y0 -= mean<double>(y0);
      double denominator = ::sqrt(sum<double>(x0 * x0) * sum<double>(y0 * y0));
      if(denominator == static_cast<double>(0)) {
        referenceValue = static_cast<double>(1.0);
      } else {
        referenceValue = sum<double>(x0 * y0) / denominator;
      }

      // Test computation of normalized correlation.
      double testValue = normalizedCorrelation<double>(
        inputArray0, inputArray1);

      BRICK_TEST_ASSERT(
        approximatelyEqual(testValue, referenceValue, m_defaultTolerance));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testSum0()
    {
      typedef typename ArithmeticTraits<Type, Type>::SumType SumType;

      Array2D<Type> testArray0("[[-2, -1, -10],"
                               " [5, 2, 3],"
                               " [10, 0, -2],"
                               " [1, 2, 2]]");
      SumType sum0_2_0_3 = -3;
      SumType sum0_3_0_2 = 14;
      SumType sum2_4_1_3 = 2;

      SumType computedSum0_2_0_3 =
        sum<SumType>(testArray0, Index2D(0, 0), Index2D(2, 3));
      SumType computedSum0_3_0_2 =
        sum<SumType>(testArray0, Index2D(0, 0), Index2D(3, 2));
      SumType computedSum2_4_1_3 =
        sum<SumType>(testArray0, Index2D(2, 1), Index2D(4, 3));

      BRICK_TEST_ASSERT(
        approximatelyEqual(computedSum0_2_0_3, sum0_2_0_3,
                           static_cast<SumType>(m_defaultTolerance)));
      BRICK_TEST_ASSERT(
        approximatelyEqual(computedSum0_3_0_2, sum0_3_0_2,
                           static_cast<SumType>(m_defaultTolerance)));
      BRICK_TEST_ASSERT(
        approximatelyEqual(computedSum2_4_1_3, sum2_4_1_3,
                           static_cast<SumType>(m_defaultTolerance)));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testTake_Array1D_Array1D()
    {
      Array1D<Type> dataArray("[0, 1, 2, 3, 4, 5]");
      Array1D<unsigned int> indexArray("[2, 1, 1, 4]");
      Array1D<Type> resultArray = take(dataArray, indexArray);

      Array1D<Type> referenceArray("[2, 1, 1, 4]");
      BRICK_TEST_ASSERT(this->equivalent(resultArray, referenceArray,
                                         static_cast<Type>(0)));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testTake_Array2D_Array1D()
    {
      Array2D<Type> dataArray("[[0, 1, 2], [3, 4, 5]]");
      Array1D<unsigned int> indexArray("[2, 1, 1, 4]");
      Array1D<Type> resultArray = take(dataArray, indexArray);

      Array1D<Type> referenceArray("[2, 1, 1, 4]");
      BRICK_TEST_ASSERT(this->equivalent(resultArray, referenceArray,
                                         static_cast<Type>(0)));
    }


    template <class Type>
    void
    UtilitiesTest<Type>::
    testTake_Array2D_Array1D_unsignedInt()
    {
      Array2D<Type> dataArray("[[0, 1, 2], "
                              " [3, 4, 5], "
                              " [6, 7, 8], "
                              " [9, 10, 11]]");
      Array1D<unsigned int> indexArray("[0, 2]");
      Array2D<Type> resultArray0 = take(dataArray, indexArray, 0);
      Array2D<Type> resultArray1 = take(dataArray, indexArray, 1);

      Array2D<Type> referenceArray0("[[0, 1, 2],"
                                    " [6, 7, 8]]");
      Array2D<Type> referenceArray1("[[0, 2],"
                                    " [3, 5],"
                                    " [6, 8],"
                                    " [9, 11]]");
      BRICK_TEST_ASSERT(this->equivalent(resultArray0, referenceArray0,
                                         static_cast<Type>(0)));
      BRICK_TEST_ASSERT(this->equivalent(resultArray1, referenceArray1,
                                         static_cast<Type>(0)));
    }

    template <class Type>
    template <class Type2>
    bool
    UtilitiesTest<Type>::
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


    template <class Type>
    template <class Type2>
    bool
    UtilitiesTest<Type>::
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


    template<class Type>
    Array1D<double>
    UtilitiesTest<Type>::
    getRandomSample(size_t dimensionality)
    {
      Array1D<double> returnValue(dimensionality);
      double randMax = static_cast<double>(RAND_MAX);
      for(size_t elementIndex = 0; elementIndex < dimensionality;
          ++elementIndex) {
        returnValue[elementIndex] = rand() / randMax;
      }
      return returnValue;
    }

  } // namespace numeric

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::numeric::UtilitiesTest<double> currentTest0("double");
  brick::numeric::UtilitiesTest<float> currentTest1("float");
  brick::numeric::UtilitiesTest<int> currentTest2("int");
  bool result = (currentTest0.run()
                 && currentTest1.run()
                 && currentTest2.run());
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::UtilitiesTest<double> currentTest0("double");
  brick::numeric::UtilitiesTest<float> currentTest1("float");
  brick::numeric::UtilitiesTest<int> currentTest2("int");

}

#endif
