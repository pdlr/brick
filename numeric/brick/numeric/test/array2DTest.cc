/**
***************************************************************************
* @file brick/numeric/test/array2DTest.cpp
*
* Source file defining Array2DTest class.
*
* Copyright (C) 2004-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <math.h>
#include <iomanip>
#include <sstream>
// #include <brick/numeric/utilities.hh>
#include <brick/numeric/array2D.hh>
#include <brick/numeric/test/arrayTestCommon.hh>
#include <brick/test/functors.hh>

namespace brick {

  namespace numeric {

    template <class Type>
    class Array2DTest
      : public ArrayTestCommon< Array2DTest<Type>, Array2D<Type>,
                                Array2D<bool> >
    {

    public:

      typedef Array2DTest<Type> TestFixtureType;


      Array2DTest(const std::string& typeName);
      virtual ~Array2DTest();

      // Inherited methods from TestFixture
      void setUp(const std::string&) {}
      void tearDown(const std::string&) {}

      // Inherited methods from ArrayTestCommon
      virtual void
      checkShapeEquality(const Array2D<Type>& array0,
                         const Array2D<Type>& array1);


      virtual void
      checkShapeEquality(const Array2D<Type>& array0,
                         const Array2D<bool>& array1);


      virtual void
      checkValueEquality(const Array2D<Type>& array0,
                         const Array2D<Type>& array1,
                         typename Array2D<Type>::value_type tolerance);


      virtual typename Array2D<Type>::value_type
      getComparisonOperatorThreshold() {
        return m_fibonacciCArray[m_defaultArraySize / 2];
      }


      virtual typename Array2D<Type>::value_type
      getEqualityOperatorTarget() {
        return m_squaresCArray[m_defaultArraySize / 2];
      }


      virtual Array2D<Type>
      getFibonacciArray() {return Array2D<Type>(m_fibonacciString);}


      virtual typename Array2D<Type>::value_type
      getIncrementOperatorArgument() {return static_cast<Type>(4);}


      virtual typename Array2D<Type>::value_type
      getMultiplicationOperatorArgument() {return static_cast<Type>(2);}


      virtual Array2D<Type>
      getSquaresArray() {return Array2D<Type>(m_squaresString);}


      // Tests of member functions.
      void testConstructor__void();
      void testConstructor__size_t__size_t();
      void testConstructor__size_t__size_t__size_t();
      void testConstructor__string();
      void testConstructor__Array2D();
      void testConstructor__size_t__size_t__TypePtr();
      void testDestructor();
      void testBegin();
      void testBeginConst();
      void testClear();
      void testColumns();
      void testCopy();
      void testCopy__Array2D();
      void testCopy__Type2Ptr();
      void testData();
      void testDataConst();
      void testData__size_t();
      void testDataConst__size_t();
      void testData__size_t__size_t();
      void testDataConst__size_t__size_t();
      void testEnd();
      void testEndConst();
      void testGetRegion();
      void testRavel();
      void testRavelConst();
      void testReadFromStream();
      void testReinit();
      void testReshape();
      void testRow();
      void testRowConst();
      void testRows();
      void testShape();
      void testShape__size_t();
      void testSize();
      void testTranspose();
      void testAssignmentOperator__Array2D();
      void testAssignmentOperator__Type();
      void testApplicationOperator__size_t();
      void testApplicationOperatorConst__size_t();
      void testApplicationOperator__size_t__size_t();
      void testApplicationOperatorConst__size_t__size_t();
      void testIndexOperator();
      void testIndexOperatorConst();

      // C++11 tests.
      void testInitializerList();

      // Tests of non-member functions.
      void testSquareRoot__Array2D();
      void testSqrt__Array2D();


    private:

      size_t m_defaultArrayColumns;
      size_t m_defaultArrayRows;
      size_t m_defaultArraySize;
      size_t m_defaultArrayValue;
      size_t m_defaultIncrement;
      size_t m_defaultMultiplier;
      Type* m_fibonacciCArray;
      std::string m_fibonacciString;
      std::string m_illegalString;
      Type* m_squaresCArray;
      std::string m_squaresString;
      double m_testEpsilon;

    }; // class Array2DTest


    /* ============== Member Function Definititions ============== */

    template <class Type>
    Array2DTest<Type>::
    Array2DTest(const std::string& typeName)
      : ArrayTestCommon< Array2DTest<Type>, Array2D<Type>, Array2D<bool> >(
        std::string("Array2DTest<") + typeName + ">"),
        m_defaultArrayColumns(4),
        m_defaultArrayRows(3),
        m_defaultArraySize(12),
        m_defaultArrayValue(9),
        m_defaultIncrement(4),
        m_defaultMultiplier(2),
        m_fibonacciCArray(0),
        m_fibonacciString(""),
        m_illegalString(),
        m_squaresCArray(0),
        m_squaresString(""),
        m_testEpsilon(1.0e-8)
    {
      // // Register all tests.
      // Tests of member functions.
      BRICK_TEST_REGISTER_MEMBER(testConstructor__void);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__size_t__size_t);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__size_t__size_t__size_t);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__string);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__Array2D);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__size_t__size_t__TypePtr);
      BRICK_TEST_REGISTER_MEMBER(testDestructor);
      BRICK_TEST_REGISTER_MEMBER(testBegin);
      BRICK_TEST_REGISTER_MEMBER(testBeginConst);
      BRICK_TEST_REGISTER_MEMBER(testClear);
      BRICK_TEST_REGISTER_MEMBER(testColumns);
      BRICK_TEST_REGISTER_MEMBER(testCopy);
      BRICK_TEST_REGISTER_MEMBER(testCopy__Array2D);
      BRICK_TEST_REGISTER_MEMBER(testCopy__Type2Ptr);
      BRICK_TEST_REGISTER_MEMBER(testData);
      BRICK_TEST_REGISTER_MEMBER(testDataConst);
      BRICK_TEST_REGISTER_MEMBER(testData__size_t);
      BRICK_TEST_REGISTER_MEMBER(testDataConst__size_t);
      BRICK_TEST_REGISTER_MEMBER(testData__size_t__size_t);
      BRICK_TEST_REGISTER_MEMBER(testDataConst__size_t__size_t);
      BRICK_TEST_REGISTER_MEMBER(testEnd);
      BRICK_TEST_REGISTER_MEMBER(testEndConst);
      BRICK_TEST_REGISTER_MEMBER(testGetRegion);
      BRICK_TEST_REGISTER_MEMBER(testRavel);
      BRICK_TEST_REGISTER_MEMBER(testRavelConst);
      BRICK_TEST_REGISTER_MEMBER(testReadFromStream);
      BRICK_TEST_REGISTER_MEMBER(testReinit);
      BRICK_TEST_REGISTER_MEMBER(testReshape);
      BRICK_TEST_REGISTER_MEMBER(testRow);
      BRICK_TEST_REGISTER_MEMBER(testRowConst);
      BRICK_TEST_REGISTER_MEMBER(testRows);
      BRICK_TEST_REGISTER_MEMBER(testShape);
      BRICK_TEST_REGISTER_MEMBER(testShape__size_t);
      BRICK_TEST_REGISTER_MEMBER(testSize);
      BRICK_TEST_REGISTER_MEMBER(testTranspose);
      BRICK_TEST_REGISTER_MEMBER(testAssignmentOperator__Array2D);
      BRICK_TEST_REGISTER_MEMBER(testAssignmentOperator__Type);
      BRICK_TEST_REGISTER_MEMBER(testApplicationOperator__size_t);
      BRICK_TEST_REGISTER_MEMBER(testApplicationOperatorConst__size_t);
      BRICK_TEST_REGISTER_MEMBER(testApplicationOperator__size_t__size_t);
      BRICK_TEST_REGISTER_MEMBER(testApplicationOperatorConst__size_t__size_t);
      BRICK_TEST_REGISTER_MEMBER(testIndexOperator);
      BRICK_TEST_REGISTER_MEMBER(testIndexOperatorConst);
      BRICK_TEST_REGISTER_MEMBER(testInitializerList);

      // Tests of non-member functions.
      BRICK_TEST_REGISTER_MEMBER(testSquareRoot__Array2D);
      BRICK_TEST_REGISTER_MEMBER(testSqrt__Array2D);


      // Set up fibonacci data for tests.
      m_fibonacciCArray = new Type[m_defaultArraySize];
      m_fibonacciCArray[0] = static_cast<Type>(1);
      m_fibonacciCArray[1] = static_cast<Type>(2);
      for(size_t index = 2; index < m_defaultArraySize; ++index) {
        m_fibonacciCArray[index] = (m_fibonacciCArray[index - 1]
                                    + m_fibonacciCArray[index - 2]);
      }

      // Set up squares data for tests.
      m_squaresCArray = new Type[m_defaultArraySize];
      for(size_t index = 0; index < m_defaultArraySize; ++index) {
        m_squaresCArray[index] = static_cast<Type>(index * index);
      }

      // Set up input strings for tests.
      std::ostringstream fibonacciBuffer;
      std::ostringstream squaresBuffer;
      size_t index0 = 0;
      fibonacciBuffer << std::fixed << "[";
      squaresBuffer << std::fixed << "[";
      for(size_t rowIndex = 0; rowIndex < m_defaultArrayRows;
          ++rowIndex) {
        fibonacciBuffer << "[" << static_cast<Type>(m_fibonacciCArray[index0]);
        squaresBuffer << "[" << static_cast<Type>(m_squaresCArray[index0]);
        ++index0;
        for(size_t column = 1; column < m_defaultArrayColumns; ++column) {
          fibonacciBuffer << ", "
                          << static_cast<Type>(m_fibonacciCArray[index0]);
          squaresBuffer << ", "
                        << static_cast<Type>(m_squaresCArray[index0]);
          ++index0;
        }
        fibonacciBuffer << "]";
        squaresBuffer << "]";
        if(rowIndex != m_defaultArrayRows - 1) {
          fibonacciBuffer << ",\n";
          squaresBuffer << ",\n";
        }
      }
      fibonacciBuffer << "]";
      squaresBuffer << "]";
      m_fibonacciString = fibonacciBuffer.str();
      m_squaresString = squaresBuffer.str();

      // Set up malformed input string for tests.
      m_illegalString = "[[1, 2, 3 4, 5, 6]\n[1, 2, 3, 4, 5, 6]]";
    }


    template <class Type>
    Array2DTest<Type>::
    ~Array2DTest()
    {
      delete[] m_fibonacciCArray;
      delete[] m_squaresCArray;
    }


    template <class Type>
    void
    Array2DTest<Type>::
    checkShapeEquality(const Array2D<Type>& array0,
                       const Array2D<Type>& array1)
    {
      BRICK_TEST_ASSERT(array0.rows() == array1.rows());
      BRICK_TEST_ASSERT(array0.columns() == array1.columns());
      BRICK_TEST_ASSERT(array0.size() == array1.size());
      BRICK_TEST_ASSERT(array0.data() != array1.data());
    }


    template <class Type>
    void
    Array2DTest<Type>::
    checkShapeEquality(const Array2D<Type>& array0,
                       const Array2D<bool>& array1)
    {
      BRICK_TEST_ASSERT(array0.rows() == array1.rows());
      BRICK_TEST_ASSERT(array0.columns() == array1.columns());
      BRICK_TEST_ASSERT(array0.size() == array1.size());
    }


    template <class Type>
    void
    Array2DTest<Type>::
    checkValueEquality(const Array2D<Type>& array0,
                       const Array2D<Type>& array1,
                       typename Array2D<Type>::value_type tolerance)
    {
      this->checkShapeEquality(array0, array1);
      BRICK_TEST_ASSERT(std::equal(array0.begin(), array0.end(), array1.begin(),
                                   ApproximatelyEqualFunctor<Type>(tolerance)));
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testConstructor__void()
    {
      // Default constructor should initialize to zero size.
      Array2D<Type> array0;
      BRICK_TEST_ASSERT(array0.rows() == 0);
      BRICK_TEST_ASSERT(array0.columns() == 0);
      BRICK_TEST_ASSERT(array0.size() == 0);
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testConstructor__size_t__size_t()
    {
      // This constructor should initialize to the specified size.
      Array2D<Type> array0(m_defaultArrayRows, m_defaultArrayColumns);
      BRICK_TEST_ASSERT(array0.rows() == m_defaultArrayRows);
      BRICK_TEST_ASSERT(array0.columns() == m_defaultArrayColumns);
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize);
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testConstructor__size_t__size_t__size_t()
    {
      // This constructor should initialize to the specified size.
      Array2D<Type> array0(m_defaultArrayRows, m_defaultArrayColumns,
                           m_defaultArrayColumns + 2);
      BRICK_TEST_ASSERT(array0.rows() == m_defaultArrayRows);
      BRICK_TEST_ASSERT(array0.columns() == m_defaultArrayColumns);
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize);

      array0 = Type(1);
      array0(1, 2) = Type(2);
      array0(2, 1) = Type(3);
      unsigned int rowStep = m_defaultArrayColumns + 2;
      for(unsigned int rr = 0; rr < m_defaultArrayRows; ++rr) {
        for(unsigned int cc = 0; cc < m_defaultArrayColumns; ++cc) {
          BRICK_TEST_ASSERT(
            &(array0(rr, cc)) == array0.getData() + cc + rr * rowStep);
          if((1 == rr) && (2 == cc)) {
            BRICK_TEST_ASSERT(2 == array0(rr, cc))
          } else if((2 == rr) && (1 == cc)) {
            BRICK_TEST_ASSERT(3 == array0(rr, cc))
          } else {
            BRICK_TEST_ASSERT(1 == array0(rr, cc))
          }
        }
      }
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testConstructor__string()
    {
      // This constructor should initialize according to the specified
      // string, which describes an array of m_defaultArraySize elements
      // having values matching those in m_fibonacciCArray.
      Array2D<Type> array0(m_fibonacciString);
      BRICK_TEST_ASSERT(array0.rows() == m_defaultArrayRows);
      BRICK_TEST_ASSERT(array0.columns() == m_defaultArrayColumns);
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize);
      BRICK_TEST_ASSERT(std::equal(array0.begin(), array0.end(),
                                   m_fibonacciCArray));

      // Construct using a malformed string.
      BRICK_TEST_ASSERT_EXCEPTION(common::ValueException,
                                  Array2D<Type> array1(m_illegalString));

    }


    template <class Type>
    void
    Array2DTest<Type>::
    testConstructor__Array2D()
    {
      // This constructor should make an array which references the same
      // data as its argument.
      Array2D<Type>* array0Ptr = new Array2D<Type>(m_fibonacciString);
      Array2D<Type> array1(*array0Ptr);
      BRICK_TEST_ASSERT(array1.data() == array0Ptr->data());
      BRICK_TEST_ASSERT(array1.rows() == array0Ptr->rows());
      BRICK_TEST_ASSERT(array1.columns() == array0Ptr->columns());
      BRICK_TEST_ASSERT(array1.size() == array0Ptr->size());
      delete array0Ptr;

      // Normally, deleting array0Ptr would invalidate the data it
      // points to.  If reference counting is working, however, the data
      // should remain valid until after array1 goes out of scope.
      BRICK_TEST_ASSERT(std::equal(array1.begin(), array1.end(),
                                   m_fibonacciCArray));
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testConstructor__size_t__size_t__TypePtr()
    {
      // This constructor initializes the array to point to the
      // specified data.
      // Start by making an array we can play with.
      Type* cArray = new Type[m_defaultArraySize];
      std::copy(m_fibonacciCArray, m_fibonacciCArray + m_defaultArraySize,
                cArray);

      // Create an array using the constructor under test.
      Array2D<Type>* array0Ptr = new Array2D<Type>(
        m_defaultArrayRows, m_defaultArrayColumns, cArray);
      BRICK_TEST_ASSERT(array0Ptr->data() == cArray);
      BRICK_TEST_ASSERT(array0Ptr->rows() == m_defaultArrayRows);
      BRICK_TEST_ASSERT(array0Ptr->columns() == m_defaultArrayColumns);
      BRICK_TEST_ASSERT(array0Ptr->size() == m_defaultArraySize);
      delete array0Ptr;

      // In this case, no reference counting should be done, and the
      // data should not be deleted on deletion of array0Ptr.
      BRICK_TEST_ASSERT(std::equal(cArray, cArray + m_defaultArraySize,
                                   m_fibonacciCArray));

      // Clean up.
      delete[] cArray;
    }


//   testConstructor__size_t__TypePtr__size_tPtr()
//   {
//     // This constructor should reference the provided test data _and_
//     // maintain the reference count.

//     // First set up some data to reference.
//     Array2D<Type> array0(m_fibonacciString);

//     // Create a size_t to hold the reference count.
//     size_t refCount = 1;

//     // Create the array.  It should now reference the same data as the
//     // first array, and the reference count should be incremented.
//     Array2D<Type>* array2Ptr = new Array2D<Type>(
//       array0.size(), array0.data(), &refCount);
//     BRICK_TEST_ASSERT(array2Ptr->data() == array0.data());
//     BRICK_TEST_ASSERT(array2Ptr->size() == array0.size());
//     BRICK_TEST_ASSERT(refCount == 2);
//     BRICK_TEST_ASSERT(std::equal(array2Ptr->begin(), array2Ptr->end(),
//                                m_fibonacciCArray));

//     // Following deletion of the new array, the reference count should
//     // be back down to 1.
//     delete array2Ptr;
//     BRICK_TEST_ASSERT(refCount == 1);
//   }



    template <class Type>
    void
    Array2DTest<Type>::
    testDestructor()
    {
      // No independent test for destructor.
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testBegin()
    {
      // Member function begin() should return an iterator pointing to
      // the first element of the array.
      Array2D<Type> array0(
        m_defaultArrayRows, m_defaultArrayColumns, m_fibonacciCArray);
      BRICK_TEST_ASSERT(&(*(array0.begin())) == m_fibonacciCArray);
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testBeginConst()
    {
      // Member function begin() should return an iterator pointing to
      // the first element of the array.
      const Array2D<Type> array0(
        m_defaultArrayRows, m_defaultArrayColumns, m_fibonacciCArray);
      BRICK_TEST_ASSERT(&(*(array0.begin())) == m_fibonacciCArray);
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testClear()
    {
      // Member function clear should reset the array to zero size,
      // abandoning all contents.
      // First set up some data to reference.
      Array2D<Type> array0(m_fibonacciString);

      // Clearing the array should release the reference.
      array0.clear();
      BRICK_TEST_ASSERT(array0.size() == 0);

//     // Member function clear should reset the array to zero size,
//     // abandoning all contents.
//     // First set up some data to reference.
//     Array2D<Type> array0(m_fibonacciString);

//     // Create a size_t to hold a reference count.
//     size_t refCount = 1;

//     // Create the test array.  It should now reference the same data
//     // as the first array, and the reference count should be
//     // incremented.
//     Array2D<Type>array1(array0.size(), array0.data(), &refCount);
//     BRICK_TEST_ASSERT(refCount == 2);

//     // Clearing the array should release the reference.
//     array1.clear();
//     BRICK_TEST_ASSERT(array1.size() == 0);
//     BRICK_TEST_ASSERT(refCount == 1);
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testColumns()
    {
      // No explicit test.
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testCopy()
    {
      // Member function copy() should allocate a new array and deep
      // copy the contents of *this.
      Array2D<Type> array0(m_fibonacciString);
      Array2D<Type> array1 = array0.copy();
      BRICK_TEST_ASSERT(array1.size() == array0.size());
      BRICK_TEST_ASSERT(array1.rows() == array0.rows());
      BRICK_TEST_ASSERT(array1.columns() == array0.columns());
      BRICK_TEST_ASSERT(std::equal(array1.begin(), array1.end(), array0.begin()));
      BRICK_TEST_ASSERT(array1.data() != array0.data());
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testCopy__Array2D()
    {
      // This member function deep copies the contents of source.  It is
      // an error if source does not have the same size as *this.
      Array2D<Type> array0(m_fibonacciString);
      Array2D<Type> array1(array0.rows(), array0.columns());
      Array2D<Type> array2(array0.rows(), array0.columns() + 1);
      Array2D<Type> array3(array0.rows() + 1, array0.columns());

      // Test that the copy goes OK.
      array1.copy(array0);
      BRICK_TEST_ASSERT(array1.rows() == array0.rows());
      BRICK_TEST_ASSERT(array1.columns() == array0.columns());
      BRICK_TEST_ASSERT(array1.size() == array0.size());
      BRICK_TEST_ASSERT(std::equal(array1.begin(), array1.end(), array0.begin()));
      BRICK_TEST_ASSERT(array1.data() != array0.data());

      // Test that size is checked.
      BRICK_TEST_ASSERT_EXCEPTION(common::ValueException, array2.copy(array0));
      BRICK_TEST_ASSERT_EXCEPTION(common::ValueException, array3.copy(array0));
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testCopy__Type2Ptr()
    {
      // This method simply copies elements from its argument.
      // Begin by making an appropriately sized array with known contents.
      Array2D<Type> array0(m_defaultArrayRows, m_defaultArrayColumns);
      array0 = 0;

      array0.copy(m_fibonacciCArray);
      BRICK_TEST_ASSERT(array0.rows() == m_defaultArrayRows);
      BRICK_TEST_ASSERT(array0.columns() == m_defaultArrayColumns);
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize);
      BRICK_TEST_ASSERT(std::equal(array0.begin(), array0.end(),
                                   m_fibonacciCArray));
      BRICK_TEST_ASSERT(array0.data() != m_fibonacciCArray);
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testData()
    {
      // This method returns a pointer to the internal data store.
      Array2D<Type> array0(
        m_defaultArrayRows, m_defaultArrayColumns, m_fibonacciCArray);
      BRICK_TEST_ASSERT(array0.data() == m_fibonacciCArray);
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testDataConst()
    {
      // This method returns a pointer to the internal data store.
      const Array2D<Type> array0(
        m_defaultArrayRows, m_defaultArrayColumns, m_fibonacciCArray);
      BRICK_TEST_ASSERT(array0.data() == m_fibonacciCArray);
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testData__size_t()
    {
      // This method returns a pointer to the index-th element of the
      // internal data store.
      Array2D<Type> array0(
        m_defaultArrayRows, m_defaultArrayColumns, m_fibonacciCArray);
      BRICK_TEST_ASSERT(array0.data(m_defaultArraySize - 1)
                        == (m_fibonacciCArray + m_defaultArraySize - 1));
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testDataConst__size_t()
    {
      // This method returns a pointer to the index-th element of the
      // internal data store.
      Array2D<Type> array0(
        m_defaultArrayRows, m_defaultArrayColumns, m_fibonacciCArray);
      BRICK_TEST_ASSERT(array0.data(m_defaultArraySize - 1)
                        == (m_fibonacciCArray + m_defaultArraySize - 1));
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testData__size_t__size_t()
    {
      // This method returns a pointer to the index-th element of the
      // internal data store.
      Array2D<Type> array0(
        m_defaultArrayRows, m_defaultArrayColumns, m_fibonacciCArray);
      BRICK_TEST_ASSERT(array0.data(0, 0) == m_fibonacciCArray);
      BRICK_TEST_ASSERT(array0.data(0, 1) == (m_fibonacciCArray + 1));
      BRICK_TEST_ASSERT(array0.data(1, 0)
                        == (m_fibonacciCArray + m_defaultArrayColumns));
      BRICK_TEST_ASSERT(array0.data(m_defaultArrayRows - 1,
                                    m_defaultArrayColumns - 1)
                        == (m_fibonacciCArray + m_defaultArraySize - 1));
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testDataConst__size_t__size_t()
    {
      // This method returns a pointer to the index-th element of the
      // internal data store.
      const Array2D<Type> array0(
        m_defaultArrayRows, m_defaultArrayColumns, m_fibonacciCArray);
      BRICK_TEST_ASSERT(array0.data(0, 0) == m_fibonacciCArray);
      BRICK_TEST_ASSERT(array0.data(0, 1) == (m_fibonacciCArray + 1));
      BRICK_TEST_ASSERT(array0.data(1, 0)
                        == (m_fibonacciCArray + m_defaultArrayColumns));
      BRICK_TEST_ASSERT(array0.data(m_defaultArrayRows - 1,
                                    m_defaultArrayColumns - 1)
                        == (m_fibonacciCArray + m_defaultArraySize - 1));
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testEnd()
    {
      // Member function end() should return an iterator pointing to
      // the last element of the array.
      Array2D<Type> array0(
        m_defaultArrayRows, m_defaultArrayColumns, m_fibonacciCArray);
      Type* finalElementPtr = m_fibonacciCArray + m_defaultArraySize;
      BRICK_TEST_ASSERT(&(*(array0.end())) == finalElementPtr);
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testEndConst()
    {
      // Member function end() should return an iterator pointing to
      // the last element of the array.
      const Array2D<Type> array0(
        m_defaultArrayRows, m_defaultArrayColumns, m_fibonacciCArray);
      Type* finalElementPtr = m_fibonacciCArray + m_defaultArraySize;
      BRICK_TEST_ASSERT(&(*(array0.end())) == finalElementPtr);
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testGetRegion()
    {
      // Create an array in which every element is set to 1.
      Array2D<Type> fullArray(10, 20);
      fullArray = Type(1);

      // Get an array that refers to a 4x4 subset of fullArray.
      Index2D corner0(5, 10);
      Index2D corner1(9, 14);
      Array2D<Type> subRegion = fullArray.getRegion(corner0, corner1);

      // Set just the elements within that subset to 100.
      subRegion = Type(100);

      // And set one more element to 50.
      subRegion(1, 2) = Type(50);

      // Check that all is as expected.
      for(unsigned int rr = 0; rr < fullArray.rows(); ++rr) {
        for(unsigned int cc = 0; cc < fullArray.columns(); ++cc) {

          if((static_cast<int>(rr) == corner0.getRow()    + 1) &&
             (static_cast<int>(cc) == corner0.getColumn() + 2)) {

            BRICK_TEST_ASSERT(Type(50) == fullArray(rr, cc));

          } else if((static_cast<int>(rr) >= corner0.getRow()) &&
                    (static_cast<int>(rr) <  corner1.getRow()) &&
                    (static_cast<int>(cc) >= corner0.getColumn()) &&
                    (static_cast<int>(cc) <  corner1.getColumn())) {

            BRICK_TEST_ASSERT(Type(100) == fullArray(rr, cc));

          } else {

            BRICK_TEST_ASSERT(Type(1) == fullArray(rr, cc));

          }
        }
      }
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testRavel()
    {
      // Member function ravel() should return an Array1D instance referencing
      // the same data.
      Array2D<Type>* array0Ptr = new Array2D<Type>(
        m_defaultArrayRows, m_defaultArrayColumns, m_fibonacciCArray);
      Array1D<Type> array1 = array0Ptr->ravel();

      // Check for correct value.
      BRICK_TEST_ASSERT(array1.size() == array0Ptr->size());
      BRICK_TEST_ASSERT(array1.data() == array0Ptr->data());
      BRICK_TEST_ASSERT(std::equal(array1.begin(), array1.end(),
                                   m_fibonacciCArray));

      // Delete array0Ptr and hope that array1 is still valid.
      delete array0Ptr;
      BRICK_TEST_ASSERT(std::equal(array1.begin(), array1.end(),
                                   m_fibonacciCArray));
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testRavelConst()
    {
      // Member function ravel() should return an Array1D instance referencing
      // the same data.
      const Array2D<Type>* array0Ptr = new Array2D<Type>(
        m_defaultArrayRows, m_defaultArrayColumns, m_fibonacciCArray);
      const Array1D<Type> array1 = array0Ptr->ravel();

      // Check for correct value.
      BRICK_TEST_ASSERT(array1.size() == array0Ptr->size());
      BRICK_TEST_ASSERT(array1.data() == array0Ptr->data());
      BRICK_TEST_ASSERT(std::equal(array1.begin(), array1.end(),
                                   m_fibonacciCArray));

      // Delete array0Ptr and hope that array1 is still valid.
      delete array0Ptr;
      BRICK_TEST_ASSERT(std::equal(array1.begin(), array1.end(),
                                   m_fibonacciCArray));
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testReadFromStream()
    {
      // This member function sets the value of the array from an input
      // stream.
      // Set up input streams
      std::istringstream inputStream0(m_fibonacciString);
      std::istringstream inputStream1(m_illegalString);

      // Try to parse the good stream.
      Array2D<Type> array0;
      array0.readFromStream(inputStream0);
      BRICK_TEST_ASSERT(inputStream0);
      BRICK_TEST_ASSERT(array0.rows() == m_defaultArrayRows);
      BRICK_TEST_ASSERT(array0.columns() == m_defaultArrayColumns);
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize);
      BRICK_TEST_ASSERT(std::equal(array0.begin(), array0.end(),
                                   m_fibonacciCArray));

      // Try to parse the bad stream.  This should not affect array0 in
      // any way, but it should change the stream state.
      array0.readFromStream(inputStream1);
      BRICK_TEST_ASSERT(!inputStream1);
      BRICK_TEST_ASSERT(array0.rows() == m_defaultArrayRows);
      BRICK_TEST_ASSERT(array0.columns() == m_defaultArrayColumns);
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize);
      BRICK_TEST_ASSERT(std::equal(array0.begin(), array0.end(),
                                   m_fibonacciCArray));
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testReinit()
    {
      // This method changes the shape of the array and reallocates
      // storage.
      Array2D<Type> array0(
        m_defaultArrayRows, m_defaultArrayColumns, m_fibonacciCArray);
      array0.reinit(m_defaultArrayRows + 1, m_defaultArrayColumns + 2);
      BRICK_TEST_ASSERT(array0.rows() == m_defaultArrayRows + 1);
      BRICK_TEST_ASSERT(array0.columns() == m_defaultArrayColumns + 2);
      BRICK_TEST_ASSERT(array0.size() ==
                        (m_defaultArrayRows + 1) * (m_defaultArrayColumns + 2));
      BRICK_TEST_ASSERT(array0.data() != m_fibonacciCArray);
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testReshape()
    {
      // This method changes the shape of the array withot reallocating
      // storage.
      Array2D<Type> array0(
        m_defaultArrayRows, m_defaultArrayColumns, m_fibonacciCArray);

      // Test that size is checked.
      BRICK_TEST_ASSERT_EXCEPTION(
        common::ValueException,
        array0.reshape(static_cast<int>(m_defaultArrayRows) + 1,
                       static_cast<int>(m_defaultArrayColumns) + 2));

      // Check that correctly sized reshapes work.
      array0.reshape(1, static_cast<int>(m_defaultArraySize));
      BRICK_TEST_ASSERT(array0.rows() == 1);
      BRICK_TEST_ASSERT(array0.columns() == m_defaultArraySize);
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize);
      BRICK_TEST_ASSERT(array0.data() == m_fibonacciCArray);
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testRow()
    {
      Array2D<Type> array0(m_fibonacciString);
      Type* dataPtr = array0.data();
      Array1D<Type> array1 = array0.row(1);

      BRICK_TEST_ASSERT(array0.data() == dataPtr);
      BRICK_TEST_ASSERT(array1.size() == m_defaultArrayColumns);
      BRICK_TEST_ASSERT(array1.data() == dataPtr + m_defaultArrayColumns);
      for(size_t index0 = 0; index0 < m_defaultArrayColumns; ++index0) {
        BRICK_TEST_ASSERT(array1[index0] ==
                          m_fibonacciCArray[m_defaultArrayColumns + index0]);
      }
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testRowConst()
    {
      const Array2D<Type> array0(m_fibonacciString);
      const Type* dataPtr = array0.data();
      const Array1D<Type> array1 = array0.row(1);

      BRICK_TEST_ASSERT(array0.data() == dataPtr);
      BRICK_TEST_ASSERT(array1.size() == m_defaultArrayColumns);
      BRICK_TEST_ASSERT(array1.data() == dataPtr + m_defaultArrayColumns);
      for(size_t index0 = 0; index0 < m_defaultArrayColumns; ++index0) {
        BRICK_TEST_ASSERT(array1[index0] ==
                          m_fibonacciCArray[m_defaultArrayColumns + index0]);
      }
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testRows()
    {
      // No explicit test.
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testShape()
    {
      Array2D<Type> array0(
        m_defaultArrayRows, m_defaultArrayColumns, m_fibonacciCArray);
      Array1D<size_t> shapeArray = array0.shape();
      BRICK_TEST_ASSERT(shapeArray.size() == 2);
      BRICK_TEST_ASSERT(shapeArray(0) == m_defaultArrayRows);
      BRICK_TEST_ASSERT(shapeArray(1) == m_defaultArrayColumns);
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testShape__size_t()
    {
      Array2D<Type> array0(
        m_defaultArrayRows, m_defaultArrayColumns, m_fibonacciCArray);
      BRICK_TEST_ASSERT(array0.shape(0) == m_defaultArrayRows);
      BRICK_TEST_ASSERT(array0.shape(1) == m_defaultArrayColumns);
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testSize()
    {
      // No explicit test.
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testTranspose()
    {
      const Array2D<Type> array0(m_fibonacciString);
      const Array2D<Type> array1 = array0.transpose();

      BRICK_TEST_ASSERT(array0.rows() == m_defaultArrayRows);
      BRICK_TEST_ASSERT(array0.columns() == m_defaultArrayColumns);
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize);
      BRICK_TEST_ASSERT(array1.rows() == m_defaultArrayColumns);
      BRICK_TEST_ASSERT(array1.columns() == m_defaultArrayRows);
      BRICK_TEST_ASSERT(array1.size() == m_defaultArraySize);
      size_t index0 = 0;
      for(size_t row = 0; row < m_defaultArrayRows;
          ++row) {
        for(size_t column = 0; column < m_defaultArrayColumns; ++column) {
          BRICK_TEST_ASSERT(array0(row, column) == m_fibonacciCArray[index0]);
          BRICK_TEST_ASSERT(array1(column, row) == m_fibonacciCArray[index0]);
          ++index0;
        }
      }
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testAssignmentOperator__Array2D()
    {
      // The assignment operator shallow copies the contents of source.
      Array2D<Type> array0(m_fibonacciString);
      Array2D<Type> array1;
      array1 = array0;
      BRICK_TEST_ASSERT(array0.rows() == m_defaultArrayRows);
      BRICK_TEST_ASSERT(array0.columns() == m_defaultArrayColumns);
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize);
      BRICK_TEST_ASSERT(array1.rows() == m_defaultArrayRows);
      BRICK_TEST_ASSERT(array1.columns() == m_defaultArrayColumns);
      BRICK_TEST_ASSERT(array1.size() == array0.size());
      BRICK_TEST_ASSERT(array1.data() == array0.data());

      // Make sure the data wasn't changed by the operation.
      BRICK_TEST_ASSERT(std::equal(array1.begin(), array1.end(),
                                   m_fibonacciCArray));
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testAssignmentOperator__Type()
    {
      // This member function assigns the same value to every element in
      // the array.
      Array2D<Type> array0(m_defaultArrayRows, m_defaultArrayColumns);
      Type* dataPtr = array0.data();
      array0 = static_cast<Type>(m_defaultArrayValue);

      // Check that the size and location of the array wasn't changed.
      BRICK_TEST_ASSERT(array0.rows() == m_defaultArrayRows);
      BRICK_TEST_ASSERT(array0.columns() == m_defaultArrayColumns);
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize);
      BRICK_TEST_ASSERT(array0.data() == dataPtr);

      // Check that all values were set by getting a pointer to the
      // first one that's not equal to m_defaultArrayValue
      typename Array2D<Type>::iterator firstRenegade =
        std::find_if(array0.begin(), array0.end(),
                     std::bind2nd(std::not_equal_to<Type>(),
                                  static_cast<Type>(m_defaultArrayValue)));
      BRICK_TEST_ASSERT(firstRenegade == array0.end());
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testApplicationOperator__size_t()
    {
      // This member function returns the (index)th element of the array
      // by reference.
      Array2D<Type> array0(m_fibonacciString);
      for(size_t index = 0; index < array0.size(); ++index) {
        BRICK_TEST_ASSERT(array0(index) == m_fibonacciCArray[index]);
      }
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testApplicationOperatorConst__size_t()
    {
      // This member function returns the (index)th element of the array
      // by reference.
      const Array2D<Type> array0(m_fibonacciString);
      for(size_t index = 0; index < array0.size(); ++index) {
        BRICK_TEST_ASSERT(array0(index) == m_fibonacciCArray[index]);
      }
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testApplicationOperator__size_t__size_t()
    {
      Array2D<Type> array0(m_fibonacciString);

      size_t index0 = 0;
      for(size_t row = 0; row < m_defaultArrayRows;
          ++row) {
        for(size_t column = 0; column < m_defaultArrayColumns; ++column) {
          BRICK_TEST_ASSERT(array0(row, column) == m_fibonacciCArray[index0]);
          ++index0;
        }
      }
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testApplicationOperatorConst__size_t__size_t()
    {
      const Array2D<Type> array0(m_fibonacciString);

      size_t index0 = 0;
      for(size_t row = 0; row < m_defaultArrayRows; ++row) {
        for(size_t column = 0; column < m_defaultArrayColumns; ++column) {
          BRICK_TEST_ASSERT(array0(row, column) == m_fibonacciCArray[index0]);
          ++index0;
        }
      }
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testIndexOperator()
    {
      // This member function returns the (index)th element of the array
      // by reference.
      Array2D<Type> array0(m_fibonacciString);
      for(size_t index = 0; index < array0.size(); ++index) {
        BRICK_TEST_ASSERT(array0[index] == m_fibonacciCArray[index]);
      }
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testIndexOperatorConst()
    {
      // This member function returns the (index)th element of the array
      // by reference.
      const Array2D<Type> array0(m_fibonacciString);
      for(size_t index = 0; index < array0.size(); ++index) {
        BRICK_TEST_ASSERT(array0[index] == m_fibonacciCArray[index]);
      }
    }


    template <class Type>
    void
    Array2DTest<Type>::
    testInitializerList()
    {
      Array2D<Type> array0 {
        {Type(0), Type(1), Type(2), Type(3)},
        {Type(4), Type(5), Type(6), Type(7)},
        {Type(8), Type(9), Type(10), Type(11)}};
      BRICK_TEST_ASSERT(array0.rows() == 3);
      BRICK_TEST_ASSERT(array0.columns() == 4);

      std::size_t cc = 0;
      for(std::size_t ii = 0; ii < array0.rows(); ++ii) {
        for(std::size_t jj = 0; jj < array0.columns(); ++jj) {
          BRICK_TEST_ASSERT(array0(ii, jj) == Type(cc));
          ++cc;
        }
      }
    }


    // General test of square root is not performed, since square root
    // only makes sense for certain types.
    template <class Type>
    void
    Array2DTest<Type>::
    testSquareRoot__Array2D()
    {
      // Empty.
    }


    // Specialization tests square root for Array2D<double>.
    template<>
    void
    Array2DTest<double>::
    testSquareRoot__Array2D()
    {
      const Array2D<double> array0(m_fibonacciString);
      const Array2D<double> array1 = squareRoot(array0);

      BRICK_TEST_ASSERT(array1.rows() == m_defaultArrayRows);
      BRICK_TEST_ASSERT(array1.columns() == m_defaultArrayColumns);
      BRICK_TEST_ASSERT(array1.size() == m_defaultArraySize);
      for(size_t index0 = 0; index0 < array1.size(); ++index0) {
        double targetValue = std::sqrt(
          static_cast<double>(m_fibonacciCArray[index0]));
        // Grr. std::sqrt(x) return value appears to vary slightly from
        // call to call.  Possibly this is due to pentium register
        // precision confusion.
        //
        // BRICK_TEST_ASSERT(array1(index0) == targetValue);
        BRICK_TEST_ASSERT(
          approximatelyEqual(array1(index0), targetValue, 1.0E-14));
      }
    }



    // General test of square root is not performed, since square root
    // only makes sense for certain types.
    template <class Type>
    void
    Array2DTest<Type>::
    testSqrt__Array2D()
    {
      // Empty.
    }


    // Specialization tests square root for Array2D<double>.
    template<>
    void
    Array2DTest<double>::
    testSqrt__Array2D()
    {
      const Array2D<double> array0(m_fibonacciString);
      const Array2D<double> array1 = sqrt(array0);

      BRICK_TEST_ASSERT(array1.rows() == m_defaultArrayRows);
      BRICK_TEST_ASSERT(array1.columns() == m_defaultArrayColumns);
      BRICK_TEST_ASSERT(array1.size() == m_defaultArraySize);
      for(size_t index0 = 0; index0 < array1.size(); ++index0) {
        double targetValue = std::sqrt(
          static_cast<double>(m_fibonacciCArray[index0]));
        // Grr. std::sqrt(x) return value appears to vary slightly from
        // call to call.  Possibly this is due to pentium register
        // precision confusion.
        //
        // BRICK_TEST_ASSERT(array1(index0) == targetValue);
        BRICK_TEST_ASSERT(
          approximatelyEqual(array1(index0), targetValue, 1.0E-14));
      }
    }

  } // namespace numeric

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::Array2DTest<double> currentTest0("double");
  brick::Array2DTest<float> currentTest1("float");
  brick::Array2DTest<int> currentTest2("int");
  brick::Array2DTest<size_t> currentTest3("size_t");
  bool result = (currentTest0.run()
                 && currentTest1.run()
                 && currentTest2.run()
                 && currentTest3.run());
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::Array2DTest<double> currentTest0("double");
  brick::numeric::Array2DTest<float> currentTest1("float");
  brick::numeric::Array2DTest<int> currentTest2("int");
  // brick::numeric::Array2DTest<size_t> currentTest3("size_t");

}

#endif
