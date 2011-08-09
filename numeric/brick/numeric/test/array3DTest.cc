/**
***************************************************************************
* @file array3DTest.cpp
* 
* Source file defining Array3DTest class.
*
* Copyright (C) 2004,2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <math.h>
#include <iomanip>
#include <sstream>
#include <brick/common/functional.hh>
#include <brick/numeric/utilities.hh>
#include <brick/numeric/array3D.hh>
#include <brick/numeric/test/arrayTestCommon.hh>

namespace brick {

  namespace numeric {
    
    template <class Type>
    class Array3DTest
      : public ArrayTestCommon< Array3DTest<Type>, Array3D<Type>,
                                Array3D<bool> >
    {

    public:

      typedef Array3DTest<Type> TestFixtureType;


      Array3DTest(const std::string& typeName);
      virtual ~Array3DTest();

      // Inherited methods from TestFixture
      void setUp(const std::string&) {}
      void tearDown(const std::string&) {}

      // Inherited methods from ArrayTestCommon
      virtual void
      checkShapeEquality(const Array3D<Type>& array0,
                         const Array3D<Type>& array1);

    
      virtual void
      checkShapeEquality(const Array3D<Type>& array0,
                         const Array3D<bool>& array1);

    
      virtual void
      checkValueEquality(const Array3D<Type>& array0,
                         const Array3D<Type>& array1,
                         typename Array3D<Type>::value_type tolerance);

    
      virtual typename Array3D<Type>::value_type
      getComparisonOperatorThreshold() {
        return m_fibonacciCArray[m_defaultArraySize / 2];
      }

    
      virtual typename Array3D<Type>::value_type
      getEqualityOperatorTarget() {
        return m_squaresCArray[m_defaultArraySize / 2];
      }

    
      virtual Array3D<Type>
      getFibonacciArray() {return Array3D<Type>(m_fibonacciString);}

    
      virtual typename Array3D<Type>::value_type
      getIncrementOperatorArgument() {return static_cast<Type>(4);}

    
      virtual typename Array3D<Type>::value_type
      getMultiplicationOperatorArgument() {return static_cast<Type>(2);}

    
      virtual Array3D<Type>
      getSquaresArray() {return Array3D<Type>(m_fibonacciString);}

    
      // Tests of member functions.
      void testConstructor__void();
      void testConstructor__size_t__size_t__size_t();
      void testConstructor__string();
      void testConstructor__Array3D();
      void testConstructor__size_t__size_t__size_t__TypePtr();
      void testDestructor();
      void testBegin();
      void testBeginConst();
      void testClear();
      void testCopy();
      void testCopy__Array3D();
      void testCopy__Type2Ptr();
      void testData();
      void testDataConst();
      void testData__size_t();
      void testDataConst__size_t();
      void testData__size_t__size_t__size_t();
      void testDataConst__size_t__size_t__size_t();
      void testEnd();
      void testEndConst();
      void testReadFromStream();
      void testReinit();
      void testReshape();
      void testShape();
      void testShape__size_t();
      void testShape0();
      void testShape1();
      void testShape2();
      void testSize();
      void testSlice();
      void testSliceConst();
      void testAssignmentOperator__Array3D();
      void testAssignmentOperator__Type();
      void testApplicationOperator__size_t();
      void testApplicationOperatorConst__size_t();
      void testApplicationOperator__size_t__size_t__size_t();
      void testApplicationOperatorConst__size_t__size_t__size_t();
      void testIndexOperator();
      void testIndexOperatorConst();

      // Tests of non-member functions.

    private:
    
      size_t m_defaultArrayShape0;
      size_t m_defaultArrayShape1;
      size_t m_defaultArrayShape2;
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
    
    }; // class Array3DTest


    /* ============== Member Function Definititions ============== */

    template <class Type>
    Array3DTest<Type>::
    Array3DTest(const std::string& typeName)
      : ArrayTestCommon< Array3DTest<Type>, Array3D<Type>, Array3D<bool> >(
        std::string("Array3DTest<") + typeName + ">"),
        m_defaultArrayShape0(5),
        m_defaultArrayShape1(4),
        m_defaultArrayShape2(3),
        m_defaultArraySize(60),
        m_defaultArrayValue(9),
        m_defaultIncrement(4),
        m_defaultMultiplier(2),
        m_fibonacciCArray(0),
        m_fibonacciString(""),
        m_illegalString(""),      
        m_squaresCArray(0),
        m_squaresString(""),
        m_testEpsilon(1.0e-8)
    {
      // // Register all tests.
      // Tests of member functions.
      BRICK_TEST_REGISTER_MEMBER(testConstructor__void);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__size_t__size_t__size_t);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__string);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__Array3D);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__size_t__size_t__size_t__TypePtr);
      BRICK_TEST_REGISTER_MEMBER(testDestructor);
      BRICK_TEST_REGISTER_MEMBER(testBegin);
      BRICK_TEST_REGISTER_MEMBER(testBeginConst);
      BRICK_TEST_REGISTER_MEMBER(testClear);
      BRICK_TEST_REGISTER_MEMBER(testCopy);
      BRICK_TEST_REGISTER_MEMBER(testCopy__Array3D);
      BRICK_TEST_REGISTER_MEMBER(testCopy__Type2Ptr);
      BRICK_TEST_REGISTER_MEMBER(testData);
      BRICK_TEST_REGISTER_MEMBER(testDataConst);
      BRICK_TEST_REGISTER_MEMBER(testData__size_t);
      BRICK_TEST_REGISTER_MEMBER(testDataConst__size_t);
      BRICK_TEST_REGISTER_MEMBER(testData__size_t__size_t__size_t);
      BRICK_TEST_REGISTER_MEMBER(testDataConst__size_t__size_t__size_t);
      BRICK_TEST_REGISTER_MEMBER(testEnd);
      BRICK_TEST_REGISTER_MEMBER(testEndConst);
      BRICK_TEST_REGISTER_MEMBER(testReadFromStream);
      BRICK_TEST_REGISTER_MEMBER(testReinit);
      BRICK_TEST_REGISTER_MEMBER(testReshape);
      BRICK_TEST_REGISTER_MEMBER(testShape);
      BRICK_TEST_REGISTER_MEMBER(testShape__size_t);
      BRICK_TEST_REGISTER_MEMBER(testShape0);
      BRICK_TEST_REGISTER_MEMBER(testShape1);
      BRICK_TEST_REGISTER_MEMBER(testShape2);
      BRICK_TEST_REGISTER_MEMBER(testSize);
      BRICK_TEST_REGISTER_MEMBER(testSlice);
      BRICK_TEST_REGISTER_MEMBER(testSliceConst);
      BRICK_TEST_REGISTER_MEMBER(testAssignmentOperator__Array3D);
      BRICK_TEST_REGISTER_MEMBER(testAssignmentOperator__Type);
      BRICK_TEST_REGISTER_MEMBER(testApplicationOperator__size_t);
      BRICK_TEST_REGISTER_MEMBER(testApplicationOperatorConst__size_t);
      BRICK_TEST_REGISTER_MEMBER(testApplicationOperator__size_t__size_t__size_t);
      BRICK_TEST_REGISTER_MEMBER(
        testApplicationOperatorConst__size_t__size_t__size_t);
      BRICK_TEST_REGISTER_MEMBER(testIndexOperator);
      BRICK_TEST_REGISTER_MEMBER(testIndexOperatorConst);

    
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
      for(size_t shape0Index = 0; shape0Index < m_defaultArrayShape0;
          ++shape0Index) {
        fibonacciBuffer << "[";
        squaresBuffer << "[";
        for(size_t shape1Index = 0; shape1Index < m_defaultArrayShape1;
            ++shape1Index) {
          fibonacciBuffer << "[" << static_cast<Type>(m_fibonacciCArray[index0]);
          squaresBuffer << "[" << static_cast<Type>(m_squaresCArray[index0]);
          ++index0;
          for(size_t shape2Index = 1; shape2Index < m_defaultArrayShape2;
              ++shape2Index) {
            fibonacciBuffer << ", "
                            << static_cast<Type>(m_fibonacciCArray[index0]);
            squaresBuffer << ", "
                          << static_cast<Type>(m_squaresCArray[index0]);
            ++index0;
          }
          fibonacciBuffer << "]";
          squaresBuffer << "]";
          if(shape1Index != m_defaultArrayShape1 - 1) {
            fibonacciBuffer << ",\n";
            squaresBuffer << ",\n";
          }
        }
        fibonacciBuffer << "]";
        squaresBuffer << "]";
        if(shape0Index != m_defaultArrayShape0 - 1) {
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
    Array3DTest<Type>::
    ~Array3DTest()
    {
      delete[] m_fibonacciCArray;
      delete[] m_squaresCArray;
    }


    template <class Type>
    void
    Array3DTest<Type>::
    checkShapeEquality(const Array3D<Type>& array0,
                       const Array3D<Type>& array1)
    {
      BRICK_TEST_ASSERT(array0.shape0() == array1.shape0());
      BRICK_TEST_ASSERT(array0.shape1() == array1.shape1());
      BRICK_TEST_ASSERT(array0.shape2() == array1.shape2());
      BRICK_TEST_ASSERT(array0.size() == array1.size());
      BRICK_TEST_ASSERT(array0.data() != array1.data());
    }


    template <class Type>
    void
    Array3DTest<Type>::
    checkShapeEquality(const Array3D<Type>& array0,
                       const Array3D<bool>& array1)
    {
      BRICK_TEST_ASSERT(array0.shape0() == array1.shape0());
      BRICK_TEST_ASSERT(array0.shape1() == array1.shape1());
      BRICK_TEST_ASSERT(array0.shape2() == array1.shape2());
      BRICK_TEST_ASSERT(array0.size() == array1.size());
    }


    template <class Type>
    void
    Array3DTest<Type>::
    checkValueEquality(const Array3D<Type>& array0,
                       const Array3D<Type>& array1,
                       typename Array3D<Type>::value_type tolerance)
    {
      this->checkShapeEquality(array0, array1);
      BRICK_TEST_ASSERT(std::equal(array0.begin(), array0.end(), array1.begin(),
                                   ApproximatelyEqualFunctor<Type>(tolerance)));
    }

  
    template <class Type>
    void
    Array3DTest<Type>::
    testConstructor__void()
    {
      // Default constructor should initialize to zero size.
      Array3D<Type> array0;
      BRICK_TEST_ASSERT(array0.shape0() == 0);
      BRICK_TEST_ASSERT(array0.shape1() == 0);
      BRICK_TEST_ASSERT(array0.shape2() == 0);
      BRICK_TEST_ASSERT(array0.size() == 0);
    }
  

    template <class Type>
    void
    Array3DTest<Type>::
    testConstructor__size_t__size_t__size_t()
    {
      // This constructor should initialize to the specified size.
      Array3D<Type> array0(m_defaultArrayShape0,
                           m_defaultArrayShape1,
                           m_defaultArrayShape2);
      BRICK_TEST_ASSERT(array0.shape0() == m_defaultArrayShape0);
      BRICK_TEST_ASSERT(array0.shape1() == m_defaultArrayShape1);
      BRICK_TEST_ASSERT(array0.shape2() == m_defaultArrayShape2);
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize);
    }
  

    template <class Type>
    void
    Array3DTest<Type>::
    testConstructor__string()
    {
      // This constructor should initialize according to the specified
      // string, which describes an array of m_defaultArraySize elements
      // having values matching those in m_fibonacciCArray.
      Array3D<Type> array0(m_fibonacciString);
      BRICK_TEST_ASSERT(array0.shape0() == m_defaultArrayShape0);
      BRICK_TEST_ASSERT(array0.shape1() == m_defaultArrayShape1);
      BRICK_TEST_ASSERT(array0.shape2() == m_defaultArrayShape2);
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize);
      BRICK_TEST_ASSERT(std::equal(array0.begin(), array0.end(),
                                   m_fibonacciCArray));

      // Construct using a malformed string.
      BRICK_TEST_ASSERT_EXCEPTION(common::ValueException,
                                  Array3D<Type> array1(m_illegalString));

    }
  

    template <class Type>
    void
    Array3DTest<Type>::
    testConstructor__Array3D()
    {
      // This constructor should make an array which references the same
      // data as its argument.
      Array3D<Type>* array0Ptr = new Array3D<Type>(m_fibonacciString);
      Array3D<Type> array1(*array0Ptr);
      BRICK_TEST_ASSERT(array1.data() == array0Ptr->data());
      BRICK_TEST_ASSERT(array1.shape0() == array0Ptr->shape0());
      BRICK_TEST_ASSERT(array1.shape1() == array0Ptr->shape1());
      BRICK_TEST_ASSERT(array1.shape2() == array0Ptr->shape2());
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
    Array3DTest<Type>::
    testConstructor__size_t__size_t__size_t__TypePtr()
    {
      // This constructor initializes the array to point to the
      // specified data.
      // Start by making an array we can play with.
      Type* cArray = new Type[m_defaultArraySize];
      std::copy(m_fibonacciCArray, m_fibonacciCArray + m_defaultArraySize,
                cArray);
    
      // Create an array using the constructor under test.
      Array3D<Type>* array0Ptr = new Array3D<Type>(
        m_defaultArrayShape0, m_defaultArrayShape1, m_defaultArrayShape2,
        cArray);
      BRICK_TEST_ASSERT(array0Ptr->data() == cArray);
      BRICK_TEST_ASSERT(array0Ptr->shape0() == m_defaultArrayShape0);
      BRICK_TEST_ASSERT(array0Ptr->shape1() == m_defaultArrayShape1);
      BRICK_TEST_ASSERT(array0Ptr->shape2() == m_defaultArrayShape2);
      BRICK_TEST_ASSERT(array0Ptr->size() == m_defaultArraySize);
      delete array0Ptr;

      // In this case, no reference counting should be done, and the
      // data should not be deleted on deletion of array0Ptr.
      BRICK_TEST_ASSERT(std::equal(cArray, cArray + m_defaultArraySize,
                                   m_fibonacciCArray));

      // Clean up.
      delete[] cArray;
    }
  

    template <class Type>
    void
    Array3DTest<Type>::
    testDestructor()
    {
      // No independent test for destructor.
    }


    template <class Type>
    void
    Array3DTest<Type>::
    testBegin()
    {
      // Member function begin() should return an iterator pointing to
      // the first element of the array.
      Array3D<Type> array0(
        m_defaultArrayShape0, m_defaultArrayShape1, m_defaultArrayShape2,
        m_fibonacciCArray);
      BRICK_TEST_ASSERT(&(*(array0.begin())) == m_fibonacciCArray);
    }


    template <class Type>
    void
    Array3DTest<Type>::
    testBeginConst()
    {
      // Member function begin() should return an iterator pointing to
      // the first element of the array.
      const Array3D<Type> array0(
        m_defaultArrayShape0, m_defaultArrayShape1, m_defaultArrayShape2,
        m_fibonacciCArray);
      BRICK_TEST_ASSERT(&(*(array0.begin())) == m_fibonacciCArray);
    }


    template <class Type>
    void
    Array3DTest<Type>::
    testClear()
    {
      // Member function clear should reset the array to zero size,
      // abandoning all contents.
      // First set up some data to reference.
      Array3D<Type> array0(m_fibonacciString);

      // Clearing the array should release the reference.
      array0.clear();
      BRICK_TEST_ASSERT(array0.size() == 0);

//     // Member function clear should reset the array to zero size,
//     // abandoning all contents.
//     // First set up some data to reference.
//     Array3D<Type> array0(m_fibonacciString);

//     // Create a size_t to hold a reference count.
//     size_t refCount = 1;

//     // Create the test array.  It should now reference the same data
//     // as the first array, and the reference count should be
//     // incremented.
//     Array3D<Type>array1(array0.size(), array0.data(), &refCount);
//     BRICK_TEST_ASSERT(refCount == 2);

//     // Clearing the array should release the reference.
//     array1.clear();
//     BRICK_TEST_ASSERT(array1.size() == 0);
//     BRICK_TEST_ASSERT(refCount == 1);
    }


    template <class Type>
    void
    Array3DTest<Type>::
    testCopy()
    {
      // Member function copy() should allocate a new array and deep
      // copy the contents of *this.
      Array3D<Type> array0(m_fibonacciString);
      Array3D<Type> array1 = array0.copy();
      BRICK_TEST_ASSERT(array1.size() == array0.size());
      BRICK_TEST_ASSERT(array1.shape0() == array0.shape0());
      BRICK_TEST_ASSERT(array1.shape1() == array0.shape1());
      BRICK_TEST_ASSERT(array1.shape2() == array0.shape2());
      BRICK_TEST_ASSERT(std::equal(array1.begin(), array1.end(), array0.begin()));
      BRICK_TEST_ASSERT(array1.data() != array0.data());    
    }


    template <class Type>
    void
    Array3DTest<Type>::
    testCopy__Array3D()
    {
      // This member function deep copies the contents of source.  It is
      // an error if source does not have the same size as *this.
      Array3D<Type> array0(m_fibonacciString);
      Array3D<Type> array1(array0.shape0(), array0.shape1(), array0.shape2());
      Array3D<Type> array2(array0.shape0(), array0.shape1(),
                           array0.shape2() + 1);
      Array3D<Type> array3(array0.shape0(), array0.shape1() + 1,
                           array0.shape2());
      Array3D<Type> array4(array0.shape0() + 1, array0.shape1(),
                           array0.shape2());

      // Test that the copy goes OK.
      array1.copy(array0);
      BRICK_TEST_ASSERT(array1.shape0() == array0.shape0());
      BRICK_TEST_ASSERT(array1.shape1() == array0.shape1());
      BRICK_TEST_ASSERT(array1.shape2() == array0.shape2());
      BRICK_TEST_ASSERT(array1.size() == array0.size());
      BRICK_TEST_ASSERT(std::equal(array1.begin(), array1.end(), array0.begin()));
      BRICK_TEST_ASSERT(array1.data() != array0.data());

      // Test that size is checked.
      BRICK_TEST_ASSERT_EXCEPTION(common::ValueException, array2.copy(array0));
      BRICK_TEST_ASSERT_EXCEPTION(common::ValueException, array3.copy(array0));
      BRICK_TEST_ASSERT_EXCEPTION(common::ValueException, array4.copy(array0));
    }


    template <class Type>
    void
    Array3DTest<Type>::
    testCopy__Type2Ptr()
    {
      // This method simply copies elements from its argument.
      // Begin by making an appropriately sized array with known contents.
      Array3D<Type> array0(m_defaultArrayShape0, m_defaultArrayShape1,
                           m_defaultArrayShape2);
      array0 = 0;

      array0.copy(m_fibonacciCArray);
      BRICK_TEST_ASSERT(array0.shape0() == m_defaultArrayShape0);
      BRICK_TEST_ASSERT(array0.shape1() == m_defaultArrayShape1);
      BRICK_TEST_ASSERT(array0.shape2() == m_defaultArrayShape2);
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize);
      BRICK_TEST_ASSERT(std::equal(array0.begin(), array0.end(),
                                   m_fibonacciCArray));
      BRICK_TEST_ASSERT(array0.data() != m_fibonacciCArray);
    }


    template <class Type>
    void
    Array3DTest<Type>::
    testData()
    {
      // This method returns a pointer to the internal data store.
      Array3D<Type> array0(
        m_defaultArrayShape0, m_defaultArrayShape1, m_defaultArrayShape2,
        m_fibonacciCArray);
      BRICK_TEST_ASSERT(array0.data() == m_fibonacciCArray);
    }


    template <class Type>
    void
    Array3DTest<Type>::
    testDataConst()
    {
      // This method returns a pointer to the internal data store.
      const Array3D<Type> array0(
        m_defaultArrayShape0, m_defaultArrayShape1, m_defaultArrayShape2,
        m_fibonacciCArray);
      BRICK_TEST_ASSERT(array0.data() == m_fibonacciCArray);
    }


    template <class Type>
    void
    Array3DTest<Type>::
    testData__size_t()
    {
      // This method returns a pointer to the index-th element of the
      // internal data store.
      Array3D<Type> array0(
        m_defaultArrayShape0, m_defaultArrayShape1, m_defaultArrayShape2,
        m_fibonacciCArray);
      BRICK_TEST_ASSERT(array0.data(m_defaultArraySize - 1)
                        == (m_fibonacciCArray + m_defaultArraySize - 1));
    }


    template <class Type>
    void
    Array3DTest<Type>::
    testDataConst__size_t()
    {
      // This method returns a pointer to the index-th element of the
      // internal data store.
      Array3D<Type> array0(
        m_defaultArrayShape0, m_defaultArrayShape1, m_defaultArrayShape2,
        m_fibonacciCArray);
      BRICK_TEST_ASSERT(array0.data(m_defaultArraySize - 1)
                        == (m_fibonacciCArray + m_defaultArraySize - 1));
    }


    template <class Type>
    void
    Array3DTest<Type>::
    testData__size_t__size_t__size_t()
    {
      // This method returns a pointer to the index-th element of the
      // internal data store.
      Array3D<Type> array0(
        m_defaultArrayShape0, m_defaultArrayShape1, m_defaultArrayShape2,
        m_fibonacciCArray);
      BRICK_TEST_ASSERT(array0.data(0, 0, 0) == m_fibonacciCArray);
      BRICK_TEST_ASSERT(array0.data(0, 0, 1) == (m_fibonacciCArray + 1));
      BRICK_TEST_ASSERT(array0.data(0, 1, 0) 
                        == (m_fibonacciCArray + m_defaultArrayShape2));
      BRICK_TEST_ASSERT(array0.data(0, 1, 1) 
                        == (m_fibonacciCArray + m_defaultArrayShape2 + 1));
      BRICK_TEST_ASSERT(array0.data(1, 0, 0) 
                        == (m_fibonacciCArray
                            + m_defaultArrayShape1 * m_defaultArrayShape2));
      BRICK_TEST_ASSERT(array0.data(1, 0, 1) 
                        == (m_fibonacciCArray
                            + m_defaultArrayShape1 * m_defaultArrayShape2 + 1));
      BRICK_TEST_ASSERT(array0.data(1, 1, 0) 
                        == (m_fibonacciCArray
                            + m_defaultArrayShape1 * m_defaultArrayShape2
                            + m_defaultArrayShape2));
      BRICK_TEST_ASSERT(array0.data(1, 1, 1) 
                        == (m_fibonacciCArray
                            + m_defaultArrayShape1 * m_defaultArrayShape2
                            + m_defaultArrayShape2 + 1));
      BRICK_TEST_ASSERT(array0.data(m_defaultArrayShape0 - 1,
                                    m_defaultArrayShape1 - 1,
                                    m_defaultArrayShape2 - 1)
                        == (m_fibonacciCArray + m_defaultArraySize - 1));
    }
  

    template <class Type>
    void
    Array3DTest<Type>::
    testDataConst__size_t__size_t__size_t()
    {
      // This method returns a pointer to the index-th element of the
      // internal data store.
      const Array3D<Type> array0(
        m_defaultArrayShape0, m_defaultArrayShape1, m_defaultArrayShape2,
        m_fibonacciCArray);
      BRICK_TEST_ASSERT(array0.data(0, 0, 0) == m_fibonacciCArray);
      BRICK_TEST_ASSERT(array0.data(0, 0, 1) == (m_fibonacciCArray + 1));
      BRICK_TEST_ASSERT(array0.data(0, 1, 0) 
                        == (m_fibonacciCArray + m_defaultArrayShape2));
      BRICK_TEST_ASSERT(array0.data(0, 1, 1) 
                        == (m_fibonacciCArray + m_defaultArrayShape2 + 1));
      BRICK_TEST_ASSERT(array0.data(1, 0, 0) 
                        == (m_fibonacciCArray
                            + m_defaultArrayShape1 * m_defaultArrayShape2));
      BRICK_TEST_ASSERT(array0.data(1, 0, 1) 
                        == (m_fibonacciCArray
                            + m_defaultArrayShape1 * m_defaultArrayShape2 + 1));
      BRICK_TEST_ASSERT(array0.data(1, 1, 0) 
                        == (m_fibonacciCArray
                            + m_defaultArrayShape1 * m_defaultArrayShape2
                            + m_defaultArrayShape2));
      BRICK_TEST_ASSERT(array0.data(1, 1, 1) 
                        == (m_fibonacciCArray
                            + m_defaultArrayShape1 * m_defaultArrayShape2
                            + m_defaultArrayShape2 + 1));
      BRICK_TEST_ASSERT(array0.data(m_defaultArrayShape0 - 1,
                                    m_defaultArrayShape1 - 1,
                                    m_defaultArrayShape2 - 1)
                        == (m_fibonacciCArray + m_defaultArraySize - 1));
    }


    template <class Type>
    void
    Array3DTest<Type>::
    testEnd()
    {
      // Member function end() should return an iterator pointing to
      // the last element of the array.
      Array3D<Type> array0(
        m_defaultArrayShape0, m_defaultArrayShape1, m_defaultArrayShape2,
        m_fibonacciCArray);
      Type* finalElementPtr = m_fibonacciCArray + m_defaultArraySize;
      BRICK_TEST_ASSERT(&(*(array0.end())) == finalElementPtr);
    }
    

    template <class Type>
    void
    Array3DTest<Type>::
    testEndConst()
    {
      // Member function end() should return an iterator pointing to
      // the last element of the array.
      const Array3D<Type> array0(
        m_defaultArrayShape0, m_defaultArrayShape1, m_defaultArrayShape2,
        m_fibonacciCArray);
      Type* finalElementPtr = m_fibonacciCArray + m_defaultArraySize;
      BRICK_TEST_ASSERT(&(*(array0.end())) == finalElementPtr);
    }


    template <class Type>
    void
    Array3DTest<Type>::
    testReadFromStream()
    {
      // This member function sets the value of the array from an input
      // stream.
      // Set up input streams
      std::istringstream inputStream0(m_fibonacciString);
      std::istringstream inputStream1(m_illegalString);

      // Try to parse the good stream.
      Array3D<Type> array0;
      array0.readFromStream(inputStream0);
      BRICK_TEST_ASSERT(inputStream0);
      BRICK_TEST_ASSERT(array0.shape0() == m_defaultArrayShape0);
      BRICK_TEST_ASSERT(array0.shape1() == m_defaultArrayShape1);
      BRICK_TEST_ASSERT(array0.shape2() == m_defaultArrayShape2);
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize);
      BRICK_TEST_ASSERT(std::equal(array0.begin(), array0.end(),
                                   m_fibonacciCArray));

      // Try to parse the bad stream.  This should not affect array0 in
      // any way, but it should change the stream state.
      array0.readFromStream(inputStream1);
      BRICK_TEST_ASSERT(!inputStream1);
      BRICK_TEST_ASSERT(array0.shape0() == m_defaultArrayShape0);
      BRICK_TEST_ASSERT(array0.shape1() == m_defaultArrayShape1);
      BRICK_TEST_ASSERT(array0.shape2() == m_defaultArrayShape2);
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize);
      BRICK_TEST_ASSERT(std::equal(array0.begin(), array0.end(),
                                   m_fibonacciCArray));
    }


    template <class Type>
    void
    Array3DTest<Type>::
    testReinit()
    {
      // This method changes the shape of the array and reallocates
      // storage.
      Array3D<Type> array0(
        m_defaultArrayShape0, m_defaultArrayShape1, m_defaultArrayShape2,
        m_fibonacciCArray);
      array0.reinit(m_defaultArrayShape0 + 1, m_defaultArrayShape1 + 2,
                    m_defaultArrayShape2 + 3);
      BRICK_TEST_ASSERT(array0.shape0() == m_defaultArrayShape0 + 1);
      BRICK_TEST_ASSERT(array0.shape1() == m_defaultArrayShape1 + 2);
      BRICK_TEST_ASSERT(array0.shape2() == m_defaultArrayShape2 + 3);
      BRICK_TEST_ASSERT(array0.size() ==
                        (m_defaultArrayShape0 + 1)
                        * (m_defaultArrayShape1 + 2)
                        * (m_defaultArrayShape2 + 3));
      BRICK_TEST_ASSERT(array0.data() != m_fibonacciCArray);
    }

  
    template <class Type>
    void
    Array3DTest<Type>::
    testReshape()
    {
      // This method changes the shape of the array withot reallocating
      // storage.
      Array3D<Type> array0(
        m_defaultArrayShape0, m_defaultArrayShape1, m_defaultArrayShape2,
        m_fibonacciCArray);

      // Test that size is checked.
      BRICK_TEST_ASSERT_EXCEPTION(
        common::ValueException, 
        array0.reshape(static_cast<int>(m_defaultArrayShape0) + 1,
                       static_cast<int>(m_defaultArrayShape1) + 2,
                       static_cast<int>(m_defaultArrayShape2) + 3));

      // Check that correctly sized reshapes work.
      array0.reshape(1, 1, static_cast<int>(m_defaultArraySize));
      BRICK_TEST_ASSERT(array0.shape0() == 1);
      BRICK_TEST_ASSERT(array0.shape1() == 1);
      BRICK_TEST_ASSERT(array0.shape2() == m_defaultArraySize);
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize);
      BRICK_TEST_ASSERT(array0.data() == m_fibonacciCArray);
    }


    template <class Type>
    void
    Array3DTest<Type>::
    testShape()
    {
      Array3D<Type> array0(
        m_defaultArrayShape0, m_defaultArrayShape1, m_defaultArrayShape2,
        m_fibonacciCArray);
      Array1D<size_t> shapeArray = array0.shape();
      BRICK_TEST_ASSERT(shapeArray.size() == 3);
      BRICK_TEST_ASSERT(shapeArray(0) == m_defaultArrayShape0);
      BRICK_TEST_ASSERT(shapeArray(1) == m_defaultArrayShape1);
      BRICK_TEST_ASSERT(shapeArray(2) == m_defaultArrayShape2);
    }


    template <class Type>
    void
    Array3DTest<Type>::
    testShape__size_t()
    {
      Array3D<Type> array0(
        m_defaultArrayShape0, m_defaultArrayShape1, m_defaultArrayShape2,
        m_fibonacciCArray);
      BRICK_TEST_ASSERT(array0.shape(0) == m_defaultArrayShape0);
      BRICK_TEST_ASSERT(array0.shape(1) == m_defaultArrayShape1);
      BRICK_TEST_ASSERT(array0.shape(2) == m_defaultArrayShape2);
    }


    template <class Type>
    void
    Array3DTest<Type>::
    testShape0()
    {
      // No explicit test.
    }

  
    template <class Type>
    void
    Array3DTest<Type>::
    testShape1()
    {
      // No explicit test.
    }

  
    template <class Type>
    void
    Array3DTest<Type>::
    testShape2()
    {
      // No explicit test.
    }
  

    template <class Type>
    void
    Array3DTest<Type>::
    testSize()
    {
      // No explicit test.
    }


    template <class Type>
    void
    Array3DTest<Type>::
    testSlice()
    {
      Array3D<Type> array0(m_fibonacciString);
      Type* dataPtr = array0.data();
      Array2D<Type> array1 = array0.slice(1);

      BRICK_TEST_ASSERT(array0.data() == dataPtr);
      BRICK_TEST_ASSERT(
        array1.size() == m_defaultArrayShape1 * m_defaultArrayShape2);
      BRICK_TEST_ASSERT(
        array1.data() == dataPtr + m_defaultArrayShape1 * m_defaultArrayShape2);
      for(size_t index0 = 0;
          index0 < m_defaultArrayShape1 * m_defaultArrayShape2;
          ++index0) {
        BRICK_TEST_ASSERT(
          array1[index0] ==
          m_fibonacciCArray[m_defaultArrayShape1 * m_defaultArrayShape2
                            + index0]);
      }
    }

  
    template <class Type>
    void
    Array3DTest<Type>::
    testSliceConst()
    {
      const Array3D<Type> array0(m_fibonacciString);
      const Type* dataPtr = array0.data();
      const Array2D<Type> array1 = array0.slice(1);

      BRICK_TEST_ASSERT(array0.data() == dataPtr);
      BRICK_TEST_ASSERT(
        array1.size() == m_defaultArrayShape1 * m_defaultArrayShape2);
      BRICK_TEST_ASSERT(
        array1.data() == dataPtr + m_defaultArrayShape1 * m_defaultArrayShape2);
      for(size_t index0 = 0;
          index0 < m_defaultArrayShape1 * m_defaultArrayShape2;
          ++index0) {
        BRICK_TEST_ASSERT(
          array1[index0] ==
          m_fibonacciCArray[m_defaultArrayShape1 * m_defaultArrayShape2
                            + index0]);
      }
    }

  
    template <class Type>
    void
    Array3DTest<Type>::
    testAssignmentOperator__Array3D()
    {
      // The assignment operator shallow copies the contents of source.
      Array3D<Type> array0(m_fibonacciString);
      Array3D<Type> array1;
      array1 = array0;
      BRICK_TEST_ASSERT(array0.shape0() == m_defaultArrayShape0);
      BRICK_TEST_ASSERT(array0.shape1() == m_defaultArrayShape1);
      BRICK_TEST_ASSERT(array0.shape2() == m_defaultArrayShape2);
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize);
      BRICK_TEST_ASSERT(array1.shape0() == m_defaultArrayShape0);
      BRICK_TEST_ASSERT(array1.shape1() == m_defaultArrayShape1);
      BRICK_TEST_ASSERT(array1.shape2() == m_defaultArrayShape2);
      BRICK_TEST_ASSERT(array1.size() == array0.size());
      BRICK_TEST_ASSERT(array1.data() == array0.data());

      // Make sure the data wasn't changed by the operation.
      BRICK_TEST_ASSERT(std::equal(array1.begin(), array1.end(),
                                   m_fibonacciCArray));
    }

  
    template <class Type>
    void
    Array3DTest<Type>::
    testAssignmentOperator__Type()
    {
      // This member function assigns the same value to every element in
      // the array.
      Array3D<Type> array0(m_defaultArrayShape0, m_defaultArrayShape1,
                           m_defaultArrayShape2);
      Type* dataPtr = array0.data();
      array0 = static_cast<Type>(m_defaultArrayValue);

      // Check that the size and location of the array wasn't changed.
      BRICK_TEST_ASSERT(array0.shape0() == m_defaultArrayShape0);
      BRICK_TEST_ASSERT(array0.shape1() == m_defaultArrayShape1);
      BRICK_TEST_ASSERT(array0.shape2() == m_defaultArrayShape2);
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize);
      BRICK_TEST_ASSERT(array0.data() == dataPtr);

      // Check that all values were set by getting a pointer to the
      // first one that's not equal to m_defaultArrayValue
      typename Array3D<Type>::iterator firstRenegade =
        std::find_if(array0.begin(), array0.end(),
                     std::bind2nd(std::not_equal_to<Type>(),
                                  static_cast<Type>(m_defaultArrayValue)));
      BRICK_TEST_ASSERT(firstRenegade == array0.end());
    }

  
    template <class Type>
    void
    Array3DTest<Type>::
    testApplicationOperator__size_t()
    {
      // This member function returns the (index)th element of the array
      // by reference.
      Array3D<Type> array0(m_fibonacciString);
      for(size_t index = 0; index < array0.size(); ++index) {
        BRICK_TEST_ASSERT(array0(index) == m_fibonacciCArray[index]);
      }
    }


    template <class Type>
    void
    Array3DTest<Type>::
    testApplicationOperatorConst__size_t()
    {
      // This member function returns the (index)th element of the array
      // by reference.
      const Array3D<Type> array0(m_fibonacciString);
      for(size_t index = 0; index < array0.size(); ++index) {
        BRICK_TEST_ASSERT(array0(index) == m_fibonacciCArray[index]);
      }
    }


    template <class Type>
    void
    Array3DTest<Type>::
    testApplicationOperator__size_t__size_t__size_t()
    {
      Array3D<Type> array0(m_fibonacciString);

      size_t flatIndex = 0;
      for(size_t index0 = 0; index0 < m_defaultArrayShape0; ++index0) {
        for(size_t index1 = 0; index1 < m_defaultArrayShape1; ++index1) {
          for(size_t index2 = 0; index2 < m_defaultArrayShape2; ++index2) {
            BRICK_TEST_ASSERT(array0(index0, index1, index2) ==
                              m_fibonacciCArray[flatIndex]);
            ++flatIndex;
          }
        }
      }
    }


    template <class Type>
    void
    Array3DTest<Type>::
    testApplicationOperatorConst__size_t__size_t__size_t()
    {
      const Array3D<Type> array0(m_fibonacciString);

      size_t flatIndex = 0;
      for(size_t index0 = 0; index0 < m_defaultArrayShape0; ++index0) {
        for(size_t index1 = 0; index1 < m_defaultArrayShape1; ++index1) {
          for(size_t index2 = 0; index2 < m_defaultArrayShape2; ++index2) {
            BRICK_TEST_ASSERT(array0(index0, index1, index2) ==
                              m_fibonacciCArray[flatIndex]);
            ++flatIndex;
          }
        }
      }
    }


    template <class Type>
    void
    Array3DTest<Type>::
    testIndexOperator()
    {
      // This member function returns the (index)th element of the array
      // by reference.
      Array3D<Type> array0(m_fibonacciString);
      for(size_t index = 0; index < array0.size(); ++index) {
        BRICK_TEST_ASSERT(array0[index] == m_fibonacciCArray[index]);
      }
    }


    template <class Type>
    void
    Array3DTest<Type>::
    testIndexOperatorConst()
    {
      // This member function returns the (index)th element of the array
      // by reference.
      const Array3D<Type> array0(m_fibonacciString);
      for(size_t index = 0; index < array0.size(); ++index) {
        BRICK_TEST_ASSERT(array0[index] == m_fibonacciCArray[index]);
      }
    }

  } // namespace numeric
  
} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::numeric::Array3DTest<double> currentTest0("double");
  brick::numeric::Array3DTest<float> currentTest1("float");
  brick::numeric::Array3DTest<int> currentTest2("int");
  brick::numeric::Array3DTest<size_t> currentTest3("size_t");
  bool result = (currentTest0.run()
                 && currentTest1.run()
                 && currentTest2.run()
                 && currentTest3.run());
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::Array3DTest<double> currentTest0("double");
  brick::numeric::Array3DTest<float> currentTest1("float");
  brick::numeric::Array3DTest<int> currentTest2("int");
  // brick::numeric::Array3DTest<size_t> currentTest3("size_t");
  
}

#endif
