/**
***************************************************************************
* @file array1DTest.cpp
* 
* Source file defining Array1DTest class.
*
* Copyright (C) 2004 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <math.h>
#include <iomanip>
#include <sstream>
// #include <brick/common/functional.hh>
// #include <brick/numeric/utilities.hh>
#include <brick/numeric/array1D.hh>
#include <brick/test/functors.hh>

#include <brick/numeric/test/arrayTestCommon.hh>

namespace brick {

  namespace numeric {

    template <class Type>
    class Array1DTest
      : public ArrayTestCommon< Array1DTest<Type>, Array1D<Type>,
                                Array1D<bool> > {

    public:

      typedef Array1DTest<Type> TestFixtureType;

    
      Array1DTest(const std::string& typeName);
      ~Array1DTest();

      // Inherited methods from TestFixture
      void setUp(const std::string&) {}
      void tearDown(const std::string&) {}

      // Inherited methods from ArrayTestCommon
      virtual void
      checkShapeEquality(const Array1D<Type>& array0,
                         const Array1D<Type>& array1);

    
      virtual void
      checkShapeEquality(const Array1D<Type>& array0,
                         const Array1D<bool>& array1);

    
      virtual void
      checkValueEquality(const Array1D<Type>& array0,
                         const Array1D<Type>& array1,
                         typename Array1D<Type>::value_type tolerance);

    
      virtual typename Array1D<Type>::value_type
      getComparisonOperatorThreshold() {
        return m_fibonacciCArray[m_defaultArraySize / 2];
      }

    
      virtual typename Array1D<Type>::value_type
      getEqualityOperatorTarget() {
        return m_squaresCArray[m_defaultArraySize / 2];
      }

    
      virtual Array1D<Type>
      getFibonacciArray() {return Array1D<Type>(m_fibonacciString);}

    
      virtual typename Array1D<Type>::value_type
      getIncrementOperatorArgument() {return static_cast<Type>(4);}

    
      virtual typename Array1D<Type>::value_type
      getMultiplicationOperatorArgument() {return static_cast<Type>(2);}

    
      virtual Array1D<Type>
      getSquaresArray() {return Array1D<Type>(m_squaresString);}

      // Tests of member functions.
      void testConstructor__void();
      void testConstructor__size_t();
      void testConstructor__string();
      void testConstructor__Array1D();
      void testConstructor__size_t__TypePtr();
      void testConstructor__size_t__TypePtr__ReferenceCount();
      void testDestructor();
      void testBegin();
      void testBeginConst();
      void testClear();
      void testCopy();
      void testCopy__Array1D();
      void testCopy__Type2Ptr();
      void testData();
      void testDataConst();
      void testData__size_t();
      void testDataConst__size_t();
      void testEnd();
      void testEndConst();
      void testLength();
      void testReadFromStream();
      void testReinit();
      void testSize();
      void testAssignmentOperator__Array1D();
      void testAssignmentOperator__Type();
      void testApplicationOperator();
      void testApplicationOperatorConst();
      void testIndexOperator();
      void testIndexOperatorConst();


    private:
      size_t m_defaultArraySize;
      size_t m_defaultArrayValue;
      size_t m_defaultIncrementSize;
      size_t m_defaultMultiplier;
      Type* m_fibonacciCArray;
      std::string m_fibonacciString;
      std::string m_illegalString;
      Type* m_squaresCArray;
      std::string m_squaresString;
      double m_testEpsilon;
    
    }; // class Array1DTest


    /* ============== Member Function Definititions ============== */

    template <class Type>
    Array1DTest<Type>::
    Array1DTest(const std::string& typeName)
      : ArrayTestCommon< Array1DTest<Type>, Array1D<Type>, Array1D<bool> >(
        std::string("Array1DTest<") + typeName + ">"),
        m_defaultArraySize(11),
        m_defaultArrayValue(9),
        m_defaultIncrementSize(4),
        m_defaultMultiplier(2),
        m_fibonacciCArray(0),
        m_fibonacciString(""),
        m_illegalString(),      
        m_squaresCArray(0),
        m_squaresString(""),
        m_testEpsilon(1.0e-8)
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testConstructor__void);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__size_t);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__string);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__Array1D);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__size_t__TypePtr);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__size_t__TypePtr__ReferenceCount);
      BRICK_TEST_REGISTER_MEMBER(testDestructor);
      BRICK_TEST_REGISTER_MEMBER(testBegin);
      BRICK_TEST_REGISTER_MEMBER(testBeginConst);
      BRICK_TEST_REGISTER_MEMBER(testClear);
      BRICK_TEST_REGISTER_MEMBER(testCopy);
      BRICK_TEST_REGISTER_MEMBER(testCopy__Array1D);
      BRICK_TEST_REGISTER_MEMBER(testCopy__Type2Ptr);
      BRICK_TEST_REGISTER_MEMBER(testData);
      BRICK_TEST_REGISTER_MEMBER(testDataConst);
      BRICK_TEST_REGISTER_MEMBER(testData__size_t);
      BRICK_TEST_REGISTER_MEMBER(testDataConst__size_t);
      BRICK_TEST_REGISTER_MEMBER(testEnd);
      BRICK_TEST_REGISTER_MEMBER(testEndConst);
      BRICK_TEST_REGISTER_MEMBER(testLength);
      BRICK_TEST_REGISTER_MEMBER(testReadFromStream);
      BRICK_TEST_REGISTER_MEMBER(testReinit);
      BRICK_TEST_REGISTER_MEMBER(testSize);
      BRICK_TEST_REGISTER_MEMBER(testAssignmentOperator__Array1D);
      BRICK_TEST_REGISTER_MEMBER(testAssignmentOperator__Type);
      BRICK_TEST_REGISTER_MEMBER(testApplicationOperator);
      BRICK_TEST_REGISTER_MEMBER(testApplicationOperatorConst);
      BRICK_TEST_REGISTER_MEMBER(testIndexOperator);
      BRICK_TEST_REGISTER_MEMBER(testIndexOperatorConst);
      BRICK_TEST_REGISTER_MEMBER(testOutputOperator);
      BRICK_TEST_REGISTER_MEMBER(testInputOperator);

    
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
        m_squaresCArray[index] = static_cast<Type>((index + 1) * (index + 1));
      }

      // Set up input strings for tests.
      std::ostringstream fibonacciBuffer;
      std::ostringstream squaresBuffer;
      fibonacciBuffer << std::fixed << "["
                      << static_cast<Type>(m_fibonacciCArray[0]);
      squaresBuffer << std::fixed
                    << "[" << static_cast<Type>(m_squaresCArray[0]);
      for(size_t index = 1; index < m_defaultArraySize; ++index) {
        fibonacciBuffer << ", " << static_cast<Type>(m_fibonacciCArray[index]);
        squaresBuffer << ", " << static_cast<Type>(m_squaresCArray[index]);
      }
      fibonacciBuffer << "]";
      squaresBuffer << "]";
      m_fibonacciString = fibonacciBuffer.str();
      m_squaresString = squaresBuffer.str();

      // Set up malformed input string for tests.
      m_illegalString = "[1, 2, 3 4, 5, 6]";
    }


    template <class Type>
    Array1DTest<Type>::
    ~Array1DTest()
    {
      delete[] m_fibonacciCArray;
      delete[] m_squaresCArray;
    }


    template <class Type>
    void
    Array1DTest<Type>::
    checkShapeEquality(const Array1D<Type>& array0,
                       const Array1D<Type>& array1)
    {
      BRICK_TEST_ASSERT(array0.size() == array1.size());
      BRICK_TEST_ASSERT(array0.data() != array1.data());
    }


    template <class Type>
    void
    Array1DTest<Type>::
    checkShapeEquality(const Array1D<Type>& array0,
                       const Array1D<bool>& array1)
    {
      BRICK_TEST_ASSERT(array0.size() == array1.size());
    }


    template <class Type>
    void
    Array1DTest<Type>::
    checkValueEquality(const Array1D<Type>& array0,
                       const Array1D<Type>& array1,
                       typename Array1D<Type>::value_type tolerance)
    {
      this->checkShapeEquality(array0, array1);
      BRICK_TEST_ASSERT(
        std::equal(array0.begin(), array0.end(), array1.begin(),
                   test::ApproximatelyEqualFunctor<Type>(tolerance)));
    }
  
  
    template <class Type>
    void
    Array1DTest<Type>::
    testConstructor__void()
    {
      // Default constructor should initialize to zero size.
      Array1D<Type> array0;
      BRICK_TEST_ASSERT(array0.size() == 0);
    }
  

    template <class Type>
    void
    Array1DTest<Type>::
    testConstructor__size_t()
    {
      // This constructor should initialize to the specified size.
      Array1D<Type> array0(m_defaultArraySize);
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize);
    }
  

    template <class Type>
    void
    Array1DTest<Type>::
    testConstructor__string()
    {
      // This constructor should initialize according to the specified
      // string, which describes an array of m_defaultArraySize elements
      // having values matching those in m_fibonacciCArray.
      Array1D<Type> array0(m_fibonacciString);
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize);
      BRICK_TEST_ASSERT(std::equal(array0.begin(), array0.end(),
                                   m_fibonacciCArray));

      // Construct using a malformed string.
      BRICK_TEST_ASSERT_EXCEPTION(common::ValueException,
                                  Array1D<Type> array1(m_illegalString));

    }
  

    template <class Type>
    void
    Array1DTest<Type>::
    testConstructor__Array1D()
    {
      // This constructor should make an array which references the same
      // data as its argument.
      Array1D<Type>* array0Ptr = new Array1D<Type>(m_fibonacciString);
      Array1D<Type> array1(*array0Ptr);
      BRICK_TEST_ASSERT(array1.data() == array0Ptr->data());
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
    Array1DTest<Type>::
    testConstructor__size_t__TypePtr()
    {
      // This constructor initializes the array to point to the
      // specified data.
      // Start by making an array we can play with.
      Type* cArray = new Type[m_defaultArraySize];
      std::copy(m_fibonacciCArray, m_fibonacciCArray + m_defaultArraySize,
                cArray);
    
      // Create an array using the constructor under test.
      Array1D<Type>* array0Ptr = new Array1D<Type>(m_defaultArraySize, cArray);
      BRICK_TEST_ASSERT(array0Ptr->data() == cArray);
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
    Array1DTest<Type>::
    testConstructor__size_t__TypePtr__ReferenceCount()
    {
      // This constructor should reference the provided test data _and_
      // maintain the reference count.

      // First set up some data to reference.
      Array1D<Type> array0(m_fibonacciString);
      BRICK_TEST_ASSERT(array0.getReferenceCount().getCount() == 1);

      // Create the array.  It should now reference the same data as the
      // first array, and the reference count should be incremented.
      Array1D<Type>* array1Ptr = new Array1D<Type>(
        array0.size(), array0.data(), array0.getReferenceCount());
      BRICK_TEST_ASSERT(array1Ptr->data() == array0.data());
      BRICK_TEST_ASSERT(array1Ptr->size() == array0.size());
      BRICK_TEST_ASSERT(array0.getReferenceCount().getCount() == 2);
      BRICK_TEST_ASSERT(array1Ptr->getReferenceCount().getCount() == 2);
      BRICK_TEST_ASSERT(std::equal(array1Ptr->begin(), array1Ptr->end(),
                                   m_fibonacciCArray));

      // Following deletion of the new array, the reference count should
      // be back down to 1.
      delete array1Ptr;
      BRICK_TEST_ASSERT(array0.getReferenceCount().getCount() == 1);
    }


    template <class Type>
    void
    Array1DTest<Type>::
    testDestructor()
    {
      // No independent test for destructor.
    }


    template <class Type>
    void
    Array1DTest<Type>::
    testBegin()
    {
      // Member function begin() should return an iterator pointing to
      // the first element of the array.
      Array1D<Type> array0(m_defaultArraySize, m_fibonacciCArray);
      BRICK_TEST_ASSERT(&(*(array0.begin())) == m_fibonacciCArray);
    }


    template <class Type>
    void
    Array1DTest<Type>::
    testBeginConst()
    {
      // Member function begin() should return an iterator pointing to
      // the first element of the array.
      const Array1D<Type> array0(m_defaultArraySize, m_fibonacciCArray);
      BRICK_TEST_ASSERT(&(*(array0.begin())) == m_fibonacciCArray);
    }


    template <class Type>
    void
    Array1DTest<Type>::
    testClear()
    {
      // Member function clear should reset the array to zero size,
      // abandoning all contents.
      // First set up some data to reference.
      Array1D<Type> array0(m_fibonacciString);

      // Create a reference count to communicate to the Array1D we're
      // about to create that we're referencing this data.
      common::ReferenceCount refCount(1);

      // Create the test array.  It should now reference the same data
      // as the first array, and the reference count should be
      // incremented.
      Array1D<Type>array1(array0.size(), array0.data(), refCount);
      BRICK_TEST_ASSERT(refCount.getCount() == 2);

      // Clearing the array should release the reference.
      array1.clear();
      BRICK_TEST_ASSERT(array1.size() == 0);
      BRICK_TEST_ASSERT(refCount.getCount() == 1);
    }


    template <class Type>
    void
    Array1DTest<Type>::
    testCopy()
    {
      // Member function copy() should allocate a new array and deep
      // copy the contents of *this.
      Array1D<Type> array0(m_fibonacciString);
      Array1D<Type> array1 = array0.copy();
      BRICK_TEST_ASSERT(array1.size() == array0.size());
      BRICK_TEST_ASSERT(std::equal(array1.begin(), array1.end(), array0.begin()));
      BRICK_TEST_ASSERT(array1.data() != array0.data());    
    }


    template <class Type>
    void
    Array1DTest<Type>::
    testCopy__Array1D()
    {
      // This member function deep copies the contents of source.  It is
      // an error if source does not have the same size as *this.
      Array1D<Type> array0(m_fibonacciString);
      Array1D<Type> array1(array0.size());
      Array1D<Type> array2(array0.size() + 1);

      // Test that the copy goes OK.
      array1.copy(array0);
      BRICK_TEST_ASSERT(array1.size() == array0.size());
      BRICK_TEST_ASSERT(std::equal(array1.begin(), array1.end(), array0.begin()));
      BRICK_TEST_ASSERT(array1.data() != array0.data());

      // Test that size is checked.
      BRICK_TEST_ASSERT_EXCEPTION(common::ValueException, array2.copy(array0));
    }


    template <class Type>
    void
    Array1DTest<Type>::
    testCopy__Type2Ptr()
    {
      // This method simply copies elements from its argument.
      // Begin by making an appropriately sized array with known contents.
      Array1D<Type> array0(m_defaultArraySize);
      array0 = 0;

      array0.copy(m_fibonacciCArray);
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize);
      BRICK_TEST_ASSERT(std::equal(array0.begin(), array0.end(),
                                   m_fibonacciCArray));
      BRICK_TEST_ASSERT(array0.data() != m_fibonacciCArray);
    }


    template <class Type>
    void
    Array1DTest<Type>::
    testData()
    {
      // This method returns a pointer to the internal data store.
      Array1D<Type> array0(m_defaultArraySize, m_fibonacciCArray);
      BRICK_TEST_ASSERT(array0.data() == m_fibonacciCArray);
    }


    template <class Type>
    void
    Array1DTest<Type>::
    testDataConst()
    {
      // This method returns a pointer to the internal data store.
      const Array1D<Type> array0(m_defaultArraySize, m_fibonacciCArray);
      BRICK_TEST_ASSERT(array0.data() == m_fibonacciCArray);
    }


    template <class Type>
    void
    Array1DTest<Type>::
    testData__size_t()
    {
      // This method returns a pointer to the index-th element of the
      // internal data store.
      Array1D<Type> array0(m_defaultArraySize, m_fibonacciCArray);
      BRICK_TEST_ASSERT(array0.data(m_defaultArraySize - 1)
                        == (m_fibonacciCArray + m_defaultArraySize - 1));
    }


    template <class Type>
    void
    Array1DTest<Type>::
    testDataConst__size_t()
    {
      // This method returns a pointer to the index-th element of the
      // internal data store.
      Array1D<Type> array0(m_defaultArraySize, m_fibonacciCArray);
      BRICK_TEST_ASSERT(array0.data(m_defaultArraySize - 1)
                        == (m_fibonacciCArray + m_defaultArraySize - 1));
    }

  
    template <class Type>
    void
    Array1DTest<Type>::
    testEnd()
    {
      // Member function end() should return an iterator pointing to
      // the last element of the array.
      Array1D<Type> array0(m_defaultArraySize, m_fibonacciCArray);
      Type* finalElementPtr = m_fibonacciCArray + m_defaultArraySize;
      BRICK_TEST_ASSERT(&(*(array0.end())) == finalElementPtr);
    }
    

    template <class Type>
    void
    Array1DTest<Type>::
    testEndConst()
    {
      // Member function end() should return an iterator pointing to
      // the last element of the array.
      const Array1D<Type> array0(m_defaultArraySize, m_fibonacciCArray);
      Type* finalElementPtr = m_fibonacciCArray + m_defaultArraySize;
      BRICK_TEST_ASSERT(&(*(array0.end())) == finalElementPtr);
    }


    template <class Type>
    void
    Array1DTest<Type>::
    testLength()
    {
      Array1D<Type> array0;
      Array1D<Type> array1(m_defaultArraySize);
      BRICK_TEST_ASSERT(array0.length() == 0);
      BRICK_TEST_ASSERT(array1.length() == m_defaultArraySize);
    }

  
    template <class Type>
    void
    Array1DTest<Type>::
    testReadFromStream()
    {
      // This member function sets the value of the array from an input
      // stream.
      // Set up input streams
      std::istringstream inputStream0(m_fibonacciString);
      std::istringstream inputStream1(m_illegalString);

      // Try to parse the good stream.
      Array1D<Type> array0;
      array0.readFromStream(inputStream0);
      BRICK_TEST_ASSERT(inputStream0);
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize);
      BRICK_TEST_ASSERT(std::equal(array0.begin(), array0.end(),
                                   m_fibonacciCArray));

      // Try to parse the bad stream.  This should not affect array0 in
      // any way, but it should change the stream state.
      array0.readFromStream(inputStream1);
      BRICK_TEST_ASSERT(!inputStream1);
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize);
      BRICK_TEST_ASSERT(std::equal(array0.begin(), array0.end(),
                                   m_fibonacciCArray));
    }


    template <class Type>
    void
    Array1DTest<Type>::
    testReinit()
    {
      // This method changes the shape of the array and reallocates
      // storage.
      Array1D<Type> array0(m_defaultArraySize, m_fibonacciCArray);
      array0.reinit(m_defaultArraySize + 1);
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize + 1);
      BRICK_TEST_ASSERT(array0.data() != m_fibonacciCArray);
    }

  
    template <class Type>
    void
    Array1DTest<Type>::
    testSize()
    {
      // No explicit test.
    }


    template <class Type>
    void
    Array1DTest<Type>::
    testAssignmentOperator__Array1D()
    {
      // The assignment operator shallow copies the contents of source.
      Array1D<Type> array0(m_fibonacciString);
      Array1D<Type> array1;
      array1 = array0;
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize);
      BRICK_TEST_ASSERT(array1.size() == array0.size());
      BRICK_TEST_ASSERT(array1.data() == array0.data());

      // Make sure the data wasn't changed by the operation.
      BRICK_TEST_ASSERT(std::equal(array1.begin(), array1.end(),
                                   m_fibonacciCArray));
    }

  
    template <class Type>
    void
    Array1DTest<Type>::
    testAssignmentOperator__Type()
    {
      // This member function assigns the same value to every element in
      // the array.
      Array1D<Type> array0(m_defaultArraySize);
      Type* dataPtr = array0.data();
      array0 = static_cast<Type>(m_defaultArrayValue);

      // Check that the size and location of the array wasn't changed.
      BRICK_TEST_ASSERT(array0.size() == m_defaultArraySize);
      BRICK_TEST_ASSERT(array0.data() == dataPtr);

      // Check that all values were set by getting a pointer to the
      //
      // Strangely, std::find_if segfaults for me under g++-4.1.
      // 
      // first one that's not equal to m_defaultArrayValue
      //   typename Array1D<Type>::iterator firstRenegade =
      //     std::find_if(array0.begin(), array0.end(),
      //                  std::bind2nd(std::not_equal_to<Type>(),
      //                               m_defaultArrayValue));
    
      typename Array1D<Type>::iterator firstRenegade = array0.begin();
      while(firstRenegade != array0.end()) {
        if(*firstRenegade != static_cast<Type>(m_defaultArrayValue)) {
          break;
        }
        ++firstRenegade;
      }

      BRICK_TEST_ASSERT(firstRenegade == array0.end());
    }

  
    template <class Type>
    void
    Array1DTest<Type>::
    testApplicationOperator()
    {
      // This member function returns the (index)th element of the array
      // by reference.
      Array1D<Type> array0(m_fibonacciString);
      for(size_t index = 0; index < array0.size(); ++index) {
        BRICK_TEST_ASSERT(array0(index) == m_fibonacciCArray[index]);
      }
    }


    template <class Type>
    void
    Array1DTest<Type>::
    testApplicationOperatorConst()
    {
      // This member function returns the (index)th element of the array
      // by reference.
      const Array1D<Type> array0(m_fibonacciString);
      for(size_t index = 0; index < array0.size(); ++index) {
        BRICK_TEST_ASSERT(array0(index) == m_fibonacciCArray[index]);
      }
    }


    template <class Type>
    void
    Array1DTest<Type>::
    testIndexOperator()
    {
      // This member function returns the (index)th element of the array
      // by reference.
      Array1D<Type> array0(m_fibonacciString);
      for(size_t index = 0; index < array0.size(); ++index) {
        BRICK_TEST_ASSERT(array0[index] == m_fibonacciCArray[index]);
      }
    }


    template <class Type>
    void
    Array1DTest<Type>::
    testIndexOperatorConst()
    {
      // This member function returns the (index)th element of the array
      // by reference.
      const Array1D<Type> array0(m_fibonacciString);
      for(size_t index = 0; index < array0.size(); ++index) {
        BRICK_TEST_ASSERT(array0[index] == m_fibonacciCArray[index]);
      }
    }

  } // namespace numeric
  
} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::Array1DTest<double> currentTest0("double");
  brick::Array1DTest<float> currentTest1("float");
  brick::Array1DTest<int> currentTest2("int");
  brick::Array1DTest<size_t> currentTest3("size_t");
  bool result = (currentTest0.run()
                 && currentTest1.run()
                 && currentTest2.run()
                 && currentTest3.run());
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::Array1DTest<double> currentTest0("double");
  brick::numeric::Array1DTest<float> currentTest1("float");
  brick::numeric::Array1DTest<int> currentTest2("int");
  // brick::numeric::Array1DTest<size_t> currentTest3("size_t"); /* Test assumes signed elements */

}

#endif

