/**
***************************************************************************
* @file brick/numeric/test/arrayNDTest.cpp
*
* Source file defining ArrayNDTest class.
*
* Copyright (C) 2004-2008, 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <math.h>
#include <iomanip>
#include <sstream>
#include <brick/common/functional.hh>
#include <brick/numeric/utilities.hh>
#include <brick/numeric/arrayND.hh>
#include <brick/numeric/test/arrayTestCommon.hh>

namespace brick {

  namespace numeric {

    template <class Type>
    class ArrayNDTest
      : public ArrayTestCommon< ArrayNDTest<Type>, ArrayND<3, Type>,
                                ArrayND<3, bool> > {

    public:

      typedef ArrayNDTest<Type> TestFixtureType;


      ArrayNDTest(const std::string& typeName);
      ~ArrayNDTest();

      // Inherited methods from TestFixture
      void setUp(const std::string&) {}
      void tearDown(const std::string&) {}

      // Inherited methods from ArrayTestCommon
      virtual void
      checkShapeEquality(const ArrayND<3, Type>& array0,
                         const ArrayND<3, Type>& array1);


      virtual void
      checkShapeEquality(const ArrayND<3, Type>& array0,
                         const ArrayND<3, bool>& array1);


      virtual void
      checkValueEquality(const ArrayND<3, Type>& array0,
                         const ArrayND<3, Type>& array1,
                         typename ArrayND<3, Type>::value_type tolerance);


      virtual typename ArrayND<3, Type>::value_type
      getComparisonOperatorThreshold() {
        return m_fibonacciCArray[m_defaultArraySize / 2];
      }


      virtual typename ArrayND<3, Type>::value_type
      getEqualityOperatorTarget() {
        return m_squaresCArray[m_defaultArraySize / 2];
      }


      virtual ArrayND<3, Type>
      getFibonacciArray() {return ArrayND<3, Type>(m_fibonacciString);}


      virtual typename ArrayND<3, Type>::value_type
      getIncrementOperatorArgument() {return static_cast<Type>(4);}


      virtual typename ArrayND<3, Type>::value_type
      getMultiplicationOperatorArgument() {return static_cast<Type>(2);}


      virtual ArrayND<3, Type>
      getSquaresArray() {return ArrayND<3, Type>(m_squaresString);}

      // Tests of member functions.
      void testConstructor__void();
      void testConstructor__Array1D();
      void testConstructor__size_t();
      void testConstructor__string();
      void testConstructor__ArrayND();
      void testConstructor__size_t__TypePtr();
      void testConstructor__size_t__TypePtr__size_tPtr();
      void testDestructor();
      void testBegin();
      void testBeginConst();
      void testClear();
      void testCopy();
      void testCopy__ArrayND();
      void testCopy__Type2Ptr();
      void testData();
      void testDataConst();
      void testEnd();
      void testEndConst();
      void testReinit();
      void testSize();
      void testAssignmentOperator__ArrayND();
      void testAssignmentOperator__Type();
      void testApplicationOperator();
      void testApplicationOperatorConst();
      void testIndexOperator();
      void testIndexOperatorConst();


    private:

      bool
      isShapeMatch(Array1D<size_t> const& shape0, Array1D<size_t> const& shape1);


      Array1D<size_t> m_defaultArrayShape;
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

    }; // class ArrayNDTest


    /* ============== Member Function Definititions ============== */

    template <class Type>
    ArrayNDTest<Type>::
    ArrayNDTest(const std::string& typeName)
      : ArrayTestCommon< ArrayNDTest<Type>, ArrayND<3, Type>, ArrayND<3, bool> >(
        std::string("ArrayNDTest<") + typeName + ">"),
        m_defaultArrayShape(3),
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
      BRICK_TEST_REGISTER_MEMBER(testConstructor__Array1D);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__string);
      BRICK_TEST_REGISTER_MEMBER(testConstructor__ArrayND);
      BRICK_TEST_REGISTER_MEMBER(testDestructor);
      BRICK_TEST_REGISTER_MEMBER(testBegin);
      BRICK_TEST_REGISTER_MEMBER(testBeginConst);
      BRICK_TEST_REGISTER_MEMBER(testClear);
      BRICK_TEST_REGISTER_MEMBER(testCopy);
      BRICK_TEST_REGISTER_MEMBER(testCopy__ArrayND);
      BRICK_TEST_REGISTER_MEMBER(testCopy__Type2Ptr);
      BRICK_TEST_REGISTER_MEMBER(testData);
      BRICK_TEST_REGISTER_MEMBER(testDataConst);
      BRICK_TEST_REGISTER_MEMBER(testEnd);
      BRICK_TEST_REGISTER_MEMBER(testEndConst);
      BRICK_TEST_REGISTER_MEMBER(testReinit);
      BRICK_TEST_REGISTER_MEMBER(testSize);
      BRICK_TEST_REGISTER_MEMBER(testAssignmentOperator__ArrayND);
      BRICK_TEST_REGISTER_MEMBER(testAssignmentOperator__Type);
      BRICK_TEST_REGISTER_MEMBER(testApplicationOperator);
      BRICK_TEST_REGISTER_MEMBER(testApplicationOperatorConst);
      BRICK_TEST_REGISTER_MEMBER(testIndexOperator);
      BRICK_TEST_REGISTER_MEMBER(testIndexOperatorConst);
      BRICK_TEST_REGISTER_MEMBER(testOutputOperator);
      BRICK_TEST_REGISTER_MEMBER(testInputOperator);

      // Define default array shape.
      m_defaultArrayShape[0] = 2;
      m_defaultArrayShape[1] = 3;
      m_defaultArrayShape[2] = 5;

      // Figure out how many elements there are in array of this
      // shape.
      size_t dataSize = std::accumulate(m_defaultArrayShape.begin(),
                                        m_defaultArrayShape.end(),
                                        1.0,
                                        std::multiplies<size_t>());

      // Set up fibonacci data for tests.
      m_fibonacciCArray = new Type[dataSize];
      m_fibonacciCArray[0] = static_cast<Type>(1);
      m_fibonacciCArray[1] = static_cast<Type>(2);
      for(size_t index = 2; index < dataSize; ++index) {
        m_fibonacciCArray[index] = (m_fibonacciCArray[index - 1]
                                    + m_fibonacciCArray[index - 2]);
      }

      // Set up squares data for tests.
      m_squaresCArray = new Type[dataSize];
      for(size_t index = 0; index < dataSize; ++index) {
        m_squaresCArray[index] = static_cast<Type>((index + 1) * (index + 1));
      }

      // Set up input strings for tests.
      std::ostringstream fibonacciBuffer;
      std::ostringstream squaresBuffer;
      fibonacciBuffer
        << std::fixed << "ArrayND {\n"
        << "  shape: " << m_defaultArrayShape << "\n"
        << "  data: "
        << Array1D<Type>(dataSize, m_fibonacciCArray) << "\n"
        << "}";
      squaresBuffer
        << std::fixed << "ArrayND {\n"
        << "  shape: " << m_defaultArrayShape << "\n"
        << "  data: "
        << Array1D<Type>(dataSize, m_squaresCArray) << "\n"
        << "}";
      m_fibonacciString = fibonacciBuffer.str();
      m_squaresString = squaresBuffer.str();

      // Set up malformed input string for tests.
      m_illegalString = "[1, 2, 3 4, 5, 6]";
    }


    template <class Type>
    ArrayNDTest<Type>::
    ~ArrayNDTest()
    {
      delete[] m_fibonacciCArray;
      delete[] m_squaresCArray;
    }


    template <class Type>
    void
    ArrayNDTest<Type>::
    checkShapeEquality(const ArrayND<3, Type>& array0,
                       const ArrayND<3, Type>& array1)
    {
      Array1D<size_t> shape0 = array0.getShape();
      Array1D<size_t> shape1 = array1.getShape();
      BRICK_TEST_ASSERT(shape0.size() == shape1.size());
      BRICK_TEST_ASSERT(std::equal(shape0.begin(), shape0.end(), shape1.begin()));
      BRICK_TEST_ASSERT(array0.data() != array1.data());
    }


    template <class Type>
    void
    ArrayNDTest<Type>::
    checkShapeEquality(const ArrayND<3, Type>& array0,
                       const ArrayND<3, bool>& array1)
    {
      Array1D<size_t> shape0 = array0.getShape();
      Array1D<size_t> shape1 = array1.getShape();
      BRICK_TEST_ASSERT(shape0.size() == shape1.size());
      BRICK_TEST_ASSERT(std::equal(shape0.begin(), shape0.end(), shape1.begin()));
    }


    template <class Type>
    void
    ArrayNDTest<Type>::
    checkValueEquality(const ArrayND<3, Type>& array0,
                       const ArrayND<3, Type>& array1,
                       typename ArrayND<3, Type>::value_type tolerance)
    {
      this->checkShapeEquality(array0, array1);
      BRICK_TEST_ASSERT(std::equal(array0.begin(), array0.end(), array1.begin(),
                                 ApproximatelyEqualFunctor<Type>(tolerance)));
    }


    template <class Type>
    void
    ArrayNDTest<Type>::
    testConstructor__void()
    {
      // Default constructor should initialize to zero size.
      ArrayND<3, Type> array0;
      BRICK_TEST_ASSERT(array0.size() == 0);
    }


    template <class Type>
    void
    ArrayNDTest<Type>::
    testConstructor__Array1D()
    {
      // This constructor should initialize to the specified size.
      ArrayND<3, Type> array0(m_defaultArrayShape);
      BRICK_TEST_ASSERT(this->isShapeMatch(array0.getShape(), m_defaultArrayShape));
    }


    template <class Type>
    void
    ArrayNDTest<Type>::
    testConstructor__string()
    {
      // This constructor should initialize according to the specified
      // string, which describes an array of m_defaultArraySize elements
      // having values matching those in m_fibonacciCArray.
      ArrayND<3, Type> array0(m_fibonacciString);
      BRICK_TEST_ASSERT(this->isShapeMatch(array0.getShape(),
                                           m_defaultArrayShape));
      BRICK_TEST_ASSERT(std::equal(array0.begin(), array0.end(),
                                 m_fibonacciCArray));

      // Construct using a malformed string.
      try {
        ArrayND<3, Type> array1(m_illegalString);
        BRICK_TEST_ASSERT(false);
      } catch(brick::common::ValueException const&) {
        // pass.
      }
    }


    template <class Type>
    void
    ArrayNDTest<Type>::
    testConstructor__ArrayND()
    {
      // This constructor should make an array which references the same
      // data as its argument.
      ArrayND<3, Type>* array0Ptr = new ArrayND<3, Type>(m_fibonacciString);
      ArrayND<3, Type> array1(*array0Ptr);
      BRICK_TEST_ASSERT(array1.data() == array0Ptr->data());
      BRICK_TEST_ASSERT(
        this->isShapeMatch(array1.getShape(), array0Ptr->getShape()));
      BRICK_TEST_ASSERT(
        array1.getShape().data() != array0Ptr->getShape().data());
      delete array0Ptr;

      // Normally, deleting array0Ptr would invalidate the data it
      // points to.  If reference counting is working, however, the data
      // should remain valid until after array1 goes out of scope.
      BRICK_TEST_ASSERT(std::equal(array1.begin(), array1.end(),
                                 m_fibonacciCArray));
    }


    template <class Type>
    void
    ArrayNDTest<Type>::
    testDestructor()
    {
      // No independent test for destructor.
    }


    template <class Type>
    void
    ArrayNDTest<Type>::
    testBegin()
    {
      // Member function begin() should return an iterator pointing to
      // the first element of the array.
      ArrayND<3, Type> array0(m_defaultArrayShape);
      BRICK_TEST_ASSERT(&(*(array0.begin())) == array0.data());
    }


    template <class Type>
    void
    ArrayNDTest<Type>::
    testBeginConst()
    {
      // Member function begin() should return an iterator pointing to
      // the first element of the array.
      const ArrayND<3, Type> array0(m_defaultArrayShape);
      BRICK_TEST_ASSERT(&(*(array0.begin())) == array0.data());
    }


    template <class Type>
    void
    ArrayNDTest<Type>::
    testClear()
    {
      ArrayND<3, Type> array0(m_fibonacciString);
      array0.clear();
      BRICK_TEST_ASSERT(array0.size() == 0);
    }


    template <class Type>
    void
    ArrayNDTest<Type>::
    testCopy()
    {
      // Member function copy() should allocate a new array and deep
      // copy the contents of *this.
      ArrayND<3, Type> array0(m_fibonacciString);
      ArrayND<3, Type> array1 = array0.copy();
      BRICK_TEST_ASSERT(this->isShapeMatch(array1.getShape(), array0.getShape()));
      BRICK_TEST_ASSERT(std::equal(array1.begin(), array1.end(), array0.begin()));
      BRICK_TEST_ASSERT(array1.data() != array0.data());
    }


    template <class Type>
    void
    ArrayNDTest<Type>::
    testCopy__ArrayND()
    {
      ArrayND<3, Type> array0(m_fibonacciString);

      // Deep copy the shape so we can modify it with impunity.
      Array1D<size_t> arrayShape = array0.getShape().copy();
      ArrayND<3, Type> array1(arrayShape);
      arrayShape[0] += 1;
      ArrayND<3, Type> array2(arrayShape);

      // Test that the copy goes OK.
      array1.copy(array0);
      BRICK_TEST_ASSERT(
        this->isShapeMatch(array1.getShape(), array0.getShape()));
      BRICK_TEST_ASSERT(
        std::equal(array1.begin(), array1.end(), array0.begin()));
      BRICK_TEST_ASSERT(array1.data() != array0.data());

      // Test that size is checked.
      BRICK_TEST_ASSERT_EXCEPTION(brick::common::IndexException,
                                  array2.copy(array0));
    }


    template <class Type>
    void
    ArrayNDTest<Type>::
    testCopy__Type2Ptr()
    {
      // This method simply copies elements from its argument.
      // Begin by making an appropriately sized array with known contents.
      ArrayND<3, Type> array0(m_defaultArrayShape);
      array0 = 0;

      array0.copy(m_fibonacciCArray);
      BRICK_TEST_ASSERT(this->isShapeMatch(array0.getShape(), m_defaultArrayShape));
      BRICK_TEST_ASSERT(std::equal(array0.begin(), array0.end(),
                                 m_fibonacciCArray));
      BRICK_TEST_ASSERT(array0.data() != m_fibonacciCArray);
    }


    template <class Type>
    void
    ArrayNDTest<Type>::
    testData()
    {
      // This method returns a pointer to the internal data store.
      ArrayND<3, Type> array0(m_defaultArrayShape);
      BRICK_TEST_ASSERT(array0.data() == &(array0[0]));
    }


    template <class Type>
    void
    ArrayNDTest<Type>::
    testDataConst()
    {
      // This method returns a pointer to the internal data store.
      ArrayND<3, Type> array0(m_defaultArrayShape);
      for(size_t index0 = 0; index0 < array0.size(); ++index0) {
        array0[index0] = static_cast<Type>(index0 + 15);
      }
      const ArrayND<3, Type> array1(array0);
      BRICK_TEST_ASSERT(*array1.data() == array1[0]);
    }


    template <class Type>
    void
    ArrayNDTest<Type>::
    testEnd()
    {
      // Member function end() should return an iterator pointing to
      // the last element of the array.
      ArrayND<3, Type> array0(m_defaultArrayShape);
      size_t total = std::accumulate(m_defaultArrayShape.begin(),
                                     m_defaultArrayShape.end(),
                                     1.0,
                                     std::multiplies<size_t>());
      Type* finalElementPtr = &(array0[total - 1]) + 1;
      BRICK_TEST_ASSERT(&(*(array0.end())) == finalElementPtr);
    }


    template <class Type>
    void
    ArrayNDTest<Type>::
    testEndConst()
    {
      // Member function end() should return an iterator pointing to
      // the last element of the array.
      const ArrayND<3, Type> array0(m_defaultArrayShape);
      size_t total = std::accumulate(m_defaultArrayShape.begin(),
                                     m_defaultArrayShape.end(),
                                     1.0,
                                     std::multiplies<size_t>());
      BRICK_TEST_ASSERT(array0.end() == array0.begin() + total);
    }


    template <class Type>
    void
    ArrayNDTest<Type>::
    testReinit()
    {
      // This method changes the shape of the array and reallocates
      // storage.
      ArrayND<3, Type> array0(m_defaultArrayShape);
      ArrayND<3, Type> array1 = array0;
      Type* dataPtr0 = array0.data();
      Type* dataPtr1 = array1.data();

      BRICK_TEST_ASSERT(dataPtr0 == dataPtr1);

      array0.reinit(m_defaultArrayShape);
      BRICK_TEST_ASSERT(this->isShapeMatch(
                        array0.getShape(), m_defaultArrayShape));
      BRICK_TEST_ASSERT(array0.data() != dataPtr0);
    }


    template <class Type>
    void
    ArrayNDTest<Type>::
    testSize()
    {
      ArrayND<3, Type> array0;
      ArrayND<3, Type> array1(m_defaultArrayShape);
      size_t total = std::accumulate(m_defaultArrayShape.begin(),
                                     m_defaultArrayShape.end(),
                                     1.0,
                                     std::multiplies<size_t>());
      BRICK_TEST_ASSERT(array0.size() == 0);
      BRICK_TEST_ASSERT(array1.size() == total);
    }


    template <class Type>
    void
    ArrayNDTest<Type>::
    testAssignmentOperator__ArrayND()
    {
      // The assignment operator shallow copies the contents of source.
      ArrayND<3, Type> array0(m_fibonacciString);
      ArrayND<3, Type> array1;
      array1 = array0;
      BRICK_TEST_ASSERT(this->isShapeMatch(array0.getShape(), m_defaultArrayShape));
      BRICK_TEST_ASSERT(this->isShapeMatch(array0.getShape(), array1.getShape()));
      BRICK_TEST_ASSERT(array1.data() == array0.data());

      // Make sure the data wasn't changed by the operation.
      BRICK_TEST_ASSERT(std::equal(array1.begin(), array1.end(),
                                 m_fibonacciCArray));
    }


    template <class Type>
    void
    ArrayNDTest<Type>::
    testAssignmentOperator__Type()
    {
      // This member function assigns the same value to every element in
      // the array.
      ArrayND<3, Type> array0(m_defaultArrayShape);
      Type* dataPtr = array0.data();
      array0 = static_cast<Type>(m_defaultArrayValue);

      // Check that the size and location of the array wasn't changed.
      BRICK_TEST_ASSERT(this->isShapeMatch(array0.getShape(), m_defaultArrayShape));
      BRICK_TEST_ASSERT(array0.data() == dataPtr);

      // Check that all values were set by getting a pointer to the
      //
      // Strangely, std::find_if segfaults for me under g++-4.1.
      //
      // first one that's not equal to m_defaultArrayValue
      //   typename ArrayND<3, Type>::iterator firstRenegade =
      //     std::find_if(array0.begin(), array0.end(),
      //                  std::bind2nd(std::not_equal_to<Type>(),
      //                               m_defaultArrayValue));

      typename ArrayND<3, Type>::iterator firstRenegade = array0.begin();
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
    ArrayNDTest<Type>::
    testApplicationOperator()
    {
      // This member function returns the (index)th element of the array
      // by reference.
      ArrayND<3, Type> array0(m_fibonacciString);
      for(size_t index = 0; index < array0.size(); ++index) {
        BRICK_TEST_ASSERT(array0(index) == m_fibonacciCArray[index]);
      }
    }


    template <class Type>
    void
    ArrayNDTest<Type>::
    testApplicationOperatorConst()
    {
      // This member function returns the (index)th element of the array
      // by reference.
      const ArrayND<3, Type> array0(m_fibonacciString);
      for(size_t index = 0; index < array0.size(); ++index) {
        BRICK_TEST_ASSERT(array0(index) == m_fibonacciCArray[index]);
      }
    }


    template <class Type>
    void
    ArrayNDTest<Type>::
    testIndexOperator()
    {
      // This member function returns the (index)th element of the array
      // by reference.
      ArrayND<3, Type> array0(m_fibonacciString);
      for(size_t index = 0; index < array0.size(); ++index) {
        BRICK_TEST_ASSERT(array0[index] == m_fibonacciCArray[index]);
      }
    }


    template <class Type>
    void
    ArrayNDTest<Type>::
    testIndexOperatorConst()
    {
      // This member function returns the (index)th element of the array
      // by reference.
      const ArrayND<3, Type> array0(m_fibonacciString);
      for(size_t index = 0; index < array0.size(); ++index) {
        BRICK_TEST_ASSERT(array0[index] == m_fibonacciCArray[index]);
      }
    }


    template<class Type>
    bool
    ArrayNDTest<Type>::
    isShapeMatch(Array1D<size_t> const& shape0, Array1D<size_t> const& shape1)
    {
      if(shape0.size() != shape1.size()) {
        return false;
      }
      return std::equal(shape0.begin(), shape0.end(), shape1.begin());
    }

  } // namespace numeric

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::numeric::ArrayNDTest<double> currentTest0("double");
  brick::numeric::ArrayNDTest<float> currentTest1("float");
  brick::numeric::ArrayNDTest<int> currentTest2("int");
  // brick::numeric::ArrayNDTest<size_t> currentTest3("size_t");
  bool result = (currentTest0.run()
                 && currentTest1.run()
                 && currentTest2.run()
                 && currentTest3.run());
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::ArrayNDTest<double> currentTest0("double");
  brick::numeric::ArrayNDTest<float> currentTest1("float");
  brick::numeric::ArrayNDTest<int> currentTest2("int");
  // brick::numeric::ArrayNDTest<size_t> currentTest3("size_t"); /* Test assumes signed elements */

}

#endif
