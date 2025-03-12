/**
***************************************************************************
* @file arrayTestCommon.h
* Source file defining ArrayTestCommon class template.
*
* Copyright (C) 2004-2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <cmath>
#include <sstream>
#include <iomanip>
// xxx #include <brick/common/functional.hh>
// xxx #include <brick/numeric/utilities.hh>
#include <brick/test/functors.hh>
#include <brick/test/testFixture.hh>

namespace brick {

  /**
   ** This class template automates many of the routine tests of the
   ** dlrLibs array classes.
   **/
  template <class FixtureType, class ArrayType, class ComparisonResultType>
  class ArrayTestCommon : public test::TestFixture< FixtureType > {

  public:

    typedef ArrayType array_type;
    typedef FixtureType TestFixtureType;

    // Constructor
    ArrayTestCommon(const std::string& testFixtureName);

    // Destructor
    virtual ~ArrayTestCommon() {}

    // Pure virtual members must be overridden by a child class.
    virtual void
    checkShapeEquality(const ArrayType& array0,
                       const ArrayType& array1) = 0;

    virtual void
    checkShapeEquality(const ArrayType& array0,
                       const ComparisonResultType& array1) = 0;

    virtual void
    checkValueEquality(const ArrayType& array0,
                       const ArrayType& array1,
                       typename ArrayType::value_type tolerance) = 0;

    virtual typename ArrayType::value_type
    getComparisonOperatorThreshold() = 0;

    virtual typename ArrayType::value_type
    getEqualityOperatorTarget() = 0;

    virtual ArrayType
    getFibonacciArray() = 0;

    virtual typename ArrayType::value_type
    getIncrementOperatorArgument() = 0;

    virtual typename ArrayType::value_type
    getMultiplicationOperatorArgument() = 0;

    virtual ArrayType
    getSquaresArray() = 0;


    // Tests of member functions.
    void testOperatorPlusEquals__ArrayND();
    void testOperatorPlusEquals__Type();
    void testOperatorMinusEquals__ArrayND();
    void testOperatorMinusEquals__Type();
    void testOperatorTimesEquals__ArrayND();
    void testOperatorTimesEquals__Type();
    void testOperatorDividedByEquals__ArrayND();
    void testOperatorDividedByEquals__Type();


    // Tests of non-member functions.
    void testOperatorEqualEqual__Type();
    void testOperatorEqualEqual__ArrayND();
    void testOperatorGreaterThan__Type();
    void testOperatorGreaterThanOrEqualTo__Type();
    void testOperatorLessThan__Type();
    void testOperatorLessThanOrEqualTo__Type();
    void testOperatorPlus__ArrayND__ArrayND();
    void testOperatorMinus__ArrayND__ArrayND();
    void testOperatorTimes__ArrayND__ArrayND();
    void testOperatorDividedBy__ArrayND__ArrayND();
    void testOperatorPlus__ArrayND__Type();
    void testOperatorMinus__ArrayND__Type();
    void testOperatorTimes__ArrayND__Type();
    void testOperatorDividedBy__ArrayND__Type();
    void testOperatorPlus__Type__ArrayND();
    // void testOperatorMinus__Type__ArrayND();
    void testOperatorTimes__Type__ArrayND();
    // void testOperatorDividedBy__Type__ArrayND();
    void testOutputOperator();
    void testInputOperator();

  protected:

    // Currently no protected members.

  }; // class ArrayTestCommon


  /* ============== Locally defined classes ============== */

  namespace {

    template<class Type0, class Type1, class Type2>
    class DividedByFunctor {
    public:
      using result_type = Type2;
      Type2 operator()(Type0& arg0, const Type1& arg1) {
        return arg0 / arg1;
      }
    };


    template<class Type0, class Type1>
    class DividedByEqualsFunctor {
    public:
      using result_type = Type0;
      Type0 operator()(Type0& arg0, const Type1& arg1) {
        return arg0 /= arg1;
      }
    };


    template<class Type0, class Type1, class Type2>
    class EqualEqualFunctor {
    public:
      using result_type = Type2;
      Type2 operator()(const Type0& arg0, const Type1& arg1) {
        return arg0 == arg1;
      }
    };


    template<class Type0, class Type1, class Type2>
    class GreaterThanFunctor {
    public:
      using result_type = Type2;
      Type2 operator()(const Type0& arg0, const Type1& arg1) {
        return arg0 > arg1;
      }
    };


    template<class Type0, class Type1, class Type2>
    class GreaterThanOrEqualToFunctor {
    public:
      using result_type = Type2;
      Type2 operator()(const Type0& arg0, const Type1& arg1) {
        return arg0 >= arg1;
      }
    };


    template<class Type0, class Type1, class Type2>
    class LessThanFunctor {
    public:
      using result_type = Type2;
      Type2 operator()(const Type0& arg0, const Type1& arg1) {
        return arg0 < arg1;
      }
    };


    template<class Type0, class Type1, class Type2>
    class LessThanOrEqualToFunctor {
    public:
      using result_type = Type2;
      Type2 operator()(const Type0& arg0, const Type1& arg1) {
        return arg0 <= arg1;
      }
    };


    template<class Type0, class Type1, class Type2>
    class MinusFunctor {
    public:
      using result_type = Type2;
      Type2 operator()(Type0& arg0, const Type1& arg1) {
        return arg0 - arg1;
      }
    };


    template<class Type0, class Type1>
    class MinusEqualsFunctor {
    public:
      using result_type = Type0;
      Type0 operator()(Type0& arg0, const Type1& arg1) {
        return arg0 -= arg1;
      }
    };


    template<class Type0, class Type1, class Type2>
    class PlusFunctor {
    public:
      using result_type = Type2;
      Type2 operator()(Type0& arg0, const Type1& arg1) {
        return arg0 + arg1;
      }
    };


    template<class Type0, class Type1>
    class PlusEqualsFunctor {
    public:
      using result_type = Type0;
      Type0 operator()(Type0& arg0, const Type1& arg1) {
        return arg0 += arg1;
      }
    };


    template<class Type0, class Type1, class Type2>
    class TimesFunctor {
    public:
      using result_type = Type2;
      Type2 operator()(Type0& arg0, const Type1& arg1) {
        return arg0 * arg1;
      }
    };


    template<class Type0, class Type1>
    class TimesEqualsFunctor {
    public:
      using result_type = Type0;
      Type0 operator()(Type0& arg0, const Type1& arg1) {
        return arg0 *= arg1;
      }
    };


    /* ============== Locally defined functions ============== */

    template <class NumericType>
    NumericType
    getTolerance(NumericType) {
      return static_cast<NumericType>(0);
    }

    // Adaptive tolerance is necessary to avoid trouble with float
    // comparisons below.  Note that IEEE 32-bit float has 23-bit
    // mantissa.
    template<>
    common::Float32
    getTolerance(common::Float32 scale)
    {
      if(scale < 1000.0) {
        return 1.0E-6;
      }
      common::Float32 exponent = static_cast<common::Float32>(
        std::log(scale) / std::log(2.0));
      return static_cast<common::Float32>(
        10.0 * std::pow((double)2.0, (double)(exponent - 23)));
    }

    // Adaptive tolerance is necessary to avoid trouble with float
    // comparisons below.  Note that IEEE 64-bit float has 52-bit
    // mantissa.
    template<>
    common::Float64
    getTolerance(common::Float64 scale)
    {
      if(scale < 1000.0) {
        return 1.0E-12;
      }
      common::Float64 exponent = static_cast<common::Float32>(
        std::log(scale) / std::log(2.0));
      return static_cast<common::Float32>(
        10.0 * std::pow((double)2.0, (double)(exponent - 52)));
    }



    /**
     * This function template simply dispatches to
     * brick::test::approximatelyEqual().  We use it instead of using
     * brick::test::approximatelyEqual() directly because we can specialize
     * this function to handle funky cases, such as comparing two
     * bools.
     *
     * @param argument0 This argument is the first element to be compared.
     *
     * @param argument1 This argument is the second element to be compared.
     *
     * @return The return value is true if the two elements are just
     * about equal, false otherwise.
     */
    template <class Type>
    bool
    testApproximatelyEqual(Type argument0, Type argument1) {
      return test::approximatelyEqual(
        argument0, argument1, static_cast<Type>(1.0E-14));
    }


    /**
     * This function specializes testApproximatelyEqual() for
     * arguments of type bool.
     *
     * @param argument0 This argument is the first element to be compared.
     *
     * @param argument1 This argument is the second element to be compared.
     *
     * @return The return value is true if the two elements are equal,
     * false otherwise.
     */
    template <>
    bool
    testApproximatelyEqual(bool argument0, bool argument1) __attribute__((__unused__));
    template <>
    bool
    testApproximatelyEqual(bool argument0, bool argument1) {
      return argument0 == argument1;
    }


    /**
     * This function is a helper to avoid code duplication in operator
     * tests.  We'd like to make this function be a member of
     * ArrayTestCommon, but instead we make it external because this is a
     * little easier on the compiler.
     *
     * @param arrayFunctor This argument will be applied to an ArrayND
     * instance and an element instance.
     *
     * @param elementFunctor This argument will be applied to element
     * instances, mimicing the effect of arrayFunctor.
     */
    template <class TestClass, class Functor0, class Functor1>
    void
    testOperatorX__ArrayND__ArrayND(TestClass& testInstance,
                                    Functor0 arrayFunctor,
                                    Functor1 elementFunctor) {
      // Create the arguments.
      typedef typename TestClass::array_type ArrayType;
      ArrayType array0 = testInstance.getSquaresArray();
      ArrayType array1 = testInstance.getFibonacciArray();

      // Apply the operator.
      typename Functor0::result_type result = arrayFunctor(array0, array1);

      // Get the data again, just in case it was screwed up somehow
      // by the operator.
      array0 = testInstance.getSquaresArray();
      array1 = testInstance.getFibonacciArray();

      // Verify result.
      // Note: it's very easy to get bitten by the x86 register
      // precision snafu here.  Basically, when compiled with
      // optimizations, the reference and computed values may have
      // different precisions, leading to unexpected answers!  We
      // solve this by having an adaptive tolerance for 32bit floats.
      testInstance.checkShapeEquality(array0, result);
      for(size_t index0 = 0; index0 < result.size(); ++index0) {
        typename Functor1::result_type targetValue = elementFunctor(
          array0[index0], array1[index0]);
        typename Functor1::result_type tolerance =
          getTolerance<typename Functor1::result_type>(targetValue);

	bool equalFlag = test::approximatelyEqual(result[index0], targetValue,
					    tolerance);
	if(!equalFlag) {
	  std::cout << "testOperatorX__ArrayND__ArrayND() fail: "
		    << result[index0] << " vs. " << targetValue
		    << " with tolerance of " << tolerance << std::endl;
	}
	BRICK_TEST_ASSERT(equalFlag);
      }
    }


    /**
     * This member function is a helper to avoid code duplication in
     * operator tests.  We'd like to make this function be a member of
     * ArrayTestCommon, but instead we make it external because this is a
     * little easier on the compiler.
     *
     * @param arrayFunctor This argument will be applied to an ArrayND
     * instance and an element instance.
     *
     * @param elementFunctor This argument will be applied to element
     * instances, mimicing the effect of arrayFunctor.
     */
    template <class TestClass, class Functor0, class Functor1,
              class ScalarType>
    void
    testOperatorX__ArrayND__Type(TestClass& testInstance,
                                 Functor0 arrayFunctor,
                                 Functor1 elementFunctor,
                                 ScalarType argument) {
      // Get the arguments.
      typedef typename TestClass::array_type ArrayType;
      ArrayType array0 = testInstance.getSquaresArray();

      // Apply the operator.
      typename Functor0::result_type result = arrayFunctor(array0, argument);

      // Get the data again, just in case it was screwed up somehow
      // by the operator.
      array0 = testInstance.getSquaresArray();

      // Verify result.
      // Note: it's very easy to get bitten by the x86 register
      // precision snafu here.  Basically, when compiled with
      // optimizations, the reference and computed values may have
      // different precisions, leading to unexpected answers!  We
      // solve this by having an adaptive tolerance for 32bit floats.
      testInstance.checkShapeEquality(array0, result);
      for(size_t index0 = 0; index0 < result.size(); ++index0) {
        typename Functor1::result_type targetValue = elementFunctor(
          array0[index0], argument);
        typename Functor1::result_type tolerance =
          getTolerance<typename Functor1::result_type>(targetValue);

        BRICK_TEST_ASSERT(
          test::approximatelyEqual(result[index0], targetValue, tolerance));
        // BRICK_TEST_ASSERT(result[index0] == targetValue);
      }
    }


    /**
     * This member function is a helper to avoid code duplication in
     * operator tests.  We'd like to make this function be a member of
     * ArrayTestCommon, but instead we make it external because this is a
     * little easier on the compiler.
     *
     * @param arrayFunctor This argument will be applied to an ArrayND
     * instance and an element instance.
     *
     * @param elementFunctor This argument will be applied to element
     * instances, mimicing the effect of arrayFunctor.
     */
    template <class TestClass, class Functor0, class Functor1,
              class ScalarType>
    void
    testOperatorX__Type__ArrayND(TestClass& testInstance,
                                 Functor0 arrayFunctor,
                                 Functor1 elementFunctor,
                                 ScalarType argument) {
      // Get the arguments.
      typedef typename TestClass::array_type ArrayType;
      ArrayType array0 = testInstance.getSquaresArray();

      // Apply the operator.
      typename Functor0::result_type result = arrayFunctor(argument, array0);

      // Get the data again, just in case it was screwed up somehow
      // by the operator.
      array0 = testInstance.getSquaresArray();

      // Verify result.
      testInstance.checkShapeEquality(array0, result);
      typename ArrayType::value_type tolerance =
        static_cast<typename ArrayType::value_type>(1.0e-14);
      for(size_t index0 = 0; index0 < result.size(); ++index0) {
        typename Functor1::result_type targetValue = elementFunctor(
          argument, array0[index0]);

	bool equalFlag = test::approximatelyEqual(result[index0], targetValue,
					    tolerance);
	if(!equalFlag) {
	  std::cout << "testOperatorX__Type__ArrayND() fail: "
		    << result[index0] << " vs. " << targetValue
		    << " with tolerance of " << tolerance << std::endl;
	}
	BRICK_TEST_ASSERT(equalFlag);
      }
    }



    /**
     * This member function is a helper to avoid code duplication in
     * operator tests.  We'd like to make this function be a member of
     * ArrayTestCommon, but instead we make it external because this is a
     * little easier on the compiler.
     *
     * @param arrayFunctor This argument will be applied to ArrayND
     * instances.
     *
     * @param elementFunctor This argument will be applied to element
     * instances, mimicing the effect of arrayFunctor.
     */
    template <class TestClass, class Functor0, class Functor1>
    void
    testOperatorXEquals__ArrayND(TestClass& testInstance,
                                 Functor0 arrayFunctor,
                                 Functor1 elementFunctor) {
      // Create the arguments.
      typedef typename TestClass::array_type ArrayType;
      ArrayType array0 = testInstance.getSquaresArray();
      ArrayType array0Copy = testInstance.getSquaresArray();
      ArrayType array1 = testInstance.getFibonacciArray();

      // Preliminary recordkeeping.
      typename ArrayType::value_type* dataPtr = array0.data();

      // Apply the operator.  This should modify array0.
      arrayFunctor(array0, array1);

      // Get the data again, just in case it was screwed up somehow
      // by the operator.
      array1 = testInstance.getFibonacciArray();

      // Verify result.  Shape and storage of array0 shouldn't change.
      testInstance.checkShapeEquality(array0, array1);
      BRICK_TEST_ASSERT(array0.data() == dataPtr);

      // Note: it's very easy to get bitten by the x86 register
      // precision snafu here.  Basically, when compiled with
      // optimizations, the reference and computed values may have
      // different precisions, leading to unexpected answers!  We
      // solve this by having an adaptive tolerance for 32bit floats.
      for(size_t index0 = 0; index0 < array0.size(); ++index0) {
        // typename Functor1::result_type targetValue = elementFunctor(
        //   array0Copy[index0], array1[index0]);
        elementFunctor(array0Copy[index0], array1[index0]);
        typename ArrayType::value_type tolerance =
          getTolerance<typename ArrayType::value_type>(array0Copy[index0]);
        BRICK_TEST_ASSERT(test::approximatelyEqual(array0[index0], array0Copy[index0],
                                           tolerance));
      }
    }


    /**
     * This member function is a helper to avoid code duplication in
     * operator tests.  We'd like to make this function be a member of
     * ArrayTestCommon, but instead we make it external because this is a
     * little easier on the compiler.
     *
     * @param arrayFunctor This argument will be applied to an ArrayND
     * instance and an element instance.
     *
     * @param elementFunctor This argument will be applied to element
     * instances, mimicing the effect of arrayFunctor.
     */
    template <class TestClass, class Functor0, class Functor1,
              class ScalarType>
    void
    testOperatorXEquals__Type(TestClass& testInstance,
                              Functor0 arrayFunctor,
                              Functor1 elementFunctor,
                              ScalarType argument) {
      // Create the arguments.
      typedef typename TestClass::array_type ArrayType;
      ArrayType array0 = testInstance.getFibonacciArray();

      // Preliminary recordkeeping.
      typename ArrayType::value_type* dataPtr = array0.data();

      // Apply the operator.
      arrayFunctor(array0, argument);

      // Get another array, for size checking.
      ArrayType array1 = testInstance.getFibonacciArray();

      // Verify result.
      testInstance.checkShapeEquality(array0, array1);
      BRICK_TEST_ASSERT(array0.data() == dataPtr);

      // Note: it's very easy to get bitten by the x86 register
      // precision snafu here.  Basically, when compiled with
      // optimizations, the reference and computed values may have
      // different precisions, leading to unexpected answers!  We
      // solve this by having an adaptive tolerance for 32bit floats.
      for(size_t index0 = 0; index0 < array0.size(); ++index0) {
        typename Functor1::result_type targetValue = (
          elementFunctor(array0[index0], argument));
        typename Functor1::result_type tolerance = (
          getTolerance<typename ArrayType::value_type>(targetValue));
        BRICK_TEST_ASSERT(
          test::approximatelyEqual(array0[index0], targetValue, tolerance));
      }
    }


  }


  /* ============== Member Function Definititions ============== */

  template <class FixtureType, class ArrayType, class ComparisonResultType>
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  ArrayTestCommon(const std::string& testFixtureName)
    : test::TestFixture<FixtureType>(testFixtureName)
  {
    // // Register all tests.
    // Tests of member functions.
    BRICK_TEST_REGISTER_MEMBER(testOperatorPlusEquals__ArrayND);
    BRICK_TEST_REGISTER_MEMBER(testOperatorPlusEquals__Type);
    BRICK_TEST_REGISTER_MEMBER(testOperatorMinusEquals__ArrayND);
    BRICK_TEST_REGISTER_MEMBER(testOperatorMinusEquals__Type);
    BRICK_TEST_REGISTER_MEMBER(testOperatorTimesEquals__ArrayND);
    BRICK_TEST_REGISTER_MEMBER(testOperatorTimesEquals__Type);
    BRICK_TEST_REGISTER_MEMBER(testOperatorDividedByEquals__ArrayND);
    BRICK_TEST_REGISTER_MEMBER(testOperatorDividedByEquals__Type);

    // Tests of non-member functions.
    BRICK_TEST_REGISTER_MEMBER(testOperatorEqualEqual__Type);
    BRICK_TEST_REGISTER_MEMBER(testOperatorEqualEqual__ArrayND);
    BRICK_TEST_REGISTER_MEMBER(testOperatorGreaterThan__Type);
    BRICK_TEST_REGISTER_MEMBER(testOperatorGreaterThanOrEqualTo__Type);
    BRICK_TEST_REGISTER_MEMBER(testOperatorLessThan__Type);
    BRICK_TEST_REGISTER_MEMBER(testOperatorLessThanOrEqualTo__Type);
    BRICK_TEST_REGISTER_MEMBER(testOperatorPlus__ArrayND__ArrayND);
    BRICK_TEST_REGISTER_MEMBER(testOperatorMinus__ArrayND__ArrayND);
    BRICK_TEST_REGISTER_MEMBER(testOperatorTimes__ArrayND__ArrayND);
    BRICK_TEST_REGISTER_MEMBER(testOperatorDividedBy__ArrayND__ArrayND);
    BRICK_TEST_REGISTER_MEMBER(testOperatorPlus__ArrayND__Type);
    BRICK_TEST_REGISTER_MEMBER(testOperatorMinus__ArrayND__Type);
    BRICK_TEST_REGISTER_MEMBER(testOperatorTimes__ArrayND__Type);
    BRICK_TEST_REGISTER_MEMBER(testOperatorDividedBy__ArrayND__Type);
    BRICK_TEST_REGISTER_MEMBER(testOperatorPlus__Type__ArrayND);
    // BRICK_TEST_REGISTER_MEMBER(testOperatorMinus__Type__ArrayND);
    BRICK_TEST_REGISTER_MEMBER(testOperatorTimes__Type__ArrayND);
    // BRICK_TEST_REGISTER_MEMBER(testOperatorDividedBy__Type__ArrayND);
    BRICK_TEST_REGISTER_MEMBER(testOutputOperator);
    BRICK_TEST_REGISTER_MEMBER(testInputOperator);
  }


  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testOperatorPlusEquals__ArrayND()
  {
    typedef typename ArrayType::value_type ElementType;

    typedef
      PlusEqualsFunctor< ArrayType, ArrayType >
      ArrayFunctor;

    typedef
      PlusEqualsFunctor<
        ElementType,
        ElementType >
      ElementFunctor;

    testOperatorXEquals__ArrayND(
      *this, ArrayFunctor(), ElementFunctor());
  }


  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testOperatorPlusEquals__Type()
  {
    typedef typename ArrayType::value_type ElementType;
    ElementType increment = this->getIncrementOperatorArgument();

    typedef
      PlusEqualsFunctor<
        ArrayType,
        ElementType >
      ArrayFunctor;

    typedef
      PlusEqualsFunctor<
        ElementType,
        ElementType >
      ElementFunctor;

    testOperatorXEquals__Type(
      *this, ArrayFunctor(), ElementFunctor(), increment);
  }


  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testOperatorMinusEquals__ArrayND()
  {
    typedef typename ArrayType::value_type ElementType;

    typedef
      MinusEqualsFunctor< ArrayType, ArrayType >
      ArrayFunctor;

    typedef
      MinusEqualsFunctor<
        ElementType,
        ElementType >
      ElementFunctor;

    testOperatorXEquals__ArrayND(
      *this, ArrayFunctor(), ElementFunctor());
  }


  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testOperatorMinusEquals__Type()
  {
    // This member function decrements each element by the value of the
    // corresponding element of its argument.
    typedef typename ArrayType::value_type ElementType;
    ElementType increment = this->getIncrementOperatorArgument();

    typedef
      MinusEqualsFunctor<
        ArrayType,
        ElementType >
      ArrayFunctor;

    typedef
      MinusEqualsFunctor<
        ElementType,
        ElementType >
      ElementFunctor;

    testOperatorXEquals__Type(
      *this, ArrayFunctor(), ElementFunctor(), increment);
  }


  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testOperatorTimesEquals__ArrayND()
  {
    // This member function increments each element by the value of the
    // corresponding element of its argument.

    typedef typename ArrayType::value_type ElementType;

    typedef
      TimesEqualsFunctor< ArrayType, ArrayType >
      ArrayFunctor;

    typedef
      TimesEqualsFunctor<
        ElementType,
        ElementType >
      ElementFunctor;

    testOperatorXEquals__ArrayND(
      *this, ArrayFunctor(), ElementFunctor());
  }


  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testOperatorTimesEquals__Type()
  {
    // This member function increments each element by the value of the
    // corresponding element of its argument.
    typedef typename ArrayType::value_type ElementType;
    ElementType multiplier = this->getMultiplicationOperatorArgument();

    typedef
      TimesEqualsFunctor<
        ArrayType,
        ElementType >
      ArrayFunctor;

    typedef
      TimesEqualsFunctor<
        ElementType,
        ElementType >
      ElementFunctor;

    testOperatorXEquals__Type(
      *this, ArrayFunctor(), ElementFunctor(), multiplier);
  }


  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testOperatorDividedByEquals__ArrayND()
  {
    // This member function divides each element by the value of the
    // corresponding element of its argument.

    typedef typename ArrayType::value_type ElementType;

    typedef
      DividedByEqualsFunctor< ArrayType, ArrayType >
      ArrayFunctor;

    typedef
      DividedByEqualsFunctor<
        ElementType,
        ElementType >
      ElementFunctor;

    testOperatorXEquals__ArrayND(
      *this, ArrayFunctor(), ElementFunctor());
  }


  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testOperatorDividedByEquals__Type()
  {
    // This member function divides each element by the value of the
    // corresponding element of its argument.
    typedef typename ArrayType::value_type ElementType;
    ElementType multiplier = this->getMultiplicationOperatorArgument();

    typedef
      DividedByEqualsFunctor<
        ArrayType,
        ElementType >
      ArrayFunctor;

    typedef
      DividedByEqualsFunctor<
        ElementType,
        ElementType >
      ElementFunctor;

    testOperatorXEquals__Type(
      *this, ArrayFunctor(), ElementFunctor(), multiplier);
  }


  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testOperatorEqualEqual__Type()
  {
    // This member compares the elements of the array with a constant.
    typedef typename ArrayType::value_type ElementType;

    ElementType target = this->getEqualityOperatorTarget();

    typedef
      EqualEqualFunctor< ArrayType,
                         ElementType,
                         ComparisonResultType >
      ArrayFunctor;

    typedef EqualEqualFunctor< ElementType,
                               ElementType,
                               typename ComparisonResultType::value_type >
      ElementFunctor;

    testOperatorX__ArrayND__Type(
      *this, ArrayFunctor(), ElementFunctor(), target);
  }


  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testOperatorEqualEqual__ArrayND()
  {
    typedef typename ArrayType::value_type ElementType;

    typedef
      EqualEqualFunctor< ArrayType,
                         ArrayType,
                         ComparisonResultType >
      ArrayFunctor;

    typedef EqualEqualFunctor< ElementType,
                               ElementType,
                               typename ComparisonResultType::value_type >
      ElementFunctor;

    testOperatorX__ArrayND__ArrayND(
      *this, ArrayFunctor(), ElementFunctor());
  }


  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testOperatorGreaterThan__Type()
  {
    // This member compares the elements of the array with a constant.
    typedef typename ArrayType::value_type ElementType;

    ElementType threshold = this->getComparisonOperatorThreshold();

    typedef
      GreaterThanFunctor< ArrayType,
                          ElementType,
                          ComparisonResultType >
      ArrayFunctor;

    typedef GreaterThanFunctor< ElementType,
                                ElementType,
                                typename ComparisonResultType::value_type >
      ElementFunctor;

    testOperatorX__ArrayND__Type(
      *this, ArrayFunctor(), ElementFunctor(), threshold);
  }



  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testOperatorGreaterThanOrEqualTo__Type()
  {
    // This member compares the elements of the array with a constant.
    typedef typename ArrayType::value_type ElementType;

    ElementType threshold = this->getComparisonOperatorThreshold();

    typedef
      GreaterThanOrEqualToFunctor<
        ArrayType,
        ElementType,
        ComparisonResultType >
      ArrayFunctor;

    typedef
      GreaterThanOrEqualToFunctor<
        ElementType,
        ElementType,
        typename ComparisonResultType::value_type >
      ElementFunctor;

    testOperatorX__ArrayND__Type(
      *this, ArrayFunctor(), ElementFunctor(), threshold);
  }


  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testOperatorLessThan__Type()
  {
    // This member compares the elements of the array with a constant.
    typedef typename ArrayType::value_type ElementType;

    ElementType threshold = this->getComparisonOperatorThreshold();

    typedef
      LessThanFunctor< ArrayType,
                       ElementType,
                       ComparisonResultType >
      ArrayFunctor;

    typedef LessThanFunctor< ElementType,
                             ElementType,
                             typename ComparisonResultType::value_type >
      ElementFunctor;

    testOperatorX__ArrayND__Type(
      *this, ArrayFunctor(), ElementFunctor(), threshold);
  }


  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testOperatorLessThanOrEqualTo__Type()
  {
    // This member compares the elements of the array with a constant.
    typedef typename ArrayType::value_type ElementType;

    ElementType threshold = this->getComparisonOperatorThreshold();

    typedef
      LessThanOrEqualToFunctor<
        ArrayType,
        ElementType,
        ComparisonResultType >
      ArrayFunctor;

    typedef
      LessThanOrEqualToFunctor<
        ElementType,
        ElementType,
        typename ComparisonResultType::value_type >
      ElementFunctor;

    testOperatorX__ArrayND__Type(
      *this, ArrayFunctor(), ElementFunctor(), threshold);
  }


  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testOperatorPlus__ArrayND__ArrayND()
  {
    // This non-member function returns an array whose elements are
    // the sums of the corresponding elements of the two input arrays.

    typedef typename ArrayType::value_type ElementType;

    typedef
      PlusFunctor< ArrayType, ArrayType, ArrayType >
      ArrayFunctor;

    typedef
      PlusFunctor<
        ElementType,
        ElementType,
        ElementType >
      ElementFunctor;

    testOperatorX__ArrayND__ArrayND(
      *this, ArrayFunctor(), ElementFunctor());
  }


  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testOperatorMinus__ArrayND__ArrayND()
  {
    // This non-member function returns an array whose elements are
    // the differences of the corresponding elements of the two input
    // arrays.

    typedef typename ArrayType::value_type ElementType;

    typedef
      MinusFunctor< ArrayType, ArrayType, ArrayType >
      ArrayFunctor;

    typedef
      MinusFunctor<
        ElementType,
        ElementType,
        ElementType >
      ElementFunctor;

    testOperatorX__ArrayND__ArrayND(
      *this, ArrayFunctor(), ElementFunctor());
  }


  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testOperatorTimes__ArrayND__ArrayND()
  {
    // This non-member function returns an array whose elements are
    // the products of the corresponding elements of the two input
    // arrays.

    typedef typename ArrayType::value_type ElementType;

    typedef
      TimesFunctor< ArrayType, ArrayType, ArrayType >
      ArrayFunctor;

    typedef
      TimesFunctor<
        ElementType,
        ElementType,
        ElementType >
      ElementFunctor;

    testOperatorX__ArrayND__ArrayND(
      *this, ArrayFunctor(), ElementFunctor());
  }


  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testOperatorDividedBy__ArrayND__ArrayND()
  {
    // This non-member function returns an array whose elements are
    // the result of dividing corresponding element of the first input
    // argument by the corresponding element of the second input
    // argument.

    typedef typename ArrayType::value_type ElementType;

    typedef
      DividedByFunctor< ArrayType, ArrayType, ArrayType >
      ArrayFunctor;

    typedef
      DividedByFunctor<
        ElementType,
        ElementType,
        ElementType >
      ElementFunctor;

    testOperatorX__ArrayND__ArrayND(
      *this, ArrayFunctor(), ElementFunctor());
  }


  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testOperatorPlus__ArrayND__Type()
  {
    // This non-member function returns an array whose elements are
    // equal to the sums of the corresponding elements of the first
    // argument and value of the second argument.
    typedef typename ArrayType::value_type ElementType;
    ElementType increment = this->getIncrementOperatorArgument();

    typedef
      PlusFunctor<
        ArrayType,
        ElementType,
        ArrayType >
      ArrayFunctor;

    typedef
      PlusFunctor<
        ElementType,
        ElementType,
        ElementType >
      ElementFunctor;

    testOperatorX__ArrayND__Type(
      *this, ArrayFunctor(), ElementFunctor(), increment);
  }


  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testOperatorMinus__ArrayND__Type()
  {
    // This non-member function returns an array whose elements are
    // equal to the differences of the corresponding elements of the
    // first argument and value of the second argument.
    typedef typename ArrayType::value_type ElementType;
    ElementType increment = this->getIncrementOperatorArgument();

    typedef
      MinusFunctor<
        ArrayType,
        ElementType,
        ArrayType >
      ArrayFunctor;

    typedef
      MinusFunctor<
        ElementType,
        ElementType,
        ElementType >
      ElementFunctor;

    testOperatorX__ArrayND__Type(
      *this, ArrayFunctor(), ElementFunctor(), increment);
  }


  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testOperatorTimes__ArrayND__Type()
  {
    // This non-member function returns an array whose elements are
    // equal to the products of the corresponding elements of the first
    // argument and value of the second argument.
    typedef typename ArrayType::value_type ElementType;
    ElementType multiplier = this->getMultiplicationOperatorArgument();

    typedef
      TimesFunctor<
        ArrayType,
        ElementType,
        ArrayType >
      ArrayFunctor;

    typedef
      TimesFunctor<
        ElementType,
        ElementType,
        ElementType >
      ElementFunctor;

    testOperatorX__ArrayND__Type(
      *this, ArrayFunctor(), ElementFunctor(), multiplier);
  }


  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testOperatorDividedBy__ArrayND__Type()
  {
    // This non-member function returns an array whose elements are
    // equal to the dividends of the corresponding elements of the
    // first argument and value of the second argument.
    typedef typename ArrayType::value_type ElementType;
    ElementType multiplier = this->getMultiplicationOperatorArgument();

    typedef
      DividedByFunctor<
        ArrayType,
        ElementType,
        ArrayType >
      ArrayFunctor;

    typedef
      DividedByFunctor<
        ElementType,
        ElementType,
        ElementType >
      ElementFunctor;

    testOperatorX__ArrayND__Type(
      *this, ArrayFunctor(), ElementFunctor(), multiplier);
  }


  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testOperatorPlus__Type__ArrayND()
  {
    // This non-member function returns an array whose elements are
    // equal to the sums of the second argument with the corresponding
    // elements of the first argument.
    typedef typename ArrayType::value_type ElementType;
    ElementType increment = this->getIncrementOperatorArgument();

    typedef
      PlusFunctor<
        ElementType,
        ArrayType,
        ArrayType >
      ArrayFunctor;

    typedef
      PlusFunctor<
        ElementType,
        ElementType,
        ElementType >
      ElementFunctor;

    testOperatorX__Type__ArrayND(
      *this, ArrayFunctor(), ElementFunctor(), increment);
  }


//   template <class FixtureType, class ArrayType, class ComparisonResultType>
//   void
//   ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
//   testOperatorMinus__Type__ArrayND()
//   {
//     // This non-member function returns an array whose elements are
//     // equal to the differences of the second argument with the
//     // corresponding elements of the first argument.
//     typedef typename ArrayType::value_type ElementType;
//     ElementType increment = this->getIncrementOperatorArgument();

//     typedef
//       MinusFunctor<
//         ElementType,
//         ArrayType,
//         ArrayType >
//       ArrayFunctor;

//     typedef
//       MinusFunctor<
//         ElementType,
//         ElementType,
//         ElementType >
//       ElementFunctor;

//     testOperatorX__Type__ArrayND(
//       *this, ArrayFunctor(), ElementFunctor(), increment);
//   }


  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testOperatorTimes__Type__ArrayND()
  {
    // This non-member function returns an array whose elements are
    // equal to the products of the second argument with the
    // corresponding elements of the first argument.
    typedef typename ArrayType::value_type ElementType;
    ElementType multiplier = this->getMultiplicationOperatorArgument();

    typedef
      TimesFunctor<
        ElementType,
        ArrayType,
        ArrayType >
      ArrayFunctor;

    typedef
      TimesFunctor<
        ElementType,
        ElementType,
        ElementType >
      ElementFunctor;

    testOperatorX__Type__ArrayND(
      *this, ArrayFunctor(), ElementFunctor(), multiplier);
  }


//   template <class FixtureType, class ArrayType, class ComparisonResultType>
//   void
//   ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
//   testOperatorDividedBy__Type__ArrayND()
//   {
//     // This non-member function returns an array whose elements are
//     // equal to the dividends of the second argument with the
//     // corresponding elements of the first argument.
//     typedef typename ArrayType::value_type ElementType;
//     ElementType multiplier = this->getMultiplicationOperatorArgument();

//     typedef
//       DividedByFunctor<
//         ElementType,
//         ArrayType,
//         ArrayType >
//       ArrayFunctor;

//     typedef
//       DividedByFunctor<
//         ElementType,
//         ElementType,
//         ElementType >
//       ElementFunctor;

//     testOperatorX__Type__ArrayND(
//       *this, ArrayFunctor(), ElementFunctor(), multiplier);
//   }


  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testOutputOperator()
  {
    ArrayType array0 = this->getFibonacciArray();
    std::ostringstream outputStream;
    outputStream << std::fixed << array0;
    std::istringstream inputStream(outputStream.str());
    ArrayType array1;
    inputStream >> array1;
    this->checkValueEquality(
      array0, array1, static_cast<typename ArrayType::value_type>(1.0E-10));
  }


  template <class FixtureType, class ArrayType, class ComparisonResultType>
  void
  ArrayTestCommon<FixtureType, ArrayType, ComparisonResultType>::
  testInputOperator()
  {
    // No explicit test.
  }


} // namespace brick
