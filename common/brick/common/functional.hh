/**
***************************************************************************
* @file brick/common/functional.hh
*
* Header file declaring functors that are not part of the C++ standard
* library.
*
* Copyright (C) 2003-2011, David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_COMMON_FUNCTIONAL_HH
#define BRICK_COMMON_FUNCTIONAL_HH

#include <algorithm>
#include <functional>
#include <utility>

namespace brick {

  namespace common {

    /**
     ** Functor template for composing one binary function,
     ** functor0(x, y), and two unary functions functor1(x) and functor2(x)
     ** so that the result is functor0(functor1(x), functor2(y)).
     **/
    template <class Functor0, class Functor1, class Functor2>
    class BinaryComposeFunctor
    {
    public:
      /**
       * The constructor accepts instances of the three functor types to
       * be composed, and makes local copies of them for use in
       * operator()(...).
       *
       * @param functor0 During evaluation of operator()(...), the
       * return value of a copy of functor1 will be passed as the first
       * argument to this functor (or rather, a copy of this functor),
       * the return value of a copy of functor2 will be passed as the
       * second argument, and the return value of the copy of this
       * functor will be passed to the calling context.
       *
       * @param functor1 During evaluation of operator()(...), the
       * return value of this functor (or rather, a copy of this
       * functor) will be passed as the first argument to a copy of
       * functor0, and the return value of the functor0 copy will be
       * passed to the calling context.
       *
       * @param functor2 During evaluation of operator()(...), the
       * return value of this functor (or rather, a copy of this
       * functor) will be passed as the second argument to a copy of
       * functor0, and the return value of the functor0 copy will be
       * passed to the calling context.
       */
      BinaryComposeFunctor(const Functor0& functor0, const Functor1& functor1,
                           const Functor2& functor2)
        : m_functor0(functor0), m_functor1(functor1), m_functor2(functor2) {}

      /**
       * This operator passes its first argument to the operator()(...)
       * method of a copy of constructor argument functor1, passes its
       * second argument to the operator()(...) method of a copy of
       * constructor argument functor2, and then passes the result of
       * these calls as the first and second arguments of the
       * operator()(...) method of functor0, and returns the result.
       *
       * @param argument0 This argument will be passed as input to the
       * functor1 copy.
       *
       * @param argument1 This argument will be passed as input to the
       * functor2 copy.
       *
       * @return The result of processing the argument with the three
       * composed functors:
       *   returnValue = functor0(functor1(argument0), functor2(argument1)).
       */
      typename Functor0::result_type
      operator()(const typename Functor1::argument_type& argument0,
                 const typename Functor2::argument_type& argument1) {
        return m_functor0(m_functor1(argument0), m_functor2(argument1));
      }

    protected:

      /// This protected member stores a copy of functor0.
      Functor0 m_functor0;

      /// This protected member stores a copy of functor1.
      Functor1 m_functor1;

      /// This protected member stores a copy of functor2.
      Functor1 m_functor2;
    };


    /**
     ** Functor template for composing one unary function,
     ** functor0(x) and one binary function functor1(x, y)
     ** so that the result is functor0(functor1(x, y)).
     **/
    template <class Functor0, class Functor1>
    class ComposeFunctor_1_2
    {
    public:
      /**
       * The constructor accepts instances of the two functor types to
       * be composed, and makes local copies of them for use in
       * operator()(...).
       *
       * @param functor0 During evaluation of operator()(...), the
       * return value of a copy of functor1 will be passed as the
       * argument to this functor (or rather, a copy of this functor),
       * and the return value of the copy of this functor will be passed
       * to the calling context.
       *
       * @param functor1 During evaluation of operator()(...), the
       * return value of this functor (or rather, a copy of this
       * functor) will be passed as the first argument to a copy of
       * functor0, and the return value of the functor0 copy will be
       * passed to the calling context.
       */
      ComposeFunctor_1_2(const Functor0& functor0, const Functor1& functor1)
        : m_functor0(functor0), m_functor1(functor1) {}


      /**
       * This operator passes its arguments to the operator()(...)
       * method of a copy of constructor argument functor1, then passes
       * the result of this call as the argument of the operator()(...)
       * method of functor0, and returns the result.
       *
       * @param argument0 This argument will be passed as the first
       * argument of the functor1 copy.
       *
       * @param argument1 This argument will be passed as the second
       * argument of the functor1 copy.
       *
       * @return The result of processing the argument with the composed
       * functors:
       *   returnValue = functor0(functor1(argument0, argument1)).
       */
      typename Functor0::result_type
      operator()(const typename Functor1::first_argument_type& argument0,
                 const typename Functor1::second_argument_type& argument1) {
        return m_functor0(m_functor1(argument0, argument1));
      }

    protected:

      /// This protected member stores a copy of functor0.
      Functor0 m_functor0;

      /// This protected member stores a copy of functor1.
      Functor1 m_functor1;
    };


    /**
     ** Functor template for extracting the first element of a
     ** std::pair.
     **/
    template <class Type0, class Type1>
    class ExtractFirstFunctor
    {
    public:
      /**
       * This is a very lightweight class.  The constructor does nothing.
       */
      ExtractFirstFunctor() {}


      /**
       * This version of the application operator returns the value of
       * the first component of its argument.
       *
       * @param argument This argument is a std::pair instance.
       *
       * @return The return value is a copy of argument.first.
       */
      Type0
      operator()(const std::pair<Type0, Type1>& argument) {
        return argument.first;
      }


      /**
       * This version of the application operator returns a reference to
       * the first component of its argument.
       *
       * @param argument This argument is a std::pair instance.
       *
       * @return The return value is a reference to argument.first.
       */
      Type0&
      operator()(std::pair<Type0, Type1>& argument) {
        return argument.first;
      }
    };


    /**
     ** Functor template for extracting the second element of a
     ** std::pair.
     **/
    template <class Type0, class Type1>
    class ExtractSecondFunctor
    {
    public:
      /**
       * This is a very lightweight class.  The constructor does nothing.
       */
      ExtractSecondFunctor() {}


      /**
       * This version of the application operator returns the value of
       * the second component of its argument.
       *
       * @param argument This argument is a std::pair instance.
       *
       * @return The return value is a copy of argument.second.
       */
      Type1
      operator()(const std::pair<Type0, Type1>& argument) {
        return argument.second;
      }


      /**
       * This version of the application operator returns a reference to
       * the second component of its argument.
       *
       * @param argument This argument is a std::pair instance.
       *
       * @return The return value is a reference to argument.second.
       */
      Type1&
      operator()(std::pair<Type0, Type1>& argument) {
        return argument.second;
      }
    };


    /**
     ** Functor template for comparing two values to determine
     ** equivalence within a specified precision.  This is just like
     ** std::equal<>, except that it has the additional capability of
     ** allowing for numerical precision in its comparisons.
     **/
    template <class Type>
    class ApproximatelyEqualFunctor
    {
    public:
      /**
       * The constructor sets the threshold for what is considered
       * approximately equal.  For example, if the constructor argument
       * is 1.0E-6, then two values will be considered equal if the
       * absolute value of their difference is less than or equal to
       * 1.0E-6.
       *
       * @param epsilon This argument sets the largest difference that
       * will be considered equivalent.
       */
      ApproximatelyEqualFunctor(const Type& epsilon=static_cast<Type>(0))
        : m_epsilon(epsilon) {}

      /**
       * The application operator returns true if the difference between
       * its two arguments is less than epsilon and greater than
       * -(epsilon), where epsilon is specified in the constructor.
       *
       * @param argument0 This argument will be compared with the second
       * argument.
       *
       * @param argument1 This argument will be compared with the first
       * argument.
       *
       * @return The return value is true if the values are
       * approximately equal, false otherwise.
       */
      inline bool
      operator()(const Type& argument0, const Type& argument1) {
        Type difference = argument0 - argument1;
        return ((difference <= m_epsilon) && (difference >= (-m_epsilon)));
      }

    protected:

      /// This protected member function stores the comparison tolerance.
      Type m_epsilon;
    };


    /**
     ** Functor template much like std::pointer_to_binary_function, but
     ** specifically for functions that take const reference arguments.
     ** Currently, calling std::mem_fun() for such a function pointer
     ** causes "forming reference to reference" compile errors.
     **/
    template <class ArgumentType0, class ArgumentType1, class ResultType>
    class PointerToBinaryFunctionRA
    {
    public:
      /// Typedef describing what type of function is to be wrapped by this
      /// class.
      typedef ResultType (*FunctionPtrType)(const ArgumentType0&,
                                            const ArgumentType1&);

      /**
       * Constructor requires a function pointer to wrap.
       */
      explicit
      PointerToBinaryFunctionRA(FunctionPtrType functionPtr)
        : m_functionPtr(functionPtr) {}

      /**
       * Call the wrapped function pointer and return the result.
       *
       * @param argument0 Will be passed as an argument to the function pointer.
       *
       * @param argument1 Will be passed as an argument to the function pointer.
       *
       * @return The result of the function call.
       */
      inline ResultType
      operator()(const ArgumentType0& argument0,
                 const ArgumentType0& argument1) const {
        return m_functionPtr(argument0, argument1);
      }

    private:
      FunctionPtrType m_functionPtr;
    };


    /**
     ** Functor template that uses static_cast to convert instances
     ** of one type into instances of another.
     **/
    template <class TypeIn, class TypeOut>
    struct StaticCastFunctor
    {

      /**
       * Static cast the input argument to TypeOut, and return the result.
       *
       * @param input Will be cast to type TypeOut.
       *
       * @return The result of the cast.
       */
      inline TypeOut
      operator()(const TypeIn& input) {
        return static_cast<TypeOut>(input);
      }
    };


    /**
     ** Functor template for composing two unary functions functor0(x)
     ** and functor1(x) so that the result is functor0(functor1(x)).
     **/
    template <class Functor0, class Functor1>
    class UnaryComposeFunctor
    {
    public:
      /**
       * Constructor accepts instances of the two functor types to be
       * composed, and makes local copies of them for use in
       * operator()(...).
       *
       * @param functor0 During evaluation of operator()(...), the
       * return value of a copy of functor1 will be passed to this
       * functor (or rather, a copy of this functor) and the return
       * value of the copy of this functor will be passed to the calling
       * context.
       *
       * @param functor1 During evaluation of operator()(...), the
       * return value of this functor (or rather, a copy of this
       * functor) will be passed to a copy of functor0, and the return
       * value of the functor0 copy will be passed to the calling
       * context.
       */
      UnaryComposeFunctor(const Functor0& functor0, const Functor1& functor1)
        : m_functor0(functor0), m_functor1(functor1) {}

      /**
       * This operator passes its argument to the operator()(...) method
       * of a copy of constructor argument functor1, then passes the
       * result to the operator()(...) method of a copy of constructor
       * argument functor0, and returns the result.
       *
       * @param argument This argument will be passed as input to the
       * functor1 copy.  @return The result of processing the argument
       * with the two composed functors:
       * returnValue = functor0(functor1(argument)).
       */
      typename Functor0::result_type
      operator()(const typename Functor1::argument_type& argument) {
        return m_functor0(m_functor1(argument));
      }

    protected:

      /// This protected member stores a copy of functor0.
      Functor0 m_functor0;

      /// This protected member stores a copy of functor1.
      Functor1 m_functor1;
    };


    // =================== Helper functions =================== //

    /**
     * This convenience function constructs an ApproximatelyEqualFunctor<Type>
     * functor and applies it two the first two arguments.  The third
     * argument sets the threshold for approximate equality.  For
     * example,
     * approximatelyEqual(1.0003, 1.0005, 1.0E-6) will return false, while
     * approximatelyEqual(1.0003, 1.0005, 1.0E-3) will return true.
     *
     * @param argument0 This argument is the first value to be compared.
     *
     * @param argument1 This argument is the second value to be compared.
     *
     * @param epsilon This argument is passed directly to the
     * constructor of ApproximatelyEqual, and specifies the threshold
     * for approximate equality.
     *
     * @return The return value is true if the values are approximately
     * equal, false otherwise.
     */
    template<class Type>
    inline bool
    approximatelyEqual(const Type& argument0, const Type& argument1,
                       const Type& epsilon) {
      return (ApproximatelyEqualFunctor<Type>(epsilon))(argument0, argument1);
    }


    /**
     * This is a convenience function that makes it easy to create
     * BinaryComposeFunctor instances.
     *
     * @param functor0 This will be passed as the first constructor
     * argument of the returned BinaryComposeFunctor instance.
     *
     * @param functor1 This will be passed as the second constructor
     * argument of the returned BinaryComposeFunctor instance.
     *
     * @param functor2 This will be passed as the third constructor
     * argument of the returned BinaryComposeFunctor instance.
     *
     * @return A BinaryComposeFunctor template instantiation of the
     * appropriate type.
     */
    template <class Functor0, class Functor1, class Functor2>
    inline BinaryComposeFunctor<Functor0, Functor1, Functor2>
    binaryComposeFunctor(const Functor0& functor0,
                         const Functor1& functor1,
                         const Functor2& functor2) {
      return BinaryComposeFunctor<Functor0, Functor1, Functor2>(
        functor0, functor1, functor2);
    }


    /**
     * This convenience function clips its first argument to lie
     * within the closed range defined by it's remaining arguments.
     *
     * @param value This argument is the value to be clipped.
     *
     * @param lowerBound The value will be clipped so that it is
     * greater than or equal to this argument.
     *
     * @param upperBound The value will be clipped so that it is less
     * than or equal to this argument.
     *
     * @return The clipped value is returned.
     */
    template <class Type>
    Type
    clip(Type value, Type lowerBound, Type upperBound) {
#ifdef WIN32
      return std::_cpp_max(std::_cpp_min(value, upperBound), lowerBound);
#else
      return std::max(std::min(value, upperBound), lowerBound);
#endif
    }


    /**
     * This is a convenience function that makes it easy to create
     * ComposeFunctor_1_2 instances.
     *
     * @param functor0 This will be passed as the first constructor
     * argument of the returned functor.
     *
     * @param functor1 This will be passed as the second constructor
     * argument of the returned functor.
     *
     * @return A ComposeFunctor_1_2 template instantiation of the
     * appropriate type.
     */
    template <class Functor0, class Functor1>
    inline ComposeFunctor_1_2<Functor0, Functor1>
    composeFunctor_1_2(const Functor0& functor0,
                       const Functor1& functor1) {
      return ComposeFunctor_1_2<Functor0, Functor1>(functor0, functor1);
    }


    /**
     * This is a convenience function that makes it easy to create
     * UnaryComposeFunctor instances.
     *
     * @param functor0 This will be passed as the first constructor
     * argument of the returned UnaryComposeFunctor instance.
     *
     * @param functor1 This will be passed as the second constructor
     * argument of the returned UnaryComposeFunctor instance.
     *
     * @return A UnaryComposeFunctor template instantiation of the
     * appropriate type.
     */
    template <class Functor0, class Functor1>
    inline UnaryComposeFunctor<Functor0, Functor1>
    unaryComposeFunctor(const Functor0& functor0,
                        const Functor1& functor1) {
      return UnaryComposeFunctor<Functor0, Functor1>(functor0, functor1);
    }


  } // namespace common

}  // namespace brick


/* ======= Declarations to maintain compatibility with legacy code. ======= */

namespace brick {

  using common::BinaryComposeFunctor;
  using common::ComposeFunctor_1_2;
  using common::ExtractFirstFunctor;
  using common::ExtractSecondFunctor;
  using common::ApproximatelyEqualFunctor;
  using common::PointerToBinaryFunctionRA;
  using common::StaticCastFunctor;
  using common::UnaryComposeFunctor;
  using common::approximatelyEqual;
  using common::binaryComposeFunctor;
  using common::composeFunctor_1_2;
  using common::unaryComposeFunctor;

} // namespace brick

/* ================= Specializations =================== */

#include <brick/common/types.hh>

namespace brick {

  namespace common {

    template <>
    inline bool
	ApproximatelyEqualFunctor<bool>::
    operator()(const bool& argument0, const bool& argument1) {
        return argument0 == argument1;
	}


	template<>
	inline bool
	ApproximatelyEqualFunctor<size_t>::
	operator()(const size_t& argument0, const size_t& argument1) {
        Int64 difference = static_cast<Int64>(argument0) - static_cast<Int64>(argument1);
        return ((difference <= static_cast<Int64>(m_epsilon))
			&& (difference >= (-static_cast<Int64>(m_epsilon))));
      }

  } // namespace common

}  // namespace brick

#endif /* #ifndef BRICK_COMMON_FUNCTIONAL_HH */
