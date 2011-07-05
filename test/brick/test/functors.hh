/**
***************************************************************************
* @file brick/test/functors.hh
* 
* Header file declaring macros for brick::test library.
*
* Copyright (C) 2004-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_TEST_FUNCTORS_HH
#define BRICK_TEST_FUNCTORS_HH

namespace brick {

  namespace test {

    /**
     ** Functor template for comparing two values to determine
     ** equivalence within a specified precision.  This is just like
     ** std::equal<>, except that it has the additional capability of
     ** allowing for numerical precision in its comparisons.
     **/
    template <class Type>
    class ApproximatelyEqualFunctor
      : public std::binary_function<Type, Type, bool>
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

  } // namespace test

} // namespace brick

#endif /* #ifndef BRICK_TEST_FUNCTORS_HH */
