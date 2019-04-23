/**
***************************************************************************
* @file brick/numeric/polynomial.hh
*
* Header file declaring a class for representing simple polynomials.
*
* Copyright (C) 2006,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_NUMERIC_POLYNOMIAL_HH
#define BRICK_NUMERIC_POLYNOMIAL_HH

#include <brick/numeric/array1D.hh>

namespace brick {

  namespace numeric {


    /**
     ** This class represents polynomials of the form
     **
     **   p(x) = k0 + (k1 * x) + (k2 * x^2) + (k3 * x^3) ...
     **
     ** It provides operators for evaluating the polynomial, as well
     ** as for combining the polynomial with other polynomials.
     **/
    template <class Type>
    class Polynomial {
    public:

      /**
       * The default constructor makes a constant polynomial: p(x) = 1.
       */
      Polynomial();


      /**
       * This constructor makes a constant polynomial:
       *
       *   p(x) = coefficient0.
       *
       * @param coefficient0 This argument specifies the (constant)
       * value of the polynomial.
       */
      explicit
      Polynomial(Type coefficient0);


      /**
       * This constructor makes a first order polynomial:
       *
       *   p(x) = coefficient0 + (coefficient1 * x)
       *
       * @param coefficient1 This argument specifies one of the
       * coefficients of the polynomial.
       *
       * @param coefficient0 This argument specifies one of the
       * coefficients of the polynomial.
       */
      Polynomial(Type coefficient1, Type coefficient0);


      /**
       * This constructor makes a second order polynomial:
       *
       *   p(x) = coefficient0 + (coefficient1 * x) + (coefficient2 * x^2)
       *
       * @param coefficient2 This argument specifies one of the
       * coefficients of the polynomial.
       *
       * @param coefficient1 This argument specifies one of the
       * coefficients of the polynomial.
       *
       * @param coefficient0 This argument specifies one of the
       * coefficients of the polynomial.
       */
      Polynomial(Type coefficient2, Type coefficient1, Type coefficient0);


      /**
       * This constructor makes a polynomial of arbitrary order.
       *
       *   p(x) = coefficients[0] + (coefficients[1] * x)
       *          + (coefficients[2] * x2) + ...
       *          + (coefficients[N - 1] * x^(N-1)
       *
       * where N = coefficients.size()
       *
       * @param coefficients This argument is an array of coefficients
       * for the polynomial, as described above.
       */
      explicit
      Polynomial(Array1D<Type> const& coefficients);


// Apparently, this construction doesn't work with g++, 2012-02
#define BRICK_NUMERIC_POLYNOMIAL_TEMPLATED_CONSTRUCTORS 0
#if BRICK_NUMERIC_POLYNOMIAL_TEMPLATED_CONSTRUCTORS
      /**
       * This constructor makes a polynomial of arbitrary order.
       *
       *   p(x) = *beginIter + (*(beginIter + 1) * x)
       *          + ((*beginIter + 2) * x2) + ...
       *
       * where p(x) has as many terms as there are elements in the
       * input sequence.
       *
       * @param beginIter This argument defines the beginning of the
       * input sequence.
       *
       * @param endIter This argument defines the endning of the
       * input sequence.
       *
       * @param isSequence This argument is ignored.  Its role is to
       * disambiguate this constructor from the Polynomial<Type>(Type,
       * Type) constructor.
       */
      template <class Iter>
      Polynomial(Iter beginIter, Iter endIter, bool isSequence);
#endif /* #if BRICK_NUMERIC_POLYNOMIAL_TEMPLATED_CONSTRUCTORS */


      /**
       * The copy constructor does a deep copy.
       *
       * @param other This argument is the Polynomial instance to copy.
       */
      Polynomial(Polynomial<Type> const& other);


      /**
       * This member function returns the coefficients of the
       * polynomial.  The array values are arranged in the same way as
       * described in the documentation for constructor
       * Polynomial(const Array1D<Type>&).
       *
       * @return The return value is an array of coefficients.
       */
      Array1D<Type>
      getCoefficientArray() const;


      /**
       * This member function returns the order of the polynomial: 0
       * for constant; 1 for linear; 2 for quadratic; and so on.
       *
       * @return The return value is the order of the polynomial.
       */
      size_t
      getOrder() const {return m_coefficientArray.size() - 1;}


      /**
       * The assigment operator does a deep copy.
       *
       * @param other This argument is the Polynomial instance to copy.
       *
       * @return The return value is a reference to *this.
       */
      Polynomial<Type>&
      operator=(Polynomial<Type> const& other);


      /**
       * This operator evaluates the polynomial.
       *
       * @param xValue This argument specifies the value of x at which
       * to evaluate the polynomial.
       *
       * @return The return value is the result of the evaluation.
       */
      Type
      operator()(Type xValue) const;


      /**
       * This operator multiplies the polynomial by another
       * polynomial.  After the multipliation, the result of calling
       * this->operator()(xValue) will be equal to the previous
       * (before the call to this->operator*=()) result of
       * this->operator()(xValue) multiplied by the result of
       * other.operator()(xValue).
       *
       * @param other This argument is the the polynomial by which to
       * multiply *this.
       *
       * @return The return value is a reference to *this.
       */
      Polynomial<Type>&
      operator*=(Polynomial<Type> const& other);


      /**
       * This operator adds another polynomial to *this.  After the
       * addition, the result of calling this->operator()(xValue)
       * will be equal to the previous (before the call to
       * this->operator+=()) result of this->operator()(xValue) plus
       * the result of other.operator()(xValue).
       *
       * @param other This argument is the the polynomial to be added
       * to *this.
       *
       * @return The return value is a reference to *this.
       */
      Polynomial<Type>&
      operator+=(Polynomial<Type> const& other);


      /**
       * This operator subtracts another polynomial from *this.  After
       * the subtraction, the result of calling
       * this->operator()(xValue) will be equal to the previous
       * (before the call to this->operator+=()) result of
       * this->operator()(xValue) minus the result of
       * other.operator()(xValue).
       *
       * @param other This argument is the the polynomial to be subtracted
       * from *this.
       *
       * @return The return value is a reference to *this.
       */
      Polynomial<Type>&
      operator-=(Polynomial<Type> const& other);

    private:

      Array1D<Type> m_coefficientArray;

    };


    /* ============ Non-member function declarations ============ */


    /**
     * This operator multiplies two Polynomial instances.  Values
     * returned by operator()(xValue) of the result of this
     * multiplication will be equal to the product of the values
     * returned by the operator()(xValue) of arg0 and arg1.
     *
     * @param arg0 This argument is one of the Polynomial instances to
     * be multiplied.
     *
     * @param arg1 This argument is one of the Polynomial instances to
     * be multiplied.
     *
     * @return The return value is the result of the multiplication.
     */
    template<class Type>
    Polynomial<Type>
    operator*(Polynomial<Type> const& arg0, Polynomial<Type> const& arg1);


    /**
     * This operator adds two Polynomial instances.  Values returned by
     * operator()(xValue) of the result of this multiplication will be
     * equal to the sum of the values returned by the operator()(xValue)
     * of arg0 and arg1.
     *
     * @param arg0 This argument is one of the Polynomial instances to
     * be added.
     *
     * @param arg1 This argument is one of the Polynomial instances to
     * be added.
     *
     * @return The return value is the result of the addition.
     */
    template<class Type>
    Polynomial<Type>
    operator+(Polynomial<Type> const& arg0, Polynomial<Type> const& arg1);


    /**
     * This operator subtracts two Polynomial instances.  Values
     * returned by operator()(xValue) of the result of this
     * multiplication will be equal to arg0.operator(xValue) -
     * arg1.operator()(xValue).
     *
     * @param arg0 This argument is the Polynomial instance from which
     * arg1 is to be subtracted.
     *
     * @param arg1 This argument is the Polynomial instance to be
     * subtracted from arg0.
     *
     * @return The return value is the result of the subtraction.
     */
    template<class Type>
    Polynomial<Type>
    operator-(Polynomial<Type> const& arg0, Polynomial<Type> const& arg1);


  } // namespace numeric

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/polynomial_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_POLYNOMIAL_HH */
