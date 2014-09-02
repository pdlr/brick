/**
***************************************************************************
* @file brick/numeric/differentiableScalar.hh
*
* Header file declaring a class for representing scalar quantities and
* their first derivatives.
*
* Copyright (C) 2014 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_NUMERIC_DIFFERENTIABLESCALAR_HH
#define BRICK_NUMERIC_DIFFERENTIABLESCALAR_HH

#include <inttypes.h>

namespace brick {

  namespace numeric {


    /**
     ** This class represents a scalar value and its first derivative
     ** with respect to one or more parameters.  Use it as follows:
     **
     ** @code
     **   // Compute the derivative with respect to x of the function
     **   // f(x) = x^2 + sin(x), evaluated at x == 5.0;
     **   DifferentiableScalar xx(5.0);
     **   DifferentiableScalar result = xx * xx + sin(xx);
     **   std::cout << "df/dx = " << result.getDerivative() << std::endl;
     ** @endcode
     **
     ** The first template argument of this class determines the
     ** precision with which calculations are to be carried out, while
     ** the second template argument specifies the dimensionality of
     ** the parameter vector with respect to which derivatives should
     ** be tracked.  For example:
     **
     ** @code
     **   // Compute the derivative with respect to x and y of the function
     **   // f(x) = x^2 + 2xy + sin(y), evaluated at x == 5.0 and y = 3.0.
     **   Array1D<double> xDerivatives("[1.0, 0.0]");
     **   Array1D<double> yDerivatives("[0.0, 1.0]");
     **   DifferentiableScalar<double, 2> xx(5.0, xDerivatives.begin());
     **   DifferentiableScalar<double, 2> yy(3.0, yDerivatives.begin());
     **   DifferentiableScalar result = xx * xx + 2.8 * xx * yy + sin(yy);
     **   std::cout << "df/dx = " << result.getPartialDerivative(0)
     **             << std::endl;
     **   std::cout << "df/dy = " << result.getPartialDerivative(1)
     **             << std::endl;
     ** @endcode
     **
     **/
    template <class Type, uint32_t Dimension = 1>
    class DifferentiableScalar {
    public:

      /** 
       * The default constructor makes a differentiableScalar with
       * value 0, first partial derivative equal to 1, and all other
       * partials equal to 0.
       */
      DifferentiableScalar();


      /** 
       * This default constructor makes a differentiableScalar with a
       * user specified value, first partial derivative equal to 1,
       * and all other partials equal to 0.
       *
       * @param value This argument specifies the scalar value of
       * *this.
       */
      explicit
      DifferentiableScalar(Type value);


      /** 
       * This default constructor makes a differentiableScalar with a
       * user specified value and partial derivatives explicitly set
       * by the second constructor argument.
       * 
       * @param value This argument specifies the scalar value of
       * *this.  @param partialsBegin This iterator points to a
       * sequence of Type specifying a vector of partial derivatives.
       * The number of values in this sequence must be at least equal
       * to the value of class template argument Dimension.
       */
      template <class Iter>
      explicit
      DifferentiableScalar(Type value, Iter partialsBegin);

      
      /** 
       * The copy constructor does a deep copy.
       * 
       * @param other This argument is the DifferentiableScalar
       * instance to copy.
       */
      DifferentiableScalar(DifferentiableScalar<Type, Dimension> const& other);


      /** 
       * The assigment operator does a deep copy.
       * 
       * @param other This argument is the DifferentiableScalar
       * instance to copy.
       *
       * @return The return value is a reference to *this.
       */
      DifferentiableScalar<Type, Dimension>&
      operator=(DifferentiableScalar<Type, Dimension> const& other);
      

      /** 
       * This operator multiplies the differentiableScalar by another
       * differentiableScalar.  After the multipliation, the partial
       * derivatives of *this will be updated to reflect the partial
       * derivatives of the product.
       * 
       * @param other This argument is the differentiableScalar
       * instance by which to multiply *this.
       * 
       * @return The return value is a reference to *this.
       */
      DifferentiableScalar<Type, Dimension>&
      operator*=(DifferentiableScalar<Type, Dimension> const& other);


      /** 
       * This operator divides the differentiableScalar by another
       * differentiableScalar.  After the division, the partial
       * derivatives of *this will be updated to reflect the partial
       * derivatives of the result.
       * 
       * @param other This argument is the differentiableScalar
       * instance by which to divide *this.
       * 
       * @return The return value is a reference to *this.
       */
      DifferentiableScalar<Type, Dimension>&
      operator/=(DifferentiableScalar<Type, Dimension> const& other);


      /** 
       * This operator adds another differentiableScalar to *this.
       * After the addition, the partial derivatives of *this will be
       * updated to reflect the partial derivatives of the sum.
       * 
       * @param other This argument is the the differentiableScalar to
       * be added to *this.
       * 
       * @return The return value is a reference to *this.
       */
      DifferentiableScalar<Type, Dimension>&
      operator+=(DifferentiableScalar<Type, Dimension> const& other);


      /** 
       * This operator subtracts another differentiableScalar from *this.
       * After the addition, the partial derivatives of *this will be
       * updated to reflect the partial derivatives of the difference.
       * 
       * @param other This argument is the the differentiableScalar to
       * be subtracted from *this.
       * 
       * @return The return value is a reference to *this.
       */
      DifferentiableScalar<Type, Dimension>&
      operator-=(DifferentiableScalar<Type, Dimension> const& other);


      /** 
       * This member function is equivalent to
       * this->getPartialDerivative(0).
       * 
       * @return The return value is the first partial derivative of
       * *this.
       */
      Type const&
      getDerivative() const {return m_partials[0];}


      /** 
       * This member function returns the first derivative of *this
       * with respect to the "ii"-th parameter.
       * 
       * @return The return value is the requested partial derivative.
       */
      Type const&
      getPartialDerivative(uint32_t ii) const {return m_partials[ii];}


      /** 
       * This member returns the current scalar value.
       * 
       * @return The return value is the value of *this.
       */
      Type const&
      getValue() const {return m_value;}


      /** 
       * This member function is equivalent to
       * this->setPartialDerivative(0, value).
       * 
       * @return The return value is the first partial derivative of
       * *this after the call to setDerivative.
       */
      Type const&
      setDerivative(Type const& derivative);

      
      /** 
       * This member function sets the first derivative of *this
       * with respect to the "ii"-th parameter.
       * 
       * @return The return value is the updated partial derivative.
       */
      Type const&
      setPartialDerivative(uint32_t ii, Type const& derivative);

      
      /** 
       * This member returns the current scalar value.  It does not
       * affect derivatives.
       * 
       * @return The return value is the value of *this after setting.
       */
      Type const&
      setValue(Type const& value) {return m_value = value;}
      
    private:

      Type m_value;
      Type m_partials[Dimension];
      
    };

    
    /* ============ Non-member function declarations ============ */

  
    /** 
     * This operator multiplies two DifferentiableScalar instances,
     * applying the chain rule to compute the derivatives of the
     * result.
     * 
     * @param arg0 This argument is one of the DifferentiableScalar
     * instances to be multiplied.
     * 
     * @param arg1 This argument is one of the DifferentiableScalar
     * instances to be multiplied.
     * 
     * @return The return value is the result of the multiplication.
     */
    template<class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>
    operator*(DifferentiableScalar<Type, Dimension> const& arg0,
              DifferentiableScalar<Type, Dimension> const& arg1);


    /** 
     * This operator divides two DifferentiableScalar instances,
     * applying the quotient rule to compute the derivatives of the
     * result.
     * 
     * @param arg0 This argument is the DifferentiableScalar instance
     * to be divided.
     * 
     * @param arg1 This argument is the DifferentiableScalar instance
     * by which to divide.
     * 
     * @return The return value is the result of the division.
     */
    template<class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>
    operator/(DifferentiableScalar<Type, Dimension> const& arg0,
              DifferentiableScalar<Type, Dimension> const& arg1);


    /** 
     * This operator adds two DifferentiableScalar instances.
     * 
     * @param arg0 This argument is one of the DifferentiableScalar
     * instances to be added.
     * 
     * @param arg1 This argument is one of the DifferentiableScalar
     * instances to be added.
     * 
     * @return The return value is the result of the addition.
     */
    template<class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>
    operator+(DifferentiableScalar<Type, Dimension> const& arg0,
              DifferentiableScalar<Type, Dimension> const& arg1);
  

    /** 
     * This operator subtracts two DifferentiableScalar instances.  Values
     * returned by operator()(xValue) of the result of this
     * multiplication will be equal to arg0.operator(xValue) -
     * arg1.operator()(xValue).
     * 
     * @param arg0 This argument is the DifferentiableScalar instance
     * from which arg1 is to be subtracted.
     * 
     * @param arg1 This argument is the DifferentiableScalar instance
     * to be subtracted from arg0.
     * 
     * @return The return value is the result of the subtraction.
     */
    template<class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>
    operator-(DifferentiableScalar<Type, Dimension> const& arg0,
              DifferentiableScalar<Type, Dimension> const& arg1);


    /** 
     * This operator computes the sine of a DifferentiableScalar
     * instance, with partial derivatives.
     * 
     * @param arg0 The sine of this argument will be computed.
     * 
     * @return The return value is the result of the calculation.
     */
    template<class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>
    sine(DifferentiableScalar<Type, Dimension> const& arg0);
    

    /** 
     * This operator computes the cosine of a DifferentiableScalar
     * instance, with partial derivatives.
     * 
     * @param arg0 The cosine of this argument will be computed.
     * 
     * @return The return value is the result of the calculation.
     */
    template<class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>
    cosine(DifferentiableScalar<Type, Dimension> const& arg0);

  } // namespace numeric

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/differentiableScalar_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_DIFFERENTIABLESCALAR_HH */
