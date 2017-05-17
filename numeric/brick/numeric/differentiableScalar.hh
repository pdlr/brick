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
#include <iostream>

#include <brick/numeric/mathFunctions.hh>
#include <brick/numeric/numericTraits.hh>

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
     **   xx.setDerivative(1.0);
     **   DifferentiableScalar result = xx * xx + sin(xx);
     **   std::cout << "df/dx = " << result.getDerivative() << std::endl;
     ** @endcode
     **
     ** In the previous example, we initialize xx with value 5.  Since
     ** we are computing d/dx, and and the quantity represented by
     ** variable xx is the value of x, the derivative of the value
     ** represented by xx is (trivially) 1.0.  After the value and
     ** derivative of xx have been set, the subsequent mathematical
     ** operations transparently update the derivative of the result.
     ** 
     ** The first template argument of this class determines the
     ** precision with which calculations are to be carried out, while
     ** the second template argument specifies the dimensionality of
     ** the parameter vector with respect to which derivatives should
     ** be tracked.  For example:
     **
     ** @code
     **   // Compute partial derivatives of the function
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
       * value 0, and all partial derivatives equal to 0.
       */
      DifferentiableScalar();


      /** 
       * This default constructor makes a differentiableScalar with a
       * user specified value, and all partial derivatives equal to 0.
       *
       * @param value This argument specifies the scalar value of
       * *this.
       */
      DifferentiableScalar(Type value);


      /** 
       * This default constructor makes a differentiableScalar with a
       * user specified value and partial derivatives explicitly set
       * by the second constructor argument.
       * 
       * @param value This argument specifies the scalar value of
       * *this.
       *
       * @param partialsBegin This iterator points to a
       * sequence of Type specifying a vector of partial derivatives.
       * The number of values in this sequence must be at least equal
       * to the value of class template argument Dimension.
       */
      template <class Iter>
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
       * This type conversion operator allows implicit casting to a
       * scalar value, discarding derivatives.
       * 
       * @return The return value is the scalar value represented by
       * *this.
       */
      template <class OtherType>
      // Explicit conversion operators require C++11. 
#if __cplusplus > 199711L 
      explicit
#endif
      operator OtherType() const {
        return static_cast<OtherType>(this->getValue());
      }


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
      getDerivative() const {return this->getPartialDerivative(0);}


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


    /* ============  Associated class declarations   ============ */

    template <class Type>
    struct NumericTraits< DifferentiableScalar<Type> >
      : public NumericTraitsBase< DifferentiableScalar<Type> > {
    public:

      static inline DifferentiableScalar<Type>
      epsilon() {
        return DifferentiableScalar<Type>(NumericTraits<Type>::epsilon());
      }
      
      static inline bool
      isIntegral() {return NumericTraits<Type>::isIntegral();}
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
     * This operator multiplies a DifferentiableScalar instance by a
     * scalar.
     * 
     * @param arg0 This argument is the DifferentiableScalar instance
     * by which arg1 is to be multiplied.
     * 
     * @param arg1 This argument is the scalar to be multiplied by
     * arg0.
     * 
     * @return The return value is the result of the multiplication.
     */
    template<class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>
    operator*(DifferentiableScalar<Type, Dimension> const& arg0,
              Type const& arg1);


    /** 
     * This operator multiplies a scalar by a DifferentiableScalar
     * instance.
     * 
     * @param arg0 This argument is the scalar by which arg1 is to
     * be multiplied.
     * 
     * @param arg1 This argument is the DifferentiableScalar instance
     * to be multiplied by arg0.
     * 
     * @return The return value is the result of the multiplication.
     */
    template<class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>
    operator*(Type const& arg0,
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
     * This operator divides a DifferentiableScalar instance by a
     * scalar.
     * 
     * @param arg0 This argument is the DifferentiableScalar instance
     * by which arg1 is to be divided.
     * 
     * @param arg1 This argument is the scalar to be divided by
     * arg0.
     * 
     * @return The return value is the result of the division.
     */
    template<class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>
    operator/(DifferentiableScalar<Type, Dimension> const& arg0,
              Type const& arg1);


    /** 
     * This operator divides a scalar by a DifferentiableScalar
     * instance.
     * 
     * @param arg0 This argument is the scalar by which arg1 is to
     * be divided.
     * 
     * @param arg1 This argument is the DifferentiableScalar instance
     * to be divided by arg0.
     * 
     * @return The return value is the result of the division.
     */
    template<class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>
    operator/(Type const& arg0,
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
     * This operator adds a DifferentiableScalar instance to a scalar.
     * 
     * @param arg0 This argument is the DifferentiableScalar instance
     * to which arg1 is to be added.
     * 
     * @param arg1 This argument is the scalar to be added to
     * arg0.
     * 
     * @return The return value is the result of the addition.
     */
    template<class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>
    operator+(DifferentiableScalar<Type, Dimension> const& arg0,
              Type const& arg1);


    /** 
     * This operator adds a scalar to a DifferentiableScalar instance.
     * 
     * @param arg0 This argument is the scalar to which arg1 is to
     * be added.
     * 
     * @param arg1 This argument is the DifferentiableScalar instance
     * to be added to arg0.
     * 
     * @return The return value is the result of the addition.
     */
    template<class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>
    operator+(Type const& arg0,
              DifferentiableScalar<Type, Dimension> const& arg1);


    /** 
     * This operator subtracts two DifferentiableScalar instances.
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
     * This operator subtracts a scalar from a DifferentiableScalar
     * instance.
     * 
     * @param arg0 This argument is the DifferentiableScalar instance
     * from which arg1 is to be subtracted.
     * 
     * @param arg1 This argument is the scalar to be subtracted from
     * arg0.
     * 
     * @return The return value is the result of the subtraction.
     */
    template<class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>
    operator-(DifferentiableScalar<Type, Dimension> const& arg0,
              Type const& arg1);


    /** 
     * This operator subtracts a DifferentiableScalar instance from a
     * scalar.
     * 
     * @param arg0 This argument is the scalar from which arg1 is to
     * be subtracted.
     * 
     * @param arg1 This argument is the DifferentiableScalar instance
     * to be subtracted from arg0.
     * 
     * @return The return value is the result of the subtraction.
     */
    template<class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>
    operator-(Type const& arg0,
              DifferentiableScalar<Type, Dimension> const& arg1);


    /** 
     * The additive inverse operator.
     * 
     * @param arg0 This argument is the the differentiableScalar to
     * be negated.
     * 
     * @return The return value is -arg0.
     */
    template<class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>
    operator-(DifferentiableScalar<Type, Dimension> const& arg0);

      
    /** 
     * Compares the values of two DifferentiableScalar instances.
     * 
     * @param arg0 This argument is the first value to be compared.
     * 
     * @param arg1 This argument is the second value to be compared.
     * 
     * @return The return value is true if arg0.getValue() ==
     * arg1.getValue(), false otherwise.
     */
    template<class Type, uint32_t Dimension>
    bool
    operator==(DifferentiableScalar<Type, Dimension> const& arg0,
               DifferentiableScalar<Type, Dimension> const& arg1)
    {return arg0.getValue() == arg1.getValue();}


    /** 
     * Compares the values of two DifferentiableScalar instances.
     * 
     * @param arg0 This argument is the first value to be compared.
     * 
     * @param arg1 This argument is the second value to be compared.
     * 
     * @return The return value is true if arg0.getValue() !=
     * arg1.getValue(), false otherwise.
     */
    template<class Type, uint32_t Dimension>
    bool
    operator!=(DifferentiableScalar<Type, Dimension> const& arg0,
              DifferentiableScalar<Type, Dimension> const& arg1)
    {return arg0.getValue() != arg1.getValue();}


    /** 
     * Compares the values of two DifferentiableScalar instances.
     * 
     * @param arg0 This argument is the first value to be compared.
     * 
     * @param arg1 This argument is the second value to be compared.
     * 
     * @return The return value is true if arg0.getValue() >
     * arg1.getValue(), false otherwise.
     */
    template<class Type, uint32_t Dimension>
    bool
    operator>(DifferentiableScalar<Type, Dimension> const& arg0,
              DifferentiableScalar<Type, Dimension> const& arg1)
    {return arg0.getValue() > arg1.getValue();}


    /** 
     * Compares the values of two DifferentiableScalar instances.
     * 
     * @param arg0 This argument is the first value to be compared.
     * 
     * @param arg1 This argument is the second value to be compared.
     * 
     * @return The return value is true if arg0.getValue() >=
     * arg1.getValue(), false otherwise.
     */
    template<class Type, uint32_t Dimension>
    bool
    operator>=(DifferentiableScalar<Type, Dimension> const& arg0,
               DifferentiableScalar<Type, Dimension> const& arg1)
    {return arg0.getValue() >= arg1.getValue();}


    /** 
     * Compares the values of two DifferentiableScalar instances.
     * 
     * @param arg0 This argument is the first value to be compared.
     * 
     * @param arg1 This argument is the second value to be compared.
     * 
     * @return The return value is true if arg0.getValue() <
     * arg1.getValue(), false otherwise.
     */
    template<class Type, uint32_t Dimension>
    bool
    operator<(DifferentiableScalar<Type, Dimension> const& arg0,
              DifferentiableScalar<Type, Dimension> const& arg1)
    {return arg0.getValue() < arg1.getValue();}


    /** 
     * Compares the values of two DifferentiableScalar instances.
     * 
     * @param arg0 This argument is the first value to be compared.
     * 
     * @param arg1 This argument is the second value to be compared.
     * 
     * @return The return value is true if arg0.getValue() <=
     * arg1.getValue(), false otherwise.
     */
    template<class Type, uint32_t Dimension>
    bool
    operator<=(DifferentiableScalar<Type, Dimension> const& arg0,
               DifferentiableScalar<Type, Dimension> const& arg1)
    {return arg0.getValue() <= arg1.getValue();}
    

    /** 
     * This function computes the absolute value of a
     * DifferentiableScalar instance, with partial derivatives.
     * 
     * @param arg0 The absolute value of this argument will be
     * computed.
     * 
     * @return The return value is the result of the calculation.
     */
    template<class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>
    absoluteValue(DifferentiableScalar<Type, Dimension> const& arg0);
    
    /** 
     * This function computes the cosine of a DifferentiableScalar
     * instance, with partial derivatives.
     * 
     * @param arg0 The cosine of this argument will be computed.
     * 
     * @return The return value is the result of the calculation.
     */
    template<class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>
    cosine(DifferentiableScalar<Type, Dimension> const& arg0);


    /** 
     * This function computes the sine of a DifferentiableScalar
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
     * This function computes the square root of a
     * DifferentiableScalar instance, with partial derivatives.
     * 
     * @param arg0 The square of this argument will be computed.
     * 
     * @return The return value is the result of the calculation.
     */
    template<class Type, uint32_t Dimension>
    DifferentiableScalar<Type, Dimension>
    squareRoot(DifferentiableScalar<Type, Dimension> const& arg0);


    /** 
     * Stream output operator.
     * 
     * @param stream This argument is the output stream to which arg1
     * will be written.
     * 
     * @param arg1 This argument will be formatted to the output stream.
     * 
     * @return The return value is the output stream after writing.
     */
    template<class Type, uint32_t Dimension>
    std::ostream&
    operator<<(std::ostream& stream,
               DifferentiableScalar<Type, Dimension> const& arg1);
    
  } // namespace numeric

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/differentiableScalar_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_DIFFERENTIABLESCALAR_HH */
