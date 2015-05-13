/**
***************************************************************************
* @file brick/common/complexNumber.hh
*
* Header file declaring ComplexNumber class template.
*
* Copyright (C) 2015 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_COMMON_COMPILETIMESTAMP_HH
#define BRICK_COMMON_COMPILETIMESTAMP_HH

#include <brick/common/constants.hh>
#include <brick/common/mathFunctions.hh>

namespace brick {

  namespace common {

    /**
     ** The ComplexNumber class template is a temporary stand-in for
     ** std::complex, until the std::complex API stabilizes, and until
     ** std::complex is expanded to permit arbitrary types as template
     ** arguments.
     **/
    template <class Type>
    class ComplexNumber
    {
    public:
      /** 
       * The default constructor sets value to 0 + 0i.
       */
      ComplexNumber()
        : m_realPart(0), m_imaginaryPart(0) {}

      
      /** 
       * This constructor allows real and imaginary parts to be
       * explicitly set.
       * 
       * @param rr This specifies what the real part of the complex number
       * will be set to.
       * 
       * @param ii This specifies what the imaginary part of the complex number
       * will be set to.
       */
      ComplexNumber(Type const& rr, Type const& ii)
        : m_realPart(rr), m_imaginaryPart(ii) {}
      
     
      /** 
       * Destroys the ComplexNumber instance.
       */
      ~ComplexNumber() {}


      /** 
       * Returns the real part of the complex number by value.
       * 
       * @return The return value is the real component of the complex
       * number.
       */
      Type const&
      getRealPart() const {return this->m_realPart;}
    

      /** 
       * Returns the imaginary part of the complex number by value.
       * 
       * @return The return value is the imaginary component of the
       * complex number.
       */
      Type const&
      getImaginaryPart() const {return this->m_imaginaryPart;}
    

      /** 
       * Sets the real part of the complex number.
       *
       * @param rr This specifies what the real part of the complex
       * number will be set to.
       *
       * @return The return value is the real component of the complex
       * number.
       */
      Type const&
      setRealPart(Type const& rr) {
        this->m_realPart = rr;
        return this->m_realPart;
      }


      /** 
       * Sets the imaginary part of the complex number.
       *
       * @param ii This specifies what the imaginary part of the
       * complex number will be set to.
       *
       * @return The return value is the imaginary component of the
       * complex number.
       */
      Type const&
      setImaginaryPart(Type const& ii) {
        this->m_imaginaryPart = ii;
        return this->m_imaginaryPart;
      }
      

      /** 
       * Sets the real and imaginary parts of the complex number.
       *
       * @param rr This specifies what the real part of the complex
       * number will be set to.
       *
       * @param ii This specifies what the imaginary part of the
       * complex number will be set to.
       *
       * @return The return value is a const reference to *this.
       */
      ComplexNumber const&
      setValue(Type const& rr, Type const& ii) {
        this->m_realPart = rr;
        this->m_imaginaryPart = ii;
        return *this;
      }
      
    private:

      Type m_realPart;
      Type m_imaginaryPart;
      
    };


    template <class Type>
    ComplexNumber<Type>
    squareRoot(ComplexNumber<Type> const& arg0) {
      // We solve this using De Moivre's Identity, which states that
      // for any real number x and integer n:
      //
      // @code
      //   (cos(x) + i*sin(x))^n = cos(n*x) + i*sin(n*x)
      // @endcode
      if(!arg0.getRealPart() && !arg0.getImaginaryPart()) {
        // Square root of 0 + 0i is still 0 + 0i.
        return arg0;
      }

      // Convert into the form r*(cos(x) + i*sin(x)).
      Type xx = arctangent2(arg0.getImaginaryPart(), arg0.getRealPart());
      while(xx >= constants::pi) {xx -= constants::twoPi;}
      while(xx < -(constants::pi)) {xx += constants::twoPi;}
      Type rr(0);
      if((xx < constants::piOverFour && xx >= -(constants::piOverFour))
         || xx >= constants::threePiOverFour
         || xx < -(constants::threePiOverFour)) {
        rr = arg0.getRealPart() / cosine(xx);
      } else {
        rr = arg0.getImaginaryPart() / sine(xx);
      }

      // OK.  Now return sqrt(r) * (cos(x/2) + i*sin(x/2)).
      Type rootR = squareRoot(rr);
      Type xOverTwo = xx / Type(2);
      return ComplexNumber<Type>(rootR * cosine(xOverTwo),
                                 rootR * sine(xOverTwo)); 
    }
    
    
    template <class Type>
    ComplexNumber<Type>
    operator-(ComplexNumber<Type> const& arg0) {
      return ComplexNumber<Type>(
        -arg0.getRealPart(), -arg0.getImaginaryPart());
    }

    
    template <class Type>
    ComplexNumber<Type>
    operator+(ComplexNumber<Type> const& arg0,
              ComplexNumber<Type> const& arg1) {
      return ComplexNumber<Type>(
        arg0.getRealPart() + arg1.getRealPart(),
        arg0.getImaginaryPart() + arg1.getImaginaryPart());
    }

    
    template <class Type>
    ComplexNumber<Type>
    operator-(ComplexNumber<Type> const& arg0,
              ComplexNumber<Type> const& arg1) {
      return ComplexNumber<Type>(
        arg0.getRealPart() - arg1.getRealPart(),
        arg0.getImaginaryPart() - arg1.getImaginaryPart());
    }

    
    template <class Type>
    ComplexNumber<Type>
    operator*(ComplexNumber<Type> const& arg0,
              ComplexNumber<Type> const& arg1) {
      return ComplexNumber<Type>(
        arg0.getRealPart() * arg1.getRealPart()
        - arg0.getImaginaryPart() * arg1.getImaginaryPart(),
        arg0.getRealPart() * arg1.getImaginaryPart()
        + arg0.getImaginaryPart() * arg1.getRealPart());
    }

    
    template <class Type>
    ComplexNumber<Type>
    operator/(ComplexNumber<Type> const& arg0,
              ComplexNumber<Type> const& arg1) {
      Type denominator = (arg1.getRealPart() * arg1.getRealPart()
                          + arg1.getImaginaryPart() * arg1.getImaginaryPart());
      Type realPart = ((arg0.getRealPart() * arg1.getRealPart()
                        + arg0.getImaginaryPart() * arg1.getImaginaryPart())
                       / denominator);
      Type imaginaryPart = ((arg0.getImaginaryPart() * arg1.getRealPart()
                             - arg0.getRealPart() * arg1.getImaginaryPart())
                            / denominator);
      return ComplexNumber<Type>(realPart, imaginaryPart);
    }

    
    template <class Type>
    ComplexNumber<Type>
    operator+(Type const& arg0, ComplexNumber<Type> const& arg1) {
      return ComplexNumber<Type>(
        arg0 + arg1.getRealPart(), arg1.getImaginaryPart());
    }

    
    template <class Type>
    ComplexNumber<Type>
    operator-(Type const& arg0, ComplexNumber<Type> const& arg1) {
      return ComplexNumber<Type>(
        arg0 - arg1.getRealPart(), -(arg1.getImaginaryPart()));
    }

    
    template <class Type>
    ComplexNumber<Type>
    operator*(Type const& arg0, ComplexNumber<Type> const& arg1) {
      return ComplexNumber<Type>(
        arg0 * arg1.getRealPart(), arg0 * arg1.getImaginaryPart());
    }

    
    template <class Type>
    ComplexNumber<Type>
    operator/(Type const& arg0, ComplexNumber<Type> const& arg1) {
      Type denominator = (arg1.getRealPart() * arg1.getRealPart()
                          + arg1.getImaginaryPart() * arg1.getImaginaryPart());
      Type realPart = (arg0 * arg1.getRealPart()) / denominator;
      Type imaginaryPart = ((-arg0) * arg1.getImaginaryPart()) / denominator;
      return ComplexNumber<Type>(realPart, imaginaryPart);
    }

    
    template <class Type>
    ComplexNumber<Type>
    operator+(ComplexNumber<Type> const& arg0, Type const& arg1) {
      return ComplexNumber<Type>(
        arg0.getRealPart() + arg1, arg0.getImaginaryPart());
    }

    
    template <class Type>
    ComplexNumber<Type>
    operator-(ComplexNumber<Type> const& arg0, Type const& arg1) {
      return ComplexNumber<Type>(
        arg0.getRealPart() - arg1, arg0.getImaginaryPart());
    }

    
    template <class Type>
    ComplexNumber<Type>
    operator*(ComplexNumber<Type> const& arg0, Type const& arg1) {
      return ComplexNumber<Type>(
        arg0.getRealPart() * arg1, arg0.getImaginaryPart() * arg1);
    }

    
    template <class Type>
    ComplexNumber<Type>
    operator/(ComplexNumber<Type> const& arg0, Type const& arg1) {
      return ComplexNumber<Type>(
        arg0.getRealPart() / arg1, arg0.getImaginaryPart() / arg1);
    }

  } // namespace common
    
} // namespace brick


#endif /* #ifndef BRICK_COMMON_COMPILETIMESTAMP_HH */
