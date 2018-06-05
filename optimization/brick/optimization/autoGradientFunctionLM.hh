/**
**********************************************************************
* @file brick/optimization/autoGradientFunctionLM.hh
*
* Header file declaring AutoGradientFunctionLM class template.
*
* Copyright (C) 2018 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
**********************************************************************
**/

#ifndef BRICK_OPTIMIZATION_AUTOGRADIENTFUNCTIONLM_HH
#define BRICK_OPTIMIZATION_AUTOGRADIENTFUNCTIONLM_HH

#include <functional>
#include <brick/numeric/array1D.hh>
#include <brick/numeric/array2D.hh>


namespace brick {

  namespace optimization {

    /**
     ** The AutoGradientFunctionLM class template is derived from
     ** std::unary_function, and adds one additional member function
     ** for computing the gradient and approximate Hessian matrix of
     ** sum-of-squares error functions.  This is primarily useful when
     ** you need to minimize a function using the OptimizerLM class
     ** template.  The gradient and Hessian values are computed using
     ** automatic differentiation, as described below.
     **
     ** Template argument SSDFunction must implement three public
     ** member functions as illustrated in the following code snippet:
     **
     ** @code
     **   struct MySSDFunction {
     **     // The input sequence of this->apply() must have at least
     **     // this many valid elements.  This function must return
     **     // a number that is no greater than AutoGradientFunctionLM
     **     // template argument NumberOfArguments.
     **     std::size_t getNumberOfArguments();
     **
     **     // The output sequence of this->apply() must be able to accept
     **     // at least this many elements.
     **     std::size_t getNumberOfErrorTerms();
     **
     **     // Takes input arguments from the sequence *argsBegin,
     **     // *(argsBegin + 1), *(argsBegin + 2), etc., and computes
     **     // a sequence of residual (error) terms that are stored in
     **     // the sequence *resultBegin, *(resultBegin + 1), etc.
     **     template <class InputIter, class OutputIter>
     **     void apply(InputIter argsBegin, OutputIter resultBegin);
     **     };
     ** @endcode
     **
     ** OptimizerLM will attempt to find the input argument values
     ** that minimize the sum of squares of the output residuals.
     **
     ** AutoGradientFunctionLM will make a copy of its SSDFunction
     ** constructor argument, and call the apply() member function of
     ** this copy in two different ways:
     **
     ** - It will be called with iterater arguments that dereference
     **   to elements of type Scalar.  When called in this way, the
     **   output sequence supplies residual values that are squared
     **   and summed to obtain the total sum-squared error.
     **
     ** - It will be called with iterater arguments that dereference
     **   to elements of type DifferentiableScalar<Scalar, NumberOfArguments>.
     **   When called in this way, the output sequence supplies
     **   automatically computed partial derivatives for each residual
     **   value.  These partial derivatives are then used to compute the
     **   gradient vector and approximate Hessian matrix.
     **
     ** The second way that apply() is called implies that template
     ** argument NumberOfArguments must be greater than or equal to
     ** the number of arguments required by SSDFunction.  If template
     ** argument NumberOfArguments is smaller than this number, it is
     ** an error.  If NumberOfArguments is larger than the number of
     ** arguments required by SSDFunction, it is not an error, but
     ** there is a performance cost associated with computing extra
     ** partial derivatives.
     **
     ** In order for partial derivatives to be automatically computed,
     ** the implementation of SSDFunction::apply() must know the type
     ** of the sequence elements and make sure partial derivatives are
     ** propagated through all residual calculations.  The following
     ** example shows how to do this.
     **
     ** @code
     **   template <class InputIter, class OutputIter>
     **   void MySSDFunction::apply(InputIter argsBegin, OutputIter resultBegin)
     **   {
     **     // Deduce the type of scalar we're working with.  After this
     **     // line, FloatType is typedef'd to AutoGradientFunctionLM::Scalar
     **     // if apply() is being called in the first way (to compute
     **     // residuals only), or
     **     // DifferentiableScalar<AutoGradientFunctionLM::Scalar,
     **     //                      AutoGradientFunctionLM::NumberOfArguments>
     **     // if apply() is being called in the second way.
     **     typedef typename std::remove_reference<decltype(*resultBegin)>::type
     **     FloatType;
     **
     **     // Through the rest of the function, any value that depends
     **     // on an input argument is declared FloatType (instead of,
     **     // say, double.  In this example, the numbers 5.0, etc. don't
     **     // have to be declared FloatType because they don't depend on
     **     // an input argument.  This example works as long as there
     **     // are operator*(double, DifferentiableScalar<...>) and
     **     // DifferentiableScalar<...>::operator=(double) operators
     **     // defined.  If there aren't, it may be necessary to explicitly
     **     // construct FloatType(5.0), FloatType(0.0), and FloatType(-2.0)
     **     // instances.
     **     FloatType result0 = 5.0 * (*argsBegin) * (*argsBegin);
     **     FloatType result1 = 0.0;
     **     FloatType result2 = -2.0 * (*argsBegin);
     **     ++argsBegin;
     **     [...]
     **     *(resultBegin++) = result0;
     **     *(resultBegin++) = result1;
     **     [...]
     **   }
     ** @endcode
     **
     ** Template argument Scalar specifies the precision with which
     ** internal calculations will be conducted and the type that will
     ** be used to represent the sum-of-squares error.  Reasonable choices
     ** are double and float.
     **
     ** Here's a usage example:
     **
     ** @code
     **   constexpr int numberOfArguments = 3;
     **   typedef AutoGradientFunctionLM<MySSDFunction, numberOfArguments>
     **     GradientFunction;
     **   
     **   MySSDFunction ssdFunction;
     **   assert(numberOfArguments >= ssdFunction.getNumberOfArguments());
     **   GradientFunction gradientFunction(ssdFunction);
     **   OptimizerLM<GradientFunction> optimizer(gradientFunction);
     **   optimizer.setStartPoint(myStartPoint);
     **   myResult = optimizer.optimum();
     ** @endcode
     **/
    template <class SSDFunction, int NumberOfArguments,
              class Scalar = brick::common::Float64>
    class AutoGradientFunctionLM
      : public std::unary_function<brick::numeric::Array1D<Scalar>, Scalar>
    {
    public:
      /** 
       * The default constructor simply uses the default SSDFunction.
       */
      AutoGradientFunctionLM();
      

      /** 
       * Constructor.
       *
       * @param ssdFunction This argument is the function object to be
       * adapted.
       */
      AutoGradientFunctionLM(SSDFunction const& ssdFunction);
      

      /** 
       * Destructor.
       */
      virtual ~AutoGradientFunctionLM() {}

      
      /** 
       * This operator evaluates the sum-of-squares error at the
       * specified point.
       * 
       * @param theta The point at which to evaluate the function.
       * @return The function value at theta.
       */
      Scalar
      operator()(brick::numeric::Array1D<Scalar> const& theta);

      
      /** 
       * This method approximates the gradient and Hessian matrix of
       * this->operator().  The Jacobian of SSDFunction::operator()()
       * is computed by automatic differentiation, the gradient is
       * computed directly from the Jacobian, and the Hessian matrix
       * is computed from the Jacobian using the Levenberg-Marquardt
       * approximation.
       *
       * @param theta The point around which to compute the gradient
       * and hessian.
       *
       * @param dEdX This argument is used to return the computed
       * gradient to the calling context.  Note that the type if this
       * argument should be controlled by a template parameter.  This
       * will be fixed when we finally get around to overhauling the
       * Optimizer* template parameters.
       * 
       * @param d2EdX2 This argument is used to return the computed
       * Hessian matrix to the calling context.  Note that the type if
       * this argument should be controlled by a template parameter.
       * This will be fixed when we finally get around to overhauling
       * the Optimizer* template parameters.
       */
      void
      computeGradientAndHessian(brick::numeric::Array1D<Scalar> const& theta,
                                brick::numeric::Array1D<Scalar>& dEdX,
                                brick::numeric::Array2D<Scalar>& d2EdX2);
      
    private:
      SSDFunction m_ssdFunction;
      Scalar m_epsilon;

    }; // class AutoGradientFunctionLM

  } // namespace optimization

} // namespace brick


/*******************************************************************
 * Member function definitions follow.  This would be a .C file
 * if it weren't templated.
 *******************************************************************/

#include <brick/numeric/differentiableScalar.hh>
#include <brick/numeric/utilities.hh>

namespace brick {

  namespace optimization {

    // Default constructor.
    template <class SSDFunction, int NumberOfArguments, class Scalar>
    AutoGradientFunctionLM<SSDFunction, NumberOfArguments, Scalar>::
    AutoGradientFunctionLM()
      : m_ssdFunction()
    {
      // Empty.
    }

    
    // Constructor.
    template <class SSDFunction, int NumberOfArguments, class Scalar>
    AutoGradientFunctionLM<SSDFunction, NumberOfArguments, Scalar>::
    AutoGradientFunctionLM(SSDFunction const& ssdFunction)
      : m_ssdFunction(ssdFunction)
    {
      // Empty.
    }

    
    // This operator evaluates the sum-of-squares error at the
    // specified point.
    template <class SSDFunction, int NumberOfArguments, class Scalar>
    Scalar
    AutoGradientFunctionLM<SSDFunction, NumberOfArguments, Scalar>::
    operator()(brick::numeric::Array1D<Scalar> const& theta)
    {
      Scalar result = Scalar(0.0);

      if(theta.size() < this->m_ssdFunction.getNumberOfArguments()) {
        std::ostringstream message;
        message << "SSDFunction requires "
                << this->m_ssdFunction.getNumberOfArguments()
                << " arguments, but input array has only "
                << theta.size() << " elements.";
        BRICK_THROW(brick::common::IndexException,
                    "AutoGradientFunctionLM::operator()()",
                    message.str().c_str());
      }
      
      unsigned int const numTerms = this->m_ssdFunction.getNumberOfErrorTerms();
      brick::numeric::Array1D<Scalar> errorTerms(numTerms);
      this->m_ssdFunction.apply(theta.begin(), errorTerms.begin());
      for(unsigned int ii = 0; ii < numTerms; ++ii) {
        result += errorTerms[ii] * errorTerms[ii];
      }
      return result;
    }


    // This method computes the gradient and Hessian matrix of
    // this->operator().
    template <class SSDFunction, int NumberOfArguments, class Scalar>
    void
    AutoGradientFunctionLM<SSDFunction, NumberOfArguments, Scalar>::
    computeGradientAndHessian(brick::numeric::Array1D<Scalar> const& theta,
                              brick::numeric::Array1D<Scalar>& dEdX,
                              brick::numeric::Array2D<Scalar>& d2EdX2)
    {
      typedef brick::numeric::DifferentiableScalar<Scalar, NumberOfArguments>
        DiffScalar;

      if(theta.size() < this->m_ssdFunction.getNumberOfArguments()) {
        std::ostringstream message;
        message << "SSDFunction requires "
                << this->m_ssdFunction.getNumberOfArguments()
                << " arguments, but input array has only "
                << theta.size() << " elements.";
        BRICK_THROW(brick::common::IndexException,
                    "AutoGradientFunctionLM::computeGradientAndHessian()",
                    message.str().c_str());
      }
      if(theta.size() < NumberOfArguments) {
        std::ostringstream message;
        message << "SSDFunction requires "
                << this->m_ssdFunction.getNumberOfArguments()
                << " arguments, but AutoGradientFunctionLM template argument "
                << "NumberOfArguments is set to only " << NumberOfArguments
                << ".";
        BRICK_THROW(brick::common::IndexException,
                    "AutoGradientFunctionLM::computeGradientAndHessian()",
                    message.str().c_str());
      }
      
      // Get oriented.
      unsigned int const numTerms = this->m_ssdFunction.getNumberOfErrorTerms();

      // Copy parameters to a type that tracks its own first derivatives.
      brick::numeric::Array1D<DiffScalar> arguments(theta.size());
      brick::numeric::Array1D<Scalar> derivatives(theta.size());
      derivatives = Scalar(0);
      for(std::size_t ii = 0; ii < arguments.size(); ++ii) {

        // Set the partial derivative of this particular argument with respect
        // to itself to 1.0, while the partial with respect to all other
        // arguments stays at 0.0.
        derivatives[ii] = Scalar(1);
        arguments[ii] = DiffScalar(theta[ii], derivatives.begin());

        // Reset partial derivatives to zero in preparation for the next
        // argument.
        derivatives[ii] = Scalar(0);
      }

      // Compute residuals (and their first derivatives).
      brick::numeric::Array1D<DiffScalar> errorTerms(numTerms);
      this->m_ssdFunction.apply(arguments.begin(), errorTerms.begin());

      // Construct an array to hold 1st derivatives.
      brick::numeric::Array2D<Scalar> jacobian(NumberOfArguments, numTerms);

      for(unsigned int rr = 0; rr < NumberOfArguments; ++rr) {
        for(unsigned int cc = 0; cc < numTerms; ++cc) {
          jacobian(rr, cc) = errorTerms[cc].getPartialDerivative(rr);
        }
      }

      // We want to do a matrix multiplication, but we can't multiply
      // matrices of dissimilar types without a contortion.  For now,
      // we just make a copy of errorTerms.
      brick::numeric::Array1D<Scalar> errorTermsNoDerivative(numTerms);
      for(unsigned int cc = 0; cc < numTerms; ++cc) {
        errorTermsNoDerivative[cc] = errorTerms[cc].getValue();
      }
      
      dEdX = brick::numeric::matrixMultiply<Scalar>(
        jacobian, errorTermsNoDerivative);
      dEdX *= 2.0;
      
      // Compute Hession estimate.
      d2EdX2 = brick::numeric::matrixMultiply<Scalar>(
        jacobian, jacobian.transpose());
      d2EdX2 *= 2.0;
    }
    
  } // namespace optimization

} // namespace brick

#endif /* #ifndef BRICK_OPTIMIZATION_AUTOGRADIENTFUNCTION_HH */
