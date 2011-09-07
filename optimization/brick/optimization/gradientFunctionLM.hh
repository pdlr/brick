/**
**********************************************************************
* @file brick/optimization/gradientFunctionLM.hh
*
* Header file declaring GradientFunctionLM class template.
*
* Copyright (C) 2010-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
**********************************************************************
**/

#ifndef BRICK_OPTIMIZATION_GRADIENTFUNCTIONLM_HH
#define BRICK_OPTIMIZATION_GRADIENTFUNCTIONLM_HH

#include <functional>

namespace brick {

  namespace optimization {

    /**
     ** The GradientFunctionLM class template is derived from
     ** std::unary_function, and adds one additional member function
     ** for computing the gradient and Hessian matrix of
     ** sum-of-squares error functions.  This is primarily useful when
     ** you need to minimize a function using the OptimizerLM class
     ** template.
     **
     ** WARNING: This interface is not stable, and is likely to
     ** change.  In particular, GradientFunctionLM currently forces
     ** you to compute your gradient numerically.  Future updates will
     ** allow you to do this symbolically if you like, but may slightl
     ** break existing code.
     **
     ** Template argument SSDFunctor is assumed to be a subclass of
     ** std::unary_function.  SSDFunctor::argument_type is assumed to be
     ** a vector or 1D array type which supports the following
     ** interface:
     **
     **   argument_type(size_t N): construct an N-element vector.
     **   size_t size(): return the number of elements in the vector.
     **   argument_type::element_type& operator[](size_t i): return a
     **                        reference to the (i)th element of the array.
     **   const argument_type::element_type& operator[](size_t i) const:
     **                        return a const reference to the (i)th element
     **                        of the array.
     **
     ** It is further assumed that element type of argument_type is a
     ** continuous scalar, and can be implicitly cast to and from
     ** the type specified by template argument Scalar.
     **
     ** SSDFunctor::result_type is _also_ assumed to be a vector or 1D
     ** array type supporting the same interface as
     ** SSDFunctor::argument_type.
     **
     ** The result of SSDFunctor::operator()(argument_type const&) is an
     ** array of error terms.  Generally, the sum of squares of these
     ** error terms is what you want to minimize using OptimizerLM.
     **
     ** Template argument Scalar specifies the precision with which
     ** internal calculations will be conducted and the type that will
     ** be used to represent the sum-of-squares error.
     **
     ** Here's a usage example:
     **
     ** @code
     **   typedef GradientFunctionLM<MySSDFunction> GradFunctor;
     **   MySSDFunction functor();
     **   GradFunctor gradFunctor(functor);
     **   OptimizerLM<GradFunctor> optimizer(gradFunctor);
     **   optimizer.setStartPoint(myStartPoint);
     **   myResult = optimizer.optimum();
     ** @endcode
     **/
    template <class SSDFunctor, class Scalar = common::Float64>
    class GradientFunctionLM
      : public std::unary_function<typename SSDFunctor::argument_type, Scalar>
    {
    public:
      /** 
       * Constructor.
       *
       * @param functor This argument is the function object to be
       * adapted.
       *
       * @param epsilon If the computeGradientAndHessian() method is
       * not overridden in a subclass, the gradient will be computed
       * by using symmetric divided differences with a total step size
       * of 2 * epsilon.
       */
      GradientFunctionLM(const SSDFunctor& functor, Scalar epsilon=1.0e-6);
      

      /** 
       * Destructor.
       */
      virtual ~GradientFunctionLM() {}

      
      /** 
       * This operator evaluates the sum-of-squares error at the
       * specified point.
       * 
       * @param theta The point at which to evaluate the function.
       * @return The function value at theta.
       */
      Scalar
      operator()(const typename SSDFunctor::argument_type& theta);

      
      /** 
       * This method approximates the gradient and Hessian matrix of
       * this->operator().  The Jacobian of SSDFunctor::operator()()
       * is computed by divided differences, the gradient is computed
       * directly from the Jacobian, and the Hessian matrix is
       * computed from the Jacobian using the Levenberg-Marquardt
       * approximation.  This method should often be overridden by a
       * subclass.
       *
       * This function will throw ValueException if you set
       * constructor argument epsilon small enough that, when added to
       * the elements of theta, it gets completely rounded away.
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
      computeGradientAndHessian(
        typename SSDFunctor::argument_type const& theta,
        numeric::Array1D<common::Float64>& dEdX,
        numeric::Array2D<common::Float64>& d2EdX2);
      
    private:
      SSDFunctor m_functor;
      Scalar m_epsilon;

    }; // class GradientFunctionLM

  } // namespace optimization

} // namespace brick


/*******************************************************************
 * Member function definitions follow.  This would be a .C file
 * if it weren't templated.
 *******************************************************************/

namespace brick {

  namespace optimization {

    // Constructor.
    template <class SSDFunctor, class Scalar>
    GradientFunctionLM<SSDFunctor, Scalar>::
    GradientFunctionLM(const SSDFunctor& functor, Scalar epsilon)
      : m_functor(functor),
        m_epsilon(epsilon)
    {
      if(epsilon == 0.0) {
        BRICK_THROW(brick::common::ValueException, 
		    "GradientFunctionLM::GradientFunctionLM()",
		    "Invalid value (0.0) for argument epsilon.");
      }
    }

    
    // This operator evaluates the sum-of-squares error at the
    // specified point.
    template <class SSDFunctor, class Scalar>
    Scalar
    GradientFunctionLM<SSDFunctor, Scalar>::
    operator()(const typename SSDFunctor::argument_type& theta)
    {
      Scalar result = 0.0;
      typename SSDFunctor::result_type errorTerms = this->m_functor(theta);
      unsigned int numTerms = errorTerms.size();
      for(unsigned int ii = 0; ii < numTerms; ++ii) {
        result += errorTerms[ii] * errorTerms[ii];
      }
      return result;
    }


    // This method approximates the gradient and Hessian matrix of
    // this->operator().
    template <class SSDFunctor, class Scalar>
    void
    GradientFunctionLM<SSDFunctor, Scalar>::
    computeGradientAndHessian(
      typename SSDFunctor::argument_type const& theta,
      numeric::Array1D<common::Float64>& dEdX,
      numeric::Array2D<common::Float64>& d2EdX2)
    {
      // Note(xxx): Move all of these Array?D constructors out of this
      // function so as to no to gobs of extra new/delete cycles.
      
      // Get oriented.
      typename SSDFunctor::result_type errorTerms = this->m_functor(theta);

      // Construct an array to hold 1st derivatives.
      unsigned int numParameters = theta.size();
      unsigned int numTerms = errorTerms.size();
      numeric::Array2D<Float64> jacobian(numParameters, numTerms);

      // Compute Jacobian of error terms by divided differences.
      typename SSDFunctor::argument_type thetaPlus(numParameters);
      typename SSDFunctor::argument_type thetaMinus(numParameters);
      numeric::Array1D<Scalar> derivative(numTerms);
      for(unsigned int ii = 0; ii < numParameters; ++ii) {
        thetaPlus[ii] = theta[ii];
        thetaMinus[ii] = theta[ii];
      }
      for(unsigned int ii = 0; ii < numParameters; ++ii) {
        // Set up divided differences.
        thetaPlus[ii] += this->m_epsilon;
        thetaMinus[ii] -= this->m_epsilon;
        Scalar delta = thetaPlus[ii] - thetaMinus[ii];

        // Compute divided differences.
        typename SSDFunctor::result_type errorTermsPlus =
          this->m_functor(thetaPlus);
        typename SSDFunctor::result_type errorTermsMinus =
          this->m_functor(thetaMinus);
        for(unsigned int jj = 0; jj < numTerms; ++jj) {
          derivative[jj] = (errorTermsPlus[jj] - errorTermsMinus[jj]) / delta;
        }
        jacobian.getRow(ii).copy(derivative);

        // Reset our parameter changes.
        thetaPlus[ii] = theta[ii];
        thetaMinus[ii] = theta[ii];
      }

      // Compute gradient estimate from Jacobian.
      numeric::Array1D<Float64> errorTermsCopy(numTerms);
      for(unsigned int jj = 0; jj < numTerms; ++jj) {
        errorTermsCopy[jj] = errorTerms[jj];
      }
      dEdX = numeric::matrixMultiply(jacobian, errorTerms);
      dEdX *= 2.0;
      
      // Compute Hession estimate.
      d2EdX2 = numeric::matrixMultiply(jacobian, jacobian.transpose());
      d2EdX2 *= 2.0;
    }
    
  } // namespace optimization

} // namespace brick

#endif /* #ifndef BRICK_OPTIMIZATION_GRADIENTFUNCTION_HH */
