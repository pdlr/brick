/**
***************************************************************************
* @file brick/optimization/optimizerLM.hh
*
* Header file declaring OptimizerLM class.
*
* (C) Copyright 2003-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying LICENSE file for details.
*
***************************************************************************
**/

#ifndef BRICK_OPTIMIZATION_OPTIMIZERLM_HH
#define BRICK_OPTIMIZATION_OPTIMIZERLM_HH

#include <vector>
#include <brick/common/types.hh>
#include <brick/linearAlgebra/linearAlgebra.hh>
#include <brick/numeric/utilities.hh>
#include <brick/optimization/optimizer.hh>
#include <brick/optimization/optimizerLineSearch.hh>

namespace brick {

  namespace optimization {

    /**
     ** OptimizerLM implements the Levenberg-Marquardt nonlinear
     ** least-squares minimization algorithm, as described in [1].  This
     ** algorithm seeks the parameter value which minimizes the
     ** objective function.  The template parameter (Functor) defines
     ** the type to use as the objective function of the minimization,
     ** and must support the GradientFunction interface.
     **
     ** [1] W. H. Press et al., Numerical Recipes in C The Art of
     ** Scientific Computing, Cambridge University Press, 1988.
     **/
    template <class Functor>
    class OptimizerLM
      : public Optimizer<Functor>
    {
    public:
      // Typedefs for convenience
      typedef typename Functor::argument_type argument_type;
      typedef typename Functor::result_type result_type;


      /**
       * The default constructor sets parameters to reasonable values
       * for functions which take values and arguments in the "normal"
       * range of 0 to 100 or so.
       */
      OptimizerLM();


      /**
       * This constructor specifies the specific Functor instance to
       * use.  Using this constructor exclusively avoids the danger of
       * calling optimalValue() or optimum() before a Functor instance
       * has been specified.
       *
       * @param functor A copy of this argument will be stored
       * internally for use in optimization.
       */
      explicit OptimizerLM(const Functor& functor);


      /** 
       * Copy constructor.  This constructor deep copies its argument.
       * 
       * @param source The OptimizerLM instance to be copied.
       */
      OptimizerLM(const OptimizerLM& source);


      /**
       * The destructor destroys the class instance and deallocates any
       * associated storage.
       */
      virtual
      ~OptimizerLM();

    
//     /** 
//      * This method returns the number of function calls required to
//      * complete the previous minimization.  If the minimization
//      * parameter "restarts" is 0, there will be only one element in
//      * the returned vector.  If restarts is greater than 0, the first
//      * element of the return value will reflect the number of function
//      * calls in the initial minimization, and subsequent numbers will
//      * reflect the number of function calls in the following restarted
//      * minimizations.  If a valid minimization has not been performed
//      * since the last update to startPoint, parameters, etc., then the
//      * return value will be an empty vector.
//      *
//      * @return a vector of function call counts.
//      */
//     virtual std::vector<size_t>
//     getNumberOfFunctionCalls() {return this->m_functionCallCount;}

    
//     /** 
//      * This method returns the number of gradient calls required to
//      * complete the previous minimization.  If the minimization
//      * parameter "restarts" is 0, there will be only one element in
//      * the returned vector.  If restarts is greater than 0, the first
//      * element of the return value will reflect the number of gradient
//      * calls in the initial minimization, and subsequent numbers will
//      * reflect the number of gradient calls in the following restarted
//      * minimizations.  If a valid minimization has not been performed
//      * since the last update to startPoint, parameters, etc., then the
//      * return value will be an empty vector.
//      *
//      * @return a vector of gradient call counts.
//      */
//     virtual std::vector<size_t>
//     getNumberOfGradientCalls() {return this->m_gradientCallCount;}
    
//     /** 
//      * This method returns the number of iterations required to
//      * complete the previous minimization.  If the minimization
//      * parameter "restarts" is 0, there will be only one element in
//      * the returned vector.  If restarts is greater than 0, the first
//      * element of the return value will reflect the number of
//      * iterations in the initial minimization, and subsequent numbers
//      * will reflect the number of iterations in the following
//      * restarted minimizations.  If a valid minimization has not been
//      * performed since the last update to startPoint, parameters,
//      * etc., then the return value will be an empty vector.
//      *
//      * @return a vector of iteration counts.
//      */
//     virtual std::vector<size_t>
//     getNumberOfIterations() {return this->m_iterationCount;}


      /**
       * This method sets one of the termination criteria of the
       * optimization.  
       * 
       * @param maxIterations Each minimization will terminate after
       * this many iterations.
       */
      virtual void
      setMaxIterations(size_t maxIterations) {m_maxIterations = maxIterations;}

      
      /**
       * This method sets one of the termination criteria of the
       * optimization.  
       * 
       * @param maxLambda Iteration will terminate if LM parameter
       * "lambda" increases beyound this amount.
       */
      virtual void
      setMaxLambda(double maxLambda) {m_maxLambda = maxLambda;}

      
      /**
       * This method sets one of the termination criteria of the
       * optimization.  
       * 
       * @param minDrop Iteration will terminate if error fails to
       * decrease by at least this proportion for "strikes"
       * consecutive iterations.
       */
      virtual void
      setMinDrop(double minDrop) {m_minDrop = minDrop;}

    
      /**
       * This method sets one of the termination criteria of the
       * optimization.  Iteration will stop if the magnitude of the
       * gradient of the objective function at the current location is
       * less than the specified value.
       *
       * @param minimumGradientMagnitude The value at which the
       * magnitude of the objective function gradient will be
       * considered small enough to terminate iteration.
       */
      virtual void
      setMinimumGradientMagnitude(double minimumGradientMagnitude);

    
      /** 
       * This method sets minimization parameters.  Default values are
       * reasonable for functions which take values and arguments in the
       * "normal" range of 0 to 100 or so.
       * 
       * @param initialLambda This argument 
       * 
       * @param maxIterations Each minimization will terminate after
       * this many iterations.
       * 
       * @param maxLambda Iteration will terminate if LM parameter
       * "lambda" increases beyound this amount.
       * 
       * @param minLambda This argument sets a limit on how small LM
       * parameter "lambda" is allowed to get.
       * 
       * @param minError Iteration will terminate if the objective
       * function value goes below this level.
       * 
       * @param minimumGradientMagnitude The value at which the
       * magnitude of the objective function gradient will be
       * considered small enough to terminate iteration.
       * 
       * @param minDrop Iteration will terminate if error fails to
       * decrease by at least this proportion for "strikes"
       * consecutive iterations.
       * 
       * @param strikes Iteration will terminate if error fails to
       * decrease by at least the proportion specified by "minDrip"
       * for this many consecutive iterations.
       * 
       * @param maxBackSteps Iteration will terminate if the value of
       * LM parameter "lambda" must be increased this many times in a
       * row.
       * 
       * @param verbosity This argument indicates the desired output
       * level.  Setting verbosity to zero mean that no standard output
       * should be generated.  Higher numbers indicate increasingly more
       * output.
       */
      virtual void
      setParameters(double initialLambda = 1.0,
                    size_t maxIterations = 40,
                    double maxLambda = 1.0E7,
                    double minLambda = 1.0E-13,
                    double minError = 0.0,
                    double minimumGradientMagnitude = 1.0E-5,
                    double minDrop = 1.0E-4,
                    size_t strikes = 3,
                    int maxBackSteps = -1,
                    int verbosity = 0);
    
    
      /** 
       * This method sets the initial conditions for the minimization.
       * Gradient based search will start at this location in parameter
       * space.
       *
       * @param startPoint Indicates a point in the parameter space of
       * the objective function.
       */
      virtual void
      setStartPoint(const typename Functor::argument_type& startPoint);

    
      /**
       * This method sets the amount of text printed to the standard
       * output during the optimization.  Currently this method does
       * nothing, since OptimizerBFGS never generates any standard
       * output.
       *
       * @param verbosity This argument indicates the desired output
       * level.  Setting verbosity to zero mean that no standard output
       * should be generated.  Higher numbers indicate increasingly more
       * output.
       */
      virtual void
      setVerbosity(int verbosity) {m_verbosity = verbosity;}

    
      /**
       * The assignment operator deep copies its argument.
       * 
       * @param source The OptimizerLM instance to be copied.
       * 
       * @return Reference to *this.
       */
      virtual OptimizerLM&
      operator=(const OptimizerLM& source);
    
    protected:

      /** 
       * This protected member function is used to asses whether the
       * algorithm has reached convergence.
       * 
       * @param theta This argument specifies the parameter values
       * (arguments to the objective function) being assessed.
       * 
       * @param value This argument specifies the function value at the
       * point described by theta.
       * 
       * @param gradient This argument specifies the function gradient
       * at the point described by theta.
       * 
       * @return The return value gets progressively smaller as we
       * approach a local minimum.
       */
      double
      gradientConvergenceMetric(const argument_type& theta,
                                const result_type& value,
                                const argument_type& gradient);

    
      /** 
       * Perform the optimization.  This virtual function overrides the 
       * definition in Optimizer.
       * 
       * @return A std::pair of the vector parameter which brings the
       * specified Functor to an optimum, and the corresponding optimal
       * Functor value.
       */
      virtual
      std::pair<typename Functor::argument_type, typename Functor::result_type>
      run();


      inline virtual void
      verboseWrite(const char* message, int verbosity);


      template <class Type>
      inline void
      verboseWrite(const char* intro, const Type& subject, int verbosity);
    
      // Data members.
      double m_initialLambda;
      int m_maxBackSteps;
      size_t m_maxIterations;
      double m_maxLambda;
      double m_minDrop;
      double m_minError;
      double m_minGrad;
      double m_minLambda;
      argument_type m_startPoint;
      size_t m_strikes;
      int m_verbosity;

    }; // class OptimizerLM

  } // namespace optimization

} // namespace brick


/* ======= Declarations to maintain compatibility with legacy code. ======= */

namespace brick {

  using optimization::OptimizerLM;

} // namespace brick


/*******************************************************************
 * Member function definitions follow.  This would be a .cpp file
 * if it weren't templated.
 *******************************************************************/

#include <iostream>
#include <cmath>

#include <brick/numeric/utilities.hh>
#include <brick/optimization/optimizerCommon.hh>

namespace brick {

  namespace optimization {
  
    template <class Functor>
    OptimizerLM<Functor>::
    OptimizerLM()
      : Optimizer<Functor>(),
        m_initialLambda(),
        m_maxBackSteps(),
        m_maxIterations(),
        m_maxLambda(),
        m_minDrop(),
        m_minError(),
        m_minGrad(),
        m_minLambda(),
        m_startPoint(),
        m_strikes(),
        m_verbosity()
    {
      this->setParameters();
    }


    template <class Functor>
    OptimizerLM<Functor>::
    OptimizerLM(const Functor& functor)
      : Optimizer<Functor>(functor),
        m_initialLambda(),
        m_maxBackSteps(),
        m_maxIterations(),
        m_maxLambda(),
        m_minDrop(),
        m_minError(),
        m_minGrad(),
        m_minLambda(),
        m_startPoint(),
        m_strikes(),
        m_verbosity()
    {
      this->setParameters();
    }
  

    template<class Functor>
    OptimizerLM<Functor>::
    OptimizerLM(const OptimizerLM& source)
      : Optimizer<Functor>(source),
        m_initialLambda(source.m_initialLambda),
        m_maxBackSteps(source.m_maxBackSteps),
        m_maxIterations(source.m_maxIterations),
        m_maxLambda(source.m_maxLambda),
        m_minDrop(source.m_minDrop),
        m_minError(source.m_minError),
        m_minGrad(source.m_minGrad),
        m_minLambda(source.m_minLambda),
        m_startPoint(source.m_startPoint.size()),
        m_strikes(source.m_strikes),
        m_verbosity(source.m_verbosity)
    {
      copyArgumentType(source.m_startPoint, this->m_startPoint);
    }

  
    template <class Functor>
    OptimizerLM<Functor>::
    ~OptimizerLM()
    {
      // Empty
    }


    // This method sets one of the termination criteria of the
    // optimization.
    template<class Functor>
    void
    OptimizerLM<Functor>::
    setMinimumGradientMagnitude(double minimumGradientMagnitude)
    {
      this->m_minGrad = minimumGradientMagnitude * minimumGradientMagnitude;
    }

  
    template<class Functor>
    void
    OptimizerLM<Functor>::
    setParameters(double initialLambda,
                  size_t maxIterations,
                  double maxLambda,
                  double minLambda,
                  double minError,
                  double minimumGradientMagnitude,
                  double minDrop,
                  size_t strikes,
                  int maxBackSteps,
                  int verbosity)
    {
      // Copy input arguments.
      this->m_initialLambda = initialLambda;
      this->m_maxIterations = maxIterations;
      this->m_maxLambda = maxLambda;
      this->m_minLambda = minLambda;
      this->m_minError = minError;
      this->m_minGrad = minimumGradientMagnitude * minimumGradientMagnitude;
      this->m_minDrop = minDrop;
      this->m_strikes = strikes;
      this->m_maxBackSteps = maxBackSteps;
      this->m_verbosity = verbosity;

//     // Reset memory of previous minimization.
//     this->m_functionCallCount.clear();
//     this->m_gradientCallCount.clear();
//     this->m_iterationCount.clear();
    
      // We've changed the parameters, so we'll have to rerun the 
      // optimization.  Indicate this by setting the inherited member
      // m_needsOptimization.    
      Optimizer<Functor>::m_needsOptimization = true;
    }    


    template<class Functor>
    void
    OptimizerLM<Functor>::
    setStartPoint(const typename Functor::argument_type& startPoint)
    {
      copyArgumentType(startPoint, this->m_startPoint);

//     // Reset memory of previous minimization.
//     this->m_functionCallCount.clear();
//     this->m_gradientCallCount.clear();
//     this->m_iterationCount.clear();

      // We've changed the parameters, so we'll have to rerun the 
      // optimization.  Indicate this by setting the inherited member
      // m_needsOptimization.    
      Optimizer<Functor>::m_needsOptimization = true;
    }


    template<class Functor>
    OptimizerLM<Functor>&
    OptimizerLM<Functor>::
    operator=(const OptimizerLM<Functor>& source)
    {
      if(&source != this) {
        Optimizer<Functor>::operator=(source);
        m_initialLambda = source.m_initialLambda;
        m_maxBackSteps = source.m_maxBackSteps;
        m_maxIterations = source.m_maxIterations;
        m_maxLambda = source.m_maxLambda;
        m_minDrop = source.m_minDrop;
        m_minError = source.m_minError;
        m_minGrad = source.m_minGrad;
        m_minLambda = source.m_minLambda;
        copyArgumentType(source.m_startPoint, this->m_startPoint);
        m_strikes = source.m_strikes;
        m_verbosity = source.m_verbosity;
      }
      return *this;
    }


    // =============== Protected member functions below =============== //

    template <class Functor>
    std::pair<typename Functor::argument_type, typename Functor::result_type>
    OptimizerLM<Functor>::
    run()
    {
      // Check that we have a valid startPoint.
      if(this->m_startPoint.size() == 0) {
        BRICK_THROW(brick::common::StateException, 
		    "OptimizerLM<Functor>::run()",
		    "startPoint has not been initialized.");
      }
    
      // Initialize working location so that we start at the right place.
      argument_type theta(this->m_startPoint.size());
      copyArgumentType(this->m_startPoint, theta);

      // Initialize variables relating to convergence and convergence
      // failure.
      size_t strikes = 0;
      int backtrackCount = 0;

      // Initialize variables which will take function return values and
      // derivatives.
      result_type errorValue;
      argument_type dEdX(theta.size());
      brick::numeric::Array2D<brick::common::Float64> d2EdX2(theta.size(), theta.size());

      // Initialize intermediate values used by the minimization.
      brick::numeric::Array2D<brick::common::Float64> BMatrix(theta.size(), theta.size());
      brick::numeric::Array1D<brick::common::Float64> deltaX(theta.size());
      argument_type xCond(theta.size());
      brick::numeric::Array1D<result_type> errorHistory =
        brick::numeric::zeros<result_type>(m_maxIterations + 1);
      brick::common::Float64 lambda = m_initialLambda;

      // Get initial value of error function.
      errorValue = this->m_functor(theta);
      errorHistory[0] = errorValue;
      this->verboseWrite("Error History:\n", errorHistory, 1);

      // Loop until termination.
      size_t iterationIndex = 0;
      for(; iterationIndex < m_maxIterations; ++iterationIndex) {

        // Compute gradient.
        this->m_functor.computeGradientAndHessian(theta, dEdX, d2EdX2);

        // Gradient almost zero?
        if(dotArgumentType(dEdX,dEdX) <= m_minGrad) {
          this->verboseWrite("Tiny gradient, terminating iteration.\n", 1);
          break;
        }

        // Adjust lambda.
        while(lambda <= m_maxLambda) {

          // Precondition the Hessian matrix.
          std::copy(d2EdX2.begin(), d2EdX2.end(), BMatrix.begin());
          for(size_t diagIndex = 0; diagIndex < theta.size(); ++diagIndex) {
            BMatrix(diagIndex, diagIndex) += lambda;
          }

          // Solve B * deltaX = dEdX'
          copyArgumentType(dEdX, deltaX);
	  brick::linearAlgebra::linearSolveInPlace(BMatrix, deltaX);
        
          // We have a new candidate location in the error space.
          for(size_t elementIndex = 0; elementIndex < theta.size();
              ++elementIndex) {
            xCond[elementIndex] = theta[elementIndex] - deltaX[elementIndex];
          }

          // Do we have a decrease in the error function at the candidate
          // location?
          errorValue = m_functor(xCond);
          if(errorValue < errorHistory[iterationIndex]) {
            // Yes. Go on to the next iteration.
            backtrackCount = 0;
            errorHistory[iterationIndex + 1] = errorValue;
            copyArgumentType(xCond, theta);
            lambda /= 10.0;
            if(lambda < m_minLambda) {
              lambda = m_minLambda;
            }
            this->verboseWrite("Lambda = ", lambda, 1);
            this->verboseWrite("Theta = ", theta, 2);
            break;
          } else {
            // Error did not decrease.  Try a bigger lambda.
            ++backtrackCount;
            lambda *= 10.0;
            if(lambda > m_maxLambda) {
              break;
            }
            this->verboseWrite("Lambda = ", lambda, 1);

            // Make sure we haven't exceeded the maxBackSteps
            // termination criterion.
            if(m_maxBackSteps >= 0
               && backtrackCount > m_maxBackSteps) {
              break;
            }

          }
        }

        if(m_verbosity >= 1 ) {
          std::cout << "Error History:\n" << errorHistory << std::endl;
        }
// Note(xxx):  Not sure if we still want this functionality.
//       if(m_iterationFunctorPtr != 0) {
//         m_iterationFunctorPtr(theta);
//       }
    
        // Test termination conditions.
        brick::common::Float64 drop =
          (errorHistory[iterationIndex] - errorValue)
          / errorHistory[iterationIndex];
        if(drop < m_minDrop) {
          ++strikes;
          if(m_verbosity >= 2) {
            std::cout << "strikes = " << strikes << std::endl;
          }
        } else {
          strikes = 0;
        }

        if(lambda >= m_maxLambda 
           || strikes == m_strikes 
           || errorHistory[iterationIndex] <= m_minError 
           || (m_maxBackSteps >= 0 && backtrackCount >= m_maxBackSteps)) {
          if(m_verbosity >= 1) {
            std::cout << "Stopping with lambda = " << lambda
                      << " (" << m_maxLambda << ")\n"
                      << "              strikes = " << strikes
                      << " (" << m_strikes << ")\n"
                      << "              error = " << errorHistory[iterationIndex]
                      << " (" << m_minError << ")\n"
                      << "              backTrackCount = " << backtrackCount
                      << " (" << m_maxBackSteps << ")" << std::endl;
          }
          break;
        }
      }
      return std::make_pair(theta, errorHistory[iterationIndex]);
    }


    template <class Functor>
    inline void
    OptimizerLM<Functor>::
    verboseWrite(const char* message, int verbosity)
    {
      if(verbosity <= this->m_verbosity) {
        std::cout << message << std::flush;
      }
    }


    template <class Functor> template <class Type>
    inline void
    OptimizerLM<Functor>::
    verboseWrite(const char* intro, const Type& subject, int verbosity)
    {
      if(verbosity <= this->m_verbosity) {
        std::cout << intro << subject << std::endl;
      }
    }

  } // namespace optimization

} // namespace brick

#endif /* #ifndef BRICK_OPTIMIZATION_OPTIMIZERLM_HH */
