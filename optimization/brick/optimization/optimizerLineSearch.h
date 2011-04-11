/**
***************************************************************************
* @file brick/optimization/optimizerLineSearch.hh
*
* Header file declaring OptimizerLineSearch class.
*
* Copyright (C) 2003-2007 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_OPTIMIZATION_OPTIMIZERLINESEARCH_HH
#define BRICK_OPTIMIZATION_OPTIMIZERLINESEARCH_HH

#include <brick/optimization/optimizer.hh>>
#include <brick/optimization/optimizerCommon.hh>>

namespace brick {

  namespace optimization {

    /**
     ** OptimizerLineSearch implements the line search algorithm
     ** described in [1]. It takes as input an n dimensional point, the
     ** value of the function and gradient at that point, a direction,
     ** and an input quantity that limits the length of the steps. It
     ** uses these to find a new point within the length requested and
     ** along the direction specified at which the function has
     ** significantly decreased, and returns that n-dimensional point.
     **
     ** WARNING: The line search algorithm is not really an
     ** optimization.  If you want to find the minimum function value
     ** along a line in parameter space, this is _not_ the class you're
     ** looking for.
     **
     ** [1] W. H. Press et al., Numerical Recipes in C The Art of
     ** Scientific Computing, Cambridge University Press, 1988.
     **/
    template <class Functor>
    class OptimizerLineSearch
      : public Optimizer<Functor>
    {
    public:
      // Typedefs for convenience
      typedef typename Functor::argument_type argument_type;
      typedef typename Functor::result_type result_type;

      /**
       * Default constructor sets parameters to reasonable values for
       * functions which take values and arguments in the "normal" range
       * of 0 to 100 or so.
       */
      OptimizerLineSearch();

      /**
       * Constructor which specifies the specific Functor instance to
       * use.  Using this constructor exclusively avoids the danger of
       * calling optimalValue() or optimum() before a Functor instance
       * has been specified.
       *
       * @param functor A copy of this argument will be stored
       * internally for use in optimization.
       */
      explicit OptimizerLineSearch(const Functor& functor);

      /**
       * Overloading of constructor 
       *
       * @param functor A copy of this argument will be stored
       * internally for use in optimization.
       * 
       * @param startPoint The initial search point on the function.
       * 
       * @param startGradient The initial gradient at startPoint.
       *
       * @param startValue The function value at startPoint.     
       */
      OptimizerLineSearch(const Functor& functor, 
                          const argument_type& startPoint,
                          const argument_type& startGradient,
                          const result_type& startValue);

      /** 
       * Copy constructor.
       * 
       * @param source The OptimizerLineSearch instance to be copied.
       */
      OptimizerLineSearch(const OptimizerLineSearch& source);

      /**
       * Destructor.
       */
      virtual
      ~OptimizerLineSearch();

      /** 
       * If a valid minimization result is available, this method returns
       * the number of function calls required to produce that result.  If
       * no valid minimization result is available, the return value is 0.
       * 
       * @return The number of function calls spent in the last
       * minimization, or 0.
       */
      size_t
      getNumberOfFunctionCalls() {return this->m_functionCallCount;}
    
      /** 
       * Sets part of the initial conditions for the minimization.
       * Search will start at this location in parameter space.  The
       * member function setInitialStep() must also be called before
       * the line search can be run.
       *
       * @param startPoint Indicates a point in the parameter space of
       * the objective function.
       */
      virtual void
      setStartPoint(const argument_type& startPoint);

      /**
       * Sets part of the initial conditions for the minimization.
       * Search will start at this location in parameter space.  Use
       * this function to avoid recomputing the function value and
       * gradient at this point.  The member function setInitialStep()
       * must also be called before the line search can be run.
       *
       * @param startPoint Indicates a point in the parameter space of
       * the objective function.
       * @param startValue The value of the objective function, when
       * evaluated at startPoint.
       * @param startGradient The gradient of the objective function,
       * when evaluated at startPoint.
       */
      virtual void
      setStartPoint(const argument_type& startPoint,
                    const result_type& startValue,
                    const argument_type& startGradient);

      /** 
       * Sets part of the initial conditions for the minimization.  This
       * specifies both the initial search step and the direction of
       * search in parameter space.  The member function setStartPoint()
       * must also be called before the line search can be run.
       *
       * @param initialStep Indicates the direction in space where
       * the objective function will be going.
       */
      virtual void
      setInitialStep(const argument_type& initialStep);

      /** 
       * Sets the line search parameters.  Default values are reasonable
       * for functions which take values and arguments in the "normal"
       * range of 0 to 100 or so.
       *
       * @param argumentTolerance Iteration will terminate when a
       * minimization step moves, along every axis, a distance less than
       * a threshold which is linearly related to this factor.
       * @param alpha Iteration will terminate when the function value
       * has decreased by an amount which is linearly related to this
       * factor.
       * @param maximumStepMagnitude Sets the maximum size of the initial
       * step of the line search.
       */
      void
      setParameters(double argumentTolerance=1.0e-7,
                    double alpha=1.0e-4,
                    double maximumStepMagnitude=100.0);
    
      /**
       * Assignment operator.
       * 
       * @param source The OptimizerLineSearch instance to be copied.
       * @return Reference to *this.
       */
      OptimizerLineSearch&
      operator=(const OptimizerLineSearch& source);
    
    protected:

      /** 
       * This protected member function verifies that all necessary
       * values have been set prior to running the optimization.  If
       * anything is missing or inconsistent, it throws an exception.
       */
      void
      checkState() {
        // Warning(xxx): currently no state checking is implemented.
      }
    
      /** 
       * Perform the minimization.  This overrides Optimizer<Functor>::run().
       * 
       * @return A std::pair of the vector parameter which brings the
       * specified Functor to a point of sufficient decrease, and the
       * corresponding Functor value.
       */
      std::pair<typename Functor::argument_type, typename Functor::result_type>
      run();
    
      // Data members
      double m_alpha;
      double m_argumentTolerance;
      size_t m_functionCallCount;
      argument_type m_initialStep;
      double m_initialStepMagnitude;
      double m_maximumStepMagnitude;
      argument_type m_startGradient;
      argument_type m_startPoint;
      result_type m_startValue;

    }; // class OptimizerLineSearch

  } // namespace optimization

} // namespace brick


/* ======= Declarations to maintain compatibility with legacy code. ======= */

namespace brick {

  using optimization::OptimizerLineSearch;

} // namespace brick


/*******************************************************************
 * Member function definitions follow.  This would be a .cpp file
 * if it weren't templated.
 *******************************************************************/

#include <cmath>
#include <sstream>
#include <brick/common/exception.hh>>
#include <brick/optimization/optimizerCommon.hh>>

namespace brick {

  namespace optimization {

    template <class Functor>
    OptimizerLineSearch<Functor>::
    OptimizerLineSearch()
      : Optimizer<Functor>(),
        m_alpha(0.0),
        m_argumentTolerance(0.0),
        m_functionCallCount(0),
        m_initialStep(),
        m_initialStepMagnitude(0.0),
        m_maximumStepMagnitude(0.0),
        m_startGradient(),
        m_startPoint(),
        m_startValue(0.0)
    {
      this->setParameters();
    }

    template <class Functor>
    OptimizerLineSearch<Functor>::
    OptimizerLineSearch(const Functor& functor)
      : Optimizer<Functor>(functor),
        m_alpha(0.0),
        m_argumentTolerance(0.0),
        m_functionCallCount(0),
        m_initialStep(),
        m_initialStepMagnitude(0.0),
        m_maximumStepMagnitude(0.0),
        m_startGradient(),
        m_startPoint(),
        m_startValue(0.0)
    {
      this->setParameters();
    }
  
    template<class Functor>
    OptimizerLineSearch<Functor>::
    OptimizerLineSearch(const OptimizerLineSearch& source)
      : Optimizer<Functor>(source),
        m_alpha(source.m_alpha),
        m_argumentTolerance(source.argumentTolerance),
        m_functionCallCount(source.m_functionCallCount),
        m_initialStep(source.m_initialStep.size()),
        m_initialStepMagnitude(source.initialStepMagnitude),
        m_maximumStepMagnitude(source.m_maximumStepMagnitude),
        m_startGradient(source.m_startGradient.size()),
        m_startPoint(source.m_startPoint.size()),
        m_startValue(source.m_startValue)
    {
      copyArgumentType(source.m_initialStep, this->m_initialStep);
      copyArgumentType(source.m_startGradient, this->m_startGradient);
      copyArgumentType(source.m_startPoint, this->m_startPoint);
    }

    template <class Functor>
    OptimizerLineSearch<Functor>::
    ~OptimizerLineSearch()
    {
      // Empty
    }

    template <class Functor>
    void OptimizerLineSearch<Functor>::
    setStartPoint(const argument_type& startPoint)
    {
      // Copy input argument.
      copyArgumentType(startPoint, this->m_startPoint);

      // Compute other starting info.  Note that we don't need to use
      // copyArgumentType(), because we have sole ownership of the
      // return value from m_functor.gradient().
      this->m_startValue = this->m_functor(startPoint);
      this->m_startGradient = this->m_functor.gradient(startPoint);
    
      // Reset the count of function calls.  We ignore the calls used 
      // to initialize m_startValue, etc. above.
      this->m_functionCallCount = 0;

      // We've changed the starting point, so we'll have to rerun the 
      // optimization.  Indicate this by setting the inherited member
      // m_needsOptimization.
      this->m_needsOptimization = true;
    }

    template<class Functor>
    void
    OptimizerLineSearch<Functor>::
    setStartPoint(const argument_type& startPoint,
                  const result_type& startValue,
                  const argument_type& startGradient)
    {
      // Copy input arguments.
      copyArgumentType(startPoint, this->m_startPoint);
      this->m_startValue = startValue;
      copyArgumentType(startGradient, this->m_startGradient);
    
      // Reset the count of function calls.
      this->m_functionCallCount = 0;

      // We've changed the starting point, so we'll have to rerun the 
      // optimization.  Indicate this by setting the inherited member
      // m_needsOptimization.
      this->m_needsOptimization = true;
    }

    template<class Functor>
    void OptimizerLineSearch<Functor>::
    setInitialStep(const argument_type& initialStep)
    {
      // First, make sure this is a valid initialStep.  
      double initialStepMagnitude =
        std::sqrt(dotArgumentType(initialStep, initialStep));
      if(initialStepMagnitude == 0.0) {
        BRICK_THROW3(ValueException, "OptimizerLineSearch::setInitialStep()",
                   "initialStep magnitude is equal to 0.0.");
      }

      // OK, valid initialStep, save it for later.  Since we'll also
      // be needing the magnitude of initialStep, save it too.
      this->m_initialStepMagnitude = initialStepMagnitude;
      copyArgumentType(initialStep, this->m_initialStep);
    
      // Reset the count of function calls.
      this->m_functionCallCount = 0;

      // We've changed the line direction, so we'll have to rerun the 
      // optimization.  Indicate this by setting the inherited member
      // m_needsOptimization.
      this->m_needsOptimization = true;
    }
  
    template<class Functor>
    void OptimizerLineSearch<Functor>::
    setParameters(double argumentTolerance,
                  double alpha,
                  double maximumStepMagnitude)
    {
      this->m_argumentTolerance = argumentTolerance;
      this->m_alpha = alpha;
      this->m_maximumStepMagnitude = maximumStepMagnitude;

      // Reset the count of function calls.
      this->m_functionCallCount = 0;

      // We've changed the parameters, so we'll have to rerun the 
      // optimization.  Indicate this by setting the inherited member
      // m_needsOptimization.    
      this->m_needsOptimization = true;
    }    
  
    template <class Functor>
    OptimizerLineSearch<Functor>& OptimizerLineSearch<Functor>::
    operator=(const OptimizerLineSearch& source)
    {
      // First, call the assignment operator of the parent class.
      Optimizer<Functor>::operator=(source);

      // Then copy all of the variables which are specific to
      // OptimizerLineSearch.
      this->m_alpha = source.m_alpha;
      this->m_argumentTolerance = source.m_argumentTolerance;
      this->m_functionCallCount = source.m_functionCallCount;
      this->m_initialStepMagnitude = source.m_initialStepMagnitude;
      this->m_maximumStepMagnitude = source.m_maximumStepMagnitude;
      this->m_startValue = source.m_startValue;
      copyArgumentType(source.m_initialStep, this->m_initialStep);
      copyArgumentType(source.m_startGradient, this->m_startGradient);
      copyArgumentType(source.m_startPoint, this->m_startPoint);
      return *this;
    }

    template <class Functor>
    std::pair<typename Functor::argument_type, typename Functor::result_type>
    OptimizerLineSearch<Functor>::
    run()
    {
      // Make sure all required starting values have been specified.
      this->checkState();

      // Reset the count of function calls.
      this->m_functionCallCount = 0;
    
      // Initialize a few variables.
      argument_type initialStep;
      copyArgumentType(this->m_initialStep, initialStep);
      double initialStepMagnitude = this->m_initialStepMagnitude;

      // Restrict length of initial step.
      if(this->m_initialStepMagnitude > this->m_maximumStepMagnitude) {
        initialStepMagnitude = this->m_maximumStepMagnitude;
        for(size_t index = 0; index < initialStep.size(); ++index) {
          initialStep[index] *=
            this->m_maximumStepMagnitude / this->m_initialStepMagnitude;
        }
      }

      // Compute degree of similarity between startGradient and initialStep.
      double slope = dotArgumentType(this->m_startGradient, initialStep);
      if(slope >= 0.0) {
        std::ostringstream message;
        message << "Initial search direction: " << initialStep << " "
                << "is not downhill with respect to initial gradient: "
                << this->m_startGradient << ".";
        BRICK_THROW3(StateException, "OptimizerLineSearch::run()",
                   message.str().c_str());
                   
      }

      // Compute smallest allowable step, allowing for numerical issues.
      double scale = contextSensitiveScale(initialStep, this->m_startPoint);
      if(scale == 0.0) {
        BRICK_THROW3(RunTimeException, "OptimizerLineSearch::run()",
                   "Invalid initial scale.  "
                   "Perhaps initialStep is the zero vector.");
      }
      double minimumLambda = this->m_argumentTolerance / scale;

      // Now do the search.
      double lambdaValue = 1.0;
      argument_type nextSamplePoint(initialStep.size());
      double previousLambdaValue = 1.0;
      double previousNextFunctionValue = this->m_startValue;
      while(1) {
        // Check first termination criterion.
        // Note: this differs from NRC in that NRC first computes
        // nextFunctionValue before doing this check, and consequently
        // returns it to the calling context if the following conditional
        // passes.
        if(lambdaValue < minimumLambda) {
          // In this case, root finding programs should double check
          // convergence!
          copyArgumentType(this->m_startPoint, nextSamplePoint);
          // Note: This return value differs from NRC, in that we return
          // startValue instead of the most recently computed value.
          // return std::make_pair(nextSamplePoint, nextFunctionValue);
          return std::make_pair(nextSamplePoint, this->m_startValue);
        }

        // Compute next place to evaluate function.
        for(size_t index = 0; index < initialStep.size(); ++index) {
          nextSamplePoint[index] =
            this->m_startPoint[index] + (lambdaValue * initialStep[index]);
        }

        // Do the function evaluation.
        result_type nextFunctionValue = this->m_functor(nextSamplePoint);
        ++(this->m_functionCallCount);

        // Check second termination criterion.
        if(nextFunctionValue
           < this->m_startValue + (this->m_alpha * lambdaValue * slope)) {
          return std::make_pair(nextSamplePoint, nextFunctionValue);
        }

        // Compute next value of lambda.
        double lambdaHat;
        if(lambdaValue == 1.0) {
          lambdaHat =
            (-slope / (2.0 * (nextFunctionValue - this->m_startValue
                              - slope)));
        } else {
          // Buld vector.
          double b0 = (nextFunctionValue - this->m_startValue
                       - (lambdaValue * slope));
          double b1 = (previousNextFunctionValue - this->m_startValue
                       - (previousLambdaValue * slope));
          // Build matrix.
          double deltaLambda = lambdaValue - previousLambdaValue;
          double a00 = 1.0 / (lambdaValue * lambdaValue * deltaLambda);
          double a01 = -1.0 / (previousLambdaValue * previousLambdaValue
                               * deltaLambda);
          double a10 = -previousLambdaValue / (lambdaValue * lambdaValue
                                               * deltaLambda);
          double a11 = lambdaValue / (previousLambdaValue * previousLambdaValue
                                      * deltaLambda);
          // Do matrix multiplication.
          double ab0 = a00 * b0 + a01 * b1;
          double ab1 = a10 * b0 + a11 * b1;

          // Use result to choose next value of lambda.
          if(ab0 == 0.0) {
            if(ab1 == 0.0) {
              BRICK_THROW3(RunTimeException, "OptimizerLineSearch::run()",
                         "Invalid value for internal variable ab1.");
            }
            lambdaHat = -slope / (2.0 * ab1);
          } else {
            double discriminant = ab1 * ab1 - (3.0 * ab0 * slope);
            if(discriminant < 0.0) {
              BRICK_THROW3(RunTimeException, "OptimizerLineSearch::run()",
                         "Roundoff error, discriminant < 0.0");
            }
            lambdaHat = (-ab1 + std::sqrt(discriminant)) / (3.0 * ab0);
          }
          if(lambdaHat > (0.5 * lambdaValue)) {
            lambdaHat = 0.5 * lambdaValue;
          }
        }
        previousLambdaValue = lambdaValue;
        previousNextFunctionValue = nextFunctionValue;
        lambdaValue = std::max(lambdaHat, 0.1 * lambdaValue);
      }
    }

  } // namespace optimization

} // namespace brick

#endif /* #ifndef BRICK_OPTIMIZATION_OPTIMIZERLINESEARCH_HH */
