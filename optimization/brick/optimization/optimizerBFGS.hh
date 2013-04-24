/**
***************************************************************************
* @file brick/optimization/optimizerBFGS.hh
*
* Header file declaring OptimizerBFGS class.
*
* (C) Copyright 2003-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying LICENSE file for details.
*
***************************************************************************
**/

#ifndef BRICK_OPTIMIZATION_OPTIMIZERBFGS_HH
#define BRICK_OPTIMIZATION_OPTIMIZERBFGS_HH

#include <limits>
#include <vector>
#include <brick/common/triple.hh>
#include <brick/optimization/optimizer.hh>
#include <brick/optimization/optimizerLineSearch.hh>

namespace brick {

  namespace optimization {

    /**
     ** OptimizerBFGS implements the Quasi-Newton method of Broyden,
     ** Fletcher, Goldfarb, and Shanno, as described in [1] (and
     ** possibly in [2]).  This algorithm seeks the parameter value
     ** which minimizes the objective function.  The template parameter
     ** (Functor) defines the type to use as the objective function of
     ** the minimization, and must support the GradientFunction
     ** interface.
     **
     ** [1] W. H. Press et al., Numerical Recipes in C The Art of
     ** Scientific Computing, Cambridge University Press, 1988.
     **
     ** [2] C. G. Broyden, R. Fletcher, M. Goldfarb and Shanno, The
     ** convergence of a class of double rank minimization algorithms,
     ** J. Inst. Math. Appl., 6:222-231, 1970.
     **/
    template <class Functor>
    class OptimizerBFGS
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
      OptimizerBFGS();

      /**
       * This constructor specifies the specific Functor instance to
       * use.  Using this constructor exclusively avoids the danger of
       * calling optimalValue() or optimum() before a Functor instance
       * has been specified.
       *
       * @param functor A copy of this argument will be stored
       * internally for use in optimization.
       */
      explicit OptimizerBFGS(const Functor& functor);

      /** 
       * Copy constructor.  This constructor simply copies the source
       * argument.
       * 
       * @param source The OptimizerBFGS instance to be copied.
       */
      OptimizerBFGS(const OptimizerBFGS& source);

      /**
       * The destructor destroys the class instance and deallocates any
       * associated storage.
       */
      virtual
      ~OptimizerBFGS();

      /** 
       * This method returns the number of function calls required to
       * complete the previous minimization.  If the minimization
       * parameter "restarts" is 0, there will be only one element in
       * the returned vector.  If restarts is greater than 0, the first
       * element of the return value will reflect the number of function
       * calls in the initial minimization, and subsequent numbers will
       * reflect the number of function calls in the following restarted
       * minimizations.  If a valid minimization has not been performed
       * since the last update to startPoint, parameters, etc., then the
       * return value will be an empty vector.
       *
       * @return a vector of function call counts.
       */
      virtual std::vector<size_t>
      getNumberOfFunctionCalls() {return this->m_functionCallCount;}

      /** 
       * This method returns the number of gradient calls required to
       * complete the previous minimization.  If the minimization
       * parameter "restarts" is 0, there will be only one element in
       * the returned vector.  If restarts is greater than 0, the first
       * element of the return value will reflect the number of gradient
       * calls in the initial minimization, and subsequent numbers will
       * reflect the number of gradient calls in the following restarted
       * minimizations.  If a valid minimization has not been performed
       * since the last update to startPoint, parameters, etc., then the
       * return value will be an empty vector.
       *
       * @return a vector of gradient call counts.
       */
      virtual std::vector<size_t>
      getNumberOfGradientCalls() {return this->m_gradientCallCount;}
    
      /** 
       * This method returns the number of iterations required to
       * complete the previous minimization.  If the minimization
       * parameter "restarts" is 0, there will be only one element in
       * the returned vector.  If restarts is greater than 0, the first
       * element of the return value will reflect the number of
       * iterations in the initial minimization, and subsequent numbers
       * will reflect the number of iterations in the following
       * restarted minimizations.  If a valid minimization has not been
       * performed since the last update to startPoint, parameters,
       * etc., then the return value will be an empty vector.
       *
       * @return a vector of iteration counts.
       */
      virtual std::vector<size_t>
      getNumberOfIterations() {return this->m_iterationCount;}
                  
      /** 
       * This method sets minimization parameters.  Default values are
       * reasonable for functions which take values and arguments in the
       * "normal" range of 0 to 100 or so.
       *
       * @param iterationLimit Each minimization will terminate after
       * this many iterations.
       *
       * @param numberOfRestarts Following successful termination, the
       * minimization will be re-run this many times to refine the
       * result accuracy.  Generally you should leave this at its
       * default value.
       *
       * @param gradientTolerance Iteration will terminate when the
       * magnitude of the gradient times the magnitude of the parameter
       * vector becomes smaller than the function value by this factor.
       *
       * @param argumentTolerance Iteration will terminate when a
       * minimization step moves, along every axis, a distance less than
       * this factor times the corresponding element of the argument
       * vector.
       *
       * @param numericEpsilon Sets the internal epsilon value of the
       * gradient update routine.
       *
       * @param maximumStepMagnitudeFactor Sets the maximum step
       * distance for each line minimization in the algorithm.
       *
       * @param minimumFunctionValue Iteration will terminate if the
       * objective function value falls to or below this value.
       */
      virtual void
      setParameters(size_t iterationLimit = 500,
                    size_t numberOfRestarts = 1,
                    double argumentTolerance = 1.2E-7, // Generally 4*(num. Eps.)
                    double gradientTolerance = 0.00001,
                    double lineSearchAlpha = 1.0E-4,
                    double lineSearchArgumentTolerance = 1.0e-7,
                    double numericEpsilon = 3.0E-8,
                    double maximumStepMagnitudeFactor = 100.0,
                    double minimumFunctionValue =
                      -std::numeric_limits<double>::max());
    
    
      /** 
       * This method sets the optimization parameter controlling the
       * maximum number of iterations, without affecting any other
       * optimization parameters.
       * 
       * @param iterationLimit Each minimization will terminate after
       * this many iterations.
       */
      virtual void
      setIterationLimit(size_t iterationLimit) {
        this->m_iterationLimit = iterationLimit;
      }

    
      /** 
       * This method sets the optimization parameter controlling the
       * number of restarts, without affecting any other optimization
       * parameters.
       * 
       * @param numberOfRestarts Following successful termination, the
       * minimization will be re-run this many times to refine the
       * result accuracy.
       */
      virtual void
      setNumberOfRestarts(size_t numberOfRestarts) {
        this->m_numberOfRestarts = numberOfRestarts;
      }

      
      /** 
       * This method sets the optimization parameter controlling the
       * function value at which the optimization will be considered
       * "close enough." 
       * 
       * @param minimumFunctionValue Iteration will terminate if the
       * objective function value falls to or below this value.
       */
      virtual void
      setMinimumFunctionValue(double minimumFunctionValue) {
        this->m_minimumFunctionValue = minimumFunctionValue;
      }

      
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
      setVerbosity(int verbosity) {this->m_verbosity = verbosity;}

    
      /**
       * Assignment operator.
       * 
       * @param source The OptimizerBFGS instance to be copied.
       * 
       * @return Reference to *this.
       */
      virtual OptimizerBFGS&
      operator=(const OptimizerBFGS& source);
    
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

      /** 
       * Perform one complete BFGS minimization, starting from the
       * specified position.
       * 
       * @param theta This argument specifies the point at which to
       * start the minimization.
       * 
       * @param startValue This argument should be set to the value of
       * the objective function evaluated at theta.
       * 
       * @param startGradient This argument should be set to the objective
       * function gradient evaluated at theta.
       * 
       * @param numberOfFunctionCalls This parameter is used to return
       * the number of function calls required to perform the
       * minimization.
       * 
       * @param numberOfGradientCalls This parameter is used to return
       * the number of gradient calls required to perform the
       * minimization.
       * 
       * @param numberOfIterations This parameter is used to return the
       * number of BFGS iterations required to perform the minimization.
       * 
       * @return A brick::triple of the vector parameter which brings the
       * specified Functor to an optimum, and the corresponding optimal
       * Functor value, and the corresponding gradient.
       */
      brick::common::Triple<typename Functor::argument_type,
                            typename Functor::result_type,
                            typename Functor::argument_type>
      doBfgs(const argument_type& theta,
             const result_type& startValue,
             const argument_type& startGradient,
             size_t& numberOfFunctionCalls,
             size_t& numberOfGradientCalls,
             size_t& numberOfIterations);

      // Data members.
      double m_argumentTolerance;
      size_t m_iterationLimit;
      double m_gradientTolerance;
      double m_lineSearchAlpha;
      double m_lineSearchArgumentTolerance;
      double m_maximumStepMagnitudeFactor;
      double m_minimumFunctionValue;
      size_t m_numberOfRestarts;
      double m_numericEpsilon;
      OptimizerLineSearch<Functor> m_optimizerLineSearch;
      argument_type m_startPoint;
      int m_verbosity;
      
      // Data members used for bookkeeping.
      std::vector<size_t> m_functionCallCount;
      std::vector<size_t> m_gradientCallCount;
      std::vector<size_t> m_iterationCount;
    
    }; // class OptimizerBFGS

  } // namespace optimization

} // namespace brick


/* ======= Declarations to maintain compatibility with legacy code. ======= */

namespace brick {

  using optimization::OptimizerBFGS;

} // namespace brick


/*******************************************************************
 * Member function definitions follow.  This would be a .cpp file
 * if it weren't templated.
 *******************************************************************/

#include <cmath>
#include <iomanip>
#include <iostream>
#include <brick/numeric/array1D.hh>
#include <brick/numeric/array2D.hh>
#include <brick/numeric/utilities.hh>
#include <brick/optimization/optimizerCommon.hh>

namespace brick {

  namespace optimization {
  
    template <class Functor>
    OptimizerBFGS<Functor>::
    OptimizerBFGS()
      : Optimizer<Functor>(),
        m_argumentTolerance(),
        m_iterationLimit(),
        m_gradientTolerance(),
        m_lineSearchAlpha(),
        m_lineSearchArgumentTolerance(),
        m_maximumStepMagnitudeFactor(),
        m_minimumFunctionValue(),
        m_numberOfRestarts(),
        m_numericEpsilon(),
        m_optimizerLineSearch(),
        m_startPoint(),
        m_verbosity(0),
        m_functionCallCount(),
        m_gradientCallCount(),
        m_iterationCount()
    {
      this->setParameters();
    }


    template <class Functor>
    OptimizerBFGS<Functor>::
    OptimizerBFGS(const Functor& functor)
      : Optimizer<Functor>(functor),
        m_argumentTolerance(),
        m_iterationLimit(),
        m_gradientTolerance(),
        m_lineSearchAlpha(),
        m_lineSearchArgumentTolerance(),
        m_maximumStepMagnitudeFactor(),
        m_minimumFunctionValue(),
        m_numberOfRestarts(),
        m_numericEpsilon(),
        m_optimizerLineSearch(functor),
        m_startPoint(),
        m_verbosity(0),
        m_functionCallCount(),
        m_gradientCallCount(),
        m_iterationCount()    
    {
      this->setParameters();
    }
  

    template<class Functor>
    OptimizerBFGS<Functor>::
    OptimizerBFGS(const OptimizerBFGS& source)
      : Optimizer<Functor>(source),
        m_argumentTolerance(source.m_argumentTolerance),
        m_iterationLimit(source.m_iterationLimit),
        m_gradientTolerance(source.m_gradientTolerance),
        m_lineSearchAlpha(source.m_lineSearchAlpha),
        m_lineSearchArgumentTolerance(source.m_lineSearchArgumentTolerance),
        m_maximumStepMagnitudeFactor(source.m_maximumStepMagnitudeFactor),
        m_minimumFunctionValue(source.m_minimumFunctionValue),
        m_numberOfRestarts(source.m_numberOfRestarts),
        m_numericEpsilon(source.m_numericEpsilon),
        m_optimizerLineSearch(source.m_optimizerLineSearch),
        m_startPoint(source.m_startPoint.size()),
        m_verbosity(source.m_verbosity),
        m_functionCallCount(source.m_functionCallCount),
        m_gradientCallCount(source.m_gradientCallCount),
        m_iterationCount(source.m_iterationCount)
    {
      copyArgumentType(source.m_startPoint, this->m_startPoint);
    }

    template <class Functor>
    OptimizerBFGS<Functor>::
    ~OptimizerBFGS()
    {
      // Empty
    }

    template<class Functor>
    void
    OptimizerBFGS<Functor>::
    setParameters(size_t iterationLimit,
                  size_t numberOfRestarts,
                  double argumentTolerance,
                  double gradientTolerance,
                  double lineSearchAlpha,
                  double lineSearchArgumentTolerance,
                  double numericEpsilon,
                  double maximumStepMagnitudeFactor,
                  double minimumFunctionValue)
    {
      // Copy input arguments.
      this->m_iterationLimit = iterationLimit;
      this->m_numberOfRestarts = numberOfRestarts;
      this->m_argumentTolerance = argumentTolerance;
      this->m_gradientTolerance = gradientTolerance;
      this->m_lineSearchAlpha = lineSearchAlpha;
      this->m_lineSearchArgumentTolerance = lineSearchArgumentTolerance;
      this->m_numericEpsilon = numericEpsilon;
      this->m_maximumStepMagnitudeFactor = maximumStepMagnitudeFactor;
      this->m_minimumFunctionValue = minimumFunctionValue;

      // Reset memory of previous minimization.
      this->m_functionCallCount.clear();
      this->m_gradientCallCount.clear();
      this->m_iterationCount.clear();
    
      // We've changed the parameters, so we'll have to rerun the 
      // optimization.  Indicate this by setting the inherited member
      // m_needsOptimization.    
      Optimizer<Functor>::m_needsOptimization = true;
    }    

    
    template<class Functor>
    void
    OptimizerBFGS<Functor>::
    setStartPoint(const typename Functor::argument_type& startPoint)
    {
      copyArgumentType(startPoint, this->m_startPoint);

      // Reset memory of previous minimization.
      this->m_functionCallCount.clear();
      this->m_gradientCallCount.clear();
      this->m_iterationCount.clear();

      // We've changed the parameters, so we'll have to rerun the 
      // optimization.  Indicate this by setting the inherited member
      // m_needsOptimization.    
      Optimizer<Functor>::m_needsOptimization = true;
    }


    template<class Functor>
    OptimizerBFGS<Functor>&
    OptimizerBFGS<Functor>::
    operator=(const OptimizerBFGS<Functor>& source)
    {
      Optimizer<Functor>::operator=(source);
      this->m_argumentTolerance = source.m_argumentTolerance;
      this->m_iterationLimit = source.m_iterationLimit;
      this->m_gradientTolerance = source.m_gradientTolerance;
      this->m_lineSearchAlpha = source.m_lineSearchAlpha;
      this->m_lineSearchArgumentTolerance = source.m_lineSearchArgumentTolerance;
      this->m_maximumStepMagnitudeFactor = source.m_maximumStepMagnitudeFactor;
      this->m_minimumFunctionValue = source.m_minimumFunctionValue;
      this->m_numberOfRestarts = source.m_numberOfRestarts;
      this->m_numericEpsilon = source.m_numericEpsilon;
      this->m_optimizerLineSearch = source.m_optimizerLineSearch;
      copyArgumentType(source.m_startPoint, this->m_startPoint);

      this->m_functionCallCount = source.m_functionCallCount;
      this->m_gradientCallCount = source.m_gradientCallCount;
      this->m_iterationCount = source.m_iterationCount;
    
      return *this;
    }


    // =============== Protected member functions below =============== //

    template <class Functor>
    double
    OptimizerBFGS<Functor>::
    gradientConvergenceMetric(const argument_type& theta,
                              const result_type& value,
                              const argument_type& gradient)
    {
      double returnValue = 0.0;
      double denominator = std::max(static_cast<double>(value), 1.0);
      for(size_t index = 0; index < theta.size(); ++index) {
        double thetaAbsValue = std::fabs(theta[index]);
        double gradientAbsValue = std::fabs(gradient[index]);
        double candidate =
          gradientAbsValue * std::max(thetaAbsValue, 1.0) / denominator;
        if(candidate > returnValue) {
          returnValue = candidate;
        }
      }
      return returnValue;
    }

    template <class Functor>
    std::pair<typename Functor::argument_type, typename Functor::result_type>
    OptimizerBFGS<Functor>::
    run()
    {
      // Check that we have a valid startPoint.
      if(this->m_startPoint.size() == 0) {
        BRICK_THROW(brick::common::StateException,
		    "OptimizerBFGS<Functor>::run()",
		    "startPoint has not been initialized.");
      }
    
      // Initialize working location so that we start at the right place.
      argument_type theta(this->m_startPoint.size());
      copyArgumentType(this->m_startPoint, theta);

      // Compute initial values of function and its gradient.
      result_type startValue = this->m_functor(theta);
      argument_type startGradient = this->m_functor.gradient(theta);

      // Now run the optimization.
      this->m_functionCallCount.clear();
      this->m_gradientCallCount.clear();
      this->m_iterationCount.clear();
      size_t functionCallCount;
      size_t gradientCallCount;
      size_t iterationCount;
      brick::common::Triple<argument_type, result_type, argument_type>
        optimum_optimalValue_gradient =
        this->doBfgs(theta, startValue, startGradient, functionCallCount,
                     gradientCallCount, iterationCount);
      // When accounting function calls, don't forget to add the initial
      // function and gradient evaluations above.
      this->m_functionCallCount.push_back(functionCallCount + 1);
      this->m_gradientCallCount.push_back(gradientCallCount + 1);
      this->m_iterationCount.push_back(iterationCount);

      // Restart as many times as requested.
      for(size_t index = 0; index < this->m_numberOfRestarts; ++index) {
        copyArgumentType(optimum_optimalValue_gradient.first, theta);
        startValue = optimum_optimalValue_gradient.second;
        copyArgumentType(optimum_optimalValue_gradient.third,
                         startGradient);
        optimum_optimalValue_gradient =
          this->doBfgs(theta, startValue, startGradient, functionCallCount,
                       gradientCallCount, iterationCount);
        this->m_functionCallCount.push_back(functionCallCount);
        this->m_gradientCallCount.push_back(gradientCallCount);
        this->m_iterationCount.push_back(iterationCount);
      }

      return std::make_pair(optimum_optimalValue_gradient.first,
                            optimum_optimalValue_gradient.second);
    }
    
    template <class Functor>
    brick::common::Triple<typename Functor::argument_type,
                          typename Functor::result_type,
                          typename Functor::argument_type>
    OptimizerBFGS<Functor>::
    doBfgs(const argument_type& theta,
           const result_type& startValue,
           const argument_type& startGradient,
           size_t& numberOfFunctionCalls,
           size_t& numberOfGradientCalls,
           size_t& numberOfIterations)
    {
      // Basic initializations.
      size_t dimensionality = theta.size();
      numberOfFunctionCalls = 0;
      numberOfGradientCalls = 0;
      numberOfIterations = 0;

      // We'll need non-const versions of the input arguments.
      argument_type thetaLocal(theta.size());
      copyArgumentType(theta, thetaLocal);
      result_type currentValue = startValue;
      if(currentValue <= this->m_minimumFunctionValue) {
        argument_type mutableGradient(theta.size());
        copyArgumentType(startGradient, mutableGradient);
        return brick::common::makeTriple(
          thetaLocal, currentValue, mutableGradient);
      }
      
      // Of course, we should only copy startGradient if it's actually
      // been initialized.
      argument_type currentGradient;
      if(startGradient.size() != 0) {
        copyArgumentType(startGradient, currentGradient);
      } else {
        currentGradient = this->m_functor.gradient(thetaLocal);
        ++numberOfGradientCalls;
      }

      // If gradient magnitude is zero, we're already at an extremum or
      // a saddle point.  In either case, we don't know which way to go.
      // Instead, just terminate.
      if(dotArgumentType(currentGradient, currentGradient) == 0) {
        if(this->m_verbosity > 0) {
          std::cout << "\nTerminating OptimizerBFGS::doBfgs() with "
                    << "zero gradient" << std::endl;
        }
        return brick::common::makeTriple(
          thetaLocal, currentValue, currentGradient);
      }
    
      // Check that gradient dimension is correct.
      if(currentGradient.size() != dimensionality) {
        std::ostringstream message;
        message << "startPoint has dimensionality " << dimensionality
                << " but objective function returns gradient with "
                << "dimensionality " << currentGradient.size() << ".";
        BRICK_THROW(brick::common::ValueException,
		    "OptimizerBFGS<Functor>::doBfgs()",
		    message.str().c_str());
      }
    
      // Set up inverse hessian estimate.
      brick::numeric::Array2D<double> inverseHessian(dimensionality, dimensionality);
      inverseHessian = 0.0;
      for(size_t index = 0; index < dimensionality; ++index) {
        // initialize to identity.
        inverseHessian(index, index) = 1.0;
      }

      // Set up initial search step.
      argument_type searchStep(dimensionality);
      for(size_t index = 0; index < dimensionality; ++index) {
        searchStep[index] = -(currentGradient[index]);
      }

      // Compute maximum allowable step size, allowing for numerical issues.
      double thetaMagnitudeSquared = 0.0;
      for(size_t index = 0; index < dimensionality; ++index) {
        thetaMagnitudeSquared += thetaLocal[index] * thetaLocal[index];
      }
      double thetaMagnitude = std::sqrt(thetaMagnitudeSquared);
      double maximumStepMagnitude =
        (this->m_maximumStepMagnitudeFactor
         * std::max(thetaMagnitude, static_cast<double>(dimensionality)));

      // Now set line search parameters
      // Make sure line search optimizer has the right objective function.
      // Since optimizer::setObjectiveFunction() doesn't update
      // m_optimizerLineSearch.
      this->m_optimizerLineSearch.setObjectiveFunction(this->m_functor);
      // Set line search parameters.
      this->m_optimizerLineSearch.setParameters(
        this->m_lineSearchArgumentTolerance, this->m_lineSearchAlpha,
        maximumStepMagnitude);
    
      // Perform the bfgs iteration.
      while(1) {
        // Perform line search.
        this->m_optimizerLineSearch.setStartPoint(thetaLocal, currentValue,
                                                  currentGradient);
        this->m_optimizerLineSearch.setInitialStep(searchStep);
        argument_type thetaNew = this->m_optimizerLineSearch.optimum();
        currentValue = this->m_optimizerLineSearch.optimalValue();
        numberOfFunctionCalls +=
          this->m_optimizerLineSearch.getNumberOfFunctionCalls();

        if(this->m_verbosity > 1) {
          std::cout << "\rCalls: " << std::setw(10) << numberOfFunctionCalls
                    << ", " << std::setw(10) << numberOfGradientCalls
                    << "     Current value: "
                    << std::setw(15) << currentValue << std::flush;
        }
        
        // Revise our initial step to match the data, and update current
        // position in parameter space.
        for(size_t index = 0; index < dimensionality; ++index) {
          searchStep[index] = thetaNew[index] - thetaLocal[index];
          thetaLocal[index] = thetaNew[index];
        }

        // Test for adequately small objective value.
        if(currentValue <= this->m_minimumFunctionValue) {
          if(this->m_verbosity > 0) {
            std::cout << "\nTerminating OptimizerBFGS::doBfgs() with "
                      << "objective value (" << currentValue 
                      << ") less than or equal to threshold ("
                      << this->m_minimumFunctionValue
                      << ")." << std::endl;
          }

          // Gradient at thetaLocal has not been computed yet.  Return
          // empty gradient.
          return brick::common::makeTriple(
            thetaLocal, currentValue, argument_type());
        }

        // Test for "insufficient parameter change" convergence.
        if(contextSensitiveScale(searchStep, thetaLocal)
           < this->m_argumentTolerance) {
          if(this->m_verbosity > 0) {
            std::cout << "\nTerminating OptimizerBFGS::doBfgs() with "
                      << "search step magnitude ("
                      << dotArgumentType(searchStep, searchStep)
                      << ") small compared to argument magnitude ("
                      << dotArgumentType(thetaLocal, thetaLocal)
                      << ")." << std::endl;
          }

          // Gradient at thetaLocal has not been computed yet.  Return
          // empty gradient.
          return brick::common::makeTriple(
            thetaLocal, currentValue, argument_type());
        }

        // Temporarily save gradient value.
        brick::numeric::Array1D<double> deltaGradient(dimensionality);
        for(size_t index = 0; index < dimensionality; ++index) {
          deltaGradient[index] = -currentGradient[index];
        }

        // Update gradient.
        currentGradient = this->m_functor.gradient(thetaLocal);
        ++numberOfGradientCalls;

        // Sanity check
        if(currentGradient.size() != dimensionality) {
          std::ostringstream message;
          message << "dimensionality of gradient changed mid-stream from "
                  << dimensionality << " to " << currentGradient.size() << ".";
          BRICK_THROW(brick::common::RunTimeException,
		      "OptimizerBFGS<Functor>::doBfgs()",
		      message.str().c_str());
        }

        // And compute change in gradient.  Remember that deltaGradient was
        // previously set to -currentGradient.
        for(size_t index = 0; index < dimensionality; ++index) {
          deltaGradient[index] += currentGradient[index];
        }

        // Test for "small gradient" convergence.
        if(this->gradientConvergenceMetric(thetaLocal, currentValue,
                                           currentGradient)
           < this->m_gradientTolerance) {
          if(this->m_verbosity > 0) {
            std::cout << "\nTerminating OptimizerBFGS::doBfgs() with "
                      << "small gradient." << std::endl;
          }
          return brick::common::makeTriple(
            thetaLocal, currentValue, currentGradient);
        }

        // Now prepare to update estimate of the inverse hessian.
        argument_type inverseHessianTimesDeltaGradient(dimensionality);
        matrixMultiplyArgumentType(inverseHessian, deltaGradient,
                                   inverseHessianTimesDeltaGradient);
        double fac = dotArgumentType(deltaGradient, searchStep);
        if(fac < 0.0) {
          fac = 0.0;
        }
        double fae = dotArgumentType(
          deltaGradient, inverseHessianTimesDeltaGradient);
        double mag2DeltaGradient = dotArgumentType(deltaGradient, deltaGradient);
        double mag2SearchStep = dotArgumentType(searchStep, searchStep);

        // Are gradient change and search direction sufficiently aligned?
        if((fac * fac)
           > (this->m_numericEpsilon * mag2DeltaGradient * mag2SearchStep)) {
          // Some useful scalars.
          double oneOverFac = 1.0 / fac;
          double fad = 1.0 / fae;

          // Use deltaGradient as temporary storage
          for(size_t index = 0; index < dimensionality; ++index) {
            deltaGradient[index] =
              (oneOverFac * searchStep[index]
               - fad * inverseHessianTimesDeltaGradient[index]);
          }

          // Now update inverse Hession estimate.
          for(size_t row = 0; row < inverseHessian.rows(); ++row) {
            for(size_t column = 0; column < inverseHessian.columns(); ++column) {
              // First DFP term.
              // 
              // Note(xxx): redundant calculation?  this multiplication
              // is done immediately above.
              double increment0 = (oneOverFac * searchStep[row]
                                   * searchStep[column]);
              // Second DFP term.
              double decrement1 = (fad * inverseHessianTimesDeltaGradient[row]
                                   * inverseHessianTimesDeltaGradient[column]);
              // BFGS term.  Remember that deltaGradient has been preempted for
              // temporary storage at this point.
              double increment2 = (fae * deltaGradient[row]
                                   * deltaGradient[column]);
              // Now do the update.
              inverseHessian(row, column) += (increment0 - decrement1
                                              + increment2);
            }
          }

          // Array2D<double> increment0 =
          //   outerProduct(oneOverFac * searchStep, searchStep);

          // Second DFP term.
          // Array2D<double> decrement1 =
          //   outerProduct(fad * inverseHessianTimesDeltaGradient,
          //                inverseHessianTimesDeltaGradient);

          // BFGS term.  Remember that deltaGradient has been preempted for
          // temporary storage at this point.
          // Array2D<double> increment2 =
          //   outerProduct(fae * deltaGradient, deltaGradient);

          // Actually do the update.
          // inverseHessian += increment0;
          // inverseHessian -= decrement1;
          // inverseHessian += increment2;
        }

        // Use Newton's method to choose the next update.
        matrixMultiplyArgumentType(inverseHessian, currentGradient, searchStep);
        for(size_t index = 0; index < dimensionality; ++index) {
          searchStep[index] = -(searchStep[index]);
        }

        // Check for convergence failure.
        ++numberOfIterations;
        if(numberOfIterations > this->m_iterationLimit) {
          // We're going to bail out, but save the result anyway, just
          // in case someone cares.
          this->setOptimum(thetaLocal, currentValue, false);

          // Now throw the exception.
          std::ostringstream message;
          message << "Iteration limit of " << this->m_iterationLimit
                  << " exceeded.";
          BRICK_THROW(brick::common::RunTimeException, 
		      "OptimizerBFGS<Functor>::doBfgs()",
		      message.str().c_str());
        }
      }
    }

  } // namespace optimization

} // namespace brick

#endif /* #ifndef BRICK_OPTIMIZATION_OPTIMIZERBFGS_HH */
