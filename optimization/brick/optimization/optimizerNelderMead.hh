/**
***************************************************************************
* @file brick/optimization/optimizerNelderMead.hh
* Header file declaring OptimizerNelderMead class.
*
* Copyright (C) 2003-2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
***************************************************************************
**/

#ifndef BRICK_OPTIMIZATION_OPTIMIZERNELDERMEAD_HH
#define BRICK_OPTIMIZATION_OPTIMIZERNELDERMEAD_HH

#include <functional>
#include <vector>
#include <brick/optimization/optimizer.hh>
#include <brick/common/exception.hh>

namespace brick {

  namespace optimization {

    /**
     ** OptimizerNelderMead implements the non-gradient based downhill
     ** simplex optimization method of Nelder and Mead, as described in
     ** [1] and [2].  This algorithm seeks the parameter value which
     ** minimizes the objective function.  The template parameter
     ** (Functor) defines the type to use as the objective function of
     ** the minimization.
     **
     ** [1] W. H. Press et al., Numerical Recipes in C The Art of
     ** Scientific Computing, Cambridge University Press, 1988.
     **
     ** [2] J.A. Nelder and R. Mead. A simplex method for function
     ** minimization. Computer Journal, 7:303--313, 1965.
     **/
    // template <std::unary_function Functor>
    template <class Functor, class FloatType = double>
    class OptimizerNelderMead
      : public Optimizer<Functor> {
    public:

      // Typedefs for convenience
      typedef typename Functor::argument_type argument_type;
      typedef typename Functor::result_type result_type;


      /**
       * Default constructor sets parameters to reasonable values for
       * functions which take values and arguments in the "normal" range
       * of 0 to 100 or so.
       */
      OptimizerNelderMead();


      /**
       * Constructor which specifies the specific Functor instance to
       * use.  Using this constructor exclusively avoids the danger of
       * calling optimalValue() or optimum() before a Functor instance
       * has been specified.
       *
       * @param functor A copy of this argument will be stored
       * internally for use in optimization.
       */
      explicit OptimizerNelderMead(const Functor& functor);


      /**
       * Copy constructor.
       *
       * @param source The OptimizerNelderMead instance to be copied.
       */
      OptimizerNelderMead(const OptimizerNelderMead& source);


      /**
       * Destructor.
       */
      ~OptimizerNelderMead();


      /**
       * Queries the number of iterations required to complete the
       * previous minimization.  If the minimization parameter
       * "restarts" is 0, there will be only one number to report.  If
       * restarts is greater than 0, the first element of the return
       * value will reflect the number of iterations in the initial
       * minimization, and subsequent numbers will reflect the number of
       * iterations in the following restarted minimizations.
       *
       * @return a vector of iteration counts.
       */
      virtual std::vector<size_t>
      getNumberOfFunctionCalls();


      /**
       * This method sets the spacing of the initial points used in the
       * nonlinear optimization, without affecting any other
       * optimization parameters.
       *
       * @param delta The method of Nelder and Mead requires not only a
       * point in parameter space at which to start the minimization,
       * but also N additional initial points, where N is the dimensionality of
       * the parameter space.  These will be automatically generated
       * according to: p[i+1] = p[0] + delta[i]*e[i] where e[i] is all
       * zeros except for a one in the i(th) position.
       */
      void
      setDelta(const argument_type& delta);


      /**
       * Sets how many times the optimization will be restarted.  No smaller
       * than 1 unless your objective function is very simple, since
       * downhill simplex search often benefits from a restart or two.
       *
       * @param numberOfRestarts
       */
      void
      setNumberOfRestarts(size_t numberOfRestarts) {
        this->m_numberOfRestarts = numberOfRestarts;
      }


      /**
       * Sets minimization parameters.  Default values are reasonable
       * for functions which take values and arguments in the "normal"
       * range of 0 to 100 or so.
       *
       * Since it is difficult to provide a default value for the
       * template parameter type Functory::argument_type, no default is
       * provided.  If a default were provided, it would be a vector of
       * all ones, with length equal to the dimensionality of the space
       * over which the optimization is to be performed.
       *
       * @param delta The method of Nelder and Mead requires not only a
       * point in parameter space at which to start the minimization,
       * but also N additional initial points, where N is the dimensionality of
       * the parameter space.  These will be automatically generated
       * according to: p[i+1] = p[0] + delta[i]*e[i] where e[i] is all
       * zeros except for a one in the i(th) position.
       *
       * @param functionCallLimit Each minimization will terminate after
       * this many iterations, even if a minimum has not been found.
       *
       * @param numberOfRestarts Following successful termination, the
       * minimization will be re-run this many times to refine the
       * result accuracy.  Generally you should set this to at least one
       * for this algorithm, since restarting the algorithm frequently
       * leads to significantly better final results.
       *
       * @param alpha Specifies the size of the initial reflection step.
       *
       * @param beta Specifies the size of the second reflection step.
       *
       * @param gamma Specifies the "shrink" factor if neither
       * reflection step is successful.
       *
       * @param minimumSimplexValueSpan Terminate if, among the N + 1
       * current points,
       *
       *   (2*(abs(f(pMax) - f(pMin)) / (abs(f(pMax)) + abs(f(pMin)))))
       *   < minSimplexValueSpan
       *
       * where pMax is the point generating the biggest value, and pMin
       * is the point generating the smallest value.
       *
       * @param verbosity Specifies how much standard output should be
       * generated.  Higher numbers mean more output.
       */
      void
      setParameters(argument_type delta,
                    size_t functionCallLimit=5000,
                    size_t numberOfRestarts=1,
                    FloatType alpha=1.0,
                    FloatType beta=0.5,
                    FloatType gamma=2.0,
                    FloatType minimumSimplexValueSpan=0.0001,
                    size_t verbosity=0);


      /**
       * Sets the initial conditions for the minimization.  Gradient
       * based search will start at this location in parameter space.
       *
       * @param startPoint Indicates a point in the parameter space of
       * the objective function.
       */
      virtual void
      setStartPoint(argument_type startPoint);

      /**
       * This function sets the level of text sent to standard
       * output by the class.
       *
       * @param verbosity This argument specifies how verbose to be.  0
       * means no output, 3 means lots.
       */
      virtual void
      setVerbosity(int verbosity) {this->m_verbosity = verbosity;}

      /**
       * Assignment operator.
       *
       * @param source The OptimizerNelderMead instance to be copied.
       * @return Reference to *this.
       */
      OptimizerNelderMead&
      operator=(const OptimizerNelderMead& source);

    protected:

      /**
       * This protected member function collapses a vector of points by
       * summing the corresponding elements of each point (useful for
       * averaging a bunch of locations in parameter space).
       *
       * @param currentPoints This argument specifies the points to sum.
       *
       * @param axisSums This argument is used to return the result.
       */
      void
      computeAxisSums(
        const std::vector<argument_type>& currentPoints,
        argument_type& axisSums);


      /**
       * This protected member function runs the actual simplex search.
       * It modifies all arguments.
       *
       * @param currentPoints This argument specifies the initial set of
       * points, and is used to return the final set of points.
       *
       * @param currentValues This argument specifies the initial set of
       * function values, and is used to return the final set of
       * function values.
       *
       * @param numberOfFunctionCalls This argument returns the total
       * number of function calls required by the minimization.
       */
      void
      doNelderMead(std::vector<argument_type>& currentPoints,
                   std::vector<result_type>& currentValues,
                   size_t& numberOfFunctionCalls);


      /**
       * This protected member function is used to decide whether a
       * proposed step is in fact a good one, and also updates its
       * arguments if the move is accepted.
       *
       * @param currentPoints This argument specifies the current set of
       * points in parameter space.
       *
       * @param currentValues This argument specifies the function value
       * at each of the points in currentPoints.
       *
       * @param axisSums This argument passes in the sums of the
       * elements of currentPoints.  See protected member function
       * computeAxisSums().
       *
       * @param factor This argument specifies the size of the step.
       *
       * @return This function returns the objective function value at
       * the evaluated point.
       */
      result_type
      evaluateMove(std::vector<argument_type>& currentPoints,
                   std::vector<result_type>& currentValues,
                   const argument_type& axisSums,
                   FloatType factor);


      /**
       * This protected member function performs the minimization.
       *
       * @return A std::pair of the vector parameter which brings the
       * specified Functor to a minimum, and the corresponding Functor
       * value.
       */
      std::pair<typename Functor::argument_type, typename Functor::result_type>
      run();

      argument_type m_delta;
      size_t m_functionCallLimit;
      size_t m_numberOfRestarts;
      FloatType m_alpha;
      FloatType m_beta;
      FloatType m_gamma;
      FloatType m_minimumSimplexValueSpan;

      std::vector<size_t> m_functionCallCount;
      argument_type m_theta0;
      size_t m_verbosity;

    }; // class OptimizerNelderMead

  } // namespace optimization

} // namespace brick


/* ======= Declarations to maintain compatibility with legacy code. ======= */

namespace brick {

  using optimization::OptimizerNelderMead;

} // namespace brick


/*******************************************************************
 * Member function definitions follow.  This would be a .cpp file
 * if it weren't templated.
 *******************************************************************/

#include <cmath>
#include <iomanip>

namespace brick {

  namespace optimization {

    /// @cond privateCode
    namespace privateCode {

      template <class Type0, class Type1>
      class extractSecond
      {
      public:
          using argument_type = std::pair<Type0, Type1>;
          using return_type = Type1;

        Type1 operator()(const std::pair<Type0, Type1>& inputPair) {
          return inputPair.second;
        }
      };

      template <class Type>
      std::vector<size_t>
      argsort(const std::vector<Type>& inputVector)
      {
        std::vector< std::pair<Type, size_t> >
          sortVector(inputVector.size());
        std::vector<size_t> resultVector(inputVector.size());

        for(size_t index = 0; index < inputVector.size(); ++index) {
          sortVector[index] =
            std::pair<Type, size_t>(inputVector[index], index);
        }
        std::sort(sortVector.begin(), sortVector.end());
        std::transform(sortVector.begin(), sortVector.end(), resultVector.begin(),
                       extractSecond<Type, size_t>());
        return resultVector;
      }

      template <class Type>
      std::vector<Type>
      take(const std::vector<Type>& inputVector,
            const std::vector<size_t>& indexVector)
      {
        if(inputVector.size() != indexVector.size()) {
          std::ostringstream message;
          message << "Incorrectly Sized Arguments:\n"
                  << " First Argument: " <<inputVector.size()
                  << ", Second Argument: " << indexVector.size() << ".";
          BRICK_THROW(brick::common::ValueException, "take()",
                    message.str().c_str());
        }
        std::vector<Type> resultVector(inputVector.size());
        for(size_t index = 0; index < resultVector.size(); ++index) {
          if(indexVector[index] >= inputVector.size()) {
            std::ostringstream message;
            message << "Out of Range Index:"
                    << indexVector[index] << ".";
            BRICK_THROW(brick::common::ValueException, "take()",
                      message.str().c_str());
          }
          resultVector[index] = inputVector[indexVector[index]];
        }
        return resultVector;
      }

    } // namespace privateCode
    /// @endcond

  } // namespace optimization

} // namespace brick


namespace brick {

  namespace optimization {

    // Sets parameters to reasonable values for functions which take
    // values and arguments in the "normal" range of 0 to 100 or so.
    // template <std::unary_function Functor>
    template <class Functor, class FloatType>
    OptimizerNelderMead<Functor, FloatType>::
    OptimizerNelderMead()
      : Optimizer<Functor>(),
        m_functionCallCount(),
        m_verbosity(0)
    {
      this->setParameters(argument_type());
    }


    // Constructor which specifies the specific Functor instance to
    // use.
    // template <std::unary_function Functor>
    template <class Functor, class FloatType>
    OptimizerNelderMead<Functor, FloatType>::
    OptimizerNelderMead(const Functor& functor)
      : Optimizer<Functor>(functor),
        m_functionCallCount(),
        m_verbosity(0)
    {
      this->setParameters(argument_type());
    }

    // Copy constructor.
    // template <std::unary_function Functor>
    template <class Functor, class FloatType>
    OptimizerNelderMead<Functor, FloatType>::
    OptimizerNelderMead(const OptimizerNelderMead& source)
      : Optimizer<Functor>(source),
        m_functionCallLimit(source.m_functionCallLimit),
        m_numberOfRestarts(source.m_numberOfRestarts),
        m_alpha(source.m_alpha),
        m_beta(source.m_beta),
        m_gamma(source.m_gamma),
        m_minimumSimplexValueSpan(source.m_minimumSimplexValueSpan),
        m_functionCallCount(source.m_functionCallCount),
        m_verbosity(source.m_verbosity)
    {
      copyArgumentType(source.m_delta, this->m_delta);
    }

    // Destructor.
    // template <std::unary_function Functor>
    template <class Functor, class FloatType>
    OptimizerNelderMead<Functor, FloatType>::
    ~OptimizerNelderMead()
    {
      // Empty
    }

    // Queries the number of functionCalls required to complete the
    // previous minimization.
    // template <std::unary_function Functor>
    template <class Functor, class FloatType>
    std::vector<size_t>
    OptimizerNelderMead<Functor, FloatType>::
    getNumberOfFunctionCalls()
    {
      return this->m_functionCallCount;
    }


    // This method sets the spacing of the initial points used in the
    // nonlinear optimization, without affecting any other optimization
    // parameters.
    template <class Functor, class FloatType>
    void
    OptimizerNelderMead<Functor, FloatType>::
    setDelta(const typename OptimizerNelderMead<Functor, FloatType>::argument_type& delta)
    {
      copyArgumentType(delta, this->m_delta);
    }


    // Sets minimization parameters.
    // template <std::unary_function Functor>
    template <class Functor, class FloatType>
    void
    OptimizerNelderMead<Functor, FloatType>::
    setParameters(typename OptimizerNelderMead<Functor, FloatType>::argument_type delta,
                  size_t functionCallLimit,
                  size_t numberOfRestarts,
                  FloatType alpha,
                  FloatType beta,
                  FloatType gamma,
                  FloatType minimumSimplexValueSpan,
                  size_t verbosity)
    {
      copyArgumentType(delta, this->m_delta);
      this->m_functionCallLimit = functionCallLimit;
      this->m_numberOfRestarts = numberOfRestarts;
      this->m_alpha = alpha;
      this->m_beta = beta;
      this->m_gamma = gamma;
      this->m_minimumSimplexValueSpan = minimumSimplexValueSpan;
      this->m_verbosity = verbosity;

      // Inherited member
      Optimizer<Functor>::m_needsOptimization = true;
    }

    // Sets the initial conditions for the minimization.
    // template <std::unary_function Functor>
    template <class Functor, class FloatType>
    void
    OptimizerNelderMead<Functor, FloatType>::
    setStartPoint(argument_type startPoint)
    {
      this->m_theta0 = startPoint;

      // Inherited member
      Optimizer<Functor>::m_needsOptimization = true;
    }

    // Assignment operator.
    // template <std::unary_function Functor>
    template <class Functor, class FloatType>
    OptimizerNelderMead<Functor, FloatType>&
    OptimizerNelderMead<Functor, FloatType>::
    operator=(const OptimizerNelderMead<Functor, FloatType>& source)
    {
      Optimizer<Functor>::operator=(source);
      this->m_delta = source.m_delta;
      this->m_functionCallLimit = source.m_functionCallLimit;
      this->m_numberOfRestarts = source.m_numberOfRestarts;
      this->m_alpha = source.m_alpha;
      this->m_beta = source.m_beta;
      this->m_gamma = source.m_gamma;
      this->m_minimumSimplexValueSpan = source.m_minimumSimplexValueSpan;
      this->m_functionCallCount = source.m_functionCallCount;
    }

    // template <std::unary_function Functor>
    template <class Functor, class FloatType>
    void
    OptimizerNelderMead<Functor, FloatType>::
    computeAxisSums(
      const std::vector<typename OptimizerNelderMead<Functor, FloatType>::argument_type>&
      currentPoints,
      typename OptimizerNelderMead<Functor, FloatType>::argument_type& axisSums)
    {
      if(currentPoints.size() == 0) {
        BRICK_THROW(brick::common::LogicException,
		    "OptimizerNelderMead::computeAxisSums()",
		    "Vector Size of Current Points is 0... "
		    "something's not right.");
        return;
      }
      copyArgumentType(currentPoints[0], axisSums);
      for(size_t index = 1; index < currentPoints.size(); ++index) {
        axisSums += currentPoints[index];
      }
    }

    // Run the actual simplex search.  Modifies all arguments.
    // template <std::unary_function Functor>
    template <class Functor, class FloatType>
    void
    OptimizerNelderMead<Functor, FloatType>::
    doNelderMead(
      std::vector<typename OptimizerNelderMead<Functor, FloatType>::argument_type>&
      currentPoints,
      std::vector<typename OptimizerNelderMead<Functor, FloatType>::result_type>&
      currentValues,
      size_t& numberOfFunctionCalls)
    {
      size_t dimension = currentValues.size() - 1;
      argument_type axisSums;
      this->computeAxisSums(currentPoints, axisSums);
      while(1) {
        std::vector<size_t> indices = privateCode::argsort(currentValues);
        currentValues = privateCode::take(currentValues, indices);
        currentPoints = privateCode::take(currentPoints, indices);

        if(this->m_verbosity >= 3) {
          std::cout << "\rCurrent best: " << std::setw(15) << currentValues[0]
                    << std::flush;
        }

        FloatType simplexValueSpan = 0;
        if((std::fabs(currentValues[dimension])
            + std::fabs(currentValues[0])) > 0) {
          simplexValueSpan =
            (2.0 * std::fabs(currentValues[dimension] - currentValues[0])
             / (std::fabs(currentValues[dimension])
                + std::fabs(currentValues[0])));
        }
        if(this->m_verbosity >= 3) {
          std::cout << "         simplexValueSpan: "
                    << std::setw(15) << simplexValueSpan << std::flush;
        }
        if((simplexValueSpan < this->m_minimumSimplexValueSpan)
           || (numberOfFunctionCalls >= this->m_functionCallLimit)) {
          if(this->m_verbosity >= 3) {
            std::cout << "\n";
          }
          if(this->m_verbosity >= 1) {
            std::cout << "Terminating with:\n"
                      << "  simplexValueSpan = " << simplexValueSpan
                      << " (" << this->m_minimumSimplexValueSpan << ")\n";
            std::cout << "  numberOfFunctionCalls = " << numberOfFunctionCalls
                      << " (" << this->m_functionCallLimit << ")" << std::endl;
          }
          break;
        }
        result_type newValue =
          this->evaluateMove(currentPoints, currentValues, axisSums,
                             (-1.0 * this->m_alpha));
        numberOfFunctionCalls += 1;
        this->computeAxisSums(currentPoints, axisSums);
        if(newValue <= currentValues[0]) {
          newValue = evaluateMove(currentPoints, currentValues, axisSums,
                                  this->m_gamma);
          numberOfFunctionCalls += 1;
          this->computeAxisSums(currentPoints, axisSums);
        } else if (newValue >= currentValues[dimension - 1]) {
          result_type oldMaxValue = currentValues[dimension];
          newValue = evaluateMove(currentPoints, currentValues, axisSums,
                                  this->m_beta);
          numberOfFunctionCalls += 1;
          if(newValue >= oldMaxValue) {
            for(size_t i = 1; i < dimension + 1; ++i) {
              currentPoints[i] = 0.5 * (currentPoints[i] + currentPoints[0]);
              currentValues[i] = this->m_functor(currentPoints[i]);
              numberOfFunctionCalls += 1;
            }
          }
          this->computeAxisSums(currentPoints, axisSums);
        }
      }
    }


    // Evaluate, and maybe accept, a new point.  You'll want to believe
    // the following math before reading this code.
    //
    // xMean + factor*(xMax - xMean)
    // = xMean + factor*xMax - factor*xMean
    // = (1 - factor)*xMean + factor*xMax
    // = ((1 - factor)/n)*sum(xSubi, xSubi != xMax) + factor*xMax
    // = ((1 - factor)/n)*sum(xSubi) - ((1 - factor)/n)*xMax + factor*xMax
    // = ((1 - factor)/n)*sum(xSubi) + (factor - ((1 - factor)/n))*xMax
    // template <std::unary_function Functor>
    template <class Functor, class FloatType>
    typename OptimizerNelderMead<Functor, FloatType>::result_type
    OptimizerNelderMead<Functor, FloatType>::
    evaluateMove(std::vector<argument_type>& currentPoints,
                 std::vector<result_type>& currentValues,
                 const argument_type& axisSums,
                 FloatType factor)
    {
      size_t dimension = currentPoints.size() - 1;
      FloatType centroidFactor = (1.0 - factor) / dimension;
      FloatType extrapolationFactor = factor - centroidFactor;
      argument_type newPoint =
        ((centroidFactor * axisSums)
         + (extrapolationFactor * currentPoints[dimension]));
      result_type newValue = this->m_functor(newPoint);
      if(newValue < currentValues[dimension]) {
        currentValues[dimension] = newValue;
        currentPoints[dimension] = newPoint;
      }
      return newValue;
    }

    // Perform the minimization (top level).
    // template <std::unary_function Functor>
    template <class Functor, class FloatType>
    std::pair<typename Functor::argument_type, typename Functor::result_type>
    OptimizerNelderMead<Functor, FloatType>::
    run()
    {
      size_t dimension = this->m_theta0.size();
      if(dimension == 0) {
        BRICK_THROW(brick::common::ValueException, "OptimizerNelderMead::run()",
                  "invalid starting point has zero size.");
      }
      if(this->m_delta.size() != dimension) {
        // Initialize delta using only operations that we expect to be
        // supported by argtype, and which we don't expect to be
        // defined in surprising ways.
        this->m_delta = (this->m_theta0 * 0.0) + 1.0;
      }

      // The algorithm requires dimension + 1 initial points and values.
      std::vector<argument_type> currentPoints(dimension + 1);
      std::vector<result_type> currentValues(dimension + 1);
      // First point is the one specified by the user.
      copyArgumentType(this->m_theta0, currentPoints[0]);
      currentValues[0] = this->m_functor(currentPoints[0]);

      // We'll repeat the whole minimization several times.
      this->m_functionCallCount.clear();
      for(size_t i = 0; i < this->m_numberOfRestarts + 1; ++i) {
        // The remaining points will be built using this->m_delta
        for(size_t j = 0; j < dimension; ++j) {
          copyArgumentType(currentPoints[0], currentPoints[j + 1]);
          (currentPoints[j + 1])[j] += this->m_delta[j];
          currentValues[j + 1] = this->m_functor(currentPoints[j + 1]);
        }
        size_t numberOfFunctionCalls = (i == 0) ? 1 : 0;
        numberOfFunctionCalls += dimension;
        // Run the minimization algorithm.
        this->doNelderMead(currentPoints, currentValues, numberOfFunctionCalls);
        // Have to add dimension to numberOfFunctionCalls because of the
        // initialization steps above.
        this->m_functionCallCount.push_back(numberOfFunctionCalls + dimension);
      }
      return std::make_pair(currentPoints[0], currentValues[0]);
    }

  } // namespace optimization

} // namespace brick


#endif /* #ifndef BRICK_OPTIMIZATION_OPTIMIZERNELDERMEAD_HH */
