/**
***************************************************************************
* @file brick/optimization/optimizer.hh
*
* Header file declaring Optimizer class.
*
* Copyright (C) 2003-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_OPTIMIZATION_OPTIMIZER_HH
#define BRICK_OPTIMIZATION_OPTIMIZER_HH

#include <functional>

namespace brick {

  /**
   ** This namespace contains classes for gradient- and non-gradient-
   ** based optimization.
   **/
  namespace optimization {

    /**
     ** Optimizer is an abstract base class defining the interface for
     ** general optimization classes.  The template parameter (Functor)
     ** defines the type to use as the objective function of the
     ** optimization.  The specific requirements for Functor are to be
     ** determined by the classes derived from Optimizer, however
     ** Functor will generally be a scalar-valued function with an array
     ** or vector type argument.
     **/
    // template <std::unary_function Functor>
    template<class Functor>
    class Optimizer {
    public:
      // Typedefs for convenience
      /// This is the Type of the objective function argument.
      typedef typename Functor::argument_type argument_type;

    
      /// This is the Type of the objective function return value.
      typedef typename Functor::result_type result_type;

    
      /**
       * This is the default constructor.  Note that the default
       * constructor in any derived classes must initialize
       * configuration parameters to reasonable values.
       */
      Optimizer();

    
      /**
       * This constructor specifies the specific Functor instance to
       * use.  Using this constructor exclusively avoids the danger of
       * calling optimalValue() or optimum() before a Functor instance
       * has been specified.
       *
       * @param functor A copy of this argument will be stored
       * internally for use in optimization.
       */
      explicit Optimizer(const Functor& functor);

    
      /** 
       * Copy constructor. This constructor simply copies the source
       * argument.
       * 
       * @param source The Optimizer instance to be copied.
       */
      Optimizer(const Optimizer& source);

    
      /**
       * Destructor.
       */
      virtual ~Optimizer();

    
      /** 
       * This method returns a copy of the Functor instance used for
       * optimization.
       * 
       * @return A Functor instance.
       */
      Functor objectiveFunction() {return m_functor;}

    
      /**
       * This is the assignment operator. It simply copies its input
       * argument.
       * 
       * @param source The Optimizer instance to be copied.
       * @return Reference to *this.
       */
      Optimizer&
      operator=(const Optimizer& source);

    
      /** 
       * This method finds the optimum of the current Functor, if
       * necessary, and returns the Functor value at that point.  Note
       * that you must have specified an objective function (Functor)
       * before calling this method.
       *
       * @return The Functor value at it's optimum.
       */
      result_type
      optimalValue();

    
      /** 
       * This method finds the optimum of the current Functor, if
       * necessary, and returns the Functor argument which produces that
       * optimum.  Note that you must have specified an objective
       * function (Functor) before calling this method.
       *
       * @return The Functor arguments which produce the optimal value.
       */
      argument_type
      optimum();

    
      /** 
       * This method specifies the Functor instance to use for the
       * optimization.  If this function is overridden by the base
       * class, it should normally either call
       * Optimizer::setObjectiveFunction(), or explicitly set the member
       * variable m_needsOptimization to true.
       *
       * @param functor A copy of this argument will be stored
       * internally for use in optimization.
       */
      void
      setObjectiveFunction(const Functor& functor);
    
    protected:

      /** 
       * Perform the optimization.  This pure virtual function must be
       * overridden by the base class.
       * 
       * @return A std::pair of the vector parameter which brings the
       * specified Functor to an optimum, and the corresponding optimal
       * Functor value.
       */
      virtual
      std::pair<typename Functor::argument_type, typename Functor::result_type>
      run() = 0;


      /** 
       * This protected member function provides a way for subclasses to
       * communicate intermediate optimization results outside of the
       * normal "return value of this->run()" method.
       * 
       * @param optimum This argument will be saved as the current optimum.
       * 
       * @param optimalValue This argument will be saved as the function
       * value a the current optimum.
       * 
       * @param needsFurtherOptimization This argument indicates whether
       * or not further refinement is necessary.
       */
      virtual
      void
      setOptimum(const typename Functor::argument_type& optimum,
                 const typename Functor::result_type& optimalValue,
                 bool needsFurtherOptimization) {
        m_optimum = optimum;
        m_optimalValue = optimalValue;
        m_needsOptimization = needsFurtherOptimization;
      }
    
    
      /// m_functor->operator()() should compute the objective function.
      Functor m_functor;

      /// Set to false if m_optimum contains a valid optimum, true otherwise.
      bool m_needsOptimization;

      /// Caches the result of the most recent optimization.
      argument_type m_optimum;

      /// Caches the result of the most recent optimization.
      result_type m_optimalValue;

    }; // class Optimizer

  } // namespace optimization

} // namespace brick


/* ======= Declarations to maintain compatibility with legacy code. ======= */

namespace brick {

  using optimization::Optimizer;

} // namespace brick


/*******************************************************************
 * Member function definitions follow.  This would be a .cpp file
 * if it weren't templated.
 *******************************************************************/

#include <dlrNumeric/array1D.hh>>
#include <brick/optimization/optimizerCommon.hh>>

namespace brick {

  namespace optimization {

    // Default constructor in derived classes must initialize
    // configuration parameters to reasonable values.
    template <class Functor>
    Optimizer<Functor>::
    Optimizer()
      : m_functor(),
        m_needsOptimization(true),
        m_optimum()
    {
      // Empty
    }
    
    // Constructor which specifies the specific Functor instance to use.
    template <class Functor>
    Optimizer<Functor>::
    Optimizer(const Functor& functor)
      : m_functor(functor),
        m_needsOptimization(true),
        m_optimum()
    {
      // Empty
    }

    // Copy constructor.
    template <class Functor>
    Optimizer<Functor>::
    Optimizer(const Optimizer& source)
      : m_functor(source.m_functor),
        m_needsOptimization(source.m_needsOptimization),
        m_optimum()
    {
      copyArgumentType(source.m_optimum, m_optimum);
    }

    // Destructor.
    template <class Functor>
    Optimizer<Functor>::
    ~Optimizer()
    {
      // Empty
    }

    // Assignment operator.
    template <class Functor>
    Optimizer<Functor>&
    Optimizer<Functor>::
    operator=(const Optimizer& source)
    {
      if(&source != this) {
        m_functor = source.m_functor;
        m_needsOptimization = source.m_needsOptimization;
        copyArgumentType(source.m_optimum, m_optimum);
      }
      return *this;
    }

    // Find the optimum of the current Functor, if necessary, and
    // return the Functor value at that point.
    template <class Functor>
    typename Optimizer<Functor>::result_type
    Optimizer<Functor>::
    optimalValue()
    {
      if(m_needsOptimization==true) {
        std::pair<argument_type, result_type> optimum_optimalValue = this->run();
        this->setOptimum(
          optimum_optimalValue.first, optimum_optimalValue.second, false);
      }
      return m_optimalValue;
    }

    // Find the optimum of the current Functor, if necessary, and
    // return the Functor argument which produces that optimum.
    template <class Functor>
    typename Optimizer<Functor>::argument_type
    Optimizer<Functor>::
    optimum()
    {
      if(m_needsOptimization==true) {
        std::pair<argument_type, result_type> optimum_optimalValue = this->run();
        this->setOptimum(
          optimum_optimalValue.first, optimum_optimalValue.second, false);
      }
      return m_optimum;
    }

    // Specify the Functor instance to use for the optimization.
    template <class Functor>
    void
    Optimizer<Functor>::
    setObjectiveFunction(const Functor& functor)
    {
      m_functor = functor;
      m_needsOptimization = true;
    }

  } // namespace optimization

} // namespace brick

#endif /* #ifndef BRICK_OPTIMIZATION_OPTIMIZER_HH */
