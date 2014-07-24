/**
***************************************************************************
* @file brick/computerVision/ransac.hh
*
* Header file declaring an implementation of the RANSAC algorithm.
*
* Copyright (C) 2008 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_RANSACCLASS_HH
#define BRICK_COMPUTERVISION_RANSACCLASS_HH

#include <vector>
#include <brick/computerVision/randomSampleSelector.hh>

namespace brick {

  namespace computerVision {

    /**
     ** This enum is used by the RANSAC algorithm to select between
     ** the various ways of deciding whether a particular sample
     ** matches the current model.  For now, you only have one choice.
     ** See RansacProblem::getNaiveNaiveErrorThreshold().
     **/
    enum RansacInlierStrategy {
      BRICK_CV_NAIVE_ERROR_THRESHOLD
    };
      
    
    /**
     ** This class template implements the RANSAC algorithm[1].
     **
     ** The template argument, Problem, provides all of the
     ** user-supplied problem-specific code, such as model estimation
     ** and error computation.  The easiest way to do make an
     ** appropriate Problem class is to derive from class
     ** RansacProblem, below.  For an example, see the file
     ** test/ransacTest.cpp.
     **
     ** [1] M. Fischler and R. Bolles. Random Sample Consensus: A
     ** Paradigm for Model Fitting with Applications to Image Analysis
     ** and Automated Cartography. Graphics and Image Processing,
     ** 24(6):381--395, 1981.
     **/
    template <class Problem>
    class Ransac {
    public:

      /**
       ** This typedef simply shadows template argument Problem.
       **/
      typedef Problem ProblemType;

      
      /**
       ** This typedef indicates the type of model that will be
       ** estimated by the RANSAC algorithm.  Its value is controlled
       ** by the user specified RansacProblem class.
       **/
      typedef typename Problem::ModelType ResultType;

      
      /** 
       * This constructor sets up the Ransac instance so that it is
       * ready to solve the model fitting problem, but does not run
       * the RANSAC algorithm.
       * 
       * @param problem This argument is a class instance implementing
       * the RansacProblem interface, which provides all of the
       * problem-specific code.
       * 
       * @param minimumConsensusSize This argument specifies the
       * smallest set of "agreeing" samples that should be taken as
       * proof that the correct model has been found (and grounds for
       * terminating the algorithm).  Setting this value to zero will
       * make the Ransac constructor compute an appropriate value
       * using arguments requiredConfidence and inlierProbability.
       * 
       * @param requiredConfidence This argument indicates how
       * confident we need to be that one run of the RANSAC algorithm
       * will find the correct model.  It affects the number of
       * iterations that the RANSAC algorithm will be allowed to run,
       * as well as the automatically computed value for
       * minimumConsensusSize (see above).
       * 
       * @param inlierProbability This argument indicates the
       * likelihood that any particular input value (data point,
       * sample, whatever) is an "inlier" for the purpose of
       * estimating a model.  If you had 10 sample points, and thought
       * that 3 of them had excessive noise, then you would supply a
       * value of 0.7 for this argument. It affects the number of
       * iterations that the RANSAC algorithm will be allowed to run.
       */
      Ransac(ProblemType const& problem,
             size_t minimumConsensusSize = 0,
             double requiredConfidence = 0.99,
             double inlierProbability = 0.5,
             unsigned int verbosity = 0);


      /**
       * The destructor cleans up any system resources and destroys *this.
       */
      virtual
      ~Ransac() {}


      /** 
       * This member function runs the RANSAC algorithm and returns
       * the computed model.
       * 
       * @return The return value is the best model estimate returned
       * by the RANSAC algorithm.
       */
      virtual ResultType
      getResult() {
        // xxx
        ResultType result;
        this->estimate(result);
        return result;
      }


      /** 
       * Controls how many times a model may be refined on each rasac
       * iteration.  Normally, in each iteration, the model is
       * repeatedly recomputed using the consensus set until it stops
       * improving.  This function allow the user to limit how many
       * refinements can happen during each iteration.
       * 
       * @param numberOfRefinements This argument specifies the
       * maximum number of allowable refinements per iteration.
       * Setting this to a negative number allows the refinements to
       * continue until convergence.
       */
      void
      setNumberOfRefinements(int numberOfRefinements) {
        m_numberOfRefinements = numberOfRefinements;
      }

    protected:

      void
      computeConsensusSet(ResultType& model, std::vector<bool>& consensusFlags);

      bool
      estimate(ResultType& model);

      bool
      isConverged(std::vector<bool> const& consensusFlags,
                  std::vector<bool>& previousConsensusFlags,
                  size_t& consensusSetSize,
                  size_t& previousConsensusSetSize,
                  size_t& strikes,
                  int refinementCount);


      size_t m_minimumConsensusSize;
      size_t m_numberOfRandomSampleSets;
      int m_numberOfRefinements;
      ProblemType m_problem;
      unsigned int m_verbosity;
    };


    /**
     ** This class template implements the "Problem" interface
     ** required by the Ransac class, above.  If you derive your
     ** Problem class from RansacProblem, you avoid some some of the
     ** hassles involved in random sampling, etc.  If you like, it's
     ** also OK to not derive your Problem class from RansacProblem,
     ** in which case you should be sure to provide both the interface
     ** shown here, and the RandomSampleSelector interface
     ** (RandomSampleSelector is a parent class of RansacProblem).
     **
     ** Template argument Sample specifies what type of sample (2D
     ** points, 3D points, etc.) are used to estimate the model.
     ** Template argument Model specifies the actual type of a model.
     **
     ** For example, you might set Sample to brick::numeric::Vector2D
     ** (representing 2D points), and Model to std::pair<double,
     ** double> (representing the slope and intercept of a line in 2D
     ** space).  For an example that does exactly this, see the
     ** LineFittingProblem class in file test/ransacTest.cpp.
     **/
    template <class Sample, class Model>
    class RansacProblem
      : public RandomSampleSelector<Sample>
    {
    public:
      
      // ========= Typedefs that must be present in order   =========
      // ========= to work with the Ransac class template.  =========
      // ========= You probably don't need to change these. =========

      /**
       ** This typedef simply mirrors the "Model" template argument
       ** described in the documentation for this class.
       **/
      typedef Model ModelType;


      /**
       ** This typedef simply mirrors the "Sample" template argument
       ** described in the documentation for this class.
       **/
      typedef Sample SampleType;


      /**
       ** This typedef specifies what type will be used to represent
       ** sequences of SampleType.  You can just use the default here,
       ** or if you need to implement your own random sample
       ** selection, you might use something simple like this:
       **
       ** @code
       **   typedef std::vector<SampleType> SampleSequenceType;
       ** @endcode
       **
       ** Or, if you're worried about the cost of copying vectors, you
       ** could pass a pair of iterators pointing to a sequence that
       ** you know will remain valid until the next call to
       ** this->getRandomSample() or this->getSubset().  The default
       ** typedef provide below does exactly this:
       **
       ** @code
       **   typedef std::pair<std::vector<SampleType>::const_iterator,
       **                     std::vector<SampleType>::const_iterator>
       **           SampleSequenceType;
       ** @endcode
       **
       ** Note(xxx): we might want to change this "remain valid"
       ** restriction so that RandomSampleSelector is more useful to
       ** non-Ransac algorithms.
       **/
      typedef typename RandomSampleSelector<SampleType>::SampleSequenceType
        SampleSequenceType;
      

      // ========= Member functions that must be provided by user =========
      // ========= You definitely do need to change these.        =========

      /** 
       * Subclasses may override this member function to reset
       * internal state prior to the beginning of each iteration of
       * the RANSAC algorithm.  Consider a problem that uses nonlinear
       * optimization inside estimateModel().  Such a problem might
       * remember the most recent model estimate, and use that
       * estimate as initial conditions for iterative refinement of
       * the model (with progressively changing consensus sets) during
       * a single RANSAC iteration.  This member function allows such
       * a problem to reset its initial conditions after each
       * iteration so that the results of a previous random sample
       * don't influence the current one.
       *
       * @param iterationNumber This argument starts at zero and
       * increments with each call to beginIteration(), indicating
       * which iteration is beginning.
       */
      virtual void
      beginIteration(size_t /* iterationNumber */) {}

      
      /** 
       * This member function should take a sequence of samples, and
       * compute the best fit model based on the sample values.  The
       * number of sample values contained in argument sampleSequence
       * will be equal to the value of constructor argument
       * sampleSize.
       * 
       * @param sampleSequence This argument specifies the sequence of
       * samples.  See the documentation for typedef
       * SampleSequenceType, above.
       * 
       * @return The return value is a model instance based on the
       * input data.
       */
      virtual ModelType
      estimateModel(SampleSequenceType const& sampleSequence) = 0;


      /** 
       * This function should take a model and a sequence of samples,
       * and fill in the output sequence with error values reflecting
       * how well each member of the sample sequence matches the
       * model.
       * 
       * @param model This argument specifies the model against which
       * to test.
       * 
       * @param sampleSequence This argument specifies the sequence of
       * samples to be tested.
       * 
       * @param ouputIter This argument is an iterator pointing to the
       * first element of the sequence of doubles that should be
       * filled with error values.  The sequence will have at least
       * this->getPoolSize() valid elements whenever this function is
       * called by Ransac.
       */
      template <class IterType>
      void
      computeError(ModelType const& model,
                   SampleSequenceType const& sampleSequence,
                   IterType ouputIter) {};


      /** 
       * This member function should return a threshold against which
       * error values (computed by this->computeError()) should be
       * compared in order to determine whether a particular sample is
       * an inlier or an outlier.  Sophisticated RANSAC
       * implementations sometimes compute thresholds for the
       * inlier/outlier decision based on the sensitivity of the model
       * estimation algorithm to noise, but we haven't implemented
       * this.  For now, your only real option is to return a single
       * threshold that will be applied to each sample error value.
       * Samples with error larger than this threshold will be
       * considered to be outliers.
       * 
       * @return The return value is the inlier/outlier threshold.
       */
      virtual double 
      getNaiveErrorThreshold() = 0;
      

      // ========= Public constructors, destructors, and    =========
      // ========= predefined member functions.             =========
      // ========= You probably don't need to change these. =========

      /** 
       * This constructor instantiate the RansacProblem instance using
       * a sequence of samples from which to randomly draw populations
       * and compute models.
       * 
       * @param sampleSize This argument specifies how many individual
       * samples are required by member function estimateModel() to
       * compute a model, and is used to control the length of the
       * sample sequence passed to estimateModel().
       * 
       * @param beginIter This argument is an iterator pointing to the
       * first element of a sequence of samples against which the
       * RANSAC algorithm will be run.
       * 
       * @param endIter This argument points, in the normal STL way,
       * one element past the end of the sequence started by argument
       * beginIter.
       */
      template <class IterType>
      RansacProblem(size_t sampleSize, IterType beginIter, IterType endIter)
        : RandomSampleSelector<SampleType>(beginIter, endIter),
          m_sampleSize(sampleSize) {}

      
      /**
       * The destructor cleans up any system resources and destroys *this.
       */
      virtual
      ~RansacProblem() {}


      /** 
       * This member function returns how many individual samples are
       * required by member function estimateModel() to compute a
       * model, and is used by the Ransac class to control the length
       * of the sample sequence passed to estimateModel().
       * 
       * @return The return value specifies the required length of the
       * sample sequence.
       */
      virtual size_t
      getSampleSize() {return m_sampleSize;}


      /** 
       * This member function controls how the inlier/outlier decision
       * is made during RANSAC operation.  For now, the only thing you
       * can return is Ransac::BRICK_CV_NAIVE_ERROR_THRESHOLD.
       * 
       * @return The return value 
       */
      RansacInlierStrategy
      getInlierStrategy() {
        return BRICK_CV_NAIVE_ERROR_THRESHOLD;
      }

    protected:

      size_t m_sampleSize;
      
    }; // class RansacProblem
    
  } // namespace computerVision
  
} // namespace brick



/* ============ Definitions of inline & template functions ============ */


#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <brick/common/exception.hh>
#include <brick/numeric/maxRecorder.hh>

namespace brick {

  namespace computerVision {

    // The default constructor currently does nothing.
    template <class Problem>
    Ransac<Problem>::
    Ransac(Problem const& problem,
           size_t minimumConsensusSize,
           double requiredConfidence,
           double inlierProbability,
           unsigned int verbosity)
      : m_minimumConsensusSize(minimumConsensusSize),
        m_numberOfRandomSampleSets(),
        m_numberOfRefinements(-1),
        m_problem(problem),
        m_verbosity(verbosity)
    {
      size_t sampleSize = m_problem.getSampleSize();

      if((requiredConfidence < 0.0) || (requiredConfidence >= 1.0)) {
        BRICK_THROW(common::ValueException, "Ransac::Ransac()",
                  "Probability value requiredConfidence is out of range.");
      }

      if((inlierProbability < 0.0) || (inlierProbability >= 1.0)) {
        BRICK_THROW(common::ValueException, "Ransac::Ransac()",
                  "Probability value inlierProbability is out of range.");
      }

      // The probability that any particular random sample is all inliers.
      double singlePickConfidence =
        std::pow(inlierProbability, static_cast<double>(sampleSize));

      // The probability that any particular random sample contains at
      // least one outlier.
      double singlePickDisconfidence = 1.0 - singlePickConfidence;

      // Number of times we have to sample in order to get a set of
      // all inliers with probability requiredConfidence.  This
      // follows Fischler's paper exactly.
      m_numberOfRandomSampleSets =
        std::ceil(std::log(1.0 - requiredConfidence)
                  / std::log(singlePickDisconfidence)) + 0.5;
      
      if(m_minimumConsensusSize == 0) {
        // User didn't specify a minimum consensus size.  Compute one
        // assuming that the probability of a sample matching an
        // incorrect model is <= 0.5, and using requiredConfidence as
        // a rough measure of how confident we have to be that we've
        // found the right model.  This follows Fischler's paper
        // approximately.
        int extraSamples = static_cast<int>(
          std::log(1.0 - requiredConfidence) / std::log(0.5) + 0.5);
        if(extraSamples < 0) {
          extraSamples = 0;
        }
        m_minimumConsensusSize = sampleSize + extraSamples;
      }
      
    }


    template <class Problem>
    void
    Ransac<Problem>::
    computeConsensusSet(typename Ransac<Problem>::ResultType& model,
                        std::vector<bool>& consensusFlags)
    {
      if(m_problem.getInlierStrategy() != BRICK_CV_NAIVE_ERROR_THRESHOLD) {
        BRICK_THROW(brick::common::NotImplementedException,
                    "Ransac::computeConsensusSet()",
                    "Currently only naive error thresholding is supported.");
      }
      if(consensusFlags.size() != m_problem.getPoolSize()) {
        consensusFlags.resize(m_problem.getPoolSize());
      }
      
      // Apply error function to entire set.
      typename Problem::SampleSequenceType testSet = m_problem.getPool();
      std::vector<double> errorMetrics(m_problem.getPoolSize());
      m_problem.computeError(model, testSet, errorMetrics.begin());

      // Find out which samples are within tolerance.
      double threshold = m_problem.getNaiveErrorThreshold();
      std::transform(
        errorMetrics.begin(), errorMetrics.end(), consensusFlags.begin(),
        std::bind2nd(std::less<double>(), threshold));
    }
    
    
    template <class Problem>
    bool
    Ransac<Problem>::
    estimate(typename Ransac<Problem>::ResultType& model)
    {
      brick::numeric::MaxRecorder<size_t, ResultType> maxRecorder;
      for(size_t iteration = 0; iteration < m_numberOfRandomSampleSets;
          ++iteration) {
        if(m_verbosity >= 3) {
          std::cout << "Ransac: running sample #" << iteration
                    << " of " << m_numberOfRandomSampleSets << std::endl;
        }
        
        // Select samples
        typename ProblemType::SampleSequenceType trialSet =
          m_problem.getRandomSample(m_problem.getSampleSize());

        std::vector<bool> consensusFlags(m_problem.getPoolSize());
        std::vector<bool> previousConsensusFlags(m_problem.getPoolSize(),
                                                 false);

        // Some problem classes may, for example, retain internal
        // state during the iterative refinement loop below.  This
        // call allows those problems to reset that state prior to
        // starting over with a new random sample.
        m_problem.beginIteration(iteration);
        
        size_t consensusSetSize = 0;
        size_t previousConsensusSetSize = 0;
        size_t strikes = 0;
        int refinementCount = 0;
        while(1) {
          // Fit the model to the reduced (randomly sampled) set.
          model = m_problem.estimateModel(trialSet);

          // Identify the consensus set, made up of samples that are
          // sufficiently consistent with the model estimate.
          this->computeConsensusSet(model, consensusFlags);

          if(m_verbosity >= 3) {
            std::cout
              << "Ransac:   consensus set size is "
              << std::count(consensusFlags.begin(), consensusFlags.end(), true)
              << " (vs. " << m_minimumConsensusSize << ")" << std::endl;
          }
          
          // See if this iteration has converged yet.
          if(this->isConverged(consensusFlags, previousConsensusFlags,
                               consensusSetSize, previousConsensusSetSize,
                               strikes, refinementCount)) {
            break;
          }

          // Not converged yet... loop so we can recompute the model
          // using the new consensus set.
          trialSet = m_problem.getSubset(
            consensusFlags.begin(), consensusFlags.end());
          ++refinementCount;
        }

        // OK, we've converged to a "best" result for this iteration.
        // Is it good enough to terminate?
        if(consensusSetSize > m_minimumConsensusSize) {
          // Reference argument model is already set appropriately.
          return true;
        }

        // Not ready to terminate yet, but remember this model (if
        // it's the best so far) in case we don't find any better.
        maxRecorder.test(consensusSetSize, model);
      }

      // Looks like we never found a gold plated correct answer.  Just
      // return report the best we found, and return false to indicate
      // our frustration.
      model = maxRecorder.getPayload();

      if(m_verbosity >= 3) {
        std::cout
          << "Ransac: terminating with best consensus set size of "
          << maxRecorder.getMaximum() << std::endl;
      }
      return false;
    }


    template <class Problem>
    bool
    Ransac<Problem>::
    isConverged(std::vector<bool> const& consensusFlags,
                std::vector<bool>& previousConsensusFlags,
                size_t& consensusSetSize,
                size_t& previousConsensusSetSize,
                size_t& strikes,
                int refinementCount)
    {
      // Do we even have enough matching points to continue
      // iteration?
      consensusSetSize = std::count(
        consensusFlags.begin(), consensusFlags.end(), true);
      if(consensusSetSize < m_problem.getSampleSize()) {
        // No. The model is so bad that not even the points we used to
        // estimate it are within the consensus set!  We're converged,
        // after a fashion.
        return true;
      }

      // Are we permitted to refine the model again?
      if(m_numberOfRefinements >= 0
         && refinementCount >= m_numberOfRefinements) {
        // No.  We're converged, after a fashion.
        return true;
      }
      
      // If previousConsensusFlags isn't initialized yet, then we're
      // clearly not converged, and can skip the remaining tests.
      if(consensusFlags.size() == previousConsensusFlags.size()) {

        // Does it look like we're in a cycle of adding/subtracting the
        // same points?
        if(consensusSetSize <= previousConsensusSetSize) {
          ++strikes;
        }
        if(strikes > 10) {
          return true;
        }

        // Hmm.  Are we selecting the same set as last time?  If so,
        // we're converged.
        if(std::equal(consensusFlags.begin(), consensusFlags.end(),
                      previousConsensusFlags.begin())) {
          return true;
        }
      }
      
      // Finally, do the bookkeeping so that previousConsensusFlags
      // gets updated.
      previousConsensusFlags = consensusFlags;
      previousConsensusSetSize = consensusSetSize;

      return false;
    }
    
  } // namespace computerVision
  
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_RANSACCLASS_HH */
