/**
***************************************************************************
* @file brick/computerVision/ransacClassInterface_impl.hh
*
* Header file declaring an implementation of the RANSAC algorithm.
*
* Copyright (C) 2008-2014 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_RANSACCLASSINTERFACE_IMPL_HH
#define BRICK_COMPUTERVISION_RANSACCLASSINTERFACE_IMPL_HH

// This file is included by ransacClassInterface.hh, and should not be
// directly included by user code, so no need to include
// ransacClassInterface.hh here.
//
// #include <brick/computerVision/ransacClassInterface.hh>

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



    // Calculate which input samples are consistent with the
    // specified model, and return them to the calling context as a
    // SampleSequenceType instance.
    template <class Problem>
    typename Problem::SampleSequenceType
    Ransac<Problem>::
    getConsensusSet(ResultType model)
    {
      // Identify the consensus set, made up of samples that are
      // sufficiently consistent with the model estimate.
      std::vector<bool> consensusFlags(m_problem.getPoolSize());
      this->computeConsensusSet(model, consensusFlags);

      // Select those elements corresponding to the flags we just
      // computed.
      return this->m_problem.getSubset(consensusFlags.begin(),
                                       consensusFlags.end());
    }


    // This member function runs the RANSAC algorithm and returns
    // the computed model.
    template <class Problem>
    typename Ransac<Problem>::ResultType
    Ransac<Problem>::
    getResult() {
      // xxx
      ResultType result;
      this->estimate(result);
      return result;
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

#endif /* #ifndef BRICK_COMPUTERVISION_RANSACCLASSINTERFACE_IMPL_HH */
