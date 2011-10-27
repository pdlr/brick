/**
***************************************************************************
* @file dlrComputerVision/ransac.h
*
* Header file defining helper functions for implementing Fischler's
* and Bolles's RANSAC algorithm.
*
* Copyright (C) 2008, 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
***************************************************************************
*/

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <dlrCommon/exception.h>
#include <dlrComputerVision/randomSampleSelector.h>
#include <dlrNumeric/maxRecorder.h>

namespace dlr {

  namespace computerVision {

    template <class InIter, class OutIter, Functor>
    void
    ransacGetConsensusSet(
      InIter inBegin, InIter inEnd, OutIter outBegin, Functor functor)
    {
      // while(inBegin != ineEnd) {
      //   if(functor(*inBegin)) {
      //     *outBegin = *inBegin;
      //   }
      // }
      std::copy_if(inBegin, inEnd, outBegin, functor);
    }

    
    template <class Type, class Functor>
    brick::numeric::Array2D<Type>
    ransacGetConsensusSetRows(
      brick::numeric::Array2D<Type> const& candidates,
      Functor& functor)
    {
      brick::numeric::Array1D<bool> indicatorArray(candidates.rows());
      for(unsigned int ii = 0; ii < candidates.rows(); ++ii) {
        if(functor(candidates[ii])) {
          indicatorArray[ii] = true;
          ++count;
        } else {
          indicatorArray[ii] = false;
        }
      }

      brick::numeric::Array2D<Type> result(count, candidates.columns());
      unsigned int outputRow = 0;
      for(unsigned int ii = 0; ii < candidates.rows(); ++ii) {
        if(indicatorArray[ii]) {
          result.getRow(outputRow).copy(candidates.getRow(ii));
        }
      }
      return result;
    }

    
    template <class InIterator, class OutIterator, 
    unsigned int
    ransacGetConsensusSetByThreshold(
    unsigned int
    ransacGetRequiredIterations(unsigned int sampleSize,
                                double requiredConfidence,
                                double inlierProbability)
    {
      if((requiredConfidence < 0.0) || (requiredConfidence >= 1.0)) {
        DLR_THROW(common::ValueException, "ransacGetRequiredIterations()",
                  "Probability value requiredConfidence is out of range.");
      }
      if((inlierProbability < 0.0) || (inlierProbability >= 1.0)) {
        DLR_THROW(common::ValueException, "ransacGetRequiredIterations()",
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
      double numberOfRandomSampleSets =
        std::ceil(std::log(1.0 - requiredConfidence)
                  / std::log(singlePickDisconfidence)) + 0.5;

      return numberOfRandomSampleSets;
    }
    

    template <class Type>
    brick::numeric::Array2D<Type>
    ransacSelectRows(brick::numeric::Array1D<Type> const& sampleArray,
                     unsigned int numberOfSamplesRequired)
    {
      brick::random::PseudoRandom pseudoRandom();
      return ransacSelectRows(sampleArray, numberOfSamplesRequired,
                              pseudoRandom);
    }
    

    template <class Type>
    brick::numeric::Array2D<Type>
    ransacSelectRows(brick::numeric::Array2D<Type> const& sampleArray,
                     unsigned int numberOfSamplesRequired,
                     brick::common::Int64& seed)
    {
      brick::random::PseudoRandom pseudoRandom(seed);
      brick::numeric::Array1D<Type> result =
        ransacSelectRows(sampleArray, numberOfSamplesRequired, pseudoRandom);
      seed = pseudoRandom.getCurrentSeed();
      return result;
    }
    
    
    template <class Type>
    brick::numeric::Array2D<Type>
    ransacSelectRows(brick::numeric::Array2D<Type> const& sampleArray,
                     unsigned int numberOfSamplesRequired,
                     brick::random::PseudoRandom& pseudoRandom)
    {
      brick::numeric::Array1D<Type> result(numberOfSamplesRequired,
                                           sampleArray.columns());
      brick::numeric::Array2D<Type> shuffleBuffer(sampleArray.rows(),
                                                  sampleArray.columns());
      brick::numeric::Array1D<bool> indicators(sampleArray.rows());
      indicators = false;

      // Select each sample in turn.  We will sample without replacement.
      for(unsigned int ii = 0; ii < numberOfSamplesRequired; ++ii) {
        // Easy enough: choose from among the remaining samples.
        unsigned int selectedRowIndex = pseudoRandom.getUniformInt(
          ii, sampleArray.rows());
        // Copy from the input array, unless this row has already been
        // selected.
        if(indicators[selectedRowIndex]) {
          // This row _has_ already been selected.  Fortunately, last
          // time it was selected, we had the forsight to copy an
          // un-selected row into shufflebuffer.  Use that one
          // instead.
          result.getRow(ii).copy(shuffleBuffer.getRow(selectedRowIndex));
        } else {
          // This row _has not_ already been selected.  Carry on.
          result.getRow(ii).copy(sampleArray.getRow(selectedRowIndex));
        }

        // The row that was just selected might get selected again
        // later, so we remove it from circulation by copying a
        // not-selected row into its place in shuffleBuffer.  We do
        // this unless selectedRow == ii, in which case the lower
        // bound of the getUniformInt() call above prevents
        // reselection of the row, so we do nothing.
        if(selectedRow != ii) {
          shuffleBuffer.getRow(selectedRowIndex).copy(sampleArray.getRow(ii));
          indicators[selectedRowIndex] = true;
        }
      }
      return result;
    }
    
  } // namespace computerVision
  
} // namespace dlr

