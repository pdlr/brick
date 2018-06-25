/**
***************************************************************************
* @file brick/computerVision/nChooseKSampleSelector_impl.hh
*
* Header file defining inline and template functions declared in
* nChooseKSampleSelector.hh.
*
* Copyright (C) 2008,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_NCHOOSEKSAMPLESELECTOR_IMPL_HH
#define BRICK_COMPUTERVISION_NCHOOSEKSAMPLESELECTOR_IMPL_HH

// This file is included by nChooseKSampleSelector.hh, and should not
// be directly included by user code, so no need to include
// nChooseKSampleSelector.hh here.
//
// #include <brick/computerVision/nChooseKSampleSelector.hh>

#include <limits>
#include <brick/common/exception.hh>

namespace brick {

  namespace computerVision {

    template <class Sample>
    template <class IterType>
    NChooseKSampleSelector<Sample>::
    NChooseKSampleSelector(size_t sampleSize,
                           IterType beginIter,
                           IterType endIter)
      : m_numberOfSamples(0),
        m_poolVector(beginIter, endIter),
        m_previousSampleNumber(std::numeric_limits<size_t>::max() - 1),
        m_sampleIndices(sampleSize),
        m_sampleVector(sampleSize)
    {
      if(m_sampleIndices.size() == 0) {
        BRICK_THROW(brick::common::ValueException,
                  "NChooseKSampleSelector::NChooseKSampleSelector()",
                  "Argument sampleSize must be nonzero.");
      }
      if(m_sampleIndices.size() > m_poolVector.size()) {
        BRICK_THROW(brick::common::ValueException,
                  "NChooseKSampleSelector::NChooseKSampleSelector()",
                  "Argument sampleSize must be less than or equal to "
                  "the number of elements in the input sequence.");
      }
    }


    template <class Sample>
    size_t
    NChooseKSampleSelector<Sample>::
    getNumberOfSamples()
    {
      if(m_numberOfSamples == 0) {
        // Compute binomial coefficient "n choose k".
        size_t nFactorial = this->getFactorial(m_poolVector.size());
        size_t kFactorial = this->getFactorial(m_sampleIndices.size());
        size_t nMinusKFactorial = this->getFactorial(
          m_poolVector.size() - m_sampleIndices.size());
        m_numberOfSamples = nFactorial / (kFactorial * nMinusKFactorial);
      }
      return m_numberOfSamples;
    }


    template <class Sample>
    typename NChooseKSampleSelector<Sample>::SampleSequenceType
    NChooseKSampleSelector<Sample>::
    getSample(size_t sampleNumber)
    {
      if(sampleNumber != m_previousSampleNumber + 1) {
        m_previousSampleNumber = sampleNumber;
        return this->statelessGetSample(sampleNumber);
      }

      // This code ought to do something more efficient than this.
      this->incrementSampleIndices();
      for(size_t kk = 0; kk < m_sampleIndices.size(); ++kk) {
        m_sampleVector[kk] = m_poolVector[m_sampleIndices[kk]];
      }
      ++m_previousSampleNumber;
      return std::make_pair(m_sampleVector.begin(), m_sampleVector.end());
    }


    template <class Sample>
    size_t
    NChooseKSampleSelector<Sample>::
    getFactorial(size_t argument)
    {
      size_t result = 1;
      for(size_t ii = 2; ii <= argument; ++ii) {
        result *= ii;
      }
      return result;
    }


    template <class Sample>
    void
    NChooseKSampleSelector<Sample>::
    incrementSampleIndices()
    {
      size_t const finalIndex = m_sampleIndices.size() - 1;

      // Try for the naive increment, ignoring overflow.
      ++(m_sampleIndices[finalIndex]);

      // Now propagate any overflow up the chain.
      size_t ii = finalIndex;
      while(m_sampleIndices[finalIndex] >= m_poolVector.size()) {
        // Remember size_t wraps around: (0 - 1) is a very big number.
        if((--ii) > m_sampleIndices.size()) {
          BRICK_THROW(brick::common::IndexException,
                    "NChooseKSampleSelector::incrementSampleIndices()",
                    "No more samples available.");
        }
        ++(m_sampleIndices[ii]);
        for(size_t jj = ii + 1; jj < m_sampleIndices.size(); ++jj) {
          m_sampleIndices[jj] = m_sampleIndices[jj - 1] + 1;
        }
      }
    }


    template <class Sample>
    typename NChooseKSampleSelector<Sample>::SampleSequenceType
    NChooseKSampleSelector<Sample>::
    statelessGetSample(size_t sampleNumber)
    {
      for(size_t ii = 0; ii < m_sampleIndices.size(); ++ii) {
        m_sampleIndices[ii] = ii;
      }
      for(size_t jj = 0; jj < sampleNumber; ++jj) {
        this->incrementSampleIndices();
      }
      for(size_t kk = 0; kk < m_sampleIndices.size(); ++kk) {
        m_sampleVector[kk] = m_poolVector[m_sampleIndices[kk]];
      }
      return std::make_pair(m_sampleVector.begin(), m_sampleVector.end());
    }

  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_NCHOOSEKSAMPLESELECTOR_IMPL_HH */
