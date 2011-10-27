/**
***************************************************************************
* @file brick/computerVision/randomSampleSelector_impl.hh
*
* Source file defining a class template for randomly sampling
* populations of things.
*
* Copyright (C) 2008, 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_RANDOMSAMPLESELECTOR_IMPL_HH
#define BRICK_COMPUTERVISION_RANDOMSAMPLESELECTOR_IMPL_HH

// This file is included by randomSampleSelector.hh, and should not be
// directly included by user code, so no need to include
// randomSampleSelector.hh here.
// 
// #include <brick/numeric/randomSampleSelector.hh>

namespace brick {

  namespace computerVision {

    // The constructor specifies the full population of samples from
    // which to randomly select.
    template <class Sample>
    template <class IterType>
    RandomSampleSelector<Sample>::
    RandomSampleSelector(IterType beginIter, IterType endIter)
      : m_sampleVector(beginIter, endIter)
    {
      // Empty.
    }


    // This member function returns a SampleSequenceType instance
    // containing the entire population passed to the constructor.
    template <class Sample>
    typename RandomSampleSelector<Sample>::SampleSequenceType
    RandomSampleSelector<Sample>::
    getPool()
    {
      return std::make_pair(m_sampleVector.begin(), m_sampleVector.end());
    }


    // This member function returns a the number of samples in the
    // entire population passed to the constructor.
    template <class Sample>
    size_t
    RandomSampleSelector<Sample>::
    getPoolSize()
    {
      return m_sampleVector.size();
    }


    // This member function returns a SampleSequenceType instance
    // drawn randomly (without replacement) from the sample
    // population.
    template <class Sample>
    typename RandomSampleSelector<Sample>::SampleSequenceType
    RandomSampleSelector<Sample>::
    getRandomSample(size_t sampleSize)
    {
      for(size_t ii = 0; ii < sampleSize; ++ii) {
        int jj = m_pseudoRandom.uniformInt(ii, m_sampleVector.size());
        std::swap(m_sampleVector[ii], m_sampleVector[jj]);
      }
      return std::make_pair(
        m_sampleVector.begin(), m_sampleVector.begin() + sampleSize);
    }


    // This member function has very limited usefulness.
    template <class Sample>
    template<class IterType>
    typename RandomSampleSelector<Sample>::SampleSequenceType
    RandomSampleSelector<Sample>::
    getSubset(IterType beginIter, IterType endIter)
    {
      typedef typename std::vector<SampleType>::iterator SampleIter;

      SampleIter poolIter = m_sampleVector.begin();
      SampleIter candidateIter = m_sampleVector.begin();
      while(beginIter != endIter && candidateIter != m_sampleVector.end()) {
        if(*beginIter) {
          if(candidateIter != poolIter) {
            std::swap(*poolIter, *candidateIter);
          }
          ++poolIter;
        }
        ++beginIter;
        ++candidateIter;
      }
      return std::make_pair(m_sampleVector.begin(), poolIter);
    }
        
  } // namespace computerVision
  
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_RANDOMSAMPLESELECTOR_IMPL_HH */
