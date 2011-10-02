/**
***************************************************************************
* @file brick/computerVision/randomSampleSelector.hh
*
* Header file declaring a class template for randomly sampling
* populations of things.
*
* Copyright (C) 2008, 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_RANDOMSAMPLESELECTOR_HH
#define BRICK_COMPUTERVISION_RANDOMSAMPLESELECTOR_HH

#include <vector>
#include <brick/random/pseudoRandom.h>

namespace brick {

  namespace computerVision {

    /**
     ** This class template provides capabilities to randomly select
     ** sequences of samples from a pool of candidates.  It is useful
     ** for implementing robust statistics algorithms such as RANSAC.
     **
     ** Template argument Sample specifies what type of thing will
     ** make up the population from which random samples will be
     ** drawn.
     **/
    template <class Sample>
    class RandomSampleSelector {
    public:

      // ========= Public typedefs. =========

      /**
       ** This typedef simply mirrors template argument Sample. 
       **/
      typedef Sample SampleType;

      /**
       ** This typedef specifies what type will be used to represent
       ** sequences of samples.  See getRandomSample(size_t) for more
       ** information.
       **/
      typedef std::pair<typename std::vector<SampleType>::const_iterator,
                        typename std::vector<SampleType>::const_iterator>
        SampleSequenceType;


      // ========= Public member functions. =========

      /** 
       * The constructor specifies the full population of samples from
       * which to randomly select.
       * 
       * @param beginIter This argument and the next specify a
       * sequence from which to copy the sample population.
       * 
       * @param endIter This argument and the previous specify a
       * sequence from which to copy the sample population.
       */
      template <class IterType>
      RandomSampleSelector(IterType beginIter, IterType endIter);


      /** 
       * This member function returns a SampleSequenceType instance
       * containing the entire population passed to the constructor.
       * This sequence will remain valid at least until the next call
       * to getRandomSample() or getSubset().
       * 
       * @return The return value is a sequence containing the entire
       * population.
       */
      SampleSequenceType
      getPool();


      /** 
       * This member function returns a the number of samples in the
       * entire population passed to the constructor.
       * 
       * @return The return value is the size of the entire
       * population.
       */
      size_t
      getPoolSize();


      /** 
       * This member function returns a SampleSequenceType instance
       * drawn randomly (without replacement) from the sample
       * population.  This sequence will remain valid at least until
       * the next call to getRandomSample() or getSubset().
       * 
       * @param sampleSize This argument specifies how many elements
       * should be in the returned sequence.
       * 
       * @return The return value is a sequence containing the the
       * requested number of randomly selected samples.
       */
      SampleSequenceType
      getRandomSample(size_t sampleSize);


      /** 
       * This member function has very limited usefulness.  Following
       * a call to getPool(), you can compute a sequence of bools (or
       * values that will implicitly cast to bools) indicating which
       * samples you like, pass this sequence to getSubset(), and get
       * back a SampleSequenceType instance containing just the
       * samples you requested.  In the current implementation, the
       * (internal) ordering of the sample population changes with
       * each call to getRandomSample() or getSubset(), so if there's
       * been a call to either getRandomSample() or getSubset() more
       * recently than your last call to getPool(), then this function
       * will return unexpected results.
       * 
       * @param beginIter This argument and the next are the
       * pair of indicators specifying which elements should be in
       * the output sequence.  Elements for which the indicator is
       * true will be included in the sequence, while elements for
       * which the indicator is false will not.
       * 
       * @param endIter This argument and the previous are the
       * pair of indicators specifying which elements should be in
       * the output sequence.
       * 
       * @return The return value is a sequence containing the
       * requested samples.
       */
      template<class IterType>
      SampleSequenceType
      getSubset(IterType beginIter, IterType endIter);
      
    private:
      
      brick::random::PseudoRandom m_pseudoRandom;
      std::vector<SampleType> m_sampleVector;

    };

  } // namespace computerVision
  
} // namespace brick

#include <brick/computerVision/randomSampleSelector_impl.hh>

/* ============ Definitions of inline & template functions ============ */


namespace brick {

  namespace computerVision {

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

#endif /* #ifndef BRICK_COMPUTERVISION_RANDOMSAMPLESELECTOR_HH */
