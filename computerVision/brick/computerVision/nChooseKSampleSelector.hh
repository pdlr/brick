/**
***************************************************************************
* @file brick/computerVision/nChooseKSampleSelector.hh
*
* Header file declaring a class for exhaustively sampling populations of
* things.
*
* Copyright (C) 2008,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_NCHOOSEKSAMPLESELECTOR_HH
#define BRICK_COMPUTERVISION_NCHOOSEKSAMPLESELECTOR_HH

#include <vector>

namespace brick {

  namespace computerVision {

    /**
     ** This class template provides capabilities to exhaustively select
     ** sequences of samples from a pool of candidates.  It is useful
     ** for implementing robust statistics algorithms such as RANSAC.
     **
     ** Template argument Sample specifies what type of thing will
     ** make up the population from which samples will be drawn.
     **/
    template <class Sample>
    class NChooseKSampleSelector {
    public:

      // ========= Public typedefs. =========

      /**
       ** This typedef simply mirrors template argument Sample. 
       **/
      typedef Sample SampleType;

      /**
       ** This typedef specifies what type will be used to represent
       ** sequences of samples.  See getSample(size_t) for more
       ** information.
       **/
      typedef std::pair<typename std::vector<SampleType>::const_iterator,
                        typename std::vector<SampleType>::const_iterator>
      SampleSequenceType;


      // ========= Public member functions. =========

      /** 
       * The constructor specifies the full population of samples from
       * which to exhaustively select.
       *
       * @param sampleSize This argument specifies the size of the
       * samples that will be drawn from the input sequence.  The
       * number of elements in the sequence [beginIter, endIter] must
       * be at least this large.
       *
       * @param beginIter This argument and the next specify a
       * sequence from which to copy the sample population.
       * 
       * @param endIter This argument and the previous specify a
       * sequence from which to copy the sample population.
       */
      template <class IterType>
      NChooseKSampleSelector(size_t sampleSize,
                             IterType beginIter,
                             IterType endIter);


      /** 
       * This member function returns the number of distinct
       * combinations of sampleSize elements that can be drawn from
       * the sequence passed to the constructor.
       * 
       * @return The return value is "N choose M," where N is the size
       * of the sequence passed to the constructor, and M is the value
       * of constructor argument sampleSize.
       */
      size_t
      getNumberOfSamples();
      

      /** 
       * This member function returns a SampleSequenceType instance
       * containing the entire population passed to the constructor.
       * 
       * @return The return value is a sequence containing the entire
       * population.
       */
      SampleSequenceType
      getPool() {
        return std::make_pair(m_poolVector.begin(), m_poolVector.end());
      }


      /** 
       * This member function returns a the number of samples in the
       * entire population passed to the constructor.
       * 
       * @return The return value is the size of the entire
       * population.
       */
      size_t
      getPoolSize() {return m_poolVector.size();}


      /** 
       * This member function returns a SampleSequenceType instance
       * drawn from the sample population.  This sequence will remain
       * valid at least until the next call to getSample().
       * 
       * @param sampleNumber This argument specifies which of the
       * available samples should be returned.  Its value should be in
       * the range 0 <= sampleNumber < this->getNumberOfSamples().
       * This function is more efficient if the value of sampleNumber
       * starts at 0 and is incremented with each cal to getSample().
       * 
       * @return The return value is a sequence containing the
       * requested sample.
       */
      SampleSequenceType
      getSample(size_t sampleNumber);

    private:

      size_t
      getFactorial(size_t argument);

      void
      incrementSampleIndices();
      
      SampleSequenceType
      statelessGetSample(size_t sampleNumber);

      size_t m_numberOfSamples;
      std::vector<SampleType> m_poolVector;
      size_t m_previousSampleNumber;
      std::vector<size_t> m_sampleIndices;
      std::vector<SampleType> m_sampleVector;

    };

  } // namespace computerVision
  
} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/nChooseKSampleSelector_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_NCHOOSEKSAMPLESELECTOR_HH */
