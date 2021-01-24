/**
***************************************************************************
* @file brick/numeric/normalizedCorrelator.hh
*
* Implementation of inline and template functions for the
* NormalizedCorrelator class.
*
* Copyright (C) 1999-2007,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_NORMALIZEDCORRELATOR_IMPL_HH
#define BRICK_NUMERIC_NORMALIZEDCORRELATOR_IMPL_HH

// This file is included by normalizedCorrelator.hh, and should not be
// directly included by user code, so no need to include
// normalizedCorrelator.hh here.
//
// #include <brick/numeric/normalizedCorrelator.hh>

namespace brick {

  namespace numeric {

    // This constructor initializes the NormalizedCorrelator
    // instance, but doesn't add any samples.
    template <class Type>
    NormalizedCorrelator<Type>::
    NormalizedCorrelator(bool trackInput)
      : m_count(0),
        m_inputTracker0Ptr(0),
        m_inputTracker1Ptr(0),
        m_sum0(static_cast<Type>(0)),
        m_sum1(static_cast<Type>(0)),
        m_sum00(static_cast<Type>(0)),
        m_sum01(static_cast<Type>(0)),
        m_sum11(static_cast<Type>(0))
    {
      if(trackInput) {
        this->enableInputTracking();
      }
    }


    // This constructor initializes the NormalizedCorrelator
    // instance using sequences of samples from the two signals to
    // be correlated.
    template <class Type>
    template <class IterType0, class IterType1>
    NormalizedCorrelator<Type>::
    NormalizedCorrelator(IterType0 begin0, IterType0 end0, IterType1 begin1,
                         bool trackInput)
      : m_count(0),
        m_inputTracker0Ptr(0),
        m_inputTracker1Ptr(0),
        m_sum0(static_cast<Type>(0)),
        m_sum1(static_cast<Type>(0)),
        m_sum00(static_cast<Type>(0)),
        m_sum01(static_cast<Type>(0)),
        m_sum11(static_cast<Type>(0))
    {
      if(trackInput) {
        this->enableInputTracking();
      }
      this->addSamples(begin0, end0, begin1);
    }


    // The destructor destroys the NormalizedCorrelator instance and
    // cleans up any associated resources.
    template <class Type>
    NormalizedCorrelator<Type>::
    ~NormalizedCorrelator()
    {
      this->enableInputTracking(false);
    }


    // This member function adds a single pair of samples (one from
    // each of the two signals to be correlated) to the normalized
    // correlation calculation.
    template <class Type>
    void
    NormalizedCorrelator<Type>::
    addSample(Type sample0, Type sample1)
    {
      // Copy input, if required to do so.
      if(this->isInputTrackingEnabled()) {
        m_inputTracker0Ptr->push_back(sample0);
        m_inputTracker1Ptr->push_back(sample1);
      }
      this->addSampleWithoutTracking(sample0, sample1);
    }


    // This member function works identically to addSample(), with
    // the exception that input tracking is never updated.
    template <class Type>
    inline void
    NormalizedCorrelator<Type>::
    addSampleWithoutTracking(Type sample0, Type sample1)
    {
      m_sum0 += sample0;
      m_sum1 += sample1;
      m_sum00 += sample0 * sample0;
      m_sum01 += sample0 * sample1;
      m_sum11 += sample1 * sample1;
      ++m_count;
    }


    // This member function adds a sequence of pairs of samples
    // (each pair containing one sample from each of the two signals
    // to be correlated) to the normalized correlation calculation.
    template <class Type>
    template <class IterType0, class IterType1>
    void
    NormalizedCorrelator<Type>::
    addSamples(IterType0 begin0, IterType0 end0, IterType1 begin1)
    {
      // Copy input, if required to do so.
      if(this->isInputTrackingEnabled()) {
        IterType0 begin0Copy = begin0;
        IterType1 begin1Copy = begin1;
        while(begin0Copy != end0) {
          m_inputTracker0Ptr->push_back(*begin0Copy);
          m_inputTracker1Ptr->push_back(*begin1Copy);
          ++begin0Copy;
          ++begin1Copy;
        }
      }

      // Update statistics.
      while(begin0 != end0) {
        this->addSampleWithoutTracking(*begin0, *begin1);
        ++begin0;
        ++begin1;
      }
    }


    // This member function removes all samples from the
    // NormalizedCorrelator instance.
    template <class Type>
    void
    NormalizedCorrelator<Type>::
    clear()
    {
      if(this->isInputTrackingEnabled) {
        // Clear input tracking cache by re-enabling.  The
        // documentation above says this has undefined result, but
        // because we control the implementation, it's ok for us to
        // abuse it.
        this->enableInputTracking(true);
      }
      m_count = 0;
      m_sum0 = static_cast<Type>(0);
      m_sum1 = static_cast<Type>(0);
      m_sum00 = static_cast<Type>(0);
      m_sum01 = static_cast<Type>(0);
      m_sum11 = static_cast<Type>(0);
    }


    // This member function enables (or disables) internal
    // recordkeeping that allows samples to be automatically removed
    // from the normalized correlation calculation following the
    // order in which they were added.
    template <class Type>
    void
    NormalizedCorrelator<Type>::
    enableInputTracking(bool trackInput)
    {
      if(trackInput) {
        if(m_inputTracker0Ptr == 0) {
          m_inputTracker0Ptr = new std::deque<Type>;
        } else {
          m_inputTracker0Ptr->clear();
        }
        if(m_inputTracker1Ptr == 0) {
          m_inputTracker1Ptr = new std::deque<Type>;
        } else {
          m_inputTracker1Ptr->clear();
        }
      } else {
        if(m_inputTracker0Ptr != 0) {
          delete m_inputTracker0Ptr;
          m_inputTracker0Ptr = 0;
        }
        if(m_inputTracker1Ptr != 0) {
          delete m_inputTracker1Ptr;
          m_inputTracker1Ptr = 0;
        }
      }
    }


    // This member function returns the number of sample pairs
    // contributing to the normalized correlation.
    template <class Type>
    inline size_t
    NormalizedCorrelator<Type>::
    getCount() const
    {
      return m_count;
    }


    // This member function returns the normalized correlation of
    // all the currently added sample pairs.
    template <class Type>
    Type
    NormalizedCorrelator<Type>::
    getNormalizedCorrelation() const
    {
      if(m_count == 0) {
        return static_cast<Type>(1);
      }
      Type oneOverN = static_cast<Type>(1) / static_cast<Type>(m_count);
      Type numerator = m_sum01 - oneOverN * m_sum0 * m_sum1;
      Type denominator = std::sqrt((m_sum00 - oneOverN * m_sum0 * m_sum0)
                                   * (m_sum11 - oneOverN * m_sum1 * m_sum1));
      return numerator / denominator;
    }


    // If input tracking is enabled, this member function removes
    // pairs of samples from the normalized correlation calculation,
    // following the order in which they were added.
    template <class Type>
    void
    NormalizedCorrelator<Type>::
    removeOldestSamples(size_t count)
    {
      if(!(this->isInputTrackingEnabled())) {
        BRICK_THROW(brick::common::LogicException,
                    "NormalizedCorrelator::removeInputSamples()",
                    "Attempt to call removeOldestSamples() method of a "
                    "NormalizedCorrelator instance that does not have "
                    "input tracking enabled.");
      }
      if(count > (m_inputTracker0Ptr->size())) {
        BRICK_THROW(brick::common::ValueException,
                    "NormalizedCorrelator::removeInputSamples()",
                    "Trying to remove more samples than have been added.");
      }
      while(count != 0) {
        // Warning(xxx): if this call is changed to removeSample(),
        // then no tests fail, but calls to removeOldestSamples() will
        // throw.
        this->removeSampleWithoutTracking(m_inputTracker0Ptr->front(),
                                          m_inputTracker1Ptr->front());
        m_inputTracker0Ptr->pop_front();
        m_inputTracker1Ptr->pop_front();
        --count;
      }
    }


    // This member function removes a pair of sample values from the
    // normalized correlation calculation.
    template <class Type>
    void
    NormalizedCorrelator<Type>::
    removeSample(Type sample0, Type sample1)
    {
      if(this->isInputTrackingEnabled()) {
        BRICK_THROW(brick::common::NotImplementedException,
                    "NormalizedCorrelator::removeSample()",
                    "Currently, removeSample() can only be called if "
                    "input tracking is not enabled.");
      }
      this->removeSampleWithoutTracking(sample0, sample1);
    }


    // This member function works identically to removeSample(),
    // with the exception that input tracking is never updated.
    template <class Type>
    inline void
    NormalizedCorrelator<Type>::
    removeSampleWithoutTracking(Type sample0, Type sample1)
    {
      m_sum0 -= sample0;
      m_sum1 -= sample1;
      m_sum00 -= sample0 * sample0;
      m_sum01 -= sample0 * sample1;
      m_sum11 -= sample1 * sample1;
      --m_count;
    }


    // This member function removes a sequence of pairs of samples
    // (each pair containing one sample from each of the two signals
    // to be correlated) from the normalized correlation
    // calculation.
    template <class Type>
    template <class IterType0, class IterType1>
    void
    NormalizedCorrelator<Type>::
    removeSamples(IterType0 begin0, IterType0 end0, IterType1 begin1)
    {
      if(this->isInputTrackingEnabled()) {
        BRICK_THROW(brick::common::NotImplementedException,
                    "NormalizedCorrelator::removeSamples()",
                    "Currently, removeSamples() can only be called if "
                    "input tracking is not enabled.");
      }

      while(begin0 != end0) {
        this->removeSampleWithoutTracking(*begin0, *begin1);
        ++begin0;
        ++begin1;
      }
    }


  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_NORMALIZEDCORRELATOR_IMPL_HH */
