/**
***************************************************************************
* @file brick/numeric/normalizedCorrelator.hh
*
* Header file declaring NormalizedCorrelator class.
*
* Copyright (C) 1999-2007,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_NORMALIZEDCORRELATOR_HH
#define BRICK_NUMERIC_NORMALIZEDCORRELATOR_HH

#include <cmath>
#include <deque>
#include <brick/numeric/array2D.hh>

namespace brick {

  namespace numeric {


    /**
     ** This class implements 1D normalized correlation, which is
     ** sometimes also called the Correlation Coefficient.  That is,
     ** given two signals, f[x] and g[x], it computes C(f, g), where
     **
     **   C(f, g) = sum_x(p[x] * q[x]),
     **
     ** and
     **
     **   p[x] = (f[x] - mean(f)) / variance(f)
     **
     ** and
     **
     **   q[x] = (g[x] - mean(g)) / variance(g)
     **
     ** For increased efficiency, the actual computation doesn't
     ** explicitly represent the mean and variance of the input
     ** signals.
     **
     ** NormalizedCorrelator is designed to support efficient
     ** incremental updates of the signals being correlated, and makes
     ** it easy to add and remove samples from the correlation
     ** on-the-fly.
     **
     ** Here are some examples of how you might use this class.
     **
     ** @code
     **
     **   std::vector<double> signalF;
     **   std::vector<double> signalG;
     **
     **   // Real code would set the contents of signalF and signalG
     **   // here.  They must have the same size.
     **
     **   // Here, we compute the normalized correlation of the two
     **   // signals.
     **   NormalizedCorrelator<double> correlator0(
     **     signalF.begin(), signalF.end(), signalG.begin());
     **   double correlation0 = correlator0.getNormalizedCorrelation();
     **
     **   // Here, we repeat the computation using a slightly
     **   // different sequence of calls.
     **   NormalizedCorrelator<double> correlator1;
     **   correlator1.addSamples(
     **     signalF.begin(), signalF.end(), signalG.begin());
     **   double correlation1 = correlator1.getNormalizedCorrelation();
     **
     **   // These lines compute the normalized correlation of the
     **   // first 100 elements of the two signals.
     **   NormalizedCorrelator<double> correlator2;
     **   correlator2.addSamples(
     **     signalF.begin(), signalF.begin() + 100, signalG.begin());
     **   double correlation2 = correlator2.getNormalizedCorrelation();
     **
     **   // These lines illustrate incremental update of the signals
     **   // being correlated.  They compute the normalized
     **   // correlation of elements 2 through 101 (note zero-based
     **   // indexing) of the two signals, but use the result of the
     **   // previous computation to save time.
     **   correlator2.removeSamples(
     **     signalF.begin(), signalF.begin() + 2, signalG.begin());
     **   correlator2.addSamples(
     **     signalF.begin() + 100, signalF.begin() + 102,
     **     signalG.begin() + 100);
     **   double correlation3 = correlator1.getNormalizedCorrelation();
     **
     **   // These lines repeat the computation of correlation2 and
     **   // correlation3, but request that the NormalizedCorrelator
     **   // instance remember the order in which samples were added
     **   // so that it can automatically remove the oldest samples.
     **   // The call to removeOldestSamples() below has the same
     **   // effect as the call to removeSamples() above, but in most
     **   // situations requires less bookkeeping in the calling
     **   // context.
     **   NormalizedCorrelator<double> correlator4(true);
     **   correlator4.addSamples(
     **     signalF.begin(), signalF.begin() + 100, signalG.begin());
     **   double correlation4 = correlator4.getNormalizedCorrelation();
     **
     **   correlator4.removeOldestSamples(2);
     **   correlator4.addSamples(
     **     signalF.begin() + 100, signalF.begin() + 102,
     **     signalG.begin() + 100);
     **   double correlation5 = correlator1.getNormalizedCorrelation();
     **
     ** @endcode
     **
     **/
    template <class Type>
    class NormalizedCorrelator
    {
    public:

      /**
       * This constructor initializes the NormalizedCorrelator
       * instance, but doesn't add any samples.
       *
       * @param trackInput This argument indicates whether the
       * NormalizedCorrelator instance should keep a record of samples
       * as they're added so that it can automatically remove them in
       * order using member function removeSamples(size_t).
       */
      NormalizedCorrelator(bool trackInput = false);


      /**
       * This constructor initializes the NormalizedCorrelator
       * instance using sequences of samples from the two signals to
       * be correlated.  After calling this constructor, the
       * normalized correlation of the two input signals is available
       * via member function getNormalizedCorrelation().  If the
       * addSamples() method is called after calling this constructor,
       * the effect will be as if the input sequences from the
       * constructor and addSamples() calls were simply concatenated.
       *
       * @param begin0 This argument is an iterator pointing to the
       * beginning of the sequence of samples from the first of the
       * two signals to be correlated.
       *
       * @param end0 This argument is an iterator pointing to the
       * end of the sequence of samples from the first of the two
       * signals to be correlated.  Just as with standard library
       * algorithms, the final element of the input sequence is the
       * one _before_ *end0.
       *
       * @param begin1 This argument is an iterator pointing to the
       * beginning of the sequence of samples from the second of the
       * two signals to be correlated.
       *
       * @param trackInput This argument indicates whether the
       * NormalizedCorrelator instance should keep a record of samples
       * as they're added so that it can automatically remove them in
       * order using member function removeSamples(size_t).
       */
      template <class IterType0, class IterType1>
      NormalizedCorrelator(IterType0 begin0,
                           IterType0 end0,
                           IterType1 begin1,
                           bool trackInput = false);


      /**
       * The destructor destroys the NormalizedCorrelator instance and
       * cleans up any associated resources.
       */
      ~NormalizedCorrelator();


      /**
       * This member function adds a single pair of samples (one from
       * each of the two signals to be correlated) to the normalized
       * correlation calculation.  You might call this repeatedly,
       * once for each of many different samples.  If input tracking
       * is enabled (via constructor argument) then this function will
       * update the internal record keeping so that pairs of samples
       * can be automatically removed by a call to member function
       * removeOldestSample() and removeOldestSamples().
       *
       * @param sample0 This argument is the sample value from the
       * first of the two input signals.
       *
       * @param sample1 This argument  is the sample value from the
       * second of the two input signals.
       */
      inline void
      addSample(Type sample0, Type sample1);


      /**
       * This member function works identically to addSample(), with
       * the exception that input tracking is never updated.  Use this
       * member function instead of addSample() if you know that you
       * will never have input tracking enabled, you can't pass your
       * samples in en masse using addSamples(), and you're in such a
       * hurry that the run-time cost of one conditional branch is
       * worth avoiding.
       *
       * @param sample0 This argument is the sample value from the
       * first of the two input signals.
       *
       * @param sample1 This argument  is the sample value from the
       * second of the two input signals.
       */
      inline void
      addSampleWithoutTracking(Type sample0, Type sample1);


      /**
       * This member function adds a sequence of pairs of samples
       * (each pair containing one sample from each of the two signals
       * to be correlated) to the normalized correlation calculation.
       * If input tracking is enabled (via constructor argument) then
       * this function will update the internal record keeping so that
       * pairs of samples can be automatically removed by a call to
       * member function removeOldestSamples().  Note that for the
       * purposes of this automatic removal, *begin0 is considered to
       * be added before *(begin0 + 1).
       *
       * @param begin0 This argument is an iterator pointing to the
       * beginning of the sequence of samples from the first of the
       * two signals to be correlated.
       *
       * @param end0 This argument is an iterator pointing to the
       * end of the sequence of samples from the first of the two
       * signals to be correlated.  Just as with standard library
       * algorithms, the final element of the input sequence is the
       * one _before_ *end0.
       *
       * @param begin1 This argument is an iterator pointing to the
       * beginning of the sequence of samples from the second of the
       * two signals to be correlated.
       */
      template <class IterType0, class IterType1>
      void
      addSamples(IterType0 begin0, IterType0 end0, IterType1 begin1);


      /**
       * This member function removes all samples from the
       * NormalizedCorrelator instance.
       */
      void
      clear();


      /**
       * This member function enables (or disables) internal
       * recordkeeping that allows samples to be automatically removed
       * from the normalized correlation calculation following the
       * order in which they were added.  When input tracking is
       * disabled, execution and storage requirements are reduced.
       * When input tracking is enabled, the NormalizedCorrelator
       * instance maintains a record of all added samples and the
       * order in which they were added.  Subsequent calls to
       * removeOldestSamples() will discard the oldest samples.  This
       * is useful if you need to maintain a running normalized
       * correlation of, for example, the last 100 sample pairs: each
       * time you get a new set of samples, you can add them using
       * addSample() or addSamples(), discard the corresponding number
       * of sample pairs using removeOldestSamples(), and recover the
       * updated normalized correlation by calling
       * getNormalizedCorrelation().
       *
       * Note that even if input tracking is disabled, you can still
       * remove samples explicitly by calling removeSample() or
       * removeSamples().
       *
       * Calling enableInputTracking(true) when tracking is already
       * enabled has undefined result.
       *
       * Currently, you cannot call removeSample() or removeSamples()
       * on a NormalizedCorrelator instance for which input tracking
       * is enabled (we expect this to change), and you cannot call
       * removeOldestSamples() on a NormalizedCorrelator instance for
       * which input tracking is not enabled.
       *
       * @param trackInput Setting this argument to true enables input
       * tracking. Setting this argument false disabled input
       * tracking.
       */
      void
      enableInputTracking(bool trackInput = true);


      /**
       * This member function returns the number of sample pairs
       * contributing to the normalized correlation.  It indicates the
       * total number of sample pairs added by calls to the
       * constructor, addSample(), and addSamples(), less the number
       * of sample pairs removed by removeSample(), removeSamples(),
       * and removeOldestSamples().
       *
       * @return The return value is number of sample pairs.
       */
      inline size_t
      getCount() const;


      /**
       * This member function returns the normalized correlation of
       * all the currently added sample pairs.  If no sample pairs
       * have been added, the return value is 1.0.
       *
       * @return The return value is the computed normalized
       * correlation.
       */
      Type
      getNormalizedCorrelation() const;


      /**
       * This member function returns a bool indicating whether or not
       * input tracking is enabled (see member function
       * enableInputTracking()).
       *
       * @return The return value is true if input tracking is
       * enabled, false otherwise.
       */
      inline bool
      isInputTrackingEnabled() const {return m_inputTracker0Ptr != 0;}


      /**
       * If input tracking is enabled, this member function removes
       * pairs of samples from the normalized correlation calculation,
       * following the order in which they were added.  It has the
       * same effect as the removeSamples() member function, but does
       * not require the calling context to explicitly specify the
       * sample values to be removed.  Calling this member function
       * when input tracking is disabled is an error.
       *
       * @param count This argument specifies how many sample pairs to
       * remove.
       */
      void
      removeOldestSamples(size_t count);


      /**
       * This member function removes a pair of sample values from the
       * normalized correlation calculation.  Note that this function
       * works regardless of whether input tracking is enabled.
       *
       * @param sample0 This argument is the sample value from the
       * first of the two input signals.
       *
       * @param sample1 This argument  is the sample value from the
       * second of the two input signals.
       */
      void
      removeSample(Type sample0, Type sample1);


      /**
       * This member function works identically to removeSample(),
       * with the exception that input tracking is never updated.  Use
       * this member function instead of removeSample() if you know
       * that you will never have input tracking enabled, you can't
       * pass your samples in en masse using removeSamples(), and
       * you're in such a hurry that the run-time cost of one
       * conditional branch is worth avoiding.
       *
       * @param sample0 This argument is the sample value from the
       * first of the two input signals.
       *
       * @param sample1 This argument  is the sample value from the
       * second of the two input signals.
       */
      inline void
      removeSampleWithoutTracking(Type sample0, Type sample1);


      /**
       * This member function removes a sequence of pairs of samples
       * (each pair containing one sample from each of the two signals
       * to be correlated) from the normalized correlation
       * calculation.  Note that this function works regardless of
       * whether input tracking is enabled (see member function
       * enableInputTracking()).
       *
       * @param begin0 This argument is an iterator pointing to the
       * beginning of the sequence of samples from the first of the
       * two signals.
       *
       * @param end0 This argument is an iterator pointing to the end
       * of the sequence of samples from the first of the two signals.
       * Just as with standard library algorithms, the final element
       * of the input sequence is the one _before_ *end0.
       *
       * @param begin1 This argument is an iterator pointing to the
       * beginning of the sequence of samples from the second of the
       * two signals.
       */
      template <class IterType0, class IterType1>
      void
      removeSamples(IterType0 begin0, IterType0 end0, IterType1 begin1);

    private:

      size_t m_count;
      std::deque<Type>* m_inputTracker0Ptr;
      std::deque<Type>* m_inputTracker1Ptr;
      Type m_sum0;
      Type m_sum1;
      Type m_sum00;
      Type m_sum01;
      Type m_sum11;

    };

  } // namespace numeric

} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/normalizedCorrelator_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_NORMALIZEDCORRELATOR_HH */
