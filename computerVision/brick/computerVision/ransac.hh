/**
***************************************************************************
* @file brick/computerVision/ransac.hh
*
* Header file declaring helper functions for implementing Fischler's and
* Bolles's RANSAC algorithm.
*
* Copyright (C) 2008, 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_RANSAC_HH
#define BRICK_COMPUTERVISION_RANSAC_HH

#include <brick/numeric/array1D.hh>
#include <brick/numeric/array2D.hh>
#include <brick/random/pseudoRandom.hh>

namespace brick {

  namespace computerVision {

    /** 
     * Selects only those elements of the input sequence that, when
     * passed as arguments to functor.operator()(), result in a true
     * return value.
     * 
     * @param inBegin This argument is an STL-style iterator for the
     * beginning of the input sequence.
     * 
     * @param inEnd This argument is an STL-style iterator for the
     * end of the input sequence.
     * 
     * @param outBegin This argument is an iterator for the beginning
     * of the output sequence.  Elements of the input sequence that
     * satisfy functor will be copied to the output sequence.
     * 
     * @param functor This argument will determine which elements are
     * copied from the input sequence to the output sequence.  Each
     * element of the input sequence will be passed as the argument of
     * a call to this functor.  Those arguments for which the functor
     * returns true will be copied to the output sequence.
     *
     * @return The return value is the number of input elements copied.
     */
    template <class InIter, class OutIter, class Functor>
    unsigned int
    ransacGetConsensusSet(
      InIter inBegin, InIter inEnd, OutIter outBegin, Functor functor);

    
    /** 
     * This is a convenience function that functions just like
     * ransacGetConsensusSet, except that the output of the functor
     * argument is passed to a second functor for evaluation.  For
     * example, this call:
     *
     * @code
     *   ransacGetConsensusSetByComparison(
     *     candidates.begin(), candidates.end(),
     *     std::back_inserter(consensusSet),
     *     myFunctorThatReturnsDouble,
     *     std::bind2nd(std::less<double>(), myThreshold));
     * @endcode
     *
     * ...will select only those elements of the input sequence for
     * which myFunctorThatReturnsDouble() returns a value less than
     * myThreshold.   This is equivalent to
     * 
     * @code
     *   ransacGetConsensusSet(
     *     candidates.begin(), candidates.end(),
     *     std::back_inserter(consensusSet),
     *     std::compose1(std::bind2nd(std::less<double>(), myThreshold),
     *                   myFunctorThatReturnsDouble));
     * @endcode
     *
     * ...except that it works for functors that do not implement the
     * std::unary_function interface.
     * 
     * @param inBegin This argument is an STL-style iterator for the
     * beginning of the input sequence.
     * 
     * @param inEnd This argument is an STL-style iterator for the
     * end of the input sequence.
     * 
     * @param outBegin This argument is an iterator for the beginning
     * of the output sequence.  Elements of the input sequence that
     * satisfy functor will be copied to the output sequence.
     * 
     * @param functor This argument, and the next, determine which
     * elements are copied from the input sequence to the output
     * sequence.  Each element of the input sequence will be passed as
     * the argument of a call to this functor.  The result of this
     * functor will be passed as an argument to functor "criterion,"
     * and those arguments for which criterion returns true will be
     * copied to the output sequence.
     *
     * @param criterion Please see the documentation for argument
     * functor.
     *
     * @return The return value is the number of input elements copied.
     */
    template <class InIter, class OutIter, class Functor, class Criterion>
    unsigned int
    ransacGetConsensusSetByComparison(
      InIter inBegin, InIter inEnd, OutIter outBegin, Functor functor,
      Criterion criterion);


    /** 
     * This functions just like ransacGetConsensusSet, except that the
     * input and output sequences are replaced by 2D arrays.  Each row
     * of the the input array is passed to functor as an Array1D<Type>
     * instance, and the output Array2D will contain only those rows
     * for which functor returns true.  You might use it like this:
     *
     * @code
     *   struct MyFunctor {
     *     bool testSample(Array1D<double> const& candidate);
     *   };
     *   
     *   Array2D<double> candidates(numberOfCandidates,
     *                              numberOfElementsInACandidate);
     *   for(unsigned int ii = 0; ii < numberOfCandidates; ++ii) {
     *     // Copy candidate into row ii of candidates array.
     *   }
     *
     *   Array2D<double> selectedRows = ransacGetConsensusSetRows(
     *     candidates, MyFunctor());
     * @endcode
     *
     * @param candidates This argument is a 2D array in which each row
     * represents one candidate to be evaluated for inclusion in the
     * consensus set.
     * 
     * @param functor This argument is a functor that accepts Array1D
     * arguments and returns a bool, indicating whether or not the
     * argument should be included in the consensus set.
     * 
     * @return The return value is an Array2D instance containing only
     * those rows for which functor returned true.
     */
    template <class Type, class Functor>
    brick::numeric::Array2D<Type>
    ransacGetConsensusSetRows(
      brick::numeric::Array2D<Type> const& candidates,
      Functor functor);

    
    /** 
     * This is is a convenience function that functions just like
     * ransacGetConsensusSetRows, except that the output of the
     * functor argument is passed to a second functor for evaluation.
     * For example, this call:
     *
     * @code
     *   ransacGetConsensusSetRowsByComparison(
     *     candidates, MyFunctorThatReturnsDouble(),
     *     std::bind2nd(std::less<double>(), myThreshold));
     * @endcode
     *
     * ...will select only those rows of the input array for which
     * myFunctorThatReturnsDouble() returns a value less than
     * myThreshold.  This is equivalent to
     *
     * @code
     *   ransacGetConsensusSetRows(
     *     candidates, 
     *     std::compose1(std::bind2nd(std::less<double>(), myThreshold),
     *                   myFunctorThatReturnsDouble));
     * @endcode
     *
     * ...except that it works for functors that do not implement the
     * std::unary_function interface.
     * 
     * @param candidates This argument is a 2D array in which each row
     * represents one candidate to be evaluated for inclusion in the
     * consensus set.
     * 
     * @param functor This argument, and the next, determine which
     * rows are copied from the input array to the output array.  Each
     * row of the input sequence will be passed as the argument of
     * a call to this functor.  The result of this functor will be
     * passed as an argument to functor "criterion," and those
     * arguments for which criterion returns true will be copied to
     * the output sequence.
     *
     * @param criterion Please see the documentation for argument
     * functor.
     *
     * @return The return value is an Array2D instance containing only
     * those rows for which functor returned true.
     */
    template <class Type, class Functor, class Criterion>
    brick::numeric::Array2D<Type>
    ransacGetConsensusSetRowsByComparison(
      brick::numeric::Array2D<Type> const& candidates,
      Functor functor,
      Criterion criterion);


    /** 
     * Computes how many RANSAC iterations are required to be sure
     * there's at least one iteration in which the model was estimated
     * using all inliers.
     * 
     * @param sampleSize This argument indicates how many observations
     * are required to estimate a model.
     * 
     * @param requiredConfidence indicates the how confident we must
     * be that at least one iteration will generate a good model.  The
     * larger this number is, the more iterations will be required.
     * Setting this to 1.0 (or higher) is an error, because (as long
     * as there is at least one outlier in the input set) there's
     * always a chance that every iteration will by bad luck include
     * that outlier.
     * 
     * @param inlierProbability This argument indicates what
     * proportion of the input sample population are believed to be
     * inliers.  A value of 0.0 means no inliers, a value of 1.0 means
     * 100% inliers.
     * 
     * @return The return value is the number 
     */
    unsigned int
    ransacGetRequiredIterations(unsigned int sampleSize,
                                double requiredConfidence,
                                double inlierProbability);

    
    /** 
     * Randomly (or rather, pseudo-randomly) selects elements of the
     * input sequence for use in RANSAC estimation.
     * 
     * @param sampleArray This argument is an vector in which each
     * element represents one observation, or sample, for use in the
     * RANSAC algorithm.
     * 
     * @param numberOfSamplesRequired This argument indicates how many
     * elements must be selected and copied.
     * 
     * @return The return value is a vector containing only those
     * elements that were selected.
     */
    template <class Type>
    std::vector<Type>
    ransacSelectElements(std::vector<Type> const& inputSequence,
                         unsigned int numberOfSamplesRequired);
    

    /** 
     * Randomly (or rather, pseudo-randomly) selects rows of the input
     * array for use in RANSAC estimation.
     * 
     * @param sampleArray This argument is an array in which each row
     * represents one observation, or sample, for use in the RANSAC
     * algorithm.
     * 
     * @param numberOfSamplesRequired This argument indicates how many
     * rows must be selected and copied.
     * 
     * @return The return value is an Array2D instance containing only
     * those rows that were selected.  It will always have
     * numberOfSamplesRequired rows, and the same number of columns as
     * sampleArray.
     */
    template <class Type>
    brick::numeric::Array2D<Type>
    ransacSelectRows(brick::numeric::Array2D<Type> const& sampleArray,
                     unsigned int numberOfSamplesRequired);
    

    /** 
     * This function is just like ransacSelectRows, except that a
     * third argument allows the calling context to control the
     * pseudo-random sequence.
     * 
     * @param sampleArray This argument is an array in which each row
     * represents one observation, or sample, for use in the RANSAC
     * algorithm.
     * 
     * @param numberOfSamplesRequired This argument indicates how many
     * rows must be selected and copied.
     * 
     * @param seed This argument passes in a seed that will be used to
     * initialize the pseudo-random number generator.  On exit, this
     * argument will be updated to contain a new seed, which can be
     * passed to the next call of ransacSelectRows() so that the
     * psuedoRandom sequence is repeatable.
     * 
     * @return The return value is an Array2D instance containing only
     * those rows that were selected.  It will always have
     * numberOfSamplesRequired rows, and the same number of columns as
     * sampleArray.
     */
    template <class Type>
    brick::numeric::Array2D<Type>
    ransacSelectRows(brick::numeric::Array2D<Type> const& sampleArray,
                     unsigned int numberOfSamplesRequired,
                     brick::common::Int64& seed);

    
    /** 
     * This function is just like ransacSelectRows, except that a
     * third argument allows the calling context to control the
     * pseudo-random sequence.
     * 
     * @param sampleArray This argument is an array in which each row
     * represents one observation, or sample, for use in the RANSAC
     * algorithm.
     * 
     * @param numberOfSamplesRequired This argument indicates how many
     * rows must be selected and copied.
     * 
     * @param pseudoRandom This argument passes in a pseudo-random
     * number generator that will be used to select rows for inclusion
     * in the return array.
     * 
     * @return The return value is an Array2D instance containing only
     * those rows that were selected.  It will always have
     * numberOfSamplesRequired rows, and the same number of columns as
     * sampleArray.
     */
    template <class Type>
    brick::numeric::Array2D<Type>
    ransacSelectRows(brick::numeric::Array2D<Type> const& sampleArray,
                     unsigned int numberOfSamplesRequired,
                     brick::random::PseudoRandom& pseudoRandom);


    
  } // namespace computerVision
  
} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/ransac_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_RANSAC_HH */
