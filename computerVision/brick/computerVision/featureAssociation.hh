/**
***************************************************************************
* @file brick/computerVision/featureAssociation.hh
*
* Header file declaring functions that pair features from two disjoint
* sets based on a user supplied similarity measure.
*
* Copyright (C) 2008,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_FEATUREASSOCIATION_HH
#define BRICK_COMPUTERVISION_FEATUREASSOCIATION_HH

#include <brick/numeric/array2D.hh>

namespace brick {

  namespace computerVision {

    /**
     * This function template implements the feature association
     * algorithm of Guy Scott and H. Christopher Longuet-Higgins, as
     * described in [1].
     *
     * Template argument Functor must implement a similarity measure
     * for comparing pairs of features.  For example, calling
     *
     *   similarityFunctor(*sequence0Begin, *sequence1Begin)
     *
     * Template argument FloatType describes the return value of this
     * functor, as well as the type used to represent real numbers
     * within the function.  The functor should return a similarity
     * measure between 0.0 and 1.0.  Scott's suggested similarity
     * measure is:
     *
     *   g = exp(-r_ij^2 / 2*sigma^2),
     *
     * where g is the similarity result, r_ij is the distance between
     * point i in the first feature set and point j in the second
     * feature set, and sigma is a parameter setting how sensitive the
     * algorithm is to increasing distance between matched features
     * (large alpha gives low sensitivity, small alpha gives high
     * sensitivity).
     *
     * Note that this routine does not employ robust statistics, and
     * requires decomposing (via SVD) an N by M matrix, where N and M
     * are the lengths of the two input sequenceds.  This limits its
     * usefulness for large feature sets.
     *
     * [1] G. L. Scott and H. C. Longuet Higgins, "An Algorithm for
     * Associating the Features of Two Images," Proceedings of
     * Biological Sciences, Vol. 244, No. 1309, pp. 21-26, April,
     * 1991.
     *
     *
     * @param sequence0Begin This argument is an STL style iterator
     * pointing to the first element of the first feature sequence.
     *
     * @param sequence0End This argument is an STL style iterator
     * pointing one-past-the-last element of the first feature
     * sequence.
     *
     * @param sequence1Begin This argument is an STL style iterator
     * pointing to the first element of the second feature sequence.
     *
     * @param sequence1End This argument is an STL style iterator
     * pointing one-past-the-last element of the second feature
     * sequence.
     *
     * @param similarityFunctor This argument specifies a functor that
     * takes two features and computes their similarity.
     *
     * @return The return value is a vector of pairs of indices in
     * which each pair of indices indicates a pair of corresponding
     * features.  The first element of the pair is an index into
     * features sequence 0, while the second is an index into feature
     * sequence 1.  The pairs will be arranged so that the sequence 0
     * indices are in ascending order.
     *
     * Note(xxx): Provide a way of returning similarity values.
     */
    template<class FloatType, class Iterator0, class Iterator1, class Functor>
    std::vector< std::pair<size_t, size_t> >
    associateFeaturesScott91(Iterator0 sequence0Begin, Iterator0 sequence0End,
                             Iterator1 sequence1Begin, Iterator1 sequence1End,
                             Functor similarityFunctor);

  } // namespace computerVision

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/featureAssociation_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_FEATUREASSOCIATION_HH */
