/**
***************************************************************************
* @file brick/computerVision/eightPointAlgorithm.hh
*
* Header file declaring the eightPointAlgorithm() function template.
*
* Copyright (C) 2008,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_EIGHTPOINTALGORITHM_HH
#define BRICK_COMPUTERVISION_EIGHTPOINTALGORITHM_HH

#include <brick/numeric/array1D.hh>
#include <brick/numeric/array2D.hh>
#include <brick/numeric/vector2D.hh>

namespace brick {

  namespace computerVision {

    /**
     * This function implements the "eight point algorithm"[1] for
     * recovering the fundamental matrix of a pair of cameras from a
     * sequence of at least eight pairs of corresponding image points.
     * That is, it recovers the matrix F, such that
     *
     *   transpose(u') * F * u = 0
     *
     * where u is a homogeneous 2D point in the first image, u' is a
     * homogeneous 2D point in the second image, and the symbol "*"
     * indicates matrix multiplication.
     *
     * This implementation implements the input transformation
     * described in [2] to improve the numerical stability of the
     * result.
     *
     * [1] Longuet-Higgins, H.C., "A Computer Algorithm for
     * Reconstructing a Scene From Two Projections," Nature, vol. 293,
     * pp. 133­135, Sept 1981.
     *
     * [2] Hartley, R. I., "In Defense of the Eight Point Algorithm."
     * IEEE Transactions on Pattern Analysis and Machine Intelligence,
     * vol. 19, No. 6, pp. 580-593, June 1997.
     *
     * @param sequence0Begin This argument is the beginning (in the
     * STL sense) of a sequence of feature points, represented as
     * brick::numeric::Vector2D<FloatType> instances, from the first image (the u
     * points in the equation above).
     *
     * @param sequence0End This argument is the end (in the STL sense)
     * of a sequence of feature points, represented as
     * brick::numeric::Vector2D<FloatType> instances, from the first image (the u
     * points in the equation above).
     *
     * @param sequence1Begin This argument is the beginning (in the
     * STL sense) of a sequence of feature points, represented as
     * brick::numeric::Vector2D<FloatType> instances, from the second image (the u'
     * points in the equation above).
     *
     * @return The return value is the recovered fundamental matrix.
     */
    template<class FloatType, class Iterator>
    brick::numeric::Array2D<FloatType>
    eightPointAlgorithm(Iterator sequence0Begin, Iterator sequence0End,
                        Iterator sequence1Begin);


    /**
     * WARNING: This function may go away at some point, or be
     * replaced with a slightly different interface.
     *
     * This function is just like the other version of
     * eightPointAlgorithm(), except that it returns by reference
     *
     * @param sequence0Begin This argument matches the corresponding
     * argument of the other version of eightPointAlgorithm().
     *
     * @param sequence0End This argument matches the corresponding
     * argument of the other version of eightPointAlgorithm().
     *
     * @param sequence1Begin This argument matches the corresponding
     * argument of the other version of eightPointAlgorithm().
     *
     * @param eigenvalues This argument returns by reference a vector
     * of eigenvalues used in computing the result of the function
     * call.  Only the last one should be close to zero.
     *
     * @return The return value is the recovered fundamental matrix.
     */
    template<class FloatType, class Iterator>
    brick::numeric::Array2D<FloatType>
    eightPointAlgorithm(Iterator sequence0Begin, Iterator sequence0End,
                        Iterator sequence1Begin,
                        brick::numeric::Array1D<FloatType>& eigenvalues);



    /**
     * This function is used internally by eightPointAlgorithm() to
     * translate and scale input points so that their mean lies at the
     * origin and they have isotropic unit variance.  It is exposed
     * here to facilitate testing, and in case it's useful.
     *
     * @param inputPoints This argument is an Nx3 Array2D<FloatType>
     * instance in which each row` represents one 2D point of the
     * input set using homogeneous coordinates.  The last column of
     * this array will normally contain only ones (although this is
     * not required to be true).`
     *
     * @param outputPoints This argument is an Nx3 Array2D<FloatType>
     * instance used to return the transformed points to the calling
     * context.  The elements of the third column of this array will
     * almost certainly not be equal to 1.0.
     *
     * @param transform this argument is return the affine transform
     * that takes points from inputPoints to outputPoints via left
     * multiplication.
     */
    template<class FloatType>
    void
    normalizePointSequence(
      brick::numeric::Array2D<FloatType> const& inputPoints,
      brick::numeric::Array2D<FloatType>& outputPoints,
      brick::numeric::Array2D<FloatType>& transform);

  } // namespace computerVision

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/eightPointAlgorithm_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_EIGHTPOINTALGORITHM_HH */
