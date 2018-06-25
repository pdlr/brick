/**
***************************************************************************
* @file brick/computerVision/threePointAlgorithm.hh
*
* Header file declaring the threePointAlgorithm() function template.
*
* Copyright (C) 2009,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_THREEPOINTALGORITHM_HH
#define BRICK_COMPUTERVISION_THREEPOINTALGORITHM_HH

#include <brick/computerVision/cameraIntrinsicsPinhole.hh>
#include <brick/numeric/transform3D.hh>
#include <brick/numeric/vector2D.hh>
#include <brick/numeric/vector3D.hh>
#include <brick/random/pseudoRandom.hh>

namespace brick {

  namespace computerVision {

    /**
     * This function implements the "three point perspective pose
     * estimation algorithm" of Grunert[1][2] for recovering the
     * camera-frame coordinates of the corners of a triangle of known
     * size, given the projections of those corners in the camera
     * image.  We follow the derivation in [2].  That is, given w_0,
     * w_1, and w_2, all 3D positions in world coordinates, u_0, u_1,
     * u_2, the projections of those coordinates in the image, and
     * pinhole camera parameters, the algorithm recovers p_0, p_1, and
     * p_2, the positions of those points in camera coordinates.
     *
     * [1] J. A. Grunert, "Das Pothenotische Problem in Erweiterter
     * Gestalt Nebst Uber Seine Anwendungen in der Geodisie," Grunerts
     * Archiv fur Mathematik und Physik, Band 1, 1841, pp. 238-248.
     *
     * [2] R. M. Haralick, C. Lee, K. Ottenberg, and M. Nolle, "Review
     * and Analysis of Solutions of the Three Point Perspective Pose
     * Estimation Problem," International Journal of Computer Vision,
     * 13, 3, 331-356 (1994).
     *
     * @param w0 This argument is the first of the three 3D points in
     * world coordinates.
     *
     * @param w1 This argument is the second of the three 3D points in
     * world coordinates.
     *
     * @param w2 This argument is the third  of the three 3D points in
     * world coordinates.
     *
     * @param u0 This argument is the image location of the projection
     * of world point w0.
     *
     * @param u1 This argument is the image location of the projection
     * of world point w1.
     *
     * @param u2 This argument is the image location of the projection
     * of world point w2.
     *
     * @param intrinsics This argument describes the intrinsic
     * calibration of the camera that generated the image from which
     * u0, u1, u2 were drawn.
     *
     * @param p0OutputIter This argument is used to return estimates
     * of the point in camera coordinates corresponding to w0.  It is
     * an iterator pointing to the beginning of an output sequence of
     * Vector3D<FloatType>, and must be able to accept at least four
     * values.
     *
     * @param p1OutputIter This argument is used to return estimates
     * of the point in camera coordinates corresponding to w1.  It is
     * an iterator pointing to the beginning of an output sequence of
     * Vector3D<FloatType>, and must be able to accept at least four
     * values.
     *
     * @param p2OutputIter This argument is used to return estimates
     * of the point in camera coordinates corresponding to w2.  It is
     * an iterator pointing to the beginning of an output sequence of
     * Vector3D<FloatType>, and must be able to accept at least four
     * values.
     *
     * @param epsilon This argument sets some internal tolerances of
     * the algorithm, and should be left at its default value for now.
     *
     * @return The return value indicates how many solutions for the
     * camera position and orientation were found.  If returnValue >=
     * 1, then the first elements of the three output sequences
     * correspond to the first solution found.  If return value >= 2,
     * then the second elements of the three output sequences
     * correspond to the second solution found, and so on.
     */
    template <class FloatType, class IterType>
    unsigned int
    threePointAlgorithm(brick::numeric::Vector3D<FloatType> const& w0,
                        brick::numeric::Vector3D<FloatType> const& w1,
                        brick::numeric::Vector3D<FloatType> const& w2,
                        brick::numeric::Vector2D<FloatType> const& u0,
                        brick::numeric::Vector2D<FloatType> const& u1,
                        brick::numeric::Vector2D<FloatType> const& u2,
                        CameraIntrinsicsPinhole<FloatType> const& intrinsics,
                        IterType p0OutputIter,
                        IterType p1OutputIter,
                        IterType p2OutputIter,
                        FloatType epsilon = 1.0E-8);

    /**
     * This function implements the "robust" version of
     * threePointAlgorithm().  Multiple solutions for camera pose are
     * computed using randomly selected sets of three input points, an
     * error value is computed (based on all of the input points) for
     * each of the potential solutions, and the solution having the
     * best error value is retained.
     *
     * Template argument InIter3D is an iterator type describing a
     * sequence of points in 3D space.  Template argument InIter2D is
     * an iterator type describing the corresponding sequence of 2D
     * points.
     *
     * @param worldPointsBegin This argument is the beginning (in the
     * STL sense) of a sequence of 3D points expressed in world
     * coordinates, and represented as brick::numeric::Vector3D<FloatType>
     * instances.
     *
     * @param worldPointsEnd This argument is the end (in the STL
     * sense) of the sequence begun by worldPointsBegin.
     *
     * @param imagePointsBegin This argument is the beginning (in the
     * STL sense) of a sequence of 2D points corresponding to the
     * elements of [worlPointsBegin, worldPointsEnd], and expressed in
     * image coordinates.
     *
     * @param intrinsics This argument describes the intrinsic
     * calibration of the camera that generated the input points.
     *
     * @param iterations This argument specifies how many random
     * samples of three input points should be processed to generate
     * solution hypotheses.
     *
     * @param inlierProportion This argument specifies what proportion
     * of the input points are expected to be "inliers" and conform to
     * the correct solution (once we find it).  It is used to tune the
     * error value computation.
     *
     * @param score This argument is a projection residual indicating
     * the goodness of the final solution.
     *
     * @param pRandom This argument is a pseudorandom number generator
     * used by the algorithm to select sets of three input points.
     *
     * @return The return value is a coordinate tranformation that
     * takes points in world coordinates and converts them to camera
     * coordinates.
     */
    template <class FloatType, class InIter3D, class InIter2D>
    brick::numeric::Transform3D<FloatType>
    threePointAlgorithmRobust(
      InIter3D worldPointsBegin,
      InIter3D worldPointsEnd,
      InIter2D imagePointsBegin,
      CameraIntrinsicsPinhole<FloatType> const& intrinsics,
      size_t iterations,
      FloatType inlierProportion,
      FloatType& score,
      brick::random::PseudoRandom& pRandom = brick::random::PseudoRandom());


    /**
     * This sort-of-private function solves a quartic function that
     * shows up repeatedly in the threePointAlgorithm() problem.
     * Template argument OutIter is an iterator type used to pass the
     * multiple results back to the calling context.
     *
     * Documentation for the actual equation goes here.
     *
     * @param cosAlpha This argument is one of six used to pass the
     * constant coefficients of the quartic equation.
     *
     * @param cosBeta This argument is one of six used to pass the
     * constant coefficients of the quartic equation.
     *
     * @param cosGamma This argument is one of six used to pass the
     * constant coefficients of the quartic equation.
     *
     * @param a2 This argument is one of six used to pass the constant
     * coefficients of the quartic equation.
     *
     * @param b2 This argument is one of six used to pass the constant
     * coefficients of the quartic equation.
     *
     * @param c2 This argument is one of six used to pass the constant
     * coefficients of the quartic equation.
     *
     * @param epsilon This argument specifies the threshold at which
     * quantities should be considered "almost zero" for numerical
     * precision.
     *
     * @param s0Iter This argument passes results back to the calling
     * context.
     *
     * @param s1Iter This argument passes results back to the calling
     * context.
     *
     * @param s2Iter This argument passes results back to the calling
     * context.
     *
     * @param condition This argument returns a number indicating the
     * stability of the result.  Lower is better.
     *
     * @return The return value indicates how many potential solutions
     * were found.
     */
    template <class FloatType, class OutIter>
    unsigned int
    solveThreePointAlgorithmQuarticSystem(
      FloatType cosAlpha, FloatType cosBeta, FloatType cosGamma,
      FloatType a2, FloatType b2, FloatType c2, FloatType epsilon,
      OutIter s0Iter, OutIter s1Iter, OutIter s2Iter,
      FloatType& condition);

  } // namespace computerVision

} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/threePointAlgorithm_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_THREEPOINTALGORITHM_HH */
