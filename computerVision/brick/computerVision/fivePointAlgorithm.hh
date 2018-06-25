/**
***************************************************************************
* @file brick/computerVision/fivePointAlgorithm.hh
*
* Header file declaring the fivePointAlgorithm() function template.
*
* Copyright (C) 2009,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_FIVEPOINTALGORITHM_HH
#define BRICK_COMPUTERVISION_FIVEPOINTALGORITHM_HH

#include <brick/numeric/array1D.hh>
#include <brick/numeric/array2D.hh>
#include <brick/numeric/transform3D.hh>
#include <brick/numeric/vector2D.hh>
#include <brick/random/pseudoRandom.hh>


namespace brick {

  namespace computerVision {

    /**
     * WARNING(xxx): The essential matrix returned by this function is
     * currently not normalized to reasonable magnitude.  You can
     * easily get back a matrix in which the average element magnitude
     * is on the order of 1.0E12.
     *
     * This function implements the "five point algorithm"[1] for
     * recovering the essential matrix of a pair of cameras from a
     * sequence of at least five pairs of corresponding image points.
     * That is, it recovers the matrix E, such that
     *
     * @code
     *   transpose(q') * E * q = 0
     * @endcode
     *
     * where q is a homogeneous 2D point in the "calibrated coordinate
     * system" of the first image (see below), q' is a homogeneous 2D
     * point in calibrated coordinate system of the second image, and
     * the symbol "*" indicates matrix multiplication.
     *
     * This algorithm differs from the eight point algorithm in that
     * the intrinsic parameters of both cameras must be known.  For
     * pinhole cameras, points in the calibrated coordinate system of
     * the camera are related to points in the pixel coordinate system
     * by the equations
     *
     * @code
     *   q = inverse(K) * u
     *   q' = inverse(K') * u'
     * @endcode
     *
     * where K and K' are 3x3 matrices encoding the camera intrinsic
     * parameters.  Typically, a K matrix will look something like
     * this:
     *
     * @code
     *       |f/k_u,     0, C_u|
     *   K = |    0, f/k_v, C_v|
     *       |    0,     0,   1|
     * @endcode
     *
     * where f is the focal length of the camera, k_u & k_v describe
     * the physical size of camera pixels in the horizontal and
     * vertical directions, respectively, and C_u & C_v are the image
     * coordinates of the camera projection center.
     *
     * [1] Henrik Stewénius, Christopher Engels, and David Nistér,
     * "Recent Developments on Direct Relative Orientation."  ISPRS
     * Journal of Photogrammetry and Remote Sensing, vol. 60, no. 4,
     * pp. 284-294, January 2006.
     *
     * @param sequence0Begin This argument is the beginning (in the
     * STL sense) of a sequence of calibrated feature points,
     * represented as brick::numeric::Vector2D instances, from the first
     * image (the q points in the equation above).  You might generate
     * this sequence of points by taking raw image coordinates and
     * then left-multiplying them by inverse(K).
     *
     * @param sequence0End This argument is the end (in the STL sense)
     * of a sequence begun by sequence0Begin.
     *
     * @param sequence1Begin This argument is the beginning (in the STL
     * sense) of a sequence of calibrated feature points, represented
     * as brick::numeric::Vector2D instances, from the second image (the
     * q' points in the equation above).  You might generate this
     * sequence of points by taking raw image coordinates and then
     * left-multiplying them by inverse(K').
     *
     * @return The return value is the recovered essential matrix.
     */
    template<class FloatType, class Iterator>
    std::vector< brick::numeric::Array2D<FloatType> >
    fivePointAlgorithm(Iterator sequence0Begin, Iterator sequence0End,
                       Iterator sequence1Begin);


    /**
     * WARNING(xxx): The essential matrix returned by this funtion is
     * currently not normalized to reasonable magnitude.  You can
     * easily get back a matrix in which the average element magnitude
     * is on the order of 1.0E12.
     *
     * Warning: this interface may change.
     *
     * This function implements the two-view robust five-point
     * algorithm described in section 5 of [2].  Note that our method
     * of computing the goodness of a potential solution for the
     * essential matrix involves testing feature point pairs against
     * the (image-space) epipolar constraint.  We do not project
     * features into 3D space and compute errors there.
     *
     * [2] David Nister, "An Efficient Solution to the Five-Point
     * Relative Pose Problem."  IEEE Transactions on Pattern Analysis
     * and Machine Intelligence (PAMI), vol. 26, no. 6, pp. 756-770,
     * June 2004.
     *
     * @param sequence0Begin This argument is the beginning (in the
     * STL sense) of a sequence of calibrated feature points,
     * represented as brick::numeric::Vector2D instances, from the first
     * image (the q points in the equation above).  You might generate
     * this sequence of points by taking raw image coordinates and
     * then left-multiplying them by inverse(K).
     *
     * @param sequence0End This argument is the end (in the STL sense)
     * of a sequence begun by sequence0Begin.
     *
     * @param sequence1Begin This argument is the beginning (in the
     * STL sense) of a sequence of calibrated feature points,
     * represented as brick::numeric::Vector2D instances, from the
     * second image (the q' points in the equation above).  You might
     * generate this sequence of points by taking raw image
     * coordinates and then left-multiplying them by inverse(K').
     *
     * @param iterations This argument specifies how many random
     * samples of five points each should be used to generate
     * hypothetical essential matrices.  Setting this number high
     * increases the chances that you'll find the right answer,
     * however runtime is directly proportional to the number of
     * iterations.
     *
     * @param inlierProportion This argument specifies what proportion
     * of the input point pairs you expect to be correctly matched.
     * Set this a little low, so that you're confident there are at
     * least inlierProportion * (sequence0End - sequence0Begin)
     * correct feature matches in the input data.
     *
     * @param score This argument returns a value indicating how good
     * the returned essential matrix is.  0.0 is perfect.  1.0 is
     * terrible.
     *
     * @return The return value is the recovered "best-fit" essential
     * matrix.
     */
    template<class FloatType, class Iterator>
    brick::numeric::Array2D<FloatType>
    fivePointAlgorithmRobust(Iterator sequence0Begin, Iterator sequence0End,
                             Iterator sequence1Begin,
                             size_t iterations,
                             FloatType inlierProportion,
                             FloatType& score,
                             brick::random::PseudoRandom pRandom
                             = brick::random::PseudoRandom());


    /**
     * Warning: this interface may change.
     *
     * This function implements the three-view robust five-point
     * algorithm described in section 5 of [2].
     *
     * @param sequence0Begin This argument is the beginning (in the
     * STL sense) of a sequence of calibrated feature points,
     * represented as brick::numeric::Vector2D instances, from the first
     * image (the q points in the equation above).  You might generate
     * this sequence of points by taking raw image coordinates and
     * then left-multiplying them by inverse(K).
     *
     * @param sequence0End This argument is the end (in the STL sense)
     * of a sequence begun by sequence0Begin.
     *
     * @param sequence1Begin This argument is the beginning (in the
     * STL sense) of a sequence of calibrated feature points,
     * represented as brick::numeric::Vector2D instances, from the
     * second image (the q' points in the equation above).  You might
     * generate this sequence of points by taking raw image
     * coordinates and then left-multiplying them by inverse(K').
     *
     * @param sequence2Begin This argument is the beginning (in the STL
     * sense) of a sequence of calibrated feature points, represented
     * as brick::numeric::Vector2D instances, from the third image.
     *
     * @param iterations This argument specifies how many random
     * samples of fiv points each should be used to generate
     * hypothetical essential matrices.  Setting this number high
     * increases the chances that you'll find the right answer,
     * however runtime is directly proportional to the number of
     * iterations.
     *
     * @param inlierProportion This argument specifies what proportion
     * of the input point pairs you expect to be correctly matched.
     * Set this a little low, so that you're confident there are at
     * least inlierProportion * (sequence0End - sequence0Begin)
     * correct feature matches in the input data.
     *
     * @param cam2Ecam0 This argument returns by reference the
     * recovered Essential matrix between the third image and the
     * first image, so that transpose(q'') * cam2Ecam0 * q = 0.
     *
     * @param cam0Tcam2 This argument returns by reference a recovered
     * coordinate transformation taking 3D points from the coordinate
     * frame of third camera and converting them to the coordinate
     * frame of first camera.
     *
     * @param cam1Tcam2 This argument returns by reference a recovered
     * coordinate transformation taking 3D points from the coordinate
     * frame of third camera and converting them to the coordinate
     * frame of second camera.
     *
     * @param score This argument returns a value indicating how good
     * the returned parameters are.  0.0 is perfect.  1.0 is
     * terrible.
     *
     * @param pRandom This argument is a pseudorandom number generator
     * used by the algorithm to select sets of three input points.
     */
    template<class FloatType, class Iterator>
    void
    fivePointAlgorithmRobust(Iterator sequence0Begin, Iterator sequence0End,
                             Iterator sequence1Begin,
                             Iterator sequence2Begin,
                             size_t iterations,
                             FloatType inlierProportion,
                             brick::numeric::Array2D<FloatType>& cam2Ecam0,
                             brick::numeric::Transform3D<FloatType>& cam0Tcam2,
                             brick::numeric::Transform3D<FloatType>& cam1Tcam2,
                             FloatType& score,
                             brick::random::PseudoRandom pRandom
                             = brick::random::PseudoRandom());


    // Return value is residual in pix^2.
    template <class FloatType>
    FloatType
    checkEpipolarConstraint(
      brick::numeric::Array2D<FloatType> const& fundamentalMx,
      brick::numeric::Vector2D<FloatType>& point0,
      brick::numeric::Vector2D<FloatType>& point1);


    // WARNING(xxx): We have observed at least one case in which the
    // return value from this function is fishy.  We currently do not
    // trust it.
    //
    // Warning Think of this function as being "private."  Gets cam1Tcam0.
    // Might throw if test points are on parallel rays.
    template <class FloatType>
    brick::numeric::Transform3D<FloatType>
    getCameraMotionFromEssentialMatrix(
      brick::numeric::Array2D<FloatType> const& EE,
      brick::numeric::Vector2D<FloatType> const& testPointCamera0,
      brick::numeric::Vector2D<FloatType> const& testPointCamera1);


    // Returns point in coordinate system of camera 0.  Note that
    // argument is the inverse of what's returned by
    // getCameraMotionFromEssentialMatrix().  Might throw if test
    // points are on parallel rays.
    template <class FloatType>
    brick::numeric::Vector3D<FloatType>
    triangulateCalibratedImagePoint(
      brick::numeric::Transform3D<FloatType> const& c0Tc1,
      brick::numeric::Vector2D<FloatType> const& testPointCamera0,
      brick::numeric::Vector2D<FloatType> const& testPointCamera1);


    /**
     * This function is used internally by fivePointAlgorithm() to
     * generate a 10x20 matrix of coefficients of polynomial
     * constraints.  It will normally not be useful unless you happen
     * to be implementing a similar five point algorithm solver (such
     * as the one described in [2], which is faster, but slightly less
     * accurate.
     *
     * @param E0Array This argument is the first basis element of the
     * null space of the system of linear constraints on the essential
     * matrix.
     *
     * @param E1Array This argument is the first basis element of the
     * null space of the system of linear constraints on the essential
     * matrix.
     *
     * @param E2Array This argument is the first basis element of the
     * null space of the system of linear constraints on the essential
     * matrix.
     *
     * @param E3Array This argument is the first basis element of the
     * null space of the system of linear constraints on the essential
     * matrix.
     *
     * @return The return value is the coefficient matrix, M,
     * described in [1].
     */
    template <class FloatType>
    brick::numeric::Array2D<FloatType>
    generateFivePointConstraintMatrix(
      brick::numeric::Array2D<FloatType> const& E0Array,
      brick::numeric::Array2D<FloatType> const& E1Array,
      brick::numeric::Array2D<FloatType> const& E2Array,
      brick::numeric::Array2D<FloatType> const& E3Array);

  } // namespace computerVision

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/fivePointAlgorithm_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_FIVEPOINTALGORITHM_HH */
