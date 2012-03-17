/**
***************************************************************************
* @file brick/computerVision/registerPoints3D.hh
*
* Header file declaring an implementation of Horn's method for
* registering 3D point sets.
*
* Copyright (C) 1998-2007,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_COMPUTERVISION_REGISTERPOINTS3D_HH
#define BRICK_COMPUTERVISION_REGISTERPOINTS3D_HH

#include <brick/numeric/transform3D.hh>

namespace brick {

  namespace computerVision {
  
    /**
     * Using Horn's method [1], find the Rigid body tranform which
     * takes one set of points (represented as
     * brick::numeric::Vector3D<FloatType>) and most nearly registers
     * them with a second set.
     *
     * [1] Horn, B, "Closed-form solution of absolute orientation
     * using unit quaternions", J. Opt. Soc. Am., Vol.  4(4), April,
     * 1987.
     * 
     * @param fromPointsBegin This iterator, along with fromPointsEnd,
     * defines one of the two sets of points to be registered.  The
     * type of (*fromPointsBegin) must be Vector3D<FloatType>.  The
     * returned Transform3D<FloatType> instance will take points in
     * the range specified by this iterator pair and transform them so
     * that they match the corresponding elements of the "to" range as
     * closely as possible in the least-squares sense.
     * 
     * @param fromPointsEnd See the documentation for argument
     * fromPointsBegin.
     * 
     * @param toPointsBegin This argument defines the beginning of one
     * of the two sets of points to be registered.  The returned
     * Transform3D<FloatType> instance will take points in the range
     * specified by arguments fromPointsBegin and fromPointsEnd, and
     * transform them so that they match as closely as possible (in
     * the least squares sense) the corresponding elements of the
     * range starting from toPointsBegin and extending (fromPointsEnd
     * - fromPointsBegin) elements.  possible in the least-squares
     * sense.
     * 
     * @return The return value is a Transform3D<FloatType> instance
     * which takes the points in the "from" range and matches them to
     * the points in the "to" range.
     */
    template <class FloatType, class InIter0, class InIter1>
    brick::numeric::Transform3D<FloatType>
    registerPoints3D(InIter0 fromPointsBegin, InIter0 fromPointsEnd,
                     InIter1 toPointsBegin);


    /**
     * This function works just like the three-argument form of
     * registerPoints3D(), except that only selected points are
     * considered in the registration.
     *
     * @param fromPointsBegin This iterator, along with fromPointsEnd,
     * defines one of the two sets of points to be registered.  The
     * type of (*fromPointsBegin) must be Vector3D<FloatType>.  The
     * returned Transform3D<FloatType> instance will take points in
     * the range specified by this iterator pair and transform them so
     * that they match the corresponding elements of the "to" range as
     * closely as possible in the least-squares sense.
     * 
     * @param fromPointsEnd See the documentation for argument
     * fromPointsBegin.
     * 
     * @param toPointsBegin This argument defines the beginning of one
     * of the two sets of points to be registered.  The returned
     * Transform3D<FloatType> instance will take points in the range
     * specified by arguments fromPointsBegin and fromPointsEnd, and
     * transform them so that they match as closely as possible (in
     * the least squares sense) the corresponding elements of the
     * range starting from toPointsBegin and extending (fromPointsEnd
     * - fromPointsBegin) elements.
     *
     * @param flagsBegin This argument defines the beginning of a
     * sequence of boolean values indicating which points should be
     * included in the registration.  That is, flagsBegin is an iterator
     * which must dereference to a bool, and the sequence from
     * flagsBegin to flagsBegin + (fromPointsEnd - fromPointsBegin) must
     * be valid.  Points will be included in the registration iff the
     * corresponding bool is true.
     *
     * @return The return value is a Transform3D<FloatType> instance
     * which takes the points in the "from" range and matches them to
     * the points in the "to" range.
     */
    template <class FloatType, class InIter0, class InIter1, class InIter2>
    brick::numeric::Transform3D<FloatType>
    registerPoints3D(InIter0 fromPointsBegin, InIter0 fromPointsEnd,
                     InIter1 toPointsBegin, InIter2 selectedFlagsBegin);



    /**
     * This function works just like the three-argument form of
     * registerPoints3D(), except that points can be weighted,
     * indicating how much influence they should have on the resulting
     * solution.
     *
     * @param fromPointsBegin This iterator, along with fromPointsEnd,
     * defines one of the two sets of points to be registered.  The
     * type of (*fromPointsBegin) must be Vector3D<FloatType>.  The
     * returned Transform3D<FloatType> instance will take points in
     * the range specified by this iterator pair and transform them so
     * that they match the corresponding elements of the "to" range as
     * closely as possible in the least-squares sense.
     * 
     * @param fromPointsEnd See the documentation for argument
     * fromPointsBegin.
     * 
     * @param toPointsBegin This argument defines the beginning of one
     * of the two sets of points to be registered.  The returned
     * Transform3D<FloatType> instance will take points in the range
     * specified by arguments fromPointsBegin and fromPointsEnd, and
     * transform them so that they match as closely as possible (in
     * the least squares sense) the corresponding elements of the
     * range starting from toPointsBegin and extending (fromPointsEnd
     * - fromPointsBegin) elements.
     *
     * @param weightsBegin This argument defines the beginning of a
     * sequence of FloatType values indicating the weight to be given
     * to each point.  That is, weightsBegin is an iterator which must
     * dereference to FloatType, and the sequence from flagsBegin to
     * flagsBegin + (fromPointsEnd - fromPointsBegin) must be valid.
     * The influence of each point on the result of the registration
     * will be proportional to its weight.
     *
     * @param dummy This argument distambiguates the weighted version
     * of registerPoints3D from the version that accepts boolean
     * flags.
     * 
     * @return The return value is a Transform3D<FloatType> instance
     * that takes into account the input points and their weights.
     */
    template <class FloatType, class InIter0, class InIter1, class InIter2>
    brick::numeric::Transform3D<FloatType>
    registerPoints3D(InIter0 fromPointsBegin, InIter0 fromPointsEnd,
                     InIter1 toPointsBegin, InIter2 weightsBegin,
                     bool dummy);



    /**
     * This function calls the four-argument form of registerPoints3D()
     * repeatedly while trying to identify and ignore outliers.
     * Iteration terminates when a call to the four-argument form of
     * registerPoints3D() does not change the selected outlier list.
     * Calling this function with argument inclusion set to 1.0 and
     * argument maximumResidual less than zero is just the same as
     * calling the four-argument form of of registerPoints3D with
     * flagsBegin pointing to a sequence of all true.
     *
     * @param fromPointsBegin This iterator, along with fromPointsEnd,
     * defines one of the two sets of points to be registered.  The
     * type of (*fromPointsBegin) must be Vector3D<FloatType>.  The
     * returned Transform3D<FloatType> instance will take points in
     * the range specified by this iterator pair and transform them so
     * that they match the corresponding elements of the "to" range as
     * closely as possible in the least-squares sense.
     * 
     * @param fromPointsEnd See the documentation for argument
     * fromPointsBegin.
     * 
     * @param toPointsBegin This argument defines the beginning of one
     * of the two sets of points to be registered.  The returned
     * Transform3D<FloatType> instance will take points in the range
     * specified by arguments fromPointsBegin and fromPointsEnd, and
     * transform them so that they match as closely as possible (in
     * the least squares sense) the corresponding elements of the
     * range starting from toPointsBegin and extending (fromPointsEnd
     * - fromPointsBegin) elements.  possible in the least-squares
     * sense.
     *
     * @param flagsBegin This output iterator will be used to return a
     * sequence of bools indicating which points were used in the
     * registration.  For points which were included in the
     * registration, the corresponding bool will be true.
     *
     * @param inclusion This argument specifies the proportion of the
     * dataset that is expected to be inliers.  Setting this value to
     * 1.0 or greater indicates that all of the points should be
     * included in the registration, unless their inclusion is
     * countermanded by argument maximumResidual.  Setting this value
     * less than 0.0 has the same effect as settint it to 1.0.  To
     * clarify(?): at each call to the three-argument form of
     * registerPoints3D(), at most floor(inclusion *
     * (int)(fromPointsEnd - fromPointsBegin)) pairs of points will be
     * used, and the points used will be those with the smallest
     * residual prior to the call to the three-argument form of
     * registerPoints3D().
     *
     * @param maximumResidual This argument specifies the largest
     * expected residual between corresponding points after the
     * registration.  Pairs of points that differ by more than this
     * amount will be assumed to be outliers and ignored during the
     * next registration.  Setting this argument less than zero
     * indicates that all pairs of points should be included in the
     * registration, unless some points are countermanded by argument
     * inclusion.
     *
     * @param maximumIterations This argument how many iterations are
     * permissible.  The loop (register -> compute outliers -> repeat)
     * will terminate after this many iterations even if the set of
     * points chosen as outliers is still changing.
     *
     * @return The return value is a Transform3D<FloatType> instance
     * which takes the points in the "from" range and matches them to
     * the points in the "to" range.
     */
    template <class FloatType, class InIter0, class InIter1, class OutIter0>
    brick::numeric::Transform3D<FloatType>
    registerPoints3D(InIter0 fromPointsBegin, InIter0 fromPointsEnd,
                     InIter1 toPointsBegin, OutIter0 selectedFlagsBegin,
                     FloatType inclusion, FloatType maximumResidual = -1.0,
                     size_t maximumIterations = 5);
  
  
  } // namespace computerVision    

} // namespace brick

// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/registerPoints3D_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_REGISTERPOINTS3D_HH */
