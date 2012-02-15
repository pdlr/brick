/**
***************************************************************************
* @file brick/numeric/rotations.hh
*
* Header file declaring functions which convert between different
* representations of 3D rotation.
*
* Copyright (C) 2005-2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
* The functions in this file deal with several different ways of
* representing 3D rotation, which are described here:
*
* - Angle-Axis representation is specified by a scalar and a Vector3D
* instance.  The scalar specifies the size of the rotation in radians.
* Positive numbers mean clockwise rotation as you look from the origin
* along the axis of rotation.  Negative numbers mean counter-clockwise
* rotation.  The Vector3D instance specifies the axis of rotation.
*
* - Unit Quaternion representation uses a single Quaternion< instance.
* For a given rotation, alpha, around the unit axis [x, y, z], the
* corresponding unit quaternion representation is
* Quaternion(cos(alpha/2.0), sin(alpha/2.0) * x, sin(alpha/2.0) * y,
* sin(alpha/2.0) * z).
*
* - Roll-Pitch-Yaw representation specifies a 3D rotation as
* consecutive rotations around the three axes assuming ISO
* coordinates.  Roll, pitch, and yaw are all specified in radians.
* The 3D rotation is obtained by first rotating to the yaw angle
* around the (vertical) Z axis, then rotating to the pitch angle
* around the new (transverse) Y axis, and finally rotating to the roll
* angle around the new (fore-aft) X axis.
*
* - Transform3D representation specifies a 3D rotation using a
* coordinate transformation matrix.  Only the upper-left 3x3
* sub-matrix is considered.  In order to be a valid rotation matrix,
* the rows (and therefore columns) of this 3x3 sub-matrix must be
* orthonormal.  Note that no attempt is made to check for this
* condition, or to correct for non-orthonormal matrices.  If you pass
* an invalid rotation matrix to any of these routines, the results are
* undefined.
*
* - Euler angle representation is a generalization of Roll-Pitch-Yaw
* representation in which the three rotations need not be around the
* Z, Y, and X axes.  Euler angle routines allow you to explicitly set
* the rotation axes.
*
* - Rodrigues parameters specify rotation using a 3D vector pointing
* along the axis of rotation.  The length of the vector specifies the
* rotation angle in radians.  Rotation follows the right-hand rule
* around the rotation axis.  Note the similarity to Angle-Axis
* representation.  Rodrigues parameters are useful for interacting
* with OpenCV.
* 
***************************************************************************
**/

#ifndef BRICK_NUMERIC_ROTATIONS_HH
#define BRICK_NUMERIC_ROTATIONS_HH

#include <brick/numeric/quaternion.hh>
#include <brick/numeric/transform3D.hh>
#include <brick/numeric/vector3D.hh>

namespace brick {

  namespace numeric {

    /**
     ** Enum for representing rotation axes.
     **/
    enum Axis {
      BRICK_AXIS_X,
      BRICK_AXIS_Y,
      BRICK_AXIS_Z
    };
    
      
    /** 
     * This function converts a rotation from angle-axis representation
     * to unit quaternion representation.
     * 
     * @param angle This argument specifies the size of the rotation in
     * radians.  For further information, please see the rotations.hh
     * file comment.
     * 
     * @param axis This argument specifies the axis of rotation.  For
     * further information, please see the rotations.hh file comment.
     *
     * @param isNormalized If set to true, this argument disables
     * normalization of the axis prior to use.  You should only set it
     * to true if you know that the magnitude of argument axis is equal
     * to 1.0.
     * 
     * @return The return value is a unit quaternion representing the
     * specified rotation.  For further information, please see the
     * rotations.hh file comment.
     */
    template <class Type>
    Quaternion<Type>
    angleAxisToQuaternion(Type const& angle, const Vector3D<Type>& axis,
                          bool isNormalized=false);
  

    /** 
     * This function converts a rotation from angle-axis representation
     * to roll-pitch-yaw representation.
     * 
     * @param angle This argument specifies the size of the rotation in
     * radians.  For further information, please see the rotations.hh
     * file comment.
     * 
     * @param axis This argument specifies the axis of rotation.  For
     * further information, please see the rotations.hh file comment.
     * 
     * @param isNormalized If set to true, this argument disables
     * normalization of the axis prior to use.  You should only set it
     * to true if you know that the magnitude of argument axis is equal
     * to 1.0.
     *
     * @return The return value is a Vector3D instance containing roll,
     * pitch, and yaw.  For further information, please see the
     * rotations.hh file comment.
     */
    template <class Type>
    Vector3D<Type>
    angleAxisToRollPitchYaw(Type const& angle, const Vector3D<Type>& axis,
                            bool isNormalized=false);
  

    /** 
     * This function converts a rotation from angle-axis representation
     * to Transform3D representation.
     * 
     * @param angle This argument specifies the size of the rotation in
     * radians.  For further information, please see the rotations.hh
     * file comment.
     * 
     * @param axis This argument specifies the axis of rotation.  For
     * further information, please see the rotations.hh file comment.
     * 
     * @param isNormalized If set to true, this argument disables
     * normalization of the axis prior to use.  You should only set it
     * to true if you know that the magnitude of argument axis is equal
     * to 1.0.
     *
     * @return The return value is a Transform3D instance representing
     * the rotation.  For further information, please see the
     * rotations.hh file comment.
     */
    template <class Type>
    Transform3D<Type>
    angleAxisToTransform3D(Type const& angle, const Vector3D<Type>& axis,
                           bool isNormalized=false);
  

    /** 
     * This function converts a rotation from general euler angle
     * representation to Transform3D representation.  The resulting
     * transform will take points in the "from" coordinate system,
     * rotate them around axis0, then rotate around the new (rotated)
     * axis1, then rotate around the new (doubly rotated) axis2.
     *
     * One way to clarify this rotation order is to note that the
     * following three Transform3D instances are equivalent (up to
     * numerical precision):
     *
     * @code
     *   Axis anyAxis;
     *
     *   Transform3D<double> xf0 = eulerToTransform3D(
     *     angle0, axis0, angle1, axis1, angle2, axis2);
     *
     *   Transform3D<double> xf1 = (
     *     eulerToTransform3D(angle2, axis2, 0.0, anyAxis, 0.0, anyAxis)
     *     * eulerToTransform3D(angle1, axis1, 0.0, anyAxis, 0.0, anyAxis)
     *     * eulerToTransform3D(angle0, axis0, 0.0, anyAxis, 0.0, anyAxis));
     *
     *   Transform3D<double> xf2 = (
     *     eulerToTransform3D(0.0, anyAxis, 0.0, anyAxis, angle2, axis2)
     *     * eulerToTransform3D(0.0, anyAxis, angle1, axis1, 0.0, anyAxis)
     *     * eulerToTransform3D(angle0, axis0, 0.0, anyAxis, 0.0, anyAxis));
     * @code
     * 
     * @param angle0 This argument indicates the angle of rotation, in
     * radians, around the first axis.
     * 
     * @param axis0 This argument specifies which axis should be
     * rotated around first.
     * 
     * @param angle1 This argument indicates the angle of rotation, in
     * radians, around the second axis.
     * 
     * @param axis1 This argument specifies which axis should be
     * rotated around second.
     * 
     * @param angle2 This argument indicates the angle of rotation, in
     * radians, around the third axis.
     * 
     * @param axis2 This argument specifies which axis should be
     * rotated around third.
     * 
     * @return The return value is a rotation matrix equivalent to the
     * composition of the three euler components.
     */
    template <class Type>
    Transform3D<Type>
    eulerToTransform3D(Type const& angle0, Axis axis0,
                       Type const& angle1, Axis axis1,
                       Type const& angle2, Axis axis2);
  

    /** 
     * This function converts a rotation from Quaternion representation
     * to Angle-Axis representation.
     *
     * @param quaternion This argument specifies the rotation using a
     * Quaternion instance.  If the Quaternion is not normalized, the
     * normalization will be done inside the function. For further
     * information, please see the rotations.hh file comment.
     * 
     * @return The return value is a pair which specifies the Angle-Axis
     * representation of the input rotation.  For further information,
     * please see the rotations.hh file comment.
     */
    template <class Type>
    std::pair< Type, Vector3D<Type> >
    quaternionToAngleAxis(const Quaternion<Type>& quaternion);
  

    /** 
     * This function converts a rotation from Quaternion representation
     * to Roll-Pitch-Yaw representation.
     * 
     * @param quaternion This argument specifies the rotation using a
     * Quaternion instance.  If the Quaternion is not normalized, the
     * normalization will be done inside the function. For further
     * information, please see the rotations.hh file comment.
     * 
     * @return The return value is a Vector3D instance containing roll,
     * pitch, and yaw.  For further information, please see the
     * rotations.hh file comment.
     */
    template <class Type>
    Vector3D<Type>
    quaternionToRollPitchYaw(const Quaternion<Type>& quaternion);


    /** 
     * This function converts a rotation from Quaternion representation
     * to Transform3D representation.
     * 
     * @param quaternion This argument specifies the rotation using a
     * Quaternion instance.  If the Quaternion is not normalized, the
     * normalization will be done inside the function. For further
     * information, please see the rotations.hh file comment.
     * 
     * @return The return value is a Transform3D instance representing
     * the rotation.  For further information, please see the
     * rotations.hh file comment.
     */
    template <class Type>
    Transform3D<Type>
    quaternionToTransform3D(const Quaternion<Type>& quaternion);

    
    /** 
     * This function converts a rotation from Rodrigues representation
     * to Transform3D representation.
     * 
     * @param rodrigues This argument is a Vector3D instance in which
     * the x, y, and z components specify the first, second, and third
     * Rodrigues parameters, respectively.
     * 
     * @return The return value is a Transform3D instance representing
     * the rotation.  For further information, please see the
     * rotations.hh file comment.
     */
    template <class Type>
    Transform3D<Type>
    rodriguesToTransform3D(Vector3D<Type> const& rodrigues);
  

    /** 
     * This function converts a rotation from Roll-Pitch-Yaw representation
     * to Angle-Axis representation.
     *
     * @param rollPitchYaw This argument specifies the rotation using a
     * Vector3D instance representing Roll, Pitch, and Yaw.  For further
     * information, please see the rotations.hh file comment.
     *
     * @return The return value is a pair which specifies the Angle-Axis
     * representation of the input rotation.  For further information,
     * please see the rotations.hh file comment.
     */
    template <class Type>
    std::pair< Type, Vector3D<Type> >
    rollPitchYawToAngleAxis(const Vector3D<Type>& rollPitchYaw);

  
    /** 
     * This function converts a rotation from Roll-Pitch-Yaw representation
     * to Quaternion representation.
     * 
     * @param rollPitchYaw This argument specifies the rotation using a
     * Vector3D instance representing Roll, Pitch, and Yaw.  For further
     * information, please see the rotations.hh file comment.
     * 
     * @return The return value is a unit quaternion representing the
     * specified rotation.  For further information, please see the
     * rotations.hh file comment.
     */
    template <class Type>
    Quaternion<Type>
    rollPitchYawToQuaternion(const Vector3D<Type>& rollPitchYaw);

  
    /** 
     * This function converts a rotation from Roll-Pitch-Yaw representation
     * to Transform3D representation.
     * 
     * @param rollPitchYaw This argument specifies the rotation using a
     * Vector3D instance representing Roll, Pitch, and Yaw.  For further
     * information, please see the rotations.hh file comment.
     * 
     * @return The return value is a Transform3D instance representing
     * the rotation.  For further information, please see the
     * rotations.hh file comment.
     */
    template <class Type>
    Transform3D<Type>
    rollPitchYawToTransform3D(const Vector3D<Type>& rollPitchYaw);


    /** 
     * This function converts a rotation from Transform3D representation
     * to Angle-Axis representation.
     *
     * @param transform3D This argument specifies the rotation using a
     * Transform3D instance.  For further information, please see the
     * rotations.hh file comment.
     *
     * @return The return value is a pair which specifies the Angle-Axis
     * representation of the input rotation.  For further information,
     * please see the rotations.hh file comment.
     */
    template <class Type>
    std::pair< Type, Vector3D<Type> >
    transform3DToAngleAxis(const Transform3D<Type>& transform3D);


    /** 
     * This function converts a rotation from Transform3D representation
     * to Quaternion representation.
     * 
     * This routine draws from _On Homogeneous Transforms, Quaternions,
     * and Computation Efficiency_, by Funda, Taylor, and Paul, IEEE R&A,
     * June 1990.
     *
     * @param transform3D This argument specifies the rotation using a
     * Transform3D instance.  For further information, please see the
     * rotations.hh file comment.
     * 
     * @return The return value is a unit quaternion representing the
     * specified rotation.  For further information, please see the
     * rotations.hh file comment.
     */
    template <class Type>
    Quaternion<Type>
    transform3DToQuaternion(const Transform3D<Type>& transform3D);


    /** 
     * This function converts a rotation from Transform3D representation
     * to Rodrigues representation.
     * 
     * @param transform3D This argument specifies the rotation using a
     * Transform3D instance.  For further information, please see the
     * rotations.hh file comment.
     * 
     * @return The return value is a Vector3D instance in which the x,
     * y, and z components represent the first, second, and third
     * Rodrigues parameters, respectively.  For further information,
     * please see the rotations.hh file comment.
     */
    template <class Type>
    Vector3D<Type>
    transform3DToRodrigues(const Transform3D<Type>& transform3D);


    /** 
     * This function converts a rotation from Transform3D representation
     * to Roll-Pitch-Yaw representation.
     *
     * @param transform3D This argument specifies the rotation using a
     * Transform3D instance.  For further information, please see the
     * rotations.hh file comment.
     * 
     * @return The return value is a Vector3D instance containing roll,
     * pitch, and yaw.  For further information, please see the
     * rotations.hh file comment.
     */
    template <class Type>
    Vector3D<Type>
    transform3DToRollPitchYaw(const Transform3D<Type>& transform3D);

  } // namespace numeric

} // namespace brick

#include <brick/numeric/rotations_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_ROTATIONS_HH */
