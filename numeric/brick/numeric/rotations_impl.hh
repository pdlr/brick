/**
***************************************************************************
* @file rotations_impl.cc
*
* Source file declaring functions which convert between different
* representations of 3D rotation.
*
* Copyright (C) 2005-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <cmath>
#include <brick/common/exception.hh>
#include <brick/common/functional.hh>
#include <brick/common/mathFunctions.hh>
#include <brick/numeric/rotations.hh>
#include <brick/numeric/utilities.hh>

namespace {

  using namespace brick::numeric;

  // This templated local function is a quick hack to let us test
  // whether things are approximately equal without digging too deeply
  // into numerical issues.  Ultimately, these functions should be
  // removed, and brick::numeric::NumericTraits<> should be extended
  // to provide a better version of this functionality.
  template <class Type>
  Type l_getRotationsEpsilon() {
    BRICK_THROW(brick::common::NotImplementedException,
                "getRotationsEpsilon()",
                "Rotation conversions are not implemented for this scalar "
                "type.");
    return Type(0.0);
  }

  template <>
  brick::common::Float64 l_getRotationsEpsilon<brick::common::Float64>() {
    return 1.0E-10;
  }
  template <>
  brick::common::Float32 l_getRotationsEpsilon<brick::common::Float32>() {
    return 1.0E-5;
  }
  

  template <class Type>
  Transform3D<Type>
  l_getEulerComponent(const Type& angle, Axis axis)
  {
    Type cosineAngle = brick::common::cosine(angle);
    Type sineAngle = brick::common::sine(angle);

    switch(axis) {
    case BRICK_AXIS_X:
      return Transform3D<Type>(1.0, 0.0, 0.0, 0.0,
                         0.0, cosineAngle, -sineAngle, 0.0,
                         0.0, sineAngle, cosineAngle, 0.0,
                         0.0, 0.0, 0.0, 1.0);
      break;
    case BRICK_AXIS_Y:
      return Transform3D<Type>(cosineAngle, 0.0, sineAngle, 0.0,
                         0.0, 1.0, 0.0, 0.0,
                         -sineAngle, 0.0, cosineAngle, 0.0,
                         0.0, 0.0, 0.0, 1.0);
      break;
    case BRICK_AXIS_Z:
      return Transform3D<Type>(cosineAngle, -sineAngle, 0.0, 0.0,
                         sineAngle, cosineAngle, 0.0, 0.0,
                         0.0, 0.0, 1.0, 0.0, 
                         0.0, 0.0, 0.0, 1.0);
      break;
    default:
      // Should never get here.
      BRICK_THROW(brick::common::LogicException, "l_getEulerComponent()",
                  "Unrecognized value for argument axis.");
      break;
    }
    // Should never get here either.
    return Transform3D<Type>();
  }
                    
}


namespace brick {

  namespace numeric {
    
    template <class Type>
    Quaternion<Type>
    angleAxisToQuaternion(const Type& angle, const Vector3D<Type>& axis, bool isNormalized)
    {
      // Deal with the angle.
      Type angleOverTwo = angle / 2.0;
      Type cosineValue = brick::common::cosine(angleOverTwo);
      Type sineValue = brick::common::sine(angleOverTwo);

      // If the axis is already unit length, we're done!
      if(isNormalized) {
        return Quaternion<Type>(cosineValue, sineValue * axis.x(),
                          sineValue * axis.y(), sineValue * axis.z());
      }

      // Axis is not known to be unit length, so we have more work to do.
      Vector3D<Type> axisCopy = axis;
      Type axisMagnitude = magnitude(axis);
      if(axisMagnitude != 0.0) {
        axisCopy /= axisMagnitude;
      } else {
        // Axis is too small to observe.  Pick an arbitrary axis.
        axisCopy.setValue(1.0, 0.0, 0.0);
      }
      return Quaternion<Type>(cosineValue, sineValue * axisCopy.x(),
                        sineValue * axisCopy.y(), sineValue * axisCopy.z());
    }

  
    template <class Type>
    Vector3D<Type>
    angleAxisToRollPitchYaw(const Type& angle, const Vector3D<Type>& axis,
                            bool isNormalized)
    {
      Quaternion<Type> quaternion = angleAxisToQuaternion(angle, axis, isNormalized);
      return quaternionToRollPitchYaw(quaternion);
    }    

  
    template <class Type>
    Transform3D<Type>
    angleAxisToTransform3D(const Type& angle, const Vector3D<Type>& axis, bool isNormalized)
    {
      Quaternion<Type> quaternion = angleAxisToQuaternion(angle, axis, isNormalized);
      return quaternionToTransform3D(quaternion);
    }

  
    // This function converts a rotation from general euler angle
    // representation to Transform3D representation.
    template <class Type>
    Transform3D<Type>
    eulerToTransform3D(const Type& angle0, Axis axis0,
                       const Type& angle1, Axis axis1,
                       const Type& angle2, Axis axis2)
    {
      return (l_getEulerComponent(angle2, axis2)
              * l_getEulerComponent(angle1, axis1)
              * l_getEulerComponent(angle0, axis0));
    }

    
    template <class Type>
    std::pair< Type, Vector3D<Type> >
    quaternionToAngleAxis(const Quaternion<Type>& quaternion)
    {
      // Normalize the quaternion, if necessary.
      Quaternion<Type> quaternionCopy(quaternion);
      quaternionCopy.normalize();

      // Recover the obvious quantities.
      Type cosineHalfTheta = quaternionCopy.s();
      Vector3D<Type> sineTimesAxis(
        quaternionCopy.i(), quaternionCopy.j(), quaternionCopy.k());
      Type sineHalfTheta = magnitude(sineTimesAxis);

      // Recover angle and axis.
      Type angle = 2 * std::atan2(sineHalfTheta, cosineHalfTheta);
      if(sineHalfTheta == 0.0) {
        return std::make_pair(angle, Vector3D<Type>(1.0, 0.0, 0.0));
      }
      return std::make_pair(angle, sineTimesAxis / sineHalfTheta);
    }

  
    template <class Type>
    Vector3D<Type>
    quaternionToRollPitchYaw(const Quaternion<Type>& quaternion)
    {
      // The quaternion will be normalized inside quaternionToTransform3D().
      Transform3D<Type> transform3D = quaternionToTransform3D(quaternion);
      return transform3DToRollPitchYaw(transform3D);
    }


    template <class Type>
    Transform3D<Type>
    quaternionToTransform3D(const Quaternion<Type>& quaternion)
    {
      // Normalize the quaternion, if necessary.
      Quaternion<Type> quaternionCopy(quaternion);
      quaternionCopy.normalize();

      // Some convenience variables.
      Type ii = 2.0 * quaternionCopy.i() * quaternionCopy.i();
      Type jj = 2.0 * quaternionCopy.j() * quaternionCopy.j();
      Type kk = 2.0 * quaternionCopy.k() * quaternionCopy.k();
      Type si = 2.0 * quaternionCopy.s() * quaternionCopy.i();
      Type sj = 2.0 * quaternionCopy.s() * quaternionCopy.j();
      Type sk = 2.0 * quaternionCopy.s() * quaternionCopy.k();
      Type ij = 2.0 * quaternionCopy.i() * quaternionCopy.j();
      Type ik = 2.0 * quaternionCopy.i() * quaternionCopy.k();
      Type jk = 2.0 * quaternionCopy.j() * quaternionCopy.k();

      return Transform3D<Type>(1 - jj - kk, ij - sk, ik + sj, 0.0,
                         ij + sk, 1 - ii - kk, jk - si, 0.0,
                         ik - sj, jk + si, 1 - ii - jj, 0.0,
                         0.0, 0.0, 0.0, 1);
    }

  
    template <class Type>
    std::pair< Type, Vector3D<Type> >
    rollPitchYawToAngleAxis(const Vector3D<Type>& rollPitchYaw)
    {
      Quaternion<Type> quaternion = rollPitchYawToQuaternion(rollPitchYaw);
      return quaternionToAngleAxis(quaternion);
    }

  
    template <class Type>
    Quaternion<Type>
    rollPitchYawToQuaternion(const Vector3D<Type>& rollPitchYaw)
    {
      Transform3D<Type> transform3D = rollPitchYawToTransform3D(rollPitchYaw);
      return transform3DToQuaternion(transform3D);
    }
  

    template <class Type>
    Transform3D<Type>
    rollPitchYawToTransform3D(const Vector3D<Type>& rollPitchYaw)
    {
      // First compute some convenience values.
      Type cosineRoll = brick::common::cosine(rollPitchYaw.x());
      Type cosinePitch = brick::common::cosine(rollPitchYaw.y());
      Type cosineYaw = brick::common::cosine(rollPitchYaw.z());

      Type sineRoll = brick::common::sine(rollPitchYaw.x());
      Type sinePitch = brick::common::sine(rollPitchYaw.y());
      Type sineYaw = brick::common::sine(rollPitchYaw.z());

      // Each of roll, pitch, yaw, correspond to a rotation about one
      // axis.
      Transform3D<Type> rollTransform(1.0, 0.0, 0.0, 0.0,
                                0.0, cosineRoll, -sineRoll, 0.0,
                                0.0, sineRoll, cosineRoll, 0.0,
                                0.0, 0.0, 0.0, 1.0);
      Transform3D<Type> pitchTransform(cosinePitch, 0.0, sinePitch, 0.0,
                                 0.0, 1.0, 0.0, 0.0,
                                 -sinePitch, 0.0, cosinePitch, 0.0,
                                 0.0, 0.0, 0.0, 1.0);
      Transform3D<Type> yawTransform(cosineYaw, -sineYaw, 0.0, 0.0,
                               sineYaw, cosineYaw, 0.0, 0.0,
                               0.0, 0.0, 1.0, 0.0, 
                               0.0, 0.0, 0.0, 1.0);

      // Compose the three rotations to get the result.
      return rollTransform * (pitchTransform * yawTransform);
    }

  
    template <class Type>
    std::pair< Type, Vector3D<Type> >
    transform3DToAngleAxis(const Transform3D<Type>& transform3D)
    {
      Quaternion<Type> quaternion = transform3DToQuaternion(transform3D);
      return quaternionToAngleAxis(quaternion);
    }

  
    // This routine draws from _On Homogeneous Transforms, Quaternions,
    // and Computation Efficiency_, by Funda, Taylor, and Paul, IEEE R&A,
    // June 1990.
    template <class Type>
    Quaternion<Type>
    transform3DToQuaternion(const Transform3D<Type>& transform3D)
    {
      // For convenience, we make temporary variables representing the
      // elements of transform3D.
      // 
      // Type t00 = transform3D.value<0, 0>();
      // Type t01 = transform3D.value<0, 1>();
      // Type t02 = transform3D.value<0, 2>();
      // Type t10 = transform3D.value<1, 0>();
      // Type t11 = transform3D.value<1, 1>();
      // Type t12 = transform3D.value<1, 2>();
      // Type t20 = transform3D.value<2, 0>();
      // Type t21 = transform3D.value<2, 1>();
      // Type t22 = transform3D.value<2, 2>();
      Type t00 = transform3D(0, 0);
      Type t01 = transform3D(0, 1);
      Type t02 = transform3D(0, 2);
      Type t10 = transform3D(1, 0);
      Type t11 = transform3D(1, 1);
      Type t12 = transform3D(1, 2);
      Type t20 = transform3D(2, 0);
      Type t21 = transform3D(2, 1);
      Type t22 = transform3D(2, 2);
    
      // First compute s.
      Type sSquaredTimesFour = (1.0 + t00 + t11 + t22);
      // Allow for numerical errors.
      if(sSquaredTimesFour < 0.0) {
        sSquaredTimesFour = 0.0;
      }
      Type sValue = brick::common::squareRoot(sSquaredTimesFour) / 2.0;
      // Allow for numerical errors.
      if(sValue > 1.0) {
        sValue = 1.0;
      }

      // This vector points in the direction of the axis of rotation, but
      // unfortunately goes to zero for theta approaching 0 degrees and
      // 180 degrees.
      Vector3D<Type> axis0(t21 - t12, t02 - t20, t10 - t01);
    
      // We need to find another parallel vector using the elements of
      // transform3D. Start by noting which axis dominates.
      int axisOfLargestRotation = 0;
      if(t00 > t11) {
        if(t00 > t22) {
          axisOfLargestRotation = 0;
        } else {
          axisOfLargestRotation = 2;
        }
      } else {
        if(t11 > t22) {
          axisOfLargestRotation = 1;
        } else {
          axisOfLargestRotation = 2;
        }
      }

      // Now compute the parallel vector and add it to axis0.
      if(axisOfLargestRotation == 0) {
        Vector3D<Type> axis1(1.0 + t00 - t11 - t22, t10 + t01, t20 + t02);
        if(axis0.x() >= 0) {
          axis0 += axis1;
        } else {
          axis0 -= axis1;
        }
      } else if(axisOfLargestRotation == 1) {
        Vector3D<Type> axis1(t10 + t01, 1.0 + t11 - t00 - t22, t21 + t12);
        if(axis0.y() >= 0) {
          axis0 += axis1;
        } else {
          axis0 -= axis1;
        }
      } else if(axisOfLargestRotation == 2) {
        Vector3D<Type> axis1(t20 + t02, t21 + t12, 1.0 + t22 - t00 - t11);
        if(axis0.z() >= 0) {
          axis0 += axis1;
        } else {
          axis0 -= axis1;
        }
      }

      // Now see about normalizing the unit quaternion.
      Type axisMagnitudeSquared = dot(axis0, axis0);
      if(approximatelyEqual(axisMagnitudeSquared,
                            Type(0.0), l_getRotationsEpsilon<Type>())) {
        // Hmm, we still have a very small axis.  Assume this means that
        // the input rotation is nearly zero.
        if(sValue >= 0) {
          return Quaternion<Type>(1.0, 0.0, 0.0, 0.0);
        } else {
          return Quaternion<Type>(-1.0, 0.0, 0.0, 0.0);
        }
      }
      axis0 *= brick::common::squareRoot((1 - (sValue * sValue)) / axisMagnitudeSquared);

      // Done.
      return Quaternion<Type>( sValue, axis0.x(), axis0.y(), axis0.z());
    }

  
    template <class Type>
    Vector3D<Type>
    transform3DToRollPitchYaw(const Transform3D<Type>& transform3D)
    {
      // There must be a better way to get this value.
      const Type piOverTwo = 1.57079632679;

      // Start by recovering pitch.
      // Type sinePitch = transform3D.value<0, 2>();
      Type sinePitch = transform3D();
      // Type cosinePitch = brick::common::squareRoot(
      //   transform3D.value<0, 0>() * transform3D.value<0, 0>()
      //   + transform3D.value<0, 1>() * transform3D.value<0, 1>());
      Type cosinePitch = brick::common::squareRoot(
        transform3D(0, 0) * transform3D(0, 0)
        + transform3D(0, 1) * transform3D(0, 1));
      Type pitch = std::atan2(sinePitch, cosinePitch);

      // Now recover roll and yaw.
      Type roll;
      Type yaw;
      if((!approximatelyEqual(pitch, piOverTwo,
                              l_getRotationsEpsilon<Type>()))
         && (!approximatelyEqual(pitch, -piOverTwo,
                                 l_getRotationsEpsilon<Type>()))) {
        // Choose a more numerically stable version of cos(pitch).
        cosinePitch = brick::common::cosine(pitch);
        // roll = std::atan2(-(transform3D.value<1, 2>() / cosinePitch),
        //                   (transform3D.value<2, 2>() / cosinePitch));
        // yaw = std::atan2(-(transform3D.value<0, 1>() / cosinePitch),
        //                  (transform3D.value<0, 0>() / cosinePitch));
        roll = std::atan2(-(transform3D(1, 2) / cosinePitch),
                          (transform3D(2, 2) / cosinePitch));
        yaw = std::atan2(-(transform3D(0, 1) / cosinePitch),
                         (transform3D(0, 0) / cosinePitch));

      } else {
        // roll = std::atan2(transform3D.value<2, 1>(),
        //                   transform3D.value<1, 1>());
        roll = std::atan2(transform3D(2, 1),
                          transform3D(1, 1));
        yaw = 0.0;
      }
      return Vector3D<Type>(roll, pitch, yaw);
    }

  
  } // namespace numeric

} // namespace brick
