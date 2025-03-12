/**
***************************************************************************
* @file brick/computerVision/cameraIntrinsicsDistortedPinhole_impl.hh
*
* Definitions of inline and template functions for
* CameraIntrinsicsDistortedPinhole.
*
* Copyright (C) 2007-2014 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_CAMERAINTRINSICSDISTORTEDPINHOLE_IMPL_HH
#define BRICK_COMPUTERVISION_CAMERAINTRINSICSDISTORTEDPINHOLE_IMPL_HH

// This file is included by cameraIntrinsicsDistortedPinhole.hh, and should
// not be directly included by user code, so no need to include
// cameraIntrinsicsDistortedPinhole.hh here.
//
// #include <brick/numeric/cameraIntrinsicsDistortedPinhole.hh>

#include <iomanip>
#include <brick/common/expect.hh>
#include <brick/optimization/optimizerBFGS.hh>
#include <brick/optimization/optimizerNelderMead.hh>

namespace brick {

  namespace computerVision {

    namespace privateCode {

      /**
       ** This functor is called by
       ** CameraIntrinsicsDistortedPinhole::reverseProject() during
       ** iterative approximation of reverse projection.
       **/
      template <class FloatType>
      class ReverseProjectionObjective
      {
      public:
        using argument_type = brick::numeric::Array1D<FloatType>;
        using result_type = FloatType;

        ReverseProjectionObjective();

        ReverseProjectionObjective(
          CameraIntrinsicsDistortedPinhole<FloatType> const& intrinsics,
          const brick::numeric::Vector2D<FloatType>& uvTarget,
          FloatType const& maxAzimuthTangent,
          FloatType const& maxElevationTangent);

        FloatType
        operator()(const brick::numeric::Array1D<FloatType>& theta);

        FloatType
        getOffset();

        brick::numeric::Array1D<FloatType>
        gradient(const brick::numeric::Array1D<FloatType>& theta);

      private:

        FloatType
        computeBoundsPenalty(
          const brick::numeric::Vector2D<FloatType>& uvPosition);

        void
        computeBoundsPenaltyGradient(FloatType uValue, FloatType vValue,
                                     FloatType dUdX, FloatType dUdY,
                                     FloatType dVdX, FloatType dVdY,
                                     FloatType& dPdX, FloatType& dPdY);

        FloatType
        computeFieldOfViewPenalty(
          const brick::numeric::Vector3D<FloatType>& candidate);

        void
        computeFieldOfViewPenaltyGradient(
          const FloatType& xValue, const FloatType& yValue,
          FloatType& dPdX, FloatType& dPdY);

        CameraIntrinsicsDistortedPinhole<FloatType> const* m_intrinsicsPtr;
        FloatType m_offset;
        brick::numeric::Vector2D<FloatType> m_uvTarget;
        FloatType m_maxAzimuthTangent;
        FloatType m_maxElevationTangent;
        FloatType m_azimuthViolationScale;
        FloatType m_elevationViolationScale;
      };

      // Implementation of ReverseProjectionObjective is at the bottom of this
      // file.

    } // namespace privateCode;



    // Default constructor.
    template <class FloatType>
    CameraIntrinsicsDistortedPinhole<FloatType>::
    CameraIntrinsicsDistortedPinhole()
      : CameraIntrinsics<FloatType>(),
        m_pinholeIntrinsics()
    {
      // Empty.
    }


    // This constructor allows the caller to explicitly set the
    // camera pinhole intrinsic parameters.
    template <class FloatType>
    CameraIntrinsicsDistortedPinhole<FloatType>::
    CameraIntrinsicsDistortedPinhole(unsigned int numPixelsX,
                                     unsigned int numPixelsY,
                                     FloatType focalLengthX,
                                     FloatType focalLengthY,
                                     FloatType centerU,
                                     FloatType centerV)
      : CameraIntrinsics<FloatType>(),
        m_pinholeIntrinsics(numPixelsX, numPixelsY,
                            focalLengthX, focalLengthY,
                            centerU, centerV)
    {
      // Empty.
    }


    // This member function takes a point in 2D pixel coordinates
    // and returns a ray in 3D camera coordinates passing through
    // all of the 3D points that project to the specified 2D
    // position.
    template <class FloatType>
    geometry::Ray3D<FloatType>
    CameraIntrinsicsDistortedPinhole<FloatType>::
    reverseProject(const brick::numeric::Vector2D<FloatType>& pixelPosition,
                   bool normalize,
                   FloatType const& maxAzimuthTangent,
                   FloatType const& maxElevationTangent) const
    {
      // We need a starting point for our optimization.  Ostensibly,
      // this is a 3D point that gets projected into the image.  The
      // optimization is then to tweak the point until it projects to
      // pixelPosition.  A full 3D point would be an
      // overparameterization, however, so we assume Z = 1, and
      // represent the point with a Vector2D instance.  We want to
      // start the optimization with a point we're _sure_ will project
      // into the image, so we choose the safest start point we can
      // think of... the center of projection.
      brick::numeric::Array1D<FloatType> startPoint(2);
      startPoint = FloatType(0.0);

      // Now optimize to find a better answer.
      typedef privateCode::ReverseProjectionObjective<FloatType> LocalObjective;
      LocalObjective objective(*this, pixelPosition,
                               maxAzimuthTangent, maxElevationTangent);
      OptimizerBFGS<LocalObjective, FloatType> optimizer(objective);
      optimizer.setStartPoint(startPoint);
      brick::numeric::Array1D<FloatType> endPoint = optimizer.optimum();
      FloatType residual = optimizer.optimalValue();

      // This would be the place to check convergence and try a
      // different start point, if we were so inclined.

      if(residual > (this->getMaximumReverseProjectionResidual()
                     + objective.getOffset())) {
        BRICK_THROW(
          brick::common::ValueException,
          "CameraIntrinsicsDistortedPinhole<FloatType>::reverseProject()",
          "Reverse projection failed to converge.");
      }

      return brick::geometry::Ray3D<FloatType>(
        brick::numeric::Vector3D<FloatType>(0.0, 0.0, 0.0),
        brick::numeric::Vector3D<FloatType>(endPoint[0], endPoint[1], 1.0),
        normalize);
    }


    // Implementation of ReverseProjectionObjective.
    namespace privateCode {

      template <class FloatType>
      ReverseProjectionObjective<FloatType>::
      ReverseProjectionObjective()
        : m_intrinsicsPtr(0),
          m_offset(0.0),
          m_uvTarget(0.0, 0.0),
          m_maxAzimuthTangent(0.0),
          m_maxElevationTangent(0.0),
          m_azimuthViolationScale(0.0),
          m_elevationViolationScale(0.0)
      {
        // Empty.
      }


      template <class FloatType>
      ReverseProjectionObjective<FloatType>::
      ReverseProjectionObjective(
        const CameraIntrinsicsDistortedPinhole<FloatType>& intrinsics,
        const brick::numeric::Vector2D<FloatType>& uvTarget,
        FloatType const& maxAzimuthTangent,
        FloatType const& maxElevationTangent)
        : m_intrinsicsPtr(&intrinsics),
          m_offset(0.0001),
          m_uvTarget(uvTarget),
          m_maxAzimuthTangent(maxAzimuthTangent),
          m_maxElevationTangent(maxElevationTangent),
          m_azimuthViolationScale(0.0),
          m_elevationViolationScale(0.0)
      {
        // See comments in memberFunction getFieldOfViewPenalty() for
        // an explantation of these two lines.
        this->m_azimuthViolationScale =
          (FloatType(10.0) * FloatType(intrinsics.getImageWidth())
           / maxAzimuthTangent);
        this->m_elevationViolationScale =
          (FloatType(10.0) * FloatType(intrinsics.getImageHeight())
           / maxElevationTangent);
      }


      template <class FloatType>
      FloatType
      ReverseProjectionObjective<FloatType>::
      operator()(const brick::numeric::Array1D<FloatType>& theta)
      {
        if(m_intrinsicsPtr == 0) {
          return FloatType(-1.0);
        }
        brick::numeric::Vector3D<FloatType> candidate(theta[0], theta[1], 1.0);
        brick::numeric::Vector2D<FloatType> projection =
          m_intrinsicsPtr->project(candidate);

        // The high order terms of the distortion model can lead to
        // false minima outside the image boundaries.  Depending on
        // what information is available, we use one of two methods to
        // mitigate this.
        FloatType boundsPenalty(0.0);
        if((m_maxAzimuthTangent > FloatType(0.0))
           || (m_maxElevationTangent > FloatType(0.0))) {
          boundsPenalty = this->computeFieldOfViewPenalty(candidate);
        } else {
          boundsPenalty = this->computeBoundsPenalty(projection);
        }

        return (brick::numeric::magnitudeSquared<FloatType>(
                  projection - m_uvTarget) + boundsPenalty + m_offset);
      }


      template <class FloatType>
      FloatType
      ReverseProjectionObjective<FloatType>::
      getOffset()
      {
        return m_offset;
      }


      template <class FloatType>
      brick::numeric::Array1D<FloatType>
      ReverseProjectionObjective<FloatType>::
      gradient(const brick::numeric::Array1D<FloatType>& theta)
      {
        brick::numeric::Array1D<FloatType> gradientArray(2);
        if(m_intrinsicsPtr == 0) {
          gradientArray = FloatType(0.0);
          return gradientArray;
        }

        FloatType uValue;
        FloatType vValue;
        FloatType dUdX;
        FloatType dUdY;
        FloatType dVdX;
        FloatType dVdY;
        m_intrinsicsPtr->projectWithPartialDerivatives(
          theta[0], theta[1], uValue, vValue, dUdX, dUdY, dVdX, dVdY);

        FloatType twoTimesDeltaU = 2.0 * (uValue - m_uvTarget.x());
        FloatType twoTimesDeltaV = 2.0 * (vValue - m_uvTarget.y());

        // Using chain rule,
        //
        // (d/dx)(deltaU^2 + deltaV^2) =
        //    2*deltaU*(d/dx)deltaU + 2*deltaV*(d/dx)deltaV
        //
        // (d/dy)(deltaU^2 + deltaV^2) =
        //    2*deltaU*(d/dy)deltaU + 2*deltaV*(d/dy)deltaV
        gradientArray[0] = twoTimesDeltaU * dUdX + twoTimesDeltaV * dVdX;
        gradientArray[1] = twoTimesDeltaU * dUdY + twoTimesDeltaV * dVdY;

        // Don't forget to add gradient for bounds penalty.
        FloatType dPdX;
        FloatType dPdY;
        if((m_maxAzimuthTangent > FloatType(0.0))
           || (m_maxElevationTangent > FloatType(0.0))) {
          this->computeFieldOfViewPenaltyGradient(
            theta[0], theta[1], dPdX, dPdY);
        } else {
          this->computeBoundsPenaltyGradient(
            uValue, vValue, dUdX, dUdY, dVdX, dVdY, dPdX, dPdY);
        }

        gradientArray[0] += dPdX;
        gradientArray[1] += dPdY;

        return gradientArray;
      }


      template <class FloatType>
      FloatType
      ReverseProjectionObjective<FloatType>::
      computeBoundsPenalty(const brick::numeric::Vector2D<FloatType>& uvPosition)
      {
        if(m_intrinsicsPtr == 0) {
          return 0.0;
        }

        // Start by computing a vaguely normalized measure of how far
        // out-of-bounds uvPosition is.  We'll call this "violation."
        FloatType uViolation = 0.0;
        FloatType scaleU = m_intrinsicsPtr->getNumPixelsX() / 10.0;
        if(uvPosition.x() < static_cast<FloatType>(0)) {
          uViolation = -uvPosition.x() / scaleU;
        } else if(uvPosition.x()
                  >= static_cast<FloatType>(m_intrinsicsPtr->getNumPixelsX())) {
          uViolation =
            ((uvPosition.x()
              - static_cast<FloatType>(m_intrinsicsPtr->getNumPixelsX()))
             / scaleU);
        }
        FloatType vViolation = 0.0;
        FloatType scaleV = m_intrinsicsPtr->getNumPixelsY() / 10.0;
        if(uvPosition.y() < static_cast<FloatType>(0)) {
          vViolation = -uvPosition.y() / scaleV;
        } else if(uvPosition.y()
                  >= static_cast<FloatType>(m_intrinsicsPtr->getNumPixelsY())) {
          vViolation =
            ((uvPosition.y()
              - static_cast<FloatType>(m_intrinsicsPtr->getNumPixelsY()))
             / scaleV);
        }

        // Bounds penaly is violation**8 so as to dominate the r**6
        // term in the distortion model.
        uViolation *= uViolation;
        uViolation *= uViolation;
        uViolation *= uViolation;
        vViolation *= vViolation;
        vViolation *= vViolation;
        vViolation *= vViolation;
        return uViolation + vViolation;
      }


      template <class FloatType>
      void
      ReverseProjectionObjective<FloatType>::
      computeBoundsPenaltyGradient(FloatType uValue, FloatType vValue,
                                   FloatType dUdX, FloatType dUdY,
                                   FloatType dVdX, FloatType dVdY,
                                   FloatType& dPdX, FloatType& dPdY)
      {
        if(m_intrinsicsPtr == 0) {
          dPdX = 0.0;
          dPdY = 0.0;
          return;
        }

        // Start by computing a vaguely normalized measure of how far
        // out-of-bounds uvPosition is.  We'll call this "violation."
        // Simultaneously compute the first derivative of violation
        // wrt x and y.
        FloatType uViolation = 0.0;
        FloatType scaleU = m_intrinsicsPtr->getNumPixelsX() / 10.0;
        FloatType dUVdX = 0.0;  // <-- derivative of uViolation wrt X.
        FloatType dUVdY = 0.0;  // <-- derivative of uViolation wrt Y.
        if(uValue < static_cast<FloatType>(0)) {
          uViolation = -uValue / scaleU;
          dUVdX = -dUdX / scaleU;
          dUVdY = -dUdY / scaleU;
        } else if(uValue
                  >= static_cast<FloatType>(m_intrinsicsPtr->getNumPixelsX())) {
          uViolation =
            (uValue - static_cast<FloatType>(m_intrinsicsPtr->getNumPixelsX()))
            / scaleU;
          dUVdX = dUdX / scaleU;
          dUVdY = dUdY / scaleU;
        }
        FloatType vViolation = 0.0;
        FloatType scaleV = m_intrinsicsPtr->getNumPixelsY() / 10.0;
        FloatType dVVdX = 0.0;  // <-- derivative of vViolation wrt X.
        FloatType dVVdY = 0.0;  // <-- derivative of vViolation wrt Y.
        if(vValue < static_cast<FloatType>(0)) {
          vViolation = -vValue / scaleV;
          dVVdX = -dVdX / scaleV;
          dVVdY = -dVdY / scaleV;
        } else if(vValue
                  >= static_cast<FloatType>(m_intrinsicsPtr->getNumPixelsY())) {
          vViolation =
            (vValue - static_cast<FloatType>(m_intrinsicsPtr->getNumPixelsY()))
            / scaleV;
          dVVdX = dVdX / scaleV;
          dVVdY = dVdY / scaleV;
        }
        // Bounds penaly is violation**8 so as to dominate the r**6
        // term in the distortion model.  Using chain rule, then first
        // derivative is (8 * violation^7 * derivative-of-violation).
        FloatType uViolationTo2 = uViolation * uViolation;
        FloatType uViolationTo4 = uViolationTo2 * uViolationTo2;
        FloatType uViolationTo7 = uViolationTo4 * uViolationTo2 * uViolation;
        FloatType vViolationTo2 = vViolation * vViolation;
        FloatType vViolationTo4 = vViolationTo2 * vViolationTo2;
        FloatType vViolationTo7 = vViolationTo4 * vViolationTo2 * vViolation;

        dPdX = 8.0 * uViolationTo7 * dUVdX + 8.0 * vViolationTo7 * dVVdX;
        dPdY = 8.0 * uViolationTo7 * dUVdY + 8.0 * vViolationTo7 * dVVdY;
      }


      template <class FloatType>
      FloatType
      ReverseProjectionObjective<FloatType>::
      computeFieldOfViewPenalty(
        const brick::numeric::Vector3D<FloatType>& candidate)
      {
        FloatType penalty(0.0);

        // Penalize unrectified points that should project outside the
        // image.  We make this term be zero at the edge of the field
        // of view, ramping quadratically to the cost equivalent of
        // being a whole image off when the azimuth exceeds the max by
        // 10%.

        FloatType azimuthTangent = brick::numeric::absoluteValue(
          candidate.x() / candidate.z());
        FloatType elevationTangent = brick::numeric::absoluteValue(
          candidate.y() / candidate.z());

        // Are we in violation horizontally?
        if(azimuthTangent > this->m_maxAzimuthTangent) {

          // This is how far out of bounds we are.
          FloatType azimuthViolation =
            azimuthTangent - this->m_maxAzimuthTangent;

          // We define our cost to have the form (k * azimuthViolation)^2,
          // where k is a scaling constant. We want to choose k so
          // that the cost is equal to imageWidth^2 when
          // azimuthVioltion is equal to 0.1 * maxAzimuthTangent.
          // Below, t_max is maxAzimuthTangent and w is imageWidth.
          // We write:
          //
          // @verbatim
          //   (k * 0.1 * t_max)^2 = w^2
          //   k * 0.1 * t_max = w
          //   k = 10 * w / t_max
          // @endverbatim
          //
          // We precomputed k in the constructor.  We called it
          // m_azimuthViolationScale.
          azimuthViolation *= this->m_azimuthViolationScale;
          azimuthViolation *= azimuthViolation;
          penalty += azimuthViolation;
        }

        // Are we in violation vertically?
        if(elevationTangent > this->m_maxElevationTangent) {
          // See the comment above regarding azimuthViolation.
          FloatType elevationViolation =
            elevationTangent - this->m_maxElevationTangent;
          elevationViolation *= this->m_elevationViolationScale;
          elevationViolation *= elevationViolation;
          penalty += elevationViolation;
        }

        return penalty;
      }


      template <class FloatType>
      void
      ReverseProjectionObjective<FloatType>::
      computeFieldOfViewPenaltyGradient(
        const FloatType& xValue, const FloatType& yValue,
        FloatType& dPdX, FloatType& dPdY)
      {
        // Initialize penalty gradients to zero.
        dPdX = FloatType(0.0);
        dPdY = FloatType(0.0);

        // Penalize unrectified points that should project outside the
        // image.
        FloatType azimuthTangent = xValue;
        FloatType elevationTangent = yValue;
        FloatType dAzdX(1.0);
        FloatType dEldY(1.0);

        // We do the absoluteValue operation explicitly.
        if(azimuthTangent < FloatType(0.0)) {
          azimuthTangent = -azimuthTangent;
          dAzdX = -dAzdX;
        }
        if(elevationTangent < FloatType(0.0)) {
          elevationTangent = -elevationTangent;
          dEldY = -dEldY;
        }

        // Are we in violation horizontally?
        if(azimuthTangent > this->m_maxAzimuthTangent) {
          FloatType azimuthViolation =
            ((azimuthTangent - this->m_maxAzimuthTangent)
             * this->m_azimuthViolationScale);
          dPdX = (FloatType(2) * azimuthViolation
                  * this->m_azimuthViolationScale * dAzdX);
        }

        // Are we in violation vertically?
        if(elevationTangent > this->m_maxElevationTangent) {
          FloatType elevationViolation =
            ((elevationTangent - this->m_maxElevationTangent)
             * this->m_elevationViolationScale);
          dPdY = (FloatType(2) * elevationViolation
                  * this->m_elevationViolationScale * dEldY);
        }
      }

    } // namespace privateCode;

  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_CAMERAINTRINSICSDISTORTEDPINHOLE_IMPL_HH */
