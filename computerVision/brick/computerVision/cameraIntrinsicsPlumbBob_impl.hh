/**
***************************************************************************
* @file brick/computerVision/cameraIntrinsicsPlumbBob_impl.hh
*
* Definitions of inline and template functions for CameraIntrinsicsPlumbBob.
*
* Copyright (C) 2007-2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_CAMERAINTRINSICSPLUMBBOB_IMPL_HH
#define BRICK_COMPUTERVISION_CAMERAINTRINSICSPLUMBBOB_IMPL_HH

// This file is included by array1D.hh, and should not be directly included
// by user code, so no need to include cameraIntrinsicsPlumbBob.hh here.
// 
// #include <brick/numeric/cameraIntrinsicsPlumbBob.hh>

#include <iomanip>
#include <brick/common/expect.hh>
#include <brick/computerVision/cameraIntrinsicsPlumbBob.hh>
#include <brick/optimization/optimizerBFGS.hh>
#include <brick/optimization/optimizerNelderMead.hh>

namespace brick {

  namespace computerVision {

    namespace privateCode {

      // We require reverse projection accuracy of better than 1/100th of
      // a pixel.
      template <class FloatType>
      const FloatType&
      getPlumbBobMaximumReverseProjectionResidual() {
        static FloatType maxResidual(0.01 * 0.01);
        return maxResidual;
      }

      
      /**
       ** This functor is called by
       ** CameraIntrinsicsPlumbBob::reverseProject() during iterative
       ** approximation of reverse projection.
       **/
      template <class FloatType>
      class PlumbBobObjective
        : public std::unary_function<brick::numeric::Array1D<FloatType>, FloatType>
      {
      public:
        PlumbBobObjective(const CameraIntrinsicsPlumbBob<FloatType>& intrinsics,
                          const brick::numeric::Vector2D<FloatType>& uvTarget);
        
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

        CameraIntrinsicsPlumbBob<FloatType> m_intrinsics;
        FloatType m_offset;
        brick::numeric::Vector2D<FloatType> m_uvTarget;
      };

      // Implementation of PlumbBobObjective is at the bottom of this
      // file.
      
    } // namespace privateCode;


    // The default constructor initializes the CameraIntrinsicsPlumbBob
    // instance to a consistent (but not terribly useful) state.
    template <class FloatType>
    CameraIntrinsicsPlumbBob<FloatType>::
    CameraIntrinsicsPlumbBob()
      : CameraIntrinsicsDistorted<FloatType>(),
        m_allowSixthOrderRadial(true),
        m_allowSkew(true),
        m_centerU(50),
        m_centerV(50),
        m_kX(1.0),
        m_kY(1.0),
        m_numPixelsX(100),
        m_numPixelsY(100),
        m_radialCoefficient0(0.0),
        m_radialCoefficient1(0.0),
        m_radialCoefficient2(0.0),
        m_skewCoefficient(0.0),
        m_tangentialCoefficient0(0.0),
        m_tangentialCoefficient1(0.0)
    {
      // Empty.
    }
      

    // This constructor allows the caller to explicitly set the
    // camera intrinsic parameters.
    template <class FloatType>
    CameraIntrinsicsPlumbBob<FloatType>::
    CameraIntrinsicsPlumbBob(unsigned int numPixelsX,
                             unsigned int numPixelsY,
                             FloatType focalLengthX,
                             FloatType focalLengthY,
                             FloatType centerU,
                             FloatType centerV,
                             FloatType skewCoefficient,
                             FloatType radialCoefficient0,
                             FloatType radialCoefficient1,
                             FloatType radialCoefficient2,
                             FloatType tangentialCoefficient0,
                             FloatType tangentialCoefficient1)
      : CameraIntrinsicsDistorted<FloatType>(),
        m_allowSixthOrderRadial(true),
        m_allowSkew(true),
        m_centerU(centerU),
        m_centerV(centerV),
        m_kX(focalLengthX),
        m_kY(focalLengthY),
        m_numPixelsX(numPixelsX),
        m_numPixelsY(numPixelsY),
        m_radialCoefficient0(radialCoefficient0),
        m_radialCoefficient1(radialCoefficient1),
        m_radialCoefficient2(radialCoefficient2),
        m_skewCoefficient(skewCoefficient),
        m_tangentialCoefficient0(tangentialCoefficient0),
        m_tangentialCoefficient1(tangentialCoefficient1)
    {
      // Empty.
    }
    

    // Call this function before running calibration routines to
    // specify whether or not you want to allow a nonzero third
    // radial distortion coefficient.
    template <class FloatType>
    void
    CameraIntrinsicsPlumbBob<FloatType>::
    allowSixthOrderRadial(bool flag)
    {
      m_allowSixthOrderRadial = flag;
      if(flag == false) {
        m_radialCoefficient2 = 0.0;
      }
    }


    // Call this function before running calibration routines to
    // specify whether or not you want to allow a nonzero skew
    // coefficient.
    template <class FloatType>
    void
    CameraIntrinsicsPlumbBob<FloatType>::
    allowSkew(bool flag)
    {
      m_allowSkew = flag;
      if(flag == false) {
        m_skewCoefficient = 0.0;
      }
    }
    

    // This function exposes a subset of the intrinsic parameters
    // for use in calibration routines.
    template <class FloatType>
    typename CameraIntrinsicsPlumbBob<FloatType>::ParameterVectorType
    CameraIntrinsicsPlumbBob<FloatType>::
    getFreeParameters() const
    {
      unsigned int numParameters = 4;
      if(m_allowSixthOrderRadial) {
        ++numParameters;
      }
      if(m_allowSkew) {
        ++numParameters;
      }
      typename CameraIntrinsicsPlumbBob<FloatType>::ParameterVectorType
        result(numParameters);
      unsigned int resultIndex = 0;
      result[resultIndex] = m_radialCoefficient0; ++resultIndex;
      result[resultIndex] = m_radialCoefficient1; ++resultIndex;
      if(m_allowSixthOrderRadial) {
        result[resultIndex] = m_radialCoefficient2; ++resultIndex;
      }
      if(m_allowSkew) {
        result[resultIndex] = m_skewCoefficient; ++resultIndex;
      }
      result[resultIndex] = m_tangentialCoefficient0; ++resultIndex;
      result[resultIndex] = m_tangentialCoefficient1; ++resultIndex;
      return result;        
    }
      
      
    // This function provides a reasonable starting point for
    // intrinsic parameters that are generally estimated by
    // nonlinear optimization.
    template <class FloatType>
    typename CameraIntrinsicsPlumbBob<FloatType>::ParameterVectorType
    CameraIntrinsicsPlumbBob<FloatType>::
    getNominalFreeParameters() const
    {
      unsigned int numParameters = 4;
      if(m_allowSixthOrderRadial) {
        ++numParameters;
      }
      if(m_allowSkew) {
        ++numParameters;
      }
      typename CameraIntrinsicsPlumbBob<FloatType>::ParameterVectorType
        result(numParameters);
      result = 0.0;
      return result;
    }
      
      
    // This member function takes a point in 3D camera coordinates
    // and projects it into pixel coordinates.
    template <class FloatType>
    brick::numeric::Vector2D<FloatType>
    CameraIntrinsicsPlumbBob<FloatType>::
    project(const brick::numeric::Vector3D<FloatType>& point) const
    {
      brick::numeric::Vector3D<FloatType> skewedPoint = this->projectThroughDistortion(point);
      return brick::numeric::Vector2D<FloatType>(m_kX * skewedPoint.x() + m_centerU,
                      m_kY * skewedPoint.y() + m_centerV);
    }
    

    // This member function sets the calibration from an input
    // stream.
    template <class FloatType>
    std::istream&
    CameraIntrinsicsPlumbBob<FloatType>::
    readFromStream(std::istream& stream)
    {
      // If stream is in a bad state, we can't read from it.
      if (!stream){
        return stream;
      }
    
      // We'll silently skip whitespace.
      common::Expect::FormatFlag flags = common::Expect::SkipWhitespace;

      // Read input data into temporary variables.
      FloatType centerU;
      FloatType centerV;
      FloatType kX;
      FloatType kY;
      unsigned int numpixelsX;
      unsigned int numpixelsY;
      FloatType radialCoefficient0;
      FloatType radialCoefficient1;
      FloatType radialCoefficient2;
      FloatType skewCoefficient;
      FloatType tangentialCoefficient0;
      FloatType tangentialCoefficient1;
      
      stream >> common::Expect("CameraIntrinsicsPlumbBob", flags);
      stream >> common::Expect("{", flags);
      stream >> numpixelsX;
      stream >> common::Expect(",", flags);
      stream >> numpixelsY;
      stream >> common::Expect(",", flags);
      stream >> kX;
      stream >> common::Expect(",", flags);
      stream >> kY;
      stream >> common::Expect(",", flags);
      stream >> centerU;
      stream >> common::Expect(",", flags);
      stream >> centerV;
      stream >> common::Expect(",", flags);
      stream >> skewCoefficient;
      stream >> common::Expect(",", flags);
      stream >> radialCoefficient0;
      stream >> common::Expect(",", flags);
      stream >> radialCoefficient1;
      stream >> common::Expect(",", flags);
      stream >> radialCoefficient2;
      stream >> common::Expect(",", flags);
      stream >> tangentialCoefficient0;
      stream >> common::Expect(",", flags);
      stream >> tangentialCoefficient1;
      stream >> common::Expect("}", flags);

      if(stream) {
        m_centerU = centerU;
        m_centerV = centerV;
        m_kX = kX;
        m_kY = kY;
        m_numPixelsX = numpixelsX;
        m_numPixelsY = numpixelsY;
        m_radialCoefficient0 = radialCoefficient0;
        m_radialCoefficient1 = radialCoefficient1;
        m_radialCoefficient2 = radialCoefficient2;
        m_skewCoefficient = skewCoefficient;
        m_tangentialCoefficient0 = tangentialCoefficient0;
        m_tangentialCoefficient1 = tangentialCoefficient1;
      }
      return stream;
    }


    // This member function takes a point in 2D pixel coordinates
    // and returns a ray in 3D camera coordinates passing through
    // all of the 3D points that project to the specified 2D
    // position.
    template <class FloatType>
    geometry::Ray3D<FloatType>
    CameraIntrinsicsPlumbBob<FloatType>::
    reverseProject(const brick::numeric::Vector2D<FloatType>& pixelPosition,
                   bool normalize) const
    {
      // TODO(xxx): Collapse this optimization to 1D when tangential
      // distortion is zero.

      // Ignoring distortion, we calculate an initial guess using the
      // pinhole camera model.  See
      // CameraIntrinsicsPinhole::reverseProject() for more detailed
      // comments.
      //
      // Array1D<FloatType> startPoint(2);
      // startPoint[0] = (pixelPosition.x() - m_centerU) / m_kX;
      // startPoint[1] = (pixelPosition.y() - m_centerV) / m_kY;

      // Actually, for extreme distortions, the estimate above can be
      // a point that projects quite far outside the image.  This can
      // cause numerical problems because we have an eighth power
      // out-of-bounds error term that blows up quickly.  Instead of
      // using the pinhole projection above, we choose the safest
      // start point we can think of... the center of projection.
      brick::numeric::Array1D<FloatType> startPoint(2);
      startPoint = 0.0;
      
      // Now optimize to find a better answer.
      typedef privateCode::PlumbBobObjective<FloatType> LocalObjective;
      LocalObjective objective(*this, pixelPosition);
      OptimizerBFGS<LocalObjective> optimizer(objective);
      optimizer.setStartPoint(startPoint);
      brick::numeric::Array1D<FloatType> endPoint = optimizer.optimum();
      FloatType residual = optimizer.optimalValue();

      // This would be the place to check convergence and try a
      // different start point, if we were so inclined.
      
      if(residual >
         (privateCode::getPlumbBobMaximumReverseProjectionResidual<FloatType>()
          + objective.getOffset())) {
        BRICK_THROW(brick::common::ValueException,
                    "CameraIntrinsicsPlumbBob<FloatType>::reverseProject()",
                    "Reverse projection failed to converge.");
      }
      
      return brick::geometry::Ray3D<FloatType>(
        brick::numeric::Vector3D<FloatType>(0.0, 0.0, 0.0),
        brick::numeric::Vector3D<FloatType>(endPoint[0], endPoint[1], 1.0),
        normalize);
    }


    // This function exposes a subset of the intrinsic parameters
    // for use in calibration routines.
    template <class FloatType>
    void
    CameraIntrinsicsPlumbBob<FloatType>::
    setFreeParameters(
      typename CameraIntrinsicsPlumbBob<FloatType>::ParameterVectorType
      const& parameterVector)
    {
      unsigned int ii = 0;
      m_radialCoefficient0 = parameterVector[ii]; ++ii;
      m_radialCoefficient1 = parameterVector[ii]; ++ii;
      if(m_allowSixthOrderRadial) {
        m_radialCoefficient2 = parameterVector[ii]; ++ii;
      } else {
        m_radialCoefficient2 = 0.0;
      }
      if(m_allowSkew) {
        m_skewCoefficient = parameterVector[ii]; ++ii;
      } else {
        m_skewCoefficient = 0.0;
      }
      m_tangentialCoefficient0 = parameterVector[ii]; ++ii;
      m_tangentialCoefficient1 = parameterVector[ii]; ++ii;
    }

    
    // This member function writes the calibration to an
    // outputstream in a format which is compatible with member
    // function readFromStream().
    template <class FloatType>
    std::ostream&
    CameraIntrinsicsPlumbBob<FloatType>::
    writeToStream(std::ostream& stream) const
    {
      std::ios::fmtflags flags = stream.flags();
      std::streamsize precision = stream.precision();
      stream.precision(15);
      stream << "CameraIntrinsicsPlumbBob {"
             << std::fixed << std::setw(15)
             << m_numPixelsX << ", "
             << m_numPixelsY << ", "
             << m_kX << ", "
             << m_kY << ", "
             << m_centerU << ", "
             << m_centerV << ", "
             << m_skewCoefficient << ", "
             << m_radialCoefficient0 << ", "
             << m_radialCoefficient1 << ", "
             << m_radialCoefficient2 << ", "
             << m_tangentialCoefficient0 << ", "
             << m_tangentialCoefficient1 << "}";
      stream.precision(precision);
      stream.flags(flags);
      return stream;
    }


    template <class FloatType>
    void
    CameraIntrinsicsPlumbBob<FloatType>::
    projectWithPartialDerivatives(FloatType xNorm,
                                  FloatType yNorm,
                                  FloatType& uValue,
                                  FloatType& vValue,
                                  FloatType& dUdX,
                                  FloatType& dUdY,
                                  FloatType& dVdX,
                                  FloatType& dVdY) const
    {
      FloatType xSquared = xNorm * xNorm;
      FloatType ySquared = yNorm * yNorm;
      FloatType rSquared = xSquared + ySquared;
      FloatType rFourth = rSquared * rSquared;
      FloatType rSixth = rSquared * rFourth;

      FloatType dRSquaredDX = 2.0 * xNorm;
      FloatType dRSquaredDY = 2.0 * yNorm;
      FloatType dRFourthDX = 2.0 * rSquared * dRSquaredDX;
      FloatType dRFourthDY = 2.0 * rSquared * dRSquaredDY;
      FloatType dRSixthDX = 3.0 * rFourth * dRSquaredDX;
      FloatType dRSixthDY = 3.0 * rFourth * dRSquaredDY;

      // Compute distortion terms according to plumb bob model.
      FloatType radialDistortion = (1.0 + m_radialCoefficient0 * rSquared
                                 + m_radialCoefficient1 * rFourth
                                 + m_radialCoefficient2 * rSixth);
      FloatType dRadialDistortionDX =
        (m_radialCoefficient0 * dRSquaredDX
         + m_radialCoefficient1 * dRFourthDX
         + m_radialCoefficient2 * dRSixthDX);
      FloatType dRadialDistortionDY =
        (m_radialCoefficient0 * dRSquaredDY
         + m_radialCoefficient1 * dRFourthDY
         + m_radialCoefficient2 * dRSixthDY);

      FloatType crossTerm = xNorm * yNorm;
      FloatType dCrossTermDX = yNorm;
      FloatType dCrossTermDY = xNorm;
      
      brick::numeric::Vector2D<FloatType> tangentialDistortion(
        (2.0 * m_tangentialCoefficient0 * crossTerm
         + m_tangentialCoefficient1 * (rSquared + 2.0 * xSquared)),
        (m_tangentialCoefficient0 * (rSquared + 2.0 * ySquared)
         + 2.0 * m_tangentialCoefficient1 * crossTerm));
      FloatType dTangXDX =
        (2.0 * m_tangentialCoefficient0 * dCrossTermDX
         + m_tangentialCoefficient1 * (dRSquaredDX + 4.0 * xNorm));
      FloatType dTangXDY =
        (2.0 * m_tangentialCoefficient0 * dCrossTermDY
         + m_tangentialCoefficient1 * dRSquaredDY);
      FloatType dTangYDX =
        (m_tangentialCoefficient0 * dRSquaredDX
         + 2.0 * m_tangentialCoefficient1 * dCrossTermDX);
      FloatType dTangYDY =
        (m_tangentialCoefficient0 * (dRSquaredDY + 4.0 * yNorm)
         + 2.0 * m_tangentialCoefficient1 * dCrossTermDY);

      // Apply distortion and skew, then project into pixel coordinates.
      brick::numeric::Vector2D<FloatType> distortedPoint(
        radialDistortion * xNorm + tangentialDistortion.x(),
        radialDistortion * yNorm + tangentialDistortion.y());

      FloatType dDistortedXDX =
        dRadialDistortionDX * xNorm + radialDistortion + dTangXDX;
      FloatType dDistortedXDY = dRadialDistortionDY * xNorm + dTangXDY;
      FloatType dDistortedYDX = dRadialDistortionDX * yNorm + dTangYDX;
      FloatType dDistortedYDY =
        dRadialDistortionDY * yNorm + radialDistortion + dTangYDY;
      
      brick::numeric::Vector2D<FloatType> skewedPoint(
        distortedPoint.x() + m_skewCoefficient * distortedPoint.y(),
        distortedPoint.y());
      FloatType dSkewedXDX = dDistortedXDX + m_skewCoefficient * dDistortedYDX;
      FloatType dSkewedXDY = dDistortedXDY + m_skewCoefficient * dDistortedYDY;
      FloatType dSkewedYDX = dDistortedYDX;
      FloatType dSkewedYDY = dDistortedYDY;
      
      uValue = m_kX * skewedPoint.x() + m_centerU;
      vValue = m_kY * skewedPoint.y() + m_centerV;
      dUdX = m_kX * dSkewedXDX;
      dUdY = m_kX * dSkewedXDY;
      dVdX = m_kY * dSkewedYDX;
      dVdY = m_kY * dSkewedYDY;
    }
    
    
    
    template <class FloatType>
    inline brick::numeric::Vector3D<FloatType>
    CameraIntrinsicsPlumbBob<FloatType>::
    projectThroughDistortion(const brick::numeric::Vector3D<FloatType>& point)
      const
    {
      // The Vector2D constructor will divide by point.z(), which is will
      // throw if point.z() == 0.0.
      brick::numeric::Vector2D<FloatType> normalizedPoint(point.x(), point.y(), point.z());
    
      FloatType xSquared = normalizedPoint.x() * normalizedPoint.x();
      FloatType ySquared = normalizedPoint.y() * normalizedPoint.y();
      FloatType rSquared = xSquared + ySquared;
      FloatType rFourth = rSquared * rSquared;
      FloatType rSixth = rSquared * rFourth;
    
      // Compute distortion terms according to plumb bob model.
      FloatType radialDistortion = (1.0 + m_radialCoefficient0 * rSquared
                                 + m_radialCoefficient1 * rFourth
                                 + m_radialCoefficient2 * rSixth);
    
      FloatType crossTerm = normalizedPoint.x() * normalizedPoint.y();
      brick::numeric::Vector2D<FloatType> tangentialDistortion(
        (2.0 * m_tangentialCoefficient0 * crossTerm
         + m_tangentialCoefficient1 * (rSquared + 2.0 * xSquared)),
        (m_tangentialCoefficient0 * (rSquared + 2.0 * ySquared)
         + 2.0 * m_tangentialCoefficient1 * crossTerm));
    
      // Apply distortion and skew, then project into pixel coordinates.
      brick::numeric::Vector2D<FloatType> distortedPoint(
        radialDistortion * normalizedPoint + tangentialDistortion);
      return brick::numeric::Vector3D<FloatType>(
        distortedPoint.x() + m_skewCoefficient * distortedPoint.y(),
        distortedPoint.y(),
        1.0);
    }



    // Implementation of PlumbBobObjective.
    namespace privateCode {

      template <class FloatType>
      PlumbBobObjective<FloatType>::
      PlumbBobObjective(const CameraIntrinsicsPlumbBob<FloatType>& intrinsics,
                        const brick::numeric::Vector2D<FloatType>& uvTarget)
        : m_intrinsics(intrinsics),
          m_offset(0.0001),
          m_uvTarget(uvTarget)
      {
        // Empty.
      }

        
      template <class FloatType>
      FloatType
      PlumbBobObjective<FloatType>::
      operator()(const brick::numeric::Array1D<FloatType>& theta)
      {
        brick::numeric::Vector3D<FloatType> candidate(theta[0], theta[1], 1.0);
        brick::numeric::Vector2D<FloatType> projection = m_intrinsics.project(candidate);

        // The high order terms of the distortion model can lead to
        // false minima outside the image boundaries.  We add an even
        // higher order penalty for projecting outside the image,
        // hopefully reducing this problem.
        FloatType boundsPenalty = this->computeBoundsPenalty(projection);

        return (brick::numeric::magnitudeSquared<FloatType>(
                  projection - m_uvTarget) + boundsPenalty + m_offset);
      }


      template <class FloatType>
      FloatType
      PlumbBobObjective<FloatType>::
      getOffset()
      {
        return m_offset;
      }

      
      template <class FloatType>
      brick::numeric::Array1D<FloatType>
      PlumbBobObjective<FloatType>::
      gradient(const brick::numeric::Array1D<FloatType>& theta)
      {
        FloatType uValue;
        FloatType vValue;
        FloatType dUdX;
        FloatType dUdY;
        FloatType dVdX;
        FloatType dVdY;
        m_intrinsics.projectWithPartialDerivatives(
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
        brick::numeric::Array1D<FloatType> gradientArray(2);
        gradientArray[0] = twoTimesDeltaU * dUdX + twoTimesDeltaV * dVdX;
        gradientArray[1] = twoTimesDeltaU * dUdY + twoTimesDeltaV * dVdY;

        // Don't forget to add gradient for bounds penalty.
        FloatType dPdX;
        FloatType dPdY;
        this->computeBoundsPenaltyGradient(
          uValue, vValue, dUdX, dUdY, dVdX, dVdY, dPdX, dPdY);
        gradientArray[0] += dPdX;
        gradientArray[1] += dPdY;
        
        return gradientArray;
      }


      template <class FloatType>
      FloatType
      PlumbBobObjective<FloatType>::
      computeBoundsPenalty(const brick::numeric::Vector2D<FloatType>& uvPosition)
      {
        // Start by computing a vaguely normalized measure of how far
        // out-of-bounds uvPosition is.  We'll call this "violation."
        FloatType uViolation = 0.0;
        FloatType scaleU = m_intrinsics.getNumPixelsX() / 10.0;
        if(uvPosition.x() < 0) {
          uViolation = -uvPosition.x() / scaleU;
        } else if(uvPosition.x() >= m_intrinsics.getNumPixelsX()) {
          uViolation = ((uvPosition.x() - m_intrinsics.getNumPixelsX())
                        / scaleU);
        }
        FloatType vViolation = 0.0;
        FloatType scaleV = m_intrinsics.getNumPixelsY() / 10.0;
        if(uvPosition.y() < 0) {
          vViolation = -uvPosition.y() / scaleV;
        } else if(uvPosition.y() >= m_intrinsics.getNumPixelsY()) {
          vViolation = ((uvPosition.y() - m_intrinsics.getNumPixelsY())
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
      PlumbBobObjective<FloatType>::
      computeBoundsPenaltyGradient(FloatType uValue, FloatType vValue,
                                   FloatType dUdX, FloatType dUdY,
                                   FloatType dVdX, FloatType dVdY,
                                   FloatType& dPdX, FloatType& dPdY)
      {
        // Start by computing a vaguely normalized measure of how far
        // out-of-bounds uvPosition is.  We'll call this "violation."
        // Simultaneously compute the first derivative of violation
        // wrt x and y.
        FloatType uViolation = 0.0;
        FloatType scaleU = m_intrinsics.getNumPixelsX() / 10.0;
        FloatType dUVdX = 0.0;  // <-- derivative of uViolation wrt X.
        FloatType dUVdY = 0.0;  // <-- derivative of uViolation wrt Y.
        if(uValue < 0) {
          uViolation = -uValue / scaleU;
          dUVdX = -dUdX / scaleU;
          dUVdY = -dUdY / scaleU;
        } else if(uValue >= m_intrinsics.getNumPixelsX()) {
          uViolation = (uValue - m_intrinsics.getNumPixelsX()) / scaleU;
          dUVdX = dUdX / scaleU;
          dUVdY = dUdY / scaleU;
        }
        FloatType vViolation = 0.0;
        FloatType scaleV = m_intrinsics.getNumPixelsY() / 10.0;
        FloatType dVVdX = 0.0;  // <-- derivative of vViolation wrt X.
        FloatType dVVdY = 0.0;  // <-- derivative of vViolation wrt Y.
        if(vValue < 0) {
          vViolation = -vValue / scaleV;
          dVVdX = -dVdX / scaleV;
          dVVdY = -dVdY / scaleV;
        } else if(vValue >= m_intrinsics.getNumPixelsY()) {
          vViolation = (vValue - m_intrinsics.getNumPixelsY()) / scaleV;
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
        

    } // namespace privateCode;
    
  } // namespace computerVision
  
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_CAMERAINTRINSICSPLUMBBOB_IMPL_HH */
