/**
***************************************************************************
* @file brick/computerVision/cameraIntrinsicsPlumbBob.cc
*
* Source file defining a CameraIntrinsics subclass for pinhole
* cameras.
*
* Copyright (C) 2007-2008 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#include <iomanip>
#include <brick/common/expect.hh>
#include <brick/computerVision/cameraIntrinsicsPlumbBob.hh>
#include <brick/optimization/optimizerBFGS.hh>
#include <brick/optimization/optimizerNelderMead.hh>

using namespace brick::numeric;
using namespace brick::geometry;

namespace {

  // We require reverse projection accuracy of better than 1/100th of
  // a pixel.
  const double l_maximumReverseProjectionResidual = 0.01 * 0.01;
  
}


namespace brick {

  namespace computerVision {


    // The default constructor initializes the CameraIntrinsicsPlumbBob
    // instance to a consistent (but not terribly useful) state.
    CameraIntrinsicsPlumbBob::
    CameraIntrinsicsPlumbBob()
      : CameraIntrinsicsDistorted(),
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
    CameraIntrinsicsPlumbBob::
    CameraIntrinsicsPlumbBob(size_t numPixelsX,
                             size_t numPixelsY,
                             double focalLengthX,
                             double focalLengthY,
                             double centerU,
                             double centerV,
                             double skewCoefficient,
                             double radialCoefficient0,
                             double radialCoefficient1,
                             double radialCoefficient2,
                             double tangentialCoefficient0,
                             double tangentialCoefficient1)
      : CameraIntrinsicsDistorted(),
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
    void
    CameraIntrinsicsPlumbBob::
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
    void
    CameraIntrinsicsPlumbBob::
    allowSkew(bool flag)
    {
      m_allowSkew = flag;
      if(flag == false) {
        m_skewCoefficient = 0.0;
      }
    }
    

    // This function exposes a subset of the intrinsic parameters
    // for use in calibration routines.
    CameraIntrinsicsPlumbBob::ParameterVectorType
    CameraIntrinsicsPlumbBob::
    getFreeParameters() const
    {
      size_t numParameters = 4;
      if(m_allowSixthOrderRadial) {
        ++numParameters;
      }
      if(m_allowSkew) {
        ++numParameters;
      }
      ParameterVectorType result(numParameters);
      size_t resultIndex = 0;
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
    CameraIntrinsicsPlumbBob::ParameterVectorType
    CameraIntrinsicsPlumbBob::
    getNominalFreeParameters() const
    {
      size_t numParameters = 4;
      if(m_allowSixthOrderRadial) {
        ++numParameters;
      }
      if(m_allowSkew) {
        ++numParameters;
      }
      ParameterVectorType result(numParameters);
      result = 0.0;
      return result;
    }
      
      
    // This member function takes a point in 3D camera coordinates
    // and projects it into pixel coordinates.
    Vector2D<double>
    CameraIntrinsicsPlumbBob::
    project(const Vector3D<double>& point) const
    {
      Vector3D<double> skewedPoint = this->projectThroughDistortion(point);
      return Vector2D<double>(m_kX * skewedPoint.x() + m_centerU,
                      m_kY * skewedPoint.y() + m_centerV);
    }
    

    // This member function sets the calibration from an input
    // stream.
    std::istream&
    CameraIntrinsicsPlumbBob::
    readFromStream(std::istream& stream)
    {
      // If stream is in a bad state, we can't read from it.
      if (!stream){
        return stream;
      }
    
      // We'll silently skip whitespace.
      common::Expect::FormatFlag flags = common::Expect::SkipWhitespace;

      // Read input data into temporary variables.
      double centerU;
      double centerV;
      double kX;
      double kY;
      size_t numpixelsX;
      size_t numpixelsY;
      double radialCoefficient0;
      double radialCoefficient1;
      double radialCoefficient2;
      double skewCoefficient;
      double tangentialCoefficient0;
      double tangentialCoefficient1;
      
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
    geometry::Ray3D<double>
    CameraIntrinsicsPlumbBob::
    reverseProject(const Vector2D<double>& pixelPosition,
                   bool normalize) const
    {
      // TODO(xxx): Collapse this optimization to 1D when tangential
      // distortion is zero.

      // Ignoring distortion, we calculate an initial guess using the
      // pinhole camera model.  See
      // CameraIntrinsicsPinhole::reverseProject() for more detailed
      // comments.
      //
      // Array1D<double> startPoint(2);
      // startPoint[0] = (pixelPosition.x() - m_centerU) / m_kX;
      // startPoint[1] = (pixelPosition.y() - m_centerV) / m_kY;

      // Actually, for extreme distortions, the estimate above can be
      // a point that projects quite far outside the image.  This can
      // cause numerical problems because we have an eighth power
      // out-of-bounds error term that blows up quickly.  Instead of
      // using the pinhole projection above, we choose the safest
      // start point we can think of... the center of projection.
      Array1D<double> startPoint(2);
      startPoint = 0.0;
      
      // Now optimize to find a better answer.
      privateCode::PlumbBobObjective objective(*this, pixelPosition);
      OptimizerBFGS<privateCode::PlumbBobObjective> optimizer(objective);
      optimizer.setStartPoint(startPoint);
      Array1D<double> endPoint = optimizer.optimum();
      double residual = optimizer.optimalValue();

      // This would be the place to check convergence and try a
      // different start point, if we were so inclined.
      
      if(residual
         > l_maximumReverseProjectionResidual + objective.getOffset()) {
        BRICK_THROW(brick::common::ValueException,
                    "CameraIntrinsicsPlumbBob::reverseProject()",
                    "Reverse projection failed to converge.");
      }
      
      return Ray3D<double>(Vector3D<double>(0.0, 0.0, 0.0),
                           Vector3D<double>(endPoint[0], endPoint[1], 1.0),
                           normalize);
    }


    // This function exposes a subset of the intrinsic parameters
    // for use in calibration routines.
    void
    CameraIntrinsicsPlumbBob::
    setFreeParameters(ParameterVectorType const& parameterVector)
    {
      size_t ii = 0;
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
    std::ostream&
    CameraIntrinsicsPlumbBob::
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


    void
    CameraIntrinsicsPlumbBob::
    projectWithPartialDerivatives(double xNorm,
                                  double yNorm,
                                  double& uValue,
                                  double& vValue,
                                  double& dUdX,
                                  double& dUdY,
                                  double& dVdX,
                                  double& dVdY) const
    {
      double xSquared = xNorm * xNorm;
      double ySquared = yNorm * yNorm;
      double rSquared = xSquared + ySquared;
      double rFourth = rSquared * rSquared;
      double rSixth = rSquared * rFourth;

      double dRSquaredDX = 2.0 * xNorm;
      double dRSquaredDY = 2.0 * yNorm;
      double dRFourthDX = 2.0 * rSquared * dRSquaredDX;
      double dRFourthDY = 2.0 * rSquared * dRSquaredDY;
      double dRSixthDX = 3.0 * rFourth * dRSquaredDX;
      double dRSixthDY = 3.0 * rFourth * dRSquaredDY;

      // Compute distortion terms according to plumb bob model.
      double radialDistortion = (1.0 + m_radialCoefficient0 * rSquared
                                 + m_radialCoefficient1 * rFourth
                                 + m_radialCoefficient2 * rSixth);
      double dRadialDistortionDX =
        (m_radialCoefficient0 * dRSquaredDX
         + m_radialCoefficient1 * dRFourthDX
         + m_radialCoefficient2 * dRSixthDX);
      double dRadialDistortionDY =
        (m_radialCoefficient0 * dRSquaredDY
         + m_radialCoefficient1 * dRFourthDY
         + m_radialCoefficient2 * dRSixthDY);

      double crossTerm = xNorm * yNorm;
      double dCrossTermDX = yNorm;
      double dCrossTermDY = xNorm;
      
      Vector2D<double> tangentialDistortion(
        (2.0 * m_tangentialCoefficient0 * crossTerm
         + m_tangentialCoefficient1 * (rSquared + 2.0 * xSquared)),
        (m_tangentialCoefficient0 * (rSquared + 2.0 * ySquared)
         + 2.0 * m_tangentialCoefficient1 * crossTerm));
      double dTangXDX =
        (2.0 * m_tangentialCoefficient0 * dCrossTermDX
         + m_tangentialCoefficient1 * (dRSquaredDX + 4.0 * xNorm));
      double dTangXDY =
        (2.0 * m_tangentialCoefficient0 * dCrossTermDY
         + m_tangentialCoefficient1 * dRSquaredDY);
      double dTangYDX =
        (m_tangentialCoefficient0 * dRSquaredDX
         + 2.0 * m_tangentialCoefficient1 * dCrossTermDX);
      double dTangYDY =
        (m_tangentialCoefficient0 * (dRSquaredDY + 4.0 * yNorm)
         + 2.0 * m_tangentialCoefficient1 * dCrossTermDY);

      // Apply distortion and skew, then project into pixel coordinates.
      Vector2D<double> distortedPoint(
        radialDistortion * xNorm + tangentialDistortion.x(),
        radialDistortion * yNorm + tangentialDistortion.y());

      double dDistortedXDX =
        dRadialDistortionDX * xNorm + radialDistortion + dTangXDX;
      double dDistortedXDY = dRadialDistortionDY * xNorm + dTangXDY;
      double dDistortedYDX = dRadialDistortionDX * yNorm + dTangYDX;
      double dDistortedYDY =
        dRadialDistortionDY * yNorm + radialDistortion + dTangYDY;
      
      Vector2D<double> skewedPoint(
        distortedPoint.x() + m_skewCoefficient * distortedPoint.y(),
        distortedPoint.y());
      double dSkewedXDX = dDistortedXDX + m_skewCoefficient * dDistortedYDX;
      double dSkewedXDY = dDistortedXDY + m_skewCoefficient * dDistortedYDY;
      double dSkewedYDX = dDistortedYDX;
      double dSkewedYDY = dDistortedYDY;
      
      uValue = m_kX * skewedPoint.x() + m_centerU;
      vValue = m_kY * skewedPoint.y() + m_centerV;
      dUdX = m_kX * dSkewedXDX;
      dUdY = m_kX * dSkewedXDY;
      dVdX = m_kY * dSkewedYDX;
      dVdY = m_kY * dSkewedYDY;
    }
    


    namespace privateCode {

      PlumbBobObjective::
      PlumbBobObjective(const CameraIntrinsicsPlumbBob& intrinsics,
                        const Vector2D<double>& uvTarget)
        : m_intrinsics(intrinsics),
          m_offset(0.0001),
          m_uvTarget(uvTarget)
      {
        // Empty.
      }

        
      double
      PlumbBobObjective::
      operator()(const Array1D<double>& theta)
      {
        Vector3D<double> candidate(theta[0], theta[1], 1.0);
        Vector2D<double> projection = m_intrinsics.project(candidate);

        // The high order terms of the distortion model can lead to
        // false minima outside the image boundaries.  We add an even
        // higher order penalty for projecting outside the image,
        // hopefully reducing this problem.
        double boundsPenalty = this->computeBoundsPenalty(projection);

        return (magnitudeSquared<double>(projection - m_uvTarget)
                + boundsPenalty + m_offset);
      }


      double
      PlumbBobObjective::
      getOffset()
      {
        return m_offset;
      }

      
      Array1D<double>
      PlumbBobObjective::
      gradient(const Array1D<double>& theta)
      {
        double uValue;
        double vValue;
        double dUdX;
        double dUdY;
        double dVdX;
        double dVdY;
        m_intrinsics.projectWithPartialDerivatives(
          theta[0], theta[1], uValue, vValue, dUdX, dUdY, dVdX, dVdY);

        double twoTimesDeltaU = 2.0 * (uValue - m_uvTarget.x());
        double twoTimesDeltaV = 2.0 * (vValue - m_uvTarget.y());

        // Using chain rule,
        //
        // (d/dx)(deltaU^2 + deltaV^2) =
        //    2*deltaU*(d/dx)deltaU + 2*deltaV*(d/dx)deltaV
        //
        // (d/dy)(deltaU^2 + deltaV^2) =
        //    2*deltaU*(d/dy)deltaU + 2*deltaV*(d/dy)deltaV 
        Array1D<double> gradientArray(2);
        gradientArray[0] = twoTimesDeltaU * dUdX + twoTimesDeltaV * dVdX;
        gradientArray[1] = twoTimesDeltaU * dUdY + twoTimesDeltaV * dVdY;

        // Don't forget to add gradient for bounds penalty.
        double dPdX;
        double dPdY;
        this->computeBoundsPenaltyGradient(
          uValue, vValue, dUdX, dUdY, dVdX, dVdY, dPdX, dPdY);
        gradientArray[0] += dPdX;
        gradientArray[1] += dPdY;
        
        return gradientArray;
      }


      double
      PlumbBobObjective::
      computeBoundsPenalty(const Vector2D<double>& uvPosition)
      {
        // Start by computing a vaguely normalized measure of how far
        // out-of-bounds uvPosition is.  We'll call this "violation."
        double uViolation = 0.0;
        double scaleU = m_intrinsics.getNumPixelsX() / 10.0;
        if(uvPosition.x() < 0) {
          uViolation = -uvPosition.x() / scaleU;
        } else if(uvPosition.x() >= m_intrinsics.getNumPixelsX()) {
          uViolation = ((uvPosition.x() - m_intrinsics.getNumPixelsX())
                        / scaleU);
        }
        double vViolation = 0.0;
        double scaleV = m_intrinsics.getNumPixelsY() / 10.0;
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
        

      void
      PlumbBobObjective::
      computeBoundsPenaltyGradient(double uValue, double vValue,
                                   double dUdX, double dUdY,
                                   double dVdX, double dVdY,
                                   double& dPdX, double& dPdY)
      {
        // Start by computing a vaguely normalized measure of how far
        // out-of-bounds uvPosition is.  We'll call this "violation."
        // Simultaneously compute the first derivative of violation
        // wrt x and y.
        double uViolation = 0.0;
        double scaleU = m_intrinsics.getNumPixelsX() / 10.0;
        double dUVdX = 0.0;  // <-- derivative of uViolation wrt X.
        double dUVdY = 0.0;  // <-- derivative of uViolation wrt Y.
        if(uValue < 0) {
          uViolation = -uValue / scaleU;
          dUVdX = -dUdX / scaleU;
          dUVdY = -dUdY / scaleU;
        } else if(uValue >= m_intrinsics.getNumPixelsX()) {
          uViolation = (uValue - m_intrinsics.getNumPixelsX()) / scaleU;
          dUVdX = dUdX / scaleU;
          dUVdY = dUdY / scaleU;
        }
        double vViolation = 0.0;
        double scaleV = m_intrinsics.getNumPixelsY() / 10.0;
        double dVVdX = 0.0;  // <-- derivative of vViolation wrt X.
        double dVVdY = 0.0;  // <-- derivative of vViolation wrt Y.
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
        double uViolationTo2 = uViolation * uViolation;
        double uViolationTo4 = uViolationTo2 * uViolationTo2;
        double uViolationTo7 = uViolationTo4 * uViolationTo2 * uViolation;
        double vViolationTo2 = vViolation * vViolation;
        double vViolationTo4 = vViolationTo2 * vViolationTo2;
        double vViolationTo7 = vViolationTo4 * vViolationTo2 * vViolation;

        dPdX = 8.0 * uViolationTo7 * dUVdX + 8.0 * vViolationTo7 * dVVdX;
        dPdY = 8.0 * uViolationTo7 * dUVdY + 8.0 * vViolationTo7 * dVVdY;
      }
        

    } // namespace privateCode;
    
  } // namespace computerVision
  
} // namespace brick
