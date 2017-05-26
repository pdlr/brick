/**
***************************************************************************
* @file brick/computerVision/cameraIntrinsicsRational_impl.hh
*
* Definitions of inline and template functions for CameraIntrinsicsRational.
*
* Copyright (C) 2014 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_CAMERAINTRINSICSRATIONAL_IMPL_HH
#define BRICK_COMPUTERVISION_CAMERAINTRINSICSRATIONAL_IMPL_HH

// This file is included by cameraIntrinsicsRational.hh, and should
// not be directly included by user code, so no need to include
// cameraIntrinsicsRational.hh here.
// 
// #include <brick/numeric/cameraIntrinsicsRational.hh>

#include <iomanip>
#include <brick/common/expect.hh>
#include <brick/common/types.hh>
#include <brick/computerVision/cameraIntrinsicsRational.hh>
#include <brick/numeric/differentiableScalar.hh>

namespace brick {

  namespace computerVision {

    // The default constructor initializes the CameraIntrinsicsRational
    // instance to a consistent (but not terribly useful) state.
    template <class FloatType>
    CameraIntrinsicsRational<FloatType>::
    CameraIntrinsicsRational()
      : CameraIntrinsicsDistortedPinhole<FloatType>(),
        m_radialCoefficient0(0.0),
        m_radialCoefficient1(0.0),
        m_radialCoefficient2(0.0),
        m_radialCoefficient3(0.0),
        m_radialCoefficient4(0.0),
        m_radialCoefficient5(0.0),
        m_tangentialCoefficient0(0.0),
        m_tangentialCoefficient1(0.0)
    {
      // Empty.
    }
      

    // This constructor allows the caller to explicitly set the
    // camera intrinsic parameters.
    template <class FloatType>
    CameraIntrinsicsRational<FloatType>::
    CameraIntrinsicsRational(unsigned int numPixelsX,
                             unsigned int numPixelsY,
                             FloatType focalLengthX,
                             FloatType focalLengthY,
                             FloatType centerU,
                             FloatType centerV,
                             FloatType radialCoefficient0,
                             FloatType radialCoefficient1,
                             FloatType radialCoefficient2,
                             FloatType radialCoefficient3,
                             FloatType radialCoefficient4,
                             FloatType radialCoefficient5,
                             FloatType tangentialCoefficient0,
                             FloatType tangentialCoefficient1)
      : CameraIntrinsicsDistortedPinhole<FloatType>(numPixelsX, numPixelsY,
                                                    focalLengthX, focalLengthY,
                                                    centerU, centerV),
        m_radialCoefficient0(radialCoefficient0),
        m_radialCoefficient1(radialCoefficient1),
        m_radialCoefficient2(radialCoefficient2),
        m_radialCoefficient3(radialCoefficient3),
        m_radialCoefficient4(radialCoefficient4),
        m_radialCoefficient5(radialCoefficient5),
        m_tangentialCoefficient0(tangentialCoefficient0),
        m_tangentialCoefficient1(tangentialCoefficient1)
    {
      // Empty.
    }
    

    // This function exposes the distortion parameters of the camera
    // model.
    template <class FloatType>
    typename CameraIntrinsicsRational<FloatType>::ParameterVectorType
    CameraIntrinsicsRational<FloatType>::
    getDistortionCoefficients() const
    {
      brick::common::UInt32 const numParameters = 8;
      typename CameraIntrinsicsRational<FloatType>::ParameterVectorType
        result(numParameters);
      unsigned int resultIndex = 0;
      result[resultIndex] = m_radialCoefficient0; ++resultIndex;
      result[resultIndex] = m_radialCoefficient1; ++resultIndex;
      result[resultIndex] = m_radialCoefficient2; ++resultIndex;
      result[resultIndex] = m_radialCoefficient3; ++resultIndex;
      result[resultIndex] = m_radialCoefficient4; ++resultIndex;
      result[resultIndex] = m_radialCoefficient5; ++resultIndex;
      result[resultIndex] = m_tangentialCoefficient0; ++resultIndex;
      result[resultIndex] = m_tangentialCoefficient1; ++resultIndex;
      if(resultIndex != numParameters) {
        BRICK_THROW(brick::common::LogicException,
                    "CameraIntrinsicsRational::getDistortionCoefficients()",
                    "Wrong number of parameters.");
      }
      return result;        
    }
      
      
    // This function provides a reasonable starting point for
    // intrinsic parameters that are generally estimated by
    // nonlinear optimization.
    template <class FloatType>
    typename CameraIntrinsicsRational<FloatType>::ParameterVectorType
    CameraIntrinsicsRational<FloatType>::
    getNominalFreeParameters() const
    {
      brick::common::UInt32 const numParameters = 8;
      typename CameraIntrinsicsRational<FloatType>::ParameterVectorType
        result(numParameters);
      result = 0.0;
      return result;
    }
      

    // Returns a vector of all parameters of the class.
    template <class FloatType>
    typename CameraIntrinsicsRational<FloatType>::ParameterVectorType
    CameraIntrinsicsRational<FloatType>::
    getParameters() const
    {
      brick::common::UInt32 constexpr numParameters = 4 + 8;
      typename CameraIntrinsicsRational<FloatType>::ParameterVectorType
        result(numParameters);
      unsigned int resultIndex = 0;
      result[resultIndex] = this->getFocalLengthX(); ++resultIndex;
      result[resultIndex] = this->getFocalLengthY(); ++resultIndex;
      result[resultIndex] = this->getCenterU(); ++resultIndex;
      result[resultIndex] = this->getCenterV(); ++resultIndex;
      result[resultIndex] = m_radialCoefficient0; ++resultIndex;
      result[resultIndex] = m_radialCoefficient1; ++resultIndex;
      result[resultIndex] = m_radialCoefficient2; ++resultIndex;
      result[resultIndex] = m_radialCoefficient3; ++resultIndex;
      result[resultIndex] = m_radialCoefficient4; ++resultIndex;
      result[resultIndex] = m_radialCoefficient5; ++resultIndex;
      result[resultIndex] = m_tangentialCoefficient0; ++resultIndex;
      result[resultIndex] = m_tangentialCoefficient1; ++resultIndex;
      if(resultIndex != numParameters) {
        BRICK_THROW(brick::common::LogicException,
                    "CameraIntrinsicsRational::getFreeParameters()",
                    "Wrong number of parameters.");
      }
      return result;        
    }
    

    // This member function takes a point in 3D camera coordinates
    // and projects it into pixel coordinates.
    template <class FloatType>
    brick::numeric::Vector2D<FloatType>
    CameraIntrinsicsRational<FloatType>::
    project(const brick::numeric::Vector3D<FloatType>& point) const
    {
      brick::numeric::Vector3D<FloatType> distortedPoint =
        this->projectThroughDistortion(point);
      return brick::numeric::Vector2D<FloatType>(
        this->getFocalLengthX() * distortedPoint.x() + this->getCenterU(),
        this->getFocalLengthY() * distortedPoint.y() + this->getCenterV());
    }
    

    // This member function sets the calibration from an input
    // stream.
    template <class FloatType>
    std::istream&
    CameraIntrinsicsRational<FloatType>::
    readFromStream(std::istream& stream)
    {
      // If stream is in a bad state, we can't read from it.
      if (!stream){
        return stream;
      }
    
      // We'll silently skip whitespace.
      common::Expect::FormatFlag flags = common::Expect::SkipWhitespace();

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
      FloatType radialCoefficient3;
      FloatType radialCoefficient4;
      FloatType radialCoefficient5;
      FloatType tangentialCoefficient0;
      FloatType tangentialCoefficient1;
      
      stream >> common::Expect("CameraIntrinsicsRational", flags);
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
      stream >> radialCoefficient0;
      stream >> common::Expect(",", flags);
      stream >> radialCoefficient1;
      stream >> common::Expect(",", flags);
      stream >> radialCoefficient2;
      stream >> common::Expect(",", flags);
      stream >> radialCoefficient3;
      stream >> common::Expect(",", flags);
      stream >> radialCoefficient4;
      stream >> common::Expect(",", flags);
      stream >> radialCoefficient5;
      stream >> common::Expect(",", flags);
      stream >> tangentialCoefficient0;
      stream >> common::Expect(",", flags);
      stream >> tangentialCoefficient1;
      stream >> common::Expect("}", flags);

      if(stream) {
        this->setDependentParameters(
          numpixelsX, numpixelsY, kX, kY, centerU, centerV);
        m_radialCoefficient0 = radialCoefficient0;
        m_radialCoefficient1 = radialCoefficient1;
        m_radialCoefficient2 = radialCoefficient2;
        m_radialCoefficient3 = radialCoefficient3;
        m_radialCoefficient4 = radialCoefficient4;
        m_radialCoefficient5 = radialCoefficient5;
        m_tangentialCoefficient0 = tangentialCoefficient0;
        m_tangentialCoefficient1 = tangentialCoefficient1;
      }
      return stream;
    }


    // This function iteratively computes and returns a ray in 3D
    // camera coordinates starting at the camera focus and passing
    // through the specified pixel position.
    template <class FloatType>
    geometry::Ray3D<FloatType>
    CameraIntrinsicsRational<FloatType>::
    reverseProjectEM(
      const brick::numeric::Vector2D<FloatType>& pixelPosition,
      bool normalize,
      FloatType requiredPrecision,
      std::size_t maximumIterations,
      std::size_t minimumIterations) const
    {
      // We'll want this later to avoid computing a bunch of square roots.
      FloatType requiredPrecisionSquared =
        requiredPrecision * requiredPrecision;
      
      // First project back through pinhole parameters to get a point
      // we can reverse project through the distortion model.  See
      // CameraIntrinsicsPinhole::reverseProject() for an explanation
      // of this line.
      brick::numeric::Vector2D<FloatType> x0(
        (pixelPosition.x() - this->getCenterU()) / this->getFocalLengthX(),
        (pixelPosition.y() - this->getCenterV()) / this->getFocalLengthY());

      // Iterate to find the undistorted point.
      brick::numeric::Vector2D<FloatType> xHat = x0;
      std::size_t ii = 1;
      while(1) {
        // We have x0 = a(x) * x + b(x),
        // where a(x) is radial distortion, and b(x) is tangential distortion.
        // 
        // Approximate that as x0 ~= a(xHat) * x + b(xHat)
        // Gives us x ~= (x0 - b(xHat)) / a(xHat)
        // 
        // Iterate on this.  This converges as long as the gradient of
        // the distortion field has magnitude less than 1.0.
    
        FloatType xSquared = xHat.x() * xHat.x();
        FloatType ySquared = xHat.y() * xHat.y();
        FloatType rSquared = xSquared + ySquared;
        FloatType rFourth = rSquared * rSquared;
        FloatType rSixth = rSquared * rFourth;
    
        // Compute radial distortion terms.
        FloatType radialDistortionNumerator =
          (1.0 + m_radialCoefficient0 * rSquared
           + m_radialCoefficient1 * rFourth
           + m_radialCoefficient2 * rSixth);
        FloatType radialDistortionDenominator = 
          (1.0 + m_radialCoefficient3 * rSquared
           + m_radialCoefficient4 * rFourth
           + m_radialCoefficient5 * rSixth);
        FloatType oneOverRadialDistortion =
          radialDistortionDenominator/ radialDistortionNumerator;

        // Compute tangential distortion.
        FloatType crossTerm = xHat.x() * xHat.y();
        brick::numeric::Vector2D<FloatType> tangentialDistortion(
          (2.0 * m_tangentialCoefficient0 * crossTerm
           + m_tangentialCoefficient1 * (rSquared + 2.0 * xSquared)),
          (m_tangentialCoefficient0 * (rSquared + 2.0 * ySquared)
           + 2.0 * m_tangentialCoefficient1 * crossTerm));
    
        // Apply distortion.
        brick::numeric::Vector2D<FloatType> xNext =
          (x0 - tangentialDistortion) * oneOverRadialDistortion;

        // Check for termination criteria.
        if(ii >= minimumIterations) {
          FloatType incrementSquared =
            brick::numeric::magnitudeSquared<FloatType>(xNext - xHat);
          if(incrementSquared < requiredPrecisionSquared) {
            xHat = xNext;
            break;
          }
        }
        if(ii > maximumIterations) {
          BRICK_THROW(
            brick::common::ValueException,
            "CameraIntrinsicsRational<FloatType>::reverseProjectEM()",
            "Reverse projection failed to converge.");
        }

        // Prepare for the next iteration.
        xHat = xNext;
        ++ii;
      }

      return geometry::Ray3D<FloatType>(
        brick::numeric::Vector3D<FloatType>(0.0, 0.0, 0.0),
        brick::numeric::Vector3D<FloatType>(xHat.x(), xHat.y(), 1.0),
        normalize);
    }

    // This sets the value of a subset of the intrinsic parameters,
    // and is commonly used by in calibration routines.
    template <class FloatType>
    void
    CameraIntrinsicsRational<FloatType>::
    setFreeParameters(
      typename CameraIntrinsicsRational<FloatType>::ParameterVectorType
      const& parameterVector)
    {
      unsigned int ii = 0;
      m_radialCoefficient0 = parameterVector[ii]; ++ii;
      m_radialCoefficient1 = parameterVector[ii]; ++ii;
      m_radialCoefficient2 = parameterVector[ii]; ++ii;
      m_radialCoefficient3 = parameterVector[ii]; ++ii;
      m_radialCoefficient4 = parameterVector[ii]; ++ii;
      m_radialCoefficient5 = parameterVector[ii]; ++ii;
      m_tangentialCoefficient0 = parameterVector[ii]; ++ii;
      m_tangentialCoefficient1 = parameterVector[ii]; ++ii;
    }


    // Sets the internal state of *this based on a parameter vector,
    // such as the one described in member function getParameters().
    template <class FloatType>
    void
    CameraIntrinsicsRational<FloatType>::
    setParameters(
      typename CameraIntrinsicsRational<FloatType>::ParameterVectorType
      const& parameterVector)
    {
      brick::common::UInt32 constexpr numParameters = 4 + 8;

      unsigned int ii = 0;
      this->setFocalLengthX(parameterVector[ii]); ++ii;
      this->setFocalLengthY(parameterVector[ii]); ++ii;
      this->setCenterU(parameterVector[ii]); ++ii;
      this->setCenterV(parameterVector[ii]); ++ii;
      m_radialCoefficient0 = parameterVector[ii]; ++ii;
      m_radialCoefficient1 = parameterVector[ii]; ++ii;
      m_radialCoefficient2 = parameterVector[ii]; ++ii;
      m_radialCoefficient3 = parameterVector[ii]; ++ii;
      m_radialCoefficient4 = parameterVector[ii]; ++ii;
      m_radialCoefficient5 = parameterVector[ii]; ++ii;
      m_tangentialCoefficient0 = parameterVector[ii]; ++ii;
      m_tangentialCoefficient1 = parameterVector[ii]; ++ii;
      if(ii != numParameters) {
        BRICK_THROW(brick::common::LogicException,
                    "CameraIntrinsicsRational::setParameters()",
                    "Wrong number of parameters.");
      }
    }
    
    
    // This member function writes the calibration to an
    // outputstream in a format that is compatible with member
    // function readFromStream().
    template <class FloatType>
    std::ostream&
    CameraIntrinsicsRational<FloatType>::
    writeToStream(std::ostream& stream) const
    {
      std::ios::fmtflags flags = stream.flags();
      std::streamsize precision = stream.precision();
      stream.precision(15);
      stream << "CameraIntrinsicsRational {"
             << std::fixed << std::setw(15)
             << this->getNumPixelsX() << ", "
             << this->getNumPixelsY() << ", "
             << this->getFocalLengthX() << ", "
             << this->getFocalLengthY() << ", "
             << this->getCenterU() << ", "
             << this->getCenterV() << ", "
             << m_radialCoefficient0 << ", "
             << m_radialCoefficient1 << ", "
             << m_radialCoefficient2 << ", "
             << m_radialCoefficient3 << ", "
             << m_radialCoefficient4 << ", "
             << m_radialCoefficient5 << ", "
             << m_tangentialCoefficient0 << ", "
             << m_tangentialCoefficient1 << "}";
      stream.precision(precision);
      stream.flags(flags);
      return stream;
    }


    template <class FloatType>
    void
    CameraIntrinsicsRational<FloatType>::
    projectWithPartialDerivatives(FloatType xNorm,
                                  FloatType yNorm,
                                  FloatType& uValue,
                                  FloatType& vValue,
                                  FloatType& dUdX,
                                  FloatType& dUdY,
                                  FloatType& dVdX,
                                  FloatType& dVdY) const
    {
      // First handle distortion.
      FloatType xDistorted;
      FloatType yDistorted;
      FloatType dXDdX;
      FloatType dXDdY;
      FloatType dYDdX;
      FloatType dYDdY;
      this->projectThroughDistortionWithPartialDerivatives(
        xNorm, yNorm, xDistorted, yDistorted, dXDdX, dXDdY, dYDdX, dYDdY);

      // Next handle pinhole projection.
      FloatType kX = this->getFocalLengthX();
      FloatType kY = this->getFocalLengthY();
      uValue = kX * xDistorted + this->getCenterU();
      vValue = kY * yDistorted + this->getCenterV();
      dUdX = kX * dXDdX;
      dUdY = kX * dXDdY;
      dVdX = kY * dYDdX;
      dVdY = kY * dYDdY;
    }
    
    
    // This member function takes a 2D point in the Z==1 plane of
    // camera coordinates, and returns an "distorted" version of
    // that 2D point.
    template <class FloatType>
    inline brick::numeric::Vector2D<FloatType>
    CameraIntrinsicsRational<FloatType>::
    projectThroughDistortion(const brick::numeric::Vector2D<FloatType>& point)
      const
    {
      FloatType xSquared = point.x() * point.x();
      FloatType ySquared = point.y() * point.y();
      FloatType rSquared = xSquared + ySquared;
      FloatType rFourth = rSquared * rSquared;
      FloatType rSixth = rSquared * rFourth;
    
      // Compute radial distortion terms.
      FloatType radialDistortionNumerator =
        (1.0 + m_radialCoefficient0 * rSquared
         + m_radialCoefficient1 * rFourth
         + m_radialCoefficient2 * rSixth);
      FloatType radialDistortionDenominator = 
        (1.0 + m_radialCoefficient3 * rSquared
         + m_radialCoefficient4 * rFourth
         + m_radialCoefficient5 * rSixth);
      FloatType radialDistortion =
        radialDistortionNumerator / radialDistortionDenominator;

      // Compute tangential distortion.
      FloatType crossTerm = point.x() * point.y();
      brick::numeric::Vector2D<FloatType> tangentialDistortion(
        (2.0 * m_tangentialCoefficient0 * crossTerm
         + m_tangentialCoefficient1 * (rSquared + 2.0 * xSquared)),
        (m_tangentialCoefficient0 * (rSquared + 2.0 * ySquared)
         + 2.0 * m_tangentialCoefficient1 * crossTerm));
    
      // Apply distortion and return.
      return brick::numeric::Vector2D<FloatType> (
        radialDistortion * point + tangentialDistortion);
    }

    
    template <class FloatType>
    inline brick::numeric::Vector3D<FloatType>
    CameraIntrinsicsRational<FloatType>::
    projectThroughDistortion(const brick::numeric::Vector3D<FloatType>& point)
      const
    {
      // The if clause here is to avoid Vector2D throwing an exception
      // if point.z() == 0.
      brick::numeric::Vector3D<FloatType> result;
      if(point.z() != static_cast<FloatType>(0.0)) {
        brick::numeric::Vector2D<FloatType> normalizedPoint;
        normalizedPoint.setValue(point.x(), point.y(), point.z());
        brick::numeric::Vector2D<FloatType> distortedPoint =
          this->projectThroughDistortion(normalizedPoint);
        result.setValue(distortedPoint.x(), distortedPoint.y(), 1.0);
      }
      return result;
    }


    // This member function takes a 2D point in the Z==1 plane of
    // camera coordinates, and returns an "distorted" version of
    // that 2D point.
    template <class FloatType>
    void
    CameraIntrinsicsRational<FloatType>::
    projectThroughDistortionWithPartialDerivatives(
      FloatType xNorm,
      FloatType yNorm,
      FloatType& xDistorted,
      FloatType& yDistorted,
      FloatType& dXDdX,
      FloatType& dXDdY,
      FloatType& dYDdX,
      FloatType& dYDdY) const
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

      // Compute distortion.
      FloatType radialDistortionNumerator =
        (1.0 + m_radialCoefficient0 * rSquared
         + m_radialCoefficient1 * rFourth
         + m_radialCoefficient2 * rSixth);
      FloatType radialDistortionDenominator =
        (1.0 + m_radialCoefficient3 * rSquared
         + m_radialCoefficient4 * rFourth
         + m_radialCoefficient5 * rSixth);
      FloatType radialDistortion =
        radialDistortionNumerator / radialDistortionDenominator;

      // Compute distortion derivatives.
      FloatType dRadialDistortionNumeratorDX =
        (m_radialCoefficient0 * dRSquaredDX
         + m_radialCoefficient1 * dRFourthDX
         + m_radialCoefficient2 * dRSixthDX);
      FloatType dRadialDistortionDenominatorDX =
        (m_radialCoefficient3 * dRSquaredDX
         + m_radialCoefficient4 * dRFourthDX
         + m_radialCoefficient5 * dRSixthDX);
      FloatType dRadialDistortionNumeratorDY =
        (m_radialCoefficient0 * dRSquaredDY
         + m_radialCoefficient1 * dRFourthDY
         + m_radialCoefficient2 * dRSixthDY);
      FloatType dRadialDistortionDenominatorDY =
        (m_radialCoefficient3 * dRSquaredDY
         + m_radialCoefficient4 * dRFourthDY
         + m_radialCoefficient5 * dRSixthDY);

      // Using the quotient rule for differentiation:
      // @code
      //   (f(x) / g(x))' = (g(x)f'(x) - f(x)g'(x)) / (g(x)g(x))
      // @endcode
      FloatType dRadialDistortionDX =
        (radialDistortionDenominator * dRadialDistortionNumeratorDX
         - radialDistortionNumerator * dRadialDistortionDenominatorDX)
        / (radialDistortionDenominator * radialDistortionDenominator);
      FloatType dRadialDistortionDY =
        (radialDistortionDenominator * dRadialDistortionNumeratorDY
         - radialDistortionNumerator * dRadialDistortionDenominatorDY)
        / (radialDistortionDenominator * radialDistortionDenominator);

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

      // Apply distortion.
      xDistorted = radialDistortion * xNorm + tangentialDistortion.x();
      yDistorted = radialDistortion * yNorm + tangentialDistortion.y();

      dXDdX = dRadialDistortionDX * xNorm + radialDistortion + dTangXDX;
      dXDdY = dRadialDistortionDY * xNorm + dTangXDY;
      dYDdX = dRadialDistortionDX * yNorm + dTangYDX;
      dYDdY = dRadialDistortionDY * yNorm + radialDistortion + dTangYDY;
    }


    template <class FloatType>
    bool
    reverseProjectWithJacobian(
      brick::numeric::Vector2D<FloatType>& rectifiedPoint,
      brick::numeric::Array2D<FloatType>& jacobian,
      brick::numeric::Vector2D<FloatType> imagePoint,
      CameraIntrinsicsRational<FloatType> intrinsics,
      FloatType requiredPrecision,
      std::size_t maximumIterations)
    {
      // Four pinhole parameters plus image U and V == 6 parameters
      // for initial reprojection.
      std::size_t constexpr numParameters0 = 6;

      // Eight distortion parameters plus X, Y = 10 parameters for
      // reprojection through distortion.
      std::size_t constexpr numParameters1 = 10;

      // 4 pinhole parameters + 8 distortion parameters + 2 image coordinates
      // == 14 total.
      std::size_t constexpr numParameters = 14;

      // We'll use automatic differentiation for some of this work.
      // These typedefs give us easy access to the necessary types.
      typedef brick::numeric::DifferentiableScalar<FloatType, numParameters0>
        DiffScalar0;
      typedef brick::numeric::DifferentiableScalar<FloatType, numParameters1>
        DiffScalar1;
      typedef CameraIntrinsicsRational<DiffScalar1> DiffIntrinsics;

      // Make a CameraIntrinsics instance that matches intrinsics, but
      // has machinery for doing automatic differentiation.
      typename CameraIntrinsicsRational<FloatType>::ParameterVectorType
        parameters = intrinsics.getDistortionCoefficients();
      typename DiffIntrinsics::ParameterVectorType diffParameters(
        parameters.size());
      for(std::size_t ii = 0; ii < parameters.size(); ++ii) {
        diffParameters[ii].setValue(parameters[ii]);
        diffParameters[ii].setPartialDerivative(ii, FloatType(1.0));
      }
      DiffIntrinsics diffIntrinsics;
      diffIntrinsics.setNumPixelsX(intrinsics.getNumPixelsX());
      diffIntrinsics.setNumPixelsY(intrinsics.getNumPixelsY());
      diffIntrinsics.setFreeParameters(diffParameters);

      // Start by reverse-projecting through the pinhole parameters,
      // tracking derivatives, but leaving the distortion model for
      // later.
      DiffScalar0 fX = intrinsics.getFocalLengthX();
      DiffScalar0 fY = intrinsics.getFocalLengthY();
      DiffScalar0 cU = intrinsics.getCenterU();
      DiffScalar0 cV = intrinsics.getCenterV();
      DiffScalar0 imageX = imagePoint.x();
      DiffScalar0 imageY = imagePoint.y();
      fX.setPartialDerivative(0, 1.0);
      fY.setPartialDerivative(1, 1.0);
      cU.setPartialDerivative(2, 1.0);
      cV.setPartialDerivative(3, 1.0);
      imageX.setPartialDerivative(4, 1.0);
      imageY.setPartialDerivative(5, 1.0);
      DiffScalar0 distortedX = (imageX - cU) / fX;
      DiffScalar0 distortedY = (imageY - cV) / fY;
      
      // Reverse projection is done iteratively, in a way that doesn't
      // lend itself to automatic differentiation.  Here we do the
      // reverse projection _without_ derivatives.
      brick::geometry::Ray3D<FloatType> reprojection =
        intrinsics.reverseProjectEM(imagePoint, true, requiredPrecision,
                                    maximumIterations);

      // Now _forward_ project, tracking derivatives.  This should
      // give us our original image point back, plus a bunch of
      // jacobian information.  Here we use the last two partial
      // derivatives to refer to partials with respect to rectified X
      // and Y.
      DiffScalar1 rectifiedX(reprojection.getDirectionVector().x()/
                             reprojection.getDirectionVector().z());
      DiffScalar1 rectifiedY(reprojection.getDirectionVector().y()/
                             reprojection.getDirectionVector().z());
      rectifiedX.setPartialDerivative(numParameters1 - 2, FloatType(1.0));
      rectifiedY.setPartialDerivative(numParameters1 - 1, FloatType(1.0));
      brick::numeric::Vector2D<DiffScalar1> rectifiedVector2D(
        rectifiedX, rectifiedY);
      brick::numeric::Vector2D<DiffScalar1> distortedPoint =
        diffIntrinsics.projectThroughDistortion(rectifiedVector2D);
        
      // Now that we know the partial derivatives of the forward
      // projection, Use the inverse function theorem to recover
      // partials of the rectified point (the point we just found) wrt
      // the distorted point (the point that came directly from our
      // input argument (via pinhole parameters)).  Let matrix A be
      // the jacobian of distorted point x, y wrt rectified point.
      // Let aij be the element of matrix A at row i, column j.
      FloatType a00 =
        distortedPoint.x().getPartialDerivative(numParameters1 - 2);
      FloatType a01 =
        distortedPoint.x().getPartialDerivative(numParameters1 - 1);
      FloatType a10 =
        distortedPoint.y().getPartialDerivative(numParameters1 - 2);
      FloatType a11 =
        distortedPoint.y().getPartialDerivative(numParameters1 - 1);

      // Here we do a quick cofactor inverse to get the jacobian of
      // the rectified point with respect to the distorted point.
      FloatType determinant = a00 * a11 - a01 * a10;
      if(determinant < std::numeric_limits<FloatType>::epsilon()) {
        return false;
      }
      FloatType b00 = a11 / determinant;
      FloatType b01 = -a01 / determinant;
      FloatType b10 = -a10 / determinant;
      FloatType b11 = a00 / determinant;
      
      // Now we use this jacobian to propagate all our other partial
      // derivatives (including those we got from the forward
      // projection) through the reverse projection.
      brick::numeric::Array2D<FloatType> dRectdDist(2, 2);
      dRectdDist(0, 0) = b00;
      dRectdDist(0, 1) = b01;
      dRectdDist(1, 0) = b10;
      dRectdDist(1, 1) = b11;

      // We also need the jacobian of the distorted point wrt all of
      // our parameters.  We allocate space for this here.
      brick::numeric::Array2D<FloatType> dDistdParams(2, numParameters);
      
      // The first four columns describe derivatives with respect to
      // pinhole parameters.  For these parameters, the jacobian we're
      // after is a straightforward composition of [the partial
      // derivative of the distorted point wrt parameters] and [the
      // partial derivative of rectified point wrt distorted point].
      // Accordingly, we copy partials of the distorted point directly
      // into dDistdParams, where it will be multiplied with
      // dRectdDist, below.
      std::size_t ii = 0;
      for(std::size_t jj = 0; jj < 4; ++jj) {
        dDistdParams(0, ii) = distortedX.getPartialDerivative(jj);
        dDistdParams(1, ii) = distortedY.getPartialDerivative(jj);
        ++ii;
      }

      // The next 8 columns describe derivatives wrt distortion
      // parameters.  There's a little cheat here... we negate the
      // partial derivatives we copy into dDistdParams.  This is
      // because the partials reported by distortedPoint describe the
      // same function as dRectdDist.  Here's one way to think about it:
      //   - Consider a differential change to distorted coordinates due
      //     to differential changes in distortion parameters.
      //   - We're looking for the change in rectified coordinates that
      //     exactly cancels this change.
      //   - Following this through, we're looking for the negative of the
      //     composition of the partials from distortedPoint and the
      //     partials in dRectdDist.  We make this happen by including a
      //     factor of -1 here.
      for(std::size_t jj = 0; jj < 8; ++jj) {
        dDistdParams(0, ii) = -distortedPoint.x().getPartialDerivative(jj);
        dDistdParams(1, ii) = -distortedPoint.y().getPartialDerivative(jj);
        ++ii;
      }

      // The last two columns describe derivatives with respect to the
      // input image point, and are reflected in our original reverse
      // projection through pinhole parameters.
      for(std::size_t jj = 0; jj < 2; ++jj) {
        dDistdParams(0, ii) = distortedX.getPartialDerivative(4 + jj);
        dDistdParams(1, ii) = distortedY.getPartialDerivative(4 + jj);
        ++ii;
      }

      jacobian = brick::numeric::matrixMultiply<FloatType>(
        dRectdDist, dDistdParams);
      rectifiedPoint.setValue(rectifiedX.getValue(), rectifiedY.getValue());

      return true;
    }
    
  } // namespace computerVision
  
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_CAMERAINTRINSICSRATIONAL_IMPL_HH */
