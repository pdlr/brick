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
#include <brick/optimization/optimizerBFGS.hh>
#include <brick/optimization/optimizerNelderMead.hh>

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
                    "CameraIntrinsicsRational::getFreeParameters()",
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

      // Apply distortion, then project into pixel coordinates.
      brick::numeric::Vector2D<FloatType> distortedPoint(
        radialDistortion * xNorm + tangentialDistortion.x(),
        radialDistortion * yNorm + tangentialDistortion.y());

      FloatType dDistortedXDX =
        dRadialDistortionDX * xNorm + radialDistortion + dTangXDX;
      FloatType dDistortedXDY = dRadialDistortionDY * xNorm + dTangXDY;
      FloatType dDistortedYDX = dRadialDistortionDX * yNorm + dTangYDX;
      FloatType dDistortedYDY =
        dRadialDistortionDY * yNorm + radialDistortion + dTangYDY;

      FloatType kX = this->getFocalLengthX();
      FloatType kY = this->getFocalLengthY();

      uValue = kX * distortedPoint.x() + this->getCenterU();
      vValue = kY * distortedPoint.y() + this->getCenterV();
      dUdX = kX * dDistortedXDX;
      dUdY = kX * dDistortedXDY;
      dVdX = kY * dDistortedYDX;
      dVdY = kY * dDistortedYDY;
    }
    
    
    
    template <class FloatType>
    inline brick::numeric::Vector3D<FloatType>
    CameraIntrinsicsRational<FloatType>::
    projectThroughDistortion(const brick::numeric::Vector3D<FloatType>& point)
      const
    {
      // The if clause here is to avoid Vector2D throwing an exception
      // if point.z() == 0.
      brick::numeric::Vector2D<FloatType> normalizedPoint;
      if(point.z() != 0.0) {
        normalizedPoint.setValue(point.x(), point.y(), point.z());
      }
    
      FloatType xSquared = normalizedPoint.x() * normalizedPoint.x();
      FloatType ySquared = normalizedPoint.y() * normalizedPoint.y();
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
      FloatType crossTerm = normalizedPoint.x() * normalizedPoint.y();
      brick::numeric::Vector2D<FloatType> tangentialDistortion(
        (2.0 * m_tangentialCoefficient0 * crossTerm
         + m_tangentialCoefficient1 * (rSquared + 2.0 * xSquared)),
        (m_tangentialCoefficient0 * (rSquared + 2.0 * ySquared)
         + 2.0 * m_tangentialCoefficient1 * crossTerm));
    
      // Apply distortion and return.
      brick::numeric::Vector2D<FloatType> distortedPoint(
        radialDistortion * normalizedPoint + tangentialDistortion);
      return brick::numeric::Vector3D<FloatType>(
        distortedPoint.x(), distortedPoint.y(), 1.0);
    }
    
  } // namespace computerVision
  
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_CAMERAINTRINSICSRATIONAL_IMPL_HH */
