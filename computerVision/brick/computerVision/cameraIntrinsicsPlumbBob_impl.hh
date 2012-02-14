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
    
namespace brick {

  namespace computerVision {

    namespace privateCode {

      /**
       ** This functor is called by
       ** CameraIntrinsicsPlumbBob::reverseProject() during iterative
       ** approximation of reverse projection.
       **/
      class PlumbBobObjective
        : public std::unary_function<brick::numeric::Array1D<double>, double>
      {
      public:
        PlumbBobObjective(const CameraIntrinsicsPlumbBob& intrinsics,
                          const brick::numeric::Vector2D<double>& uvTarget);
        
        double
        operator()(const brick::numeric::Array1D<double>& theta);

        double
        getOffset();
        
        brick::numeric::Array1D<double>
        gradient(const brick::numeric::Array1D<double>& theta);

      private:

        double
        computeBoundsPenalty(
          const brick::numeric::Vector2D<double>& uvPosition);

        void
        computeBoundsPenaltyGradient(double uValue, double vValue,
                                     double dUdX, double dUdY,
                                     double dVdX, double dVdY,
                                     double& dPdX, double& dPdY);

        CameraIntrinsicsPlumbBob m_intrinsics;
        double m_offset;
        brick::numeric::Vector2D<double> m_uvTarget;
      };
      
    } // namespace privateCode;
    
    
    inline brick::numeric::Vector3D<double>
    CameraIntrinsicsPlumbBob::
    projectThroughDistortion(const brick::numeric::Vector3D<double>& point)
      const
    {
      // The Vector2D constructor will divide by point.z(), which is will
      // throw if point.z() == 0.0.
      brick::numeric::Vector2D<double> normalizedPoint(point.x(), point.y(), point.z());
    
      double xSquared = normalizedPoint.x() * normalizedPoint.x();
      double ySquared = normalizedPoint.y() * normalizedPoint.y();
      double rSquared = xSquared + ySquared;
      double rFourth = rSquared * rSquared;
      double rSixth = rSquared * rFourth;
    
      // Compute distortion terms according to plumb bob model.
      double radialDistortion = (1.0 + m_radialCoefficient0 * rSquared
                                 + m_radialCoefficient1 * rFourth
                                 + m_radialCoefficient2 * rSixth);
    
      double crossTerm = normalizedPoint.x() * normalizedPoint.y();
      brick::numeric::Vector2D<double> tangentialDistortion(
        (2.0 * m_tangentialCoefficient0 * crossTerm
         + m_tangentialCoefficient1 * (rSquared + 2.0 * xSquared)),
        (m_tangentialCoefficient0 * (rSquared + 2.0 * ySquared)
         + 2.0 * m_tangentialCoefficient1 * crossTerm));
    
      // Apply distortion and skew, then project into pixel coordinates.
      brick::numeric::Vector2D<double> distortedPoint(
        radialDistortion * normalizedPoint + tangentialDistortion);
      return brick::numeric::Vector3D<double>(
        distortedPoint.x() + m_skewCoefficient * distortedPoint.y(),
        distortedPoint.y(),
        1.0);
    }

  } // namespace computerVision
  
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_CAMERAINTRINSICSPLUMBBOB_IMPL_HH */
