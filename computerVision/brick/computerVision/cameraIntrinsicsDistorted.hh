/**
***************************************************************************
* @file brick/computerVision/cameraIntrinsicsDistorted.hh
*
* Header file declaring a parent class from which to derive classes which
* represent "distorted" camera intrinsic parameters.
*
* Copyright (C) 2007-2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_CAMERAINTRINSICSDISTORTED_HH
#define BRICK_COMPUTERVISION_CAMERAINTRINSICSDISTORTED_HH

#include <brick/computerVision/cameraIntrinsics.hh>
#include <brick/numeric/array1D.hh>

namespace brick {

  namespace computerVision {

    /**
     ** This abstract base class defines an interface for classes that
     ** describe camera projection parameters involving distortion
     ** models.  It requires the designer to explicitly distinguish
     ** between "free" parameters (those that must be found using
     ** nonlinear optimization) and "dependent" parameters (such as
     ** pinhole projection parameters, that can be determined
     ** closed-form once the free parameters have been set).  The
     ** advantage of this approach is that camera calibration
     ** routines, such as those declared in calibrationTools.h, can
     ** estimate distorted camera intrinsics by optimizing over only
     ** the free parameters.
     **/
    class CameraIntrinsicsDistorted : public CameraIntrinsics {
    public:

      typedef numeric::Array1D<double> ParameterVectorType;

      
      /** 
       * Default constructor.
       */
      CameraIntrinsicsDistorted() {};

      
      /** 
       * Destructor.
       */
      virtual
      ~CameraIntrinsicsDistorted() {}


      virtual ParameterVectorType
      getFreeParameters() const = 0;

      
      virtual ParameterVectorType
      getNominalFreeParameters() const = 0;


      virtual void
      setFreeParameters(ParameterVectorType const& parameterVector) = 0;

      
      virtual numeric::Vector3D<double>
      projectThroughDistortion(numeric::Vector3D<double> const& inputPoint)
        const = 0;

    protected:

    };


  } // namespace computerVision
  
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_CAMERAINTRINSICSDISTORTED_HH */
