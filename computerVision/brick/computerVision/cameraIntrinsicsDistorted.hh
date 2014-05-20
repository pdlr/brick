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

    namespace privateCode {

      // Forward declaration of class to help with reverse projection.
      template <class FloatType>
      class ReverseProjectionObjective;

    } // namespace privateCode;

    
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
     **
     ** In addition, because many distortion models are not easily
     ** invertable, this class implements machinery for doing
     ** iterative reverse-projection and un-distortion of image
     ** points.  If a derived class has a better way to do these
     ** things, the relevant member functions can be overridden.
     ** There are, however, a few pure virtual member functions that
     ** support the iterative implementation, and these must
     ** overridden regardless of whether the derived class plans to
     ** make use of them or not.
     **/
    template <class FloatType = double>
    class CameraIntrinsicsDistorted : public CameraIntrinsics<FloatType> {
    public:

      typedef numeric::Array1D<FloatType> ParameterVectorType;

      // This class is defined externally so that it can be more
      // easily unit tested.
      friend class privateCode::ReverseProjectionObjective<FloatType>;
      
      /** 
       * Default constructor.
       */
      CameraIntrinsicsDistorted() {};

      
      /** 
       * Destructor.
       */
      virtual
      ~CameraIntrinsicsDistorted() {}


      /** 
       * This function exposes a subset of the intrinsic parameters
       * for use in calibration routines.  Parameters that can
       * generally be calculated closed-form are omitted from this
       * return vector, leaving only those that are normally estimated
       * using nonlinear optimization.
       * 
       * @return The return value is a vector of free parameters.  For
       * example, if the subclass implements the Brown-Conrady model,
       * it might be a vector containing radial coefficients and
       * tangential coefficients.
       */
      virtual ParameterVectorType
      getFreeParameters() const = 0;


      /** 
       * This function provides a reasonable starting point for
       * intrinsic parameters that are generally estimated by
       * nonlinear optimization.  See getFreeParameters().
       *
       * @return The return value is a vector of free parameters
       * suitable for passing to setFreeParameters().
       */
      virtual ParameterVectorType
      getNominalFreeParameters() const = 0;


      /** 
       * This member function returns the number of columns in images
       * produced by the camera corresponding to *this.
       * 
       * @return The return value is the number of pixels in each
       * image row.
       */
      virtual unsigned int
      getNumPixelsX() const = 0;


      /** 
       * This member function returns the number of rows in images
       * produced by the camera corresponding to *this.
       * 
       * @return The return value is the number of pixels in each
       * image column.
       */
      virtual unsigned int
      getNumPixelsY() const = 0;

      
      /** 
       * This member function takes a point in camera coordinates, and
       * returns an "undistorted" version of that 3D point.  The
       * undistorted point is not guaranteed to be similar to the
       * input point at all, but will project through the idealized
       * pinhole parameters associated with *this in such a way that
       * its projection is coincident with the projection of the input
       * point.  This member function is generally not useful for user
       * code.  It is provided here to help with camera calibration
       * algorithms.
       * 
       * @param point This argument is the point to be projected,
       * represented in world coordinates.
       * 
       * @return The return value is represented in a fictional
       * undistorted 3D world coordinate system, and is one of the
       * infinitely many points that lie on the ray projecting to the
       * 2D image point that corresponds to the input argument.
       */
      virtual numeric::Vector3D<FloatType>
      projectThroughDistortion(numeric::Vector3D<FloatType> const& inputPoint)
        const = 0;

      
      /** 
       * This function iteratively computes and returns a ray in 3D
       * camera coordinates starting at the camera focus and passing
       * through the specified pixel position.  Derived classes may
       * choose to implement closed-form (or simply more efficient)
       * versions of this function.
       * 
       * @param pixelPosition This argument is the point in pixel
       * coordinates through which the returned ray should pass.
       * 
       * @param normalize This argument indicates whether the ray
       * should be normalized to unit length before being returned.
       * 
       * @return The return value is the resulting ray.
       */
      virtual geometry::Ray3D<FloatType>
      reverseProject(const brick::numeric::Vector2D<FloatType>& pixelPosition,
                     bool normalize = true) const;
      

      /** 
       * This sets the value of a subset of the intrinsic parameters,
       * and is commonly used by in calibration routines.  Parameters
       * that can generally be calculated closed-form are omitted from
       * this return vector, leaving only those that are normally
       * estimated using nonlinear optimization.  The omitted
       * parameters are generally pinhole projection parameters.
       * 
       * @param parameterVector This argument specifies values for the
       * free parameters as a vector of free parameters containing, in
       * order, the three radial coefficients, the skew coefficient,
       * and the two tangential coefficients.  If the third radial
       * coefficent and/or the skew coefficient have been disallowed
       * (by this->allow*()), then they must be omitted from this
       * vector.
       */
      virtual void
      setFreeParameters(ParameterVectorType const& parameterVector) = 0;

    protected:

      // By default, we require reverse projection accuracy of better
      // than 1/10th of a pixel.
      virtual FloatType
      getMaximumReverseProjectionResidual() const {
        return static_cast<FloatType>(0.01); // 0.1 * 0.1 == 0.01.
      }


      // This pure virtual function must be overridden to provide both
      // the projection of the 3D point (xNorm, yNorm, 1.0) and the
      // derivatives of that projection with respect to xNorm and
      // yNorm.  See CameraIntrinsicsPlumbBob for an example.
      virtual void
      projectWithPartialDerivatives(FloatType xNorm,
                                    FloatType yNorm,
                                    FloatType& uValue,
                                    FloatType& vValue,
                                    FloatType& dUdX,
                                    FloatType& dUdY,
                                    FloatType& dVdX,
                                    FloatType& dVdY) const = 0;

    };


  } // namespace computerVision
  
} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/cameraIntrinsicsDistorted_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_CAMERAINTRINSICSDISTORTED_HH */
