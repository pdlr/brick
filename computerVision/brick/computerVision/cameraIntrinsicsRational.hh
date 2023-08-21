/**
***************************************************************************
* @file brick/computerVision/cameraIntrinsicsRational.hh
*
* Header file declaring a CameraIntrinsics subclass for cameras
* conforming to the "full" camera model of OpenCV 2.4.
*
* Copyright (C) 2014 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_CAMERAINTRINSICSRATIONAL_HH
#define BRICK_COMPUTERVISION_CAMERAINTRINSICSRATIONAL_HH

#include <iostream>
#include <brick/computerVision/cameraIntrinsicsDistortedPinhole.hh>
#include <brick/numeric/array1D.hh>
#include <brick/numeric/utilities.hh>

namespace brick {

  namespace computerVision {

    /**
     ** This class represents calibration parameters for cameras
     ** conforming to the "full" camera model of OpenCV 2.4.
     **
     ** Our representation differs from the OpenCV convention in that
     ** OpenCV places the pixel coordinate (0, 0) in the center of the
     ** upper left pixel of the image, whereas we place (0, 0) at the
     ** upper left corner of the upper left pixel of the image.  In
     ** practice, you can often neglect this second difference because
     ** it is dominated by uncertainty in the principal point
     ** (centerU, centerV).  If you want to be absolutely correct,
     ** just modify constructor arguments as described in the
     ** constructor documentation.
     **
     ** The distortion model is as follows:
     **
     ** @code
     **   x_d = x * ((1 + k_0*r^2 + k_1*r^4 + k_2*r^6)
     **              / (1 + k_3*r^2 + k_4*r^4 + k_5*r^6))
     **         + 2 * p_0 * x * y
     **         + p_1 * (r^2 + 2 * x^2)
     **
     **   y_d = y * ((1 + k_0*r^2 + k_0*r^4 + k_2*r^6)
     **              / (1 + k_3*r^2 + k_4*r^4 + k_5*r^6))
     **         + p_0 * (r^2 + 2 * y^2)
     **         + 2 * p_1 * x * y
     ** @endcode
     **
     ** where x and y are undistorted coordinates, x_d and y_d are
     ** distorted coefficients, k_0 ... k_5 are radial distortion
     ** coefficients, and p0 and p1 are tangential coefficients.
     **
     ** In our implementation, we continue to use the two coordinate
     ** systems described in the CameraIntrinsicsPinhole class
     ** documentation: the 3D camera coordinate system, and the 2D
     ** pixel coordinate system.  Note that distortion parameters are
     ** applied in a third coordinate system: a physical 2D coordinate
     ** system coincident with the image plane.  The
     ** CameraIntrinsicsRational interface provides no access to this
     ** third coordinate system, so you can just ignore it.
     **/
    template <class FloatType>
    class CameraIntrinsicsRational
      : public CameraIntrinsicsDistortedPinhole<FloatType> {

    public:

      // Public member functions inherited from
      // CameraIntrinsicsDistortedPinhole.
      //
      // FloatType getCenterU();
      // FloatType getCenterV();
      // FloatType getFocalLengthX();
      // FloatType getFocalLengthY();
      // unsigned int getImageHeight();
      // unsigned int getImageWidth();
      // unsigned int getNumPixelsX();
      // unsigned int getNumPixelsY();
      // geometry::Ray3D<FloatType> reverseProject(...);
      // void setDependentParameters(...);
      // void setNumPixelsX(unsigned int);
      // void setNumPixelsY(unsigned int);


      /**
       * The default constructor initializes the
       * CameraIntrinsicsRational instance to a consistent (but not
       * terribly useful) state.
       */
      CameraIntrinsicsRational();


      /**
       * This constructor allows the caller to explicitly set the
       * camera intrinsic parameters.
       *
       * @param numPixelsX This argument specifies how many columns
       * there are in the camera images.
       *
       * @param numPixelsY This argument specifies how many rows there
       * are in the camera images.
       *
       * @param focalLengthX This argument the distance from the
       * camera focus to the image plane, expressed in pixel-width
       * sized units.  Generally this number should be positive,
       * indicating that the the image plane lies at a positive Z
       * coordinate in the 3D camera coordinate frame.
       *
       * @param focalLengthY This argument the distance from the
       * camera focus to the image plane, expressed in pixel-height
       * sized units.  Generally this number should be positive,
       * indicating that the the image plane lies at a positive Z
       * coordinate in the 3D camera coordinate frame.
       *
       * @param centerU This argument and the next specify the
       * position in pixel coordinates at which the Z axis passes
       * through the image plane.  If you are calling the constructor
       * using parameters computed by OpenCV, add 0.5 to OpenCV's
       * first principal point value to get the correct value for this
       * argument (see the documentation for CameraIntrinsicsRational
       * for more information).
       *
       * @param centerV This argument and the previous specify the
       * position in pixel coordinates at which the Z axis passes
       * through the image plane.  If you are calling the construtor
       * using parameters computed by OpenCV, add 0.5 to OpenCV's
       * second principal point value to get the correct value for this
       * argument (see the documentation for CameraIntrinsicsRational
       * for more information).
       *
       * @param radialCoefficient0 This argument specifies the
       * first quadratic term of the radial distortion model.  It
       * corresponds to k_0 in class documentation.
       *
       * @param radialCoefficient1 This argument specifies the first
       * 4th power term of the radial distortion model.  It
       * corresponds to k_1 in the class documentation.
       *
       * @param radialCoefficient2 This argument specifies the first
       * 6th power term of the radial distortion model.  It
       * corresponds to k_2 in the class documentation.
       *
       * @param radialCoefficient3 This argument specifies the second
       * quadratic term of the radial distortion model.  It
       * corresponds to k_3 in the class documentation.
       *
       * @param radialCoefficient4 This argument specifies the second
       * 6th power term of the radial distortion model.  It
       * corresponds to k_4 in the class documentation.
       *
       * @param radialCoefficient5 This argument specifies the second
       * 6th power term of the radial distortion model.  It
       * corresponds to k_5 in the class documentation.
       *
       * @param tangentialCoefficient0 This argument specifies the
       * first term of the tangential distortion model.  It
       * corresponds to p_0 in the class documentation.
       *
       * @param tangentialCoefficient1 This argument specifies the
       * second term of the tangential distortion model.  It
       * corresponds to p_1 in the class documentation.
       */
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
                               FloatType tangentialCoefficient1);


      /**
       * Destructor.
       */
      virtual
      ~CameraIntrinsicsRational() {}


      /**
       * This function exposes the distortion parameters of the camera model.
       *
       * @return The return value is a vector of free parameters
       * containing, in order, the six radial coefficients and the two
       * tangential coefficients.
       */
      typename CameraIntrinsicsRational<FloatType>::ParameterVectorType
      getDistortionCoefficients() const;


      /**
       * This function exposes a subset of the intrinsic parameters
       * for use in calibration routines.  Parameters that can
       * generally be calculated closed-form are omitted from this
       * return vector, leaving only those that are normally estimated
       * using nonlinear optimization.  Normally, this means leaving
       * out the pinhole projection parameters.
       *
       * @return The return value is a vector of free parameters
       * containing, in order, the six radial coefficients and the two
       * tangential coefficients.
       */
      virtual typename CameraIntrinsicsRational<FloatType>::ParameterVectorType
      getFreeParameters() const {return this->getDistortionCoefficients();}


      /**
       * This function provides a reasonable starting point for
       * intrinsic parameters that are generally estimated by
       * nonlinear optimization.  See getFreeParameters().
       *
       * @return The return value is a vector of free parameters
       * suitable for passing to setFreeParameters().
       */
      virtual typename CameraIntrinsicsRational<FloatType>::ParameterVectorType
      getNominalFreeParameters() const;


      /**
       * Returns a vector of all continuous parameters of the
       * class. Note that image dimensions (in pixels) are not
       * included in this vector.  The parameters (in order) are:
       *
       * @verbatim
       *   // Pinhole parameters
       *   - FocalLengthX
       *   - FocalLengthY
       *   - CenterU
       *   - CenterV
       *
       *   // Distortion coefficients.
       *   - RadialCoefficient0
       *   - RadialCoefficient1
       *   - RadialCoefficient2
       *   - RadialCoefficient3
       *   - RadialCoefficient4
       *   - RadialCoefficient5
       *   - TangentialCoefficient0
       *   - TangentialCoefficient1
       * @endverbatim
       *
       * @return The return value is the parameter vector described above.
       */
      virtual typename CameraIntrinsicsRational<FloatType>::ParameterVectorType
      getParameters() const;


      /**
       * This member function takes a point in 3D camera coordinates
       * and projects it into pixel coordinates.
       *
       * @param point This argument specifies the 3D point to be projected.
       *
       * @return The return value gives the point in pixel coordinates
       * to which the input point will project.
       */
      virtual brick::numeric::Vector2D<FloatType>
      project(const brick::numeric::Vector3D<FloatType>& point) const;


      /**
       * This member function takes a 2D point in the Z==1 plane of
       * camera coordinates, and returns an "distorted" version of
       * that 2D point.  The distorted point is not guaranteed to be
       * similar to the input point at all, but will project through
       * the idealized pinhole parameters associated with *this in
       * such a way that its projection is coincident with the
       * projection of the input point.  This member function is
       * generally not useful for user code.  It is provided here to
       * help with camera calibration algorithms.
       *
       * @param point This argument is the point to be projected,
       * represented in world coordinates.
       *
       * @return The return value is represented in a fictional
       * undistorted 3D world coordinate system, and is one of the
       * infinitely many points that lie on the ray projecting to the
       * 2D image point that corresponds to the input argument.
       */
      inline virtual numeric::Vector2D<FloatType>
      projectThroughDistortion(numeric::Vector2D<FloatType> const& point) const;


      /**
       * This member function takes a point in camera coordinates, and
       * returns an "distorted" version of that 3D point.  The
       * distorted point is not guaranteed to be similar to the input
       * point at all, but will project through the idealized pinhole
       * parameters associated with *this in such a way that its
       * projection is coincident with the projection of the input
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
      inline virtual numeric::Vector3D<FloatType>
      projectThroughDistortion(numeric::Vector3D<FloatType> const& point) const;


      /**
       * This member function takes a 2D point in the Z==1 plane of
       * camera coordinates, and returns an "distorted" version of
       * that 2D point, along with partial derivatives of the result
       * with respect to the input X and Y position.
       *
       * @param xNorm This argument is the X coordinate of the input
       * point.  Z is implicitly set to 1.0.
       *
       * @param yNorm This argument is the Y coordinate of the input
       * point.  Z is implicitly set to 1.0.
       *
       * @param xDistorted This argument returns the distorted X
       * coordinate by reference.  The distorted Z coordinate is
       * always equal to 1.0.
       *
       * @param yDistorted This argument returns the distorted Y
       * coordinate by reference.  The distorted Z coordinate is
       * always equal to 1.0.
       *
       * @param dXDdX This argument returns by reference the partial
       * derivative of the returned X coordinate with respect to the
       * input X coordinate.
       *
       * @param dXDdY This argument returns by reference the partial
       * derivative of the returned X coordinate with respect to the
       * input Y coordinate.
       *
       * @param dYDdX This argument returns by reference the partial
       * derivative of the returned Y coordinate with respect to the
       * input X coordinate.
       *
       * @param dYDdY This argument returns by reference the partial
       * derivative of the returned Y coordinate with respect to the
       * input Y coordinate.
       */
      void
      projectThroughDistortionWithPartialDerivatives(
        FloatType xNorm, FloatType yNorm,
        FloatType& xDistorted, FloatType& yDistorted,
        FloatType& dXDdX, FloatType& dXDdY,
        FloatType& dYDdX, FloatType& dYDdY) const;


      /**
       * This member function sets the calibration from an input
       * stream.  *this is modified only if the read was successful,
       * otherwise it is not modified, and failbit is set in the
       * stream state.
       *
       * @param inputStream This is the stream from which to read the
       * data.
       *
       * @return The return value is a reference to inputStream.
       */
      std::istream&
      readFromStream(std::istream& inputStream);


      /**
       * This function iteratively computes and returns a ray in 3D
       * camera coordinates starting at the camera focus and passing
       * through the specified pixel position.  It differs from member
       * function reverseProject() in that it uses a less general
       * iterative algorithm that requires the gradient of the
       * distortion function (from R2 --> R2) to have magnitude less
       * than 1.0.  This requirement is frequently met by real
       * cameras.  If you're used to using the OpenCV
       * undistortPoints() function, then this is what you're looking
       * for.
       *
       * @param pixelPosition This argument is the point in pixel
       * coordinates through which the returned ray should pass.
       *
       * @param normalize This argument indicates whether the ray
       * should be normalized to unit length before being returned.
       *
       * @param requiredPrecision This specifies the termination
       * criterion for the iterative algorithm.  After
       * minimumIterations has been reached, the algorithm will
       * terminate when an iteration improves the estimate by less
       * than this amount.
       *
       * @param maximumIterations Limits how long the algorithm can
       * run.  If requiredPrecision hasn't been satisfied within this
       * many iterations, an exception is thrown.
       *
       * @param minimumIterations Prevents the algorithm from
       * terminating prematurely.
       *
       * @return The return value is the resulting ray.
       */
      virtual geometry::Ray3D<FloatType>
      reverseProjectEM(
        const brick::numeric::Vector2D<FloatType>& pixelPosition,
        bool normalize = true,
        FloatType requiredPrecision = FloatType(1.0E-5),
        std::size_t maximumIterations = 25,
        std::size_t minimumIterations = 5) const;


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
       * order, the six radial coefficients and the two tangential
       * coefficients.
       */
      virtual void
      setFreeParameters(
        typename CameraIntrinsicsRational<FloatType>::ParameterVectorType
        const& parameterVector);


      /**
       * Sets the internal state of *this based on a parameter vector,
       * such as the one described in member function getParameters().
       * Note that image dimensions (in pixels) are not included in
       * this vector.
       *
       * @verbatim
       *   // Pinhole parameters
       *   - FocalLengthX
       *   - FocalLengthY
       *   - CenterU
       *   - CenterV
       *
       *   // Distortion coefficients.
       *   - RadialCoefficient0
       *   - RadialCoefficient1
       *   - RadialCoefficient2
       *   - RadialCoefficient3
       *   - RadialCoefficient4
       *   - RadialCoefficient5
       *   - TangentialCoefficient0
       *   - TangentialCoefficient1
       * @endverbatim
       */
      virtual void
      setParameters(
        typename CameraIntrinsicsRational<FloatType>::ParameterVectorType
        const& parameterVector);


      /**
       * This member function writes the calibration to an
       * outputstream in a format which is compatible with member
       * function readFromStream().
       *
       * @param outputStream This is the stream to which to write the
       * data.
       *
       * @return The return value is a reference to outputStream.
       */
      std::ostream&
      writeToStream(std::ostream& outputStream) const;


    protected:

      // Protected member function used during iterative approximation
      // in CameraIntrinsicsRational::reverseProject().
      void
      projectWithPartialDerivatives(FloatType xNorm,
                                    FloatType yNorm,
                                    FloatType& uValue,
                                    FloatType& vValue,
                                    FloatType& dUdX,
                                    FloatType& dUdY,
                                    FloatType& dVdX,
                                    FloatType& dVdY) const;


      // Protected data members.
      FloatType m_radialCoefficient0;
      FloatType m_radialCoefficient1;
      FloatType m_radialCoefficient2;
      FloatType m_radialCoefficient3;
      FloatType m_radialCoefficient4;
      FloatType m_radialCoefficient5;
      FloatType m_tangentialCoefficient0;
      FloatType m_tangentialCoefficient1;
    };


    /**
     * This function outputs a text representation of a
     * CameraIntrinsicsRational instance to a std::ostream.  The output
     * format looks like this:
     *
     * CameraIntrinsicsRational {240.0, 320.0, 2000.0, 3000.0, 640, 480, ...}
     *
     * @param stream This argument is a reference to the the output
     * stream.
     *
     * @param intrinsics This argument is a const reference to the
     * CameraIntrinsicsRational instance to be output.
     *
     * @return The return value is a reference to the input stream after
     * the write has taken place.
     */
    template <class FloatType>
    inline std::ostream&
    operator<<(std::ostream& stream,
               const CameraIntrinsicsRational<FloatType>& intrinsics)
    {
      return intrinsics.writeToStream(stream);
    }


    /**
     * This function sets the value of a CameraIntrinsicsRational
     * instance from a std::istream.  The input format is as described
     * for operator<<(std::ostream&, const CameraIntrinsicsRational&)
     * above.
     *
     * @param stream This argument is a reference to the the input
     * stream from which to read.
     *
     * @param intrinsics This argument is a reference to the
     * CameraIntrinsicsRational which will take the input.
     *
     * @return The return value is a reference to the input stream after
     * the read has taken place.
     */
    template <class FloatType>
    inline std::istream&
    operator>>(std::istream& stream,
               CameraIntrinsicsRational<FloatType>& intrinsics)
    {
      return intrinsics.readFromStream(stream);
    }


    template <class FloatType>
    bool
    reverseProjectWithJacobian(
      brick::numeric::Vector2D<FloatType>& rectifiedPoint,
      brick::numeric::Array2D<FloatType>& jacobian,
      brick::numeric::Vector2D<FloatType> imagePoint,
      CameraIntrinsicsRational<FloatType> intrinsics,
      FloatType requiredPrecision,
      std::size_t maximumIterations = 25);

  } // namespace computerVision

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/cameraIntrinsicsRational_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_CAMERAINTRINSICSRATIONAL_HH */
