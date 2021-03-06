/**
***************************************************************************
* @file brick/computerVision/cameraIntrinsicsPlumbBob.hh
*
* Header file declaring a CameraIntrinsics subclass for cameras
* conforming to Brown-Conrady "plumb bob" camera model.
*
* Copyright (C) 2007-2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_CAMERAINTRINSICSPLUMBBOB_HH
#define BRICK_COMPUTERVISION_CAMERAINTRINSICSPLUMBBOB_HH

#include <iostream>
#include <brick/computerVision/cameraIntrinsicsDistortedPinhole.hh>
#include <brick/numeric/array1D.hh>
#include <brick/numeric/utilities.hh>

namespace brick {

  namespace computerVision {


    /**
     ** This class represents calibration parameters for cameras
     ** conforming to Brown-Conrady "plumb bob" camera model, as
     ** described in [1].  This is the intrinsic model used by the
     ** Matlab camera calibration toolbox.
     **
     ** Our representation differs from the matlab convention in two
     ** respects.  The first difference is that the order of
     ** distortion coefficients in the constructor arguments differs
     ** slightly from the order used by the toolbox.  This first
     ** difference is described in the constructor documentation.  The
     ** second difference is that matlab places the pixel coordinate
     ** (0, 0) in the center of the upper left pixel of the image,
     ** whereas we place (0, 0) at the upper left corner of the upper
     ** left pixel of the image.  In practice, you can often neglect
     ** this second difference because it is dominated by uncertainty
     ** in the principal point (centerU, centerV).  If you want to be
     ** absolutely correct, just modify constructor arguments as
     ** described in the constructor documentation.
     **
     ** In our implementation, we continue to use the two coordinate
     ** systems described in the CameraIntrinsicsPinhole class
     ** documentation: the 3D camera coordinate system, and the 2D
     ** pixel coordinate system.  Note that distortion parameters are
     ** applied in a third coordinate system: a physical 2D coordinate
     ** system coincident with the image plane.  The
     ** CameraIntrinsicsPlumbBob interface provides no access to this
     ** third coordinate system, so you can just ignore it.
     **
     ** [1] D.C. Brown, Decentering Distortion of Lenses.  Photometric
     ** Engineering, pp. 444-462, Vol. 32, No. 3, 1966.
     **/
    template <class FloatType>
    class CameraIntrinsicsPlumbBob
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
       * CameraIntrinsicsPlumbBob instance to a consistent (but not
       * terribly useful) state.
       */
      CameraIntrinsicsPlumbBob();


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
       * using parameters computed by matlab, add 0.5 to matlab's
       * first principal point value to get the correct value for this
       * argument (see the documentation for CameraIntrinsicsPlumbBob
       * for more information).
       *
       * @param centerV This argument and the previous specify the
       * position in pixel coordinates at which the Z axis passes
       * through the image plane.  If you are calling the construtor
       * using parameters computed by matlab, add 0.5 to matlab's
       * second principal point value to get the correct value for this
       * argument (see the documentation for CameraIntrinsicsPlumbBob
       * for more information).
       *
       * @param skewCoefficient This argument specifies any deviation
       * from squareness in the pixel grid.  For most modern cameras,
       * this parameter will be zero.  Specifically, skewCoefficient is
       * the cosine of the angle between the x and y axes, where the
       * nominal angle is positive pi/2 radians.
       *
       * @param radialCoefficient0 This argument specifies the
       * quadratic term of the radial distortion model.  It
       * corresponds to kc(1) in the matlab representation.
       *
       * @param radialCoefficient1 This argument specifies the 4th
       * power term of the radial distortion model.  It corresponds to
       * kc(2) in the matlab representation.
       *
       * @param radialCoefficient2 This argument specifies the 6th
       * power term of the radial distortion model.  It corresponds to
       * kc(5) in the matlab representation.
       *
       * @param tangentialCoefficient0 This argument specifies the
       * first term of the tangential distortion model.  It
       * corresponds to kc(3) in the matlab representation.
       *
       * @param tangentialCoefficient1 This argument specifies the
       * second term of the tangential distortion model.  It
       * corresponds to kc(4) in the matlab representation.
       */
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
                               FloatType tangentialCoefficient1);


      /**
       * Destructor.
       */
      virtual
      ~CameraIntrinsicsPlumbBob() {}


      /**
       * Call this function before running calibration routines to
       * specify whether or not you want to allow a nonzero third
       * radial distortion coefficient.
       *
       * @param flag This argument indicates whether the coefficient
       * should be estimate (true) or zeroed (false).
       */
      void
      allowSixthOrderRadial(bool flag = true);


      /**
       * Call this function before running calibration routines to
       * specify whether or not you want to allow a nonzero skew
       * coefficient.
       *
       * @param flag This argument indicates whether the coefficient
       * should be estimate (true) or zeroed (false).
       */
      void
      allowSkew(bool flag = true);


      /**
       * Indicates whether or not a nonzero third radial distortion
       * coefficient is allowed.
       *
       * @return The return value indicates whether a nonzero
       * coefficient is allowed (true) or not (false).
       */
      bool
      getAllowSixthOrderRadial() const {return m_allowSixthOrderRadial;}


      /**
       * Indicates whether or not a nonzero skew coefficient is
       * allowed.
       *
       * @return The return value indicates whether a nonzero
       * coefficient is allowed (true) or not (false).
       */
      bool
      getAllowSkew() const {return m_allowSkew;}


      /**
       * This function exposes the distortion parameters of the camera model.
       *
       * @return The return value is a vector of free parameters
       * containing, in order, the three radial coefficients, the skew
       * coefficient, and the two tangential coefficients.  If the
       * third radial coefficent and/or the skew coefficient have been
       * disallowed (by this->allow*()), then they are omitted from
       * the return vector.
       */
      typename CameraIntrinsicsPlumbBob<FloatType>::ParameterVectorType
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
       * containing, in order, the three radial coefficients, the skew
       * coefficient, and the two tangential coefficients.  If the
       * third radial coefficent and/or the skew coefficient have been
       * disallowed (by this->allow*()), then they are omitted from
       * the return vector.
       */
      virtual typename CameraIntrinsicsPlumbBob<FloatType>::ParameterVectorType
      getFreeParameters() const {return this->getDistortionCoefficients();}


      /**
       * This function provides a reasonable starting point for
       * intrinsic parameters that are generally estimated by
       * nonlinear optimization.  See getFreeParameters().
       *
       * @return The return value is a vector of free parameters
       * suitable for passing to setFreeParameters().
       */
      virtual typename CameraIntrinsicsPlumbBob<FloatType>::ParameterVectorType
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
       *   - RadialCoefficient2  // Only if allowed by
       *                         // this->allowSixthOrderRadial()
       *   - SkewCoefficient     // Only if allowed by this->allowSkew()
       *   - TangentialCoefficient0
       *   - TangentialCoefficient1
       * @endverbatim
       *
       * @return The return value is the parameter vector described above.
       */
      virtual typename CameraIntrinsicsPlumbBob<FloatType>::ParameterVectorType
      getParameters() const;


      /**
       * This member function provides access to the value of
       * parameter radialCoefficient0, as described in the constructor
       * documentation.
       *
       * @return The return value is the requested parameter.
       */
      FloatType
      getRadialCoefficient0() const {return m_radialCoefficient0;}


      /**
       * This member function provides access to the value of
       * parameter radialCoefficient1, as described in the constructor
       * documentation.
       *
       * @return The return value is the requested parameter.
       */
      FloatType
      getRadialCoefficient1() const {return m_radialCoefficient1;}


      /**
       * This member function provides access to the value of
       * parameter radialCoefficient2, as described in the constructor
       * documentation.
       *
       * @return The return value is the requested parameter.
       */
      FloatType
      getRadialCoefficient2() const {return m_radialCoefficient2;}


      /**
       * This member function provides access to the value of
       * parameter skewCoefficient, as described in the constructor
       * documentation.
       *
       * @return The return value is the requested parameter.
       */
      FloatType
      getSkewCoefficient() const {return m_skewCoefficient;}


      /**
       * This member function provides access to the value of
       * parameter tangentialCoefficient0, as described in the
       * constructor documentation.
       *
       * @return The return value is the requested parameter.
       */
      FloatType
      getTangentialCoefficient0() const {return m_tangentialCoefficient0;}


      /**
       * This member function provides access to the value of
       * parameter tangentialCoefficient1, as described in the
       * constructor documentation.
       *
       * @return The return value is the requested parameter.
       */
      FloatType
      getTangentialCoefficient1() const {return m_tangentialCoefficient1;}


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
       * order, the three radial coefficients, the skew coefficient,
       * and the two tangential coefficients.  If the third radial
       * coefficent and/or the skew coefficient have been disallowed
       * (by this->allow*()), then they must be omitted from this
       * vector.
       */
      virtual void
      setFreeParameters(
        typename CameraIntrinsicsPlumbBob<FloatType>::ParameterVectorType
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
       *   - RadialCoefficient2  // Only if allowed by
       *                         // this->allowSixthOrderRadial()
       *   - SkewCoefficient     // Only if allowed by this->allowSkew()
       *   - TangentialCoefficient0
       *   - TangentialCoefficient1
       * @endverbatim
       */
      virtual void
      setParameters(
        typename CameraIntrinsicsPlumbBob<FloatType>::ParameterVectorType
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
      // in CameraIntrinsicsPlumbBob::reverseProject().
      void
      projectWithPartialDerivatives(FloatType xNorm,
                                    FloatType yNorm,
                                    FloatType& uValue,
                                    FloatType& vValue,
                                    FloatType& dUdX,
                                    FloatType& dUdY,
                                    FloatType& dVdX,
                                    FloatType& dVdY) const;

      bool m_allowSixthOrderRadial;
      bool m_allowSkew;
      FloatType m_radialCoefficient0;
      FloatType m_radialCoefficient1;
      FloatType m_radialCoefficient2;
      FloatType m_skewCoefficient;
      FloatType m_tangentialCoefficient0;
      FloatType m_tangentialCoefficient1;
    };


    /**
     * This function outputs a text representation of a
     * CameraIntrinsicsPlumbBob instance to a std::ostream.  The output
     * format looks like this:
     *
     * CameraIntrinsicsPlumbBob {240.0, 320.0, 2000.0, 3000.0, 640, 480, ...}
     *
     * @param stream This argument is a reference to the the output
     * stream.
     *
     * @param intrinsics This argument is a const reference to the
     * CameraIntrinsicsPlumbBob instance to be output.
     *
     * @return The return value is a reference to the input stream after
     * the write has taken place.
     */
    template <class FloatType>
    inline std::ostream&
    operator<<(std::ostream& stream,
               const CameraIntrinsicsPlumbBob<FloatType>& intrinsics)
    {
      return intrinsics.writeToStream(stream);
    }


    /**
     * This function sets the value of a CameraIntrinsicsPlumbBob
     * instance from a std::istream.  The input format is as described
     * for operator<<(std::ostream&, const CameraIntrinsicsPlumbBob&)
     * above.
     *
     * @param stream This argument is a reference to the the input
     * stream from which to read.
     *
     * @param intrinsics This argument is a reference to the
     * CameraIntrinsicsPlumbBob which will take the input.
     *
     * @return The return value is a reference to the input stream after
     * the read has taken place.
     */
    template <class FloatType>
    inline std::istream&
    operator>>(std::istream& stream,
               CameraIntrinsicsPlumbBob<FloatType>& intrinsics)
    {
      return intrinsics.readFromStream(stream);
    }


    template <class FloatType>
    bool
    reverseProjectWithJacobian(
      brick::numeric::Vector2D<FloatType>& rectifiedPoint,
      brick::numeric::Array2D<FloatType>& jacobian,
      brick::numeric::Vector2D<FloatType> imagePoint,
      CameraIntrinsicsPlumbBob<FloatType> intrinsics,
      FloatType requiredPrecision,
      std::size_t maximumIterations = 25);

  } // namespace computerVision

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/cameraIntrinsicsPlumbBob_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_CAMERAINTRINSICSPLUMBBOB_HH */
