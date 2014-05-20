/**
***************************************************************************
* @file brick/computerVision/cameraIntrinsicsPlumbBob.hh
*
* Header file declaring a CameraIntrinsics subclass for cameras
* conforming to Brown-Conrady "plumb bob" camera model.
*
* Copyright (C) 2007-2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_CAMERAINTRINSICSPLUMBBOB_HH
#define BRICK_COMPUTERVISION_CAMERAINTRINSICSPLUMBBOB_HH

#include <iostream>
#include <brick/computerVision/cameraIntrinsicsDistorted.hh>
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
      : public CameraIntrinsicsDistorted<FloatType> {
      
    public:

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
       * This member function provides access to the value of
       * parameter centerU, as described in the constructor
       * documentation.
       * 
       * @return The return value is the requested parameter.
       */
      FloatType
      getCenterU() const {return m_centerU;}


      /** 
       * This member function provides access to the value of
       * parameter centerV, as described in the constructor
       * documentation.
       * 
       * @return The return value is the requested parameter.
       */
      FloatType
      getCenterV() const {return m_centerV;}


      /** 
       * This member function provides access to the value of
       * parameter focalLengthX, as described in the constructor
       * documentation.
       * 
       * @return The return value is the requested parameter.
       */
      FloatType
      getFocalLengthX() const {return m_kX;}


      /** 
       * This member function provides access to the value of
       * parameter focalLengthY, as described in the constructor
       * documentation.
       * 
       * @return The return value is the requested parameter.
       */
      FloatType
      getFocalLengthY() const {return m_kY;}


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
      getFreeParameters() const;
      

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
       * This member function returns the number of columns in images
       * produced by the camera corresponding to *this.
       * 
       * @return The return value is the number of pixels in each
       * image row.
       */
      virtual unsigned int
      getNumPixelsX() const {return m_numPixelsX;}


      /** 
       * This member function returns the number of rows in images
       * produced by the camera corresponding to *this.
       * 
       * @return The return value is the number of pixels in each
       * image column.
       */
      virtual unsigned int
      getNumPixelsY() const {return m_numPixelsY;}

      
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
      inline virtual numeric::Vector3D<FloatType>
      projectThroughDistortion(numeric::Vector3D<FloatType> const& point) const;
      

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
       * This function is the counterpart to setFreeParameters().  It
       * allows the calling context to specify pinhole parameters that
       * are normally estimated closed-form.
       * 
       * @param numPixelsX This argument specifies how many columns
       * there are in the camera images.
       * 
       * @param numPixelsY This argument specifies how many rows there
       * are in the camera images.
       * 
       * @param focalLengthX This argument the distance from the
       * camera focus to the image plane in units of "pixel width."
       * Generally this number should be positive, indicating that the
       * the image plane lies at a positive Z coordinate in the 3D
       * camera coordinate frame.
       * 
       * @param focalLengthY This argument the distance from the
       * camera focus to the image plane in units of "pixel height."
       * Generally this number should be positive, indicating that the
       * the image plane lies at a positive Z coordinate in the 3D
       * camera coordinate frame.
       * 
       * @param centerU This argument and the next specify the
       * position in pixel coordinates at which the Z axis passes
       * through the image plane.
       * 
       * @param centerV This argument and the previous specify the
       * position in pixel coordinates at which the Z axis passes
       * through the image plane.
       */
      virtual void
      setDependentParameters(unsigned int numPixelsX,
                             unsigned int numPixelsY,
                             FloatType focalLengthX,
                             FloatType focalLengthY,
                             FloatType centerU,
                             FloatType centerV) {
        m_numPixelsX = numPixelsX;
        m_numPixelsY = numPixelsY;
        m_kX = focalLengthX;
        m_kY = focalLengthY;
        m_centerU = centerU;
        m_centerV = centerV;
      }

      
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
       * This member function specifies the width of the image in
       * pixels.
       * 
       * @param numPixelsX The width of the image.
       */
      virtual void
      setNumPixelsX(unsigned int numPixelsX) {m_numPixelsX = numPixelsX;}


      /** 
       * This member function specifies the height of the image in
       * pixels.
       * 
       * @param numPixelsX The height of the image.
       */
      virtual void
      setNumPixelsY(unsigned int numPixelsY) {m_numPixelsY = numPixelsY;}
      
      
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
      FloatType m_centerU;
      FloatType m_centerV;
      FloatType m_kX;
      FloatType m_kY;
      unsigned int m_numPixelsX;
      unsigned int m_numPixelsY;
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


  } // namespace computerVision
  
} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/cameraIntrinsicsPlumbBob_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_CAMERAINTRINSICSPLUMBBOB_HH */
