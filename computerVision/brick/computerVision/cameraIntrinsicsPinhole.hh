/**
***************************************************************************
* @file brick/computerVision/cameraIntrinsicsPinhole.hh
*
* Header file declaring a CameraIntrinsics subclass for pinhole
* cameras.
*
* Copyright (C) 2007 - 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_CAMERAINTRINSICSPINHOLE_HH
#define BRICK_COMPUTERVISION_CAMERAINTRINSICSPINHOLE_HH

#include <iostream>
#include <brick/computerVision/cameraIntrinsics.hh>
#include <brick/geometry/ray3D.hh>
#include <brick/numeric/array2D.hh>
#include <brick/numeric/vector2D.hh>


namespace brick {

  namespace computerVision {

    /**
     ** This class represents calibration parameters for a simple
     ** pinhole camera model, as described in [1].  There are two
     ** coordinate systems associated with the pinhole camera: the 3D
     ** camera coordinate system, and the 2D pixel coordinate system.
     **
     ** The camera coordinate system is right-handed, and has its
     ** origin at the focus of the camera.  The Z axis is
     ** perpendicular to the image plane, points in the direction of
     ** the camera field of view, and intersects the image plane at
     ** (0, 0, focalLength).  The X axis runs parallel to the rows of
     ** the image, increasing from the left of the image to the right.
     ** The Y axis runs parallel to the columns of the image,
     ** increasing from the top of the image to the bottom.  We
     ** generally consider X, Y, Z to have units of meters, although
     ** nothing in the implementation of this class requires this.
     **
     ** The pixel coordinate system has its origin at one corner of
     ** the image.  We call the two axes U and V.  The U axis points
     ** along the first row of the image.  The U axis is parallel to,
     ** and increases along, the X axis of the camera coordinate
     ** system.  The V axis points along the first column, parallel to
     ** (and increasing along) the Y axis of the camera coordinate
     ** system.  U and V have units of pixels.  The pixel at indices
     ** (r, c) in the image (that is, the pixel at row r and column c)
     ** has its four corners at (u, v) coordinates (c, r), (c + 1, r),
     ** (c, r + 1), and (c + 1, r + 1).  For example, the upper left
     ** pixel in the image (row 0, column 0) has its upper left corner
     ** at (u, v) = (0, 0) and its center at (u, v) = (0.5, 0.5).
     **
     ** When displaying an image, our convention is to orient the
     ** display so that the origin of the pixel coordinate system is
     ** at the top left of the display, with the U axis pointing to
     ** the right and the V axis pointing down.
     **/
    template <class FloatType = double>
    class CameraIntrinsicsPinhole : public CameraIntrinsics<FloatType> {
    public:

      /** 
       * The default constructor initializes the
       * CameraIntrinsicsPinhole instance to a consistent (but not
       * terribly useful) state.
       */
      CameraIntrinsicsPinhole();

      
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
       * @param focalLength This argument the distance from the camera
       * focus to the image plane.  Generally this number should be
       * positive, indicating that the the image plane lies at a
       * positive Z coordinate in the 3D camera coordinate frame.  Our
       * convention is to specify it in units of meters, but you can
       * use whatever unit you like, provided you use the same unit
       * for arguments pixelSizeX and pixelSizeY.
       * 
       * @param pixelSizeX This argument specifies the width of an
       * individual pixel.  In other words, pixelSizeX is the distance
       * between the centers of adjacent pixels in the same row.  Our
       * convention is to specify it in units of meters, but you can
       * use whatever unit you like, provided you use the same unit
       * for arguments focalLength and pixelSizeY.
       * 
       * @param pixelSizeY This argument specifies the height of an
       * individual pixel.  In other words, pixelSizeX is the distance
       * between the centers of adjacent pixels in the same column.
       * Our convention is to specify it in units of meters, but you
       * can use whatever unit you like, provided you use the same
       * unit for arguments focalLength and pixelSizeX.
       * 
       * @param centerU This argument and the next specify the
       * position in pixel coordinates at which the Z axis passes
       * through the image plane.
       * 
       * @param centerV This argument and the previous specify the
       * position in pixel coordinates at which the Z axis passes
       * through the image plane.
       */
      CameraIntrinsicsPinhole(unsigned int numPixelsX,
                              unsigned int numPixelsY,
                              FloatType focalLength,
                              FloatType pixelSizeX,
                              FloatType pixelSizeY,
                              FloatType centerU,
                              FloatType centerV);


      /** 
       * Destructor.
       */
      virtual
      ~CameraIntrinsicsPinhole() {}


      /** 
       * This member function returns the U coordinate of the
       * projection center (principle point) of the image.
       * 
       * @return The return value is the U coordinate (pixel
       * coordinate) at which the optical axis intersects the image
       * plane.
       */
      FloatType
      getCenterU() const {return m_centerU;}
          

      /** 
       * This member function returns the V coordinate of the
       * projection center (principle point) of the image.
       * 
       * @return The return value is the V coordinate (pixel
       * coordinate) at which the optical axis intersects the image
       * plane.
       */
      FloatType
      getCenterV() const {return m_centerV;}
          

      /** 
       * This member function returns the focal length of the camera,
       * as passed to the CameraIntrinsicsPinhole constructor.
       * 
       * @return The return value is the originally specified focal length.
       */
      FloatType
      getFocalLength() const {return m_focalLength;}
          

      /** 
       * This member function returns the camera focal length
       * expressed in units of pixel width.
       * 
       * @return The return value is focalLength / xPixelSize.
       */
      FloatType
      getKx() const {return m_kX;}
          

      /** 
       * This member function returns the camera focal length
       * expressed in units of pixel height.
       * 
       * @return The return value is focalLength / yPixelSize.
       */
      FloatType
      getKy() const {return m_kY;}
          

      /** 
       * This member function returns the width of the image in
       * pixels.
       * 
       * @return The return value is the number of pixels in the X
       * direction.
       */
      unsigned int
      getNumPixelsX() const {return m_numPixelsX;}
          

      /** 
       * This member function returns the height of the image in
       * pixels.
       * 
       * @return The return value is the number of pixels in the Y
       * direction.
       */
      unsigned int
      getNumPixelsY() const {return m_numPixelsY;}
          

      /** 
       * This member function returns the width of the image pixels in
       * whatever units were used in the constructor call.
       * 
       * @return The return value is the size of each pixel in the X
       * direction.
       */
      FloatType
      getPixelSizeX() const {return m_focalLength / m_kX;}
          

      /** 
       * This member function returns the height of the image pixels in
       * whatever units were used in the constructor call.
       * 
       * @return The return value is the size of each pixel in the Y
       * direction.
       */
      FloatType
      getPixelSizeY() const {return m_focalLength / m_kY;}
          

      /** 
       * This member function returns a coordinate transform that
       * "matches" *this.  That is, projecting a world point using the
       * returned coordinate transform has the same effect as passing
       * that world point as an argument to this->project().
       * 
       * @return The return value is the relevant projection matrix.
       */
      brick::numeric::Array2D<FloatType>
      getProjectionMatrix() const;

      
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
       * This member function takes a point in 2D pixel coordinates
       * and returns a ray in 3D camera coordinates passing through
       * all of the 3D points that project to the specified 2D
       * position.
       * 
       * @param pixelPosition This argument is the point to be
       * projected out into 3D camera coordinates.
       * 
       * @param normalize This argument indicates whether the returned
       * vector should be normalized to unit length.  Setting this to
       * false saves a few arithmetic operations.
       * 
       * @return The return value is the ray in 3D camera coordinates
       * corresponding to the input 2D point.
       */
      virtual brick::geometry::Ray3D<FloatType>
      reverseProject(const brick::numeric::Vector2D<FloatType>& pixelPosition,
                     bool normalize = true) const;


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

      FloatType m_centerU;
      FloatType m_centerV;
      FloatType m_focalLength;
      FloatType m_kX;
      FloatType m_kY;
      unsigned int m_numPixelsX;
      unsigned int m_numPixelsY;

    };


    /** 
     * This function outputs a text representation of a
     * CameraIntrinsicsPinhole instance to a std::ostream.  The output
     * format looks like this:
     *
     * CameraIntrinsicsPinhole {240.0, 320.0, 2000.0, 3000.0, 640, 480}
     *
     * @param stream This argument is a reference to the the output
     * stream.
     *
     * @param intrinsics This argument is a const reference to the
     * CameraIntrinsicsPinhole instance to be output.
     *
     * @return The return value is a reference to the input stream after
     * the write has taken place.
     */
    template <class FloatType>
    inline std::ostream&
    operator<<(std::ostream& stream,
               const CameraIntrinsicsPinhole<FloatType>& intrinsics)
    {
      return intrinsics.writeToStream(stream);
    }


    /** 
     * This function sets the value of a CameraIntrinsicsPinhole instance from a
     * std::istream.  The input format is as described for
     * operator<<(std::ostream&, const CameraIntrinsicsPinhole&) above.
     * 
     * @param stream This argument is a reference to the the input
     * stream from which to read.
     *
     * @param intrinsics This argument is a reference to the
     * CameraIntrinsicsPinhole which will take the input.
     *
     * @return The return value is a reference to the input stream after
     * the read has taken place.
     */
    template <class FloatType>
    inline std::istream&
    operator>>(std::istream& stream,
               CameraIntrinsicsPinhole<FloatType>& intrinsics)
    {
      return intrinsics.readFromStream(stream);
    }

    
  } // namespace computerVision
  
} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/cameraIntrinsicsPinhole_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_CAMERAINTRINSICSPINHOLE_HH */
