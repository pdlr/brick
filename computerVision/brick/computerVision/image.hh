/**
***************************************************************************
* @file brick/computerVision/image.hh
*
* Header file declaring Image class.
*
* Copyright (C) 2005-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_IMAGE_HH
#define BRICK_COMPUTERVISION_IMAGE_HH

#include <brick/computerVision/imageFormatTraits.hh>
#include <brick/numeric/array2D.hh>

namespace brick {

  /**
   ** This namespace is still in the early stages of development. Feel
   ** free to play around with it, but please bear in mind that its
   ** interface is not stable.
   **/
  namespace computerVision {
    
    /**
     ** This class template represents a 2D image.  The template
     ** parameter indicates the format of the image, such as GRAY8,
     ** RGB8, YUV420, etc.  Please see the ImageFormat enum for legal
     ** values.  Please see the ImageFormatTraits class template for the
     ** characteristics of the available image formats.
     **
     ** Here are some examples of how you might use this class.
     **
     ** @code
     **   // Create an image of 480 rows and 640 columns.  The
     **   // resulting image will be reference counted, and will
     **   // manage its own memory allocation/deallocation.
     **   Image<RGB8> image0(480, 640);
     **
     **   // Set every pixel to red.
     **   image0 = PixelRGB8(255, 0, 0);
     **
     **   // Shallow copy the image.  The copy remains valid even
     **   // after image0 is destroyed.
     **   Image<RGB8> image1 = image0;
     **
     **   // Deep copy the image.
     **   Image<RGB8> image2 = image0.copy();
     **
     **   // Copy a region of the image.
     **   const unsigned int startRow = 10;
     **   const unsigned int startColum = 20;
     **   const unsigned int stopRow = 10;
     **   const unsigned int stopColum = 20;
     **   for(unsigned int rr = startRow; rr < stopRow; ++rr) {
     **     for(unsigned int cc = startColumn; cc < stopColumn; ++cc) {
     **       image2(rr, cc) = image1(rr, cc);
     **     }
     **   }
     **
     **   // Convert to a different format.
     **   Image<GRAY16> grayImage0 = convertColorspace<GRAY16>(image0);
     ** @endcode
     **
     ** More examples TBD.
     **/
    template <ImageFormat FORMAT>
    class Image
      : public brick::numeric::Array2D<
          typename ImageFormatTraits<FORMAT>::PixelType>
    {
    public:

      /* ******** Public typedefs ******** */

      typedef typename ImageFormatTraits<FORMAT>::PixelType PixelType;


      /* ******** Public member functions ******** */

      /** 
       * Default constructor initializes to zero size.
       */
      Image()
        : brick::numeric::Array2D<PixelType>() {}

    
      /** 
       * Constructs a "rows x columns" element image.
       * 
       * @param numRows Number of rows in the image after successful
       * construction.
       *
       * @param numColumns Number of columns in the image after successful
       * construction.
       */
      Image(size_t numRows, size_t numColumns)
        : brick::numeric::Array2D<PixelType>(numRows, numColumns) {}
    

      /** 
       * The copy constructor does a shallow copy.  The newly created
       * image points to the same data as copied image.
       * 
       * @param source The Image instance to be copied.
       */
      Image(const Image<FORMAT> &source)
        : brick::numeric::Array2D<PixelType>(source) {}

    
      /** 
       * This constructor allows us to implicitly make an Image instance
       * from an Array2D.  As with the copy constructor, the newly
       * created image points to the same data as copied array.
       * 
       * @param source The Array2D instance to be copied.
       */
      Image(const brick::numeric::Array2D<PixelType> &source)
        : brick::numeric::Array2D<PixelType>(source) {}

    
      /**
       * Construct an image around external data.  Images constructed in
       * this way will not implement reference counting, and will not
       * delete dataPtr when done.  The elements of the Image are
       * generally organized in row-major order, however special case
       * formats such as YUV420 may specify their own ordering.
       * 
       * @param numRows Number of rows in the image after successful
       * construction.
       *
       * @param numColumns Number of columns in the image after successful
       * construction.
       *
       * @param dataPtr A C-style array of PixelType into which the newly
       * constructed Image should index.
       */
      Image(size_t numRows, size_t numColumns, PixelType* const dataPtr)
        : brick::numeric::Array2D<PixelType>(numRows, numColumns, dataPtr) {}      


      /**
       * Construct an image around external data with reference
       * counting.  Images constructed in this way will implement
       * reference counting, and will delete dataPtr when done.  This
       * constructor is provide for ease of interaction with Array2D
       * classes.
       * 
       * @param numRows Number of rows in the image after successful
       * construction.
       *
       * @param numColumns Number of columns in the image after successful
       * construction.
       *
       * @param dataPtr A C-style array of PixelType into which the newly
       * constructed Image should index.
       *
       * @param referenceCountPtr A pointer to the associated reference
       * count.
       */
      Image(size_t numRows, size_t numColumns, PixelType* const dataPtr,
            size_t* referenceCountPtr)
        : brick::numeric::Array2D<PixelType>(
          numRows, numColumns, dataPtr, referenceCountPtr) {}      


      /**
       * Destroys the Image instance and deletes the internal data
       * store if no remaining images point to it.
       */
      virtual
      ~Image() {}


      /** 
       * This assignment operator copies its argument into each pixel of
       * the image.  It is provided avoid an implicit cast when using
       * the corresponding Array2D operator.
       * 
       * @param value This argument is the value to be copied.
       * 
       * @return The return value is a reference to *this.
       */
      virtual
      Image<FORMAT>
      operator=(const PixelType& value) {
        return brick::numeric::Array2D<PixelType>::operator=(value);
      }

    
    private:


    };

  } // namespace computerVision
    
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_IMAGE_HH */
