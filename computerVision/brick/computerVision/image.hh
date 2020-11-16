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
#include <brick/numeric/index2D.hh>

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
       * The copy assignment operator does a shallow copy.  The newly
       * created image points to the same data as copied image.
       *
       * @param source The Image instance to be copied.
       */
      Image&
      operator=(const Image<FORMAT> &source)
      {
          brick::numeric::Array2D<PixelType>::operator=(source);
          return *this;
      }


      /**
       * This copy assignment operator allows us to implicitly make an
       * Image instance from an Array2D.  As with the copy assignment
       * operator, the newly created image points to the same data as
       * copied array.
       *
       * @param source The Array2D instance to be copied.
       */
      Image&
      operator=(const brick::numeric::Array2D<PixelType> &source)
      {
          brick::numeric::Array2D<PixelType>::operator=(source);
          return *this;
      }


      /**
       * Returns an image that references only a rectangular region of
       * interest drawn from *this.  The returned image will reference
       * the same memory, but may have different start, end, and
       * rowstep.  The returned image _will not_ be reference counted,
       * so that its contents are valid only as long as the contents
       * of the original image are valid.  The operation of this
       * function is illustrated in the following example:
       *
       * @code
       *   // Create an array in which every pixel is set to 1.
       *   Image<GRAY8>* fullImageP = new Image<GRAY8>(10, 20);
       *   *fullImageP = 1;
       *
       *   // Get an image that refers to a 4x4 subset of fullImage.
       *   Index2D corner0(5, 10);
       *   Index2D corner1(9, 14);
       *   Image<GRAY8> subRegion = fullImageP->getROI(corner0, corner1);
       *
       *   // Set just the pixels within that region of interest to
       *   // 100.
       *   subRegion = 100;
       *
       *   for(unsigned int rr = 0, rr < fullImage.rows(); ++rr) {
       *     for(unsigned int cc = 0, cc < fullImage.columns(); ++cc) {
       *       if((rr >= corner0.getRow()) &&
       *          (rr <  corner1.getRow()) &&
       *          (cc >= corner0.getColumn()) &&
       *          (cc <  corner1.getColumn())) {
       *         std::assert(100 == (*fullImageP)(rr, cc));
       *       } else {
       *         std::assert(1   == (*fullImageP)(rr, cc));
       *       }
       *     }
       *   }
       *
       *   // Contents of *fullImageP are deleted here.
       *   delete fullImageP;
       *
       *   // ERROR! The data to which subRegion refers was deleted
       *   // along with fullImageP.
       *   subRegion(2, 3) = 90;
       * @endCode
       *
       * @param corner0 This argument and the next define the
       * subregion.  Corner0 is included in the region, but corner1 is
       * not.  It is an error if corner0 is not above and to the left
       * of corner1.  It is an error if corner0.getRow() is less than
       * zero, corner0.getColumn() is less than zero, corner0.getRow()
       * is greater than this->getRows(), or corner0.getColumn() is
       * greater than this->getColumns().
       *
       * @param corner1 This argument and the previous define the
       * subregion.  It is an error if corner1.getRow() is less than
       * zero, corner1.getColumn() is less than zero, corner1.getRow()
       * is greater than this->getRows(), or corner1.getColumn() is
       * greater than this->getColumns().
       *
       * @return The return value is a shallow copy of the selected
       * region of the original image.
       */
      Image<FORMAT>
      getROI(brick::numeric::Index2D const& corner0,
             brick::numeric::Index2D const& corner1) {
        return Image<FORMAT>(this->getRegion(corner0, corner1));
      }


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
