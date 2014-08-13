/**
***************************************************************************
* @file brick/computerVision/imageIO_impl.hh
*
* Header file defining inline and template functions for reading and
* writing images.
*
* Copyright (C) 2004-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_IMAGEIO_IMPL_HH
#define BRICK_COMPUTERVISION_IMAGEIO_IMPL_HH

#include <cmath>
#include <fstream>
#include <sstream>
#include <brick/common/byteOrder.hh>
#include <brick/common/exception.hh>
#include <brick/computerVision/utilities.hh>

namespace brick {

  namespace computerVision {

    /// @cond privateCode
    namespace privateCode {

      template<class Type0, class Type1>
      const Type1*
      copyPNMData(const Type0* sourcePtr, size_t numberOfElements);

  
      template<class Type0, class Type1>
      void
      deletePNMData(const Type0* outputBuffer, const Type1* sourcePtr);

  
      template<class Type0, class Type1>
      const Type1*
      normalizePNMData(const Type0* sourcePtr, size_t numberOfElements);

  
      template<class Type>
      void
      writeRawPGM(const std::string& fileName,
                  const Type* imageData,
                  const size_t rows,
                  const size_t columns);

  
      template<class Type>
      void
      writeRawPPM(const std::string& fileName,
                  const Type* imageData,
                  const size_t rows,
                  const size_t columns);

  
      template<class Type>
      void
      writePlainPGM(const std::string& fileName,
                    const Type* imageData,
                    const size_t rows,
                    const size_t columns);

      
      template<class Type>
      void
      writePlainPPM(const std::string& fileName,
                    const Type* imageData,
                    const size_t rows,
                    const size_t columns);

    } // namespace privateCode
    /// @endcond
    

    template<class Type>
    void
    writePGM(const std::string& fileName,
             const Type* imageData,
             const size_t rows,
             const size_t columns,
             bool normalize,
             bool rawFormat,
             int bitsPerPixel)
    {
      // Check arguments.
      if((bitsPerPixel != 8) && (bitsPerPixel != 16)) {
        BRICK_THROW(brick::common::ValueException, "writePGM(...)",
                  "Argument bitsPerPixel must be either 8 or 16.");
      }

      size_t numberOfElements = rows * columns;
      if(bitsPerPixel == 8) {
        const unsigned char* charBuffer;

        // Normalize if necessary.
        if(normalize) {
          charBuffer = privateCode::normalizePNMData<Type, brick::common::UnsignedInt8>(
            imageData, numberOfElements);

        } else {
          charBuffer = privateCode::copyPNMData<Type, brick::common::UnsignedInt8>(
            imageData, numberOfElements);
        }
      
        // Dispatch to the appropriate IO routine.
        if(rawFormat) {
          privateCode::writeRawPGM(fileName, charBuffer, rows, columns);
        } else {
          privateCode::writePlainPGM(fileName, charBuffer, rows, columns);
        }

        // And clean up.
        privateCode::deletePNMData(charBuffer, imageData);
      } else if(bitsPerPixel == 16) {
        const unsigned short* shortBuffer;

        // Normalize if necessary.
        if(normalize) {
          shortBuffer = privateCode::normalizePNMData<Type, brick::common::UnsignedInt16>(
            imageData, numberOfElements);
        } else {
          shortBuffer = privateCode::copyPNMData<Type, brick::common::UnsignedInt16>(
            imageData, numberOfElements);
        }
      
        // Dispatch to the appropriate IO routine.
        if(rawFormat) {
          privateCode::writeRawPGM(fileName, shortBuffer, rows, columns);
        } else {
          privateCode::writePlainPGM(fileName, shortBuffer, rows, columns);
        }

        // And clean up.
        privateCode::deletePNMData(shortBuffer, imageData);
      }
    }


    template<class Type>
    void
    writePPM(const std::string& fileName,
             const Type* imageData,
             const size_t rows,
             const size_t columns,
             bool normalize,
             bool rawFormat,
             int bitsPerPixel)
    {
      // Check arguments.
      if((bitsPerPixel != 8) && (bitsPerPixel != 16)) {
        BRICK_THROW(brick::common::ValueException, "writePGM(...)",
                  "Argument bitsPerPixel must be either 8 or 16.");
      }

      size_t numberOfElements = rows * columns * 3;
      if(bitsPerPixel == 8) {
        const unsigned char* charBuffer;

        // Normalize if necessary.
        if(normalize) {
          charBuffer = privateCode::normalizePNMData<Type, brick::common::UnsignedInt8>(
            imageData, numberOfElements);
        } else {
          charBuffer = privateCode::copyPNMData<Type, brick::common::UnsignedInt8>(
            imageData, numberOfElements);
        }
      
        // Dispatch to the appropriate IO routine.
        if(rawFormat) {
          privateCode::writeRawPPM(fileName, charBuffer, rows, columns);
        } else {
          privateCode::writePlainPPM(fileName, charBuffer, rows, columns);
        }

        // And clean up.
        privateCode::deletePNMData(charBuffer, imageData);
      } else {
        const unsigned short* shortBuffer;

        // Normalize if necessary.
        if(normalize) {
          shortBuffer = privateCode::normalizePNMData<Type, brick::common::UnsignedInt16>(
            imageData, numberOfElements);
        } else {
          shortBuffer = privateCode::copyPNMData<Type, brick::common::UnsignedInt16>(
            imageData, numberOfElements);
        }
      
        // Dispatch to the appropriate IO routine.
        if(rawFormat) {
          privateCode::writeRawPPM(fileName, shortBuffer, rows, columns);
        } else {
          privateCode::writePlainPPM(fileName, shortBuffer, rows, columns);
        }

        // And clean up.
        privateCode::deletePNMData(shortBuffer, imageData);
      }
    }

  } // namespace computerVision

} // namespace brick


#if HAVE_LIBPNG

#include <brick/computerVision/pngReader.hh>

namespace brick {

  namespace computerVision {

    namespace privateCode {

      inline Image<GRAY8>
      fixPngEndianness(Image<GRAY8> const& outputImage) {
        return outputImage;
      }

      inline Image<GRAY16>
      fixPngEndianness(Image<GRAY16> const& outputImage) {
        if(brick::common::getByteOrder() != brick::common::BRICK_BIG_ENDIAN) {
          Image<GRAY16> swappedImage = outputImage.copy();
          brick::common::switchByteOrder(
            swappedImage.data(), swappedImage.size(),
            brick::common::getByteOrder(), brick::common::BRICK_BIG_ENDIAN);
          return swappedImage;
        }
        return outputImage;
      }

      inline Image<RGB8>
      fixPngEndianness(Image<RGB8> const& outputImage) {
        return outputImage;
      }

      inline Image<RGB16>
      fixPngEndianness(Image<RGB16> const& outputImage) {

        // Set up a plausible return value, just in case no swapping
        // is necessary.  This is a shallow copy.
        Image<RGB16> swappedImage = outputImage;

        // Do we need to swap bytes?
        if(brick::common::getByteOrder() != brick::common::BRICK_BIG_ENDIAN) {

          // Yes.  Get access to an array of UInt16 (rather than RGB16
          // structs).
          typedef ImageFormatTraits<RGB16>::ComponentType ComponentType;
          brick::numeric::Array2D<ComponentType> componentArray;
          if(brick::computerVision::dissociateColorComponents(
               const_cast<Image<RGB16>&>(outputImage), componentArray)) {

            // The array of UInt16 shares data with outputImage.  Copy
            // it before swapping bytes.
            componentArray = componentArray.copy();

          }

          // Swap bytes.
          brick::common::switchByteOrder(
            componentArray.data(), componentArray.size(),
            brick::common::getByteOrder(), brick::common::BRICK_BIG_ENDIAN);

          // Repack into an RGB image.
          // 
          // TBD(xxx): We could skip this step (and save a copy) by
          // changing the return value of fixPngEndianness to
          // Array2D<ComponentType>.
          
          if(brick::computerVision::associateColorComponents(
               componentArray, swappedImage)) {

            // If swappedImage shares data with componentArray, and
            // this data will become invalid as soon as componentArray
            // goes out of scope.  Make a copy to avoid this.
            swappedImage = swappedImage.copy();

          }
          
        } // if(getByteOrder...)
        
        return swappedImage;
      }

      
      inline void
      setPngHeaderInfo(png_structp pngPtr,
                       png_infop infoPtr,
                       Image<GRAY8> const& outputImage) {
        png_set_IHDR(
          pngPtr, infoPtr, outputImage.columns(), outputImage.rows(),
          sizeof(ImageFormatTraits<GRAY8>::PixelType) * 8,
          PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
          PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
      }
      
      inline void
      setPngHeaderInfo(png_structp pngPtr,
                       png_infop infoPtr,
                       Image<GRAY16> const& outputImage) {

#if 0
        // Seems like this should be a sane way to handle endianness,
        // but it doesn't appear to work.
        
        // PNG files are natively big-endian, but can be coerced to
        // store little-endian data without a swap.
        if(brick::common::getByteOrder()
           == brick::common::BRICK_LITTLE_ENDIAN) {
          png_set_swap(pngPtr);
        }
#endif /* #if 0 */

        png_set_IHDR(
          pngPtr, infoPtr, outputImage.columns(), outputImage.rows(),
          sizeof(ImageFormatTraits<GRAY16>::PixelType) * 8,
          PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
          PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
      }
      
      inline void
      setPngHeaderInfo(png_structp pngPtr,
                       png_infop infoPtr,
                       Image<RGB8> const& outputImage) {
        if(!Image<RGB8>::value_type::isContiguous()) {
          BRICK_THROW(
            brick::common::NotImplementedException, "writePNG()",
            "Your compiler appears to pack RGB pixels in an unusual way.  "
            "Please update brick::computerVision::writePNG() to handle this.");
        }
        png_set_IHDR(
          pngPtr, infoPtr, outputImage.columns(), outputImage.rows(),
          sizeof(ImageFormatTraits<GRAY8>::PixelType) * 8,
          PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
          PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
      }

      inline void
      setPngHeaderInfo(png_structp pngPtr,
                       png_infop infoPtr,
                       Image<RGB16> const& outputImage) {
        if(!Image<RGB16>::value_type::isContiguous()) {
          BRICK_THROW(
            brick::common::NotImplementedException, "writePNG()",
            "Your compiler appears to pack RGB pixels in an unusual way.  "
            "Please update brick::computerVision::writePNG() to handle this.");
        }
        
#if 0
        // Seems like this should be a sane way to handle endianness,
        // but it doesn't appear to work.
        
        // PNG files are natively big-endian, but can be coerced to
        // store little-endian data without a swap.
        if(brick::common::getByteOrder()
           == brick::common::BRICK_LITTLE_ENDIAN) {
          png_set_swap(pngPtr);
        }
#endif /* #if 0 */

        // Pick reasonable defaults for the rest of the header
        // configuration.
        png_set_IHDR(
          pngPtr, infoPtr, outputImage.columns(), outputImage.rows(),
          sizeof(ImageFormatTraits<GRAY16>::PixelType) * 8,
          PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
          PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
      }

    } // namespace privateCode
      
    
    template <ImageFormat Format>
    Image<Format>
    readPNG(const std::string& fileName,
             std::string& /* commentString */)
    {
      PngReader pngReader(fileName);
      return pngReader.getImage<Format>();
    }


    template<ImageFormat Format>
    void
    writePNG(const std::string& fileName,
             const Image<Format>& outputImage,
             const std::string& /* comment */)
    {
      // This code is heavily in debt to example.c from the libpng 1.2.1
      // distribition, which carries the following header comment:
      // /* example.c - an example of using libpng
      //  * Last changed in libpng 1.2.1 December 7, 2001.
      //  * This file has been placed in the public domain by the authors.
      //  * Maintained 1998-2001 Glenn Randers-Pehrson
      //  * Maintained 1996, 1997 Andreas Dilger)
      //  * Written 1995, 1996 Guy Eric Schalnat, Group 42, Inc.)
      //  */

      FILE* fp = fopen(fileName.c_str(), "wb");
      if (fp == 0) {
        BRICK_THROW(brick::common::IOException, "ImageIO::writePNG()",
                    "Couldn't open output file.");
      }

      // Make sure we clean up open FILE.
      try {

        // Create and initialize the png_struct with the default stderr and
        // longjump error handler functions.
        png_structp pngPtr = png_create_write_struct(
          PNG_LIBPNG_VER_STRING, 0, 0, 0);
        if (pngPtr == 0) {
          BRICK_THROW(brick::common::RunTimeException, "ImageIO::writePng()",
                      "Couldn't initialize png_structp.");
        }

        // This variable is just to let us choreograph nicely with libpng
        // cleanup functions.
        png_infop* infoPtrPtrForCleanup = png_infopp_NULL;

        // Make sure we clean up pngPtr (and eventually infoPtr).
        try { 

          // We'll need a place to record details (colorspace, etc.)
          // about the image.
          png_infop infoPtr = png_create_info_struct(pngPtr);
          if (infoPtr == 0) {
            BRICK_THROW(brick::common::RunTimeException, "ImageIO::writePng()",
                        "Couldn't initialize png_infop.");
          }
          infoPtrPtrForCleanup = &infoPtr;
     
          // Set error handling in case libpng calls longjmp().
          if(setjmp(png_jmpbuf(pngPtr))) {
            // If we get here, we had a problem reading the file
            BRICK_THROW(brick::common::IOException, "ImageIO::writePng()",
                        "Trouble reading from file.");
          }

          // Set up the output control.
          png_init_io(pngPtr, fp);

          // Note(xxx): Just in case we need it...
          // png_set_compression_level(pngPtr, Z_BEST_COMPRESSION);

          // Seems like endianness should be handled by setting pngPtr
          // options, but that hasn't worked so far, so we brute force
          // it.
          Image<Format> savableImage =
            privateCode::fixPngEndianness(outputImage);
            
          privateCode::setPngHeaderInfo(pngPtr, infoPtr, savableImage);
          
          // Fill in png structure here.
          png_bytep* rowPointers = (png_bytep*)png_malloc(
            pngPtr, outputImage.rows() * png_sizeof(png_bytep));
          if(rowPointers == 0) {
            BRICK_THROW(brick::common::RunTimeException, "ImageIO::writePng()",
                        "Couldn't allocate row pointers.");
          }

          // Be sure to free rowPointers.
          try {
            
            for(size_t rowIndex = 0; rowIndex < outputImage.rows();
                ++rowIndex) {
              rowPointers[rowIndex] =
                (png_bytep)(savableImage.getRow(rowIndex).data());
            }

            png_set_rows(pngPtr, infoPtr, rowPointers);
            png_write_png(
              pngPtr, infoPtr, PNG_TRANSFORM_IDENTITY, png_voidp_NULL);
          } catch(...) {
            png_free(pngPtr, rowPointers);
          }
          png_free(pngPtr, rowPointers);

          /* If you png_malloced a palette, free it here (don't free
             infoPtr->palette, as recommended in versions 1.0.5m and
             earlier of this example; if libpng mallocs infoPtr->palette,
             libpng will free it).  If you allocated it with malloc()
             instead of png_malloc(), use free() instead of
             png_free(). */
          // png_free(pngPtr, palette); xxx;
          // palette=NULL;

        } catch(...) {
          png_destroy_write_struct(&pngPtr, infoPtrPtrForCleanup);
          throw;
        }
     
        /* clean up after the write, and free any memory allocated */
        png_destroy_write_struct(&pngPtr, infoPtrPtrForCleanup);

      } catch(...) {
        fclose(fp);
        throw;
      }

      fclose(fp);
    }


    
  } // namespace computerVision

} // namespace brick


# endif /* #if HAVE_LIBPNG */


namespace brick {

  namespace computerVision {

    /// @cond privateCode
    namespace privateCode {

      template<class Type0, class Type1>
      const Type1*
      copyPNMData(const Type0* sourcePtr, size_t numberOfElements)
      {
        // Make a buffer to hold the normalized image data.
        Type1* outputBuffer = new Type1[numberOfElements];

        // Convert to Type1.
        for(size_t index0 = 0; index0 < numberOfElements; ++index0) {
          outputBuffer[index0] = static_cast<Type1>(sourcePtr[index0]);
        }

        // All done.
        return outputBuffer;
      }
  

      template<class Type0, class Type1>
      void
      deletePNMData(const Type0* outputBuffer, const Type1* sourcePtr)
      {
        // Eventually, we might make copyPNMData() simply return the input
        // pointer if the input pointer has the right type.  Check to make
        // sure that didn't happen.
        if(reinterpret_cast<const Type1*>(outputBuffer) == sourcePtr) {
          return;
        }

        // On with the deletion!
        delete[] outputBuffer;
      }

  
      template<class Type0, class Type1>
      const Type1*
      normalizePNMData(const Type0* sourcePtr, size_t numberOfElements)
      {
        // Compute the largest acceptable value for the target type.
        const int imageMax =
          static_cast<int>(std::pow(2.0, 8.0 * sizeof(Type1)) - 0.5);

        // Make a buffer to hold the normalized image data.
        Type1* outputBuffer = new Type1[numberOfElements];

        // Figure out the range of pixel values.
        Type0 maxValue =
          *std::max_element(sourcePtr, sourcePtr + numberOfElements);
        Type0 minValue =
          *std::min_element(sourcePtr, sourcePtr + numberOfElements);

        // Try to avoid dividing by zero.
        double scaleFactor = 0.0;
        if(maxValue != minValue) {
          scaleFactor = imageMax / (maxValue - minValue);      
        }

        // Do the actual normalization.
        for(size_t index0 = 0; index0 < numberOfElements; ++index0) {
          outputBuffer[index0] = static_cast<Type1>(
            (sourcePtr[index0] - minValue) * scaleFactor);
        }

        // All done.
        return outputBuffer;
      }

  
      template<class Type>
      void
      writeRawPGM(const std::string& fileName,
                  const Type* imageData,
                  const size_t rows,
                  const size_t columns)
      {
        // Compute the largest acceptable value for the target type.
        const int imageMax =
          static_cast<int>(std::pow(2.0, 8.0 * sizeof(Type)) - 0.5);
    
        // Open the file.
        std::ofstream outputStream(fileName.c_str(), std::ios::binary);
        if(!outputStream) {
          std::ostringstream message;
          message << "Couldn't open output file: " << fileName;
          BRICK_THROW(brick::common::IOException, "writeRawPGM()", message.str().c_str());
        }

        // Write the header.
        outputStream << "P5\n"
                     << columns << " " << rows << "\n"
                     << imageMax << "\n";

        // Write the image data.
        size_t numberOfElements = rows * columns;
        if(sizeof(Type) == 1) {
          outputStream.write(
            (char*)imageData, static_cast<std::streamsize>(numberOfElements));
        } else if(sizeof(Type) == 2) {
          size_t numberOfBytes = 2 * numberOfElements;
          if(brick::common::getByteOrder() != brick::common::BRICK_BIG_ENDIAN) {
            char* swabbedData = new char[numberOfBytes];
            switchByteOrder(
              imageData, numberOfElements,
              reinterpret_cast<Type*>(swabbedData),
              brick::common::getByteOrder(), brick::common::BRICK_BIG_ENDIAN);
            outputStream.write(
              swabbedData, static_cast<std::streamsize>(numberOfBytes));
            delete[] swabbedData;
          } else {
            outputStream.write(
              (char*)imageData, static_cast<std::streamsize>(numberOfBytes));
          }
        } else {
          outputStream.close();
          std::ostringstream message;
          message << "Can't write a pgm with " << sizeof(Type)
                  << "-byte values.";
          BRICK_THROW(brick::common::ValueException, "writeRawPGM()", message.str().c_str());
        }      
        
        // Check for errors.
        if(!outputStream) {
          outputStream.close();
          std::ostringstream message;
          message << "Error writing to open output file: " << fileName;
          BRICK_THROW(brick::common::IOException, "writeRawPGM()", message.str().c_str());
        }
    
        // All done!
        outputStream.close();
      }

  
      template<class Type>
      void
      writeRawPPM(const std::string& fileName,
                  const Type* imageData,
                  const size_t rows,
                  const size_t columns)
      {
        // Compute the largest acceptable value for the target type.
        const int imageMax =
          static_cast<int>(std::pow(2.0, 8.0 * sizeof(Type)) - 0.5);
    
        // Open the file.
        std::ofstream outputStream(fileName.c_str(), std::ios::binary);
        if(!outputStream) {
          std::ostringstream message;
          message << "Couldn't open output file: " << fileName;
          BRICK_THROW(brick::common::IOException, "writeRawPPM()", message.str().c_str());
        }

        // Write the header.
        outputStream << "P6\n"
                     << columns << " " << rows << "\n"
                     << imageMax << "\n";

        // Write the image data.
        size_t numberOfElements = rows * columns * 3;
        if(sizeof(Type) == 1) {
          outputStream.write(
            (char*)imageData, static_cast<std::streamsize>(numberOfElements));
        } else if(sizeof(Type) == 2) {
          size_t numberOfBytes = 2 * numberOfElements;
          if(brick::common::getByteOrder() != brick::common::BRICK_BIG_ENDIAN) {
            char* swabbedData = new char[numberOfBytes];
            switchByteOrder(
              imageData, numberOfElements,
              reinterpret_cast<Type*>(swabbedData),
              brick::common::getByteOrder(), brick::common::BRICK_BIG_ENDIAN);
            outputStream.write(
              swabbedData, static_cast<std::streamsize>(numberOfBytes));
            delete[] swabbedData;
          } else {
            outputStream.write(
              (char*)imageData, static_cast<std::streamsize>(numberOfBytes));
          }
        } else {
          outputStream.close();
          std::ostringstream message;
          message << "Can't write a ppm with " << sizeof(Type)
                  << "-byte values.";
          BRICK_THROW(brick::common::ValueException, "writeRawPPM()", message.str().c_str());
        }      
        
        // Check for errors.
        if(!outputStream) {
          std::ostringstream message;
          message << "Error writing to open output file: " << fileName;
          BRICK_THROW(brick::common::IOException, "writeRawPPM()", message.str().c_str());
        }
    
        // All done!
        outputStream.close();
      }

  
      template<class Type>
      void
      writePlainPGM(const std::string& fileName,
                    const Type* imageData,
                    const size_t rows,
                    const size_t columns)
      {
        // Compute the largest acceptable value for the target type.
        const int imageMax =
          static_cast<int>(std::pow(2.0, 8.0 * sizeof(Type)) - 0.5);
    
        // Open the file.
        std::ofstream outputStream(fileName.c_str());
        if(!outputStream) {
          std::ostringstream message;
          message << "Couldn't open output file: " << fileName;
          BRICK_THROW(brick::common::IOException, "writePlainPGM()", message.str().c_str());
        }

        // Write the header.
        outputStream << "P2\n"
                     << columns << " " << rows << "\n"
                     << imageMax << "\n";

        // Write the image data.
        for(size_t row = 0; row < rows; ++row) {
          for(size_t column = 0; column < columns; ++column) {
            outputStream << static_cast<int>(*imageData) << " ";
            ++imageData;
          }
          outputStream << "\n";
        }

        // Check for errors.
        if(!outputStream) {
          std::ostringstream message;
          message << "Error writing to open output file: " << fileName;
          BRICK_THROW(brick::common::IOException, "writePlainPGM()", message.str().c_str());
        }

        // All done!
        outputStream.close();
      }
    

      template<class Type>
      void
      writePlainPPM(const std::string& fileName,
                    const Type* imageData,
                    const size_t rows,
                    const size_t columns)
      {
        // Compute the largest acceptable value for the target type.
        const int imageMax =
          static_cast<int>(std::pow(2.0, 8.0 * sizeof(Type)) - 0.5);
    
        // Open the file.
        std::ofstream outputStream(fileName.c_str());
        if(!outputStream) {
          std::ostringstream message;
          message << "Couldn't open output file: " << fileName;
          BRICK_THROW(brick::common::IOException, "writePlainPPM()", message.str().c_str());
        }

        // Write the header.
        outputStream << "P3\n"
                     << columns << " " << rows << "\n"
                     << imageMax << "\n";

        // Write the image data.
        for(size_t row = 0; row < rows; ++row) {
          for(size_t column = 0; column < columns; ++column) {
            outputStream << static_cast<int>(*imageData) << " "
                         << static_cast<int>(*(imageData + 1)) << " "
                         << static_cast<int>(*(imageData + 2)) << " ";
            imageData += 3;
          }
          outputStream << "\n";
        }

        // Check for errors.
        if(!outputStream) {
          std::ostringstream message;
          message << "Error writing to open output file: " << fileName;
          BRICK_THROW(brick::common::IOException, "writePlainPPM()", message.str().c_str());
        }

        // All done!
        outputStream.close();
      }

    } // namespace privateCode
    /// @endcond
      
  } // namespace computerVision

} // namespace brick


#endif /* #ifndef BRICK_COMPUTERVISION_IMAGEIO_IMPL_HH */
