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
    

#if HAVE_LIBPNG

#include <png.h>

    template <ImageFormat FORMAT>
    Image<FORMAT>
    readPNG8(const std::string& /* fileName */,
             std::string& /* commentString */)
    {
      BRICK_THROW(brick::common::NotImplementedException, "readPNG8()",
                  "This function is not yet implemented for the specified "
                  "image format.");
      return Image<FORMAT>();
    }


    template <>
    Image<GRAY8>
    readPNG8(const std::string& fileName,
             std::string& /* commentString */)
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

      Image<GRAY8> result;

      // We'll check eight bytes of magic at the beginning of the file
      // to make sure it's actually a png image.
      const size_t pngSignatureSize = 8;

      FILE *fp = fopen(fileName.c_str(), "rb");
      if (fp == 0) {
        BRICK_THROW(brick::common::IOException,
                    "dlr::computerVision::readPng8()",
                    "Couldn't open input file.");
      }

      // Be sure to clean up the open file.
      try {

        // Read and check the png magic to see if we have an actual
        // png image.
        unsigned char header[pngSignatureSize + 1];
        if(fread(header, 1, pngSignatureSize, fp) != pngSignatureSize) {
          BRICK_THROW(brick::common::IOException,
                      "brick::computerVision::readPng8()",
                      "Couldn't read png signature from file.");
        }
        if(png_sig_cmp(header, 0, pngSignatureSize) != 0) {
          BRICK_THROW(brick::common::IOException,
                      "brick::computerVision::readPng8()",
                      "File doesn't seem to be a PNG image.");
        }
        
        // Create and initialize the png_struct with the default stderr
        // and longjump error functions.
        png_structp pngPtr = png_create_read_struct(
          PNG_LIBPNG_VER_STRING, 0, 0, 0);
        if(pngPtr == 0) {
          BRICK_THROW(brick::common::RunTimeException,
                      "dlr::computerVision::readPng8()",
                      "Couldn't initialize png_structp.");
        }

        // This variable is just to let us choreograph nicely with libpng
        // cleanup functions.
        png_infop* infoPtrPtrForCleanup = png_infopp_NULL;

        // Be sure to clean up pngPtr (and later, infoPtr).
        try {

          // Allocate/initialize the memory for image information.
          png_infop infoPtr = png_create_info_struct(pngPtr);
          if (infoPtr == 0) {
            BRICK_THROW(brick::common::RunTimeException,
                        "brick::computerVision::readPng8()",
                        "Couldn't initialize png_infop.");
          }
          infoPtrPtrForCleanup = &infoPtr;
       
          // Set error handling in case libpng calls longjmp().
          if(setjmp(png_jmpbuf(pngPtr))) {
            BRICK_THROW(brick::common::IOException,
                        "dlr::computerVision::readPng8()",
                        "Trouble reading from file.");
          }
        
          // Set up the input control.
          png_init_io(pngPtr, fp);

          // Let libpng know that we've already checked some magic.
          png_set_sig_bytes(pngPtr, pngSignatureSize);
          
          // Read entire image into pngPtr.
          png_read_png(pngPtr, infoPtr, PNG_TRANSFORM_IDENTITY, png_voidp_NULL);

          // Find out about our image.
          png_uint_32 width, height;
          int bitDepth, colorType, interlaceType;
          int compressionType, filterMethod;
          png_get_IHDR(pngPtr, infoPtr, &width, &height, &bitDepth, &colorType,
                       &interlaceType, &compressionType, &filterMethod);
          if(bitDepth != 8) {
            BRICK_THROW(brick::common::IOException,
                        "brick::computerVision::readPng8()",
                        "Image file is not 8 bit.");
          }
          if(colorType != PNG_COLOR_TYPE_GRAY) {
            BRICK_THROW(brick::common::IOException,
                        "brick::computerVision::readPng8()",
                        "Can currently only handle grayscale images.");
          }
          if(interlaceType != PNG_INTERLACE_NONE) {
            BRICK_THROW(brick::common::IOException,
                        "brick::computerVision::readPng8()",
                        "Can currently only handle non-interlaced images.");
          }

          // Copy image and return.
          result.reinit(height, width);
          png_bytep* rowPointers = png_get_rows(pngPtr, infoPtr);
          for(size_t rowIndex = 0; rowIndex < height; ++rowIndex) {
            result.getRow(rowIndex).copy(rowPointers[rowIndex]);
          }

        } catch(...) {
          png_destroy_read_struct(
            &pngPtr, infoPtrPtrForCleanup, png_infopp_NULL);
          throw;
        }
        
        // clean up after the read, and free any memory allocated.
        png_destroy_read_struct(&pngPtr, infoPtrPtrForCleanup, png_infopp_NULL);

      } catch(...) {
        fclose(fp);
        throw;
      }
      fclose(fp);

      return result;
    }


    template <>
    Image<RGB8>
    readPNG8(const std::string& fileName,
             std::string& /* commentString */)
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

      if(!PixelRGB8::isContiguous()) {
        BRICK_THROW(brick::common::NotImplementedException,
                    "readPNG8RGB()",
                    "This function currently only works with compilers that "
                    "don't add padding to the PixelRGB8 memory layout.");
      }

      Image<RGB8> result;

      // We'll check eight bytes of magic at the beginning of the file
      // to make sure it's actually a png image.
      const size_t pngSignatureSize = 8;

      FILE *fp = fopen(fileName.c_str(), "rb");
      if (fp == 0) {
        BRICK_THROW(brick::common::IOException,
                    "dlr::computerVision::readPng8()",
                    "Couldn't open input file.");
      }

      // Be sure to clean up the open file.
      try {

        // Read and check the png magic to see if we have an actual
        // png image.
        unsigned char header[pngSignatureSize + 1];
        if(fread(header, 1, pngSignatureSize, fp) != pngSignatureSize) {
          BRICK_THROW(brick::common::IOException,
                      "brick::computerVision::readPng8()",
                      "Couldn't read png signature from file.");
        }
        if(png_sig_cmp(header, 0, pngSignatureSize) != 0) {
          BRICK_THROW(brick::common::IOException,
                      "brick::computerVision::readPng8()",
                      "File doesn't seem to be a PNG image.");
        }
        
        // Create and initialize the png_struct with the default stderr
        // and longjump error functions.
        png_structp pngPtr = png_create_read_struct(
          PNG_LIBPNG_VER_STRING, 0, 0, 0);
        if(pngPtr == 0) {
          BRICK_THROW(brick::common::RunTimeException,
                      "dlr::computerVision::readPng8()",
                      "Couldn't initialize png_structp.");
        }

        // This variable is just to let us choreograph nicely with libpng
        // cleanup functions.
        png_infop* infoPtrPtrForCleanup = png_infopp_NULL;

        // Be sure to clean up pngPtr (and later, infoPtr).
        try {

          // Allocate/initialize the memory for image information.
          png_infop infoPtr = png_create_info_struct(pngPtr);
          if (infoPtr == 0) {
            BRICK_THROW(brick::common::RunTimeException,
                        "brick::computerVision::readPng8()",
                        "Couldn't initialize png_infop.");
          }
          infoPtrPtrForCleanup = &infoPtr;
       
          // Set error handling in case libpng calls longjmp().
          if(setjmp(png_jmpbuf(pngPtr))) {
            BRICK_THROW(brick::common::IOException,
                        "dlr::computerVision::readPng8()",
                        "Trouble reading from file.");
          }
        
          // Set up the input control.
          png_init_io(pngPtr, fp);

          // Let libpng know that we've already checked some magic.
          png_set_sig_bytes(pngPtr, pngSignatureSize);
          
          // Read entire image into pngPtr.
          png_read_png(pngPtr, infoPtr, PNG_TRANSFORM_IDENTITY, png_voidp_NULL);

          // Find out about our image.
          png_uint_32 width, height;
          int bitDepth, colorType, interlaceType;
          int compressionType, filterMethod;
          png_get_IHDR(pngPtr, infoPtr, &width, &height, &bitDepth, &colorType,
                       &interlaceType, &compressionType, &filterMethod);
          if(bitDepth != 8) {
            BRICK_THROW(brick::common::IOException,
                        "brick::computerVision::readPng8()",
                        "Image file is not 8 bit.");
          }
          /*if(colorType != PNG_COLOR_TYPE_GRAY) {
            BRICK_THROW(brick::common::IOException,
            "brick::computerVision::readPng8()",
            "Can currently only handle grayscale images.");
            }*/
          if(interlaceType != PNG_INTERLACE_NONE) {
            BRICK_THROW(brick::common::IOException,
                        "brick::computerVision::readPng8()",
                        "Can currently only handle non-interlaced images.");
          }

          // Copy image and return.
          result.reinit(height, width);
          png_bytep* rowPointers = png_get_rows(pngPtr, infoPtr);
          for(size_t rowIndex = 0; rowIndex < height; ++rowIndex) {
            std::copy(rowPointers[rowIndex], rowPointers[rowIndex] + width * 3, 
                      reinterpret_cast<common::UInt8*>(
                        result.getRow(rowIndex).data()));
          }

        } catch(...) {
          png_destroy_read_struct(
            &pngPtr, infoPtrPtrForCleanup, png_infopp_NULL);
          throw;
        }
        
        // clean up after the read, and free any memory allocated.
        png_destroy_read_struct(&pngPtr, infoPtrPtrForCleanup, png_infopp_NULL);

      } catch(...) {
        fclose(fp);
        throw;
      }
      fclose(fp);

      return result;
    }
    
#endif /* #if HAVE_LIBPNG */

    
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
