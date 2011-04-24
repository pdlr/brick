/**
***************************************************************************
* @file brick/utilities/imageIO.hh
*
* Header file declaring routines for importing & exporting images.
* These are bare-minimum utility functions; see brick::computerVision
* for better image IO routines.
*
* Copyright (C) 2004-2011, David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_UTILITIES_IMAGEIO_HH
#define BRICK_UTILITIES_IMAGEIO_HH

#include <string>
#include <brick/common/byteOrder.hh>

namespace brick {

  namespace utilities {

    /** 
     * Write a flat array of data to a PGM image file.
     * 
     * @param fileName This argument is the path of the file to be written.
     * 
     * @param imageData This argument is an array of at least (rows *
     * columns) data elements, one for each pixel in the image,
     * arranged in raster order, left to right, top to bottom.
     * 
     * @param rows This argument indicates how many rows are in the image.
     * 
     * @param columns This argument indicates how many columns are in the image.
     * 
     * @param normalize This argument indicates whether the dynamic
     * range of the input data should be stretched/compressed to match
     * the dynamic range of the target PGM format.
     * 
     * @param rawFormat This argument indicates whether the PGM file
     * should be written in ASCII or raw format.  True indicates raw,
     * false indicates ASCII.
     * 
     * @param bitsPerPixel This argument must be either 8 or 16, and
     * indicates the bit depth of the target file.
     */
    template<class Type>
    void
    writePGM(std::string const& fileName,
             const Type* imageData,
             const size_t rows,
             const size_t columns,
             bool normalize=false,
             bool rawFormat=true,
             int bitsPerPixel=8);


    /** 
     * Write a flat array of data to a PPM image file.
     * 
     * @param fileName This argument is the path of the file to be written.
     * 
     * @param imageData This argument is an array of at least (rows *
     * columns * 3) data elements, one for each color channel of each
     * pixel in the image, arranged in raster order, left to right,
     * top to bottom.  Color chanels are interleaved [red, green,
     * blue, red, green, blue, etc.].
     * 
     * @param rows This argument indicates how many rows are in the image.
     * 
     * @param columns This argument indicates how many columns are in the image.
     * 
     * @param normalize This argument indicates whether the dynamic
     * range of the input data should be stretched/compressed to match
     * the dynamic range of the target PGM format.  The three color
     * channels will be normalized together, so that each channel is
     * scaled and offset by the same amount.
     * 
     * @param rawFormat This argument indicates whether the PGM file
     * should be written in ASCII or raw format.  True indicates raw,
     * false indicates ASCII.
     * 
     * @param bitsPerChannel This argument must be either 8 or 16, and
     * indicates the bit depth of the target file.
     */
    template<class Type>
    void
    writePPM(std::string const& fileName,
             const Type* imageData,
             const size_t rows,
             const size_t columns,
             bool normalize=false,
             bool rawFormat=true,
             int bitsPerChannel=8);


  } // namespace utilities

} // namespace brick


/* ================ Function definitions follow ================ */

#include <cmath>
#include <fstream>
#include <sstream>
#include <brick/common/exception.hh>

namespace brick {

  namespace utilities {

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
      writeRawPGM(std::string const& fileName,
                  const Type* imageData,
                  const size_t rows,
                  const size_t columns);

  
      template<class Type>
      void
      writeRawPPM(std::string const& fileName,
                  const Type* imageData,
                  const size_t rows,
                  const size_t columns);

  
      template<class Type>
      void
      writePlainPGM(std::string const& fileName,
                    const Type* imageData,
                    const size_t rows,
                    const size_t columns);

      
      template<class Type>
      void
      writePlainPPM(std::string const& fileName,
                    const Type* imageData,
                    const size_t rows,
                    const size_t columns);

    } // namespace privateCode
    /// @endcond
    

    template<class Type>
    void
    writePGM(std::string const& fileName,
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
    writePPM(std::string const& fileName,
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
      writeRawPGM(std::string const& fileName,
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
          BRICK_THROW(brick::common::IOException, "writeRawPGM()",
                      message.str().c_str());
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
            brick::common::switchByteOrder(
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
          BRICK_THROW(brick::common::ValueException, "writeRawPGM()",
                      message.str().c_str());
        }      
        
        // Check for errors.
        if(!outputStream) {
          outputStream.close();
          std::ostringstream message;
          message << "Error writing to open output file: " << fileName;
          BRICK_THROW(brick::common::IOException, "writeRawPGM()",
                      message.str().c_str());
        }
    
        // All done!
        outputStream.close();
      }

  
      template<class Type>
      void
      writeRawPPM(std::string const& fileName,
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
          BRICK_THROW(brick::common::IOException, "writeRawPPM()",
                      message.str().c_str());
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
            brick::common::switchByteOrder(
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
          BRICK_THROW(brick::common::ValueException, "writeRawPPM()",
                      message.str().c_str());
        }      
        
        // Check for errors.
        if(!outputStream) {
          std::ostringstream message;
          message << "Error writing to open output file: " << fileName;
          BRICK_THROW(brick::common::IOException, "writeRawPPM()",
                      message.str().c_str());
        }
    
        // All done!
        outputStream.close();
      }

  
      template<class Type>
      void
      writePlainPGM(std::string const& fileName,
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
          BRICK_THROW(brick::common::IOException, "writePlainPGM()",
                      message.str().c_str());
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
          BRICK_THROW(brick::common::IOException, "writePlainPGM()",
                      message.str().c_str());
        }

        // All done!
        outputStream.close();
      }
    

      template<class Type>
      void
      writePlainPPM(std::string const& fileName,
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
          BRICK_THROW(brick::common::IOException, "writePlainPPM()",
                      message.str().c_str());
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
          BRICK_THROW(brick::common::IOException, "writePlainPPM()",
                      message.str().c_str());
        }

        // All done!
        outputStream.close();
      }

    } // namespace privateCode
    /// @endcond

  } // namespace utilities
  
} // namespace brick

#endif /* #ifndef BRICK_UTILITIES_IMAGEIO_HH */
