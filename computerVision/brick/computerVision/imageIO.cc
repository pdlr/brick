/**
***************************************************************************
* @file brick/computerVision/imageIO.cc
*
* Source file defining functions for reading and writing images.
*
* Copyright (C) 2005, 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#include <fstream>

#include <brick/common/byteOrder.hh>
#include <brick/computerVision/imageIO.hh>

namespace {

  void
  discardComments(std::istream& inputStream,
                  std::istream::char_type commentCharacter)
  {
    std::istream::char_type inputCharacter;
    while(1) {
      inputStream >> inputCharacter;
      if(inputCharacter == commentCharacter) {
        std::string dummyString;
        getline(inputStream, dummyString);
      } else {
        inputStream.putback(inputCharacter);
        break;
      }
    }
  }


  std::string
  readComments(std::istream& inputStream,
               std::istream::char_type commentCharacter)
  {
    std::istream::char_type inputCharacter;
    std::ostringstream commentStream;
    while(1) {
      inputStream >> inputCharacter;
      if(inputCharacter == commentCharacter) {
        std::string dummyString;
        getline(inputStream, dummyString);
        commentStream << dummyString;
      } else {
        inputStream.putback(inputCharacter);
        break;
      }
    }
    return commentStream.str();
  }

} // Anonymous namespace


namespace brick {

  namespace computerVision {


#if 0
    
    // Constructor.
    ImageReader::
    ImageReader()
    {
      // Empty.
    }
      
      
    // Destructor.
    ImageReader::
    ~ImageReader()
    {
      // Empty.
    }


    void
    ImageReader::
    readFile(std::string const& fileName)
    {

    }


    ImageFormat
    ImageReader::
    getNativeFormat()
    {
      return GRAY8;
    }

#endif
    
    Image<GRAY8>
    readPGM8(const std::string& fileName)
    {
      std::string commentString;
      return readPGM8(fileName, commentString);
    }

  
    Image<GRAY8>
    readPGM8(const std::string& fileName, std::string& commentString)
    {
      std::ifstream inputStream(fileName.c_str(), std::ios::binary);
      if(!inputStream) {
        std::ostringstream message;
        message << "Couldn't open input file: " << fileName;
        BRICK_THROW(brick::common::IOException,
                    "readPGM8()", message.str().c_str());
      }

      // Read the header.
      commentString.clear();
      std::string magic;
      size_t columns;
      size_t rows;
      long long int imageMax;
      inputStream >> magic;
      commentString += readComments(inputStream, '#');
      inputStream >> columns >> rows;
      commentString += readComments(inputStream, '#');
      inputStream >> imageMax;

      // Image data starts after the next newline.
      std::string dummy;
      std::getline(inputStream, dummy);

      // Check for errors.
      if(!inputStream) {
        std::ostringstream message;
        message << "Couldn't read image header from file: " << fileName;
        BRICK_THROW(brick::common::IOException,
                    "readPGM8()", message.str().c_str());
      }

      // Check that this image will fit in an 8-bit pixel array.
      if(imageMax > 255LL) {
        std::ostringstream message;
        message << "File " << fileName << " has max value of " << imageMax
                << ", which is too big for an 8-bit image.";
        BRICK_THROW(brick::common::IOException,
                    "readPGM8()", message.str().c_str());
      }

      // Allocate storage space.
      Image<GRAY8> newImage(rows, columns);

      // And read the pixels.
      if(magic == "P5") {
        // Looks like a raw image.
        inputStream.read(reinterpret_cast<char*>(newImage.data()),
                         newImage.size());
      } else if(magic == "P2") {
        // Looks like a plain image.
        unsigned short interpreter;
        for(size_t pixelIndex = 0; pixelIndex < newImage.size(); ++pixelIndex) {
          // We can't simply read into the image, since the compiler
          // will assume we want the ascii values of the characters in
          // the file.
          inputStream >> interpreter;
          newImage(pixelIndex) =
            static_cast<Image<GRAY8>::PixelType>(interpreter);
        }
      } else {
        std::ostringstream message;
        message << "Incorrect magic, " << magic << " in file " << fileName
                << ".";
        BRICK_THROW(brick::common::IOException,
                    "readPGM8()", message.str().c_str());
      }

      // Check for errors.
      if(!inputStream) {
        std::ostringstream message;
        message << "Error reading image data from input file: " << fileName;
        BRICK_THROW(brick::common::IOException,
                    "readPGM8()", message.str().c_str());
      }
    
      // All done!
      inputStream.close();
      return newImage;
    }


    Image<GRAY16>
    readPGM16(const std::string& fileName)
    {
      std::ifstream inputStream(fileName.c_str(), std::ios::binary);
      if(!inputStream) {
        std::ostringstream message;
        message << "Couldn't open input file: " << fileName;
        BRICK_THROW(brick::common::IOException,
                    "readPGM16()", message.str().c_str());
      }

      // Read the header.
      std::string magic;
      size_t columns;
      size_t rows;
      long long int imageMax;
      inputStream >> magic;
      discardComments(inputStream, '#');
      inputStream >> columns >> rows;
      discardComments(inputStream, '#');
      inputStream >> imageMax;

      // Image data starts after the next newline.
      std::string dummy;
      std::getline(inputStream, dummy);

      // Check for errors.
      if(!inputStream) {
        std::ostringstream message;
        message << "Couldn't read image header from file: " << fileName;
        BRICK_THROW(brick::common::IOException,
                    "readPGM16()", message.str().c_str());
      }

      // Check that this image will fit in an 16-bit pixel array.
      if(imageMax > 65535LL) {
        std::ostringstream message;
        message << "File " << fileName << " has max value of " << imageMax
                << ", which is too big for an 16-bit image.";
        BRICK_THROW(brick::common::IOException,
                    "readPGM16()", message.str().c_str());
      }

      // Allocate storage space.
      Image<GRAY16> newImage(rows, columns);

      // And read the pixels.
      if(magic == "P5") {
        // Looks like a raw image.
        if(imageMax <= 255LL) {
          // Looks like an 8 bit image.
          BRICK_THROW(brick::common::NotImplementedException, "readPGM16()",
                      "This routine currently cannot read 8-bit PGMs");
        } else {
          size_t numberOfBytes = 2 * newImage.size();
          inputStream.read(
            reinterpret_cast<char*>(newImage.data()), numberOfBytes);
          switchByteOrder(newImage.data(), newImage.size(),
                          brick::common::BRICK_BIG_ENDIAN,
                          brick::common::getByteOrder());
        }
      } else if(magic == "P2") {
        // Looks like a plain image.
        unsigned short interpreter;
        for(size_t pixelIndex = 0; pixelIndex < newImage.size(); ++pixelIndex) {
          // We can't simply read into the image, since the compiler
          // will assume we want the ascii values of the characters in
          // the file.
          inputStream >> interpreter;
          newImage(pixelIndex) =
            static_cast<Image<GRAY16>::PixelType>(interpreter);
        }
      } else {
        std::ostringstream message;
        message << "Incorrect magic, " << magic << " in file " << fileName
                << ".";
        BRICK_THROW(brick::common::IOException,
                    "readPGM16()", message.str().c_str());
      }

      // Check for errors.
      if(!inputStream) {
        std::ostringstream message;
        message << "Error reading image data from input file: " << fileName;
        BRICK_THROW(brick::common::IOException,
                    "readPGM16()", message.str().c_str());
      }
    
      // All done!
      inputStream.close();
      return newImage;
    }


    Image<RGB8>
    readPPM8(const std::string& fileName)
    {
      std::string commentString;
      return readPPM8(fileName, commentString);
    }
      

    Image<RGB8>
    readPPM8(const std::string& fileName,
             std::string& commentString)
    {
      std::ifstream inputStream(fileName.c_str(), std::ios::binary);
      if(!inputStream) {
        std::ostringstream message;
        message << "Couldn't open input file: " << fileName;
        BRICK_THROW(brick::common::IOException,
                    "readPPM8()", message.str().c_str());
      }

      // Read the header.
      commentString.clear();
      std::string magic;
      size_t columns;
      size_t rows;
      long long int imageMax;
      inputStream >> magic;
      commentString += readComments(inputStream, '#');
      inputStream >> columns >> rows;
      commentString += readComments(inputStream, '#');
      inputStream >> imageMax;

      // Image data starts after the next newline.
      std::string dummy;
      std::getline(inputStream, dummy);

      // Check for errors.
      if(!inputStream) {
        std::ostringstream message;
        message << "Couldn't read image header from file: " << fileName;
        BRICK_THROW(brick::common::IOException,
                    "readPPM8()", message.str().c_str());
      }

      // Check that this image will fit in an 8-bit pixel array.
      if(imageMax > 255LL) {
        std::ostringstream message;
        message << "File " << fileName << " has max value of " << imageMax
                << ", which is too big for an 8-bit image.";
        BRICK_THROW(brick::common::IOException,
                    "readPPM8()", message.str().c_str());
      }

      // Allocate storage space.
      Image<RGB8> newImage(rows, columns);

      // And read the pixels.
      if(magic == "P6") {
        // Looks like a raw image.
        inputStream.read(reinterpret_cast<char*>(newImage.data()),
                         newImage.size() * 3);
      } else if(magic == "P3") {
        // Looks like a plain image.
        unsigned short interpreterRed;
        unsigned short interpreterGreen;
        unsigned short interpreterBlue;
        for(size_t pixelIndex = 0; pixelIndex < newImage.size(); ++pixelIndex) {
          // We can't simply read into the image, since the compiler
          // will assume we want the ascii values of the characters in
          // the file.
          inputStream >> interpreterRed >> interpreterGreen >> interpreterBlue;
          newImage(pixelIndex).red =
            static_cast<brick::common::UnsignedInt8>(interpreterRed);
          newImage(pixelIndex).green =
            static_cast<brick::common::UnsignedInt8>(interpreterGreen);
          newImage(pixelIndex).blue =
            static_cast<brick::common::UnsignedInt8>(interpreterBlue);
        }
      } else {
        std::ostringstream message;
        message << "Incorrect magic, " << magic << " in file " << fileName
                << ".";
        BRICK_THROW(brick::common::IOException,
                    "readPPM8()", message.str().c_str());
      }

      // Check for errors.
      if(!inputStream) {
        std::ostringstream message;
        message << "Error reading image data from input file: " << fileName;
        BRICK_THROW(brick::common::IOException,
                    "readPPM8()", message.str().c_str());
      }
    
      // All done!
      inputStream.close();
      return newImage;
    }


    void
    writePGM8(const std::string& fileName,
              const Image<GRAY8>& outputImage,
              const std::string& comment)
    {
      std::ofstream outputStream(fileName.c_str(), std::ios::binary);
      if(!outputStream) {
        std::ostringstream message;
        message << "Couldn't open output file: " << fileName;
        BRICK_THROW(brick::common::IOException,
                    "writePGM8()", message.str().c_str());
      }

      // Write the header.
      outputStream << "P5\n";
      if(comment != "") {
        outputStream << "# " << comment << "\n";
      }
      outputStream << outputImage.columns() << " " << outputImage.rows() << "\n"
                   << "255\n";

      // Check for errors.
      if(!outputStream) {
        std::ostringstream message;
        message << "Couldn't write image header to file: " << fileName;
        BRICK_THROW(brick::common::IOException,
                    "writePGM8()", message.str().c_str());
      }

      // And write the pixels.
      outputStream.write(reinterpret_cast<const char*>(outputImage.data()),
                         outputImage.size());

      // Check for errors.
      if(!outputStream) {
        std::ostringstream message;
        message << "Error writing image data to output file: " << fileName;
        BRICK_THROW(brick::common::IOException,
                    "writePGM8()", message.str().c_str());
      }
    
      // All done!
      outputStream.close();
    }
  

    void
    writePGM16(const std::string& fileName,
               const Image<GRAY16>& outputImage,
               const std::string& comment)
    {
      std::ofstream outputStream(fileName.c_str(), std::ios::binary);
      if(!outputStream) {
        std::ostringstream message;
        message << "Couldn't open output file: " << fileName;
        BRICK_THROW(brick::common::IOException,
                    "writePGM16()", message.str().c_str());
      }

      // Write the header.
      outputStream << "P5\n";
      if(comment != "") {
        outputStream << "# " << comment << "\n";
      }
      outputStream << outputImage.columns() << " " << outputImage.rows() << "\n"
                   << "65535\n";

      // Check for errors.
      if(!outputStream) {
        std::ostringstream message;
        message << "Couldn't write image header to file: " << fileName;
        BRICK_THROW(brick::common::IOException,
                    "writePGM16()", message.str().c_str());
      }

      // And write the pixels.
      size_t numberOfElements = outputImage.size();
      size_t numberOfBytes =
        numberOfElements * sizeof(Image<GRAY16>::value_type);
      if(brick::common::getByteOrder() == brick::common::BRICK_BIG_ENDIAN) {
        outputStream.write(
          reinterpret_cast<const char*>(outputImage.data()), numberOfBytes);
      } else {
        char* swabbedData = new char[numberOfBytes];
        switchByteOrder(
          outputImage.data(), numberOfElements,
          reinterpret_cast<Image<GRAY16>::value_type*>(swabbedData),
          brick::common::getByteOrder(), brick::common::BRICK_BIG_ENDIAN);
        outputStream.write(swabbedData, numberOfBytes);
        delete[] swabbedData;
      }

      // Check for errors.
      if(!outputStream) {
        outputStream.close();
        std::ostringstream message;
        message << "Error writing image data to output file: " << fileName;
        BRICK_THROW(brick::common::IOException,
                    "writePGM16()", message.str().c_str());
      }
    
      // All done!
      outputStream.close();
    }
  

    void
    writePPM8(const std::string& fileName,
              const Image<RGB8>& outputImage,
              const std::string& comment)
    {
      if(!Image<RGB8>::value_type::isContiguous()) {
        BRICK_THROW(
          brick::common::NotImplementedException, "writePPM8()",
          "Your compiler appears to pack RGB pixels in an unusual way.  "
          "Please update brick::computerVision::writePPM8() to handle this.");
      }
      
      std::ofstream outputStream(fileName.c_str(), std::ios::binary);
      if(!outputStream) {
        std::ostringstream message;
        message << "Couldn't open output file: " << fileName;
        BRICK_THROW(brick::common::IOException,
                    "writePPM8()", message.str().c_str());
      }

      // Write the header.
      outputStream << "P6\n";
      if(comment != "") {
        outputStream << "# " << comment << "\n";
      }
      outputStream << outputImage.columns() << " " << outputImage.rows() << "\n"
                   << "255\n";

      // Check for errors.
      if(!outputStream) {
        std::ostringstream message;
        message << "Couldn't write image header to file: " << fileName;
        BRICK_THROW(brick::common::IOException,
                    "writePPM8()", message.str().c_str());
      }

      // And write the pixels.
      outputStream.write(reinterpret_cast<const char*>(outputImage.data()),
                         3 * outputImage.size());

      // Check for errors.
      if(!outputStream) {
        std::ostringstream message;
        message << "Error writing image data to output file: " << fileName;
        BRICK_THROW(brick::common::IOException,
                    "writePPM8()", message.str().c_str());
      }
    
      // All done!
      outputStream.close();
    }


    void
    writePPM16(const std::string& fileName,
               const Image<RGB16>& outputImage,
               const std::string& comment)
    {
      if(!Image<RGB16>::value_type::isContiguous()) {
        BRICK_THROW(
          brick::common::NotImplementedException, "writePPM16()",
          "Your compiler appears to pack RGB pixels in an unusual way.  "
          "Please update brick::computerVision::writePPM16() to handle this.");
      }
      
      std::ofstream outputStream(fileName.c_str(), std::ios::binary);
      if(!outputStream) {
        std::ostringstream message;
        message << "Couldn't open output file: " << fileName;
        BRICK_THROW(brick::common::IOException,
                    "writePPM16()", message.str().c_str());
      }

      // Write the header.
      outputStream << "P6\n";
      if(comment != "") {
        outputStream << "# " << comment << "\n";
      }
      outputStream << outputImage.columns() << " " << outputImage.rows() << "\n"
                   << "65535\n";

      // Check for errors.
      if(!outputStream) {
        std::ostringstream message;
        message << "Couldn't write image header to file: " << fileName;
        BRICK_THROW(brick::common::IOException,
                    "writePPM16()", message.str().c_str());
      }

      // And write the pixels.
      size_t numberOfElements = outputImage.size();
      size_t numberOfBytes =
        numberOfElements * sizeof(Image<RGB16>::value_type);
      if(brick::common::getByteOrder() == brick::common::BRICK_BIG_ENDIAN) {
        outputStream.write(
          reinterpret_cast<const char*>(outputImage.data()), numberOfBytes);
      } else {
        char* swabbedData = new char[numberOfBytes];
        switchByteOrder(
          reinterpret_cast<brick::common::UInt16 const*>(outputImage.data()),
          numberOfBytes / 2,
          reinterpret_cast<brick::common::UInt16*>(swabbedData),
          brick::common::getByteOrder(), brick::common::BRICK_BIG_ENDIAN);
        outputStream.write(swabbedData, numberOfBytes);
        delete[] swabbedData;
      }

      // Check for errors.
      if(!outputStream) {
        outputStream.close();
        std::ostringstream message;
        message << "Error writing image data to output file: " << fileName;
        BRICK_THROW(brick::common::IOException,
                    "writePPM16()", message.str().c_str());
      }
    
      // All done!
      outputStream.close();
    }

  } // namespace computerVision

} // namespace brick

    
#if HAVE_LIBPNG

#include <png.h>

namespace brick {

  namespace computerVision {
    
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


    void
    writePNG8(const std::string& fileName,
              const Image<GRAY8>& outputImage,
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
        BRICK_THROW(brick::common::IOException, "ImageIO::writePng()",
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

          // Prepare for writing.
          png_set_IHDR(
            pngPtr, infoPtr, outputImage.columns(), outputImage.rows(),
            sizeof(ImageFormatTraits<GRAY8>::PixelType) * 8,
            PNG_COLOR_TYPE_GRAY, PNG_INTERLACE_NONE,
            PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
         
     
          // xxx Fill in png structure here.
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
              rowPointers[rowIndex] = const_cast<unsigned char*>(
                outputImage.getRow(rowIndex).data());
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


    void
    writePNG8(const std::string& fileName,
              const Image<RGB8>& outputImage,
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
        BRICK_THROW(brick::common::IOException, "ImageIO::writePng()",
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

          // Prepare for writing.
          png_set_IHDR(
            pngPtr, infoPtr, outputImage.columns(), outputImage.rows(),
            sizeof(brick::computerVision::ImageFormatTraits<brick::computerVision::GRAY8>::PixelType) * 8,
            PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
            PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);           

          // xxx Fill in png structure here.
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
              rowPointers[rowIndex] = (png_bytep) const_cast<brick::computerVision::PixelRGB<unsigned char>*>(
                outputImage.getRow(rowIndex).data());
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

      // End writing PNG in color
    }

  } // namespace computerVision

} // namespace brick

#endif /* #if HAVE_LIBPNG */
