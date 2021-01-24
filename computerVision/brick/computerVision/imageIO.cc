/**
***************************************************************************
* @file brick/computerVision/imageIO.cc
*
* Source file defining functions for reading and writing images.
*
* Copyright (C) 2005, 2011 David LaRose, dlr@davidlarose.com
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

namespace brick {

  namespace computerVision {

    template<>
    void
    writePNG<GRAY8>(const std::string& fileName,
                    const Image<GRAY8>& outputImage,
                    const std::string& /* comment */)
    {
      try {
        png::image<png::gray_pixel> pngImage(outputImage.columns(),
                                          outputImage.rows());
        for (png::uint_32 yy = 0; yy < pngImage.get_height(); ++yy) {
          for (png::uint_32 xx = 0; xx < pngImage.get_width(); ++xx) {
            common::UInt8 const& pixel = outputImage(yy, xx);
            pngImage[yy][xx] = png::gray_pixel(pixel);
          }
        }
        pngImage.write(fileName);
      } catch(...) {
        std::ostringstream message;
        message << "Failed to write png image to " << fileName;
        BRICK_THROW(brick::common::IOException, "writePNG()",
                    message.str().c_str());
      }
    }

    template<>
    void
    writePNG<GRAY16>(const std::string& fileName,
                     const Image<GRAY16>& outputImage,
                     const std::string& /* comment */)
    {
      try {
        png::image<png::gray_pixel_16> pngImage(outputImage.columns(),
                                              outputImage.rows());
        for (png::uint_32 yy = 0; yy < pngImage.get_height(); ++yy) {
          for (png::uint_32 xx = 0; xx < pngImage.get_width(); ++xx) {
            common::UInt16 const& pixel = outputImage(yy, xx);
            pngImage[yy][xx] = png::gray_pixel_16(pixel);
          }
        }
        pngImage.write(fileName);
      } catch(...) {
        std::ostringstream message;
        message << "Failed to write png image to " << fileName;
        BRICK_THROW(brick::common::IOException, "writePNG()",
                    message.str().c_str());
      }
    }

    template<>
    void
    writePNG<RGB8>(const std::string& fileName,
                   const Image<RGB8>& outputImage,
                   const std::string& /* comment */)
    {
      try {
        png::image<png::rgb_pixel> pngImage(outputImage.columns(),
                                            outputImage.rows());
        for (png::uint_32 yy = 0; yy < pngImage.get_height(); ++yy) {
          for (png::uint_32 xx = 0; xx < pngImage.get_width(); ++xx) {
            PixelRGB8 const& pixel = outputImage(yy, xx);
            pngImage[yy][xx] = png::rgb_pixel(
              pixel.red, pixel.green, pixel.blue);
          }
        }
        pngImage.write(fileName);
      } catch(...) {
        std::ostringstream message;
        message << "Failed to write png image to " << fileName;
        BRICK_THROW(brick::common::IOException, "writePNG()",
                    message.str().c_str());
      }
    }

    template<>
    void
    writePNG<RGB16>(const std::string& fileName,
                    const Image<RGB16>& outputImage,
                    const std::string& /* comment */)
    {
      try {
        png::image<png::rgb_pixel_16> pngImage(outputImage.columns(),
                                               outputImage.rows());
        for (png::uint_32 yy = 0; yy < pngImage.get_height(); ++yy) {
          for (png::uint_32 xx = 0; xx < pngImage.get_width(); ++xx) {
            PixelRGB16 const& pixel = outputImage(yy, xx);
            pngImage[yy][xx] = png::rgb_pixel_16(
              pixel.red, pixel.green, pixel.blue);
          }
        }
        pngImage.write(fileName);
      } catch(...) {
        std::ostringstream message;
        message << "Failed to write png image to " << fileName;
        BRICK_THROW(brick::common::IOException, "writePNG()",
                    message.str().c_str());
      }
    }

  } // namespace computerVision

} // namespace brick


# endif /* #if HAVE_LIBPNG */
