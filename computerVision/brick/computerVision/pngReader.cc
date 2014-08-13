/**
***************************************************************************
* @file brick/computerVision/pngReader.cc
*
* Source file defining a class for interfacing with libpng.
*
* Copyright (C) 2014 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#include <brick/common/byteOrder.hh>
#include <brick/computerVision/pngReader.hh>
#include <brick/computerVision/utilities.hh>

#if HAVE_LIBPNG

namespace brick {

  namespace computerVision {

    // The constructor opens a png image file and reads its contents.
    PngReader::
    PngReader(std::string const& fileName)
      : m_pngPtr(png_structp_NULL),
        m_infoPtr(0),
        m_width(0),
        m_height(0),
        m_bitDepth(0),
        m_colorType(0),
        m_interlaceType(0),
        m_compressionType(0),
        m_filterMethod(0)
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
        this->m_pngPtr = png_create_read_struct(
          PNG_LIBPNG_VER_STRING, 0, 0, 0);
        if(this->m_pngPtr == 0) {
          BRICK_THROW(brick::common::RunTimeException,
                      "dlr::computerVision::readPng8()",
                      "Couldn't initialize png_structp.");
        }

        // Be sure to clean up pngPtr (and later, infoPtr).
        try {

          // Allocate/initialize the memory for image information.
          this->m_infoPtr = png_create_info_struct(this->m_pngPtr);
          if (this->m_infoPtr == 0) {
            BRICK_THROW(brick::common::RunTimeException,
                        "brick::computerVision::readPng8()",
                        "Couldn't initialize png_infop.");
          }
       
          // Set error handling in case libpng calls longjmp().
          if(setjmp(png_jmpbuf(this->m_pngPtr))) {
            BRICK_THROW(brick::common::IOException,
                        "dlr::computerVision::readPng8()",
                        "Trouble reading from file.");
          }
        
          // Set up the input control.
          png_init_io(this->m_pngPtr, fp);

          // Let libpng know that we've already checked some magic.
          png_set_sig_bytes(this->m_pngPtr, pngSignatureSize);

#if 0
          // Seems like this should work, but it doesn't.
          
          // PNG files are natively big-endian, but can be coerced to
          // provide little-endian data without a swap.
          if(brick::common::getByteOrder()
             == brick::common::BRICK_LITTLE_ENDIAN) {
            png_set_swap(this->m_pngPtr);
          }
#endif /* #if 0 */
          
          // Read entire image into pngPtr.
          png_read_png(this->m_pngPtr, this->m_infoPtr,
                       PNG_TRANSFORM_IDENTITY, png_voidp_NULL);

          // Find out about our image.
          png_get_IHDR(this->m_pngPtr, this->m_infoPtr,
                       &(this->m_width), &(this->m_height),
                       &(this->m_bitDepth), &(this->m_colorType),
                       &(this->m_interlaceType), &(this->m_compressionType),
                       &(this->m_filterMethod));

          // Here we brute force the endianness fix that didn't work above.
          if(this->m_bitDepth == 16) {

            // TBD(xxx): Handle RGBA, etc.
            brick::common::UInt32 numberOfComponents = 1;
            if(this->m_colorType == PNG_COLOR_TYPE_RGB) {
              numberOfComponents = 3;
            }
            
            // Swap each row individually.
            png_bytep* rowPointers = png_get_rows(
              this->m_pngPtr, this->m_infoPtr);

            for(size_t rowIndex = 0; rowIndex < this->m_height; ++rowIndex) {
              brick::common::switchByteOrder(
                reinterpret_cast<common::UInt16*>(rowPointers[rowIndex]),
                this->m_width * numberOfComponents,
                brick::common::BRICK_BIG_ENDIAN,
                brick::common::getByteOrder());
            }

          } // if(this->m_bitDepth == 16)
            
        } catch(...) {
          png_destroy_read_struct(
            &(this->m_pngPtr), &(this->m_infoPtr), png_infopp_NULL);
          throw;
        }
        
      } catch(...) {
        fclose(fp);
        throw;
      }
      fclose(fp);
    }
    
    
    // Destructor.
    PngReader::
    ~PngReader()
    {
      png_destroy_read_struct(&(this->m_pngPtr), &(this->m_infoPtr),
                              png_infopp_NULL);
    }
    

    // Returns the contents of the the image file in the requested
    // image format.
    template <ImageFormat Format>
    Image<Format>
    PngReader::
    getImage()
    {
      BRICK_THROW(brick::common::NotImplementedException,
                  "PngReader::getImage()",
                  "Unsupported image format.");
    }


    template <>
    Image<GRAY8>
    PngReader::
    getImage<GRAY8>()
    {
      if(this->getNativeImageFormat() == GRAY8) {
        Image<GRAY8> result(this->m_height, this->m_width);
        png_bytep* rowPointers = png_get_rows(
          this->m_pngPtr, this->m_infoPtr);
        for(size_t rowIndex = 0; rowIndex < this->m_height; ++rowIndex) {
          result.getRow(rowIndex).copy(rowPointers[rowIndex]);
        }

        return result;
      }
    
      return this->convertImage<GRAY8>();
    }


    template <>
    Image<GRAY16>
    PngReader::
    getImage<GRAY16>()
    {
      if(this->getNativeImageFormat() == GRAY16) {
        Image<GRAY16> result(this->m_height, this->m_width);
        png_bytep* rowPointers = png_get_rows(
          this->m_pngPtr, this->m_infoPtr);
        for(size_t rowIndex = 0; rowIndex < this->m_height; ++rowIndex) {
          std::copy(rowPointers[rowIndex],
                    rowPointers[rowIndex] + this->m_width * 2, 
                    reinterpret_cast<common::UInt8*>(
                      result.getRow(rowIndex).data()));
        }

        return result;
      }

      return this->convertImage<GRAY16>();
    }
  

    template <>
    Image<RGB8>
    PngReader::
    getImage<RGB8>()
    {
      if(!PixelRGB8::isContiguous()) {
        BRICK_THROW(brick::common::NotImplementedException,
                    "PngReader::getImage<RGB8>()",
                    "This function currently only works with compilers that "
                    "don't add padding to the PixelRGB8 memory layout.");
      }
        
      if(this->getNativeImageFormat() == RGB8) {
        Image<RGB8> result(this->m_height, this->m_width);
        png_bytep* rowPointers = png_get_rows(
          this->m_pngPtr, this->m_infoPtr);
        for(size_t rowIndex = 0; rowIndex < this->m_height; ++rowIndex) {
          std::copy(rowPointers[rowIndex],
                    rowPointers[rowIndex] + this->m_width * 3, 
                    reinterpret_cast<common::UInt8*>(
                      result.getRow(rowIndex).data()));
        }

        return result;
      }

      return this->convertImage<RGB8>();
    }


    template <>
    Image<RGB16>
    PngReader::
    getImage<RGB16>()
    {
      if(!PixelRGB16::isContiguous()) {
        BRICK_THROW(brick::common::NotImplementedException,
                    "PngReader::getImage<RGB16>()",
                    "This function currently only works with compilers that "
                    "don't add padding to the PixelRGB16 memory layout.");
      }
        
      if(this->getNativeImageFormat() == RGB16) {
        Image<RGB16> result(this->m_height, this->m_width);
        png_bytep* rowPointers = png_get_rows(
          this->m_pngPtr, this->m_infoPtr);
        for(size_t rowIndex = 0; rowIndex < this->m_height; ++rowIndex) {
          std::copy(rowPointers[rowIndex],
                    rowPointers[rowIndex] + this->m_width * 6, 
                    reinterpret_cast<common::UInt8*>(
                      result.getRow(rowIndex).data()));
        }

        return result;
      }

      return this->convertImage<RGB16>();
    }
      

    // Returns the native image format of the .png file.
    ImageFormat
    PngReader::
    getNativeImageFormat()
    {
      if(this->m_bitDepth == 8) {
        if(this->m_colorType == PNG_COLOR_TYPE_GRAY) {
          return GRAY8;
        } else if(this->m_colorType == PNG_COLOR_TYPE_RGB) {
          return RGB8;
        }
      } else if(this->m_bitDepth == 16) {
        if(this->m_colorType == PNG_COLOR_TYPE_GRAY) {
          return GRAY16;
        } else if(this->m_colorType == PNG_COLOR_TYPE_RGB) {
          return RGB16;
        }
      }

      BRICK_THROW(brick::common::NotImplementedException,
                  "PngReader::getNativeImageFormat()",
                  "Unsupported image format.");
      return GRAY8;  // Keep the compiler happy.
    }
      

    // Indicates whether the image data in the file is interlaced.
    bool
    PngReader::
    isInterlaced()
    {
      return this->m_interlaceType != PNG_INTERLACE_NONE;
    }


    // ---- Private members below this line. ----
      
    // Helper function for getImage().
    template <ImageFormat Format>
    Image<Format>
    PngReader::
    convertImage()
    {
      Image<Format> result;
      switch(this->getNativeImageFormat()) {
      case GRAY8:
      {
        Image<GRAY8> tempImage = this->getImage<GRAY8>();
        result = convertColorspace<Format>(tempImage);
        break;
      }
      case GRAY16:
      {
        Image<GRAY16> tempImage = this->getImage<GRAY16>();
        result = convertColorspace<Format>(tempImage);
        break;
      }
      case RGB8:
      {
        Image<RGB8> tempImage = this->getImage<RGB8>();
        result = convertColorspace<Format>(tempImage);
        break;
      }
      case RGB16:
      {
        Image<RGB16> tempImage = this->getImage<RGB16>();
        result = convertColorspace<Format>(tempImage);
        break;
      }
      default:
        BRICK_THROW(brick::common::NotImplementedException,
                    "PngReader::convertImage()",
                    "Unsupported image format.");
        break;
      }

      return result;
    }
  
  } // namespace computerVision

} // namespace brick

#endif /* #if HAVE_LIBPNG */
