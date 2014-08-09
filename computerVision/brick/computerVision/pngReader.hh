/**
***************************************************************************
* @file brick/computerVision/pngReader.hh
*
* Header file declaring a class for interfacing with libpng.
*
* Copyright (C) 2014 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_PNGREADER_HH

#ifndef HAVE_LIBPNG
#define HAVE_LIBPNG 1
#endif

#if HAVE_LIBPNG

#include <png.h>

#include <brick/computerVision/image.hh>

namespace brick {

  namespace computerVision {

    /**
     ** Wrapper class to make it easy to interact with libpng.
     **/
    class PngReader {
    public:

      /** 
       * The constructor opens a png image file and reads its contents.
       * 
       * @param fileName This argument is the name of the file to be
       * opened.
       */
      PngReader(std::string const& fileName);

        
      /**
       * Destructor.
       */
      virtual
      ~PngReader();


      /** 
       * Returns the contents of the the image file in the requested
       * image format.
       * 
       * @return The return value is an image reflecting the
       * contents of the .png file.
       */
      template <ImageFormat Format>
      Image<Format>
      getImage();


      /** 
       * Returns the native image format of the .png file.
       * 
       * @return The return value indicates how the data is stored
       * in the file.
       */
      ImageFormat
      getNativeImageFormat();
      

      /** 
       * Indicates whether the image data in the file is interlaced.
       * 
       * @return The return value is true if the .png file contains
       * interlaced data, false otherwise.
       */
      bool
      isInterlaced();

    private:

      // Helper function for getImage().
      template <ImageFormat Format>
      Image<Format>
      convertImage();


      // ---- Data members below this line ----

      png_structp m_pngPtr;
      png_infop m_infoPtr;
        
      png_uint_32 m_width;
      png_uint_32 m_height;
      int m_bitDepth;
      int m_colorType;
      int m_interlaceType;
      int m_compressionType;
      int m_filterMethod;
        
    };


    // Declare template specializations.
    template <>
    Image<GRAY8>
    PngReader::
    getImage<GRAY8>();


    template <>
    Image<GRAY16>
    PngReader::
    getImage<GRAY16>();
  

    template <>
    Image<RGB8>
    PngReader::
    getImage<RGB8>();


    template <>
    Image<RGB16>
    PngReader::
    getImage<RGB16>();

  } // namespace computerVision

} // namespace brick

#endif /* #if HAVE_LIBPNG */

#endif /* #ifndef BRICK_COMPUTERVISION_PNGREADER_HH */
