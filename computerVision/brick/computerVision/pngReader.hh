/**
***************************************************************************
* @file brick/computerVision/pngInterface.hh
*
* Header file declaring a class for interfacing with libpng.
*
* Copyright (C) 2014 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_PNGINTERFACE_HH

#ifndef HAVE_LIBPNG
#define HAVE_LIBPNG 1
#endif

#if HAVE_LIBPNG

#include <png.h>

#include <brick/computerVision/image.hh>

namespace brick {

  namespace computerVision {

    /// @cond privateCode
    namespace privateCode {

      /**
       ** Wrapper class to make it easy to interact with libpng.
       **/
      class PngInterface {
      public:

        /** 
         * The constructor opens a png image file and reads its contents.
         * 
         * @param fileName This argument is the name of the file to be
         * opened.
         */
        PngInterface(std::string const& fileName);

        
        /**
         * Destructor.
         */
        virtual
        ~PngInterface();


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
      PngInterface::
      getImage<GRAY8>();


      template <>
      Image<GRAY16>
      PngInterface::
      getImage<GRAY16>();
  

      template <>
      Image<RGB8>
      PngInterface::
      getImage<RGB8>();


      template <>
      Image<RGB16>
      PngInterface::
      getImage<RGB16>();

      
    } // namespace privateCode

  } // namespace computerVision

} // namespace brick

#endif /* #if HAVE_LIBPNG */

#endif /* #ifndef BRICK_COMPUTERVISION_PNGINTERFACE_HH */
