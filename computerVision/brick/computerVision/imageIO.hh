/**
***************************************************************************
* @file brick/computerVision/imageIO.hh
*
* Header file declaring functions for reading and writing images.
*
* Copyright (C) 2005-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_IMAGEIO_HH
#define BRICK_COMPUTERVISION_IMAGEIO_HH

#include <string>
#include <brick/computerVision/image.hh>

namespace brick {

  namespace computerVision {

#if 0  /* Here's the interface I think we ultimately want. */
    
    class ImageReader {
    public:
      
      /**
       * Constructor.
       */
      ImageReader();
      
      
      /**
       * Destructor.
       */
      virtual
      ~ImageReader();


      void
      readFile(std::string const& fileName);


      ImageFormat
      getNativeFormat();


      template<ImageFormat format>
      Image<format>
      getImage() {return Image<format>();}
      
    private:

      bool
      readPGM(std::string const& fileName);

      bool
      readPNG(std::string const& fileName);

      bool
      readPPM(std::string const& fileName);
      
    };

#endif

    /* ================ Non-member function declarations. ================ */
  
    Image<GRAY8>
    readPGM8(const std::string& fileName);

  
    Image<GRAY8>
    readPGM8(const std::string& fileName,
             std::string& commentString);

  
    Image<GRAY16>
    readPGM16(const std::string& fileName);

  
    Image<RGB8>
    readPPM8(const std::string& fileName);


    Image<RGB8>
    readPPM8(const std::string& fileName,
             std::string& commentString);
    
  
    void
    writePGM8(const std::string& fileName,
              const Image<GRAY8>& outputImage,
              const std::string& comment = "");
  
  
    void
    writePGM16(const std::string& fileName,
               const Image<GRAY16>& outputImage,
               const std::string& comment = "");
  

    void
    writePPM8(const std::string& fileName,
              const Image<RGB8>& outputImage,
              const std::string& comment = "");


    void
    writePPM16(const std::string& fileName,
               const Image<RGB16>& outputImage,
               const std::string& comment = "");
    
    
#ifndef HAVE_LIBPNG
#define HAVE_LIBPNG 1
#endif
    
#if HAVE_LIBPNG

    /** 
     * WARNING: This routine may not stick around for long.
     * 
     * @param fileName This argument ...
     * 
     * @param commentString This argument is currently ignored.
     * 
     * @return The return value ...
     */
    Image<GRAY8>
    readPNG8(const std::string& fileName,
             std::string& commentString);


    /** 
     * WARNING: This routine may not stick around for long.
     * 
     * @param fileName This argument ...
     * 
     * @param outputImage This argument ...
     * 
     * @param comment This argument is currently ignored.
     */
    void
    writePNG8(const std::string& fileName,
              const Image<GRAY8>& outputImage,
              const std::string& comment = "");

#endif /* #if HAVE_LIBPNG */



    /* ========== Lowest common denominator image IO functions =========== */

    template<class Type>
    void
    writePGM(const std::string& fileName,
             const Type* imageData,
             const size_t rows,
             const size_t columns,
             bool normalize=false,
             bool rawFormat=true,
             int bitsPerPixel=8);

  
    template<class Type>
    void
    writePPM(const std::string& fileName,
             const Type* imageData,
             const size_t rows,
             const size_t columns,
             bool normalize=false,
             bool rawFormat=true,
             int bitsPerPixel=8);
    
      
  } // namespace computerVision

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/imageIO_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_IMAGEIO_HH */
