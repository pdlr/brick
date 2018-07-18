/**
***************************************************************************
* @file brick/computerVision/imageFormat.hh
*
* Header file declaring ImageFormat enum.
*
* Copyright (C) 2005-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_IMAGEFORMAT_HH
#define BRICK_COMPUTERVISION_IMAGEFORMAT_HH

namespace brick {

  namespace computerVision {

    /**
     ** This enum indicates the acceptable image format values.  Please
     ** see the ImageFormatTraits class template for the characteristics
     ** of these image formats.
     **/
    enum ImageFormat {
      GRAY1,
      GRAY8,
      GRAY16,
      GRAY32,
      GRAY64,
      GRAY_SIGNED16,
      GRAY_SIGNED32,
      GRAY_FLOAT32,
      GRAY_FLOAT64,
      RGB8,
      RGB16,
      RGB_SIGNED16,
      RGB_SIGNED32,
      RGB_FLOAT32,
      RGB_FLOAT64,
      HSV_FLOAT64,
      YIQ_FLOAT64,
      BGRA8,
      RGBA8,
      YUV420,
      NO_FORMAT
    };

  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_IMAGEFORMAT_HH */
