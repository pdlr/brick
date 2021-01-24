/**
***************************************************************************
* @file brick/numeric/convolutionStrategy.hh
*
* Header file declaring an enum which represents different ways of
* handling border effects when performing convolution and correlation.
*
* Copyright (C) 2006-2007, 2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_CONVOLUTIONSTRATEGY_HH
#define BRICK_NUMERIC_CONVOLUTIONSTRATEGY_HH

namespace brick {

  namespace numeric {

    enum ConvolutionStrategy {
      BRICK_CONVOLVE_TRUNCATE_RESULT,
      BRICK_CONVOLVE_PAD_RESULT,
      BRICK_CONVOLVE_PAD_SIGNAL,
      BRICK_CONVOLVE_ZERO_PAD_SIGNAL,
      BRICK_CONVOLVE_REFLECT_SIGNAL,
      BRICK_CONVOLVE_WRAP_SIGNAL
    };


    enum ConvolutionROI {
      BRICK_CONVOLVE_ROI_SAME,
      BRICK_CONVOLVE_ROI_VALID,
      BRICK_CONVOLVE_ROI_FULL
    };

  } // namespace numeric

} // namespace brick

#endif /* #ifndef BRICK_NUMERIC_CONVOLUTIONSTRATEGY_HH */
