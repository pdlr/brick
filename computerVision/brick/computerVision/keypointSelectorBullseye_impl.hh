/**
***************************************************************************
* @file brick/computerVision/keypointSelectorBullseye_impl.hh
*
* Header file defining a class template for selecting stable keypoints
* from an image.
*
* Copyright (C) 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_KEYPOINTSELECTORBULLSEYE_IMPL_HH
#define BRICK_COMPUTERVISION_KEYPOINTSELECTORBULLSEYE_IMPL_HH

// This file is included by keypointSelectorBullseye.hh, and should not be directly included
// by user code, so no need to include keypointSelectorBullseye.hh here.
// 
// #include <brick/computerVision/keypointSelectorBullseye.hh>

namespace brick {

  namespace computerVision {

    // Currently empty, but code will go here as we optimize.

    // ============== Private member functions below this line ==============


    accumulateAsymmetrySums(Int32 pixel0, Int32 pixel1,
                            UnsignedInt32& pixelSum,
                            UnsignedInt32& pixelSquaredSum,
                            UnsignedInt32& asymmetrySum)
    {
      pixelSum += pixel0 + pixel1;
      pixelSquaredSum += pixel0 * pixel0 + pixel1 * pixel1;
      asymmetrySum += common::absoluteValue(pixel0 - pixel1);
    }
    
  } // namespace computerVision
  
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_KEYPOINTSELECTORBULLSEYE_IMPL_HH */
