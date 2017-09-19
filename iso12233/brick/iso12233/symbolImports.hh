/**
***************************************************************************
* @file brick/iso12233/symbolImports.hh
*
* Header file that includes other brick library symbols in the
* iso12233 namespace to avoid lots of explicit namespace qualifiers.
*
* Copyright (C) 2017 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_ISO12233_SYMBOLIMPORTS_HH
#define BRICK_ISO12233_SYMBOLIMPORTS_HH

#include <brick/computerVision/image.hh>
#include <brick/computerVision/imageFormat.hh>
#include <brick/computerVision/imageFormatTraits.hh>
#include <brick/numeric/array1D.hh>
#include <brick/numeric/array2D.hh>

namespace brick {

  namespace iso12233 {

    using brick::computerVision::Image;
    using brick::computerVision::ImageFormat;
    using brick::computerVision::ImageFormatTraits;

    using brick::numeric::Array1D;
    using brick::numeric::Array2D;
    
  } // namespace iso12233
  
} // namespace brick

#endif /* #ifndef BRICK_ISO12233_SYMBOLIMPORTS_HH */
