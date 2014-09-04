/**
***************************************************************************
* @file brick/pixelGraphics/draw2D.hh
*
* Header file declaring routines for drawing simple 2D graphics in
* user-defined arrays.
*
* Copyright (C) 2014 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_PIXELGRAPHICS_DRAW2D_HH
#define BRICK_PIXELGRAPHICS_DRAW2D_HH

#include <brick/geometry/lineSegment2D.hh>

namespace brick {

  namespace pixelGraphics {

    template <class ArrayType, class CoordinateType>
    void
    draw2D(ArrayType& canvas,
           brick::geometry::LineSegment2D<CoordinateType> const& lineSegment,
           typename ArrayType::value_type const& color,
           brick::numeric::Transform2D<double> const& pixelFromWorld =
             brick::numeric::Transform2D<double>());
    
  } // namespace pixelGraphics

} // namespace brick

#include <brick/pixelGraphics/draw2D_impl.hh>

#endif /* #ifdef BRICK_PIXELGRAPHICS_DRAW2D_HH */
