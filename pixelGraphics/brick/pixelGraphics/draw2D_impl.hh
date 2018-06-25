/**
***************************************************************************
* @file brick/pixelGraphics/draw2D_impl.hh
*
* Implementation file defining templates for drawing simple 2D
* graphics in user-defined arrays.
*
* Copyright (C) 2014 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_PIXELGRAPHICS_DRAW2D_IMPL_HH
#define BRICK_PIXELGRAPHICS_DRAW2D_IMPL_HH

// This file is included by draw2D.hh, and should not be
// directly included by user code, so no need to include
// draw2D.hh here.
//
// #include <brick/pixelGraphics/draw2D.hh>

#include <brick/numeric/amanatidesWoo2D.hh>

namespace brick {

  namespace pixelGraphics {

    template <class ArrayType, class CoordinateType>
    void
    draw2D(ArrayType& canvas,
           brick::geometry::LineSegment2D<CoordinateType> const& lineSegment,
           typename ArrayType::value_type const& color,
           brick::numeric::Transform2D<double> const& pixelFromWorld)
    {
      typedef typename AmanatidesWoo2D<ArrayType>::iterator awIterator;

      // Figure out in which direction to draw the line.
      brick::numeric::Vector2D<CoordinateType> direction =
        lineSegment.getVertex1() - lineSegment.getVertex0();

      // Use the fast voxel traversal algorithm of Amanatides & Woo to
      // draw the line.
      AmanatidesWoo2D<ArrayType> rayTracer(
        canvas, pixelFromWorld, lineSegment.getVertex0(), direction, true);
      awIterator iterator0 = rayTracer.begin();
      awIterator iterator1 = rayTracer.end();

      // This loop terminates if we iterate off the edge of the image,
      // and breaks (see the if() clause, below) when we reach the end
      // of the line segment.
      for(; iterator0 != iterator1; ++iterator0) {

        // Do the actual drawing.
        *iterator0 = color;

        // The AmanatidesWoo class is tracing the line defined by
        //
        //   pointOnLine = startPoint + t * direction
        //
        // where t is a free parameters.  Argument _direction_ is
        // defined above so that we reach the end of the line when
        // t == 1.0.  This makes the termination criterion easy.
        if(iterator0.tExit() >= 1.0) {
          break;
        }
      }
    }


  } // namespace pixelGraphics

} // namespace brick

#endif /* #ifndef BRICK_PIXELGRAPHICS_DRAW2D_IMPL_HH */
