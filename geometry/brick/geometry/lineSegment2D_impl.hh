/**
***************************************************************************
* @file brick/geometry/lineSegment2D_impl.hh
*
* Source file defining the LineSegment2D class.
*
* Copyright (C) 2007-2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_GEOMETRY_LINESEGMENT2D_IMPL_HH
#define BRICK_GEOMETRY_LINESEGMENT2D_IMPL_HH

// This file is included by lineSegment2D.hh, and should not be
// directly included by user code, so no need to include
// lineSegment2D.hh here.
//
// #include <brick/geometry/lineSegment2D.hh>

#include <brick/common/expect.hh>
// #include <brick/numeric/utilities.hh>

namespace brick {

  namespace geometry {

    // The default constructor initializes to the line segment that starts
    // at the origin and ends at (1, 0).
    template <class Type>
    LineSegment2D<Type>::
    LineSegment2D()
	: m_startPoint(0.0, 0.0), m_endPoint(1.0, 0.0)
    {
      // Empty.
    }


    // This constructor initializes the line segment using a pair of points.
    template <class Type>
    LineSegment2D<Type>::
    LineSegment2D(brick::numeric::Vector2D<Type> const& startPoint,
                  brick::numeric::Vector2D<Type> const& endPoint)
      : m_startPoint(startPoint), m_endPoint(endPoint)
    {
      // Empty.
    }


    // This constructor initializes the line segment using four
    // coordinates representing a pair of points.
    template <class Type>
    LineSegment2D<Type>::
    LineSegment2D(Type const& x0, Type const& y0,
                  Type const& x1, Type const& y1)
      : m_startPoint(x0, y0), m_endPoint(x1, y1)
    {
      // Empty.
    }


    // The copy constructor deep copies its argument.
    template <class Type>
    LineSegment2D<Type>::
    LineSegment2D(LineSegment2D<Type> const& source)
      : m_startPoint(source.m_startPoint), m_endPoint(source.m_endPoint)
    {
      // Empty.
    }


    // Destructor.
    template <class Type>
    LineSegment2D<Type>::
    ~LineSegment2D()
    {
      // Empty.
    }


    // The assignment operator deep copies its argument.
    template <class Type>
    LineSegment2D<Type>&
    LineSegment2D<Type>::
    operator=(const LineSegment2D<Type>& source)
    {
      if(&source != this) {
        m_startPoint = source.m_startPoint;
        m_endPoint = source.m_endPoint;
      }
      return *this;
    }


    // This member function returns the end point of the line segment.
    template <class Type>
    brick::numeric::Vector2D<Type> const&
    LineSegment2D<Type>::
    getEndPoint() const
    {
      return m_endPoint;
    }


    // This member function returns the start point of the line segment.
    template <class Type>
    brick::numeric::Vector2D<Type> const&
    LineSegment2D<Type>::
    getStartPoint() const {
      return m_startPoint;
    }


    // This member function returns the start point of the line segment.
    template <class Type>
    brick::numeric::Vector2D<Type> const&
    LineSegment2D<Type>::
    getVertex0() const
    {
      return m_startPoint;
    }


    // This member function returns the end point of the line segment.
    template <class Type>
    brick::numeric::Vector2D<Type> const&
    LineSegment2D<Type>::
    getVertex1() const
    {
      return m_endPoint;
    }


    // This member function changes the start point and end point of
    // the line segment.
    template <class Type>
    LineSegment2D<Type>&
    LineSegment2D<Type>::
    setValue(Type const& startPointX, Type const& startPointY,
             Type const& endPointX, Type const& endPointY)
    {
      m_startPoint.setValue(startPointX, startPointY);
      m_endPoint.setValue(endPointX, endPointY);
      return *this;
    }


    // This member function changes the start point and end point of
    // the line segment.
    template <class Type>
    LineSegment2D<Type>&
    LineSegment2D<Type>::
    setValue(brick::numeric::Vector2D<Type> const& startPoint,
             brick::numeric::Vector2D<Type> const& endPoint)
    {
      m_startPoint = startPoint; m_endPoint = endPoint; return *this;
    }


    /* ======= Non-member functions. ======= */

    template <class Type>
    std::istream&
    operator>>(std::istream& stream, LineSegment2D<Type>& lineSegment)
    {
      // If stream is in a bad state, we can't read from it.
      if (!stream){
        return stream;
      }

      double startPointX;
      double startPointY;
      double endPointX;
      double endPointY;

      brick::common::Expect::FormatFlag flags =
        brick::common::Expect::SkipWhitespace;
      stream >> brick::common::Expect("LineSegment2D", flags)
             >> brick::common::Expect("{", flags)
             >> startPointX
             >> brick::common::Expect(",", flags)
             >> startPointY
             >> brick::common::Expect(",", flags)
             >> endPointX
             >> brick::common::Expect(",", flags)
             >> endPointY
             >> brick::common::Expect("}", flags);

      // If reading failed, don't change lineSegment.
      if (!stream){
        return stream;
      }

      lineSegment.setValue(startPointX, startPointY,
                           endPointX, endPointY);
      return stream;
    }


    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream, const LineSegment2D<Type>& lineSegment)
    {
      stream << "LineSegment2D { "
             << lineSegment.getStartPoint().x() << ", "
             << lineSegment.getStartPoint().y() << ", "
             << lineSegment.getEndPoint().x() << ", "
             << lineSegment.getEndPoint().y() << " }";
      return stream;
    }


  } // namespace geometry

} // namespace brick

#endif /* #ifndef BRICK_GEOMETRY_LINESEGMENT2D_IMPL_HH */
