/**
***************************************************************************
* @file brick/geometry/lineSegment2D.hh
*
* Header file declaring LineSegment2D class template.
*
* Copyright (C) 2007-2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_GEOMETRY_LINESEGMENT2D_HH
#define BRICK_GEOMETRY_LINESEGMENT2D_HH

#include <iostream>
#include <brick/numeric/vector2D.hh>

namespace brick {

  namespace geometry {

    /**
     ** The LineSegment2D class represents a line segment in 2D space.
     **/
    template <class Type>
    class LineSegment2D {
    public:

      /**
       * The default constructor initializes to the line segment that starts
       * at the origin and ends at (1, 0).
       */
      LineSegment2D();


      /**
       * This constructor initializes the line segment using a pair of points.
       *
       * @param point This argument specifies the start point of the
       * line segment.
       *
       * @param endPoint This argument specifies the end point of the
       * line segment.
       */
      LineSegment2D(brick::numeric::Vector2D<Type> const& startPoint,
                    brick::numeric::Vector2D<Type> const& endPoint);


      /**
       * This constructor initializes the line segment using four
       * coordinates representing a pair of points.
       *
       * @param x0 This argument is the x coordinate of the start
       * point of the line segment.
       *
       * @param y0 This argument  is the y coordinate of the start
       * point of the line segment.
       *
       * @param x1 This argument  is the x coordinate of the end
       * point of the line segment.
       *
       * @param y1 This argument  is the y coordinate of the end
       * point of the line segment.
       */
      LineSegment2D(Type const& x0, Type const& y0,
                    Type const& x1, Type const& y1);


      /**
       * The copy constructor deep copies its argument.
       *
       * @param source This argument is the class instance to be
       * copied.
       */
      LineSegment2D(LineSegment2D<Type> const& source);


      /**
       * Destructor.
       */
      ~LineSegment2D();


      /**
       * The assignment operator deep copies its argument.
       *
       * @param source This argument is the class instance to be
       * copied.
       *
       * @return The return value is a reference to *this.
       */
      LineSegment2D<Type>&
      operator=(LineSegment2D<Type> const& source);


      /**
       * This member function returns the end point of the line segment.
       *
       * @return The return value is a const reference to Vector2D
       * representing the end point of the line segment.
       */
      brick::numeric::Vector2D<Type> const&
      getEndPoint() const;


      /**
       * This member function returns the start point of the line segment.
       *
       * @return The return value is a const reference to Vector2D
       * representing the start point of the line segment.
       */
      brick::numeric::Vector2D<Type> const&
      getStartPoint() const;


      /**
       * This member function returns the start point of the line segment.
       *
       * @return The return value is a const reference to Vector2D
       * representing the start point of the line segment.
       */
      brick::numeric::Vector2D<Type> const&
      getVertex0() const;


      /**
       * This member function returns the end point of the line segment.
       *
       * @return The return value is a const reference to Vector2D
       * representing the end point of the line segment.
       */
      brick::numeric::Vector2D<Type> const&
      getVertex1() const;


      /**
       * This member function changes the start point and end point of
       * the line segment.
       *
       * @param startPointX This argument is the new start point X coordinate.
       *
       * @param startPointY This argument is the new start point Y coordinate.
       *
       * @param endPointX This argument is the new end point X coordinate.
       *
       * @param endPointY This argument is the new end point Y coordinate.
       *
       * @return The return value is a reference to *this after the
       * change has been applied.
       */
      LineSegment2D<Type>&
      setValue(Type const& startPointX, Type const& startPointY,
               Type const& endPointX, Type const& endPointY);


      /**
       * This member function changes the start point and end point of
       * the line segment.
       *
       * @param startPoint This argument is the new start point.
       *
       * @param endPoint This argument is the new end point.
       *
       * @return The return value is a reference to *this after the
       * change has been applied.
       */
      LineSegment2D<Type>&
      setValue(brick::numeric::Vector2D<Type> const& startPoint,
               brick::numeric::Vector2D<Type> const& endPoint);


   private:
      // Private member functions.

      // Private data members.
      brick::numeric::Vector2D<Type> m_startPoint;
      brick::numeric::Vector2D<Type> m_endPoint;

    }; // class LineSegment2D



    /* ======= Non-member functions. ======= */

    template <class Type>
    std::istream&
    operator>>(std::istream& stream, LineSegment2D<Type>& lineSegment);


    template <class Type>
    std::ostream&
    operator<<(std::ostream& stream, LineSegment2D<Type> const& lineSegment);


  } // namespace utilities

} // namespace brick


// Include definitions of inline and template functions.
#include <brick/geometry/lineSegment2D_impl.hh>

#endif /* #ifndef BRICK_GEOMETRY_LINESEGMENT2D_HH */
