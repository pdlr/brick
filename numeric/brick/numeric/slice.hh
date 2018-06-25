/**
***************************************************************************
* @file brick/numeric/slice.hh
*
* Header file declaring Slice class.
*
* Copyright (C) 2002-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_NUMERIC_SLICE_HH
#define BRICK_NUMERIC_SLICE_HH

namespace brick {

  namespace numeric {

    /**
     ** A simple Slice class to work with SubArrays.  This is modeled
     ** after std::slice, but has "stop" instead of "size," and permits
     ** negative indexing.
     **/
    class Slice {
    public:
      /**
       * Default constructor initializes everything to zero.  Some
       * classes interpret this as meaning "every row" or "every
       * column."
       */
      Slice() : m_start(0), m_stop(0), m_stride(1) {}

      /**
       * This constructor sets start and stop as specified, and defaults
       * stride to 1.
       *
       * @param start This argument specifies the desired first element.
       * @param stop This argument specifies an element one beyond the
       * desired last element.
       */
      Slice(int startIndex, int stopIndex)
        : m_start(startIndex), m_stop(stopIndex), m_stride(1) {}

      /**
       * This constructor lets the caller specify all three values.
       *
       * @param start This argument specifies the desired first element.
       * @param stop This argument specifies an element one beyond the
       * desired last element.
       * @param stride0 This argument specifies how many elements to skip
       * with each step.
       */
      Slice(int startIndex, int stopIndex, int stride0)
        : m_start(startIndex), m_stop(stopIndex), m_stride(stride0) {}

      /**
       * Accessor function for the value of start.
       *
       * @return The index of the first element in the slice.
       */
      inline int getStart() const {return m_start;}
      inline int start() const {return m_start;}

      /**
       * Accessor function for the value of stop.
       *
       * @return The index of an element one beyond the last element in
       * the slice.
       */
      inline int getStop() const {return m_stop;}
      inline int stop() const {return m_stop;}

      /**
       * Accessor function for the value of stride.
       *
       * @return The value of stride, where the slice accesses every
       * stride'th element from start to just before stop.
       */
      inline int getStride() const {return m_stride;}
      inline int stride() const {return m_stride;}

    private:
      int m_start;
      int m_stop;
      int m_stride;
    };

  } // namespace numeric

} // namespace brick

#endif /* #ifdef BRICK_NUMERIC_SLICE_HH */
