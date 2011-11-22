/**
***************************************************************************
* @file brick/numeric/stencil2D.hh
*
* Header file declaring the Stencil2D class template.
*
* Copyright (C) 2006, 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_NUMERIC_STENCIL2D_HH
#define BRICK_NUMERIC_STENCIL2D_HH

#include <vector>
#include <brick/numeric/array2D.hh>
#include <brick/numeric/index2D.hh>

namespace brick {

  namespace numeric {
    
    /** 
     * WARNING: This class is still under development and quite
     * unstable.  Use at your own risk.
     *
     * This iterator class allows STL-style interaction with the array
     * elements slected by a stencil.  For more information, please see
     * the documentation for class Stencil2D.
     */
    template <class Type, int Size>
    class StencilIterator
    {
    public:

      /**
       * This constructor is called by Stencil2D when creating a
       * StencilIterator instance to return from Stencil2D::begin() or
       * Stencil2D::end().
       *
       * @param ptr This argument is a pointer to the data of the
       * underlying array.
       * 
       * @param incrementPtr This argument is a pointer to a C-style
       * array of spacings between subsequent stencil elements.  Adding
       * these numbers to ptr will cause ptr to point to each stencil
       * element in turn.
       * 
       * @param position This argument 
       */
      StencilIterator(Type* ptr, int* incrementPtr, size_t position=0)
        : m_incrementPtr(incrementPtr),
          m_position(position),
          m_ptr(ptr) {}

      /** 
       * This operator returns true if its argument refers to a
       * different memory location *this.
       * 
       * @param other This argument is the StencilIterator instance to
       * which *this should be compared.
       * 
       * @return The return value is true if the two StencilIterator
       * instances refer to diferent referents, false otherwise.
       */
      bool
      operator!=(const StencilIterator<Type, Size>& other) {
        return m_ptr != other.m_ptr;
      }

    
      /** 
       * This operator returns true if its argument refers to the same
       * memory location *this.
       * 
       * @param other This argument is the StencilIterator instance to
       * which *this should be compared.
       * 
       * @return The return value is true if the two StencilIterator
       * instances refer to the same referents, false otherwise.
       */
      bool
      operator==(const StencilIterator<Type, Size>& other) {
        return m_ptr == other.m_ptr;
      }

    
      /** 
       * The pre-increment operator moves the iterator to the next
       * element in the sequence.
       * 
       * @return A reference to the incremented StencilIterator instance.
       */
      StencilIterator<Type, Size>&
      operator++() {m_ptr += m_incrementPtr[m_position++]; return *this;}

    
      /** 
       * The post-increment operator is just like the pre-increment
       * operator, above, but returns a copy of the StencilIterator which
       * has not been incremented.
       * 
       * @return An un-incremented copy of *this.
       */
      StencilIterator<Type, Size>
      operator++(int) {
        StencilIterator<Type, Size> result(*this);
        ++(*this);
        return result;
      }


      /** 
       * The pre-decrement operator moves the iterator to the previous
       * element in the sequence.
       * 
       * @return A reference to the decremented StencilIterator instance.
       */
      StencilIterator<Type, Size>&
      operator--() {m_ptr += m_incrementPtr[--m_position]; return *this;}

    
      /** 
       * The post-decrement operator is just like the pre-decrement
       * operator, above, but returns a copy of the StencilIterator which
       * has not been decremented.
       * 
       * @return An un-decremented copy of *this.
       */
      StencilIterator<Type, Size>
      operator--(int) {
        StencilIterator<Type, Size> result(*this);
        --(*this);
        return result;
      }


      /** 
       * The dereference operator permits reading and writing the value
       * to which this StencilIterator instance currently points.  For
       * example, you might write:
       *
       *   StencilIterator<char, 10> ptr(...);
       *   char ch = *ptr;
       *
       * just as you would with a native pointer type.
       * 
       * @return Reference to the value to which this StencilIterator
       * instance currently points.
       */
      Type&
      operator*() {return *m_ptr;}

    
      /** 
       * The access operator allows you to access the member functions
       * and member variables of the value to which this StencilIterator
       * instance currently points, just as you would with a normal
       * pointer.  Use it like this:
       *
       *   StencilIterator<MyClass, 10> ptr(...);
       *   ptr->someMethod();
       * 
       * @return C-style pointer to the value pointed to by this
       * StencilIterator instance.
       */
      Type* operator->() {return m_ptr;}

    private:
      int* m_incrementPtr;
      size_t m_position;
      Type* m_ptr;
    };


    /** 
     * WARNING: This class is still under development and quite
     * unstable.  Use at your own risk.
     *
     * This class makes it easy to access Array2D elements in a
     * repeating pattern which moves around the array.  It is useful
     * for implementing things like 2D convolution or Bayer filtering.
     *
     * Template argument Type specifies the element type of the
     * Array2D instance.
     *
     * Template argument Size specifies a maximum number of elements
     * in the stencil pattern.  For example, you must set Size greater
     * than or equal to the number of elements in your convolution
     * kernel or Bayer filtering pattern.  Setting Size too big simply
     * increases memory usage.  Setting Size too small will very
     * likely cause a segfault.
     *
     * For an example of how to use Stencil2D, please refer to
     * brick/numeric/convolve2D.hh.
     */
    template <class Type, int Size>
    class Stencil2D {
    public:
      /** 
       * The default constructor creates a blank Stencil2D instance.
       */
      Stencil2D();


      /** 
       * This constructor creates a Stencil2D instance which accesses
       * each element of a (rows x columns) pattern in raster order
       * from top left to bottom right.
       * 
       * @param rows This argument specifies the height of the stencil.
       * 
       * @param columns This argument specifies the width of the stencil.
       */
      Stencil2D(size_t rows, size_t columns);
    

      /** 
       * This constructor creates a Stencil2D instance which accesses
       * a pattern of elements that has the same shape as its
       * argument, but includes only those elements for which the
       * corresponding element of argument pattern is true.  Elements
       * wil be accessed in raster order from top left to bottom
       * right.
       * 
       * @param pattern This argument is a 2D array of bools
       * indicating which elements should be accessed as part of the
       * stencil pattern.
       */
      Stencil2D(const Array2D<bool>& pattern);


      /** 
       * This constructor is like Stencil2D(const Array2D<bool>&),
       * except that the elements which make up the stencil are
       * specified explicitly by passing in a vector of their
       * coordinates with respect to the pattern origin (top left
       * corner).  Column numbers increase from left to right, and row
       * numbers increase from top to bottom.
       * 
       * @param pattern This argument is a vector of element
       * coordinates from which to construct the stencil.
       */
      Stencil2D(const std::vector<Index2D>& pattern);
    

      /** 
       * Destructor.
       */
      ~Stencil2D() {}


      /** 
       * This member function moves the stencil one column to the
       * right, but does not check to see whether such a move is
       * valid.  You can easily advance off right edge of the array if
       * you're not careful.
       */
      inline void
      advance() {++m_ptr;}
    

      /** 
       * This member functio refers an iterator pointing to the first
       * element in the stencil sequence.  In order for the result to
       * be valid, member function setTarget() must have been called.
       * 
       * @return The return value is a Stencil2DIterator instance
       * pointing to the begining of the stencil sequence.
       */
      StencilIterator<Type, Size>
      begin() {return StencilIterator<Type, Size>(
                 m_ptr + m_offsetArray[0], m_incrementArray, 0);}


      /** 
       * This member functio refers an iterator pointing one element
       * past the final element of the stencil sequence.  In order for the
       * result to be valid, member function setTarget() must have
       * been called.
       * 
       * @return The return value is an "end" iterator to use along
       * with member function begin() in stl-style algorithms.
       */
      StencilIterator<Type, Size>
      end() {return StencilIterator<Type, Size>(
               m_ptr + m_targetSize, m_incrementArray, m_numberOfElements);}
    

      /** 
       * This member function returns a reference to a single stencil
       * element.  In order for the result to be valid, member
       * function setTarget() must have been called.
       * 
       * @param elementNumber This argument specifies which element is
       * to be returned.
       * 
       * @return The return value is a reference to the requested
       * element.
       */
      inline Type&
      getReference(size_t elementNumber);


      /** 
       * This member function returns the value of a single stencil
       * element.  In order for the result to be valid, member
       * function setTarget() must have been called.
       * 
       * @param elementNumber This argument specifies which element is
       * to be returned.
       * 
       * @return The value of the requested element is returned.
       */
      inline Type
      getValue(size_t elementNumber) const;


      /** 
       * This member function moves the upper left corner of the
       * stencil to the specified location in the target Array2D
       * instance.
       * 
       * @param index This argument specifies the requested location
       * as described for Array2D::operator()(size_t).
       */
      inline void
      goTo(size_t index);


      /** 
       * This member function moves the upper left corner of the
       * stencil to the specified location in the target Array2D
       * instance.
       * 
       * @param row This argument specifies row coordinate of the
       * requested location.
       *
       * @param column This argument specifies column coordinate of the
       * requested location.
       */
      inline void
      goTo(size_t row, size_t column);


      /** 
       * This member function sets the stencil pattern as described in
       * the documentation of Stencil2D<Type>::Stencil2D(const
       * Array2D<bool>&).
       * 
       * @param pattern This argument is a 2D array of bools
       * indicating which elements should be accessed as part of the
       * stencil pattern.
       */
      void
      setPattern(const Array2D<bool>& pattern);


      /** 
       * This member function sets the Array2D instance to which the
       * stencil will be applied.
       * 
       * @param target This argument is the Array2D instance whose
       * elements should be accessed.
       */
      template <class Type2>
      void
      setTarget(const Array2D<Type2>& target);


      /** 
       * This member function sets the Array2D instance to which the
       * stencil will be applied.
       * 
       * @param target This argument is the Array2D instance whose
       * elements should be accessed.
       */
      template <class Type2>
      void
      setTarget(Array2D<Type2>& target);
    
      // private:

      inline void
      checkBounds(size_t elementNumber) const;


      // Data members to allow moving around in the target array.
      Type* m_basePtr;
      size_t m_rows;
      size_t m_columns;

      // Data members to allow bounds checking.
      int m_targetSize;

      // Data members which actually do the work.
      size_t m_numberOfElements;
      Type* m_ptr;
      int m_offsetArray[Size];
      int m_incrementArray[Size];
      Index2D m_patternArray[Size];
    };

  } // namespace numeric

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/numeric/stencil2D_impl.hh>

#endif /* #ifndef BRICK_NUMERIC_STENCIL2D_HH */
