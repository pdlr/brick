/**
***************************************************************************
* @file brick/common/stridedPointer.hh> 
* Header file declaring StridedPointer class.
*
* Copyright (C) 2002-2003 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_COMMON_STRIDEDPOINTER_HH
#define BRICK_COMMON_STRIDEDPOINTER_HH

#include <iterator>

namespace brick {

  namespace common {
    
    /**
     ** The StridedPointer class permits convenient iterator-style
     ** access to regularly spaced elements within a C-style array.
     ** That is, a StridedPointer instance with stride == 1 acts a lot
     ** like an ordinary pointer, while a StridedPointer instance with
     ** stride == 2 "sees" only every other element of the array, and so
     ** on.
     **/
    template <class Type>
    class StridedPointer
      : public std::iterator<std::random_access_iterator_tag, Type, ptrdiff_t>
    {
    public:

      /** 
       * The default constructor creates a null pointer with stride ==
       * 1.
       */
      StridedPointer() : m_pointer(0), m_stride(1) {}

      
      /** 
       * This constructor creates a StridedPointer which points to the
       * same memory as the argument, and has stride == 1.
       * 
       * @param pointer The pointer from which to take the initial value
       * of *this.
       */
      StridedPointer(Type* pointer) : m_pointer(pointer), m_stride(1) {}

      
      /** 
       * This constructor creates a StridedPointer which has the same
       * value as the argument, but has the specified stride.
       * 
       * @param pointer The pointer from which to take the initial value
       * of *this.
       * @param stride The StridedPointer instance will be able to
       * access only every stride'th element in the array.
       */
      StridedPointer(Type* pointer, int stride)
        : m_pointer(pointer), m_stride(stride) {}

      
      /** 
       * This is the copy constructor.  It copies both address and
       * stride from its argument.
       * 
       * @param source The StridedPointer instance to be copied.
       */
      StridedPointer(const StridedPointer<Type>& source)
        : m_pointer(source.m_pointer), m_stride(source.m_stride) {}

      
      /** 
       * Destroys the StridedPointer instance.
       */
      ~StridedPointer() {}

      
      /** 
       * The dereference operator permits reading and writing the value
       * to which this StridedPointer instance currently points.  For
       * example, you might write:
       *
       * @code
       *   StridedPointer<char> ptr(...);
       *   char ch = *ptr;
       * @endcode
       *
       * just as you would with a native pointer type.
       * 
       * @return Reference to the value to which this StridedPointer
       * instance currently points.
       */
      Type&
      operator*() {return *m_pointer;}

      
      /** 
       * The access operator allows you to access the member functions
       * and member variables of the value to which this StridedPointer
       * instance currently points, just as you would with a normal
       * pointer.  Use it like this:
       *
       * @code
       *   StridedPointer<myClass> ptr(...);
       *   ptr->someMethod();
       * @endcode
       * 
       * @return C-style pointer to the value pointed to by this
       * StridedPointer instance.
       */
      Type* operator->() {return m_pointer;}

      
      /** 
       * The pre-increment operator actually adds the value of stride
       * (which was set at construction time) to *this.  That is, if the
       * value of stride is N, then incrementing a StridedPointer
       * instance makes it point to an array element N spaces away.
       * 
       * @return A Reference to the incremented StridedPointer instance.
       */
      StridedPointer<Type>& operator++() {m_pointer += m_stride; return *this;}

      
      /** 
       * The post-increment operator is just like the pre-increment
       * operator, above, but returns a copy of the StridedPointer which
       * has not been incremented.
       * 
       * @return An un-incremented copy of *this.
       */
      StridedPointer<Type> operator++(int) {
        StridedPointer<Type> result(*this);
        ++(*this);
        return result;
      }

      
      /** 
       * The pre-decrement operator is just like the pre-increment
       * operator, above, but decrements instead of incrementing.
       * 
       * @return A Reference to the decremented StridedPointer instance.
       */
      StridedPointer<Type>& operator--() {m_pointer -= m_stride; return *this;}


      /** 
       * The post-decrement operator is just like the post-increment
       * operator, above, but decrements instead of incrementing.
       * 
       * @return An un-decremented copy of *this.
       */
      StridedPointer<Type> operator--(int) {
        StridedPointer<Type> result(*this);
        --(*this);
        return result;
      }

      
      /** 
       * This operator has the same effect as incrementing a copy of
       * *this offset times, but runs in O(1) time.  That is, a call
       * to operator+() returns a StridedPointer instance which
       * references a place in memory that is (stride * offset) items
       * advanced from the value referenced by *this.
       * 
       * @param offset This argument specifies how many steps to
       * advance.
       * @return A copy of *this that references the new location.
       */
      StridedPointer<Type> operator+(ptrdiff_t offset) {
        StridedPointer<Type> result(*this);
        result += offset;
        return result;
      }

      
      /** 
       * This operator works just like operator+() above, but decreases
       * the StridedPointer value instead of increasing it.
       * 
       * @param offset This argument specifies how many steps to
       * move back.
       * @return A copy of *this which references the new location.
       */
      StridedPointer<Type> operator-(ptrdiff_t offset) {
        StridedPointer<Type> result(*this);
        result -= offset;
        return result;
      }

      
      /** 
       * This operator attempts to return a value N so that *(*this)
       * == *(other + N).  This value is meaningless if the two
       * StridedPointer instances have different stride values.  Note
       * that if stride is not equal to 1, it is very possible to have
       * two StridedPointer instances which are "out of phase" so that
       * there is no value of N which satisfies the above equation.
       * 
       * @param other The StridedPointer instance to subtract from
       * *this.
       * @return A difference value as described above.
       */
      ptrdiff_t operator-(const StridedPointer<Type>& other) {
        ptrdiff_t distance = m_pointer - other.m_pointer;
        return distance / m_stride;
      }

      
      /** 
       * This operator has the same effect as incrementing *this offset
       * times, but runs in O(1) time.
       * 
       * @param offset This argument specifies how many steps to
       * advance the StridedPointer.
       * @return A reference to *this.
       */
      StridedPointer<Type>& operator+=(ptrdiff_t offset) {
        m_pointer += offset * m_stride;
        return *this;
      }

      
      /** 
       * This operator has the same effect as decrementing *this offset
       * times, but runs in O(1) time.
       * 
       * @param offset This argument specifies how many steps to
       * decrement the StridedPointer.
       * @return A reference to *this.
       */
      StridedPointer<Type>& operator-=(ptrdiff_t offset) {
        m_pointer -= offset * m_stride;
        return *this;
      }

      
      /** 
       * This operator returns true if its argument references the
       * same memory location as *this, false otherwise.  Note that it
       * does not require the two StridedPointer instances to have the
       * same stride.  It is sufficient that they reference the same
       * memory location.
       * 
       * @param other A StridedPointer instance to be compared with *this.
       * @return true if other references same memory location as *this,
       * false otherwise.
       * 
       * @return The return value is true if the two StridedPointer
       * instances point to the same memory location, false otherwise.
       */
      bool operator==(const StridedPointer<Type>& other) const {
        return m_pointer == other.m_pointer;
      }

      
      /** 
       * This operator returns false if its argument references the same
       * memory location as *this, true otherwise.
       * 
       * @param other A StridedPointer instance to be compared with *this.
       * @return false if other references same memory location as *this,
       * true otherwise.
       */
      bool operator!=(const StridedPointer<Type>& other) const {
        return m_pointer != other.m_pointer;
      }

      
      /** 
       * This operator return true if the address of the memory location
       * referenced by *this is less than the address of the memory
       * location referenced by other, false otherwise.
       * 
       * @param other The StridedPointer instance to be compared with
       * *this.
       * @return true if the address of the memory location referenced
       * by *this is less than the address of the memory location
       * referenced by other, false otherwise.
       */
      bool operator<(const StridedPointer<Type>& other) const {
        return m_pointer < other.m_pointer;
      }

      
      /** 
       * This operator return true if the address of the memory location
       * referenced by *this is greater than the address of the memory
       * location referenced by other, false otherwise.
       * 
       * @param other The StridedPointer instance to be compared with
       * *this.
       * @return true if the address of the memory location referenced
       * by *this is greater than the address of the memory location
       * referenced by other, false otherwise.
       */
      bool operator>(const StridedPointer<Type>& other) const {
        return m_pointer > other.m_pointer;
      }

      
      /** 
       * This operator return true if the address of the memory location
       * referenced by *this is less than or equal to the address of the
       * memory location referenced by other, false otherwise.
       * 
       * @param other The StridedPointer instance to be compared with
       * *this.
       * @return true if the address of the memory location referenced
       * by *this is less than or equal to the address of the memory
       * location referenced by other, false otherwise.
       */
      bool operator<=(const StridedPointer<Type>& other) const {
        return m_pointer <= other.m_pointer;
      }

      
      /** 
       * This operator return true if the address of the memory location
       * referenced by *this is greater than or equal to the address of
       * the memory location referenced by other, false otherwise.
       * 
       * @param other The StridedPointer instance to be compared with
       * *this.
       * @return true if the address of the memory location referenced
       * by *this is greater than or equal to the address of the memory
       * location referenced by other, false otherwise.
       */
      bool operator>=(const StridedPointer<Type>& other) const {
        return m_pointer >= other.m_pointer;
      }

      
      /** 
       * The assignment operator copies both address and stride from its
       * argument.
       * 
       * @param source The StridedPointer instance to be copied.
       * @return A reference to *this.
       */
      StridedPointer<Type>& operator=(const StridedPointer<Type>& source) {
        m_pointer = source.m_pointer;
        m_stride = source.m_stride;
        return *this;
      }

    private:
      Type* m_pointer;
      int m_stride;
    };

  } // namespace common
  
} // namespace brick

#endif /* #ifndef BRICK_COMMON_STRIDEDPOINTER_HH */
