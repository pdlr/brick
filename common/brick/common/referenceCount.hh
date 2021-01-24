/**
***************************************************************************
* @file brick/common/referenceCount.hh
*
* Header file declaring ReferenceCount class.
*
* Copyright (C) 2003-2010 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#ifndef BRICK_COMMON_REFERENCECOUNT_HH
#define BRICK_COMMON_REFERENCECOUNT_HH

#include <cstddef>

namespace brick {

  namespace common {

    /**
     ** The ReferenceCount class provides a convenient way to track a
     ** shared resource so you know when to delete it.  ReferenceCount
     ** has two states: counted and uncounted.  It defaults to, and is
     ** normally used in, the counted state.
     **
     ** In the counted state, ReferenceCount maintains an internal
     ** count, which defaults to 1 at construction time.  Copying the
     ** ReferenceCount instance increments the count, and destroying
     ** the instance or any of its copies decrements the count.  If
     ** the user initializes the count to 1 when constructing the
     ** ReferenceCount instance, and then copies the instance whenever
     ** the pointer to the shared resource is copied, then the count
     ** provides a way of tracking how many copies of the shared
     ** pointer are in use.  Increment and decrement operators are
     ** provided so that the user can explicitly change the count if
     ** necessary, but be very careful if you do this, as you can (for
     ** example) fool the ReferenceCount instance into thinking that
     ** no other ReferenceCount instances share the same count,
     ** leading to memory errors when you destroy it.  Explicitly
     ** decrementing the count past 1 is legal, but it's not clear why
     ** you would want to do this.  Destroying the ReferenceCount with
     ** a count of 1 or less will delete the internal count.
     **
     ** Calling the reset() method while in the counted state
     ** decrements the count (affecting any copies), then abandons the
     ** count and starts a new one. Calling reset() with argument 0
     ** puts the ReferenceCount instance into the uncounted state.
     ** Calling the reset() method with an argument greater than 1
     ** starts the count at the specified value, indicating that other
     ** references to the shared resource already exist.
     **
     ** In the uncounted state, the increment and decrement operators
     ** have no effect, copying the ReferenceCount instance simply
     ** creates another uncounted ReferenceCount instance, and
     ** destroying copies has no effect on the count.  This is useful
     ** if you need to create a ReferenceCount instance now, but only
     ** initialize the count (using reset()) later when the shared
     ** resource is allocated.  The only way to move from the
     ** uncounted state to the counted state is to call the reset()
     ** method with no argument, or with a nonzero argument.
     **
     ** Here is a simple example of how you might use ReferenceCount
     ** to implement an vector class with automatically managed
     ** shallow copy semantics:
     **
     ** @code
     **   template <class ElementType>
     **   class MyVector {
     **     MyVector()
     **       : m_vectorPtr(new std::vector<ElementType>),
     **         m_referenceCount() {}
     **
     **     ~MyVector() {
     **       if(!m_referenceCount.isShared()) {
     **         delete m_vectorPtr;
     **       }
     **     }
     **
     **   // ** Public interface for element access, etc. goes here. **
     **
     **   private:
     **     std::vector* m_vectorPtr;
     **     ReferenceCount m_referenceCount;
     **   };
     ** @endcode
     **
     ** In this example, the compiler-generated copy constructor and
     ** assignment operator will directly copy the pointer to
     ** m_vectorPtr.  The will also call the ReferenceCount copy
     ** constructor and assignment operator, which will update the
     ** internal state of m_referenceCount so that both the copied
     ** instance of MyVector and the copying instance of MyVector know
     ** that m_vectorPtr is shared.  Subsequent copies will further
     ** increment the count inside m_referenceCount.  When the copies
     ** are destroyed, the ReferenceCount destructor will decrement
     ** this count until, when the very last copy is destroyed,
     ** m_vectorPtr will be deleted.
     **
     ** One note regarding thread safety: if you use ReferenceCount as
     ** a member of a class that is shared between threads, and this
     ** class manages its own mutex locking (such as
     ** brick::thread::Monitor), construction and destruction of the
     ** members of the containing class will most likely happen
     ** outside of the locked sections of code.  This means that the
     ** ReferenceCount copy constructor and destructor will be called
     ** without protection, which is not thread-safe.  To solve this,
     ** you can dynamically allocate the ReferenceCount instance using
     ** new, keep a pointer to the ReferenceCount as a member of the
     ** class, and explicitly manage copying, incrementing, and
     ** decrementing of the ReferenceCount inside the locked sections
     ** of code.  Alternatively, you could do the locking outside the
     ** class (perhaps using brick::thread::Monitor) so that
     ** constructor and destructor calls are protected.
     **/
    class ReferenceCount
    {
    public:
      /**
       * The default constructor sets the reference count to 1,
       * indicating a counted and unshared condition.  To indicate
       * that there are existing references, pass a constructor
       * argument greater than 1.  Use an argument of zero to indicate
       * that the shared data has not been initialized, and leave
       * *this in the uncounted state.
       *
       * @param count This argument specifies to what value the count
       * should be initialized.  For most applications, this argument
       * should be set to 1.  Setting this argument to zero indicates
       * that no reference counting should be done (until a subsequent
       * call to the reset() method).
       */
      ReferenceCount(size_t count=1)
        : m_countPtr(0) {
        this->reset(count);
      }


      /**
       * This is the copy constructor.  If *this is in the counted
       * state, then after copying both ReferenceCount instances share
       * the same count, and the count is incremented by one.  If
       * *this is in the uncounted state, then the copy will also be
       * uncounted.
       *
       * @param other The ReferenceCount instance to be copied.
       */
      ReferenceCount(const ReferenceCount& other)
        : m_countPtr(other.m_countPtr) {
        ++(*this);
      }


      /**
       * Decrements the count (if the ReferenceCount instance is in
       * the counted state) and destroys the ReferenceCount instance.
       */
      ~ReferenceCount() {
        --(*this);
        this->deleteIfNecessary();
      }


      /**
       * The pre-increment operator increments the count by one, if
       * the ReferenceCount instance is in the counted state, and has
       * no impact in the uncounted state.
       *
       * @return A Reference to the incremented ReferenceCount instance.
       */
      ReferenceCount&
      operator++() {
        if(m_countPtr != 0) {++(*m_countPtr);}
        return *this;
      }


      /**
       * Due to the semantics of the ReferenceCount class, the
       * post-increment operator is identical to the the pre-increment
       * operator.
       *
       * @return A Reference to the incremented ReferenceCount instance.
       */
      ReferenceCount
      operator++(int) {return ++(*this);}


      /**
       * The pre-decrement operator decrements the count by one, if
       * the ReferenceCount instance is in the counted state, and has
       * no impact in the uncounted state.
       *
       * @return A Reference to the decremented ReferenceCount instance.
       */
      ReferenceCount&
      operator--() {
        if(m_countPtr != 0) {--(*m_countPtr);}
        return *this;
      }


      /**
       * Due to the semantics of the ReferenceCount class, the
       * post-decrement operator is identical to the the pre-decrement
       * operator.
       *
       * @return A Reference to the decremented ReferenceCount instance.
       */
      ReferenceCount
      operator--(int) {return --(*this);}


      /**
       * This operator has the same effect as incrementing *this offset
       * times, but runs in O(1) time.
       *
       * @param offset This argument specifies how many times to increment
       * @return A reference to *this.
       */
      ReferenceCount&
      operator+=(size_t offset) {
        if(m_countPtr != 0) {(*m_countPtr) += offset;}
        return *this;
      }


      /**
       * This operator has the same effect as decrementing *this offset
       * times, but runs in O(1) time.
       *
       * @param offset This argument specifies how many times to decrement
       * @return A reference to *this.
       */
      ReferenceCount&
      operator-=(size_t offset) {
        if(m_countPtr != 0) {(*m_countPtr) -= offset;}
        return *this;
      }


      /**
       * The assignment operator copies its argument.  After copying a
       * ReferenceCount instance that is in the counted state, both
       * ReferenceCount instances share the same count, and the count
       * is incremented by one.  After copying a ReferenceCount
       * instance that is in the uncounted state, both ReferenceCount
       * instances are uncounted.
       *
       * @param source The ReferenceCount instance to be copied.
       * @return A reference to *this.
       */
      ReferenceCount&
      operator=(const ReferenceCount& source) {
        // Check for self-assignment.
        if (this != &source) {
          // Release the count that was previously tracked by *this.
          --(*this);
          this->deleteIfNecessary();

          // Adopt the new count and increment it.
          m_countPtr = source.m_countPtr;
          ++(*this);
        }
        return *this;
      }


      /**
       * This member function returns the current count if the
       * ReferenceCount instance is in the counted state, or 0
       * otherwise.
       *
       * @return The internal reference count.
       */
      int
      getCount() const {
        if(this->isCounted()) {return *m_countPtr;}
        return 0;
      }


      /**
       * This member function returns true if the ReferenceCount
       * instance is in the counted state (see class documentation for
       * ReferenceCount).
       *
       * @return true if *this is in the counted state.
       */
      bool
      isCounted() const {return this->m_countPtr != 0;}


      /**
       * This member function returns true if more than one
       * ReferenceCount object is sharing the count with *this.  That
       * is, it returns true if *this is in the counted state, and the
       * internal count is greater than 1.
       *
       * @return true if the internal count is greater than 1.
       */
      bool
      isShared() const {return (this->getCount() > 1);}


      /**
       * This member function decrements the count and releases the
       * reference, then reinitializes with a fresh count.  Use this
       * if you want to release the shared resource and create a new
       * one.
       *
       * @param count This argument specifies to what value the count
       * should be reinitialized.  For most applications, this argument
       * should be set to 1.
       */
      void
      reset(size_t count=1) {
        --(*this);
        this->deleteIfNecessary();
        if(count != 0) {
          m_countPtr = new int;
          *m_countPtr = count;
        }
      }


    private:

      /**
       * This member function deletes the internal count pointer
       * pointer if no references remain, and (always) resets the
       * pointer to 0.
       */
      void deleteIfNecessary() {
        if(m_countPtr != 0) {
          if((*m_countPtr) <= 0) {
            delete m_countPtr;
          }
          m_countPtr = 0;
        }
      }


      int* m_countPtr;
    };

  } // namespace common

} // namespace brick

#endif /* #ifndef BRICK_COMMON_REFERENCECOUNT_HH */
