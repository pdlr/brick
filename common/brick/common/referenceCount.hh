/**
***************************************************************************
* @file brick/common/referenceCount.hh
*
* Header file declaring ReferenceCount class.
*
* Copyright (C) 2003-2010 David LaRose, dlr@cs.cmu.edu
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
     ** shared resource so you know when to delete it.  Here is a
     ** simple example of how you might use ReferenceCount to
     ** implement an vector class with automatically managed shallow
     ** copy semantics:
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
     ** are deleted, the ReferenceCount destructor will decrement this
     ** count until when the very last copy is destroyed, m_vectorPtr
     ** will be deleted.
     **
     ** One note regarding thread safety: if you use ReferenceCount as
     ** a member of a class that is shared between threads, and that
     ** manages its own mutex locking (such as dlrThread::Monitor)
     ** construction and destruction of the members of the containing
     ** class will most likely happen outside of the locked sections
     ** of code.  This means that the ReferenceCount copy constructor
     ** and destructor will be called without protection, which is not
     ** thread-safe.  To solve this, you can dynamically allocate the
     ** ReferenceCount instance using new, keep a pointer to the
     ** ReferenceCount as a member of the class, and explicitly manage
     ** copying, incrementing, and decrementing of the ReferenceCount
     ** inside the locked sections of code.
     **/
    class ReferenceCount
    {
    public:
      /** 
       * The default constructor sets the reference count to 1, indicating
       * a new and unshared condition.  To indicate that there are existing
       * references, pass a constructor argument greater than 1.  Use an
       * argument of zero to indicate that the shared data has not been
       * initialized.
       * 
       * @param count This argument specifies to what value the count
       * should be initialized.  For most applications, this argument
       * should be set to 1.
       */
      ReferenceCount(size_t count=1)
        : m_countPtr(0) {
        if(count != 0) {
          m_countPtr = new size_t;
          *m_countPtr = count;
        }
      }

    
      /** 
       * This is the copy constructor.  After copying, both ReferenceCount
       * instances share the same count, and the count is incremented by one.
       * 
       * @param other The ReferenceCount instance to be copied.
       */
      ReferenceCount(const ReferenceCount& other)
        : m_countPtr(other.m_countPtr) {
        if(m_countPtr != 0) {++(*m_countPtr);}
      }

    
      /** 
       * Decrements the count and destroys the ReferenceCount instance.
       */
      ~ReferenceCount() {this->deleteIfNecessary();}

    
      /** 
       * The pre-increment operator increments the count by one.
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
       * The pre-decrement operator decrements the count by one.
       * 
       * @return A Reference to the decremented ReferenceCount instance.
       */
      ReferenceCount&
      operator--() {
        if(m_countPtr != 0) {
          if((*m_countPtr) != 0) {--(*m_countPtr);}
        }
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
        if(m_countPtr != 0) {
          if((*m_countPtr) > offset) {(*m_countPtr) -= offset;}
          else {(*m_countPtr) = 0;}
        }
        return *this;
      }

    
      /** 
       * The assignment operator copies its argument.  After copying,
       * both ReferenceCount instances share the same count, and the count
       * is incremented by one.
       * 
       * @param source The ReferenceCount instance to be copied.
       * @return A reference to *this.
       */
      ReferenceCount&
      operator=(const ReferenceCount& source) {
        // Check for self-assignment.
        if (this != &source) {
          this->deleteIfNecessary();
          m_countPtr = source.m_countPtr;
          if(m_countPtr != 0) {++(*m_countPtr);}
        }
        return *this;
      }

    
      /** 
       * This member function returns the current count.
       * 
       * @return The internal reference count.
       */
      size_t
      count() const {return *m_countPtr;}

    
      /** 
       * This member function returns true if more than one
       * ReferenceCount object is sharing the same data.  That is, it
       * returns true if the internal count is greater than 1.
       * 
       * @return true if the internal count is greater than 1.
       */
      bool
      isShared() const {
        if(m_countPtr == 0) {return false;}
        return ((*m_countPtr) > 1);
      }

    
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
        this->deleteIfNecessary();
        if(count != 0) {
          m_countPtr = new size_t;
          *m_countPtr = count;
        }
      }

      
      /** 
       * This member function is identical to member function
       * isShared().  It is included for backwards compatibility with
       * previous versions of ReferenceCount.
       * 
       * @return true if the internal count is greater than 1.
       */
      bool
      shared() const {
        return this->isShared();
      }

    
    private:

      /** 
       * This member function decrements the count and deletes the
       * corresponding pointer if no references remain.
       */
      void deleteIfNecessary() {
        if(m_countPtr != 0) {
          switch(*m_countPtr) {
          case 0:
          case 1:
            delete m_countPtr;
            m_countPtr = 0;
            break;
          default:
            --(*m_countPtr);
          }
        }
      }
    
      size_t* m_countPtr;
    };

  } // namespace common
  
} // namespace brick

#endif /* #ifndef BRICK_COMMON_REFERENCECOUNT_HH */
