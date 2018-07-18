/**
***************************************************************************
* @file brick/computerVision/disjointSet_impl.hh
*
* Header file defining inline and template functions declared in
* disjointSet.hh.
*
* Copyright (C) 2008,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_DISJOINTSET_IMPL_HH
#define BRICK_COMPUTERVISION_DISJOINTSET_IMPL_HH

// This file is included by disjointSet.hh, and should not be directly included
// by user code, so no need to include disjointSet.hh here.
//
// #include <brick/computerVision/disjointSet.hh>


namespace brick {

  namespace computerVision {

    // The default constructor creates a leaf in a disjoint set
    // tree.
    template <class Type>
    DisjointSet<Type>::
    DisjointSet()
      : m_parentPtr(this), m_rank(0), m_size(1), m_payload()
    {
      // Empty.
    }


    // This constructor allows you to explicitly set the payload
    // value.
    template <class Type>
    DisjointSet<Type>::
    DisjointSet(const Type& payload)
      : m_parentPtr(this), m_rank(0), m_size(1), m_payload(payload)
    {
      // Empty.
    }


    // Destructor.
    template <class Type>
    DisjointSet<Type>::
    ~DisjointSet()
    {
      // Empty.
    }


    // This member function returns a reference to the head of the
    // set to which *this belongs.
    template <class Type>
    DisjointSet<Type>&
    DisjointSet<Type>::
    find()
    {
      // Follow parent pointers to find the head of the tree.
      DisjointSet<Type>* headPtr = this;
      while(headPtr->m_parentPtr != headPtr) {
        headPtr = headPtr->m_parentPtr;
      }

      // Now do the same thing again, but this time update all of
      // the parent pointers to point directly at the head.
      DisjointSet<Type>* updatePtr = this;
      while(updatePtr->m_parentPtr != updatePtr) {
        DisjointSet* tmpPtr = updatePtr->m_parentPtr;
        updatePtr->m_parentPtr = headPtr;
        updatePtr = tmpPtr;
      }
      return *m_parentPtr;
    }


    // This member function implements a crippled version of find
    // that doesn't fully update parent pointers.
    template <class Type>
    DisjointSet<Type>&
    DisjointSet<Type>::
    findNoUpdate()
    {
      DisjointSet<Type>* headPtr = this;
      while(headPtr->m_parentPtr != headPtr) {
        headPtr = headPtr->m_parentPtr;
      }
      m_parentPtr = headPtr;
      return *m_parentPtr;
    }


    // This member function returns a reference to the head of the
    // set to which *this belongs.
    template <class Type>
    DisjointSet<Type>&
    DisjointSet<Type>::
    findRecursive()
    {
      if(m_parentPtr == this) {
        return *this;
      }
      m_parentPtr = &(m_parentPtr->findRecursive());
      return *m_parentPtr;
    }


    // This member function merges two sets.
    template <class Type>
    void
    DisjointSet<Type>::
    merge(DisjointSet& other)
    {
      DisjointSet& myParent = this->find();
      DisjointSet& otherParent = other.find();
      if(&myParent == &otherParent) {
        return;
      }
      if(myParent.m_rank > otherParent.m_rank) {
        otherParent.m_parentPtr = &myParent;
        myParent.m_size += otherParent.m_size;
      } else if(myParent.m_rank < otherParent.m_rank) {
        myParent.m_parentPtr = &otherParent;
        otherParent.m_size += myParent.m_size;
      } else {
        otherParent.m_parentPtr = &myParent;
        myParent.m_size += otherParent.m_size;
        ++(myParent.m_rank);
      }
    }


    // This member function returns the payload associated with
    // *this.
    template <class Type>
    const Type&
    DisjointSet<Type>::
    getPayload() const
    {
      return m_payload;
    }


    // This member function returns the number of members in the set
    // of which *this is a member.
    template <class Type>
    size_t
    DisjointSet<Type>::
    getSize()
    {
      return (this->find()).m_size;
    }


    // This member function associates a payload vaule with *this.
    template <class Type>
    void
    DisjointSet<Type>::
    setPayload(const Type& payload)
    {
      m_payload = payload;
    }

  } // namespace computerVision

} // namespace brick


#endif /* #ifndef BRICK_COMPUTERVISION_DISJOINTSET_IMPL_HH */
