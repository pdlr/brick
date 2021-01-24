/**
***************************************************************************
* @file brick/computerVision/disjointSet.hh
*
* Header file declaring a class to be used in implementing algorithms
* that need a "forest of disjoint sets" data structure.
*
* Copyright (C) 2008,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_DISJOINTSET_HH
#define BRICK_COMPUTERVISION_DISJOINTSET_HH


namespace brick {

  namespace computerVision {

    /**
     ** This class implements one tree in the "forest of disjoint
     ** sets" data structure first described by Bernard Galler and
     ** Michael Fischer[1].  This implementation includes union by
     ** rank and path compression, yielding "nearly" constant time
     ** amortized complexity for find/merge operations.
     **
     ** Our implementation includes a size member that allows the
     ** calling context to easily see how many members belong to the
     ** set.  This adds a small constant to the execution time of the
     ** merge() operation.
     **
     ** To make a forest of disjoint sets, you'll need to put a
     ** bunch of DisjointSet instances in some kind of container,
     ** such as an Array1D instance (brick::numeric::Array1D<
     ** DisjointSet<Type> >).  Unfortunately, most of the C++
     ** standard library containers, such as std::vector, use the
     ** copy constructor to fill in the array, which breaks the
     ** parent pointers.  So don't use a std::vector (or std::list,
     ** or...).
     **
     ** An Improved Equivalence Algorithm.  Bernard Galler, and
     ** Michael Fischer. Communications of the ACM, Volume 7, pages
     ** 301-303, 1964.
     **/
    template <class Type>
    class DisjointSet {
    public:

      /**
       * The default constructor creates a leaf in a disjoint set
       * tree.
       */
      DisjointSet();


      /**
       * This constructor allows you to explicitly set the payload
       * value.
       *
       * @param payload This argument will be copied into the
       * payload, and then ignored until you call
       * this->getPayload().
       */
      DisjointSet(const Type& payload);


      /**
       * Copying is disallowed because copying a DisjointSet instance
       * violates assumptions that each DisjointSet instance uniquely
       * represents a tree in a forest of disjoint sets.  Furthermore,
       * if you copy a DisjointSet instance, and then the original
       * instance goes out of scope, internal bookkeeping gets messed
       * up.
       *
       * @param other This function is disallowed.
       */
      DisjointSet(const DisjointSet<Type>& other) = delete;


      /**
       * Moving is disallowed (see copy constructor documentation).
       *
       * @param other This function is disallowed.
       */
      DisjointSet(const DisjointSet<Type>&& other) = delete;


      /**
       * Destructor.
       */
      virtual
      ~DisjointSet();


      /**
       * Copying is disallowed (see copy constructor documentation).
       *
       * @param other This function is disallowed.
       * @return Nothing, because this function is disallowed.
       */
      DisjointSet<Type>&
      operator=(const DisjointSet<Type>& other) = delete;


      /**
       * Moving is disallowed (see copy constructor documentation).
       *
       * @param other This function is disallowed.
       * @return Nothing, because this function is disallowed.
       */
      DisjointSet<Type>&
      operator=(const DisjointSet<Type>&& other) = delete;


      /**
       * This member function returns a reference to the head of the
       * set to which *this belongs.  Although which member of a set
       * gets to be "head" is a little arbitrary, all members of the
       * set will report the same head (until the set is merged with
       * another set, when the head may change).  This
       * implementation does not use recursion.  This makes it less
       * compact than a recursive implementation, but possibly
       * faster.
       *
       * @return The return value is the head of the set.
       */
      DisjointSet&
      find();


      /**
       * This member function implements a crippled version of find
       * that doesn't fully update parent pointers.  It follows the
       * implementation in Pedro Felzenszwalb's segmentation example
       * code[2].  It turns out to be a little faster for that
       * particular application.
       *
       * [2] Efficient Graph-Based Image Segmentation.  Pedro
       * F. Felzenszwalb and Daniel P. Huttenlocher. International
       * Journal of Computer Vision, Volume 59, Number 2, September
       * 2004.
       *
       * @return The return value is a reference to the head of the
       * tree to which *this belongs.
       */
      DisjointSet&
      findNoUpdate();


      /**
       * This member function returns a reference to the head of the
       * set to which *this belongs.  It is equivalent to member
       * function find().  Although which member of a set gets to be
       * "head" is a little arbitrary, all members of the set will
       * report the same head (until the set is merged with another
       * set, when the head may change).  This implementation is
       * taken directly from the Wikipedia page on disjoint sets,
       * and uses recursion.
       *
       * @return The return value is the head of the set.
       */
      DisjointSet&
      findRecursive();


      /**
       * This member function merges two sets.  It doesn't matter
       * whether you call x.merge(y), or y.merge(x); both will get
       * the job done.  After merging, all members of both sets will
       * report the same head, which will be the head reported by
       * one of the two sets before merging.  The head reported by
       * the other of the two sets before merging will stop being a
       * head.
       *
       * @param other This argument is the set with which to merge.
       */
      void
      merge(DisjointSet& other);


      /**
       * This member function returns the payload associated with
       * *this.  The payload has no affect on the functioning of
       * *this.  It's just along for the ride.
       *
       * @return The return value is a const reference to the copy
       * that was made when the payload was specified.
       */
      const Type&
      getPayload() const;


      /**
       * This member function returns the number of members in the set
       * of which *this is a member.
       *
       * @return The return value is a count of the number of members
       * in the set.
       */
      size_t
      getSize();


      /**
       * This member function associates a payload vaule with *this.
       *
       * @param payload This argument is copied into *this as a
       * payload, and then ignored until you call member function
       * getPayload().
       */
      void
      setPayload(const Type& payload);

    protected:

      DisjointSet* m_parentPtr;
      size_t m_rank;
      size_t m_size;

      Type m_payload;
    };

  } // namespace computerVision

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/disjointSet_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_DISJOINTSET_HH */
