/**
***************************************************************************
* @file brick/computerVision/nTree.hh
*
* Header file declaring a class implementing a KD-Tree data structure.
*
* Copyright (C) 2009,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_NTREE_HH
#define BRICK_COMPUTERVISION_NTREE_HH

#include <cstdlib>
#include <functional>
#include <vector>

// For KDComparator.
#include <brick/computerVision/kdTree.hh>

namespace brick {

  namespace computerVision {

    /**
     ** This class template is used by the NTree class template to
     ** interact with the data points to be stored in the N-Tree.  It
     ** provides functions for comparing points during a sort, for
     ** computing the distance between points, and for finding the
     ** lower bound on the distance between a particular point and the
     ** region of space represented by a particular branch of the
     ** N-Tree.  Template argument Dimension specifies how many
     ** dimension the tree spans.  Template argument Type specifies
     ** what kind of element is contained in the tree.  Logically,
     ** Type represents a point in multi-dimensional space.  In order
     ** to work with this class template, Type must allow access to
     ** individual coordinates via Type::operator[](size_t), which
     ** will be called with arguments in the range [0 ... (Dimension -
     ** 1)].  Also, operator==(Type const&, Type const&) must return
     ** true if the two arguments represent the same point.
     **
     ** If you need to build a N-Tree using a Type that doesn't
     ** support this interface, you can simply specialize NComparator
     ** (or specific member functions of NComparator) for your type.
     ** For an example of such a specialization, see the file
     ** brick/computerVision/test/nTreeTest.cc.
     **/
    template <unsigned int Dimension, class Type, class FloatType = double>
    class NComparator
      : public std::binary_function<Type, Type, bool>
    {
      public:

      /**
       * Constructor.  Each NComparator is associated with a specific
       * level of the N-Tree, and therefor with a specific axis in
       * multi-dimensional space.  This constructor specifies whith
       * which axis *this is associated.
       *
       * @param axis This argument must be an integer between 0 and
       * (Dimension - 1), inclusive.
       */
      NComparator(unsigned int axis = 0)
        : m_axis(axis) {};


      /**
       * Destructor.
       */
      virtual
      ~NComparator() {}


      /**
       * This member function computes the square of the Euclidean
       * distance between two Type instances.
       *
       * @param arg0 This argument is the first of the two Type
       * instances.
       *
       * @param arg1 This argument is the second of the two Type
       * instances.
       *
       * @return The return value is the Euclidean distance between
       * arg0 and arg1.
       */
      FloatType
      computeDistance(Type const& arg0, Type const& arg1) const {
        FloatType distance = 0.0;
        for(size_t ii = 0; ii < Dimension; ++ii) {
          FloatType newTerm = arg0[ii] - arg1[ii];
          distance += newTerm * newTerm;
        }
        return distance;
      }


      unsigned int
      getAxis {
        return this->m_axis;
      }


      FloatType
      getElement(Type const& arg0, std::size_t ii) {
        return arg0[ii];
      }

      /**
       * This member function computes a lower bound on the distance
       * between the specified Type instance and Type instances
       * contained in the "far" branch of a N-Tree.  The dividing
       * hyperplane between the two brances is assumed to run
       * perpendicular to the axis specified by this->m_axis, and
       * intersect the point specified by argument arg1.
       *
       * @param arg0 This argument is the point for which to compute
       * the lower bound distance.
       *
       * @param arg1 This argument is the point through which the
       * separating hyperplane passes.
       *
       * @return The return value is the minimum distance between arg0
       * and the separating hyperplane.
       */
      FloatType
      getPrimarySeparation(Type const& arg0, Type const& arg1) const {
        FloatType difference = arg0[m_axis] - arg1[m_axis];
        return difference * difference;
      }


      /**
       * This member function returns true if its two arguments
       * represent the same point.
       *
       * @param arg0 This argument is the first Type instance to be
       * compared.
       *
       * @param arg1 This argument is the second Type instance to be
       * compared.
       *
       * @return The return value is true if the two arguments
       * represent the same point, false otherwise.
       */
      bool isEqual(Type const& arg0, Type const& arg1) const {
        return arg0 == arg1;
      }


      /**
       * This operator implements "arg0 < arg1" so that sequences of
       * Type instances can be sorted using std::sort().  NComparator
       * instances will be passed as the final (functor) argument to
       * std::sort().  For normal N-Tree operation, it should begin
       * by comparing the locations of arg1 and arg0 along the axis
       * associated with *this.  In the current implementation, Type
       * instances that have identical coordinates on this one axis
       * will be compared along other axes.
       *
       * @return The return value is true if arg0 is "less than" arg1,
       * false otherwise.
       */
      bool
      operator()(Type const& arg0, Type const& arg1) const {
        if(arg0[m_axis] < arg1[m_axis]) {
          return true;
        }
        if(arg0[m_axis] > arg1[m_axis]) {
          return false;
        }
        for(size_t ii = m_axis + 1; ii < Dimension; ++ii) {
          if(arg0[ii] < arg1[ii]) {
            return true;
          }
          if(arg0[ii] > arg1[ii]) {
            return false;
          }
        }
        for(size_t ii = 0; ii < m_axis; ++ii) {
          if(arg0[ii] < arg1[ii]) {
            return true;
          }
          if(arg0[ii] > arg1[ii]) {
            return false;
          }
        }
        return false;
      }


    private:

      unsigned int m_axis;

    };
 

    /**
     ** This class implements a modified KD-Tree-like data structure that
     ** partitions a bounded space, with a preset partitioning scheme,
     ** rather than choosing its partitions based on input data.  One
     ** disadvantage of this approach is that the partitioning doesn't
     ** match the data as well.  If your data are concentrated in a small
     ** part of the space, you can wind up with a very unbalanced tree.
     ** It has the advantage, though, that the order of point presentation
     ** is much less important, so you can add points one-by-one without
     ** compromising performance.  A second disadvantage is that there's
     ** bookkeeping information to track, so memory footprint is bigger.
     **
     ** Template argument Dimension specifies how many dimensions the
     ** tree will span.  Template argument Type specifies what kind of
     ** element will be contained in the tree.  Logically, Type
     ** represents a point in multi-dimensional space.  It must
     ** support default construction, copying and assignment.  It must
     ** also fulfill the requirements of the NComparator class
     ** template.  Note that if you need to use a Type that doesn't
     ** support the NComparator requirements, you can always
     ** specialize NComparator for your specific Type.
     **
     ** Here's an example of how to use the NTree class template:
     **
     ** @code
     **   namespce cv = brick::computerVision;
     **   namespce num = brick::numeric;
     **
     **   std::vector<num::Vector3D> myPoints;
     **   myPoints.push_back(num::Vector3D(2.1, 3.5, 2.0));
     **   myPoints.push_back(num::Vector3D(5.0, 3.2, 2.5));
     **   myPoints.push_back(num::Vector3D(2.4, 1.6, 1.3));
     **   myPoints.push_back(num::Vector3D(7.7, 4.7, 1.1));
     **   myPoints.push_back(num::Vector3D(-2.0, 6.3, 5.0));
     **   myPoints.push_back(num::Vector3D(0.0, 0.0, 0.0));
     **   myPoints.push_back(num::Vector3D(3.1, 4.7, 5.4));
     **
     **   // If you're doing this, you may as well use a KDTree.
     **   // cv::NTree<3, num::Vector3D> nTree(
     **   //     myPoints.begin(), myPoints.end());
     **   cv::NTree<3, num::Vector3D> nTree(
     **          num::Vector3D(0.0, 0.0, 0.0),
     **          num::Vector3D(10.0, 10.0, 10.0));
     **   for(auto point : myPoints) {
     **       nTree.add(point);
     **   }
     **
     **   FloatType distance;
     **   num::Vector3D testPoint(3.5, 6.9, 4.4);
     **   num::Vector3D nearestPoint = nTree.findNearest(testPoint, distance);
     **
     **   std::cout << "The closest point was " << nearestPoint << ", "
     **             << "which was " << distance << " distance from "
     **             << testPoint << std::endl;
     ** @endcode
     **/
    template <unsigned int Dimension, class Type, class FloatType = double>
    class NTree {
    public:

      /**
       * The default constructor creates an empty tree.
       */
      NTree();


      /**
       * This constructor is like calling the default constructor, followed
       * by the setBounds() member function.
       */
      NTree(Type const& corner0, Type const& corner1);


      /**
       * This constructor creates a tree and populates it with the
       * specified sample points.  It has complexity O(N*log(N)),
       * where N is the number of elements to be inserted into the
       * tree.
       *
       * @param beginIter This argument is an iterator pointing to the
       * beginning of a sequence of Type instances that is to be
       * inserted into the tree.
       *
       * @param endIter This argument is an interator pointing one
       * element past the last Type instance in the sequence that is
       * to be inserted into the tree.
       */
      template <class Iter>
      NTree(Iter beginIter, Iter endIter,
            Type const& corner0, Type const& corner1);


      /**
       * The destructor cleans up any system resources during
       * destruction.
       */
      virtual
      ~NTree();


      /**
       * This member function adds a single point to the tree without first
       * clearing it.  It has worst-case complexity O(N), but
       * for non-pathological datasets performance is more like O(log(N)),
       * where N is the number of elements in the tree.
       *
       * @param sample This argument is the point to be added.  It 
       * must lie within the bounds specified by constructor or setBounds.
       */
      void
      add(Type const& sample);


      /**
       * This member function adds points to the tree without first
       * clearing it.  It has worst-case complexity O(N * (N + M)), but
       * for non-pathological datasets performance is more like
       * O(N*log(N + M)), where N is the number of elements to be
       * inserted into the tree, and M is the number of elements already
       * in the tree.
       *
       * @param beginIter This argument is an iterator pointing to the
       * beginning of a sequence of Type instances that is to be
       * inserted into the tree.
       *
       * @param endIter This argument is an interator pointing one
       * element past the last Type instance in the sequence that is
       * to be inserted into the tree.
       */
      template <class Iter>
      void
      addSamples(Iter beginIter, Iter endIter);


      /**
       * This member function removes all samples, leaving an empty tree.
       */
      void
      clear();


      /**
       * This member function returns true if the specified Type
       * instance has already been inserted into the tree.  It has
       * complexity O(log(N), where N is the number of points
       * contained in the tree.
       *
       * @param point This argument is the Type instance to search
       * for.  It will be compared to elements in the tree using
       * NComparator<Dimension, Type>::isEqual().
       *
       * @return The return value is true if a matching point is found
       * in the tree, false otherwise.
       */
      bool
      find(Type const& point) const;


      /**
       * This member function returns a const reference to the tree
       * element that is closest (Euclidean distance) to the specified
       * point.  It has complexity O(log(N), where N is the number of
       * points contained in the tree.
       *
       * @param point This argument is the Type instance to search
       * for.  It will be compared to elements in the tree using
       * NComparator<Dimension, Type>::computeDistance(Type const&,
       * Type const&).
       *
       * @param distance This argument is used to return the distance
       * between the point for which we're searching and the closest
       * point in the tree.  It will be computed using
       * NComparator<Dimension, Type>::computeDistance(Type const&,
       * Type const&).
       *
       * @return The return value is a const reference to the closest
       * point in the tree.
       */
      Type const&
      findNearest(Type const& point, FloatType& distance) const;


      // void
      // rebalance();


    protected:

      template <class Iter>
      NTree(Iter beginIter, Iter endIter, size_t vectorSize,
                    size_t level);


      template <class Iter>
      void
      construct(Iter beginIter, Iter endIter, size_t vectorSize, size_t level);

      void
      findNearestIterative(Type const& point,
                           Type const*& bestPointPtr,
                           FloatType& bestDistance) const;


      void
      findNearestRecursive(Type const& point,
                           Type const*& bestPointPtr,
                           FloatType& bestDistance) const;


      NComparator<Dimension, Type> m_comparator;
      std::array<FloatType, Dimension> m_corner0;
      std::array<FloatType, Dimension> m_corner1;
      Type m_point;
      NTree* m_leftChild;
      NTree* m_rightChild;
    };

  } // namespace computerVision

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/nTree_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_NTREE_HH */
