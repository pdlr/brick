/**
***************************************************************************
* @file brick/computerVision/kdTree.hh
*
* Header file declaring a class implementing a KD-Tree data structure.
*
* Copyright (C) 2009,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_KDTREE_HH
#define BRICK_COMPUTERVISION_KDTREE_HH

#include <cstdlib>
#include <functional>
#include <vector>


namespace brick {

  namespace computerVision {

    /**
     ** This class template is used by the KDTree class template to
     ** interact with the data points to be stored in the KD-Tree.  It
     ** provides functions for comparing points during a sort, for
     ** computing the distance between points, and for finding the
     ** lower bound on the distance between a particular point and the
     ** region of space represented by a particular branch of the
     ** KD-Tree.  Template argument Dimension specifies how many
     ** dimension the tree spans.  Template argument Type specifies
     ** what kind of element is contained in the tree.  Logically,
     ** Type represents a point in multi-dimensional space.  In order
     ** to work with this class template, Type must allow access to
     ** individual coordinates via Type::operator[](size_t), which
     ** will be called with arguments in the range [0 ... (Dimension -
     ** 1)].  Also, operator==(Type const&, Type const&) must return
     ** true if the two arguments represent the same point.
     **
     ** If you need to build a KD-Tree using a Type that doesn't
     ** support this interface, you can simply specialize KDComparator
     ** (or specific member functions of KDComparator) for your type.
     ** For an example of such a specialization, see the file
     ** brick/computerVision/test/kdTreeTest.cc.
     **/
    template <unsigned int Dimension, class Type, class FloatType = double>
    class KDComparator
      : public std::binary_function<Type, Type, bool>
    {
      public:

      /**
       * Constructor.  Each KDComparator is associated with a specific
       * level of the KD-Tree, and therefor with a specific axis in
       * multi-dimensional space.  This constructor specifies whith
       * which axis *this is associated.
       *
       * @param axis This argument must be an integer between 0 and
       * (Dimension - 1), inclusive.
       */
      KDComparator(unsigned int axis = 0)
        : m_axis(axis) {};


      /**
       * Destructor.
       */
      virtual
      ~KDComparator() {}


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


      /**
       * This member function computes a lower bound on the distance
       * between the specified Type instance and Type instances
       * contained in the "far" branch of a KD-Tree.  The dividing
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
       * Type instances can be sorted using std::sort().  KDComparator
       * instances will be passed as the final (functor) argument to
       * std::sort().  For normal KD-Tree operation, it should begin
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
     ** This class implements a basic KD-Tree data structure.
     ** Template argument Dimension specifies how many dimensions the
     ** tree will span.  Template argument Type specifies what kind of
     ** element will be contained in the tree.  Logically, Type
     ** represents a point in multi-dimensional space.  It must
     ** support default construction, copying and assignment.  It must
     ** also fulfill the requirements of the KDComparator class
     ** template.  Note that if you need to use a Type that doesn't
     ** support the KDComparator requirements, you can always
     ** specialize KDComparator for your specific Type.  This class
     ** currently does not support adding points after construction or
     ** rebalancing.
     **
     ** Here's an example of how to use the KDTree class template:
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
     **   cv::KDTree<3, num::Vector3D> kdTree(myPoints.begin(), myPoints.end());
     **
     **   FloatType distance;
     **   num::Vector3D testPoint(3.5, 6.9, 4.4);
     **   num::Vector3D nearestPoint = kdTree.findNearest(testPoint, distance);
     **
     **   std::cout << "The closest point was " << nearestPoint << ", "
     **             << "which was " << distance << " distance from "
     **             << testPoint << std::endl;
     ** @endcode
     **/
    template <unsigned int Dimension, class Type, class FloatType = double>
    class KDTree {
    public:

      /**
       * The default constructor creates an empty tree.
       */
      KDTree();


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
      KDTree(Iter beginIter, Iter endIter);


      /**
       * The destructor cleans up any system resources during
       * destruction.
       */
      virtual
      ~KDTree();


      /**
       * This member function clears a KDTree instance, and then
       * populates it with the specified sample points.  It has
       * complexity O(N*log(N)), where N is the number of elements to
       * be inserted into the tree.
       *
       * Future implementations may add points to the existing tree,
       * without clearing the tree first, so please explicitly
       * call the clear() member function prior to adding samples.
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
       * KDComparator<Dimension, Type>::isEqual().
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
       * KDComparator<Dimension, Type>::computeDistance(Type const&,
       * Type const&).
       *
       * @param distance This argument is used to return the distance
       * between the point for which we're searching and the closest
       * point in the tree.  It will be computed using
       * KDComparator<Dimension, Type>::computeDistance(Type const&,
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
      KDTree(Iter beginIter, Iter endIter, size_t vectorSize, size_t level);


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


      KDComparator<Dimension, Type> m_comparator;
      Type m_point;
      KDTree* m_leftChild;
      KDTree* m_rightChild;
    };

  } // namespace computerVision

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/kdTree_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_KDTREE_HH */
