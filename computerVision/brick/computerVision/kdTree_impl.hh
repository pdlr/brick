/**
***************************************************************************
* @file brick/computerVision/kdTree_impl.hh
*
* Header file defining inline and template functions from kdTree.hh.
*
* Copyright (C) 2009,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_KDTREE_IMPL_HH
#define BRICK_COMPUTERVISION_KDTREE_IMPL_HH

// This file is included by kdTree.hh, and should not be directly included
// by user code, so no need to include kdTree.hh here.
//
// #include <brick/computerVision/kdTree.hh>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stack>

namespace brick {

  namespace computerVision {


    template <unsigned int Dimension, class Type, class FloatType>
    KDTree<Dimension, Type, FloatType>::
    KDTree()
      : m_comparator(0),
        m_point(),
        m_leftChild(0),
        m_rightChild(0)
    {}


    template <unsigned int Dimension, class Type, class FloatType>
    template <class Iter>
    KDTree<Dimension, Type, FloatType>::
    KDTree(Iter beginIter, Iter endIter)
      : m_comparator(0),
        m_point(),
        m_leftChild(0),
        m_rightChild(0)
    {
      this->addSamples(beginIter, endIter);
    }


    // The destructor cleans up any system resources during destruction.
    template <unsigned int Dimension, class Type, class FloatType>
    KDTree<Dimension, Type, FloatType>::
    ~KDTree()
    {
      this->clear();
    }


    // This member function populate a KDTree instance with the
    // specified sample points.
    template <unsigned int Dimension, class Type, class FloatType>
    template <class Iter>
    void
    KDTree<Dimension, Type, FloatType>::
    addSamples(Iter beginIter, Iter endIter)
    {
      // It would be nice to do away with this, so that points can be
      // added in stages, if desired.
      this->clear();

      if(beginIter == endIter) {
        return;
      }

      std::vector<Type> pointVector;
      std::copy(beginIter, endIter, std::back_inserter(pointVector));

      this->construct(pointVector.begin(), pointVector.end(),
                      pointVector.size(), 0);
    }


    // This member function removes all samples, leaving an empty tree.
    template <unsigned int Dimension, class Type, class FloatType>
    void
    KDTree<Dimension, Type, FloatType>::
    clear()
    {
      if(m_leftChild != 0) {
        delete m_leftChild;
      }
      if(m_rightChild != 0) {
        delete m_rightChild;
      }
    }


    template <unsigned int Dimension, class Type, class FloatType>
    bool
    KDTree<Dimension, Type, FloatType>::
    find(Type const& point) const
    {
      if(m_comparator.isEqual(m_point, point)) {
        return true;
      }

      if(m_comparator(point, m_point)) {
        if(m_leftChild == 0) {
          return false;
        }
        return m_leftChild->find(point);
      } else {
        if(m_rightChild == 0) {
          return false;
        }
        return m_rightChild->find(point);
      }
    }


    template <unsigned int Dimension, class Type, class FloatType>
    Type const&
    KDTree<Dimension, Type, FloatType>::
    findNearest(Type const& point, FloatType& distance) const
    {
      distance = std::numeric_limits<FloatType>::max();
      Type const* bestPointPtr = &m_point;
      // this->findNearestRecursive(point, bestPointPtr, distance);
      this->findNearestIterative(point, bestPointPtr, distance);
      return *bestPointPtr;
    }


    /* ================ Protected ================= */

    template <unsigned int Dimension, class Type, class FloatType>
    template <class Iter>
    KDTree<Dimension, Type, FloatType>::
    KDTree(Iter beginIter, Iter endIter, size_t vectorSize, size_t level)
      : m_comparator(level % Dimension),
        m_point(),
        m_leftChild(0),
        m_rightChild(0)
    {
      this->construct(beginIter, endIter, vectorSize, level);
    }


    template <unsigned int Dimension, class Type, class FloatType>
    template <class Iter>
    void
    KDTree<Dimension, Type, FloatType>::
    construct(Iter beginIter, Iter endIter, size_t vectorSize, size_t level)
    {
      std::sort(beginIter, endIter, m_comparator);
      size_t partitionIndex = vectorSize / 2;

      m_point = *(beginIter + partitionIndex);

      if(vectorSize == 1) {
        return;
      }
      m_leftChild = new KDTree(
        beginIter, beginIter + partitionIndex,
        partitionIndex, level + 1);

      if(vectorSize == 2) {
        return;
      }
      m_rightChild = new KDTree(
        beginIter + (partitionIndex + 1), endIter,
        vectorSize - (partitionIndex + 1), level + 1);
    }


    template <unsigned int Dimension, class Type, class FloatType>
    void
    KDTree<Dimension, Type, FloatType>::
    findNearestIterative(Type const& point,
                         Type const*& bestPointPtr,
                         FloatType& bestDistance) const
    {
      // Contents of this stack are std::pairs in which the first
      // element points to an un-searched tree, and the second element
      // is a lower bound on the distance from argument point to the
      // points contained in the un-searched tree.
      std::stack< std::pair< KDTree<Dimension, Type, FloatType> const*, FloatType> >
        kdTreeStack;
      kdTreeStack.push(std::make_pair(this, 0.0));

      while(!(kdTreeStack.empty())) {
        KDTree<Dimension, Type, FloatType> const* currentTree = kdTreeStack.top().first;
        FloatType bound = kdTreeStack.top().second;
        kdTreeStack.pop();

        if(bestDistance < bound) {
          continue;
        }

        FloatType myDistance = currentTree->m_comparator.computeDistance(
          point, currentTree->m_point);
        if(myDistance < bestDistance) {
          bestDistance = myDistance;
          bestPointPtr = &currentTree->m_point;
        }

        if(currentTree->m_leftChild == 0 && currentTree->m_rightChild == 0) {
          continue;
        }

        KDTree* nearChildPtr;
        KDTree* farChildPtr;
        bool isLeft = currentTree->m_comparator(point, currentTree->m_point);
        if(isLeft) {
          nearChildPtr = currentTree->m_leftChild;
          farChildPtr = currentTree->m_rightChild;
        } else {
          nearChildPtr = currentTree->m_rightChild;
          farChildPtr = currentTree->m_leftChild;
        }

        // Push remote child first, and near child second, so that
        // near child will be popped first, increasing the chance that
        // remote child will be eliminated.

        // Only push the remote child if it's plausible that it
        // contains a closer point than our best so far.  This
        // duplicates a similar test above. Not sure if it's worth
        // duplicating the test to avoid the occasional extra
        // push/pop.
        if(farChildPtr) {
          FloatType remoteDistanceLowerBound =
            currentTree->m_comparator.getPrimarySeparation(
              point, currentTree->m_point);
          if(remoteDistanceLowerBound < bestDistance) {
            kdTreeStack.push(
              std::make_pair(farChildPtr, remoteDistanceLowerBound));
          }
        }

        if(nearChildPtr) {
          kdTreeStack.push(std::make_pair(nearChildPtr, 0.0));
        }
      }
    }


    template <unsigned int Dimension, class Type, class FloatType>
    void
    KDTree<Dimension, Type, FloatType>::
    findNearestRecursive(Type const& point,
                         Type const*& bestPointPtr,
                         FloatType& bestDistance) const
    {
      FloatType myDistance = m_comparator.computeDistance(point, m_point);
      if(myDistance < bestDistance) {
        bestDistance = myDistance;
        bestPointPtr = &m_point;
      }

      if(m_leftChild == 0 && m_rightChild == 0) {
        return;
      }

      bool isLeft = m_comparator(point, m_point);
      FloatType remoteDistanceLowerBound = m_comparator.getPrimarySeparation(
        point, m_point);

      KDTree* nearChildPtr;
      KDTree* farChildPtr;
      if(isLeft) {
        nearChildPtr = m_leftChild;
        farChildPtr = m_rightChild;
      } else {
        nearChildPtr = m_rightChild;
        farChildPtr = m_leftChild;
      }

      if(nearChildPtr) {
        nearChildPtr->findNearestRecursive(point, bestPointPtr, bestDistance);
      }

      if((remoteDistanceLowerBound < bestDistance) && farChildPtr) {
        farChildPtr->findNearestRecursive(point, bestPointPtr, bestDistance);
      }
    }

  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_KDTREE_IMPL_HH */
