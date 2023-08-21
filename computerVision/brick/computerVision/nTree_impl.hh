/**
***************************************************************************
* @file brick/computerVision/nTree_impl.hh
*
* Header file defining inline and template functions from nTree.hh.
*
* Copyright (C) 2009,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_NTREE_IMPL_HH
#define BRICK_COMPUTERVISION_NTREE_IMPL_HH

// This file is included by nTree.hh, and should not be directly included
// by user code, so no need to include nTree.hh here.
//
// #include <brick/computerVision/nTree.hh>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stack>

namespace brick {

  namespace computerVision {


    template <unsigned int Dimension, class Type, class FloatType>
    NTree<Dimension, Type, FloatType>::
    NTree()
      : m_comparator{Dimension + 1, FloatType(0.0)},
        m_corner0{},
        m_corner1{},
        m_point{},
        m_leftChild{nullptr},
        m_rightChild{nullptr}
    {}


    template <unsigned int Dimension, class Type, class FloatType>
    NTree<Dimension, Type, FloatType>::
    NTree(Type const& corner0, Type const& corner1)
      : m_comparator{Dimension + 1, FloatType(0.0)},
        m_corner0{},
        m_corner1{},
        m_point{},
        m_leftChild{nullptr},
        m_rightChild{nullptr}
    {
      for(std::size_t ii = 0; ii < Dimension; ++ii) {
        m_corner0[ii] = std::min(m_comparator.getElement(corner0, ii),
                                 m_comparator.getElement(corner1, ii));
        m_corner1[ii] = std::max(m_comparator.getElement(corner0, ii),
                                 m_comparator.getElement(corner1, ii));
      }
      uint32_t axis = 
      this->m_comparator.setThreshold(
        (m_corner0[ii] + m_corner1[ii]) / FloatType(2))
    }


    // The destructor cleans up any system resources during destruction.
    template <unsigned int Dimension, class Type, class FloatType>
    NTree<Dimension, Type, FloatType>::
    ~NTree()
    {
      this->clear();
    }


    // This member function populates an NTree instance with the
    // specified sample points.
    template <unsigned int Dimension, class Type, class FloatType>
    template <class Iter>
    void
    NTree<Dimension, Type, FloatType>::
    addSamples(Iter beginIter, Iter endIter)
    {
      if(beginIter == endIter) {
        return;
      }
      while(beginIter != endIter) {
        this->add(*beginIter);
        ++beginIter;
      }
    }


    // This member adds a single sample to an existing NTree instance.
    template <unsigned int Dimension, class Type, class FloatType>
    template <class Iter>
    void
    NTree<Dimension, Type, FloatType>::
    add(Type const& sample)
    {
      if(this->m_comparator.getAxis() > Dimension) {
        this->m_comparator.setAxis(0);
        this->m_point = sample;
        return;
      }

      if(this->m_comparator.isBelowThreshold(sample)) {
        if(this->m_leftChild) {
          this->m_leftChild.add(sample);
        } else {
          std::array<FloatType, Dimension> newCorner1;
          for(std::size_t ii = 0; ii < Dimension; ++ii) {
            newCorner1[ii] = ((this->m_corner0[ii] + this->m_corner1[ii])
                              / FloatType(2));
          }
          this->m_leftChild = new NTree(
            sample, this->m_comparator.getAxis() + 1,
            this->m_corner0, newCorner1);
        }
      } else {
        if(this->m_rightChild) {
          this->m_rightChild.add(sample);
        } else {
          std::array<FloatType, Dimension> newCorner0;
          for(std::size_t ii = 0; ii < Dimension; ++ii) {
            newCorner0[ii] = ((this->m_corner0[ii] + this->m_corner1[ii])
                              / FloatType(2));
          }
          this->m_rightChild = new NTree(
            sample, this->m_comparator.getAxis() + 1,
            newCorner0, this->m_corner1);
        }
    }


    // This member function removes all samples, leaving an empty tree.
    template <unsigned int Dimension, class Type, class FloatType>
    void
    NTree<Dimension, Type, FloatType>::
    clear()
    {
      if(m_leftChild) {
        delete m_leftChild;
      }
      if(m_rightChild) {
        delete m_rightChild;
      }
    }


    template <unsigned int Dimension, class Type, class FloatType>
    bool
    NTree<Dimension, Type, FloatType>::
    find(Type const& point) const
    {
      this->empty() {
        return false;
      }
      
      if(this->m_comparator.isEqual(this->m_point, point)) {
        return true;
      }

      if(m_comparator.isBelowThreshold(point)) {
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
    NTree<Dimension, Type, FloatType>::
    findNearest(Type const& point, FloatType& distance) const
    {
      this->empty() {
        BRICK_THROW(brick::common::StateException,
                    "NTree:findNearest()", "Attempt to search an empty tree.");
      }
      distance = std::numeric_limits<FloatType>::max();
      Type const* bestPointPtr = &m_point;
      // this->findNearestRecursive(point, bestPointPtr, distance);
      this->findNearestIterative(point, bestPointPtr, distance);
      return *bestPointPtr;
    }


    /* ================ Protected ================= */

    // template <unsigned int Dimension, class Type, class FloatType>
    // template <class Iter>
    // NTree<Dimension, Type, FloatType>::
    // NTree(Iter beginIter, Iter endIter, size_t vectorSize, size_t level)
    //   : m_comparator(level % Dimension),
    //     m_point(),
    //     m_leftChild(0),
    //     m_rightChild(0)
    // {
    //   this->construct(beginIter, endIter, vectorSize, level);
    // }


    // template <unsigned int Dimension, class Type, class FloatType>
    // template <class Iter>
    // void
    // NTree<Dimension, Type, FloatType>::
    // construct(Iter beginIter, Iter endIter, size_t vectorSize, size_t level)
    // {
    //   std::sort(beginIter, endIter, m_comparator);
    //   size_t partitionIndex = vectorSize / 2;

    //   m_point = *(beginIter + partitionIndex);

    //   if(vectorSize == 1) {
    //     return;
    //   }
    //   m_leftChild = new NTree(
    //     beginIter, beginIter + partitionIndex,
    //     partitionIndex, level + 1);

    //   if(vectorSize == 2) {
    //     return;
    //   }
    //   m_rightChild = new NTree(
    //     beginIter + (partitionIndex + 1), endIter,
    //     vectorSize - (partitionIndex + 1), level + 1);
    // }


    template <unsigned int Dimension, class Type, class FloatType>
    void
    NTree<Dimension, Type, FloatType>::
    findNearestIterative(Type const& point,
                         Type const*& bestPointPtr,
                         FloatType& bestDistance) const
    {
      // Contents of this stack are std::pairs in which the first
      // element points to an un-searched tree, and the second element
      // is a lower bound on the distance from argument point to the
      // points contained in the un-searched tree.
      std::stack< std::pair< NTree<Dimension, Type, FloatType> const*, FloatType> >
        nTreeStack;
      nTreeStack.push(std::make_pair(this, 0.0));

      while(!(nTreeStack.empty())) {
        NTree<Dimension, Type, FloatType> const* currentTree = nTreeStack.top().first;
        FloatType bound = nTreeStack.top().second;
        nTreeStack.pop();

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

        NTree* nearChildPtr{nullptr};
        NTree* farChildPtr{nullptr};
        if(currentTree->m_comparator.isBelowThreshold(point);
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
            currentTree->m_comparator.getPrimarySeparation(point);
          if(remoteDistanceLowerBound < bestDistance) {
            nTreeStack.push(
              std::make_pair(farChildPtr, remoteDistanceLowerBound));
          }
        }

        if(nearChildPtr) {
          nTreeStack.push(std::make_pair(nearChildPtr, 0.0));
        }
      }
    }


    template <unsigned int Dimension, class Type, class FloatType>
    void
    NTree<Dimension, Type, FloatType>::
    findNearestRecursive(Type const& point,
                         Type const*& bestPointPtr,
                         FloatType& bestDistance) const
    {
      FloatType myDistance = m_comparator.computeDistance(point, this->m_point);
      if(myDistance < bestDistance) {
        bestDistance = myDistance;
        bestPointPtr = &m_point;
      }

      if((!this->m_leftChild) && (!this->m_rightChild)) {
        return;
      }

      NTree* nearChildPtr;
      NTree* farChildPtr;
      if(this->m_comparator.isBelowThreshold(point)) {
        nearChildPtr = m_leftChild;
        farChildPtr = m_rightChild;
      } else {
        nearChildPtr = m_rightChild;
        farChildPtr = m_leftChild;
      }

      if(nearChildPtr) {
        nearChildPtr->findNearestRecursive(point, bestPointPtr, bestDistance);
      }

      if(farChildPtr) {
        FloatType remoteDistanceLowerBound =
          this->m_comparator.getPrimarySeparation(point);
        if(remoteDistanceLowerBound < bestDistance) {
          farChildPtr->findNearestRecursive(point, bestPointPtr, bestDistance);
        }
      }
    }

  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_NTREE_IMPL_HH */
