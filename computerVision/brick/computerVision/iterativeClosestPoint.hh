/**
***************************************************************************
* @file brick/computerVision/iterativeClosestPoint.hh
*
* Header file declaring a class template implementing a derivative of
* Besl's and McKay's Iterative Closest Point algorithm.
*
* Copyright (C) 2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_ITERATIVECLOSESTPOINT_HH
#define BRICK_COMPUTERVISION_ITERATIVECLOSESTPOINT_HH

#include<brick/computerVision/kdTree.hh>
#include<brick/numeric/transform3D.hh>

namespace brick {

  namespace computerVision {

    /**
     ** This class implements a basic ICP algorithm of Besl and McKay
     ** [1].  Template argument Dimension specifies how many
     ** dimensions the point sets span (for example, set this to 3 for
     ** 3D points.  Template argument Type specifies what kind of
     ** element will make up the point sets (for example,
     ** Vector2D<double>).  Logically, Type represents a point in
     ** multi-dimensional space.  It must support default
     ** construction, copying and assignment.  It must also fulfill
     ** the requirements of the KDComparator class template (see
     ** brick/computerVision/kdTree.hh).  Note that if you need to use
     ** a Type that doesn't support the KDComparator requirements, you
     ** can always specialize KDComparator for your specific Type.
     **
     ** Here's an example of how to use the ICP class template:
     **
     ** @code
     ** @endcode
     **
     ** [1] P. J. Besl and N. D. McKay, A Method for Registration of
     ** 3-D Shapes, IEEE Transactions on Pattern Analysis and Machine
     ** Intelligence, Vol 14(2), pp 239-256, February, 1992.
     **/
    template <unsigned int Dimension, class Type, class FloatType = double>
    class IterativeClosestPoint {
    public:

      /** 
       * The default constructor.
       */
      IterativeClosestPoint();


      /**
       * The destructor cleans up any system resources.
       */
      virtual
      ~IterativeClosestPoint();


      brick::numeric::Transform3D<FloatType>
      getTransform();


      template <class Iter>
      brick::numeric::Transform3D<FloatType>
      registerPoints(
        Iter beginIter, Iter endIter,
        brick::numeric::Transform3D<FloatType> const& modelFromQueryEstimate 
        = brick::numeric::Transform3D<FloatType>());
      
      
      template <class Iter>
      void
      setModelPoints(Iter beginIter, Iter endIter);


      void
      setInitialTransform(brick::numeric::Transform3D<FloatType> const&
                          modelFromQueryEstimate);


    protected:

      brick::numeric::Transform3D<FloatType>
      estimateTransformModelFromQuery(
        std::vector<Type> const& selectedQueryPoints,
        std::vector<Type> const& matchingModelPoints,
        std::vector<FloatType> const& weights);


      bool
      findMatches(std::vector<Type>& matchingModelPoints,
                  std::vector<FloatType>& weights,
		  unsigned int& count,
		  FloatType& rmsError,
                  std::vector<Type> const& queryPoints,
                  brick::numeric::Transform3D<FloatType> const& modelFromQuery);


      void
      selectQueryPoints(std::vector<Type>& selectedQueryPoints,
                        std::vector<Type> const& allQueryPoints);


      FloatType                          m_convergenceThreshold;
      FloatType                          m_distanceThreshold;
      KDTree<Dimension, Type, FloatType> m_modelTree;

    };

  } // namespace computerVision
  
} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/iterativeClosestPoint_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_ITERATIVECLOSESTPOINT_HH */
