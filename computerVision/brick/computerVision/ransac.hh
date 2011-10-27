/**
***************************************************************************
* @file brick/computerVision/ransac.hh
*
* Header file declaring helper functions for implementing Fischler's and
* Bolles's RANSAC algorithm.
*
* Copyright (C) 2008, 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
***************************************************************************
*/

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <brick/common/exception.hh>
#include <brick/numeric/array2D.hh>
#include <brick/random/pseudoRandom.hh>

namespace brick {

  namespace computerVision {

    template <class InIter, class OutIter, class Functor>
    void
    ransacGetConsensusSet(
      InIter inBegin, InIter inEnd, OutIter outBegin, Functor functor);

    
    template <class Type, class Functor>
    brick::numeric::Array2D<Type>
    ransacGetConsensusSetRows(
      brick::numeric::Array2D<Type> const& candidates,
      Functor& functor);

    
    unsigned int
    ransacGetRequiredIterations(unsigned int sampleSize,
                                double requiredConfidence,
                                double inlierProbability);
    

    template <class Type>
    brick::numeric::Array2D<Type>
    ransacSelectRows(brick::numeric::Array2D<Type> const& sampleArray,
                     unsigned int numberOfSamplesRequired);
    

    template <class Type>
    brick::numeric::Array2D<Type>
    ransacSelectRows(brick::numeric::Array2D<Type> const& sampleArray,
                     unsigned int numberOfSamplesRequired,
                     brick::common::Int64& seed);

    
    template <class Type>
    brick::numeric::Array1D<Type>
    ransacSelectRows(brick::numeric::Array2D<Type> const& sampleArray,
                     unsigned int numberOfSamplesRequired,
                     brick::random::PseudoRandom& pseudoRandom);


    
  } // namespace computerVision
  
} // namespace brick

