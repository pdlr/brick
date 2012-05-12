/**
***************************************************************************
* @file brick/computerVision/featureAssociation_impl.hh
*
* Header file defining inline and template functions declared in
* featureAssociation.hh.
*
* Copyright (C) 2008,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_FEATUREASSOCIATION_IMPL_HH
#define BRICK_COMPUTERVISION_FEATUREASSOCIATION_IMPL_HH

#include <brick/numeric/array2D.hh>

// This file is included by featureAssociation.hh, and should not be
// directly included by user code, so no need to include
// featureAssociation.hh here.
// 
// #include <brick/computerVision/featureAssociation.hh>

#include <cmath>
#include <brick/linearAlgebra/linearAlgebra.hh>
#include <brick/numeric/maxRecorder.hh>
#include <brick/numeric/utilities.hh>


namespace brick {

  namespace computerVision {


    // This function template implements the feature association
    // algorithm of Guy Scott and H. Christopher Longuet-Higgins.
    template<class FloatType, class Iterator0, class Iterator1, class Functor>
    std::vector< std::pair<size_t, size_t> >
    associateFeaturesScott91(Iterator0 sequence0Begin, Iterator0 sequence0End,
                             Iterator1 sequence1Begin, Iterator1 sequence1End,
                             Functor similarityFunctor)
    {
      // Count the number of features in each sequence.  Even if this
      // operation is O(n), it will be dominated by the more expensive
      // steps below.
      size_t sequence0Length = sequence0End - sequence0Begin;
      size_t sequence1Length = sequence1End - sequence1Begin;

      // Compute the similarity matrix.
      brick::numeric::Array2D<double> GMatrix(sequence0Length, sequence1Length);
      Iterator0 begin0 = sequence0Begin;
      for(size_t rr = 0; rr < sequence0Length; ++rr) {
        Iterator1 begin1 = sequence1Begin;
        for(size_t cc = 0; cc < sequence1Length; ++cc) {
          GMatrix(rr, cc) = similarityFunctor(*begin0, *begin1);
          ++begin1;
        }
        ++begin0;
      }

      // Compute a new similarity matrix, comprising only rotations
      // and reflections, that is "similar" to G (see the paper for
      // more details).
      brick::numeric::Array2D<double> uMatrix;
      brick::numeric::Array1D<double> sigmaArray;
      brick::numeric::Array2D<double> vTransposeMatrix;
      brick::linearAlgebra::singularValueDecomposition(
        GMatrix, uMatrix, sigmaArray, vTransposeMatrix);
      brick::numeric::Array2D<double> PMatrix =
        brick::numeric::matrixMultiply<FloatType>(uMatrix, vTransposeMatrix);

      // Find the max elements for each row and column of P.
      std::vector< brick::numeric::MaxRecorder<double, size_t> > rowMaxes(
        sequence0Length);
      std::vector< brick::numeric::MaxRecorder<double, size_t> > columnMaxes(
        sequence1Length);
      for(size_t rr = 0; rr < sequence0Length; ++rr) {
        for(size_t cc = 0; cc < sequence1Length; ++cc) {
          double similarity = PMatrix(rr, cc);
          rowMaxes[rr].test(similarity, cc);
          columnMaxes[cc].test(similarity, rr);
        }
      }

      // Valid correspondences are those for which the similarity is
      // max of both row and column.
      std::vector< std::pair<size_t, size_t> > result;
      for(size_t rr = 0; rr < sequence0Length; ++rr) {
        size_t bestColumn = rowMaxes[rr].getPayload();
        if(columnMaxes[bestColumn].getPayload() == rr) {
          result.push_back(std::make_pair(rr, bestColumn));
        }
      }

      return result;
    }

  } // namespace computerVision
    
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_FEATUREASSOCIATION_IMPL_HH */
