/**
***************************************************************************
* @file dlrComputerVision/ransac.cc
*
* Source file defining helper functions for implementing Fischler's
* and Bolles's RANSAC algorithm.
*
* Copyright (C) 2008, 2011 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
***************************************************************************
*/

#include <cmath>
#include <brick/common/exception.hh>
#include <brick/computerVision/ransac.hh>

namespace brick {

  namespace computerVision {

    // Computes how many RANSAC iterations are required to be sure
    // there's at least one iteration in which the model was estimated
    // using all inliers.
    unsigned int
    ransacGetRequiredIterations(unsigned int sampleSize,
                                double requiredConfidence,
                                double inlierProbability)
    {
      if((requiredConfidence < 0.0) || (requiredConfidence >= 1.0)) {
        BRICK_THROW(brick::common::ValueException,
                  "ransacGetRequiredIterations()",
                  "Probability value requiredConfidence is out of range.");
      }
      if((inlierProbability < 0.0) || (inlierProbability >= 1.0)) {
        BRICK_THROW(brick::common::ValueException,
                  "ransacGetRequiredIterations()",
                  "Probability value inlierProbability is out of range.");
      }

      // The probability that any particular random sample is all inliers.
      double singlePickConfidence =
        std::pow(inlierProbability, static_cast<double>(sampleSize));

      // The probability that any particular random sample contains at
      // least one outlier.
      double singlePickDisconfidence = 1.0 - singlePickConfidence;

      // Number of times we have to sample in order to get a set of
      // all inliers with probability requiredConfidence.  This
      // follows Fischler's paper exactly.
      double numberOfRandomSampleSets =
        std::ceil(std::log(1.0 - requiredConfidence)
                  / std::log(singlePickDisconfidence)) + 0.5;

      return std::max(int(1),
                      static_cast<int>(numberOfRandomSampleSets + 0.5));
    }

  } // namespace computerVision

} // namespace brick
