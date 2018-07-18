/**
***************************************************************************
* @file brick/computerVision/keypointSelectorLepetit.cc
*
* Source file defining a class for selecting stable keypoints from an
* image.
*
* Copyright (C) 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#include <brick/computerVision/keypointSelectorLepetit.hh>

namespace brick {

  namespace computerVision {

    KeypointSelectorLepetit::
    KeypointSelectorLepetit()
      : m_levelVector(),
        m_locationVector(),
        m_thresholdLaplacianMagnitude(0.0),
        m_thresholdPixelSimilarity(0.0)
    {
      // Empty.
    }


    std::vector<numeric::Index2D>
    KeypointSelectorLepetit::
    getKeypoints()
    {
      return m_locationVector;
    }


    float
    KeypointSelectorLepetit::
    getLaplacianMagnitudeThreshold()
    {
      return m_thresholdLaplacianMagnitude;
    }


    float
    KeypointSelectorLepetit::
    getPixelSimilarityThreshold()
    {
      return m_thresholdPixelSimilarity;
    }

  } // namespace brick

} // namespace computerVision
