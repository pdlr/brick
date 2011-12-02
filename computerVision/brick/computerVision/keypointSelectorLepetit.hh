/**
***************************************************************************
* @file brick/computerVision/keypointSelectorLepetit.hh
*
* Header file declaring a class template for selecting stable
* keypoints from an image.
*
* Copyright (C) 2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_KEYPOINTSELECTORLEPETIT_HH
#define BRICK_COMPUTERVISION_KEYPOINTSELECTORLEPETIT_HH

#include <vector>
#include <brick/computerVision/image.hh>
#include <brick/numeric/index2D.hh>

namespace brick {

  namespace computerVision {

    /**
     ** This class template selects stable keypoints from a grayscale
     ** input image.  It works by constructing a coarse scale-space
     ** pyramid, then selecting extrema from the scale space using a
     ** computationally inexpensive heuristic to discard extrema that
     ** have low contrast, or are poorly localized due to lying on
     ** intensity edges in the image.
     **
     ** Template argument Format specifies the format of the input
     ** image.
     **/
    class KeypointSelectorLepetit {
    public:

      // ========= Public member functions. =========

      KeypointSelectorLepetit();

      std::vector<numeric::Index2D> getKeypoints();


      Image<GRAY_FLOAT32> getImageLevel(unsigned int level);


      std::vector<unsigned int> getKeypointLevels();
      

      float getLaplacianMagnitudeThreshold();


      float getPixelSimilarityThreshold();


      template <ImageFormat Format>
      void
      setImage(Image<Format> const& image);

    private:

      void
      estimateThresholds(Image<GRAY_FLOAT32> const& image,
                         float& pixelSimilarityThreshold,
                         float& laplacianMagnitudeThreshold);

      void
      measurePixelThresholds(unsigned int row, unsigned int column,
                             Image<GRAY_FLOAT32> const& image,
                             float& pixelSimilarity,
                             float& laplacianMagnitude);

      bool
      testPixel(unsigned int row, unsigned int column,
                Image<GRAY_FLOAT32> const& image,
                float pixelSimilarityThreshold,
                float laplacianMagnitudeThreshold);



      std::vector<float> m_levelVector;
      std::vector<brick::numeric::Index2D> m_locationVector;
      float m_thresholdLaplacianMagnitude;
      float m_thresholdPixelSimilarity;
    };

  } // namespace computerVision
  
} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/keypointSelectorLepetit_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_KEYPOINTSELECTORLEPETIT_HH */
