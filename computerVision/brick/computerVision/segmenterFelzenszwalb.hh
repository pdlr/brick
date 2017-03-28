/**
***************************************************************************
* @file brick/computerVision/segmenterFelzenszwalb.hh
*
* Header file declaring an image segmentation class.
*
* Copyright (C) 2008,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_SEGMENTERFELZENSZWALB_HH
#define BRICK_COMPUTERVISION_SEGMENTERFELZENSZWALB_HH

#include <vector>
#include <brick/computerVision/disjointSet.hh>
#include <brick/computerVision/imageFilter.hh>
#include <brick/computerVision/image.hh>
#include <brick/computerVision/kernels.hh>
#include <brick/computerVision/utilities.hh>

namespace brick {

  namespace computerVision {

    template <class FloatType>
    class EdgeDefaultFunctor {
    public:
      template <ImageFormat FORMAT>
      FloatType operator()(Image<FORMAT> const& inImage,
                        size_t index0, size_t index1) {
        FloatType result = inImage[index1] - inImage[index0];
        return (result < 0.0) ? -result : result;
      }
    };
      

    template <class FloatType>
    struct Edge {
      size_t end0;
      size_t end1;
      FloatType weight;
    };


    /**
     ** This class implements the image segmentation algorithm
     ** described [1].  Essentially grouping pixels based on local
     ** differences so that segmented regions have similar local
     ** variances.
     **
     ** Here is an example of how to use this class:
     **
     ** @code
     **   Image<GRAY8> inputImage0 = readPGM8(getTestImageFileNamePGM0());
     **   SegmenterFelzenszwalb< EdgeDefaultFunctor<Float32> > segmenter(
     **     kappa, sigma, minSegmentSize);
     **   segmenter.segment(inputImage0);
     **   Array2D<UnsignedInt32> labelArray = segmenter.getLabelArray();
     ** @endCode
     **
     ** If you need more control over the edge weights that are used
     ** to do the (graph based) segmentation, you can do the
     ** following:
     **
     ** @code
     **   Image<GRAY8> inputImage0 = readPGM8(getTestImageFileNamePGM0());
     **   SegmenterFelzenszwalb<MyEdgeFunctor> segmenter(
     **     kappa, sigma, minSegmentSize);
     **   
     **   // Generate graph edges for segmentation.
     **   Image<HSV8> hsvImage = myPreprocessingRoutine(inputImage0);
     **   std::vector< Edge<Float32> > edges = segmenter.getEdges(hsvImage);
     **
     **   // Actually do the segmentation.
     **   segmenter.segmentFromEdges(hsvImage.rows(), hsvImage.columns(),
     **                              edges.begin(), edges.end());
     **   Array2D<UnsignedInt32> labelArray = segmenter.getLabelArray();
     ** @endcode
     **
     ** [1] Felzenszwalb, P., and Huttenlocher, D., Efficient
     ** Graph-Based Image Segmentation, International Journal of
     ** Computer Vision, Volume 59, Number 2, September 2004.
     **/
    template <class EdgeFunctor, class FloatType>
    class SegmenterFelzenszwalb {
    public:

      SegmenterFelzenszwalb(
        float k = 200.0f,
        float sigma = 0.8f,
        size_t minSegmentSize = 20,
        EdgeFunctor const& edgeFunctor = EdgeFunctor());

      
      virtual
      ~SegmenterFelzenszwalb() {}


      template <ImageFormat FORMAT>
      std::vector< Edge<FloatType> >
      getEdges(const Image<FORMAT>& inImage);

      
      template <ImageFormat FORMAT>
      std::vector< Edge<FloatType> >
      getEdges4Connected(const Image<FORMAT>& inImage);


      template <ImageFormat FORMAT>
      std::vector< Edge<FloatType> >
      getEdges8Connected(const Image<FORMAT>& inImage);


      virtual brick::numeric::Array2D<brick::common::UnsignedInt32>
      getLabelArray();

      
      virtual brick::numeric::Array2D<brick::common::UnsignedInt32>
      getLabelArray(brick::common::UnsignedInt32& numberOfSegments,
                    std::vector<size_t>& segmentSizes);


      template <ImageFormat FORMAT>
      void
      segment(const Image<FORMAT>& inputImage);


      template <class ITER>
      void
      segmentFromEdges(size_t imageRows, size_t imageColumns,
                       ITER edgeBegin, ITER edgeEnd);

      
    protected:

      typedef DisjointSet<float> Segment;


      inline float
      getCost(const Segment& C_i, const Segment& C_j);


      template <ImageFormat FORMAT>
      inline void
      setEdge(Edge<FloatType>& edge, size_t index0, size_t index1,
              Image<FORMAT> inImage);

      
      inline void
      updateCost(Segment& C_i, float weight);


      EdgeFunctor m_edgeFunctor;
      numeric::Index2D m_imageSize;
      float m_k;
      size_t m_minimumSegmentSize;
      brick::numeric::Array1D<Segment> m_segmentation;
      float m_sigma;
      size_t m_smoothSize;
      
    };


    template <class FloatType>
    inline bool
    operator<(Edge<FloatType> const& arg0, Edge<FloatType> const& arg1);
    
  } // namespace computerVision
  
} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/segmenterFelzenszwalb_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_SEGMENTERFELZENSZWALB_HH */
