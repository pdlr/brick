/**
***************************************************************************
* @file brick/computerVision/segmenterFelzenszwalb_impl.hh
*
* Header file defining inline and template functions declared in
* segmenterFelzenszwalb.hh.
*
* Copyright (C) 2008,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_SEGMENTERFELZENSZWALB_IMPL_HH
#define BRICK_COMPUTERVISION_SEGMENTERFELZENSZWALB_IMPL_HH

// This file is included by segmenterFelzenszwalb.hh, and should not be
// directly included by user code, so no need to include
// segmenterFelzenszwalb.hh here.
//
// #include <brick/computerVision/segmenterFelzenszwalb.hh>

namespace brick {

  namespace computerVision {

  } // namespace computerVision

} // namespace brick


/* ============ Definitions of inline & template functions ============ */


#include <cmath>
#include <limits>

namespace brick {

  namespace computerVision {

    template <class EdgeFunctor, class FloatType>
    SegmenterFelzenszwalb<EdgeFunctor, FloatType>::
    SegmenterFelzenszwalb(float k, float sigma, size_t minSegmentSize,
                          EdgeFunctor const& edgeFunctor)
      : m_edgeFunctor(edgeFunctor),
        m_imageSize(0, 0),
        m_k(k),
        m_minimumSegmentSize(minSegmentSize),
        m_segmentation(),
        m_sigma(sigma),
        m_smoothSize(static_cast<size_t>(std::fabs(6 * sigma + 1)))
    {
      // Smoothing kernel must not have even size or filter2D() will
      // complain later.
      if(m_smoothSize % 2 == 0) {
        m_smoothSize += 1;
      }
    }


    template <class EdgeFunctor, class FloatType>
    template <ImageFormat FORMAT>
    std::vector< Edge<FloatType> >
    SegmenterFelzenszwalb<EdgeFunctor, FloatType>::
    getEdges(const Image<FORMAT>& inImage)
    {
      return this->getEdges8Connected(inImage);
    }


    template <class EdgeFunctor, class FloatType>
    template <ImageFormat FORMAT>
    std::vector< Edge<FloatType> >
    SegmenterFelzenszwalb<EdgeFunctor, FloatType>::
    getEdges4Connected(const Image<FORMAT>& inImage)
    {
      // 4-connected means 2 undirected edges per pixel.
      size_t numEdges = 2 * inImage.size();

      // Except that some point off the side of the image.
      numEdges -= inImage.rows();

      // And some point off the bottom/top of the image.
      numEdges -= inImage.columns();

      std::vector< Edge<FloatType> > edges(numEdges);
      size_t interiorRows = inImage.rows() - 1;
      size_t interiorColumns = inImage.columns() - 1;
      size_t edgeNumber = 0;
      size_t pixelIndex0 = 0;
      size_t pixelIndex1;
      for(size_t row = 0; row < interiorRows; ++row) {
        // Get the first pixel of the row and the interior pixels of the row.
        for(size_t column = 0; column < interiorColumns; ++column) {
          pixelIndex1 = pixelIndex0 + 1;
          this->setEdge(edges[edgeNumber++], pixelIndex0, pixelIndex1, inImage);

          pixelIndex1 += (inImage.columns() - 1);
          this->setEdge(edges[edgeNumber++], pixelIndex0, pixelIndex1, inImage);
          ++pixelIndex0;
        }

        // Get the last pixel of the row.
        pixelIndex1 = pixelIndex0 + inImage.columns();
        this->setEdge(edges[edgeNumber++], pixelIndex0, pixelIndex1, inImage);
        ++pixelIndex0;
      }

      // Sanity check.
      if(pixelIndex0 != (inImage.rows() - 1) * inImage.columns()) {
        BRICK_THROW(brick::common::LogicException,
                    "SegmenterFelzenszwalb::getEdges4Connected()",
                    "Indexing error.");
      }

      // Get the last row of pixels;
      for(size_t column = 0; column < interiorColumns; ++column) {
        pixelIndex1 = pixelIndex0 + 1;
        this->setEdge(edges[edgeNumber++], pixelIndex0, pixelIndex1, inImage);
        ++pixelIndex0;
      }

      // Sanity check.
      if(edgeNumber != numEdges) {
        BRICK_THROW(brick::common::LogicException,
                    "SegmenterFelzenszwalb::getEdges4Connected()",
                    "Edge counting error.");
      }
      return edges;
    }


    template <class EdgeFunctor, class FloatType>
    template <ImageFormat FORMAT>
    std::vector< Edge<FloatType> >
    SegmenterFelzenszwalb<EdgeFunctor, FloatType>::
    getEdges8Connected(const Image<FORMAT>& inImage)
    {
      // 8-connected means 4 undirected edges per pixel.
      size_t numEdges = 4 * inImage.size();

      // Except that some point off the side of the image.
      numEdges -= 3 * inImage.rows();

      // And some point off the bottom/top of the image.
      numEdges -= 3 * inImage.columns();

      // Oops! We subtracted the bottom-right and bottom-left corners twice.
      numEdges += 2;

      std::vector< Edge<FloatType> > edges(numEdges);
      size_t interiorRows = inImage.rows() - 1;
      size_t interiorColumns = inImage.columns() - 1;
      size_t edgeNumber = 0;
      size_t pixelIndex0 = 0;
      size_t pixelIndex1;
      for(size_t row = 0; row < interiorRows; ++row) {
        // Get the first pixel of the row;
        pixelIndex1 = pixelIndex0 + 1;
        this->setEdge(edges[edgeNumber++], pixelIndex0, pixelIndex1, inImage);

        pixelIndex1 += inImage.columns();
        this->setEdge(edges[edgeNumber++], pixelIndex0, pixelIndex1, inImage);

        pixelIndex1 -= 1;
        this->setEdge(edges[edgeNumber++], pixelIndex0, pixelIndex1, inImage);
        ++pixelIndex0;

        // Get the interior pixels of the row.
        for(size_t column = 1; column < interiorColumns; ++column) {
          pixelIndex1 = pixelIndex0 + 1;
          this->setEdge(edges[edgeNumber++], pixelIndex0, pixelIndex1, inImage);

          pixelIndex1 += inImage.columns();
          this->setEdge(edges[edgeNumber++], pixelIndex0, pixelIndex1, inImage);

          pixelIndex1 -= 1;
          this->setEdge(edges[edgeNumber++], pixelIndex0, pixelIndex1, inImage);

          pixelIndex1 -= 1;
          this->setEdge(edges[edgeNumber++], pixelIndex0, pixelIndex1, inImage);
          ++pixelIndex0;
        }

        // Get the last pixel of the row.
        pixelIndex1 = pixelIndex0 + inImage.columns();
        this->setEdge(edges[edgeNumber++], pixelIndex0, pixelIndex1, inImage);

        pixelIndex1 -= 1;
        this->setEdge(edges[edgeNumber++], pixelIndex0, pixelIndex1, inImage);
        ++pixelIndex0;
      }

      // Sanity check.
      if(pixelIndex0 != (inImage.rows() - 1) * inImage.columns()) {
        BRICK_THROW(brick::common::LogicException,
                    "SegmenterFelzenszwalb::getEdges8Connected()",
                    "Indexing error.");
      }

      // Get the last row of pixels;
      for(size_t column = 0; column < interiorColumns; ++column) {
        pixelIndex1 = pixelIndex0 + 1;
        this->setEdge(edges[edgeNumber++], pixelIndex0, pixelIndex1, inImage);
        ++pixelIndex0;
      }

      // Sanity check.
      if(edgeNumber != numEdges) {
        BRICK_THROW(brick::common::LogicException,
                    "SegmenterFelzenszwalb::getEdges8Connected()",
                    "Edge counting error.");
      }
      return edges;
    }


    template <class EdgeFunctor, class FloatType>
    brick::numeric::Array2D<brick::common::UnsignedInt32>
    SegmenterFelzenszwalb<EdgeFunctor, FloatType>::
    getLabelArray()
    {
      brick::numeric::Array2D<brick::common::UnsignedInt32> labelArray(
        m_imageSize.getRow(), m_imageSize.getColumn());

      brick::numeric::Array1D<Segment>::iterator setIter =
        m_segmentation.begin();
      brick::numeric::Array2D<brick::common::UnsignedInt32>::iterator
        labelIter = labelArray.begin();
      while(setIter != m_segmentation.end()) {
        Segment& head = setIter->find();

        // Warning(xxx): Assuming we know something about how both
        // vectors and DisjointSets are implemented.
        brick::common::UnsignedInt32 labelIndex = static_cast<brick::common::UnsignedInt32>(
          &head - &(m_segmentation[0]));
        *labelIter = labelIndex;

        ++labelIter;
        ++setIter;
      }

      return labelArray;
    }


    template <class EdgeFunctor, class FloatType>
    brick::numeric::Array2D<brick::common::UnsignedInt32>
    SegmenterFelzenszwalb<EdgeFunctor, FloatType>::
    getLabelArray(brick::common::UnsignedInt32& numberOfSegments,
                  std::vector<size_t>& segmentSizes)
    {
      brick::numeric::Array2D<brick::common::UnsignedInt32> labelArray(
        m_imageSize.getRow(), m_imageSize.getColumn());
      std::vector<brick::common::UnsignedInt32> labelMap(
        m_segmentation.size(),
        std::numeric_limits<brick::common::UnsignedInt32>::max());
      brick::common::UnsignedInt32 currentLabel = 0;
      segmentSizes.clear();

      // Iterate over each pixel.
      brick::numeric::Array1D<Segment>::iterator setIter =
        m_segmentation.begin();
      brick::numeric::Array2D<brick::common::UnsignedInt32>::iterator
        labelIter = labelArray.begin();
      while(setIter != m_segmentation.end()) {

        // Figure out to which segment the current pixel belongs.
        //
        // Warning(xxx): Assuming we know something about how both
        // vectors and DisjointSets are implemented.
        Segment& head = setIter->find();
        brick::common::UnsignedInt32 labelIndex =
          static_cast<brick::common::UnsignedInt32>(
            &head - &(m_segmentation[0]));

        // Have we labeled this segment yet?
        if(labelMap[labelIndex] > currentLabel) {
          // No.  Label it now and remember how big the segment is.
          labelMap[labelIndex] = currentLabel;
          segmentSizes.push_back(head.getSize());
          ++currentLabel;
        }
        // Record the label in our output label image.
        *labelIter = labelMap[labelIndex];

        // Move on to next pixel.
        ++labelIter;
        ++setIter;
      }

      numberOfSegments = currentLabel;
      return labelArray;
    }


    template <class EdgeFunctor, class FloatType>
    template <ImageFormat FORMAT>
    void
    SegmenterFelzenszwalb<EdgeFunctor, FloatType>::
    segment(const Image<FORMAT>& inputImage)
    {
      m_imageSize.setValue(inputImage.rows(), inputImage.columns());

      // Smooth the image slightly to reduce artifacts.
      Image<GRAY_FLOAT32> smoothedImage;
      if(m_sigma == 0.0) {
        smoothedImage = convertColorspace<GRAY_FLOAT32>(inputImage);
      } else {
        Kernel<brick::common::Float32> gaussian =
          getGaussianKernelBySize<brick::common::Float32>(
            m_smoothSize, m_smoothSize,
            m_sigma, m_sigma);
        smoothedImage =
          filter2D<GRAY_FLOAT32, FORMAT>(
            gaussian, inputImage, brick::common::Float32(0));
      }

      // Get a vector of the edges in the image, sorted in ascending
      // order.
      std::vector< Edge<FloatType> > edges = this->getEdges(smoothedImage);

      this->segmentFromEdges(smoothedImage.rows(), smoothedImage.columns(),
                             edges.begin(), edges.end());
    }


    template <class EdgeFunctor, class FloatType>
    template <class ITER>
    void
    SegmenterFelzenszwalb<EdgeFunctor, FloatType>::
    segmentFromEdges(size_t imageRows, size_t imageColumns,
                     ITER edgeBegin, ITER edgeEnd)
    {
      size_t numPixels = imageRows * imageColumns;
      m_imageSize.setValue(imageRows, imageColumns);
      std::sort(edgeBegin, edgeEnd);

      // Start with segmentation S^0, where every vertex is its own component.
      typedef brick::numeric::Array1D<Segment>::iterator SegmentIter;
      m_segmentation.reinit(numPixels);
      for(SegmentIter segmentIter = m_segmentation.begin();
          segmentIter != m_segmentation.end(); ++segmentIter) {
        segmentIter->setPayload(m_k);
      }

      // Iteratively merge segments, as described in the paper.
      ITER edgeIter = edgeBegin;
      while(edgeIter != edgeEnd) {
        Segment& C_i = m_segmentation[edgeIter->end0].find();
        Segment& C_j = m_segmentation[edgeIter->end1].find();
        if(&C_i != &C_j) {
          float threshold = this->getCost(C_i, C_j);
          if(edgeIter->weight <= threshold) {
            C_i.merge(C_j);
            this->updateCost(C_i, edgeIter->weight);
          }
        }
        ++edgeIter;
      }

      // Merge any undersize segments, merging weak edges first.
      edgeIter = edgeBegin;
      while(edgeIter != edgeEnd) {
        Segment& C_i = m_segmentation[edgeIter->end0].find();
        Segment& C_j = m_segmentation[edgeIter->end1].find();
        if(C_i.getSize() < m_minimumSegmentSize
           || C_i.getSize() < m_minimumSegmentSize) {
          C_i.merge(C_j);
        }
        ++edgeIter;
      }
    }


    template <class EdgeFunctor, class FloatType>
    inline float
    SegmenterFelzenszwalb<EdgeFunctor, FloatType>::
    getCost(const Segment& C_i, const Segment& C_j)
    {
      return std::min(C_i.getPayload(), C_j.getPayload());
    }


    template <class EdgeFunctor, class FloatType>
    template <ImageFormat FORMAT>
    inline void
    SegmenterFelzenszwalb<EdgeFunctor, FloatType>::
    setEdge(Edge<FloatType>& edge, size_t index0, size_t index1,
            Image<FORMAT> inImage)
    {
      edge.end0 = index0;
      edge.end1 = index1;
      edge.weight = m_edgeFunctor(inImage, index0, index1);
    }


    template <class EdgeFunctor, class FloatType>
    inline void
    SegmenterFelzenszwalb<EdgeFunctor, FloatType>::
    updateCost(Segment& C_i, float weight)
    {
      Segment& head = C_i.find();
      head.setPayload(weight + m_k / head.getSize());
    }


    template <class FloatType>
    inline bool
    operator<(Edge<FloatType> const& arg0, Edge<FloatType> const& arg1) {
      return arg0.weight < arg1.weight;
    }


  } // namespace computerVision

} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_SEGMENTERFELZENSZWALB_IMPL_HH */
