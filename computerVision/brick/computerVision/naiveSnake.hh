/**
***************************************************************************
* @file brick/computerVision/naiveSnake.hh
*
* Header file declaring a naive snakes implementation.
*
* Copyright (C) 2006,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_NAIVESNAKE_HH
#define BRICK_COMPUTERVISION_NAIVESNAKE_HH

#include <vector>
#include <brick/computerVision/image.hh>
#include <brick/numeric/vector2D.hh>

namespace brick {

  namespace computerVision {

    enum SnakeStrategy {
      SNAKE_EUCLIDEAN_DISTANCE
    };


    template <class FloatType>
    class Snake {
    public:

      Snake();


      ~Snake() {}


      void
      enableCornerAdditionAndDeletion(bool enableFlag=true);


      std::vector<FloatType>
      getBetaVector() {return m_betaVector;}


      bool
      isConverged() {return m_isConverged;}


      std::vector< brick::numeric::Vector2D<FloatType> >
      run();


      std::vector< brick::numeric::Vector2D<FloatType> >
      runOneIteration();


      void
      setBendingConstant(FloatType beta);


      void
      setCornerAdditionAngle(FloatType theta);


      void
      setCornerDeletionAngle(FloatType theta);


      void
      setClosedCurve(bool isClosed=true)
        {m_isClosed = isClosed; m_isConverged = false;}


      void
      setExternalForceConstant(FloatType kappa)
        {m_kappa = kappa; m_isConverged = false;}


      // void
      // setFixedEndpoints(bool isFixed=true)
      //   {m_isFixed = isFixed; m_isConverged = false;}


      void
      setInterestImage(const Image<GRAY1>& interestImage);


      void
      setMaxIterations(size_t maxIterations) {m_maxIterations = maxIterations;}


      void
      setMaxSnakeSize(size_t maxSize) {m_maxSnakeSize = maxSize;}


      void
      setMaximumSpanLength(size_t spanLength)
        {m_maxSpanLength = spanLength; m_isConverged = false;}


      void
      setMinimumSpanLength(size_t spanLength)
        {m_minSpanLength = spanLength; m_isConverged = false;}


      void
      setSeedPoints(const std::vector< brick::numeric::Vector2D<FloatType> >& seedPoints);


      void
      setSeedPoints(const std::vector< brick::numeric::Vector2D<FloatType> >& seedPoints,
                    const std::vector<bool> cornerFlags);


      void
      setStepsPerIteration(size_t numSteps) {m_contourIterations = numSteps;}


      void
      setStretchingConstant(FloatType alpha)
        {m_alpha = alpha; m_isConverged = false;}


      void
      setViscosityConstant(FloatType gamma)
        {m_gamma = gamma; m_isConverged = false;}


    private:

      void
      addResampledSpan(std::vector< brick::numeric::Vector2D<FloatType> >& snake,
                       std::vector<FloatType>& betaVector,
                       const brick::numeric::Vector2D<FloatType>& newPoint,
                       FloatType newBeta,
                       bool isLast=false);


      void
      adjustBetas(const std::vector< brick::numeric::Vector2D<FloatType> >& snake,
                  std::vector<FloatType>& betaVector);


      brick::numeric::Array2D<FloatType>
      buildForceBalanceMatrix(size_t numberOfSnakePoints);


      void
      buildForceBalanceRHS(const std::vector< brick::numeric::Vector2D<FloatType> >& snake,
                           brick::numeric::Array1D<FloatType>& xRHS,
                           brick::numeric::Array1D<FloatType>& yRHS);


      bool
      isConverged(const std::vector< brick::numeric::Vector2D<FloatType> >& snake,
                  const std::vector< brick::numeric::Vector2D<FloatType> >& oldSnake);


      std::pair< std::vector< brick::numeric::Vector2D<FloatType> >, std::vector<FloatType> >
      resampleSnake();

      void
      updateSnakePosition(std::vector< brick::numeric::Vector2D<FloatType> >& snake);


      FloatType m_alpha;
      FloatType m_beta;
      std::vector<FloatType> m_betaVector;
      size_t m_contourIterations;
      FloatType m_cornerAdditionThreshold;
      FloatType m_cornerDeletionThreshold;
      brick::numeric::Array2D<FloatType> m_externalForceGradientX;
      brick::numeric::Array2D<FloatType> m_externalForceGradientY;
      brick::numeric::Array2D<FloatType> m_forceBalanceMatrix;
      FloatType m_gamma;
      bool m_isClosed;
      bool m_isConverged;
      bool m_isFixed;
      FloatType m_kappa;
      size_t m_maxIterations;
      size_t m_maxSnakeSize;
      FloatType m_maxSpanLength;
      FloatType m_minSpanLength;
      std::vector< brick::numeric::Vector2D<FloatType> > m_snake;
    };

  } // namespace computerVision

} // namespace brick


// Include file containing definitions of inline and template
// functions.
#include <brick/computerVision/naiveSnake_impl.hh>

#endif /* #ifndef BRICK_COMPUTERVISION_NAIVESNAKE_HH */
