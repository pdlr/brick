/**
***************************************************************************
* @file brick/computerVision/naiveSnake_impl.hh
*
* Header file defining inline and template functions from
* naiveSnake.hh.
*
* Copyright (C) 2006,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
*/

#ifndef BRICK_COMPUTERVISION_NAIVESNAKE_IMPL_HH
#define BRICK_COMPUTERVISION_NAIVESNAKE_IMPL_HH

// This file is included by naiveSnake.hh, and should not be directly included
// by user code, so no need to include naiveSnake.hh here.
// 
// #include <brick/computerVision/naiveSnake.hh>

#include <brick/computerVision/getEuclideanDistance.hh>
#include <brick/computerVision/naiveSnake.hh>
#include <brick/computerVision/sobel.hh>
#include <brick/linearAlgebra/linearAlgebra.hh>
#include <brick/numeric/bilinearInterpolator.hh>
#include <brick/numeric/utilities.hh>

namespace brick {

  namespace computerVision {
    
    template <class FloatType>
    Snake<FloatType>::
    Snake()
      : m_alpha(0.1),
        m_beta(0.0),   // 10.0 is reasonable, but slows things down a lot!
        m_betaVector(),
        m_contourIterations(10),
        m_cornerAdditionThreshold(-2.0),
        m_cornerDeletionThreshold(2.0),
        m_externalForceGradientX(),
        m_externalForceGradientY(),
        m_forceBalanceMatrix(),
        m_gamma(1.0),
        m_isClosed(true),
        m_isConverged(false),
        m_isFixed(true),
        m_kappa(1.0),
        m_maxIterations(100),
        m_maxSnakeSize(500),
        m_maxSpanLength(3.0),
        m_minSpanLength(0.5),
        m_snake()
    {
      // Empty.
    }


    template <class FloatType>
    void
    Snake<FloatType>::
    enableCornerAdditionAndDeletion(bool enableFlag)
    {
      if(!enableFlag) {
        // Disable corner addition and deletion by setting thresholds
        // to out-of-bounds values.
        m_cornerAdditionThreshold = -2.0;
        m_cornerDeletionThreshold = 2.0;
      }
    }

    
    template <class FloatType>
    std::vector< brick::numeric::Vector2D<FloatType> >
    Snake<FloatType>::
    run()
    {
      // Check state.
      if(m_externalForceGradientX.size() == 0) {
        BRICK_THROW(brick::common::StateException, "Snake::run()",
                  "Interest image has not been set.  Please call "
                  "setInterestImage() before calling run().");
      }
      
      // If snake has not been seeded, pick a reasonable default.
      if(m_snake.size() == 0) {
        size_t rows = m_externalForceGradientX.rows();
        size_t columns = m_externalForceGradientX.columns();
        m_snake.push_back(brick::numeric::Vector2D<FloatType>(columns / 4.0, rows / 4.0));
        m_betaVector.push_back(m_beta);
        m_snake.push_back(brick::numeric::Vector2D<FloatType>(columns / 4.0, 3 * rows / 4.0));
        m_betaVector.push_back(m_beta);
        m_snake.push_back(brick::numeric::Vector2D<FloatType>(3 * columns / 4.0, 3 * rows / 4.0));
        m_betaVector.push_back(m_beta);
        m_snake.push_back(brick::numeric::Vector2D<FloatType>(3 * columns / 4.0, rows / 4.0));
        m_betaVector.push_back(m_beta);
      }

      size_t iterationNumber = 0;
      while(!this->isConverged()) {
        if(iterationNumber >= m_maxIterations) {
          break;
        }
        this->runOneIteration();
        ++iterationNumber;
      }
      return m_snake;
    }


    template <class FloatType>
    std::vector< brick::numeric::Vector2D<FloatType> >
    Snake<FloatType>::
    runOneIteration()
    {
      // No need to do anything if we've already converged.
      if(m_isConverged) {
        return m_snake;
      }
      
      // Check state.
      if(m_externalForceGradientX.size() == 0) {
        BRICK_THROW(brick::common::StateException, "Snake::runOneIteration()",
                  "Interest image has not been set.  Please call "
                  "setInterestImage() before calling runOneIteration().");
      }

      // If snake has not been seeded, pick a reasonable default.
      if(m_snake.size() == 0) {
        BRICK_THROW(brick::common::StateException, "Snake::runOneIteration()",
                  "Snake has not been initialized.  Please call "
                  "setSeedPoints() before calling runOneIteration().");
      }
    
      if(m_snake.size() < 4) {
        BRICK_THROW(brick::common::StateException, "Snake::runOneIteration()",
                  "Snake is too short.");
      }

      if(m_snake.size() > m_maxSnakeSize) {
        BRICK_THROW(brick::common::StateException, "Snake::runOneIteration()",
                  "Snake is too long.");
      }

      this->adjustBetas(m_snake, m_betaVector);

      std::pair< std::vector< brick::numeric::Vector2D<FloatType> >, std::vector<FloatType> >
        snake_betaVector = this->resampleSnake();
      m_snake = snake_betaVector.first;
      m_betaVector = snake_betaVector.second;

      for(size_t updateNumber = 0; updateNumber < m_contourIterations;
          ++updateNumber) {
        std::vector< brick::numeric::Vector2D<FloatType> > oldSnake = m_snake;
        this->updateSnakePosition(m_snake);
        if(this->isConverged(m_snake, oldSnake)) {
          m_isConverged = true;
          break;
        }
      }
      return m_snake;
    }


    template <class FloatType>
    void
    Snake<FloatType>::
    setBendingConstant(FloatType beta)
    {
      m_beta = beta;
      for(size_t index0 = 0; index0 < m_betaVector.size(); ++index0) {
        if(m_betaVector[index0] != 0.0) {
          m_betaVector[index0] = m_beta;
        }
      }
      m_isConverged = false;
    }


    template <class FloatType>
    void
    Snake<FloatType>::
    setCornerAdditionAngle(FloatType theta)
    {
      if(theta < 0.0 || theta > 3.15) {
        BRICK_THROW(brick::common::ValueException, "Snake::setCornerAdditionAngle",
                  "Argument theta should be in the range [0.0, Pi].");
      }
      m_cornerAdditionThreshold = std::cos(theta);
    }


    template <class FloatType>
    void
    Snake<FloatType>::
    setCornerDeletionAngle(FloatType theta)
    {
      if(theta < 0.0 || theta > 3.15) {
        BRICK_THROW(brick::common::ValueException, "Snake::setCornerDeletionAngle",
                  "Argument theta should be in the range [0.0, Pi].");
      }
      m_cornerDeletionThreshold = std::cos(theta);
    }
    
    
    template <class FloatType>
    void
    Snake<FloatType>::
    setInterestImage(const Image<GRAY1>& interestImage)
    {
      size_t numberOfPassesUsed;
      brick::numeric::Array2D<FloatType> distanceMatrix =
        getEuclideanDistance<FloatType>(interestImage, 10, numberOfPassesUsed);
      std::cout << "Euclidean distance calculation took " << numberOfPassesUsed
                << " passes." << std::endl;
      m_externalForceGradientX =
        applySobelX(Image<GRAY_FLOAT64>(distanceMatrix)) / -8.0;
      m_externalForceGradientY =
        applySobelY(Image<GRAY_FLOAT64>(distanceMatrix)) / -8.0;
      m_isConverged = false;
    }


    template <class FloatType>
    void
    Snake<FloatType>::
    setSeedPoints(const std::vector< brick::numeric::Vector2D<FloatType> >& seedPoints)
    {
      std::vector<bool> cornerFlags(seedPoints.size(), false);
      this->setSeedPoints(seedPoints, cornerFlags);
    }


    template <class FloatType>
    void
    Snake<FloatType>::
    setSeedPoints(const std::vector< brick::numeric::Vector2D<FloatType> >& seedPoints,
                  const std::vector<bool> cornerFlags)
    {
      m_snake = seedPoints;
      m_betaVector = std::vector<FloatType>(cornerFlags.size());
      for(size_t index0 = 0; index0 < m_betaVector.size(); ++index0) {
        if(cornerFlags[index0]) {
          m_betaVector[index0] = 0.0;
        } else {
          m_betaVector[index0] = m_beta;
        }
      }
      m_isConverged = false;
    }

    
    template <class FloatType>
    void
    Snake<FloatType>::
    addResampledSpan(std::vector< brick::numeric::Vector2D<FloatType> >& snake,
                     std::vector<FloatType>& betaVector,
                     const brick::numeric::Vector2D<FloatType>& newPoint,
                     FloatType newBeta,
                     bool isLast)
    {
      // This loop pops elements off the back of the snake until the
      // gap between the new point and the previous point is bigger
      // than m_minSpanLength.  If, during this popping loop, the
      // length of the snake reaches zero, then the new point is added
      // and we return.  If we encounter a corner (small beta) point
      // in the course of our popping, then leave the corner and
      // ignore newPoint instead.
      size_t snakeSize = snake.size();
      FloatType spanLength = 0.0;
      while(1) {
        if(snakeSize == 0) {
          snake.push_back(newPoint);
          betaVector.push_back(newBeta);
          return;
        }
      
        // If the new point is sufficiently far from the previous point,
        // we can continue (leave the while loop).
        spanLength = brick::numeric::magnitude<FloatType>(newPoint - snake[snakeSize - 1]);
        if(spanLength >= m_minSpanLength) {
          break;
        }
        
        // If we get here, then the new point is very close to the
        // previous point.  Delete the previous pont to make the space
        // is bigger.  We assume that smaller beta values are more
        // unusual and should be preserved, so we delete the more
        // unusual of newPoint and previousPoint.  Being the final
        // point of an open snake is most unusual of all.
        if((newBeta <= betaVector[betaVector.size() - 1])
           || (isLast && !m_isClosed)) {
          snake.pop_back();
          betaVector.pop_back();
          --snakeSize;
        } else {
          // looks like the previous point was more important than
          // the one were about to add.  Since the previous point
          // and the one we're adding are very close together, we
          // simply won't add the new point.
          return;
        }
      }
      
      // If the new point is very far from the previous point, add
      // some intermediate points to fill up the gap.
      if(spanLength >= m_maxSpanLength) {
        int numberToAdd = int(spanLength / m_maxSpanLength);
        brick::numeric::Vector2D<FloatType> previousPoint = snake[snakeSize - 1];
        for(int extraPointIndex = 0; extraPointIndex < numberToAdd;
            ++extraPointIndex) {
          FloatType fraction =
            FloatType(extraPointIndex + 1) / FloatType(numberToAdd + 1);
          FloatType extraX = (fraction * newPoint.x()
                           + (1.0 - fraction) * previousPoint.x());
          FloatType extraY = (fraction * newPoint.y()
                           + (1.0 - fraction) * previousPoint.y());
          snake.push_back(brick::numeric::Vector2D<FloatType>(extraX, extraY));
          betaVector.push_back(m_beta);
        }
      }

      // Add the new point to the snake, unless this is the last
      // point of a closed snake (which would overlap the first
      // point).
      if(!(isLast && m_isClosed)) {
        snake.push_back(newPoint);
        betaVector.push_back(newBeta);
      }
    }


    template <class FloatType>
    void
    Snake<FloatType>::
    adjustBetas(const std::vector< brick::numeric::Vector2D<FloatType> >& snake,
                std::vector<FloatType>& betaVector)
    {
      if((m_cornerAdditionThreshold < -1.0)
         && (m_cornerDeletionThreshold >= 1.0)) {
        return;
      }

      // Note(xxx): Figure out what these were for.
      // size_t startIndex = 0;
      // size_t stopIndex = snake.size();
      if(!m_isClosed) {
        betaVector[0] = 0.0;
        betaVector[snake.size() - 1] = 0.0;
        // startIndex = 1;
        // stopIndex = snake.size() - 1;
      }
      
      for(size_t snakeIndex = 1; snakeIndex < snake.size() - 1; ++snakeIndex) {
        size_t previousIndex =
          (snakeIndex == 0) ? snake.size() - 1 : snakeIndex - 1;
        size_t nextIndex =
          (snakeIndex == snake.size() - 1) ? 0 : snakeIndex + 1;
          
        const brick::numeric::Vector2D<FloatType>& currentPoint = snake[snakeIndex];
        const brick::numeric::Vector2D<FloatType>& previousPoint = snake[previousIndex];
        const brick::numeric::Vector2D<FloatType>& nextPoint = snake[nextIndex];

        brick::numeric::Vector2D<FloatType> v0 = currentPoint - previousPoint;
        brick::numeric::Vector2D<FloatType> v1 = nextPoint - currentPoint;
        v0 /= brick::numeric::magnitude<FloatType>(v0);
        v1 /= brick::numeric::magnitude<FloatType>(v1);
        FloatType dotProduct = brick::numeric::dot<FloatType>(v0, v1);
        if(dotProduct <= m_cornerAdditionThreshold) {
          betaVector[snakeIndex] = 0.0;
        } else if(dotProduct > m_cornerDeletionThreshold) {
          betaVector[snakeIndex] = m_beta;
        }
      }
    }

      
    template <class FloatType>
    brick::numeric::Array2D<FloatType>
    Snake<FloatType>::
    buildForceBalanceMatrix(size_t numberOfSnakePoints)
    {
      // For a traditional snake, we want to minimize the total
      // energy:
      // 
      //   E_tot = \alpha * E_stretch + \beta * E_bend + \kappa * E_ext
      // 
      // Where
      // 
      //   E_stretch = \sum_j((x_j - x_j-1)^2 + (y_j - y_j-1)^2)
      //   E_bend = \sum_j(((x_j+1 - x_j) - (x_j - x_j-1))^2
      //                  + ((y_j+1 - y_j) - (y_j - y_j-1))^2)
      //   E_ext = \sum_j(-1 * InterestArray[x_j, y_j])
      // 
      // To find the minimum, we do the usual trick of setting the
      // gradient to zero.
      // 
      //   For all j,
      //   0 = \alpha * \partial(E_stretch) / \partial(x_j)
      //       + \beta * \partial(E_bend) / \partial(x_j)
      //       + \kappa * \partial(E_external) / \partial(x_j)
      //   And similarly for y_j.
      //
      // Once after taking the derivative, this gives us
      //
      //   For all j,
      //   0 = \alpha * (4 * x_j - 2 * x_j-1 - 2 * x_j+1)
      //       + \beta * (12 * x_j - 8 * x_j-1 - 8 * x_j+1 + 2 * x_j-2
      //                  + 2 * x_j+2)
      //       - \kappa * \partial(InterestArray) / \partial(x_j)
      //
      // This is a linear equation, which we can solve directly for the
      // x_j (and y_j).  Of course, we still have to iterate because the
      // external forces vary spacially.  To stabiize the system, we add
      // a "viscosity" term... a force proportional to the change in x_j
      // & y_j at each iteration.
      //
      //   For all j,
      //   0 = \alpha * (4 * x_j - 2 * x_j-1 - 2 * x_j+1)
      //       + \beta * (12 * x_j - 8 * x_j-1 - 8 * x_j+1 + 2 * x_j-2
      //                  + 2 * x_j+2)
      //       - \kappa * \partial(InterestArray) / \partial(x_j)
      //       + \gamma * (x_j - prevousX_j)
      //
      //   (\gamma * previousX_j
      //    + \kappa * \partial(InterestArray) / \partial(x_j))
      //   = (\alpha * (4 * x_j - 2 * x_j-1 - 2 * x_j+1)
      //      + \beta * (12 * x_j - 8 * x_j-1 - 8 * x_j+1 + 2 * x_j-2
      //                 + 2 * x_j+2)
      //      + \gamma * x_j)
      //
      // For our implementation, we'd like to be able to have "corner"
      // points which bend freely, so we maintain a different beta for
      // each value of j.  In our implementation, the betas all have
      // the same value, except at corner points, which have beta set
      // to 0.0.  Redoing the above derivation with unique betas for
      // each value of j, and rearranging, we get:
      // 
      //   For all j,
      //   (\kappa * \partial(InterestArray) / \partial(x_j)
      //    + \gamma * previousX_j)
      //   = (x_j * (4 * \alpha + 8 * \beta_j 
      //             + 2 * \beta_j-1 + 2 * \beta_j+1 + \gamma)
      //      - x_j-1 * (2 * \alpha + 4 * \beta_j + 4 * \beta_j-1)
      //      - x_j+1 * (2 * \alpha + 4 * \beta_j + 4 * \beta_j+1)
      //      + x_j-2 * (2 * \beta_j-1) + x_j+2 * (2 * \beta_j+1))
      //      
      // We still have a special case to consider: what if we have a
      // non-closed snake with fixed endpoints?  In this case, \beta_0
      // and \beta_N-1 are set to 0.0, and the equations for j == 0
      // and j == N-1 are trivially changed to keep x_j and x_N-1
      // fixed.
      //
      //   For j == 0, j == N - 1,
      //   x_j = previousX_j
      //
      // For a closed snake, we simply imagine that the tail (j = N-1)
      // is connected to the head (j = 0), so that indices j == 0, j
      // == -1, j == -2, etc. are equivalent to j == N, j == N - 1, j
      // = N - 2, etc.
      // 
      // These equations, which we'll call the force balance equations,
      // are reflected in AMatrix and the vectors returned by
      // buildForceBalanceRHS().  We invert AMatrix because we want to
      // solve for x_j, y_j.

      // We use the single letter variable N because we'll be indexing
      // with it many times below.
      size_t N = numberOfSnakePoints;
      
      // Argh!  It hurts to make (& invert) such a big matrix, knowing
      // that it'll be mostly zeros.
      brick::numeric::Array2D<FloatType> AMatrix(N, N);
      AMatrix = 0.0;
    
      // Handle the special case of a non-closed, fixed endpoints snake.
      size_t loopStartRow = 0;
      size_t loopStopRow = N;
      if((!m_isClosed) && m_isFixed) {
        AMatrix(0, 0) = 1.0;
        AMatrix(N - 1, N - 1) = 1.0;
        loopStartRow = 1;
        loopStopRow = N - 1;
      }

      // Fill in the matrix.
      for(size_t rowIndex = loopStartRow; rowIndex < loopStopRow;
          ++rowIndex) {
        // Sort out the indices.
        size_t currentJ = rowIndex;
        size_t previousJ = (rowIndex == 0) ? (N - 1) : (rowIndex - 1);
        size_t prepreviousJ = (previousJ == 0) ? (N - 1) : (previousJ - 1);
        size_t nextJ = (rowIndex == N - 1) ? 0 : (rowIndex + 1);
        size_t nextnextJ = (nextJ == N - 1) ? 0 : (nextJ + 1);

        // We allow different betas for each node so that some nodes
        // can be stiff and others can be corners.
        FloatType currentBeta = m_betaVector[currentJ];
        FloatType previousBeta = m_betaVector[previousJ];
        FloatType nextBeta = m_betaVector[nextJ];

        // Assign matrix elements as described in the comments above.
        AMatrix(currentJ, prepreviousJ) = 2.0 * previousBeta;
        AMatrix(currentJ, previousJ) = (
          -2.0 * m_alpha - 4.0 * currentBeta - 4.0 * previousBeta);
        AMatrix(currentJ, currentJ) = (
          4.0 * m_alpha + 8.0 * currentBeta + 2.0 * previousBeta
          + 2.0 * nextBeta + m_gamma);
        AMatrix(currentJ, nextJ) = (
          -2.0 * m_alpha - 4.0 * currentBeta - 4.0 * nextBeta);
        AMatrix(currentJ, nextnextJ) = 2.0 * nextBeta;
      }

      // Invert AMatrix and return.
      return brick::linearAlgebra::inverse(AMatrix);
    }


    template <class FloatType>
    void
    Snake<FloatType>::
    buildForceBalanceRHS(const std::vector< brick::numeric::Vector2D<FloatType> >& snake,
                         brick::numeric::Array1D<FloatType>& xRHS,
                         brick::numeric::Array1D<FloatType>& yRHS)
    {
      // Please read the comments in buildForceBalanceMatrix() in order
      // to understand this function.  For unimportant reasons, the
      // "Right Hand Side" quantities which we're building here show up
      // on the left hand side of the equations in those comments.
      if(xRHS.size() != snake.size()) {
        xRHS.reinit(snake.size());
      }
      if(yRHS.size() != snake.size()) {
        yRHS.reinit(snake.size());
      }

      // We use the single letter variable N because we'll be indexing
      // with it many times below.
      size_t N = snake.size();
      
      brick::numeric::BilinearInterpolator<FloatType> xInterp(m_externalForceGradientX);
      brick::numeric::BilinearInterpolator<FloatType> yInterp(m_externalForceGradientY);
      if((!m_isClosed) && m_isFixed) {
        xRHS[0] = snake[0].x();
        yRHS[0] = snake[0].y();
        xRHS[N - 1] = snake[N - 1].x();
        yRHS[N - 1] = snake[N - 1].y();
      } else {
        xRHS[0] = (m_gamma * snake[0].x()
                   + m_kappa * xInterp(snake[0].y(), snake[0].x()));
        yRHS[0] = (m_gamma * snake[0].y()
                   + m_kappa * yInterp(snake[0].y(), snake[0].x()));
        xRHS[N - 1] = (
          m_gamma * snake[N - 1].x()
          + m_kappa * xInterp(snake[N - 1].y(), snake[N - 1].x()));
        yRHS[N - 1] = (
          m_gamma * snake[N - 1].y()
          + m_kappa * yInterp(snake[N - 1].y(), snake[N - 1].x()));
      }
    
      for(size_t rowIndex = 1; rowIndex < N - 1; ++rowIndex) {
        const brick::numeric::Vector2D<FloatType>& snakePoint = snake[rowIndex];
        xRHS[rowIndex] = (m_gamma * snakePoint.x()
                          + m_kappa * xInterp(snakePoint.y(), snakePoint.x()));
        yRHS[rowIndex] = (m_gamma * snakePoint.y()
                          + m_kappa * yInterp(snakePoint.y(), snakePoint.x()));
      }
    }
  
  
    template <class FloatType>
    bool
    Snake<FloatType>::
    isConverged(const std::vector< brick::numeric::Vector2D<FloatType> >& snake,
                const std::vector< brick::numeric::Vector2D<FloatType> >& oldSnake)
    {
      bool convergedFlag = true;
      for(size_t pointNumber = 0; pointNumber < snake.size();
          ++pointNumber) {
        if(brick::numeric::magnitude<FloatType>(snake[pointNumber] - oldSnake[pointNumber]) >= 0.5) {
          convergedFlag = false;
          break;
        }
      }
      return convergedFlag;
    }
  
    
    template <class FloatType>
    std::pair< std::vector< brick::numeric::Vector2D<FloatType> >, std::vector<FloatType> >
    Snake<FloatType>::
    resampleSnake()
    {
      std::vector< brick::numeric::Vector2D<FloatType> > newSnake;
      std::vector<FloatType> newBetaVector;
      
      size_t snakeIndex = 0;
      for(; snakeIndex < m_snake.size() - 1; ++snakeIndex) {
        this->addResampledSpan(newSnake, newBetaVector,
                               m_snake[snakeIndex], m_betaVector[snakeIndex]);
      }

      if(m_isClosed) {
        this->addResampledSpan(newSnake, newBetaVector,
                               m_snake[snakeIndex], m_betaVector[snakeIndex]);
        this->addResampledSpan(newSnake, newBetaVector,
                               m_snake[0], m_betaVector[0], true);
      } else {
        this->addResampledSpan(newSnake, newBetaVector, m_snake[snakeIndex],
                               m_betaVector[snakeIndex], true);
      }

      return std::make_pair(newSnake, newBetaVector);
    }


    template <class FloatType>
    void
    Snake<FloatType>::
    updateSnakePosition(std::vector< brick::numeric::Vector2D<FloatType> >& snake)
    {
      if(snake.size() != m_forceBalanceMatrix.size()) {
        m_forceBalanceMatrix = this->buildForceBalanceMatrix(snake.size());
      }
      brick::numeric::Array1D<FloatType> xRHS;
      brick::numeric::Array1D<FloatType> yRHS;
      this->buildForceBalanceRHS(snake, xRHS, yRHS);
      brick::numeric::Array1D<FloatType> newX = brick::numeric::matrixMultiply<FloatType>(m_forceBalanceMatrix, xRHS);
      brick::numeric::Array1D<FloatType> newY = brick::numeric::matrixMultiply<FloatType>(m_forceBalanceMatrix, yRHS);
      for(size_t pointNum = 0; pointNum < snake.size(); ++pointNum) {
        snake[pointNum].setValue(newX[pointNum], newY[pointNum]);
      }
    }

  } // namespace computerVision
    
} // namespace brick

#endif /* #ifndef BRICK_COMPUTERVISION_NAIVESNAKE_IMPL_HH */
