/**
***************************************************************************
* @file brick/computerVision/naiveSnakeTest.cc
*
* Source file defining tests for naive Snake class.
*
* Copyright (C) 2006,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/computerVision/naiveSnake.hh>
#include <brick/numeric/subArray2D.hh>
#include <brick/test/testFixture.hh>

#include <sstream>
#include <iomanip>
#include <brick/computerVision/utilities.hh>
#include <brick/computerVision/imageIO.hh>

namespace brick {

  namespace computerVision {

    class NaiveSnakeTest
      : public brick::test::TestFixture<NaiveSnakeTest> {

    public:

      NaiveSnakeTest();
      ~NaiveSnakeTest() {}
      void setVerboseOutput() {m_verbose = true;}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testExternalForce();
      void testStretchingAndBendingForces();

    private:

      void writeOutputFile(const std::vector< numeric::Vector2D<double> >& snakePoints,
                           const std::string& fileNameBase,
                           size_t iterationCount);


      std::vector<bool> m_cornerFlags1;
      Image<GRAY1> m_interestImage0;
      std::vector< numeric::Vector2D<double> > m_seedPoints0;
      std::vector< numeric::Vector2D<double> > m_seedPoints1;
      std::vector< numeric::Vector2D<double> > m_seedPoints2;
      size_t m_testImageColumns;
      size_t m_testImageRows;
      bool m_verbose;

    }; // class NaiveSnakeTest


    /* ============== Member Function Definititions ============== */

    NaiveSnakeTest::
    NaiveSnakeTest()
      : brick::test::TestFixture<NaiveSnakeTest>("NaiveSnakeTest"),
        m_cornerFlags1(),
        m_interestImage0(),
        m_seedPoints0(),
        m_seedPoints1(),
        m_testImageColumns(320),
        m_testImageRows(240),
        m_verbose(false)
    {
      BRICK_TEST_REGISTER_MEMBER(testExternalForce);
      BRICK_TEST_REGISTER_MEMBER(testStretchingAndBendingForces);

      // Build an interest image.
      m_interestImage0.reinit(m_testImageRows, m_testImageColumns);
      m_interestImage0 = false;
      numeric::subArray(m_interestImage0, numeric::Slice(70, 150), 100) = true;
      numeric::subArray(m_interestImage0, numeric::Slice(70, 150), 220) = true;
      numeric::subArray(m_interestImage0, 70, numeric::Slice(100, 221)) = true;
      numeric::subArray(m_interestImage0, 149, numeric::Slice(100, 221)) = true;

      // Build some initial snake points.
      m_seedPoints0.push_back(numeric::Vector2D<double>(80, 60));
      m_seedPoints0.push_back(numeric::Vector2D<double>(240, 60));
      m_seedPoints0.push_back(numeric::Vector2D<double>(240, 180));
      m_seedPoints0.push_back(numeric::Vector2D<double>(80, 180));

      m_seedPoints1.push_back(numeric::Vector2D<double>(36, 69));
      m_seedPoints1.push_back(numeric::Vector2D<double>(145, 6));
      m_seedPoints1.push_back(numeric::Vector2D<double>(311, 82));
      m_seedPoints1.push_back(numeric::Vector2D<double>(176, 237));
      m_cornerFlags1.push_back(true);
      m_cornerFlags1.push_back(true);
      m_cornerFlags1.push_back(true);
      m_cornerFlags1.push_back(true);

      m_seedPoints2.push_back(numeric::Vector2D<double>(145, 22));
      m_seedPoints2.push_back(numeric::Vector2D<double>(137, 217));
      m_seedPoints2.push_back(numeric::Vector2D<double>(166, 222));
      m_seedPoints2.push_back(numeric::Vector2D<double>(182, 20));
    }


    void
    NaiveSnakeTest::
    testExternalForce()
    {
      Snake<double> snake;
      snake.setInterestImage(m_interestImage0.copy());
      snake.setSeedPoints(m_seedPoints1);

      // Convigure snake to use only external force.
      snake.setStretchingConstant(0.0);
      snake.setBendingConstant(0.0);

      std::vector< numeric::Vector2D<double> > snakePoints;
      if(m_verbose) {
        // snake.setStepsPerIteration(1);
        // snake.setMinimumSpanLength(10);
        // snake.setMaximumSpanLength(30);
        size_t iterationCount = 0;
        while(!snake.isConverged()) {
          snakePoints = snake.runOneIteration();
          this->writeOutputFile(
            snakePoints, "externalForceTest", iterationCount);
          ++iterationCount;
        }
      } else {
        snakePoints = snake.run();
      }

      for(size_t index0 = 0; index0 < snakePoints.size(); ++index0) {
        size_t row = size_t(snakePoints[index0].y() + 0.5);
        size_t column = size_t(snakePoints[index0].x() + 0.5);
        BRICK_TEST_ASSERT(row < m_interestImage0.rows());
        BRICK_TEST_ASSERT(column < m_interestImage0.columns());
        BRICK_TEST_ASSERT(m_interestImage0(row, column) == true);
      }
    }



    void
    NaiveSnakeTest::
    testStretchingAndBendingForces()
    {
      Snake<double> snake;
      snake.setInterestImage(m_interestImage0.copy());
      snake.setSeedPoints(m_seedPoints1, m_cornerFlags1);

      // Convigure snake to use only external force.
      snake.setStretchingConstant(0.0);
      snake.setBendingConstant(10.0);
      snake.enableCornerAdditionAndDeletion();
      snake.setCornerAdditionAngle(0.25);
      snake.setCornerDeletionAngle(0.20);

      std::vector< numeric::Vector2D<double> > snakePoints;
      if(m_verbose) {
        // snake.setStepsPerIteration(1);
        // snake.setMinimumSpanLength(10);
        // snake.setMaximumSpanLength(30);
        size_t iterationCount = 0;
        while(!snake.isConverged()) {
          snakePoints = snake.runOneIteration();
          this->writeOutputFile(
            snakePoints, "stretchingAndBendingForceTest", iterationCount);
          ++iterationCount;
        }
      } else {
        snakePoints = snake.run();
      }

      for(size_t index0 = 0; index0 < snakePoints.size(); ++index0) {
        size_t row = size_t(snakePoints[index0].y() + 0.5);
        size_t column = size_t(snakePoints[index0].x() + 0.5);
        BRICK_TEST_ASSERT(row < m_interestImage0.rows());
        BRICK_TEST_ASSERT(column < m_interestImage0.columns());
        BRICK_TEST_ASSERT(m_interestImage0(row, column) == true);
      }
    }


    void
    NaiveSnakeTest::
    writeOutputFile(const std::vector< numeric::Vector2D<double> >& snakePoints,
                    const std::string& fileNameBase,
                    size_t iterationCount)
    {
      Image<GRAY1> resultImage(m_testImageRows, m_testImageColumns);
      resultImage = false;
      for(size_t index0 = 0; index0 < snakePoints.size(); ++index0) {
        size_t row = size_t(snakePoints[index0].y() + 0.5);
        size_t column = size_t(snakePoints[index0].x() + 0.5);
        BRICK_TEST_ASSERT(row < resultImage.rows());
        BRICK_TEST_ASSERT(column < resultImage.columns());
        resultImage(row, column) = true;
      }
      BRICK_TEST_ASSERT(resultImage.size() != 0);
      std::ostringstream fileNameStream;
      fileNameStream << fileNameBase << std::setw(3) << std::setfill('0')
                     << iterationCount << ".pgm";
      writePGM8(fileNameStream.str(), convertColorspace<GRAY8>(resultImage));
    }

  } // namespace computervision

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::computerVision::NaiveSnakeTest currentTest;
  if(argc > 1) {
    currentTest.setVerboseOutput();
  }
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::NaiveSnakeTest currentTest;

}

#endif
