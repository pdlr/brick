/**
***************************************************************************
* @file brick/numeric/test/amanatidesWoo2DTest.cc
* 
* Source file defining AmanatidesWoo2DTest class.
*
* Copyright (C) 2004-2007,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <math.h>
#include <brick/numeric/array2D.hh>
#include <brick/numeric/vector2D.hh>
#include <brick/numeric/utilities.hh>
#include <brick/numeric/amanatidesWoo2D.hh>

#include <brick/test/testFixture.hh>

namespace brick {

  namespace numeric {

    class AmanatidesWoo2DTest
      : public brick::test::TestFixture<AmanatidesWoo2DTest> {

    public:

      AmanatidesWoo2DTest();
      ~AmanatidesWoo2DTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      void testConstructor();
      void testCopyConstructor();
      void testDestructor();
      void testBegin0();
      void testBegin1();
      void testBegin2();
      void testBegin3();
      void testBegin4();
      void testEnd();
      void testGetData();
      void testValidIntersection0();
      void testValidIntersection1();
      void testAssignmentOperator();

    private:
      void
      compareImages(const Array2D<double>& testImage,
                    const Array2D<double>& groundTruthImage);
      void
      traverseImage(AmanatidesWoo2D< Array2D<double> >& aw2D);

    
      Vector2D<double> m_directionXY0;
      Vector2D<double> m_directionXY1;
      Vector2D<double> m_directionXY2;
      Array2D<double> m_downStreamImage;
      Array2D<double> m_fullImage;
      Transform2D<double> m_pixelTworld;
      Vector2D<double> m_startXY0;
      Vector2D<double> m_startXY1;
      Vector2D<double> m_startXY2;
      double m_testEpsilon;
      Array2D<double> m_zeroImage;
    
    }; // class AmanatidesWoo2DTest


    /* ============== Member Function Definititions ============== */

    AmanatidesWoo2DTest::
    AmanatidesWoo2DTest()
      : TestFixture<AmanatidesWoo2DTest>("AmanatidesWoo2DTest"),
        m_directionXY0(),
        m_directionXY1(),
        m_directionXY2(),
        m_downStreamImage(),
        m_fullImage(),
        m_pixelTworld(),
        m_startXY0(),
        m_startXY1(),
        m_startXY2(),
        m_testEpsilon(1.0e-8),
        m_zeroImage()
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testConstructor);
      BRICK_TEST_REGISTER_MEMBER(testCopyConstructor);
      BRICK_TEST_REGISTER_MEMBER(testDestructor);
      BRICK_TEST_REGISTER_MEMBER(testBegin0);
      BRICK_TEST_REGISTER_MEMBER(testBegin1);
      BRICK_TEST_REGISTER_MEMBER(testBegin2);
      BRICK_TEST_REGISTER_MEMBER(testBegin3);
      BRICK_TEST_REGISTER_MEMBER(testBegin4);
      BRICK_TEST_REGISTER_MEMBER(testEnd);
      BRICK_TEST_REGISTER_MEMBER(testGetData);
      BRICK_TEST_REGISTER_MEMBER(testValidIntersection0);
      BRICK_TEST_REGISTER_MEMBER(testValidIntersection1);
      BRICK_TEST_REGISTER_MEMBER(testAssignmentOperator);

      // Set up ground truth.  Each element of these arrays represents
      // ray traversal distance through the corresponding region of
      // space.  See the documentation file AmanatidesWoo2DTest.fig for
      // more details.
      m_downStreamImage = Array2D<double>(
        "[[0.0,          0.0,          1.1180339887, 2.2360679775],"
        " [0.0,          1.1180339887, 1.1180339887, 0.0],"
        " [0.0,          0.0,          0.0,          0.0]]");
      m_fullImage = Array2D<double>(
        "[[0.0,          0.0,          1.1180339887, 2.2360679775],"
        " [1.1180339887, 2.2360679775, 1.1180339887, 0.0],"
        " [1.1180339887, 0.0,          0.0,          0.0]]");
      m_zeroImage = zeros<double>(3, 4);

      // The "XY0" data members describe a ray which originates inside
      // the pixel array.
      m_startXY0 = Vector2D<double>(11.0, 8.0);
      Vector2D<double> endXY0(17.0, 5.0);
      m_directionXY0 = endXY0 - m_startXY0;
      m_directionXY0 /= magnitude<double>(m_directionXY0);

      // The "XY1" data members describe a ray which originates outside
      // the pixel array, but _does_ intersect the pixel array.
      // m_directionXY1 points toward the pixel array.
      m_startXY1 = Vector2D<double>(5.0, 11.0);
      Vector2D<double> endXY1(17.0, 5.0);
      m_directionXY1 = endXY1 - m_startXY1;
      m_directionXY1 /= magnitude<double>(m_directionXY1);

      // The "XY2" data members describe a ray which originates outside
      // the pixel array, but _does not_ intersect the pixel array.
      m_startXY2 = Vector2D<double>(5.0, 11.0);
      Vector2D<double> endXY2(6.0, 6.0);
      m_directionXY2 = endXY2 - m_startXY2;
      m_directionXY2 /= magnitude<double>(m_directionXY2);
    
      // Coordinate transformation which converts world coordinates to pixel
      // coordinates.
      m_pixelTworld = Transform2D<double>(0.5, 0.0, -4.0,
                                  0.0, 0.5, -2.5,
                                  0.0, 0.0, 1.0);
    }


    void
    AmanatidesWoo2DTest::
    testConstructor()
    {
      // Error case from the real world.
      Array2D<double> traceImage(480, 640);
      Vector2D<double> startPoint(88.799742766559362, 488.52994772105603);
      Vector2D<double> direction(10.0, 0.0);

      // This should not throw.
      AmanatidesWoo2D< Array2D<double> >
        aw2D(traceImage, numeric::Transform2D<double>(), startPoint, direction, true);
    }


    void
    AmanatidesWoo2DTest::
    testCopyConstructor()
    {
      // Set up the image to be traversed.
      Array2D<double> traceImage = m_zeroImage.copy();

      // Construct the AmanatidesWoo2D instance.  Start point is inside
      // of the pixel array.  DownStreamOnly flag is true.
      AmanatidesWoo2D< Array2D<double> >
        aw2D(traceImage, m_pixelTworld, m_startXY0, m_directionXY0, true);

      // Construct a second AmanatidesWoo2D instance using the copy
      // constructor.
      AmanatidesWoo2D< Array2D<double> > aw2DCopy(aw2D);

      // Actual voxel traversal and image comparison is done in this
      // since it's identical for all of the testBegin*() functions.
      this->traverseImage(aw2DCopy);
      this->compareImages(traceImage, m_downStreamImage);
    }


    void
    AmanatidesWoo2DTest::
    testDestructor()
    {
      // No independent test for destructor.
    }

  
    void
    AmanatidesWoo2DTest::
    testBegin0()
    {
      // Set up the image to be traversed.
      Array2D<double> traceImage = m_zeroImage.copy();

      // Construct the AmanatidesWoo2D instance.  Start point is inside
      // of the pixel array.  DownStreamOnly flag is true.
      AmanatidesWoo2D< Array2D<double> >
        aw2D(traceImage, m_pixelTworld, m_startXY0, m_directionXY0, true);

      // Actual voxel traversal and image comparison is done in this
      // since it's identical for all of the testBegin*() functions.
      this->traverseImage(aw2D);
      this->compareImages(traceImage, m_downStreamImage);
    }

    void
    AmanatidesWoo2DTest::
    testBegin1()
    {
      // Set up the image to be traversed.
      Array2D<double> traceImage = m_zeroImage.copy();

      // Construct the AmanatidesWoo2D instance.  Start point is inside
      // of the pixel array.  DownStreamOnly flag is false.
      AmanatidesWoo2D< Array2D<double> >
        aw2D(traceImage, m_pixelTworld, m_startXY0, m_directionXY0, false);

      // Actual voxel traversal and image comparison is done in this
      // since it's identical for all of the testBegin*() functions.
      this->traverseImage(aw2D);
      this->compareImages(traceImage, m_fullImage);
    }

    void
    AmanatidesWoo2DTest::
    testBegin2()
    {
      // Set up the image to be traversed.
      Array2D<double> traceImage = m_zeroImage.copy();

      // Construct the AmanatidesWoo2D instance.  Start point is outside
      // of the pixel array.  DownStreamOnly flag is true.
      AmanatidesWoo2D< Array2D<double> >
        aw2D(traceImage, m_pixelTworld, m_startXY1, m_directionXY1, true);

      // Actual voxel traversal and image comparison is done in this
      // since it's identical for all of the testBegin*() functions.
      this->traverseImage(aw2D);
      this->compareImages(traceImage, m_fullImage);
    }

    void
    AmanatidesWoo2DTest::
    testBegin3()
    {
      // Set up the image to be traversed.
      Array2D<double> traceImage = m_zeroImage.copy();

      // Construct the AmanatidesWoo2D instance.  Start point is outside
      // of the pixel array.  DownStreamOnly flag is true, and direction
      // is away from the pixel array.
      AmanatidesWoo2D< Array2D<double> >
        aw2D(traceImage, m_pixelTworld, m_startXY1, -1.0 * m_directionXY1,
             true);

      // Actual voxel traversal and image comparison is done in this
      // since it's identical for all of the testBegin*() functions.
      this->traverseImage(aw2D);
      this->compareImages(traceImage, m_zeroImage);
    }

    void
    AmanatidesWoo2DTest::
    testBegin4()
    {
      // Set up the image to be traversed.
      Array2D<double> traceImage = m_zeroImage.copy();

      // Construct the AmanatidesWoo2D instance.  Start point is outside
      // of the pixel array.  Ray does not intersect the pixel array.
      AmanatidesWoo2D< Array2D<double> >
        aw2D(traceImage, m_pixelTworld, m_startXY2, m_directionXY2, false);

      // Actual voxel traversal and image comparison is done in this
      // since it's identical for all of the testBegin*() functions.
      this->traverseImage(aw2D);
      this->compareImages(traceImage, m_zeroImage);
    }


    void
    AmanatidesWoo2DTest::
    testEnd()
    {
      // Passing testBegin*() methods is sufficient.
    }


    void
    AmanatidesWoo2DTest::
    testGetData()
    {
      // Set up the image to be traversed.
      Array2D<double> traceImage = m_zeroImage.copy();

      // Construct the AmanatidesWoo2D instance.  Start point is inside
      // of the pixel array.  DownStreamOnly flag is true.
      AmanatidesWoo2D< Array2D<double> >
        aw2D(traceImage, m_pixelTworld, m_startXY0, m_directionXY0, true);

      // Actual voxel traversal and image comparison is done here.
      this->traverseImage(aw2D);
      this->compareImages(aw2D.getData(), traceImage);
    }

  
    void
    AmanatidesWoo2DTest::
    testValidIntersection0()
    {
      // Set up the image to be traversed.
      Array2D<double> traceImage = m_zeroImage.copy();

      // Construct the AmanatidesWoo2D instance.  Start point is inside
      // of the pixel array.  Ray does intersect the pixel array.
      AmanatidesWoo2D< Array2D<double> >
        aw2D(traceImage, m_pixelTworld, m_startXY0, m_directionXY0, true);

      // The validIntersection() member function should return true.
      BRICK_TEST_ASSERT(aw2D.validIntersection());
    }


    void
    AmanatidesWoo2DTest::
    testValidIntersection1()
    {
      // Set up the image to be traversed.
      Array2D<double> traceImage = m_zeroImage.copy();

      // Construct the AmanatidesWoo2D instance.  Start point is outside
      // of the pixel array.  DownStreamOnly flag is true, and direction
      // is away from the pixel array.
      AmanatidesWoo2D< Array2D<double> >
        aw2D(traceImage, m_pixelTworld, m_startXY1, -1.0 * m_directionXY1,
             true);

      // The validIntersection() member function should return false.
      BRICK_TEST_ASSERT(!(aw2D.validIntersection()));
    }

  
    void
    AmanatidesWoo2DTest::
    testAssignmentOperator()
    {
      // Set up the images to be traversed.
      Array2D<double> traceImage = m_zeroImage.copy();
      Array2D<double> traceImageCopy = m_zeroImage.copy();

      // Construct the AmanatidesWoo2D instance.  Start point is inside
      // of the pixel array.  DownStreamOnly flag is true.
      AmanatidesWoo2D< Array2D<double> >
        aw2D(traceImage, m_pixelTworld, m_startXY0, m_directionXY0, true);

      // Construct a second AmanatidesWoo2D having different contents.
      AmanatidesWoo2D< Array2D<double> >
        aw2DCopy(traceImageCopy, m_pixelTworld, m_startXY2, m_directionXY2,
                 false);
    
      // Copy using the assignment operator.
      aw2DCopy = aw2D;

      // Actual voxel traversal and image comparison is done in this
      // since it's identical for all of the testBegin*() functions.
      this->traverseImage(aw2DCopy);
      this->compareImages(traceImage, m_downStreamImage);
    }

  
    void
    AmanatidesWoo2DTest::
    compareImages(const Array2D<double>& testImage,
                  const Array2D<double>& groundTruthImage)
    {
      // Verify result.
      BRICK_TEST_ASSERT(testImage.rows() == groundTruthImage.rows());
      BRICK_TEST_ASSERT(testImage.columns() == groundTruthImage.columns());
      for(size_t index = 0; index < groundTruthImage.size(); ++index) {
        double residual = std::fabs(testImage(index) - groundTruthImage(index));
        BRICK_TEST_ASSERT(residual < m_testEpsilon);
      }
    }


    void
    AmanatidesWoo2DTest::
    traverseImage(AmanatidesWoo2D< Array2D<double> >& aw2D)
    {
      // Do the voxel traversal.
      typedef AmanatidesWoo2D< Array2D<double> >::iterator aw2DIterator;
      aw2DIterator endIterator = aw2D.end();
      for(aw2DIterator iter = aw2D.begin(); iter != endIterator; ++iter) {
        double parametricDistance = iter.tExit() - iter.tEntry();
        *iter = parametricDistance;
      }
    }

  } // namespace numeric

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::numeric::AmanatidesWoo2DTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::AmanatidesWoo2DTest currentTest;

}

#endif
