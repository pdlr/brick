/**
***************************************************************************
* @file amanatidesWoo3DTest.cc
* 
* Source file defining AmanatidesWoo3DTest class.
*
* Copyright (C) 2004-2007,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <math.h>
#include <brick/numeric/array3D.hh>
#include <brick/numeric/vector3D.hh>
#include <brick/numeric/utilities.hh>
#include <brick/numeric/amanatidesWoo3D.hh>

#include <brick/test/testFixture.hh>

namespace brick {

  namespace numeric {

    class AmanatidesWoo3DTest
      : public brick::test::TestFixture<AmanatidesWoo3DTest> {

    public:

      AmanatidesWoo3DTest();
      ~AmanatidesWoo3DTest() {}

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
      compareImages(const Array3D<double>& testImage,
                    const Array3D<double>& groundTruthImage);

    
      void
      traverseImage(AmanatidesWoo3D< Array3D<double> >& aw3D);

    
      Vector3D<double> m_directionXYZ0;
      Vector3D<double> m_directionXYZ1;
      Vector3D<double> m_directionXYZ2;
      Array3D<double> m_downStreamImage;
      Array3D<double> m_fullImage;
      Transform3D<double> m_pixelTworld;
      Vector3D<double> m_startXYZ0;
      Vector3D<double> m_startXYZ1;
      Vector3D<double> m_startXYZ2;
      double m_testEpsilon;
      Array3D<double> m_zeroImage;
    
    }; // class AmanatidesWoo3DTest


    /* ============== Member Function Definititions ============== */

    AmanatidesWoo3DTest::
    AmanatidesWoo3DTest()
      : brick::test::TestFixture<AmanatidesWoo3DTest>("AmanatidesWoo3DTest"),
        m_directionXYZ0(),
        m_directionXYZ1(),
        m_directionXYZ2(),
        m_downStreamImage(),
        m_fullImage(),
        m_pixelTworld(),
        m_startXYZ0(),
        m_startXYZ1(),
        m_startXYZ2(),
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
      // space.  See the documentation file AmanatidesWoo3DTest.fig for
      // more details.
      m_downStreamImage = Array3D<double>(
        "[[[0.0,          0.0,          0.0,          0.0],"
        "  [0.0,          0.0,          0.0,          0.0],"
        "  [0.0,          0.0,          0.0,          0.0]],"
        " [[0.0,          0.0,          0.0,          0.0],"
        "  [0.0,          1.3017082793, 0.6508541397, 0.0],"
        "  [0.0,          0.0,          0.0,          0.0]],"
        " [[0.0,          0.0,          1.3017082793, 1.9525624189],"
        "  [0.0,          0.0,          0.6508541397, 0.0],"
        "  [0.0,          0.0,          0.0,          0.0]]]");

      m_fullImage = Array3D<double>(
        "[[[0.0,          0.0,          0.0,          0.0],"
        "  [0.6508541397, 0.0,          0.0,          0.0],"
        "  [1.3017082793, 0.0,          0.0,          0.0]],"
        " [[0.0,          0.0,          0.0,          0.0],"
        "  [0.6508541397, 2.6034165586, 0.6508541397, 0.0],"
        "  [0.0,          0.0,          0.0,          0.0]],"
        " [[0.0,          0.0,          1.3017082793, 1.9525624189],"
        "  [0.0,          0.0,          0.6508541397, 0.0],"
        "  [0.0,          0.0,          0.0,          0.0]]]");

      m_zeroImage = zeros<double>(3, 3, 4);

      // The "XYZ0" data members describe a ray which originates inside
      // the pixel array.
      m_startXYZ0 = Vector3D<double>(11.0, 8.0, 2.0);
      Vector3D<double> endXYZ0(17.0, 5.0, -2.0);
      m_directionXYZ0 = endXYZ0 - m_startXYZ0;
      m_directionXYZ0 /= magnitude<double>(m_directionXYZ0);

      // The "XYZ1" data members describe a ray which originates outside
      // the pixel array, but _does_ intersect the pixel array.
      // m_directionXYZ1 points toward the pixel array.
      m_startXYZ1 = Vector3D<double>(5.0, 11.0, 6.0);
      Vector3D<double> endXYZ1(17.0, 5.0, -2.0);
      m_directionXYZ1 = endXYZ1 - m_startXYZ1;
      m_directionXYZ1 /= magnitude<double>(m_directionXYZ1);

      // The "XYZ2" data members describe a ray which originates outside
      // the pixel array, but _does not_ intersect the pixel array.
      m_startXYZ2 = Vector3D<double>(5, 11.0, 6.0);
      Vector3D<double> endXYZ2(6.0, 6.0, 5.0);
      m_directionXYZ2 = endXYZ2 - m_startXYZ2;
      m_directionXYZ2 /= magnitude<double>(m_directionXYZ2);
    
      // Coordinate transformation which converts world coordinates to pixel
      // coordinates.
      m_pixelTworld = Transform3D<double>(0.5, 0.0, 0.0, -4.0,
                                  0.0, 0.5, 0.0, -2.5,
                                  0.0, 0.0, -0.5, 2.5,
                                  0.0, 0.0, 0.0, 1.0);
    }


    void
    AmanatidesWoo3DTest::
    testConstructor()
    {
      // Passing testBegin*() methods is sufficient.
    }


    void
    AmanatidesWoo3DTest::
    testCopyConstructor()
    {
      // Set up the image to be traversed.
      Array3D<double> traceImage = m_zeroImage.copy();

      // Construct the AmanatidesWoo3D instance.  Start point is inside
      // of the pixel array.  DownStreamOnly flag is true.
      AmanatidesWoo3D< Array3D<double> >
        aw3D(traceImage, m_pixelTworld, m_startXYZ0, m_directionXYZ0, true);

      // Construct a second AmanatidesWoo3D instance using the copy
      // constructor.
      AmanatidesWoo3D< Array3D<double> > aw3DCopy(aw3D);

      // Actual voxel traversal and image comparison is done in this
      // since it's identical for all of the testBegin*() functions.
      this->traverseImage(aw3DCopy);
      this->compareImages(traceImage, m_downStreamImage);
    }


    void
    AmanatidesWoo3DTest::
    testDestructor()
    {
      // No independent test for destructor.
    }

  
    void
    AmanatidesWoo3DTest::
    testBegin0()
    {
      // Set up the image to be traversed.
      Array3D<double> traceImage = m_zeroImage.copy();

      // Construct the AmanatidesWoo3D instance.  Start point is inside
      // of the pixel array.  DownStreamOnly flag is true.
      AmanatidesWoo3D< Array3D<double> >
        aw3D(traceImage, m_pixelTworld, m_startXYZ0, m_directionXYZ0, true);

      // Actual voxel traversal and image comparison is done in this
      // since it's identical for all of the testBegin*() functions.
      this->traverseImage(aw3D);
      this->compareImages(traceImage, m_downStreamImage);
    }

    void
    AmanatidesWoo3DTest::
    testBegin1()
    {
      // Set up the image to be traversed.
      Array3D<double> traceImage = m_zeroImage.copy();

      // Construct the AmanatidesWoo3D instance.  Start point is inside
      // of the pixel array.  DownStreamOnly flag is false.
      AmanatidesWoo3D< Array3D<double> >
        aw3D(traceImage, m_pixelTworld, m_startXYZ0, m_directionXYZ0, false);

      // Actual voxel traversal and image comparison is done in this
      // since it's identical for all of the testBegin*() functions.
      this->traverseImage(aw3D);
      this->compareImages(traceImage, m_fullImage);
    }

    void
    AmanatidesWoo3DTest::
    testBegin2()
    {
      // Set up the image to be traversed.
      Array3D<double> traceImage = m_zeroImage.copy();

      // Construct the AmanatidesWoo3D instance.  Start point is outside
      // of the pixel array.  DownStreamOnly flag is true.
      AmanatidesWoo3D< Array3D<double> >
        aw3D(traceImage, m_pixelTworld, m_startXYZ1, m_directionXYZ1, true);

      // Actual voxel traversal and image comparison is done in this
      // since it's identical for all of the testBegin*() functions.
      this->traverseImage(aw3D);
      this->compareImages(traceImage, m_fullImage);
    }

    void
    AmanatidesWoo3DTest::
    testBegin3()
    {
      // Set up the image to be traversed.
      Array3D<double> traceImage = m_zeroImage.copy();

      // Construct the AmanatidesWoo3D instance.  Start point is outside
      // of the pixel array.  DownStreamOnly flag is true, and direction
      // is away from the pixel array.
      AmanatidesWoo3D< Array3D<double> >
        aw3D(traceImage, m_pixelTworld, m_startXYZ1,
             -1.0 * m_directionXYZ1, true);

      // Actual voxel traversal and image comparison is done in this
      // since it's identical for all of the testBegin*() functions.
      this->traverseImage(aw3D);
      this->compareImages(traceImage, m_zeroImage);
    }

    void
    AmanatidesWoo3DTest::
    testBegin4()
    {
      // Set up the image to be traversed.
      Array3D<double> traceImage = m_zeroImage.copy();

      // Construct the AmanatidesWoo3D instance.  Start point is outside
      // of the pixel array.  Ray does not intersect the pixel array.
      AmanatidesWoo3D< Array3D<double> >
        aw3D(traceImage, m_pixelTworld, m_startXYZ2, m_directionXYZ2, false);

      // Actual voxel traversal and image comparison is done in this
      // since it's identical for all of the testBegin*() functions.
      this->traverseImage(aw3D);
      this->compareImages(traceImage, m_zeroImage);
    }


    void
    AmanatidesWoo3DTest::
    testEnd()
    {
      // Passing testBegin*() methods is sufficient.
    }


    void
    AmanatidesWoo3DTest::
    testGetData()
    {
      // Set up the image to be traversed.
      Array3D<double> traceImage = m_zeroImage.copy();

      // Construct the AmanatidesWoo3D instance.  Start point is inside
      // of the pixel array.  DownStreamOnly flag is true.
      AmanatidesWoo3D< Array3D<double> >
        aw3D(traceImage, m_pixelTworld, m_startXYZ0, m_directionXYZ0, true);

      // Actual voxel traversal and image comparison is done here.
      this->traverseImage(aw3D);
      this->compareImages(aw3D.getData(), traceImage);
    }

  
    void
    AmanatidesWoo3DTest::
    testValidIntersection0()
    {
      // Set up the image to be traversed.
      Array3D<double> traceImage = m_zeroImage.copy();

      // Construct the AmanatidesWoo3D instance.  Start point is inside
      // of the pixel array.  Ray does intersect the pixel array.
      AmanatidesWoo3D< Array3D<double> >
        aw3D(traceImage, m_pixelTworld, m_startXYZ0, m_directionXYZ0, true);

      // The validIntersection() member function should return true.
      BRICK_TEST_ASSERT(aw3D.validIntersection());
    }


    void
    AmanatidesWoo3DTest::
    testValidIntersection1()
    {
      // Set up the image to be traversed.
      Array3D<double> traceImage = m_zeroImage.copy();

      // Construct the AmanatidesWoo3D instance.  Start point is outside
      // of the pixel array.  DownStreamOnly flag is true, and direction
      // is away from the pixel array.
      AmanatidesWoo3D< Array3D<double> >
        aw3D(traceImage, m_pixelTworld, m_startXYZ1,
             -1.0 * m_directionXYZ1, true);

      // The validIntersection() member function should return false.
      BRICK_TEST_ASSERT(!(aw3D.validIntersection()));
    }

  
    void
    AmanatidesWoo3DTest::
    testAssignmentOperator()
    {
      // Set up the images to be traversed.
      Array3D<double> traceImage = m_zeroImage.copy();
      Array3D<double> traceImageCopy = m_zeroImage.copy();

      // Construct the AmanatidesWoo3D instance.  Start point is inside
      // of the pixel array.  DownStreamOnly flag is true.
      AmanatidesWoo3D< Array3D<double> >
        aw3D(traceImage, m_pixelTworld, m_startXYZ0, m_directionXYZ0, true);

      // Construct a second AmanatidesWoo3D having different contents.
      AmanatidesWoo3D< Array3D<double> >
        aw3DCopy(traceImageCopy, m_pixelTworld, m_startXYZ2, m_directionXYZ2,
                 false);
    
      // Copy using the assignment operator.
      aw3DCopy = aw3D;

      // Actual voxel traversal and image comparison is done in this
      // since it's identical for all of the testBegin*() functions.
      this->traverseImage(aw3DCopy);
      this->compareImages(traceImage, m_downStreamImage);
    }

  
    void
    AmanatidesWoo3DTest::
    compareImages(const Array3D<double>& testImage,
                  const Array3D<double>& groundTruthImage)
    {
      // Verify result.
      BRICK_TEST_ASSERT(testImage.shape0() == groundTruthImage.shape0());
      BRICK_TEST_ASSERT(testImage.shape1() == groundTruthImage.shape1());
      BRICK_TEST_ASSERT(testImage.shape2() == groundTruthImage.shape2());
      for(size_t index = 0; index < groundTruthImage.size(); ++index) {
        double residual = std::fabs(testImage(index) - groundTruthImage(index));
        BRICK_TEST_ASSERT(residual < m_testEpsilon);
      }
    }


    void
    AmanatidesWoo3DTest::
    traverseImage(AmanatidesWoo3D< Array3D<double> >& aw3D)
    {
      // Do the voxel traversal.
      typedef AmanatidesWoo3D< Array3D<double> >::iterator aw3DIterator;
      aw3DIterator endIterator = aw3D.end();
      for(aw3DIterator iter = aw3D.begin(); iter != endIterator; ++iter) {
        double parametricDistance = iter.tExit() - iter.tEntry();
        *iter = parametricDistance;
      }
    }

  } // namespace numeric

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::numeric::AmanatidesWoo3DTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::numeric::AmanatidesWoo3DTest currentTest;

}

#endif
