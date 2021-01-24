/**
***************************************************************************
* @file brick/computerVision/featureAssociationTest.cc
*
* Source file defining tests for featureAssociation functions.
*
* Copyright (C) 2008,2012 David LaRose, dlr@davidlarose.com
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <algorithm>
#include <brick/computerVision/featureAssociation.hh>
#include <brick/test/testFixture.hh>

namespace {

  const double s_sigma = 4.0;

  class SimilarityFunctor {
  public:

    SimilarityFunctor(double sigma) : m_sigma(sigma) {}
    ~SimilarityFunctor() {}

    double operator()(brick::numeric::Vector2D<double> const& feature0,
                      brick::numeric::Vector2D<double> const& feature1) {
      double distance = brick::numeric::magnitude<double>(feature1 - feature0);
      return std::exp(-(distance * distance) / (2.0 * m_sigma * m_sigma));
    }

  private:
    double m_sigma;
  };

} // namespace


namespace brick {

  namespace computerVision {

    class FeatureAssociationTest
      : public brick::test::TestFixture<FeatureAssociationTest> {

    public:

      FeatureAssociationTest();
      ~FeatureAssociationTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testFeatureAssociation();

    private:

    }; // class FeatureAssociationTest


    /* ============== Member Function Definititions ============== */

    FeatureAssociationTest::
    FeatureAssociationTest()
      : brick::test::TestFixture<FeatureAssociationTest>("FeatureAssociationTest")
    {
      BRICK_TEST_REGISTER_MEMBER(testFeatureAssociation);
    }


    void
    FeatureAssociationTest::
    testFeatureAssociation()
    {
      // Make two sequences of features.
      std::vector< numeric::Vector2D<double> > features0;
      for(size_t ii = 0; ii < 4; ++ii) {
        features0.push_back(
          numeric::Vector2D<double>(static_cast<double>(ii),
                                    static_cast<double>(ii)));
      }
      std::vector< numeric::Vector2D<double> > features1;
      for(size_t ii = 0; ii < 4; ++ii) {
        features1.push_back(
          features0[ii] + numeric::Vector2D<double>(2.75, 2.0));
      }

      // Do the association, but put strong emphasis on proximity..
      std::vector< std::pair<size_t, size_t> > correspondences =
        associateFeaturesScott91<double>(features0.begin(), features0.end(),
                                         features1.begin(), features1.end(),
                                         SimilarityFunctor(1.0));

      // Correspondences should form based on shortest distance.
      bool twoGoesToZero = false;
      bool threeGoesToOne = false;
      for(size_t ii = 0; ii < correspondences.size(); ++ii) {
        if((correspondences[ii].first == 2)
           && (correspondences[ii].second == 0)) {
          twoGoesToZero = true;
        }
        if((correspondences[ii].first == 3)
           && (correspondences[ii].second == 1)) {
          threeGoesToOne = true;
        }
      }
      BRICK_TEST_ASSERT(twoGoesToZero);
      BRICK_TEST_ASSERT(threeGoesToOne);

      // Do the association, but put allow long-distance relationships.
      correspondences =
        associateFeaturesScott91<double>(features0.begin(), features0.end(),
                                         features1.begin(), features1.end(),
                                         SimilarityFunctor(4.0));

      // Now correspondence should form between the first, second,
      // third pairs in the two sequences.
      BRICK_TEST_ASSERT(correspondences.size() == features0.size());
      for(size_t ii = 0; ii < 4; ++ii) {
        BRICK_TEST_ASSERT(correspondences[ii].first == ii);
        BRICK_TEST_ASSERT(correspondences[ii].second == ii);
      }

      // Repeat the above, but reverse the order of features1.
      std::reverse(features1.begin(), features1.end());
      correspondences =
        associateFeaturesScott91<double>(features0.begin(), features0.end(),
                                         features1.begin(), features1.end(),
                                         SimilarityFunctor(4.0));

      // Now correspondence should form between the first, second,
      // third pairs in the two sequences.
      BRICK_TEST_ASSERT(correspondences.size() == features0.size());
      for(size_t ii = 0; ii < 4; ++ii) {
        BRICK_TEST_ASSERT(correspondences[ii].first == ii);
        BRICK_TEST_ASSERT(correspondences[ii].second
                        == features0.size() - ii - 1);
      }

    }

  } // namespace computerVision

} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::computerVision::FeatureAssociationTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::FeatureAssociationTest currentTest;

}

#endif
