/**
***************************************************************************
* @file brick/computerVision/test/extendedKalmanFilterTest.cc
*
* Source file defining tests for ExtendedKalmanFilter class.
*
* Copyright (C) 2009,2012 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <brick/computerVision/extendedKalmanFilter.hh>
#include <brick/linearAlgebra/linearAlgebra.hh>
#include <brick/numeric/utilities.hh>
#include <brick/test/testFixture.hh>


namespace brick {

  namespace computerVision {
    
    class ExtendedKalmanFilterTest
      : public brick::test::TestFixture<ExtendedKalmanFilterTest> {

    public:

      ExtendedKalmanFilterTest();
      ~ExtendedKalmanFilterTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      // Tests.
      void testExtendedKalmanFilter1D();
      void testExtendedKalmanFilter2D();
      void testCheckMeasurementJacobians();
      void testCheckProcessJacobians();

    private:

      // Internal class will be tested to make sure it tracks the process
      // correctly.
      class LinearProcessFilter
        : public ExtendedKalmanFilter<double> {
      public:
        
        LinearProcessFilter(
          numeric::Array2D<double> const& AMatrix,
          numeric::Array2D<double> const& BMatrix,
          numeric::Array2D<double> const& HMatrix,
          numeric::Array2D<double> const& processCovariance,
          numeric::Array2D<double> const& measurementCovariance)
          : ExtendedKalmanFilter<double>(),
            m_AMatrix(AMatrix),
            m_BMatrix(BMatrix),
            m_HHMatrix(HMatrix),
            m_processCovariance(processCovariance),
            m_measurementCovariance(measurementCovariance) {}

        virtual
        ~LinearProcessFilter() {}

      protected:

        virtual
        numeric::Array1D<double>
        applyMeasurementModel(unsigned int /* measurementID */,
                              double /* currentTime */,
                              double /* previousTime */,
                              numeric::Array1D<double> const& currentState) {
          return brick::numeric::matrixMultiply<double>(m_HHMatrix, currentState);
        }

        virtual
        numeric::Array1D<double>
        applyProcessModel(double /* currentTime */,
                          double /* previousTime */,
                          numeric::Array1D<double> const& previousState,
                          numeric::Array1D<double> const& controlInput) {
          return (brick::numeric::matrixMultiply<double>(m_AMatrix, previousState)
                  + brick::numeric::matrixMultiply<double>(m_BMatrix, controlInput));
        }

        virtual
        void
        getMeasurementJacobians(unsigned int /* measurementID */,
                                double /* currentTime */,
                                double /* previousTime */,
                                numeric::Array1D<double> const& /* state */,
                                numeric::Array2D<double>& stateJacobian,
                                numeric::Array2D<double>& noiseJacobian) {
          stateJacobian = m_HHMatrix.copy();
          noiseJacobian.reinit(m_HHMatrix.rows(), m_HHMatrix.rows());
          noiseJacobian = 0.0;
          for(unsigned int ii = 0; ii < m_HHMatrix.rows(); ++ii) {
            noiseJacobian(ii, ii) = 1.0;
          }
        }

        
        virtual
        void
        getProcessJacobians(double /* currentTime */,
                            double /* previousTime */,
                            numeric::Array1D<double> const& /* state */,
                            numeric::Array2D<double>& stateJacobian,
                            numeric::Array2D<double>& noiseJacobian) {
          stateJacobian = m_AMatrix.copy();
          noiseJacobian.reinit(m_AMatrix.rows(), m_AMatrix.rows());
          noiseJacobian = 0.0;
          for(unsigned int ii = 0; ii < m_AMatrix.rows(); ++ii) {
            noiseJacobian(ii, ii) = 1.0;
          }
        }


        virtual
        numeric::Array2D<double>
        getMeasurementNoiseCovariance(unsigned int /* measurementID */,
                                      double /* currentTime */,
                                      double /* previousTime */) {
          return m_measurementCovariance.copy();
        }

      
        virtual
        numeric::Array2D<double>
        getProcessNoiseCovariance(double /* currentTime */,
                                  double /* previousTime */) {

          return m_processCovariance.copy();
        }


      private:
        
        numeric::Array2D<double> m_AMatrix;
        numeric::Array2D<double> m_BMatrix;
        numeric::Array2D<double> m_HHMatrix;
        numeric::Array2D<double> m_processCovariance;
        numeric::Array2D<double> m_measurementCovariance;
      };


      // Create a 2D LinearProcessFilter instance for use in tests.
      LinearProcessFilter
      initialize2DFilter(numeric::Array1D<double>& initialState,
                         numeric::Array2D<double>& AMatrix,
                         numeric::Array2D<double>& HMatrix,
                         numeric::Array2D<double>& measurementNoiseArray);

      // Member variables.
      double m_estimateTolerance;
      numeric::Array2D<double> m_whiteNoiseArray;
      
    }; // class ExtendedKalmanFilterTest


    /* ============== Member Function Definititions ============== */

    ExtendedKalmanFilterTest::
    ExtendedKalmanFilterTest()
      : brick::test::TestFixture<ExtendedKalmanFilterTest>("ExtendedKalmanFilterTest"),
        m_estimateTolerance(0.1),
        m_whiteNoiseArray(
          "[[-0.34777470766943475, -2.5167002812462926, 0.29152891873455761, "
          "  0.24396752157073318, 0.32853635502875528, -0.62132836079277953, "
          "  0.52368697096200245, -0.51782573311816438, 0.4192428557019065, "
          "  1.1606517722474001, 0.80637113174256814, 0.28310761899516684, "
          "  -0.47826667807836459, -0.10507894125617172, -0.68110182690922216, "
          "  0.60364480345798088, 0.36271643140613952, -0.55467596241092398, "
          "  -1.3765244118896232, 0.54461300103674692, -0.94246670778733199, "
          "  0.65200364339242733, 0.66989746132746342, -0.62732446956750354, "
          "  1.2274013388721239, -0.024248580579372594, 0.37110010015522504, "
          "  0.020526535468818993, 0.83843571392376537, 0.46188430115330326, "
          "  -0.18384973535906426, -0.61866559129442311, 1.5244441262776489, "
          "  0.25787709909967271, 0.25886566027293562, 0.84953821841605681, "
          "  1.4881831176125766, -1.0060271426894807, 1.4244968344155386, "
          "  -0.38795069000407595, 0.73827056103971389, 0.21609474185760788, "
          "  0.97595299425837656, 0.3879175915569596, -0.50447601452334123, "
          "  -1.7641301338102178, 0.47949636657185152, 0.58184287554918135, "
          "  1.0276249945222393, 0.027941670250590848, 2.1275492895438988, "
          "  -1.4553905750013381, 0.25269044669218049, 1.2559982597002974, "
          "  -0.33920966003844238, -0.53482175824681788, -1.9048889705807066, "
          "  0.57630942376919458, -0.7481389490641317, 0.2915494788049442, "
          "  0.992382154727939, -0.021835641171445742, -1.1675319658918746, "
          "  -2.1457907925912316, 0.058897368370674268, 1.2145386026940457, "
          "  -0.10695745372508686, -1.0831169431521384, 0.24697964299620273, "
          "  0.067598288991645286, -0.29329253154169782, -0.05136016939704467, "
          "  -0.50367509460100357, -0.18147759707353134, 0.88762647432794772, "
          "  -2.1513482352671591, 0.64670859601119191, 0.64687476170010949, "
          "  -0.86735137804769602, -0.79352051833799608, 0.91992158216691833, "
          "  -0.4162532408332742, -0.63855038173405243, 0.47858372194082172, "
          "  -1.9603028846547463, -2.0826280114288029, -1.3843484308933014, "
          "  0.10897286155833757, -0.11813368610509058, -0.6764840166360635, "
          "  0.46740866370859813, 0.055324063525083594, 0.060453302741196716, "
          "  -1.4316303321020896, -0.51543971186200843, 0.62477520368391182, "
          "  -0.8907335721675721, 1.5763340404363035, 0.74813617964686385, "
          "  0.15863633423218287],\n"
          " [0.19152037021367921, 0.51723475112534856, 1.3697103309467167, "
          "  0.58693424957880325, -0.96055404770804431, 0.32733478821244799, "
          "  1.234038403951923, -0.17063766977183581, 0.21436767390214836, "
          "  0.60831019118203755, 1.1711968810158098, 0.0034351992924098406, "
          "  0.1432223716435139, -0.1160251240365664, -0.8231304346488012, "
          "  0.93471405595394752, 2.3022570456706832, -0.26341967019828461, "
          "  0.41633699088272208, 0.60429988023256775, 0.61150138708632051, "
          "  0.70694049388113533, 0.059994940135596392, -1.3536931912881016, "
          "  0.86551191108138037, -0.014754809660149869, 0.61431357642385576, "
          "  -2.473493834070517, -1.3710385099847242, 0.53832548412384973, "
          "  -1.290945309899201, -1.9678217691072832, 0.83337261541372754, "
          "  0.628005611923379, -0.21153575573739417, -1.0136828527776782, "
          "  0.26235090551533552, 0.34089461359767825, -0.7778706390505068, "
          "  -0.58738511003358296, -0.95337528774661306, -1.0438980744847162, "
          "  0.66986824933520661, -0.89167941414896146, -1.1696765054855494, "
          "  1.4523177530841171, 0.27322672643687435, 0.17871325067265892, "
          "  0.24247894198109171, 0.49240900239066104, -0.16537802926198711, "
          "  1.5192716454890383, -2.4207782371148467, 1.1607306804820809, "
          "  -0.40788804171170634, 0.86373855301478375, 1.2321737893762994, "
          "  -0.92333020953652578, -0.052086158017863045, 0.81359520662332208, "
          "  -0.19757308631674106, -0.12274574960765659, -0.43160885691996748, "
          "  -2.4663493829574903, 0.94939389527639062, -0.1646899611827298, "
          "  -0.19229972263328882, 0.049532378760731473, 1.7443379859282322, "
          "  1.0288827452322873, -1.4338342518761851, 0.1507231212681088, "
          "  0.57290700168905251, 0.75579521047691312, -0.29139156820661227, "
          "  0.61447406215995115, -0.26008747536453186, -0.63205663852031002, "
          "  1.6190934136633608, -1.3714413766234732, -1.8977630964850098, "
          "  0.28709008249167794, -0.59910897799463358, 0.55555620615171641, "
          "  0.3597678013167262, 0.61103209472202746, 1.1230001808674461, "
          "  0.57688844295474284, -0.11983423621864027, 0.33797755092727599, "
          "  1.1781804236683893, 0.87866645997628923, -1.8322910442993776, "
          "  1.3972482839627325, -0.086548914824670795, 1.0361363720602994, "
          "  -0.090028476206786348, 0.89160201952527796, -1.4651344893933715, "
          "  0.75194168019873209]]")
    {
      BRICK_TEST_REGISTER_MEMBER(testExtendedKalmanFilter1D);
      BRICK_TEST_REGISTER_MEMBER(testExtendedKalmanFilter2D);
      BRICK_TEST_REGISTER_MEMBER(testCheckMeasurementJacobians);
      BRICK_TEST_REGISTER_MEMBER(testCheckProcessJacobians);
    }


    void
    ExtendedKalmanFilterTest::
    testExtendedKalmanFilter1D()
    {
      const double xActual = -0.3;
      const double x0 = 0.0; // Initial x estimate.
      const double p0 = 1.0; // Estimated variance of our initial guess.

      const double processVariance = 1.0e-5;
      const double measurementSigma = 1.0e-1;

      // x_k = x_(k-1) + w_k
      const numeric::Array2D<double> AMatrix("[[1.0]]");
      const numeric::Array2D<double> BMatrix("[[0.0]]");

      // z_k = 2.5 * x_k + v_k
      const numeric::Array2D<double> HMatrix("[[2.5]]");

      numeric::Array2D<double> processCovariance(1, 1);
      numeric::Array2D<double> measurementCovariance(1, 1);
      processCovariance(0, 0) = processVariance;
      measurementCovariance(0, 0) = measurementSigma * measurementSigma;
      
      // Put initial state in appropriate form for the filter.
      numeric::Array1D<double> initialState(1);
      initialState[0] = x0;
      numeric::Array2D<double> initialVariance(1, 1);
      initialVariance(0, 0) = p0;

      // Pre-generate our noise.
      numeric::Array1D<double> measurementNoiseArray =
        m_whiteNoiseArray.getRow(0) * measurementSigma;
      numeric::Array1D<double> measurementArray =
        measurementNoiseArray + HMatrix[0] * xActual;
      
      LinearProcessFilter filter(AMatrix, BMatrix, HMatrix,
                                 processCovariance, measurementCovariance);
      filter.setStateEstimate(0, initialState, initialVariance);

      // Run the filter.
      numeric::Array1D<double> timestampArray(measurementArray.size());
      numeric::Array1D<double> estimateArray(measurementArray.size());
      numeric::Array1D<double> varianceArray(measurementArray.size());
      numeric::Array1D<double> measurement(1);
      numeric::Array1D<double> dummyControlInput(1);
      dummyControlInput = 0.0;
      for(unsigned int count = 0; count < measurementArray.size(); ++count) {
        // Update the filter
        measurement[0] = measurementArray[count];
        filter.addMeasurement(0, count + 1, measurement, dummyControlInput);

        // Query and record the resulting state estimate.
        numeric::Array1D<double> xHat(1);
        numeric::Array2D<double> xVarianceHat(1, 1);
        filter.getStateEstimate(timestampArray[count], xHat, xVarianceHat);
        estimateArray[count] = xHat[0];
        varianceArray[count] = xVarianceHat[0];
      }

      // Make sure we converged.
      for(unsigned int ii = estimateArray.size() / 2;
          ii < estimateArray.size(); ++ii) {
        BRICK_TEST_ASSERT(timestampArray[ii] == ii + 1);
        BRICK_TEST_ASSERT(std::fabs(estimateArray[ii] - xActual)
                        < m_estimateTolerance);
        BRICK_TEST_ASSERT(varianceArray[ii] > 0.0);
        BRICK_TEST_ASSERT(varianceArray[ii] < m_estimateTolerance);
      }

      // Make sure variance is monotonically decreasing.
      for(unsigned int ii = 1; ii < estimateArray.size(); ++ii) {
        BRICK_TEST_ASSERT(varianceArray[ii] > 0.0);
        BRICK_TEST_ASSERT(varianceArray[ii] < varianceArray[ii - 1]);
      }

      // Make sure variance estimate is appropriate.  We'll pick
      // thresholds here that are likely to pass with, high confidence
      // if the filter is working right.  Hopefully this particular
      // test will be one of the cases where it passes.  If so, we're
      // in good shape, because we're using the same pseudorandom data
      // each time.
      bool isReasonableFlag = false;
      for(unsigned int ii = 1; ii < estimateArray.size(); ++ii) {
        double residual = std::fabs(xActual - estimateArray[ii]);
        double sigmaHat = std::sqrt(varianceArray[ii]);
        double fourSigma = 4.0 * sigmaHat;
        BRICK_TEST_ASSERT(residual < fourSigma);
        isReasonableFlag |= (residual > sigmaHat);
      }
      BRICK_TEST_ASSERT(isReasonableFlag);
    }


    void
    ExtendedKalmanFilterTest::
    testExtendedKalmanFilter2D()
    {
      // Get a filter to test.
      numeric::Array1D<double> initialState;
      numeric::Array2D<double> AMatrix;
      numeric::Array2D<double> HMatrix;
      numeric::Array2D<double> measurementNoiseArray;
      LinearProcessFilter filter = this->initialize2DFilter(
        initialState, AMatrix, HMatrix, measurementNoiseArray);
      
      // Run the filter.
      numeric::Array1D<double> timestampArray(m_whiteNoiseArray.columns());
      numeric::Array2D<double> referenceArray(
        AMatrix.rows(), m_whiteNoiseArray.columns());
      numeric::Array2D<double> estimateArray(
        AMatrix.rows(), m_whiteNoiseArray.columns());
      numeric::Array2D<double> covarianceArray(
        AMatrix.size(), m_whiteNoiseArray.columns());
      numeric::Array1D<double> sqErrorArray(estimateArray.columns());
      numeric::Array1D<double> actualState = initialState.copy();
      numeric::Array1D<double> dummyControlInput(1);
      dummyControlInput = 0.0;
      for(unsigned int count = 0; count < m_whiteNoiseArray.columns();
          ++count) {
        // Update actual state, and record it for later reference.
        actualState = brick::numeric::matrixMultiply<double>(AMatrix, actualState);
        for(unsigned int jj = 0; jj < actualState.size(); ++jj) {
          referenceArray(jj, count) = actualState[jj];
        }

        // Synthesize a measurement.
        numeric::Array1D<double> measurement = brick::numeric::matrixMultiply<double>(HMatrix, actualState);
        for(unsigned int jj = 0; jj < measurement.size(); ++jj) {
          measurement[jj] += measurementNoiseArray(jj, count);
        }
        
        // Update the filter
        filter.addMeasurement(0, count + 1, measurement, dummyControlInput);

        // Query and record the resulting state estimate.
        numeric::Array1D<double> xHat;
        numeric::Array2D<double> xCovarianceHat;
        filter.getStateEstimate(timestampArray[count], xHat, xCovarianceHat);
        for(unsigned int row = 0; row < xHat.size(); ++row) {
          estimateArray(row, count) = xHat(row);
          for(unsigned int column = 0; column < xHat.size(); ++column) {
            covarianceArray(column + row * xHat.size(), count) =
              xCovarianceHat(row, column);
          }
        }

        // Record residual error.
        sqErrorArray[count] = 0.0;
        for(size_t jj = 0; jj < estimateArray.rows(); ++jj) {
          double residual =
            referenceArray(jj, count) - estimateArray(jj, count);
          sqErrorArray[count] += residual * residual;
        }
        
      }

      // Make sure we converged.
      for(unsigned int ii = estimateArray.columns() / 2;
          ii < estimateArray.columns(); ++ii) {
        BRICK_TEST_ASSERT(timestampArray[ii] == ii + 1);
        BRICK_TEST_ASSERT(std::sqrt(sqErrorArray[ii]) < m_estimateTolerance);

        for(size_t jj = 0; jj < estimateArray.rows(); ++jj) {
          size_t tmpIndex = jj + jj * estimateArray.rows();
          BRICK_TEST_ASSERT(covarianceArray(tmpIndex, ii) > 0.0);
          BRICK_TEST_ASSERT(covarianceArray(tmpIndex, ii) < m_estimateTolerance);
        }
      }

      // Make sure variance is generally decreasing.
      numeric::Array1D<double> maxEigenvalueArray(estimateArray.columns());
      for(unsigned int ii = 0; ii < estimateArray.columns(); ++ii) {
        numeric::Array2D<double> covariance(estimateArray.rows(), estimateArray.rows());
        for(unsigned int jj = 0; jj < covariance.size(); ++jj) {
          covariance[jj] = covarianceArray(jj, ii);
        }
        maxEigenvalueArray[ii] =
          linearAlgebra::eigenvaluesSymmetric(covariance)[0];
        BRICK_TEST_ASSERT(maxEigenvalueArray[ii] > 0.0);
        if(ii > 0) {
          // Should be non-increasing, not counting noise.
          double decrease = maxEigenvalueArray[ii - 1] - maxEigenvalueArray[ii];
          BRICK_TEST_ASSERT(decrease >= 0.0);
        }
      }

      // Make sure variance estimate is appropriate.  We'll pick
      // thresholds here that are likely to pass with, high confidence
      // if the filter is working right.  Hopefully this particular
      // test will be one of the cases where it passes.  If so, we're
      // in good shape, because we're using the same pseudorandom data
      // each time.
      bool isReasonableFlag = false;
      for(unsigned int ii = 10; ii < estimateArray.columns(); ++ii) {
        double sigmaHat = std::sqrt(maxEigenvalueArray[ii]);
        double fourSigma = 4.0 * sigmaHat;
        BRICK_TEST_ASSERT(std::sqrt(sqErrorArray[ii]) < fourSigma);
        isReasonableFlag |= (std::sqrt(sqErrorArray[ii]) > sigmaHat);
      }

      BRICK_TEST_ASSERT(isReasonableFlag);
    }


    void
    ExtendedKalmanFilterTest::
    testCheckMeasurementJacobians()
    {
      // Get a filter to test.
      numeric::Array1D<double> initialState;
      numeric::Array2D<double> AMatrix;
      numeric::Array2D<double> HMatrix;
      numeric::Array2D<double> measurementNoiseArray;
      LinearProcessFilter filter = this->initialize2DFilter(
        initialState, AMatrix, HMatrix, measurementNoiseArray);

      // Call the function under test.
      numeric::Array1D<double> epsilonArray(2);
      epsilonArray = 1.0E-6;
      numeric::Array2D<double> residualArray;
      filter.checkMeasurementJacobians(0, epsilonArray, residualArray);

      // Because the process is actually linear, we expect the
      // residual to be negligible.
      BRICK_TEST_ASSERT(residualArray.rows() == 2);
      BRICK_TEST_ASSERT(residualArray.columns() == 2);
      for(unsigned int ii = 0; ii < residualArray.size(); ++ii) {
        BRICK_TEST_ASSERT(std::fabs(residualArray[ii]) < 1.0E-8);
      }
    }


    void
    ExtendedKalmanFilterTest::
    testCheckProcessJacobians()
    {
      // Get a filter to test.
      numeric::Array1D<double> initialState;
      numeric::Array2D<double> AMatrix;
      numeric::Array2D<double> HMatrix;
      numeric::Array2D<double> measurementNoiseArray;
      LinearProcessFilter filter = this->initialize2DFilter(
        initialState, AMatrix, HMatrix, measurementNoiseArray);

      // Call the function under test.
      numeric::Array1D<double> epsilonArray(2);
      epsilonArray = 1.0E-6;
      numeric::Array1D<double> controlInput(1);
      controlInput = 0.0;
      numeric::Array2D<double> residualArray;
      filter.checkProcessJacobians(epsilonArray, controlInput, residualArray);

      // Because the process is actually linear, we expect the
      // residual to be negligible.
      BRICK_TEST_ASSERT(residualArray.rows() == 2);
      BRICK_TEST_ASSERT(residualArray.columns() == 2);
      for(unsigned int ii = 0; ii < residualArray.size(); ++ii) {
        BRICK_TEST_ASSERT(std::fabs(residualArray[ii]) < 1.0E-8);
      }
    }


    ExtendedKalmanFilterTest::LinearProcessFilter
    ExtendedKalmanFilterTest::
    initialize2DFilter(numeric::Array1D<double>& initialState,
                       numeric::Array2D<double>& AMatrix,
                       numeric::Array2D<double>& HMatrix,
                       numeric::Array2D<double>& measurementNoiseArray)
    {
      // Actual start position.
      const numeric::Array1D<double> xStart("[1.0, 0.0]");

      // Initial state estimate.
      initialState = numeric::Array1D<double>("[-0.5, 1.0]");
      const numeric::Array2D<double> initialVariance("[[1.0, 0.0], [0.0, 1.0]]");

      // Process model: rotate by 10 degrees every timestep.
      // 
      //   x_k = A * x_(k-1) + w_k
      AMatrix = numeric::Array2D<double>(
        "[[0.99619469809174555, 0.087155742747658166], "
        " [-0.087155742747658166, 0.99619469809174555]]");
      const numeric::Array2D<double> BMatrix("[[0.0], [0.0]]");
      const numeric::Array2D<double> processCovariance("[[1.0e-5, 0.0], [0.0, 1.0e-5]]");

      // z_k = H * x_k + v_k
      HMatrix = numeric::Array2D<double>("[[-1.0, 1.5], [0.0, 0.5]]");
      const numeric::Array2D<double> measurementCovariance(
        "[[1.0e-2, 0.0], [0.0, 1.0e-2]]");

      // Pre-generate our noise.  For now, we can only tolerate
      // diagonal covariances.
      for(unsigned int row = 0; row < measurementCovariance.rows(); ++row) {
        for(unsigned int column = 0; column < measurementCovariance.columns();
            ++column) {
          if(row != column) {
            if(measurementCovariance(row, column) != 0.0) {
              BRICK_THROW(
                common::LogicException,
                "ExtendedKalmanFilterTest::testExtendedKalmanFilter2D()",
                "Measurement noise covariance must be diagonal.");
            }
          }
        }
      }
      if(measurementCovariance.rows() > m_whiteNoiseArray.rows()) {
        BRICK_THROW(
          common::LogicException,
          "ExtendedKalmanFilterTest::testExtendedKalmanFilter2D()",
          "Member variable m_witeNoiseArray needs more rows.");
      }
      measurementNoiseArray = numeric::Array2D<double>(
        measurementCovariance.rows(), m_whiteNoiseArray.columns());
      for(unsigned int ii = 0; ii < measurementCovariance.rows(); ++ii) {
        double sigma = std::sqrt(measurementCovariance(ii, ii));
        measurementNoiseArray.getRow(ii).copy(
          m_whiteNoiseArray.getRow(ii) * sigma);
      }

      // Set up the filter.
      LinearProcessFilter filter(AMatrix, BMatrix, HMatrix,
                                 processCovariance, measurementCovariance);
      filter.setStateEstimate(0, initialState, initialVariance);
      return filter;
    }
    
  } // namespace computerVision
  
} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::computerVision::ExtendedKalmanFilterTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::computerVision::ExtendedKalmanFilterTest currentTest;

}

#endif
