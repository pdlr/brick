/**
***************************************************************************
* @file brick/linearAlgebra/test/linearAlgebraTest.cc
* Source file defining linearAlgebraTest class.
*
* Copyright (C) 2005-2011 David LaRose, dlr@cs.cmu.edu
* See accompanying file, LICENSE.TXT, for details.
*
***************************************************************************
**/

#include <algorithm>
#include <complex>
#include <brick/common/functional.hh>
#include <brick/linearAlgebra/linearAlgebra.hh>
#include <brick/numeric/utilities.hh>

#include <brick/test/testFixture.hh>

namespace brick {

  namespace linearAlgebra {

    class LinearAlgebraTest : public test::TestFixture<LinearAlgebraTest> {

    public:

      LinearAlgebraTest();
      ~LinearAlgebraTest() {}

      void setUp(const std::string& /* testName */) {}
      void tearDown(const std::string& /* testName */) {}

      void testCholeskyFactorization();
      void testEigenvaluesSymmetric();
      void testEigenvectors();
      void testEigenvectorsSymmetric();
      void testInverse();
      void testLinearLeastSquares();
      void testLinearSolveInPlace();
      void testQrFactorization();
      void testSingularValueDecomposition();
      void testSingularValues();

    private:

      bool
      approximatelyEqual(const numeric::Array1D<common::Float64>& array0,
                         const numeric::Array1D<common::Float64>& array1);

      bool
      approximatelyEqual(const numeric::Array2D<common::Float64>& array0,
                         const numeric::Array2D<common::Float64>& array1);


      bool
      isUpperTriangular(const numeric::Array2D<common::Float64>& array0);


      static const size_t numberOfTestMatrixSets = 2;
    
      std::vector< numeric::Array2D<common::Float64> > m_aMatrices;
      std::vector< numeric::Array1D<common::Float64> > m_bVectors;
      std::vector< numeric::Array2D<common::Float64> > m_inverseMatrices;
      std::vector< numeric::Array1D<common::Float64> > m_squareBVectors;
      std::vector< numeric::Array1D< std::complex<common::Float64> > > m_squareEigenvalues;
      std::vector< numeric::Array2D< std::complex<common::Float64> > > m_squareEigenvectors;
      std::vector< numeric::Array2D<common::Float64> > m_squareMatrices;
      std::vector< numeric::Array1D<common::Float64> > m_squareXVectors;
      std::vector< numeric::Array1D<common::Float64> > m_symmetricEigenvalues;
      std::vector< numeric::Array2D<common::Float64> > m_symmetricEigenvectors;
      std::vector< numeric::Array2D<common::Float64> > m_symmetricMatrices;
      std::vector< numeric::Array2D<common::Float64> > m_upperTriangularMatrices;
      std::vector< numeric::Array1D<common::Float64> > m_sVectors;
      std::vector< numeric::Array2D<common::Float64> > m_testMatrices;
      std::vector< numeric::Array2D<common::Float64> > m_uMatrices;
      std::vector< numeric::Array2D<common::Float64> > m_vtMatrices;
      std::vector< numeric::Array1D<common::Float64> > m_xVectors;
    
    }; // class LinearAlgebraTest


    /* ============== Member Function Definititions ============== */

    LinearAlgebraTest::
    LinearAlgebraTest()
      : brick::test::TestFixture<LinearAlgebraTest>("LinearAlgebraTest"),
        m_aMatrices(LinearAlgebraTest::numberOfTestMatrixSets),
        m_bVectors(LinearAlgebraTest::numberOfTestMatrixSets),
        m_inverseMatrices(LinearAlgebraTest::numberOfTestMatrixSets),
        m_squareBVectors(LinearAlgebraTest::numberOfTestMatrixSets),
        m_squareEigenvalues(LinearAlgebraTest::numberOfTestMatrixSets),
        m_squareEigenvectors(LinearAlgebraTest::numberOfTestMatrixSets),
        m_squareMatrices(LinearAlgebraTest::numberOfTestMatrixSets),
        m_squareXVectors(LinearAlgebraTest::numberOfTestMatrixSets),
        m_symmetricEigenvalues(LinearAlgebraTest::numberOfTestMatrixSets),
        m_symmetricEigenvectors(LinearAlgebraTest::numberOfTestMatrixSets),
        m_symmetricMatrices(LinearAlgebraTest::numberOfTestMatrixSets),
        m_upperTriangularMatrices(LinearAlgebraTest::numberOfTestMatrixSets),
        m_sVectors(LinearAlgebraTest::numberOfTestMatrixSets),
        m_testMatrices(LinearAlgebraTest::numberOfTestMatrixSets),
        m_uMatrices(LinearAlgebraTest::numberOfTestMatrixSets),
        m_vtMatrices(LinearAlgebraTest::numberOfTestMatrixSets),
        m_xVectors(LinearAlgebraTest::numberOfTestMatrixSets)      
    {
      // Register all tests.
      BRICK_TEST_REGISTER_MEMBER(testCholeskyFactorization);
      BRICK_TEST_REGISTER_MEMBER(testEigenvaluesSymmetric);
      BRICK_TEST_REGISTER_MEMBER(testEigenvectors);
      BRICK_TEST_REGISTER_MEMBER(testEigenvectorsSymmetric);
      BRICK_TEST_REGISTER_MEMBER(testInverse);
      BRICK_TEST_REGISTER_MEMBER(testLinearLeastSquares);
      BRICK_TEST_REGISTER_MEMBER(testLinearSolveInPlace);
      BRICK_TEST_REGISTER_MEMBER(testQrFactorization);
      BRICK_TEST_REGISTER_MEMBER(testSingularValueDecomposition);
      BRICK_TEST_REGISTER_MEMBER(testSingularValues);

      // Set up test matrices.
      if(LinearAlgebraTest::numberOfTestMatrixSets != 2) {
        BRICK_THROW(common::LogicException,
                    "LinearAlgebraTest::LinearAlgebraTest()",
                    "Incorrect number of test sets.");
      }

      m_aMatrices[0] = numeric::Array2D<common::Float64>(
        "[[1.0, 2.0, 3.0],"
        " [4.0, 5.0, 6.0],"
        " [10.0, 3.0, 1.0],"
        " [1.0, 5.0, 9.0],"
        " [6.0, 6.0, 6.0]]");
      
      m_aMatrices[1] = numeric::Array2D<common::Float64>(
        "[[1.0, 2.0, 3.0, 4.0],"
        " [4.0, 2.0, 6.0, 1.0],"
        " [10.0, 3.0, 1.0, 1.0]]");

    
      m_bVectors[0] = numeric::Array1D<common::Float64>(
        "[14.0, 32.0, 19.0, 38.0, 36.0]");

      m_bVectors[1] = numeric::Array1D<common::Float64>(
        "[23.5, 20.5, 21.0]");

    
      m_inverseMatrices[0] = numeric::Array2D<common::Float64>(
        "[[-1.47619048,  0.61904762, -0.14285714],"
        " [ 1.28571429, -0.57142857,  0.28571429],"
        " [-0.14285714,  0.28571429, -0.14285714]]");

      m_inverseMatrices[1] = numeric::Array2D<common::Float64>(
        "[[ 0.06666667,  0.33333333, -0.06666667, -2.73333333],"
        " [ 1.26666667, -0.66666667, -0.26666667, -1.93333333],"
        " [-1.53333333,  0.33333333,  0.53333333,  4.86666667],"
        " [ 0.        ,  0.        ,  0.        ,  1.        ]]");


      m_squareBVectors[0] = numeric::Array1D<common::Float64>(
        "[4., 2., 7.]");

      m_squareBVectors[1] = numeric::Array1D<common::Float64>(
        "[1., 2., 7., 4.]");

      m_squareEigenvalues[0] = numeric::Array1D< std::complex<common::Float64> >(
        "[(9.87406, 0.0), (-4.38955, 0.0), (-0.48451, 0.0)]");

      m_squareEigenvalues[1] = numeric::Array1D< std::complex<common::Float64> >(
        "[(11.80772, 0.0), (1.22733, 0.0), (-1.03505, 0.0), (1.00000, 0.0)]");

      m_squareEigenvectors[0] = numeric::Array2D< std::complex<common::Float64> >(
        "[[(-0.20982, 0.0), (-0.30042, 0.0), (-0.73405, 0.0)],"
        " [(-0.69368, 0.0), (-0.40674, 0.0), ( 0.66163, 0.0)],"
        " [(-0.68905, 0.0), ( 0.86273, 0.0), (-0.15299, 0.0)]]");

      m_squareEigenvectors[1] = numeric::Array2D< std::complex<common::Float64> >(
        "[[(-0.31529, 0.0), (-0.24503, 0.0), ( 0.29668, 0.0), (-0.27156, 0.0)],"
        " [(-0.20248, 0.0), (-0.37076, 0.0), (-0.82200, 0.0), (-0.37717, 0.0)],"
        " [(-0.92714, 0.0), ( 0.89582, 0.0), ( 0.48611, 0.0), ( 0.88508, 0.0)],"
        " [( 0.00000, 0.0), ( 0.00000, 0.0), ( 0.00000, 0.0), ( 0.02514, 0.0)]]"
        );

      m_squareMatrices[0] = numeric::Array2D<common::Float64>(
        "[[ 0.,  1.,  2.],"
        " [ 3.,  4.,  5.],"
        " [ 6.,  7.,  1.]]");
      
      m_squareMatrices[1] = numeric::Array2D<common::Float64>(
        "[[ 4.,  3.,  2.,  7.],"
        " [ 4.,  1.,  1.,  8.],"
        " [ 9.,  8.,  7.,  6.],"
        " [ 0.,  0.,  0.,  1.]]");


      m_squareXVectors[0] = numeric::Array1D<common::Float64>(
        "[-5.66666667, 6., -1.]");


      m_squareXVectors[1] = numeric::Array1D<common::Float64>(
        "[-10.66666667, -9.66666667, 22.33333333, 4.]");

    
      m_symmetricEigenvalues[0] = numeric::Array1D<common::Float64>(
        "[24.06253512, 0.5580362, 0.18491359, -0.80548492]");

      m_symmetricEigenvalues[1] = numeric::Array1D<common::Float64>(
        "[22.25237164, 10.12484443, -2.53641946, -12.84079662]");


      m_symmetricEigenvectors[0] = numeric::Array2D<common::Float64>(
        "[[ -0.22593827, 0.56047033, -0.32976505, -0.72531367],"
        " [ -0.44322186, -0.77633831, -0.3153256, -0.31846973],"
        " [ -0.57278788, 0.05844028, 0.805111, -0.14246073],"
        " [ -0.6514755, 0.28241204, -0.37897369, 0.59346614]]");

      m_symmetricEigenvectors[1] = numeric::Array2D<common::Float64>(
        "[[ 0.57388773, -0.40927089, 0.65613403, -0.26951503],"
        " [ 0.65948776, -0.25973728, -0.69857373, 0.09801622],"
        " [ -0.48540238, -0.84196315, -0.16831305, -0.16478259],"
        " [ 0.01064409, -0.2369218, 0.23055064, 0.94371668]]");

    
      m_symmetricMatrices[0] = numeric::Array2D<common::Float64>(
        "[[  1.,   2.,   3.,   4.],"
        " [  2.,   5.,   6.,   7.],"
        " [  3.,   6.,   8.,   9.],"
        " [  4.,   7.,   9.,  10.]]");

      m_symmetricMatrices[1] = numeric::Array2D<common::Float64>(
        "[[  7.,  11.,  -3.,   4.],"
        " [ 11.,   9.,  -5.,   0.],"
        " [ -3.,  -5.,  12.,   4.],"
        " [  4.,   0.,   4., -11.]]");
    

      m_upperTriangularMatrices[0] = numeric::Array2D<common::Float64>(
        "[[  1.,   2.,   3.,   4.],"
        " [  0.,   5.,   6.,   7.],"
        " [  0.,   0.,   8.,   9.],"
        " [  0.,   0.,   0.,  10.]]");

      m_upperTriangularMatrices[1] = numeric::Array2D<common::Float64>(
        "[[  7.,  11.,  -3.,   4.],"
        " [  0.,   9.,  -5.,   0.],"
        " [  0.,   0.,  12.,   4.],"
        " [  0.,   0.,   0., -11.]]");
    

      m_sVectors[0] = numeric::Array1D<common::Float64>(
        "[ 18.52607955,   6.51255861,   0.60906229]");

      m_sVectors[1] = numeric::Array1D<common::Float64>(
        "[ 18.67643567,   6.03617265,   0.86912045]");

    
      m_testMatrices[0] = numeric::Array2D<common::Float64>(
        "[[  0.,   1.,   2.],"
        " [  3.,   4.,   5.],"
        " [  6.,   7.,   8.],"
        " [  9.,  10.,  1.]]");

      m_testMatrices[1] = numeric::Array2D<common::Float64>(
        "[[  0.,   1.,   2.,   3.],"
        " [  4.,   5.,   6.,   7.],"
        " [  8.,   9.,  10.,  1.]]");


      m_uMatrices[0] = numeric::Array2D<common::Float64>(
        "[[-0.08198544, -0.23824706,  0.877411  ],"
        " [-0.35697841, -0.38359746,  0.24238963],"
        " [-0.63197137, -0.52894786, -0.39263173],"
        " [-0.68297656,  0.718544  ,  0.13129179]]");

      m_uMatrices[1] = numeric::Array2D<common::Float64>(
        "[[-0.14130401,  0.42003822, -0.89643799],"
        " [-0.55423802,  0.71674227,  0.42320295],"
        " [-0.82027642, -0.55664029, -0.13152259]]");

    
      m_vtMatrices[0] = numeric::Array2D<common::Float64>(
        "[[-0.5942732 , -0.68894578, -0.41496154],"
        " [ 0.32896694,  0.26259542, -0.90709669],"
        " [-0.73390743,  0.67557188, -0.07058696]]");

      m_vtMatrices[1] = numeric::Array2D<common::Float64>(
        "[[-0.47006632, -0.5512284 , -0.63239049, -0.27434863],"
        " [-0.26277466, -0.16666405, -0.07055345,  0.94773139],"
        " [ 0.73710275,  0.04127557, -0.65455161,  0.16290502]]");

    
      m_xVectors[0] = numeric::Array1D<common::Float64>(
        "[1.0, 2.0, 3.0]");

      m_xVectors[1] = numeric::Array1D<common::Float64>(
        "[1.070794, 1.72098834, 1.52915047, 3.59994448]");

    }


    void
    LinearAlgebraTest::
    testCholeskyFactorization()
    {
      for(size_t index0 = 0;
          index0 < LinearAlgebraTest::numberOfTestMatrixSets;
          ++index0) {

        // Construct a symmetric positive definite matrix for which we
        // know the Cholesky factorization.
        numeric::Array2D<common::Float64> aMatrix =
          numeric::matrixMultiply<common::Float64>(
            (m_upperTriangularMatrices[index0]).transpose(),
            m_upperTriangularMatrices[index0]);
        numeric::Array2D<common::Float64> aMatrixCopy = aMatrix.copy();
      
      
        // Do the computation once for upper triangular.
        numeric::Array2D<common::Float64> kMatrix;
        choleskyFactorization(aMatrix, kMatrix, true);

        // Since the Cholesky factorization is not unique, we can't
        // just compare kMatrix with m_upperTriangularMatrices[index0].
        // Instead we reconstruct aMatrix using kMatrix.
        numeric::Array2D<common::Float64> kTk = numeric::matrixMultiply<common::Float64>(
          kMatrix.transpose(), kMatrix);
        
        // Check that input is unchanged.
        BRICK_TEST_ASSERT(
          this->approximatelyEqual(aMatrix, aMatrixCopy));

        // Check that results are correct.
        BRICK_TEST_ASSERT(this->isUpperTriangular(kMatrix));
        BRICK_TEST_ASSERT(this->approximatelyEqual(kTk, aMatrix));

        // Do the computation once for lower triangular.
        choleskyFactorization(aMatrix, kMatrix, false);
        numeric::Array2D<common::Float64> kkT = numeric::matrixMultiply<common::Float64>(
          kMatrix, kMatrix.transpose());
        
        // Check that input is unchanged.
        BRICK_TEST_ASSERT(
          this->approximatelyEqual(aMatrix, aMatrixCopy));

        // Check that results are correct.
        BRICK_TEST_ASSERT(this->isUpperTriangular(kMatrix.transpose()));
        BRICK_TEST_ASSERT(this->approximatelyEqual(kkT, aMatrix));
      }
    }


    void
    LinearAlgebraTest::
    testEigenvaluesSymmetric()
    {
      for(size_t index0 = 0;
          index0 < LinearAlgebraTest::numberOfTestMatrixSets;
          ++index0) {
        // Do the computation.
        numeric::Array2D<common::Float64> aMatrix = m_symmetricMatrices[index0].copy();
        numeric::Array1D<common::Float64> eigenvalueArray = eigenvaluesSymmetric(aMatrix);

        // Check that input is unchanged.
        BRICK_TEST_ASSERT(
          this->approximatelyEqual(aMatrix, m_symmetricMatrices[index0]));

        // Check that results are correct.
        BRICK_TEST_ASSERT(
          this->approximatelyEqual(
            eigenvalueArray, m_symmetricEigenvalues[index0]));
      }
    }

  
    void
    LinearAlgebraTest::
    testEigenvectors()
    {
      for(size_t index0 = 0;
          index0 < LinearAlgebraTest::numberOfTestMatrixSets;
          ++index0) {
        // Do the computation.
        numeric::Array2D<common::Float64> aMatrix = m_squareMatrices[index0].copy();
        numeric::Array1D< std::complex<common::Float64> > eigenvalueArray;
        numeric::Array2D< std::complex<common::Float64> > eigenvectorArray;
        eigenvectors(aMatrix, eigenvalueArray, eigenvectorArray, true);

        // Check that input is unchanged.
        BRICK_TEST_ASSERT(
          this->approximatelyEqual(aMatrix, m_squareMatrices[index0]));

        // Check that results are correct.
        BRICK_TEST_ASSERT(
          eigenvalueArray.size() == m_squareEigenvalues[index0].size());
        BRICK_TEST_ASSERT(
          eigenvectorArray.rows() == m_squareEigenvectors[index0].rows());
        BRICK_TEST_ASSERT(
          eigenvectorArray.columns() == m_squareEigenvectors[index0].columns());

        // Check eigenvalues.
        for(size_t ii = 0; ii < eigenvalueArray.size(); ++ii) {
          BRICK_TEST_ASSERT(
            brick::approximatelyEqual(
              eigenvalueArray[ii].real(),
              m_squareEigenvalues[index0][ii].real(), 1.0E-4));
          BRICK_TEST_ASSERT(
            brick::approximatelyEqual(
              eigenvalueArray[ii].imag(),
              m_squareEigenvalues[index0][ii].imag(), 1.0E-4));
        }

        // Remember that -1 times an eigenvector is still an
        // eigenvector.  We correct any sign errors here.
        for(size_t ii = 0; ii < eigenvectorArray.columns(); ++ii) {
          // This works because we know (from checking) that the first
          // element of all of our eigenvectors has nonzero real part.
          if((eigenvectorArray(0, ii).real()
              * m_squareEigenvectors[index0](0, ii).real())
             < 0.0) {
            for(size_t jj = 0; jj < eigenvectorArray.rows(); ++jj) {
              eigenvectorArray(jj, ii) *= -1.0;
            }
          }
        }
            

        // Check eigenvectors
        for(size_t ii = 0; ii < eigenvectorArray.size(); ++ii) {
          BRICK_TEST_ASSERT(
            brick::approximatelyEqual(
              eigenvectorArray[ii].real(),
              m_squareEigenvectors[index0][ii].real(), 1.0E-4));
          BRICK_TEST_ASSERT(
            brick::approximatelyEqual(
              eigenvectorArray[ii].imag(),
              m_squareEigenvectors[index0][ii].imag(), 1.0E-4));
        }
      }
    }
    

    void
    LinearAlgebraTest::
    testEigenvectorsSymmetric()
    {
      for(size_t index0 = 0;
          index0 < LinearAlgebraTest::numberOfTestMatrixSets;
          ++index0) {
        // Do the computation.
        numeric::Array2D<common::Float64> aMatrix = m_symmetricMatrices[index0].copy();
        numeric::Array1D<common::Float64> eigenvalueArray;
        numeric::Array2D<common::Float64> eigenvectorArray;
        eigenvectorsSymmetric(aMatrix, eigenvalueArray, eigenvectorArray);

        // Check that input is unchanged.
        BRICK_TEST_ASSERT(
          this->approximatelyEqual(aMatrix, m_symmetricMatrices[index0]));

        // Check that results are correct.
        BRICK_TEST_ASSERT(
          this->approximatelyEqual(
            eigenvalueArray, m_symmetricEigenvalues[index0]));
        BRICK_TEST_ASSERT(
          this->approximatelyEqual(
            eigenvectorArray, m_symmetricEigenvectors[index0]));
      }
    }


    void
    LinearAlgebraTest::
    testInverse()
    {
      for(size_t index0 = 0;
          index0 < LinearAlgebraTest::numberOfTestMatrixSets;
          ++index0) {
        // Do the inverse.
        numeric::Array2D<common::Float64> squareMatrix = m_squareMatrices[index0].copy();
        numeric::Array2D<common::Float64> inverseMatrix = inverse(squareMatrix);

        // Check that input is unchanged.
        BRICK_TEST_ASSERT(
          this->approximatelyEqual(squareMatrix, m_squareMatrices[index0]));

        // Check that results are correct.
        BRICK_TEST_ASSERT(
          this->approximatelyEqual(inverseMatrix, m_inverseMatrices[index0]));
      }
    }


    void
    LinearAlgebraTest::
    testLinearLeastSquares()
    {
      for(size_t index0 = 0;
          index0 < LinearAlgebraTest::numberOfTestMatrixSets;
          ++index0) {
        // Compute the solution.
        numeric::Array2D<common::Float64> aMatrix = m_aMatrices[index0].copy();
        numeric::Array1D<common::Float64> bVector = m_bVectors[index0].copy();
        numeric::Array1D<common::Float64> xVector = linearLeastSquares(aMatrix, bVector);
      
        // Check that input is unchanged.
        BRICK_TEST_ASSERT(
          this->approximatelyEqual(aMatrix, m_aMatrices[index0]));
        BRICK_TEST_ASSERT(
          this->approximatelyEqual(bVector, m_bVectors[index0]));

        // Check that results are correct.
        BRICK_TEST_ASSERT(
          this->approximatelyEqual(xVector, m_xVectors[index0]));
      }
    }
  

    void
    LinearAlgebraTest::
    testLinearSolveInPlace()
    {
      for(size_t index0 = 0;
          index0 < LinearAlgebraTest::numberOfTestMatrixSets;
          ++index0) {
        // Compute the solution.
        numeric::Array2D<common::Float64> aMatrix = m_squareMatrices[index0].copy();
        numeric::Array1D<common::Float64> bVector = m_squareBVectors[index0].copy();
        linearSolveInPlace(aMatrix, bVector);
      
        // Check that results are correct.
        BRICK_TEST_ASSERT(
          this->approximatelyEqual(bVector, m_squareXVectors[index0]));
      }
    }
  

    void
    LinearAlgebraTest::
    testQrFactorization()
    {
      // Create some test matrices.
      std::vector< numeric::Array2D<common::Float64> > testMatrices;
      for(size_t index0 = 0;
          index0 < LinearAlgebraTest::numberOfTestMatrixSets;
          ++index0) {
        testMatrices.push_back(m_testMatrices[index0]);
        testMatrices.push_back((m_testMatrices[index0]).transpose());
      }

      for(size_t index0 = 0; index0 < testMatrices.size(); ++index0) {

        // Start with a general matrix
        numeric::Array2D<common::Float64> aMatrix = testMatrices[index0];
        numeric::Array2D<common::Float64> aMatrixCopy = aMatrix.copy();
      
        // Do the factorization.
        numeric::Array2D<common::Float64> qMatrix;
        numeric::Array2D<common::Float64> rMatrix;
        qrFactorization(aMatrix, qMatrix, rMatrix);

        // Since q is orthonormal, q^T * q == I.  Reconstruct the
        // identity matrix, and make a reference identity against
        // which to compare.
        numeric::Array2D<common::Float64> qTq = numeric::matrixMultiply<common::Float64>(
          qMatrix.transpose(), qMatrix);
        
        // numeric::Array2D<common::Float64> iMatrix = num::identity<common::Float64>(
        //   aMatrix.rows(), aMatrix.rows());
        numeric::Array2D<common::Float64> iMatrix =
          numeric::identity<common::Float64>(aMatrix.rows(), aMatrix.rows());

        // Reconstruct aMatrix using qMatrix and rMatrix.
        numeric::Array2D<common::Float64> qrProduct =
          numeric::matrixMultiply<common::Float64>(qMatrix, rMatrix);

        // Check the sizes of the returned matrices.
        BRICK_TEST_ASSERT(qMatrix.rows() == aMatrix.rows());
        BRICK_TEST_ASSERT(qMatrix.columns() == aMatrix.rows());
        BRICK_TEST_ASSERT(rMatrix.rows() == rMatrix.rows());
        BRICK_TEST_ASSERT(rMatrix.columns() == rMatrix.columns());
        
        // Check that rMatrix is upper triangular with non-negative diagonal.
        BRICK_TEST_ASSERT(this->isUpperTriangular(rMatrix));
        for(size_t ii = 0; ii < std::min(rMatrix.rows(), rMatrix.columns());
            ++ii) {
          BRICK_TEST_ASSERT(rMatrix(ii, ii) >= 0.0);
        }
        
        // Check that input is unchanged.
        BRICK_TEST_ASSERT(
          this->approximatelyEqual(aMatrix, aMatrixCopy));

        // Check that results are correct.
        BRICK_TEST_ASSERT(this->approximatelyEqual(qTq, iMatrix));
        BRICK_TEST_ASSERT(this->approximatelyEqual(qrProduct, aMatrix));
      }
    }


    void
    LinearAlgebraTest::
    testSingularValueDecomposition()
    {
      numeric::Array2D<common::Float64> aMatrix;
      numeric::Array2D<common::Float64> uMatrix;
      numeric::Array1D<common::Float64> sVector;
      numeric::Array2D<common::Float64> vtMatrix;

      for(size_t index0 = 0;
          index0 < LinearAlgebraTest::numberOfTestMatrixSets;
          ++index0) {
        // Do the SVD.
        aMatrix = m_testMatrices[index0].copy();
        singularValueDecomposition(aMatrix, uMatrix, sVector, vtMatrix);

        // Check that input is unchanged.
        BRICK_TEST_ASSERT(
          this->approximatelyEqual(aMatrix, m_testMatrices[index0]));

        // Check that results are correct.
        BRICK_TEST_ASSERT(this->approximatelyEqual(uMatrix, m_uMatrices[index0]));
        BRICK_TEST_ASSERT(this->approximatelyEqual(sVector, m_sVectors[index0]));
        BRICK_TEST_ASSERT(
          this->approximatelyEqual(vtMatrix, m_vtMatrices[index0]));
      }
    }


    void
    LinearAlgebraTest::
    testSingularValues()
    {
      for(size_t index0 = 0;
          index0 < LinearAlgebraTest::numberOfTestMatrixSets;
          ++index0) {
        // Do the SVD.
        numeric::Array2D<common::Float64> aMatrix = m_testMatrices[index0].copy();
        numeric::Array1D<common::Float64> sVector = singularValues(aMatrix);

        // Check that input is unchanged.
        BRICK_TEST_ASSERT(
          this->approximatelyEqual(aMatrix, m_testMatrices[index0]));

        // Check that results are correct.
        BRICK_TEST_ASSERT(this->approximatelyEqual(sVector, m_sVectors[index0]));
      }
    }


    bool
    LinearAlgebraTest::
    approximatelyEqual(const numeric::Array1D<common::Float64>& array0,
                       const numeric::Array1D<common::Float64>& array1)
    {
      if(array0.size() != array1.size()) {
        return false;
      }
      return std::equal(array0.begin(), array0.end(), array1.begin(),
                        ApproximatelyEqualFunctor<common::Float64>(1.0E-5));
    }

  
    bool
    LinearAlgebraTest::
    approximatelyEqual(const numeric::Array2D<common::Float64>& array0,
                       const numeric::Array2D<common::Float64>& array1)
    {
      if(array0.rows() != array1.rows()) {
        return false;
      }
      if(array0.columns() != array1.columns()) {
        return false;
      }
      return std::equal(array0.begin(), array0.end(), array1.begin(),
                        ApproximatelyEqualFunctor<common::Float64>(1.0E-5));
    }

    
    bool
    LinearAlgebraTest::
    isUpperTriangular(const numeric::Array2D<common::Float64>& array0)
    {
      for(size_t rowIndex = 0; rowIndex < array0.rows(); ++rowIndex) {
        for(size_t columnIndex = 0; columnIndex < rowIndex; ++columnIndex) {
          if(array0(rowIndex, columnIndex) != 0.0) {
            return false;
          }
        }
      }
      return true;
    }
    

  } // namespace linearAlgebra
    
} // namespace brick


#if 0

int main(int argc, char** argv)
{
  brick::LinearAlgebraTest currentTest;
  bool result = currentTest.run();
  return (result ? 0 : 1);
}

#else

namespace {

  brick::linearAlgebra::LinearAlgebraTest currentTest;

}

#endif
