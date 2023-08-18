#include "../include/sle_solver_gauss_seidel.h"

#include "../../../Source/RegressionHelper.h"
#include "../../eigen/Eigen/IterativeLinearSolvers"


#include <cmath>
#include <assert.h>

atg_scs::GaussSeidelSleSolver::GaussSeidelSleSolver()
    : atg_scs::SleSolver(true)
{
    m_maxIterations = 128;
    m_minDelta = 1E-1;

    m_M.initialize(1, 1);

    Eigen::initParallel();
}

atg_scs::GaussSeidelSleSolver::~GaussSeidelSleSolver() {
    m_M.destroy();
    m_reg.destroy();
}

bool atg_scs::GaussSeidelSleSolver::solve(
        SparseMatrix<3> &J, // J_sparse (3n x m_f)
        Matrix &W,          // M_inv   
        Matrix &right,      // right
        Matrix *previous,
        Matrix *result)
{
    const int n = right.getHeight();

    result->resize(1, n);

    if (previous != nullptr && previous->getHeight() == n) {
        result->set(previous);
    }

    J.rightScale(W, &m_reg);
    m_reg.multiplyTranspose(J, &m_M);

    for (int i = 0; i < m_maxIterations; ++i) {
        const double maxDelta = solveIteration(
                m_M,
                right,
                result, 
                result);
 
        if (maxDelta < m_minDelta) {
            return true; 
        }
    }

    return false;
}

bool atg_scs::GaussSeidelSleSolver::solveWithLimits(
    SparseMatrix<3> &J,
    Matrix &W,
    Matrix &right,
    Matrix &limits,
    Matrix *result,
    Matrix *previous)
{
    const int n = right.getHeight();
    if (result->getHeight() != n) {
        result->initialize(1, n);
    }

    if (previous != nullptr && previous->getHeight() == n) {
        result->set(previous);
    }

    // J: 33 x 36
    // W: 36 x 1
    // (J .* W) * J^T 
    J.rightScale(W, &m_reg);
    m_reg.multiplyTranspose(J, &m_M);

    /*

    // Solves Ax = b
    Eigen::ConjugateGradient<Eigen::MatrixXf> cg;
    cg.setMaxIterations(m_maxIterations);
    cg.compute(m_M.m_matrix);
    cg.setTolerance(m_minDelta);
    // result->m_matrix = cg.solveWithGuess(right.m_matrix, previous->m_matrix);
    result->m_matrix = cg.solve(right.m_matrix);
    return true;
    //*/

    // M: (33 x 33)
    // Solve: M * x = right (33 x 1)
    for (int i = 0; i < m_maxIterations; ++i) {
        const double maxDelta = solveIteration(
            m_M,
            right,
            limits,
            result,
            result);

        if (maxDelta < m_minDelta) {
            return true;
        }
    }

    return false;
}

double atg_scs::GaussSeidelSleSolver::solveIteration(
        Matrix &left,
        Matrix &right,
        Matrix *k_next,
        Matrix *k)
{
    double maxDifference = 0.0;
    const int n = k->getHeight();

    for (int i = 0; i < n; ++i) {

        // Dot product
        #ifdef ATG_S2C_USE_VECTORIZE
        #define M(x) ((x).m_matrix)
        double s0 = (
            M(left).row(i).segment(0, i) *
            M(*k_next).col(0).segment(0, i)
        ).sum();
        
       double s1 = (
            M(left).row(i).segment(i + 1, n - i - 1) * 
            M(*k).col(0).segment(i + 1, n - i - 1)
        ).sum();

        #else
        double s0 = 0.0, s1 = 0.0;
        for (int j = 0; j < i; ++j) {
            s0 += left.get(j, i) * k_next->get(0, j);
        }

        for (int j = i + 1; j < n; ++j) {
            s1 += left.get(j, i) * k->get(0, j);
        }
        #endif

        const double k_next_i =
            (1 / left.get(i, i)) * (right.get(0, i) - s0 - s1);

        const double min_k = 1E-3 > k->get(0, i) ? 1E-3 : k->get(0, i);
        const double delta = (std::abs(k_next_i) - min_k) / min_k;
        maxDifference = (delta > maxDifference)
            ? delta
            : maxDifference;

        k_next->set(0, i, k_next_i);
    }

    return maxDifference;
}

template<typename T>
T cond(int c, T a, T b) {
    return c * a + (1 - c) * b;
}

template<typename T>
T fastmax(const T& a, const T& b) {
    if (a > b) {
        return a;
    }
    return b;
    return std::max(a, b);
    return cond(a > b, a, b);
}

template<typename T>
T fastmin(const T& a, const T& b) {
    if (a < b) {
        return a;
    }
    return b;
    return std::min(a, b);
    return cond(a < b, a, b);
}

template<typename T>
T fastabs(const T& a) {
    if (a > 0) {
        return a;
    }
    return -a;
    return std::abs(a);
    return cond(a > 0, a, -a);
}


double atg_scs::GaussSeidelSleSolver::solveIteration(
        Matrix &left,
        Matrix &right,
        Matrix &limits,
        Matrix *k_next,
        Matrix *k)
{
#define ATG_S2C_USE_JACOBI
#ifdef ATG_S2C_USE_JACOBI
    using Mat = Eigen::MatrixXf;

    // This is the jacobi method now not GaussSeidel
    // GaussSeidel converges faster but jacobi is more easily vectorized
    Mat& A       = left.m_matrix;
    Mat& B       = right.m_matrix;
    Mat& MLimits = limits.m_matrix;
    Mat& Mk      = k->m_matrix;
    Mat& Mk_next = k_next->m_matrix;

    // Mk == Mk_next but this does not work for the jacobi method
    // so we make a copy that we will update later
    static Mat Previous = Mk;

    const int rows = A.rows();
    const int cols = A.cols();

    static Mat Mask = (1 - Mat::Identity(rows, cols).array()).matrix();
    
    // Mask: 33x33
    // k   : 33x1
    //   S : 33x1
    Mat S = ((A.array() * Mask.array()).matrix() * Mk);

    //      (bi       -   ij * xj) / a ii    where i != j
    Mat u = (B.array() - S.array()) / A.diagonal().array();

    Previous = Mk_next;
    Mk_next = u.cwiseMax(MLimits.col(0)).cwiseMin(MLimits.col(1));

    return ((Mk_next - Previous).cwiseAbs().array() / Previous.cwiseAbs().cwiseMax(1e-3).array()).maxCoeff();

#else
    double maxDifference = 0.0;
    const int n = k->getHeight();

    for (int i = 0; i < n; ++i) {
        
        // Dot product
        #ifdef ATG_S2C_USE_VECTORIZE
        #define M(x) ((x).m_matrix)
        double s0 = (
            M(left).row(i).segment(0, i) *
            M(*k_next).col(0).segment(0, i)
        ).sum();
        
       double s1 = (
            M(left).row(i).segment(i + 1, n - i - 1) * 
            M(*k).col(0).segment(i + 1, n - i - 1)
        ).sum();

        #else
        double s0 = 0.0, s1 = 0.0;
        for (int j = 0; j < i; ++j) {
            s0 += left.get(j, i) * k_next->get(0, j); // [0:i, i]  k_next[0, 0:i]
        }
        for (int j = i + 1; j < n; ++j) {
            s1 += left.get(j, i) * k->get(0, j);    //    [i + 1:n, i] k[0, i + 1:n]
        }
        #endif

        // regression::fprint("gaus_seidel", s0, s1);

        const double k_next_i =(1 / left.get(i, i)) * (right.get(0, i) - s0 - s1);
        const double limitMin = limits.get(0, i);
        const double limitMax = limits.get(1, i);
        const double x = fastmax(limitMin, fastmin(limitMax, k_next_i));

        const double min_k = fastmax(1E-3, fastabs(k->get(0, i)));
        const double delta = fastabs(x - k->get(0, i)) / min_k;
        maxDifference = (delta > maxDifference)
            ? delta
            : maxDifference;

        k_next->set(0, i, x);
    }

    return maxDifference;
#endif
}
