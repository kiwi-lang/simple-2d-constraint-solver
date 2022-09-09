#include "../include/llt_sle_solver.h"

#include <cmath>
#include <assert.h>

atg_scs::LLTSolver::LLTSolver()
    : atg_scs::SleSolver(false)
{
}

atg_scs::LLTSolver::~LLTSolver()
{
    m_mreg0.destroy();
    m_mreg1.destroy();
    m_Ap.destroy();
    m_x.destroy();
    m_r.destroy();
    m_p.destroy();
}

bool atg_scs::LLTSolver::solve(
    SparseMatrix<3> &J,
    Matrix &W,
    Matrix &right,
    Matrix *previous,
    Matrix *result)
{
    const int n = right.getHeight();

    // A = J * W * J_T
    Eigen::MatrixXd A = (J.m_matrix) * (W.m_matrix) * (J.m_matrix.transpose());

    result->m_matrix = A.llt().solve(right);
    return true;
}
