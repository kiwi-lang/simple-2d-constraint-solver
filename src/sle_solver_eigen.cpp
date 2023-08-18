#include "../include/sle_eigen.h"

atg_scs::SleEigenSolver::SleEigenSolver(bool supportsLimits) {
    m_supportsLimits = supportsLimits;
}

atg_scs::SleEigenSolver::~SleEigenSolver() {
    /* void */
}

bool atg_scs::SleEigenSolver::solve(
        SparseMatrix<3> &J,
        Matrix &W,
        Matrix &right,
        Matrix *result,
        Matrix *previous)
{
    return false;
}

bool atg_scs::SleEigenSolver::solveWithLimits(
        SparseMatrix<3> &J,
        Matrix &W,
        Matrix &right,
        Matrix &limits,
        Matrix *result,
        Matrix *previous)
{
    return false;
}
