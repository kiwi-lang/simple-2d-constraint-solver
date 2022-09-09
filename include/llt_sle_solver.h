#ifndef ATG_SIMPLE_2D_CONSTRAINT_SOLVER_CONJUGATE_GRADIENT_SLE_SOLVER_H
#define ATG_SIMPLE_2D_CONSTRAINT_SOLVER_CONJUGATE_GRADIENT_SLE_SOLVER_H

#include "sle_solver.h"

namespace atg_scs
{
    class LLTSolver : public SleSolver
    {
    public:
        LLTSolver();
        virtual ~LLTSolver();

        virtual bool solve(
            SparseMatrix<3> &J,
            Matrix &W,
            Matrix &right,
            Matrix *result,
            Matrix *previous);

    protected:
        void multiply(SparseMatrix<3> &J, Matrix &W, Matrix &x, Matrix *target);
        bool sufficientlySmall(Matrix &x, Matrix &target) const;

        Matrix
            m_mreg0,
            m_mreg1,
            m_Ap,
            m_x,
            m_r,
            m_p,
            m_A;

        int m_maxIterations;
        double m_maxError;
        double m_minError;
    };
} /* namespace atg_scs */

#endif /* ATG_SIMPLE_2D_CONSTRAINT_SOLVER_CONJUGATE_GRADIENT_SLE_SOLVER_H */
