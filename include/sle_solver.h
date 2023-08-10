#ifndef ATG_SIMPLE_2D_CONSTRAINT_SOLVER_SLE_SOLVER_H
#define ATG_SIMPLE_2D_CONSTRAINT_SOLVER_SLE_SOLVER_H

#include "matrix.h"
#include "matrix_sparse.h"

namespace atg_scs {
    // Simple Linear Equation Solver
    class SleSolver {
        public:
            SleSolver(bool supportsLimits);
            virtual ~SleSolver();
            
            // Solves for (x) J * X + W = right
            virtual bool solve(
                    SparseMatrix<3> &J,
                    Matrix &W,
                    Matrix &right,
                    Matrix *result,
                    Matrix *previous);
            virtual bool solveWithLimits(
                    SparseMatrix<3> &J,
                    Matrix &W,
                    Matrix &right,
                    Matrix &limits,
                    Matrix *result,
                    Matrix *previous);

            bool supportsLimits() const { return m_supportsLimits; }

        private:
            bool m_supportsLimits;
    };
} /* namespace atg_scs */

#endif /* ATG_SIMPLE_2D_CONSTRAINT_SOLVER_SLE_SOLVER_H */
