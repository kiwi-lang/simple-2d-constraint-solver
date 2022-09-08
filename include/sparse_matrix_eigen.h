#ifndef ATG_SIMPLE_2D_CONSTRAINT_SOLVER_SPARSE_EIGEN_MATRIX_H
#define ATG_SIMPLE_2D_CONSTRAINT_SOLVER_SPARSE_EIGEN_MATRIX_H

#include "matrix_eigen.h"
#include "Eigen/Core"

namespace atg_scs
{
    template <int T_Stride = 3, int T_Entries = 2>
    class SparseMatrix
    {
    public:
        void initialize(int width, int height)
        {
            resize(width, height);
        }

        void dump()
        {
            Eigen::IOFormat fmt;
            std::cout << m_matrix << std::endl;
        }

        void resize(int width, int height)
        {
            m_matrix = Eigen::MatrixXd::Zero(height, width);
        }

        void destroy()
        {
        }

        void expand(Matrix *matrix)
        {
            matrix->m_matrix = m_matrix;
        }

        void expandTransposed(Matrix *matrix)
        {
            matrix->m_matrix = m_matrix.transpose();
        }

        int getWidth() const { return int(m_matrix.cols()); }
        int getHeight() const { return int(m_matrix.rows()); }

        inline void setBlock(int row, int entry, uint8_t index)
        {
            // this is problematic, it can move blocks everywhere
        }

        inline void set(int row, int entry, int slice, double v)
        {
            assert(row >= 0 && row < getHeight());
            assert(entry >= 0 && entry < T_Entries);
            assert(slice < T_Stride);

            m_matrix(row, entry * T_Stride + slice) = v;
        }

        inline double get(int row, int entry, int slice)
        {
            assert(row >= 0 && row < getHeight());
            assert(entry >= 0 && entry < T_Entries);
            assert(slice < T_Stride);

            return m_matrix(row, entry * T_Stride + slice);
        }

        inline void setEmpty(int row, int col)
        {
            assert(row >= 0 && row < getHeight());
            assert(col >= 0 && col < T_Entries);

            // Removes a block
            m_matrix.block(row, col * T_Stride, 0, T_Stride) = 0;
        }

        void multiply(Matrix &b, Matrix *target) const
        {
            target->m_matrix = m_matrix * b.m_matrix;
        }

        void multiplyTranspose(const SparseMatrix<T_Stride, T_Entries> &b_T, Matrix *target) const
        {
            // a: (H x W)
            // b: (p x W)
            // R: (H x p)
            target->m_matrix = m_matrix * b_T.m_matrix.transpose();
        }

        void transposeMultiplyVector(Matrix &b, Matrix *target) const
        {
            // b : (H x 1)
            const int b_w = b.getWidth();
            const int b_h = b.getHeight();

            assert(b_w == 1);
            assert(getHeight() == b_h);

            // (self: H x W)T * (b: H x 1)     : W x 1
            target->m_matrix = m_matrix.transpose() * b.m_matrix;
        }

        void rightScale(Matrix &scale, SparseMatrix<T_Stride> *target)
        {
            auto vector = scale.m_matrix.col(0).transpose().array();
            target->m_matrix = (m_matrix.array().rowwise() * vector).matrix();
        }
        void leftScale(Matrix &scale, SparseMatrix<T_Stride> *target)
        {
            auto vector = scale.m_matrix.col(0).array();
            target->m_matrix = (m_matrix.array().colwise() * vector).matrix();
        }

    private:
        Eigen::MatrixXd m_matrix;
    };

}

#endif
