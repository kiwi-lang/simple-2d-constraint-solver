#ifndef ATG_SIMPLE_2D_CONSTRAINT_SOLVER_SPARSE_EIGEN_MATRIX_H
#define ATG_SIMPLE_2D_CONSTRAINT_SOLVER_SPARSE_EIGEN_MATRIX_H

#include "matrix_eigen.h"
#include "Eigen/Core"

namespace atg_scs
{
    /* Fake sparse matrix: this is actually much slower than the sparce matrix
     *
     * A.llt().solve(b)
     */
    template <int T_Stride = 3, int T_Entries = 2>
    class SparseMatrix
    {
    public:
        typedef Eigen::Matrix<uint8_t, -1, -1> BlockMatrix;

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
            m_block = BlockMatrix::Zero(height, T_Entries);
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
            int previous = m_block(row, entry);

            if (previous != 0)
            {
                auto old_block = m_matrix.block(row, previous * T_Stride, 0, T_Stride);
                m_matrix.block(row, index * T_Stride, 0, T_Stride) = old_block;
            }

            m_block(row, entry) = index;
        }

        inline void set(int row, int entry, int slice, double v)
        {
            assert(row >= 0 && row < getHeight());
            assert(entry >= 0 && entry < T_Entries);
            assert(slice < T_Stride);

            int offset = m_block(row, entry);
            m_matrix(row, offset * T_Stride + slice) = v;
        }

        inline double get(int row, int entry, int slice)
        {
            assert(row >= 0 && row < getHeight());
            assert(entry >= 0 && entry < T_Entries);
            assert(slice < T_Stride);

            int offset = m_block(row, entry);
            return m_matrix(row, offset * T_Stride + slice);
        }

        inline void setEmpty(int row, int col)
        {
            assert(row >= 0 && row < getHeight());
            assert(col >= 0 && col < T_Entries);

            // Removes a block
            int offset = m_block(row, entry);
            m_matrix.block(row, offset * T_Stride, 0, T_Stride) = 0;
        }

        void multiply(Matrix &b, Matrix *target) const
        {
            target->initialize(b.getWidth(), getHeight());
            target->m_matrix = m_matrix * b.m_matrix;
        }

        void multiplyTranspose(const SparseMatrix<T_Stride, T_Entries> &b_T, Matrix *target) const
        {
            target->initialize(b_T.getHeight(), getHeight());
            // a: (H x W)
            // b: (p x W)
            // R: (H x p)
            target->m_matrix = m_matrix * b_T.m_matrix.transpose();
        }

        void transposeMultiplyVector(Matrix &b, Matrix *target) const
        {
            target->initialize(1, getWidth());

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
            target->initialize(getWidth(), getHeight());

            auto vector = scale.m_matrix.col(0).transpose().array();
            target->m_matrix = (m_matrix.array().rowwise() * vector).matrix();
        }
        void leftScale(Matrix &scale, SparseMatrix<T_Stride> *target)
        {
            target->initialize(getWidth(), getHeight());

            auto vector = scale.m_matrix.col(0).array();
            target->m_matrix = (m_matrix.array().colwise() * vector).matrix();
        }

    private:
        BlockMatrix m_block;
        Eigen::MatrixXd m_matrix;
    };

}

#endif
