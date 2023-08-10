#ifndef ATG_SIMPLE_2D_CONSTRAINT_SOLVER_SPARSE_EIGEN_MATRIX_H
#define ATG_SIMPLE_2D_CONSTRAINT_SOLVER_SPARSE_EIGEN_MATRIX_H

#include <assert.h>
#include <stdint.h>
#include <string.h>

#ifdef ATG_S2C_USE_EIGEN

#define EIGEN_NO_DEBUG
#define EIGEN_NO_STATIC_ASSERT

#include "Eigen/Dense"
#include "Eigen/SparseCore"

namespace atg_scs {
    class Matrix;

    template <int T_Stride = 3, int T_Entries = 2>
    class SparseMatrix {
        // Use dense type for now
        typedef Eigen::MatrixXf MatrixType;
        // typedef Eigen::SparseMatrix<double, Eigen::RowMajor> MatrixType;

        public:
            SparseMatrix() {
            }

            ~SparseMatrix() {
            }

            void initialize(int width, int height) {
                m_matrix = MatrixType::Zero(height, width);
            }

            void resize(int width, int height) {
                m_matrix.conservativeResize(height, width);
            }

            void destroy() {}

            void expand(Matrix *matrix) {
                matrix->m_matrix = m_matrix;
            }

            void expandTransposed(Matrix *matrix) {
                matrix->m_matrix = m_matrix.transpose();
            }

            inline void setBlock(int row, int col, double v1, double v2, double v3) {
                m_matrix(row, col + 0) = v1;
                m_matrix(row, col + 1) = v2;
                m_matrix(row, col + 2) = v3;
            }

            // inline void setBlock(int row, int entry, uint8_t index) {
            //     assert(row >= 0 && row < m_height);
            //     assert(entry >= 0 && entry < T_Entries);
            //     assert(index < m_width);

            //     m_blockData[row * T_Entries + entry] = index;
            // }

            // inline void set(int row, int entry, int slice, double v) {
            //     assert(row >= 0 && row < m_height);
            //     assert(entry >= 0 && entry < T_Entries);
            //     assert(slice < T_Stride);

            //     m_matrix(row, entry * T_Stride + slice) = v;
            // }

            inline void set(int row, int col, double value) {
                m_matrix(row, col) = value;
            }

            // inline double get(int row, int entry, int slice) {
            //     assert(row >= 0 && row < m_height);
            //     assert(entry >= 0 && entry < T_Entries);
            //     assert(slice < T_Stride);

            //     return m_matrix(row, entry * T_Stride + slice);
            // }

            inline double get(int row, int col) {
                return m_matrix(row, col);
            }

            inline void setEmpty(int row, int col) {
                m_matrix = MatrixType::Zero(row, col);
            }

            void multiplyTranspose(const SparseMatrix<T_Stride, T_Entries> &b_T, Matrix *target) const {
                //   self: h x w
                //    b_T: n x w
                // target: h x n
                target->m_matrix = m_matrix * b_T.m_matrix.transpose();

                assert(target->getWidth() == b_T.getHeight());
                assert(target->getHeight() == getHeight());
            }

            void transposeMultiplyVector(Matrix &b, Matrix *target) const {
                //   self: h x w
                //      b: h x 1
                // target: w x 1
                target->m_matrix = m_matrix.transpose() * b.m_matrix;

                assert(target->getWidth() == 1);
                assert(target->getHeight() == getWidth());
            }

            void multiply(Matrix &b, Matrix *target) const {
                //   self: h x w
                //      b: w x n
                // target: h x n

                target->m_matrix = m_matrix * b.m_matrix;

                assert(target->getHeight() == getHeight());
                assert(target->getWidth() == b.getWidth());
            }

            void rightScale(Matrix &scale, SparseMatrix<T_Stride> *target) {
                // Rowwise multiply
                auto vector = scale.m_matrix.col(0).transpose().array();
                target->m_matrix = (m_matrix.array().rowwise() * vector).matrix();

                assert(target->m_matrix.getHeight() == getHeight());
                assert(target->m_matrix.getWidth() == getWidth());
            }

            void leftScale(Matrix &scale, SparseMatrix<T_Stride> *target) {
                // Column wise multiply
                auto vector = scale.m_matrix.col(0).array();
                target->m_matrix = (m_matrix.array().colwise() * vector).matrix();
            
                assert(target->m_matrix.getHeight() == getHeight());
                assert(target->m_matrix.getWidth() == getWidth());
            }

            __forceinline int getWidth() const { return m_matrix.cols(); }
            __forceinline int getHeight() const { return m_matrix.rows(); }

        public:
            MatrixType m_matrix;
    };
} /* namespace atg_scs */

#endif

#endif /* ATG_SIMPLE_2D_CONSTRAINT_SOLVER_SPARSE_MATRIX_H */
