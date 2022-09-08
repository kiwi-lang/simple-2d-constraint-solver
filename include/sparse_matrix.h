#ifndef ATG_SIMPLE_2D_CONSTRAINT_SOLVER_SPARSE_MATRIX_H
#define ATG_SIMPLE_2D_CONSTRAINT_SOLVER_SPARSE_MATRIX_H

#include <assert.h>
#include <stdint.h>
#include <string.h>

#define USE_TRUE_SPARSE

#ifndef USE_TRUE_SPARSE
#include "sparce_matrix_eigen.h"
#else

#include "Eigen/Dense"

namespace atg_scs {
    class Matrix;

    // Stride is the size of a block
    /*
     */
    template <int T_Stride = 3, int T_Entries = 2>
    class SparseMatrix {
        typedef Eigen::Matrix<uint8_t, Eigen::Dynamic, T_Entries> BlockMatrix;
        typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> DataMatrix;

        public:
            SparseMatrix() {
                m_width = m_height = 0;
                m_capacityHeight = 0;
            }

            ~SparseMatrix() {
                // assert(m_matrix == nullptr);
                // assert(m_data == nullptr);
                // assert(m_blockData == nullptr);
            }

            void initialize(int width, int height) {
                resize(width, height);
                // memset(m_blockData, 0xFFFFFF, sizeof(uint8_t) * T_Entries * m_height);
            }

            void resize(int width, int height) {
                if (width == m_width && height == m_height) return;
                else if (height > m_capacityHeight) {
                    destroy();

                    m_capacityHeight = (height > m_capacityHeight)
                        ? height
                        : m_capacityHeight;

                    // size_t data_size = T_Stride * T_Entries * m_capacityHeight;
                    // size_t block_size = m_capacityHeight * T_Entries;

                    // m_data = new double[data_size];
                    // m_matrix = new double *[m_capacityHeight];
                    m_data = DataMatrix::Constant(m_capacityHeight, T_Stride * T_Entries);
                    m_blockData = BlockMatrix::Constant(m_capacityHeight, 0xFFFFFF);
                    // m_blockData = new uint8_t[block_size];
                }

                m_height = height;
                m_width = width;

                // for (int i = 0; i < height; ++i) {
                //     m_matrix[i] = &m_data[i * T_Entries * T_Stride];
                // }
            }

            void destroy() {
                // if (m_matrix == nullptr) {
                //     return;
                // }

                // delete[] m_matrix;
                // delete[] m_data;
                // delete[] m_blockData;

                // m_matrix = nullptr;
                // m_data = nullptr;
                // m_blockData = nullptr;

                m_width = m_height = 0;
            }

            void expand(Matrix *matrix) {
                matrix->initialize(m_width, m_height);

                for (int i = 0; i < m_height; ++i) {
                    for (int j = 0; j < T_Entries; ++j) {
                        // const uint8_t block = m_blockData[i * T_Entries + j];
                        const uint8_t block = m_blockData(i, j);

                        if (block == 0xFF) continue;
                        else {
                            for (int k = 0; k < T_Stride; ++k) {
                                double value = m_data(i, j * T_Stride + k);
                                matrix->set(block * T_Stride + k, value);
                            }
                        }
                    }
                }
            }

            void expandTransposed(Matrix *matrix) {
                matrix->initialize(m_height, m_width);

                for (int i = 0; i < m_height; ++i) {
                    for (int j = 0; j < T_Entries; ++j) {
                        // const uint8_t block = m_blockData[i * T_Entries + j];
                        const uint8_t block = m_blockData(i, j);
                        if (block == 0xFF) continue;
                        else {
                            for (int k = 0; k < T_Stride; ++k) {
                                matrix->set(i, block * T_Stride + k, m_data(i, j * T_Stride + k));
                            }
                        }
                    }
                }
            }

            inline void setBlock(int row, int entry, uint8_t index) {
                assert(row >= 0 && row < m_height);
                assert(entry >= 0 && entry < T_Entries);
                assert(index * (T_Stride - 1) < m_width);

                m_blockData(row, entry) = index;
            }

            inline void set(int row, int entry, int slice, double v) {
                assert(row >= 0 && row < m_height);
                assert(entry >= 0 && entry < T_Entries);
                assert(slice < T_Stride);

                m_data(row, entry * T_Stride + slice) = v;
            }

            inline double get(int row, int entry, int slice) {
                assert(row >= 0 && row < m_height);
                assert(entry >= 0 && entry < T_Entries);
                assert(slice < T_Stride);

                return m_data(row, entry * T_Stride + slice);
            }

            inline void setEmpty(int row, int col) {
                assert(row >= 0 && row < m_height);
                assert(col >= 0 && col < T_Entries);

                m_blockData(row, col) = 0xFF;
                for (int i = 0; i < T_Stride; ++i) {
                    m_data(row, col * T_Stride + i) = 0;
                }
            }

            void multiplyTranspose(const SparseMatrix<T_Stride, T_Entries> &b_T, Matrix *target) const {
                assert(m_width == b_T.m_width);

                target->initialize(b_T.m_height, m_height);

                for (int i = 0; i < m_height; ++i) {
                    for (int j = 0; j < b_T.m_height; ++j) {
                        double dot = 0;
                        for (int k = 0; k < T_Entries; ++k) {
                            const uint8_t block0 = m_blockData(i, k);
                            if (block0 == 0xFF) continue;

                            for (int l = 0; l < T_Entries; ++l) {
                                const uint8_t block1 = b_T.m_blockData(j, l);
                                if (block0 == block1) {
                                    for (int m = 0; m < T_Stride; ++m) {
                                        dot +=
                                            m_data(i, k * T_Stride + m)
                                            * b_T.m_data(j, l * T_Stride + m);
                                    }
                                }
                            }
                        }

                        target->set(j, i, dot);
                    }
                }
            }

            void transposeMultiplyVector(Matrix &b, Matrix *target) const {
                const int b_w = b.getWidth();
                const int b_h = b.getHeight();

                assert(b_w == 1);
                assert(m_height == b_h);

                target->initialize(1, m_width);

                for (int i = 0; i < m_height; ++i) {
                    double v = 0.0;
                    for (int k = 0; k < T_Entries; ++k) {
                        const int offset = k * T_Stride;
                        const uint8_t block = m_blockData(i, k);
                        if (block == 0xFF) continue;

                        for (int l = 0; l < T_Stride; ++l) {
                            const int j = block * T_Stride + l;
                            target->add(0, j, m_data(i, offset + l) * b.get(0, i));
                        }
                    }
                }
            }

            void multiply(Matrix &b, Matrix *target) const {
                const int b_w = b.getWidth();
                const int b_h = b.getHeight();

                assert(m_width == b_h);

                target->initialize(b.getWidth(), m_height);

                for (int i = 0; i < m_height; ++i) {
                    for (int j = 0; j < b_w; ++j) {
                        double v = 0.0;
                        for (int k = 0; k < T_Entries; ++k) {
                            const int offset = k * T_Stride;

                            const uint8_t block = m_blockData(i, j);
                            // const uint8_t block = m_blockData[i * T_Entries + k];
                            if (block == 0xFF) continue;

                            for (int l = 0; l < T_Stride; ++l) {
                                v += m_data(i, offset + l) * b.get(j, block * T_Stride + l);
                            }
                        }

                        target->set(j, i, v);
                    }
                }
            }

            void rightScale(Matrix &scale, SparseMatrix<T_Stride> *target) {
                assert(scale.getWidth() == 1);
                assert(scale.getHeight() == m_width);

                target->initialize(m_width, m_height);

                for (int i = 0; i < m_height; ++i) {
                    for (int j = 0; j < T_Entries; ++j) {
                        const uint8_t index = m_blockData(i, j);
                        //const uint8_t index = m_blockData[i * T_Entries + j];
                        if (index == 0xFF) continue;

                        target->setBlock(i, j, index);

                        for (int k = 0; k < T_Stride; ++k) {
                            target->set(
                                i,
                                j,
                                k,
                                scale.get(0, index * T_Stride + k) * m_data(i, j * T_Stride + k));
                        }
                    }
                }
            }

            void leftScale(Matrix &scale, SparseMatrix<T_Stride> *target) {
                assert(scale.getWidth() == 1 || m_height == 0);
                assert(scale.getHeight() == m_height);

                target->initialize(m_width, m_height);

                for (int i = 0; i < m_height; ++i) {
                    for (int j = 0; j < T_Entries; ++j) {
                        const uint8_t index = m_blockData(i, j);
                        //const uint8_t index = m_blockData[i * T_Entries + j];

                        if (index == 0xFF) continue;

                        target->setBlock(i, j, index);

                        for (int k = 0; k < T_Stride; ++k) {
                            target->set(
                                i,
                                j,
                                k,
                                scale.get(0, i) *m_data(i, j * T_Stride + k));
                        }
                    }
                }
            }

            __forceinline int getWidth() const { return m_width; }
            __forceinline int getHeight() const { return m_height; }

        protected:
            // uint8_t *m_blockData;
            BlockMatrix m_blockData;
            DataMatrix m_data;

            int m_width;
            int m_height;
            int m_capacityHeight;
    };
} /* namespace atg_scs */
#endif

#endif /* ATG_SIMPLE_2D_CONSTRAINT_SOLVER_SPARSE_MATRIX_H */
