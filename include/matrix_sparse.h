#ifndef ATG_SIMPLE_2D_CONSTRAINT_SOLVER_MATRIX_SPARSE_H
#define ATG_SIMPLE_2D_CONSTRAINT_SOLVER_MATRIX_SPARSE_H

#define ATG_S2C_USE_EIGEN_SPARSE

#ifdef ATG_S2C_USE_EIGEN_SPARSE
#include "matrix_sparse_eigen.h"
#else
#include "matrix_sparse_custom.h"
#endif



#endif /* ATG_SIMPLE_2D_CONSTRAINT_SOLVER_MATRIX_H */
