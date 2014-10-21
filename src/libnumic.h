#ifndef __NUM_LIBNUMIC_H__
#define __NUM_LIBNUMIC_H__

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/* Incomplete types. */
typedef double scalar;
typedef struct matrix matrix;
typedef matrix vector;          /* As a 1D matrix */

/* Cache tiling block, defined by user. */
#define BLOCK_ROWS 64
#define BLOCK_COLS 48

/* Interface. */
matrix *create_matrix(int rows, int cols);
void destroy_matrix(matrix *mp);
void print_matrix(matrix *mp);
void copy_matrix(matrix *src, matrix *dst);

/* Matrix operations. */
inline scalar get_element(matrix *mp, int i, int j);
inline void set_element(matrix *mp, int i, int j, scalar val);
void get_block(matrix *mat, int i, int j, matrix *blk);
void set_block(matrix *mat, int i, int j, matrix *blk);
inline int get_rows(matrix *mp);
inline int get_cols(matrix *mp);

void transpose(matrix *src, matrix* dst);
void zero_matrix(matrix *mp);
void scalar_matrix_mul(matrix *dst, scalar alpha, matrix *src);
void matrix_mul(matrix *dst, matrix *src1, matrix *src2);
void subtract_matrix(matrix *dst, matrix *src);

void qr_decompose_cgs(matrix *A, matrix *Q, matrix *R);

/**
 * Map the matrix methods to vector.
 */
inline vector *create_col_vector(int dim);
inline vector *create_row_vector(int dim);
inline void destroy_vector(vector *v);
inline void print_vector(vector *v);
inline void copy_vector(vector *src, vector *dst);

/* Vector operations. */
inline scalar get_vector_element(vector *v, int k);
inline void set_vector_element(vector *v, int k, scalar val);
inline int get_dim(vector *v);

inline void transpose_vector(vector *src, vector* dst);
inline void zero_vector(vector *v);
inline void scalar_vector_mul(vector *dst, scalar alpha, vector *src);

void saxpy(vector *y, scalar a, vector *x);
void gaxpy(vector *y, matrix *A, vector *x);
scalar dot_product(vector *v1, vector *v2);
scalar vector_norm(vector *vp);
void out_product(matrix *mp, vector *v1, vector *v2);

void get_col_vector(matrix *mp, int k, vector *vp);
void set_col_vector(matrix *mp, int k, vector *vp);

void householder_vector(vector *x, vector *v, scalar *beta, int k);

#endif  /* __NUM_LIBNUMIC_H__ */
