#ifndef LIBMATRIX_H
#define LIBMATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


/* Types */
typedef double        Scalar;
typedef struct Vector Vector;
typedef struct Matrix Matrix;

typedef enum bTypeV {           /* Vector type: column vector is default. */
    COL,
    ROW
} bTypeV;


/* Types operations */
Vector *vector_new(int dim, bTypeV t);
void    vector_del(Vector *v);
void    vector_set(Vector *v, int k, Scalar val);
Scalar  vector_get(Vector *v, int k);
int     vector_dim(Vector *v);
bTypeV  vector_type(Vector *v);
void    vector_print(Vector *v);
void    vector_cpy(Vector *dst, Vector *src);

Matrix *matrix_new(int row, int col);
void    matrix_del(Matrix *m);
void    matrix_set(Matrix *m, int i, int j, Scalar val);
Scalar  matrix_get(Matrix *m, int i, int j);
int     matrix_row(Matrix *m);
int     matrix_col(Matrix *m);
void    matrix_print(Matrix *m);
void    matrix_cpy(Matrix *dst, Matrix *src);


/* Core operations */
void   blas_saxpy(Vector *y, Scalar  a, Vector *x); /* y = y + ax: Scalar A X Plus Y */
void   blas_gaxpy(Vector *y, Matrix *A, Vector *x); /* y = y + Ax: Generalized sAXPY */
void   blas_gemul(Matrix *C, Matrix *A, Matrix *B); /* C = C + aAB: GEneral matrix-matrix MULtiplication */

void   vector_transpose(Vector *y, Vector *x);
Scalar vector_dot(Vector *x, Vector *y);
void   vector_outer(Matrix *A, Vector *x, Vector *y);
void   vector_scalar_mul(Vector *y, Scalar a, Vector *x);
void   vector_scalar_add(Vector *y, Scalar a, Vector *x);
void   vector_add(Vector *z, Vector *x, Vector *y);
void   vector_sub(Vector *z, Vector *x, Vector *y);
void   vector_point_mul(Vector *z, Vector *x, Vector *y);
void   vector_point_div(Vector *z, Vector *x, Vector *y);
int    vector_any(Vector *x);
int    vector_all(Vector *x);
void   vector_zero(Vector *x);

void   matrix_transpose(Matrix *C, Matrix *A);
void   matrix_scalar_mul(Matrix *C, Scalar a, Matrix *A);
void   matrix_scalar_add(Matrix *C, Scalar a, Matrix *A);
void   matrix_add(Matrix *C, Matrix *A, Matrix *B);
void   matrix_sub(Matrix *C, Matrix *A, Matrix *B);
void   matrix_mul(Matrix *C, Matrix *A, Matrix *B);
void   matrix_point_mul(Matrix *C, Matrix *A, Matrix *B);
void   matrix_point_div(Matrix *C, Matrix *A, Matrix *B);
void   matrix_get_col_vector(Vector *y, Matrix *A, int j);
void   matrix_get_row_vector(Vector *y, Matrix *A, int i);
void   matrix_set_col_vector(Matrix *A, int j, Vector *x);
void   matrix_set_row_vector(Matrix *A, int i, Vector *x);
void   matrix_col_any(Vector *y, Matrix *A);
void   matrix_col_all(Vector *y, Matrix *A);
void   matrix_row_any(Vector *y, Matrix *A);
void   matrix_row_all(Vector *y, Matrix *A);
int    matrix_any(Matrix *A);
int    matrix_all(Matrix *A);
void   matrix_zero(Matrix *A);
void   matrix_eye(Matrix *A);


/* Self operations */


/* Linear algebra */
Scalar vector_L2_norm(Vector *v);
void   vector_housh(Vector *v, Scalar *beta, Vector *x, int j);
void   vector_givens(Scalar *c, Scalar *s, Scalar y, Scalar z);

void   matrix_full_qr_housh(Matrix *A, Matrix *Q, Matrix *R);
void   matrix_full_qr_givens(Matrix *A, Matrix *Q, Matrix *R);
void   matrix_full_qr_cgs(Matrix *A, Matrix *Q, Matrix *R);
void   matrix_full_qr_mgs(Matrix *A, Matrix *Q, Matrix *R);
void   matrix_bidiagonal(Matrix *A, Matrix *P, Matrix *B, Matrix * Q);

#endif /* LIBMATRIX_H */
