#ifndef __NUM_LIBNUMIC_H__
#define __NUM_LIBNUMIC_H__

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/* Incomplete type. */
typedef double scalar;
typedef struct matrix matrix;

/* Cache tiling size, defined by user. */
#define BLOCK_ROWS 64
#define BLOCK_COLS 48

/* Interface. */
matrix *create_matrix(int rows, int cols);
void destroy_matrix(matrix *mp);
void print_matrix(matrix *mp);
void copy_matrix(matrix *src, matrix *dst);
inline scalar get_element(matrix *mp, int i, int j);
inline void set_element(matrix *mp, int i, int j, scalar val);

/* Basic operations. */
void transpose(matrix *src, matrix* dst);
void zero_matrix(matrix *mp);

/**
 * Vector's part: definition and operations. We treat vector as a kind of
 * special matrix. Then all the methods are "syntactic sugars".
 */
/* Define vector as matrix. */
typedef matrix vector;

/* Map the matrix methods to vector. */
inline vector *create_col_vector(int dim);
inline vector *create_row_vector(int dim);
inline void destroy_vector(vector *v);
inline void print_vector(vector *v);
inline void copy_vector(vector *src, vector *dst);
inline scalar get_vector_element(vector *v, int k);
inline void set_vector_element(vector *v, int k, scalar val);

inline void transpose_vector(vector *src, vector* dst);
inline void zero_vector(vector *v);

#endif  /* __NUM_LIBNUMIC_H__ */
