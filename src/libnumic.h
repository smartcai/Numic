#ifndef __NUM_LIBNUMIC_H__
#define __NUM_LIBNUMIC_H__

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

/* Incomplete type. */
typedef struct matrix matrix;

/* Interface. */
matrix *create_matrix(int rows, int cols);
void destroy_matrix(matrix *mp);
double get_element(matrix *mp, int i, int j);
void set_element(matrix *mp, int i, int j, double val);

/* Basic operations. */
void transpose(matrix *src, matrix* dst);

#endif  /* __NUM_LIBNUMIC_H__ */
