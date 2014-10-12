#ifndef __NUM_LIBNUMIC_H__
#define __NUM_LIBNUMIC_H__

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>


/* Incomplete type. */
typedef struct matrix matrix;

/* Interface. */
matrix *create_matrix(int rows, int cols);
void destroy_matrix(matrix *mp);
void print_matrix(matrix *mp);
inline double get_element(matrix *mp, int i, int j);
inline void set_element(matrix *mp, int i, int j, double val);

/* Basic operations. */
void transpose(matrix *src, matrix* dst);
void zero_matrix(matrix *mp);

#endif  /* __NUM_LIBNUMIC_H__ */
