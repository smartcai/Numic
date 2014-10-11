#ifndef __NUM_LIBNUMIC_H__
#define __NUM_LIBNUMIC_H__

#include <stddef.h>
#include <stdlib.h>

typedef struct Matrix_Object MObject;

/**
 * @Matrix_Object: a wrapper of matrix.
 * - Put all the interfaces here.
 *
 * @matrix: an incomplete type.
 *
 * - Fully declared in "libnumic.c". So the matrix data is not accessible
 * directly in user space, just like "private class member" in C++. The matrix
 * store format is critical to performance. So design it seperately.
 */
struct Matrix_Object {
	struct matrix *m;
};

MObject *matrix_object_new(int rows, int cols);
void matrix_object_del(MObject *mobj);

#endif  /* __NUM_LIBNUMIC_H__ */
