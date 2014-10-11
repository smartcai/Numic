#include "libnumic.h"

#define PTR_FREE(p) (free(p), p = NULL)

/**
 * @matrix: defined as math.
 */
typedef struct matrix {
	int rows;
	int cols;
	double *data;
} mat;

static mat* mat_new(int rows, int cols);
static void mat_del(mat *m);

MObject *matrix_object_new(int rows, int cols)
{
	MObject *mobj = (MObject*)malloc(sizeof(MObject));
	mobj->m = mat_new(rows, cols);
	return mobj;
}
void matrix_object_del(MObject *mobj)
{
	if (mobj != NULL) {
		mat_del(mobj->m);
		PTR_FREE(mobj);
	}
}

static mat* mat_new(int rows, int cols)
{
	mat* m = (mat*)malloc(sizeof(mat));
	m->rows = rows;
	m->cols = cols;
	m->data = (double*)malloc(sizeof(double) * rows * cols);
	return m;
}

static void mat_del(mat *m)
{
	if (m != NULL) {
		if (m->data != NULL) {
			PTR_FREE(m->data);
		}
		PTR_FREE(m);
	}
}
