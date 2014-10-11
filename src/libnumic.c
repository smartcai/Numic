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
static void mat_set(mat *m, int i, int j, double e);
static double mat_get(mat *m, int i, int j);

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

void matrix_object_set(MObject *mobj, int i, int j, double elt)
{
	mat_set(mobj->m, i, j, elt);
}
double matrix_object_get(MObject *mobj, int i, int j)
{
	return mat_get(mobj->m, i, j);
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

static void mat_set(mat *m, int i, int j, double e)
{
	int step = m->cols;
	int addr = i*step + j;
	m->data[addr] = e;
}

static double mat_get(mat *m, int i , int j)
{
	double e;
	int step = m->cols;
	int addr = i*step + j;
	e = m->data[addr];
	return e;
}
