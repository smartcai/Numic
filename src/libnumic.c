#include "libnumic.h"

#define PTR_FREE(p) (free(p), p = NULL)

/* Define a matrix as math */
typedef struct Matrix {
	int rows;
	int cols;
	double *data;
} Mat;

static Mat* Matrix_New(int rows, int cols);
static void Matrix_Del(Mat *m);

MObject *Matrix_Object_New(int rows, int cols)
{
	MObject *mobj = (MObject*)malloc(sizeof(MObject));
	mobj->m = Matrix_New(rows, cols);
	return mobj;
}
void Matrix_Object_Del(MObject *mobj)
{
	if (mobj != NULL) {
		Matrix_Del(mobj->m);
		PTR_FREE(mobj);
	}
}

static Mat* Matrix_New(int rows, int cols)
{
	Mat* m = (Mat*)malloc(sizeof(Mat));
	m->rows = rows;
	m->cols = cols;
	m->data = (double*)malloc(sizeof(double) * rows * cols);
	return m;
}

static void Matrix_Del(Mat *m)
{
	if (m != NULL) {
		if (m->data != NULL) {
			PTR_FREE(m->data);
		}
		PTR_FREE(m);
	}
}
