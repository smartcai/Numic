#include "libnumic.h"

#define FREE_SET_NULL(p) (free(p), p = NULL)

/**
 * @matrix: defined as math.
 */
struct matrix {
	int rows;
	int cols;
	double *data;
};

matrix *create_matrix(int rows, int cols)
{
	matrix *mp = (matrix*) malloc (sizeof(matrix));
	mp->rows = rows;
	mp->cols = cols;
	mp->data = (double*) malloc (sizeof(double) * rows * cols);
	if (mp->data == NULL)
	{
		printf("Allocating matrix data failed.\n");
		exit(-1);
	}
	return mp;
}

void destroy_matrix(matrix *mp)
{
	if (mp != NULL)
		FREE_SET_NULL(mp->data);
	FREE_SET_NULL(mp);
}

double get_element(matrix *mp, int i, int j)
{
	double *p = mp->data;
	int step = mp->rows;
	int offset = i * step + j;
	return *(p+offset);
}

void set_element(matrix *mp, int i, int j, double val)
{
	double *p = mp->data;
	int step = mp->rows;
	int offset = i * step + j;
	*(p+offset) = val;
}
