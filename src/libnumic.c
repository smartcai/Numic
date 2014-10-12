#include "libnumic.h"

#define FREE_SET_NULL(p) (free(p), p = NULL)

/**
 * @matrix: defined as math.
 */
struct matrix {
	int cols;
	int rows;
	double *array;
};

matrix *create_matrix(int rows, int cols)
{
	matrix *mp = (matrix*) malloc (sizeof(matrix));
	mp->rows = rows;
	mp->cols = cols;
	mp->array = (double*) calloc (cols * rows, sizeof(double));
	if (mp->array == NULL)
	{
		printf("%s: %d, out of memory.\n", __FILE__, __LINE__);
		exit(-1);
	}
	return mp;
}

void destroy_matrix(matrix *mp)
{
	if (mp != NULL)
		FREE_SET_NULL(mp->array);
	FREE_SET_NULL(mp);
}

double get_element(matrix *mp, int i, int j)
{
	double *p = mp->array;
	int step = mp->rows;
	int offset = i + j * step;
	return *(p+offset);
}

void set_element(matrix *mp, int i, int j, double val)
{
	double *p = mp->array;
	int step = mp->rows;
	int offset = i + j * step;
	*(p+offset) = val;
}
