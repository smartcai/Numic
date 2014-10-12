#include "libnumic.h"

#define FREE_SET_NULL(p) (free(p), p = NULL)                            \

#define ERR_MSG(Str)                                                    \
	printf("ERROR: %s: L%d, %s " Str, __FILE__, __LINE__, __FUNCTION__)

#define ASSERT(stmt, S) do {                    \
		if (!stmt) {                            \
			ERR_MSG(S);                         \
			exit(-1);                           \
		}                                       \
	} while(0)

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

	ASSERT((mp->array != NULL), ", out of memory.");

	return mp;
}

void destroy_matrix(matrix *mp)
{
	if (mp != NULL)
		FREE_SET_NULL(mp->array);
	FREE_SET_NULL(mp);
}

void print_matrix(matrix *mp)
{
	int i, j;
	for (i = 0; i < mp->rows; i++) {
		for (j = 0; j < mp->cols; j++)
			printf("%f ", get_element(mp, i, j));
		printf("\n");
	}
}

inline double get_element(matrix *mp, int i, int j)
{
	return *(mp->array + (i + j * mp->rows));
}

inline void set_element(matrix *mp, int i, int j, double val)
{
	*(mp->array + (i + j * mp->rows)) = val;
}

void transpose(matrix *src, matrix* dst)
{
	int i, j;
	int n = src->cols;
	int m = src->rows;

	ASSERT(m == dst->cols && n == dst->rows, ", mismatching size.");

	for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
			set_element(dst, j, i, get_element(src, i, j));
		}
	}
}
