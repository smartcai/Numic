#include "libnumic.h"

#define FREE_SET_NULL(p) (free(p), p = NULL)                            \

#define ERR_MSG(Str)                                                    \
	printf("ERROR: %s: L%d, %s " Str "\n", __FILE__, __LINE__, __FUNCTION__)

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
	scalar *array;
};

matrix *create_matrix(int rows, int cols)
{
	matrix *mp = (matrix*) malloc (sizeof(matrix));

	ASSERT((rows > 0 && cols > 0), ", invalid dimensions.");

	mp->rows = rows;
	mp->cols = cols;
	mp->array = (scalar*) calloc (cols * rows, sizeof(scalar));

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

void copy_matrix(matrix *src, matrix *dst)
{
	int n = src->cols;
	int m = src->rows;

	ASSERT((n == dst->cols && m == dst->rows),  \
	       ", mismatching size.");

	memcpy(dst->array, src->array, n * m * sizeof(scalar));
}

inline scalar get_element(matrix *mp, int i, int j)
{
	return *(mp->array + (i + j * mp->rows));
}

inline void set_element(matrix *mp, int i, int j, scalar val)
{
	*(mp->array + (i + j * mp->rows)) = val;
}

void transpose(matrix *src, matrix* dst)
{
	int i, j;
	int n = src->cols;
	int m = src->rows;

	ASSERT((m == dst->cols && n == dst->rows),  \
	       ", mismatching size.");

	zero_matrix(dst);

	for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
			set_element(dst, j, i, get_element(src, i, j));
		}
	}
}

void zero_matrix(matrix *mp)
{
	memset(mp->array, 0, mp->cols * mp->rows * sizeof(scalar));
}
