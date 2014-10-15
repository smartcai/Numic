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

#define ERR_OUT_OF_MEMORY     ", out of memory."
#define ERR_MISMATCH_SIZE     ", mismatching size."
#define ERR_INVALID_DIMENSION ", invalid dimension."

/**
 * @matrix: defined as math.
 */
struct matrix {
	int     cols;
	int     rows;
	scalar *array;
};

#define isSameSize(m1, m2)                              \
	((m1->rows == m2->rows) && (m1->cols == m2->cols))

#define isTransposedSize(m1, m2)                        \
	((m1->rows == m2->cols) && (m1->cols == m2->rows))

#define isSquareMatrix(m)                       \
	(m->rows == m->cols)

/**
 * @vector: a special kind of matrix.
 */
#define isRowVector(v) (v->rows == 1)
#define isColVector(v) (v->cols == 1)
#define isVector(v)    ((isRowVector(v)) ^ (isColVector(v)))

matrix *create_matrix(int rows, int cols)
{
	matrix *mp = (matrix*) malloc (sizeof(matrix));

	ASSERT((rows > 0 && cols > 0), ERR_INVALID_DIMENSION);

	mp->rows = rows;
	mp->cols = cols;
	mp->array = (scalar*) calloc (cols * rows, sizeof(scalar));

	ASSERT((mp->array != NULL), ERR_OUT_OF_MEMORY);

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
	       ERR_MISMATCH_SIZE);

	memcpy(dst->array, src->array, n * m * sizeof(scalar));
}

inline scalar get_element(matrix *mp, int i, int j)
{
	return mp->array[i + j * mp->rows];
}

inline void set_element(matrix *mp, int i, int j, scalar val)
{
	mp->array[i + j * mp->rows] = val;
}

inline int get_rows(matrix *mp)
{
	return mp->rows;
}

inline int get_cols(matrix *mp)
{
	return mp->cols;
}

void transpose(matrix *src, matrix* dst)
{
	int i, j;
	int n = src->cols;
	int m = src->rows;

	ASSERT((m == dst->cols && n == dst->rows),  \
	       ERR_MISMATCH_SIZE);

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

/* void qr_decompose_cgs(matrix *src, matrix *Q, matrix *R)
 * {
 * 	/\* Classical Gram-Schmidt Algorithm. *\/
 * 
 * 	ASSERT(isSameSize(src, Q), ERR_MISMATCH_SIZE);
 * 	ASSERT(isSquareMatrix(R), ERR_MISMATCH_SIZE);
 * 	ASSERT((Q->cols == R->rows), ERR_MISMATCH_SIZE);
 * } */

inline vector *create_col_vector(int dim)
{
	return create_matrix(dim, 1);
}

inline vector *create_row_vector(int dim)
{
	return create_matrix(1, dim);
}

inline void destroy_vector(vector *v)
{
	destroy_matrix(v);
}

inline void print_vector(vector *v)
{
	ASSERT((isVector(v)), ERR_INVALID_DIMENSION);
	print_matrix(v);
}

inline void copy_vector(vector *src, vector *dst)
{
	ASSERT((isVector(src) && isVector(dst)), ERR_INVALID_DIMENSION);
	copy_matrix(src, dst);
}

inline scalar get_vector_element(vector *v, int k)
{
	ASSERT((isVector(v)), ERR_INVALID_DIMENSION);
	return get_element(v, k, 0);
}

inline void set_vector_element(vector *v, int k, scalar val)
{
	ASSERT((isVector(v)), ERR_INVALID_DIMENSION);
	set_element(v, k, 0, val);
}

inline int get_dim(vector *v)
{
	return isColVector(v) ? v->rows : v->cols;
}

inline void transpose_vector(vector *src, vector* dst)
{
	ASSERT((isVector(src) && isVector(dst)), ERR_INVALID_DIMENSION);
	ASSERT((isTransposedSize(src, dst)), ERR_MISMATCH_SIZE);
	transpose(src, dst);
}

inline void zero_vector(vector *v)
{
	ASSERT((isVector(v)), ERR_INVALID_DIMENSION);
	zero_matrix(v);
}

scalar dot_product(vector *v1, vector *v2)
{
	double c;
	int n, k;

	ASSERT((isVector(v1) && isVector(v2)), ERR_INVALID_DIMENSION);
	ASSERT((isSameSize(v1, v2)), ERR_MISMATCH_SIZE);

	n = get_dim(v1);

	c = 0;
	for (k = 0; k < n; k++) {
		c += v1->array[k] * v2->array[k];
	}

	return c;
}
