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

#define isSameBufSize(m1, m2)                       \
	(m1->rows * m1->cols == m2->rows * m2->cols)

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
	ASSERT(isSameBufSize(src, dst), ERR_MISMATCH_SIZE);

	dst->rows = src->rows;
	dst->cols = src->cols;

	memcpy(dst->array, src->array, src->rows * src->cols * sizeof(scalar));
}

void get_block(matrix *mat, int i, int j, matrix *blk)
{
	/* Get a block from matrix, starting from point (i, j). */

	int k;
	scalar *p1, *p2;

	ASSERT( (i + blk->rows <= mat->rows && j + blk->cols <= mat->cols)  \
	       , ERR_INVALID_DIMENSION);

	p1 = mat->array + (i + j * mat->rows); /* Point (i, j) */
	p2 = blk->array;

	/* Copy by columns. */
	for (k = 0; k < blk->cols; k++) {
		memcpy(p2, p1, blk->rows * sizeof(scalar));
		p1 += mat->rows;
		p2 += blk->rows;
	}
}

void set_block(matrix *mat, int i, int j, matrix *blk)
{
	/* Set a block to matrix, starting from point (i, j). */

	int k;
	scalar *p1, *p2;

	ASSERT( (i + blk->rows <= mat->rows && j + blk->cols <= mat->cols)  \
	       , ERR_INVALID_DIMENSION);

	p1 = mat->array + (i + j * mat->rows); /* Point (i, j) */
	p2 = blk->array;

	/* Copy by columns. */
	for (k = 0; k < blk->cols; k++) {
		memcpy(p1, p2, blk->rows * sizeof(scalar));
		p1 += mat->rows;
		p2 += blk->rows;
	}
}

inline scalar get_element(matrix *mp, int i, int j)
{
	/* Column-Major Order */
	return mp->array[i + j * mp->rows];
}

inline void set_element(matrix *mp, int i, int j, scalar val)
{
	/* Column-Major Order */
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

	ASSERT(isSameBufSize(src, dst), ERR_MISMATCH_SIZE);

	matrix *tmp = create_matrix(src->cols, src->rows);

	for (j = 0; j < src->cols; j++) {
		for (i = 0; i < src->rows; i++) {
			set_element(tmp, j, i, get_element(src, i, j));
		}
	}

	copy_matrix(tmp, dst);
	destroy_matrix(tmp);
}

void zero_matrix(matrix *mp)
{
	memset(mp->array, 0, mp->cols * mp->rows * sizeof(scalar));
}

void scalar_matrix_mul(matrix *dst, scalar alpha, matrix *src)
{
	/* dst = alpha * src */

	ASSERT(isSameSize(dst, src), ERR_MISMATCH_SIZE);

	int k;
	int m = src->rows;
	int n = src->cols;

	for (k = 0; k < m * n; k++) {
		dst->array[k] = alpha * src->array[k];
	}
}

/* void matrix_mul(matrix *dst, matrix *src1, matrix *src2)
 * {
 * 	ASSERT((src1->cols == src2->rows), ERR_MISMATCH_SIZE);
 * 	ASSERT((dst->rows == src1->rows) && (dst->cols == src2->cols),  \
 * 	       ERR_MISMATCH_SIZE);
 * } */

/* void qr_decompose_cgs(matrix *src, matrix *Q, matrix *R)
 * {
 * 	/\* Classical Gram-Schmidt Algorithm. *\/
 * 
 * 	vector *vp;
 * 
 * 	ASSERT(isSameSize(src, Q), ERR_MISMATCH_SIZE);
 * 	ASSERT(isSquareMatrix(R), ERR_MISMATCH_SIZE);
 * 	ASSERT((Q->cols == R->rows), ERR_MISMATCH_SIZE);
 * 
 * 	vp = create_col_vector(get_rows(src));
 * 
 * 	/\* Initialization. *\/
 * 	copy_matrix(src, Q);
 * 	zero_matrix(R);
 * 
 * 	get_col_vector(src, 0, vp);
 * 	set_element(R, 0, 0, vector_norm(vp));
 * 
 * 	free(vp);
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

inline void scalar_vector_mul(vector *dst, scalar alpha, vector *src)
{
	ASSERT(((isVector(src) && isVector(dst))), ERR_INVALID_DIMENSION);
	scalar_matrix_mul(dst, alpha, src);
}

void saxpy(vector *y, scalar a, vector *x)
{
	/* Saxpy: y = y + ax */
	int k, n = get_dim(y);

	ASSERT(isSameSize(y, x), ERR_MISMATCH_SIZE);

	for (k = 0; k < n; k++) {
		y->array[k] = y->array[k] + a * x->array[k];
	}
}

void gaxpy(vector *y, matrix *A, vector *x)
{
	/* Gaxpy: y = y + Ax */
	int m, n, j;
	vector *vp;

	ASSERT((get_rows(A) == get_dim(y)), ERR_MISMATCH_SIZE);
	ASSERT((get_cols(A) == get_dim(x)), ERR_MISMATCH_SIZE);

	m = get_rows(A);
	n = get_cols(A);
	vp = create_col_vector(m);

	for (j = 0; j < n; j++) {
		get_col_vector(A, j, vp);
		saxpy(y, x->array[j], vp);
	}

	destroy_vector(vp);
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

scalar vector_norm(vector *vp)
{
	return sqrt(dot_product(vp, vp));
}

void get_col_vector(matrix *mp, int k, vector *vp)
{
	/* Get the k-th column vector. */

	ASSERT(isColVector(vp), ERR_INVALID_DIMENSION);
	ASSERT((get_rows(mp) == get_dim(vp)), ERR_MISMATCH_SIZE);
	ASSERT((k < get_cols(mp)), ERR_INVALID_DIMENSION);

	get_block(mp, 0, k, vp);
}

void set_col_vector(matrix *mp, int k, vector *vp)
{
	/* Get the k-th column vector. */

	ASSERT(isColVector(vp), ERR_INVALID_DIMENSION);
	ASSERT((get_rows(mp) == get_dim(vp)), ERR_MISMATCH_SIZE);
	ASSERT((k < get_cols(mp)), ERR_INVALID_DIMENSION);

	set_block(mp, 0, k, vp);
}
