#include "libnumic.h"

#define EPSILON (pow(2.0, -50))
#define ABS_F(x) (x < 0 ? -x : x)
#define isNotZero(x) (ABS_F(x) > EPSILON)
#define isZero(x) (ABS_F(x) <= EPSILON)

#define FREE_SET_NULL(p) (free(p), p = NULL)

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
#define ERR_DIVIDE_BY_ZERO    ", divide by zero."

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
#define isVector(v)    ((isRowVector(v)) || (isColVector(v)))

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

void matrix_mul(matrix *dst, matrix *src1, matrix *src2)
{
	/* dst = src1 * src2 */

	int m, r, n;
	int i, j, k;
	vector *y, *x;

	ASSERT((src1->cols == src2->rows), ERR_MISMATCH_SIZE);
	ASSERT((dst->rows == src1->rows) && (dst->cols == src2->cols),  \
	       ERR_MISMATCH_SIZE);

	m = get_rows(src1);
	r = get_cols(src1);
	n = get_cols(src2);

	zero_matrix(dst);
	y = create_col_vector(m);
	x = create_col_vector(r);

	for (j = 0; j < n; j++) {
		zero_vector(y);
		get_col_vector(src2, j, x);
		gaxpy(y, src1, x);
		set_col_vector(dst, j, y);
	}

	destroy_matrix(y);
	destroy_matrix(x);
}

void subtract_matrix(matrix *dst, matrix *src)
{
	int m, n, k;

	ASSERT(isSameSize(dst, src), ERR_MISMATCH_SIZE);

	m = get_rows(dst);
	n = get_cols(src);

	for (k = 0; k < m * n; k++) {
		dst->array[k] -= src->array[k];
	}
}

void add_matrix(matrix *dst, matrix *src)
{
	int m, n, k;

	ASSERT(isSameSize(dst, src), ERR_MISMATCH_SIZE);

	m = get_rows(dst);
	n = get_cols(src);

	for (k = 0; k < m * n; k++) {
		dst->array[k] += src->array[k];
	}
}

void qr_decompose_cgs(matrix *A, matrix *Q, matrix *R)
{
	/* CGS - Classical Gram-Schmidt Algorithm. */

	vector *Ak, *Qk;
	scalar L2, r;
	int m, n;
	int k, kk;

	ASSERT(isSameSize(A, Q), ERR_MISMATCH_SIZE);
	ASSERT(isSquareMatrix(R), ERR_MISMATCH_SIZE);
	ASSERT((Q->cols == R->rows), ERR_MISMATCH_SIZE);

	m = get_rows(A);
	n = get_cols(A);
	Ak = create_col_vector(m);
	Qk = create_col_vector(m);

	/* Initialization. */
	copy_matrix(A, Q);
	zero_matrix(R);

	get_col_vector(A, 0, Ak);
	L2 = vector_norm(Ak);
	set_element(R, 0, 0, L2);

	ASSERT(isNotZero(L2), ERR_DIVIDE_BY_ZERO);
	scalar_vector_mul(Ak, 1.0 / L2, Ak);
	set_col_vector(Q, 0, Ak);

	for (k = 1; k < n; k++) {
		get_col_vector(A, k, Ak);
		for (kk = 0; kk < k; kk++) {
			get_col_vector(Q, kk, Qk);
			r = dot_product(Ak, Qk);
			set_element(R, kk, k, r);
			r = -r;
			saxpy(Ak, r, Qk);
		}
		L2 = vector_norm(Ak);
		set_element(R, k, k, L2);
		ASSERT(isNotZero(L2), ERR_DIVIDE_BY_ZERO);
		scalar_vector_mul(Ak, 1.0 / L2, Ak);
		set_col_vector(Q, k, Ak);
	}

	destroy_matrix(Ak);
	destroy_matrix(Qk);
}

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

void out_product(matrix *mp, vector *v1, vector *v2)
{
	/* mp = v1 * transpose(v2) */
	int m = get_dim(v1);
	int n = get_dim(v2);
	int i, j;
	scalar product;

	ASSERT((m == get_rows(mp)), ERR_MISMATCH_SIZE);
	ASSERT((n == get_cols(mp)), ERR_MISMATCH_SIZE);

	for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
			product = get_vector_element(v1, i) * get_vector_element(v2, j);
			set_element(mp, i, j, product);
		}
	}
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

void householder_vector(vector *x, vector *v, scalar *beta, int k)
{
	/**
	 * Rotate x to the dimension of k.
	 * (I - beta*v*vt)x =   norm(x) * e(k) if x(k) <  0
	 * (I - beta*v*vt)x = - norm(x) * e(k) if x(k) >= 0
	 */
	int m = get_dim(x);
	scalar xk = get_vector_element(x, k), vk;
	scalar sigma = dot_product(x, x) - xk * xk;
	scalar mu;
	copy_vector(x, v);
	set_vector_element(v, k, 1);

	if (isZero(sigma) && xk >= 0)
		*beta = 0;
	else if (isZero(sigma) && xk < 0)
		*beta = -2;
	else {
		mu = sqrt(xk * xk + sigma);
		if (xk <= 0)
			vk = xk - mu;
		else
			vk = -sigma / (xk + mu);
		set_vector_element( v, k, vk );
		*beta = 2 * vk * vk / (sigma + vk * vk);
		scalar_vector_mul(v, 1.0 / vk, v);
	}
}

void house_matrix_columns(matrix *src, matrix *dst)
{
	/* Factored-form Householder QR
	 * Store the Householder vectors of src to the lower triangular dst
	 * Store the R to the diag and upper triangular. */

	int m, n;
	int i, j;
	matrix *tmp;

	ASSERT(isSameSize(src, dst), ERR_MISMATCH_SIZE);

	m = get_rows(src);
	n = get_cols(src);

	copy_matrix(src, dst);
	tmp = create_matrix(m, n);

	for (j = 0; j < n; j++) {
		vector *vj = create_col_vector(m - j);
		vector *vh = create_col_vector(m - j);
		vector *vt = create_row_vector(m - j);
		matrix *Qj = create_matrix(m - j, n - j);
		matrix *Qp = create_matrix(m - j, n - j);
		vector *vn = create_row_vector(n - j);
		scalar beta;

		get_block(dst, j, j, vj);
		get_block(dst, j, j, Qj);

		householder_vector(vj, vh, &beta, 0);
		set_block(tmp, j, j, vj);

		transpose_vector(vh, vt);
		matrix_mul(vn, vt, Qj);
		scalar_vector_mul(vh, beta, vh);

		out_product(Qp, vh, vn);
		subtract_matrix(Qj, Qp);

		set_block(dst, j, j, Qj);

		destroy_vector(vj);
		destroy_vector(vh);
		destroy_matrix(Qj);
		destroy_vector(vn);
		destroy_vector(vt);
		destroy_matrix(Qp);
	}

	for (i = 0; i < n; i++) {
		set_element(tmp, i, i, 0);
	}

	add_matrix(dst, tmp);

	destroy_matrix(tmp);
}
