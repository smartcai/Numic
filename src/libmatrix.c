#include "libmatrix.h"

/* Utilities */
#define FREE_SET_NULL(p) (free(p), p = NULL)

#define ERR_MSG(str)                                                    \
    printf("ERROR: %s: L%d, %s " str "\n", __FILE__, __LINE__, __FUNCTION__)

#define ASSERT(stmt, msg)                       \
do {                                            \
    if (!(stmt)) {                              \
        ERR_MSG(msg);                           \
        exit(-1);                               \
    }                                           \
} while (0)

#define ERR_MISMATCHED_TYPE ", mismatched type."
#define ERR_MISMATCHED_SIZE ", mismatched size."
#define ERR_DIVIDED_BY_ZERO ", divided by zero."


/* Machine EPSILON -- IEEE754 */
#define SCALAR_EPSILON    (pow(2.0, -52))
#define SCALAR_ABS(x)     (x < 0 ? -x : x)
#define SCALAR_EQLZERO(x) ((SCALAR_ABS(x)) <= (SCALAR_EPSILON))
#define SCALAR_NOTZERO(x) ((SCALAR_ABS(x)) >  (SCALAR_EPSILON))
#define SCALAR_MIN(a, b)  (a > b ? b : a)
#define SCALAR_MAX(a, b)  (a < b ? b : a)

/* Store in column-major order */
#define COL_MAJOR_ORDER(i, j, row) (i + j * row)
#define ROW_MAJOR_ORDER(i, j, col) (i * col + j)

struct Vector
{
    bTypeV  t;                  /* Vector type */
    int     dim;                /* Dimention */
    Scalar *arr;                /* Elements */
};

struct Matrix
{
    int     row;                /* Row number */
    int     col;                /* Column number */
    Scalar *arr;                /* Elements */
};

Vector *vector_new(int dim, bTypeV t)
{
    Vector *v = (Vector*) malloc (sizeof(Vector));

    v->t   = t;
    v->dim = dim;
    v->arr = (Scalar*) calloc (dim, sizeof(Scalar));

    return v;
}

void vector_del(Vector *v)
{
    if (v != NULL)
        FREE_SET_NULL(v);
    FREE_SET_NULL(v);
}

void vector_set(Vector *v, int k, Scalar val)
{
    v->arr[k] = val;
}

Scalar vector_get(Vector *v, int k)
{
    return v->arr[k];
}

int vector_dim(Vector *v)
{
    return v->dim;
}

bTypeV vector_type(Vector *v)
{
    return v->t;
}

void vector_print(Vector *v)
{
    int k;
    if (v->t == COL) {
        for (k = 0; k < v->dim; k++) {
            printf("%lf\n", vector_get(v, k));
        }
    }
    else {
        for (k = 0; k < v->dim; k++) {
            printf("%lf ", vector_get(v, k));
        }
        printf("\n");
    }
    printf("\n");
}

void vector_cpy(Vector *dst, Vector *src)
{
    int n;

    ASSERT(dst->dim == src->dim, ERR_MISMATCHED_SIZE);

    n = dst->dim; dst->t = src->t;
    memcpy(dst->arr, src->arr, n * sizeof(Scalar));
}


Matrix *matrix_new(int row, int col)
{
    Matrix *m = (Matrix *) malloc (sizeof(Matrix));

    m->row = row; m->col = col;
    m->arr = (Scalar*) calloc (row * col, sizeof(Scalar));

    return m;
}

void matrix_del(Matrix *m)
{
    if (m != NULL)
        FREE_SET_NULL(m->arr);
    FREE_SET_NULL(m->arr);
}

int matrix_row(Matrix *m)
{
    return m->row;
}

int matrix_col(Matrix *m)
{
    return m->col;
}

void matrix_set(Matrix *m, int i, int j, Scalar val)
{
    int idx     = COL_MAJOR_ORDER(i, j, m->row);
    m->arr[idx] = val;
}

Scalar matrix_get(Matrix *m, int i, int j)
{
    int idx = COL_MAJOR_ORDER(i, j, m->row);
    return m->arr[idx];
}

void matrix_print(Matrix *m)
{
    int i, j;
    for (i = 0; i < m->row; i++) {
        for (j = 0; j < m->col; j++) {
            printf("%lf ", matrix_get(m, i, j));
        }
        printf("\n");
    }
    printf("\n");
}

void matrix_cpy(Matrix *dst, Matrix *src)
{
    int m, n;

    ASSERT(dst->row == src->row, ERR_MISMATCHED_SIZE);
    ASSERT(dst->col == src->col, ERR_MISMATCHED_SIZE);

    m = dst->row; n = dst->col;
    memcpy(dst->arr, src->arr, m * n * sizeof(Scalar));
}

void matrix_save(Matrix *m, char *fn)
{
    FILE *f;
    int i, j;

    if ((f = fopen(fn, "w")) == NULL) {
        printf("Cannot open file. %s, %d\n", fn, __LINE__);
        exit(-1);
    }

    for (i = 0; i < m->row; i++) {
        for (j = 0; j < m->col; j++) {
            Scalar val = matrix_get(m, i, j);
            fprintf(f, "%.15lf,  ", val);
        }
        fprintf(f, "\n");
    }

    fclose(f);
}

void matrix_load(Matrix *m, char *fn)
{
    FILE *f;
    int i, j;

    if ((f = fopen(fn, "r")) == NULL) {
        printf("Cannot open file. %s, %d\n", fn, __LINE__);
        exit(-1);
    }

    for (i = 0; i < m->row; i++) {
        for (j = 0; j < m->col; j++) {
            Scalar val;
            fscanf(f, "%lf,", &val);
            matrix_set(m, i, j, val);
        }
    }

    fclose(f);
}


void blas_saxpy(Vector *y, Scalar  a, Vector *x)
/* y = y + ax: Scalar A X Plus Y */
{
    int k, n;

    ASSERT(y->t == x->t, ERR_MISMATCHED_TYPE);
    ASSERT(y->dim == x->dim, ERR_MISMATCHED_SIZE);

    n = y->dim;
    for (k = 0; k < n; k++) {
        y->arr[k] = y->arr[k] + a * x->arr[k];
    }
}

void blas_gaxpy(Vector *y, Matrix *A, Vector *x)
/* y = y + Ax: Generalized sAXPY */
{
    int k, m, n; Vector *v;

    ASSERT(x->t == COL, ERR_MISMATCHED_TYPE);
    ASSERT(A->col == x->dim, ERR_MISMATCHED_SIZE);

    m = A->row; n =A->col;
    v = vector_new(m, COL);
    for (k = 0; k < n; k++) {
        matrix_get_col_vector(v, A, k);
        blas_saxpy(y, x->arr[k], v);
    }
    vector_del(v);
}

void blas_gemul(Matrix *C, Matrix *A, Matrix *B)
/* C = C + AB: GEneral matrix-matrix MULtiplication */
{
    int j, m, n; Vector *cj, *bj;

    ASSERT(C->row == A->row, ERR_MISMATCHED_SIZE);
    ASSERT(C->col == B->col, ERR_MISMATCHED_SIZE);
    ASSERT(A->col == B->row, ERR_MISMATCHED_SIZE);

    m = C->row; n = C->col;
    cj = vector_new(m, COL);
    bj = vector_new(B->row, COL);
    for (j = 0; j < n; j++) {
        matrix_get_col_vector(cj, C, j);
        matrix_get_col_vector(bj, B, j);
        blas_gaxpy(cj, A, bj);
        matrix_set_col_vector(C, j, cj);
    }
    vector_del(cj);
    vector_del(bj);
}

void vector_transpose(Vector *y, Vector *x)
/* y = transpose(x), a.k.a., y = x^t */
{
    int  n;

    ASSERT(y->t != x->t, ERR_MISMATCHED_TYPE);
    ASSERT(y->dim == x->dim, ERR_MISMATCHED_SIZE);

    n = y->dim;
    memcpy(y->arr, x->arr, n * sizeof(Scalar));
}

Scalar vector_dot(Vector *x, Vector *y)
/* return dot(x, y) */
{
    Scalar s;
    int k, n;

    ASSERT(x->t == y->t, ERR_MISMATCHED_TYPE);
    ASSERT(x->dim == y->dim, ERR_MISMATCHED_SIZE);

    s = 0; n = x->dim;
    for (k = 0; k < n; k++) {
        s = s + x->arr[k] * y->arr[k];
    }

    return s;
}

void vector_outer(Matrix *A, Vector *x, Vector *y)
/* A = xy^t */
{
    int m, n, i, j;

    ASSERT(A->row == x->dim, ERR_MISMATCHED_SIZE);
    ASSERT(A->col == y->dim, ERR_MISMATCHED_SIZE);
    ASSERT(x->t == COL, ERR_MISMATCHED_TYPE);
    ASSERT(y->t == COL, ERR_MISMATCHED_TYPE);

    m = A->row; n = A->col;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            Scalar val = vector_get(x, i) * vector_get(y, j);
            matrix_set(A, i, j, val);
        }
    }
}

void vector_scalar_mul(Vector *y, Scalar a, Vector *x)
/* y = ax */
{
    int k, n;

    ASSERT(x->t == y->t, ERR_MISMATCHED_TYPE);
    ASSERT(x->dim == y->dim, ERR_MISMATCHED_SIZE);

    n = x->dim;
    for (k = 0; k < n; k++) {
        y->arr[k] = a * x->arr[k];
    }
}

void vector_scalar_add(Vector *y, Scalar a, Vector *x)
/* y(k) = x(k) + a */
{
    int k, n;

    ASSERT(x->t == y->t, ERR_MISMATCHED_TYPE);
    ASSERT(x->dim == y->dim, ERR_MISMATCHED_SIZE);

    n = x->dim;
    for (k = 0; k < n; k++) {
        y->arr[k] = x->arr[k] + a;
    }
}

void vector_add(Vector *z, Vector *x, Vector *y)
/* z = x + y */
{
    int k, n;

    ASSERT(z->t == x->t, ERR_MISMATCHED_TYPE);
    ASSERT(z->t == y->t, ERR_MISMATCHED_TYPE);
    ASSERT(z->dim == x->dim, ERR_MISMATCHED_SIZE);
    ASSERT(z->dim == y->dim, ERR_MISMATCHED_SIZE);

    n = x->dim;
    for (k = 0; k < n; k++) {
        z->arr[k] = x->arr[k] + y->arr[k];
    }
}

void vector_sub(Vector *z, Vector *x, Vector *y)
/* z = x - y */
{
    int k, n;

    ASSERT(z->t == x->t, ERR_MISMATCHED_TYPE);
    ASSERT(z->t == y->t, ERR_MISMATCHED_TYPE);
    ASSERT(z->dim == x->dim, ERR_MISMATCHED_SIZE);
    ASSERT(z->dim == y->dim, ERR_MISMATCHED_SIZE);

    n = x->dim;
    for (k = 0; k < n; k++) {
        z->arr[k] = x->arr[k] - y->arr[k];
    }
}

void vector_point_mul(Vector *z, Vector *x, Vector *y)
/* z(k) = x(k) * y(k) */
{
    int k, n;

    ASSERT(z->t == x->t, ERR_MISMATCHED_TYPE);
    ASSERT(z->t == y->t, ERR_MISMATCHED_TYPE);
    ASSERT(z->dim == x->dim, ERR_MISMATCHED_SIZE);
    ASSERT(z->dim == y->dim, ERR_MISMATCHED_SIZE);

    n = x->dim;
    for (k = 0; k < n; k++) {
        z->arr[k] = x->arr[k] * y->arr[k];
    }
}

void vector_point_div(Vector *z, Vector *x, Vector *y)
/* z(k) = x(k) / y(k), y(k) != 0 */
{
    int k, n;

    ASSERT(z->t == x->t, ERR_MISMATCHED_TYPE);
    ASSERT(z->t == y->t, ERR_MISMATCHED_TYPE);
    ASSERT(z->dim == x->dim, ERR_MISMATCHED_SIZE);
    ASSERT(z->dim == y->dim, ERR_MISMATCHED_SIZE);
    ASSERT(vector_all(y), ERR_DIVIDED_BY_ZERO);

    n = x->dim;
    for (k = 0; k < n; k++) {
        z->arr[k] = x->arr[k] / y->arr[k];
    }
}

int vector_any(Vector *x)
/* return true, when x has any nonzero. */
{
    int k, n, bool;

    bool = 0; n = x->dim;
    for (k = 0; k < n; k++) {
        if (SCALAR_NOTZERO(x->arr[k])) {
            bool = 1;
            break;
        }
    }

    return bool;
}

int vector_all(Vector *x)
/* return true, when x has all nonzero. */
{
    int k, n, bool;

    bool = 1; n = x->dim;
    for (k = 0; k < n; k++) {
        if (SCALAR_EQLZERO(x->arr[k])) {
            bool = 0;
            break;
        }
    }

    return bool;
}

void vector_zero(Vector *x)
{
    int k = x->dim;
    memset(x->arr, 0, k * sizeof(Scalar));
}

void matrix_transpose(Matrix *C, Matrix *A)
/* C = transpose(A), a.k.a., C = A^t */
{
    int m, n, i, j;

    ASSERT(C->row == A->col, ERR_MISMATCHED_SIZE);
    ASSERT(C->col == A->row, ERR_MISMATCHED_SIZE);

    m = C->row; n = C->col;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            Scalar val = matrix_get(A, j, i);
            matrix_set(C, i, j, val);
        }
    }
}

void matrix_scalar_mul(Matrix *C, Scalar a, Matrix *A)
{
    int m, n, i, j;

    ASSERT(C->row == A->row, ERR_MISMATCHED_SIZE);
    ASSERT(C->col == A->col, ERR_MISMATCHED_SIZE);

    m = C->row; n = C->col;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            Scalar val = a * matrix_get(A, i, j);
            matrix_set(C, i, j, val);
        }
    }
}

void matrix_scalar_add(Matrix *C, Scalar a, Matrix *A)
{
    int m, n, i, j;

    ASSERT(C->row == A->row, ERR_MISMATCHED_SIZE);
    ASSERT(C->col == A->col, ERR_MISMATCHED_SIZE);

    m = C->row; n = C->col;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            Scalar val = a + matrix_get(A, i, j);
            matrix_set(C, i, j, val);
        }
    }
}

void matrix_add(Matrix *C, Matrix *A, Matrix *B)
{
    int m, n, i, j;

    ASSERT(C->row == A->row, ERR_MISMATCHED_SIZE);
    ASSERT(C->col == A->col, ERR_MISMATCHED_SIZE);
    ASSERT(C->row == B->row, ERR_MISMATCHED_SIZE);
    ASSERT(C->col == B->col, ERR_MISMATCHED_SIZE);

    m = C->row; n = C->col;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            Scalar val = matrix_get(A, i, j) + matrix_get(B, i, j);
            matrix_set(C, i, j, val);
        }
    }
}

void matrix_sub(Matrix *C, Matrix *A, Matrix *B)
{
    int m, n, i, j;

    ASSERT(C->row == A->row, ERR_MISMATCHED_SIZE);
    ASSERT(C->col == A->col, ERR_MISMATCHED_SIZE);
    ASSERT(C->row == B->row, ERR_MISMATCHED_SIZE);
    ASSERT(C->col == B->col, ERR_MISMATCHED_SIZE);

    m = C->row; n = C->col;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            Scalar val = matrix_get(A, i, j) - matrix_get(B, i, j);
            matrix_set(C, i, j, val);
        }
    }
}

void matrix_mul(Matrix *C, Matrix *A, Matrix *B)
{
    matrix_zero(C);
    blas_gemul(C, A, B);
}

void matrix_point_mul(Matrix *C, Matrix *A, Matrix *B)
{
    int m, n, i, j;

    ASSERT(C->row == A->row, ERR_MISMATCHED_SIZE);
    ASSERT(C->col == A->col, ERR_MISMATCHED_SIZE);
    ASSERT(C->row == B->row, ERR_MISMATCHED_SIZE);
    ASSERT(C->col == B->col, ERR_MISMATCHED_SIZE);

    m = C->row; n = C->col;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            Scalar val = matrix_get(A, i, j) * matrix_get(B, i, j);
            matrix_set(C, i, j, val);
        }
    }
}

void matrix_point_div(Matrix *C, Matrix *A, Matrix *B)
{
    int m, n, i, j;

    ASSERT(C->row == A->row, ERR_MISMATCHED_SIZE);
    ASSERT(C->col == A->col, ERR_MISMATCHED_SIZE);
    ASSERT(C->row == B->row, ERR_MISMATCHED_SIZE);
    ASSERT(C->col == B->col, ERR_MISMATCHED_SIZE);
    ASSERT(matrix_any(B), ERR_DIVIDED_BY_ZERO);

    m = C->row; n = C->col;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            Scalar val = matrix_get(A, i, j) / matrix_get(B, i, j);
            matrix_set(C, i, j, val);
        }
    }
}

void matrix_get_col_vector(Vector *y, Matrix *A, int j)
/* y = A(:, j) */
{
    int k, n;

    ASSERT(y->t == COL, ERR_MISMATCHED_TYPE);
    ASSERT(y->dim == A->row, ERR_MISMATCHED_SIZE);

    n = y->dim;
    for (k = 0; k < n; k++) {
        y->arr[k] = matrix_get(A, k, j);
    }
}

void matrix_get_row_vector(Vector *y, Matrix *A, int i)
/* y = A(i, :) */
{
    int k, n;

    ASSERT(y->t == ROW, ERR_MISMATCHED_TYPE);
    ASSERT(y->dim == A->col, ERR_MISMATCHED_SIZE);

    n = y->dim;
    for (k = 0; k < n; k++) {
        y->arr[k] = matrix_get(A, i, k);
    }
}

void matrix_set_col_vector(Matrix *A, int j, Vector *x)
/* A(:, j) = x */
{
    int k, n;

    ASSERT(x->t == COL, ERR_MISMATCHED_TYPE);
    ASSERT(x->dim == A->row, ERR_MISMATCHED_SIZE);

    n = x->dim;
    for (k = 0; k < n; k++) {
        matrix_set(A, k, j, x->arr[k]);
    }
}

void matrix_set_row_vector(Matrix *A, int i, Vector *x)
/* A(i, :) = x */
{
    int k, n;

    ASSERT(x->t == ROW, ERR_MISMATCHED_TYPE);
    ASSERT(x->dim == A->col, ERR_MISMATCHED_SIZE);

    n = x->dim;
    for (k = 0; k < n; k++) {
        matrix_set(A, i, k, x->arr[k]);
    }
}

void matrix_col_any(Vector *y, Matrix *A)
{
    int k, m, n; Vector *v;

    ASSERT(y->t == ROW, ERR_MISMATCHED_TYPE);
    ASSERT(y->dim == A->col, ERR_MISMATCHED_SIZE);

    m = A->row; n = A->col;
    v = vector_new(m, COL);

    for (k = 0; k < n; k++) {
        matrix_get_col_vector(v, A, k);
        y->arr[k] = vector_any(v);
    }

    vector_del(v);
}

void matrix_col_all(Vector *y, Matrix *A)
{
    int k, m, n; Vector *v;

    ASSERT(y->t == ROW, ERR_MISMATCHED_TYPE);
    ASSERT(y->dim == A->col, ERR_MISMATCHED_SIZE);

    m = A->row; n = A->col;
    v = vector_new(m, COL);

    for (k = 0; k < n; k++) {
        matrix_get_col_vector(v, A, k);
        y->arr[k] = vector_all(v);
    }

    vector_del(v);
}

void matrix_row_any(Vector *y, Matrix *A)
{
    int k, m, n; Vector *v;

    ASSERT(y->t == COL, ERR_MISMATCHED_TYPE);
    ASSERT(y->dim == A->row, ERR_MISMATCHED_SIZE);

    m = A->row; n = A->col;
    v = vector_new(n, ROW);

    for (k = 0; k < m; k++) {
        matrix_get_row_vector(v, A, k);
        y->arr[k] = vector_any(v);
    }

    vector_del(v);
}

void matrix_row_all(Vector *y, Matrix *A)
{
    int k, m, n; Vector *v;

    ASSERT(y->t == COL, ERR_MISMATCHED_TYPE);
    ASSERT(y->dim == A->row, ERR_MISMATCHED_SIZE);

    m = A->row; n = A->col;
    v = vector_new(n, ROW);

    for (k = 0; k < m; k++) {
        matrix_get_row_vector(v, A, k);
        y->arr[k] = vector_all(v);
    }

    vector_del(v);
}

int matrix_any(Matrix *A)
{
    int k, m, n, bool;

    bool = 0; m = A->row; n = A->col;
    for (k = 0; k < m * n; k++) {
        if (SCALAR_NOTZERO(A->arr[k])) {
            bool = 1;
            break;
        }
    }

    return bool;
}

int matrix_all(Matrix *A)
{
    int k, m, n, bool;

    bool = 1; m = A->row; n = A->col;
    for (k = 0; k < m * n; k++) {
        if (SCALAR_EQLZERO(A->arr[k])) {
            bool = 0;
            break;
        }
    }

    return bool;
}

void matrix_zero(Matrix *A)
/* A(i, j) = 0 */
{
    int m, n;
    m = A->row; n = A->col;
    memset(A->arr, 0, m * n * sizeof(Scalar));
}

void matrix_eye(Matrix *A)
/* A(k, k) = 1 */
{
    int k, m, n;

    matrix_zero(A);
    m = A->row; n = A->col;
    for (k = 0; k < SCALAR_MIN(m, n); k++) {
        matrix_set(A, k, k, 1);
    }
}

Scalar vector_L2_norm(Vector *v)
{
    Scalar s;
    s = sqrt(vector_dot(v, v));
    return s;
}

void vector_housh(Vector *v, Scalar *beta, Vector *x, int j)
/**
 * Rotate x to the dimension of j.
 * (I - beta*v*v^t)x =   norm(x) * e(j) if x(j) <  0
 * (I - beta*v*v^t)x = - norm(x) * e(j) if x(j) >= 0
 */
{
    Scalar sigma, L2, xj, vj;

    ASSERT(v->t == x->t, ERR_MISMATCHED_TYPE);
    ASSERT(v->dim == x->dim, ERR_MISMATCHED_SIZE);

    xj = vector_get(x, j);
    L2 = vector_L2_norm(x);
    sigma = L2 * L2 - xj * xj;
    vector_cpy(v, x);
    vector_set(v, j, 1);

    if (SCALAR_EQLZERO(sigma)) {
        if (xj >= 0)
            *beta = 0;
        else
            *beta = 2;
    }
    else {
        if (xj <= 0)
            vj = xj - L2;
        else
            vj = - sigma / (xj + L2);
        vector_set(v, j, vj);
        *beta = 2 * vj * vj / (sigma + vj * vj);
        vector_scalar_mul(v, 1.0 / vj, v);
    }
}

void vector_givens(Scalar *c, Scalar *s, Scalar y, Scalar z)
/**
 * Compute a 2x2 matrix G to rotate angle theta such that:
 * G^t * [y; z] = [r; 0],
 * where G = [c s; -s c]
 *       c = cos(theta)
 *       s = sin(theta)
 * G * G^t = I
 */
{
    Scalar r;

    if (SCALAR_EQLZERO(y)) {
        *c = 1;
        *s = 0;
    }
    else {
        if (SCALAR_ABS(z) > SCALAR_ABS(y)) {
            r = -y / z;
            *s = 1 / sqrt(1 + r * r);
            *c = (*s) * r;
        }
        else {
            r = -z / y;
            *c = 1 / sqrt(1 + r * r);
            *s = (*c) * r;
        }
    }
}

void matrix_full_qr_housh(Matrix *A, Matrix *Q, Matrix *R)
/* A = Q * R */
{
    int m, n, i, j, k;
    Matrix *Hj, *Qt, *Rt;
    Vector *xj, *vj;
    Scalar beta;

    ASSERT(A->row == R->row, ERR_MISMATCHED_SIZE);
    ASSERT(A->col == R->col, ERR_MISMATCHED_SIZE);
    ASSERT(Q->col == R->row, ERR_MISMATCHED_SIZE);
    ASSERT(Q->row == Q->col, ERR_MISMATCHED_SIZE);

    m = A->row; n = A->col;
    k = m < n ? m : n;
    Hj = matrix_new(m, m);
    Qt = matrix_new(m, m);
    Rt = matrix_new(m, n);
    xj = vector_new(m, COL);
    vj = vector_new(m, COL);
    matrix_cpy(R, A);
    matrix_eye(Q);
    for (j = 0; j < k; j++) {
        vector_zero(xj);
        for (i = j; i < m; i++) {
            xj->arr[i] = matrix_get(R, i, j);
        }

        vector_housh(vj, &beta, xj, j);
        vector_outer(Hj, vj, vj);
        matrix_scalar_mul(Hj, -beta, Hj);
        for (i = 0; i < m; i++) {
            Scalar val = matrix_get(Hj, i, i) + 1;
            matrix_set(Hj, i, i, val);
        }

        matrix_mul(Qt, Q, Hj);
        matrix_cpy(Q, Qt);
        matrix_mul(Rt, Hj, R);
        matrix_cpy(R, Rt);
    }

    matrix_del(Hj);
    matrix_del(Qt);
    matrix_del(Rt);
    vector_del(xj);
    vector_del(vj);
}

void matrix_full_qr_givens(Matrix *A, Matrix *Q, Matrix *R)
{
    ;
}

void matrix_full_qr_cgs(Matrix *A, Matrix *Q, Matrix *R)
{
    ;
}

void matrix_full_qr_mgs(Matrix *A, Matrix *Q, Matrix *R)
{
    ;
}

void matrix_bidiagonal(Matrix *A, Matrix *P, Matrix *B, Matrix *Q)
/* A = PBQ */
{
    int m, n, i, j, k, ii, jj;
    Matrix *Hj, *Hi, *Pt, *Bt, *Qt;
    Vector *xj, *vj, *xi, *vi;
    Scalar beta, c, s, y, z;

    ASSERT(A->row == B->row, ERR_MISMATCHED_SIZE);
    ASSERT(A->col == B->col, ERR_MISMATCHED_SIZE);
    ASSERT(P->row == P->col, ERR_MISMATCHED_SIZE);
    ASSERT(Q->row == Q->col, ERR_MISMATCHED_SIZE);
    ASSERT(B->row == P->col, ERR_MISMATCHED_SIZE);
    ASSERT(B->col == Q->row, ERR_MISMATCHED_SIZE);

    m = A->row; n = A->col;
    k = m < n ? m : n;
    Hj = matrix_new(m, m);
    Hi = matrix_new(n, n);
    Pt = matrix_new(m, m);
    Bt = matrix_new(m, n);
    Qt = matrix_new(n, n);
    xj = vector_new(m, COL);
    vj = vector_new(m, COL);
    xi = vector_new(n, COL);
    vi = vector_new(n, COL);
    matrix_cpy(B, A);
    matrix_eye(P);
    matrix_eye(Q);
    for (j = 0; j < k; j++) {
        if (1) {
            vector_zero(xj);
            for (i = j; i < m; i++) {
                xj->arr[i] = matrix_get(B, i, j);
            }
            vector_housh(vj, &beta, xj, j);
            vector_outer(Hj, vj, vj);
            matrix_scalar_mul(Hj, -beta, Hj);
            for (i = 0; i < m; i++) {
                Scalar val = matrix_get(Hj, i, i) + 1;
                matrix_set(Hj, i, i, val);
            }
            matrix_mul(Pt, P, Hj);
            matrix_cpy(P, Pt);
            matrix_mul(Bt, Hj, B);
            matrix_cpy(B, Bt);
        }

        i = j;
        if (j + 1 < n) {
            vector_zero(xi);
            for (jj = j + 1; jj < n; jj++) {
                xi->arr[jj] = matrix_get(B, i, jj);
            }
            vector_housh(vi, &beta, xi, j + 1);
            vector_outer(Hi, vi, vi);
            matrix_scalar_mul(Hi, -beta, Hi);
            for (ii = 0; ii < n; ii++) {
                Scalar val = matrix_get(Hi, ii, ii) + 1;
                matrix_set(Hi, ii, ii, val);
            }
            matrix_mul(Qt, Hi, Q);
            matrix_cpy(Q, Qt);
            matrix_mul(Bt, B, Hi);
            matrix_cpy(B, Bt);
        }
    }

    matrix_del(Hj);
    matrix_del(Hi);
    matrix_del(Pt);
    matrix_del(Bt);
    matrix_del(Qt);
    vector_del(xj);
    vector_del(vj);
    vector_del(xi);
    vector_del(vi);
}

void matrix_svd(Matrix *A, Matrix *U, Matrix *S, Matrix *V)
/* A = U * S * V^t */
{
    int m, n, i, j, k, r, iter, flag, f1, f2;
    Scalar c, s, y, z, y1, z1, sigma;
    Matrix *Vt;

    ASSERT(A->row == S->row, ERR_MISMATCHED_SIZE);
    ASSERT(A->col == S->col, ERR_MISMATCHED_SIZE);
    ASSERT(U->row == U->col, ERR_MISMATCHED_SIZE);
    ASSERT(V->row == V->col, ERR_MISMATCHED_SIZE);
    ASSERT(S->row == U->col, ERR_MISMATCHED_SIZE);
    ASSERT(S->col == V->row, ERR_MISMATCHED_SIZE);

    matrix_bidiagonal(A, U, S, V);

    m = A->row; n = A->col;
    r = m < n ? m : n;
    Vt = matrix_new(n, n);
    flag = 1;
    while (flag) {
        flag = 0;
        for (k = 0; k < r - 1; k++) {
            /* Right multiplication */
            y = matrix_get(S, k, k);
            z = matrix_get(S, k, k + 1);
            if (SCALAR_EQLZERO(z)) {
                f1 = 0;
                flag = flag || f1;
            }
            else {
                f1 = 1;
                flag = flag || f1;
            }
            vector_givens(&c, &s, y, z);
            for (i = 0; i < m; i++) {
                y = matrix_get(S, i, k);
                z = matrix_get(S, i, k + 1);
                y1 = c * y - s * z;
                z1 = s * y + c * z;
                matrix_set(S, i, k, y1);
                matrix_set(S, i, k + 1, z1);
            }
            for (j = 0; j < n; j++) {
                y = matrix_get(V, k, j);
                z = matrix_get(V, k + 1, j);
                y1 = c * y - s * z;
                z1 = s * y + c * z;
                matrix_set(V, k, j, y1);
                matrix_set(V, k + 1, j, z1);
            }

            /* Left multiplication */
            y = matrix_get(S, k, k);
            z = matrix_get(S, k + 1, k);
            if (SCALAR_EQLZERO(z)) {
                f1 = 0;
                flag = flag || f1;
            }
            else {
                f1 = 1;
                flag = flag || f1;
            }
            vector_givens(&c, &s, y, z);
            for (j = 0; j < n; j++) {
                y = matrix_get(S, k, j);
                z = matrix_get(S, k + 1, j);
                y1 = c * y - s * z;
                z1 = s * y + c * z;
                matrix_set(S, k, j, y1);
                matrix_set(S, k + 1, j, z1);
            }
            for (i = 0; i < m; i++) {
                y = matrix_get(U, i, k);
                z = matrix_get(U, i, k + 1);
                y1 = c * y - s * z;
                z1 = s * y + c * z;
                matrix_set(U, i, k, y1);
                matrix_set(U, i, k + 1, z1);
            }
        }

        /* Remove the tail element */
        if ((n > m) && (iter == 0)) {
            /* Right multiplication */
            y = matrix_get(S, r - 1, r - 1);
            z = matrix_get(S, r - 1, r);
            vector_givens(&c, &s, y, z);
            for (i = 0; i < m; i++) {
                y = matrix_get(S, i, r - 1);
                z = matrix_get(S, i, r);
                y1 = c * y - s * z;
                z1 = s * y + c * z;
                matrix_set(S, i, r - 1, y1);
                matrix_set(S, i, r, z1);
            }
            for (j = 0; j < n; j++) {
                y = matrix_get(V, k, j);
                z = matrix_get(V, k + 1, j);
                y1 = c * y - s * z;
                z1 = s * y + c * z;
                matrix_set(V, k, j, y1);
                matrix_set(V, k + 1, j, z1);
            }
        }
    }

    /* Change signs of eigenvalues */
    for (k = 0; k < r; k++) {
        sigma = matrix_get(S, k, k);
        if (sigma < 0) {
            matrix_set(S, k, k, -sigma);
            Scalar val;
            for (i = 0; i < m; i++) {
                val = matrix_get(U, i, k);
                matrix_set(U, i, k, -val);
            }
        }
    }

    matrix_transpose(Vt, V);
    matrix_cpy(V, Vt);

    matrix_del(Vt);
}
