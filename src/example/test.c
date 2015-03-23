#include <stdio.h>
#include <stdlib.h>

#include "libmatrix.h"

Scalar arr1[15] = {
     0.3 ,  0.6,  0.9,  1.2,  1.5,
    -0.4 , -0.7, -1.0, -1.3, -1.6,
     0.2 ,  0.8, -0.4,  1.4, -1.0
};

Scalar arr2[15] = {
     2.23,  2.87, 3.96,  1.25, -1.11,
    -9.78, -1.29, 3.14, -0.78,  0.05,
     1.00,  2.71, 2.00, -4.80, -1.11
};


void test_types(void)
{
    int i, j, k, cnt;
    Vector *v;
    Matrix *m;

    printf("Testing basic types.\n");

    v = vector_new(15, COL);
    m = matrix_new(3, 5);

    printf("Testing new, get and print.\n");
    vector_print(v);
    matrix_print(m);

    printf("Testing set a vector.\n");
    cnt = 0;
    for (k = 0; k < vector_dim(v); k++) {
        vector_set(v, k, arr1[cnt++]);
    }
    vector_print(v);

    printf("Testing set a matrix.\n");
    cnt = 0;
    for (i = 0; i < matrix_row(m); i++) {
        for (j = 0; j < matrix_col(m); j++) {
            matrix_set(m, i, j, arr1[cnt++]);
        }
    }
    matrix_print(m);

    vector_del(v);
    matrix_del(m);
}

void test_optns(void)
{
    int k;
    Scalar  s1;
    Vector *v1, *v2;

    printf("Testing operations.\n");

    v1 = vector_new(15, COL);
    v2 = vector_new(15, COL);
    for (k = 0; k < 15; k++) {
        vector_set(v1, k, arr1[k]);
        vector_set(v2, k, arr2[k]);
    }

    s1 = vector_dot(v1, v2);
    printf("dot(v1, v2) = %lf\n", s1);

    vector_del(v1);
    vector_del(v2);
}

void test_vector_housh(void)
{
    Vector *x, *v, *e;
    Matrix *A, *I;
    Scalar beta;
    int k, n;

    printf("Testing vector_housh.\n");

    x = vector_new(15, COL);
    v = vector_new(15, COL);
    e = vector_new(15, COL);
    A = matrix_new(15, 15);
    I = matrix_new(15, 15);

    n = 15;
    for (k = 0; k < n; k++) {
        vector_set(x, k, arr1[k]);
    }

    vector_housh(v, &beta, x, 3);
    vector_outer(A, v, v);
    matrix_scalar_mul(A, beta, A);
    matrix_eye(I);
    matrix_sub(I, I, A);
    blas_gaxpy(e, I, x);
    printf("beta:%lf\n\n", beta);
    vector_print(x);
    vector_print(v);
    vector_print(e);
    printf("L2 norm:%lf\n\n", vector_L2_norm(x));

    vector_del(x);
    vector_del(v);
    vector_del(e);
    matrix_del(A);
    matrix_del(I);
}

void test_vector_givens(void)
{
    Matrix *G;
    Vector *x, *y;
    Scalar c, s;

    x = vector_new(2, COL);
    y = vector_new(2, COL);
    G = matrix_new(2, 2);

    vector_set(x, 0, 1);
    vector_set(x, 1, 1);
    vector_givens(&c, &s, 1, 1);
    matrix_set(G, 0, 0, c);
    matrix_set(G, 0, 1, -s);
    matrix_set(G, 1, 0, s);
    matrix_set(G, 1, 1, c);
    blas_gaxpy(y, G, x);
    matrix_print(G);
    vector_print(y);

    matrix_del(G);
    vector_del(x);
    vector_del(y);
}

void test_matrix_full_qr_house(void)
{
    Matrix *A, *Q, *R, *B, *I, *Qt;
    int i, j, k;

    A  = matrix_new(5, 3);
    B  = matrix_new(5, 3);
    Q  = matrix_new(5, 5);
    R  = matrix_new(5, 3);
    I  = matrix_new(5, 5);
    Qt = matrix_new(5, 5);

    k = 0;
    for (i = 0; i < 5; i++) {
        for (j = 0; j < 3; j++) {
            matrix_set(A, i, j, arr1[k++]);
        }
    }

    matrix_full_qr_housh(A, Q, R);
    matrix_mul(B, Q, R);

    matrix_print(A);
    matrix_print(B);
    matrix_print(Q);
    matrix_print(R);

    matrix_transpose(Qt, Q);
    matrix_mul(I, Qt, Q);
    matrix_print(I);
    matrix_mul(I, Q, Qt);
    matrix_print(I);

    matrix_del(A);
    matrix_del(B);
    matrix_del(Q);
    matrix_del(R);
    matrix_del(I);
    matrix_del(Qt);
}

void test_matrix_bidiagonal_housh(void)
{
    Matrix *A, *P, *B, *Q, *Bt;
    int i, j, k;

    A  = matrix_new(3, 5);
    B  = matrix_new(3, 5);
    Bt = matrix_new(3, 5);
    P  = matrix_new(3, 3);
    Q  = matrix_new(5, 5);

    k = 0;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 5; j++) {
            matrix_set(A, i, j, arr1[k++]);
        }
    }
    matrix_bidiagonal(A, P, B, Q);
    matrix_print(P);
    matrix_print(B);
    matrix_print(Q);

    matrix_mul(Bt, B, Q);
    matrix_cpy(B, Bt);
    matrix_mul(Bt, P, B);
    matrix_print(A);
    matrix_print(Bt);

    matrix_del(A);
    matrix_del(P);
    matrix_del(B);
    matrix_del(Q);
    matrix_del(Bt);
}

int main(int argc, char *argv[])
{
    printf("Hello libmatrix!\n");
    printf("This is the test.\n");

    /* test_types(); */
    /* test_optns(); */
    /* test_vector_housh(); */
    /* test_vector_givens(); */
    /* test_matrix_full_qr_house(); */
    test_matrix_bidiagonal_housh();

    return 0;
}
