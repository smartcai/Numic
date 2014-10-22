#include <stdio.h>
#include "test.h"

void test_matrix(void);
void test_vector(void);
void test_saxpy(void);
void test_gaxpy(void);
void test_matrix_mul(void);
void test_qr_decompose_cgs(void);
void test_householder_vector(void);
void test_house_matrix_columns(void);

int main(int argc, char *argv[])
{
	test_matrix();
	test_vector();
	test_saxpy();
	test_gaxpy();
	test_matrix_mul();
	test_qr_decompose_cgs();
	test_householder_vector();
	test_house_matrix_columns();
	return 0;
}

void test_matrix(void)
{
	int i, j;
	int rows = 4;
	int cols = 3;
	matrix *mp1 = create_matrix(rows, cols);
	matrix *mp2 = create_matrix(rows, cols);
	matrix *mpt = create_matrix(cols, rows);
	matrix *sub = create_matrix(rows-1, cols-1);
	vector *vec = create_col_vector(rows);

	printf("****test_matrix start****\n");
	printf("mp1\n");
	print_matrix(mp1);
	printf("mp2\n");
	print_matrix(mp2);
	printf("mpt\n");
	print_matrix(mpt);

	/* Pretend to use columns-major store */
	for (j = 0; j < cols ; j++) {
		for (i = 0; i < rows; i++) {
			set_element(mp1, i, j , i + j * rows);
		}
	}
	printf("mp1\n");
	print_matrix(mp1);

	copy_matrix(mp1, mp2);
	printf("mp2\n");
	print_matrix(mp2);

	get_block(mp1, 1, 1, sub);
	printf("sub\n");
	print_matrix(sub);

	set_block(mp2, 0, 1, sub);
	printf("mp2\n");
	print_matrix(mp2);

	printf("vec\n");
	print_vector(vec);
	get_col_vector(mp1, 2, vec);
	printf("vec\n");
	print_vector(vec);
	printf("vec\n");
	scalar_vector_mul(vec, 3.14, vec);
	print_vector(vec);

	set_col_vector(mp2, 0, vec);
	printf("mp2\n");
	print_matrix(mp2);

	printf("mp2\n");
	scalar_matrix_mul(mp2, 0.5, mp2);
	print_matrix(mp2);

	transpose(mp1, mpt);
	printf("mpt\n");
	print_matrix(mpt);

	zero_matrix(mp1);
	printf("mp1\n");
	print_matrix(mp1);

	copy_matrix(mp1, mp2);
	printf("mp2\n");
	print_matrix(mp2);

	printf("%d %d\n", get_rows(mp1), get_cols(mp1));

	printf("mpt\n");
	print_matrix(mpt);
	transpose(mpt, mpt);
	printf("mpt\n");
	print_matrix(mpt);

	destroy_matrix(mp1);
	destroy_matrix(mp2);
	destroy_matrix(mpt);
	destroy_matrix(sub);
	destroy_vector(vec);

	printf("****test_matrix end****\n");
}

void test_vector(void)
{
	int k;
	int dim = 5;
	scalar a = 3.14;
	vector *v1 = create_col_vector(dim);
	vector *v2 = create_col_vector(dim);
	vector *vt = create_row_vector(dim);

	printf("****test_vector start****\n");
	print_vector(v1);
	print_vector(v2);
	print_vector(vt);

	for (k = 0; k < dim; k++) {
		set_vector_element(v1, k, k);
	}
	printf("v1\n");
	print_vector(v1);

	printf("Test saxpy:\n");
	saxpy(v1, a, v1);
	print_vector(v1);
	saxpy(v2, 2, v1);
	print_vector(v2);

	copy_vector(v1, v2);
	print_vector(v2);

	printf("%f\n", dot_product(v1, v2));
	printf("%f\n", vector_norm(v1));
	printf("%f\n", vector_norm(v2));

	transpose_vector(v1, vt);
	print_vector(vt);

	zero_vector(v1);
	copy_vector(v1, v2);

	print_vector(v2);

	printf("%d\n", get_dim(v1));

	destroy_vector(v1);
	destroy_vector(v2);
	destroy_vector(vt);

	printf("****test_vector end****\n");
}

void test_saxpy(void)
{
	printf("**** Test saxpy ****\n");
	vector *y, *x;
	scalar a = 3.14;

	int n = 10, k;

	y = create_col_vector(n);
	x = create_col_vector(n);
	for (k = 0; k < n; k++) {
		set_vector_element(y, k, 1);
		set_vector_element(x, k, k);
	}

	printf("y\n");
	print_vector(y);
	printf("x\n");
	print_vector(x);
	saxpy(y, a, x);
	printf("y after saxpy\n");
	print_vector(y);

	destroy_vector(y);
	destroy_vector(x);
	printf("**** Test saxpy end ****\n");
}

void test_gaxpy(void)
{
	printf("**** Test gaxpy ****\n");
	vector *y, *x;
	matrix *A;

	int n = 10, m = 5, i, j, k;

	y = create_col_vector(m);
	x = create_col_vector(n);
	A = create_matrix(m, n);

	for (k = 0; k < m; k++) {
		set_vector_element(y, k, 1);
	}

	for (k = 0; k < n; k++) {
		set_vector_element(x, k, k);
	}

	for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
			set_element( A, i, j, (i + j) / 10.0 );
		}
	}


	printf("A\n");
	print_matrix(A);
	printf("y\n");
	print_vector(y);
	printf("x\n");
	print_vector(x);
	gaxpy(y, A, x);
	printf("y after gaxpy\n");
	print_vector(y);

	destroy_vector(y);
	destroy_vector(x);
	destroy_matrix(A);
	printf("**** Test gaxpy end ****\n");
}

void test_matrix_mul(void)
{
	printf("**** Test matrix_mul ****\n");
	int m = 5, r = 3, n = 4;
	matrix *A = create_matrix(m, r);
	matrix *B = create_matrix(r, n);
	matrix *C = create_matrix(m, n);

	int i, j, k;

	for (j = 0; j < r; j++) {
		for (i = 0; i < m; i++) {
			set_element( A, i, j, (i + j) / 10.0 );
		}
	}
	printf("A\n");
	print_matrix(A);

	for (j = 0; j < n; j++) {
		for (i = 0; i < r; i++) {
			set_element( B, i, j, (2 * i + j) / 10.0 );
		}
	}
	printf("B\n");
	print_matrix(B);

	matrix_mul(C, A, B);
	printf("C\n");
	print_matrix(C);

	destroy_matrix(A);
	destroy_matrix(B);
	destroy_matrix(C);
	printf("**** Test matrix_mul end ****\n");
}

void test_qr_decompose_cgs(void)
{
	printf("**** Test qr_decompose_cgs ****\n");
	int m = 3, r = 2, n = 2;
	matrix *A = create_matrix(m, n);
	matrix *Q = create_matrix(m, n);
	matrix *R = create_matrix(n, n);
	vector *q1 = create_col_vector(m);
	vector *q2 = create_col_vector(m);
	scalar c = 0;

	int i, j, k;

	for (j = 0; j < r; j++) {
		for (i = 0; i < m; i++) {
			set_element( A, i, j, (i + j) );
		}
	}
	printf("A\n");
	print_matrix(A);

	qr_decompose_cgs(A, Q, R);
	printf("Q\n");
	print_matrix(Q);
	printf("R\n");
	print_matrix(R);

	matrix_mul(A, Q, R);
	printf("A\n");
	print_matrix(A);

	get_col_vector(Q, 0, q1);
	get_col_vector(Q, 1, q2);

	c = dot_product(q1, q2);
	printf("dot(q1, q2): %f\n", c);

	destroy_matrix(A);
	destroy_matrix(Q);
	destroy_matrix(R);
	destroy_vector(q1);
	destroy_vector(q2);
	printf("**** Test qr_decompose_cgs end ****\n");
}

void test_householder_vector(void)
{
	printf("**** Test householder_vector ****\n");
	vector *x, *v, *y;
	scalar beta;
	int k;
	int m = 7, i;

	matrix *I;
	matrix *prod;
	x = create_col_vector(m);
	v = create_col_vector(m);
	y = create_col_vector(m);
	I = create_matrix(m, m);
	prod = create_matrix(m, m);

	for (i = 0; i < m; i++) {
		set_element(I, i, i, 1);
	}

	printf("I\n");
	print_matrix(I);

	for (i = 0; i < m; i++) {
		set_vector_element(x, i, pow(-1, i) * sqrt(i));
	}

	printf("x\n");
	print_vector(x);
	householder_vector(x, v, &beta, 3);
	printf("v\n");
	print_vector(v);
	printf("beta: %f\n", beta);

	out_product(prod, v, v);
	printf("prod\n");
	print_matrix(prod);

	scalar_matrix_mul(prod, beta, prod);
	subtract_matrix(I, prod);
	printf("I\n");
	print_matrix(I);

	gaxpy(y, I, x);
	printf("y\n");
	print_vector(y);

	printf("norm x:%f\n", vector_norm(x));

	destroy_vector(x);
	destroy_vector(v);
	destroy_matrix(I);
	destroy_matrix(prod);
	destroy_matrix(y);
	printf("**** Test householder_vector end ****\n");
}

void test_house_matrix_columns(void)
{
	printf("**** Test house_matrix_columns ****\n");
	int m = 5, n = 3;
	matrix *A, *H;
	int i, j;
	A = create_matrix(m, n);
	H = create_matrix(m, n);

	for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
			set_element(A, i, j, i + j);
		}
	}

	house_matrix_columns(A, H);

	printf("A\n");
	print_matrix(A);
	printf("H\n");
	print_matrix(H);

	destroy_matrix(A);
	destroy_matrix(H);
	printf("**** Test house_matrix_columns end ****\n");
}
