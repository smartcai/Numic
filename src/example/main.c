#include <stdio.h>
#include "test.h"

void test_matrix(void);
void test_vector(void);

int main(int argc, char *argv[])
{
	test_matrix();
	test_vector();
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

	printf("****test_matrix start****\n");
	print_matrix(mp1);
	print_matrix(mp2);
	print_matrix(mpt);

	/* Pretend to use columns-major store */
	for (j = 0; j < cols ; j++) {
		for (i = 0; i < rows; i++) {
			set_element(mp1, i, j , i + j * rows);
		}
	}
	print_matrix(mp1);

	copy_matrix(mp1, mp2);
	print_matrix(mp2);

	transpose(mp1, mpt);
	print_matrix(mpt);

	zero_matrix(mp1);
	print_matrix(mp1);

	copy_matrix(mp1, mp2);
	print_matrix(mp2);

	printf("%d %d\n", get_rows(mp1), get_cols(mp1));

	destroy_matrix(mp1);
	destroy_matrix(mp2);
	destroy_matrix(mpt);

	printf("****test_matrix end****\n");
}

void test_vector(void)
{
	int k;
	int dim = 5;
	vector *v1 = create_col_vector(dim);
	vector *v2 = create_col_vector(dim);
	vector *vt = create_row_vector(dim);

	printf("****test_vector start****\n");
	print_vector(v1);
	print_vector(v2);
	print_vector(vt);

	for (k = 0; k < dim; k++) {
		set_vector_element(v1, k, k*k);
	}
	print_vector(v1);

	copy_vector(v1, v2);
	print_vector(v2);

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
