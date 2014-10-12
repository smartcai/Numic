#include <stdio.h>
#include "test.h"

int main(int argc, char *argv[])
{
	int i, j;
	int rows = 100;
	int cols = 200;
	matrix *mp = create_matrix(rows, cols);

	/* Pretend to use columns-major store */
	for (j = 0; j < cols ; j++) {
		for (i = 0; i < rows; i++) {
			set_element(mp, i, j , i + j * rows);
		}
	}

	for (j = 0; j < cols ; j++) {
		for (i = 0; i < rows; i++) {
			printf("%f ", get_element(mp, i, j));
		}
		printf("\n");
	}

	destroy_matrix(mp);
	return 0;
}
