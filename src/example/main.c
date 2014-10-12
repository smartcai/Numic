#include <stdio.h>
#include "test.h"

int main(int argc, char *argv[])
{
	int i, j;
	int rows = 3;
	int cols = 4;
	matrix *mp = create_matrix(rows, cols);
	matrix *tmp = create_matrix(cols, rows);

	/* Pretend to use columns-major store */
	for (j = 0; j < cols ; j++) {
		for (i = 0; i < rows; i++) {
			set_element(mp, i, j , i + j * rows);
		}
	}

	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols ; j++) {
			printf("%f ", get_element(mp, i, j));
		}
		printf("\n");
	}

	transpose(mp, tmp);

	for (i = 0; i < cols; i++) {
		for (j = 0; j < rows ; j++) {
			printf("%f ", get_element(tmp, i, j));
		}
		printf("\n");
	}

	destroy_matrix(mp);
	destroy_matrix(tmp);
	return 0;
}
