#include <stdio.h>
#include "test.h"

int main(int argc, char *argv[])
{
	int i, j;
	int rows = 100;
	int cols = 200;
	matrix *mp = create_matrix(rows, cols);

	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols ; j++) {
			set_element(mp, i, j, i * cols + j);
		}
	}

	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols ; j++) {
			printf("%f ", get_element(mp, i, j));
		}
		printf("\n");
	}

	destroy_matrix(mp);
	return 0;
}

