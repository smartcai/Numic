#include <stdio.h>
#include "test.h"

int main(int argc, char *argv[])
{
	int i, j;
	MObject *mobj = matrix_object_new(3, 5);
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 5; j++) {
			matrix_object_set(mobj, i, j, i*5+j);
		}
	}

	for (i = 0; i < 3; i++) {
		for (j = 0; j < 5; j++) {
			printf("%f\n", matrix_object_get(mobj, i, j));
		}
	}

	matrix_object_del(mobj);
	return 0;
}

