#include "test.h"

int main(int argc, char *argv[])
{
	MObject *m = matrix_object_new(3, 5);
	matrix_object_del(m);
	return 0;
}

