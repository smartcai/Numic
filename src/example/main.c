#include "test.h"

int main(int argc, char *argv[])
{
	MObject *m = Matrix_Object_New(3, 5);
	Matrix_Object_Del(m);
	return 0;
}

