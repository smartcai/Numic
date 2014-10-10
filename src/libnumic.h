#ifndef __NUM_LIBNUMIC_H__
#define __NUM_LIBNUMIC_H__

#include <stddef.h>
#include <stdlib.h>

typedef struct Matrix_Object MObject;
struct Matrix_Object {
	struct Matrix *m;
};

MObject *Matrix_Object_New(int rows, int cols);
void Matrix_Object_Del(MObject *mobj);

#endif  /* __NUM_LIBNUMIC_H__ */
