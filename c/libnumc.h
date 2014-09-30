#ifndef _LIB_NUM_C_H_
#define _LIB_NUM_C_H_

#include <stddef.h>

typedef struct mat
{
	void *elt;
	size_t type;
	unsigned dim;
} mat;

#endif  /* _LIB_NUM_C_H_ */
