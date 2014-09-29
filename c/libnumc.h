#ifndef _LIB_NUM_C_H_
#define _LIB_NUM_C_H_

#include <stddef.h>

typedef struct vector
{
	void *elt;
	size_t type;
	const unsigned *dim;
} vector;

#endif  /* _LIB_NUM_C_H_ */
