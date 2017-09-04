#include "../param.h"
#include "../tdef.h"

#ifndef LLMSORT_INCLUDED
#define LLMSORT_INCLUDED

void *sort_linked_list(void *p, unsigned index,
   int (*compare)(void *, void *, void *), void *pointer,
   unsigned long *pcount);

#endif

