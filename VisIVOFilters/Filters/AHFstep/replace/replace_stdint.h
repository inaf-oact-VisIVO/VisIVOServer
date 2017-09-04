#ifndef STDINT_REPLACE_H
#define STDINT_REPLACE_H

#include <sys/types.h>

#ifndef INT32_C
#	define INT32_C(a) (a)
#endif

#ifndef INT8_C
#	define INT8_C(a) (a)
#endif

#ifndef PRIu32
#	define PRIu32 "u"
#endif
#ifndef PRIi32
#	define PRIi32 "i"
#endif
#ifndef PRIi8
#	define PRIi8 "i"
#endif

#endif /* STDINT_REPLACE_H */
