/**
 * \file replace_stdbool.h 
 *
 * This file takes care of not sufficiently C99-compliant compilers and
 * environments concerning stdbool.h
 */

#ifndef REPLACE_STDBOOL_H
#define REPLACE_STDBOOL_H

/* Assuming we don't have stdbool.h at all */
#ifndef bool
#	define bool    short
#endif

#ifndef true
#	define true    1
#endif

#ifndef false
#	define false   0
#endif

#endif /* STDBOOL_REPLACE_H */
