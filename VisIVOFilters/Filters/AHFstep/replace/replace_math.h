/**
 * \file replace_math.h 
 *
 * This file takes care of not sufficiently C99-compliant compilers and
 * environments concerning uncomplete math.h.
 */

#if (NONC99 == 1)
/* This is for SunOS 2.9, tweaked to work on lomond.epcc.ed.ac.uk */

	/* Try to get as much as possible from math.h anyway */
#	include <math.h>

	/* This'll give us log2. NOTE: You must then also link with -lsunmath! */
#	include <sunmath.h>

	/* round() is still unkown, but check before defining */
#	ifndef round
#		define round(a) \
		   ((a) < 0.0 ? ceil((a) - 0.5) : floor((a) + 0.5))
#	else
		/* Tell the user that round() apparently is defined. */
#		warning No need to define round(), you seem to have that.
#	endif

	/* Comparing floating points */
	/* XXX These defines are inferior to the proper C99 ones. I might do
	 *     something about that
	 */
#	ifndef isless
#		define isless(a,b) \
		        (a<b)
#	endif
#	ifndef isgreater
#		define isgreater(a,b) \
		        (a>b)
#	endif
#	ifndef islessequal
#		define islessequal(a,b) \
		        (a<=b)
#	endif

#endif
