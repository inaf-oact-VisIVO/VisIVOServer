#ifndef STARTRUN_INCLUDED
#define STARTRUN_INCLUDED

#include "param.h"
#include "tdef.h"

#ifndef NEWSTARTRUN
void startrun(char AMIGA_input[MAXSTRING], double *timecounter, double *timestep, int *no_first_timestep);
#else

#ifdef HAVE_STDINT_H
#	include <stdint.h>
#else
#	include <replace_stdint.h>
#endif

/**
 * \brief Sets up the simulation.
 *
 * \param *paramfile          Filename of the parameter file. This
 *                            can be NULL in which case the parameters
 *                            will be asked from the commandline.
 * \param *timecounter        TODO
 * \param *timestep           TODO
 * \param *no_first_timestep  TODO
 *
 * \return Nothing.
 */
extern void
startrun(char *paramfile,
         double *timecounter,
         double *timestep,
         int32_t *no_first_timestep);

/**
 * \brief Small function doing some clean-up work.
 *
 * \return Nothing.
 */
extern void
stoprun(void);

#endif /* NEWSTARTRUN */

#endif /* STARTRUN_INCLUDED */

