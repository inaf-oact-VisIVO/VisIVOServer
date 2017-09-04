#ifndef IO_PARAMETER_H
#define IO_PARAMETER_H

/* $Id: io_parameter.h,v 1.3 2006/08/30 15:14:44 knolli Exp $ */

/**
 * \file io_parameter.h
 *
 * Provides functions for reading in AMIGA parameter files.
 */


/***********************************************************************\
 *    Includes                                                         * 
\***********************************************************************/
#include "io_parameter_def.h"


/***********************************************************************\
 *    Global defines, structure definitions and typedefs               * 
\***********************************************************************/

/** The standard name for parameter file */
#define IO_PARAMETER_FNAME "amiga.input"

/** Defines from where to read the parameters */
typedef enum {
	/* Get the parameters from a file */
	IO_PARAMETER_FROM_FNAME = 0,
	/* Get the parameters from stdin */
	IO_PARAMETER_FROM_STDIN = 1
} io_parameter_source_t;


/***********************************************************************\
 *    Prototypes of global functions                                   * 
\***********************************************************************/

/**
 * \brief Gets the parameters.
 *
 * When in MPI mode this function will check if the parameters
 * are going to be read from a file and only then proceed. Keep
 * in  mind that even though you can thus call the function from
 * each MPI task, the parameter files needs to be readable from
 * every process. And you might produce a lot of traffic.
 *
 * \param where   Selects from where to read the file. Only if
 *                reading from file is requested the next
 *                parameter is evaluated.
 * \param *fname  The filename of the parameter file. May be
 *                NULL, in which case a standarized filename is
 *                used.
 *
 * \return Returns a structure holding the information present
 *         in the file. Might be NULL in the case of errors.
 */
extern io_parameter_t
io_parameter_get(io_parameter_source_t where, char *fname);

/**
 * \brief Disposes a parameter object.
 *
 * \param *param  Pointer to the variable holding the parameter
 *                object. This will be set to NULL.
 *
 * \return Nothing.
 */
extern void
io_parameter_del(io_parameter_t *params);

#endif /* IO_PARAMETER_H */
