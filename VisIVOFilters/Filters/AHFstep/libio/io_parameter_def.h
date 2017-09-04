#ifndef IO_PARAMETER_DEF_H
#define IO_PARAMETER_DEF_H

/* $Id: io_parameter_def.h,v 1.4 2007/11/01 09:23:49 knolli Exp $ */

/**
 * \file io_parameter_def.h
 *
 * Provides functions for reading in AMIGA parameter files.
 */


/***********************************************************************\
 *    Includes                                                         * 
\***********************************************************************/
#include "io_file.h"
#ifdef HAVE_STDINT_H
#	include <stdint.h>
#else
#	include <replace_stdint.h>
#endif


/***********************************************************************\
 *    Global defines, structure definitions and typedefs               * 
\***********************************************************************/

#define IO_PARAMETER_MAXSTRING 1024

/**
 * The header structure itself
 */
struct io_parameter_struct {
	char *icfile_name;
	io_file_type_t ic_filetype;
	uint32_t reader;
	char *outfile_prefix;
	io_file_type_t outfile_filetype;
#	ifdef LIGHTCONE
	char *lightcone_prefix;
#	endif
	int NGRID_DOM;
	double Nth_dom;
	double Nth_ref;
	double final_z;
#	ifdef LIGHTCONE
	double lightcone_z;
	int lightcone_type;
	double patchangle[3];
	double patchsize[2];
#	endif
	int out_dumps;
	int no_outputs;
	double *a_out;
#	ifdef MOND
	double g0;
	double h0;
#	endif
};

/** Convenient typedef */
typedef struct io_parameter_struct io_parameter_struct_t;

/** Convenient typedef */
typedef io_parameter_struct_t *io_parameter_t;


#endif /* IO_PARAMETER_DEF_H */
