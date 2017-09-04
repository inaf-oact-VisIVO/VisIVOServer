/* $Id: io_parameter.c,v 1.9 2007/11/01 09:23:49 knolli Exp $ */

/**
 * \file io_parameter.c
 *
 * Provides functions for reading in AMIGA parameter files.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "io_parameter.h"
#include "io_util.h"


/**********************************************************************\
 *    Local defines, structure definitions and typedefs               * 
\**********************************************************************/


/**********************************************************************\
 *    Prototypes of local functions                                   * 
\**********************************************************************/
inline static io_parameter_t
local_getparams(FILE *in);


/**********************************************************************\
 *    Implementation of global functions                              * 
\**********************************************************************/
extern io_parameter_t
io_parameter_get(io_parameter_source_t where, char *fname)
{
	FILE *in;
	io_parameter_t dummy;

#	ifdef WITH_MPI
	/* Make sure that we only proceed it we read from a        *\
	 * parameter file, reading from stdin in MPI mode is a bad *
	\* idea.                                                   */
	if (where != IO_PARAMETER_FROM_FNAME)
		return NULL;
#	endif

	/* Check from where to read */
	if (where == IO_PARAMETER_FROM_FNAME) {
		/* Okay, reading from file, open it then */
		if (fname == NULL)
			in = fopen(IO_PARAMETER_FNAME, "r");
		else
			in = fopen(fname, "r");
		/* Did the open work? */
		if (in == NULL) {
			/* Nope, nothing we can do here, return */
			fprintf(stderr,
			        "FATAL: Could not open parameter file (%s)\n",
			        ((fname == NULL) ? IO_PARAMETER_FNAME : fname));
			return NULL;
		}
	} else {
		in = stdin;
	}

	dummy = local_getparams(in);

	/* Make sure to properly close the parameter file */
	if (where == IO_PARAMETER_FROM_FNAME) {
		fclose(in);
	}

	return dummy;
}

extern void
io_parameter_del(io_parameter_t *params)
{
	if (params == NULL)
		return;

	if (*params == NULL);
		return;

	free((*params)->icfile_name);
	free((*params)->outfile_prefix);
#	ifdef LIGHTCONE
	free((*params)->lightcone_prefix);
#	endif
	free((*params)->a_out);

	free(*params);
	
	*params = NULL;

	return;
}

/**********************************************************************\
 *    Implementation of local functions                               * 
\**********************************************************************/
inline static io_parameter_t                                            
local_getparams(FILE *in)
{
	io_parameter_t dummy;
	int i;
	char str[IO_PARAMETER_MAXSTRING];
	char *tmp;

	dummy = (io_parameter_t)malloc(sizeof(io_parameter_struct_t));
	if (dummy == NULL) {
		fprintf(stderr,
		        "Could not allocate memory for parameters.\n");
		return NULL;
	}

#	ifndef WITH_MPI
	printf("please give name of file with initial conditions: ");
#	endif
	io_util_readline(in, str, IO_PARAMETER_MAXSTRING-1);
#	ifndef WITH_MPI
	printf("%s\n", str);
#	endif

	/* See if there is something besides a filename, which will be
	 * interpreted as the file type
	 */
	tmp = strstr(str, " ");
	if (tmp == NULL) {
		dummy->ic_filetype = IO_FILE_AMIGA;
	} else {
		*tmp = '\0';
		do {
			tmp++;
		} while ( (*tmp != '\0') && !isdigit(*tmp) );
		if (sscanf(tmp, "%d", (int *)&(dummy->ic_filetype)) != 1) {
			fprintf(stderr,
			        "Could not gather the file type! "
			        "Assuming AMIGA format");
			dummy->ic_filetype = IO_FILE_AMIGA;
		}
#		ifdef WITH_MPI
		do {
			tmp++;
		} while ( (*tmp != '\0') && isdigit(*tmp) );
		if (sscanf(tmp, "%"SCNu32, &(dummy->reader)) != 1) {
			fprintf(stderr,
			        "Number of reading processes not given. "
			        "Assuming 1.\n");
			dummy->reader = UINT32_C(1);
		}
#		else
		dummy->reader = UINT32_C(1);
#		endif
	} 

	/* Get the filename now */
	dummy->icfile_name = io_util_strdup(str);
	if (dummy->icfile_name == NULL) {
		fprintf(stderr, "Could not store filename!\n");
		free(dummy);
		return NULL;
	}

#	ifndef WITH_MPI
	printf("please give prefix for output file names:         ");
#	endif
	fscanf(in, "%s", str);
#	ifndef WITH_MPI
	printf("%s\n", str);
#	endif
	dummy->outfile_prefix = io_util_strdup(str);
	if (dummy->outfile_prefix == NULL) {
		fprintf(stderr, "Could not store filename!\n");
		free(dummy->icfile_name);
		free(dummy);
		return NULL;
	}

#	ifdef LIGHTCONE
#		ifndef WITH_MPI
	printf("please give prefix for lightcone file names:      ");
#		endif
	fscanf("%s", str);
#		ifndef WITH_MPI
	printf("%s\n", str);
#		endif
	dummy->lightcone_prefix = io_util_strdup(str);
	if (dummy->lightcone_prefix == NULL) {
		fprintf(stderr, "Could not store filename!\n");
		free(dummy->outfile_prefix);
		free(dummy->icfile_name);
		free(dummy);
		return NULL;
	}

#	endif /* LIGHTCONE */

#	ifndef WITH_MPI
	printf("please give number of domain grid cells (1D):     ");
#	endif
	fscanf(in, "%d", &(dummy->NGRID_DOM));
#	ifndef WITH_MPI
	printf("%d\n", dummy->NGRID_DOM);
#	endif

#	ifndef WITH_MPI
	printf("please give Nth for domain grid:                  ");
#	endif
	fscanf(in, "%lf", &(dummy->Nth_dom));
#	ifndef WITH_MPI
	printf("%g\n", dummy->Nth_dom);
#	endif

#	ifndef WITH_MPI
	printf("please give Nth for refinements:                  ");
#	endif
	fscanf(in, "%lf", &(dummy->Nth_ref));
#	ifndef WITH_MPI
	printf("%g\n", dummy->Nth_ref);
#	endif

#	ifndef WITH_MPI
	printf("please give final redshift:                       ");
#	endif
	fscanf(in, "%lf", &(dummy->final_z));
#	ifndef WITH_MPI
	printf("%g\n", dummy->final_z);
#	endif

#	ifdef LIGHTCONE
#		ifndef WITH_MPI
	printf("please give lightcone redshift limit:             ");
#	endif
	fscanf(in, "%lf", &(dummy->lightcone_z));
#		ifndef WITH_MPI
	printf("%g\n", dummy->lightcone_z);
#	endif

#		ifndef WITH_MPI
	printf("please give lightcone type:                       ");
#	endif
	fscanf(in, "%d", &(dummy->lightcone_type));
#		ifndef WITH_MPI
	printf("%d\n", dummy->lightcone_type);
#	endif

	if(lightcone_type==-1) {
#		ifndef WITH_MPI
		printf("please give patch orientation angles (degrees):   ");
#		endif
		fscanf(in, "%lf %lf %lf",
		       dummy->patchangle,
		       dummy->patchangle+1,
		       dummy->patchangle+2);
#		ifndef WITH_MPI
		printf("%g %g %g\n",
		       dummy->patchangle[0],
		       dummy->patchangle[1],
		       dummy->patchangle[2]);
#		endif

#		ifndef WITH_MPI
		printf("please give patch sizes (degrees):                ");
#		endif
		fscanf(in, "%lf %lf",
		       dummy->patcharea, dummy->patcharea+1);
#		ifndef WITH_MPI
		printf("%g %g\n",
		       dummy->patcharea[0],
		       dummy->patcharea[2]);
#		endif
	}
#	endif /* LIGHTCONE */ 

#	ifndef WITH_MPI
	printf("please give dump step:                            ");
#	endif
	fscanf(in, "%d", &(dummy->out_dumps));
#	ifndef WITH_MPI
	printf("%d\n", dummy->out_dumps);
#	endif

#	ifndef WITH_MPI
	printf("please give total number of output files:         ");
#	endif
	fscanf(in, "%d", &(dummy->no_outputs));
#	ifndef WITH_MPI
	printf("%d\n", dummy->no_outputs);
#	endif
   
	/* allocate the a_out array and get the redshift for which to
	 * write outputs */
	dummy->a_out = (double *)malloc(dummy->no_outputs * sizeof(double));
	if (dummy->a_out == NULL) {
		fprintf(stderr, "Could create outa array!\n");
#		ifdef LIGHTCONE
		free(dummy->lightcone_prefix);
#		endif
		free(dummy->outfile_prefix);
		free(dummy->icfile_name);
		free(dummy);
		return NULL;
	}

#	ifndef ISOLATED
	/* for -DISOLATED the output times are *not* provided by the user! */
	for(i = 0; i < dummy->no_outputs; i++) {
#		ifndef WITH_MPI
		printf("please give redshift for %05d. output file:      ",
		       i+1);
#		endif
		fscanf(in, "%lf", dummy->a_out+i);
#		ifndef WITH_MPI
		printf("%g\n", dummy->a_out[i]);
#		endif
		dummy->a_out[i] = (double)1.0/((double)1.0+dummy->a_out[i]);
	}
#	endif                    
                                                                
#	ifdef MOND 
	/* for -DMOND two additional parameters need to be provided */
#		ifndef WITH_MPI
	printf("please give MOND acceleration g0:                 ");
#		endif
	fscanf(in, "%lf", &(dummy->g0));
#		ifndef WITH_MPI
	printf("%g\n", dummy->g0);
#		endif

#		ifndef WITH_MPI
	printf("please give Hubble parameter H0:                  ");
#		endif
	fscanf(in, "%lf", &(dummy->h0));
#		ifndef WITH_MPI
	printf("%g\n", dummy->h0);
#		endif
#	endif

	/* TODO This is not standard for input files, hence define it  *\
	\*      here for the time being.                               */
	dummy->outfile_filetype = IO_FILE_AMIGA;

	return dummy;
}
