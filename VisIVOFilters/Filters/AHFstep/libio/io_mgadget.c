/* $Id: io_mgadget.c,v 1.20 2007/12/10 13:13:08 knolli Exp $ */

/**
 * \file io_gadget.c
 *
 * Provides functions for reading and writing Gadget files.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#if (NONC99 == 0)
#	include <math.h>
#else
#	include <replace_math.h>
#endif
#include <stddef.h>

#include "io_mgadget.h"
#include "io_gadget.h"
#include "io_util.h"


/**********************************************************************\
 *    Local defines, structure definitions and typedefs               * 
\**********************************************************************/


/**********************************************************************\
 *    Prototypes of local functions                                   * 
\**********************************************************************/


/**********************************************************************\
 *    Implementation of global functions                              * 
\**********************************************************************/
extern io_mgadget_t
io_mgadget_open(io_logging_t log,
                char *fname,
                io_file_swap_t swapped,
                io_file_mode_t mode,
                uint32_t reader)
{
	int32_t i;
	io_mgadget_t f;
	char **fnames;

	/* XXX THIS IS CURRENTLY ONLY FOR READING! */
	if (mode != IO_FILE_READ)
		return NULL;

	/* Get memory for the structure */
	f = (io_mgadget_t)malloc(sizeof(io_mgadget_struct_t));
	if (f == NULL) {
		io_logging_memfatal(log,  "io_mgadget structure");
		return NULL;
	}

	/* Start filling the structure */

	/* Okay, we are a Multiple Gadget file */
	f->ftype = IO_FILE_MGADGET;

	/* And we can just copy in the parallel information */
#	ifdef WITH_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &(f->rank));
	MPI_Comm_size(MPI_COMM_WORLD, &(f->size));
	MPI_Comm_split(MPI_COMM_WORLD, 1,
	               f->rank, &(f->mycomm));
	MPI_Comm_size(f->mycomm, &(f->size_mycomm));
	MPI_Comm_rank(f->mycomm, &(f->rank_mycomm));
#	endif

	/* Split the filename in path and stem */
	f->path = NULL;
	f->stem = NULL;
	if (  io_util_split_pathfname(fname, &(f->path), &(f->stem))
	     == NULL) {
		io_logging_fatal(log,
		                 "Could not split %s in path and filename.",
		                 fname);
		free(f);
		return NULL;
	}
	io_logging_msg(log, INT32_C(1),
	               "Will look in %s for %s.",
	               f->path, f->stem);

	/* Get the filenames */
	f->numfiles = io_util_findfiles(f->path, f->stem, &fnames);
	if (f->numfiles <= 0) {
		io_logging_fatal(log,
		                 "Could not open anything starting with %s "
		                 "in %s.",
		                 f->stem, f->path);
		free(f->stem);
		free(f->path);
		free(f);
		return NULL;
	}

	/* Glue the files into the MGadget structure */
	f->files = (io_gadget_t *)malloc( sizeof(io_gadget_t)*(f->numfiles));
	if (f->files == NULL) {
		for (i=0; i<f->numfiles; i++)
			free(fnames[i]);
		free(fnames);
		free(f->stem);
		free(f->path);
		free(f);
		return NULL;
	}
	for (i=0; i<f->numfiles; i++) {
#		ifdef WITH_MPI
		/* TODO 
		 * THIS	IS JUST A NASTY HACK TO PREVENT io_gadget.c FROM REDOING
		 * THE MPI-SPLIT. CALLED FROM HERE IT IS ONLY SUPPOSED TO ACT AS
		 * A DUMMY INTERFACE TO A GADGET FILE.
		 * TODO
		 */
		(f->files)[i] = io_gadget_open(log, fnames[i], swapped, mode,
		                               f->size+1);
#		else
		(f->files)[i] = io_gadget_open(log, fnames[i], swapped, mode,
		                               reader);
#		endif
		if ((f->files)[i] == NULL) {
			int32_t j;
			for (j=i; i<f->numfiles; j++)
				free(fnames[j]);
			free(fnames);
			while (i>0) {
				i--;
				io_gadget_close(log, &((f->files)[i]));
			}
			free(f->stem);
			free(f->path);
			free(f);
			return NULL;
		}
		free(fnames[i]);
	}
	free(fnames);

	/* Set initial values */
	f->no_part = UINT64_C(0);
	f->multimass = INT8_C(0);
	f->mmass = 1e40;
	f->minweight = 1e40;
	f->maxweight = 0.0;
	f->sumweight = 0.0;
	f->no_species = INT32_C(0);
	f->minpos[0] = f->minpos[1] = f->minpos[2] = 1e40;
	f->maxpos[0] = f->maxpos[1] = f->maxpos[2] = -1e40;
	f->posscale = 1.0;
	f->weightscale = 1.0;
	

	return f;
}

extern void
io_mgadget_close(io_logging_t log,
                 io_mgadget_t *f)
{
	int32_t i;

	/* Catch NULLs */
	if (f == NULL || *f == NULL)
		return;

	/* Put header to the file if necessary */
	/* XXX Only relevant for writing */

	/* Close */
	for (i=0; i<(*f)->numfiles; i++)
		io_gadget_close(log, &(((*f)->files)[i]));
	free((*f)->files);
	free(*f);
	*f = NULL;

	return;
}

extern void
io_mgadget_init(io_logging_t log,
                io_mgadget_t f)
{
	int32_t i;

	if (f == NULL)
		return;

	if (f->files[0]->mode != IO_FILE_READ) {
		io_logging_warn(log, INT32_C(1),
		                "%s (first file of %" PRIi32 ") is not opened "
		                "for reading. Will do nothing.",
		                f->files[0]->fname, f->numfiles);
		return;
	}

	/* Check for multimass and sum up the particle count */
	f->no_part = UINT64_C(0);
	f->multimass = INT8_C(0);
	for (i=0; i<f->numfiles; i++) {
		io_gadget_init(log, (f->files)[i]);
		f->multimass |= f->files[i]->multimass;
		f->no_part += f->files[i]->no_part;
	}

	return;
}

extern uint64_t
io_mgadget_readpart(io_logging_t log,
                    io_mgadget_t f,
                    uint64_t pskip,
                    uint64_t pread,
                    io_file_strg_struct_t strg)
{
	uint64_t part_read, tmp;
	double box[3], shift[3];
	double scale_pos, scale_mom, scale_weight;

	/*
	 * First read the particles unscaled.
	 */
	part_read = io_mgadget_readpart_raw(log, f, pskip, pread, strg);
	if (part_read != pread) {
		return UINT64_C(0);
	}

	/* And do the scaling */
#ifdef WITH_MPI
	io_gadget_scale_global(log, f->mycomm,  f->maxpos,
	                       f->minpos, &(f->mmass));
#endif
	tmp = io_gadget_scale_particles(log, f->maxpos, f->minpos, 
	                                &(f->files[0]->header->boxsize),
	                                f->files[0]->header->expansion,
	                                f->posscale, f->mmass,
	                                part_read, strg);
	if (tmp != part_read) {
		return tmp;
	}

	/* Wow, we are done! */
	return part_read;
}


extern uint64_t
io_mgadget_readpart_raw(io_logging_t log,
                        io_mgadget_t f,
                        uint64_t pskip,
                        uint64_t pread,
                        io_file_strg_struct_t strg)
{
	long tmp;
	uint64_t partread, pread_file, pread_done, partinfile;
	uint64_t pskip_file, pskip_done;
	bool something_to_read = false;
	int32_t i;

	/* See if there is anything to do */
	if ( (f == NULL) || (f->files == NULL) )
		return UINT64_C(0);

	/* Initialize accounting of skipping and reading */
	pskip_done = pread_done = UINT64_C(0);

	/* Read the particles from the different files */
	for (i=0; i<f->numfiles; i++) {
		/*
		 * First figure out how many particles are in the file (use the
		 * tmporary long variable for that and copy it over to the 64
		 * bit integer afterwards.
		 */
		if (    io_gadget_get(log, f->files[i],
		                      IO_FILE_GET_NOPART, &tmp)
		     != true) {
			io_logging_fatal(log,
			                 "Could not get number of particles from %s",
			                 f->files[i]->fname);
		}
		partinfile = (uint64_t)tmp;

		/* Then do some arithmetic and set the skipping and reading
		 * numbers correctly for the file 
		 */
		pskip_file = (pskip_done<pskip) ? pskip-pskip_done : UINT64_C(0);
		pread_file = (pread_done<pread) ? pread-pread_done : UINT64_C(0);
		if (pskip_file > partinfile) {
			pskip_file = partinfile;
			pread_file = UINT64_C(0);
		} else {
			if (pread_file > partinfile - pskip_file)
				pread_file = partinfile - pskip_file;
		}

		/* Now read the particles */
		partread = io_gadget_readpart_raw(log, f->files[i],
		                                  pskip_file, pread_file,
		                                  strg);
		if (pread_file != partread) {
			io_logging_fatal(log,
			                 "Something went wrong. Wanted to read %"
			                 PRIu64 " particles, but got %" PRIu64
			                 ". Aborting.", pread_file, partread);
			return pread_done + partread;
		}

		/* Update the skipping arithmetic */
		pskip_done += pskip_file;
		pread_done += pread_file;

		/* Move the particle pointers */
		strg.posx.val = (void *)(((char *)strg.posx.val) 
		                         + strg.posx.stride*partread);
		strg.posy.val = (void *)(((char *)strg.posy.val)
		                         + strg.posy.stride*partread);
		strg.posz.val = (void *)(((char *)strg.posz.val)
		                         + strg.posz.stride*partread);
		strg.momx.val = (void *)(((char *)strg.momx.val)
		                         + strg.momx.stride*partread);
		strg.momy.val = (void *)(((char *)strg.momy.val)
		                         + strg.momy.stride*partread);
		strg.momz.val = (void *)(((char *)strg.momz.val)
		                         + strg.momz.stride*partread);
		if (strg.weight.val != NULL)
			strg.weight.val = (void *)(((char *)strg.weight.val) 
			                           + strg.weight.stride*partread);
		if (strg.id.val != NULL)
			strg.id.val = (void *)(((char *)strg.id.val)
			                       + strg.id.stride*partread);
		if (strg.u.val != NULL)
			strg.u.val = (void *)(((char *)strg.u.val)
			                      + strg.u.stride*partread);

		/* Update the important things needed for scaling */
		if (isless(f->files[i]->mmass, f->mmass))
			f->mmass = f->files[i]->mmass;
		if (isless(f->files[i]->minweight, f->minweight))
			f->minweight = f->files[i]->minweight;
		if (isgreater(f->files[i]->maxweight, f->maxweight))
			f->maxweight = f->files[i]->maxweight;
		f->sumweight += f->files[i]->sumweight;
		/* TODO Something indeed needs be done with no_species... */
		if (isgreater(f->files[i]->maxpos[0], f->maxpos[0]))
			f->maxpos[0] = f->files[i]->maxpos[0];
		if (isgreater(f->files[i]->maxpos[1], f->maxpos[1]))
			f->maxpos[1] = f->files[i]->maxpos[1];
		if (isgreater(f->files[i]->maxpos[2], f->maxpos[2]))
			f->maxpos[2] = f->files[i]->maxpos[2];
		if (isless(f->files[i]->minpos[0], f->minpos[0]))
			f->minpos[0] = f->files[i]->minpos[0];
		if (isless(f->files[i]->minpos[1], f->minpos[1]))
			f->minpos[1] = f->files[i]->minpos[1];
		if (isless(f->files[i]->minpos[2], f->minpos[2]))
			f->minpos[2] = f->files[i]->minpos[2];
	} /* End of particle reading loop*/

	return pread_done;
}

extern bool
io_mgadget_get(io_logging_t log,
               io_mgadget_t f,
               io_file_get_t what,
               void *res)
{
	if ( (f == NULL) || (f->files[0]->header == NULL) )
		return false;

	switch (what) {
		case IO_FILE_GET_NOPART_IN_FILE:
		case IO_FILE_GET_NOPART:
			*((long *)res) = (long)(f->no_part);
			break;
		case IO_FILE_GET_NOVPART:
			if (f->no_part > UINT64_C(0)) {
				if (    isgreater(f->sumweight, 0.0)
				     && isgreater(f->mmass, 0.0))
					*((double *)res) = f->sumweight / f->mmass;
				else
					*((double *)res) = (double)f->no_part;
			} else {
					io_logging_warn(log, INT32_C(0),
					                "Cannot calculate novpart yet. "
					                "You first need to read the "
					                "particles.");
					return false;
			}
			break;
		case IO_FILE_GET_NOSPECIES:
			if (isgreater(f->files[0]->header->massarr[1], 0.0))
				*((int *)res) = 1;
			else {
				*((int *)res) = 1;
				io_logging_warn(log, INT32_C(1),
				                "Not implemented yet.");
				return false;
			}
			break;
		case IO_FILE_GET_BOXSIZE:
			*((double *)res) =   f->files[0]->header->boxsize
			                   * f->posscale;
			break;
		case IO_FILE_GET_PMASS:
			*((double *)res) =   f->mmass * f->weightscale;
			break;
		case IO_FILE_GET_ZINITIAL:
			io_logging_warn(log, INT32_C(1),
			                "zinitial is not set in a Gadget file, "
			                "using current redshift");
		case IO_FILE_GET_Z:
			*((double *)res) = f->files[0]->header->redshift;
			break;
		case IO_FILE_GET_AINITIAL:
			io_logging_warn(log, INT32_C(1),
			                "ainitial is not set in a Gadget file, "
			                "using current expansion");
		case IO_FILE_GET_A:
			*((double *)res) = f->files[0]->header->expansion;
			break;
		case IO_FILE_GET_OMEGA0:
			*((double *)res) = f->files[0]->header->omega0;
			break;
		case IO_FILE_GET_OMEGAL:
			*((double *)res) = f->files[0]->header->omegalambda;
			break;
		case IO_FILE_GET_H:
			*((double *)res) = f->files[0]->header->hubbleparameter;
			break;
		case IO_FILE_GET_DOUBLE:
			io_logging_warn(log, INT32_C(1),
			                "Gadget files don't store the use of "
			                "double precision. Assuming it is not "
			                "double precision.");
			*((int *)res) = 0;
			break;
		case IO_FILE_GET_MMASS:
			if (isgreater(f->files[0]->header->massarr[1], 0.0))
				*((int *)res) = 0;
			else
				*((int *)res) = 1;
			break;
		case IO_FILE_GET_NOTSTEP:
			io_logging_warn(log, INT32_C(1),
			                "Gadget files don't store the step number. "
			                "Setting to 0.");
			*((int32_t *)res) = 0;
			break;
		case IO_FILE_GET_TSTEP:
			io_logging_warn(log, INT32_C(1),
			                "Gadget files don't store the timestep. "
			                "Setting to 0.0");
			*((double *)res) = 0.0;
			break;
		case IO_FILE_GET_HEADERSTR:
			io_logging_warn(log, INT32_C(1),
			                "Gadget files don't have a header string. "
			                "Using a dummy one.");
			*((char **)res) = "No header string.";
			break;
		case IO_FILE_GET_MINWEIGHT:
			if (isgreater(f->mmass, 0.0))
				*((double *)res) = f->minweight / f->mmass;
			else {
				io_logging_warn(log, INT32_C(1),
				                "Don't know minweight yet, setting to "
				                "0.0.");
				*((double *)res) = 0.0;
			}
			break;
		case IO_FILE_GET_MAXWEIGHT:
			if (isgreater(f->mmass, 0.0))
				*((double *)res) = f->maxweight / f->mmass; 
			else {
				io_logging_warn(log, INT32_C(1),
				                "Don't know maxweight yet, setting to "
				                "0.0.");
				*((double *)res) = 0.0;
			}
			break;
		default:
			io_logging_fatal(log, "Requesting something unkown in %s.",
			                 __func__);
			return false;
	}

	return true;
}

extern bool
io_mgadget_set(io_logging_t log,
               io_mgadget_t f,
               io_file_get_t what,
               void *res)
{
	int32_t i;

	if (f == NULL)
		return false;
	for (i=0; i<f->numfiles; i++)
		if (f->files[i]->header == NULL)
			return false;

	switch (what) {
		case IO_FILE_GET_BOXSIZE:
			for (i=0; i<f->numfiles; i++) {
				f->files[i]->header->boxsize = *((double *)res);
			}
			break;
		case IO_FILE_GET_PMASS:
			for (i=0; i<f->numfiles; i++) {
				f->files[i]->header->massarr[1] = *((double *)res);
			}
			break;
		case IO_FILE_GET_Z:
			for (i=0; i<f->numfiles; i++) {
				f->files[i]->header->redshift = *((double *)res);
			}
			break;
		case IO_FILE_GET_A:
			for (i=0; i<f->numfiles; i++) {
				f->files[i]->header->expansion = *((double *)res);
			}
			break;
		case IO_FILE_GET_OMEGA0:
			for (i=0; i<f->numfiles; i++) {
				f->files[i]->header->omega0 = *((double *)res);
			}
			break;
		case IO_FILE_GET_OMEGAL:
			for (i=0; i<f->numfiles; i++) {
				f->files[i]->header->omegalambda = *((double *)res);
			}
			break;
		case IO_FILE_GET_H:
			for (i=0; i<f->numfiles; i++) {
				f->files[i]->header->hubbleparameter = *((double *)res);
			}
			break;
		default:
			io_logging_fatal(log, "Requesting something unkown in %s.",
			                 __func__);
			return false;
	}

	return true;
}

extern void
io_mgadget_log(io_logging_t log, io_mgadget_t f)
{
	int32_t i;

	io_logging_msg(log, INT32_C(5),
	               "Fileobject information:");
	io_logging_msg(log, INT32_C(5),
	               "  Filetype:            %s",
	               io_file_typestr(f->ftype));
	io_logging_msg(log, INT32_C(5),
	               "  Path:                %s",
	               f->path);
	io_logging_msg(log, INT32_C(5),
	               "  Stem:                %s",
	               f->stem);
	io_logging_msg(log, INT32_C(5),
	               "  Number of files:     %" PRIi32,
	               f->numfiles);
	io_logging_msg(log, INT32_C(5),
	               "  Number of particles: %" PRIu64,
	               f->no_part);
	io_logging_msg(log, INT32_C(5),
	               "  Multimass         :  %" PRIi8,
	               f->multimass);
	io_logging_msg(log, INT32_C(5),
	               "  Mmass             :  %g",
	               f->mmass);
	io_logging_msg(log, INT32_C(5),
	               "  Minimal Weight    :  %g",
	               f->minweight);
	io_logging_msg(log, INT32_C(5),
	               "  Maximal Weight    :  %g",
	               f->maxweight);
	io_logging_msg(log, INT32_C(5),
	               "  Sum of weights    :  %g",
	               f->sumweight);
	for (i=0; i<f->numfiles; i++) {
		io_logging_msg(log, INT32_C(5),
		               "  ---> File %"PRIi32 ":",
		               i);
		io_gadget_log(log, (f->files)[i]);
	}

	return;
}

extern void
io_mgadget_resetscale(io_logging_t log,
                      io_mgadget_t f,
                      double posscale,
                      double weightscale) {
	int32_t i;

	if (f == NULL)
		return;

	for (i=0; i<f->numfiles; i++) {
		io_gadget_resetscale(log, f->files[i], posscale, weightscale);
	}
	f->posscale = posscale;
	f->weightscale = weightscale;

	return;
}


/**********************************************************************\
 *    Implementation of local functions                               * 
\**********************************************************************/
