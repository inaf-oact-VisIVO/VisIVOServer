/* $Id: io_amiga.c,v 1.24 2007/12/10 13:13:08 knolli Exp $ */

/**
 * \file io_amiga.c
 *
 * Provides functions for reading and writing AMIGA files.
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

#include "io_amiga.h"
#include "io_amiga_header.h"
#include "io_util.h"


/**********************************************************************\
 *    Local defines, structure definitions and typedefs               * 
\**********************************************************************/


/**********************************************************************\
 *    Prototypes of local functions                                   * 
\**********************************************************************/

/**
 * \brief Helper function to open the file
 *
 * This function is supposed to get inlined anyway. Makes
 * io_amiga_open more readable.
 *
 * \param log   The logging object.
 * \param f     The AMIGA file object sofar
 * \param mode  The mode in which to open the file.
 *
 * \returns In case of an error NULL is returned, otherwise f.
 */
inline static io_amiga_t
local_openopen(io_logging_t log, io_amiga_t f, io_file_mode_t mode);

/**
 * \brief Helper funtion to set the swap-state and write some
 *        messages
 *
 * \param log      The logging object.
 * \param f        The AMIGA file object sofar.
 * \param swapped  The swap state.
 *
 * \return Nothing.
 */
inline static void
local_openswapped(io_logging_t log,
                  io_amiga_t f,
                  io_file_swap_t swapped);

/**
 * \brief Helper function to figure out the sizeof(long) in the file.
 *
 * Additionally this is used to determine the byteswap state of the file
 * if not overwritten.
 *
 * \param log  The logging object.
 * \param f    The AMIGA file object sofar.
 *
 * \return Nothing.
 */
inline static void
local_openfsol(io_logging_t log, io_amiga_t f);

/**
 * \brief Check before writing particles and place the file pointer at
 *        the right place.
 *
 * \param log      The logging object.
 * \param f        The file object.
 * \param *weight  If NULL no particle weight will be read.
 * \param pskip    Number of particles to skip in the file.
 * \param bytes    Number of bytes per floating point value.
 *
 * \return Returns 1 if everything went fine, 0 otherwise.
 */
inline static int32_t
local_write_common(io_logging_t log,
                   io_amiga_t f,
                   void *weight,
                   uint64_t pskip,
                   int32_t bytes);


/**********************************************************************\
 *    Implementation of global functions                              * 
\**********************************************************************/
extern io_amiga_t
io_amiga_open(io_logging_t log,
              char *fname,
              io_file_swap_t swapped,
              io_file_mode_t mode,
              uint32_t reader)
{
	io_amiga_t f;

	/* Get memory for the structure */
	f = (io_amiga_t)malloc(sizeof(io_amiga_struct_t));
	if (f == NULL) {
		io_logging_memfatal(log,  "io_amiga structure");
		return NULL;
	}

	/* Start filling the structure */

	/* Store the filename */
	f->fname = (char *)malloc(sizeof(char) * (strlen(fname) + 1));
	if (f->fname == NULL) {
		io_logging_memfatal(log, "filename of AMIGAFILE");
		free(f);
		return NULL;
	}
	strncpy(f->fname, fname, strlen(fname)+1);

	/* Okay, we are an AMIGA file */
	f->ftype = IO_FILE_AMIGA;

	/* And we can just copy in the parallel information */
#	ifdef WITH_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &(f->rank));
	MPI_Comm_size(MPI_COMM_WORLD, &(f->size));
	MPI_Comm_split(MPI_COMM_WORLD, 1,
	               f->rank, &(f->mycomm));
	MPI_Comm_size(f->mycomm, &(f->size_mycomm));
	MPI_Comm_rank(f->mycomm, &(f->rank_mycomm));
#	endif

	/* Try to open the file */
	if (local_openopen(log, f, mode) == NULL)
		return NULL;

	/* Set the mode and then print a message when in READ mode about the
	 * swapping*/
	local_openswapped(log, f, swapped);

	/* Set the file_sizeof_long */
	local_openfsol(log, f);

	/* Nothing for the header for now */
	f->header = NULL;

	/* Initialise the rest to safe parameters */
	f->minweight = 1e40;
	f->maxweight = 0.0;

	return f;
}

extern void
io_amiga_close(io_logging_t log,
               io_amiga_t *f)
{
	/* Catch NULLs */
	if (f == NULL || *f == NULL)
		return;

	/* Put the header to the file if necessary */
	if (    ((*f)->mode == IO_FILE_WRITE)
	     && ((*f)->header != NULL)) {
		io_amiga_header_write(log, (*f)->header, *f);
	}

	/* Close */
	if ((*f)->file != NULL)
		fclose((*f)->file);
	
	/* Clean up */
	if ((*f)->header != NULL)
		io_amiga_header_del(log, &((*f)->header));
	if ((*f)->fname != NULL)
		free((*f)->fname);
#	ifdef WITH_MPI
	MPI_Comm_free(&((*f)->mycomm));
#	endif
	free(*f);
	*f = NULL;

	return;
}

extern void
io_amiga_init(io_logging_t log,
              io_amiga_t f)
{
	if (f->mode == IO_FILE_READ) {
		io_logging_msg(log, INT32_C(5),
		               "Starting to initialize file object from %s",
		               f->fname);
		f->header = io_amiga_header_get(log, f);
		io_logging_msg(log, INT32_C(5),
		               "Done with initializing file object from %s",
		               f->fname);
		if (f->header->multi_mass == 1) {
			f->minweight = f->maxweight = 1.0;
		}
	} else {
		io_logging_warn(log, INT32_C(1),
		                "%s is not opened for reading. "
		                "Will do nothing.",
		                f->fname);
	}

	return;
}

/* Will be used to advance to the next particle. */
#define INCR(a,b) \
	(a) = (void *)(((char *)(a)) + (b))

extern uint64_t
io_amiga_readpart(io_logging_t log,
                  io_amiga_t f,
                  uint64_t pskip,
                  uint64_t pread,
                  io_file_strg_struct_t strg)
{
	uint64_t i;
	uint32_t bytes_file;
	uint32_t partsize;
	float dummy;
	double fposx, fposy, fposz;
	double fmomx, fmomy, fmomz;
	double fweight;

	/* Check if we actually need to do something */
	if (    (f == NULL)
		 || ( f->header == NULL)
	     || (f->header->no_part <= 0)
	     || (((uint64_t)(f->header->no_part)) < pskip) )
		return UINT64_C(0);

	/* Do some ground-work */
	if (f->header->double_precision == 0)
		bytes_file = sizeof(float);
	else
		bytes_file = sizeof(double);

	/* Show what is going to happen */
#define show(bf, bs, what) {\
		if (bf > bs) {\
			io_logging_warn(log, INT32_C(3), \
			                "   " what " (read, downcast)");\
		} else if (bf == bs) {\
			io_logging_msg(log, INT32_C(3), "   " what " (read)");\
		} else {\
			io_logging_warn(log, INT32_C(3),\
			                "   " what " (read, upcast)");\
		}\
	}
	io_logging_msg(log, INT32_C(3), "Will do the following things:");
	show(bytes_file, strg.bytes_float, "posx ");
	show(bytes_file, strg.bytes_float, "posy ");
	show(bytes_file, strg.bytes_float, "posz ");
	show(bytes_file, strg.bytes_float, "momx ");
	show(bytes_file, strg.bytes_float, "momy ");
	show(bytes_file, strg.bytes_float, "momz ");
	if ( (f->header->multi_mass == 1)) {
		if (strg.weight.val != NULL) {
			show(bytes_file, strg.bytes_float, "weight");
		} else {
			io_logging_warn(log, INT32_C(3), "   weight (read, ignored)");
		}
	} else {
		if (strg.weight.val != NULL)
			io_logging_warn(log, INT32_C(3), "   weight (set to 1.0)");
	}
	if (strg.id.val != NULL)
		io_logging_warn(log, INT32_C(3),
		                "   id     (set to position in file)");
	if (strg.u.val != NULL)
		io_logging_warn(log, INT32_C(3), "   u      (set to 0.0)");
#undef show

	/* Set the number of particles to loop over correctly */
	if ( (uint64_t)(f->header->no_part) - pskip < pread)
		pread = (uint64_t)(f->header->no_part) - pskip;

	io_logging_msg(log, INT32_C(3),
	               "Starting to read %" PRIu64 " particles, "
	               "skipping %" PRIu64,
	               pread, pskip);

	/* Position the file at the right spot */
	partsize = bytes_file * 6;
	if (f->header->multi_mass == 1)
		partsize += bytes_file;
	/*
	 * From the beginning of the file, the first 4bytes or for the
	 * swapping and sizeof_long detection. Then we have to skip the
	 * header. If the header was written on a 64bit machine (sizeof_long
	 * 8), then the header structure was unaligned in the memory and we
	 * need to take 2x4b of garbage into account.
	 * Finally we have to skip the particles we don't want to read.
	 */
	fseek(f->file,
	       (long)4+(long)(AMIGA_HEADER_SIZE)
	      +(long)(pskip * partsize)
	      +(f->file_sizeof_long == 8 ? 8 : 0),
	      SEEK_SET);

	/* Loop over all particles to be read */
	for (i=0; i<pread; i++) {
		/************************************************\
		 * Step 1: Read the particle data from the file *
		\************************************************/
		if (bytes_file == 4) {
			io_util_readfloat(f->file, &dummy, f->swapped);
			fposx = (double)dummy;
			io_util_readfloat(f->file, &dummy, f->swapped);
			fmomx = (double)dummy;
			io_util_readfloat(f->file, &dummy, f->swapped);
			fposy = (double)dummy;
			io_util_readfloat(f->file, &dummy, f->swapped);
			fmomy = (double)dummy;
			io_util_readfloat(f->file, &dummy, f->swapped);
			fposz = (double)dummy;
			io_util_readfloat(f->file, &dummy, f->swapped);
			fmomz = (double)dummy;
		} else {
			io_util_readdouble(f->file, &fposx, f->swapped);
			io_util_readdouble(f->file, &fmomx, f->swapped);
			io_util_readdouble(f->file, &fposy, f->swapped);
			io_util_readdouble(f->file, &fmomy, f->swapped);
			io_util_readdouble(f->file, &fposz, f->swapped);
			io_util_readdouble(f->file, &fmomz, f->swapped);
		}
		if (f->header->multi_mass == 1) {
			if (bytes_file == 4) {
				io_util_readfloat(f->file, &dummy, f->swapped);
				fweight = (double)dummy;
			} else {
				io_util_readdouble(f->file, &fweight, f->swapped);
			}
			if (isless(fweight, f->minweight))
				f->minweight = fweight;
			else if (isgreater(fweight, f->maxweight))
				f->maxweight = fweight;
		} else {
         f->minweight = fweight;
         f->maxweight = fweight;
			fweight = 1.0;
		}

		/************************************************\
		 * Step 2: Store the particle in the array      *
		\************************************************/
		if (strg.bytes_float == sizeof(float)) {
			*((float *)strg.posx.val) = (float)fmod(fposx + 1.0, 1.0);
			*((float *)strg.posy.val) = (float)fmod(fposy + 1.0, 1.0);
			*((float *)strg.posz.val) = (float)fmod(fposz + 1.0, 1.0);
			*((float *)strg.momx.val) = (float)fmomx;
			*((float *)strg.momy.val) = (float)fmomy;
			*((float *)strg.momz.val) = (float)fmomz;
			if (strg.weight.val != NULL)
				*((float *)strg.weight.val) = (float)fweight;
		} else {
			*((double *)strg.posx.val) = fmod(fposx + 1.0, 1.0);
			*((double *)strg.posy.val) = fmod(fposy + 1.0, 1.0);
			*((double *)strg.posz.val) = fmod(fposz + 1.0, 1.0);
			*((double *)strg.momx.val) = fmomx;
			*((double *)strg.momy.val) = fmomy;
			*((double *)strg.momz.val) = fmomz;
			if (strg.weight.val != NULL)
				*((double *)strg.weight.val) = fweight;
		}

		/********************************************************\
		 * Step 3: Increment the pointers to the next particle  *
		\********************************************************/
		INCR(strg.posx.val, strg.posx.stride);
		INCR(strg.momx.val, strg.momx.stride);
		INCR(strg.posy.val, strg.posy.stride);
		INCR(strg.momy.val, strg.momy.stride);
		INCR(strg.posz.val, strg.posz.stride);
		INCR(strg.momz.val, strg.momz.stride);
		if (strg.weight.val != NULL)
			INCR(strg.weight.val, strg.weight.stride);

	} /* particle loop ends here*/

	/* Loop again if required to set the ID */
	if (strg.id.val != NULL) {
		if (strg.bytes_int == sizeof(uint64_t)) {
			for (i=0; i<pread; i++) {
				*((uint64_t *)strg.id.val) = pskip + i;
				INCR(strg.id.val, strg.id.stride);
			}
		} else {
			for (i=0; i<pread; i++) {
				*((uint32_t *)strg.id.val) = (uint32_t)(pskip + i);
				INCR(strg.id.val, strg.id.stride);
			}
		}
	}
	/* Loop again if required to set the interal energy to 0.0 */
	if (strg.u.val != NULL) {
		if (strg.bytes_int == sizeof(float)) {
			for (i=0; i<pread; i++) {
				*((float *)strg.u.val) = 0.0;
				INCR(strg.u.val, strg.u.stride);
			}
		} else {
			for (i=0; i<pread; i++) {
				*((double *)strg.u.val) = 0.0;
				INCR(strg.u.val, strg.u.stride);
			}
		}
	}

	io_logging_msg(log, INT32_C(3),
	               "Done with reading %" PRIu64 " particles from %s",
	               i, f->fname);

	return i;
}

extern uint64_t
io_amiga_writepart(io_logging_t log,
                   io_amiga_t f,
                   uint64_t pskip,
                   uint64_t pwrite,
                   io_file_strg_struct_t strg)
{
	uint64_t no_part = 0;

	if (local_write_common(log, f, strg.weight.val,
	                       pskip, strg.bytes_float) == 0)
		return UINT64_C(0);

	/* Start putting the particles to the file */
	for (no_part=0; no_part<pwrite; no_part++) {
		/* Write the current particle */
		fwrite(strg.posx.val, strg.bytes_float, 1, f->file);
		fwrite(strg.momx.val, strg.bytes_float, 1, f->file);
		fwrite(strg.posy.val, strg.bytes_float, 1, f->file);
		fwrite(strg.momy.val, strg.bytes_float, 1, f->file);
		fwrite(strg.posz.val, strg.bytes_float, 1, f->file);
		fwrite(strg.momz.val, strg.bytes_float, 1, f->file);
		if (strg.weight.val != NULL)
			fwrite(strg.weight.val, strg.bytes_float, 1, f->file);

		/* Advance to next particle */
		INCR(strg.posx.val, strg.posx.stride);
		INCR(strg.momx.val, strg.momx.stride);
		INCR(strg.posy.val, strg.posy.stride);
		INCR(strg.momy.val, strg.momy.stride);
		INCR(strg.posz.val, strg.posz.stride);
		INCR(strg.momz.val, strg.momz.stride);
		if (strg.weight.val != NULL)
			INCR(strg.weight.val, strg.weight.stride);
	}

	/* Update the header */
	f->header->no_part = pskip + no_part;
	f->header->multi_mass = (strg.weight.val == NULL ? 0 : 1);
	f->header->double_precision = (strg.bytes_float > 4 ? 1 : 0);

	/* Update logfile */
	io_logging_msg(log, INT32_C(3),
	               "Wrote %" PRIu64 " particles to file",
	               no_part);

	/* Done */
	return no_part;
}

extern uint64_t
io_amiga_writepart_ord(io_logging_t log,
                       io_amiga_t f,
                       uint64_t pskip,
                       uint64_t pwrite,
                       void *nxt_part,
                       io_file_strg_struct_t strg)
{
	uint64_t no_part = 0;
	ptrdiff_t stride;

	if (local_write_common(log, f, strg.weight.val,
	                       pskip, strg.bytes_float) == 0)
		return UINT64_C(0);

	/* Start putting the particles to the file */
	do {
		/* Write the current particle */
		fwrite(strg.posx.val, strg.bytes_float, 1, f->file);
		fwrite(strg.momx.val, strg.bytes_float, 1, f->file);
		fwrite(strg.posy.val, strg.bytes_float, 1, f->file);
		fwrite(strg.momy.val, strg.bytes_float, 1, f->file);
		fwrite(strg.posz.val, strg.bytes_float, 1, f->file);
		fwrite(strg.momz.val, strg.bytes_float, 1, f->file);
		if (strg.weight.val != NULL)
			fwrite(strg.weight.val, strg.bytes_float, 1, f->file);
		no_part++;

		/* Advance to the next particle */
		stride = (char *)*((char **)nxt_part) - (char *)nxt_part;
		nxt_part = (void *)((char *)nxt_part + stride);
		INCR(strg.posx.val, stride);
		INCR(strg.momx.val, stride);
		INCR(strg.posy.val, stride);
		INCR(strg.momy.val, stride);
		INCR(strg.posz.val, stride);
		INCR(strg.momz.val, stride);
		if (strg.weight.val != NULL)
			INCR(strg.weight.val, stride);
	} while (    (*((char **)nxt_part) != NULL)
	          && (no_part < pwrite) );

	io_logging_msg(log, INT32_C(3),
	               "Wrote %" PRIu64 " particles to file",
	               no_part);

	/* Update the header */
	f->header->no_part = (long)(pskip + no_part);
	f->header->multi_mass = (strg.weight.val == NULL ? 0 : 1);
	f->header->double_precision = (strg.bytes_float > 4 ? 1 : 0);

	/* Done */
	return no_part;
}

#undef INCR

extern bool
io_amiga_get(io_logging_t log,
             io_amiga_t f,
             io_file_get_t what,
             void *res)
{
	if ( (f == NULL) || (f->header == NULL) )
		return false;

	switch (what) {
		case IO_FILE_GET_NOPART_IN_FILE:
		case IO_FILE_GET_NOPART:
			*((long *)res) = (long)(f->header->no_part);
			break;
		case IO_FILE_GET_NOVPART:
			*((double *)res) = f->header->no_vpart;
			break;
		case IO_FILE_GET_NOSPECIES:
			*((int *)res) = f->header->no_species;
			break;
		case IO_FILE_GET_BOXSIZE:
			*((double *)res) = f->header->boxsize;
			break;
		case IO_FILE_GET_PMASS:
			*((double *)res) = f->header->pmass;
			break;
		case IO_FILE_GET_ZINITIAL:
			*((double *)res) = 1./(f->header->a_initial) - 1.;
			break;
		case IO_FILE_GET_Z:
			*((double *)res) = 1./(f->header->a_current) - 1.;
			break;
		case IO_FILE_GET_AINITIAL:
			*((double *)res) = f->header->a_initial;
			break;
		case IO_FILE_GET_A:
			*((double *)res) = f->header->a_current;
			break;
		case IO_FILE_GET_OMEGA0:
			*((double *)res) = f->header->omega0;
			break;
		case IO_FILE_GET_OMEGAL:
			*((double *)res) = f->header->lambda0;
			break;
		case IO_FILE_GET_H:
			io_logging_fatal(log,
			                 "The Hubble parameter is not available "
			                 "for AMIGA files.");
			return false;
		case IO_FILE_GET_DOUBLE:
			*((int *)res) = f->header->double_precision;
			break;
		case IO_FILE_GET_MMASS:
			*((int *)res) = f->header->multi_mass;
			break;
		case IO_FILE_GET_NOTSTEP:
			*((int32_t *)res) = (int32_t)(f->header->no_timestep);
			break;
		case IO_FILE_GET_TSTEP:
			*((double *)res) = f->header->timestep;
			break;
		case IO_FILE_GET_HEADERSTR:
			*((char **)res) = &(f->header->header[0]);
			break;
		case IO_FILE_GET_MINWEIGHT:
			*((double *)res) = f->minweight;
			break;
		case IO_FILE_GET_MAXWEIGHT:
			*((double *)res) = f->maxweight;
			break;
		default:
			io_logging_fatal(log, "Requesting something unkown in %s.",
			                 __func__);
			return false;
	}

	return true;
}

extern bool
io_amiga_set(io_logging_t log,
             io_amiga_t f,
             io_file_get_t what,
             void *res)
{
	if (f == NULL)
		return false;

	if (f->header == NULL) {
		io_logging_warn(log, INT32_C(3),
		                "File does not have a header yet, "
		                "creating one.");
		f->header = io_amiga_header_new(log);
		if (f->header == NULL)
			return false;
	}

	switch (what) {
		case IO_FILE_GET_BOXSIZE:
			f->header->boxsize = *((double *)res);
			break;
		case IO_FILE_GET_PMASS:
			f->header->pmass = *((double *)res);
			break;
		case IO_FILE_GET_Z:
			f->header->a_current = 1./(*((double *)res) + 1.);
			f->header->a_initial = 0.1/(*((double *)res) + 1.);
			break;
		case IO_FILE_GET_A:
			f->header->a_current = *((double *)res);
			f->header->a_initial = *((double *)res)/10.0;
			break;
		case IO_FILE_GET_OMEGA0:
			f->header->omega0 = *((double *)res);
			break;
		case IO_FILE_GET_OMEGAL:
			f->header->lambda0 = *((double *)res);
			break;
		case IO_FILE_GET_H:
			io_logging_fatal(log,
			                 "The Hubble parameter is not available "
			                 "for AMIGA files.");
			return false;
		default:
			io_logging_fatal(log, "Requesting something unkown in %s.",
			                 __func__);
			return false;
	}

	return true;
}

extern void
io_amiga_log(io_logging_t log, io_amiga_t f)
{
	io_logging_msg(log, INT32_C(5),
	               "Fileobject information:");
	io_logging_msg(log, INT32_C(5),
	               "  Filetype:          %s",
	               io_file_typestr(f->ftype));
	io_logging_msg(log, INT32_C(5),
	               "  Filename:          %s",
	               f->fname);
	io_logging_msg(log, INT32_C(5),
	               "  Mode:              %" PRIi8,
	               f->mode);
	io_logging_msg(log, INT32_C(5),
	               "  Swapping:          %" PRIi8,
	               f->swapped);
	io_logging_msg(log, INT32_C(5),
	               "  File sizeof(long): %" PRIi8,
	               f->file_sizeof_long);
	io_logging_msg(log, INT32_C(5),
	               "  minweight:         %g",
	               f->minweight);
	io_logging_msg(log, INT32_C(5),
	               "  maxweight:         %g",
	               f->maxweight);
	io_amiga_header_log(log, f->header);

	return;
}


/**********************************************************************\
 *    Implementation of local functions                               * 
\**********************************************************************/
inline static io_amiga_t
local_openopen(io_logging_t log, io_amiga_t f, io_file_mode_t mode)
{
	if (mode == IO_FILE_READ) {
		f->file = fopen(f->fname, IO_FILE_MODE_READ);
		if (f->file == NULL) {
			io_logging_fatal(log,
			                 "Could not open '%s' for reading.",
			                 f->fname);
			free(f);
			return NULL;
		}
	} else {
		f->file = fopen(f->fname, IO_FILE_MODE_WRITE);
		if (f->file == NULL) {
			io_logging_fatal(log,
			                 "Could not open '%s' for writing.",
			                 f->fname);
			free(f);
			return NULL;
		}
	}
	
	f->mode = mode;

	return f;
}

inline static void
local_openswapped(io_logging_t log,
                  io_amiga_t f,
                  io_file_swap_t swapped)
{
	if (f->mode == IO_FILE_READ) {
		switch (swapped) {
			case IO_FILE_ISNOT_SWAPPED:
				io_logging_msg(log, INT32_C(3),
				               "Assuming unswapped file");
				break;
			case IO_FILE_IS_SWAPPED:
				io_logging_msg(log, INT32_C(3),
				               "Assuming swapped file");
				break;
			case IO_FILE_UNKOWN_SWAPPING:
			default:
				io_logging_msg(log, INT32_C(3),
				               "Will find out swap status");
		}
	}

	/* Now set the swapping */
	f->swapped = swapped;

	return;
}

inline static void
local_openfsol(io_logging_t log, io_amiga_t f)
{
	/* In write mode, that is easy */
	if (f->mode == IO_FILE_WRITE) {
		f->file_sizeof_long = (int32_t)sizeof(uint64_t);
		return;
	}

	/* We are reading from a file, so checking for the swapping */
	if (f->swapped == IO_FILE_UNKOWN_SWAPPING) {
		/* But we only do that if the swapping is unkown */
		io_util_readint32(f->file, &(f->file_sizeof_long),
		                  IO_FILE_ISNOT_SWAPPED);
		if (f->file_sizeof_long > 8) {
			/* Wupps, byteswapped value found */
			io_logging_msg(log, INT32_C(2),
			               "Byteswapped file, will do the swapping.");
			io_util_sexchange((void *)(&(f->file_sizeof_long)), 4);
			if (f->swapped == IO_FILE_UNKOWN_SWAPPING)
				f->swapped = IO_FILE_IS_SWAPPED;
			/* Verify that the value is now making sense */
			if (f->file_sizeof_long > 8) {
				/* Still not making sense, corrupt file? */
				io_logging_warn(log, INT32_C(0),
				                "The file does not make sense. Maybe "
				                "it is corrupt or not an AMIGA binary "
				                "file?");
				/* Anyhow, continue */
			}
		} else {
			io_logging_msg(log, INT32_C(2),
			               "File is not byteswapped, everything fine.");
			f->swapped = IO_FILE_ISNOT_SWAPPED;
		}
	} else {
		/* Okay, we believe the swapping value and read the int */
		io_util_readint32(f->file, &(f->file_sizeof_long), f->swapped);
		/* But be nice and check for sanity */
		if (f->file_sizeof_long > 8) {
			io_logging_warn(log, INT32_C(0),
			                "Trying to read with the predefined "
			                "byteswapping state, but the value for "
			                "sizeof(long) is weird. Are you sure "
			                "the swapping is correct?");
		}
	}

	return;
}

inline static int32_t
local_write_common(io_logging_t log,
                   io_amiga_t f,
                   void *weight,
                   uint64_t pskip,
                   int32_t bytes)
{
	int32_t partsize;

	/* File sanity checks */
	if (f == NULL) {
		io_logging_warn(log, INT32_C(1),
		                "File object does not exist. Not writing.");
		return INT32_C(0);
	}
	if (f->mode != IO_FILE_WRITE) {
		io_logging_warn(log, INT32_C(1),
		                "File not opened for writing. Not writing.");
		return INT32_C(0);
	}
	if (f->file == NULL) {
		io_logging_warn(log, INT32_C(1),
		                "File claims to be opened, but it is not. "
		                "Not writing.");
		return INT32_C(0);
	}
	if (f->header == NULL) {
		io_logging_warn(log, INT32_C(3),
		                "File does not have a header yet, "
		                "creating one.");
		f->header = io_amiga_header_new(log);
		if (f->header == NULL)
			return INT32_C(0);
	}

	/* See if we have to write NULLs to the header part */
	if (fseek(f->file, 0L, SEEK_END) != 0) {
		io_logging_fatal(log, "Could not seek in the file!");
		abort();
	}
	if ((uint64_t)ftell(f->file) <   sizeof(io_amiga_header_struct_t)
	                                + sizeof(uint32_t)) {
		int8_t nix = 0;
		int32_t i;

		io_logging_msg(log, INT32_C(3),
		               "Need to write a dummy header");
		rewind(f->file);
		fwrite((void *)&(f->file_sizeof_long),
		       sizeof(uint32_t),
		       1,
		       f->file);
		for (i=0; i<AMIGA_HEADER_SIZE; i++)
			fwrite((void *)&nix,
			       sizeof(int8_t),
			       1,
			       f->file);
		io_logging_msg(log, INT32_C(3),
		               "Wrote a dummy header");
	}

	/* Put the file pointer to the right place */
	partsize = bytes * 6;
	if (weight != NULL)
		partsize += bytes;
	fseek(f->file,
	        4L + (long)(sizeof(io_amiga_header_struct_t))
	      + (pskip*partsize),
	      SEEK_SET);

	return INT32_C(1);
}
