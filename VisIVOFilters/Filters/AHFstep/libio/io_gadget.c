/* $Id: io_gadget.c,v 1.31 2008/07/31 08:22:59 knolli Exp $ */

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
#ifdef WITH_MPI
#	include <mpi.h>
#endif

#include "io_gadget.h"
#include "io_gadget_header.h"
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
 * io_gadget_open more readable.
 *
 * \param log   The logging object.
 * \param f     The Gadget file object sofar
 * \param mode  The mode in which to open the file.
 *
 * \returns In case of an error NULL is returned, otherwise f.
 */
inline static io_gadget_t
local_openopen(io_logging_t log, io_gadget_t f, io_file_mode_t mode);

 /**
  * \brief Helper funtion to set the swap-state and write some
  *        messages
  *
  * \param log      The logging object.
  * \param f        The Gadget file object sofar.
  * \param swapped  The swap state.
  *
  * \return Nothing.
  */
inline static void
local_openswapped(io_logging_t log,
                  io_gadget_t f,
                  io_file_swap_t swapped);

/**
 * \brief Try to find out swapping status.
 *
 * \param log  The logging object.
 * \param f    The file object.
  *
  * \return Returns the file object or NULL in case of an error.
 */
inline static io_gadget_t
local_opengetswap(io_logging_t log, io_gadget_t f);

/**
 * \brief Tries to figure out which Gadget file version is to be used.
 *
 * \param log  The logging object.
 * \param f    The file object.
 *
 * \return Nothing.
 */
inline static void
local_openversion(io_logging_t log, io_gadget_t f);

/**
 * \brief Check before writing particles and place file pointer at the
 *        right spot.
 *
 * \param log
 * \param f
 * \param pskip
 * \param bytes
 *
 * \return 
 */
inline static int32_t
local_write_common(io_logging_t log,
                   io_gadget_t f,
                   uint64_t pskip,
                   int32_t bytes);


/**********************************************************************\
 *    Implementation of global functions                              * 
\**********************************************************************/
extern io_gadget_t
io_gadget_open(io_logging_t log,
               char *fname,
               io_file_swap_t swapped,
               io_file_mode_t mode,
               uint32_t reader)
{
	io_gadget_t f;

	/* Get memory for the structure */
	f = (io_gadget_t)malloc(sizeof(io_gadget_struct_t));
	if (f == NULL) {
		io_logging_memfatal(log,  "io_gadget structure");
		return NULL;
	}

	/* Start filling the structure */

	/* Store the filename */
	f->fname = (char *)malloc(sizeof(char) * (strlen(fname) + 1));
	if (f->fname == NULL) {
		io_logging_memfatal(log, "filename of GadgetFile");
		free(f);
		return NULL;
	}
	strncpy(f->fname, fname, strlen(fname)+1);

	/* Okay, we are a Gadget file */
	f->ftype = IO_FILE_GADGET;

	/* And we can just copy in the parallel information */
#	ifdef WITH_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &(f->rank));
	MPI_Comm_size(MPI_COMM_WORLD, &(f->size));
	if (f->size >= reader) {
		/* TODO 
		 * THIS IS JUST A QUICK HACK TO PREVENT MGADGET FILES TO
		 * AGAIN TRY TO SPLIT THE COMMUNICATOR, THAT IS ALREADY DONE
		 * IN mgadget.c 
		 * TODO 
		 */
		MPI_Comm_split(MPI_COMM_WORLD, 1,
		               f->rank, &(f->mycomm));
		MPI_Comm_size(f->mycomm, &(f->size_mycomm));
		MPI_Comm_rank(f->mycomm, &(f->rank_mycomm));
	} else {
		f->mycomm = MPI_COMM_NULL;
		f->size_mycomm = -1;
		f->rank_mycomm = -1;
	}
#	endif

	/* Try to open the file and set mode */
	if (local_openopen(log, f, mode) == NULL) {
		free(f->fname);
		free(f);
		return NULL;
	}

	/* Set swapping */
	local_openswapped(log, f, swapped);
	if (    (f->mode == IO_FILE_READ)
	     && (f->swapped == IO_FILE_UNKOWN_SWAPPING)) {
		if (local_opengetswap(log, f) != f) {
			io_logging_fatal(log, "Cannot open this file.");
			free(f->fname);
			free(f);
			return NULL;
		}
	}

	/* Identify Gadget format */
	local_openversion(log, f);

	/* Nothing for the header for now */
	f->header = NULL;

	/* Set some dummy values */
	f->no_part = UINT64_C(0);
	f->no_part_with_mass = UINT64_C(0);
	f->multimass = INT8_C(0);
	f->mmass = 1e40;
	f->minweight = 1e40;
	f->maxweight = 0.0;
	f->sumweight = 0.0;
	f->no_species = INT32_C(0);
	f->posscale = 1.0;
	f->weightscale = 1.0;

	return f;
}

extern void
io_gadget_close(io_logging_t log,
                io_gadget_t *f)
{
	/* Catch NULLs */
	if (f == NULL || *f == NULL)
		return;

	/* Put header to the file if necessary */
	if (    ((*f)->mode == IO_FILE_WRITE)
	     && ((*f)->header != NULL)) {
		io_gadget_header_write(log, (*f)->header, *f);
	}

	/* Close */
	if ((*f)->header != NULL)
		io_gadget_header_del(log, &((*f)->header));
	if ((*f)->fname != NULL)
		free((*f)->fname);
#	ifdef WITH_MPI
	if ((*f)->mycomm != MPI_COMM_NULL)
		MPI_Comm_free(&((*f)->mycomm));
#	endif

	/* Actually close the file */
	if ((*f)->file != NULL)
		fclose((*f)->file);

	/* Cleaning */
	free(*f);
	*f = NULL;

	return;
}

extern void
io_gadget_init(io_logging_t log,
               io_gadget_t f)
{
	if (f == NULL)
		return;

	if (f->header != NULL) {
		io_logging_warn(log, INT32_C(1),
		                "Already have the header information! Rereading.");
		io_gadget_header_del(log, &(f->header));
	}

	if (f->mode != IO_FILE_READ) {
		io_logging_warn(log, INT32_C(1),
		                "%s is not opened for reading. "
		                "Will do nothing.",
		                f->fname);
		return;
	}

	io_logging_msg(log, INT32_C(5),
	               "Starting to initialize file object from %s",
	               f->fname);
	f->header = io_gadget_header_get(log, f);
	io_logging_msg(log, INT32_C(5),
	               "Done with initializing file object from %s",
	               f->fname);

	/* Check for multimass file and also sum up the particle count */
	f->multimass = 0;
	{ 
		int i;
		for (i=0; i<6; i++) {
			f->no_part += (uint64_t)(f->header->np[i]);
			if (    islessequal(f->header->massarr[i], 0.0)
			     && (f->header->np[i]>0)) {
				f->multimass |= (1<<i);
				io_logging_msg(log, INT32_C(5),
				               "Particle type %d requires reading "
				               "of a mass array.", i);
				f->no_part_with_mass += (uint64_t)(f->header->np[i]);
			}
		}
	}

	return;
}

extern uint64_t
io_gadget_readpart(io_logging_t log,
                   io_gadget_t f,
                   uint64_t pskip,
                   uint64_t pread,
                   io_file_strg_struct_t strg)
{
	uint64_t particles_read, tmp;
	double box[3], shift[3];
	double scale_pos, scale_mom, scale_weight;

	/* 
	 * First read the particles unscaled. This will set important
	 * scaling information
	 */
	particles_read = io_gadget_readpart_raw(log, f, pskip, pread, strg);

	if (particles_read != pread) {
		return UINT64_C(0);
	}

	/* And do the scaling */
#ifdef WITH_MPI
	io_gadget_scale_global(log, f->mycomm,  f->maxpos,
	                       f->minpos, &(f->mmass));
#endif
	tmp = io_gadget_scale_particles(log, f->maxpos, f->minpos,
	                                &(f->header->boxsize),
	                                f->header->expansion,
	                                f->posscale, f->mmass,
	                                particles_read, strg);
	if (tmp != particles_read) {
		return tmp;
	}

	/* Wow, we are done! */
	return particles_read;
}

/* No we define a bunch of macros to make life easier */
#define SKIP {io_util_readint32(f->file, &blocksize, f->swapped);}
#define SKIP2 {io_util_readint32(f->file, &blocksize2, f->swapped);}
#define CHECK_BLOCK {\
	if (blocksize != blocksize2) {\
		io_logging_fatal(log,\
		                 "The block boundaries (beginning: %" PRIi32\
		                 " end: %" PRIi32 ") are not identical. "\
		                 "Corrupt file?", blocksize, blocksize2);\
		return UINT64_C(0);\
	} else {\
		io_logging_msg(log, INT32_C(5),\
		               "Block claimed correctly to be %f MB long.", \
		               (float)(blocksize/1024./1024.));\
	}\
}
#define DESCRIBE_BLOCK {\
	if (f->ver == 2) {\
		SKIP;\
		io_util_readstring(f->file, str, (size_t)4);\
		io_util_readint32(f->file, &nextblocksize, f->swapped);\
		io_logging_msg(log, INT32_C(1),\
		               "Arrived at block %s, size of it will be %" \
		               PRIi32, str, nextblocksize);\
		SKIP2;\
		CHECK_BLOCK;\
	}\
}
#define CHECK_FLOATBYTES(bfile, bstore) {\
	if (bfile == sizeof(float)) { \
		io_logging_msg(log, INT32_C(1), \
		               "Obviously the file uses float for " \
		               "floating point values (%" PRIi32 " bytes).", \
		               bfile); \
	} else if (bfile == sizeof(double)) { \
		io_logging_msg(log, INT32_C(1), \
		               "Obviously the file uses double for " \
		               "floating point values (%" PRIi32 " bytes).", \
		               bfile); \
	} else { \
		io_logging_fatal(log, \
		                 "No clue what kind of floating point uses " \
		                 "%" PRIi32 " bytes. Aborting reading.", \
		                 bfile); \
		return UINT64_C(0); \
	}\
	if (bfile < bstore) { \
		io_logging_msg(log, INT32_C(1), \
		               "The floating point values in the file have " \
		               "less precision than the particle storage " \
		               "(%" PRIi32 " bytes vs. %" PRIi32 " bytes). " \
		               "No problem, will upcast.", \
		                bfile, bstore); \
	} else if (bfile > bstore) { \
		io_logging_warn(log, INT32_C(1), \
		                "The floating point values in the file have " \
		                "a higher precision than the particle storage " \
		                "(%" PRIi32 " bytes vs. %" PRIi32 " bytes). " \
		                "Will downcast, but precision might be lost.", \
		                bfile, bstore); \
	} \
}
extern uint64_t
io_gadget_readpart_raw(io_logging_t log,
                       io_gadget_t f,
                       uint64_t pskip,
                       uint64_t pread,
                       io_file_strg_struct_t strg)
{
	uint64_t i, j, psum;
	int ptype;
	uint32_t bytes_file, bytes_int_file;
	int32_t blocksize, blocksize2, partsize;
	int32_t nextblocksize;
	long skipsize;
	double fposx, fposy, fposz;
	double fmomx, fmomy, fmomz;
	double fweight, oldfweight;
	double fu;
	uint64_t fid;
	float dummy;
	uint32_t dummy_int;
	char str[5];
	/** Used to figure out at which particle type we are */
	uint32_t curprtt;


	/* Check if we actually have to do something */
	if ( (f == NULL) || (f->header == NULL) )
		return UINT64_C(0);

	/* Make sure that we are at the beginning of the file */
	rewind(f->file);

	/* Now go to the first block of stuff, skipping the header and in
	 * case of a Gadget-2 file, the descriptive block */
	skipsize = 2*sizeof(int)+GADGET_HEADER_SIZE;
	if (f->ver == 2)
		skipsize += (sizeof(int) * 3 + 4*sizeof(char));
	fseek(f->file, skipsize, SEEK_SET);

	/*******************************************************************\
	 *  Start with particle positions                                  *
	\*******************************************************************/
	DESCRIBE_BLOCK;

	/* Figure out how many bytes are used for float storage */
	SKIP;
	bytes_file = blocksize / (3*f->no_part);
	CHECK_FLOATBYTES(bytes_file, strg.bytes_float);
	io_logging_msg(log, INT32_C(1),
	               "A total of %" PRIu64 " particle positions "
	               "with %" PRIi32 " bytes per float (%f MB total) "
	               "are stored.",
	               f->no_part, bytes_file,
	               (float)(blocksize/1024./1024.));

	/* Set the number of particles to loop over correctly */
	io_logging_msg(log, INT32_C(3),
	               "Asked to read %" PRIu64 " and to skip %" PRIu64
	               " particles. Checking those numbers.",
	               pread, pskip);
	if (pskip > f->no_part) {
		pskip = f->no_part;
		io_logging_msg(log, INT32_C(3),
		               "Cannot skip more than there is, will now "
		               "only skip %" PRIu64 " particles.",
		               pskip);
	}
	if ( f->no_part - pskip < pread ) {
		pread = f->no_part - pskip;
		io_logging_msg(log, INT32_C(3),
		               "Cannot read more than there is left after "
		               "skipping. Will only read %" PRIu64
		               "particles.", pread);
	}
	partsize = 3*bytes_file;

	/* Go to the first particle we want to read */
	fseek(f->file, partsize*pskip, SEEK_CUR);

	/* Set extreme position detectors */
	f->minpos[0] = f->minpos[1] = f->minpos[2] = 1e10;
	f->maxpos[0] = f->maxpos[1] = f->maxpos[2] = -1e10;

	/* Loop over the particle positions */
	for (i=0; i<pread; i++) {
		/* STEP 1:  Read the particle data form the file */
		if (bytes_file == sizeof(float)) {
			io_util_readfloat(f->file, &dummy, f->swapped);
			fposx = (double)dummy;
			io_util_readfloat(f->file, &dummy, f->swapped);
			fposy = (double)dummy;
			io_util_readfloat(f->file, &dummy, f->swapped);
			fposz = (double)dummy;
		} else if (bytes_file == sizeof(double)){
			io_util_readdouble(f->file, &fposx, f->swapped);
			io_util_readdouble(f->file, &fposy, f->swapped);
			io_util_readdouble(f->file, &fposz, f->swapped);
		}
		/* STEP 2:  Detect extreme positions */
		if (isless(fposx, f->minpos[0]))
			f->minpos[0] = fposx;
		if (isless(fposy, f->minpos[1]))
			f->minpos[1] = fposy;
		if (isless(fposz, f->minpos[2]))
			f->minpos[2] = fposz;
		if (isgreater(fposx, f->maxpos[0]))
			f->maxpos[0] = fposx;
		if (isgreater(fposy, f->maxpos[1]))
			f->maxpos[1] = fposy;
		if (isgreater(fposz, f->maxpos[2]))
			f->maxpos[2] = fposz;
		/* STEP 3:  Store the particle in the array */
		if (strg.bytes_float == sizeof(float)) {
			*((float *)strg.posx.val) = (float)fposx;
			*((float *)strg.posy.val) = (float)fposy;
			*((float *)strg.posz.val) = (float)fposz;
		} else {
			*((double *)strg.posx.val) = fposx;
			*((double *)strg.posy.val) = fposy;
			*((double *)strg.posz.val) = fposz;
		}
		/* STEP 4:  Increment the pointers to the next particle */
		strg.posx.val = (void *)(((char *)strg.posx.val)
		                          + strg.posx.stride);
		strg.posy.val = (void *)(((char *)strg.posy.val)
		                          + strg.posy.stride);
		strg.posz.val = (void *)(((char *)strg.posz.val)
		                          + strg.posz.stride);
	} /* End of particle position loop */

	/* Go to the end of the particle position block */
	fseek(f->file, partsize*(f->no_part - (pread + pskip)), SEEK_CUR);
	SKIP2;
	CHECK_BLOCK;


	/*******************************************************************\
	 *  Start with particle velocities                                 *
	\*******************************************************************/
	/* Start with the next block */
	DESCRIBE_BLOCK;

	/* Figure out how many bytes are used for float storage */
	SKIP;
	bytes_file = blocksize / (3*f->no_part);
	CHECK_FLOATBYTES(bytes_file, strg.bytes_float);
	io_logging_msg(log, INT32_C(1),
	               "A total of %" PRIu64 " particle velocities "
	               "with %" PRIi32 " bytes per float (%f MB total) "
	               "are stored.",
	               f->no_part, bytes_file,
	               (float)(blocksize/1024./1024.));

	/* Go to the first particle we want to read */
	fseek(f->file, partsize*pskip, SEEK_CUR);

	/* Loop over the particle velocities */
	for (i=0; i<pread; i++) {
		/* STEP 1:  Read the particle data form the file */
		if (bytes_file == sizeof(float)) {
			io_util_readfloat(f->file, &dummy, f->swapped);
			fmomx = (double)dummy;
			io_util_readfloat(f->file, &dummy, f->swapped);
			fmomy = (double)dummy;
			io_util_readfloat(f->file, &dummy, f->swapped);
			fmomz = (double)dummy;
		} else {
			io_util_readdouble(f->file, &fmomx, f->swapped);
			io_util_readdouble(f->file, &fmomy, f->swapped);
			io_util_readdouble(f->file, &fmomz, f->swapped);
		}
		/* STEP 2:  Store the particle in the array */
		if (strg.bytes_float == sizeof(float)) {
			*((float *)strg.momx.val) = (float)fmomx;
			*((float *)strg.momy.val) = (float)fmomy;
			*((float *)strg.momz.val) = (float)fmomz;
		} else {
			*((double *)strg.momx.val) = fmomx;
			*((double *)strg.momy.val) = fmomy;
			*((double *)strg.momz.val) = fmomz;
		}
		/* STEP 3:  Increment the pointers to the next particle */
		strg.momx.val = (void *)(((char *)strg.momx.val)
		                         + strg.momx.stride);
		strg.momy.val = (void *)(((char *)strg.momy.val)
		                         + strg.momx.stride);
		strg.momz.val = (void *)(((char *)strg.momz.val)
		                         + strg.momx.stride);
	} /* End of particle velocity loop */

	/* Go to the end of the particle velocity block */
	fseek(f->file, partsize*(f->no_part - (pread + pskip)), SEEK_CUR);
	SKIP2;
	CHECK_BLOCK;

	/*******************************************************************\
	 *  Start with particle IDs                                        *
	\*******************************************************************/
	/* Start with the next block */
	DESCRIBE_BLOCK;

	/* Figure out how many bytes are used for the int storage */
	SKIP;
	bytes_int_file = blocksize / f->no_part;
	if (    (bytes_int_file != sizeof(uint32_t))
	     && (bytes_int_file != sizeof(uint64_t))) {
		io_logging_fatal(log,
		                 "Can't handle reading of integers "
		                 "with %" PRIi32 " bytes. Aborting.",
		                 bytes_int_file);
		return i;
	}
	io_logging_msg(log, INT32_C(1),
	               "A total of %" PRIu64 " particle IDs "
	               "with %" PRIi32 " bytes per int (%f MB total) "
	               "are stored.",
	               f->no_part, bytes_int_file,
	               (float)(blocksize/1024./1024.));
	if (bytes_int_file < strg.bytes_int) {
		io_logging_warn(log, INT32_C(1),
		                "File uses %" PRIi32 " bytes per integer for "
		                "the IDs, the particle storage has %" PRIi32
		                " bytes available. Will upcast.",
		                bytes_int_file, strg.bytes_int);
	}
	if (bytes_int_file > strg.bytes_int) {
		io_logging_warn(log, INT32_C(1),
		                "File uses %" PRIi32 " bytes per integer for "
		                "the IDs, the particle storage has only %" PRIi32
		                " bytes available. Will downcast, be aware that "
		                "this might lead to bogus values...",
		                bytes_int_file, strg.bytes_int);
	}

	/* Reset partsize to the bytes user for integer storage value */
	partsize = bytes_int_file;

	/* Go to the first particle we want to read */
	fseek(f->file, partsize*pskip, SEEK_CUR);

	/* See if we have to read the IDs */
	if (strg.id.val == NULL) {
		io_logging_warn(log, INT32_C(1),
	    	            "Discarding IDs as no storage for the IDs "
		                "has been specified.");
		fseek(f->file, partsize*pread, SEEK_CUR);
	} else {
		/* Loop over the particle IDs */
		for (i=0; i<pread; i++) {
			/* STEP 1:  Read the ID from the file */
			if (bytes_int_file == sizeof(uint32_t)) {
				io_util_readuint32(f->file, &dummy_int, f->swapped);
				fid = (uint64_t)dummy_int;
			} else  {
				io_util_readuint64(f->file, &fid, f->swapped);
			}
			/* STEP 2:  Store the ID in the array */
			if (strg.bytes_int == 4) {
				*((uint32_t *)strg.id.val) = (uint32_t)fid;
			} else {
				*((uint64_t *)strg.id.val) = fid;
			}
			/* STEP 3:  Increment the pointers to the next particle */
			strg.id.val = (void *)(((char *)strg.id.val)
			                       + strg.id.stride);
		} /* End of particle ID loop */
	} /* End of catch NULL-Id */

	/* Go to the end of the particle ID block */
	fseek(f->file, partsize*(f->no_part - (pread + pskip)), SEEK_CUR);
	SKIP2;
	CHECK_BLOCK;

	/*******************************************************************\
	 *  Start with particle masses                                     *
	\*******************************************************************/
	/* We are going to need that a few times for book-keeping */
#	define BOOK_KEEPING {\
		if (   isgreater(fweight, oldfweight) \
		    || isless(fweight, oldfweight)) { \
			f->no_species++; \
			oldfweight = fweight; \
			if (isless(fweight,f->minweight)) \
				f->minweight = fweight; \
			if (isgreater(fweight, f->maxweight)) \
				f->maxweight = fweight; \
			if (    (i >= f->header->np[0]) \
			     && (i < f->header->np[0]+f->header->np[1]) \
			     && isless(fweight, f->mmass) ) \
				f->mmass = fweight; \
		} \
	}

	/* First see if there is a mass block at all */
	if (f->multimass != 0) {
		/* Start with the next block */
		DESCRIBE_BLOCK;

		/* Figure out how many bytes are used for float storage */
		SKIP;
		bytes_file = blocksize / (f->no_part_with_mass);
		io_logging_msg(log, INT32_C(1), "blocksize = %i",
		               (int)blocksize);
		CHECK_FLOATBYTES(bytes_file, strg.bytes_float);
		partsize = bytes_file;
		io_logging_msg(log, INT32_C(1),
		               "A total of %" PRIu64 " particle weights "
		               "with %" PRIi32 " bytes per float (%f MB total) "
		               "are stored.",
		               f->no_part_with_mass, bytes_file,
		               (float)(blocksize/1024./1024.));
	
	} else {
		partsize = strg.bytes_float;
		io_logging_msg(log, INT32_C(1),
		               "Will construct mass information solely from"
		               " the header.");
	}

	/* Check if we can store the masses, if not they will be
	 * discarded */
	if (strg.weight.val == NULL) {
		io_logging_warn(log, INT32_C(1),
	    	            "Discarding masses as no storage for the "
		                "masses has been specified.");
	}

	/* Initialize some things */
	f->sumweight = 0.0;
	oldfweight = 0.0;
	f->no_species = 0;

	/* Set the particle type to the first type that actually occurs */
	j = 0;
	curprtt = 0;
	while (f->header->np[curprtt] == 0) {
		curprtt++;
	}

	/* Now we really read the particles */
	for (i=0; i<f->no_part; i++) {
		/* STEP 1:  Update the current particle type */
		if ((long)j >= f->header->np[curprtt]) {
			do {
				curprtt++;
			} while (f->header->np[curprtt] == 0);
			io_logging_msg(log, INT32_C(0),
			               "Detected new particle species! "
			               "i=%" PRIu64 " j=%" PRIu64
			               " curprtt=%" PRIi32 " np[curprtt]=%" PRIu64,
			               i, j, curprtt, f->header->np[curprtt]);
			j = 1;
		} else {
			j++;
		}
		/* STEP 2:  Get the particle weight */
		if (f->multimass & (1<<curprtt)) {
			/* Okay, for this particle type it is in the file */
			if (bytes_file == sizeof(float)) {
				io_util_readfloat(f->file, &dummy, f->swapped);
				fweight = (double)dummy;
			} else {
				io_util_readdouble(f->file, &fweight, f->swapped);
			}
		} else {
			/* For this particle type it is in the header */
			fweight = f->header->massarr[curprtt];
		}
		/* STEP 3:  Do the book-keeping */
		BOOK_KEEPING;
		f->sumweight += fweight;
		/* STEP 4:  Store the particle if it is requested */
		if ( (pskip <= i) && (i-pskip < pread) ) {
			if (strg.weight.val != NULL) {
				/* Okay, there is something to write to */
				if (strg.bytes_float == 4) {
					*((float *)strg.weight.val) = (float)fweight;
				} else {
					*((double *)strg.weight.val) = fweight;
				}
				strg.weight.val = (void *)(((char *)strg.weight.val)
				                           + strg.weight.stride);
			} else {
				/* WE DISCARD IT! */
				;
			}
		}
	} /* End of particle weight loop */

	/* If necessary, verify that the reading went okay */
	if (f->multimass != 0) {
		SKIP2;
		CHECK_BLOCK;
	}
#	undef BOOK_KEEPING

	/*******************************************************************\
	 *  Start with particle energies (gas)                             *
	\*******************************************************************/
	/* See if there is a gas block */
	if (f->header->np[0] > 0) {
		/* Start with the next block */
		DESCRIBE_BLOCK;

		/* Figure out how many bytes are used for float storage */
		SKIP;
		bytes_file = blocksize / (f->header->np[0]);
		CHECK_FLOATBYTES(bytes_file, strg.bytes_float);
		io_logging_msg(log, INT32_C(1),
		               "A total of %" PRIu64 " gas particle energies "
		               "with %" PRIi32 " bytes per float (%f MB total) "
		               "are stored.",
		               f->header->np[0], bytes_file,
		               (float)(blocksize/1024./1024.));
	}

	/* Only go to the gas block if required and if it is actually there */
	if (    f->header->np[0] > 0 && pskip <= (uint64_t)(f->header->np[0])
	     && strg.u.val != NULL) {
		/* Go to the first particle we want to read */
		fseek(f->file, partsize*pskip, SEEK_CUR);

		/* Loop over gas particles */
		for (i=0; i<f->header->np[0]-pskip && i<pread; i++) {
			/* STEP 1:  Read the energy from the file */
			if (bytes_file == sizeof(float)) {
				io_util_readfloat(f->file, &dummy, f->swapped);
				fu = (double)dummy;
			} else  {
				io_util_readdouble(f->file, &fu, f->swapped);
			}
			/* STEP 2:  Store the energy in the array */
			if (strg.bytes_float == sizeof(float)) {
				*((float *)strg.u.val) = (float)fu;
			} else {
				*((double *)strg.u.val) = fu;
			}
			/* STEP 3:  Increment the pointers to the next particle */
			strg.u.val = (void *)(((char *)strg.u.val) + strg.u.stride);
		}

		/* Skip to the end of the energy block */
		fseek(f->file, partsize*(f->header->np[0]-pskip-i), SEEK_CUR);
	} else {
		/* i carries the information of how many energies got read,
		 * since we read none, set it to 0 */
		i = 0;
		/* Skip to the end of the energy block */
		fseek(f->file, partsize*(f->header->np[0]), SEEK_CUR);
	}

	/* If there was an energy block, finish reading it. */
	if (f->header->np[0] > 0) {
		SKIP2;
		CHECK_BLOCK;
	}

	/* Set the energy to the (negative) particle type for the rest*/
	/* Set the sum of all particles already looped over (either skipped,
	 * or actually read) */
	psum = f->header->np[0];
	/* The sum incorporates all particles up to and including this type */
	ptype = 0;
	/* Now we do something ugly and use the energy to store the type of
	 * the particles.  However, we will not do that, if there is no
	 * energy storage provided in the particle structure. */
	if (strg.u.val != NULL) {
		for (; i<pread; i++) {
			/* Figure out if we are still at the right particle type */
			while (i+pskip>=psum && ptype<5) {
				ptype++;
				psum += f->header->np[ptype];
			}
			fu = (double)(-ptype);
			/* STEP 1:  Store the energy in the array */
			if (strg.bytes_float == sizeof(float)) {
				*((float *)strg.u.val) = (float)fu;
			} else {
					*((double *)strg.u.val) = fu;
			}
			/* STEP 2:  Increment the pointers to the next particle */
			strg.u.val = (void *)(((char *)strg.u.val) + strg.u.stride);
		}
	}

	/*******************************************************************\
	 *  Done with reading the gadget file, yeah!                       *
	\*******************************************************************/
	/* Return the number of particles read */
	return pread;
}
/* Getting rid of the macros */
#undef SKIP
#undef SKIP2
#undef CHECK_BLOCK
#undef DESCRIBE_BLOCK
#undef CHECK_FLOATBYTES


extern uint64_t
io_gadget_writepart(io_logging_t log,
                    io_gadget_t f,
                    uint64_t pskip,
                    uint64_t pwrite,
                    io_file_strg_struct_t strg)
{
	return 0;
}

extern uint64_t
io_gadget_writepart_ord(io_logging_t log,
                        io_gadget_t f,
                        uint64_t pskip,
                        uint64_t pwrite,
                        void *nxt_part,
                        io_file_strg_struct_t strg)
{
	return UINT64_C(0);
}

extern bool
io_gadget_get(io_logging_t log,
              io_gadget_t f,
              io_file_get_t what,
              void *res)
{
	if ( (f == NULL) || (f->header == NULL) )
		return false;

	switch (what) {
		case IO_FILE_GET_NOPART_IN_FILE:
		case IO_FILE_GET_NOPART:
			*((long *)res) = (long)f->no_part;
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
			if (f->multimass) {
				if (f->no_species >= 1)
					*((int *)res) = f->no_species;
				else {
					io_logging_warn(log, INT32_C(0),
					                "Don't know the number of particle "
					                "species yet. You first need to "
					                "read the particles.");
					return false;
				}
			} else {
				*((int *)res) = 1;
			}
			break;
		case IO_FILE_GET_BOXSIZE:
			*((double *)res) =   f->header->boxsize
			                   * f->posscale;
			break;
		case IO_FILE_GET_PMASS:
			if (f->multimass) {
				*((double *)res) = f->mmass * f->weightscale;
			} else {
				*((double *)res) =   f->header->massarr[1]
				                   * f->weightscale;
			}
			break;
		case IO_FILE_GET_ZINITIAL:
			io_logging_warn(log, INT32_C(1),
			                "zinitial is not set in a Gadget file, "
			                "using current redshift");
		case IO_FILE_GET_Z:
			*((double *)res) = f->header->redshift;
			break;
		case IO_FILE_GET_AINITIAL:
			io_logging_warn(log, INT32_C(1),
			                "ainitial is not set in a Gadget file, "
			                "using current expansion.");
		case IO_FILE_GET_A:
			*((double *)res) = f->header->expansion;
			break;
		case IO_FILE_GET_OMEGA0:
			*((double *)res) = f->header->omega0;
			break;
		case IO_FILE_GET_OMEGAL:
			*((double *)res) = f->header->omegalambda;
			break;
		case IO_FILE_GET_H:
			*((double *)res) = f->header->hubbleparameter;
			break;
		case IO_FILE_GET_DOUBLE:
			io_logging_warn(log, INT32_C(1),
			                "Gadget files don't store the use of "
			                "double precision. Assuming it is not "
			                "double precision.");
			*((int *)res) = 0;
			break;
		case IO_FILE_GET_MMASS:
			if (isgreater(f->header->massarr[1], 0.0))
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
io_gadget_set(io_logging_t log,
              io_gadget_t f,
              io_file_get_t what,
              void *res)
{
	if ( (f == NULL) || (f->header == NULL) )
		return false;

	switch (what) {
		case IO_FILE_GET_BOXSIZE:
			f->header->boxsize = *((double *)res);
			break;
		case IO_FILE_GET_PMASS:
			f->header->massarr[1] = *((double *)res);
			break;
		case IO_FILE_GET_Z:
			f->header->redshift = *((double *)res);
			break;
		case IO_FILE_GET_A:
			f->header->expansion = *((double *)res);
			break;
		case IO_FILE_GET_OMEGA0:
			f->header->omega0 = *((double *)res);
			break;
		case IO_FILE_GET_OMEGAL:
			f->header->omegalambda = *((double *)res);
			break;
		case IO_FILE_GET_H:
			f->header->hubbleparameter = *((double *)res);
			break;
		default:
			io_logging_fatal(log, "Requesting something unkown in %s.",
			                 __func__);
			return false;
	}

	return true;
}

extern void
io_gadget_log(io_logging_t log, io_gadget_t f)
{
	io_logging_msg(log, INT32_C(5),
	               "Fileobject information:");
	io_logging_msg(log, INT32_C(5),
	               "  Filetype:             %s",
	               io_file_typestr(f->ftype));
	io_logging_msg(log, INT32_C(5),
	               "  Filename:             %s",
	               f->fname);
	io_logging_msg(log, INT32_C(5),
	               "  Mode:                 %" PRIi8,
	               f->mode);
	io_logging_msg(log, INT32_C(5),
	               "  Swapping:             %" PRIi8,
	               f->swapped);
	io_logging_msg(log, INT32_C(5),
	               "  File version:         %" PRIi8,
	               f->ver);
	io_logging_msg(log, INT32_C(5),
	               "  Header size:          %" PRIi8,
	               f->ver);
	io_logging_msg(log, INT32_C(5),
	               "  No. particles:        %" PRIu64,
	               f->no_part);
	io_logging_msg(log, INT32_C(5),
	               "  No. particles w/mass: %" PRIu64,
	               f->no_part_with_mass);
	io_logging_msg(log, INT32_C(5),
	               "  Multimass:            %" PRIi8,
	               f->multimass);
   {
      int i;
      for (i=0; i<6; i++)
         io_logging_msg(log, INT32_C(5),
                        "      Part type %d:      %" PRIi8, i,
                        ((int8_t)(f->multimass) & (1<<i)) >> i);
   }
	io_logging_msg(log, INT32_C(5),
	               "  MMass (Halo parts):   %g",
	               f->mmass);
	io_logging_msg(log, INT32_C(5),
	               "  Minimal Weight:       %g",
	               f->minweight);
	io_logging_msg(log, INT32_C(5),
	               "  Maximal Weight:       %g",
	               f->maxweight);
	io_logging_msg(log, INT32_C(5),
	               "  Sum of all weights:   %g",
	               f->sumweight);
	io_logging_msg(log, INT32_C(5),
	               "  No. of species:       %" PRIi32,
	               f->no_species);
	io_logging_msg(log, INT32_C(5),
	               "  Position scale:       %g",
	               f->posscale);
	io_logging_msg(log, INT32_C(5),
	               "  Weight scale:         %g",
	               f->weightscale);
	io_gadget_header_log(log, f->header);

	return;
}

extern void
io_gadget_resetscale(io_logging_t log,
                     io_gadget_t f,
                     double posscale,
                     double weightscale) {
	if (f == NULL)
		return;

	io_logging_msg(log, INT32_C(8),
	               "Old posscale: %g   New posscale: %g",
	               f->posscale, posscale);
	io_logging_msg(log, INT32_C(8),
	               "Old weightscale: %g   New weightscale: %g",
	               f->weightscale, weightscale);
	f->posscale = posscale;
	f->weightscale = weightscale;

	return;
}

extern uint64_t
io_gadget_scale_particles(io_logging_t log,
                          double maxpos[],
                          double minpos[],
                          double *boxsize,
                          double expansion,
                          double posscale,
                          double mmass,
                          uint64_t particles_read,
                          io_file_strg_struct_t strg)
{
	double box[3], shift[3];
	double scale_pos, scale_mom, scale_weight;
	uint64_t i;

	/* Now we can do the scaling */
	box[0] = fabs(maxpos[0] - minpos[0]);
	box[1] = fabs(maxpos[1] - minpos[1]);
	box[2] = fabs(maxpos[2] - minpos[2]);
	if (isgreater(box[0], *boxsize)) {
		io_logging_warn(log, INT32_C(1),
		                "x-Separation of particles exceeds boxsize "
		                "(%g > %g), resetting boxsize.",
		                box[0], *boxsize);
		*boxsize = box[0];
	}
	if (isgreater(box[1], *boxsize)) {
		io_logging_warn(log, INT32_C(1),
		                "y-Separation of particles exceeds boxsize "
		                "(%g > %g), resetting boxsize.",
		                box[1], *boxsize);
		*boxsize = box[1];
	}
	if (isgreater(box[2], *boxsize)) {
		io_logging_warn(log, INT32_C(1),
		                "z-Separation of particles exceeds boxsize "
		                "(%g > %g), resetting boxsize.",
		                box[2], *boxsize);
		*boxsize = box[2];
	}
	io_logging_msg(log, INT32_C(4),
	               "Extreme positions: xmin = %g  xmax = %g",
	               minpos[0], maxpos[0]);
	io_logging_msg(log, INT32_C(4),
	               "                   ymin = %g  ymax = %g",
	               minpos[1], maxpos[1]);
	io_logging_msg(log, INT32_C(4),
	               "                   zmin = %g  zmax = %g",
	               minpos[2], maxpos[2]);
	shift[0] = (isless(minpos[0], 0.0) ? -minpos[0] : 0.0);
	shift[1] = (isless(minpos[1], 0.0) ? -minpos[1] : 0.0);
	shift[2] = (isless(minpos[2], 0.0) ? -minpos[2] : 0.0);
	io_logging_msg(log, INT32_C(4),
	               "Applying shift: (%g, %g, %g)",
	               shift[0], shift[1], shift[2]);

	/* Set scaling values */
	scale_pos = 1.0/(*boxsize);
	scale_mom =   sqrt(expansion) * expansion
	            / (*boxsize * posscale * 100.);
	scale_weight = 1.0/mmass;
	io_logging_msg(log, INT32_C(3),
	               "Scaling by:  positions:  %g", scale_pos);
	io_logging_msg(log, INT32_C(3),
	               "             velocities: %g", scale_mom);
	io_logging_msg(log, INT32_C(3),
	               "             weights:    %g", scale_weight);

	/* Define the actual scaling calls type independent */
#	define SCALE_CALL(type) {\
		*((type *)strg.posx.val) += (type)(shift[0]); \
		*((type *)strg.posx.val) *= (type)(scale_pos); \
		*((type *)strg.posy.val) += (type)(shift[1]); \
		*((type *)strg.posy.val) *= (type)(scale_pos); \
		*((type *)strg.posz.val) += (type)(shift[2]); \
		*((type *)strg.posz.val) *= (type)(scale_pos); \
		*((type *)strg.momx.val) *= (type)(scale_mom); \
		*((type *)strg.momy.val) *= (type)(scale_mom); \
		*((type *)strg.momz.val) *= (type)(scale_mom); \
		if (strg.weight.val != NULL) \
			*((type *)strg.weight.val) *= (type)(scale_weight); \
		strg.posx.val = (void *)(((char *)strg.posx.val)\
		                         + strg.posx.stride); \
		strg.posy.val = (void *)(((char *)strg.posy.val)\
		                         + strg.posy.stride); \
		strg.posz.val = (void *)(((char *)strg.posz.val)\
		                         + strg.posz.stride); \
		strg.momx.val = (void *)(((char *)strg.momx.val)\
		                         + strg.momx.stride); \
		strg.momy.val = (void *)(((char *)strg.momy.val)\
		                         + strg.momy.stride); \
		strg.momz.val = (void *)(((char *)strg.momz.val)\
		                         + strg.momz.stride); \
		if (strg.weight.val != NULL) \
			strg.weight.val = (void *)(((char *)strg.weight.val)\
			                            + strg.weight.stride); \
	}

	/* Do the scaling depending on the storage type */
	if (strg.bytes_float == sizeof(float)) {
		for (i=0; i<particles_read; i++) {
			SCALE_CALL(float);
		}
	} else if (strg.bytes_float == sizeof(double)) {
		for (i=0; i<particles_read; i++) {
			SCALE_CALL(double);
		}
	} else if (strg.bytes_float == sizeof(long double)) {
		for (i=0; i<particles_read; i++) {
			SCALE_CALL(long double);
		}
	} else {
		io_logging_fatal(log,
		                 "Don't know which floating point types "
		                 "has %" PRIi32 " bytes. Aborting read.",
		                 strg.bytes_float);
		return UINT64_C(0);
	}

	/* Clean the macro away */
#	undef SCALE_CALL

	/* And we are done */
	return particles_read;
}

#ifdef WITH_MPI
extern void
io_gadget_scale_global(io_logging_t log,
                       MPI_Comm comm,
                       double *maxpos,
                       double *minpos,
                       double *mmass)
{
	int size, rank;
	double buffer[3];

	io_logging_msg(log, INT32_C(5),
	               "Updating local scale values to global values.");
	MPI_Allreduce((void *)maxpos, (void *)buffer, 3,
	              MPI_DOUBLE, MPI_MAX, comm);
	io_logging_msg(log, INT32_C(5),
	               "local : maxpos[0] = %g \t"
	               "maxpos[1] = %g \t"
	               "maxpos[2] = %g",
	               maxpos[0], maxpos[1], maxpos[2]);
	maxpos[0] = buffer[0];
	maxpos[1] = buffer[1];
	maxpos[2] = buffer[2];
	io_logging_msg(log, INT32_C(5),
	               "global: maxpos[0] = %g \t"
	               "maxpos[1] = %g \t"
	               "maxpos[2] = %g",
	               maxpos[0], maxpos[1], maxpos[2]);

	MPI_Allreduce((void *)minpos, (void *)buffer, 3,
	              MPI_DOUBLE, MPI_MIN, comm);
	io_logging_msg(log, INT32_C(5),
	               "local : minpos[0] = %g \t"
	               "minpos[1] = %g \t"
	               "minpos[2] = %g",
	               minpos[0], minpos[1], minpos[2]);
	minpos[0] = buffer[0];
	minpos[1] = buffer[1];
	minpos[2] = buffer[2];
	io_logging_msg(log, INT32_C(5),
	               "global: minpos[0] = %g \t"
	               "minpos[1] = %g \t"
	               "minpos[2] = %g",
	               minpos[0], minpos[1], minpos[2]);

	MPI_Allreduce((void *)mmass, (void *)buffer, 1,
	              MPI_DOUBLE, MPI_MIN, comm);
	io_logging_msg(log, INT32_C(5), "local : mmass = %g", *mmass);
	*mmass = buffer[0];
	io_logging_msg(log, INT32_C(5), "global: mmass = %g", *mmass);

	return;
}
#endif /* WITH_MPI */


/**********************************************************************\
 *    Implementation of local functions                               * 
\**********************************************************************/

inline static io_gadget_t
local_openopen(io_logging_t log, io_gadget_t f, io_file_mode_t mode)
{
	if (mode == IO_FILE_READ) {
		f->file = fopen(f->fname, IO_FILE_MODE_READ);
		if (f->file == NULL) {
			io_logging_fatal(log,
			                 "Could not open '%s' for reading.",
			                 f->fname);
			return NULL;
		}
	} else {
		f->file = fopen(f->fname, IO_FILE_MODE_WRITE);
		if (f->file == NULL) {
			io_logging_fatal(log,
			                 "Could not open '%s' for writing.",
			                 f->fname);
			return NULL;
		}
	}

	f->mode = mode;

	return f;
}

inline static void
local_openswapped(io_logging_t log,
                  io_gadget_t f,
                  io_file_swap_t swapped)
{
	if (f->mode == IO_FILE_READ) {
		switch(swapped) {
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
				               "Will try to find out swap status");
		}
	}

	/* Now set the swapping */
	f->swapped = swapped;

	return;
}

inline static io_gadget_t
local_opengetswap(io_logging_t log, io_gadget_t f)
{
	int32_t bbound1, bbound2, tmp;

	/* Get the first block size from the start of the file */
	rewind(f->file);
	fread((void *)(&tmp), sizeof(int32_t), 1, f->file);

	/* Do a byteswap on a copy of that block size */
	bbound1 = tmp;
	io_util_sexchange(&bbound1, sizeof(int32_t));

	io_logging_msg(log, INT32_C(0),
	               "Boundary (file): %" PRIi32
	               " Boundary (swapped): %" PRIi32,
	               tmp, bbound1);

	/* If the bytewapped block size is smaller, then it is probable that
	 * the file is byteswapped.
	 */
	if (bbound1 > tmp) {
		bbound1 = tmp;
		f->swapped = IO_FILE_ISNOT_SWAPPED;
		io_logging_msg(log, INT32_C(0),
		               "Trying nonswapped...");
	} else {
		f->swapped = IO_FILE_IS_SWAPPED;
		io_logging_msg(log, INT32_C(0),
		               "Trying swapped...");
	}

	io_logging_msg(log, INT32_C(0),
	               "Will skip %" PRIi32 " bytes...", bbound1);

	/* Test the assumption by skipping the block and reading the end
	 * block delimiter */
	fseek(f->file, (long)bbound1, SEEK_CUR);
	if (io_util_readint32(f->file, &bbound2, f->swapped) == 0) {
		io_logging_fatal(log,
		                 "Could not read the second block delimiter. "
		                 "Corrupt file?");
		return NULL;
	}

	io_logging_msg(log, INT32_C(0),
	               "Second boundary (file): %" PRIi32,
	               bbound2);

	if (bbound1 == bbound2) {
		if (f->swapped == IO_FILE_IS_SWAPPED) {
			io_logging_msg(log, INT32_C(2),
			               "Identified a byte swapped file.");
		} else {
			io_logging_msg(log, INT32_C(2),
			               "Identified a not byte swapped file.");
		}
	} else {
		/* Swap assumption failed, trying it the other way around */
		io_util_sexchange(&bbound1, sizeof(int32_t));
		fseek(f->file, (long)(sizeof(int32_t)+bbound1), SEEK_SET);

		if (f->swapped == IO_FILE_IS_SWAPPED)
			f->swapped = IO_FILE_ISNOT_SWAPPED;
		else
			f->swapped = IO_FILE_IS_SWAPPED;

		if (io_util_readint32(f->file, &bbound2, f->swapped) == 0) {
			io_logging_fatal(log,
			                 "Could not read the second block delimiter. "
			                 "Corrupt file?");
			return NULL;
		}

		/* See if it worked now */
		if (bbound1 == bbound2) {
			if (f->swapped == IO_FILE_IS_SWAPPED) {
				io_logging_msg(log, INT32_C(2),
				               "Identified a byte swapped file (2. try).");
			} else {
				io_logging_msg(log, INT32_C(2),
				               "Identified a not byte swapped file "
				               "(2. try).");
			}
		}
		else {
			/* Bad, cannot read this file */
			io_logging_fatal(log, "Cannot identify swapping status!");
			return NULL;
		}
	}

	return f;
}

inline static void
local_openversion(io_logging_t log, io_gadget_t f)
{
	int bbound;
	char tmp[5];

	/* We always want to write Gadget2 files */
	if (f->mode == IO_FILE_WRITE) {
		f->ver = 2;
		return;
	}

	/* Read the first block boundary and the first 4 bytes of the block */
	fseek(f->file, 0L, SEEK_SET);
	io_util_readint32(f->file, &bbound, f->swapped);
	io_util_readstring(f->file, tmp, 4);

	/* If it is a Gadget 2 file, only HEAD and an integer will be
	 * stored in the block.
	 */
	if ( (bbound == 4+sizeof(int)) && (strncmp("HEAD", tmp, 4) == 0) ) {
		f->ver = 2;
		io_logging_msg(log, INT32_C(2), "Found a Gadget2 file.");
	} else {
		f->ver = 1;
		io_logging_msg(log, INT32_C(2), "Assuming a Gadget1 file.");
	}

	return;
}

inline static int32_t
local_write_common(io_logging_t log,
                   io_gadget_t f,
                   uint64_t pskip,
                   int32_t bytes)
{
	return 0;
}
