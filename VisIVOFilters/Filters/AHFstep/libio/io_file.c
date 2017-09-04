/* $Id: io_file.c,v 1.17 2007/12/10 13:13:08 knolli Exp $ */

/**
 * \file io_file.c
 *
 * Provides general definitions and function for reading and writing
 * files.
 */

/***********************************************************************\
 *    Includes                                                         * 
\***********************************************************************/
#ifdef WITH_MPI
#	include <mpi.h>
#endif
#include "io_file.h"
#include "io_amiga.h"
#include "io_ares.h"
#include "io_gadget.h"
#include "io_mgadget.h"
#include "io_deva.h"
#include "io_tipsy.h"


/***********************************************************************\
 *    Definitions of local functions                                   * 
\***********************************************************************/
#ifdef WITH_MPI
inline static void
local_recalcreadskip(uint32_t reader,
                     uint64_t pskip,
                     uint64_t pread,
                     uint64_t *pskip_parallel,
                     uint64_t *pread_parallel);
#endif


/***********************************************************************\
 *    Implemenation of global functions                                * 
\***********************************************************************/
extern const char*
io_file_typestr(io_file_type_t type)
{
	switch (type) {
		case IO_FILE_AMIGA:
			return IO_FILE_AMIGA_STR;
		case IO_FILE_PAMIGA:
			return IO_FILE_PAMIGA_STR;
		case IO_FILE_ARES:
			return IO_FILE_ARES_STR;
		case IO_FILE_MLAPM:
			return IO_FILE_MLAPM_STR;
		case IO_FILE_ASCII:
			return IO_FILE_ASCII_STR;
		case IO_FILE_GADGET:
			return IO_FILE_GADGET_STR;
		case IO_FILE_MGADGET:
			return IO_FILE_MGADGET_STR;
		case IO_FILE_ART:
			return IO_FILE_ART_STR;
      case IO_FILE_DEVA:
         return IO_FILE_DEVA_STR;
      case IO_FILE_TIPSY:
       return IO_FILE_TIPSY_STR;
		case IO_FILE_UNKOWN:
			return IO_FILE_UNKOWN_STR;
		case IO_FILE_EMPTY:
			return IO_FILE_EMPTY_STR;
		default:
			;
			/* Not used, fall through and return */
	}

	return "NOT_DEFINED";
}

extern io_file_t
io_file_open(io_logging_t log,
             char *fname,
             io_file_type_t type,
             io_file_swap_t swapped,
             io_file_mode_t mode,
             uint32_t reader)
{
	io_file_t dummy;
	int size = 1;
	int rank = 0;

#	ifdef WITH_MPI
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#	endif

	/* Do a sanity check */
	if ((mode == IO_FILE_READ) && (reader > size)) {
		io_logging_warn(log, INT32_C(1),
		                "Trying to read with %" PRIu32 " processes, "
		                "but there are only %i. "
		                "Using all those now.", reader, size);
		reader = size;
	}

	/* See if we have anything to do */
	if ( (reader>0) && (mode == IO_FILE_READ) && (rank >= reader) ) {
		io_logging_msg(log, INT32_C(2),
		               "I am not involved in reading.");
		dummy = (io_file_t)malloc(sizeof(io_file_struct_t));
		if (dummy == NULL) {
			io_logging_memfatal(log, "file object");
			return NULL;
		}
		dummy->ftype = IO_FILE_EMPTY;
#	ifdef WITH_MPI
		dummy->rank = rank;
		dummy->size = size;
		MPI_Comm_split(MPI_COMM_WORLD, 99,
		               rank, &(dummy->mycomm));
		MPI_Comm_size(dummy->mycomm, &(dummy->size_mycomm));
		MPI_Comm_rank(dummy->mycomm, &(dummy->rank_mycomm));
#	endif
		return dummy;
	} else {
		io_logging_msg(log, INT32_C(2),
		               "Opening %s, a %s file on %" PRIu32
		                " processes.",
	    	           fname, io_file_typestr(type),
		               ((reader == 0) ? 1 : reader));
	}

  switch (type) {
    case IO_FILE_AMIGA:
      dummy = (io_file_t)io_amiga_open(log, fname, swapped, mode,
                                       reader);
      break;
    case IO_FILE_ARES:
      dummy = (io_file_t)io_ares_open(log, fname, swapped, mode,
                                      reader);
      break;
    case IO_FILE_GADGET:
      dummy = (io_file_t)io_gadget_open(log, fname, swapped, mode,
                                        reader);
      break;
    case IO_FILE_MGADGET:
      dummy = (io_file_t)io_mgadget_open(log, fname, swapped, mode,
                                         reader);
      break;
    case IO_FILE_ART:
      io_logging_fatal(log,
                       "File format %s not supported for %s!",
                       io_file_typestr(type), __func__);
      dummy = NULL;
      break;
    case IO_FILE_DEVA:
      dummy = (io_file_t)io_deva_open(log, fname, swapped, mode,
                                      reader);
      break;
    case IO_FILE_TIPSY:
      dummy = (io_file_t)io_tipsy_open(log, fname, swapped, mode,
                                      reader);
      break;
    default:
      io_logging_fatal(log,
                       "File format %s not supported for %s!",
                       io_file_typestr(type), __func__);
      dummy = NULL;
  }

	return dummy;
}

extern void
io_file_close(io_logging_t log,
              io_file_t *f)
{
  if ( (f == NULL) || (*f == NULL) )
    return;
  
  io_logging_msg(log, INT32_C(2),
                 "Closing a %s file.",
                 io_file_typestr((*f)->ftype));
  switch ((*f)->ftype) {
    case IO_FILE_AMIGA:
      io_amiga_close(log, (io_amiga_t *)f);
      break;
    case IO_FILE_ARES:
      io_ares_close(log, (io_ares_t *)f);
      break;
    case IO_FILE_GADGET:
      io_gadget_close(log, (io_gadget_t *)f);
      break;
    case IO_FILE_MGADGET:
      io_mgadget_close(log, (io_mgadget_t *)f);
      break;
    case IO_FILE_ART:
      io_logging_fatal(log,
                       "File format %s not supported for %s!",
                       io_file_typestr((*f)->ftype), __func__);
      break;
    case IO_FILE_DEVA:
      io_deva_close(log, (io_deva_t *)f);
      break;
    case IO_FILE_TIPSY:
      io_tipsy_close(log, (io_tipsy_t *)f);
      break;
    case IO_FILE_EMPTY:
#			ifdef WITH_MPI
      if ((*f)->mycomm != MPI_COMM_NULL)
        MPI_Comm_free(&((*f)->mycomm));
#			endif
      free(*f);
      *f = NULL;
      break;
    default:
      io_logging_fatal(log,
                       "File format %s not supported for %s!",
                       io_file_typestr((*f)->ftype), __func__);
  }
  return;
}

extern void
io_file_init(io_logging_t log,
             io_file_t f)
{
	if ( (f == NULL) || (f->ftype == IO_FILE_EMPTY) )
		return;

	switch (f->ftype) {
		case IO_FILE_AMIGA:
			io_amiga_init(log, (io_amiga_t)f);
			break;
		case IO_FILE_ARES:
			io_ares_init(log, (io_ares_t)f);
			break;
		case IO_FILE_GADGET:
			io_gadget_init(log, (io_gadget_t)f);
			break;
		case IO_FILE_MGADGET:
			io_mgadget_init(log, (io_mgadget_t)f);
			break;
     case IO_FILE_ART:
       io_logging_fatal(log,
                        "File format %s not supported for %s!",
                        io_file_typestr(f->ftype), __func__);
       break;
     case IO_FILE_DEVA:
       io_deva_init(log, (io_deva_t)f);
       break;
     case IO_FILE_TIPSY:
       io_tipsy_init(log, (io_tipsy_t)f);
       break;
     default:
			io_logging_fatal(log,
			                 "File format %s not supported for %s!",
			                 io_file_typestr(f->ftype), __func__);
	}
  
  return;
}

extern uint64_t
io_file_readpart(io_logging_t log,
                 io_file_t f,
                 uint64_t pskip,
                 uint64_t pread,
                 io_file_strg_struct_t strg)
{
	uint64_t tmp;
	uint64_t pskip_parallel, pread_parallel;

	if ( (f == NULL) || (f->ftype == IO_FILE_EMPTY) )
		return UINT64_C(0);

#	ifdef WITH_MPI
	local_recalcreadskip(f->size_mycomm, pskip, pread,
	                     &pskip_parallel, &pread_parallel);
#	else
		pskip_parallel = pskip;
		pread_parallel = pread;
#	endif

	switch (f->ftype) {
		case IO_FILE_AMIGA:
			tmp = io_amiga_readpart(log, (io_amiga_t)f,
			                        pskip_parallel, pread_parallel,
			                        strg);
			break;
		case IO_FILE_ARES:
			tmp = io_ares_readpart(log, (io_ares_t)f,
			                       pskip_parallel, pread_parallel,
			                       strg);
			break;
		case IO_FILE_GADGET:
			tmp = io_gadget_readpart(log, (io_gadget_t)f,
			                         pskip_parallel, pread_parallel,
			                         strg);
			break;
		case IO_FILE_MGADGET:
			tmp = io_mgadget_readpart(log, (io_mgadget_t)f,
			                          pskip_parallel, pread_parallel,
			                          strg);
			break;
     case IO_FILE_DEVA:
       tmp = io_deva_readpart(log, (io_deva_t)f,
                                 pskip_parallel, pread_parallel,
                                 strg);
       break;
     case IO_FILE_TIPSY:
       tmp = io_tipsy_readpart(log, (io_tipsy_t)f,
                              pskip_parallel, pread_parallel,
                              strg);
       break;
     default:
			io_logging_fatal(log,
			                 "File format %s not supported for %s!",
			                 io_file_typestr(f->ftype), __func__);
			tmp = UINT64_C(0);
	}

	return tmp;
}

extern uint64_t
io_file_writepart(io_logging_t log,
                  io_file_t f,
                  uint64_t pskip,
                  uint64_t pread,
                  io_file_strg_struct_t strg)
{
	uint64_t tmp;

	if ( (f == NULL) || (f->ftype == IO_FILE_EMPTY) )
		return UINT64_C(0);

	switch (f->ftype) {
		case IO_FILE_AMIGA:
			tmp = io_amiga_writepart(log, (io_amiga_t)f,
			                         pskip, pread, strg);
			break;
		case IO_FILE_ARES:
			tmp = io_ares_writepart(log, (io_ares_t)f,
			                        pskip, pread, strg);
			break;
		default:
			io_logging_fatal(log,
			                 "File format %s not supported for %s!\n",
			                 io_file_typestr(f->ftype), __func__);
			tmp = UINT64_C(0);
	}

	return tmp;
}

extern uint64_t
io_file_writepart_ord(io_logging_t log,
                      io_file_t f,
                      uint64_t pskip,
                      uint64_t pread,
                      void *nxt_part,
                      io_file_strg_struct_t strg)
{
	uint64_t tmp;

	if ( (f == NULL) || (f->ftype == IO_FILE_EMPTY) )
		return UINT64_C(0);

	switch (f->ftype) {
		case IO_FILE_AMIGA:
			tmp = io_amiga_writepart_ord(log, (io_amiga_t)f,
			                             pskip, pread, nxt_part, strg);
			break;
		case IO_FILE_ARES:
			tmp = io_ares_writepart_ord(log, (io_ares_t)f,
			                            pskip, pread, nxt_part, strg);
			break;
		default:
			io_logging_fatal(log,
			                 "File format %s not supported for %s!",
			                 io_file_typestr(f->ftype), __func__);
			tmp = 0;
	}

	return tmp;
}

extern uint64_t
io_file_get_numpart(io_logging_t log,
                    io_file_t f,
                    uint64_t *pskip,
                    uint64_t *pread)
{
	uint64_t pskip_parallel, pread_parallel;
	uint64_t numpart = UINT64_C(0);
	long tmp = 0L;

	/* Sanity check */
	if ( (f == NULL) || (f->ftype == IO_FILE_EMPTY) )
		return UINT64_C(0);

	/* 
	 * Get the number of particles available in the file, to do yet
	 * another sanity check
	 */
	io_file_get(log, f, IO_FILE_GET_NOPART_IN_FILE, (void *)(&tmp));
	numpart = (uint64_t)tmp;
    
	/* See if we are reading more than available */
	if (*pskip > numpart) {
		io_logging_warn(log, INT32_C(1),
		                "Trying to skip %"PRIu64" particles, but there"
		                "are only %"PRIu64 " particles in the file. "
		                "Adjusting skipping to %"PRIu64,
		                *pskip, numpart, numpart);
		*pskip = numpart;
	}
	if (*pskip + *pread > numpart) {
		io_logging_warn(log, INT32_C(1),
		                "There are %"PRIu64" particles in the file, "
		                "but the choice of pread and pskip would need "
		                "%"PRIu64 ", adjusting to %"PRIu64,
		                numpart, *pskip + *pread, numpart - *pskip);
		*pread = numpart - *pskip;
	}

	/* Now figure out how many particles actually will be read by
	 * that choice of pskip and pread
	 */
#	ifdef WITH_MPI
	local_recalcreadskip(f->size_mycomm, *pskip, *pread,
	                     &pskip_parallel, &pread_parallel);
#	else
	pskip_parallel = *pskip;
	pread_parallel = *pread;
#	endif

	/* And return the number of particles which will be read */
	return pread_parallel;
}

extern bool
io_file_get(io_logging_t log,
            io_file_t f,
            io_file_get_t what,
            void *res)
{
	if ( (f == NULL) || (res == NULL) || (f->ftype == IO_FILE_EMPTY) )
		return false;

	switch (f->ftype) {
		case IO_FILE_AMIGA:
			return io_amiga_get(log, (io_amiga_t)f, what, res);
		case IO_FILE_ARES:
			return io_ares_get(log, (io_ares_t)f, what, res);
		case IO_FILE_GADGET:
			return io_gadget_get(log, (io_gadget_t)f, what, res);
		case IO_FILE_MGADGET:
			return io_mgadget_get(log, (io_mgadget_t)f, what, res);
      case IO_FILE_DEVA:
         return io_deva_get(log, (io_deva_t)f, what, res);
      case IO_FILE_TIPSY:
         return io_tipsy_get(log, (io_tipsy_t)f, what, res);
		default:
			io_logging_fatal(log,
			                 "File format %s not supported for %s!",
			                 io_file_typestr(f->ftype), __func__);
	}

	return false;
}

extern bool
io_file_set(io_logging_t log,
            io_file_t f,
            io_file_get_t what,
            void *res)
{
	if ( (f == NULL) || (res == NULL) || (f->ftype == IO_FILE_EMPTY) )
		return false;

	switch (f->ftype) {
		case IO_FILE_AMIGA:
			return io_amiga_set(log, (io_amiga_t)f, what, res);
		case IO_FILE_ARES:
			return io_ares_set(log, (io_ares_t)f, what, res);
		case IO_FILE_GADGET:
			return io_gadget_set(log, (io_gadget_t)f, what, res);
		case IO_FILE_MGADGET:
			return io_mgadget_set(log, (io_mgadget_t)f, what, res);
		default:
			io_logging_fatal(log,
			                 "File format %s not supported for %s!",
			                 io_file_typestr(f->ftype), __func__);
	}

	return false;
}

extern void
io_file_log(io_logging_t log,
            io_file_t f)
{
	if ( (f == NULL) || (f->ftype == IO_FILE_EMPTY) )
		return;

	switch (f->ftype) {
		case IO_FILE_AMIGA:
			io_amiga_log(log, (io_amiga_t)f);
			break;
		case IO_FILE_ARES:
			io_ares_log(log, (io_ares_t)f);
			break;
		case IO_FILE_GADGET:
			io_gadget_log(log, (io_gadget_t)f);
			break;
		case IO_FILE_MGADGET:
			io_mgadget_log(log, (io_mgadget_t)f);
			break;
     case IO_FILE_ART:
       io_logging_fatal(log,
                        "File format %s not supported for %s!",
                        io_file_typestr(f->ftype), __func__);
       break;
     case IO_FILE_DEVA:
       io_deva_log(log, (io_deva_t)f);
       break;
     case IO_FILE_TIPSY:
       io_tipsy_log(log, (io_tipsy_t)f);
       break;
     default:
			io_logging_fatal(log,
			                 "File format %s not supported for %s!",
			                 io_file_typestr(f->ftype), __func__);
	}

	return;
}


/***********************************************************************\
 *    Implementations of local functions                               * 
\***********************************************************************/
#ifdef WITH_MPI
inline static void
local_recalcreadskip(uint32_t reader,
                     uint64_t pskip,
                     uint64_t pread,
                     uint64_t *pskip_parallel,
                     uint64_t *pread_parallel)
{
	uint64_t partsperreader;
	int size = 1;
	int rank = 0;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	partsperreader = pread/reader;

	*pskip_parallel = pskip + rank*partsperreader;
	*pread_parallel = partsperreader;
	if (rank == reader-1)
		*pread_parallel = pread - (reader-1)*partsperreader;

	return;
}
#endif /* WITH_MPI */
