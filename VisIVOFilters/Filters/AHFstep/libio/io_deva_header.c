
/**
 * \file io_deva_header.c
 *
 * Provides functions for reading and writing the header of DEVA
 * files.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#ifdef HAVE_STDINT_H
#	include <stdint.h>
#else
#	include <replace_stdint.h>
#endif

#include "io_deva_header.h"
#include "io_util.h"


/**********************************************************************\
 *    Local defines, structure definitions and typedefs               * 
\**********************************************************************/
#define SKIP(f) {fseek(f, 4L, SEEK_CUR);}


/**********************************************************************\
 *    Prototypes of local functions                                   * 
\**********************************************************************/


/**********************************************************************\
 *    Implementation of global functions                              * 
\**********************************************************************/
extern io_deva_header_t
io_deva_header_get(io_logging_t log, io_deva_t f)
{
	io_deva_header_t dummy;
	long skipsize;
  int iloop;

	/* Some sanity checks */
	if ((f == NULL) || (f->file == NULL))
		return NULL;

	/* Check if there already is a header, do nothing then */
	if (f->header != NULL)
		return f->header;

	/* Create the header structure array */
	dummy = (io_deva_header_t)malloc((size_t)DEVA_HEADER_SIZE);
	if (dummy == NULL) {
		io_logging_memfatal(log, "DEVA header structure");
		return NULL;
	}

	/* just to be on the safe side: rewind the file to the start... */
  rewind(f->file);
	
  /* Now start reading the header */
  io_util_readint32(f->file, &(dummy->itime), f->swapped);
  io_util_readint32(f->file, &(dummy->itstop), f->swapped);
  io_util_readint32(f->file, &(dummy->itdump), f->swapped);
  io_util_readint32(f->file, &(dummy->iout), f->swapped);
  io_util_readint32(f->file, &(dummy->nsformed), f->swapped);
  io_util_readint32(f->file, &(dummy->nsdead), f->swapped);
  io_util_readint32(f->file, &(dummy->irun), f->swapped);
  io_util_readint32(f->file, &(dummy->nobj), f->swapped);
  io_util_readint32(f->file, &(dummy->ngas), f->swapped);
  io_util_readint32(f->file, &(dummy->ndark), f->swapped);
  io_util_readint32(f->file, &(dummy->L), f->swapped);
  io_util_readfloat(f->file, &(dummy->CHEMEVOL), f->swapped);
  io_util_readfloat(f->file, &(dummy->ANOTHER), f->swapped);
  io_util_readfloat(f->file, &(dummy->COOL), f->swapped);
  io_util_readfloat(f->file, &(dummy->REFINEMENT), f->swapped);
  io_util_readfloat(f->file, &(dummy->HYDRO), f->swapped);
  io_util_readfloat(f->file, &(dummy->GRAVITY), f->swapped);
  io_util_readint32(f->file, &(dummy->ISOLATED), f->swapped);
  io_util_readfloat(f->file, &(dummy->EXPAND), f->swapped);
  io_util_readfloat(f->file, &(dummy->COMOVING), f->swapped);
  io_util_readfloat(f->file, &(dummy->STARFORM), f->swapped);
  io_util_readfloat(f->file, &(dummy->GRADH), f->swapped);
  io_util_readint32(f->file, &(dummy->INITIALCOND), f->swapped);
  io_util_readint32(f->file, &(dummy->nstar), f->swapped);
  io_util_readint32(f->file, &(dummy->iseed1), f->swapped);
  io_util_readint32(f->file, &(dummy->ispec), f->swapped);
  io_util_readint32(f->file, &(dummy->indxsp), f->swapped);
  io_util_readint32(f->file, &(dummy->n_neigh), f->swapped);
  io_util_readint32(f->file, &(dummy->lastbar), f->swapped);
  for(iloop=0; iloop<100-29; iloop++)
    io_util_readint32(f->file, &(dummy->fill1[iloop]), f->swapped);

  io_util_readfloat(f->file, &(dummy->time), f->swapped);
  io_util_readfloat(f->file, &(dummy->atime), f->swapped);
  io_util_readfloat(f->file, &(dummy->htime), f->swapped);
  io_util_readfloat(f->file, &(dummy->dtime), f->swapped);
  io_util_readfloat(f->file, &(dummy->E_init), f->swapped);
  io_util_readfloat(f->file, &(dummy->E_kin), f->swapped);
  io_util_readfloat(f->file, &(dummy->E_ther), f->swapped);
  io_util_readfloat(f->file, &(dummy->E_pot), f->swapped);
  io_util_readfloat(f->file, &(dummy->Radiation), f->swapped);
  io_util_readfloat(f->file, &(dummy->Esum), f->swapped);
  io_util_readfloat(f->file, &(dummy->Rsum), f->swapped);
  io_util_readfloat(f->file, &(dummy->cpu), f->swapped);
  io_util_readfloat(f->file, &(dummy->time_end), f->swapped);
  io_util_readfloat(f->file, &(dummy->tout), f->swapped);
  io_util_readfloat(f->file, &(dummy->padding), f->swapped);
  io_util_readfloat(f->file, &(dummy->Tlost), f->swapped);
  io_util_readfloat(f->file, &(dummy->Qlost), f->swapped);
  io_util_readfloat(f->file, &(dummy->Ulost), f->swapped);
  io_util_readfloat(f->file, &(dummy->delta_min), f->swapped);
  io_util_readfloat(f->file, &(dummy->delta_max), f->swapped);
  io_util_readfloat(f->file, &(dummy->T_min), f->swapped);
  io_util_readfloat(f->file, &(dummy->avisc), f->swapped);
  io_util_readfloat(f->file, &(dummy->bvisc), f->swapped);
  io_util_readfloat(f->file, &(dummy->eta2), f->swapped);
  io_util_readfloat(f->file, &(dummy->rho_star), f->swapped);
  io_util_readfloat(f->file, &(dummy->c_star), f->swapped);
  io_util_readfloat(f->file, &(dummy->rmtot), f->swapped);
  io_util_readfloat(f->file, &(dummy->rmsep), f->swapped);
  io_util_readfloat(f->file, &(dummy->dnthres), f->swapped);
  io_util_readfloat(f->file, &(dummy->sft0), f->swapped);
  io_util_readfloat(f->file, &(dummy->sftmin), f->swapped);
  io_util_readfloat(f->file, &(dummy->sftmax), f->swapped);
  io_util_readfloat(f->file, &(dummy->h100), f->swapped);
  io_util_readfloat(f->file, &(dummy->box100), f->swapped);
  io_util_readfloat(f->file, &(dummy->rmgas), f->swapped);
  io_util_readfloat(f->file, &(dummy->rmdark), f->swapped);
  io_util_readfloat(f->file, &(dummy->omega0), f->swapped);
  io_util_readfloat(f->file, &(dummy->xlambda0), f->swapped);
  io_util_readfloat(f->file, &(dummy->h0t0), f->swapped);
  io_util_readfloat(f->file, &(dummy->omegab0), f->swapped);
  io_util_readfloat(f->file, &(dummy->sigma80), f->swapped);
  io_util_readfloat(f->file, &(dummy->ztime0), f->swapped);
  io_util_readfloat(f->file, &(dummy->e0), f->swapped);
  for(iloop=0; iloop<100-43; iloop++)
    io_util_readfloat(f->file, &(dummy->fill2[iloop]), f->swapped);

	f->header = dummy;

	return dummy;
}

extern void
io_deva_header_del(io_logging_t log, io_deva_header_t *header)
{
	if ( (header == NULL) || (*header == NULL) )
		return;

	free(*header);

	*header = NULL;

	return;
}

extern void
io_deva_header_write(io_logging_t log,
                       io_deva_header_t header,
                       io_deva_t f)
{
	if ( (header == NULL) || (f == NULL))
		return;

	if (f->mode != IO_FILE_WRITE)
		return;

	if (f->file == NULL)
		return;

	if (header != f->header) {
		io_logging_msg(log, INT32_C(1),
		               "Writing a different header than stored in "
		               "the file object to the file.");
	}

	/* TODO: Write the header */

	return;
}


extern void
io_deva_header_log(io_logging_t log, io_deva_header_t header)
{
  io_logging_msg(log, INT32_C(5),
                 "  itime:                         %" PRIi32,
                 header->itime);
  io_logging_msg(log, INT32_C(5),
                 "  irun:                          %" PRIi32,
                 header->irun);
  io_logging_msg(log, INT32_C(5),
                 "  nobj:                          %" PRIi32,
                 header->nobj);
  io_logging_msg(log, INT32_C(5),
                 "  ngas:                          %" PRIi32,
                 header->ngas);
  io_logging_msg(log, INT32_C(5),
                 "  ndark:                         %" PRIi32,
                 header->ndark);
  io_logging_msg(log, INT32_C(5),
                 "  ztime:                         %e",
                 1./header->atime-1.);
  io_logging_msg(log, INT32_C(5),
                 "  atime:                         %e",
                 header->atime);
  io_logging_msg(log, INT32_C(5),
                 "  b100:                          %e",
                 header->box100);
  io_logging_msg(log, INT32_C(5),
                 "  omega0:                        %e",
                 header->omega0);
  io_logging_msg(log, INT32_C(5),
                 "  xlambda0:                      %e",
                 header->xlambda0);
  io_logging_msg(log, INT32_C(5),
                 "  omegab0:                       %e",
                 header->omegab0);
  io_logging_msg(log, INT32_C(5),
                 "  h100:                          %e",
                 header->h100);
  io_logging_msg(log, INT32_C(5),
                 "  sigma80:                       %e",
                 header->sigma80);
  io_logging_msg(log, INT32_C(5),
                 "  h0t0:                          %e",
                 header->h0t0);
  io_logging_msg(log, INT32_C(5),
                 "  rmtot:                         %e",
                 header->rmtot);
  io_logging_msg(log, INT32_C(5),
                 "  rmdark:                        %e",
                 header->rmdark);
  io_logging_msg(log, INT32_C(5),
                 "  rmgas:                         %e",
                 header->rmgas);
  
  return;
}


/**********************************************************************\
 *    Implementation of local functions                               * 
\**********************************************************************/
