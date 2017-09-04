/**
 * \file loadbalance.c
 *
 * Provides functionality for loadbalancing.
 */

void dummy_loadbalance()
{
}

/************************************************************\
 * We only need that stuff, when in MPI mode, otherwise the *
 * particle structure will look different anyways.          *
\************************************************************/
#ifdef WITH_MPI

#include "loadbalance.h"

#include <mpi.h>
#include <stdlib.h>
#include <inttypes.h>
/* Required for memcpy() and round() */
#include <string.h>
#if (NONC99 == 0)
#	include <math.h>
#else
#	include <replace_math.h>
#endif

#include "../common.h"


/**********************************************************************\
 *    Prototypes of local functions                                   *
\**********************************************************************/

/**
 * \brief Does a loadbalancing based on equal computational volume per
 *        processor.
 *
 * This function gives every processor the same number of cells,
 * allowing for a mismatch of 1 cell (if the number of cells is not
 * divisable by the number of processors). 
 *
 * \param log      The logging module.
 * \param loadbal  The loadbalancing structure which will be modified to
 *                 hold the new distribution.
 *
 * \return Nothing. The loadbalance structure will be changed though.
 */
static void
local_equalvol(io_logging_t log, loadbalance_t loadbal);

/**
 * \brief Does a loadbalancing based on equals number of particles per
 *        processor.
 *
 * This function tries to distribute the cells in such a way over the
 * processors, that in the end every processor has the same number of
 * particles.
 *
 * \param log      The logging module.
 * \param loadbal  The loadbalancing structure which will be modified to
 *                 hold the new distribution.
 *
 * \return Nothing. The loadbalance structure will be changed though.
 */
static void
local_equalpart(io_logging_t log, loadbalance_t loadbal);

inline static void
local_recount(io_logging_t log,
              loadbalance_t loadbal,
              partptr part,
              uint64_t no_parts);

inline static void
local_calcDistribution(io_logging_t log,
                       loadbalance_t loadbal);

inline static void
local_recalc_localparts(io_logging_t log,
                        loadbalance_t loadbal);

inline static void
local_recalc_parts(io_logging_t log,
                   loadbalance_t loadbal);


/**********************************************************************\
 *    Implementation of global functions                              * 
\**********************************************************************/
extern const char *
loadbalance_schemestr(loadbalance_scheme_t scheme)
{
	switch (scheme) {
		case LOADBALANCE_EQUALVOL:
			return LOADBALANCE_EQUALVOL_STR;
		case LOADBALANCE_EQUALPART:
			return LOADBALANCE_EQUALPART_STR;
	}

	return "Unknown scheme";
}

extern loadbalance_t
loadbalance_new(io_logging_t log,
                loadbalance_scheme_t scheme,
                sfc_curve_t ctype,
                int level,
                int ncpu)
{
	loadbalance_t dummy;

	io_logging_msg(log, INT32_C(3),
	               "Getting memory for the loadbalance structure");
	dummy = (loadbalance_t)malloc(sizeof(struct loadbalance_struct));
	if (dummy == NULL) {
		io_logging_memfatal(log, "loadbalance structure");
		common_terminate(EXIT_FAILURE);
	}

	dummy->scheme = scheme;
	dummy->ctype = ctype;
	dummy->level = level;
	dummy->startkey = (sfc_key_t)0;
	dummy->totkeys = 1<<(level*3);
	io_logging_msg(log, INT32_C(3),
	               "Getting memory for the bf array (%fMB)",
	               sizeof(uint32_t)*dummy->totkeys/(1024.*1024.));
	dummy->bf = (uint32_t *)malloc(sizeof(uint32_t)*dummy->totkeys);
	if (dummy->bf == NULL) {
		io_logging_memfatal(log, "bf array");
		free(dummy);
		common_terminate(EXIT_FAILURE);
	}
	io_logging_msg(log, INT32_C(3),
	               "Getting memory for the local bf array (%fMB)",
	               sizeof(uint32_t)*dummy->totkeys/(1024.*1024.));
	dummy->loc_bf = (uint32_t *)malloc(sizeof(uint32_t)*dummy->totkeys);
	if (dummy->bf == NULL) {
		io_logging_memfatal(log, "local bf array");
		free(dummy->bf);
		free(dummy);
		common_terminate(EXIT_FAILURE);
	}

	dummy->ncpu = ncpu;
	io_logging_msg(log, INT32_C(3),
	               "Getting memory for fst- and lstkey arrays (2x%fkB)",
	               sizeof(sfc_key_t)*ncpu/1024.);
	dummy->fstkey = (sfc_key_t *)malloc(sizeof(sfc_key_t)*ncpu*2);
	if (dummy->fstkey == NULL) {
		io_logging_memfatal(log, "fst- and lstkey arrays");
		free(dummy->loc_bf);
		free(dummy->bf);
		free(dummy);
		common_terminate(EXIT_FAILURE);
	}
	dummy->lstkey = dummy->fstkey + ncpu;

	io_logging_msg(log, INT32_C(3),
	               "Getting memory for (boundary) particle per cpu "
	               "arrays (2x%fkB)",
	               sizeof(uint64_t)*ncpu/1024.);
	dummy->no_parts = (uint64_t *)malloc(sizeof(uint64_t)*ncpu*2);
	if (dummy->no_parts == NULL) {
		io_logging_memfatal(log, "(boundary) particle per cpu arrays");
		free(dummy->fstkey);
		free(dummy->loc_bf);
		free(dummy->bf);
		free(dummy);
		common_terminate(EXIT_FAILURE);
	}
	dummy->no_parts_loc = dummy->no_parts + ncpu;

	dummy->bound = NULL;

	io_logging_msg(log, INT32_C(3),
	               "Done with creating loadbalance structure.");

	return dummy;
}

extern void
loadbalance_del(loadbalance_t *loadbal)
{
	/* Catch wrong usage */
	if (loadbal == NULL)
		return;

	/* If it is an unfreed loadbalance structure, dispose everything */
	if (*loadbal != NULL) {
		sfc_boundary_2_del(&((*loadbal)->bound));
		free((*loadbal)->no_parts);
		free((*loadbal)->fstkey);
		free((*loadbal)->loc_bf);
		free((*loadbal)->bf);
		free(*loadbal);
		*loadbal = NULL;
	}

	/* Return */
	return;
}

extern void
loadbalance_update(io_logging_t log,
                   loadbalance_t loadbal,
                   partptr part,
                   uint64_t no_parts)
{
	uint64_t i;
	int j;
	
	local_recount(log, loadbal, part, no_parts);
	local_calcDistribution(log, loadbal);
	local_recalc_localparts(log, loadbal);

	return;
}

extern void
loadbalance_update_plain(io_logging_t log,
                         loadbalance_t loadbal,
                         partptr part,
                         uint64_t no_parts)
{
	local_recount(log, loadbal, part, no_parts);
	local_recalc_localparts(log, loadbal);
	local_recalc_parts(log, loadbal);
}

extern void
loadbalance_log(io_logging_t log, loadbalance_t lb)
{
	uint64_t i;
	int j;

	io_logging_msg(log, INT32_C(5), "Loadbalance structure:");
	io_logging_msg(log, INT32_C(5),
	               "  scheme:          %s",
	               loadbalance_schemestr(lb->scheme));
	io_logging_msg(log, INT32_C(5),
	               "  level:           %i", lb->level);
	io_logging_msg(log, INT32_C(5),
	               "  startkey:        %" SFC_PRIkey, lb->startkey);
	io_logging_msg(log, INT32_C(5),
	               "  totkeys:         %" PRIu64, lb->totkeys);
	io_logging_msg(log, INT32_C(5),
	               "  bf:              skipping %fkB",
	               lb->totkeys*sizeof(uint32_t)/1024.);
	io_logging_msg(log, INT32_C(5),
	               "  loc_bf:          skipping %fkB",
	               lb->totkeys*sizeof(uint32_t)/1024.);
	io_logging_msg(log, INT32_C(5),
	               "  ncpu:            %i", lb->ncpu);
	io_logging_msg(log, INT32_C(5), "  fstkey/lstkey:");
	for (j=0; j<lb->ncpu; j++)
		io_logging_msg(log, INT32_C(5),
	               "    %04i    %" SFC_PRIkey "/%" SFC_PRIkey,
		           j, lb->fstkey[j], lb->lstkey[j]);
	io_logging_msg(log, INT32_C(5), "  no_parts/no_parts_loc:");
	for (j=0; j<lb->ncpu; j++)
		io_logging_msg(log, INT32_C(5),
	               "    %04i:   %" PRIu64 "/%" PRIu64,
		           j, lb->no_parts[j], lb->no_parts_loc[j]);
	io_logging_msg(log, INT32_C(5), "  Inner boundary:");
	for (j=0; j<lb->bound->num_inner; j++) {
		io_logging_msg(log, INT32_C(5),
		               "    To CPU %i, Quality: %" SFC_PRIkey
		               "/%" SFC_PRIkey " = %g", j, 
		                  (lb->bound->inner+j)->maxkey
		                - (lb->bound->inner+j)->minkey + 1,
		               (lb->bound->inner+j)->num,
		                 ((double)( (lb->bound->inner+j)->maxkey
		                           -(lb->bound->inner+j)->minkey+1))
		               / ((double)((lb->bound->inner+j)->num)));
	}
	io_logging_msg(log, INT32_C(5), "  Outer boundary:");
	io_logging_msg(log, INT32_C(5),
	               "    Quality: %" SFC_PRIkey "/%" SFC_PRIkey " = %g",
	               lb->bound->outer->maxkey - lb->bound->outer->minkey + 1,
	               lb->bound->outer->num,
	                 ((double)( lb->bound->outer->maxkey
	                           -lb->bound->outer->minkey+1))
	               / ((double)(lb->bound->outer->num)));

}


/**********************************************************************\
 *    Implementation of local functions                               * 
\**********************************************************************/
static void
local_equalvol(io_logging_t log, loadbalance_t loadbal)
{
	uint64_t cellpercpu;
	uint64_t cellsleft;
	uint64_t j;
	int i;

	/*
	 * How many cells per CPU and how many are left (might not work
	 * evenly)
	 */
	cellpercpu = loadbal->totkeys/loadbal->ncpu;
	cellsleft = loadbal->totkeys%loadbal->ncpu;

	/*
	 * Set the numbers for the first CPU, and give it one more cell if
	 * there are some left. Also set the number of particle the CPU
	 * holds correctly.
	 */
	loadbal->fstkey[0] = (sfc_key_t)0;
	loadbal->lstkey[0] = (sfc_key_t)cellpercpu;
	if (cellsleft > 0) {
		loadbal->lstkey[0] += (sfc_key_t)1;
		cellsleft--;
	}
	loadbal->no_parts[0] = UINT64_C(0);
	for (j=loadbal->fstkey[0]; j<=loadbal->lstkey[0]; j++)
		loadbal->no_parts[i] += (uint64_t)(loadbal->bf[j]);

	/*
	 * Now do it for the rest of the CPU, if still required they get
	 * one cell more.  Also set the number of particles on that CPU
	 * correctly.
	 */
	for (i=1; i<loadbal->ncpu; i++) {
		loadbal->fstkey[i] = (sfc_key_t)(loadbal->lstkey[i-1]+1);
		loadbal->lstkey[i] = loadbal->fstkey[i] + cellpercpu;
		if (cellsleft > 0) {
			loadbal->lstkey[i] += (sfc_key_t)1;
			cellsleft--;
		}
		loadbal->no_parts[i] = UINT64_C(0);
		for (j=loadbal->fstkey[i]; j<=loadbal->lstkey[i]; j++)
			loadbal->no_parts[i] += (uint64_t)(loadbal->bf[j]);
	}

	return;
}

static void
local_equalpart(io_logging_t log, loadbalance_t loadbal)
{
	uint64_t i;
	uint64_t partpercpu;
	uint64_t no_part;
	uint64_t totnum;
	int32_t cpu;

	/* Count all particles */
	no_part = loadbal->bf[0];
	for (i=1; i<loadbal->totkeys; i++) {
		no_part += loadbal->bf[i];
	}

	/* Setup the loop */
	cpu = INT32_C(0);
	loadbal->no_parts[0] = UINT64_C(0);
	totnum = UINT64_C(0);
	loadbal->fstkey[0] = (sfc_key_t)0;
	partpercpu = (uint64_t)round((double)no_part / (double)loadbal->ncpu);

	/* Loop over all loadbalance domain cells */
	for (i=0; i<loadbal->totkeys; i++) {
		/* Now check if we should include this cell on the current    *\
		 * CPU, or if we have to, because there are no more CPUs to   *
		\* distribute particles to                                    */
		if (   (   labs((long)(loadbal->no_parts[cpu] - partpercpu))
		     < labs((long)(  loadbal->no_parts[cpu]
		                   + loadbal->bf[i] - partpercpu)))
		    && ( cpu < loadbal->ncpu-1 ) ) {
			/* Got as close as possible to partpercpu in step i-1 */
			loadbal->lstkey[cpu] = (sfc_key_t)(i-1);
			cpu++;
			loadbal->no_parts[cpu] = loadbal->bf[i];
			totnum += loadbal->bf[i];
			/* Now recalculate partpercpu in order to adapt *\
			\* to the new situation                         */
			partpercpu =   (int64_t)round((double)(no_part-totnum)
			             / (double)(loadbal->ncpu-cpu));
			loadbal->fstkey[cpu] = (sfc_key_t)i;
		} else {
			loadbal->no_parts[cpu] += loadbal->bf[i];
			totnum += loadbal->bf[i];
		}
	}
	loadbal->lstkey[loadbal->ncpu - 1] = loadbal->totkeys-1;

	return;
}

inline static void
local_recount(io_logging_t log,
              loadbalance_t loadbal,
              partptr part,
              uint64_t no_parts)
{
	uint64_t i;

	/* Set the counting array */
	io_logging_msgplain(log, INT32_C(2),
	                    "Setting local counting array to 0... ");
	for (i=0; i<loadbal->totkeys; i++) {
		loadbal->loc_bf[i] = UINT32_C(0);
	}
	io_logging_msgplain(log, INT32_C(2), "done\n");
	
	/* Update the SFC-Key and do the counting */
	io_logging_msgplain(log, INT32_C(2),
	                    "Updating SFC-Key and counting... ");
	for (i=0; i<no_parts; i++, part++) {
		part->sfckey = sfc_curve_calcKey(loadbal->ctype,
		                                 (double)(part->pos[0]),
		                                 (double)(part->pos[1]),
		                                 (double)(part->pos[2]),
		                                 BITS_PER_DIMENSION);
		/* Count that particle */
		loadbal->loc_bf[sfc_curve_contract(loadbal->level,
		                                   BITS_PER_DIMENSION,
		                                   loadbal->ctype,
		                                   part->sfckey)]++;
	}
	io_logging_msgplain(log, INT32_C(2), "done\n");

#	ifdef MPI_DEBUG
	for (i=0; i < loadbal->totkeys; i++) {
		io_logging_msg(log, INT32_C(11),
		               "local_bucket[%05" PRIu64 "] = %07" PRIu32,
		               i, loadbal->loc_bf[i]);
	}
#	endif

	/* Establish global bf status and cleaning up */
	MPI_Allreduce((void *)(loadbal->loc_bf), (void *)(loadbal->bf),
	              loadbal->totkeys,
	              global_mpi.dt_uint32, global_mpi.op_sum,
	              MPI_COMM_WORLD);

#	ifdef MPI_DEBUG
	for (i=0; i < loadbal->totkeys; i++) {
		io_logging_msg(log, INT32_C(10),
		               "global_bucket[%05" PRIu64 "] = %07" PRIu32,
		               i, loadbal->bf[i]);
	}
#	endif

	return;
}

inline static void
local_calcDistribution(io_logging_t log,
                       loadbalance_t loadbal)
{
	uint64_t i;

	/* Now actually calculate the distribution */
	io_logging_msg(log, INT32_C(2),
	               "Calculating the distribution of the cells.");
	switch(loadbal->scheme) {
		case LOADBALANCE_EQUALVOL:
			local_equalvol(log, loadbal);
			break;
		case LOADBALANCE_EQUALPART:
			local_equalpart(log, loadbal);
			break;
		default:
			io_logging_fatal(log, "Unkown distribution scheme.");
			common_terminate(EXIT_FAILURE);
	}
#	ifdef MPI_DEBUG
	for (i=0; i < loadbal->ncpu; i++) {
		io_logging_msg(log, INT32_C(5),
		               "cpu [%02" PRIu64 "] holds (inclusive) "
		               "keys %" SFC_PRIkey " to %" SFC_PRIkey
		               " (total of %" PRIu64 " particles).",
		               i, loadbal->fstkey[i], loadbal->lstkey[i],
		               loadbal->no_parts[i]);
	}
#	endif

	/* Throw away the old boundary information */
	if (loadbal->bound != NULL)
		sfc_boundary_2_del(&(loadbal->bound));
	/* Generate new boundary and die if that fails */
	loadbal->bound = sfc_boundary_2_get(global_io.log,
	                                    (uint32_t)global_mpi.rank,
	                                    (uint32_t)global_mpi.size,
	                                    loadbal->fstkey,
	                                    loadbal->lstkey,
	                                    LOADBALANCE_DOMAIN_LEVEL,
	                                    SFC_CURVE_HILBERT);
# ifdef MPI_DEBUG
	sfc_boundary_log(global_io.log, loadbal->bound->outer);
# endif
	if (loadbal->bound == NULL) {
		io_logging_fatal(log, "Generation of boundary failed. Aborting.");
		common_terminate(EXIT_FAILURE);
	}
	
	return;
}

inline static void
local_recalc_localparts(io_logging_t log,
                        loadbalance_t loadbal)
{
	uint64_t i;
	int j;

	io_logging_msg(log, INT32_C(2),
	               "Calculating how many particle of which process "
	               "we are holding.");
	for (j=0; j<loadbal->ncpu; j++) {
		loadbal->no_parts_loc[j] = UINT64_C(0);
		for (i=loadbal->fstkey[j]; i<=loadbal->lstkey[j]; i++) {
			loadbal->no_parts_loc[j] += (uint64_t)(loadbal->loc_bf[i]);
		}
	}
}

inline static void
local_recalc_parts(io_logging_t log,
                   loadbalance_t loadbal)
{
	uint64_t i;
	int j;

	io_logging_msg(log, INT32_C(2),
	               "Calculating how many particles are supposed"
	               "to be on which CPU");
	for (j=0; j<loadbal->ncpu; j++) {
		loadbal->no_parts[j] = UINT64_C(0);
		for (i=loadbal->fstkey[j]; i<=loadbal->lstkey[j]; i++) {
			loadbal->no_parts[j] += (uint64_t)(loadbal->bf[i]);
		}
	}
}



#ifdef VOLLMURKSIGNORE
extern partptr
loadbalance_dummy(io_logging_t log,
                  partptr fstprt,
                  uint64_t *no_part,
                  int32_t ncpu)
{
	/** Loop counter */
	uint64_t i;
	/**
	 * Stores how many particles have been send already; needed in the
	 * case that more particles than MAX_SEND_PARTICLES have to be
	 * send.
	 */
	uint64_t skip_particles = 0;
	/** Stores how many particles need to be send (regular + boundary) */
	uint64_t totalsend = 0;
	/** Holds the bucketfill status */
	uint64_t *bf;
	/** Holds the intermediate bucketfill status (local values) */
	uint64_t *bf_tmp;
	/** Particle iterator */
	partptr part;
	/** Loadbalance structure */
	loadbalance_t loadbal;

	/* Calculate how many keys are taken into account */
	global_info.totkeys = 1<<(LOADBALANCE_DOMAIN_LEVEL*3);

	/* Aquire the counting arrays */
	bf = (uint64_t *)calloc(global_info.totkeys, sizeof(uint64_t));
	if (bf == NULL) {
		/* Bad, no memory left */
		io_logging_memfatal(log,
		                    "bucket fill array, aborting");
		common_terminate(EXIT_FAILURE);
	}
	bf_tmp = (uint64_t *)calloc(global_info.totkeys, sizeof(uint64_t));
	if (bf_tmp == NULL) {
		/* Bad, no memory left */
		io_logging_memfatal(log,
		                    "bucket fill array, aborting");
		free(bf);
		common_terminate(EXIT_FAILURE);
	}

	/* First calculate the right SFC Key for every particle */
	part = fstprt;
	for (i=0; i<*no_part; i++, part++) {
		part->sfckey = sfc_curve_calcKey(global_info.ctype,
		                                 (double)(part->pos[0]),
		                                 (double)(part->pos[1]),
		                                 (double)(part->pos[2]),
		                                 BITS_PER_DIMENSION);
		/* Count that particle */
		bf_tmp[sfc_curve_contract(LOADBALANCE_DOMAIN_LEVEL,
		                          BITS_PER_DIMENSION,
		                          global_info.ctype,
		                          part->sfckey)]++;
	}

#	ifdef MPI_DEBUG
	for (i=0; i < global_info.totkeys; i++) {
		io_logging_msg(log, INT32_C(10),
		               "local_bucket[%05" PRIu64 "] = %07" PRIu64,
		               i, bf_tmp[i]);
	}
#	endif

	/* Establish global bf status */
	MPI_Allreduce((void *)bf_tmp, (void *)bf, global_info.totkeys,
	              global_mpi.dt_uint64, global_mpi.op_sum,
	              MPI_COMM_WORLD);

#	ifdef MPI_DEBUG
	for (i=0; i < global_info.totkeys; i++) {
		io_logging_msg(log, INT32_C(10),
		               "global_bucket[%05" PRIu64 "] = %07" PRIu64,
		               i, bf[i]);
	}
#	endif

	/* Clean up */
	free(bf_tmp);

	/* Do the loadbalancing */
	loadbal = loadbalance_equalpart(log, ncpu, bf);

	/* Calculate from which to which key my area belongs */
	global_info.minkey = loadbal->fstkey[global_mpi.rank];
	if (global_mpi.rank == global_mpi.size-1) {
		global_info.maxkey = loadbal->totkeys;
	} else {
		global_info.maxkey = loadbal->fstkey[global_mpi.rank+1] - 1;
	}

	/* XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX 
	   XXX     The following part only works if only the first     XXX 
	   XXX     process is reading particles from the file.         XXX 
	   XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX */
	if (global_mpi.rank == 0) {
		/** Loop counter */
		uint64_t j;
		/** Particle iterator */
		partptr part;
		/** particle array */
		partptr *newparts;
		/** Counter array */
		partptr *parts;
		sfc_key_t redkey, minkey, maxkey;
		sfc_boundary_t *bounds;

		/**
		 * Create an array to store the boundaries for all processors.
		 */
		bounds = (sfc_boundary_t *)malloc( sizeof(sfc_boundary_t)
		                                  *ncpu);
		if (bounds == NULL) {
			/* Bad, no memory left */
			io_logging_memfatal(log,
			                    "boundary array, aborting");
			common_terminate(EXIT_FAILURE);
		}
	
		/* Get the boundary */
		maxkey = 0;
		for (i=0; i<ncpu-1; i++) {
			minkey = maxkey;
			maxkey = loadbal->fstkey[i+1] - 1;
			bounds[i] = sfc_boundary_get(log,
			                             minkey,
			                             maxkey,
			                             LOADBALANCE_DOMAIN_LEVEL,
			                             global_info.ctype);
			maxkey++;
		}
		minkey = maxkey;
		maxkey = (1<<(LOADBALANCE_DOMAIN_LEVEL*3)) - 1;
		bounds[i] = sfc_boundary_get(log,
		                             minkey,
		                             maxkey,
		                             LOADBALANCE_DOMAIN_LEVEL,
		                             global_info.ctype);
		/* Loop over all boundary cells to figure out how many particles
		 * are stored there
		 */
		for (i=0; i<ncpu; i++) {
			sfc_boundary_log(log, bounds[i]);
			loadbal->no_parts_bound[i] = 0;
			for (j=0; j<((bounds[i])->num); j++) {
				loadbal->no_parts_bound[i] += bf[(bounds[i])->bound[j]];
			}
		}
	
		/******************************\ 
		 * Now split up the particles *
		\******************************/
	
		/* First create particle arrays for each processor */
		newparts = (partptr *)malloc(sizeof(partptr) * ncpu);
		if (newparts == NULL) {
			io_logging_memfatal(log,
			                    "new particle pointer array, aborting.");
			common_terminate(EXIT_FAILURE);
		}
		parts = (partptr *)malloc(sizeof(partptr) * ncpu);
		if (parts == NULL) {
			io_logging_memfatal(log,
			                    "new particle counter array, aborting.");
			common_terminate(EXIT_FAILURE);
		}
		for (i=0; i<ncpu; i++) {
#			ifdef MPI_DEBUG
			io_logging_msg(log, INT32_C(5),
			               "I want to mallocate storage for %" PRIu64
			               " particles (%" PRIu64 " bytes), plus %"
			               PRIu64 "  particles (% " PRIu64 " bytes)"
			               " in the boundary.",
			               loadbal->no_parts[i],
			               loadbal->no_parts[i]*sizeof(struct particle),
			               loadbal->no_parts_bound[i],
			                 loadbal->no_parts_bound[i]
			                *sizeof(struct particle));
			               
#			endif
			newparts[i] = (partptr)malloc( sizeof(struct particle)
			                              *( loadbal->no_parts[i]
			                                +loadbal->no_parts_bound[i]));
			if (newparts[i] == NULL) {
				io_logging_memfatal(log,
				                    "new particle array (size %010"
				                    PRIu64 "bytes), aborting.",
				                    (uint64_t)( sizeof(struct particle)
				                               *(loadbal->no_parts[i]
				                          +loadbal->no_parts_bound[i])));
				common_terminate(EXIT_FAILURE);		
			}
			/* Don't forget to set the counter array to the start */
			parts[i] = newparts[i];
		}
		
		/* Loop over all particles to distribute them */
		part = fstprt;
		for (i=0; i<*no_part; i++, part++) {
			redkey = sfc_curve_contract(LOADBALANCE_DOMAIN_LEVEL,
			                            BITS_PER_DIMENSION,
			                            global_info.ctype,
			                            part->sfckey);
			/* Find where that particle belongs to */
			j = 0;
			do {
				j++;
				if (redkey < loadbal->fstkey[j])
					break;
			} while ( (j < ncpu) );
			j--;

			/* Copy the particle */
			(void *)memcpy((void*)(parts[j]),
			               (void *)part,
			               sizeof(struct particle));
	
			/* Increment the right particle counter */
			parts[j]++;
		}
		io_logging_msg(log, INT32_C(3),
		               "Copied particles according to domain decomp.");
		/* Now we loop again to fit the boundary particles at the end */
		part = fstprt;
		for (i=0; i<*no_part; i++, part++) {
			redkey = sfc_curve_contract(LOADBALANCE_DOMAIN_LEVEL,
			                            BITS_PER_DIMENSION,
			                            global_info.ctype,
			                            part->sfckey);
			for (j=0; j<ncpu; j++) {
				if (sfc_boundary_findKeyPos(log,
				                            bounds[j],
				                            redkey,
				                            NULL)) {
					(void *)memcpy((void*)(parts[j]),
					               (void *)part,
					               sizeof(struct particle));
					parts[j]++;
				}
			}
		}

		/* Now send the other particles around */
		for (i=1; i<ncpu; i++) {
			/* FIXME Change that to MPI_Isend! FIXME*/
			/* This will send the number of particles that need to be
			 * recieved */
			totalsend =   loadbal->no_parts[i] 
			            + loadbal->no_parts_bound[i];
			MPI_Send(&totalsend,
			         1,
			         global_mpi.dt_uint64,
			         i, 0, MPI_COMM_WORLD);
			/* Send the actual particle array bytewise */
			skip_particles = 0;
			while (  (totalsend - skip_particles)
			        > MAX_SEND_PARTICLES) {
				io_logging_msg(log, INT32_C(0),
				               "Have to send too many particles. "
				               "Sending in parts.");
				MPI_Send(newparts[i]+skip_particles,
				         MAX_SEND_PARTICLES * sizeof(struct particle),
				         MPI_BYTE,
				         i, 0, MPI_COMM_WORLD);
				skip_particles += MAX_SEND_PARTICLES;
				io_logging_msg(log, INT32_C(0),
				               "Sent a pack of %" PRIu64 " particles. "
				               "A total of %" PRIu64 " particles has "
				               "been sent sofar.",
				               MAX_SEND_PARTICLES, skip_particles);
			}
			io_logging_msg(log, INT32_C(0),
			               "Now sending the rest of %" PRIu64
			               " particles.",
			               (uint64_t)(  totalsend
			                          - skip_particles));
			MPI_Send(newparts[i]+skip_particles,
			           (totalsend - skip_particles)
			         * sizeof(struct particle),
			         MPI_BYTE,
			         i, 0, MPI_COMM_WORLD);
		}
	
		/* Throw away all particles from the large array */
		free(fstprt);

		/* But I keep mine */
		fstprt = newparts[0];

		/* Remember how many particles are mine */
		*no_part = loadbal->no_parts[0] + loadbal->no_parts_bound[0];

		/* Clean up */
		sfc_boundary_del(&(bounds[0]));
		for (i=1; i<ncpu; i++) {
			free(newparts[i]);
			sfc_boundary_del(&(bounds[i]));
		}
		free(bounds);
		free(newparts);
		free(parts);
	} else {
		MPI_Status stat;

		/* Figure out how many particles need to be recieved */
		MPI_Recv(no_part, 1, global_mpi.dt_uint64,
		         0, 0, MPI_COMM_WORLD, &stat);
		io_logging_msg(log, INT32_C(2),
		               "Needing to recieve %" PRIu64 " particles.",
		               *no_part);

		/* Create a storage that is large enough for all particles */
		fstprt = (partptr)malloc(sizeof(struct particle)* (*no_part));
		if (fstprt == NULL) {
			io_logging_memfatal(log,
			                    "particle storage, aborting.");
			common_terminate(EXIT_FAILURE);
		}

		/* Recieve the particles bytewise */
		skip_particles = 0;
		while ( (*no_part - skip_particles) > MAX_SEND_PARTICLES) {
			io_logging_msg(log, INT32_C(0),
			               "Have to recieve too many particles. "
			               "Recieving in parts.");
			MPI_Recv(fstprt+skip_particles,
			         MAX_SEND_PARTICLES * sizeof(struct particle),
			         MPI_BYTE,
			         0, 0, MPI_COMM_WORLD, &stat);
			skip_particles += MAX_SEND_PARTICLES;
			io_logging_msg(log, INT32_C(0),
			               "Recieved a pack of %" PRIu64 " particles. "
			               "A total of %" PRIu64 " particles has "
			               "been recieved sofar.",
			               MAX_SEND_PARTICLES, skip_particles);
		}
		io_logging_msg(log, INT32_C(0),
		               "Now recieving the rest of %" PRIu64
		               " particles.",
		               (uint64_t)(*no_part - skip_particles));
		MPI_Recv(fstprt+skip_particles,
		         (*no_part - skip_particles) * sizeof(struct particle),
		         MPI_BYTE,
		         0, 0, MPI_COMM_WORLD, &stat);

		/* We are in the happy possession of our particles */
		;
	}

	/* Clean up */
	free(bf);
	
	return fstprt;
}

extern loadbalance_t
loadbalance_equalpart(io_logging_t log, int32_t ncpu, uint64_t *bf)
{
	loadbalance_t loadbal;
	uint64_t i;
	uint64_t partpercpu;
	uint64_t no_part;
	uint64_t totnum;
	int32_t cpu;

	loadbal = local_initloadbal(log, ncpu);

	/* Count all particles */
	no_part = bf[0];
	for (i=1; i<loadbal->totkeys; i++) {
		no_part += bf[i];
	}

	/* Setup the loop */
	cpu = INT32_C(0);
	loadbal->no_parts[0] = UINT64_C(0);
	totnum = UINT64_C(0);
	loadbal->fstkey[0] = (sfc_key_t)0;
	partpercpu = (uint64_t)round((double)no_part / (double)ncpu);

	/* Loop over all loadbalance domain cells */
	for (i=0; i<loadbal->totkeys; i++) {
		/* Now check if we should include this bucket on the current  *\
		 * CPU, or if we have to, because there are no more CPUs to   *
		\* distribute particles to                                    */
		if (   (   labs((long)(loadbal->no_parts[cpu] - partpercpu))
		     < labs((long)(loadbal->no_parts[cpu] + bf[i] - partpercpu)))
		    && ( cpu < ncpu-1 ) ) {
			/* Got as close as possible to partpercpu in step i-1 */
			cpu++;
			loadbal->no_parts[cpu] =  bf[i];
			totnum += bf[i];
			/* Now recalculate partpercpu in order to adapt *\
			\* to the new situation                         */
			partpercpu =   (int64_t)round((double)(no_part-totnum)
			             / (double)(ncpu-cpu));
			loadbal->fstkey[cpu] = (sfc_key_t)i;
		} else {
			loadbal->no_parts[cpu] +=  bf[i];
			totnum += bf[i];
		}
		
	}

#	ifdef MPI_DEBUG
	for (i=0; i < ncpu; i++) {
		io_logging_msg(log, INT32_C(5),
		               "cpu [%02" PRIu64 "] starts at key %07" PRIu64
		               " (total of %07" PRIu64 " particles).",
		               i, (uint64_t)(loadbal->fstkey[i]),
		               loadbal->no_parts[i]);
	}
#	endif

	return loadbal;
}

static loadbalance_t
local_initloadbal(io_logging_t log, int32_t ncpu)
{
	loadbalance_t loadbal;

	/* Get memory for the structure */
	loadbal = (loadbalance_t)malloc(sizeof(struct loadbalance_struct));
	if (loadbal == NULL) {
		io_logging_memfatal(log,
		                    "loadbalance structure, aborting.");
		common_terminate(EXIT_FAILURE);
	}

	/* Try to get memory for the first key array */
	loadbal->fstkey = (sfc_key_t *)calloc(ncpu, sizeof(sfc_key_t));
	if (loadbal->fstkey == NULL) {
		io_logging_memfatal(log,
		                    "first key array, aborting.");
		free(loadbal);
		common_terminate(EXIT_FAILURE);
	}

	/* Try to get memory for the no_parts array */
	loadbal->no_parts = (int64_t *)calloc(ncpu, sizeof(int64_t));
	if (loadbal->no_parts == NULL) {
		io_logging_memfatal(log,
		                    "no_parts array, aborting.");
		free(loadbal->fstkey);
		free(loadbal);
		common_terminate(EXIT_FAILURE);
	}

	/* Try to get memory for the no_parts_bound array */
	loadbal->no_parts_bound = (int64_t *)calloc(ncpu, sizeof(int64_t));
	if (loadbal->no_parts_bound == NULL) {
		io_logging_memfatal(log,
		                    "no_parts_bound array, aborting.");
		free(loadbal->no_parts);
		free(loadbal->fstkey);
		free(loadbal);
		common_terminate(EXIT_FAILURE);
	}

	/* Setup the rest */
	loadbal->ncpu = ncpu;
	loadbal->totkeys = (1 << (LOADBALANCE_DOMAIN_LEVEL * 3));

	return loadbal;
}
#endif /* VOLLMURKSIGNORE */


#endif /* WITH_MPI */
