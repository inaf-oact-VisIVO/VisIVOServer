#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "common.h"
#ifdef WITH_MPI
#	include <mpi.h>
#	include "libutility/loadbalance.h"
#	include "comm.h"
#endif
#ifdef NEWSTARTRUN
#	include "libsfc/sfc.h"
#endif

/* ...and now for the actual includes */
#include "startrun.h"
#include "step.h"
#include "libutility/utility.h"
#include "libparticles/particles.h"
#include "libgravity/gravity.h"
#include "libgrids/grids.h"
#include "libio_serial/io_serial.h"

#ifdef HYDRO
#include "libmhd/mhd.h"
#endif

#ifdef WITH_MPI
static void
local_communicate_all(void);
#endif

int main_AHFstep(char *AMIGA_input,int *argc, char ***argv)
{
   gridls  *grid_list;        /* pointer to list of grids            */
   gridls  *for_grid;         /* for looping over all grids          */
   
   int     n_grids;
   int     no_grids;          /* total number of grids               */
   int     no_timestep;       /* number of coarse grid timesteps     */
   int     no_first_timestep; /* number of initial timestep          */
      
   double  timecounter;       /* time variable                       */
   double  timestep;          /* timestep size                       */
   double  timecounter_final; /* for all sorts of tests...           */
   
//   char     command[MAXSTRING], AMIGA_input[MAXSTRING];
   char     command[MAXSTRING];
   FILE    *termfile;
  
#ifdef WITH_MPI
   uint64_t newparts;
#endif

  /* we always read the relevant parameters from an input file! */
/*  if(argc<2)
    {
      fprintf(stderr,"usage: %s AMIGA AMIGA.input\n",*argv);
      exit(1);
    }
  strcpy(AMIGA_input, argv[1]);*/
  
#if (!defined WITH_MPI)

  fprintf(stderr,"============================================================\n");
  fprintf(stderr,"\t     AA      M   M    I    GGGG      AA    \n");
  fprintf(stderr,"\t   A   A    MM MM    I    G        A   A   \n");
  fprintf(stderr,"\t  AAAAA    M M M    I    G GGG    AAAAA   \n");
  fprintf(stderr,"\t A   A    M   M    I    G   G    A   A   \n");
  fprintf(stderr,"\tA   A    M   M    I     GGG     A   A  (v%3.1f/%d)\n",VERSION,BUILT);
  fprintf(stderr,"============================================================\n\n");
#endif

#	ifdef WITH_MPI
	/* Initialize the MPI environment */
	common_initmpi(&argc, &argv); 
#		ifdef MPI_TIMING
	global_mpi.start = MPI_Wtime();
#		endif
#	endif

   /*======================================================== 
    * startrun:    input the initial data from infile 
    *========================================================*/
#	ifdef NEWSTARTRUN
	startrun(AMIGA_input,
	         &timecounter,
	         &timestep,
	         &no_first_timestep);
#		if (defined WITH_MPI && defined MPI_TIMING)
	global_mpi.stop = MPI_Wtime();
	io_logging_msg(global_io.log, INT32_C(1),
	               "Startrun done in %fs",
	               global_mpi.stop-global_mpi.start);
	global_mpi.start = global_mpi.stop;
#		endif

#		ifdef WITH_MPI
	/* Sort the particles in a particle block structure */
	io_logging_section(global_io.log,
	                   "Initial Load-Balancing and Particle Distribution");
   io_logging_subsection(global_io.log, "Loadbalancing");
   loadbalance_update(global_io.log,
                      global_info.loadbal,
                      global_info.fst_part,
                      global_info.no_part);
#			ifdef MPI_TIMING
	global_mpi.stop = MPI_Wtime();
	io_logging_msg(global_io.log, INT32_C(1),
	               "Loadbalance done in %fs",
	               global_mpi.stop-global_mpi.start);
	global_mpi.start = global_mpi.stop;
#			else
   io_logging_msg(global_io.log, INT32_C(1),
                  "Loadbalance done.");
#			endif
   loadbalance_log(global_io.log, global_info.loadbal);
#		else
	/* Generate the SFC keys for all particles */
	for (uint64_t i=0; i<global_info.no_part; i++) {
		partptr part=global_info.fst_part+i;
		part->sfckey = sfc_curve_calcKey(global_info.ctype,
		                                 (double)(part->pos[0]),
		                                 (double)(part->pos[1]),
		                                 (double)(part->pos[2]),
		                                 BITS_PER_DIMENSION);
	}
	/* Sorting all particles to have fast access later on */
	qsort(global_info.fst_part,
	      global_info.no_part,
	      sizeof(part),
	      &cmp_sfckey_part);
#		endif /* WITH_MPI*/

#		ifdef WITH_MPI
   /* Do a first sort of the particles, required for distributing */
   io_logging_subsection(global_io.log, "Sorting particles");
   qsort(global_info.fst_part,
         global_info.no_part,
         sizeof(part),
         &cmp_sfckey_part);
#			ifdef MPI_TIMING
	global_mpi.stop = MPI_Wtime();
	io_logging_msg(global_io.log, INT32_C(1),
	               "Sorting done in %fs",
	               global_mpi.stop-global_mpi.start);
	global_mpi.start = global_mpi.stop;
#			else
   io_logging_msg(global_io.log, INT32_C(1),
                  "Sorting done.");
#			endif

   /* Distribute the particles */
   io_logging_subsection(global_io.log, "Distributing particles");
   comm_dist_part(global_io.log,
                  &(global_info.fst_part),
                  &(global_info.no_part),
                  global_info.loadbal);
#			ifdef MPI_TIMING
	global_mpi.stop = MPI_Wtime();
	io_logging_msg(global_io.log, INT32_C(1),
	               "Sorting done in %fs",
	               global_mpi.stop-global_mpi.start);
	global_mpi.start = global_mpi.stop;
#			else
   io_logging_msg(global_io.log, INT32_C(1),
                  "Distributing done.");
#			endif
   io_logging_msg(global_io.log, INT32_C(0),
                  "Having %"PRIu64" particles!",
                  global_info.no_part);

#			ifdef AHFstep
	/* Do the AHF distribution*/
	io_logging_subsection(global_io.log,
	                      "AHF distribution (duplicating)");
	newparts = comm_dist_part_ahf(global_io.log,
	                              &(global_info.fst_part),
	                              &(global_info.no_part),
	                              global_info.loadbal);
	io_logging_msg(global_io.log, INT32_C(0),
	               "Recieved %"PRIu64" new particles.",
	               newparts);
	/* We need to sort the particles again */
	qsort(global_info.fst_part, global_info.no_part,
	      sizeof(part), &cmp_sfckey_part);
#				ifdef MPI_TIMING
	global_mpi.stop = MPI_Wtime();
	io_logging_msg(global_io.log, INT32_C(1),
	               "AHF distribution done in %fs",
	               global_mpi.stop-global_mpi.start);
	global_mpi.start = global_mpi.stop;
#				else
	io_logging_msg(global_io.log, INT32_C(1),
	               "AHF distribution done.");
#				endif
#			endif /* AHFstep */

#		endif /* WITH_MPI */

#	else /* NEWSTARTRUN */

#		ifdef WITH_MPI
	/* -DWITH_MPI requires -DNEWSTARTRUN */
	fprintf(stderr,
	        "Please use -DWITH_MPI only if -DNEWSTARTRUN "
	        "is defined too.\n");
	common_terminate(EXIT_FAILURE);
#		else
  startrun(AMIGA_input, &timecounter, &timestep, &no_first_timestep);
#		endif
#	endif
   

#	ifdef NEWSTARTRUN
	io_logging_msg(global_io.log, INT32_C(5),
	               "amiga_main:  running with %" PRIu64 " particles",
	               global_info.no_part);

	io_logging_part(global_io.log,
	                "Handing over logging to AMIGA");
#	else
   fprintf(stderr,"\namiga_main:     data read for %ld particles\n", simu.no_part);
   fprintf(stderr,"                using %ld particles\n\n", global.no_part);
#	endif

#if (AHFstep && AHFstep_split_only)
	{
		io_file_t dumpf;
		io_file_strg_struct_t strg;
		char *fname;

		/* Start tge section */
		io_logging_section(global_io.log, "Dumping AHF chunk to file");

		/* First generate the filename */
		fname = (char *)malloc( sizeof(char)
		                       *( strlen(global_io.params->outfile_prefix)
		                         +30));
		if (fname == NULL) {
			io_logging_memfatal(global_io.log, "filename string");
			common_terminate(EXIT_FAILURE);
		}
		sprintf(fname, "%s.chunk.%04i.dump",
		        global_io.params->outfile_prefix,
		        global_mpi.rank);
		io_logging_msg(global_io.log, UINT32_C(0),
		               "Used filename: %s", fname);

		fflush(NULL);
		MPI_Barrier(MPI_COMM_WORLD);

		/* Assign particles to structure */
		strg.posx.val = (void *)(global_info.fst_part->pos);
		strg.posx.stride =   (char *)((global_info.fst_part+1)->pos)
		                   - (char *)(global_info.fst_part->pos);
		strg.posy.val = (void *)(global_info.fst_part->pos+1);
		strg.posy.stride =   (char *)((global_info.fst_part+1)->pos+1)
		                   - (char *)(global_info.fst_part->pos+1);
		strg.posz.val = (void *)(global_info.fst_part->pos+2);
		strg.posz.stride =   (char *)((global_info.fst_part+1)->pos+2)
		                   - (char *)(global_info.fst_part->pos+2);
		strg.momx.val = (void *)(global_info.fst_part->mom);
		strg.momx.stride =   (char *)((global_info.fst_part+1)->mom)
		                   - (char *)(global_info.fst_part->mom);
		strg.momy.val = (void *)(global_info.fst_part->mom+1);
		strg.momy.stride =   (char *)((global_info.fst_part+1)->mom+1)
		                   - (char *)(global_info.fst_part->mom+1);
		strg.momz.val = (void *)(global_info.fst_part->mom+2);
		strg.momz.stride =   (char *)((global_info.fst_part+1)->mom+2)
		                   - (char *)(global_info.fst_part->mom+2);
#	ifdef MULTIMASS
		strg.weight.val = (void *)&(global_info.fst_part->weight);
		strg.weight.stride =   (char *)&((global_info.fst_part+1)->weight)
		                     - (char *)&(global_info.fst_part->weight);
#	else
		strg.weight.val = NULL;
		strg.weight.stride = (ptrdiff_t)0;
#	endif /* MULTIMASS */
#	ifdef GAS_PARTICLES
		strg.u.val = (void *)&(global_info.fst_part->u);
		strg.u.stride =   (char *)&((global_info.fst_part+1)->u)
		                - (char *)&(global_info.fst_part->u);
#	else
		strg.u.val = NULL;
		strg.u.stride = (ptrdiff_t)0;
#	endif /* GAS_PARTICLE */
		strg.id.val = &(global_info.fst_part->id);
		strg.id.stride =   (char *)&((global_info.fst_part+1)->id)
		                 - (char *)&(global_info.fst_part->id);
		strg.bytes_float = sizeof(global_info.fst_part->pos[0]);
		strg.bytes_int = sizeof(global_info.fst_part->id);

		/* Open the dump file now */
		dumpf = io_file_open(global_io.log, fname,
		                     IO_FILE_ARES, IO_FILE_UNKOWN_SWAPPING,
		                     IO_FILE_WRITE, 0);

		/* Write the particles */
		io_file_writepart(global_io.log, dumpf, 0, global_info.no_part,
		                  strg);

		/* Set the header values */
		((io_ares_t)dumpf)->header->no_part = (uint64_t)simu.no_part;
		((io_ares_t)dumpf)->header->no_species = UINT64_C(0);
		((io_ares_t)dumpf)->header->no_vpart = simu.no_vpart;
		((io_ares_t)dumpf)->header->boxsize = simu.boxsize;
		((io_ares_t)dumpf)->header->omega0 = simu.omega0;
		((io_ares_t)dumpf)->header->lambda0 = simu.lambda0;
		((io_ares_t)dumpf)->header->pmass = simu.pmass;
		((io_ares_t)dumpf)->header->minweight = simu.min_weight;
		((io_ares_t)dumpf)->header->maxweight = simu.max_weight;
		((io_ares_t)dumpf)->header->a_initial = simu.a_initial;
		((io_ares_t)dumpf)->header->a_current = global.a;
		((io_ares_t)dumpf)->header->timestep = timestep;
		((io_ares_t)dumpf)->header->minkey =
		       global_info.loadbal->fstkey[global_mpi.rank];
		((io_ares_t)dumpf)->header->maxkey =
		       global_info.loadbal->lstkey[global_mpi.rank];
		((io_ares_t)dumpf)->header->lb_level =
		       global_info.loadbal->level;
		((io_ares_t)dumpf)->header->rank = global_mpi.rank;
		((io_ares_t)dumpf)->header->size = global_mpi.size;

		/* Log the file */
		io_file_log(global_io.log, dumpf);

		/* Close the file and clean up*/		
		io_file_close(global_io.log, &dumpf);
		free(fname);
	}
	common_terminate(EXIT_SUCCESS);
#endif 

#	ifndef NEWSTARTRUN
#ifdef CONVERT
   /* write AMIGA binary of data just read in */
   output(io.outfile_prefix, 0);
#ifdef CONVERT_TERM
   fprintf(stderr,"wrote output file ... exiting now\n");
   common_terminate(0);
#endif
#endif
#	else
	/* FIXME: Writing of output file */
#	endif /* NEWSTARTRUN */
   
   /*===================================================================== 
    * generate the domain grids: simu.NGRID_MIN^3, ...., simu.NGRID_DOM^3 
    *=====================================================================*/
   grid_list = gen_domgrids(&no_grids);   
   
#if (defined HYDRO && !defined INIT_DMHYDRO)
   /* input_grid() allocates the grid from sratch while reading it in from file... */
   free_pquad(global.dom_grid->pquad);
   free(global.dom_grid->pquad);
   
   /* generate global.dom_grid from scratch and fill with values found in IC file */
   input_grid(global.dom_grid);
#endif
   
#ifdef HYDRO_TEST
#if (HYDRO_TEST != 7 && HYDRO_TEST < 100)
   /*===================================================================== 
    *                         HYDRO TESTS
    *=====================================================================*/
   hydro_test(grid_list);
#endif
#endif

   
   /*===================================================================== 
    * build initial linked list 
    *=====================================================================*/
#	if (defined WITH_MPI || defined NEWSTARTRUN)
	ll(global_info.no_part, global_info.fst_part, global.dom_grid);
	global.fst_part = global_info.fst_part;
#	else
   ll(global.no_part, global.fst_part, global.dom_grid);
#	endif

   
#ifdef HYDRO
#if (HYDRO_TEST==7)
   /* assign gas particles to grid */
   zero_udens(global.dom_grid);
   init_udens(global.dom_grid);
   
   /* gas particles are not needed anymore */
   NULL_ll(global.dom_grid);
   
   timestep = 0.1;
#else /* HYDRO_TEST==7 */
   
   
#ifdef INIT_DMHYDRO
   /*===================================================================== 
    * initialize hydro-variables accounting for baryons
    * (i.e. put a baryon down at every DM particle position...)
    *=====================================================================*/
   init_DMhydro();
#endif
   
   
#ifdef NO_DM
   NULL_ll(global.dom_grid);
#endif /* NO_DM */
#endif /* HYDRO_TEST==7 */
#endif /* HYDRO */
   
   

   /*================================================================
    * assign particles to the domain grid with simu.NGRID_DOM^3 nodes 
    *================================================================*/
   zero_dens(global.dom_grid);
   assign_dens(global.dom_grid);
#ifdef HYDRO
   add_gasdens(global.dom_grid);
#endif

#if (HYDRO_TEST==7)
   solve_gravity(grid_list, global.domgrid_no);
   calc_forces(global.dom_grid);
#endif

#ifdef POWERSPECTRUM
   /*================================================================
    * get initial power spectrum 
    *================================================================*/
   PkSpectrum.dump_Pk = TRUE;
   solve_dom_gravity(grid_list);

   PkSpectrum.Dgrowth_ini = calc_growth(global.a);
   PkSpectrum.Pk_ini      = PkSpectrum.Pk[PKMODE];
#endif
   
   /*================================================================
    * initialize some counters 
    *================================================================*/
   no_timestep         = no_first_timestep+1;  /* count total number of integration steps */
   global.total_time   = 0.;                   /* cumulative total time for simulation    */
   global.output_count = 0;                    /* count the number of outputs             */
   
   /* make *current* time step available to AHF/etc. routines */
   global.no_timestep = no_first_timestep;
   
   /*========================================================================================= 
    * all the following features of AMIGA does not require time-stepping...
    *=========================================================================================*/
#if (AHFstep || FORCETEST)
   global.ioflag    = TRUE;
   global.fst_cycle = TRUE;

   step(&grid_list, &no_grids, timestep, timecounter);
   
   update_logfile(timecounter, timestep, no_timestep);
   
   goto exit_AMIGA;
#endif



#	ifdef WITH_MPI
	/* FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME     *\
	 *       This is just a marker. Everything below here has not      *
	 *       been checked to work in MPI mode, hence noone will        *
	 *       do it. And we stop.                                       *
	\* FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME FIXME     */

	/* Ending the run */
	stoprun();

	/* And exit. */
	common_terminate(EXIT_SUCCESS);
#	endif
   
      
   /*========================================================================================= 
    *                                     TIME-STEP LOOP 
    *=========================================================================================*/
   /**********************************************************************************************************/
   fprintf(stderr,"\n\n************************************************************************************\n");
   fprintf(stderr,"*                          finished setting up simulation                          *\n");
   fprintf(stderr,"*                            STARTING TIME STEPPING LOOP!                          *\n");
   fprintf(stderr,"************************************************************************************\n");
   /**********************************************************************************************************/  
   
   fprintf(io.logfile,"\n\n");
   fprintf(io.logfile,"**********************************************************************************\n");
   fprintf(io.logfile,"*                        STARTING MAIN TIME STEPPING LOOP                        *\n");
   fprintf(io.logfile,"**********************************************************************************\n\n\n");
   fflush(io.logfile);
   
   /*======================================================================
    * step until final step is reached or terminate signal is encountered
    *======================================================================*/
#ifdef REVERSE_INTEGRATION
   while((global.a-simu.a_final) > ZERO && global.terminate == FALSE)
#else
   while((simu.a_final-global.a) > ZERO && global.terminate == FALSE)
#endif
     {
	 
      /* make no_timestep globally available */
      global.no_timestep = no_timestep;
      
      /* each grid has its own timestep */
      (global.dom_grid)->timestep = timestep;
      
      /* reset speeding flags */
      global.speeding    = FALSE;
      global.no_speeders = 0;
      global.no_ssteps   = 0;
      
      /* reset leaver counters */
      global.no_leavers = 0;
      global.no_lsteps  = 0;
      
      /* calling step from here means: first cycle starts */
      global.fst_cycle  = TRUE;
      
      /*============================================================*/
      step(&grid_list, &no_grids, timestep, timecounter);
      /*============================================================*/
      
#if (HYDRO_TEST==7)
      fprintf(stderr,"step %8d done:    time=%12.8g->%12.8g  timestep=%12.8g   cfl_timestep=%12.8g   cfl_speed=%12.8g   z=%g\n",
              no_timestep, timecounter-timestep, timecounter+timestep, timestep, 
              global.a*global.dom_grid->spacing/global.cfl_speed,global.cfl_speed,1./global.a-1.);      
#endif
#ifdef REVERSE_INTEGRATION
      /* decrement integration variable */
      timecounter -= timestep;
#else
      /* increment integration variable */
      timecounter += timestep;
#endif

      /* make integration variable globally known */
      global.super_t = timecounter;
      global.a       = calc_super_a(timecounter);
      global.t       = calc_t(global.a);
      global.z       = 1./global.a - 1.;
            
#ifdef ADAPTIVE
      /* rehash linked-list's properly ... don't rely on relink_back() too much */
       NULL_ll(global.dom_grid);
       ll(global.no_part, global.fst_part, global.dom_grid);
#endif
      
write_log:
      /*==========================================================================
       * dump information about energy, timing, sizes of grids, etc. into logfile
       *==========================================================================*/
      update_logfile(timecounter, timestep, no_timestep);
      
#ifdef AHFstep
      goto exit_AMIGA;
#endif
      
#ifdef VERBOSE
      fprintf(stderr,"\n\nAMIGAmain:      STEP %d done\n",
              no_timestep);
      fprintf(stderr,"                -------\n");
      fprintf(stderr,"\n\n");
#endif
     
      /*============================================================================
       *                           manage output files
       *============================================================================*/
      manage_outputs(timecounter, timestep);
      
      
      
      /*============================================================================
       *                     unscheduled termination of AMIGA?
       *============================================================================*/
      if((termfile = fopen(global.termfile_name,"r")) != NULL)
        {
         /* delete that file */
         fclose(termfile);
         sprintf(command,"rm -f %s",global.termfile_name);
         system(command);
         
         /* terminate AMIGA */
         global.terminate = TRUE;
        }
      
      /*============================================================================
       *                          adaptive timestepping
       *============================================================================*/
#ifdef ATS
      timestep = adjust_timestep(timecounter, timestep);
#endif /* ATS */
      
      
      
      /* now it's safe to reset the maxima values */
      global.max_dr2 = 0.0;
      global.max_dp2 = 0.0;
      
#ifdef REVERSE_INTEGRATION
      /* do not move: no_timestep REALLY needs updating HERE */
      no_timestep--;
#else
      /* do not move: no_timestep REALLY needs updating HERE */
      no_timestep++;
#endif

#ifdef LIGHTCONE 
      /* The operations below could be done from the output_lc(),
         * but doing them here guarantees the possibility to start the
         * lightcone at any time. */
      /* 1) calculate and save the lightcone radius */
      io.rcone0=r_cone(global.a);
      /* 2) save the present (old) redshift       */
      io.z0=global.z;	 
      /* 3) store the backup particle coordinates */
      store_backup();
#endif /* LIGHTCONE */      

      
      fprintf(io.logfile, 
              "-----------------------------------------------------------------------\n\n");
      fflush(io.logfile);
     }
   /*===============================================================================
    *                       ...END OF TIME-STEPPING LOOP
    *===============================================================================*/
   
   

   
   /*============================================================================
    *                        analyse final configuration
    *============================================================================*/
   /* definitely save last step done */
   no_timestep -= 1;
  
   fprintf(io.logfile," saving last step: %d (z=%8.3g)\n\n",
           no_timestep, 1.0/calc_super_a(timecounter)-1.0);
   output(io.outfile_prefix, 0);
#ifdef HYDRO
   output_grid(global.dom_grid, 0);
#endif
   global.ioflag = TRUE;
#ifdef LIGHTCONE 
   /* the lightcone data has been written already;
   * we close the file if it exists. */
   if(io.lightcone!=NULL) {
      fclose(io.lightcone);
      io.lightcone=NULL;
   }
#endif /* LIGHTCONE */
   
   /* ...do final AHF analysis */
#ifdef AHF
   step(&grid_list, &no_grids, timestep, timecounter);
#endif
   
#ifdef POWERSPECTRUM
   /* dump final P(k) */
   zero_dens(global.dom_grid);
   assign_dens(global.dom_grid);
   
   step(&grid_list, &no_grids, timestep, timecounter);
   
   /* call to power_spectrum() requires float *data rather than global.dom_grid... */
#endif
   
   /* for -DAHFstep and -DBDMstep the grids have already been free'd */
   for_grid = global.dom_grid;
   no_grids = global.domgrid_no+1;
   while(no_grids > 0)
     {
      free_grid(for_grid, &no_grids);
      for_grid--;
     }
   
   
exit_AMIGA:
   /*============================================================================
    *                                   BYE BYE!
    *============================================================================*/
      
   /* free all allocated memory... */
   free(grid_list);

   if(io.a_out != NULL)
      free(io.a_out);
   
   free(io.icfile_name);
   free(io.dumpfile_name);
   free(io.logfile_name);
   free(io.outfile_prefix);
   free(global.termfile_name);
   
   free(global.fst_part);
   if(global.fst_gas)
     free(global.fst_gas);
   if(global.fst_star)
     free(global.fst_star);
  
#ifdef POWERSPECTRUM
   free(PkSpectrum.rk);
   free(PkSpectrum.Pk);
#endif
   
   fprintf(io.logfile,
           "==========================================================\n");
   fprintf(io.logfile,
           "                       FINISHED (v%3.1f/%d)\n",VERSION,BUILT);
   fprintf(io.logfile,
           "==========================================================\n");
   fclose(io.logfile);

#	ifdef WITH_MPI
	/* Gracefully terminate MPI */
	MPI_Finalize();
#	endif
   
   return EXIT_SUCCESS;
}

#ifdef WITH_MPI
static void
local_communicate_all(void)
{
	if (global_mpi.rank == 0) {
		/*********************\
		 *      SENDING      *
		\*********************/
		/** This is where to send to */
		int target;

		for (target = 1; target < global_mpi.size; target++) {
			MPI_Send(&global, sizeof(struct info_global),
			         MPI_BYTE,
			         target, 0, MPI_COMM_WORLD);
		}
		/*********************\
		 *    DONE SENDING   *
		\*********************/
	} else {
		/*********************\
		 *     RECIEVING     *
		\*********************/
		/** The status */
		MPI_Status stat;

		MPI_Recv(&global, sizeof(struct info_global),
		         MPI_BYTE,
		         0, 0, MPI_COMM_WORLD, &stat);

		/*********************\
		 *  DONE  RECIEVING  *
		\*********************/
	}

	return;
}
#endif
