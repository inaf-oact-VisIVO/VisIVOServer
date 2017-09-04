#include "define.h"

#ifndef NEWSTARTRUN

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* the important definitions have to be included first */
#include "common.h"
#include "param.h"
#include "tdef.h"

/* ...and now for the actual includes */
#include "libutility/utility.h"
#include "libio_serial/io_serial.h"

#ifdef WITH_MPI
	/* When in MPI mode, don't try to read from stdin */
#	define WITH_AMIGA_input
#endif

void get_io_filenames   (FILE *startrun_dat);

/*==============================================================================
* set all parameters for starting a new AMIGA run 
*==============================================================================*/
void startrun(char AMIGA_input[MAXSTRING], double *timecounter, double *timestep, int *no_first_timestep)
{
  FILE   *startrun_dat;
  double  dtemp, dtemp2;
  int     i,j,k;
  
  uparam user_data;   // structure holding all relevant parameters provided by the user
  
  /* open AMIGA.input */
  if((startrun_dat = fopen(AMIGA_input,"r")) == NULL)
    {
      fprintf(stderr,"\n\nstartrun:  could not open file parameterfile %s\n",AMIGA_input);
      exit(1);
    }
  
  /*==================================================================================
   * read all relevant parameters and filenames from the user-provided parameter file
   *==================================================================================*/
  get_user_data(startrun_dat, &user_data);
  
  /*========================================================================================
   * copy information over from "user structure" to the globally available AMIGA structures
   *========================================================================================*/
  init_io_structure(user_data);
  
  /*========================================================================================
   * start logging information
   *========================================================================================*/
  init_logfile(user_data);
    
  
  /*========================================================================================
   * read dark matter particles from file
   * (input simply fills the io. structure and nothing else...)
   *========================================================================================*/
  input(io.icfile_name);  

  
  
  
  /****************************************************************************************
   *           WE ARE NOW GOING TO INITIALIZE ALL OTHER GLOBAL STRUCTURES
   *           ----------------------------------------------------------
   ****************************************************************************************/
  
  /*====================================================
   * override parameter values with actual values
   *====================================================*/
  io.header.paramNSTEPS       = NSTEPS;
  io.header.paramNGRID_DOM    = user_data.NGRID_DOM;
  io.header.paramNth_dom      = user_data.Nth_dom;
  io.header.paramNth_ref      = user_data.Nth_ref;
  
  io.header.paramE_UPDATE     = E_UPDATE;
  io.header.paramCELLFRAC_MAX = CELLFRAC_MAX;
  io.header.paramCELLFRAC_MIN = CELLFRAC_MIN;
  io.header.paramCA_CRIT      = CA_CRIT;
  io.header.paramMAX_L1DIM    = MAX_L1DIM;
#ifdef FFT
  io.header.paramDOMSWEEPS    = -1;
#else
  io.header.paramDOMSWEEPS    = DOMSWEEPS;
#endif
  io.header.paramREFSWEEPS    = REFSWEEPS;
  
  io.header.paramAHF_MINPART  = AHF_MINPART;
  io.header.paramAHF_VTUNE    = AHF_VTUNE;
  io.header.paramAHF_RISE     = AHF_RISE;
  io.header.paramAHF_SLOPE    = AHF_SLOPE;
  io.header.paramAHF_MAXNRISE = AHF_MAXNRISE;

  
  /*====================================================
   * simu.
   *====================================================*/
  init_simu_structure(user_data);
  
  /*====================================================
   * global.
   *====================================================*/
  init_global_structure(user_data);
  
  
  /*====================================================
   * energy.
   *====================================================*/
  energy.K_initial   = io.header.K_initial;
  energy.U_initial   = io.header.U_initial;
  energy.K_current   = io.header.K_current;
  energy.U_current   = io.header.U_current;
  energy.integral    = io.header.Eintegral;
  energy.econst      = io.header.Econst;
  energy.aold        = io.header.a_current;
  
  /*=====================================================================*
   * PkSpectrum.
   *=====================================================================*/
  PkSpectrum.rk = (double *) calloc(simu.NGRID_DOM, sizeof(double));
  PkSpectrum.Pk = (double *) calloc(simu.NGRID_DOM, sizeof(double));
  PkSpectrum.Dgrowth_ini = -1.0;
  PkSpectrum.Dgrowth_now = -1.0;
  PkSpectrum.Pk_ini      = -1.0;
  PkSpectrum.Pk_now      = -1.0;
  PkSpectrum.dump_Pk     = FALSE;
  PkSpectrum.time        = 0;
  /* -> "now" values will be calculated in get_Pk() prior to check in main() */
  
  
#ifdef ENERGY
  if(io.header.no_timestep > 0)
    energy.echeck = (energy.econst/(energy.U_current));
#endif
  
  /*======================================================================
   * initialize timecounter and timestep
   * (here we explicily assume integration in supercomoving coordinates!)
   *======================================================================*/
  *timecounter       = global.super_t;
  *no_first_timestep = global.no_timestep;
  *timestep          = init_timestep(user_data);
  
  /* do >>not<< allow for a larger timestep than found in the file! */
  if(*timestep > io.header.timestep && fabs(io.header.timestep) > ZERO)
    *timestep = io.header.timestep;
  
  
#ifdef ISOLATED
  /*==============================================================
   * calculate number of output files manually (dt = const!)
   *==============================================================*/
  dtemp = (simu.t_final-simu.t_initial)/(double)user.no_outputs;
  for(i=0; i<user.no_outputs; i++)
    {
      dtemp2      = simu.t_initial + (double)i * dtemp;
      io.a_out[i] = calc_a(dtemp2);
    }
#endif
    
  
  /****************************************************************************************
   *           WE ARE NOW GOING TO REMOVE SOME PARTICLES, IF NEEDED
   *           ----------------------------------------------------
   *
   *  simu.no_part    => original number of particles
   *  global.no_part  =>  actual  number of particles
   ****************************************************************************************/
  
  
#ifdef REMOVE_PARTICLES
  /*=====================================================================*
   * use only particles in a region about some X/Y/Z
   *=====================================================================*/
  {
    double        Xc, Yc, Zc, Rvir, Xrange, Yrange, Zrange, dX, dY, dZ;
    
    long unsigned no_part;
    partptr       fst_part, cur_part, new_part;
    
    
    Xc   = 3.56593968   / simu.boxsize;
    Yc   = 12.4243316   / simu.boxsize;
    Zc   = 0.17760938   / simu.boxsize;
    Rvir = 0.55910      / simu.boxsize;
    Xrange = 2.5;
    Yrange = 3.5;
    Zrange = 4.2;
    
    
    no_part  = 0;
    for(cur_part=global.fst_part; cur_part<(global.fst_part+simu.no_part); cur_part++)
      {
        dX = fabs(cur_part->pos[X]-Xc);
        dY = fabs(cur_part->pos[Y]-Yc);
        dZ = fabs(cur_part->pos[Z]-Zc);
        
        /* take care of periodic boundary conditions! */
        if(dX > 0.5) dX -= 1.0;
        if(dY > 0.5) dY -= 1.0;
        if(dZ > 0.5) dZ -= 1.0;
        
        if(dX < Xrange*Rvir && dY < Yrange*Rvir && dZ < Zrange*Rvir)
          no_part++;
      }
    
    /* allocate memory for particles */
    fst_part = c_part(no_part);
    
    /* shift particles to new array */
    new_part = fst_part;
    for(cur_part=global.fst_part; cur_part<(global.fst_part+simu.no_part); cur_part++)
      {
        
        dX = fabs(cur_part->pos[X]-Xc);
        dY = fabs(cur_part->pos[Y]-Yc);
        dZ = fabs(cur_part->pos[Z]-Zc);
        
        /* take care of periodic boundary conditions! */
        if(dX > 0.5) dX -= 1.0;
        if(dY > 0.5) dY -= 1.0;
        if(dZ > 0.5) dZ -= 1.0;
        
        if(dX < Xrange*Rvir && dY < Yrange*Rvir && dZ < Zrange*Rvir)
          {
            new_part->pos[X] = cur_part->pos[X];
            new_part->pos[Y] = cur_part->pos[Y];
            new_part->pos[Z] = cur_part->pos[Z];
            new_part->mom[X] = cur_part->mom[X];
            new_part->mom[Y] = cur_part->mom[Y];
            new_part->mom[Z] = cur_part->mom[Z];
#ifdef MULTIMASS
            new_part->weight = cur_part->weight;
#endif
            new_part++;
          }
      }
    
    /* erase old particle list and store new one */
    free(global.fst_part);
    global.fst_part = fst_part;
    
    /* update global.no_part parameter */
    global.no_part  = no_part;
  }
#endif
  
  /*=====================================================================*
   * remove all particles with cur_weight != 1
   *
   *              USE AHFmmfocus AT YOUR OWN RISK!
   * 
   * or in other words, make sure that you are removing those particles
   * you actually intend to remove by tailoring the criterion!!!!!!
   *
   * NOTE: the criterion has to be specified two times:
   *        1. when counting the number of particles to keep
   *        2. when actually throwing away the particles
   *=====================================================================*/
#ifdef MULTIMASS
#ifdef AHFmmfocus
  fprintf(stderr,"\n==================================================================\n");
  fprintf(stderr,"                          AHFmmfocus\n");
  fprintf(stderr,"               ? ARE YOU SURE ABOUT THIS FLAG ?\n");
  fprintf(stderr,"==================================================================\n");
  fprintf(stderr,"AMIGA will now remove all particles whose mass is not %g Msun/h\n",io.header.pmass);
  fprintf(stderr," simu.no_part = %12ld -> global.no_part = ",simu.no_part);
  fprintf(io.logfile,"\n==================================================================\n");
  fprintf(io.logfile,"                          AHFmmfocus\n");
  fprintf(io.logfile,"               ? ARE YOU SURE ABOUT THIS FLAG ?\n");
  fprintf(io.logfile,"==================================================================\n");
  fprintf(io.logfile,"AMIGA will now remove all particles whose mass is not %g Msun/h\n",io.header.pmass);
  fprintf(io.logfile," simu.no_part = %12ld -> global.no_part = ",simu.no_part);
  fflush(io.logfile);
  if(simu.no_species > 1)
    {
      long unsigned no_part;
      partptr       fst_part, cur_part, new_part;
      
      /* 1. count number of particles to keep */
      no_part  = 0;
      for(cur_part=global.fst_part; cur_part<(global.fst_part+simu.no_part); cur_part++)
        {
         /* we want to keep only the most common particles...
         * ...and the weight of the most common particle is 1 */
#if (defined MLAPM || defined ART)
         if(fabs(cur_part->weight-global.fst_part->weight) < ZERO)
#else
         if(fabs(cur_part->weight-1.0) < ZERO)
#endif
            no_part++;
        }
      
      /* allocate memory for particles */
      fst_part = c_part(no_part);
      
      
      /* 2. remove all other particles */
      new_part = fst_part;
      for(cur_part=global.fst_part; cur_part<(global.fst_part+simu.no_part); cur_part++)
        {
          /* we want to keep only the most common particles...
           * ...and the weight of the most common particle is 1 */
#if (defined MLAPM || defined ART)
         if(fabs(cur_part->weight-global.fst_part->weight) < ZERO)
#else
         if(fabs(cur_part->weight-1.0) < ZERO)
#endif
              {
              new_part->pos[X] = cur_part->pos[X];
              new_part->pos[Y] = cur_part->pos[Y];
              new_part->pos[Z] = cur_part->pos[Z];
              new_part->mom[X] = cur_part->mom[X];
              new_part->mom[Y] = cur_part->mom[Y];
              new_part->mom[Z] = cur_part->mom[Z];
              new_part->weight = cur_part->weight;
              new_part++;
            }
        }
      
      /* erase old particle list and store new one */
      free(global.fst_part);
      global.fst_part = fst_part;
      
      /* update global.no_part parameter */
      global.no_part  = no_part;
    }
  fprintf(stderr,"%12ld\n",global.no_part);
  fprintf(io.logfile,"%12ld\n",global.no_part);
  fflush(io.logfile);
  
#endif
#endif
  
  /*===============================================
   * the unit matrix is used with the hydro solver
   *===============================================*/
  for(i=0; i<NDIM; i++)
    {
      unit_matrix[i][X] = 0;
      unit_matrix[i][Y] = 0;
      unit_matrix[i][Z] = 0;
    }
  unit_matrix[X][X] = 1;
  unit_matrix[Y][Y] = 1;
  unit_matrix[Z][Z] = 1;
  
  
  /*=====================================================================
   * write information to logfile as well as to parameter file
   *=====================================================================*/
  log_structures();
  
  write_parameterfile();
}

#else /* NEWSTARTRUN */

#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#if (NONC99 == 0)
#	include <math.h>
#else
#	include <replace_math.h>
#endif

#include "common.h"
#include "libutility/utility.h"
#include "libio/io.h"

/**
 * \brief Helper function for reading the parameters.
 */
static void
local_startrunParams(char *paramfile);

/**
 * \brief Helper function for setting up the logging.
 */
static void
local_startrunLog(void);

/**
 * \brief Helper function for opening and initializing the IC file.
 */
static void
local_startrunFopen(void);

/**
 * \brief Helper function for reading the IC file.
 */
static void
local_startrunRead(void);

/**
 * \brief Helper function for setting the simulation parameters.
 */
static void
local_startrunSimparams();

/**
 * \brief Helper function for setting the return values.
 */
static void
local_startrunRetset(double *timecounter,
                     double *timestep,
                     int32_t *no_first_timestep);

extern void
startrun(char *paramfile,
         double *timecounter,
         double *timestep,
         int32_t *no_first_timestep)
{
	/* Read the parameters */
	local_startrunParams(paramfile);

	/* Now set up the logging */
	local_startrunLog();
	io_logging_part(global_io.log, "Setting up the run");
   
	/* FIXME Do some glueing FIXME */
	io.logfile = global_io.log->logfile;

	/* Open the file */
	local_startrunFopen();

	/* Set global_info */
	global_info.fst_part = NULL;
	global_info.no_part = UINT64_C(0);
#	ifdef WITH_MPI
	global_info.loadbal = loadbalance_new(global_io.log,
	                                      LOADBALANCE_EQUALPART,
	                                      SFC_CURVE_HILBERT,
	                                      LOADBALANCE_DOMAIN_LEVEL,
	                                      global_mpi.size);
#	endif

	/* Now read the initial conditions (sets the rest of global_info) */
	local_startrunRead();

	/* Set the simulation parameters */
	local_startrunSimparams();
   
	/* Now set the returned simulation counter stuff thingies */
	local_startrunRetset(timecounter, timestep, no_first_timestep);

	/* Set global time counter */
	global.super_t = *timecounter;                      
	global.a       = calc_super_a(global.super_t);    
	global.t       = calc_t(global.a);            
	global.z       = 1./global.a - 1.;
  
	/* And now lets be so nice and close the file... */
	io_logging_section(global_io.log, "Tidying");
	io_file_close(global_io.log, &(global_io.file));

	/* FIXME Do some glueing FIXME */

	return;
}

extern void
stoprun(void)
{	
	io_parameter_del(&(global_io.params));
	io_logging_stop(&(global_io.log));

	return;
}

static void
local_startrunParams(char *paramfile)
{
#	ifdef WITH_MPI
	global_io.params =
	    io_parameter_get(IO_PARAMETER_FROM_FNAME, paramfile);
#	else
#		ifdef WITH_AMIGA_input
	global_io.params = 
	    io_parameter_get(IO_PARAMETER_FROM_FNAME, paramfile);
#		else
	global_io.params = 
	    io_parameter_get(IO_PARAMETER_FROM_STDIN, NULL);
#		endif
#	endif
	if (global_io.params == NULL) {
		common_terminate(EXIT_FAILURE);
	}

	return;
}

static void
local_startrunLog(void)
{
	global_io.log = 
	    io_logging_start(global_io.params->outfile_prefix,
		                 INT32_C(VERBOSITY),
	                     IO_LOGGING_FLAGS_DUPLICATE_CRITICAL);
	if (global_io.log == NULL) {
		io_parameter_del(&(global_io.params));
		common_terminate(EXIT_FAILURE);
	}

	/* And since we can log now, log a bit. */
	io_logging_hello(global_io.log, "0.1");
	io_logging_msg(global_io.log, INT32_C(2), "Used filenames:");
	io_logging_msg(global_io.log, INT32_C(2),
	               "IC file  = %s (%s)",
	               global_io.params->icfile_name,
	               io_file_typestr(global_io.params->ic_filetype));
	io_logging_msg(global_io.log, INT32_C(2),
	               "log file = %s",
	               global_io.log->fname);

	return;
}

static void
local_startrunFopen(void)
{
	io_logging_section(global_io.log, "Opening the data file");

	global_io.file = io_file_open(global_io.log,
	                              global_io.params->icfile_name,
	                              global_io.params->ic_filetype,
#	ifdef BYTESWAP
	                              IO_FILE_IS_SWAPPED,
#	else
	                              IO_FILE_UNKOWN_SWAPPING,
#	endif

	                              IO_FILE_READ,
	                              global_io.params->reader);
	if (global_io.file == NULL)
		common_terminate(EXIT_FAILURE);

	/* Set the scaling */
	if (global_io.file->ftype == IO_FILE_GADGET)
		io_gadget_resetscale(global_io.log,
		                     (io_gadget_t)global_io.file,
		                     (double)GADGET_LUNIT,
		                     (double)GADGET_MUNIT);
	if (global_io.file->ftype == IO_FILE_MGADGET)
		io_mgadget_resetscale(global_io.log,
		                      (io_mgadget_t)global_io.file,
		                      (double)GADGET_LUNIT,
		                      (double)GADGET_MUNIT);
  
	/* Init the file */
	io_file_init(global_io.log, global_io.file);

  
	/* Now dump the file information */
	io_file_log(global_io.log, global_io.file);

	return;
}

static void
local_startrunRead(void)
{
	io_file_strg_struct_t strg;
	uint64_t pskip = UINT64_C(0);
	uint64_t pread = UINT64_MAX;
	partptr tmppart;

	io_logging_section(global_io.log, "Reading data from file");

	/* See if we are supposed to read anything at all */
	if (global_io.file->ftype == IO_FILE_EMPTY) {
		global_info.no_part = 0;
		global_info.fst_part = NULL;
		return;
	}

	/* First create particle storage */
	io_logging_subsection(global_io.log, "Creating Storage");
	global_info.no_part = io_file_get_numpart(global_io.log,
	                                          global_io.file,
	                                          &pskip, &pread);
	global_info.fst_part = c_part((long)global_info.no_part);

	/* Create the description of the storage */
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
#	endif
	strg.id.val = &(global_info.fst_part->id);
	strg.id.stride =   (char *)&((global_info.fst_part+1)->id)
	                 - (char *)&(global_info.fst_part->id);
#	ifdef GAS_PARTICLES
	strg.u.val = &(global_info.fst_part->u);
	strg.u.stride =   (char *)&((global_info.fst_part+1)->u)
	                - (char *)&(global_info.fst_part->u);
#	else
	strg.u.val = NULL;
	strg.u.stride = (ptrdiff_t)0;
#	endif
	strg.bytes_float = sizeof(global_info.fst_part->pos[0]);
	strg.bytes_int = sizeof(global_info.fst_part->id);

#	ifdef VERBOSE
	/* Print the description */
	io_file_strg_log(global_io.log, strg);
#	endif


  /* Now read the particles */
	io_logging_subsection(global_io.log, "Reading");
	if (io_file_readpart(global_io.log, global_io.file,
	                     pskip, pread, strg) == UINT64_C(0) ) {
		/* We read 0 particles from the file, this is an error */
		common_terminate(EXIT_FAILURE);
	}

	/* Print the first two particles to the logfile */
	io_logging_subsection(global_io.log, "Short sanity check");
	tmppart = global_info.fst_part;
	io_logging_msg(global_io.log, INT32_C(5),
	               "First particle:");
	io_logging_msg(global_io.log, INT32_C(5),
	               "    positions (x,y,z):      %g  %g  %g",
	               tmppart->pos[0],
	               tmppart->pos[1],
	               tmppart->pos[2]);
	io_logging_msg(global_io.log, INT32_C(5),
	               "    velocities (vx,vy,vz):  %g  %g  %g",
	               tmppart->mom[0],
	               tmppart->mom[1],
	               tmppart->mom[2]);
#	ifdef MULTIMASS
	io_logging_msg(global_io.log, INT32_C(5),
	               "    weight:                 %g",
	               tmppart->weight);
#	endif
	io_logging_msg(global_io.log, INT32_C(5),
	               "    ID:                     %" PRIpartid,
	               tmppart->id);
#	ifdef GAS_PARTICLES
	io_logging_msg(global_io.log, INT32_C(5),
	               "    energy:                 %g",
	               tmppart->u);
#	endif
	tmppart = global_info.fst_part+global_info.no_part-1;
	io_logging_msg(global_io.log, INT32_C(5),
	               "Last particle:");
	io_logging_msg(global_io.log, INT32_C(5),
	               "    positions (x,y,z):      %g  %g  %g",
	               (tmppart)->pos[0],
	               (tmppart)->pos[1],
	               (tmppart)->pos[2]);
	io_logging_msg(global_io.log, INT32_C(5),
	               "    velocities (vx,vy,vz):  %g  %g  %g",
	               (tmppart)->mom[0],
	               (tmppart)->mom[1],
	               (tmppart)->mom[2]);
#	ifdef MULTIMASS
	io_logging_msg(global_io.log, INT32_C(5),
	               "    weight:                 %g",
	               (tmppart)->weight);
#	endif
	io_logging_msg(global_io.log, INT32_C(5),
	               "    ID:                     %" PRIpartid,
	               (tmppart)->id);
#	ifdef GAS_PARTICLES
	io_logging_msg(global_io.log, INT32_C(5),
	               "    energy:                 %g",
	               (tmppart)->u);
#	endif

#	ifdef MPI_DEBUG
	io_logging_subsection(global_io.log, "Longer sanity check");
	io_logging_msg(global_io.log, INT32_C(0),
	               "Fileobject after reading particles (supposedly "
	               "with correct multimass information now).");
	io_file_log(global_io.log, global_io.file);
#	endif

  
  return;
}

static void
local_startrunSimparams()
{
	io_logging_section(global_io.log, "Setting simulation parameter");

	io_logging_subsection(global_io.log, "Information from file");
#	ifdef WITH_MPI
	if (global_mpi.rank != 0) {
		io_logging_msg(global_io.log, INT32_C(4),
		               "Not setting up myself, will receive.");
      fflush(NULL);
	} else {
#	else
	{
#	endif
		int32_t no_timestep;

		io_file_get(global_io.log, global_io.file,
		            IO_FILE_GET_BOXSIZE, (void *)&(simu.boxsize));
		io_file_get(global_io.log, global_io.file,
		            IO_FILE_GET_OMEGA0, (void *)&(simu.omega0));
		io_file_get(global_io.log, global_io.file,
		            IO_FILE_GET_OMEGAL, (void *)&(simu.lambda0));
		io_file_get(global_io.log, global_io.file,
		            IO_FILE_GET_PMASS, (void *)&(simu.pmass));
		io_file_get(global_io.log, global_io.file,
		            IO_FILE_GET_NOPART, (void *)&(simu.no_part));
#		ifdef MULTIMASS
		io_file_get(global_io.log, global_io.file,
		            IO_FILE_GET_NOVPART, (void *)&(simu.no_vpart));
#		else
		simu.no_vpart = (double)(simu.no_part);
#		endif
		io_file_get(global_io.log, global_io.file,
		            IO_FILE_GET_NOSPECIES, (void *)&(simu.no_species));
		io_file_get(global_io.log, global_io.file,
		            IO_FILE_GET_AINITIAL, (void *)&(simu.a_initial));
		io_file_get(global_io.log, global_io.file,
		            IO_FILE_GET_DOUBLE, (void *)&(simu.double_precision));
		io_file_get(global_io.log, global_io.file,
		            IO_FILE_GET_MMASS, (void *)&(simu.multi_mass));
		io_file_get(global_io.log, global_io.file,
		            IO_FILE_GET_MINWEIGHT, (void *)&(simu.min_weight));
		io_file_get(global_io.log, global_io.file,
		            IO_FILE_GET_MAXWEIGHT, (void *)&(simu.max_weight));
	
		/* Copy over the information contained in the parameter file */
		simu.z_final   = global_io.params->final_z;
		simu.NGRID_DOM = global_io.params->NGRID_DOM;
#		ifndef FFT
		simu.NGRID_MIN = MIN_L1DIM;
#		else
		simu.NGRID_MIN = simu.NGRID_DOM;
#		endif
		simu.Nth_dom = global_io.params->Nth_dom;
		simu.Nth_ref = global_io.params->Nth_ref;
#		ifdef MOND
		simu.g0 = global_io.params->g0;
		simu.h0 = global_io.params->h0;
#		endif

		/* Set quantities given by constants */
		simu.NGRID_MAX = MAX_L1DIM;
#		ifdef NP_LIMIT
		simu.np_limit = TRUE;
#		else
		simu.np_limit = FALSE;
#		endif
		simu.mean_dens = (double) 1.0;

#		ifdef AHFmmfocus
		fprintf(stderr, "AHFmmfocus not supported with NEWSTARTRUN\n");
		common_terminate(1);
#		else
		simu.mmfocus  = 0;
#		endif
#		ifdef HYDRO
		fprintf(stderr, "HYDRO not supported with NEWSTARTRUN\n");
		common_terminate(1);
#		else
		simu.hydro    = 0;
#		endif
#		ifdef MHD
		fprintf(stderr, "MHD not supported with NEWSTARTRUN\n");
		common_terminate(1);
#		else
		simu.magneto  = 0;
#		endif
      
		/* Set the time unit */
		simu.t_unit = 1/H0; // we assume that we only ever deal with
		                    // cosmological simulations...
      
		/* Set derived quantities */
		simu.SHIFT     = ((double)0.5000000/(double) simu.NGRID_DOM);
		simu.z_initial = (double)1.0/simu.a_initial - (double)1.0;
		simu.a_final   = (double)1.0/((double)1.0 + simu.z_final);
		simu.FourPiG   = 1.5*simu.omega0;

		/* Do some sanity checks */
		io_file_get(global_io.log, global_io.file,
		            IO_FILE_GET_NOTSTEP, (void *)&(no_timestep));

		if ( isless(fabs(simu.a_initial-simu.a_final), ZERO) ) {
			io_logging_warn(global_io.log, INT32_C(3),
			                "Since a_initial = %g is equal to "
			                "a_final = %g, create_timeline will not "
			                "function correctly, setting "
			                "a_initial = .1 * a_final = %g",
			                simu.a_initial, simu.a_final,
			                simu.a_final / 10.0);
			simu.a_initial = simu.a_final / 10.0;
		}
	} /* End of stuff done solely by process 0 */

	io_logging_subsection(global_io.log,
	                      "Gathering from reading processes");
#	ifdef WITH_MPI
	io_logging_msg(global_io.log, INT32_C(4),
	               "Broadcast of simulation parameters!");
	MPI_Bcast(&simu, sizeof(struct param_simu),
	          MPI_BYTE,
	          0,
	          MPI_COMM_WORLD);
	io_logging_msg(global_io.log, INT32_C(4),
	               "Broadcast done.");
#	endif


	/* Create timeline */
	io_logging_subsection(global_io.log, "Local setup");
	io_logging_msg(global_io.log, INT32_C(2),
	               "Creating timeline from a = %g to a = %g",
	               simu.a_initial/10., simu.a_final);
	create_timeline(simu.a_initial/10., simu.a_final, &simu.timeline);
	io_logging_msg(global_io.log, INT32_C(2),
	               "Timeline created");

	/* Set the SFC information */
	io_logging_msg(global_io.log, INT32_C(2),
	               "Setting volume boundaries");
#	ifdef AHFrestart
	global_info.minkey = (sfc_key_t)(
	   ((io_ares_t)(global_io.file))->header->minkey);
	global_info.maxkey = (sfc_key_t)(
	   ((io_ares_t)(global_io.file))->header->maxkey);
	global_info.level = ((io_ares_t)(global_io.file))->header->lb_level;
#	else
	global_info.level = LOADBALANCE_DOMAIN_LEVEL;
	global_info.minkey = (sfc_key_t)0;
	global_info.maxkey = (sfc_key_t)((1<<(3*global_info.level))-1);
#	endif
	global_info.ctype = SFC_CURVE_HILBERT;
	io_logging_msg(global_io.log, INT32_C(2),
	               "  minkey: %"SFC_PRIkey, global_info.minkey);
	io_logging_msg(global_io.log, INT32_C(2),
	               "  maxkey: %"SFC_PRIkey, global_info.maxkey);
	io_logging_msg(global_io.log, INT32_C(2),
	               "  level : %i", global_info.level);
	io_logging_msg(global_io.log, INT32_C(2),
	               "  ctype : %s", sfc_curve_typestr(global_info.ctype));

	/* Now that we have the timeline, set the time variables */
	simu.super_t_initial = calc_super_t(simu.a_initial);
	simu.super_t_final   = calc_super_t(simu.a_final);
	simu.t_initial       = calc_t(simu.a_initial);
	simu.t_final         = calc_t(simu.a_final);


	/* FIXME
	 * Not set or not properly set simu-structure members:
	 *
	 * Hydro variables:
	 * gamma
	 * omegab
	 * omegaDM
	 * f_b
	 * H_frac
	 * T_init
	 * e_init
	 * med_weight
	 * l_unit
	 * m_unit
	 *
	 * AHF variable:
	 * no_halos
	 *
	 * Unkown:
	 * ifdef LIGHTCONE: z_lightcone
	 * ifdef GAS_PARTICLES: no_gas
	 * ifdef GADGET: no_stars
    *
	 */
	/* Will use dummy values */
	simu.gamma    = 0.0;
	simu.omegab   = 0.0;
	simu.omegaDM  = simu.omega0;
	simu.f_b      = 0.0;
	simu.H_frac   = 0.0;
	simu.T_init   = 0.0;
	simu.B_init   = 0.0;
	simu.e_init   = 0.0;
	simu.no_halos = 0;
	simu.med_weight = simu.max_weight; // TODO: this is very conservative yet leads to more credible halos
	simu.l_unit = 0.0;
	simu.m_unit = 0.0;
#	ifdef LIGHTCONE
	simu.z_lightcone = 0.0;
#	endif
	simu.no_gas   = 0;
	simu.no_stars = 0;

//#	ifdef VERBOSE
	/* Be so kind and write everything to the logfile */
	io_logging_subsection(global_io.log, "Used simulation parameters");
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.omega0          :  %g", simu.omega0);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.lambda0         :  %g", simu.lambda0);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.boxsize         :  %g", simu.boxsize);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.a_initial       :  %g", simu.a_initial);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.a_final         :  %g", simu.a_final);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.z_initial       :  %g", simu.z_initial);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.z_final         :  %g", simu.z_final);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.t_initial       :  %g", simu.t_initial);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.t_final         :  %g", simu.t_final);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.super_t_initial :  %g", simu.super_t_initial);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.super_t_final   :  %g", simu.super_t_final);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.mean_dens       :  %g", simu.mean_dens);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.FourPiG         :  %g", simu.FourPiG);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.pmass           :  %g", simu.pmass);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.t_unit          :  %g", simu.t_unit);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.gamma           :  %g", simu.gamma);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.omegab          :  %g", simu.omegab);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.omegaDM         :  %g", simu.omegaDM);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.f_b             :  %g", simu.f_b);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.H_frac          :  %g", simu.H_frac);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.T_init          :  %g", simu.T_init);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.B_init          :  %g", simu.B_init);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.e_init          :  %g", simu.e_init);
#		ifdef LIGHTCONE
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.z_lightcone     :  %g", simu.z_lightcone);
#		endif
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.timeline (ptr)  :  %p", (void*)&(simu.timeline));
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.no_part         :  %lu", simu.no_part);
#		ifdef GAS_PARTICLES
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.no_gas          :  %lu", simu.no_gas);
#		endif
#		ifdef GADGET
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.no_stars        :  %lu", simu.no_stars);
#		endif
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.no_vpart        :  %g", simu.no_vpart);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.no_species      :  %i", simu.no_species);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.no_halos        :  %lu", simu.no_halos);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.NGRID_DOM       :  %i", simu.NGRID_DOM);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.NGRID_MIN       :  %i", simu.NGRID_MIN);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.NGRID_MAX       :  %i", simu.NGRID_MAX);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.Nth_dom         :  %g", simu.Nth_dom);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.Nth_ref         :  %g", simu.Nth_ref);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.SHIFT           :  %g", simu.SHIFT);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.min_weight      :  %g", simu.min_weight);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.max_weight      :  %g", simu.max_weight);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.np_limit        :  %i", simu.np_limit);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.mmfocus         :  %i", simu.mmfocus);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.multi_mass      :  %i", simu.multi_mass);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.double_precision:  %i", simu.double_precision);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.hydro           :  %i", simu.hydro);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.magneto         :  %i", simu.magneto);
#		ifdef MOND
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.g0              :  %g", simu.g0);
	io_logging_msg(global_io.log, INT32_C(5),
	               "simu.h0              :  %g", simu.h0);
#		endif
//#	endif /* VERBOSE */

	return;
}

static void
local_startrunRetset(double *timecounter,
                     double *timestep,
                     int32_t *no_first_timestep)
{
	io_logging_subsection(global_io.log, "Setting time counter");

#	ifdef WITH_MPI
	if (global_mpi.rank == 0) {
#	else
	{
#	endif
		double a_current;

		io_file_get(global_io.log, global_io.file,
		            IO_FILE_GET_NOTSTEP, (void *)no_first_timestep);
		io_file_get(global_io.log, global_io.file,
		            IO_FILE_GET_TSTEP, (void *)timestep);
		io_file_get(global_io.log, global_io.file,
		            IO_FILE_GET_A, (void *)&a_current);
		*timecounter = calc_super_t(a_current);
	}

#	ifdef WITH_MPI
	MPI_Bcast(timecounter, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(timestep, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(no_first_timestep, 1, MPI_INT, 0, MPI_COMM_WORLD);
#	endif

	return;
}

#endif /* NEWSTARTRUN */
