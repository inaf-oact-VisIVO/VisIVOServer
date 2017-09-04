#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "io_serial.h"
#include "../libutility/utility.h"

/*======================================================================
 *  structure holdig information about timing on each grid level
 * (will be calculated from information stored in each grid structure!)
 *======================================================================*/
struct  total_timing
{
   time_t potential;             /* deriving the potenital by GS sweeps      */
   time_t density;               /* deriving the proper densities            */
   time_t DK;                    /* drifting and kicking particles           */
   time_t grid;                  /* everything related to grid hierarchy     */
   time_t hydro;                 /* time spent for hydro-solver              */
   time_t ahf;
}grid_timing;


/*==============================================================================
 * here we decide what to write into the logfile to keep the user up-to-date...
 *==============================================================================*/
void update_logfile(double timecounter, double timestep, int no_timestep)
{
   
   gridls *for_grid;          /* for looping over all grids          */
   double  RAM;               /* how much RAM has been used                    */
   double  PDgrowth;
   long    dom_nopart;        /* number of particles linked to domain grid     */
   double  speeders;          /* ratio of speeding particles to all particles  */
   double  leavers;           /* ratio of leaver   particles to all particles  */
   double  da_a;
   double  total_time_step;
   
   /* write current timestep into logfile */
   fprintf(io.logfile, "\n");
   fprintf(io.logfile, "step no.        = %d  (dT = %g)\n", no_timestep,timestep);
   fprintf(io.logfile, "supercomoving T = %f -> %f\n", timecounter-timestep,timecounter);
   fprintf(io.logfile, "scale factor  a = %f -> %f (%f)\n",  calc_super_a(timecounter-timestep),calc_super_a(timecounter), 2*(calc_super_a(timecounter)-calc_super_a(timecounter-timestep))/(calc_super_a(timecounter-timestep)+calc_super_a(timecounter)));
   fprintf(io.logfile, "redshift      z = %f -> %f\n", 1.0/calc_super_a(timecounter-timestep)-1.0,1.0/calc_super_a(timecounter)-1.0);
   fflush(io.logfile);
   
   
#ifdef MOND
   fprintf(io.logfile, "maxima       dr  = %g,  dp  = %g\n", sqrt(global.max_dr2), sqrt(global.max_dp2)*simu.boxsize/simu.t_unit/global.a);
   fprintf(io.logfile, "MOND maxima  gN  = %g,  gM  = %g\n", global.max_gN, global.max_gM);
   fprintf(io.logfile, "MOND means  <gN> = %g, <gM> = %g\n", global.mean_gN/(double)global.steps, global.mean_gM/(double)global.steps);
   fprintf(io.logfile, "MOND events     = %g %%\n", 100.*global.no_MONDevents/global.no_ALLevents);
   global.steps         = 0;
   global.max_gN        = 0.0;
   global.max_gM        = 0.0;
   global.mean_gN       = 0.0;
   global.mean_gM       = 0.0;
   global.no_ALLevents  = 0;
   global.no_MONDevents = 0;
#else
   fprintf(io.logfile, "maxima       dr = %g, dp = %g\n", sqrt(global.max_dr2), sqrt(global.max_dp2)*simu.boxsize/simu.t_unit/global.a);
#endif
#ifdef ATS
   speeders = (double)global.no_speeders/(double)global.no_part*100.;
   fprintf(io.logfile,"speeders        = %5.2f %%  (%ld %ld)\n", speeders, global.no_speeders,global.no_ssteps);
#endif
   
   leavers = (double)global.no_leavers/(double)global.no_part*100.;
   fprintf(io.logfile,"leavers         = %5.2f %%  (%ld %ld)\n", leavers,global.no_leavers,global.no_lsteps);
   
#ifdef ENERGY
   fprintf(io.logfile, "energy check  C = %g, K = %e, U = %e\n", energy.econst, energy.K_current, energy.U_current);
   fprintf(io.logfile, "accuracy        = %5.2f %%\n", 100.*energy.echeck);
#endif
   
   /* reset time counter */
   grid_timing.potential = 0;
   grid_timing.density   = 0;
   grid_timing.DK        = 0;
   grid_timing.grid      = 0;
   grid_timing.hydro     = 0;
   
#ifdef POWERSPECTRUM
   /* check for linearity of fundamental mode */
   PDgrowth = (PkSpectrum.Pk_now/PkSpectrum.Pk_ini)/(pow2(PkSpectrum.Dgrowth_now/PkSpectrum.Dgrowth_ini));
   fprintf(io.logfile,"Pk mode #%1d  P/D = %g (%ld sec.)\n", PKMODE, PDgrowth, PkSpectrum.time);
   
   if(PDgrowth > PKGROWLIMIT)
      global.terminate = TRUE;
#endif
   
   
#ifdef HYDRO
   fprintf(io.logfile,"CFL speed         = %g\n", global.cfl_speed);
   fprintf(io.logfile,"CFL timestep      = %g\n", CFL_TUNE * calc_super_a(timecounter-timestep)*global.dom_grid->spacing/global.cfl_speed);
#endif
   
   /* get number of particles attached to domain grid */
   if(global.fin_l1dim > global.dom_grid->l1dim)
     {
      dom_nopart = global.no_part;
      for(for_grid=(global.dom_grid+1); for_grid->l1dim<global.fin_l1dim; for_grid++)
         dom_nopart -= for_grid->size.no_part;
      dom_nopart -= for_grid->size.no_part;
      
      global.dom_grid->size.no_part = dom_nopart;
     }
   
   /* start accumulating RAM used during this step */
   RAM = 0.0;

   
   /* write grid information to logfile */
   for(for_grid = global.dom_grid; for_grid->l1dim < global.fin_l1dim; for_grid++)
     {
      fprintf(io.logfile,"grid %6ld:  %10ld %10ld (%6ld %6ld %6ld %6ld %6ld) %ld %g\n",
              for_grid->l1dim,
              for_grid->size.no_nodes,
              for_grid->size.no_part,
              for_grid->time.potential,
              for_grid->time.density,
              for_grid->time.DK,
              for_grid->time.grid,
              for_grid->time.hydro,
              for_grid->no_sweeps,
              for_grid->cur_resid);
      
      grid_timing.potential += for_grid->time.potential;
      grid_timing.density   += for_grid->time.density;
      grid_timing.DK        += for_grid->time.DK;
      grid_timing.grid      += for_grid->time.grid;
      grid_timing.hydro     += for_grid->time.hydro;
      RAM                   += for_grid->size.no_nodes * global.bytes_node * bytes2GB;
      
      for_grid->time.potential = 0;
      for_grid->time.density   = 0;
      for_grid->time.DK        = 0;
      for_grid->time.grid      = 0;
      for_grid->time.hydro     = 0;
     }
   
   /* treat grid with for_grid->l1dim==global.fin_l1dim separately
    * in order to avoid incrementing for_grid++ too far!!! */

   /* write grid information to logfile */
   fprintf(io.logfile,"grid %6ld:  %10ld %10ld (%6ld %6ld %6ld %6ld %6ld) %ld %g\n",
           for_grid->l1dim,
           for_grid->size.no_nodes,
           for_grid->size.no_part,
           for_grid->time.potential,
           for_grid->time.density,
           for_grid->time.DK,
           for_grid->time.grid,
           for_grid->time.hydro,
           for_grid->no_sweeps,
           for_grid->cur_resid);
   
   grid_timing.potential += for_grid->time.potential;
   grid_timing.density   += for_grid->time.density;
   grid_timing.DK        += for_grid->time.DK;
   grid_timing.grid      += for_grid->time.grid;
   grid_timing.hydro     += for_grid->time.hydro;
   RAM                   += for_grid->size.no_nodes * global.bytes_node * bytes2GB;
   
   for_grid->time.potential = 0;
   for_grid->time.density   = 0;
   for_grid->time.DK        = 0;
   for_grid->time.grid      = 0;
   for_grid->time.hydro     = 0;

   
   
   
   /* add temporal FFT grid to RAM */
#ifndef AHFstep
#ifdef FFT
   RAM += 2*pow3(global.dom_grid->l1dim) * 4 * bytes2GB;
#endif
#endif
   /* add particles to RAM */
#ifndef WITH_MPI
   RAM += global.no_part * global.bytes_part * bytes2GB;
#else
   RAM += global_info.no_part * global.bytes_part * bytes2GB;
#endif /* WITH_MPI */
   
   /* finish logfile entries for this step... */
   total_time_step = (grid_timing.potential +
                      grid_timing.density   +
                      grid_timing.DK        +
                      grid_timing.grid      +
                      grid_timing.hydro     +
                      ahf.time              +
                      energy.time           +
                      PkSpectrum.time        )     /3600.;
   global.total_time += total_time_step;
   
   fprintf(io.logfile, "                                     %6ld %6ld %6ld %6ld %6ld\n", grid_timing.potential, grid_timing.density, grid_timing.DK, grid_timing.grid, grid_timing.hydro);
   fprintf(io.logfile, "force resolution    ~ %8.2g kpc/h\n", 3000.*simu.boxsize/(double)for_grid->l1dim);
   fprintf(io.logfile, "total time for step = %8.0f seconds (%8.3g hours)\n", total_time_step * 3600., total_time_step);
#ifdef AHF
   if(ahf.time > 0)
      fprintf(io.logfile, "total time for AHF  = %8ld seconds (%8.3g hours)\n",ahf.time,ahf.time/3600.);
#endif
   fprintf(io.logfile, "total time          = %8.3g hours (%g days)\n", global.total_time, global.total_time/24.);
   fprintf(io.logfile, "memory during step  = %8.3g GB\n",RAM);
#	ifdef MPI_TIMING
	global_mpi.stop = MPI_Wtime();
	io_logging_msg(global_io.log, INT32_C(0),
	               "total time (MPI measure) =  %f seconds (%f hours)",
	               (global_mpi.stop-global_mpi.start),
	               (global_mpi.stop-global_mpi.start)/3600.);
#	endif
   
   fflush(io.logfile);
   
   ahf.time        = 0;
   energy.time     = 0;
   PkSpectrum.time = 0;
}


/*==============================================================================
 * keep user up-to-date with simulation parameters initialized by startrun()
 *==============================================================================*/
void log_structures()
{
  int i;
  
  /*===========================================================
   * write information to logfile
   *===========================================================*/
  if(global.restart == TRUE)
    {
      fprintf(io.logfile,"restarting simulation with AMIGA v%3.1f/%d\n",VERSION,BUILT);
      fprintf(io.logfile,"========================================\n");
    }
  
  fprintf(io.logfile,"\n");
  fprintf(io.logfile,"simu.\n");
  fprintf(io.logfile,"=====\n");
  fprintf(io.logfile,"no_part       = %ld (about %d in 1D)\n", simu.no_part, (int)(pow((double)simu.no_part,(double)0.33333333333)+0.5));
  fprintf(io.logfile,"no_vpart      = %g (about %d in 1D)\n",  simu.no_vpart,(int)(pow((double)simu.no_vpart,(double)0.33333333333)+0.5));
  fprintf(io.logfile,"no_species    = %ld\n", simu.no_species);
  fprintf(io.logfile,"min_weight    = %g\n",  simu.min_weight);
  fprintf(io.logfile,"max_weight    = %g\n",  simu.max_weight);
  fprintf(io.logfile,"NGRID_DOM     = %d\n",  simu.NGRID_DOM);
  fprintf(io.logfile,"NGRID_MIN     = %d\n",  simu.NGRID_MIN);
  fprintf(io.logfile,"NGRID_MAX     = %d\n",  simu.NGRID_MAX);
  fprintf(io.logfile,"Nth_dom       = %f\n",  simu.Nth_dom);
  fprintf(io.logfile,"Nth_ref       = %f\n",  simu.Nth_ref);
  fprintf(io.logfile,"SHIFT         = %f\n",  simu.SHIFT);
  fprintf(io.logfile,"np_limit      = %d\n",  simu.np_limit);
  fprintf(io.logfile,"mmfocus       = %d\n",  simu.mmfocus);
  fprintf(io.logfile,"multi_mass    = %d\n",  simu.multi_mass);
  fprintf(io.logfile,"double_precsn = %d\n",  simu.double_precision);
  fprintf(io.logfile,"hydro         = %d\n",  simu.hydro);
  fprintf(io.logfile,"magneto       = %d\n",  simu.magneto);
  fprintf(io.logfile,"boxsize       = %f\n",  simu.boxsize);
  fprintf(io.logfile,"omega0        = %f\n",  simu.omega0);
  fprintf(io.logfile,"lambda0       = %f\n",  simu.lambda0);
  fprintf(io.logfile,"mean_dens     = %f\n",  simu.mean_dens);
  fprintf(io.logfile,"FourPiG       = %f\n",  simu.FourPiG);
  fprintf(io.logfile,"t_unit        = %f\n",  simu.t_unit);
  fprintf(io.logfile,"pmass         = %f\n",  simu.pmass);
  fprintf(io.logfile,"gamma         = %f\n",  simu.gamma);
  fprintf(io.logfile,"H_frac        = %f\n",  simu.H_frac);
  fprintf(io.logfile,"e_init        = %f\n",  simu.e_init);
  fprintf(io.logfile,"T_init        = %f\n",  simu.T_init);
  fprintf(io.logfile,"omegaDM       = %f\n",  simu.omegaDM);
  fprintf(io.logfile,"omegab        = %f\n",  simu.omegab);
  fprintf(io.logfile,"f_b           = %f\n",  simu.f_b);
#ifdef MOND
  fprintf(io.logfile,"g0            = %f\n",  simu.g0);
  fprintf(io.logfile,"H0            = %f\n",  simu.h0);
#endif
  
  fprintf(io.logfile,"z_initial     = %G\n",  simu.z_initial);
  fprintf(io.logfile,"z_final       = %G\n",  simu.z_final);
  fprintf(io.logfile,"a_initial     = %G\n",  simu.a_initial);
  fprintf(io.logfile,"a_final       = %G\n",  simu.a_final);
  fprintf(io.logfile,"super_t_init  = %G\n",  simu.super_t_initial);
  fprintf(io.logfile,"super_t_final = %G\n",  simu.super_t_final);
#ifdef LIGHTCONE
  fprintf(io.logfile,"z_lightcone   = %f\n",  simu.z_lightcone);
#endif
  
  
  fprintf(io.logfile," -> integration performed in supercomoving coordinates!\n");
#ifdef ATS
  fprintf(io.logfile," -> integration performed using adaptive time stepping!\n");
#endif
  
  fprintf(io.logfile,"\n");
  fprintf(io.logfile,"global.\n");
  fprintf(io.logfile,"=======\n");
  fprintf(io.logfile,"architecture  = %d  (0:little_endian, 1:big_endian)\n",  global.architecture);
  fprintf(io.logfile,"no_part       = %ld\n", global.no_part);
  fprintf(io.logfile,"no_timestep   = %d\n",  global.no_timestep);
  fprintf(io.logfile,"z             = %g\n",  global.z);
  fprintf(io.logfile,"a             = %g\n",  global.a);
  fprintf(io.logfile,"t             = %g\n",  global.t);
  fprintf(io.logfile,"super_t       = %g\n",  global.super_t);
#ifdef GADGET
  fprintf(io.logfile,"no_gadget_file= %d\n",  gadget.no_gadget_files);
#endif
  
  
  fprintf(io.logfile,"\n");
  fprintf(io.logfile,"io.\n");
  fprintf(io.logfile,"===\n");
  fprintf(io.logfile,"no_outputs    = %d\n",  io.no_outputs);
  fprintf(io.logfile,"out_dumps     = %d\n",  io.out_dumps);
  
  for(i = 0; i < io.no_outputs; i++)
    {
#ifdef ISOLATED
      fprintf(io.logfile,"a_out[%5d] = %8.4g      z=%10.8g       t=%8.4g Gyr       super_t=%8.4g\n",i,io.a_out[i],1/io.a_out[i]-1,calc_t(io.a_out[i]) * 2*simu.t_unit*Mpc/Gyr,calc_super_t(io.a_out[i]));
#else
      fprintf(io.logfile,"a_out[%5d] = %8.4g      z=%10.8g       t=%8.4g       super_t=%8.4g\n",i,io.a_out[i],1/io.a_out[i]-1,calc_t(io.a_out[i]) * simu.t_unit*Mpc/Gyr,calc_super_t(io.a_out[i]));
#endif
    }
  
  /* dump (revised) header information to logfile */
  fprintf(io.logfile,"\n");
  fprintf(io.logfile,"io.header.\n");
  fprintf(io.logfile,"==========\n");
  fprintf(io.logfile,"%s\n",io.header.header);
  fprintf(io.logfile,"header.multi_mass            = %d\n",io.header.multi_mass);
  fprintf(io.logfile,"header.double_precision      = %d\n",io.header.double_precision);
  fprintf(io.logfile,"header.hydro                 = %d\n",io.header.multi_mass);
  fprintf(io.logfile,"header.magneto               = %d\n",io.header.magneto);
  fprintf(io.logfile,"header.no_part               = %ld\n",io.header.no_part);
  fprintf(io.logfile,"header.no_species            = %ld\n",io.header.no_species);
  fprintf(io.logfile,"header.min_weight            = %g\n",io.header.min_weight);
  fprintf(io.logfile,"header.max_weight            = %g\n",io.header.max_weight);
  fprintf(io.logfile,"header.no_vpart              = %g\n",io.header.no_vpart);
  fprintf(io.logfile,"header.timestep              = %g\n",io.header.timestep);
  fprintf(io.logfile,"header.no_timestep           = %d\n",io.header.no_timestep);
  fprintf(io.logfile,"header.boxsize               = %g\n",io.header.boxsize);
  fprintf(io.logfile,"header.omega0                = %g\n",io.header.omega0);
  fprintf(io.logfile,"header.omegab                = %g\n",io.header.omegab);
  fprintf(io.logfile,"header.lambda0               = %g\n",io.header.lambda0);   
  fprintf(io.logfile,"header.gamma                 = %g\n",io.header.gamma);
  fprintf(io.logfile,"header.H_frac                = %g\n",io.header.H_frac);
  fprintf(io.logfile,"header.T_init                = %g\n",io.header.T_init);
  fprintf(io.logfile,"header.pmass                 = %g\n",io.header.pmass);
  fprintf(io.logfile,"header.t_unit                = %g\n",io.header.t_unit);
  fprintf(io.logfile,"header.cur_reflevel          = %g\n",io.header.cur_reflevel);
  fprintf(io.logfile,"header.cur_frcres            = %g\n",io.header.cur_frcres);
  fprintf(io.logfile,"header.a_initial             = %g\n",io.header.a_initial);
  fprintf(io.logfile,"header.a_current             = %g\n",io.header.a_current);
  fprintf(io.logfile,"header.K_initial             = %g\n",io.header.K_initial);
  fprintf(io.logfile,"header.K_current             = %g\n",io.header.K_current);
  fprintf(io.logfile,"header.U_initial             = %g\n",io.header.U_initial);
  fprintf(io.logfile,"header.U_current             = %g\n",io.header.U_current);
  fprintf(io.logfile,"header.Eintegral             = %g\n",io.header.Eintegral);
  fprintf(io.logfile,"header.Econst                = %g\n",io.header.Econst);
  fprintf(io.logfile,"header.paramNSTEPS           = %g\n",io.header.paramNSTEPS);
  fprintf(io.logfile,"header.paramNGRID_DOM        = %g\n",io.header.paramNGRID_DOM);   
  fprintf(io.logfile,"header.paramNth_dom          = %g\n",io.header.paramNth_dom);
  fprintf(io.logfile,"header.paramNth_ref          = %g\n",io.header.paramNth_ref);
  fprintf(io.logfile,"header.paramE_UPDATE         = %g\n",io.header.paramE_UPDATE);
  fprintf(io.logfile,"header.paramCELLFRAC_MAX     = %g\n",io.header.paramCELLFRAC_MAX);
  fprintf(io.logfile,"header.paramCELLFRAC_MIN     = %g\n",io.header.paramCELLFRAC_MIN);
  fprintf(io.logfile,"header.paramCA_CRIT          = %g\n",io.header.paramCA_CRIT);
  fprintf(io.logfile,"header.paramMAX_L1DIM        = %g\n",io.header.paramMAX_L1DIM);
  fprintf(io.logfile,"header.paramDOMSWEEPS        = %g\n",io.header.paramDOMSWEEPS);
  fprintf(io.logfile,"header.paramREFSWEEPS        = %g\n",io.header.paramREFSWEEPS);
  fprintf(io.logfile,"header.paramAHF_MINPART      = %g\n",io.header.paramAHF_MINPART);
  fprintf(io.logfile,"header.paramAHF_VTUNE        = %g\n",io.header.paramAHF_VTUNE);
  fprintf(io.logfile,"header.paramAHF_RISE         = %g\n",io.header.paramAHF_RISE);
  fprintf(io.logfile,"header.paramAHF_SLOPE        = %g\n",io.header.paramAHF_SLOPE);
  fprintf(io.logfile,"header.paramAHF_MAXNRISE     = %g\n\n",io.header.paramAHF_MAXNRISE);
  fflush(io.logfile);   
  
#ifdef POWERSPECTRUM
  fprintf(io.logfile,"\n");
  fprintf(io.logfile,"PkSpectrum.\n");
  fprintf(io.logfile,"===========\n");
  fprintf(io.logfile,"Pk_ini        = %g\n",  PkSpectrum.Pk_ini);
  fprintf(io.logfile,"Pk_now        = %g\n",  PkSpectrum.Pk_now);
  fprintf(io.logfile,"Dgrowth_ini   = %g\n",  PkSpectrum.Dgrowth_ini);
  fprintf(io.logfile,"Dgrowth_now   = %g\n",  PkSpectrum.Dgrowth_now);
  fprintf(io.logfile,"dump_Pk       = %d\n",  PkSpectrum.dump_Pk);
#endif
  
#ifdef LIGHTCONE 
  fprintf(io.logfile,"lightcone:\n");
  fprintf(io.logfile,"==========\n");
  fprintf(io.logfile,"   lightcone starts at z = %G\n",  simu.z_lightcone);
  fprintf(io.logfile,"   lightcone type = %d\n",  io.conetype);
  
  if(io.conetype==-1) {
    fprintf(io.logfile,"   patch angles: a=%g, b=%g, g=%g\n",io.sa/D2R,io.sb/D2R,io.sg/D2R);
    fprintf(io.logfile,"   patch sizes: dx=%g, dy=%g\n",io.dcx*2./D2R,io.dcy*2./D2R);
    fprintf(io.logfile,"\n   Euler rotation matrix:\n");
    fprintf(io.logfile,"   %g\t%g\t%g\n",io.R[0][0],io.R[0][1],io.R[0][2]);
    fprintf(io.logfile,"   %g\t%g\t%g\n",io.R[1][0],io.R[1][1],io.R[1][2]);
    fprintf(io.logfile,"   %g\t%g\t%g\n",io.R[2][0],io.R[2][1],io.R[2][2]);
  }
#endif /* LIGHTCONE */
  
  
  fprintf(io.logfile,"\n\n==================================================================================\n\n");
  fflush(io.logfile);
  
}


/*==============================================================================
 * just dump all relevant parameter into a file >prefix<_parameter
 *==============================================================================*/
void write_parameterfile()
{
  FILE  *fpparam;
  char   param_file[MAXSTRING];
  
#ifdef WITH_MPI
  /* make sure to only write one parameter file when running on multile CPUs */
  if(global_mpi.rank == 0)
#endif
    if(io.outfile_prefix != NULL)
      {
#	ifdef NEWSTARTRUN
        strcpy(param_file, global_io.params->outfile_prefix);
        strcat(param_file,"parameter");
#	else
        strcpy(param_file,io.outfile_prefix);
        strcat(param_file,"parameter");
#	endif
        
        /* check if there is already a paramater file */
        if( (fpparam = fopen(param_file,"r")) != NULL)
          {
            fclose(fpparam);
            fprintf(io.logfile," NOTE: parameter file %s already exists\n",param_file);
            strcat(param_file,"-restart");
            fprintf(io.logfile,"       writing parameter to %s instead!\n",param_file);
          }
        
        
        if( (fpparam = fopen(param_file,"w")) == NULL)
          {
            fprintf(io.logfile," NOTE: could not open %s to dump parameter\n",param_file);
          }
        else
          {
            /* write >>all<< relevant paramater to file */
            fprintf(fpparam,"*=========================================================================*\n");
            fprintf(fpparam,"      parameter file for AMIGA version %g (built %d)\n",VERSION,BUILT);
            fprintf(fpparam,"*=========================================================================*\n");
            
            fprintf(fpparam,"\n");
            fprintf(fpparam,"general parameter:\n");
            fprintf(fpparam,"------------------\n");
            fprintf(fpparam,"TERMINATE_AMIGA            \t\t%s\n",      global.termfile_name);
            fprintf(fpparam,"E_UPDATE                   \t\t%d\n",                  E_UPDATE);
            fprintf(fpparam,"MAX_L1DIM                  \t\t%d\n",                 MAX_L1DIM);
            fprintf(fpparam,"MIN_NNODES                 \t\t%d\n",                MIN_NNODES);
            fprintf(fpparam,"PKMODE                     \t\t%d\n",                    PKMODE);
            fprintf(fpparam,"PKGROWLIMIT                \t\t%g\n",               PKGROWLIMIT);
            fprintf(fpparam,"NSTEPS                     \t\t%d\n",                    NSTEPS);
            fprintf(fpparam,"MAXTIME                    \t\t%d\n",                   MAXTIME);
            fprintf(fpparam,"MAXSTRING                  \t\t%d\n",                 MAXSTRING);
            fprintf(fpparam,"AMIGAHEADER                \t\t%d\n",               AMIGAHEADER);
            fprintf(fpparam,"HEADERSTRING               \t\t%d\n",              HEADERSTRING);
            fprintf(fpparam,"HEADERSIZE                 \t\t%d\n",           (int)HEADERSIZE);
            fprintf(fpparam,"FILLHEADER                 \t\t%d\n",           (int)FILLHEADER);
            fprintf(fpparam,"CRITMULTI                  \t\t%g\n",                 CRITMULTI);
            fprintf(fpparam,"NP_RATIO                   \t\t%g\n",                  NP_RATIO);
            fprintf(fpparam,"ZERO                       \t\t%g\n",                      ZERO);
            fprintf(fpparam,"MACHINE_ZERO               \t\t%g\n",              MACHINE_ZERO);
            fprintf(fpparam,"MZERO                      \t\t%g\n",                     MZERO);
            
            fprintf(fpparam,"\n");
            fprintf(fpparam,"gravity solver related parameter:\n");
            fprintf(fpparam,"---------------------------------\n");
            fprintf(fpparam,"DOMSWEEPS                  \t\t%d\n",                 DOMSWEEPS);
            fprintf(fpparam,"REFSWEEPS                  \t\t%d\n",                 REFSWEEPS);
            fprintf(fpparam,"W_SOR                      \t\t%g\n",                     W_SOR);
            fprintf(fpparam,"ETA                        \t\t%g\n",                       ETA);
            fprintf(fpparam,"CONVCRIT                   \t\t%g\n",                  CONVCRIT);
            fprintf(fpparam,"DOMCORRECT                 \t\t%g\n",                DOMCORRECT);
            
            fprintf(fpparam,"\n");
            fprintf(fpparam,"hydro solver related parameter:\n");
            fprintf(fpparam,"-------------------------------\n");
            fprintf(fpparam,"CFL_TUNE                   \t\t%g\n",                  CFL_TUNE);
            fprintf(fpparam,"ETA_DUAL_ENERGY1           \t\t%g\n",           ETA_DUAL_ENERGY1);
            fprintf(fpparam,"ETA_DUAL_ENERGY2           \t\t%g\n",           ETA_DUAL_ENERGY2);
            
            fprintf(fpparam,"\n");
            fprintf(fpparam,"DKD related parameter:\n");
            fprintf(fpparam,"----------------------\n");
            fprintf(fpparam,"CELLFRAC_MAX               \t\t%g\n",              CELLFRAC_MAX);
            fprintf(fpparam,"CELLFRAC_MIN               \t\t%g\n",              CELLFRAC_MIN);
            fprintf(fpparam,"CF_MEAN                    \t\t%g\n",                   CF_MEAN);
            fprintf(fpparam,"CA_CRIT                    \t\t%g\n",                   CA_CRIT);
            fprintf(fpparam,"SPEEDFRAC                  \t\t%g\n",                 SPEEDFRAC);
            
            fprintf(fpparam,"\n");
            fprintf(fpparam,"MPI related parameter:\n");
            fprintf(fpparam,"----------------------\n");
            fprintf(fpparam,"BITS_PER_DIMENSION         \t\t%d\n",        BITS_PER_DIMENSION);
            fprintf(fpparam,"LOADBALANCE_DOMAIN_LEVEL   \t\t%d\n",  LOADBALANCE_DOMAIN_LEVEL);
            fprintf(fpparam,"MAX_SEND_PARTICLES         \t\t%d\n",        MAX_SEND_PARTICLES);
            fprintf(fpparam,"VERBOSITY                  \t\t%d\n",                 VERBOSITY);
            
            fprintf(fpparam,"\n");
            fprintf(fpparam,"AHF related parameter:\n");
            fprintf(fpparam,"----------------------\n");
            fprintf(fpparam,"AHF_MINPART                \t\t%d\n",               AHF_MINPART);
            fprintf(fpparam,"AHF_VTUNE                  \t\t%g\n",                 AHF_VTUNE);
            fprintf(fpparam,"AHF_RISE                   \t\t%g\n",                  AHF_RISE);
            fprintf(fpparam,"AHF_SLOPE                  \t\t%g\n",                 AHF_SLOPE);
            fprintf(fpparam,"AHF_MAXNRISE               \t\t%d\n",              AHF_MAXNRISE);
            fprintf(fpparam,"AHF_MASSMIX                \t\t%g\n",               AHF_MASSMIX);
            fprintf(fpparam,"AHF_MAXHALO                \t\t%d\n",               AHF_MAXHALO);
            fprintf(fpparam,"AHF_MAX_GATHER_RAD         \t\t%g\n",        AHF_MAX_GATHER_RAD);
            
            fprintf(fpparam,"\n");
            fprintf(fpparam,"GADGET related parameter:\n");
            fprintf(fpparam,"-------------------------\n");
            fprintf(fpparam,"GADGET_MUNIT               \t\t%g\n",              GADGET_MUNIT);
            fprintf(fpparam,"GADGET_LUNIT               \t\t%g\n",              GADGET_LUNIT);
            
            fprintf(fpparam,"\n");
            fprintf(fpparam,"physical constants:\n");
            fprintf(fpparam,"-------------------\n");
            fprintf(fpparam,"Gyr                        \t\t%g\n",                       Gyr);
            fprintf(fpparam,"Mpc                        \t\t%g\n",                       Mpc);
            fprintf(fpparam,"H0                         \t\t%g\n",                        H0);
            fprintf(fpparam,"rhoc0                      \t\t%g\n",                     rhoc0);
            fprintf(fpparam,"Grav                       \t\t%g\n",                      Grav);
            fprintf(fpparam,"cH0                        \t\t%g\n",                       cH0);
            fprintf(fpparam,"kB_per_mp                  \t\t%g\n",                 kB_per_mp);
            fprintf(fpparam,"kBoltzman                  \t\t%g\n",                 kBoltzman);
            
            fprintf(fpparam,"\n");
            fprintf(fpparam,"*=====================*\n");
            fprintf(fpparam,"     DEFINEFLAGS:\n");
            fprintf(fpparam,"*=====================*\n");
            fprintf(fpparam,"\n");
            fprintf(fpparam,"various flags:\n");
            fprintf(fpparam,"--------------\n");
#ifdef NONC99
            fprintf(fpparam,"NONC99=                    \t\t%d\n",NONC99);
#else
            fprintf(fpparam,"NONC99=                    \t\t undefined\n");
#endif
#ifdef WITH_MPI
            fprintf(fpparam,"WITH_MPI                   \t\t1\n");
#else
            fprintf(fpparam,"WITH_MPI                   \t\t0\n");
#endif
#ifdef NEWSTARTRUN
            fprintf(fpparam,"NEWSTARTRUN                \t\t1\n");
#else
            fprintf(fpparam,"NEWSTARTRUN                \t\t0\n");
#endif
#ifdef BYTESWAP
            fprintf(fpparam,"BYTESWAP                   \t\t1\n");
#else
            fprintf(fpparam,"BYTESWAP                   \t\t0\n");
#endif
#ifdef MULTIMASS
            fprintf(fpparam,"MULTIMASS                  \t\t1\n");
#else
            fprintf(fpparam,"MULTIMASS                  \t\t0\n");
#endif
#ifdef ATS
            fprintf(fpparam,"ATS                        \t\t1\n");
#else
            fprintf(fpparam,"ATS                        \t\t0\n");
#endif
#ifdef VERBOSE
            fprintf(fpparam,"VERBOSE                    \t\t1\n");
#else
            fprintf(fpparam,"VERBOSE                    \t\t0\n");
#endif
#ifdef CONVERT
            fprintf(fpparam,"CONVERT                    \t\t1\n");
#else
            fprintf(fpparam,"CONVERT                    \t\t0\n");
#endif
#ifdef CONVERT_TERM
            fprintf(fpparam,"CONVERT_TERM               \t\t1\n");
#else
            fprintf(fpparam,"CONVERT_TERM               \t\t0\n");
#endif
#ifdef DOUBLE
            fprintf(fpparam,"DOUBLE                     \t\t1\n");
#else
            fprintf(fpparam,"DOUBLE                     \t\t0\n");
#endif
#ifdef LIGHTCONE
            fprintf(fpparam,"LIGHTCONE                  \t\t1\n");
#else
            fprintf(fpparam,"LIGHTCONE                  \t\t0\n");
#endif
#ifdef OUTDUMPS
            fprintf(fpparam,"OUTDUMPS                   \t\t1\n");
#else
            fprintf(fpparam,"OUTDUMPS                   \t\t0\n");
#endif
#ifdef POWERSPECTRUM
            fprintf(fpparam,"POWERSPECTRUM              \t\t1\n");
#else
            fprintf(fpparam,"POWERSPECTRUM              \t\t0\n");
#endif
            
            fprintf(fpparam,"\n");
            fprintf(fpparam,"grid related flags:\n");
            fprintf(fpparam,"-------------------\n");
#ifdef ISOLATED
            fprintf(fpparam,"ISOLATED                   \t\t1\n");
#else
            fprintf(fpparam,"ISOLATED                   \t\t0\n");
#endif
#ifdef PM
            fprintf(fpparam,"PM                         \t\t1\n");
#else
            fprintf(fpparam,"PM                         \t\t0\n");
#endif
#ifdef TSC
            fprintf(fpparam,"TSC                        \t\t1\n");
#else
            fprintf(fpparam,"TSC                        \t\t0\n");
#endif
#ifdef CIC
            fprintf(fpparam,"CIC                        \t\t1\n");
#else
            fprintf(fpparam,"CIC                        \t\t0\n");
#endif
#ifdef NGP
            fprintf(fpparam,"NGP                        \t\t1\n");
#else
            fprintf(fpparam,"NGP                        \t\t0\n");
#endif
            
            
            
            fprintf(fpparam,"\n");
            fprintf(fpparam,"gravity solver related flags:\n");
            fprintf(fpparam,"-----------------------------\n");
#ifdef FFT
            fprintf(fpparam,"FFT                        \t\t1\n");
#else
            fprintf(fpparam,"FFT                        \t\t0\n");
#endif
#ifdef FIXED_SWEEPS
            fprintf(fpparam,"FIXED_SWEEPS               \t\t1\n");
#else
            fprintf(fpparam,"FIXED_SWEEPS               \t\t0\n");
#endif
#ifdef NO_SOR
            fprintf(fpparam,"NO_SOR                     \t\t1\n");
#else
            fprintf(fpparam,"NO_SOR                     \t\t0\n");
#endif
#ifdef SWAP_XY
            fprintf(fpparam,"SWAP_XY                    \t\t1\n");
#else
            fprintf(fpparam,"SWAP_XY                    \t\t0\n");
#endif
#ifdef SWAP_XZ
            fprintf(fpparam,"SWAP_XZ                    \t\t1\n");
#else
            fprintf(fpparam,"SWAP_XZ                    \t\t0\n");
#endif
            
            
            fprintf(fpparam,"\n");
            fprintf(fpparam,"MHD solver related flags:\n");
            fprintf(fpparam,"-------------------------\n");
#ifdef HYDRO
            fprintf(fpparam,"HYDRO                      \t\t1\n");
#else
            fprintf(fpparam,"HYDRO                      \t\t0\n");
#endif
#ifdef MHD
            fprintf(fpparam,"MHD                        \t\t1\n");
#else
            fprintf(fpparam,"MHD                        \t\t0\n");
#endif
#ifdef MONOATOMIC
            fprintf(fpparam,"MONOATOMIC                 \t\t1\n");
#else
            fprintf(fpparam,"MONOATOMIC                 \t\t0\n");
#endif
#ifdef DUAL_ENERGY
            fprintf(fpparam,"DUAL_ENERGY                \t\t1\n");
#else
            fprintf(fpparam,"DUAL_ENERGY                \t\t0\n");
#endif
#ifdef RYU_CRITERION
            fprintf(fpparam,"RYU_CRITERION              \t\t1\n");
#else
            fprintf(fpparam,"RYU_CRITERION              \t\t0\n");
#endif
#ifdef SLOPE_vanLeer
            fprintf(fpparam,"SLOPE_vanLeer              \t\t1\n");
#else
            fprintf(fpparam,"SLOPE_vanLeer              \t\t0\n");
#endif
#ifdef SLOPE_minmod
            fprintf(fpparam,"SLOPE_minmod               \t\t1\n");
#else
            fprintf(fpparam,"SLOPE_minmod               \t\t0\n");
#endif
#ifdef SLOPE_superbee
            fprintf(fpparam,"SLOPE_superbee             \t\t1\n");
#else
            fprintf(fpparam,"SLOPE_superbee             \t\t0\n");
#endif
#ifdef RECONSTRUCT_PRESSURE
            fprintf(fpparam,"RECONSTRUCT_PRESSURE       \t\t1\n");
#else
            fprintf(fpparam,"RECONSTRUCT_PRESSURE       \t\t0\n");
#endif
#ifdef RECONSTRUCT_VELOCITIES
            fprintf(fpparam,"RECONSTRUCT_VELOCITIES     \t\t1\n");
#else
            fprintf(fpparam,"RECONSTRUCT_VELOCITIES     \t\t0\n");
#endif
#ifdef FLAT_RECONSTRUCTION
            fprintf(fpparam,"FLAT_RECONSTRUCTION        \t\t1\n");
#else
            fprintf(fpparam,"FLAT_RECONSTRUCTION        \t\t0\n");
#endif
#ifdef CHECK_RECONSTRUCTION
            fprintf(fpparam,"CHECK_RECONSTRUCTION       \t\t1\n");
#else
            fprintf(fpparam,"CHECK_RECONSTRUCTION       \t\t0\n");
#endif
#ifdef CHECK_ADVECTION
            fprintf(fpparam,"CHECK_ADVECTION            \t\t1\n");
#else
            fprintf(fpparam,"CHECK_ADVECTION            \t\t0\n");
#endif
            
            
            fprintf(fpparam,"\n");
            fprintf(fpparam,"AHF related flags:\n");
            fprintf(fpparam,"------------------\n");
#ifdef AHF
            fprintf(fpparam,"AHF                        \t\t1\n");
#else
            fprintf(fpparam,"AHF                        \t\t0\n");
#endif
#ifdef AHFstep
            fprintf(fpparam,"AHFstep                    \t\t1\n");
#else
            fprintf(fpparam,"AHFstep                    \t\t0\n");
#endif
#ifdef AHFisodensity
            fprintf(fpparam,"AHFisodensity              \t\t1\n");
#else
            fprintf(fpparam,"AHFisodensity              \t\t0\n");
#endif
#ifdef AHFmaxdenscentre
            fprintf(fpparam,"AHFmaxdenscentre           \t\t1\n");
#else
            fprintf(fpparam,"AHFmaxdenscentre           \t\t0\n");
#endif
#ifdef AHFmaxhalo
            fprintf(fpparam,"AHFmaxhalo                 \t\t1\n");
#else
            fprintf(fpparam,"AHFmaxhalo                 \t\t0\n");
#endif
#ifdef AHFmmfocus
            fprintf(fpparam,"AHFmmfocus                 \t\t1\n");
#else
            fprintf(fpparam,"AHFmmfocus                 \t\t0\n");
#endif
#ifdef AHFnoremunbound
            fprintf(fpparam,"AHFnoremunbound            \t\t1\n");
#else
            fprintf(fpparam,"AHFnoremunbound            \t\t0\n");
#endif
#ifdef AHFgeom
            fprintf(fpparam,"AHFgeom                    \t\t1\n");
#else
            fprintf(fpparam,"AHFgeom                    \t\t0\n");
#endif
#ifdef AHFmaxdenscentre
            fprintf(fpparam,"AHFmaxdenscentre           \t\t1\n");
#else
            fprintf(fpparam,"AHFmaxdenscentre           \t\t0\n");
#endif
#ifdef AHFgeomcentre
            fprintf(fpparam,"AHFgeomcentre              \t\t1\n");
#else
            fprintf(fpparam,"AHFgeomcentre              \t\t0\n");
#endif
#ifdef AHFcomcentre
            fprintf(fpparam,"AHFcomcentre               \t\t1\n");
#else
            fprintf(fpparam,"AHFcomcentre               \t\t0\n");
#endif
#ifdef AHFpotcentre
            fprintf(fpparam,"AHFpotcentre               \t\t1\n");
#else
            fprintf(fpparam,"AHFpotcentre               \t\t0\n");
#endif
#ifdef AHFphi_infty
            fprintf(fpparam,"AHFphi_infty               \t\t1\n");
#else
            fprintf(fpparam,"AHFphi_infty               \t\t0\n");
#endif
#ifdef AHFphi_GMr
            fprintf(fpparam,"AHFphi_GMr                 \t\t1\n");
#else
            fprintf(fpparam,"AHFphi_GMr                 \t\t0\n");
#endif
#ifdef AHFsplinefit
            fprintf(fpparam,"AHFsplinefit               \t\t1\n");
#else
            fprintf(fpparam,"AHFsplinefit               \t\t0\n");
#endif
#ifdef AHFprofilerise
            fprintf(fpparam,"AHFprofilerise             \t\t1\n");
#else
            fprintf(fpparam,"AHFprofilerise             \t\t0\n");
#endif
#ifdef MANUAL_DVIR
            fprintf(fpparam,"MANUAL_DVIR                \t\t1\n");
#else
            fprintf(fpparam,"MANUAL_DVIR                \t\t0\n");
#endif
#ifdef AHFreducedinertiatensor
            fprintf(fpparam,"AHFreducedinertiatensor     \t\t1\n");
#else
            fprintf(fpparam,"AHFreducedinertiatensor     \t\t0\n");
#endif
#ifdef AHFcentrefile
            fprintf(fpparam,"AHFcentrefile              \t\t1\n");
#else
            fprintf(fpparam,"AHFcentrefile              \t\t0\n");
#endif
#ifdef AHFsubstructure
            fprintf(fpparam,"AHFsubstructure            \t\t1\n");
#else
            fprintf(fpparam,"AHFsubstructure            \t\t0\n");
#endif
#ifdef AHFgridsubstructure
            fprintf(fpparam,"AHFgridsubstructure        \t\t1\n");
#else
            fprintf(fpparam,"AHFgridsubstructure        \t\t0\n");
#endif
#ifdef AHFabsangmom
            fprintf(fpparam,"AHFabsangmom            \t\t1\n");
#else
            fprintf(fpparam,"AHFabsangmom            \t\t0\n");
#endif
#ifdef PARDAU_DISTANCE
            fprintf(fpparam,"PARDAU_DISTANCE            \t\t1\n");
#else
            fprintf(fpparam,"PARDAU_DISTANCE            \t\t0\n");
#endif
#ifdef PARDAU_NODES
            fprintf(fpparam,"PARDAU_NODES               \t\t1\n");
#else
            fprintf(fpparam,"PARDAU_NODES               \t\t0\n");
#endif
#ifdef PARDAU_PARTS
            fprintf(fpparam,"PARDAU_PARTS               \t\t1\n");
#else
            fprintf(fpparam,"PARDAU_PARTS               \t\t0\n");
#endif
#ifdef AHFgridinfofile
            fprintf(fpparam,"AHFgridinfofile          \t\t1\n");
#else
            fprintf(fpparam,"AHFgridinfofile          \t\t0\n");
#endif
            
            
            
            fprintf(fpparam,"\n");
            fprintf(fpparam,"GADGET related flags:\n");
            fprintf(fpparam,"---------------------\n");
#ifdef GADGET
            fprintf(fpparam,"GADGET                     \t\t1\n");
#else
            fprintf(fpparam,"GADGET                     \t\t0\n");
#endif
#ifdef GADGET2
            fprintf(fpparam,"GADGET2                    \t\t1\n");
#else
            fprintf(fpparam,"GADGET2                    \t\t0\n");
#endif
#ifdef GADGET_GAS_ONLY
            fprintf(fpparam,"GADGET_GAS_ONLY            \t\t1\n");
#else
            fprintf(fpparam,"GADGET_GAS_ONLY            \t\t0\n");
#endif
#ifdef GADGET_STARS_ONLY
            fprintf(fpparam,"GADGET_STARS_ONLY          \t\t1\n");
#else
            fprintf(fpparam,"GADGET_STARS_ONLY          \t\t0\n");
#endif
#ifdef PARTICLES_INFO
            fprintf(fpparam,"PARTICLES_INFO             \t\t1\n");
#else
            fprintf(fpparam,"PARTICLES_INFO             \t\t0\n");
#endif
#ifdef GADGET_LUNIT_KPC
            fprintf(fpparam,"GADGET_LUNIT_KPC           \t\t1\n");
#else
            fprintf(fpparam,"GADGET_LUNIT_KPC           \t\t0\n");
#endif
            
            
            fprintf(fpparam,"\n");
            fprintf(fpparam,"debugging flags:\n");
            fprintf(fpparam,"----------------\n");
#ifdef REF_TEST
            fprintf(fpparam,"REF_TEST                   \t\t1\n");
#else
            fprintf(fpparam,"REF_TEST                   \t\t0\n");
#endif
            
            fprintf(fpparam,"\n");
            fprintf(fpparam,"MOND flags:\n");
            fprintf(fpparam,"-----------\n");
#ifdef MOND
            fprintf(fpparam,"MOND                       \t\t1\n");
#else
            fprintf(fpparam,"MOND                       \t\t0\n");
#endif
            
            
            
            fprintf(fpparam,"\n");
            fprintf(fpparam,"other input file formats:\n");
            fprintf(fpparam,"-------------------------\n");
#ifdef ART
            fprintf(fpparam,"ART                        \t\t1\n");
#else
            fprintf(fpparam,"ART                        \t\t0\n");
#endif
#ifdef TIPSY
            fprintf(fpparam,"TIPSY                      \t\t1\n");
#else
            fprintf(fpparam,"TIPSY                      \t\t0\n");
#endif
#ifdef ASCII
            fprintf(fpparam,"ASCII                      \t\t1\n");
#else
            fprintf(fpparam,"ASCII                      \t\t0\n");
#endif
#ifdef MLAPM
            fprintf(fpparam,"MLAPM                      \t\t1\n");
#else
            fprintf(fpparam,"MLAPM                      \t\t0\n");
#endif
#ifdef AMIGA_ONE_FORMAT
            fprintf(fpparam,"AMIGA_ONE_FORMAT           \t\t1\n");
#else
            fprintf(fpparam,"AMIGA_ONE_FORMAT           \t\t0\n");
#endif
            
            
            fclose(fpparam);
          }
        
        
      }   
}   


