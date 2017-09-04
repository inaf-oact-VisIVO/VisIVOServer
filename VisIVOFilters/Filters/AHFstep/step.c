#include <stddef.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

/* the important definitions have to be included first */
#include "common.h"
#include "param.h"
#include "tdef.h"

/* ...and now for the actual includes */
#include "libutility/utility.h"
#include "libparticles/particles.h"
#include "libgravity/gravity.h"
#include "libgrids/grids.h"
#include "libmhd/mhd.h"

#ifdef STORE_REFS
#include "store_refs.h"
#endif

#ifdef AHF
#include "libahf/ahf.h"
#endif

/*================================================================================
* take a timestep on a given grid
*================================================================================*/
boolean step(gridls **grid_list, int *no_grids, double timestep, double timecounter)
{
  gridls *cur_grid;                   /* current grid pointer             */
  gridls *fin_grid;                   /* fine grid pointer                */
  gridls *for_grid;                   /* sometimes needed for debugging...*/
  int     fingrid_no, curgrid_no;     /* index of finest and current grid */
  int     grid_no;
  
  boolean fin_ms, cur_ms;             /* multi-step change                */
  
  double  KICK;                       /* KICK operator                    */
  double  DRIFT1, DRIFT2;             /* DRIFT operators                  */
  
  double  x_fac, v_fac;
  
  partptr cur_part;
  
#ifndef AHFstep
  partptr *ll_tmp;
  long unsigned i;
#endif
  
#ifdef ADAPTIVE
  boolean refined;                    /* refinement flag                  */
#endif
#ifdef CHECK_MASS
  float   masscheck=0.0;
#endif
  
#ifdef FORCETEST
  FILE    *CheckFile;
  fvect   TotFrc;
  float   AbsFrc;
  char    CheckFileName[MAXSTRING];
#endif

  /* up to now everything is TRUE ;-) */
  fin_ms = TRUE;
  cur_ms = TRUE;
  
  char file_name[50];  // for mondian test.
  
   /*=========================================================================
    * try to refine finest grid:
    *
    * 1. get (number!) density field on currently finest grid right
    * 2. try to refine grid using 'number of particles per node' criterion
    *=========================================================================*/
#ifdef ADAPTIVE
   /* currently the finest grid... */
  fin_grid = *grid_list + *no_grids - 1;
  
    if(!(((global.fst_cycle == FALSE) && (fin_grid->l1dim == global.fin_l1dim)) || 
         ((fin_grid->l1dim) == simu.NGRID_MAX)))
      {
        /* get number density: "npart/node" */
        fin_grid->time.grid -= time(NULL);
        zero_dens   (fin_grid);
        assign_npart(fin_grid);
        
#ifdef ADAPTIVE2
        if(fin_grid != global.dom_grid)
          {
            zero_dens   (fin_grid-1);
            assign_npart(fin_grid-1);
            adjust_dens(fin_grid, TRUE);
          }
#endif
        fin_grid->time.grid += time(NULL);
        
        /* try to refine grid under "npart/node" criterion */
        refined = gen_refgrid(grid_list, no_grids);
      }
  else
    {
    refined = FALSE;
    } 
  
   /*--------------------------------------------------
    * when refined we make a recursive call to step()
    *--------------------------------------------------*/
  if(refined == TRUE) 
    {
      /* recursive call to STEP */
      fin_ms = step(grid_list, no_grids, timestep/2, timecounter);
      
      if(fin_ms == TRUE) /* finest wasn't destroyed -> switch back to parent grid */
        {
          curgrid_no      = *no_grids  - 2;
          fingrid_no      = *no_grids  - 1;
          global.dom_grid = *grid_list + global.domgrid_no;
          cur_grid        = *grid_list + curgrid_no;
          fin_grid        =  cur_grid  + fingrid_no;
        }
      else               /* finest grid was destroyed... */
        {
          curgrid_no      = *no_grids - 1;
          fingrid_no      = *no_grids - 1;
          global.dom_grid = *grid_list + global.domgrid_no;
          cur_grid        = *grid_list + curgrid_no;
        }
    }
  
   /*---------------------------------------------
   * otherwise stay with the actual finest grid
   *---------------------------------------------*/
  else /* refined == FALSE */
    {
      /* no new refinement => stay with current finest grid */
      curgrid_no      = *no_grids - 1;
      fingrid_no      = *no_grids - 1;
      global.dom_grid = *grid_list + global.domgrid_no;
      cur_grid        = *grid_list + curgrid_no;
      
      /* store force resolution in header info */
      io.header.cur_reflevel = (float)(curgrid_no-global.domgrid_no);
      io.header.cur_frcres   = (float)(simu.boxsize/(double)cur_grid->l1dim);
      io.header.cur_frcres  *= 3000.;
      
      
      energy.time -= time(NULL);
#ifdef ENERGY
      /*==================================================================
       * this is the right time to calculate kinetic and potential energy
       *==================================================================*/
      if(((global.no_timestep % E_UPDATE) == 0) && (global.fst_cycle == TRUE))
        {
          layzer_irvine(timecounter, *grid_list, curgrid_no);
        }
      else if((global.no_timestep == 1) && (global.fst_cycle == TRUE))
        {
          layzer_irvine(timecounter, *grid_list, curgrid_no);
          energy.K_initial = energy.K_current;
          energy.U_initial = energy.U_current;
#ifndef ISOLATED
          energy.econst    = 0.0;
          energy.echeck    = 0.0;
          energy.integral  = 0.0;
#endif
        }
#endif /* ENERGY */
      energy.time += time(NULL);
      
      if(global.fst_cycle == TRUE)
        {
          global.fin_l1dim = cur_grid->l1dim;
          global.fst_cycle = FALSE;
        }
    } /* if(refined) */
  
#else /* ADAPTIVE */
  
  cur_grid   = global.dom_grid;
  curgrid_no = global.domgrid_no;
  fingrid_no = global.domgrid_no;

  energy.time -= time(NULL);
#ifdef ENERGY
  if( (((global.no_timestep % E_UPDATE) == 0) && (global.fst_cycle == TRUE)) || global.restart == TRUE )
    {
      layzer_irvine(timecounter, *grid_list, curgrid_no);
      
      global.restart = FALSE;
    }
  else if((global.no_timestep == 1) && (global.fst_cycle == TRUE))
    {
      layzer_irvine(timecounter, *grid_list, curgrid_no);
      energy.K_initial = energy.K_current;
      energy.U_initial = energy.U_current;
      energy.econst    = 0.0;
      energy.echeck    = 0.0;
      energy.integral  = 0.0;
    }
#endif  /* ENERGY   */
  energy.time += time(NULL);

#endif  /* ADAPTIVE */
  
  
  
   /*=========================================================================
    *
    *         WE NOW HAVE ACCESS TO THE FULL BLOWN GRID HIERACHY
    *
    *=========================================================================*/
  
   /* this is the point where to perform some test... */
#ifdef FORCETEST
  printf("\ngetting density right on all grids...");
  for(for_grid = cur_grid; for_grid >= global.dom_grid; for_grid--)
    restore_dens(for_grid);
  
  if(curgrid_no != fingrid_no) 
    refill_dens(cur_grid);
  
  if(cur_grid != global.dom_grid)
    stack_dens(cur_grid, global.dom_grid+1);
  printf("done\n\n");
  printf("\n=========== HERE STARTS THE INTERESTING PART ================================================\n");
  printf("\nsolving for potential on all grids (global.a=%f)...",global.a);
  
  /* solve for potential from next coarser grid down to current grid */
  solve_gravity(*grid_list, curgrid_no);
  
  /* difference potential to get forces on ALL grids */
  for(for_grid = cur_grid; for_grid >= global.dom_grid; for_grid--)
      calc_forces(for_grid);
  
  
  printf("\n..........forces on all grids available..............\n\n");
  
  
  for(for_grid=cur_grid; for_grid >= global.dom_grid; for_grid--)
    kick_mom(for_grid, (double)0.0);
  strcpy(CheckFileName,io.outfile_prefix);
  strcat(CheckFileName,"forces");
  fprintf(stderr,"writing to file %s\n\n",CheckFileName);
  
  CheckFile = fopen(CheckFileName,"w");
  
  TotFrc[X] = 0.0;
  TotFrc[Y] = 0.0;
  TotFrc[Z] = 0.0;
  AbsFrc    = 0.0;
  
  /* we may want to convert pos[] and mom[] to physical units?! */
  x_fac = simu.boxsize;
  v_fac = H0*simu.boxsize/global.a;
  
  for(cur_part=global.fst_part; cur_part<global.fst_part+global.no_part; cur_part++)
    {
      TotFrc[X] += cur_part->forces[X];
      TotFrc[Y] += cur_part->forces[Y];
      TotFrc[Z] += cur_part->forces[Z];
      
      AbsFrc    += sqrt(pow2(cur_part->forces[X]) + 
                        pow2(cur_part->forces[Y]) + 
                        pow2(cur_part->forces[Z]));
      fprintf(CheckFile,"%f %f %f   %f %f %f   %f %f\n",
              cur_part->pos[X],cur_part->pos[Y],cur_part->pos[Z],
              cur_part->forces[X],cur_part->forces[Y],cur_part->forces[Z],
              cur_part->pot,
              cur_part->dens);
    }
  
  fclose(CheckFile);
  
  
  fprintf(stderr," AbsFrc         = %f (AbsFrc/N = %f)\n",AbsFrc,AbsFrc/(double)simu.no_part);
  fprintf(stderr," TotFrc         = %f\n",sqrt(pow2(TotFrc[X])+pow2(TotFrc[Y])+pow2(TotFrc[Z])));
  fprintf(stderr," TotFrc/AbsFrc  = %f\n",
          sqrt(pow2(TotFrc[X])+pow2(TotFrc[Y])+pow2(TotFrc[Z]))/AbsFrc);
  
  printf("\nSTOP\n");
  exit(0);
#endif /* FORCETEST */
  

  /*=========================================================================
   * Calculating grid statistics and AHF halo centres 
   *=========================================================================*/
  ahf.time -= time(NULL);
#ifdef AHF
  if(global.ioflag == TRUE)
  {
    /* get density right on all grids */
    for(for_grid = cur_grid; for_grid >= global.dom_grid; for_grid--)
       restore_dens(for_grid);      
    if(curgrid_no != fingrid_no) 
       refill_dens(cur_grid);
    if(cur_grid != global.dom_grid)
       stack_dens(cur_grid, global.dom_grid+1);
 
#ifdef AHFpotcentre
     /* solve for potential on domain grids */
#ifdef VERBOSE
  fprintf(stderr,"AHFpotcentres:   solving on domain grid %d ... ",global.dom_grid->l1dim);
    solve_dom_gravity(*grid_list);
    fprintf(stderr,"done\n");
    
     /* solve for potential on refinement grids */
    fprintf(stderr,"                 solving on refinement grids ... ");
    if(curgrid_no > global.domgrid_no)
      for(grid_no = global.domgrid_no+1; grid_no <= curgrid_no; grid_no++)
      {
        fprintf(stderr," %d  ",(global.dom_grid+grid_no)->l1dim);
        solve_ref_gravity(*grid_list, grid_no);
      }
    fprintf(stderr,"done\n");
#endif /* VERBOSE */
    
     /* solve for potential on domain grids */
    solve_dom_gravity(*grid_list);
    
     /* solve for potential on refinement grids */
    if(curgrid_no > global.domgrid_no)
      for(grid_no = global.domgrid_no+1; grid_no <= curgrid_no; grid_no++)
      {
        solve_ref_gravity(*grid_list, grid_no);
      }
#endif /* AHFpotcentres */
    
#ifndef AHFstep
     /* backup linked-list structure */
    if((ll_tmp = (partptr *) calloc(global.no_part, sizeof(partptr))) == NULL)
    {
    fprintf(io.logfile,"step: failed to allocate memory for ll_tmp");
      exit(0);
    }
    cur_part = global.fst_part;
    for(i=0; i<global.no_part; i++)
    {
      ll_tmp[i] = cur_part->ll;
      cur_part++;
    }
#endif
    
     /* get spatially connected refinement patches */
    ahf_gridinfo(*grid_list, curgrid_no);
    
    /* get AHF halos */
    ahf_halos(*grid_list);
    
#ifndef AHFstep
     /* restore linked-list structure */
    cur_part = global.fst_part;
    for(i=0; i<global.no_part; i++)
    {
      cur_part->ll = ll_tmp[i];
      cur_part++;
    }
#endif
  }
  
  /* reset ioflag */
  global.ioflag = FALSE;
  
#endif /* AHF */
  ahf.time += time(NULL);
  
#ifdef AHFstep
  return(cur_ms);
#endif
  
  
#ifdef VERBOSE
  fprintf(stderr,"\n        doing step on %ld grid with timestep = %g\n", cur_grid->l1dim, timestep);
#endif
  
  
  
  /*=========================================================================
   * STEP on current grid:
   *
   * 1. DRIFT/2     positions     (drift_pos, move_part, and treat leavers)
   * 2. get density fields right  (restore_dens, refill_dens, adjust_dens)
   * 3. solve for the potential   (solve_gravity)
   * 4. difference  potential     (calc_forces)
   * 5. KICK        velocities    (kick_mom, treat leavers)
   * 6. DRIFT/2     positions     (drift_pos, move_part, and treat leavers)
   *=========================================================================*/
  
  
   cur_grid->timecounter = timecounter;
   cur_grid->timestep    = timestep;
  
  cur_grid->time.hydro -= time(NULL);
#ifdef HYDRO
  /*--------------------------------------------------------------------
   *                    hydro-predictor step
   *--------------------------------------------------------------------*/
  global.a = calc_super_a(timecounter); // some of the dependent routines require knowledge about the current value of a (e.g. calc_p()!)
  calc_Flux(cur_grid);
  global.a = calc_super_a(timecounter + timestep); // all routines from now on need a(T+dT) !!
  advect_hydro(cur_grid, timestep);
  ave_hydro(cur_grid);
#endif /* HYDRO */
  cur_grid->time.hydro += time(NULL);
  
  /* get DRIFT and KICK operators */
  get_DKoperators(timestep, &KICK, &DRIFT1, &DRIFT2);
  
  /*--------------------------------------------------------------------
   *                         DRIFT/2 positions
   *--------------------------------------------------------------------*/
  cur_grid->time.DK -= time(NULL);
  drift_pos(cur_grid, DRIFT1);                  /* performs NULL_newll(cur_grid) */ 
  cur_grid->time.DK += time(NULL);

  /* update DRIFT and timecounter for that grid */
  cur_grid->timecounter += timestep/2;
  cur_grid->multistep++;
  
  /* move particles to new linked-lists */
  cur_grid->time.DK -= time(NULL);
  move_part(cur_grid);                          /* performs zero_dens(cur_grid) */
  cur_grid->time.DK += time(NULL);

#ifdef ADAPTIVE
  if(cur_grid != global.dom_grid)
  {
     /* relink edge particles back to next coarser grid */
    cur_grid->time.DK -= time(NULL);
    relink_bk_edge(cur_grid-1, cur_grid);

    /* take care of particles that crossed the boundary */
    unstep_leavers(cur_grid, DRIFT1, DRIFT2, KICK);
    cur_grid->time.DK += time(NULL);
    
     /* reassign P1/P2 density on coarser grid */
    cur_grid->time.density -= time(NULL);
    restore_dens(cur_grid-1);
    cur_grid->time.density += time(NULL);
  }      
#endif
  
  /*--------------------------------------------------------------------
   *                    obtain new density field
   *--------------------------------------------------------------------*/
  cur_grid->time.density -= time(NULL);
  /*------------
   * DM density
   *------------*/
  assign_dens(cur_grid);
  
#ifdef ADAPTIVE
  /* refill density from P3 particles (re-)linked to finer grid (fin_grid) */
  if(curgrid_no != fingrid_no) 
    refill_dens(cur_grid);
  
  /* account for overlapping density with next coarser grid */
  if(cur_grid != global.dom_grid)
    adjust_dens(cur_grid, FALSE);
#endif
  cur_grid->time.density += time(NULL);
  

  cur_grid->time.hydro -= time(NULL);
#ifdef HYDRO
  /*-----------------
   * add gas density
   *-----------------*/
  add_gasdens(cur_grid);
#endif
  cur_grid->time.hydro += time(NULL);
  
  
#ifdef CHECK_MASS
  fprintf(stderr,"\n                mass conservation check before entering solve_gravity:\n");
  
#ifndef HYDRO
  for(for_grid = cur_grid; for_grid >= global.dom_grid; for_grid--)
    restore_dens(for_grid);
  
  if(curgrid_no != fingrid_no) 
    refill_dens(cur_grid);
  
  if(cur_grid != global.dom_grid)
    stack_dens(cur_grid, global.dom_grid+1);
#endif
  
  for(for_grid=global.dom_grid; for_grid <= *grid_list+curgrid_no; for_grid++)
   {
     fprintf(stderr,"                  grid with length %ld: ",for_grid->l1dim);
     masscheck = sum_dens(for_grid);
     fprintf(stderr,"%f\n",masscheck);
   }
  fprintf(stderr,"                  ...done\n");
#endif
  
  
  /*--------------------------------------------------------------------
   *                         solve gravity
   *--------------------------------------------------------------------*/
  cur_grid->time.potential -= time(NULL);
  
  solve_gravity(*grid_list, curgrid_no);
  
  /*--------------------------------------------------------------------
   *                         obtain forces
   *--------------------------------------------------------------------*/
  calc_forces(cur_grid);
  
  cur_grid->time.potential += time(NULL);
  
  
  /*--------------------------------------------------------------------
   *                         KICK momenta
   *--------------------------------------------------------------------*/
  cur_grid->time.DK -= time(NULL);
  kick_mom(cur_grid, KICK);
  cur_grid->time.DK += time(NULL);
  
#ifdef ADAPTIVE
  /* take care of particles that already crossed the boundary */
  cur_grid->time.DK -= time(NULL);
  if(cur_grid != global.dom_grid && cur_grid->multistep == 3)
    advance_leavers(cur_grid, DRIFT2, KICK);
  cur_grid->time.DK += time(NULL);
#endif
  
  
  cur_grid->time.hydro -= time(NULL);
#ifdef HYDRO
  /*--------------------------------------------------------------------
   *                  hydro-corrector step + source terms
   *--------------------------------------------------------------------*/
  //global.a = calc_super_a(timecounter+timestep); // has already been set in predictor step after first calc_Flux() call
  calc_Flux(cur_grid);
  advect_gravity_hydro(cur_grid, timecounter, timestep);
#endif /* HYDRO */
  cur_grid->time.hydro += time(NULL);
  
  /*--------------------------------------------------------------------
   *                         DRIFT/2 positions
   *--------------------------------------------------------------------*/
  cur_grid->time.DK -= time(NULL);
  drift_pos(cur_grid, DRIFT2);                  /* performs NULL_newll(cur_grid) */
  cur_grid->time.DK += time(NULL);
  
  /* update DRIFT counter and timecounter for that grid */
  cur_grid->timecounter += timestep/2;
  cur_grid->multistep++;
  
  /* move particles to new linked-lists */
  cur_grid->time.DK -= time(NULL);
  move_part(cur_grid);                          /* performs zero_dens(cur_grid) */
  cur_grid->time.DK += time(NULL);

#ifdef ADAPTIVE
  if(cur_grid != global.dom_grid)
    {
      /* relink edge particles back to next coarser grid */
      cur_grid->time.DK -= time(NULL);
      relink_bk_edge(cur_grid-1, cur_grid);
      
      /* take care of particles that crossed the boundary */
      unstep_leavers(cur_grid, DRIFT1, DRIFT2, KICK);
      cur_grid->time.DK += time(NULL);
      
      /* reassign P1/P2 density on coarser grids */
      cur_grid->time.density -= time(NULL);
      restore_dens(cur_grid-1);
      cur_grid->time.density += time(NULL);
    }      
#endif  // ADAPTIVE  
  
  /* assign updated positions to check ratio nodes/particles */
  /* cur_grid->time.density -= time(NULL);
   * assign_dens(cur_grid);
   * cur_grid->time.density += time(NULL); */
  
  /*=========================================================================
   *
   *                 STEP on current grid has finished
   *
   *=========================================================================*/
  
  
  
  
/*----------------------------------------------------------------------------
   * in case the current grid is NOT the finest grid we need to make the 
   * finishing recursive call to step()
   * (we do not try to refine the finest grid again...check usage of fin_l1dim)
   *----------------------------------------------------------------------------*/
#ifdef ADAPTIVE
  if(refined == TRUE && fin_ms == TRUE)
  {
     /* finishing recursive call to STEP */
    fin_ms = step(grid_list, no_grids, timestep/2, timecounter+timestep/2);
    
    if(fin_ms == TRUE)
    {
        /* grid_list block might have moved in memory */
      global.dom_grid = *grid_list + global.domgrid_no;
      cur_grid        = *grid_list + curgrid_no;
      fin_grid        =  cur_grid  + 1;
      
        /*  relink particles back to parent grid */
      fin_grid->time.grid -= time(NULL);
      relink_back(cur_grid, fin_grid);
      fin_grid->time.grid += time(NULL);
      
#ifdef STORE_REFS
        /* store pquad pointer to access old hydro-variables next time */
      fin_grid->time.grid -= time(NULL);
      store_pquad(fin_grid, no_grids);
      fin_grid->time.grid += time(NULL);
#else
        /* finally destroy pquad/cquad/nquad structure */
      fin_grid->time.grid -= time(NULL);
      free_grid(fin_grid, no_grids);      /* performs no_grids-- */
      fin_grid->time.grid += time(NULL);
#endif
    }
  }
#endif
  
  return(cur_ms);
} 
