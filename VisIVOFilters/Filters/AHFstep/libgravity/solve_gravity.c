#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "gravity.h"
#include "../libutility/utility.h"
#include "../libgrids/grids.h"

/*
 * there are two ways to stop the GS iteration procedure:
 *
 * 1. the residuals are smaller than the truncation error
 *
 * 2. the residuals are smaller than LIMIT given below
 *
 */
#define LIMIT     0.00001

#define ONE_THIRD 0.33333333333333333333

#include "gravity.h"

#ifdef FFT
#define COA_FFT
#endif

/*==============================================================================
* solve on coarsest grid
*==============================================================================*/
void solve_cg(gridls *cur_grid)
{
   
#ifdef COA_FFT
   flouble *dens_array;         /* density array pointer */
   pqptr    cur_pquad;          /* current pquad         */
   cqptr    cur_cquad;          /* current cquad         */
   nqptr    cur_nquad;          /* current nquad         */
   nptr     cur_node;           /* current node          */
   long     i, j, k, l1dim, FFTarray_length;
   double   FourPiGa;

   /* conversion factor to go from densito to source term */
   FourPiGa = simu.FourPiG*calc_super_a(cur_grid->timecounter);

   /* array dimension */
   l1dim           = cur_grid->l1dim;
   FFTarray_length = 2*l1dim*l1dim*l1dim;
   
   /* generate complex (!) density array for FFT */
   if((dens_array = (flouble *) calloc(FFTarray_length, sizeof(flouble))) == NULL)
     {
      fprintf(io.logfile,"solve_cg: could not allocate density array for FFT\n");
      fflush(io.logfile);
      fclose(io.logfile);
      exit(1);
     }
         
   /* fill density array for FFT ... no need for quad-ll's !! */
   cur_pquad = cur_grid->pquad;
   for(k = 0, cur_cquad = cur_pquad->loc; k < l1dim; k++, cur_cquad++)
      for(j = 0, cur_nquad = cur_cquad->loc; j < l1dim; j++, cur_nquad++)
         for(i = 0, cur_node = cur_nquad->loc; i < l1dim; i++, cur_node++)
           {
            dens_array[Re(i,j,k,l1dim)] = cur_node->dens * FourPiGa;  /* real part      */
            dens_array[Im(i,j,k,l1dim)] = 0.0;                         /* imaginary part */
           }
            
            /* solve by FFT */
            fft_potential(dens_array, l1dim);
   
   /* fill node potential values */
   for(k = 0, cur_cquad = cur_pquad->loc; k < cur_grid->l1dim; k++, cur_cquad++)
      for(j = 0, cur_nquad = cur_cquad->loc; j < cur_grid->l1dim; j++, cur_nquad++)
         for(i = 0, cur_node = cur_nquad->loc; i < cur_grid->l1dim; i++,cur_node++)
            cur_node->pot = dens_array[Re(i,j,k,l1dim)];
   
   /* destroy memory assigned to dens_array */
   free(dens_array);
   
#else   /* for this case simu.NGRID_MIN needs to be 2 !!!! */
   
   nptr    node_arr[2][2][2];       /* 2x2x2 array of pointers to nodes */ 
  double node_dens[3][3][3];
   int   i,j,k;                   /* loop indices x,y,z               */
   int   not_i, not_j, not_k;     /* indices opp to i, j & k          */
   pqptr cur_pquad;               /* current pquad                    */
   cqptr cur_cquad;               /* current cquad                    */
   nqptr cur_nquad;               /* current nquad                    */
   nptr  cur_node;                /* current node                     */
   double  FourPiGa;
   
   /* conversion factor to go from density to source term */
   FourPiGa = simu.FourPiG*calc_super_a(cur_grid->timecounter);
   
   /* find pointers to all the nodes */
   cur_pquad = cur_grid->pquad;
   cur_cquad = cur_pquad->loc;
   for(k = 0, cur_cquad = cur_pquad->loc; k < 2; k++, cur_cquad++)
      for(j = 0, cur_nquad = cur_cquad->loc; j < 2; j++, cur_nquad++)
         for(i = 0, cur_node = cur_nquad->loc; i < 2; i++, cur_node++)
         {
           node_arr[i][j][k]  = cur_node;
           node_dens[i][j][k] = cur_node->dens * FourPiGa;
         }

   /* loop over grid */
   for(i = 0; i < 2; i++)
     {
      not_i = set_opposite(i);
      for(j = 0; j < 2; j++)
        {
         not_j = set_opposite(j);
         for(k = 0; k < 2; k++)
           {
            not_k = set_opposite(k);
             node_arr[i][j][k]->pot = ((double)cur_grid->spacing2 / 48.0) *
               ((7 * ((double)node_arr[i][not_j][k]->dens +
                      (double)node_arr[not_i][j][k]->dens +
                      (double)node_arr[i][j][not_k]->dens)) +
                (9 * ((double)node_arr[not_i][not_j][k]->dens +
                      (double)node_arr[i][not_j][not_k]->dens +
                      (double)node_arr[not_i][j][not_k]->dens)) + 
                (10 * (double)node_arr[not_i][not_j][not_k]->dens));
           }
        }
     }
#endif
}


/*============================================================================
* test convergence of current grid
*============================================================================*/
boolean converged(gridls *cur_grid)
{
   double trunc_error;
   double residual;
   
   residual     =             cur_grid->cur_resid;
   trunc_error  = ONE_THIRD * cur_grid->trunc_err;
   
   if     (residual <= LIMIT)
      return TRUE;
   else if(residual <= (trunc_error * CONVCRIT)) 
      return TRUE;
   else
      return FALSE;
}

/*=============================================================================
* test for slow convergence
*=============================================================================*/
boolean slow_conv(gridls *cur_grid)
{
#ifdef NO_MULTIGRID
  return FALSE;
#endif
   if(cur_grid->cur_resid > (ETA * cur_grid->old_resid))
      return TRUE;
   else
      return FALSE;
}

/*=============================================================================
* solve for potential on domain grids
*=============================================================================*/
void solve_dom_gravity(gridls *grid_list)
{
  int     grid_no;             /* the current grid number     */
  int     lstgrid_no;          /* no_grids minus one          */
  int     i;                   /* index for GS sweep loop     */
  gridls *cur_grid;            /* pointer to the current grid */
  boolean cg_conv;             /* has current grid converged? */
  int     gs_sweeps;
  
  gs_sweeps = DOMSWEEPS;
  
#ifdef FFT
  solve_cg(global.dom_grid);
  return;
#endif

  /* last grid to be taken into account */
  lstgrid_no = global.domgrid_no;
  
   /* zero residuals on all multi-grids and set timecounter accordingly */
  for(cur_grid = grid_list; cur_grid <= grid_list + lstgrid_no; cur_grid++)
    {
     cur_grid->cur_resid   = 0.;
     
     /* all Multi-Grid grids need to be synchronized with dom_grid as conversion
      * from density to source term depends on super_t */
     cur_grid->timecounter = global.dom_grid->timecounter;
    }
  
   /* prepare to start with finest domain grid */
  cur_grid = global.dom_grid;
  grid_no  = global.domgrid_no;
  cg_conv  = FALSE;
  
  cur_grid->no_sweeps = 0;
  
   /* loop over the grid until converged on finest grid */
  while(!((cg_conv == TRUE) && (grid_no == lstgrid_no)))
  {
    if(grid_no == 0)   /* are we on the coarsest grid ? */
    {
      solve_cg(cur_grid);
      cur_grid = go_down(cur_grid);
      grid_no++;
    }
    else  /* grid_no == 0 */
    {
       /* keep old residual in mind */
       cur_grid->old_resid = cur_grid->cur_resid;   
       
       /* do GS relaxation sweep */
      for(i = 0; i < gs_sweeps; i++)
        gs(cur_grid);
      
      cur_grid->no_sweeps += gs_sweeps;
      cur_grid->cur_resid  = residual(cur_grid);     /* calc new residual */
      cur_grid->trunc_err  = trunc_err(cur_grid);         /* calc trunc. error */
            
      cg_conv = converged(cur_grid);
      
        /* converged */
      if(cg_conv == TRUE)
      {
        if(grid_no != lstgrid_no)   /* are we on the finest grid ?   */
        {
          grid_no++;
          cur_grid            = go_down(cur_grid);
          cur_grid->no_sweeps = 0;
          cg_conv             = FALSE;
          gs_sweeps           = DOMSWEEPS;
        }
      }
      /* not converged but rather stuck */
      else if(fabs(cur_grid->old_resid - cur_grid->cur_resid) < ZERO)
      {
        gs_sweeps += DOMSWEEPS;   /* just in case we are somehow stuck */
      }
      /* slow convergence */
      else if(slow_conv(cur_grid) == TRUE)
      {
        if((grid_no != 0))        /* are we on the coarsest grid ?     */
        {
          grid_no--;
          cur_grid            = go_up(cur_grid);
          cur_grid->no_sweeps = 0;
        }
      }
    }
  }
  
   /* do GS sweeps one more time */
  for(i = 0; i < DOMSWEEPS; i++)
    gs(cur_grid);

  cur_grid->cur_resid = residual(cur_grid);
  
}


/*==============================================================================
* solve for potential on refinements
*==============================================================================*/
void solve_ref_gravity(gridls *grid_list, int grid_no)
{
   int     i;                   /* index for GS sweep loop     */
   gridls *cur_grid;            /* pointer to the current grid */
   boolean cg_conv;             /* has current grid converged? */
   int     gs_sweeps;
   
   /* set cur_grid */
   cur_grid = grid_list+grid_no;
   
   /* get boundary values from next coarser level */
   if(cur_grid->multistep == 3)
      go_down(grid_list+grid_no-1);
   
   cur_grid->no_sweeps  = 0;
   cur_grid->cur_resid  = 0.;
   cg_conv              = FALSE;
   gs_sweeps            = REFSWEEPS;
  
   /* do GS sweeps until converged */
   while(cg_conv == FALSE)
   {
      
     /* do GS relaxation sweep */
     for(i = 0; i < gs_sweeps; i++)
       gs(cur_grid);
     
     cur_grid->no_sweeps += REFSWEEPS;
     
     /* calculate current residual (stored for each node in force.temp[0]) */
     cur_grid->cur_resid  = residual(cur_grid);
     
     
#ifdef FIXED_SWEEPS
     cg_conv = TRUE;
#else
     /* calculate current truncation error (uses force.temp[1] and force.temp[2])  */
     cur_grid->trunc_err = trunc_err(cur_grid);
          
     /* check for convergence, i.e. compare residual against truncation error */
     cg_conv             = converged(cur_grid);
     
     
#ifdef SWEEP_TEST
     fprintf(stderr,"grid=%10ld   resid=%16.8g  trunc_err=%16.8g   -> cg_conv=%10d\n",
             cur_grid->l1dim, cur_grid->cur_resid, CONVCRIT*ONE_THIRD*cur_grid->trunc_err, cg_conv);
#endif
     
#endif
     
     if(fabs(cur_grid->old_resid - cur_grid->cur_resid) < ZERO)
       gs_sweeps += REFSWEEPS;   /* just in case we are somehow stuck */
     
   }   /* cg_conv == FALSE */
  
  /* do GS sweeps one more time */
  for(i = 0; i < REFSWEEPS; i++)
    gs(cur_grid);
  
  cur_grid->cur_resid = residual(cur_grid);

}


/*==============================================================================
* control routine for multi grid solver
*==============================================================================*/
void solve_gravity(gridls *grid_list, int curgrid_no)
{
#ifdef FORCETEST
   int     grid_no;
   gridls *cur_grid;
   
   fprintf(stderr,"\n  solving on      %7ld grid   ",global.dom_grid->l1dim);
   solve_dom_gravity(grid_list);
   
   for(grid_no = global.domgrid_no+1; grid_no <= curgrid_no; grid_no++)
     {
      cur_grid = global.dom_grid+(grid_no-global.domgrid_no);

      fprintf(stderr,"  going down from %7ld grid\n",(cur_grid-1)->l1dim);
      go_down(cur_grid-1);
      
      fprintf(stderr,"  solving on      %7ld grid   ",cur_grid->l1dim);
      solve_ref_gravity(grid_list, grid_no);
     }
   
#else /* FORCETEST */
   if(curgrid_no == global.domgrid_no)
      solve_dom_gravity(grid_list);
   else
      solve_ref_gravity(grid_list, curgrid_no);
#endif /* FORCETEST */
}


