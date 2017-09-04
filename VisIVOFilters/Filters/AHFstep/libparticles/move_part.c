#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef WITH_OPENMP
#include <omp.h>
#endif

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "particles.h"
#include "../libutility/utility.h"
#include "../libgrids/grids.h"

/*-------------------------------------------------------------------------------
* NOTE:
*       whenever cur_part->pos[] is modified call move_part() afterwards !
*
*       but never use move_part() as a stand-alone:
*       it requires a call to NULL_newll() beforehand.
* 
*       however, drift_pos() performs a NULL_newll() automatically
*       and hence there is no need for that call within step()...
*
*-------------------------------------------------------------------------------*/

/*===============================================================================
* rearrange linked list on cur_grid
*===============================================================================*/
void move_part(gridls *cur_grid)
{
   long    ipquad;
   pqptr   cur_pquad;
   cqptr   cur_cquad, icur_cquad;
   nqptr   cur_nquad, icur_nquad;
   nptr    cur_node;
   long    x,y,z;
   nptr    sur_nodes[3][3][3];
   partptr cur_part, next_part;
   dvect   xyz_coords;
   dvect   temp_coords;
   int     new_loc[NDIM];
   int     i;
   double  pnarg_a;
   double  pnarg_b;
   double  tpnarg;
   dvect   pn_sep;
   
   double cur_shift;
   
   /* shift of cell centre as compared to edge of box [grid units] */
   cur_shift = 0.5/(double)cur_grid->l1dim;
   
   cur_grid->size.no_nodes = 0;
   cur_grid->size.no_part  = 0;
   
   /* loop over all nodes -> particles */
   for(cur_pquad=cur_grid->pquad; cur_pquad!=NULL; cur_pquad=cur_pquad->next)
       
      for(cur_cquad = cur_pquad->loc, z = cur_pquad->z; 
          cur_cquad < cur_pquad->loc + cur_pquad->length;
          cur_cquad++, z++)
         for(icur_cquad=cur_cquad; icur_cquad!=NULL; icur_cquad=icur_cquad->next)
            
            for(cur_nquad = icur_cquad->loc, y = icur_cquad->y;
                cur_nquad < icur_cquad->loc + icur_cquad->length; 
                cur_nquad++, y++)
               for(icur_nquad=cur_nquad; icur_nquad!=NULL; icur_nquad=icur_nquad->next)
                  
                  for(cur_node = icur_nquad->loc, x = icur_nquad->x;
                      cur_node < icur_nquad->loc + icur_nquad->length; 
                      cur_node++, x++, cur_grid->size.no_nodes++)
                     
                     if(cur_node->ll != NULL)
                       {
                        sur_nodes[1][1][1] = cur_node;   /* current node in centre */
                        
                        /* calculate realspace coords of cur_node */
                        temp_coords[X] = (((double)x)/(double)cur_grid->l1dim)+cur_shift;
                        temp_coords[Y] = (((double)y)/(double)cur_grid->l1dim)+cur_shift;
                        temp_coords[Z] = (((double)z)/(double)cur_grid->l1dim)+cur_shift;
                        xyz_coords[X]  = f1mod(temp_coords[X]+1.0, 1.0);
                        xyz_coords[Y]  = f1mod(temp_coords[Y]+1.0, 1.0);
                        xyz_coords[Z]  = f1mod(temp_coords[Z]+1.0, 1.0);
                        
                        /* get pointers to surrounding nodes */
                        get_TSCnodes(cur_grid, cur_pquad, icur_cquad, icur_nquad, 
                                     sur_nodes, &z, &y, &x);
                        
                        /* loop over (old) ll */
                        for(cur_part=cur_node->ll; cur_part!=NULL; cur_part=next_part)
                          {
                           /* count particles on this grid */
                           cur_grid->size.no_part++;
                           
                           /* remember next particle because add_part_to_ll modifies cur_part */
                           next_part = cur_part->ll;
                           
                           /* calc new location for cur_part */
                           for(i = X; i <= Z; i++)
                             {
                              pn_sep[i] = (double)cur_grid->l1dim * 
                              ((double)cur_part->pos[i]-xyz_coords[i]);
                              
                              if(fabs(pn_sep[i]) > 0.5*(double)cur_grid->l1dim)
                                {
                                 tpnarg    = (double)cur_part->pos[i] + 0.5;
                                 pnarg_a   = f1mod(tpnarg+1.0, 1.0);
                                 tpnarg    = xyz_coords[i]            + 0.5;
                                 pnarg_b   = f1mod(tpnarg+1.0, 1.0);
                                 pn_sep[i] = (pnarg_a-pnarg_b)*(double)cur_grid->l1dim;
                                }
                              
                              /* get index for sur_node[][][] array */
                              new_loc[i] = (int) floor((pn_sep[i] + 1.5));
                             }
                           
                           if((new_loc[X] < 0) || (new_loc[X] > 2) || 
                              (new_loc[Y] < 0) || (new_loc[Y] > 2) || 
                              (new_loc[Z] < 0) || (new_loc[Z] > 2))
                             {
                              /* 1. particle moved further than TSC neighbors */
#ifdef VERBOSELOG2
                              fprintf(io.logfile,
                                      "move_part: particle %ld moved further than TSC neighbors on %ld grid\n",
                                      cur_part-global.fst_part,cur_grid->l1dim);
                              fflush(io.logfile);
#endif
                              extended_llsearch(cur_part, cur_grid);
                             }
                           else
                             {
                              /* 2. particle moved to one of the TSC neighbors */
#ifdef VERBOSELOG2
                              fprintf(io.logfile,
                                      "move_part: particle %ld moved to one of the TSC neighbors on %ld grid.  (%d,%d,%d)\n",
                                      cur_part-global.fst_part,cur_grid->l1dim, new_loc[X], new_loc[Y], new_loc[Z]);
                              fflush(io.logfile);
#endif
                              if(sur_nodes[new_loc[Z]][new_loc[Y]][new_loc[X]]!=NULL)
                                 add_part_to_ll(sur_nodes[new_loc[Z]][new_loc[Y]][new_loc[X]], cur_part);
                              
                              /* 3. particle definitely crossed refinement boundary */
                              else
                                {
#ifdef PERIODIC
#ifdef VERBOSELOG2
                                 fprintf(io.logfile,
                                         "move_part: particle %ld is a leaver\n");    
#endif
                                 store_leaver(cur_part, cur_grid);
#else
                                 /* only store_leaver, if particle is still within computational domain... */
                                 
                                 /* (to be coded) */
                                 
                                 
                                 /* ...otherwise we need to remove the particle somehow! */
                                 
                                 /* (to be coded) */
                                 
                                 fprintf(io.logfile,"\n\n particle left grid\n\n but you are running an isolated simulation for which proper boundary conditions have not been coded yet :-(\n\n");
#endif
                                }
                             }		      
                          }
                       }
                        
                        
   /* transfer force.new_ll to ll and perform zero_dens(cur_grid) */
#ifdef WITH_OPENMP
#pragma omp parallel private(ipquad, cur_pquad, cur_cquad, icur_cquad, cur_nquad, icur_nquad, cur_node) shared(cur_grid)
#pragma omp for schedule(static)
   for(ipquad=0; ipquad<cur_grid->no_pquad; ipquad++)
     {
      cur_pquad = cur_grid->pquad_array[ipquad];
#else
   for(cur_pquad=cur_grid->pquad;cur_pquad!=NULL;cur_pquad=cur_pquad->next)
     {
#endif
      
      for(cur_cquad=cur_pquad->loc;cur_cquad<cur_pquad->loc+cur_pquad->length; 
          cur_cquad++)
         for(icur_cquad=cur_cquad;icur_cquad!=NULL;icur_cquad=icur_cquad->next)
            
            for(cur_nquad=icur_cquad->loc;cur_nquad<icur_cquad->loc+icur_cquad->length; 
                cur_nquad++)
               for(icur_nquad=cur_nquad;icur_nquad!=NULL;icur_nquad=icur_nquad->next)
                  
                  for(cur_node=icur_nquad->loc;
                      cur_node<icur_nquad->loc+icur_nquad->length;cur_node++)
                    {
                     cur_node->ll   = cur_node->force.new_ll;
                     
                     /* take the chance to reset the density */
                     cur_node->dens = -simu.mean_dens;
                     
                     /* make sure to use union force as floats from now on */
                     cur_node->force.temp[0] = 0.0;
                     cur_node->force.temp[1] = 0.0;
                     cur_node->force.temp[2] = 0.0;
                    }
     }
}



