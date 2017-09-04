#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "utility.h"
#include "../libparticles/particles.h"
#include "../libgravity/gravity.h"
#include "../libgrids/grids.h"

/*
 * NOTE:  these routine do not distinguish between FAST and NORMAL
 */

/*====================================================================
* integrate Layzer-Irvine cosmic energy equation
*====================================================================*/
void calc_econst(double a_current)
{
   double a_initial, da, Told, Tnow, Cnow;
   
   a_initial = simu.a_initial;
   da        = a_current - energy.aold;
   
   Told = energy.Kold/pow2(energy.aold);
   Tnow = energy.K_current/pow2(a_current);
   
   energy.integral += (Told+Tnow)/2 * da;
   
#ifdef ISOLATED
   Cnow          = energy.K_current + energy.U_current;
   energy.echeck = fabs(Cnow - energy.econst)/ fabs(Cnow);
   energy.econst = Cnow;
#else
   energy.econst = (energy.K_current/a_current + energy.U_current) -
      (energy.K_initial/a_initial + energy.U_initial) + energy.integral;
   energy.echeck = (energy.econst/(energy.U_current));
#endif
   
   energy.aold   = a_current;
   
   return;
}


/*====================================================================
* calculate the kinetic energy of all the particles
*====================================================================*/
double kinetic(double a_current)
{
   double  total_T;       /* total kinetic energy */
   int     i;             /* loop index           */
   partptr cur_part;      /* current particle     */
   
   
   total_T = 0.0;
   
   
   /* loop over all the particles */
   for(cur_part=global.fst_part; cur_part<global.fst_part+global.no_part; cur_part++)
      for(i = X; i <= Z; i++)
#ifdef MULTIMASS
         total_T += (double)cur_part->weight*pow2((double)cur_part->mom[i]);
#else
         total_T += pow2((double)cur_part->mom[i]);
#endif
   
   total_T /= 2;
   
   return(total_T);
}


/*====================================================================
* calculate the potential energy on given grid
*====================================================================*/
double potential(gridls *grid_list, int lstgrid_no)
{
   pqptr   cur_pquad;
   cqptr   cur_cquad, icur_cquad;
   nqptr   cur_nquad, icur_nquad;
   nptr    cur_node;
   long    x, y, z;
   
   gridls *fst_grid, *lst_grid;
   gridls *for_grid;
   
   partptr cur_part;
   nptr    tsc_nodes[3][3][3];
   dvect   weights[3];
   dvect   xyz_coords;
   dvect   temp_coords;
   dvect   pn_sep;
   int     idim,grid_no;
   int     i,j,k;
   double  pnarg_a;
   double  pnarg_b;
   double  tpnarg;
   
   double  total_U, U, Up;
   double  Ui, dVolume, pot_mean;
   
   double  cur_shift;
  
   long    ipquad;
   
   /* make use of full grid hierarchy */
   fst_grid = grid_list + global.domgrid_no;
   lst_grid = grid_list + lstgrid_no;
   
   /* get density on finest grid */
   zero_dens(lst_grid);
   assign_dens(lst_grid);
   
   /* get density on all coarser grids */
#ifdef ADAPTIVE
   for(for_grid = lst_grid-1; for_grid >= global.dom_grid; for_grid--)
      restore_dens(for_grid);      
   
   if(lst_grid != global.dom_grid)
      stack_dens(lst_grid, global.dom_grid+1);
#endif /* ADAPTIVE */
   
   /* fill timecounter with correct super-comoving time */
   for(for_grid = lst_grid; for_grid >= global.dom_grid; for_grid--)
      for_grid->timecounter = calc_super_t(global.a);         

   /* solve for potential on ALL grids */
#ifdef VERBOSE2
      fprintf(stderr, "potential:  calling solve_dom_gravity.\n");
#endif
   solve_dom_gravity(grid_list);
  
#ifdef ADAPTIVE
   if(lstgrid_no > global.domgrid_no)
      for(grid_no = global.domgrid_no+1; grid_no <= lstgrid_no; grid_no++)
      {
	  for_grid = grid_list+grid_no;
          go_down(for_grid-1);
#ifdef VERBOSE
         fprintf(stderr, "potential:  calling solve_dom_gravity.\n");
#endif
         solve_ref_gravity(grid_list, grid_no);
      }
#endif  // ADAPTIVE
  
   /* 1. zero temporary storage */
   for(for_grid = fst_grid; for_grid <= lst_grid; for_grid++)
      
      for(cur_pquad=for_grid->pquad; cur_pquad!=NULL; cur_pquad=cur_pquad->next)
         
         for(cur_cquad = cur_pquad->loc, z = cur_pquad->z; 
             cur_cquad < cur_pquad->loc + cur_pquad->length;
             cur_cquad++, z++)
            for(icur_cquad=cur_cquad; icur_cquad!=NULL; icur_cquad=icur_cquad->next)
               
               for(cur_nquad = icur_cquad->loc, y = icur_cquad->y;
                   cur_nquad < icur_cquad->loc + icur_cquad->length; 
                   cur_nquad++, y++)
                  for(icur_nquad=cur_nquad;icur_nquad!=NULL;icur_nquad=icur_nquad->next)
                     
                     for(cur_node = icur_nquad->loc, x = icur_nquad->x;
                         cur_node < icur_nquad->loc + icur_nquad->length; 
                         cur_node++, x++)
                       {
                        /* temporary storage used during volume integration */
                        cur_node->force.temp[0] = 0.0;
                        cur_node->force.temp[1] = 0.0;
                        cur_node->force.temp[2] = 0.0;
                       }
                        
                        
   /* 2. get mean potential value */
   pot_mean = (double)0.0;
   for(for_grid = lst_grid; for_grid >= fst_grid; for_grid--)
     {
      for(cur_pquad=for_grid->pquad; cur_pquad!=NULL; cur_pquad=cur_pquad->next)
         
         for(cur_cquad = cur_pquad->loc, z = cur_pquad->z; 
             cur_cquad < cur_pquad->loc + cur_pquad->length;
             cur_cquad++, z++)
            for(icur_cquad=cur_cquad; icur_cquad!=NULL; icur_cquad=icur_cquad->next)
               
               for(cur_nquad = icur_cquad->loc, y = icur_cquad->y;
                   cur_nquad < icur_cquad->loc + icur_cquad->length; 
                   cur_nquad++, y++)
                  for(icur_nquad=cur_nquad;icur_nquad!=NULL;icur_nquad=icur_nquad->next)
                     
                     for(cur_node = icur_nquad->loc, x = icur_nquad->x;
                         cur_node < icur_nquad->loc + icur_nquad->length; 
                         cur_node++, x++)
                       {
                        /* take node contribution ? */
                        if(cur_node->force.temp[0] < 1.9)
                          {
                           if(cur_node->force.temp[0] < 0.9)
                              dVolume = (double) 1.0 / for_grid->masstodens;
                           else
                              dVolume = (double) cur_node->force.temp[1];
                           
                           /* dVolume already includes the factor 
                            * simu.no_part/simu.no_vpart (via masstodens) !!! */
                           pot_mean += (double)cur_node->pot*dVolume
#ifdef MULTIMASS
                              /(double)simu.no_vpart;
#else
                              /(double)simu.no_part;
#endif
                          }
                        
                        /* zero temporary storage again... */
                        cur_node->force.temp[0] = 0.0;
                        cur_node->force.temp[1] = 0.0;
                       }
                        
      /* get proper volume for each cell on next coarser level */
      if(for_grid != fst_grid)
         f2c_volume(for_grid);  
     }
 
   /* 3. get Ui by volume integration */
   Ui = (double)0.0;
   for(for_grid = lst_grid; for_grid >= fst_grid; for_grid--)
     {
      for(cur_pquad=for_grid->pquad; cur_pquad!=NULL; cur_pquad=cur_pquad->next)
         
         for(cur_cquad = cur_pquad->loc, z = cur_pquad->z; 
             cur_cquad < cur_pquad->loc + cur_pquad->length;
             cur_cquad++, z++)
            for(icur_cquad=cur_cquad; icur_cquad!=NULL; icur_cquad=icur_cquad->next)
               
               for(cur_nquad = icur_cquad->loc, y = icur_cquad->y;
                   cur_nquad < icur_cquad->loc + icur_cquad->length; 
                   cur_nquad++, y++)
                  for(icur_nquad=cur_nquad;icur_nquad!=NULL;icur_nquad=icur_nquad->next)
                     
                     for(cur_node = icur_nquad->loc, x = icur_nquad->x;
                         cur_node < icur_nquad->loc + icur_nquad->length; 
                         cur_node++, x++)
                       {
#ifndef ISOLATED
                        /* subtract mean potential from *all* nodes */
                        cur_node->pot -= pot_mean;
#endif

                        /* take node contribution ? */
                        if(cur_node->force.temp[0] < 1.9)
                          {
                           if(cur_node->force.temp[0] < 0.9)
                              dVolume = (double) 1.0 / for_grid->masstodens;
                           else
                              dVolume = (double) cur_node->force.temp[1];
                           
#ifdef ISOLATED
                           Ui += ((double)cur_node->dens) * 
                              (double)cur_node->pot * dVolume; 
#else
                           Ui += ((double)cur_node->dens + simu.mean_dens) * 
                              (double)cur_node->pot * dVolume; 
#endif
                          }
                       }

      if(for_grid != fst_grid)
         f2c_volume(for_grid);  
     }
   
   /* 4. get Up by summing over all particles */
   Up = (double)0.0;
   for(for_grid = fst_grid; for_grid <= lst_grid; for_grid++)
     {
      
      /* shift of cell centre as compared to edge of box [grid units] */
      cur_shift = 0.5/(double)for_grid->l1dim;

      for(cur_pquad=for_grid->pquad; cur_pquad!=NULL; cur_pquad=cur_pquad->next)
         
         for(cur_cquad = cur_pquad->loc, z = cur_pquad->z; 
             cur_cquad < cur_pquad->loc + cur_pquad->length;
             cur_cquad++, z++)
            for(icur_cquad=cur_cquad; icur_cquad!=NULL; icur_cquad=icur_cquad->next)
               
               for(cur_nquad = icur_cquad->loc, y = icur_cquad->y;
                   cur_nquad < icur_cquad->loc + icur_cquad->length; 
                   cur_nquad++, y++)
                  for(icur_nquad=cur_nquad;icur_nquad!=NULL;icur_nquad=icur_nquad->next)
                     
                     for(cur_node = icur_nquad->loc, x = icur_nquad->x;
                         cur_node < icur_nquad->loc + icur_nquad->length; 
                         cur_node++, x++)
                       {
                        if(cur_node->ll != NULL)
                          {
                           tsc_nodes[1][1][1] = cur_node;
                           
                           temp_coords[X] =(((double)x)/(double)for_grid->l1dim)+cur_shift;
                           temp_coords[Y] =(((double)y)/(double)for_grid->l1dim)+cur_shift;
                           temp_coords[Z] =(((double)z)/(double)for_grid->l1dim)+cur_shift;
                           
                           xyz_coords[X]  = f1mod(temp_coords[X]+1.0, 1.0);
                           xyz_coords[Y]  = f1mod(temp_coords[Y]+1.0, 1.0);
                           xyz_coords[Z]  = f1mod(temp_coords[Z]+1.0, 1.0);
                           
                           /* get pointers to tsc nodes */
                           get_TSCnodes(for_grid,cur_pquad,icur_cquad,icur_nquad,
                                        tsc_nodes, &z, &y, &x);
                           
                           if(test_tsc(tsc_nodes) == TRUE)
                             {
                              for(cur_part  = cur_node->ll;
                                  cur_part != NULL;
                                  cur_part  = cur_part->ll)
                                {
                                 for(idim = 0; idim < NDIM; idim++)
                                   {
                                    pn_sep[idim]=
                                    ((double)cur_part->pos[idim]-xyz_coords[idim])
                                    * (double)for_grid->l1dim;
                                    if(fabs(pn_sep[idim]) > 0.5*(double)for_grid->l1dim)
                                      {
                                       tpnarg       = 
                                       (double)cur_part->pos[idim]+0.5;
                                       pnarg_a      = f1mod(tpnarg+1.0, 1.0);
                                       tpnarg       = xyz_coords[idim] + 0.5;
                                       pnarg_b      = f1mod(tpnarg+1.0, 1.0);
                                       pn_sep[idim] =(pnarg_a-pnarg_b)*
                                          (double)for_grid->l1dim;
                                      }
#ifdef TSC
                                    weights[1][idim] = 0.75-pow2(pn_sep[idim]);
                                    weights[0][idim] = pow2((0.5-pn_sep[idim]))/2;
                                    weights[2][idim] = pow2((0.5+pn_sep[idim]))/2;
#endif
#ifdef CIC
                                    if(pn_sep[idim] > 0.)
                                      {
                                       weights[0][idim] = 0.0;
                                       weights[1][idim] = 1.0 - pn_sep[idim];
                                       weights[2][idim] =       pn_sep[idim];
                                      }
                                    else
                                      {
                                       weights[0][idim] =     - pn_sep[idim];
                                       weights[1][idim] = 1.0 + pn_sep[idim];
                                       weights[2][idim] = 0.0;
                                      }
#endif
#ifdef NGP
                                    weights[0][idim] = 0.0;
                                    weights[1][idim] = 1.0;
                                    weights[2][idim] = 0.0;
#endif
                                   }
                                 
                                 for(k = 0; k < 3; k++)
                                    for(j = 0; j < 3; j++)
                                       for(i = 0; i < 3; i++)
                                         {
                                          Up += ((double)
                                                 tsc_nodes[k][j][i]->pot * 
#ifdef MULTIMASS
                                                 (double)cur_part->weight*
#endif
                                                 weights[k][Z]*
                                                 weights[j][Y]*
                                                 weights[i][X]);
                                         }
                                }
                             }
                          }
                       }
     }
                        
      /*
       * even though we calculated the potential energy U using two different 
       * approaches (volume integration and summation over all particles)
       * we only use the particle sum as this is more accurate...
       * you might want to compare both results and have the chance to do so
       * right now, right at this place of the code.
       *
       * anyway, in case you feel like removing the volume integration
       * keep in mind that we are simultaneously correcting the mean
       * value of the potential which should always be zero but gradually
       * drifts away from that mean...
       */
      
      /* take mean value of integral and summation potential energy */
#ifdef UPOTMEAN
      U = (Ui+Up)/2;
#else
      U = Up;
#endif
      
      /* get units right */
#ifdef ISOLATED
      total_U = 0.5 * U;
#else
      total_U = 0.5 * U / global.a;
#endif /* ISOLATED */
      
      return (total_U);
}

#ifdef ISOLATED
/*==============================================================================
* simple energy conservation for isolated boundary conditions
*==============================================================================*/
void layzer_irvine(double timecounter, gridls *grid_list, int fingrid_no)
{
   double a_current;
   
   a_current = 1.0;
   
   /* get kinetic energy by looping over all particles */
   energy.Kold      = energy.K_current;
   energy.K_current = kinetic(a_current);
   
   /* get potential energy using all available grids */
   energy.Uold      = energy.U_current;
   energy.U_current = potential(grid_list, fingrid_no);
   
   /* check Layzer-Irvine energy conservation */
   calc_econst(a_current);
}
#else
/*==============================================================================
* Layzer-Irvine energy check control routine
*==============================================================================*/
void layzer_irvine(double timecounter, gridls *grid_list, int fingrid_no)
{
   double a_current;
   
   a_current = calc_super_a(timecounter);
   
   global.a = a_current;
   
   /* get kinetic energy by looping over all particles */
   energy.Kold      = energy.K_current;
   energy.K_current = kinetic(a_current);
   
   /* get potential energy using all available grids */
   energy.Uold      = energy.U_current;
   energy.U_current = potential(grid_list, fingrid_no);
   
   /* check Layzer-Irvine energy conservation */
   calc_econst(a_current);
}
#endif
