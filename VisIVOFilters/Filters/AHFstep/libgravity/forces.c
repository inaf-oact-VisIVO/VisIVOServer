#include <stddef.h>
#ifdef WITH_OPENMP
#include <omp.h>
#endif

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "../libutility/utility.h"
#include "../libgrids/grids.h"

/*==============================================================================
* forces_on_boundary: 
*
* interpolating the forces from a coarse to a fine grid.
* This becomes necessary for the boundary nodes of each fine grid...
*==============================================================================*/
void forces_on_boundary(gridls *coa_grid, gridls *fin_grid)
{
   long          i, j, k, ipquad;
   
   pqptr         coa_pquad, tcoa_pquad;
   cqptr         coa_cquad, tcoa_cquad;
   nqptr         coa_nquad, tcoa_nquad;
   nptr          coa_node;
   long          coa_x, coa_y, coa_z;
   nptr          tsc_coanodes[3][3][3];
   
   pqptr         fin_pquad;
   cqptr         fin_cquad, ifin_cquad;
   nqptr         fin_nquad, ifin_nquad;
   nptr          fin_node;
   long          fin_x, fin_y, fin_z;
   nptr          tsc_finnodes[3][3][3];
   
   double        fc_sep[NDIM], slope[NDIM], func[NDIM][NDIM][NDIM];
   
   /*==============================================================
    * loop over fine grid (and simultanesouly over coarse grid...)
    *==============================================================*/
#ifdef WITH_OPENMP
#pragma omp parallel private(ipquad, i, j, k, coa_pquad, tcoa_pquad, coa_cquad, tcoa_cquad, coa_nquad, tcoa_nquad, coa_node, coa_x, coa_y, coa_z, tsc_coanodes, fin_pquad, fin_cquad, ifin_cquad, fin_nquad, ifin_nquad, fin_node, fin_x, fin_y, fin_z, tsc_finnodes, fc_sep, slope, func) shared(fin_grid)
#pragma omp for schedule(static)
   for(ipquad=0; ipquad<fin_grid->no_pquad; ipquad++)
     {
      fin_pquad = fin_grid->pquad_array[ipquad];
#else
   for(fin_pquad=fin_grid->pquad; fin_pquad != NULL; fin_pquad=fin_pquad->next)
     {
#endif
      fin_z = fin_pquad->z;
      coa_z = fin_z/2;
      
      /* find correct coa_pquad */
      for(tcoa_pquad = coa_grid->pquad; coa_z > tcoa_pquad->z + tcoa_pquad->length; tcoa_pquad = tcoa_pquad->next)
         ;
      
      /* jump to correct cquad */
      coa_cquad = tcoa_pquad->loc + (coa_z - tcoa_pquad->z);
      
      
      for(fin_cquad = fin_pquad->loc;
          fin_cquad < fin_pquad->loc + fin_pquad->length; 
          fin_cquad++, fin_z++)  
        {  
         for(ifin_cquad  = fin_cquad; 
             ifin_cquad != NULL; 
             ifin_cquad  = ifin_cquad->next)
           {
            fin_y = ifin_cquad->y;
            coa_y = fin_y/2;
            
            /* find correct coa_cquad */
            for(tcoa_cquad = coa_cquad; coa_y > tcoa_cquad->y + tcoa_cquad->length; tcoa_cquad = tcoa_cquad->next)
               ;
            
            /* jump to correct nquad */
            coa_nquad = tcoa_cquad->loc + (coa_y - tcoa_cquad->y);
            
            for(fin_nquad = ifin_cquad->loc;  
                fin_nquad < ifin_cquad->loc + ifin_cquad->length; 
                fin_nquad++, fin_y++) 
              { 
               for(ifin_nquad  = fin_nquad; 
                   ifin_nquad != NULL; 
                   ifin_nquad  = ifin_nquad->next)
                 {
                  fin_x = ifin_nquad->x;
                  coa_x = fin_x/2;
                  
                  /* find correct coarse nquad */
                  for(tcoa_nquad = coa_nquad; coa_x > tcoa_nquad->x + tcoa_nquad->length; tcoa_nquad = tcoa_nquad->next)
                     ;
                  
                  /* jump to correct coa_node */
                  coa_node = tcoa_nquad->loc + (coa_x - tcoa_nquad->x);
                  
                  for(fin_node = ifin_nquad->loc; 
                      fin_node < ifin_nquad->loc + ifin_nquad->length; 
                      fin_node++, fin_x++)
                    {
                     
                     /* find all 26 neighbouring fine nodes */
                     tsc_finnodes[1][1][1] = fin_node;
                     get_TSCnodes(fin_grid, fin_pquad, ifin_cquad, ifin_nquad, tsc_finnodes, &fin_z, &fin_y, &fin_x);
                     
                     /* fin_node == boundary node? */
                     if(test_tsc(tsc_finnodes) == FALSE)
                       {
                        /* distance of fin_node to mother coa_node (in coa_grid units!) */
                        fc_sep[X] = ((double)fin_x - (double)(2*coa_x) - 0.5) / 2.;
                        fc_sep[Y] = ((double)fin_y - (double)(2*coa_y) - 0.5) / 2.;
                        fc_sep[Z] = ((double)fin_z - (double)(2*coa_z) - 0.5) / 2.;
                        
                        /* find all 26 neighbouring coarse nodes */
                        tsc_coanodes[1][1][1] = coa_node;
                        get_TSCnodes(coa_grid, tcoa_pquad, tcoa_cquad, tcoa_nquad, tsc_coanodes, &coa_z, &coa_y, &coa_x);
                        
                        /* x-direction */
                        /* calculate slopes in all three dimensions (note: the spacing is 1 in coa_grid units!) */
                        for(k=0; k<NDIM; k++)
                           for(j=0; j<NDIM; j++)
                              for(i=0; i<NDIM; i++)
                                 func[k][j][i] = tsc_coanodes[k][j][i]->force.forces[X];
                        get_c2fslope(func, slope);
                        
                        /* interpolate to fin_node */
                        fin_node->force.forces[X] = coa_node->force.forces[X] + ( slope[X]*fc_sep[X]
                                                                                  +   slope[Y]*fc_sep[Y]
                                                                                  +   slope[Z]*fc_sep[Z] );
                        /* y-direction */
                        /* calculate slopes in all three dimensions (note: the spacing is 1 in coa_grid units!) */
                        for(k=0; k<NDIM; k++)
                           for(j=0; j<NDIM; j++)
                              for(i=0; i<NDIM; i++)
                                 func[k][j][i] = tsc_coanodes[k][j][i]->force.forces[Y];
                        get_c2fslope(func, slope);
                        
                        /* interpolate to fin_node */
                        fin_node->force.forces[Y] = coa_node->force.forces[Y] + ( slope[X]*fc_sep[X]
                                                                                  +   slope[Y]*fc_sep[Y]
                                                                                  +   slope[Z]*fc_sep[Z] );
                        /* z-direction */
                        /* calculate slopes in all three dimensions (note: the spacing is 1 in coa_grid units!) */
                        for(k=0; k<NDIM; k++)
                           for(j=0; j<NDIM; j++)
                              for(i=0; i<NDIM; i++)
                                 func[k][j][i] = tsc_coanodes[k][j][i]->force.forces[Z];
                        get_c2fslope(func, slope);
                        
                        /* interpolate to fin_node */
                        fin_node->force.forces[Z] = coa_node->force.forces[Z] + ( slope[X]*fc_sep[X]
                                                                                  +   slope[Y]*fc_sep[Y]
                                                                                  +   slope[Z]*fc_sep[Z] );
                       }
                     
                     /* move to next coa_node */
                     if(is_even(fin_x) == FALSE)
                       {
                        coa_node++;
                        coa_x++;
                       }
                     
                    }
                 }
               
               /* move to next coa_nquad */
               if(is_even(fin_y) == FALSE)
                 {
                  coa_nquad++;
                  coa_y++;
                 }
              }
           }
         
         /* move to next coa_cquad */
         if(is_even(fin_z) == FALSE)
           {
            coa_cquad++;
            coa_z++;
           }
        }
     }  
   
}



/*==============================================================================
* calc_forces: 
*
* use the cur_node->pot values to calculate the forces by mid-centre differences
*==============================================================================*/
void calc_forces(gridls *cur_grid)
{
   pqptr   cur_pquad;
   cqptr   cur_cquad, icur_cquad;
   nqptr   cur_nquad, icur_nquad;
   nptr    cur_node, tsc_nodes[3][3][3];
   long    x, y, z, ipquad;
   double  force[NDIM], dx;
   
   /* the denominator of the gradient operator */
   dx = 2.*cur_grid->spacing;
   
   /* loop over all nodes */
#ifdef WITH_OPENMP
#pragma omp parallel firstprivate(dx) private(ipquad, cur_pquad, cur_cquad, icur_cquad, cur_nquad, icur_nquad, cur_node, tsc_nodes, x, y, z, force) shared(cur_grid)
#pragma omp for schedule(static)
   for(ipquad=0; ipquad<cur_grid->no_pquad; ipquad++) 
     {
      cur_pquad = cur_grid->pquad_array[ipquad];
#else
   for(cur_pquad=cur_grid->pquad; cur_pquad!=NULL; cur_pquad=cur_pquad->next)
     {
#endif
      for(cur_cquad = cur_pquad->loc, z = cur_pquad->z; 
          cur_cquad < cur_pquad->loc + cur_pquad->length;
          cur_cquad++, z++)
        {
         for(icur_cquad=cur_cquad; icur_cquad!=NULL; icur_cquad=icur_cquad->next)
           {
            for(cur_nquad = icur_cquad->loc, y = icur_cquad->y;
                cur_nquad < icur_cquad->loc + icur_cquad->length; 
                cur_nquad++, y++)
              {
               for(icur_nquad=cur_nquad; icur_nquad!=NULL; icur_nquad=icur_nquad->next)
                 {
                  for(cur_node = icur_nquad->loc, x = icur_nquad->x;
                      cur_node < icur_nquad->loc + icur_nquad->length; 
                      cur_node++, x++)
                    {
                     /* get neighboring nodes to check ...*/
                     tsc_nodes[1][1][1] = cur_node;
                     get_TSCnodes(cur_grid, cur_pquad, icur_cquad, icur_nquad, tsc_nodes, &z, &y, &x);
                     
                     /* ..if node fully lies inside refinement */
                     if(test_tsc(tsc_nodes) == TRUE)
                       {
                        force[X] = ((double)tsc_nodes[1][1][0]->pot - (double)tsc_nodes[1][1][2]->pot)/dx;
                        force[Y] = ((double)tsc_nodes[1][0][1]->pot - (double)tsc_nodes[1][2][1]->pot)/dx;
                        force[Z] = ((double)tsc_nodes[0][1][1]->pot - (double)tsc_nodes[2][1][1]->pot)/dx;
                        
                        cur_node->force.forces[X] = force[X];
                        cur_node->force.forces[Y] = force[Y];
                        cur_node->force.forces[Z] = force[Z];
                       }
                    }
                 }
              }
           }
        }
     }

#ifdef ADAPTIVE
   /* in case we are not on the domain grid obtain forces on boundary from next coarser level */
   if(cur_grid != global.dom_grid)
     {
      calc_forces(cur_grid-1);
      forces_on_boundary(cur_grid-1, cur_grid);
     }
#endif
}
