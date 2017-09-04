#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#ifdef WITH_OPENMP
#include <omp.h>
#endif

void dummy_store_refs()
{
}

#ifdef STORE_REFS

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"


/* ...and now for the actual includes */
#include "grids.h"
#include "../libutility/utility.h"

/*===============================================================================
*
* The visible functions is:
*
* copy_grid(grid1, grid2):       copies hydro-variables from grid1 -> grid2
*
*===============================================================================*/

void copy_grid(gridls *cur_grid)
{
   pqptr         cur_pquad1;
   cqptr         cur_cquad1, icur_cquad1;
   nqptr         cur_nquad1, icur_nquad1;
   nptr          cur_node1;
   
   pqptr         cur_pquad2;
   cqptr         cur_cquad2, icur_cquad2;
   nqptr         cur_nquad2, icur_nquad2;
   nptr          cur_node2;
   
   long          x, y, z, ipquad;
   int           ivar;
   
   /*----------------------------------------------
    * loop over all nodes accessible via old_pquad
    *----------------------------------------------*/
#ifdef WITH_OPENMP
#pragma omp parallel private(ipquad, cur_pquad1, cur_cquad1, icur_cquad1, cur_nquad1, icur_nquad1, cur_node1, cur_pquad2, cur_cquad2, icur_cquad2, cur_nquad2, icur_nquad2, cur_node2, x, y, z, ivar)  shared(cur_grid)
#pragma omp for schedule(static)
   for(ipquad=0; ipquad<cur_grid->no_pquad; ipquad++)
     {
      cur_pquad1=cur_grid->pquad_array[ipquad];
#else
   for(cur_pquad1=cur_grid->old_pquad; cur_pquad1!=NULL; cur_pquad1=cur_pquad1->next)
     {
#endif
      z = cur_pquad1->z;
      for(cur_cquad1 = cur_pquad1->loc;
          cur_cquad1 < cur_pquad1->loc + cur_pquad1->length; 
          cur_cquad1++, z++)  
        {  
         for(icur_cquad1  = cur_cquad1; 
             icur_cquad1 != NULL; 
             icur_cquad1  = icur_cquad1->next)
           {
            y = icur_cquad1->y;
            for(cur_nquad1 = icur_cquad1->loc;  
                cur_nquad1 < icur_cquad1->loc + icur_cquad1->length; 
                cur_nquad1++, y++) 
              { 
               for(icur_nquad1  = cur_nquad1; 
                   icur_nquad1 != NULL; 
                   icur_nquad1  = icur_nquad1->next)
                 {
                  x = icur_nquad1->x;
                  for(cur_node1 = icur_nquad1->loc; 
                      cur_node1 < icur_nquad1->loc + icur_nquad1->length; 
                      cur_node1++, x++)
                    {
                     /*-------------------------------------------------
                     * locate node (x,y,z) in newly created refinement
                     *-------------------------------------------------*/
                     for(cur_pquad2 = cur_grid->pquad;
                         (z > cur_pquad2->z + cur_pquad2->length) && 
                         (cur_pquad2->next != NULL); 
                         cur_pquad2 = cur_pquad2->next)
                        ;
                     
                     if((z >= cur_pquad2->z) && 
                        (z < cur_pquad2->z + cur_pquad2->length))
                       {
                        for(cur_cquad2 = cur_pquad2->loc + (z-cur_pquad2->z);
                            (y > cur_cquad2->y + cur_cquad2->length) && 
                            (cur_cquad2->next != NULL); 
                            cur_cquad2 = cur_cquad2->next)
                           ;
                        
                        if((y >= cur_cquad2->y) && 
                           (y < cur_cquad2->y + cur_cquad2->length))
                          {
                           for(cur_nquad2=cur_cquad2->loc+(y-cur_cquad2->y);
                               (x > cur_nquad2->x + cur_nquad2->length) &&
                               (cur_nquad2->next != NULL);
                               cur_nquad2 = cur_nquad2->next)
                              ;
                           
                           if((x >= cur_nquad2->x) && 
                              (x < cur_nquad2->x + cur_nquad2->length))
                             {
                              /* eventually reached correct node */
                              cur_node2 = cur_nquad2->loc+(x-cur_nquad2->x);
                              
                              /* copy hydro-variables */
                              for(ivar=0; ivar<NHYDRO; ivar++)
                                 cur_node2->u[ivar] = cur_node1->u[ivar];
                             }
                           else
                             {
                              /* node not present */
                             }
                          }
                        else
                          {
                           /* node not present */
                          }
                       }
                     else
                       {
                        /* node not present */
                       }
                    }
                 }
              }
           }
        }
     }
   
}

/*==============================================================================
* store_grid: save pquad pointer
*==============================================================================*/
void store_pquad(gridls *cur_grid, int *no_grids)
{
   
   /* save pointer to whole grid structure and reset actual pointer */
   cur_grid->old_pquad = cur_grid->pquad;
   cur_grid->pquad     = NULL;
   
   /* reset all sorts of things to zero */
   cur_grid->multistep = 0;
   
   /* number of grids needs to be reduced anyways... */
   (*no_grids) -= 1;
}

#endif
