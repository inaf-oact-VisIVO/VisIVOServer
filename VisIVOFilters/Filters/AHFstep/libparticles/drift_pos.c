#include <math.h>
#include <stdlib.h>
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

/*==============================================================================
 * drift all particles linked to cur_grid by looping over all nodes
 *==============================================================================*/
void drift_pos(gridls *cur_grid, double DRIFT)
{
   long    ipquad;
   pqptr   cur_pquad;
   cqptr   cur_cquad, icur_cquad;
   nqptr   cur_nquad, icur_nquad;
   nptr    cur_node;
   partptr cur_part;
   int     i;
   double  distance;
   double  dr2, cur_part_dr2, cellfrac_max2;
  
   cellfrac_max2 = pow2(CELLFRAC_MAX);
   
   /* loop over all nodes */
#ifdef WITH_OPENMP
#pragma omp parallel firstprivate(cellfrac_max2) private(ipquad, cur_pquad, cur_cquad, icur_cquad, cur_nquad, icur_nquad, cur_node, cur_part, i, distance, dr2, cur_part_dr2) shared(cur_grid, io, global)
#pragma omp for schedule(static)
   for(ipquad=0; ipquad<cur_grid->no_pquad; ipquad++)     
     {
      cur_pquad = cur_grid->pquad_array[ipquad];
#else
   for(cur_pquad=cur_grid->pquad; cur_pquad!=NULL; cur_pquad=cur_pquad->next)
     {
#endif
      for(cur_cquad=cur_pquad->loc; cur_cquad<cur_pquad->loc+cur_pquad->length;	cur_cquad++) {
         
         for(icur_cquad=cur_cquad; icur_cquad!=NULL; icur_cquad=icur_cquad->next) {
            for(cur_nquad=icur_cquad->loc;cur_nquad<icur_cquad->loc+icur_cquad->length; cur_nquad++) {
               
               for(icur_nquad=cur_nquad; icur_nquad!=NULL; icur_nquad=icur_nquad->next) {
                  for(cur_node = icur_nquad->loc; cur_node < icur_nquad->loc + icur_nquad->length; cur_node++)
                    {
                     /* zero new_ll pointer (this erases force.temp[0] and force.temp[1]) */
                     cur_node->force.new_ll = NULL;
                     
                     /* loop over ll */
                     for(cur_part=cur_node->ll; cur_part!=NULL; cur_part=cur_part->ll)
                       {
                           cur_part_dr2 = 0.0;
                           
                           /* update pos[] of cur_part */
                            for(i = X; i <= Z; i++)
                            {
                              distance          = (double)cur_part->mom[i] * DRIFT; 
                              cur_part->pos[i]  = f1mod((double)cur_part->pos[i]+distance+1.0,1.0);
                              
                              cur_part_dr2     += pow2(distance);
                            }
                           
                           dr2 = cur_part_dr2/cur_grid->spacing2;
                           
                           /* keep track of maximum distance (in grid units) travelled by a particle */
                           if(dr2 > global.max_dr2)
                              global.max_dr2 = dr2;
                           

                           /* did particle travel farther than certain fraction of cell? */
                           if(dr2 > cellfrac_max2)
                             {
                              /* caught speeding particle */
                              global.speeding  = TRUE;
                              global.no_speeders++;
#ifndef ATS
#ifdef MULTIMASS
                              /* in case we do not run with ATS notify the user about speeder... */
                              fprintf(io.logfile,
                                      "drift_pos: %f > %f (grid: %ld, weight = %g)\n", 
                                      sqrt(cur_part_dr2)/cur_grid->spacing, CELLFRAC_MAX, 
                                      cur_grid->l1dim, cur_part->weight);
#else
                              fprintf(io.logfile,"drift_pos: %f > %f (grid: %ld)\n", 
                                      sqrt(cur_part_dr2)/cur_grid->spacing, CELLFRAC_MAX, 
                                      cur_grid->l1dim);
#endif
                              fflush(io.logfile);
                              
#ifdef TERMINATE     
                              /* ... and terminate AMIGA */
                              global.terminate = TRUE;
#endif
#endif /* !ATS */
                             }
                           
                       } /* cur_part loop */

                    } /* cur_node loop */
                  
               } } } } } /* all those quad loops */
   
   /* increase multi-step counter used as normalisation for global.no_speeders */
   if(global.speeding)
      global.no_ssteps++;

}







