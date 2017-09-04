#include <stddef.h>
#include <math.h>
#ifdef WITH_OPENMP
#include <omp.h>
#endif

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"


/* ...and now for the actual includes */
#include "../libutility/utility.h"
#include "../libgravity/gravity.h"
#include "../libgrids/grids.h"

/*==============================================================================
 * kick_mom: update particle momentum
 *==============================================================================*/
void kick_mom(gridls *cur_grid, double KICK)
{
   pqptr   cur_pquad;
   cqptr   cur_cquad, icur_cquad;
   nqptr   cur_nquad, icur_nquad;
   nptr    cur_node;
   partptr cur_part;
   long    x, y, z;
   double  cur_shift;
   long    ipquad;
   
   nptr    tsc_nodes[3][3][3];  /* nodes to assign to                 */
   dvect   xyz_coords;          /* actual double coords of node        */
   dvect   pn_sep;              /* particle-node separation / spacing */
   int     idim;                /* coord changing index               */
   int     i,j,k;               /* indices for 3D arrays              */
   int     n;                   /* index for components               */
   dvect   weights[3];          /* weights in each dimension          */
   double  pnarg_a;             /* pn sep temp arg                    */
   double  pnarg_b;             /* pn sep temp arg                    */
   double  tpnarg;              /* temp to calc pn args               */
   dvect   temp_coords;         /* double vector - calc un_mod coords  */
   double  force_comp[NDIM];    /* array for force components         */
   
   double  momchange;
   double  cur_part_dp2;
   
   int     icount;
   double  mean_gM;
   double  mean_gN; 
   double  g_0;
   double  g_N, g_Nr;
   double  g_M, g_Mr;
   double  a_current;
   
#ifdef FORCETEST
   double  ppot, pdens;
#endif

  icount    = 0;
   mean_gM   = 0.0;
   mean_gN   = 0.0;
   
   /* the MOND correction requires knowledge about the current expansion factor a */
   a_current = calc_super_a(cur_grid->timecounter);
   
   /* shift of cell centre as compared to edge of box [grid units] */
   cur_shift = 0.5/(double)cur_grid->l1dim;
   
   
   /* loop over all nodes */
#ifdef WITH_OPENMP
#pragma omp parallel reduction(+:icount) reduction(+:mean_gM) reduction(+:mean_gN) firstprivate(cur_shift, a_current) private(ipquad, cur_pquad, cur_cquad, icur_cquad, cur_nquad, icur_nquad, cur_node, cur_part, x, y, z, tsc_nodes, xyz_coords, pn_sep, idim, i, j, k, n, weights, pnarg_a, pnarg_b, tpnarg, temp_coords, force_comp, momchange, cur_part_dp2, g_0, g_N, g_Nr, g_M, g_Mr) shared(cur_grid, global)
#pragma omp for schedule(static)
   for(ipquad=0; ipquad<cur_grid->no_pquad; ipquad++)     
     {
      cur_pquad=cur_grid->pquad_array[ipquad];
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
                     /* are there any particles attached to this node? */
                     if(cur_node->ll != NULL)
                       {
                        /* current node in centre */
                        tsc_nodes[1][1][1] = cur_node;
                        
                        /* calculate realspace coords of cur_node */
                        temp_coords[X] = (((double)x) / (double)cur_grid->l1dim) + cur_shift;
                        temp_coords[Y] = (((double)y) / (double)cur_grid->l1dim) + cur_shift;
                        temp_coords[Z] = (((double)z) / (double)cur_grid->l1dim) + cur_shift;
                        
                        xyz_coords[X]  = f1mod(temp_coords[X]+1.0, 1.0);
                        xyz_coords[Y]  = f1mod(temp_coords[Y]+1.0, 1.0);
                        xyz_coords[Z]  = f1mod(temp_coords[Z]+1.0, 1.0);
                        
                        /* get pointers to tsc nodes */
                        get_TSCnodes(cur_grid, cur_pquad, icur_cquad, icur_nquad, tsc_nodes, &z, &y, &x);
                        
                        /* loop over all particles at current node */
                        if(test_tsc(tsc_nodes) == TRUE)
                          {
                           for(cur_part = cur_node->ll; cur_part != NULL; cur_part = cur_part->ll)
                             {
                              /* calc fraction to be assigned to each tsc node */
                              for(idim = 0; idim < NDIM; idim++)
                                {
                                 pn_sep[idim] = ((double)cur_part->pos[idim] - xyz_coords[idim]) * (double)cur_grid->l1dim;
                                 
                                 /* check for periodic boundaries */
                                 if(fabs(pn_sep[idim]) > 0.5*(double)cur_grid->l1dim)
                                   {
                                    tpnarg          = (double)cur_part->pos[idim] + 0.5;
                                    pnarg_a         = f1mod(tpnarg+1.0, 1.0);
                                    tpnarg          = xyz_coords[idim] + 0.5;
                                    pnarg_b         = f1mod(tpnarg+1.0, 1.0);
                                    pn_sep[idim] = (pnarg_a - pnarg_b) * (double)cur_grid->l1dim;
                                   }
#ifdef TSC
                                 weights[0][idim] = pow2((0.5 - pn_sep[idim]))/2;
                                 weights[1][idim] = 0.75 - pow2(pn_sep[idim]);
                                 weights[2][idim] = pow2((0.5 + pn_sep[idim]))/2;
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
                                    weights[0][idim] =      -pn_sep[idim];
                                    weights[1][idim] = 1.0 + pn_sep[idim];
                                    weights[2][idim] = 0.0;
                                   }
#endif
#ifdef NGP
                                 weights[0][idim] = 0.0;
                                 weights[1][idim] = 1.0;
                                 weights[2][idim] = 0.0;
#endif
                                 /* use opportunity to reset force vector */
                                 force_comp[idim] = 0.0;
                                }	
                              
#ifdef FORCETEST
                               /* interpolate potential and density off the grid to particle position */
                               ppot  = 0.0;
                               pdens = 0.0;
#endif

                               /* calculate force at particle position */
                              for(k = 0; k < 3; k++)
                                 for(j = 0; j < 3; j++)
                                    for(i = 0; i < 3; i++)
                                      {
                                       for(n = X; n <= Z; n++)
                                         {
                                          force_comp[n] += weights[k][Z]*weights[j][Y]*weights[i][X]
                                          * (double)tsc_nodes[k][j][i]->force.forces[n];
                                         }
#ifdef FORCETEST
                                        pdens += weights[k][Z]*weights[j][Y]*weights[i][X]
                                        * (double)tsc_nodes[k][j][i]->dens;
                                        ppot  += weights[k][Z]*weights[j][Y]*weights[i][X]
                                        * (double)tsc_nodes[k][j][i]->pot;
#endif
                                      }
                                       
#ifdef MOND
                              /* get peculiar acceleration */
                              g_N = sqrt(pow2(force_comp[0]) +
                                         pow2(force_comp[1]) + 
                                         pow2(force_comp[2]));
                              g_Nr = g_N / pow2(a_current);
                              
                              /* get (peculiar )MONDian acceleration */
                              g_0  = simu.g0;
                              g_Mr = get_MONDg(g_0, g_Nr);
                              
                              /* transfer to internal acceleration */
                              g_M  = g_Mr * pow2(a_current);
                              
                              /* rescale acceleration vector */
                              force_comp[0] *= g_M/g_N;
                              force_comp[1] *= g_M/g_N;
                              force_comp[2] *= g_M/g_N;
                              
                              /* some statistics */
                              global.no_ALLevents += 1;
                              if(g_Mr/g_Nr > 1.1)           /* 10% bigger g => MOND event */
                                 global.no_MONDevents += 1;
                              
                              if(g_Mr > global.max_gM)
                                 global.max_gM = g_Mr;
                              
                              mean_gM += g_Mr;
                              /* log mean (Newtonian) peculiar acceleration */
                              icount++;
                              mean_gN += g_Nr;
                              if(g_Nr > global.max_gN)
                                 global.max_gN = g_Nr;
                              
#endif   /* MOND */
                              cur_part_dp2 = (double)0.0;
                              
                              /* update particle's momentum */
                              for(n = X; n <= Z; n++)
                                {
                                 /* the actual change in momentum*/
                                 momchange        = (double)force_comp[n] * KICK;
                                 
                                 /* update particle's momentum */
                                 cur_part->mom[n] +=  momchange;

                                 /* statistics for the logfile */
                                 cur_part_dp2    += pow2(momchange);                                  
                                }
                               
#ifdef FORCETEST
                               /* store force vector, potential and density at particle position */
                               cur_part->forces[X] = force_comp[X];
                               cur_part->forces[Y] = force_comp[Y];
                               cur_part->forces[Z] = force_comp[Z];
                               cur_part->pot       = ppot;
                               cur_part->dens      = pdens;
#endif
                              
                              /* statistics for the logfile */
                              if(cur_part_dp2 > global.max_dp2)
                                 global.max_dp2 = cur_part_dp2;
                             }
                          }
                       }
                    }
                 }
              }
           }
        }
     }
   
#ifdef MOND
   global.steps   += 1;
   global.mean_gN += mean_gN/(double)icount;
   global.mean_gM += mean_gM/(double)icount;
#endif
}