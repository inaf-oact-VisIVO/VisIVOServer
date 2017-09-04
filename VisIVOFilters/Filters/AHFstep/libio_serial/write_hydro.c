#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void dummy_write_hydro()
{
}

#ifdef HYDRO

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "../libutility/utility.h"
#include "../libgrids/grids.h"
#include "../libmhd/mhd.h"

/*=============================================================================
* simply write the position of all nodes
*=============================================================================*/
void write_hydro(gridls *cur_grid)
{
   char          filename[MAXSTRING], file_no[MAXSTRING];
   FILE          *hydrofile;
   partptr       cur_part;
   pqptr         cur_pquad;
   cqptr         cur_cquad, icur_cquad;
   nqptr         cur_nquad, icur_nquad;
   nptr          cur_node;
   long          ipart, ipart_node, x, y, z;
   int           i, ivar, idim, slen;
   
   double        cur_shift;
   
   double        dist, T, u[NHYDRO], mean_dens, mean_dens_tmp;
   long          nnodes;
   
   double        output_shift;
   
   
#ifdef MHD
   nptr          tsc_nodes[3][3][3];
   nptr          MHDnodes[5][5][5];
   double        B[NDIM], divB;
#endif // MHD
   
   cur_shift     = 0.5/(double)cur_grid->l1dim;
   mean_dens     = 0.0;
   mean_dens_tmp = 0.0;
   nnodes        = 0;
   
   /* prepare filename */
   sprintf(file_no,"%05d",global.no_timestep);
   strcpy(filename,io.outfile_prefix);
   strcat(filename,"hydro-");
   strcat(filename, file_no);
   
   /* actually open file */
   if((hydrofile = fopen(filename,"w")) == NULL)
     {
      fprintf(stderr,"could not open %s\n", filename);
      exit(1);
     }
#if (HYDRO_TEST==5 || HYDRO_TEST==6 || HYDRO_TEST==7 || HYDRO_TEST==9 || HYDRO_TEST==10 || HYDRO_TEST==11 || HYDRO_TEST==12)
   
#ifdef MHD
   fprintf(hydrofile,"#           x(1)         Udens(2)     UmomdensX(3)        UEdens(4)        Untrpy(5)        Uedens(6)        Upress(7)         force(8)             T(9)            M(10)   rho+DMdens(11)           Bx(12)           By(13)           Bz(14)         divB(15)    UmomdensY(16)    UmomdensZ(17)\n");
#else // MHD   
   fprintf(hydrofile,"#           x(1)         Udens(2)     UmomdensX(3)        UEdens(4)        Untrpy(5)        Uedens(6)        Upress(7)         force(8)             T(9)            M(10)   rho+DMdens(11)\n");
#endif // MHD
   
#else // (HYDRO_TEST==5 || HYDRO_TEST==6 || HYDRO_TEST==7 || HYDRO_TEST==9 || HYDRO_TEST==10 || HYDRO_TEST==11 || HYDRO_TEST==12)
   fprintf(hydrofile,"#x(1)  y(2)  z(3)  dens(4)  Vx(5)   Vy(6)   Vz(7)  Edens(8)   S(9)  p_e(10)  p_S(11)  edens_Edens(12)  edens_S(13)    T(14)    dist(15)  phi(16)  Fx(17) Fy(18) Fz(19)\n");
#endif // (HYDRO_TEST==5 || HYDRO_TEST==6 || HYDRO_TEST==7 || HYDRO_TEST==9 || HYDRO_TEST==10 || HYDRO_TEST==11 || HYDRO_TEST==12)
   
   /* loop over all nodes */
   for(cur_pquad=cur_grid->pquad; cur_pquad!=NULL; cur_pquad=cur_pquad->next)
     {
#if (HYDRO_TEST==0 || HYDRO_TEST==1 || HYDRO_TEST==4 || HYDRO_TEST==8 || NO_DM)
      z         = cur_grid->l1dim/2;
      cur_cquad = cur_pquad->loc + z;
#else
      for(cur_cquad = cur_pquad->loc, z = cur_pquad->z; 
          cur_cquad < cur_pquad->loc + cur_pquad->length;
          cur_cquad++, z++)
#endif
        {
         for(icur_cquad=cur_cquad; icur_cquad!=NULL; icur_cquad=icur_cquad->next)
           {
#if (HYDRO_TEST==0 || HYDRO_TEST==2 || HYDRO_TEST==4 || HYDRO_TEST==8 || NO_DM)
            y         = cur_grid->l1dim/2;
            cur_nquad = icur_cquad->loc + y;
#else
            for(cur_nquad = icur_cquad->loc, y = icur_cquad->y;
                cur_nquad < icur_cquad->loc + icur_cquad->length; 
                cur_nquad++, y++)
#endif
              {
               for(icur_nquad=cur_nquad; icur_nquad!=NULL; icur_nquad=icur_nquad->next)
                 {
#if (HYDRO_TEST==1 || HYDRO_TEST==2)
                  x        = cur_grid->l1dim/2;
                  cur_node = icur_nquad->loc + x;
#else
                  for(cur_node = icur_nquad->loc, x = icur_nquad->x;
                      cur_node < icur_nquad->loc + icur_nquad->length; 
                      cur_node++, x++)
#endif
                    {
                     
                     dist = sqrt(pow2((double)x-(double)(cur_grid->l1dim/2-0.5))+
                                 pow2((double)y-(double)(cur_grid->l1dim/2-0.5))+
                                 pow2((double)z-(double)(cur_grid->l1dim/2-0.5)))/(double)cur_grid->l1dim;
#if (HYDRO_TEST<4)
                     if(x > 1 && x < cur_grid->l1dim-2 &&
                        y > 1 && y < cur_grid->l1dim-2 &&
                        z > 1 && z < cur_grid->l1dim-2  )
#endif
                       {
                        /* double check mean density */
                        mean_dens     += cur_node->u[Udens];
                        mean_dens_tmp += cur_node->u_tmp[Udens];
                        nnodes++;
                        
                        /*==============================================================
                         *            WRITE ALL HYDRO-VARIABLES TO FILE 
                         *==============================================================*/
#ifdef BLAST_WAVE
                        if(fabs(cur_node->u[Udens]-1.) > ZERO)
#else
                        if(cur_node->u[Udens] > -10.)
#endif
                          {
                           for(ivar=0; ivar<NHYDRO; ivar++)
                              u[ivar] = (double)cur_node->u[ivar];
#ifdef MHD
                           for(idim=0; idim<NDIM; idim++)
                              B[idim] = (double)cur_node->B[idim];
#endif // MHD

                           T     = calc_T(u[Uedens]/u[Udens]) * pow2(simu.boxsize/simu.t_unit/global.a);
                           
#if (HYDRO_TEST==5 || HYDRO_TEST==6 || HYDRO_TEST==7 || HYDRO_TEST==9 || HYDRO_TEST==10 || HYDRO_TEST==11 || HYDRO_TEST==12)
#ifdef MHD
                             // get neighbour nodes and calculate div B
                             tsc_nodes[1][1][1] = cur_node;
                             get_TSCnodes(cur_grid, cur_pquad, icur_cquad, icur_nquad, tsc_nodes, &z, &y, &x);
                             
                             get_MHDnodes(cur_grid, cur_pquad, z, y, x, MHDnodes);
                             if(test_mhd(MHDnodes)==TRUE)
                             //if(test_tsc(tsc_nodes)==TRUE)
                                
                             divB =  ( tsc_nodes[1][1][2]->B[X] - cur_node->B[X] )
                                   + ( tsc_nodes[1][2][1]->B[Y] - cur_node->B[Y] )
                                   + ( tsc_nodes[2][1][1]->B[Z] - cur_node->B[Z] );
                             else
                                divB = 0.;
                             
#if (defined(SHOCK_MHD) || HYDRO_TEST==12)
                             output_shift = 0.;
#else
                             output_shift = 0.25;
#endif
                             
                             
                           // exactly as normal hydro file output - plus the values of Bx, By, Bz, div B, UmomdensY, UmondensZ
                           /* this line needs to be adjusted when using -DSWAP_XY or -DSWAP_XZ */
#if (HYDRO_TEST==12)
                           if(z==cur_grid->l1dim/2) // write data for the whole xy-plane
                            fprintf(hydrofile,"%16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g\n",
                                 fmod((double)x/(double)cur_grid->l1dim+cur_shift+output_shift,1.),
                                 fmod((double)y/(double)cur_grid->l1dim+cur_shift+output_shift,1.),
#else // HYDRO_TEST==12
                           if(y==cur_grid->l1dim/2 && z==cur_grid->l1dim/2) // write data only in x-direction
                            fprintf(hydrofile,"%16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g\n",
                                        fmod((double)x/(double)cur_grid->l1dim+cur_shift+output_shift,1.),
#endif // HYDRO_TEST==12
                                        /* numerical solutions: */
                                        u[Udens],
                                        u[UmomdensX],
                                        u[UEdens],
                                        u[Untrpy],
                                        u[Uedens],
                                        (simu.gamma-1.)*u[Uedens],
                                        sqrt(pow2(cur_node->force.forces[X])+pow2(cur_node->force.forces[Y])+pow2(cur_node->force.forces[Z])),
                                        /* derived quantities: (in physical units) */
                                        pow2(simu.boxsize/simu.t_unit)*calc_T(u[Uedens]/u[Udens]),
                                        0.0,
                                        u[Udens]+cur_node->dens,
                                        B[X],
                                        B[Y],
                                        B[Z],
                                        divB,                                    
                                        u[UmomdensY],
                                        u[UmomdensZ]); 
                             
#else // MHD //
                           /* this line needs to be adjusted when using -DSWAP_XY or -DSWAP_XZ */
                           if(z==cur_grid->l1dim/2 && y==cur_grid->l1dim/2)
                            fprintf(hydrofile,"%16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g %16.8g\n",
                                    fmod((double)x/(double)cur_grid->l1dim+cur_shift+0.25,1.),
                                    /* numerical solutions: */
                                    u[Udens],
                                    u[UmomdensX],
                                    u[UEdens],
                                    u[Untrpy],
                                    u[Uedens],
                                    (simu.gamma-1.)*u[Uedens],
                                    sqrt(pow2(cur_node->force.forces[X])+pow2(cur_node->force.forces[Y])+pow2(cur_node->force.forces[Z])),
                                    /* derived quantities: (in physical units) */
                                    pow2(simu.boxsize/simu.t_unit)*calc_T(u[Uedens]/u[Udens]),
                                    0.0,
                                    u[Udens]+cur_node->dens); 
#endif // MHD //
#else /* HYDRO_TEST==5 || HYDRO_TEST==6 || HYDRO_TEST==7 || HYDRO_TEST==9  || HYDRO_TEST==10 || HYDRO_TEST==11 || HYDRO_TEST==12 */
                             
                           
#if (NO_DM || DM_GAS)
                           fprintf(hydrofile,"%g  %g %g %g %g %g\n",
                                   (double)(x+0.5)/(double)cur_grid->l1dim * simu.boxsize,
                                   (double)(y+0.5)/(double)cur_grid->l1dim * simu.boxsize,
                                   (double)(z+0.5)/(double)cur_grid->l1dim * simu.boxsize,
                                   u[UmomdensX]/u[Udens] * simu.boxsize/simu.t_unit/calc_super_a(cur_grid->timecounter),
                                   u[UmomdensY]/u[Udens] * simu.boxsize/simu.t_unit/calc_super_a(cur_grid->timecounter),
                                   u[UmomdensZ]/u[Udens] * simu.boxsize/simu.t_unit/calc_super_a(cur_grid->timecounter));
#else /* NO_DM || DM_GAS */
                           fprintf(hydrofile,"%lf %lf %lf    %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf %12.8lf     %g %g %g %g %g   %g  %g   %g %g %g\n",
                                   (double)x/(double)cur_grid->l1dim+cur_shift,
                                   (double)y/(double)cur_grid->l1dim+cur_shift,
                                   (double)z/(double)cur_grid->l1dim+cur_shift,
                                   
                                   u[Udens],
                                   u[UmomdensX],
                                   u[UmomdensY],
                                   u[UmomdensZ],
                                   u[UEdens],
                                   u[Untrpy],
                                   
                                   (simu.gamma-1.) * u[Uedens],
                                   u[Untrpy] * pow(u[Udens],simu.gamma-1),
                                   u[Uedens],
                                   u[Untrpy] * pow(u[Udens],simu.gamma-1) / (simu.gamma-1.),
                                   
                                   T,
                                   
                                   dist,
                                   
                                   cur_node->pot,
                                   cur_node->force.forces[X],
                                   cur_node->force.forces[Y],
                                   cur_node->force.forces[Z]);
#endif /* NO_DM */
#endif
                          }
                     
                       }
                    }
                 }
              }
           }
        }
     }
      
   fprintf(stderr,"write_hydro:   mean_dens     = %g\n",mean_dens/(double)nnodes);
   fprintf(stderr,"               mean_dens_tmp = %g\n",mean_dens_tmp/(double)nnodes);
   fclose(hydrofile);
}

#endif
