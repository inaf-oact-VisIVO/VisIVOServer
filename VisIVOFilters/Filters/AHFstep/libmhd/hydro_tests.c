#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "mhd.h"
#include "../libutility/utility.h"
#include "../libgrids/grids.h"
#include "../libgravity/gravity.h"
#include "../libparticles/particles.h"
#include "../libio_serial/io_serial.h"

/* make Zeldovich parameters globally available */
static int    n1dim;
static double z_caust, z_caust2, kwave, kwave2;
static double x_fac, v_fac, e_fac, T_fac;

//#define WRITE_ZELDOVICH
//#define WRITE_ZELDOVICH_TERMINATE

//#define WRITE_GAS_SPHERE

void dummy_hydro_tests()
{
}

#ifdef HYDRO
void write_ZW(double timecounter, int no_timestep)
{
   double qx, rho, phi, Fx, Amp, Amp2, a, a_dot, x, vx, Temp, rho_init, a_init;
   int    i;
   char   fname[MAXSTRING];
   FILE   *ftmp;
   
   sprintf(fname,"ZeldovichWave-%d.DAT",no_timestep);
   ftmp = fopen(fname,"w");
   
   Amp    = 1.+z_caust;
   Amp2   = 1.+z_caust2;
   
   a_init = simu.a_initial;
   
   a      = calc_super_a(timecounter);
   a_dot  = calc_a_dot(a);
   a_dot  = H0*sqrt(1./a);
   
   v_fac  = pow2(a)/(H0*simu.boxsize); /* v_fac depends on a */
   
   fprintf(ftmp,"#x(1) rho(2) vx(3) phi(4) Fx(5)  T(6)\n");

   for(i=0; i<n1dim; i++)
     {
      /* double check Zeldovich wave */
      qx       = (((double)i+0.5)/(double)n1dim) * simu.boxsize;
      x        = qx +    a     * (Amp) *cos(kwave *qx)/kwave
                    +    a     * (Amp2)*cos(kwave2*qx-TWOPI/4.)/kwave2;
      vx       =         a_dot * (Amp) *cos(kwave *qx)/kwave
                    +    a_dot * (Amp2)*cos(kwave2*qx-TWOPI/4.)/kwave2;
      rho      = 1./(1. - a    * (Amp) *sin(kwave *qx) 
                        - a    * (Amp2)*sin(kwave2*qx-TWOPI/4.)) - 1.;
      
#ifdef GREG_PANCAKE
      /* Greg Bryan's formula */
      vx       =         a_dot * (Amp) *cos(kwave *qx)
                    +    a_dot * (Amp2)*cos(kwave2*qx-TWOPI/4.);
#endif
      
      /* phi and Fx only correct for single pancake */
      phi      = -a*Amp/pow2(kwave) * sin(kwave*qx)*(1.-a*Amp/2.*sin(kwave*qx));
      Fx       =  a*Amp/(kwave)*cos(kwave*qx);
      Temp     = simu.T_init * pow(((rho_init+1.)/(rho+1.) ), (simu.gamma-1.));
      
#ifdef NO_DM
      fprintf(ftmp,"%lf   %lf\n",
              x,vx);
#else
      fprintf(ftmp,"%lf   %lf   %lf   %lf   %lf    %lf\n",
              fmod(x_fac*x+0.25,1.),rho+1.,v_fac*vx,phi,Fx,Temp);
#endif
     }
   
   fclose(ftmp);
}

/*==============================================================================
 * zero u-values on cur_grid
 *==============================================================================*/
void zero_udens(gridls *cur_grid)
{
   pqptr   cur_pquad;
   cqptr   cur_cquad, icur_cquad;
   nqptr   cur_nquad, icur_nquad;
   nptr    cur_node;
   partptr cur_part;
   int     idim, i, j, k;
   long    x, y, z;
   
   /*----------------------------------------------------------------------
    * reset hydro-variables
    *----------------------------------------------------------------------*/
   for(cur_pquad=cur_grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next)
     {
      z = cur_pquad->z;
      for(cur_cquad = cur_pquad->loc;
          cur_cquad < cur_pquad->loc + cur_pquad->length; 
          cur_cquad++, z++)  
        {  
         for(icur_cquad  = cur_cquad; 
             icur_cquad != NULL; 
             icur_cquad  = icur_cquad->next)
           {
            y = icur_cquad->y;
            for(cur_nquad = icur_cquad->loc;  
                cur_nquad < icur_cquad->loc + icur_cquad->length; 
                cur_nquad++, y++) 
              { 
               for(icur_nquad  = cur_nquad; 
                   icur_nquad != NULL; 
                   icur_nquad  = icur_nquad->next)
                 {
                  x = icur_nquad->x;
                  for(cur_node = icur_nquad->loc; 
                      cur_node < icur_nquad->loc + icur_nquad->length; 
                      cur_node++, x++)
                    {
                     cur_node->u[Udens]     = 0.0;
                     cur_node->u[UmomdensX] = 0.0;
                     cur_node->u[UmomdensY] = 0.0;
                     cur_node->u[UmomdensZ] = 0.0;
                     cur_node->u[UEdens]    = 0.0;
                     cur_node->u[Untrpy]    = 0.0;
                     cur_node->u[Uedens]    = 0.0;
#ifdef MHD
                     cur_node->B[X] = 0.0;
                     cur_node->B[Y] = 0.0;
                     cur_node->B[Z] = 0.0;
#endif
                    }
                 }
              }
           }
        }
     }
}

/*==============================================================================
 * initialize u-values on cur_grid using particles attached to linked list
 *==============================================================================*/
void init_udens(gridls *cur_grid)
{
   pqptr   cur_pquad;
   cqptr   cur_cquad, icur_cquad;
   nqptr   cur_nquad, icur_nquad;
   nptr    cur_node, tsc_nodes[3][3][3];
   partptr cur_part;
   int     idim, i, j, k;
   long    x, y, z;
   double  cur_shift;
   double  pnarg_a;
   double  pnarg_b;
   double  tpnarg;
   double  xyz_coords[NDIM];
   double  pn_sep[NDIM];
   double  weights[3][NDIM];
   double  temp_coords[NDIM];   
   
   double  e_kin, e_therm, e_tot, Bdens;
   
   cur_shift = 0.5/(double)cur_grid->l1dim;
   
   /*----------------------------------------------------------------------
    * get mass and momentum density
    *----------------------------------------------------------------------*/
   for(cur_pquad=cur_grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next)
     {
      z = cur_pquad->z;
      for(cur_cquad = cur_pquad->loc;
          cur_cquad < cur_pquad->loc + cur_pquad->length; 
          cur_cquad++, z++)  
        {  
         for(icur_cquad  = cur_cquad; 
             icur_cquad != NULL; 
             icur_cquad  = icur_cquad->next)
           {
            y = icur_cquad->y;
            for(cur_nquad = icur_cquad->loc;  
                cur_nquad < icur_cquad->loc + icur_cquad->length; 
                cur_nquad++, y++) 
              { 
               for(icur_nquad  = cur_nquad; 
                   icur_nquad != NULL; 
                   icur_nquad  = icur_nquad->next)
                 {
                  x = icur_nquad->x;
                  for(cur_node = icur_nquad->loc; 
                      cur_node < icur_nquad->loc + icur_nquad->length; 
                      cur_node++, x++)
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
                     
                     /* loop over ll */
                     for(cur_part = cur_node->ll; cur_part != NULL; cur_part = cur_part->ll)
                       {
                        /* calc fraction to be assigned to each tsc node */
                        for(idim = 0; idim < NDIM; idim++)
                          {
                           pn_sep[idim] = ((double)cur_part->pos[idim] - xyz_coords[idim]) * (double)cur_grid->l1dim;
                           
                           if(fabs(pn_sep[idim]) > 0.5*(double)cur_grid->l1dim)
                             {
                              tpnarg       = (double)cur_part->pos[idim] + 0.5;
                              pnarg_a      = f1mod(tpnarg+1.0, 1.0);
                              tpnarg       = xyz_coords[idim] + 0.5;
                              pnarg_b      = f1mod(tpnarg+1.0, 1.0);
                              pn_sep[idim] = (pnarg_a - pnarg_b) * (double)cur_grid->l1dim;
                             }
                           
#ifdef TSC
                           weights[1][idim] = 0.75 - pow2(pn_sep[idim]);
                           weights[0][idim] = pow2((0.5 - pn_sep[idim]))/2;
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
                          } /* idim */
                        
                        /* assign mass and momentum */
                        for(k = 0; k < 3; k++)
                           for(j = 0; j < 3; j++)
                              for(i = 0; i < 3; i++)
                                 if(tsc_nodes[k][j][i] != NULL)
                                   {
                                    /* mass density u[Udens] */
                                    tsc_nodes[k][j][i]->u[Udens] +=  
                                       cur_grid->masstodens * weights[k][Z]*weights[j][Y]*weights[i][X];
                                    
                                    /* mass density u_tmp[Udens] */
                                    tsc_nodes[k][j][i]->u_tmp[Udens] +=  
                                       cur_grid->masstodens * weights[k][Z]*weights[j][Y]*weights[i][X];
                                    
                                    /* momentum density */
                                    tsc_nodes[k][j][i]->u[UmomdensX] +=  
                                       cur_grid->masstodens * weights[k][Z]*weights[j][Y]*weights[i][X] *
                                       cur_part->mom[X];
                                    tsc_nodes[k][j][i]->u[UmomdensY] +=  
                                       cur_grid->masstodens * weights[k][Z]*weights[j][Y]*weights[i][X] *
                                       cur_part->mom[Y];
                                    tsc_nodes[k][j][i]->u[UmomdensZ] +=  
                                       cur_grid->masstodens * weights[k][Z]*weights[j][Y]*weights[i][X] *
                                       cur_part->mom[Z];

                                    /* check max. density for stereo2 file */
                                    if(tsc_nodes[k][j][i]->u[Udens] > global.fdummy)
                                       global.fdummy = tsc_nodes[k][j][i]->u[Udens];
                                   }
                       }
                    }
                 }
              }
           }
        }
     }
   
   
// set initial magnetic field
   
#ifdef MHD
   for(cur_pquad=cur_grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next)
   {
      z = cur_pquad->z;
      for(cur_cquad = cur_pquad->loc;
          cur_cquad < cur_pquad->loc + cur_pquad->length; 
          cur_cquad++, z++)  
      {  
         for(icur_cquad  = cur_cquad; 
             icur_cquad != NULL; 
             icur_cquad  = icur_cquad->next)
         {
            y = icur_cquad->y;
            for(cur_nquad = icur_cquad->loc;  
                cur_nquad < icur_cquad->loc + icur_cquad->length; 
                cur_nquad++, y++) 
            { 
               for(icur_nquad  = cur_nquad; 
                   icur_nquad != NULL; 
                   icur_nquad  = icur_nquad->next)
               {
                  x = icur_nquad->x;
                  for(cur_node = icur_nquad->loc; 
                      cur_node < icur_nquad->loc + icur_nquad->length; 
                      cur_node++, x++)
                  {
                     
                     // set initial B field here!
                     cur_node->B[X] = 0.0;
                     cur_node->B[Y] = 1e-8;
                     cur_node->B[Z] = 0.0;
                     
                  }
               }
            }
         }
      }
   }
#endif

         
   /*----------------------------------------------------------------------
    * get total energy density, entropy, etc.
    *----------------------------------------------------------------------*/
   for(cur_pquad=cur_grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next)
     {
      z = cur_pquad->z;
      for(cur_cquad = cur_pquad->loc;
          cur_cquad < cur_pquad->loc + cur_pquad->length; 
          cur_cquad++, z++)  
        {  
         for(icur_cquad  = cur_cquad; 
             icur_cquad != NULL; 
             icur_cquad  = icur_cquad->next)
           {
            y = icur_cquad->y;
            for(cur_nquad = icur_cquad->loc;  
                cur_nquad < icur_cquad->loc + icur_cquad->length; 
                cur_nquad++, y++) 
              { 
               for(icur_nquad  = cur_nquad; 
                   icur_nquad != NULL; 
                   icur_nquad  = icur_nquad->next)
                 {
                  x = icur_nquad->x;
                  for(cur_node = icur_nquad->loc; 
                      cur_node < icur_nquad->loc + icur_nquad->length; 
                      cur_node++, x++)
                    {
                       
                     /* internal + kinetic energy density */
                     if(cur_node->u[Udens] > MACHINE_ZERO)
                     {
                          
#ifdef MHD
                        tsc_nodes[1][1][1] = cur_node;
                        get_TSCnodes(cur_grid, cur_pquad, icur_cquad, icur_nquad, tsc_nodes, &z, &y, &x);
                          
                        // calculate Bdens = 0.5*B^2 (magnetic energy density).
                        // idea: we get the cell-averaged values by averaging over pairs of opposing interfaces.
                        // remember that the B components are saved at the interfaces!
                        
                        if(test_tsc(tsc_nodes) == TRUE)  
                        Bdens = 0.125 * (pow2(cur_node->B[X] + tsc_nodes[1][1][2]->B[X]) +
                                         pow2(cur_node->B[Y] + tsc_nodes[1][2][1]->B[Y]) +
                                         pow2(cur_node->B[Z] + tsc_nodes[2][1][1]->B[Z]) );
                        else
                        Bdens = 0.5 * ( pow2(cur_node->B[X]) + pow2(cur_node->B[Y]) + pow2(cur_node->B[Z]) );
#else
                        Bdens = 0.0;
#endif
                        
                        e_therm = cur_node->u[Udens]*simu.e_init;
                        e_kin   = 0.5 * ( pow2(cur_node->u[UmomdensX])+
                                          pow2(cur_node->u[UmomdensY])+
                                          pow2(cur_node->u[UmomdensZ]) ) / cur_node->u[Udens];
                        e_tot   = e_therm + e_kin + Bdens;
                        
                        cur_node->u[UEdens] =  e_tot;                        
                        cur_node->u[Uedens] =  e_therm;

                        /* S = (gamma-1) edens / dens^(gamma-1) */
                        cur_node->u[Untrpy] = (simu.gamma-1.)*e_therm/pow(cur_node->u[Udens],(simu.gamma-1.));
                        
                       }
                     else
                       {
                        cur_node->u[UEdens] = 0.0;
                        cur_node->u[Uedens] = 0.0;
                        cur_node->u[Untrpy] = 0.0;
                       }                     
                    }
                 }
              }
           }
        }
     }
}


/*==============================================================================
 * set up a double one dimensional pancake test (HYDRO_TEST==5, HYDRO_TEST==6)
 *==============================================================================*/
void startrun_zeldovich_test(gridls *cur_grid)
{
   pqptr   cur_pquad;
   cqptr   cur_cquad, icur_cquad;
   nqptr   cur_nquad, icur_nquad;
   nptr    cur_node;
   partptr cur_part;
   int     i, j, k;
   long    npart;
   double  qx, qy, qz, x, y, z, vx, vy, vz, rho, phi, Fx;
   double  a, a_dot, Amp, Amp2;
   double  e_init;

   double  timecounter;
   int     no_timestep;
   
   /*-------------------------------------------
    *      Zeldovich wave characteristics
    * (according to Greg Bryan's PhD thesis...)
    *-------------------------------------------*/
   simu.boxsize   = 64.;
   simu.omega0    = 1.0;
   simu.omegab    = 1.0;
   simu.lambda0   = 0.0;
   simu.T_init    = 100.;
   simu.z_initial = 100.;
   simu.z_final   = 0.0;
   
   n1dim          = 256;  /* number of initial zones */

   z_caust        = 1.0;

#if (HYDRO_TEST == 5)
   z_caust2       = -1;
   strcpy(io.header.header,"SinglePancake");
#else
   z_caust2       = 1.45;
   strcpy(io.header.header,"DoublePancake");
#endif

   kwave          = 1. * TWOPI/simu.boxsize;
   kwave2         = 4. * TWOPI/simu.boxsize;

   Amp            = 1.+z_caust;
   Amp2           = 1.+z_caust2;

   a              = 1./(1.+simu.z_initial);
   a_dot          = H0*sqrt(1./a);             /* Friedmann eq. for SCDM cosmology */
   
   
   /* supercomoving conversion factors */
   x_fac = 1./simu.boxsize;
   v_fac = pow2(a)/(H0*simu.boxsize);
   e_fac = pow2(a)/pow2(H0*simu.boxsize);

   e_init = calc_e(simu.T_init);
   
   /* initialize AMIGA parameter */
   npart                = pow3(n1dim);
   simu.no_part         = npart;
   simu.mean_dens       = 1.0;
   simu.pmass           = simu.omega0 * rhoc0 * pow3(simu.boxsize)/simu.no_part;
   simu.FourPiG         = 1.5*simu.omega0;
   simu.gamma           = 5./3.;
   simu.e_init          = e_init * e_fac;
   simu.a_initial       = 1./(1.+simu.z_initial);
   simu.a_final         = 1./(1.+simu.z_final);
   simu.t_initial       = calc_t(simu.a_initial);
   simu.t_final         = calc_t(simu.a_final);
   global.no_timestep   = 0;
   
   /* we need to create the timeline */
   create_timeline(simu.a_initial/10., simu.a_final, &simu.timeline);
   
   
   fprintf(stderr,"\n\na       = %g %g\n",simu.a_initial,simu.a_final);
   fprintf(stderr,"t       = %g %g\n",(calc_t(simu.a_initial)),(calc_t(simu.a_final)));
   fprintf(stderr,"t       = %g %g\n",simu.t_initial,simu.t_final);
   fprintf(stderr,"super_t = %g %g\n",calc_super_t(simu.a_initial),calc_super_t(simu.a_final));
   fprintf(stderr,"a_super = %g %g\n\n",calc_super_a(calc_super_t(simu.a_initial)),calc_super_a(calc_super_t(simu.a_final)));
   
   /* some technical aspects */
   global.fst_part      = c_part(npart);
   cur_part             = global.fst_part;
   cur_grid->masstodens = pow3((double)cur_grid->l1dim)/(double)npart;   
   
   
#if (HYDRO_TEST == 5)
   fprintf(stderr,"startrun_single_pancake_test:  a = %g\n\n",a);
#else
   fprintf(stderr,"startrun_double_pancake_test:  a = %g\n\n",a);
#endif
   
   /* two Zeldovich waves */
   for(k=0;k<n1dim;k++)
     {
      for(j=0;j<n1dim;j++)
        {
         for(i=0; i<n1dim; i++)
           {
            /* Lagrangian coordinates */
            qx = ((double)i+0.5)/(double)n1dim * simu.boxsize;
            qy = ((double)j+0.5)/(double)n1dim * simu.boxsize;
            qz = ((double)k+0.5)/(double)n1dim * simu.boxsize;

            /* x-coordinate */
            x  = qx + a     * (Amp) * cos(kwave*qx)/kwave + a     * (Amp2) * cos(kwave2*qx-TWOPI/4.)/kwave2;
            vx =      a_dot * (Amp) * cos(kwave*qx)/kwave + a_dot * (Amp2) * cos(kwave2*qx-TWOPI/4.)/kwave2;
            
            /* y-coordinate */
            y  = qy;
            vy = 0.0;
            
            /* z-coordinate */
            z  = qz;
            vz = 0.0;
            

            /* transfer to AMIGA's particle structure */
            cur_part->pos[X] = fmod(x_fac*x + 1.0 , 1.0);
            cur_part->mom[X] = v_fac * vx;
            
            cur_part->pos[Y] = fmod(x_fac*y + 1.0 , 1.0);
            cur_part->mom[Y] = v_fac * vy;
            
            cur_part->pos[Z] = fmod(x_fac*z + 1.0 , 1.0);
            cur_part->mom[Z] = v_fac * vz;
            
            /* move to next particle... */
            cur_part++;            
           }
        }
     }
   
   /* check for shell crossing -> particles should be ordered in x-direction! */
   cur_part = global.fst_part;
   for(i=0; i<n1dim-1; i++)
     {
      qx = (((double)i+0.5)/(double)n1dim);
      if(cur_part->pos[X] > (cur_part+1)->pos[X] && cur_part->pos[X]-(cur_part+1)->pos[X]<0.5)
         fprintf(stderr,"shell crossing: %d  %g %g  (%g)\n",i,cur_part->pos[X],(cur_part+1)->pos[X],qx);
      
      cur_part++;
     }
   
#ifdef MULTIMASS
   io.header.multi_mass       = 1;
#else
   io.header.multi_mass       = 0;
#endif
#ifdef DOUBLE
   io.header.double_precision = 1;
#else
   io.header.double_precision = 0;
#endif
   io.header.no_part          = simu.no_part;
   io.header.no_species       = 1;
   io.header.no_vpart         = (double)simu.no_part;
   io.header.timestep         = 0.1;
   io.header.no_timestep      = 0;
   io.header.boxsize          = simu.boxsize;
   io.header.omega0           = simu.omega0;
   io.header.lambda0          = simu.lambda0;
   io.header.pmass            = simu.pmass;
   io.header.cur_reflevel     = 0.0;
   io.header.cur_frcres       = simu.boxsize/(double)simu.NGRID_DOM;
   io.header.a_initial        = simu.a_initial;
   io.header.a_current        = simu.a_initial;
   io.header.K_initial        = 0.0;
   io.header.K_current        = 0.0;
   io.header.U_initial        = 0.0;
   io.header.U_current        = 0.0;
   io.header.Eintegral        = 0.0;
   io.header.Econst           = 0.0;
   /* ignore all paramXYZ values... */
   io.header.version          = VERSION;
   io.header.built            = BUILT;
   io.header.omegab           = simu.omegab;
   io.header.gamma            = simu.gamma;
   io.header.H_frac           = simu.H_frac;
   io.header.T_init           = simu.T_init;
   
   timecounter        = calc_super_t(simu.a_initial);
   no_timestep        = 0;
   global.z           = simu.z_initial;
   global.no_timestep = no_timestep;
   global.no_part       = npart;
#ifdef WRITE_ZELDOVICH
   output(io.outfile_prefix, 0);
#endif
   
#ifdef WRITE_ZELDOVICH_TERMINATE
   exit(0);
#endif
   
   /* generate linked-list */
   ll(npart, global.fst_part, cur_grid);
   
   /* assign particles to grid: store results in u[] */
   zero_udens(cur_grid);
   init_udens(cur_grid); 
   
   free(global.fst_part);
}


/*==================================================================================
 * set up the non-rotating cloud (Ziegler, 2005, A&A 345, 285)
 *==================================================================================*/
void startrun_nonrot_cloud_test(gridls *cur_grid)
{
   int     i, j, k;
   FILE   *fpin, *fpout;
   double  timecounter;
   int     no_timestep;
   pqptr   cur_pquad;
   cqptr   cur_cquad, icur_cquad;
   nqptr   cur_nquad, icur_nquad;
   nptr    cur_node;
   long    x, y, z;
   
   double  dx, dy, dz;
   
   double  R_cl, Rho_cl, tau_ff, Grav_SI;
   

   
   /* this is all non-solicited information! */
   simu.no_part         = 50000;
   simu.boxsize         = 1.;
   simu.a_initial       = 1.0;
   simu.a_final         = 1.0;
   simu.z_initial       = (double)1.0/simu.a_initial - (double)1.0;
   simu.t_initial       = calc_t(simu.a_initial);
   simu.omega0          = 1.0;
   simu.lambda0         = 0.0;
   simu.pmass           = simu.omega0 * rhoc0 * pow3(simu.boxsize)/simu.no_part;
   
   /* this is solicited information! */
   R_cl                 = 7.8E13;       // m
   Rho_cl               = 1E-12;        // kg/m^3
   tau_ff               = 6.645E10;     // s
   Grav_SI              = 6.6726E-11;   // m^3 kg^-1 s^-2
   simu.mean_dens       = 0.0;
   simu.FourPiG         = Rho_cl*pow2(tau_ff)*4*PI*Grav_SI;
   simu.gamma           = 1.001;
   simu.T_init          = 3.63;
   simu.e_init          = kB_per_mp*1.0E6 * pow2(tau_ff/R_cl) * 1.5*simu.T_init;
   
   
   /* we need to create the timeline */
   create_timeline(simu.a_initial/10., simu.a_final, &simu.timeline);
   
   /*----------------------------------------------------------------------
    * get total energy density and entropy
    *----------------------------------------------------------------------*/
   for(cur_pquad=cur_grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next)
     {
      z = cur_pquad->z;
      for(cur_cquad = cur_pquad->loc;
          cur_cquad < cur_pquad->loc + cur_pquad->length; 
          cur_cquad++, z++)  
        {  
         for(icur_cquad  = cur_cquad; 
             icur_cquad != NULL; 
             icur_cquad  = icur_cquad->next)
           {
            y = icur_cquad->y;
            for(cur_nquad = icur_cquad->loc;  
                cur_nquad < icur_cquad->loc + icur_cquad->length; 
                cur_nquad++, y++) 
              { 
               for(icur_nquad  = cur_nquad; 
                   icur_nquad != NULL; 
                   icur_nquad  = icur_nquad->next)
                 {
                  x = icur_nquad->x;
                  for(cur_node = icur_nquad->loc; 
                      cur_node < icur_nquad->loc + icur_nquad->length; 
                      cur_node++, x++)
                    {
                     
                     dx = ((double)x - ((double)cur_grid->l1dim/2 - 0.5));
                     dy = ((double)y - ((double)cur_grid->l1dim/2 - 0.5));
                     dz = ((double)z - ((double)cur_grid->l1dim/2 - 0.5));
                     
                     if((pow2(dx)+pow2(dy)+pow2(dz)) < pow2((double)cur_grid->l1dim/6.))
                        cur_node->u[Udens]     = 1.0;
                     else
                        cur_node->u[Udens]     = 0.01;
                     cur_node->u[UmomdensX] = 0.0;
                     cur_node->u[UmomdensY] = 0.0;
                     cur_node->u[UmomdensZ] = 0.0;
                     cur_node->u[Uedens]    = simu.e_init*cur_node->u[Udens];
                     cur_node->u[UEdens]    = cur_node->u[Uedens];
                     cur_node->u[Untrpy]    = (simu.gamma-1.)*cur_node->u[UEdens]/pow(cur_node->u[Udens], (simu.gamma-1.));
                     
                     /* CAREFUL: add_gasdens() requires the density in u_tmp[]! */
                     cur_node->u_tmp[Udens] = cur_node->u[Udens];
                    }
                 }
              }
           }
        }
     }
         
}


/*==============================================================================
 * set up the IC's for some simple tests incl. the BlastWave and the ShockTube
 *==============================================================================*/
void startrun_source_free_test(gridls *cur_grid)
{
   pqptr   cur_pquad;
   cqptr   cur_cquad, icur_cquad;
   nqptr   cur_nquad, icur_nquad;
   nptr    cur_node;
   long    x, y, z, left, right;
   int     idim, i;
   double  rho0_left, rho0_right, p0_left, p0_right;
   double  vx0_left, vx0_right, vy0_left, vy0_right, vz0_left, vz0_right;
   double  Bx0_left, Bx0_right, By0_left, By0_right, Bz0_left, Bz0_right;
   double  dx, dy, dz;
   
#ifdef ORSZAG_TANG
   double x_value, y_value, Bx, By, l;
#endif //ORSZAG_TANG
      
#ifdef DENSITY_CONSTANT
   simu.gamma = 5./3.;
#endif
#ifdef DENSITY_STEP
   simu.gamma = 1.4;
#endif
 
#ifdef SHOCK_TUBE
 
#ifdef SHOCK_RICKER
   /* Ricker (2000) */
   simu.gamma = 1.4;
   
   rho0_left  = 1.0;
   vx0_left   = 0.0;
   vy0_left   = 0.0;
   vz0_left   = 0.0;
   p0_left    = 1.0;
   
   rho0_right = 0.125;
   vx0_right  = 0.;
   vy0_right  = 0.;
   vz0_right  = 0.;
   p0_right   = 0.1;
#endif
 
#ifdef SHOCK_RYU
   /* Ryu et al. (1993) */
   simu.gamma = 1.4;

   rho0_left  = 1.0;
   vx0_left   = 0.0;
   vy0_left   = 0.0;
   vz0_left   = 0.0;
   p0_left    = 1.0;

   rho0_right = 1.0;
   vx0_right  = 0.;
   vy0_right  = 0.;
   vz0_right  = 0.;
   p0_right   = 0.1;
#endif
   
#ifdef SHOCK_MHD   
   
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   // what follows are different sets of 1D MHD test problems (see Ziegler 2004, Dolag/Stasyszyn 2008) /////////////
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
   // TEST Nr. 1A (Ryu 1995, Dolag/Stasyszyn 2008)
   
//   simu.gamma = 5./3.;
//   
//   rho0_left  = 1.;
//   vx0_left   = 10.;
//   vy0_left   = 0.;
//   vz0_left   = 0.;
//   p0_left    = 20.;
//   Bx0_left   = 5.0/sqrt(4.0*PI);
//   By0_left   = 5.0/sqrt(4.0*PI);
//   Bz0_left   = 0.;
//   
//   rho0_right = 1.;
//   vx0_right  = -10.;
//   vy0_right  = 0.;
//   vz0_right  = 0.;
//   p0_right   = 1.;
//   Bx0_right  = 5.0/sqrt(4.0*PI);
//   By0_right  = 5.0/sqrt(4.0*PI);
//   Bz0_right  = 0.;
   
   //set of parameters for Brio-Wu Riemann problem (Ziegler 2004)

   simu.gamma = 2.;
   
   rho0_left  = 1.;
   vx0_left   = 0.;
   vy0_left   = 0.;
   vz0_left   = 0.;
   p0_left    = 1.;
   Bx0_left   = 0.75;
   By0_left   = 1.;
   Bz0_left   = 0.;
   
   rho0_right = 0.125;
   vx0_right  = 0.;
   vy0_right  = 0.;
   vz0_right  = 0.;
   p0_right   = 0.1;
   Bx0_right  = 0.75;
   By0_right  = -1.;
   Bz0_right  = 0.;
   
   // set of parameters for SHOCK_RICKER with zero magnetic field - end time: 0.206
   
//   simu.gamma = 1.4;    
//   
//   rho0_left  = 1.0;
//   vx0_left   = 0.0;
//   vy0_left   = 0.0;
//   vz0_left   = 0.0;
//   p0_left    = 1.0;
//   Bx0_left   = 0.;
//   By0_left   = 0.;
//   Bz0_left   = 0.;
//   
//   rho0_right = 0.125;
//   vx0_right  = 0.;
//   vy0_right  = 0.;
//   vz0_right  = 0.;
//   p0_right   = 0.1;
//   Bx0_right  = 0.;
//   By0_right  = 0.;
//   Bz0_right  = 0.;  
   
   
#endif // SHOCK_MHD
 
#endif /* SHOCK_TUBE */
   
#ifdef ORSZAG_TANG  // for the funny propeller-shaped image found in the literature set final time to 0.5!
   simu.gamma = 5./3.;
   l = cur_grid->l1dim;  // box size length, needed to set up staggered B initial conditions
   // all the rest is set directly in the loop
#endif //ORSZAG_TANG
 
#ifdef BLAST_WAVE
   /* Ziegler (2004) */
   simu.gamma = 1.4;
#endif
   
   left  = 0;                      //cur_grid->l1dim/2-cur_grid->l1dim/4;
   right = cur_grid->l1dim/2;      //+cur_grid->l1dim/4;
   
   
   /*----------------------------------------------------------------------
    * loop over all nodes in cur_grid
    *----------------------------------------------------------------------*/
   for(cur_pquad=cur_grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next)
     {
      z = cur_pquad->z;
      for(cur_cquad = cur_pquad->loc;
          cur_cquad < cur_pquad->loc + cur_pquad->length; 
          cur_cquad++, z++)  
        {  
         for(icur_cquad  = cur_cquad; 
             icur_cquad != NULL; 
             icur_cquad  = icur_cquad->next)
           {
            y = icur_cquad->y;
            for(cur_nquad = icur_cquad->loc;  
                cur_nquad < icur_cquad->loc + icur_cquad->length; 
                cur_nquad++, y++) 
              { 
               for(icur_nquad  = cur_nquad; 
                   icur_nquad != NULL; 
                   icur_nquad  = icur_nquad->next)
                 {
                  x = icur_nquad->x;
                  for(cur_node = icur_nquad->loc; 
                      cur_node < icur_nquad->loc + icur_nquad->length; 
                      cur_node++, x++)
                    {
                     
#ifdef DENSITY_CONSTANT
                     cur_node->u[Udens]     = 12345.6789;
                     cur_node->u[UmomdensX] = 0.0;
                     cur_node->u[UmomdensY] = 0.0;
                     cur_node->u[UmomdensZ] = 0.0;
                     cur_node->u[UEdens]    = 0.0;
                     cur_node->u[Untrpy]    = 0.0;
                     cur_node->u[Uedens]    = 0.0;
#endif //DENSITY_CONSTANT
                     
#ifdef DENSITY_STEP
                     if(left <= x && x < right)
                       {
                        cur_node->u[Udens]  = 1.0;
                       }
                     else
                       {
                        cur_node->u[Udens]  = 0.125;
                       }
                     cur_node->u[UmomdensX] = 0.0;
                     cur_node->u[UmomdensY] = 0.0;
                     cur_node->u[UmomdensZ] = 0.0;
                     cur_node->u[UEdens]    = 0.0;
                     cur_node->u[Untrpy]    = 0.0;
                     cur_node->u[Uedens]    = 0.0;
#endif // DENSITY_STEP
                     
#ifdef SHOCK_TUBE
                     // Sod's shock tube problem
#if (HYDRO_TEST==0 || HYDRO_TEST==9)
                     if(left <= x && x < right)
#endif
#if (HYDRO_TEST==1 || HYDRO_TEST==10)
                     if(left <= y && y < right)
#endif
#if (HYDRO_TEST==2 || HYDRO_TEST==11)
                     if(left <= z && z < right)
#endif
                       {
                        cur_node->u[Udens]     = rho0_left;
                        cur_node->u[UmomdensX] = vx0_left;
                        cur_node->u[UmomdensY] = vy0_left;
                        cur_node->u[UmomdensZ] = vz0_left;
                        cur_node->u[Uedens]    = p0_left/(simu.gamma-1.0);
                        cur_node->u[Untrpy]    = (simu.gamma-1)*cur_node->u[Uedens]/pow(cur_node->u[Udens],(simu.gamma-1));
                        cur_node->u[UEdens]    = cur_node->u[Uedens] + 
                           0.5*(
                                pow2(cur_node->u[UmomdensX])+
                                pow2(cur_node->u[UmomdensY])+
                                pow2(cur_node->u[UmomdensZ])
                                )/cur_node->u[Udens];
#ifdef SHOCK_MHD
                        cur_node->B[X]         = Bx0_left;
                        cur_node->B[Y]         = By0_left;
                        cur_node->B[Z]         = Bz0_left;
                        cur_node->u[UEdens]   += (pow2(cur_node->B[X])+pow2(cur_node->B[Y])+pow2(cur_node->B[Z]))/2.;
#endif //SHOCK_MHD
                       }
                     else
                       {
                        cur_node->u[Udens]     = rho0_right;
                        cur_node->u[UmomdensX] = vx0_right;
                        cur_node->u[UmomdensY] = vy0_right;
                        cur_node->u[UmomdensZ] = vz0_right;
                        cur_node->u[Uedens]    = p0_right/(simu.gamma-1.0);
                        cur_node->u[Untrpy]    = (simu.gamma-1)*cur_node->u[Uedens]/pow(cur_node->u[Udens],(simu.gamma-1));
                        cur_node->u[UEdens]    = cur_node->u[Uedens] + 
                           0.5*(
                                pow2(cur_node->u[UmomdensX])+
                                pow2(cur_node->u[UmomdensY])+
                                pow2(cur_node->u[UmomdensZ])
                                )/cur_node->u[Udens];
#ifdef SHOCK_MHD
                        cur_node->B[X]         = Bx0_right;
                        cur_node->B[Y]         = By0_right;
                        cur_node->B[Z]         = Bz0_right;
                        cur_node->u[UEdens]   += (pow2(cur_node->B[X])+pow2(cur_node->B[Y])+pow2(cur_node->B[Z]))/2.;
#endif //SHOCK_MHD
                       }
#endif //SHOCK_TUBE
                       
// new MHD test: Orszag-Tang vortex from www.astro.princeton.edu/~jstone/tests/orszag-tang/pagesource.html ////////////////
                       
#ifdef ORSZAG_TANG
                       // most important thing here: don't forget about the fact that the B field is defined
                       // on the STAGGERED grid (= on the interfaces) when setting up initial conditions.
                       
                       // x,y,z can have values between 0 and l1dim-1; but we need (x,y) = [0,1]^2. so we do a
                       // conversion: x_value and y_value are now the "real" x and y values which lie between 0 and 1
                       
                       x_value = (double)x / l;     // l was set to cur_node->l1dim above
                       y_value = (double)y / l; 
                       Bx = -0.28209479177 * sin(2.*PI* y_value ); 
                       By = 0.28209479177 *sin(4.*PI* x_value );   //factor 1/sqrt(4*PI) = 0.28209479177                       
                       //Pairs of opposing interface B values are identical FOR THIS TEST ONLY (because Bx=f(y), By=f(x))
                  
                       cur_node->u[Udens]     = 25./(36.*PI);
                       cur_node->u[UmomdensX] = -sin(2.*PI*y_value) * cur_node->u[Udens];
                       cur_node->u[UmomdensY] = sin(2.*PI*x_value) * cur_node->u[Udens];
                       cur_node->u[UmomdensZ] = 0.;
                       cur_node->u[Uedens]    = (5./(12.*PI))/(simu.gamma-1.0);  // initial pressure p0 = 5./(12.*PI)
                       cur_node->u[Untrpy]    = (simu.gamma-1)*cur_node->u[Uedens]/pow(cur_node->u[Udens],(simu.gamma-1));
                       cur_node->u[UEdens]    = cur_node->u[Uedens] + 
                       0.5*(
                            pow2(cur_node->u[UmomdensX])+
                            pow2(cur_node->u[UmomdensY])+
                            pow2(cur_node->u[UmomdensZ])
                         )/cur_node->u[Udens];
                       cur_node->B[X]         = Bx;
                       cur_node->B[Y]         = By;
                       cur_node->B[Z]         = 0.;
                       cur_node->u[UEdens]   += 0.5 * ( Bx*Bx + By*By );                         
                       
#endif //ORSZAG_TANG
                       
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        
#ifdef BLAST_WAVE
                       /* Gaussian smoothing of "blast in the centre" over 1.5 cells */
                       dx = (x - (cur_grid->l1dim/2 - 0.5));
                       dy = (y - (cur_grid->l1dim/2 - 0.5));
                       dz = (z - (cur_grid->l1dim/2 - 0.5));
                       
#ifdef BLAST_RICKER00
                       /* Ricker (2000) */
                       cur_node->u[Udens]     = 1.0;
                       cur_node->u[UmomdensX] = 0.0;
                       cur_node->u[UmomdensY] = 0.0;
                       cur_node->u[UmomdensZ] = 0.0;
                       cur_node->u[UEdens]    = 1e-5+1.1*exp(-(pow2(dx)+pow2(dy)+pow2(dz)) / 1.5);
                       cur_node->u[Untrpy]    = (simu.gamma-1)*cur_node->u[UEdens]/pow(cur_node->u[Udens],(simu.gamma-1));
                       cur_node->u[Uedens]    = cur_node->u[UEdens];
#endif /* BLAST_RICKER00 */ 
                       
#ifdef BLAST_ZIEGLER04
                       /* Ziegler (2004) */
                       cur_node->u[Udens]     = 1.0;
                       cur_node->u[UmomdensX] = 0.0;
                       cur_node->u[UmomdensY] = 0.0;
                       cur_node->u[UmomdensZ] = 0.0;
                       if(pow2(dx)+pow2(dy)+pow2(dz) < pow2(0.1*cur_grid->l1dim))
                       {
                          cur_node->u[Uedens] = cur_node->u[Udens] * 1.0E4/(simu.gamma-1.);
                       }
                       else
                       {
                          cur_node->u[Uedens] = cur_node->u[Udens] * 1.0/(simu.gamma-1.);
                       }
                       cur_node->u[UEdens]    = cur_node->u[Uedens];
                       cur_node->u[Untrpy]    = (simu.gamma-1)*cur_node->u[UEdens]/pow(cur_node->u[Udens],(simu.gamma-1));
                       
                       
                       /* THIS IS ACTUALLY AN MHD TEST AND HENCE WE ALSO NEED TO INITIALIZE THE B-FIELD*/
                       
                       
#endif /* BLAST_ZIEGLER04 */ 
#endif // BLAST_WAVE
                     

                    }
                 }
              }
           }
        }
     }
   
   /* double check hydro-variables*/
   write_hydro(cur_grid);

}

/*==============================================================================
 *               Shock Tube and Blast Wave Hydro-Tests
 *==============================================================================*/
void test_hydro_srcfree()
{
   double timecounter, timecounter_final, timestep;
   int    no_timestep;
   
   timecounter               = 0.0;
   no_timestep               = 0;
   global.no_timestep        = no_timestep;
   global.dom_grid->timestep = timestep;
   global.cfl_speed          = 0.0;
   
#ifdef SHOCK_TUBE
 timestep                  = 0.0005;
 
#ifndef SHOCK_RYU
#ifndef SHOCK_RICKER
#ifndef SHOCK_MHD
   fprintf(stderr,"please either define SHOCK_RYU, SHOCK_RICKER, or SHOCK_MHD in define.h\n");
   exit(0);
#endif
#endif
#endif
 
#ifdef SHOCK_RICKER
   /* Ricker (2000) */
   timecounter_final         = 0.206;
#endif
 
#ifdef SHOCK_RYU
   /* Ryu et al (1993) */
   timecounter_final         = 0.5;
#endif
 
#endif /* SHOCK_TUBE */
   
#ifdef ORSZAG_TANG
   timestep                  = 0.0005;
   global.a                  = 1.0;
   timecounter_final         = 0.5;
#endif //ORSZAG_TANG
 
#ifdef BLAST_WAVE
   
#ifdef BLAST_ZIEGLER04
   /* Ziegler (2004) */
   timestep                  = 1E-5;
   timecounter_final         = 2.5E-3;
#endif
#ifdef BLAST_RICKER00
   /* Ricker (2000) */
   timestep                  = 0.00001;
   timecounter_final         = 20.0;
#endif

#endif // BLAST_WAVE
   
#ifdef SHOCK_MHD
   global.a                  = 1.0;
   timecounter_final         = 0.1;      // final time for MHD shock tube
#endif
   
   /* initialize some test scenario*/
   startrun_source_free_test(global.dom_grid);
   
   /* integrate the Euler equations in time */
   while((timecounter_final-timecounter) > ZERO)
     {
      /* obtain hydro-variables: t -> t+timestep */
      solve_dom_hydro(timestep);
      
      /* increment integration variable */
      no_timestep++;
      timecounter += timestep;
      
      /* update globally accessible integration variable */
      global.no_timestep        = no_timestep;
      global.dom_grid->timestep = timestep;
      global.super_t            = timecounter;
#ifdef MHD // should do the same as #if ( defined(SHOCK_MHD) || defined(ORSZAG_TANG) )
      global.a                  = 1.0;
#else
      global.a                  = calc_super_a(timecounter);
#endif
      global.t                  = calc_t(global.a);
      global.z                  = 1./global.a - 1.;
      
      fprintf(stderr,"step %8d done:    time=%12.8g  timestep=%12.8g   cfl_timestep=%12.8g   cfl_speed=%12.8g\n",
              no_timestep, timecounter, timestep, 
              global.dom_grid->spacing/global.cfl_speed,global.cfl_speed);
      
      /* dump newly calculated values */
//      if((no_timestep % 10) == 0)
//         write_hydro(global.dom_grid);
      
      timestep = CFL_TUNE * global.dom_grid->spacing/global.cfl_speed;
        
      // reduce the timestep to halt exactly at timecounter_final (at the end of the simulation)   
      if ( timestep > timecounter_final - timecounter)
         timestep = timecounter_final - timecounter;
        
     }
   
   write_hydro(global.dom_grid);
}

/*==============================================================================
 *                       collapse of a non-rotating cloud
 *==============================================================================*/
void test_nonrot_cloud(gridls *grid_list)
{
   double  timecounter, timecounter_final, timestep, a;
   int     no_timestep;
   gridls *cur_grid;
   
   cur_grid = global.dom_grid;
   
   fprintf(stderr,"\n\n\n\n*===================================================*\n");
   fprintf(stderr,"*               Non-Rotating Cloud Collapse\n");
   
   startrun_nonrot_cloud_test(cur_grid);
   
   fprintf(stderr,"                   (FourPiG = %g)\n",simu.FourPiG);
   fprintf(stderr,"*===================================================*\n");

   timecounter               = 0.0;   // in units of tau_ff
   timecounter_final         = 0.996; // in units of tau_ff
   timestep                  = 0.01;
   no_timestep               = 0;
   
   global.no_timestep        = no_timestep;
   global.dom_grid->timestep = timestep;
   global.cfl_speed          = 0.0;
   
   global.super_t            = timecounter;
   global.a                  = calc_super_a(timecounter);
   global.t                  = calc_t(global.a);
   global.z                  = 1./global.a - 1.;
   
   cur_grid->timecounter     = timecounter;
   
   fprintf(stderr,"\nnonrot_cloud_test from       t = %g -> %g (dt=%g)\n",timecounter,timecounter_final,timestep);
   
   
   
   zero_dens(cur_grid);
   add_gasdens(cur_grid);
   solve_gravity(grid_list, global.domgrid_no);
   calc_forces(cur_grid);
   
   write_hydro(cur_grid);
   
   
   /* integrate the Euler equations in time */
   while((timecounter_final-timecounter) > ZERO)
     {
      calc_Flux(cur_grid);
      advect_hydro(cur_grid, timestep);
      ave_hydro(cur_grid);
      
      cur_grid->timecounter += timestep/2;

      zero_dens(cur_grid);
      add_gasdens(cur_grid);
      solve_gravity(grid_list, global.domgrid_no);
      calc_forces(cur_grid);
      
      calc_Flux(cur_grid);
      advect_gravity_hydro(cur_grid, timecounter, timestep);
      
      cur_grid->timecounter += timestep/2;

      /* update timecounter and the likes... */
      no_timestep++;
      timecounter              += timestep;
      a                         = calc_super_a(timecounter);
      global.no_timestep        = no_timestep;
      global.dom_grid->timestep = timestep;
      global.super_t            = timecounter;
      global.a                  = a;
      global.t                  = calc_t(a);
      global.z                  = 1./a - 1.;
      
      
      /* be verbose */
      fprintf(stderr,"step %8d done:    time=%12.8g  timestep=%12.8g   cfl_timestep=%12.8g   cfl_speed=%12.8g\n",
              no_timestep, timecounter, timestep, 
              global.dom_grid->spacing/global.cfl_speed,global.cfl_speed);
      
      
      /* adjust timestep */
      if(global.cfl_speed > MACHINE_ZERO)
         if(CFL_TUNE * global.dom_grid->spacing/global.cfl_speed < timestep)
            timestep = CFL_TUNE * global.dom_grid->spacing/global.cfl_speed;
      
      
      /* dump newly calculated values? */
      if((no_timestep % 10) == 0)
       write_hydro(global.dom_grid);

      if(fabs(timecounter-0.945) < 0.01 || 
         fabs(timecounter-0.955) < 0.01 ||
         fabs(timecounter-0.955) < 0.01 ||
         fabs(timecounter-0.965) < 0.01 ||
         fabs(timecounter-0.975) < 0.01 ||
         fabs(timecounter-0.985) < 0.01 ||
         fabs(timecounter-0.995) < 0.01   ) 
       write_hydro(global.dom_grid);
      
     } /* while-loop */     
   
   
   
   write_hydro(global.dom_grid);
   
}


/*==============================================================================
 *                            Zeldovich Pancake
 *==============================================================================*/
void test_zeldovich(gridls *grid_list)
{
   double  timecounter, timecounter_final, timestep, a;
   int     no_timestep;
   gridls *cur_grid;
   
   cur_grid = global.dom_grid;
   
   fprintf(stderr,"\n\n\n\n*===================================================*\n");
   fprintf(stderr,"*                 Zeldovich Pancake\n");
#if (HYDRO_TEST==5)
   fprintf(stderr,"                   (single pancake)\n");
#else
   fprintf(stderr,"                   (double pancake)\n");
#endif

   
   startrun_zeldovich_test(cur_grid);
   
   fprintf(stderr,"*===================================================*\n");

   timecounter               = calc_super_t(simu.a_initial);
   timecounter_final         = calc_super_t(simu.a_final);
   timestep                  = 0.1;
   no_timestep               = 0;
   
   global.no_timestep        = no_timestep;
   global.dom_grid->timestep = timestep;
   global.cfl_speed          = 0.0;
   
   global.super_t            = timecounter;
   global.a                  = calc_super_a(timecounter);
   global.t                  = calc_t(global.a);
   global.z                  = 1./global.a - 1.;
   
   cur_grid->timecounter     = timecounter;

   a = calc_super_a(timecounter);

   fprintf(stderr,"\ncomplete_hydro_test from       a = %g -> %g\n",simu.a_initial,simu.a_final);
   fprintf(stderr,"complete_hydro_test from       t = %g -> %g\n",simu.t_initial,simu.t_final);
   fprintf(stderr,"complete_hydro_test from super_t = %g -> %g\n",timecounter,timecounter_final);
   fprintf(stderr,"complete_hydro_test from a_super = %g -> %g\n\n",
           calc_super_a(timecounter),calc_super_a(timecounter_final));
      
   
   //write_ZW(timecounter, no_timestep);
   
   
   zero_dens(cur_grid);
   add_gasdens(cur_grid);
   solve_gravity(grid_list, global.domgrid_no);
   calc_forces(cur_grid);
   
   /* dump initial conditions to file */
   write_hydro(cur_grid);
   
   /* integrate the Euler equations in time */
   while((timecounter_final-timecounter) > ZERO)
     {
      calc_Flux(cur_grid);
      advect_hydro(cur_grid, timestep);      
      ave_hydro(cur_grid);
      
      cur_grid->timecounter += timestep/2;
      
      zero_dens(cur_grid);
      add_gasdens(cur_grid);
      solve_gravity(grid_list, global.domgrid_no);
      calc_forces(cur_grid);

      calc_Flux(cur_grid);
      advect_gravity_hydro(cur_grid, timecounter, timestep);
      
      cur_grid->timecounter += timestep/2;
      
      /* update timecounter and the likes... */
      no_timestep++;
      timecounter              += timestep;
      a                         = calc_super_a(timecounter);
      global.no_timestep        = no_timestep;
      global.dom_grid->timestep = timestep;
      global.super_t            = timecounter;
      global.a                  = a;
      global.t                  = calc_t(a);
      global.z                  = 1./a - 1.;
      
      
      /* be verbose */
      fprintf(stderr,"step %8d done:    time=%12.8g->%12.8g  timestep=%12.8g   cfl_timestep=%12.8g   cfl_speed=%12.8g   z=%g\n",
              no_timestep, timecounter-timestep, timecounter, timestep, 
              a*global.dom_grid->spacing/global.cfl_speed,global.cfl_speed,1./a-1.);
      
      /* calculate new timestep according to various criteria */
      timestep = adjust_timestep(timecounter, timestep);   
      
      /* dump newly calculated values? */
      if(fabs(1./a-1. - 1.45) < 0.1 || fabs(1./a-1. - 2.) < 0.1)
        {
         write_hydro(global.dom_grid);
         //write_ZW(timecounter, no_timestep);
        }
      
     } /* while-loop */     
   
   
   
   write_hydro(global.dom_grid);
   //write_ZW(timecounter, no_timestep);   
   
}


/*==============================================================================
 * simple wrapper for all HYDRO_TESTs...
 *
 * HYDRO_TEST == 0                    X
 * HYDRO_TEST == 1    -> ShockTube in Y direction
 * HYDRO_TEST == 2                    Z
 * HYDRO_TEST == 3    -> BlastWave
 * HYDRO_TEST == 4    -> isothermal gas sphere
 * HYDRO_TEST == 5    -> ZeldovichWave (single pancake)
 * HYDRO_TEST == 6    -> ZeldovichWave (double pancake)
 * HYDRO_TEST == 7    -> testing implementation in step()
 * HYDRO_TEST == 8    -> non-rotating cloud (Ziegler U., 2005, A&A 435, 385)
 *
 *==============================================================================*/
void hydro_test(gridls *grid_list)
{
#if (HYDRO_TEST==8)
   test_nonrot_cloud(grid_list);
#endif
   
#if (HYDRO_TEST==5 || HYDRO_TEST==6)
   /*---------------------------------------------------------------------
    *               Zeldovich Pancake Tests ala Greg Bryan
    *---------------------------------------------------------------------*/
   test_zeldovich(grid_list);
#endif /* HYDRO_TEST==5 || 6 */
   
#if (HYDRO_TEST==4)
   /*---------------------------------------------------------------------
   *                        Gas Sphere Collapse
   *---------------------------------------------------------------------*/
   test_gas_sphere(grid_list);
#endif /* HYDRO_TEST==4 */

#if (HYDRO_TEST<4)
   /*---------------------------------------------------------------------
    *                 Shock Tube & Blast Wave Hydro-Tests
    *---------------------------------------------------------------------*/
   test_hydro_srcfree();
#endif /* HYDRO_TEST<4 */
   
#if (HYDRO_TEST==9 || HYDRO_TEST==10 || HYDRO_TEST==11)
   /*---------------------------------------------------------------------
    *                        MHD Shock Tube Test
    *---------------------------------------------------------------------*/
   test_hydro_srcfree();
#endif /* HYDRO_TEST==9 || HYDRO_TEST==10 || HYDRO_TEST==11 */
   
#if (HYDRO_TEST==12)
   /*---------------------------------------------------------------------
    *              3D-periodic MHD tests (e.g. Orszag-Tang)
    *---------------------------------------------------------------------*/
   test_hydro_srcfree();
#endif /* HYDRO_TEST==12 */
      
   common_terminate(0);
}
#endif /* HYDRO */


