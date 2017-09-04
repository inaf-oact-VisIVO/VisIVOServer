#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef WITH_OPENMPM
#include <omp.h>
#endif

void dummy_hydro()
{
}

#ifdef HYDRO

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "mhd.h"
#include "../libutility/utility.h"
#include "../libgrids/grids.h"

/*===============================================================================
*
* the numerical scheme to solve the equations of (magneto-)hydrodynamics
* (used with solve_hydro.c and found in hydro.c) is freely adapted from 
*
*     Ziegler U., Journal of Computational Physics 196 (2004) 393-416
*
* and based upon the seminal paper by (and hence the name KNP)
*
*    Kurganov, Noelle & Petrova, SIAM J. Sci. Comput. 23 (2001) 707
*
*
* some subroutines may seem unnecessary, but they have been coded with
* a possible extension to allow for the full blown MHD equations in mind...
*
*
*
* conventions used:
* -----------------
*
* stencil[3]            3 consecutive nodes along a given coordinate axis
*
* u[NHYDRO]             vector of hydrodynamic variables
*
* u_node[NHYDRO]        cell-averaged hydrodynamic variables
*                       - du_ltd() requires u[3] which is >>a given u<< along a 3-stencil
*
* u_infc[NHYDRO]        hydro-variables at cell interface
*
* flux_infc[NHYDRO]     fluxes at cell interface
*
*===============================================================================*/

/* the unit matrix has been initialized in startrun.c */
extern int unit_matrix[NDIM][NDIM];

/*===============================================================================
* hydro_flux:
* -----------
*
* returns the fluxes for the conserved variables in flux[NHYDRO]
* - the direction is specified by the non-zero element of idim[NDIM]!
*
* NOTE:  passing flux[] as an argument is the fail-safe version!
*        -> if we were to let hydro_flux return a pointer to a newly
*           allocated flux[] we would also need to make a call to free()!
*===============================================================================*/
void hydro_flux(double u[NHYDRO], double B[NDIM], double flux[NHYDRO], int idim)
{
   /* shortcut for hydro-variables */
   double a, mom[NDIM], dens, p, Bx, By, Bz, B2, B2over2, VdotB, V;
   int    i;
   
   dens   = u[Udens];
   mom[X] = u[UmomdensX];
   mom[Y] = u[UmomdensY];
   mom[Z] = u[UmomdensZ];
   
// we only placed an #ifdef below to avoid 3 floating point operations for HYDRO 
#ifdef MHD
   Bx       = B[X]; By = B[Y]; Bz = B[Z];
   B2       = Bx*Bx + By*By + Bz*Bz ;  // here, we can do it this way, because we are at the interface
   B2over2  = B2 / 2.;
   VdotB    = ( mom[X]*Bx + mom[Y]*By + mom[Z]*Bz ) / dens ;
#else // MHD
   Bx       = 0.; By = 0.; Bz = 0.;
   B2       = 0.;
   B2over2  = 0.;
   VdotB    = 0.;
#endif // MHD   
   
   /* the fluxes according to the hydro equations */
   if(dens > MACHINE_ZERO)  // TODO: think about the case of small densities (vacuum) !
     {
      p = calc_p(u);
      V = mom[idim]/dens;
      
      flux[Udens]     = mom[idim];
      flux[UmomdensX] = V*mom[X] + p*(double)unit_matrix[idim][X];
      flux[UmomdensY] = V*mom[Y] + p*(double)unit_matrix[idim][Y];
      flux[UmomdensZ] = V*mom[Z] + p*(double)unit_matrix[idim][Z];
      flux[UEdens]    = V*(u[UEdens]+p);
      flux[Untrpy]    = V*u[Untrpy];
      
#ifdef MHD
      // additional terms in the fluxes according to the supercomoving MHD equations
      flux[UmomdensX] += B2over2 * (double)unit_matrix[idim][X] - Bx * B[idim];
      flux[UmomdensY] += B2over2 * (double)unit_matrix[idim][Y] - By * B[idim];
      flux[UmomdensZ] += B2over2 * (double)unit_matrix[idim][Z] - Bz * B[idim];
      flux[UEdens]    += B2over2 * V - B[idim] * VdotB;
#endif
      
     }
   else
     {
      for(i=0; i<NADVECT; i++)
         flux[i] = (double)0.0;
     }
}


/*===============================================================================
* hydro_source:
* -------------
*
* calculate the right-hand-side of the Euler equations in conservative form,
* i.e. the Source terms 
*
* input:
*         current supercomoving time
*         the node pointer
*         the hydro-variables
*
* output:
*         the sources for each hydro-variable according to Euler equations
*
*===============================================================================*/
void hydro_source(double super_t, double timestep, nptr cur_node, double Source[NHYDRO], double BSource[NDIM], double B2)
{
   double grav[NDIM];
   double edens, H_term, a, a_dot;
   double *u;
#ifdef MHD
   double *B;
#endif
   int    ihydro;
   
   /* NOTE: the temporally averaged hydro variables are stored in u_tmp[] (cf. ave_hydro()!) */
   u = &(cur_node->u_tmp[0]);
#ifdef MHD
   B = &(cur_node->B_tmp[0]);  // new!!!
#endif
   
   if(u[Udens] > MACHINE_ZERO)
     {
#if (HYDRO_TEST==4 || HYDRO_TEST==8 || ISOLATED)
      a     = 1.0;
      a_dot = 0.0;
#else
      /* the Runge-Kutta time integration requires sources at "mid_a" */
      a     =  calc_super_a(super_t + timestep/2);
      a_dot =  calc_a_dot(a);
#endif
      
      /* NOTE: grav = -gradPhi !!!! (force.forces already contains the '-' sign) */
      grav[X] = cur_node->force.forces[X];
      grav[Y] = cur_node->force.forces[Y];
      grav[Z] = cur_node->force.forces[Z];
      
#ifdef MONOATOMIC
      Source[Udens]     = 0.;
      Source[UmomdensX] = u[Udens]*grav[X];
      Source[UmomdensY] = u[Udens]*grav[Y];
      Source[UmomdensZ] = u[Udens]*grav[Z];
      Source[UEdens]    = u[UmomdensX]*grav[X]+u[UmomdensY]*grav[Y]+u[UmomdensZ]*grav[Z];
      Source[Untrpy]    = 0.;
           
#ifdef MHD
      Source[UEdens]   += a * a_dot * B2 / 2.;  // additional source term from supercomoving MHD equations
      BSource[X]        = 0.5 * a * a_dot * B[X];
      BSource[Y]        = 0.5 * a * a_dot * B[Y];
      BSource[Z]        = 0.5 * a * a_dot * B[Z];
#endif

#else /* MONOATOMIC */                              // MHD only implemented for MONOATOMIC!
      edens             = u[Uedens];
      H_term            = -a * a_dot * (3*simu.gamma-5);
      
      Source[Udens]     = 0.;
      Source[UmomdensX] = u[Udens]*grav[X];
      Source[UmomdensY] = u[Udens]*grav[Y];
      Source[UmomdensZ] = u[Udens]*grav[Z];
      Source[UEdens]    = H_term*edens
         + u[UmomdensX]*grav[X]
         + u[UmomdensY]*grav[Y]
         + u[UmomdensZ]*grav[Z];
      Source[Untrpy]    = H_term*u[Untrpy];
#endif
     }
   else
     {
      Source[Udens]     = 0.;
      Source[UmomdensX] = 0.;
      Source[UmomdensY] = 0.;
      Source[UmomdensZ] = 0.;
      Source[UEdens]    = 0.;
      Source[Untrpy]    = 0.;
     }
   
}


/*===============================================================================
* min_max_speed:
* --------------
*
* calculate the minimum and maximum local speed at cell-interface.
*
* Note: Umomdens is a placeholder for UmomdensX, UmomdensY or UmomdensZ depending
* on which direction were are in. This function does not know the direction!
*===============================================================================*/
void min_max_speeds(double u_infc[3][2][NHYDRO], double B_infc[3][2][NDIM], double a[2][2], int Umomdens)
{
   double v[2], c[2], B2[3][2];
   int i, j;
   
   for (i=0;i<3;i++) {     
      for (j=0;j<2;j++) {     // calculate B^2 for passing on to calc_cf()
#ifdef MHD
         B2[i][j] = ( B_infc[i][j][X] * B_infc[i][j][X] + B_infc[i][j][Y] * B_infc[i][j][Y] + B_infc[i][j][Z] * B_infc[i][j][Z] );
#else
         B2[i][j] = 0.;
#endif   
      }
   }
   
   /* 1. "left" interface
      * ------------------- */
   /* u^E_i-1 = u_infc[0][1] speeds */
   if(u_infc[0][1][Udens] > MACHINE_ZERO)
     {
      v[0] = u_infc[0][1][Umomdens]/u_infc[0][1][Udens];
      c[0] = calc_cf(u_infc[0][1], B2[0][1]);
     }
   else
     {
      v[0] = 0.0;
      c[0] = 0.0;
     }
   
   /* u^W_i = u_infc[1][0] speeds */
   if(u_infc[1][0][Udens] > MACHINE_ZERO)
     {
      v[1] = u_infc[1][0][Umomdens]/u_infc[1][0][Udens];
      c[1] = calc_cf(u_infc[1][0], B2[1][0]);
     }
   else
     {
      v[1] = 0.0;
      c[1] = 0.0;
     }
   
   /* local min/max speeds */
   a[0][0]      = MIN(v[0]-c[0], v[1]-c[1]);
   a[0][1]      = MAX(v[0]+c[0], v[1]+c[1]);
   
   /* 2. "right" interface
      * -------------------- */
   /* u^E_i = u_infc[1][1] interface */
   if(u_infc[1][1][Udens] > MACHINE_ZERO)
     {
      v[0] = u_infc[1][1][Umomdens]/u_infc[1][1][Udens];
      c[0] = calc_cf(u_infc[1][1], B2[1][1]);
     }
   else
     {
      v[0] = 0.0;
      c[0] = 0.0;
     }
   
   /* u^W_i+1 = u_infc[2][0] interface */
   if(u_infc[2][0][Udens] > MACHINE_ZERO)
     {
      v[1] = u_infc[2][0][Umomdens]/u_infc[2][0][Udens];
      c[1] = calc_cf(u_infc[2][0], B2[2][0]);
     }
   else
     {
      v[1] = 0.0;
      c[1] = 0.0;
     }
   
   /* local min/max speeds */
   a[1][0]      = MIN(v[0]-c[0], v[1]-c[1]);
   a[1][1]      = MAX(v[0]+c[0], v[1]+c[1]);
   
   
   /* limit local speeds to 0 */
   a[0][0] = MIN(a[0][0], 0.0);
   a[0][1] = MAX(a[0][1], 0.0);
   a[1][0] = MIN(a[1][0], 0.0);
   a[1][1] = MAX(a[1][1], 0.0);
}

/*====================================================================================
* calc_fluxKNP_infc:
* ------------------
*
* actually calculate the numerical fluxes at the "left" and "right" cell interfaces 
*              -> this is equation (6) of Ziegler (2004)
*
* input:
*         flux_infc     -> analytical fluxes at "left" and "right" interfaces
*         u_infc        -> reconstructed hydro-variables at interfaces
*         a             -> local min/max (sound) speeds
*
* output:
*         Flux[2]       -> vector of "left" and "right" fluxes in all hydro-variables
*
*====================================================================================*/
void calc_fluxKNP_infc(double flux_infc[3][2][NHYDRO], double u_infc[3][2][NHYDRO], double a[2][2],
                       double Flux[2][NHYDRO])
{
   double da;
   int    ivar;
   
   /* F^x_i-0.5 */
   da = a[0][1]-a[0][0];
   if(fabs(da) > MACHINE_ZERO)
     {
      for(ivar=0; ivar<NADVECT; ivar++)
        {
         Flux[0][ivar] = (   a[0][1] * flux_infc[0][1][ivar] 
                             - a[0][0] * flux_infc[1][0][ivar]
                             + a[0][1] * a[0][0] * (u_infc[1][0][ivar] - u_infc[0][1][ivar])
                             ) / da;
        }
     }
   else
     {
      for(ivar=0; ivar<NADVECT; ivar++)
        {
         Flux[0][ivar] = 0.0;
        }
     }
   
   /* F^x_i+0.5 */
   da = a[1][1]-a[1][0];
   if(fabs(da) > MACHINE_ZERO)
     {
      for(ivar=0; ivar<NADVECT; ivar++)
        {
         Flux[1][ivar] = (   a[1][1] * flux_infc[1][1][ivar]
                             - a[1][0] * flux_infc[2][0][ivar]
                             + a[1][1] * a[1][0] * (u_infc[2][0][ivar] - u_infc[1][1][ivar])
                             ) / da;
        }
     }
   else
     {
      for(ivar=0; ivar<NADVECT; ivar++)
        {
         Flux[1][ivar] = 0.0;
        }
     }
}

/*===============================================================================
* calc_FluxKNP:
* -------------
*
* obtain the numerical fluxes at the "left" and "right" cell interfaces (only for hydro without MHD!)
*
* input:
*         a 3D vector which is a 5x5x5 stencil giving access to the nodes: MHDnodes[z][y][x]
*         examples of usage:
*           MHDnodes[2][2][2] = current cell i,j,k
*           MHDnodes[1][3][0] = cell i-2,j+1,k-1
*         etc.
*
* output:
*         the "left" and "right" numerical fluxes at nodes i-1, i, i+1  (for all three dimensions)
*
*
*===============================================================================*/
void calc_FluxKNP(nptr MHDnodes[5][5][5], double Flux[NDIM][2][NHYDRO])
{
   nptr     stencil[3];
   double    u_infc[3][2][NHYDRO];      // elements: [i-1,i,i+1] [E,W] [NHYDRO]
   double    B_infc[3][2][NDIM];        // elements: [i-1,i,i+1] [E,W] [x,y,z]
   double flux_infc[3][2][NHYDRO];      // elements: [i-1,i,i+1] [E,W] [NHYDRO]
   double a[2][2], da;                  // a[0][0] = a^-_(i-1/2); a[0][1] = a^+_(i-1/2); 
                                        // a[1][0] = a^-_(i+1/2); a[1][1] = a^+_(i+1/2)  -- min/max speeds
   int    Umomdens[NDIM];
   int    i,j,k;
   
   /* initialize required matrizes and vectors */
   Umomdens[X] = UmomdensX;
   Umomdens[Y] = UmomdensY;
   Umomdens[Z] = UmomdensZ;
   
   zero_B_infc(B_infc);    // even without MHD, we're passing a B_infc pointer, so it has to be initialized to zero!
   
   // go over all three spatial directions to obtain hydro-variables at cell interfaces ///////////////
   
   // x direction /////////////////////////////////////////////////////////////////////////////////////
   
   stencil[0] = MHDnodes[2][2][0]; stencil[1] = MHDnodes[2][2][1]; stencil[2] = MHDnodes[2][2][2];
   calc_u_infc(stencil, u_infc[0]);
   stencil[0] = MHDnodes[2][2][1]; stencil[1] = MHDnodes[2][2][2]; stencil[2] = MHDnodes[2][2][3];
   calc_u_infc(stencil, u_infc[1]);
   stencil[0] = MHDnodes[2][2][2]; stencil[1] = MHDnodes[2][2][3]; stencil[2] = MHDnodes[2][2][4];
   calc_u_infc(stencil, u_infc[2]);  
   
   hydro_flux(u_infc[0][1], B_infc[0][1], flux_infc[0][1], X);
   hydro_flux(u_infc[1][0], B_infc[1][0], flux_infc[1][0], X);
   hydro_flux(u_infc[1][1], B_infc[1][1], flux_infc[1][1], X);
   hydro_flux(u_infc[2][0], B_infc[2][0], flux_infc[2][0], X);
   min_max_speeds(u_infc, B_infc, a, Umomdens[X]);
   calc_fluxKNP_infc(flux_infc, u_infc, a, Flux[X]);
   
   // y direction /////////////////////////////////////////////////////////////////////////////////////
   
   stencil[0] = MHDnodes[2][0][2]; stencil[1] = MHDnodes[2][1][2]; stencil[2] = MHDnodes[2][2][2];
   calc_u_infc(stencil, u_infc[0]);
   stencil[0] = MHDnodes[2][1][2]; stencil[1] = MHDnodes[2][2][2]; stencil[2] = MHDnodes[2][3][2];
   calc_u_infc(stencil, u_infc[1]);
   stencil[0] = MHDnodes[2][2][2]; stencil[1] = MHDnodes[2][3][2]; stencil[2] = MHDnodes[2][4][2];
   calc_u_infc(stencil, u_infc[2]);  
   
   hydro_flux(u_infc[0][1], B_infc[0][1], flux_infc[0][1], Y);
   hydro_flux(u_infc[1][0], B_infc[1][0], flux_infc[1][0], Y);
   hydro_flux(u_infc[1][1], B_infc[1][1], flux_infc[1][1], Y);
   hydro_flux(u_infc[2][0], B_infc[2][0], flux_infc[2][0], Y);
   min_max_speeds(u_infc, B_infc, a, Umomdens[Y]);
   calc_fluxKNP_infc(flux_infc, u_infc, a, Flux[Y]);
   
   // z direction /////////////////////////////////////////////////////////////////////////////////////
   
   stencil[0] = MHDnodes[0][2][2]; stencil[1] = MHDnodes[1][2][2]; stencil[2] = MHDnodes[2][2][2];
   calc_u_infc(stencil, u_infc[0]);
   stencil[0] = MHDnodes[1][2][2]; stencil[1] = MHDnodes[2][2][2]; stencil[2] = MHDnodes[3][2][2];
   calc_u_infc(stencil, u_infc[1]);
   stencil[0] = MHDnodes[2][2][2]; stencil[1] = MHDnodes[3][2][2]; stencil[2] = MHDnodes[4][2][2];
   calc_u_infc(stencil, u_infc[2]);  
   
   hydro_flux(u_infc[0][1], B_infc[0][1], flux_infc[0][1], Z);
   hydro_flux(u_infc[1][0], B_infc[1][0], flux_infc[1][0], Z);
   hydro_flux(u_infc[1][1], B_infc[1][1], flux_infc[1][1], Z);
   hydro_flux(u_infc[2][0], B_infc[2][0], flux_infc[2][0], Z);
   min_max_speeds(u_infc, B_infc, a, Umomdens[Z]);
   calc_fluxKNP_infc(flux_infc, u_infc, a, Flux[Z]);
}


/*===============================================================================
 * calc_cf:
 * ------------
 *
 * For the standard hydrodynamic case this function simply calculates the sound speed:
 *               c_f = c_s = sqrt(gamma*p/dens)
 * 
 * For the MHD case the function instead returns an upper bound estimate for the fast magnetosonic speed:
 *               c_f = sqrt ( c_s^2 + c_a^2 )
 * where c_a is the Alfvén speed without the negative Bi term
 *               c_a = sqrt ( B^2 / density )      (mu = 1 as everywhere)
 *===============================================================================*/
double calc_cf(double u[NHYDRO], double B2)
{
   double p, cs2, ca2;
   
   /* make sure that the gas density u[Udens] is not ZERO >>before<< calling calc_cf()! */
   
   p   = calc_p(u);
   cs2 = simu.gamma * p / u[Udens];      // sound speed
   
#ifdef MHD
   ca2 = B2 / u[Udens];         // Alfvén speed   // TODO: what if the density is small???
   return sqrt ( cs2 + ca2 );    // return upper bound estimate for fast magnetosonic speed
#else
   return sqrt(cs2);             // just return the sound speed
#endif
}


/*============================================================================================
 * adjust_globalCFLspeed():
 * ------------------------
 *
 * at a given node calculate MAX(|u_x|+cs, |u_y|+cs, |u_z|+cs) 
 * and adjust global.cfl_speed accordingly
 *
 *==========================================================================================*/
void adjust_globalCFLspeed(double u[NHYDRO], double B2)
{
   double cf, v[NDIM];
   
   if(u[Udens] > MACHINE_ZERO)
     {
      cf   = calc_cf(u, B2);   // make sure that B2 = 0 ifndef MHD
        
      v[X] = fabs(u[UmomdensX]/u[Udens]) + cf;
      v[Y] = fabs(u[UmomdensY]/u[Udens]) + cf;
      v[Z] = fabs(u[UmomdensZ]/u[Udens]) + cf;
      
      global.cfl_speed = MAX(global.cfl_speed, v[X]);
      global.cfl_speed = MAX(global.cfl_speed, v[Y]);
      global.cfl_speed = MAX(global.cfl_speed, v[Z]);
     }
}


/*============================================================================================
 * du_ltd:
 * -------
 *
 * calculate a (limited) du required to reconstruct hydro-variables 
 * at cell interfaces
 *
 * input:
 *         a vector of variables along a 3-stencil i-1, i, i+1
 *
 * output:
 *         du as given by some limitation formula
 *
 * NOTE:   - the stencil can be along any coordinate axis
 *         - u_node[0,1,2] contains a certain hydro-variable along that stencil!
 *         - du = du/dx * dx/2 = u'/2 * dx => return (Delta_u) / 2
 *==========================================================================================*/
double du_ltd(double u_node[3])
{
   double du, slope_left, slope_right;
   
   slope_left  = u_node[1]-u_node[0];
   slope_right = u_node[2]-u_node[1];
   
#ifdef SLOPE_vanLeer
   du = vanLeer(slope_left, slope_right);
#endif
   
#ifdef SLOPE_MINMOD
   du = minmod(slope_left, slope_right);
#endif
   
#ifdef SLOPE_SUPERBEE
   du = superbee(slope_left, slope_right);
#endif
   
#ifdef FLAT_RECONSTRUCTION
   du = 0.0;
#endif
   
   return (du);
}

/*===============================================================================
 * update_edens:
 * -------------
 *
 * use DUAL_ENERGY formalism to update u[Uedens]:
 *
 *    - we use values stored in u[] with the DUAL_ENERGY criterion
 *      but temporarily store the new edens value in u_tmp[]
 *
 *    - DO NOT TOUCH ANYTHING ELSE IN u_tmp[] AS THESE VALUES WILL BE USED ELSEWHERE!
 *
 *===============================================================================*/
void update_edens(gridls *cur_grid)
{
   partptr       cur_part;
   pqptr         cur_pquad;
   cqptr         cur_cquad, icur_cquad;
   nqptr         cur_nquad, icur_nquad;
   nptr          cur_node;
   nptr          Tnode, Bnode, Nnode, Snode, Wnode, Enode;
   nptr          tsc_nodes[3][3][3];
   
   long          z, y, x, iVar;
      
   double        edens, Edens, Tdens, Bdens, eta2edens, gm1, cs, v[NDIM], grad_edens[NDIM];
   double       *u, *u_tmp, *Tu, *Bu, *Nu, *Su, *Wu, *Eu;
   int           use_S_system;
   
   long          icquad;
   
   gm1 = simu.gamma-1.;

   for(cur_pquad=cur_grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next)
     {
      z = cur_pquad->z;
#ifdef WITH_OPENMPM
#pragma omp parallel shared(cur_grid, cur_pquad, gm1) private(icquad, cur_cquad, icur_cquad, cur_nquad, icur_nquad, cur_node, x, y, z, tsc_nodes, Tnode, Bnode, Wnode, Enode, Snode, Nnode, u, u_tmp, use_S_system, edens, Tdens, Edens, Bdens, eta2edens, cs, v, grad_edens, Tu, Bu, Nu, Su, Wu, Eu)
#pragma omp for schedule(static)
      for(icquad=0; icquad<cur_grid->l1dim; icquad++)
        {
         z         = cur_pquad->z   + icquad;
         cur_cquad = cur_pquad->loc + icquad;          
#else
      for(cur_cquad = cur_pquad->loc; cur_cquad < cur_pquad->loc + cur_pquad->length; cur_cquad++, z++)  
        {  
#endif        
            for(icur_cquad  = cur_cquad; icur_cquad != NULL; icur_cquad  = icur_cquad->next)
           {
            y = icur_cquad->y;
            for(cur_nquad = icur_cquad->loc; cur_nquad < icur_cquad->loc + icur_cquad->length; cur_nquad++, y++) 
              { 
               for(icur_nquad  = cur_nquad; icur_nquad != NULL; icur_nquad  = icur_nquad->next)
                 {
                  x = icur_nquad->x;
                  for(cur_node = icur_nquad->loc; cur_node < icur_nquad->loc + icur_nquad->length; cur_node++, x++)
                    {
                     /* get pointers to all 26 surrounding neighbours */
                     tsc_nodes[1][1][1] = cur_node;
                     get_TSCnodes(cur_grid, cur_pquad, icur_cquad, icur_nquad, tsc_nodes, &z, &y, &x); 
                     
                     /* do not use boundary nodes */
                     if(test_tsc(tsc_nodes) == TRUE)
                       {
                        /* ...our scheme does NOT require to synchronize cur_node->u_tmp[] !!!
                        * however, we can only use u_tmp[Uedens] as a temporary variable as the
                        * other values in u_tmp[] are needed by other routines! */
                        u     = cur_node->u;
                        u_tmp = cur_node->u_tmp;

                        /* default = E-system */
                        use_S_system = FALSE;
                        
                        /* avoid division by zero */
                        if(u[Udens] > MACHINE_ZERO)
                          {
                           Tdens  = 0.5 * ( pow2(u[UmomdensX])+pow2(u[UmomdensY])+pow2(u[UmomdensZ]) ) / u[Udens];
                           Edens  = u[UEdens];
#ifdef MHD
                           // calculate Bdens = 0.5*B^2 (magnetic energy density).
                           // idea: we get the cell-averaged values by averaging over pairs of opposing interfaces.
                           // remember that the B components are saved at the interfaces!
                             
                           Bdens = 0.125 * 
                             (pow2(cur_node->B[X] + tsc_nodes[1][1][2]->B[X]) +
                              pow2(cur_node->B[Y] + tsc_nodes[1][2][1]->B[Y]) +
                              pow2(cur_node->B[Z] + tsc_nodes[2][1][1]->B[Z]) );
#else // MHD
                           Bdens = 0.0;
#endif // MHD
                           
                           
                           /* edens as given by the E-System */
                           edens  = MAX (Edens - Tdens - Bdens, MACHINE_ZERO);
#ifdef MHD_DEBUG
                           if (Bdens/Edens > 0.999999 )  // this can't be since the magnetic energy would be equal to or bigger than the total energy...
                           {
                             fprintf(io.logfile,"update_edens(): numerical error: magnetic energy >= total energy!\n");
                             fprintf(io.logfile,"   Bdens/Edens = %16.8g\n   Edens=%16.8g\n   Bdens=%16.8g\n   Tdens=%16.8g\n   edens=%16.8g\n",
                                     Bdens/Edens, Edens, Bdens, Tdens, edens);
                             fprintf(io.logfile,"coordinates: x=%i,  y=%i,  z=%i\n",x,y,z);
                             exit(1);
                           }
#endif // MHD_DEBUG                           
                           
#ifdef DUAL_ENERGY         /* should we rather use the S-System? */
#ifdef RYU_CRITERION
                           //pointers to all six direct neighbours
                           Tnode = tsc_nodes[2][1][1]; Tu = (double *)Tnode->u;
                           Bnode = tsc_nodes[0][1][1]; Bu = (double *)Bnode->u;
                           
                           Nnode = tsc_nodes[1][2][1]; Nu = (double *)Nnode->u;
                           Snode = tsc_nodes[1][0][1]; Su = (double *)Snode->u;
                           
                           Enode = tsc_nodes[1][1][2]; Eu = (double *)Enode->u;
                           Wnode = tsc_nodes[1][1][0]; Wu = (double *)Wnode->u;
                             
                           // little helper to avoid multiplicating three times the same thing in the if-condition below
                           eta2edens = ETA_DUAL_ENERGY2 * fabs (u[Uedens]) ;
                           
                           // 1. shock criterion as in Ryu et al., 1993,  ApJ 414, 1
                           //                          Gheller et al. 1998, MNRAS 295, 519
                           //    (the if-condition is TRUE if there is no shock present in X Y or Z direction => entropy is conserved)
                           grad_edens[X] = Eu[Uedens]-Wu[Uedens];
                           grad_edens[Y] = Tu[Uedens]-Bu[Uedens];
                           grad_edens[Z] = Nu[Uedens]-Su[Uedens];
                           
                           if( (pow2(grad_edens[X])+pow2(grad_edens[Y])+pow2(grad_edens[Z])) < pow2(eta2edens) )
                              use_S_system = TRUE;
                           
//                           if ( fabs(Tu[Uedens]-Bu[Uedens]) < eta2edens &&
//                                fabs(Nu[Uedens]-Su[Uedens]) < eta2edens &&
//                                fabs(Eu[Uedens]-Wu[Uedens]) < eta2edens )
//                              use_S_system = TRUE;                           
#endif /* RYU_CRITERION */
                           /* 2. dual energy formalism according to Feng, Shu & Zhang, 2004, ApJ 612, 1 */        
                           if ( (Tdens+Bdens) /Edens > 1. - ETA_DUAL_ENERGY1)
                              use_S_system = TRUE;
                           
                           if(use_S_system)
                             {
                              /* use entropy to calculate internal energy instead of Edens-Tdens */
                              edens = u[Untrpy] * pow(u[Udens], gm1) / gm1;
                             }
#endif /* DUAL_ENERGY */
                          }
                        else
                          {
                           edens = 0.0;
                          }
                        
                        /* ...store edens temporarily in u_tmp[Uedens] */
                        u_tmp[Uedens] = edens;
                       }
#ifdef VERBOSE
                     else
                       {
                        fprintf(stderr,"update_edens: x=%ld y=%ld z=%ld we should not see this message for regular grids!\n",x,y,z);
                       }
#endif
                    }
                 }
              }
           }
        }
     }
}


/*===============================================================================
 * getCFLspeed_and_syncES:
 * ------------------------
 *
 * 1. copy u_tmp[Uedens] over to u[Uedens] and simultaneously
 *
 * 2. synchronize u[UEdens] and u[Untrpy]
 *
 *===============================================================================*/
void getCFLspeed_and_syncES(gridls *cur_grid)
{
   partptr       cur_part;
   pqptr         cur_pquad;
   cqptr         cur_cquad, icur_cquad;
   nqptr         cur_nquad, icur_nquad;
   nptr          cur_node;
   nptr          tsc_nodes[3][3][3];
   
   double        *u;
   double        edens, Edens, Tdens, Sdens, gm1, B2, Bdens;
   long          z, y, x, iVar;
   
   long          icquad;

   gm1 = simu.gamma-1.;
   
   /* update_edens() updates global.cfl_speed according to the latest hydro-variables */
   global.cfl_speed = 0.0;
   
   for(cur_pquad=cur_grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next)
     {
      z = cur_pquad->z;
#ifdef WITH_OPENMPM
#pragma omp parallel shared(cur_grid, cur_pquad, gm1) private(icquad, cur_cquad, icur_cquad, cur_nquad, icur_nquad, cur_node, tsc_nodes, u, edens, Edens, Tdens, Sdens, B2, Bdens, x, y, z, iVar)
#pragma omp for schedule(static)
      for(icquad=0; icquad<cur_grid->l1dim; icquad++)
        {
         z         = cur_pquad->z   + icquad;
         cur_cquad = cur_pquad->loc + icquad;
#else
      for(cur_cquad = cur_pquad->loc; cur_cquad < cur_pquad->loc + cur_pquad->length; cur_cquad++, z++)  
        {  
#endif
            for(icur_cquad  = cur_cquad; icur_cquad != NULL; icur_cquad  = icur_cquad->next)
           {
            y = icur_cquad->y;
            for(cur_nquad = icur_cquad->loc; cur_nquad < icur_cquad->loc + icur_cquad->length; cur_nquad++, y++) 
              { 
               for(icur_nquad  = cur_nquad; icur_nquad != NULL; icur_nquad  = icur_nquad->next)
                 {
                  x = icur_nquad->x;
                  for(cur_node = icur_nquad->loc; cur_node < icur_nquad->loc + icur_nquad->length; cur_node++, x++)
                    {
                       
                     /* retrieve edens from u_tmp[] calculated in update_edens() */
                     edens = cur_node->u_tmp[Uedens];
                     
                     Tdens = 0.5* ( pow2(cur_node->u[UmomdensX]) +
                                    pow2(cur_node->u[UmomdensY]) +
                                    pow2(cur_node->u[UmomdensZ])) / cur_node->u[Udens];
                       
#ifdef MHD
                     tsc_nodes[1][1][1] = cur_node;
                     get_TSCnodes(cur_grid, cur_pquad, icur_cquad, icur_nquad, tsc_nodes, &z, &y, &x);
                       
                     // calculate 0.5*a*B^2 (cell-averaged magnetic energy density).
                     // idea: we get the cell-averaged values by averaging over pairs of opposing interfaces.
                     // remember that the B components are saved at the interfaces!
                       
                     if(test_tsc(tsc_nodes) == TRUE)  
                        B2 = 0.25 * (pow2(cur_node->B[X] + tsc_nodes[1][1][2]->B[X]) +
                                     pow2(cur_node->B[Y] + tsc_nodes[1][2][1]->B[Y]) +
                                     pow2(cur_node->B[Z] + tsc_nodes[2][1][1]->B[Z]) );
                     else
                        B2 = pow2(cur_node->B[X]) + pow2(cur_node->B[Y]) + pow2(cur_node->B[Z]);
                       
                     Bdens = 0.5 * B2;
#else // MHD
                     B2    = 0.0;
                     Bdens = 0.0;
#endif // MHD
                     
                     cur_node->u[Uedens] = edens;
                     cur_node->u[UEdens] = edens + Tdens + Bdens;
                     cur_node->u[Untrpy] = edens * gm1 / pow(cur_node->u[Udens],gm1);
                  
                     /* update global.cfl_speed by calculating maximal velocities at cur_node 
                      * using the most recent u[] values */ 
                     adjust_globalCFLspeed(cur_node->u, B2);
                     
                    }
                 }
              }
           }
        }
     }
   
}


/*===============================================================================
 * calc_e:
 * -------
 *
 * calculate internal energy for a given temperature T (in [K])
 *===============================================================================*/
double calc_e(double T)
{
   double e, molecular_weight;
   
   e = kB_per_mp * T;
   //e = T;
   
#ifdef MONOATOMIC
   /* 1/(gamma-1) = 1.5 as gamma = 5/3 */
   e *= 1.5;
#else
   e *= (1.0 / (simu.gamma-1));
   
   if(T > 1.0e4)	
      /* assuming FULL ionization */
      molecular_weight = 4 / (8 - 5 * (1 - simu.H_frac));
   else			
      /* assuming NEUTRAL gas */
      molecular_weight = 4 / (1 + 3 * simu.H_frac);
   
   e /= molecular_weight;
#endif
   
   return( e );
}

/*===============================================================================
* calc_p:
* -------
*
* calculate thermal(!!) pressure p = (gamma-1) * edens
* (the thermal pressure is the same for the MHD case)
*===============================================================================*/
double calc_p(double u[NHYDRO])
{
   double p;
   
#ifdef MONOATOMIC  
   
   p = MAX( (simu.gamma-1.) * u[Uedens] , MACHINE_ZERO );
   
#else /* MONOATOMIC */
   
   /* check formulae in Anninos & Norman (1994) */
   fprintf(stderr,"\n\ncalc_p: NON-MONOATOMIC EQUATION OF STATE NOT YET IMPLEMENTED...\n");
   exit(0);
   
#endif /* MONOATOMIC */
   
   return( p );
}

/*===============================================================================
 * calc_T:
 * -------
 *
 * calculate temperature T (in [K]) for a given internal energy
 *===============================================================================*/
double calc_T(double e)
{
   double T, molecular_weight;
   
   T = e/kB_per_mp;
   //T = e;
   
#ifdef MONOATOMIC
   /* 1/(gamma-1) = 1.5 as gamma = 5/3 */
   T /= 1.5;
#else
   T /= (1.0 / (simu.gamma-1));
   
   if(T > 1.0e4)	
      /* assuming FULL ionization */
      molecular_weight = 4 / (8 - 5 * (1 - simu.H_frac));
   else			
      /* assuming NEUTRAL gas */
      molecular_weight = 4 / (1 + 3 * simu.H_frac);
   
   T *= molecular_weight;
#endif
   
   return( T );
}


/*===========================================================================================
* calc_u_infc:
* ------------
*
* reconstruct hydro-variables at "left" and "right" cell interface of centre node = stencil[1]
*
* input:
*         a 3-stencil
*
* output:
*         values of hydro-variables at cell interfaces of centre node
*         along the coordinate axis determined by the stencil
*
* NOTE:   - the returned values are at the W/E, N/S, or T/B interface
*           depending on the orientation of the stencil!
*           the routine does not make any assumption about the spatial direction!
*         - however, stencil[0] should be the "backward", stencil[1] the "actual",
*           and stencil[2] the "forward" node...
*===========================================================================================*/
void calc_u_infc(nptr stencil[3], double u_infc[2][NHYDRO])
{
   int    inode, ivar;
   double u_node[NHYDRO][3], du;
   
   /* transfer cell-averaged hydro-variables to u_node-vector (incl. edens!) */
   for(inode=0; inode<3; inode++)
      for(ivar=0; ivar<NHYDRO; ivar++)
         u_node[ivar][inode] = (double)stencil[inode]->u[ivar];
   
#ifdef RECONSTRUCT_PRESSURE
   for(inode=0; inode<3; inode++)
      u_node[Uedens][inode] *= (simu.gamma-1.);
#endif
   
#ifdef RECONSTRUCT_VELOCITIES
   /* reconstruct velocities rather than momentum density */
   for(inode=0; inode<3; inode++)
     {
      if(u_node[Udens][inode] > MACHINE_ZERO)
        {
         u_node[UmomdensX][inode] = u_node[UmomdensX][inode]/u_node[Udens][inode];
         u_node[UmomdensY][inode] = u_node[UmomdensY][inode]/u_node[Udens][inode];
         u_node[UmomdensZ][inode] = u_node[UmomdensZ][inode]/u_node[Udens][inode];
        }
      else
        {
         u_node[UmomdensX][inode] = 0.0;
         u_node[UmomdensY][inode] = 0.0;
         u_node[UmomdensZ][inode] = 0.0;
        }
     }
#endif
   
   /* get u_infc[NHYDRO] at "left" and "right" cell interface */
   for(ivar=0; ivar<NHYDRO; ivar++)
     {
      /* du = du/dx * dx/2 */
      du = du_ltd(u_node[ivar]) / 2;
      
      /*   u_+-        =       u        +- du */
      u_infc[0][ivar]  = u_node[ivar][1] - du;
      u_infc[1][ivar]  = u_node[ivar][1] + du;
     }
   
#ifdef RECONSTRUCT_PRESSURE
   for(inode=0; inode<3; inode++)
      u_node[Uedens][inode] /= (simu.gamma-1.);
#endif

#ifdef RECONSTRUCT_VELOCITIES
   /* transfer reconstructed velocities back to momentum density */
   u_infc[0][UmomdensX] = u_infc[0][UmomdensX]*u_infc[0][Udens];
   u_infc[0][UmomdensY] = u_infc[0][UmomdensY]*u_infc[0][Udens];
   u_infc[0][UmomdensZ] = u_infc[0][UmomdensZ]*u_infc[0][Udens];
   u_infc[1][UmomdensX] = u_infc[1][UmomdensX]*u_infc[1][Udens];
   u_infc[1][UmomdensY] = u_infc[1][UmomdensY]*u_infc[1][Udens];
   u_infc[1][UmomdensZ] = u_infc[1][UmomdensZ]*u_infc[1][Udens];
#endif
   
#ifdef CHECK_RECONSTRUCTION  
   {
      double edens;
      int    i, iflat;
      
      /* are the reconstructed hydro-variables physically meaningful */
      iflat = 0;
      for(i=0; i<2; i++)
        {
         if(u_infc[i][Udens] < MACHINE_ZERO)
           {
            fprintf(stderr,"calc_u_infc: negative dens -> using flat reconstruction\n");
            iflat = 1;
           }
         
         if(u_infc[i][UEdens] < MACHINE_ZERO)
           {
            fprintf(stderr,"calc_u_infc: negative Edens -> using flat reconstruction\n");
            iflat = 1;
           }
         
         if(u_infc[i][Uedens] < MACHINE_ZERO)
           {
            fprintf(stderr,"calc_u_infc: negative edens -> using flat reconstruction\n");
            iflat = 1;
           }
         
        }
      
      /* get u_infc[NHYDRO] at "left" and "right" cell interface by "flat" reconstruction */
      if(iflat)
        {
         for(ivar=0; ivar<NHYDRO; ivar++)
           {
            /*   u_+-        =       u       */
            u_infc[0][ivar]  = u_node[ivar][1];
            u_infc[1][ivar]  = u_node[ivar][1];
           }
        }
   }
#endif
}

/*==============================================================================
 * obtain numerical fluxes according to hydro-scheme using cur_node->u[] and
 * simply store them at cur_node->F[]
 *==============================================================================*/
void calc_Flux(gridls *cur_grid)
{   
   pqptr   cur_pquad;
   cqptr   cur_cquad, icur_cquad;
   nqptr   cur_nquad, icur_nquad;
   nptr    cur_node;
   long    x, y, z;
   
   // gain access to all neighbours and next-to-neighbours in 3D (not yet completely functional!)
   nptr    MHDnodes[5][5][5];
   
   /* numerical fluxes at all six cell-interfaces */
   double Flux[NDIM][2][NHYDRO];
   double FluxX, FluxY, FluxZ;
   double FluxX_old, FluxY_old, FluxZ_old;
   double u[NHYDRO], F[NHYDRO];
   double dBdt[NDIM];
   
   int    idim, iVar, i,j,k;
   
   long          icquad;
   double        cur_shift;  // only for DEBUGMHD
   
   
   cur_shift     = 0.5/(double)cur_grid->l1dim; // only for DEBUGMHD

   
   // set all elements of MHDnodes to null-pointers to avoid un-initialized pointers...
   for (i=0;i<5;i++) {
      for (j=0;j<5;j++) {
         for (k=0;k<5;k++) {
            MHDnodes[i][j][k] = NULL; }}}
   
   /*----------------------------------------------------------------------
    * obtain Fluxes and store updated hydro-variables in u_tmp
    *----------------------------------------------------------------------*/
   for(cur_pquad=cur_grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next)
     {
#ifdef WITH_OPENMPM
#pragma omp parallel firstprivate(MHDnodes) shared(cur_grid, cur_pquad, cur_shift) private(icquad, cur_cquad, icur_cquad, cur_nquad, icur_nquad, cur_node, x, y, z, Flux, FluxX, FluxY, FluxZ, FluxX_old, FluxY_old, FluxZ_old, u, F, dBdt, idim, iVar, i, j, k)
#pragma omp for schedule(static)
      for(icquad=0; icquad<cur_grid->l1dim; icquad++)
        {
         z         = cur_pquad->z   + icquad;
         cur_cquad = cur_pquad->loc + icquad;
#else
#if (HYDRO_TEST==0 || HYDRO_TEST==1)
         z         = cur_grid->l1dim/2;
         cur_cquad = cur_pquad->loc + z;
#else
         z = cur_pquad->z;
         
         for(cur_cquad = cur_pquad->loc;
             cur_cquad < cur_pquad->loc + cur_pquad->length; 
             cur_cquad++, z++)  
#endif
           {  
#endif /* WITH_OPENMPM */
            for(icur_cquad  = cur_cquad; 
             icur_cquad != NULL; 
             icur_cquad  = icur_cquad->next)
           {
#if (HYDRO_TEST==0 || HYDRO_TEST==2)
            y         = cur_grid->l1dim/2;
            cur_nquad = icur_cquad->loc + y;
#else
            y = icur_cquad->y;
            for(cur_nquad = icur_cquad->loc;  
                cur_nquad < icur_cquad->loc + icur_cquad->length; 
                cur_nquad++, y++) 
#endif
              { 
               for(icur_nquad  = cur_nquad; 
                   icur_nquad != NULL; 
                   icur_nquad  = icur_nquad->next)
                 {
#if (HYDRO_TEST==1 || HYDRO_TEST==2)
                  x        = cur_grid->l1dim/2;
                  cur_node = icur_nquad->loc + x;
#else
                  x = icur_nquad->x;
                  for(cur_node = icur_nquad->loc; 
                      cur_node < icur_nquad->loc + icur_nquad->length; 
                      cur_node++, x++)
#endif
                    {
                     /* get pointer to all 125 neighboring nodes */
                     get_MHDnodes(cur_grid, cur_pquad, z, y, x, MHDnodes);
                     
                     /* are we at the (MHD/hydro) boundary? */
                     if(test_mhd(MHDnodes) == TRUE)
                       {

                           /* obtain numerical fluxes according to hydro-scheme */
#ifdef MHD
                           calc_FluxKNPCT (MHDnodes, Flux, dBdt);    // ...and the magnetic field change according to CT
#else
                           calc_FluxKNP   (MHDnodes, Flux);          // ...only the hydro fluxes
#endif
                           /* store fluxes (no flux for edens!) */
                           for(iVar=0; iVar<NADVECT; iVar++)
                             {
#if (HYDRO_TEST<3)
                              /* only advect in [HYDRO_TEST]-direction */
                              F[iVar] = -(Flux[HYDRO_TEST][1][iVar]-Flux[HYDRO_TEST][0][iVar]);
#else
                              /* this is the standard! */
                              FluxX = Flux[X][1][iVar]-Flux[X][0][iVar];
                              FluxY = Flux[Y][1][iVar]-Flux[Y][0][iVar];
                              FluxZ = Flux[Z][1][iVar]-Flux[Z][0][iVar];
                              
                              F[iVar] = -(FluxX + FluxY + FluxZ);     
                                
#endif /* HYDRO_TEST<3 */
                              
                              /* store the 3D flux into/out of this cell */
                              cur_node->F[iVar] = F[iVar];
                              
#ifdef HYDRO_STORE_FLUX
                              /* store all the fluxes at all interfaces */
                              cur_node->Flux[X][0][iVar] = Flux[X][0][iVar];
                              cur_node->Flux[X][1][iVar] = Flux[X][1][iVar];
                              cur_node->Flux[Y][0][iVar] = Flux[Y][0][iVar];
                              cur_node->Flux[Y][1][iVar] = Flux[Y][1][iVar];
                              cur_node->Flux[Z][0][iVar] = Flux[Z][0][iVar];
                              cur_node->Flux[Z][1][iVar] = Flux[Z][1][iVar];
#endif
                             }
                             
#ifdef MHD
                           // store the time change of B[X] at the cell
                           for(idim=0; idim<NDIM; idim++)
                           {
                              cur_node->E[idim] = dBdt[idim];
                           }
#endif // MHD
                           
                       } /* test_mhd() != NULL */
                     else
                       {
                        /* store ZERO for the fluxes */
                        for(iVar=0; iVar<NADVECT; iVar++)
                          {
                           cur_node->F[iVar] = 0.0;
#ifdef HYDRO_STORE_FLUX
                           /* store ZERO for all the fluxes at all interfaces */
                           cur_node->Flux[X][0][iVar] = 0.0;
                           cur_node->Flux[X][1][iVar] = 0.0;
                           cur_node->Flux[Y][0][iVar] = 0.0;
                           cur_node->Flux[Y][1][iVar] = 0.0;
                           cur_node->Flux[Z][0][iVar] = 0.0;
                           cur_node->Flux[Z][1][iVar] = 0.0;
#endif
                          }
                        
#ifdef MHD
                        // store ZERO for the time change of B[X] at the cell
                        for(idim=0; idim<NDIM; idim++) {
                           cur_node->E[idim] = 0.0;
                        } 
                           
#endif
                        
                       }
                    }
                 }
              }
           }
        }
     }
}


/*==============================================================================
 * advect hydro variables:
 *
 * u_new = u + dr/dx * F
 *
 * cur_node->u[]     will contain u_new
 * cur_node->u_tmp[] will contain u     (backup copy of the imput u!)
 *
 *
 * cur_node->F[] should have been calculated by a previous call to calc_Flux()
 *
 *==============================================================================*/
void advect_hydro(gridls *cur_grid, double timestep)
{   
   pqptr   cur_pquad;
   cqptr   cur_cquad, icur_cquad;
   nqptr   cur_nquad, icur_nquad;
   nptr    cur_node;
   long    x, y, z;
   
   /* gain access to nodes i-2, i-1, i, i+1, i+2 in all three dimensions */
   nptr    tsc_nodes[3][3][3];
   
   int     idim, iVar;
   double  dtdx, u_tmp;
   
   long    icquad;
   
   
   dtdx = timestep/cur_grid->spacing;
   
   
   /*----------------------------------------------------------------------
    * advect hydro variables using pre-calculated fluxes
    *----------------------------------------------------------------------*/
   for(cur_pquad=cur_grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next)
     {
#ifdef WITH_OPENMPM
#pragma omp parallel shared(cur_grid, cur_pquad, dtdx) private(icquad, cur_cquad, icur_cquad, cur_nquad, icur_nquad, cur_node, x, y, z, tsc_nodes, idim, iVar, u_tmp)
#pragma omp for schedule(static)
      for(icquad=0; icquad<cur_grid->l1dim; icquad++)
        {
         z         = cur_pquad->z   + icquad;
         cur_cquad = cur_pquad->loc + icquad;
#else /* WITH_OPENMPM */
#if (HYDRO_TEST==0 || HYDRO_TEST==1)
         z         = cur_grid->l1dim/2;
         cur_cquad = cur_pquad->loc + z;
#else
         z = cur_pquad->z;
         for(cur_cquad = cur_pquad->loc;
             cur_cquad < cur_pquad->loc + cur_pquad->length; 
             cur_cquad++, z++)  
#endif
           {  
#endif /* WITH_OPENMPM */
            for(icur_cquad  = cur_cquad; 
             icur_cquad != NULL; 
             icur_cquad  = icur_cquad->next)
           {
#if (HYDRO_TEST==0 || HYDRO_TEST==2)
            y         = cur_grid->l1dim/2;
            cur_nquad = icur_cquad->loc + y;
#else
            y = icur_cquad->y;
            for(cur_nquad = icur_cquad->loc;  
                cur_nquad < icur_cquad->loc + icur_cquad->length; 
                cur_nquad++, y++) 
#endif
              { 
               for(icur_nquad  = cur_nquad; 
                   icur_nquad != NULL; 
                   icur_nquad  = icur_nquad->next)
                 {
#if (HYDRO_TEST==1 || HYDRO_TEST==2)
                  x        = cur_grid->l1dim/2;
                  cur_node = icur_nquad->loc + x;
#else
                  x = icur_nquad->x;
                  for(cur_node = icur_nquad->loc; 
                      cur_node < icur_nquad->loc + icur_nquad->length; 
                      cur_node++, x++)
#endif
                    {                     
                     /* update hydro-variables -> but store temporarily! */
                     for(iVar=0; iVar<NADVECT; iVar++)
                       {
                        cur_node->u_tmp[iVar] = cur_node->u[iVar];
                        cur_node->u[iVar]     = cur_node->u[iVar] + dtdx*cur_node->F[iVar];
                       }
                                          
#ifdef CHECK_ADVECTION
                     if(cur_node->u[Udens] < 0.)
                        fprintf(stderr,"advect_hydro: negative dens\n");
                     if(cur_node->u[UEdens] < 0.)
                        fprintf(stderr,"advect_hydro: negative Edens\n");
                     if(cur_node->u[Uedens] < 0.)
                        fprintf(stderr,"advect_hydro: negative edens\n");
#endif

#ifdef MHD
                     for(idim=0; idim<NDIM; idim++)   // time integration of B field is completely analogous to hydro
                       {                     
                        cur_node->B_tmp[idim] = cur_node->B[idim];
                        cur_node->B[idim]     = cur_node->B[idim] + dtdx*cur_node->E[idim];
                       }
#endif
                    }
                 }
              }
           }
        }
     }

      /* refresh u[Uedens] using newly calculated hydro-variables */
   update_edens(cur_grid);
      
   /* 1. copy u_tmp[Uedens] over to u[Uedens] and simultaneously synchronize u[UEdens] and u[Untrpy]
    * 2. this loop also calculates the new global.cfl_speed using the latest u[] values! */
   getCFLspeed_and_syncES(cur_grid);  

}

/*==============================================================================
 * advect hydro variables:
 *
 * u_new = u + 0.5*dt/dx * F + dt * S
 *
 * cur_node->u[]     will contain u_new
 * cur_node->u_tmp[] will contain u     (backup copy of the imput u!)
 *
 *
 * cur_node->F[] should have been calculated by a previous call to calc_Flux()
 *
 *
 * NOTE: we are using u_tmp[] !!!!
 *==============================================================================*/
void advect_gravity_hydro(gridls *cur_grid, double timecounter, double timestep)
{   
   pqptr   cur_pquad;
   cqptr   cur_cquad, icur_cquad;
   nqptr   cur_nquad, icur_nquad;
   nptr    cur_node;
   long    x, y, z;
   
   /* gain access to nodes i-2, i-1, i, i+1, i+2 in all three dimensions */
   nptr    mhd_nodes[5][5][5];
   
   int     iVar, idim;
   double  dtdx, dtdx2;
   
   long    icquad;
   
   double  Source[NHYDRO];
   double  BSource[NDIM];
   double  B2;
   
   dtdx             = timestep/cur_grid->spacing;
   dtdx2            = dtdx/2;
   
   /*------------------------------------------------------------------------------
    * advect hydro variables using pre-calculated fluxes and include source terms 
    *------------------------------------------------------------------------------*/
   for(cur_pquad=cur_grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next)
     {
#ifdef WITH_OPENMPM
#pragma omp parallel shared(cur_grid, cur_pquad, dtdx, dtdx2) private(icquad, cur_cquad, icur_cquad, cur_nquad, icur_nquad, cur_node, x, y, z, mhd_nodes, iVar, idim, Source, BSource, B2)
#pragma omp for schedule(static)
      for(icquad=0; icquad<cur_grid->l1dim; icquad++)
        {
         z         = cur_pquad->z   + icquad;
         cur_cquad = cur_pquad->loc + icquad;
#else /* WITH_OPENMPM */
#if (HYDRO_TEST==0 || HYDRO_TEST==1)
         z         = cur_grid->l1dim/2;
         cur_cquad = cur_pquad->loc + z;
#else
         z = cur_pquad->z;
         for(cur_cquad = cur_pquad->loc;
             cur_cquad < cur_pquad->loc + cur_pquad->length; 
             cur_cquad++, z++)  
#endif
           {  
#endif /* WITH_OPENMPM */
            for(icur_cquad  = cur_cquad; 
             icur_cquad != NULL; 
             icur_cquad  = icur_cquad->next)
           {
#if (HYDRO_TEST==0 || HYDRO_TEST==2)
            y         = cur_grid->l1dim/2;
            cur_nquad = icur_cquad->loc + y;
#else
            y = icur_cquad->y;
            for(cur_nquad = icur_cquad->loc;  
                cur_nquad < icur_cquad->loc + icur_cquad->length; 
                cur_nquad++, y++) 
#endif
              { 
               for(icur_nquad  = cur_nquad; 
                   icur_nquad != NULL; 
                   icur_nquad  = icur_nquad->next)
                 {
#if (HYDRO_TEST==1 || HYDRO_TEST==2)
                  x        = cur_grid->l1dim/2;
                  cur_node = icur_nquad->loc + x;
#else
                  x = icur_nquad->x;
                  for(cur_node = icur_nquad->loc; 
                      cur_node < icur_nquad->loc + icur_nquad->length; 
                      cur_node++, x++)
#endif
                    {
                     /* get pointer to all 125 neighboring nodes */
                     get_MHDnodes(cur_grid, cur_pquad, z, y, x, mhd_nodes);
                     
                     /* check for direct boundary node */
                     if(test_mhd(mhd_nodes) == TRUE)
                       {
       
#ifdef MHD
                        // calculate cell-averaged B^2 (magnetic field squared).
                        // idea: we get the cell-averaged values by averaging over pairs of opposing interfaces.
                        // remember that the B components are saved at the interfaces!
                        
                        // important note: there is NO a factor here as this gets calculated by hydro_source()
                        
                        B2 = 0.25 * (pow2(cur_node->B[X] + mhd_nodes[2][2][3]->B[X]) +
                                     pow2(cur_node->B[Y] + mhd_nodes[2][3][2]->B[Y]) +
                                     pow2(cur_node->B[Z] + mhd_nodes[3][2][2]->B[Z]) );
#else // MHD
                        B2 = 0.;
#endif // MHD
                        
#ifndef NO_HYDRO_SOURCES
                        /* obtain source terms (uses u_tmp[]!) */
                        hydro_source(timecounter, timestep, cur_node, Source, BSource, B2);
#endif // NO_HYDRO_SOURCES
                        
                        /* update hydro-variables */
                        for(iVar=0; iVar<NADVECT; iVar++)
                          {
                           /* advect... */
                           cur_node->u[iVar]  = cur_node->u_tmp[iVar] + dtdx2 * cur_node->F[iVar];
#ifndef NO_HYDRO_SOURCES
                           /* ...and add gravity terms */
                           cur_node->u[iVar] += timestep*Source[iVar];
#endif // NO_HYDRO_SOURCES
                          }
#ifdef MHD
                        // advect magnetic field - no source terms here! (if using super_B = a^2 * proper_B)
                        for(idim=0; idim<NDIM; idim++)
                        {
                          cur_node->B[idim] = cur_node->B_tmp[idim] + dtdx2 * cur_node->E[idim];
#ifndef NO_HYDRO_SOURCES
                           cur_node->B[idim] += timestep*BSource[idim];
#endif // NO_HYDRO_SOURCES
                        }
#endif // MHD
                        
                       } /* test_mhd() != NULL */

                    }
                 }
              }
           }
        }
     }

   /* refresh u[Uedens] using newly calculated hydro-variables */
   update_edens(cur_grid);

   /* 1. copy u_tmp[Uedens] over to u[Uedens] and simultaneously synchronize u[UEdens] and u[Untrpy]
    * 2. this loop also calculates the new global.cfl_speed using the latest u[] values! */
   getCFLspeed_and_syncES(cur_grid);   
#ifdef DEBUGMHD
   // output max(divB) for checking conservation
   fprintf(stderr,"divB/|B| = %16.8g\n", calc_divB(cur_grid));
#endif
}

/*===============================================================================
 * simply take the average of u[] and u_tmp[] stored at each node
 *
 * u_tmp[] = averaged u
 * u[]     = the old u[] (unchanged!)
 *===============================================================================*/
void ave_hydro(gridls *cur_grid)
{
   pqptr         cur_pquad;
   cqptr         cur_cquad, icur_cquad;
   nqptr         cur_nquad, icur_nquad;
   nptr          cur_node;
   double        u_tmp;
   long          x, y, z;
   int           iVar, idim;
   
   long          icquad;

   /* average old and new hydro variables (store result in u_tmp[]!) */
   for(cur_pquad=cur_grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next)
     {
      z = cur_pquad->z;
#ifdef WITH_OPENMPM
#pragma omp parallel shared(cur_grid, cur_pquad) private(icquad, cur_cquad, icur_cquad, cur_nquad, icur_nquad, cur_node, x, y, z, idim, iVar)
#pragma omp for schedule(static)
      for(icquad=0; icquad<cur_grid->l1dim; icquad++)
        {
         z         = cur_pquad->z   + icquad;
         cur_cquad = cur_pquad->loc + icquad;
#else
         for(cur_cquad = cur_pquad->loc; cur_cquad < cur_pquad->loc + cur_pquad->length; cur_cquad++, z++)  
           {  
#endif
            for(icur_cquad  = cur_cquad; icur_cquad != NULL; icur_cquad  = icur_cquad->next)
           {
            y = icur_cquad->y;
            for(cur_nquad = icur_cquad->loc; cur_nquad < icur_cquad->loc + icur_cquad->length; cur_nquad++, y++) 
              { 
               for(icur_nquad  = cur_nquad; icur_nquad != NULL; icur_nquad  = icur_nquad->next)
                 {
                  x = icur_nquad->x;
                  for(cur_node = icur_nquad->loc; cur_node < icur_nquad->loc + icur_nquad->length; cur_node++, x++)
                    {
                     /* store u^n+1/2 in u_tmp[] 
                      * NOTE: we use u_tmp[Uedens] differently in update_edens()! */
                     for(iVar=0; iVar<NADVECT; iVar++)
                       {
                        cur_node->u_tmp[iVar] = 0.5*(cur_node->u[iVar]+cur_node->u_tmp[iVar]);
                       }
#ifdef MHD
                     for(idim=0; idim<NDIM; idim++)   // analogous to hydro
                       {
                        cur_node->B_tmp[idim] = 0.5*(cur_node->B[idim]+cur_node->B_tmp[idim]);
                       }
#endif
                    }
                 }
              }
           }
        }
     }
}

/*===============================================================================
 * solve_dom_hydro:
 * ----------------
 *
 * this routine is a wrapper for the 2nd order Runge-Kutta time integration
 *
 * NOTE: this wrapper >>only<< works for vanishing source terms as
 *       there are no calls to solve_gravity()!
 *===============================================================================*/
void solve_dom_hydro(double timestep)
{
   gridls *cur_grid;
   
   cur_grid = global.dom_grid;
   
   calc_Flux(cur_grid);
   advect_hydro(cur_grid, timestep);
   ave_hydro(cur_grid);
   calc_Flux(cur_grid);
   advect_gravity_hydro(cur_grid, 0.0, timestep);
}


/*===============================================================================
 * init_DMhydro:
 * -------------
 *
 * initializes the hydro-variables on a given grid; 
 *  can only be used after startrun(), generate_grids(), & ll() !
 * 
 * init_DMhydro() simply assigns mass and velocity to the grid according 
 * to the dark matter particles values...
 *
 *
 * input:
 *         n/a
 *
 * output:
 *         hydro-variables on the domain grid
 *         + re-adjusted global variables (simu.pmass, etc.)
 *
 *
 * NOTE:
 *          USE ONLY ON DOMAIN GRID AND ONLY AFTER ll()!
 *
 *===============================================================================*/
void init_DMhydro()
{
   gridls       *cur_grid;
   partptr       cur_part;
   pqptr         cur_pquad;
   cqptr         cur_cquad, icur_cquad;
   nqptr         cur_nquad, icur_nquad;
   nptr          cur_node;
   
   double        weights[3][NDIM];
   nptr          tsc_nodes[3][3][3];
   double        xyz_coords[NDIM];
   double        temp_coords[NDIM];
   double        pn_sep[NDIM];
   double        pnarg_a;
   double        pnarg_b;
   double        tpnarg;
   long          x,y,z;
   double        www, gas_dens, mean_dens;
   double        cur_shift;
   int           iVar, idim, i, j, k;
   
   double        p, edens, Tdens, Edens, Bdens, molecular_weight;
   
   /* USE THIS ROUTINE ONLY ON DOMAIN GRID */
   cur_grid = global.dom_grid;
   
#ifdef NO_GAS
   simu.f_b = 0.0;
#endif
   
#ifdef DM_GAS
   simu.omegab  = 0.5;
   simu.omegaDM = 0.5;
   simu.omega0  = simu.omegab+simu.omegaDM;
   simu.f_b     = simu.omegab/simu.omega0;
   fprintf(stderr,"DM_GAS simulation:   omgegab  = %g\n",simu.omegab);
   fprintf(stderr,"                     omgegaDM = %g\n",simu.omegaDM);
   fprintf(stderr,"                     f_b       =%g\n",simu.f_b);
#endif
   
   /*---------------------------------------------------------------------
    * as we transfer mass from the DM particles to the gas
    * we end up with a new mass unit...
    *---------------------------------------------------------------------*/  
   /* re-adjust the mass unit */
   simu.pmass *= (1.-simu.f_b);
   
#ifdef NO_DM
   simu.f_b = 1.0;
#endif
   
   fprintf(io.logfile,"-> init_DMhydro(): initializing HYDRO simulation based upon DM initial conditions file\n");
   
   cur_shift = 0.5/(double)cur_grid->l1dim;
   
   
   /* reset all hydro variables */
   for(cur_pquad=cur_grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next)
     {
      z = cur_pquad->z;
      for(cur_cquad = cur_pquad->loc; cur_cquad < cur_pquad->loc + cur_pquad->length; cur_cquad++, z++)  
        {  
         for(icur_cquad  = cur_cquad; icur_cquad != NULL; icur_cquad  = icur_cquad->next)
           {
            y = icur_cquad->y;
            for(cur_nquad = icur_cquad->loc; cur_nquad < icur_cquad->loc + icur_cquad->length; cur_nquad++, y++) 
              { 
               for(icur_nquad  = cur_nquad; icur_nquad != NULL; icur_nquad  = icur_nquad->next)
                 {
                  x = icur_nquad->x;
                  for(cur_node = icur_nquad->loc; cur_node < icur_nquad->loc + icur_nquad->length; cur_node++, x++)
                    {
                     for(iVar=0; iVar<NHYDRO; iVar++)
                        cur_node->u[iVar] = 0.0;
                    }
                 }
              }
           }
        }
     }
   
   /* initialize hydro variables */
   for(cur_pquad=cur_grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next)
     {
      z = cur_pquad->z;
      for(cur_cquad = cur_pquad->loc; cur_cquad < cur_pquad->loc + cur_pquad->length; cur_cquad++, z++)  
        {  
         for(icur_cquad  = cur_cquad; icur_cquad != NULL; icur_cquad  = icur_cquad->next)
           {
            y = icur_cquad->y;
            for(cur_nquad = icur_cquad->loc; cur_nquad < icur_cquad->loc + icur_cquad->length; cur_nquad++, y++) 
              { 
               for(icur_nquad  = cur_nquad; icur_nquad != NULL; icur_nquad  = icur_nquad->next)
                 {
                  x = icur_nquad->x;
                  for(cur_node = icur_nquad->loc; cur_node < icur_nquad->loc + icur_nquad->length; cur_node++, x++)
                    {
                     /* get TSC neighbours of cur_node */
                     tsc_nodes[1][1][1] = cur_node;
                     get_TSCnodes(cur_grid, cur_pquad, icur_cquad, icur_nquad, tsc_nodes, &z, &y, &x);
                     
                     /* calculate realspace coords of cur_node */
                     temp_coords[X] = (((double)x) / (double)cur_grid->l1dim) + cur_shift;
                     temp_coords[Y] = (((double)y) / (double)cur_grid->l1dim) + cur_shift;
                     temp_coords[Z] = (((double)z) / (double)cur_grid->l1dim) + cur_shift;
                     
                     xyz_coords[X]  = f1mod(temp_coords[X]+1.0, 1.0);
                     xyz_coords[Y]  = f1mod(temp_coords[Y]+1.0, 1.0);
                     xyz_coords[Z]  = f1mod(temp_coords[Z]+1.0, 1.0);                     
                     
                     /* loop over ll */
                     for(cur_part = cur_node->ll; cur_part != NULL; cur_part = cur_part->ll)
                       {
                        /* calc fraction to be assigned to each tsc node */
                        for(idim = 0; idim < 3; idim++)
                          {
                           pn_sep[idim] = ((double)cur_part->pos[idim] - xyz_coords[idim]) * (double)cur_grid->l1dim;
                           
                           /* take care of periodic boundaries */
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
                          }
                        
                        /* density of virtual gas particle */
#ifdef MULTIMASS
                        gas_dens = simu.f_b * (double)cur_part->weight * cur_grid->masstodens;
#else
                        gas_dens = simu.f_b * cur_grid->masstodens;
#endif
                        
                        /* assign mass */
                        for(k = 0; k < 3; k++)
                           for(j = 0; j < 3; j++)
                              for(i = 0; i < 3; i++)
                                 if(tsc_nodes[k][j][i] != NULL)
                                   {
                                    /* only a certain fraction goes to this node... */
                                    www = (weights[k][Z]*weights[j][Y]*weights[i][X]);
                                    
                                    /*-----------------------------------------------
                                     * assign baryon fraction f_b of dm mass to grid
                                     *-----------------------------------------------*/
                                    tsc_nodes[k][j][i]->u[Udens]     += www * gas_dens;
                                    
                                    /*-----------------------------------------------
                                     * assign velocity to grid
                                     *-----------------------------------------------*/
                                    tsc_nodes[k][j][i]->u[UmomdensX] += www * gas_dens * (double)cur_part->mom[X];
                                    tsc_nodes[k][j][i]->u[UmomdensY] += www * gas_dens * (double)cur_part->mom[Y];
                                    tsc_nodes[k][j][i]->u[UmomdensZ] += www * gas_dens * (double)cur_part->mom[Z];
                                   }
                       }
                              
                              
                    }
                 }
              }
           }
        }
     }
         
#ifdef MHD
   // initialise magnetic field
   init_B(cur_grid);
#endif           
         
   /* initialize energies, entropy, etc. */
   mean_dens = 0.0;
   for(cur_pquad=cur_grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next)
     {
      z = cur_pquad->z;
      for(cur_cquad = cur_pquad->loc; cur_cquad < cur_pquad->loc + cur_pquad->length; cur_cquad++, z++)  
        {  
         for(icur_cquad  = cur_cquad; icur_cquad != NULL; icur_cquad  = icur_cquad->next)
           {
            y = icur_cquad->y;
            for(cur_nquad = icur_cquad->loc; cur_nquad < icur_cquad->loc + icur_cquad->length; cur_nquad++, y++) 
              { 
               for(icur_nquad  = cur_nquad; icur_nquad != NULL; icur_nquad  = icur_nquad->next)
                 {
                  x = icur_nquad->x;
                  for(cur_node = icur_nquad->loc; cur_node < icur_nquad->loc + icur_nquad->length; cur_node++, x++)
                    {
                     /* internal energy density */
                     edens = cur_node->u[Udens]*simu.e_init;
                     Tdens = 0.5 * ( pow2(cur_node->u[UmomdensX])+
                                     pow2(cur_node->u[UmomdensY])+
                                     pow2(cur_node->u[UmomdensZ]) ) / cur_node->u[Udens];
#ifdef MHD
                     tsc_nodes[1][1][1] = cur_node;
                     get_TSCnodes(cur_grid, cur_pquad, icur_cquad, icur_nquad, tsc_nodes, &z, &y, &x);
                     
                     // calculate Bdens = 0.5*B^2 (magnetic energy density).
                     // idea: we get the cell-averaged values by averaging over pairs of opposing interfaces.
                     // remember that the B components are saved at the interfaces!         
                                          
                     if(test_tsc(tsc_nodes) == TRUE)                       
                       
                     Bdens = 0.125 *  (pow2(cur_node->B[X] + tsc_nodes[1][1][2]->B[X]) +
                                       pow2(cur_node->B[Y] + tsc_nodes[1][2][1]->B[Y]) +
                                       pow2(cur_node->B[Z] + tsc_nodes[2][1][1]->B[Z]) );
                                          
                     else
                     Bdens = 0.5 * ( pow2(cur_node->B[X]) + pow2(cur_node->B[Y]) + pow2(cur_node->B[Z]) );
#else
                     Bdens = 0.0;
#endif
                     
                     /* total energy = internal + kinetic + magnetic energy density */
                     cur_node->u[UEdens] = Tdens + edens + Bdens;
                                             
                     /* S = (gamma-1) edens / dens^(gamma-1) */
                     cur_node->u[Untrpy] = (simu.gamma-1.)*simu.e_init/pow(cur_node->u[Udens],(simu.gamma-1.));

                     /* initialize internal energy */
                     cur_node->u[Uedens] = edens;

#ifdef NO_GAS
                     for(iVar=0; iVar<NHYDRO; iVar++)
                        cur_node->u[iVar] = 0.0;
#endif
                     
                     for(iVar=0; iVar<NHYDRO; iVar++)
                        cur_node->u_tmp[iVar] = cur_node->u[iVar];
                     
                     mean_dens += cur_node->u[Udens];
                    }
                 }
              }
           }
        }
     }   
}
   
///////////////////////////////////////////////////////////////////////////////
// init_B:                                                                   //
// sets up an isotropic initial B field with average strength B_init         //
// from a vector potential A, where curl A ~ B.                              //
///////////////////////////////////////////////////////////////////////////////

#ifdef MHD   
   
void init_B(gridls* cur_grid)
{
   //set of variables needed to loop over cur_grid
   partptr       cur_part;
   pqptr         cur_pquad;
   cqptr         cur_cquad, icur_cquad;
   nqptr         cur_nquad, icur_nquad;
   nptr          cur_node;
   long          x, y, z;
   
   //variables needed for the B_init scheme
   double        B_init, B_fac;
   nptr          tsc_nodes[3][3][3];
   double        A[3][3][3][NDIM];
   int           i,j,k;
   double        sum_norm_B, new_norm_B;
   int           number_of_cells;
   
   //unit conversion to AMIGA units and supercomoving coordinates (refer to documentation)
   B_fac = 6.50778e5 * pow(simu.a_initial,2.5) / (sqrt(simu.omega0) * simu.boxsize );
   
   //set the initial magnetic field here; give the field strength in Gauss
   B_init = 1e-3 * B_fac;  // TODO: read it in somewhere else (.input file?). otherwise we have to recompile for every value change...
   
   //initialise variables needed for normalisation
   sum_norm_B = 0.;
   number_of_cells = 0;
   
   ////// loop 1: calculate initial B from vector potential A. ////////////////////////////////////////////////////
   
   for(cur_pquad=cur_grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next)
   {
      z = cur_pquad->z;
      for(cur_cquad = cur_pquad->loc; cur_cquad < cur_pquad->loc + cur_pquad->length; cur_cquad++, z++)  
      {  
         for(icur_cquad  = cur_cquad; icur_cquad != NULL; icur_cquad  = icur_cquad->next)
         {
            y = icur_cquad->y;
            for(cur_nquad = icur_cquad->loc; cur_nquad < icur_cquad->loc + icur_cquad->length; cur_nquad++, y++) 
            { 
               for(icur_nquad  = cur_nquad; icur_nquad != NULL; icur_nquad  = icur_nquad->next)
               {
                  x = icur_nquad->x;
                  for(cur_node = icur_nquad->loc; cur_node < icur_nquad->loc + icur_nquad->length; cur_node++, x++)
                  {
                     tsc_nodes[1][1][1] = cur_node;
                     get_TSCnodes(cur_grid,cur_pquad,icur_cquad,icur_nquad,tsc_nodes, &z, &y, &x);
                     if(test_tsc(tsc_nodes) == TRUE)
                     {
                        // initialise the vector potential A with the initial momentum field
                        
                        for (k=0;k<3;k++)
                        {
                           for (j=0;j<3;j++)
                           {
                              for (i=0;i<3;i++)
                              {                                 
                                 A[k][j][i][X] = tsc_nodes[k][j][i]->u[UmomdensX];
                                 A[k][j][i][Y] = tsc_nodes[k][j][i]->u[UmomdensY];
                                 A[k][j][i][Z] = tsc_nodes[k][j][i]->u[UmomdensZ];
                              }
                           }
                        }
                        
                        
                        //set the staggered magnetic field to curl A (note: the way this scheme is designed ensures div(curl A) = 0 )      
                        
                        cur_node->B[X] =  0.25*B_init * (  A[1][2][0][Z] - A[1][0][0][Z] + A[1][2][1][Z] - A[1][0][1][Z]
                                                         - A[2][1][0][Y] + A[0][1][0][Y] - A[2][1][1][Y] + A[0][1][1][Y] );

                        cur_node->B[Y] =  0.25*B_init * (  A[2][0][1][X] - A[0][0][1][X] + A[2][1][1][X] - A[0][1][1][X]
                                                         - A[1][0][2][Z] + A[1][0][0][Z] - A[1][1][2][Z] + A[1][1][0][Z] );

                        cur_node->B[Z] =  0.25*B_init * (  A[0][1][2][Y] - A[0][1][0][Y] + A[1][1][2][Y] - A[1][1][0][Y]
                                                         - A[0][2][1][X] + A[0][0][1][X] - A[1][2][1][X] + A[1][0][1][X] );
                        
                        //the following is needed to normalise B in the next loop
                        sum_norm_B += sqrt (cur_node->B[X]*cur_node->B[X] + cur_node->B[Y]*cur_node->B[Y] + cur_node->B[Z]*cur_node->B[Z]);
                        number_of_cells++;
                     }
                     else
                     {
                        fprintf(stderr,"hydro.c: function init_B: failed to set up B. neighbour cell not found.\n");
                        exit(1);
                     }
                  }
               }
            }
         }
      }
   }
   
   
   ////// loop 2: normalise B such that the average length of B vector equals B_init //////////////////////////////
   
   new_norm_B = B_init / (sum_norm_B / number_of_cells);

   for(cur_pquad=cur_grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next)
   {
      z = cur_pquad->z;
      for(cur_cquad = cur_pquad->loc; cur_cquad < cur_pquad->loc + cur_pquad->length; cur_cquad++, z++)  
      {  
         for(icur_cquad  = cur_cquad; icur_cquad != NULL; icur_cquad  = icur_cquad->next)
         {
            y = icur_cquad->y;
            for(cur_nquad = icur_cquad->loc; cur_nquad < icur_cquad->loc + icur_cquad->length; cur_nquad++, y++) 
            { 
               for(icur_nquad  = cur_nquad; icur_nquad != NULL; icur_nquad  = icur_nquad->next)
               {
                  x = icur_nquad->x;
                  for(cur_node = icur_nquad->loc; cur_node < icur_nquad->loc + icur_nquad->length; cur_node++, x++)
                  {
                     cur_node->B[X] *= new_norm_B;
                     cur_node->B[Y] *= new_norm_B;
                     cur_node->B[Z] *= new_norm_B;
                     
                     // un-comment next 3 lines to use a simple constant B instead of this scheme
//                     cur_node->B[X] = B_init;
//                     cur_node->B[Y] = 0.;
//                     cur_node->B[Z] = 0.;
                  }
               }
            }
         }
      }
   }
}
#endif // MHD //   

/*===================================================================================
 * vanLeer:  vanLeer function for flux limitation
 *===================================================================================*/
double vanLeer(double a, double b)
{
   double ab;
   
   ab = a*b;
   
   if(ab > MACHINE_ZERO)      // keep this limit in mind if dealing with small B fields...
      return( 2*ab/(a+b) );
   else
      return( 0.0 );
}

/*===================================================================================
 * minmod:   minmod function for flux limitation
 *===================================================================================*/
double minmod(double a, double b)
{
   double sgn;
   
   sgn = (sign(a)+sign(b))/2;
   
   return( sgn * MIN(fabs(a), fabs(b)) );
}

/*===================================================================================
 * superbee:   superbee function for flux limitation !!! NOT FUNCTIONAL !!!
 *===================================================================================*/
double superbee(double a, double b)
{
   if(fabs(a) > fabs(b))
      return ( minmod(a  , 2*b) );
   else
      return ( minmod(2*a,   b) );
}

//===============================================================================
// zero_B_infc sets all elements of a B_infc[NDIM][2][NDIM] array to 0.
//===============================================================================

void zero_B_infc(double B_infc[NDIM][2][NDIM])
{
   int i,j,k;
   for (i=0;i<3;i++) {
      for (j=0;j<2;j++) {
         for (k=0;k<NDIM;k++) {
            B_infc[i][j][k] = 0.;
         }
      }
   }
}

#endif /* HYDRO */
