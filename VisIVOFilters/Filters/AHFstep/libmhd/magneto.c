#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void dummy_magneto()
{
}

#ifdef MHD

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "mhd.h"
#include "../libutility/utility.h"
#include "../libgrids/grids.h"


//===============================================================================
// calc_FluxKNPCT:
//---------------
//
// does the same as calc_FluxKNP, but in addition to that it calculates the time
// change of the three B field components as stored at the central node - utilizing
// the constrained transport scheme as in Ziegler 2004.
//
// input:  as in calc_FluxKNP
// output: as in calc_FluxKNP *and* the numerical dB/dT for B[X], B[Y] and B[Z]
//           for the central node (mind that we save them on a "staggered grid"!)
//===============================================================================*/

void calc_FluxKNPCT(nptr MHDnodes[5][5][5], double Flux[NDIM][2][NHYDRO], double dBdt[NDIM])
{
   nptr     stencil[3];                        // temporary cell stencil array for calc_u_infc()
   nptr     stencil3D[3][3][3];                // temporary cell stencil array for calc_B_infc()
   
   double   u_infc[NDIM][5][3][2][NHYDRO];     // indexes: [infc direction] [row number] [cell in row: i-1,i,i+1] [E/W] [components]
   double   B_infc[NDIM][5][3][2][NDIM];        
   double   E_infc[NDIM][5][3][2][NDIM];
   double   flux_infc[NDIM][5][3][2][NHYDRO];
   
   double   a[NDIM][5][2][2];                  // indexes: [infc direction] [row number] [left or right] [a- or a+]
   double   E_tensor[3][2][NDIM];              // indexes: [i-1,i,i+1] [E,W] [NHYDRO]. same as in calc_FluxKNP() !
   double   MHDFlux[NDIM][5][2][NDIM];         // indexes: [infc direction] [row number] [left or right] [components]
   double   E_edge[9];                         // indexes: the 9 relevant cell edges (see calc_E_edge_Ziegler for explanation)
   
   int    Umomdens[NDIM];
   int    idim, row, cell;
   int    x,y,z;
   
#ifdef CT_GARDINER_STONE
   double complete_Flux[NDIM][5][2][NHYDRO];      // rhovx, rhovy and rhovz on all interfaces
   double E_ref[3][3][3][NDIM];  // reference field
#endif
   
   /* initialize required matrizes and vectors */
   Umomdens[X] = UmomdensX;
   Umomdens[Y] = UmomdensY;
   Umomdens[Z] = UmomdensZ;
   
#ifdef CT_GARDINER_STONE
   // calculate "reference field" from initial conditions
   calc_E_ref (MHDnodes, E_ref);   
#endif
   
   // calculate u_infc, B_infc and E_infc on all interfaces
   
   for ( idim = 0; idim < NDIM; idim++ ) {
      for ( row = 0; row < 5; row++ )       {
         for ( cell = 0; cell < 3; cell++ )    {
            
            transformDRCtoXYZ (idim, row, cell, &x, &y, &z);      // get correct MHDnodes[z][y][x] coordinate of current cell
            
            getStencil ( MHDnodes, x,y,z, stencil, idim );        // get array with direct neighbours in direction idim
            calc_u_infc ( stencil, u_infc[idim][row][cell] );     // calc hydro quantities on interface from cell-centered quantities
       
            getStencil3D ( MHDnodes, x,y,z, stencil3D );                // get array with direct neighbours in 3D
            calc_B_infc ( stencil3D, B_infc[idim][row][cell], idim );   // calc B field components on interface from staggered grid
            
            calc_E_infc ( u_infc[idim][row][cell], B_infc[idim][row][cell], E_infc[idim][row][cell]);   // calc E = -v x B on every interface
         }
         
         min_max_speeds ( u_infc[idim][row], B_infc[idim][row], a[idim][row], Umomdens[idim] );   // calculate min/max speeds inside each 3-stencil
      }
   }
   
   for ( idim = 0; idim < NDIM; idim++ ) {
      for ( row = 0; row < 5; row++ )       {
         hydro_flux(u_infc[idim][row][0][1], B_infc[idim][row][0][1], flux_infc[idim][row][0][1], idim);
         hydro_flux(u_infc[idim][row][1][0], B_infc[idim][row][1][0], flux_infc[idim][row][1][0], idim);
         hydro_flux(u_infc[idim][row][1][1], B_infc[idim][row][1][1], flux_infc[idim][row][1][1], idim);
         hydro_flux(u_infc[idim][row][2][0], B_infc[idim][row][2][0], flux_infc[idim][row][2][0], idim);
         
#ifdef CT_GARDINER_STONE
         // calc Flux on ALL interfaces, needed for Gardiner/stone speed...
         calc_fluxKNP_infc(flux_infc[idim][row], u_infc[idim][row], a[idim][row], complete_Flux[idim][row]);
#endif
         if (row == 0)
         calc_fluxKNP_infc(flux_infc[idim][0], u_infc[idim][0], a[idim][0], Flux[idim]);  //TODO: more effective, less redundant
      }
   }
      
   
   // at this point the hydro scheme is finished. continuing with constrained transport scheme.   
   
   // calculate the numerical MHD fluxes:

   for ( idim = 0; idim < NDIM; idim++ ) {
      for ( row = 0; row < 5; row++ )       {
         
         // calculate idim-th column of "antisymmetric flux tensor" (MHD counterpart to hydro_flux) from E_infc
         MHD_flux_tensor(E_infc[idim][row][0][1], E_tensor[0][1], idim);   
         MHD_flux_tensor(E_infc[idim][row][1][0], E_tensor[1][0], idim);
         MHD_flux_tensor(E_infc[idim][row][1][1], E_tensor[1][1], idim);
         MHD_flux_tensor(E_infc[idim][row][2][0], E_tensor[2][0], idim);
         
         // apply KNP flux formula to calcluate MHD interface fluxes (formally this is exactly the same as calc_fluxKNP_infc)
         calc_MHDfluxKNP_infc(E_tensor, B_infc[idim][row], a[idim][row], MHDFlux[idim][row]);
      }
   }
   
   // calculate edge-on electric fields from MHD interface fluxes:
   
#ifdef CT_GARDINER_STONE
   calc_E_edge_GardinerStone(complete_Flux, E_ref, MHDFlux, E_edge);
#else
   calc_E_edge_Ziegler(MHDFlux, E_edge);
#endif
   
   // calculate dB/dt from edge-on electric fields and return that:
   
   calc_dBdt(E_edge, dBdt);
}


//=======================================================================================================
// calc_E_ref:
// -----------
//
// not finished yet...
//=======================================================================================================


void calc_E_ref ( nptr MHDnodes[5][5][5], double E_ref[3][3][3][NDIM] )
{
   int i,j,k;
   double B_cell[NDIM]; //temporary variable to hold CELL-CENTERED magnetic field
   
   for (k=0;k<3;k++) {
      for (j=0;j<3;j++) {
         for (i=0;i<3;i++) {
            
            B_cell[X] = 0.5 * ( MHDnodes[k+1][j+1][i+1]->B[X] + MHDnodes[k+1][j+1][i+2]->B[X] );
            B_cell[Y] = 0.5 * ( MHDnodes[k+1][j+1][i+1]->B[Y] + MHDnodes[k+1][j+2][i+1]->B[Y] );
            B_cell[Z] = 0.5 * ( MHDnodes[k+1][j+1][i+1]->B[Z] + MHDnodes[k+2][j+1][i+1]->B[Z] );
            
            calc_E(MHDnodes[k+1][j+1][i+1]->u, B_cell, E_ref[k][j][i]);
         }
      }
   }
}



//=======================================================================================================
// calc_E_infc:
// ------------
//
// calculate E-field on all relevant interfaces (E,W,N,S,T,B) from B field and hydro variables.
//
// input:  u_infc and B_infc in one cell: left and right interfaces
// output: E_infc at the same interfaces
// notes:  - all three arrays must have matching indexes with the same definitions!
//         - the formula for E(u,B) is not important here, this is hidden in function calc_E.
//=======================================================================================================


void calc_E_infc ( double u_infc[2][NHYDRO], double B_infc[2][NDIM], double E_infc[2][NDIM] )
{
   calc_E ( u_infc[0], B_infc[0], E_infc[0] );    // left interface (W/S/B)
   calc_E ( u_infc[1], B_infc[1], E_infc[1] );    // right interface (E/N/T)
}


//=======================================================================================================
// calc_B_infc:
// ------------
//
// calculate reconstructed B-field (Bx, By, Bz) at two opposing cell interfaces of
// the central cell (at [1][1][1]) in a given stencil3D[3][3][3] array.
//
// input: the stencil3D[3][3][3] array and a spatial direction specified by idim (X, Y or Z).
//
// output: B_infc[2][NDIM]  containing Bx, By, Bz for two opposing interfaces in direction idim
//         Explanation:
//           [2]    = the two interfaces: [0] = W, S or B (backward); [1] = E, N or T (forward).
//           [NDIM] = which component of the B field? possible values: [X], [Y], [Z]
//           idim   = interfaces in which spacial direction? [X] = E or W; [Y] = N or S; [Z] = T or B.
//
// note:   the difference to calc_u_infc() is that this one can not be written in a general "1D" version
//           so we have to put the switch/case clause here.
//
// Reference for the equations used: U. Ziegler, Journal of Computational Physics 196 (2004), 393-416
//=======================================================================================================

void calc_B_infc(nptr stencil3D[3][3][3], double B_infc[2][NDIM], int idim)
{
   double B_node[3], dB_right, dB_left;
   int i;
   
   // reconstruction of the B field (chapter 3.2) - explicitly written down for every component
   // Bx, By, Bz can be calculated for all six interfaces (E,W,N,S,T,B) of cell (ijk), which is stencil3D[1][1][1]
   // However, only the two interfaces lying in the spatial direction specified by idim are calculated (switch).
   
   switch ( idim )
   {
   
   case X:
   
   B_infc[1][X] = stencil3D[1][1][2]->B[X];      // Bx^E
   B_infc[0][X] = stencil3D[1][1][1]->B[X];      // Bx^W
   
   for(i=0;i<3;i++)   B_node[i] = stencil3D[1][2][i]->B[Y]; dB_right = du_ltd(B_node) / 2.;   // dxBy_(i,j+0.5,k)
   for(i=0;i<3;i++)   B_node[i] = stencil3D[1][1][i]->B[Y]; dB_left  = du_ltd(B_node) / 2.;   // dxBy_(i,j-0,5,k)
   B_infc[1][Y] = 0.5 * ( stencil3D[1][2][1]->B[Y] + stencil3D[1][1][1]->B[Y] + dB_right + dB_left );    // By^E
   B_infc[0][Y] = 0.5 * ( stencil3D[1][2][1]->B[Y] + stencil3D[1][1][1]->B[Y] - dB_right - dB_left );    // By^W
   
   for(i=0;i<3;i++) B_node[i] = stencil3D[2][1][i]->B[Z]; dB_right = du_ltd(B_node) / 2.;   // dxBz_(i,j,k+0.5)
   for(i=0;i<3;i++) B_node[i] = stencil3D[1][1][i]->B[Z]; dB_left  = du_ltd(B_node) / 2.;   // dxBz_(i,j,k-0.5)
   B_infc[1][Z] = 0.5 * ( stencil3D[2][1][1]->B[Z] + stencil3D[1][1][1]->B[Z] + dB_right + dB_left );    // Bz^E
   B_infc[0][Z] = 0.5 * ( stencil3D[2][1][1]->B[Z] + stencil3D[1][1][1]->B[Z] - dB_right - dB_left );    // Bz^W
   
   break;
   
   case Y:
   
   for(i=0;i<3;i++)   B_node[i] = stencil3D[1][i][2]->B[X]; dB_right = du_ltd(B_node) / 2.;   // dyBx_(i+0.5,j,k)
   for(i=0;i<3;i++)   B_node[i] = stencil3D[1][i][1]->B[X]; dB_left  = du_ltd(B_node) / 2.;   // dyBx_(i-0.5,j,k))
   B_infc[1][X] = 0.5 * ( stencil3D[1][1][2]->B[X] + stencil3D[1][1][1]->B[X] + dB_right + dB_left );   // Bx^N
   B_infc[0][X] = 0.5 * ( stencil3D[1][1][2]->B[X] + stencil3D[1][1][1]->B[X] - dB_right - dB_left );   // Bx^S
   
   B_infc[1][Y] = stencil3D[1][2][1]->B[Y];     // By^N
   B_infc[0][Y] = stencil3D[1][1][1]->B[Y];     // By^S
   
   for(i=0;i<3;i++)   B_node[i] = stencil3D[2][i][1]->B[Z]; dB_right = du_ltd(B_node) / 2.;   // dyBz_(i,j,k+0.5)
   for(i=0;i<3;i++)   B_node[i] = stencil3D[1][i][1]->B[Z]; dB_left  = du_ltd(B_node) / 2.;   // dyBz_(i,j,k-0.5)
   B_infc[1][Z] = 0.5 * ( stencil3D[2][1][1]->B[Z] + stencil3D[1][1][1]->B[Z] + dB_right + dB_left );   // Bz^N
   B_infc[0][Z] = 0.5 * ( stencil3D[2][1][1]->B[Z] + stencil3D[1][1][1]->B[Z] - dB_right - dB_left );   // Bz^S
   
   break;
   
   case Z:
   
   for(i=0;i<3;i++)   B_node[i] = stencil3D[i][1][2]->B[X]; dB_right = du_ltd(B_node) / 2.;   // dzBx_(i+0.5,j,k)
   for(i=0;i<3;i++)   B_node[i] = stencil3D[i][1][1]->B[X]; dB_left  = du_ltd(B_node) / 2.;   // dzBx_(i-0.5,j,k)
   B_infc[1][X] = 0.5 * ( stencil3D[1][1][2]->B[X] + stencil3D[1][1][1]->B[X] + dB_right + dB_left );    // Bx^T
   B_infc[0][X] = 0.5 * ( stencil3D[1][1][2]->B[X] + stencil3D[1][1][1]->B[X] - dB_right - dB_left );    // Bx^B
   
   for(i=0;i<3;i++)   B_node[i] = stencil3D[i][2][1]->B[Y]; dB_right = du_ltd(B_node) / 2.;   // dzBy_(i,j+0.5,k)
   for(i=0;i<3;i++)   B_node[i] = stencil3D[i][1][1]->B[Y]; dB_left  = du_ltd(B_node) / 2.;   // dzBy_(i,j-0.5,k)
   B_infc[1][Y] = 0.5 * ( stencil3D[1][2][1]->B[Y] + stencil3D[1][1][1]->B[Y] + dB_right + dB_left );    // By^T
   B_infc[0][Y] = 0.5 * ( stencil3D[1][2][1]->B[Y] + stencil3D[1][1][1]->B[Y] - dB_right - dB_left );    // By^B  
   
   B_infc[1][Z] = stencil3D[2][1][1]->B[Z];    // Bz^T
   B_infc[0][Z] = stencil3D[1][1][1]->B[Z];    // Bz^B
   
   break;
   
   default:
   
   fprintf(stderr,"unexpected error in calc_B_infc: invalid spatial direction given\n");
   }
}


//===============================================================================
// MHD_flux_tensor: 
//   This calculates the "antisymmetric flux tensor" of the E field.
//   For the KNP flux formula (calc_fluxKNP), this quantity is the MHD analogue
//   to flux_idim and this function is the analogue to hydro_flux.
//
// input:  E_infc[NDIM]
// output: E_tensor[NDIM], set to the idim-th column of the antisymmetric
//         MHD flux tensor (Ziegler 2004, Chapter 3.3, section (i)).
//===============================================================================


void MHD_flux_tensor (double E_infc[NDIM], double E_tensor[NDIM], int idim)
{   
   switch (idim)
   {
      case X:
         E_tensor[X] = 0.;
         E_tensor[Y] = -E_infc[Z];
         E_tensor[Z] = E_infc[Y];
         break;
      case Y:
         E_tensor[X] = E_infc[Z];
         E_tensor[Y] = 0.;
         E_tensor[Z] = -E_infc[X];
         break;
      case Z:
         E_tensor[X] = -E_infc[Y];
         E_tensor[Y] = E_infc[X];
         E_tensor[Z] = 0.;
         break;
      default:
         fprintf(stderr,"unexpected error in MHD_flux_tensor: invalid spatial direction given\n");
   }
}

//===============================================================================
// calc_E_edge: 
//   This calculates the 12 E values at the cell edges by averaging over
//   the four adjacent interfaces, as in Ziegler 2004.
//
// input:  MHDflux[NDIM][5][3][2][NDIM] - the complete set of MHDflux interface values.
// output: E_edge[12], the 12 needed values of E on the cell edges.
//         For the equations, see Ziegler 2004, Chapter 3.3, section (iii).
//===============================================================================

void calc_E_edge_Ziegler (double MHDFlux[NDIM][5][2][NDIM], double E_edge[9])
{

/*
     _____________       So how do we manage these cell edges?
    /:           /|      There are 12 edges that comprise the central cell at ijk. However, as B is stored only
  8/ :          / |      at three interfaces per cell, only 9 edges touching these three interfaces are important.
  /_____2_____ /  |      We simply assign numbers to these 9 edges (see left). At each edge, only the E component
 |   7        |   |      in the same spatial direction as the edge is saved (x goes to the right, y goes inside
 |   :.....5..|...|      the screen, z goes up). So, for example, in Ziegler 2004 notation:
 3  .         1  /       
 | .6         | /4       E_edge[0] would be  Ex_(i,j-0.5,k-0.5)
 |_____0______|/         E_edge[1] would be  Ez_(i+0.5,j-0.5,k)   etc.                                    */


// Ziegler (2004), page 400, gives us three sets of edge-centered field composition
// rules for the 3D, 2D and 1D cases respectively. We implement them all:   

#ifdef CT_ZIEGLER3D
   E_edge[0] = 0.25 * ( - MHDFlux[Y][0][0][Z] - MHDFlux[Y][4][0][Z] + MHDFlux[Z][0][0][Y] + MHDFlux[Z][3][0][Y] );
   E_edge[1] = 0.25 * ( - MHDFlux[X][0][1][Y] - MHDFlux[X][4][1][Y] + MHDFlux[Y][1][0][X] + MHDFlux[Y][0][0][X] );
   E_edge[2] = 0.25 * ( - MHDFlux[Y][2][0][Z] - MHDFlux[Y][0][0][Z] + MHDFlux[Z][0][1][Y] + MHDFlux[Z][3][1][Y] );
   E_edge[3] = 0.25 * ( - MHDFlux[X][0][0][Y] - MHDFlux[X][4][0][Y] + MHDFlux[Y][0][0][X] + MHDFlux[Y][3][0][X] );
   E_edge[4] = 0.25 * ( + MHDFlux[X][0][1][Z] + MHDFlux[X][3][1][Z] - MHDFlux[Z][2][0][X] - MHDFlux[Z][0][0][X] );
   E_edge[5] = 0.25 * ( - MHDFlux[Y][0][1][Z] - MHDFlux[Y][4][1][Z] + MHDFlux[Z][1][0][Y] + MHDFlux[Z][0][0][Y] );
   E_edge[6] = 0.25 * ( + MHDFlux[X][0][0][Z] + MHDFlux[X][3][0][Z] - MHDFlux[Z][0][0][X] - MHDFlux[Z][4][0][X] );
   E_edge[7] = 0.25 * ( - MHDFlux[X][2][0][Y] - MHDFlux[X][0][0][Y] + MHDFlux[Y][0][1][X] + MHDFlux[Y][3][1][X] );
   E_edge[8] = 0.25 * ( + MHDFlux[X][1][0][Z] + MHDFlux[X][0][0][Z] - MHDFlux[Z][0][1][X] - MHDFlux[Z][4][1][X] );   
#endif   
#ifdef CT_ZIEGLER2D
   E_edge[0] = - MHDFlux[Y][0][0][Z];
   E_edge[1] = 0.25 * ( - MHDFlux[X][0][1][Y] - MHDFlux[X][4][1][Y] + MHDFlux[Y][1][0][X] + MHDFlux[Y][0][0][X] );
   E_edge[2] = - MHDFlux[Y][0][0][Z];
   E_edge[3] = 0.25 * ( - MHDFlux[X][0][0][Y] - MHDFlux[X][4][0][Y] + MHDFlux[Y][0][0][X] + MHDFlux[Y][3][0][X] );
   E_edge[4] = + MHDFlux[X][0][1][Z];
   E_edge[5] = - MHDFlux[Y][0][1][Z];
   E_edge[6] = + MHDFlux[X][0][0][Z];
   E_edge[7] = 0.25 * ( - MHDFlux[X][2][0][Y] - MHDFlux[X][0][0][Y] + MHDFlux[Y][0][1][X] + MHDFlux[Y][3][1][X] );
   E_edge[8] = + MHDFlux[X][0][0][Z];
#endif   
#ifdef CT_ZIEGLER1D
   E_edge[0] = 0.0;
   E_edge[1] = - MHDFlux[X][0][1][Y];
   E_edge[2] = 0.0;
   E_edge[3] = - MHDFlux[X][0][0][Y];
   E_edge[4] = + MHDFlux[X][0][1][Z];
   E_edge[5] = 0.0;
   E_edge[6] = + MHDFlux[X][0][0][Z];
   E_edge[7] = - MHDFlux[X][0][0][Y];
   E_edge[8] = + MHDFlux[X][0][0][Z];
#endif
}

//===============================================================================
// calc_E_edge_GardinerStone: 
//   does the same as calc_E_edge_Ziegler, but with additional corrector terms
//   as defined in the CT algorithm of Gardiner/Stone 2008
//   for this it needs additional input: the numerical hydro flux at all interfaces
//   and the reference field E_ref.
//===============================================================================

void calc_E_edge_GardinerStone (double complete_Flux[NDIM][5][2][NHYDRO], double E_ref[3][3][3][NDIM], double MHDFlux[NDIM][5][2][NDIM], double E_edge[9])
{
   E_edge[0] = 0.25 * ( - MHDFlux[Y][0][0][Z] - MHDFlux[Y][4][0][Z] + MHDFlux[Z][0][0][Y] + MHDFlux[Z][3][0][Y] 
                        + upwind_GardinerStone( complete_Flux[Y][0][0][Udens] ,  MHDFlux[Z][0][0][Y]-E_ref[0][1][1][X],  MHDFlux[Z][3][0][Y]-E_ref[0][0][1][X])
                        + upwind_GardinerStone( complete_Flux[Y][4][0][Udens] ,  MHDFlux[Z][0][0][Y]-E_ref[1][1][1][X],  MHDFlux[Z][3][0][Y]-E_ref[1][0][1][X])
                        + upwind_GardinerStone( complete_Flux[Z][0][0][Udens] , -MHDFlux[Y][0][0][Z]-E_ref[1][0][1][X], -MHDFlux[Y][4][0][Z]-E_ref[0][0][1][X])
                        + upwind_GardinerStone( complete_Flux[Z][3][0][Udens] , -MHDFlux[Y][0][0][Z]-E_ref[1][1][1][X], -MHDFlux[Y][4][0][Z]-E_ref[0][1][1][X]) );   

   E_edge[1] = 0.25 * ( - MHDFlux[X][0][1][Y] - MHDFlux[X][4][1][Y] + MHDFlux[Y][1][0][X] + MHDFlux[Y][0][0][X]
                        + upwind_GardinerStone( complete_Flux[X][0][1][Udens] ,  MHDFlux[Y][1][0][X]-E_ref[1][0][2][Z],  MHDFlux[Y][0][0][X]-E_ref[1][0][1][Z])
                        + upwind_GardinerStone( complete_Flux[X][4][1][Udens] ,  MHDFlux[Y][1][0][X]-E_ref[1][1][2][Z],  MHDFlux[Y][0][0][X]-E_ref[1][1][1][Z])
                        + upwind_GardinerStone( complete_Flux[Y][1][0][Udens] , -MHDFlux[X][0][1][Y]-E_ref[1][1][1][Z], -MHDFlux[X][4][1][Y]-E_ref[1][0][1][Z])
                        + upwind_GardinerStone( complete_Flux[Y][0][0][Udens] , -MHDFlux[X][0][1][Y]-E_ref[1][1][2][Z], -MHDFlux[X][4][1][Y]-E_ref[1][0][2][Z]) );   
   
   E_edge[2] = 0.25 * ( - MHDFlux[Y][2][0][Z] - MHDFlux[Y][0][0][Z] + MHDFlux[Z][0][1][Y] + MHDFlux[Z][3][1][Y]
                        + upwind_GardinerStone( complete_Flux[Y][2][0][Udens] ,  MHDFlux[Z][0][1][Y]-E_ref[1][1][1][X],  MHDFlux[Z][3][1][Y]-E_ref[1][0][1][X])
                        + upwind_GardinerStone( complete_Flux[Y][0][0][Udens] ,  MHDFlux[Z][0][1][Y]-E_ref[2][1][1][X],  MHDFlux[Z][3][1][Y]-E_ref[2][0][1][X])
                        + upwind_GardinerStone( complete_Flux[Z][0][1][Udens] , -MHDFlux[Y][2][0][Z]-E_ref[2][0][1][X], -MHDFlux[Y][0][0][Z]-E_ref[1][0][1][X])
                        + upwind_GardinerStone( complete_Flux[Z][3][1][Udens] , -MHDFlux[Y][2][0][Z]-E_ref[2][1][1][X], -MHDFlux[Y][0][0][Z]-E_ref[1][1][1][X]) );      
   
   E_edge[3] = 0.25 * ( - MHDFlux[X][0][0][Y] - MHDFlux[X][4][0][Y] + MHDFlux[Y][0][0][X] + MHDFlux[Y][3][0][X]
                        + upwind_GardinerStone( complete_Flux[X][0][0][Udens] ,  MHDFlux[Y][0][0][X]-E_ref[1][0][1][Z],  MHDFlux[Y][3][0][X]-E_ref[1][0][0][Z])
                        + upwind_GardinerStone( complete_Flux[X][4][0][Udens] ,  MHDFlux[Y][0][0][X]-E_ref[1][1][1][Z],  MHDFlux[Y][3][0][X]-E_ref[1][1][0][Z])
                        + upwind_GardinerStone( complete_Flux[Y][0][0][Udens] , -MHDFlux[X][0][0][Y]-E_ref[1][1][0][Z], -MHDFlux[X][4][0][Y]-E_ref[1][0][0][Z])
                        + upwind_GardinerStone( complete_Flux[Y][3][0][Udens] , -MHDFlux[X][0][0][Y]-E_ref[1][1][1][Z], -MHDFlux[X][4][0][Y]-E_ref[1][0][1][Z]) );      
   
   E_edge[4] = 0.25 * ( + MHDFlux[X][0][1][Z] + MHDFlux[X][3][1][Z] - MHDFlux[Z][2][0][X] - MHDFlux[Z][0][0][X]
                        + upwind_GardinerStone( complete_Flux[X][0][1][Udens] , -MHDFlux[Z][2][0][X]-E_ref[0][1][2][Y], -MHDFlux[Z][0][0][X]-E_ref[0][1][1][Y])
                        + upwind_GardinerStone( complete_Flux[X][3][1][Udens] , -MHDFlux[Z][2][0][X]-E_ref[1][1][2][Y], -MHDFlux[Z][0][0][X]-E_ref[1][1][1][Y])
                        + upwind_GardinerStone( complete_Flux[Z][2][0][Udens] ,  MHDFlux[X][0][1][Z]-E_ref[1][1][1][Y],  MHDFlux[X][3][1][Z]-E_ref[0][1][1][Y])
                        + upwind_GardinerStone( complete_Flux[Z][0][0][Udens] ,  MHDFlux[X][0][1][Z]-E_ref[1][1][2][Y],  MHDFlux[X][3][1][Z]-E_ref[0][1][2][Y]) );      
   
   E_edge[5] = 0.25 * ( - MHDFlux[Y][0][1][Z] - MHDFlux[Y][4][1][Z] + MHDFlux[Z][1][0][Y] + MHDFlux[Z][0][0][Y]
                        + upwind_GardinerStone( complete_Flux[Y][0][1][Udens] ,  MHDFlux[Z][1][0][Y]-E_ref[0][2][1][X],  MHDFlux[Z][0][0][Y]-E_ref[0][1][1][X])
                        + upwind_GardinerStone( complete_Flux[Y][4][1][Udens] ,  MHDFlux[Z][1][0][Y]-E_ref[1][2][1][X],  MHDFlux[Z][0][0][Y]-E_ref[1][1][1][X])
                        + upwind_GardinerStone( complete_Flux[Z][1][0][Udens] , -MHDFlux[Y][0][1][Z]-E_ref[1][1][1][X], -MHDFlux[Y][4][1][Z]-E_ref[0][1][1][X])
                        + upwind_GardinerStone( complete_Flux[Z][0][0][Udens] , -MHDFlux[Y][0][1][Z]-E_ref[1][2][1][X], -MHDFlux[Y][4][1][Z]-E_ref[0][2][1][X]) );      
   
   E_edge[6] = 0.25 * ( + MHDFlux[X][0][0][Z] + MHDFlux[X][3][0][Z] - MHDFlux[Z][0][0][X] - MHDFlux[Z][4][0][X]
                        + upwind_GardinerStone( complete_Flux[X][0][0][Udens] , -MHDFlux[Z][0][0][X]-E_ref[0][1][1][Y], -MHDFlux[Z][4][0][X]-E_ref[0][1][0][Y])
                        + upwind_GardinerStone( complete_Flux[X][3][0][Udens] , -MHDFlux[Z][0][0][X]-E_ref[1][1][1][Y], -MHDFlux[Z][4][0][X]-E_ref[1][1][0][Y])
                        + upwind_GardinerStone( complete_Flux[Z][0][0][Udens] ,  MHDFlux[X][0][0][Z]-E_ref[1][1][0][Y],  MHDFlux[X][3][0][Z]-E_ref[0][1][0][Y])
                        + upwind_GardinerStone( complete_Flux[Z][4][0][Udens] ,  MHDFlux[X][0][0][Z]-E_ref[1][1][1][Y],  MHDFlux[X][3][0][Z]-E_ref[0][1][1][Y]) );      
   
   E_edge[7] = 0.25 * ( - MHDFlux[X][2][0][Y] - MHDFlux[X][0][0][Y] + MHDFlux[Y][0][1][X] + MHDFlux[Y][3][1][X]
                        + upwind_GardinerStone( complete_Flux[X][2][0][Udens] ,  MHDFlux[Y][0][1][X]-E_ref[1][1][1][Z],  MHDFlux[Y][3][1][X]-E_ref[1][1][0][Z])
                        + upwind_GardinerStone( complete_Flux[X][0][0][Udens] ,  MHDFlux[Y][0][1][X]-E_ref[1][2][1][Z],  MHDFlux[Y][3][1][X]-E_ref[1][2][0][Z])
                        + upwind_GardinerStone( complete_Flux[Y][0][1][Udens] , -MHDFlux[X][2][0][Y]-E_ref[1][2][0][Z], -MHDFlux[X][0][0][Y]-E_ref[1][1][0][Z])
                        + upwind_GardinerStone( complete_Flux[Y][3][1][Udens] , -MHDFlux[X][2][0][Y]-E_ref[1][2][1][Z], -MHDFlux[X][0][0][Y]-E_ref[1][1][1][Z]) );      
   
   E_edge[8] = 0.25 * ( + MHDFlux[X][1][0][Z] + MHDFlux[X][0][0][Z] - MHDFlux[Z][0][1][X] - MHDFlux[Z][4][1][X]
                        + upwind_GardinerStone( complete_Flux[X][1][0][Udens] , -MHDFlux[Z][0][1][X]-E_ref[1][1][1][Y], -MHDFlux[Z][4][1][X]-E_ref[1][1][0][Y])
                        + upwind_GardinerStone( complete_Flux[X][0][0][Udens] , -MHDFlux[Z][0][1][X]-E_ref[2][1][1][Y], -MHDFlux[Z][4][1][X]-E_ref[2][1][0][Y])
                        + upwind_GardinerStone( complete_Flux[Z][0][1][Udens] ,  MHDFlux[X][1][0][Z]-E_ref[2][1][0][Y],  MHDFlux[X][0][0][Z]-E_ref[1][1][0][Y])
                        + upwind_GardinerStone( complete_Flux[Z][4][1][Udens] ,  MHDFlux[X][1][0][Z]-E_ref[2][1][1][Y],  MHDFlux[X][0][0][Z]-E_ref[1][1][1][Y]) );   
}


//===============================================================================
// upwind_GardinerStone: 
//   this is a helper function for the condition in Gardiner/Stone 2007 eq. (28)
//===============================================================================

double upwind_GardinerStone (double v, double E_left, double E_right )
{
   if ( fabs(v) < MACHINE_ZERO ) return 0.5*(E_left+E_right);
   else if ( v > 0.) return E_left;
   else return E_right;   
}

//===============================================================================
// calc_dBdt: calculates dB/dt by finite-volume-discretization of the induction
//            equation. These are the equations (7) through (9) in Ziegler 2004.
// input:  E_edge[9], the nine edge-averaged E values needed
// output: dBdt[NDIM], the time-change of B[NDIM]
//===============================================================================

void calc_dBdt (double E_edge[9], double dBdt[NDIM])
{
   dBdt[X] = - ( E_edge[7] - E_edge[3] ) + ( E_edge[8] - E_edge[6] );
   dBdt[Y] = + ( E_edge[1] - E_edge[3] ) - ( E_edge[2] - E_edge[0] );
   dBdt[Z] = - ( E_edge[4] - E_edge[6] ) + ( E_edge[5] - E_edge[0] );
}

//===============================================================================
// calc_E: calculates E = - v x B
// input:  u[NHYDRO], B[NDIM]
// output: E[NDIM]
// note:   It's not important where the quantities are defined: at the interface,
//         cell-averaged or whatever. This basically just calculates the cross product. 
//===============================================================================

void calc_E (double u[NHYDRO], double B[NDIM], double E[NDIM])
{      
   E[X] = ( - u[UmomdensY] * B[Z] + u[UmomdensZ] * B[Y] ) / u[Udens];  // TODO: small densities???
   E[Y] = ( - u[UmomdensZ] * B[X] + u[UmomdensX] * B[Z] ) / u[Udens];
   E[Z] = ( - u[UmomdensX] * B[Y] + u[UmomdensY] * B[X] ) / u[Udens];
}

//===============================================================================
// getStencil
//
// input:   MHDnodes[5][5][5] array and a cell position within specified by x,y,z.
// output:  stencil[3] array containing pointers to the cell at the specified 
//            position ( mapped to stencil[1] ) and its two neighbour cells in
//            the idim direction ( mapped to stencil[0] and stencil[2] ).
// warning: it is up to the user not to violate any array boundaries!!!
//===============================================================================

void getStencil(nptr MHDnodes[5][5][5], int x, int y, int z, nptr stencil[3], int idim)
{
   int i;
   
   switch (idim)
   {
      case X:
         for (i=0;i<3;i++) stencil[i] = MHDnodes[z][y][x+i-1]; break;
      case Y:
         for (i=0;i<3;i++) stencil[i] = MHDnodes[z][y+i-1][x]; break;
      case Z:
         for (i=0;i<3;i++) stencil[i] = MHDnodes[z+i-1][y][x]; break;
      default:
         fprintf(stderr,"unexpected error in getStencil: invalid spatial direction given\n");
   }
}

//===============================================================================
// getStencil3D
//
// input:   MHDnodes[5][5][5] array and a cell position within specified by x,y,z.
// output:  stencil3D[3][3][3] array containing pointers to the cell at the specified
//            position (mapped to stencil3D[1][1][1]) and all its neighbour cells.
// warning: it is up to the user not to violate any array boundaries!!!
//===============================================================================

void getStencil3D(nptr MHDnodes[5][5][5], int x, int y, int z, nptr stencil3D[3][3][3])
{
   int i,j,k;

   for (k=0;k<3;k++) {
      for (j=0;j<3;j++) {
         for (i=0;i<3;i++) {
            stencil3D[k][j][i] = MHDnodes[z+k-1][y+j-1][x+i-1];
         }
      }
   }
}

//===============================================================================
// transformDRCtoXYZ:
//   gives the x,y,z cell coordinates inside the MHDnodes[5][5][5] array
//   for a given set of idim/row/cell coordinates (in function calc_fluxKNPCT),
//   where the possible x,y,z values are 1-3 and [2][2][2] is the central cell.
//
// input:  idim, row, cell and pointers xp, yp, zp to where x,y,z will be written
// output: writes the corresponding x,y,z values to the address at xp,yp,zp
// note:   for further information about these coordinates refer to documentation.
//===============================================================================

void transformDRCtoXYZ(int idim, int row, int cell, int* xp, int* yp, int* zp)
{
   int x,y,z;
   
   x = cell + 1;                                  // make things faster and unreadable by using integer operations :)
   y = ( row%2 || !row ) ? 2 : ( row%4 + 1 ) ;
   z = ( row%2 ) ? ( row^2 ) : 2 ;
   
   if (idim == Y)        { x=x^y; y=x^y; x=x^y; x=x^z; z=x^z; x=x^z; } // cyclic permutation
   else if (idim == Z)   { x=x^y; y=x^y; x=x^y; y=y^z; z=y^z; y=y^z; }
   
   *xp = x; *yp = y; *zp = z;    // write the calculated coordinates to the given addresses
}

//===============================================================================
// calc_MHDfluxKNP_infc
//   this is an exact copy of calc_fluxKNP_infc, with the only difference that
//   NHYDRO/NADVECT is substituted by NDIM to fit the correct number of B field components.
//   We can use this very same function, because the MHD flux is calculated with the
//   same KNP flux formula as the hydro quantities.
//===============================================================================

void calc_MHDfluxKNP_infc(double flux_infc[3][2][NDIM], double u_infc[3][2][NDIM], double a[2][2],
                          double Flux[2][NDIM])
{
   double da;
   int    idim;
   
   /* F^x_i-0.5 */
   da = a[0][1]-a[0][0];
   if(fabs(da) > MACHINE_ZERO)
   {
      for(idim=0; idim<NDIM; idim++)
      {
         Flux[0][idim] = (  a[0][1] * flux_infc[0][1][idim] 
                          - a[0][0] * flux_infc[1][0][idim]
                          + a[0][1] * a[0][0] * (u_infc[1][0][idim] - u_infc[0][1][idim])
                          ) / da;
      }      
   }
   else
   {
      for(idim=0; idim<NDIM; idim++)
      {
         Flux[0][idim] = 0.0;
#ifdef DEBUGMHD
         fprintf(stderr,"calc_MHDfluxKNP_infc: da < 0, setting MHDflux to zero\n");
#endif
      }
   }
   
   /* F^x_i+0.5 */
   da = a[1][1]-a[1][0];
   if(fabs(da) > MACHINE_ZERO)
   {
      for(idim=0; idim<NDIM; idim++)
      {
         Flux[1][idim] = (  a[1][1] * flux_infc[1][1][idim]
                          - a[1][0] * flux_infc[2][0][idim]
                          + a[1][1] * a[1][0] * (u_infc[2][0][idim] - u_infc[1][1][idim])
                          ) / da;
      }
   }
   else
   {
      for(idim=0; idim<NDIM; idim++)
      {
         Flux[1][idim] = 0.0;
#ifdef DEBUGMHD
         fprintf(stderr,"calc_MHDfluxKNP_infc: da < 0, setting MHDflux to zero\n");
#endif
      }
   }
}

/*==============================================================================
 * calculate a mean/maximum/whatever norm of div B on a given grid
 * currently: max(divB) - maximum value of divB/|B| reached in the given grid
 *==============================================================================*/
double calc_divB(gridls *cur_grid)
{
   pqptr   cur_pquad;
   cqptr   cur_cquad, icur_cquad;
   nqptr   cur_nquad, icur_nquad;
   nptr    cur_node;
   nptr    tsc_nodes[3][3][3];
   long    x, y, z;
   double  divB, B2, maxdivB;
   long    no_nodes;
   
   divB     = 0.;
   B2       = 0.;
   maxdivB  = 0.;
   no_nodes = 0;
   
   /*----------------------------------------------------------------------
    * loop over all nodes on supplied cur_grid
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
                     no_nodes++;
                     
                     /* check for interior nodes */
                     tsc_nodes[1][1][1] = cur_node;
                     get_TSCnodes(cur_grid, cur_pquad, icur_cquad, icur_nquad, tsc_nodes, &z, &y, &x);
                     
                     if(test_tsc(tsc_nodes) == TRUE)
                     {                        
                        divB =  ( tsc_nodes[1][1][2]->B[X] - cur_node->B[X] )
                              + ( tsc_nodes[1][2][1]->B[Y] - cur_node->B[Y] )
                              + ( tsc_nodes[2][1][1]->B[Z] - cur_node->B[Z] );
                        
                        B2 = 0.25 * (pow2(cur_node->B[X] + tsc_nodes[1][1][2]->B[X]) +
                                     pow2(cur_node->B[Y] + tsc_nodes[1][2][1]->B[Y]) +
                                     pow2(cur_node->B[Z] + tsc_nodes[2][1][1]->B[Z]) );
                        
                        divB = divB / sqrt(B2);
                        
                        maxdivB = MAX(divB,maxdivB);                        
                     }
                  }
               }
            }
         }
      }
   }
   
   return(maxdivB);
   
}


#endif // MHD //
