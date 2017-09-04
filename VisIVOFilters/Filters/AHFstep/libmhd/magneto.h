#include "../param.h"
#include "../tdef.h"

#ifndef MAGNETO_INCLUDED
#define MAGNETO_INCLUDED

#ifdef MHD
void     calc_FluxKNPCT       (nptr MHDnodes[5][5][5], double Flux[NDIM][2][NHYDRO], double dBdt[NDIM]);
void     calc_E_ref           (nptr MHDnodes[5][5][5], double E_ref[3][3][3][NDIM] );
void     calc_B_infc          (nptr stencil3D[3][3][3], double B_infc[2][NDIM], int idim);
void     calc_E_infc          (double u_infc[2][NHYDRO], double B_infc[2][NDIM], double E_infc[2][NDIM]);
void     MHD_flux_tensor      (double E_infc[NDIM], double E_tensor[NDIM], int idim);
void     calc_E_edge_Ziegler  (double MHDflux[NDIM][5][2][NDIM], double E_edge[9]);
void     calc_E_edge_GardinerStone (double complete_Flux[NDIM][5][2][NHYDRO], double E_ref[3][3][3][NDIM], double MHDFlux[NDIM][5][2][NDIM], double E_edge[9]);
double   upwind_GardinerStone (double v, double E_left, double E_right );
void     calc_dBdt            (double E_edge[9], double dBdt[NDIM]);
void     calc_E               (double u[NHYDRO], double B[NDIM], double E[NDIM]);
void     getStencil           (nptr MHDnodes[5][5][5], int x, int y, int z, nptr stencil[3], int idim);
void     getStencil3D         (nptr MHDnodes[5][5][5], int x, int y, int z, nptr stencil3D[3][3][3]);
void     transformDRCtoXYZ    (int idim, int row, int cell, int* xp, int* yp, int* zp);
void     calc_MHDfluxKNP_infc (double flux_infc[3][2][NDIM], double u_infc[3][2][NDIM], double a[2][2],
                               double Flux[2][NDIM]);
double   calc_divB            (gridls *cur_grid);
#endif
#endif