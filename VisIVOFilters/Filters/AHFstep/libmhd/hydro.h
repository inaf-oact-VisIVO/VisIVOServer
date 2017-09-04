#include "../param.h"
#include "../tdef.h"

#ifndef HYDRO_INCLUDED
#define HYDRO_INCLUDED

#ifdef HYDRO

void     hydro_flux                 (double u[NHYDRO], double B[NDIM], double flux[NHYDRO], int idim);
void     hydro_source               (double super_t, double timestep, nptr cur_node, double Source[NHYDRO], double BSource[NDIM], double B2);
void     init_DMhydro               ();
void     init_B                     (gridls* cur_grid);
void     solve_dom_hydro            (double timestep);
void     calc_Flux                  (gridls *cur_grid);
void     calc_FluxKNP               (nptr MHDnodes[5][5][5], double Flux[NDIM][2][NHYDRO]);
void     calc_fluxKNP_infc          (double flux_infc[3][2][NHYDRO], double u_infc[3][2][NHYDRO], double a[2][2], double Flux[2][NHYDRO]);
void     advect_hydro               (gridls *cur_grid, double timestep);
void     advect_gravity_hydro       (gridls *cur_grid, double timecounter, double timestep);
void     ave_hydro                  (gridls *cur_grid);
void     getCFLspeed_and_syncES     (gridls *cur_grid);
void     update_edens               (gridls *cur_grid);
double   calc_T                     (double e);
double   calc_e                     (double T);
double   calc_p                     (double u[NHYDRO]);
double   calc_cf                    (double u[NHYDRO], double B2);
void     min_max_speeds             (double u_infc[3][2][NHYDRO], double B_infc[3][2][NDIM], double a[2][2], int Umomdens);
void     adjust_globalCFLspeed      (double u[NHYDRO], double B2);
void     calc_u_infc                (nptr stencil[3], double u_infc[2][NHYDRO]);
double   du_ltd                     (double u_node[3]);
double   minmod                     (double a, double b);
double   vanLeer                    (double a, double b);
double   superbee                   (double a, double b);
void     zero_B_infc                (double B_infc[3][2][NDIM]);

#endif
#endif

