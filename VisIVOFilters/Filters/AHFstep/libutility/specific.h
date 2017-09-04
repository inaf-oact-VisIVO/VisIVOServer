#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "../common.h"
#include "../param.h"
#include "../tdef.h"

#ifndef SPECIFIC_INCLUDED
#define SPECIFIC_INCLUDED

/* AMIGA specific routines */
MINMAX  MinMax               (double,double,double);
MINMAX  MinMaxBound          (double,double,double,double);
void    read_amiga_header    (FILE *infile, info_io *io, int *SWAPBYTES);
void    get_user_data        (FILE *startrun_dat, uparamptr user);
void    sanity_check         ();
void    ic_unit_conversion   ();
double  init_header_masses   ();
void    init_io_structure    (uparam user);
void    init_simu_structure  (uparam user);
void    init_global_structure(uparam user);
void    init_logfile         (uparam user);
void    binning_parameter    (HALO halo, int *nbins, double *dist_min, double *dist_max);
void    write_filename       (char *f_name, char *prefix, unsigned l1dim);
void    get_c2fslope         (double func[3][3][3], double slope[3]);
void    get_axes             (double itensor[3][3], double *axis1, double *axis2, double *axis3);
int     idx_inv              (long unsigned *idx, int numHalos, int j);
double  f1mod                (double x, double y);
void    get_DKoperators      (double timestep, double *KICK, double *DRIFT1, double *DRIFT2);
double  init_timestep        (uparam user);
double  adjust_timestep      (double timecounter, double timestep);
double  Laplace_pot          (nptr tsc_nodes[3][3][3], double spacing);
double  Laplace_temp1        (nptr tsc_nodes[3][3][3], double spacing);

#if (defined NEWSTARTRUN)
extern int
cmp_sfckey_part(const void *p1, const void *p2);
#endif

#endif
