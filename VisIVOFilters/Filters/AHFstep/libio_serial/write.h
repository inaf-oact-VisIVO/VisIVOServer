#include "../param.h"
#include "../tdef.h"

#ifndef WRITE_INCLUDED
#define WRITE_INCLUDED

/* these routines serve mainly debugging purposes and hence write information in ASCII */
void write_boundary         (gridls *cur_grid);
void write_density          (gridls *grid_ptr, char *prefix);
void write_divB             (gridls *grid_ptr, char *prefix);
void write_forces           (gridls *grid_ptr);
void write_hydro            (gridls *grid);
void write_nodepart         (gridls *grid);
void write_nodes            (gridls *grid, char *prefix);
void write_positions        (gridls *grid);
void write_residual         (gridls *grid_ptr, char *prefix);

#endif
