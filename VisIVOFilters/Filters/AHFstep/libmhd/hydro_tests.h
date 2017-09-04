#include "../param.h"
#include "../tdef.h"

#ifndef HYDRO_TESTS_INCLUDED
#define HYDRO_TESTS_INCLUDED

#ifdef HYDRO

void hydro_test (gridls  *grid_list);
void zero_udens(gridls *cur_grid);
void init_udens(gridls *cur_grid);

#endif
#endif

