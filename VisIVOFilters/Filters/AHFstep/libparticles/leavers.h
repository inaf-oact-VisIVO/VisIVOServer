#include "../param.h"
#include "../tdef.h"

#ifndef LEAVERS_INCLUDED
#define LEAVERS_INCLUDED

void    store_leaver     (partptr cur_part, gridls *cur_grid);
void    unstep_leavers   (gridls *cur_grid, double DRIFT1, double DRIFT2, double KICK);
void    advance_leavers  (gridls *cur_grid, double DRIFT, double KICK);

#endif

