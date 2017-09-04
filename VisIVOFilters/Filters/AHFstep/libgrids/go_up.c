/* the important definitions have to be included first */
#include "../param.h"
#include "../tdef.h"

/*==============================================================================
 * go up one grid in hierachy and 
 * use recently interpolated source values from fine to coarse grid
 * (assumes trunc_err calculation previously)
 *==============================================================================*/
gridls *go_up(gridls *cur_grid)
{
  /* 
   * there's nothing to be done...
   * trun_err() already did the fine to coarse interpolation of the source term
   */

  return(cur_grid-1);
}

