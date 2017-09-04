#include "../param.h"
#include "../tdef.h"

#ifndef OUTPUT_LC_INCLUDED
#define OUTPUT_LC_INCLUDED

#ifdef LIGHTCONE
void output_lc();
void store_backup();
void create_cone();
void create_slice();
int radius_check();
int vertex_check();
int face_check();
int edge_check();
void euler();
#endif /*LIGHTCONE*/

#endif
