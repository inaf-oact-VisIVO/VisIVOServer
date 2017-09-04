#ifndef GEN_CUBES_INCLUDED
#define GEN_CUBES_INCLUDED

#ifdef LIGHTCONE
void create_cone(X1,Y1,Z1,XYZ2,rsq1,rsq2);
void create_slice(N,rsq1,rsq2);
int edge_check();
void euler(x,y,z,px,py,pz);
int face_check();
int radius_check(x,y,z,rsq);
int vertex_check(x,y,z);
#endif

#endif

