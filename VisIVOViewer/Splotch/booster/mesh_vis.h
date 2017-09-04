#ifndef TREEVIS
#define TREEVIS

#include "splotch/splotchutils.h"

struct Mesh_vis {

   float Cx,Cy,Cz;
   float posx,posy,posz;
   long cell_index;
   long num_particles;
   long offset;
   bool active;
   float weight;

};

struct Mesh_dim {

   float Lx,Ly,Lz;
   float Dx,Dy,Dz;
   long Nx,Ny,Nz;
   long ncell;

};


void mesh_creator(std::vector<particle_sim> &points, Mesh_vis ** Mesh, Mesh_dim * MeshD);

void p_selector(std::vector<particle_sim> &points, Mesh_vis * Mesh, Mesh_dim MeshD, std::vector<particle_sim> &r_points);

void randomizer(std::vector<particle_sim> &points, Mesh_vis * Mesh, Mesh_dim MeshD);

void m_rotation(paramfile &params, Mesh_vis ** p, Mesh_dim MeshD, const vec3 &campos, const vec3 &lookat, vec3 sky);

#endif
