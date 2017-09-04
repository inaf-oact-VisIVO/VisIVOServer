#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "cxxsupport/arr.h"
#include "cxxsupport/paramfile.h"
#include "cxxsupport/mpi_support.h"
#include "cxxsupport/bstream.h"
#include "splotch/splotchutils.h"
#include "booster/mesh_vis.h"

using namespace std;

void randomizer(vector<particle_sim> &points, Mesh_vis * Mesh, Mesh_dim MeshD)
{

   //long npart=points.size();
   srand ( time(NULL) );

   long i,j;
   long first;
   particle_sim swap_part;

   for (j=0; j<MeshD.ncell; j++)
   {
      first = Mesh[j].offset;
      for (i=0; i<Mesh[j].num_particles; i++)
      {

// random number generation

        double random_double = double(rand())/double(RAND_MAX); 
        long random_long = long(random_double*(Mesh[j].num_particles-i));
        long index = first+random_long;

// swap particles

        swap_part = points[first];
        points[first] = points[index];
        points[index] = swap_part; 

        first++;
        

      }
   }  


}
