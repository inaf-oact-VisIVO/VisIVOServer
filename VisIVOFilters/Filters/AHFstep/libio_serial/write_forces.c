#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "../libutility/utility.h"

static double cur_shift, f_fac;

int count = 1;  // Claudio.

/*--------------------------------------------------------------------------*/
/* write_forces_nquad: nquad control routine */
/*--------------------------------------------------------------------------*/
void write_forces_nquad(gridls *cur_grid, nqptr cur_nquad, 
                        FILE *forcefile, unsigned l1dim,
                        unsigned x, unsigned y, unsigned z)
{
   nptr cur_node;          /* current node */
   
   
   /* loop over nodes */
   for(cur_node = cur_nquad->loc; cur_node < cur_nquad->loc + cur_nquad->length; 
       cur_node++, x++)
     {
      fprintf(forcefile,"%f %f %f %f %f %f  %d\n",
              (float)(x)/(float)(l1dim) + cur_shift,
              (float)(y)/(float)(l1dim) + cur_shift,
              (float)(z)/(float)(l1dim) + cur_shift,
              f_fac*cur_node->force.forces[X],
              f_fac*cur_node->force.forces[Y],
              f_fac*cur_node->force.forces[Z], count);
       ++count;
     }
   fflush(forcefile);
  
}

/*--------------------------------------------------------------------------*/
/* write_forces_cquad: cquad control routine */
/*--------------------------------------------------------------------------*/
void write_forces_cquad(gridls *cur_grid, cqptr cur_cquad, 
                        FILE *forcefile, unsigned l1dim, 
                        unsigned y, unsigned z)
{
   nqptr cur_nquad;       /* current nquad */
   nqptr icur_nquad;      /* current nquad */
   unsigned x;
   
   
   /* loop over nquads */
   for(cur_nquad = cur_cquad->loc; cur_nquad < cur_cquad->loc + cur_cquad->length; 
       cur_nquad++, y++)
     {
      for(icur_nquad = cur_nquad; icur_nquad != NULL; 
          icur_nquad = icur_nquad->next)
        {
         x = icur_nquad->x;
         write_forces_nquad(cur_grid, icur_nquad, forcefile, l1dim, x, y, z);
        }
     }
}

/*--------------------------------------------------------------------------*/
/* write_forces_pquad: pquad control routine */
/*--------------------------------------------------------------------------*/
void write_forces_pquad(gridls *cur_grid, pqptr cur_pquad, 
                        FILE *forcefile, unsigned l1dim, 
                        unsigned z)
{
   cqptr cur_cquad;      /* current nquad */
   cqptr icur_cquad;     /* cquad in linked cquad loop */
   unsigned y;
   
   /* loop over cquads */
   for(cur_cquad = cur_pquad->loc; cur_cquad < cur_pquad->loc + cur_pquad->length; 
       cur_cquad++, z++)
     {
      for(icur_cquad = cur_cquad; icur_cquad != NULL; 
          icur_cquad = icur_cquad->next)
        {
         y = icur_cquad->y;
          write_forces_cquad(cur_grid, icur_cquad, forcefile, l1dim, y, z);
        }
     }
}

/*--------------------------------------------------------------------------*/
/*-------- write_forces: loop over grid writing out density ------------------*/
/*--------------------------------------------------------------------------*/
void write_forces(gridls *cur_grid)
{
   pqptr cur_pquad;     /* pointer to the pquad */
   unsigned l1dim, z;
   
   char filename[100];
   FILE *forcefile;
   
   /* shift of cell centre as compared to edge of box [grid units] */
   cur_shift = 0.5/(double)cur_grid->l1dim;

   l1dim = cur_grid->l1dim;
   f_fac = simu.boxsize;
   
   /* prepare filename */
   write_filename(filename, "Forces.", (unsigned int)global.no_timestep);
   
   if((forcefile = fopen(filename,"w")) == NULL)
     {
      fprintf(stderr,"could not open %s\n", filename);
      exit(1);
     }
   
   /* pass cur_pquad to ass_pquad */
   for(cur_pquad = cur_grid->pquad; cur_pquad != NULL; cur_pquad = cur_pquad->next)
     {
      z = cur_pquad->z;
      write_forces_pquad(cur_grid, cur_pquad, forcefile, l1dim, z);
     }
   fclose(forcefile);
}


