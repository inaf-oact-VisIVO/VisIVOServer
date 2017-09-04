#if (NONC99 == 0)
#	include <math.h>
#else
#	include <replace_math.h>
#endif

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>


/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "../libutility/utility.h"
#include "../libgrids/grids.h"

static long unsigned totnodes, innodes;
static double        max_dens;
static double        cur_shift;

void dummy_write_divB()
{
}


#ifdef MHD

void write_divB(gridls *grid, char *prefix)
{
   char          filename[100];
   FILE          *divBfile;
   partptr       cur_part;
   pqptr         cur_pquad;
   cqptr         cur_cquad, icur_cquad;
   nqptr         cur_nquad, icur_nquad;
   nptr          cur_node;
   nptr          tsc_nodes[3][3][3];
   long          x, y, z;
   long          ncount;
   double        cur_shift;
   double        divB;
   
   /* shift of cell centre as compared to edge of box [grid units] */
   cur_shift = 0.5/(double)grid->l1dim;
   
   innodes  = 0;
   totnodes = 0;
   
   /* prepare filename */
   write_filename(filename, prefix, grid->l1dim);
   strcat(filename,"-divB");

   if((divBfile = fopen(filename,"w")) == NULL)
     {
      fprintf(stderr,"could not open %s\n", filename);
      exit(1);
     }
   
   for(cur_pquad=grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next)
     {
      z = cur_pquad->z;
      for(cur_cquad = cur_pquad->loc;
          cur_cquad < cur_pquad->loc + cur_pquad->length; 
          cur_cquad++, z++)  
        {  
         for(icur_cquad  = cur_cquad; 
             icur_cquad != NULL; 
             icur_cquad  = icur_cquad->next)
           {
            y = icur_cquad->y;
            for(cur_nquad = icur_cquad->loc;  
                cur_nquad < icur_cquad->loc + icur_cquad->length; 
                cur_nquad++, y++) 
              { 
               for(icur_nquad  = cur_nquad; 
                   icur_nquad != NULL; 
                   icur_nquad  = icur_nquad->next)
                 {
                  x = icur_nquad->x;
                  for(cur_node = icur_nquad->loc; 
                      cur_node < icur_nquad->loc + icur_nquad->length; 
                      cur_node++, x++, ncount++)
                    {
                       totnodes++;
                     
                     /* check for interior nodes */
                     tsc_nodes[1][1][1] = cur_node;
                     get_TSCnodes(grid, cur_pquad, icur_cquad, icur_nquad, tsc_nodes, &z, &y, &x);
                       
                     if(test_tsc(tsc_nodes) == TRUE)
                     {
                     
                        innodes++;
                        
                        divB =  ( tsc_nodes[1][1][2]->B[X] - cur_node->B[X] )
                              + ( tsc_nodes[1][2][1]->B[Y] - cur_node->B[Y] )
                              + ( tsc_nodes[2][1][1]->B[Z] - cur_node->B[Z] );
                        
                        fprintf(divBfile,"%f %f\n", ((float)(x)/(float)(grid->l1dim) + cur_shift), divB); 
                        fflush(divBfile);


                     }
                     
                     
                    }
                 }
              }
           }
        }
     }

   fclose(divBfile);
}

#endif // MHD
