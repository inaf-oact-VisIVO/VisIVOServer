/* 
 *   Writing and reading an existing dataset.
 */

#include <hdf5.h>
#include <stdlib.h>
#include <stdio.h>

main() {

   hid_t       file_id, dataset_id;  /* identifiers */
   herr_t      status;
   int         i, j; 
   int         npoints;
   int         nsize;
   float       *xx;
   float       *yy;
   float       *zz;
   int         iii;
   FILE        *outfile_id;
   char        file_in[80];
   char        file_out[80];
   char        dataset_name[80];

   printf("Input 1D Size\n");
   scanf ("%d",&nsize);
   printf("nsize=%d\n",nsize);

   nsize=128;   
   npoints = nsize*nsize*nsize;
   xx = (float *) malloc(npoints*sizeof(float));
   printf("x allocated\n");
   yy = (float *) malloc(npoints*sizeof(float));
   printf("y allocated\n");
   zz = (float *) malloc(npoints*sizeof(float));
   printf("z allocated\n");

   printf("Input data filename\n");
   scanf ("%s",file_in);

   printf("Input output filename\n");
   scanf ("%s",file_out);

   /* Open an existing file. */
   file_id = H5Fopen(file_in, H5F_ACC_RDWR, H5P_DEFAULT);
   printf("FILE OPENED\n");
   printf("FILEID %d \n", file_id);

   /* Open an existing dataset. */
   dataset_id = H5Dopen(file_id, "/particle_position_x"); 
   status = H5Dread(dataset_id, H5T_IEEE_F32BE, H5S_ALL, H5S_ALL, H5P_DEFAULT,xx); 
   status = H5Dclose(dataset_id);
   dataset_id = H5Dopen(file_id, "/particle_position_y"); 
   status = H5Dread(dataset_id, H5T_IEEE_F32BE, H5S_ALL, H5S_ALL, H5P_DEFAULT,yy); 
   status = H5Dclose(dataset_id);
   dataset_id = H5Dopen(file_id, "/particle_position_z"); 
   status = H5Dread(dataset_id, H5T_IEEE_F32BE, H5S_ALL, H5S_ALL, H5P_DEFAULT,zz); 
   status = H5Dclose(dataset_id);

   printf("PARTICLES READ\n");

/*
   for (iii=0; iii<256; iii++) 
	{
	printf("%d   %f \n", iii,dset_data[iii]);
	}
*/

   /* Close the file. */
   status = H5Fclose(file_id);

   outfile_id = fopen(file_out, "w");
   for(i=0; i<npoints; i++)
   {
     fwrite(&xx[i], 4, 1, outfile_id);
     fwrite(&yy[i], 4, 1, outfile_id);
     fwrite(&zz[i], 4, 1, outfile_id);
   }
   fclose(outfile_id);
   free(xx);
   free(yy);
   free(zz);
   
}

