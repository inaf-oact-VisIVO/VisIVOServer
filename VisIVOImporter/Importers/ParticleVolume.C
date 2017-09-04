#include "hdf5.h"
#include <stdio.h>
#include <iostream>
#include <string>

using namespace std;

main() {

   FILE * hdf5File;
   char hierarchyname[1000];
   char datafilename[1000];
   char outputfilename[1000];
   int ngx[3];
   long np,ng;
   int nrank = 1;
   int nghost = 3;
   long npart;
   long npartout;
   float leftside[nrank];
   float rightside[nrank];
   int nleft[nrank];
   int nright[nrank];
   int nsize[nrank];
   const int numberoffields = 6;

   float * dataarray;
   float * outarray;
   
   hsize_t * s_dims   = new hsize_t [nrank];
   hsize_t * s_maxdims   = new hsize_t [nrank];
   string fieldsnames [numberoffields];

   fieldsnames[0] = "particle_position_x";
   fieldsnames[1] = "particle_position_y";
   fieldsnames[2] = "particle_position_z";
   fieldsnames[3] = "particle_velocity_x";
   fieldsnames[4] = "particle_velocity_y";
   fieldsnames[5] = "particle_velocity_z";

   printf ("Input file name: ");
   scanf  ("%s", hierarchyname); 
   printf ("Output file name: ");
   scanf  ("%s", outputfilename); 

   float boxmin[3];
   float boxmax[3];

   printf ("Input xmin, ymin, zmin\n");
   scanf  ("%f", &boxmin[0]);
   scanf  ("%f", &boxmin[1]);
   scanf  ("%f", &boxmin[2]);
   printf ("Input xmax, ymax, zmax\n");
   scanf  ("%f", &boxmax[0]);
   scanf  ("%f", &boxmax[1]);
   scanf  ("%f", &boxmax[2]);

// create DESTINATION HDF5 file

   hid_t file_id;

   FILE * pFile1;

   pFile1 = fopen64 (outputfilename,"ab");
   if (pFile1==NULL)
   {
       fclose(pFile1);
       file_id = H5Fcreate(outputfilename, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
   } else {
       fclose(pFile1);
       file_id = H5Fopen(outputfilename, H5F_ACC_RDWR, H5P_DEFAULT);
   }

   hsize_t dims[nrank];
   hsize_t maxdims[nrank];

   printf("--> %ld\n", dims[0]);

// create all DATASETS

   hid_t obj_id;
   hid_t dataspace;
 
// open source file 

// UGO: ti puo' interessare da qua... vedi commenti "UGO" successivi 
        printf("Reading Data from: %s\n", hierarchyname);
// UGO apre il file HDF5:
	hid_t source_id = H5Fopen(hierarchyname, H5F_ACC_RDONLY, H5P_DEFAULT);

// read datasets from SOURCE

	for(int k=0; k<numberoffields; k++)
	{


// UGO apre un dataset presente nel file il cui nome e' nella variabile fieldsnames[k] (vedi sopra):
	    hid_t source_obj = H5Dopen(source_id,fieldsnames[k].c_str());
// UGO recupera la "geometria" del dataset appena aperto:
	    hid_t source_space = H5Dget_space(source_obj);
// UGO recupera la dimensionalita' del dataset (1, 2, 3...) e la mette in source_space
// UGO recupera il numero di elementi (particelle, celle...) per ogni dimensione e li mette in s_dims (dimentica s_maxdims)
   	    H5Sget_simple_extent_dims(source_space, s_dims, s_maxdims);
	    printf("Block particles = %d\n", s_dims[0]);
	    
// allocate memory for loading data

// UGO alloca l'array in base alla dimensione letta (in questo caso assumiamo un array 1D)
	    dataarray = new float[s_dims[0]];

// read data
// UGO legge i dati
	    H5Dread(source_obj, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataarray);
//	    printf("%f\n", dataarray[datasizeread-1]);

   	    npart = s_dims[0];

// allocate destination data array
            if(k==0)outarray = new float[6*npart];

// write data in proper positions
            for(long ii=0; ii<npart; ii++)outarray[(k*npart)+ii]=dataarray[ii];

// UGO chiude le cose aperte in precedenza (richiesto da HDF5)
	    H5Sclose (source_space);
	    H5Dclose (source_obj);
	    delete [] dataarray;
	}
        H5Fclose (source_id);

// check which data falls inside the box

        float * x_eff = new float[npart];
        float * y_eff = new float[npart];
        float * z_eff = new float[npart];
        float * vx_eff = new float[npart];
        float * vy_eff = new float[npart];
        float * vz_eff = new float[npart];
      
        long part_eff=0;
        for(long ii=0; ii<npart; ii++) {
           int iadd=0;
           for(int jj=0; jj<3; jj++) 
           if(outarray[(jj*npart)+ii] >= boxmin[jj] && outarray[(jj*npart)+ii] <= boxmax[jj])iadd++;
           if(iadd==3) {
                  x_eff[part_eff]=outarray[ii];
                  y_eff[part_eff]=outarray[npart+ii];
                  z_eff[part_eff]=outarray[2*npart+ii];
                  vx_eff[part_eff]=outarray[3*npart+ii];
                  vy_eff[part_eff]=outarray[4*npart+ii];
                  vz_eff[part_eff]=outarray[5*npart+ii];
                  part_eff++; 
           }
        }
        
// create output datasets


// UGO: fase di scrittura:
        printf("EFFECTIVE SIZE ----> %ld\n", part_eff-1);
        dims[0] = part_eff-1;
        maxdims[0] = part_eff-1;
// UGO: crea la geometria del dato da scrivere nel file
        dataspace = H5Screate_simple (nrank, dims, maxdims);
        for (int j=0; j<numberoffields; j++)
        {  
// UGO: crea i dataset di nome fieldsnames e geometria definita da dataspace
             obj_id = H5Dcreate(file_id, fieldsnames[j].c_str(), H5T_NATIVE_FLOAT, dataspace, H5P_DEFAULT);
             H5Dclose(obj_id);
        }
        H5Sclose (dataspace);

// write dataset on destination
	    
        float * point_data;

        for(int k=0; k<numberoffields; k++)
        {
           switch(k) {
             case 0:
                 point_data = x_eff;
                 break;
             case 1:
                 point_data = y_eff;
                 break;
             case 2:
                 point_data = z_eff;
                 break;
             case 3:
                 point_data = vx_eff;
                 break;
             case 4:
                 point_data = vy_eff;
                 break;
             case 5:
                 point_data = vz_eff;
                 break;
             }


	   obj_id  = H5Dopen(file_id, fieldsnames[k].c_str());
           printf("-------------> WRITING %s\n",fieldsnames[k].c_str());

	   dataspace = H5Dget_space(obj_id);
// UGO: scrivi i dataset nel file
	   H5Dwrite (obj_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, point_data);

	   H5Dclose (obj_id);
           H5Sclose (dataspace);
	    
	}
        
        H5Fclose (file_id);
}
