#ifdef USE_MPI
#include "mpi.h"
#endif
#include <cstdio>
#include <iostream>
#include <string>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <vector>
#include "hdf5.h"

using namespace std;

/*
long Enzo_reader (vector<float> &xpos, vector<float> &ypos, vector<float> &zpos,
                  vector<float> &scalar, vector<float> &smooth, float *maxr, float *minr)
*/
long Enzo_reader (vector<float> * xpos, vector<float> * ypos, vector<float> * zpos,
                  vector<float> * scalar, vector<float> * smooth, float *maxr, float *minr)
{

   FILE * pFile;
   FILE * auxFile;
   FILE * hdf5File;
   char hierarchyname[1000];
   char datafilename[1000];
   char outputfilename[1000];
   string groupprefix;
   string hgroup;
   string completename;
   long np,ng;
   int nrank = 3;
   int nghost = 3;
   int nfiles;
   double leftside[nrank];
   double rightside[nrank];
   int lbox[nrank];
   int rbox[nrank];
   double lxbox[nrank];
   double rxbox[nrank];
   int sbox[nrank];
   int nleft[nrank];
   int nright[nrank];
   int nsize[nrank];
   const int numberoffields = 8;
   int gridid;
   char groupid[100];
   long gcounter=0;
   float sigma;

   int ngx[nrank];
   int vngx[nrank];
   double dxvngx[nrank];
   double rdxvngx[nrank];
   double boxsize[nrank];
   double hboxsize[nrank];
   double boxcenter[nrank];
   double basesize[nrank];
   double lintersect[nrank];
   double rintersect[nrank];
   double leftcorner[nrank];
   double rightcorner[nrank];
   int lnintersect[nrank];
   int rnintersect[nrank];
   int lbintersect[nrank];
   int rbintersect[nrank];
   basesize[0] = 1.0;
   basesize[1] = 1.0;
   basesize[2] = 1.0;
   int maxlevel;
   int boxcells[nrank];

   float * dataarray;
   float * destarray;
   
   long total_size=0;
   long total_size_old=0;
   float minradius=1e30;
   float maxradius=-1e30;

   string fieldsnames [numberoffields];

   int mype=0;
   int npes=1;

#ifdef USEMPI
   MPI_Comm_rank(MPI_COMM_WORLD, &mype);
   MPI_Comm_size(MPI_COMM_WORLD, &npes);
#endif

// initialize cutout parameters


   int sf[numberoffields];
   string selfield[numberoffields];

   groupprefix = "/Grid";

   for(int i=0; i<numberoffields; i++)sf[i]=-1;

   fieldsnames[0] = "/Dark_Matter_Density";
   fieldsnames[1] = "/Density";
   fieldsnames[2] = "/Temperature";
   fieldsnames[3] = "/x-velocity";
   fieldsnames[4] = "/y-velocity";
   fieldsnames[5] = "/z-velocity";
   fieldsnames[6] = "/Gas_Energy";
   fieldsnames[7] = "/Total_Energy";
   int number_of_fields2read;

// begin of mype=0 read section

   if(mype == 0)
   {
   printf("Available fields:\n");
   for(int i=0; i<numberoffields; i++)
	printf("%d, %s\n", i, fieldsnames[i].c_str());

   char cfieldsaux[100];
   string cfields;
   printf("Input the list of required fields (e.g.: 0,3,7)\n");
   scanf("%s", cfieldsaux);
   cfields = cfieldsaux;

   int st_size = cfields.length();

   for(int i=0, j=0; i<st_size; i+=2, j++) 
   {
      selfield[j] = cfields.substr(i,1);
      sf[j] = atoi(selfield[j].c_str());
      printf("SELECTED FIELD %d %s\n", sf[j], fieldsnames[sf[j]].c_str());
      number_of_fields2read = j;
   }

   number_of_fields2read++;
   printf("NUMBER OF SELECTED FIELD %d\n", number_of_fields2read);


   printf ("Input Reduced Hierarchy file name: \n");
   scanf  ("%s", hierarchyname); 
   printf ("Input output file name: ");
   scanf  ("%s", outputfilename); 

// this is the INTEGER number of cells on the level 0 mesh

   printf ("Input number of cells in the base level (int) --> x = \n");
   scanf  ("%d", &ngx[0]); 
   printf ("Input number of cells in the base level (int) --> y = \n");
   scanf  ("%d", &ngx[1]); 
   printf ("Input number of cells in the base level (int) --> z = \n");
   scanf  ("%d", &ngx[2]); 

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// WARNING X and Z are swapped!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

// this is the position of the left corner of the cutout box in units 0-1

   printf ("Input first sub-box position (0-1) --> x = \n");
   scanf  ("%lf", &leftcorner[2]); 
   printf ("Input first sub-box position (0-1) --> y = \n");
   scanf  ("%lf", &leftcorner[1]); 
   printf ("Input first sub-box position (0-1) --> z = \n");
   scanf  ("%lf", &leftcorner[0]); 

// this is the position of the right corner of the cutout box in units 0-1

   printf ("Input second sub-box position (0-1) --> x = \n");
   scanf  ("%lf", &rightcorner[2]); 
   printf ("Input second sub-box position (0-1) --> y = \n");
   scanf  ("%lf", &rightcorner[1]); 
   printf ("Input second sub-box position (0-1) --> z = \n");
   scanf  ("%lf", &rightcorner[0]); 

   printf ("Input Max refinement level (0=Base level) MUST BE DEEPEST LEVEL:\n");
   scanf  ("%d", &maxlevel);
// end of mype=0 read section
   }
#ifdef USEMPI
   MPI_Bcast(&number_of_fields2read, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(sf, number_of_fields2read, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&hierarchyname[0], 1000, MPI_CHAR, 0, MPI_COMM_WORLD);
   MPI_Bcast(ngx, 3, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(leftcorner, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Bcast(rightcorner, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Bcast(&maxlevel, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif


// level 0 resolution:

   double dxbase[nrank];
   for (int j=0; j<nrank; j++)dxbase[j] = 1.0/float(ngx[j]);

// half box size

   float fmaxlevel = (float)maxlevel;
   for (int i=0; i<nrank; i++) {
	boxsize[i] = rightcorner[i] - leftcorner[i];
	hboxsize[i] = 0.5 * boxsize[i];
   }

// calculate virtual global grid size and resolution at the output level

   for (int i=0; i<nrank; i++) 
   {
       vngx[i] = ngx[i] * (int) pow(2.0f, fmaxlevel); 
       dxvngx[i] = 1.0 / (double)vngx[i];
       rdxvngx[i] = 1.0 / dxvngx[i];
#ifdef DEBUG   
       printf("Virtual Size: %d %d\n", i, vngx[i]);
#endif
   }

//

   double boxcoarseL[3];
   double boxcoarseR[3];
   int ndims_aux[3];

   for (int i=0; i<nrank; i++)
   {
       boxcoarseL[i] = leftcorner[i];
       boxcoarseR[i] = rightcorner[i];
       ndims_aux[i] = ((int)(boxsize[i]*(double)ngx[i]))*(int)pow(2.0f, fmaxlevel);
       //ndims_aux[i] = ((int)(boxsize[i]*(double)ngx[i])+1)*(int)pow(2.0f, fmaxlevel);
   }


// esitmate the size of the output box at the output resolution
// notice that the whole box has size basesize=1.0

   int boxstart_cell[nrank];
   int boxend_cell[nrank];
   for (int i=0; i<nrank; i++)
   {

// box references in output resolution cells

        boxstart_cell[i]  = (int)((boxcoarseL[i]/basesize[i]) * vngx[i]);
	if(boxstart_cell[i] < 0) boxstart_cell[i] = 0;
        boxend_cell[i]    = (int)((boxcoarseR[i]/basesize[i]) * vngx[i])-1;
	if(boxend_cell[i] >= vngx[i]) boxend_cell[i] = vngx[i]-1;
        boxcells[i] = boxend_cell[i]-boxstart_cell[i]+1;

// box references in universe (0-1) units

        lxbox[i] = boxcoarseL[i];
        if(lxbox[i] < 0.0) lxbox[i] = 0.0;
        rxbox[i] = boxcoarseR[i];
        if(rxbox[i] > basesize[i]) rxbox[i] = basesize[i];

#ifdef DEBUG   
        printf("Cells in the BOX: %d %d\n", i, boxcells[i]);
        printf("BOX starts at: %d %d\n", i, boxstart_cell[i]);
        printf("BOX ends at: %d %d\n", i, boxend_cell[i]);
#endif
   }  

// calculate total size and prepare HDF5 quantities and output file

   ng = 1;
   for (int j=0; j<nrank; j++) ng*=boxcells[j];

   
// create all DATASETS

   hid_t obj_id;
   hid_t dataspace;
   hid_t memoryspace;
#ifdef SP5
   long long * start = new long long [nrank];
#endif
#ifdef BCX
   hsize_t * start = new hsize_t [nrank];
#endif
#ifdef CLX
   hssize_t * start = new hssize_t [nrank];
#endif
   hsize_t * stride = new hsize_t [nrank];
   hsize_t * count  = new hsize_t [nrank];
   hsize_t * block  = new hsize_t [nrank];
   long sourcesize;
   long destinationsize;
   hid_t source_id;
   hid_t source_obj;
   hid_t source_space;

#ifdef SP5
   long long * s_start = new long long [nrank];
#endif
#ifdef BCX
   hsize_t * s_start = new hsize_t [nrank];
#endif
#ifdef CLX
   hssize_t * s_start = new hssize_t [nrank];
#endif
   hsize_t * s_stride = new hsize_t [nrank];
   hsize_t * s_count  = new hsize_t [nrank];
   hsize_t * s_block  = new hsize_t [nrank];
   hsize_t * s_dims   = new hsize_t [nrank];
   hsize_t * s_maxdims   = new hsize_t [nrank];

// create all the datasets in the HDF5 output file

// build the datasets

   pFile = fopen (hierarchyname, "r");
   fscanf (pFile, "%d", &nfiles);


   int acont=0;
   long countt=0;
   double dxxx[nrank];
   double rdxxx[nrank];
   double levelbox;
   int nlevelbox;
   int checking = 0;

   int chunk_pe = int(nfiles/npes);
   int istart_pe = mype*chunk_pe;
   int iend_pe   = istart_pe+chunk_pe;
   if(mype == npes-1 && npes > 1)iend_pe=iend_pe+(int)nfiles%npes; 

#ifdef DEBUG
   printf("PE %d STARTS FROM %d END AT %d\n", mype, istart_pe, iend_pe);
#endif

//// MAIN LOOP

   for (int i=0; i<nfiles; i++)
//   for (int i=0; i<1; i++)
     {
 
        acont++; 
        printf("Merging block %d\n",i);
        checking = 0;

        char buffer[100];
        int iaux = 10000001+i;
        sprintf(buffer,"%s%d", groupprefix.c_str(),iaux);
        hgroup = buffer;
        string aux0="0";
        hgroup.replace(5,1,aux0);
        

// read the i-th slab

        int grididaux;

        // NOTICE that due to hierarchy file box inverse mapping z and x are swapped:
        // 0 --> 2
        // 1 --> 1
        // 2 --> 0

        //fscanf(pFile, "%d", &gridid);
        for (int j=0; j<nrank; j++) fscanf(pFile, "%d", &lbox[2-j]);
        for (int j=0; j<nrank; j++) fscanf(pFile, "%d", &rbox[2-j]);
	for (int j=0; j<nrank; j++) fscanf(pFile, "%lf", &leftside[2-j]);
	for (int j=0; j<nrank; j++) fscanf(pFile, "%lf", &rightside[2-j]);
        fscanf (pFile, "%s", datafilename);

        for (int j=0; j<nrank; j++)
	{
            lbox[j] -= 3;
            rbox[j] -= 3;
	    sbox[j] = rbox[j]-lbox[j]+1; 
            dxxx[j] = (rightside[j]-leftside[j])/((double)sbox[j]);
	    rdxxx[j] = 1.0 / dxxx[j];
 
        }

        levelbox = log10(dxbase[0]/dxxx[0])/log10(2.0);
        nlevelbox = (int)levelbox;
        
        //for(int j=0; j<nrank; j++)
        //printf("LEFT and RIGHT %lf >= %lf or  %lf <= %lf --> checking=1\n", leftside[j], rxbox[j],  rightside[j], lxbox[j]);


        //for(int j=0; j<nrank; j++)printf(">>>>>>>> %f <?  %f, %f >? %f\n", rightside[j], lxbox[j], leftside[j], rxbox[j]);
        for(int j=0; j<nrank; j++)if(rightside[j] <= lxbox[j] || leftside[j] >= rxbox[j])checking=1;

#ifdef DEBUG   
	if(checking)printf("Skipping block %d, not intersecting selected region...\n", i) ;
#endif
	if(checking)continue;

// next iteration if block is not assigned to the processor

        if(i < istart_pe || i >= iend_pe)continue;

        printf("BOX = %d, LEVEL = %d\n", i, nlevelbox);
        if(nlevelbox > maxlevel) continue;

        sourcesize = 1;
        destinationsize = 1;

	for(int j=0; j<nrank; j++) 
        {

          float raux;
#ifdef DEBUG   
          printf("CHECK CORRECTNESS %f %f %f %f\n", lxbox[j],leftside[j],rxbox[j],rightside[j]);
#endif

// intersection region in world coordinates

	  lintersect[j] = (leftside[j] >= lxbox[j])?leftside[j]:lxbox[j]; 
	  rintersect[j] = (rightside[j] <= rxbox[j])?rightside[j]:rxbox[j]; 

// intersection region in slab coordinates

	  raux = ((lintersect[j]-leftside[j])*rdxxx[j]);
	  lnintersect[j] = (int)raux;
          raux = ((rintersect[j]-leftside[j])*rdxxx[j])-1;
          rnintersect[j] = (int)raux;

// intersection region in box coordinates

	  float powconv = pow(2.0f, (float)(maxlevel-nlevelbox));

	  raux = ((lintersect[j]-lxbox[j])*rdxvngx[j]);
	  lbintersect[j] = (int)raux;
          //////raux = ((rintersect[j]-lxbox[j])*rdxvngx[j])-1;
          raux = lbintersect[j] + (float)(rnintersect[j]-lnintersect[j]+1)*powconv;
	  rbintersect[j] = (int)raux - 1; 

// calculate auxiliary read array size

          destinationsize *= (long)((float)(rnintersect[j]-lnintersect[j]+1)*powconv);
          sourcesize *= (rnintersect[j]-lnintersect[j]+1);
	  

	}

#ifdef DEBUG   
        printf("SOURCESIZE = %ld\n", sourcesize);
        printf("DESTSIZE   = %ld\n", destinationsize);

        for(int j=0; j<nrank; j++)printf("------> INTERSECT LEFT ABS: %f\n", lintersect[j]);
        for(int j=0; j<nrank; j++)printf("------> INTERSECT RIGHT ABS: %f\n", rintersect[j]);
        for(int j=0; j<nrank; j++)printf("------> SLAB INTERSECT LEFT: %d\n", lnintersect[j]);
        for(int j=0; j<nrank; j++)printf("------> SLAB INTERSECT RIGHT: %d\n", rnintersect[j]);
        for(int j=0; j<nrank; j++)printf("------> BOX INTERSECT LEFT: %d\n", lbintersect[j]);
        for(int j=0; j<nrank; j++)printf("------> BOXINTERSECT RIGHT: %d\n", rbintersect[j]);

/*

        printf ("%f, %f, %f\n", leftside[0],leftside[1],leftside[2]);
        printf ("%f, %f, %f\n", rightside[0],rightside[1],rightside[2]);
        printf ("%s\n", datafilename);

*/
#endif

        for (int j=0; j<nrank; j++)
        {

	    start[j]  = (hsize_t)lbintersect[j];
	    stride[j] = 1;
	    count[j]  = (hsize_t)(rbintersect[j]-lbintersect[j]+1);
	    block[j]  = 1;


	    s_start[j]  = (hsize_t)lnintersect[j] ;
	    s_stride[j] = 1;
	    s_count[j]  = (hsize_t)(rnintersect[j]-lnintersect[j]+1);
	    s_block[j]  = 1;

#ifdef DEBUG   
            printf ("ORIGIN START %d, ORIGIN COUNT %d\n", s_start[j],s_count[j]);
            printf ("DEST  START %d, DEST   COUNT %d\n", start[j],count[j]);
#endif
        }

/*
        printf("%d, %d, %d\n", nsize[0],nsize[1],nsize[2]);
        printf("%d, %d, %d\n", start[0],start[1],start[2]);
        printf("%d, %d, %d\n", count[0],count[1],count[2]);


*/

        total_size_old = total_size;
        total_size += sourcesize;
 

// allocate memory for loading data

	dataarray = new float[sourcesize];

        scalar->resize(total_size);
	xpos->resize(total_size);
	ypos->resize(total_size);
	zpos->resize(total_size);
        smooth->resize(total_size);
	

        for(long jk=0; jk<sourcesize; jk++)dataarray[jk]=0.0;

	

// open source file 

        printf("READING DATA FROM %s\n", datafilename);
	source_id = H5Fopen(datafilename, H5F_ACC_RDONLY, H5P_DEFAULT);

// read dataset from SOURCE

        int kaux = 0;
	for(int k=0; k<number_of_fields2read; k++)
	{
            kaux = sf[k];

            completename = hgroup;
	    completename.append(fieldsnames[kaux]);

            printf("READING %s from Grid %d\n",completename.c_str(),i);
	    source_obj = H5Dopen(source_id,completename.c_str());
	    source_space = H5Dget_space(source_obj);

   	    H5Sget_simple_extent_dims(source_space, s_dims, s_maxdims);
#ifdef DEBUG   
	    printf("%d, %d, %d\n", s_dims[0], s_dims[1], s_dims[2]);
#endif

// create auxiliary memory space for reading, select region and read

            memoryspace = H5Screate_simple (nrank, s_count, s_count);
	    
	    H5Sselect_hyperslab(source_space, H5S_SELECT_SET, s_start, s_stride, s_count, s_block);

	    H5Dread(source_obj, H5T_NATIVE_FLOAT, memoryspace, source_space, H5P_DEFAULT, dataarray);

            long jaux=0;
            printf("vector size = %ld\n", (long)scalar->size());
            for(long iaux=total_size_old; iaux<total_size; iaux++)
              {
               scalar->at(iaux) = dataarray[jaux];

               long ipp = jaux;
               long ic3d = (int)((float)ipp/(float)(s_count[1]*s_count[2]));
               long naux = ipp%(s_count[1]*s_count[2]);
               long jc3d = (int)((float)naux/((float)s_count[2]));
               long kc3d = naux - s_count[2]*jc3d;


	       xpos->at(iaux) = lintersect[0] + dxxx[0] * (float)ic3d;
	       ypos->at(iaux) = lintersect[1] + dxxx[1] * (float)jc3d;
	       zpos->at(iaux) = lintersect[2] + dxxx[2] * (float)kc3d;
               smooth->at(iaux) = dxxx[0];
	       minradius = (minradius <= dxxx[0] ? minradius : dxxx[0]);
	       maxradius = (maxradius >= dxxx[0] ? maxradius : dxxx[0]);
	    
               jaux++;
              }

            H5Sclose (memoryspace);

// end of loop over fields 
       }
        
// free memory

        H5Fclose (source_id);
	delete [] dataarray;
   }


#ifdef DEBUG   
   //printf("TOTAL NUMBER OF CELLS AT LEVEL 0 = %ld\n", gcounter);
#endif

//   H5Fclose (file_id);


   fclose (pFile);

   *maxr=maxradius;
   *minr=minradius;
#ifdef USEMPI
   MPI_Allreduce(&maxradius, maxr, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
   MPI_Allreduce(&minradius, minr, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
#endif

#ifdef DEBUG
   printf("PE %d MANAGE %ld POINTS\n", mype, total_size);
#endif

   return total_size;

}