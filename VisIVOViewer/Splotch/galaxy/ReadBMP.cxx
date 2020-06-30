# include "Galaxy.h"


long ReadBMP (paramfile &params, string infile, string infile1, long numx, long numy, unsigned int Rmin, 
              unsigned int Gmin, unsigned int Bmin, unsigned int * RRR,
              unsigned int * GGG, unsigned int * BBB, unsigned int * III, float * xx, float * yy)
{

	  FILE * pFile;
	  FILE * pFile1;
	  float xc = (float)numx/2;
	  float yc = (float)numy/2;
	  float toll = 0.0;
	  long num = numx*numy;
          unsigned char * RR;
          unsigned char * GG;
          unsigned char * BB;
	  unsigned char * II;
	  int * xy;
          float * x;
	  float * y;
	  unsigned int RRaux;

	  RR = new unsigned char[num];
	  GG = new unsigned char[num];
	  BB = new unsigned char[num];
	  II = new unsigned char[num];
          xy = new int [num];
	  x = new float [num];
	  y = new float [num];

	  long i=0;
	  long counter = 0;
	  float dist = 0.0;
	  float xaux, yaux;
	  float xxx=0.0;
	  float yyy=0.0;
	  int color_depth;
	  
          color_depth = params.find<int>("ColorDepth",1);
	  if(color_depth == 1)
	  {
		Gmin = 255;
		Bmin = 255;
	  }

	  pFile = fopen(infile.c_str(), "rb");
	  pFile1 = fopen(infile1.c_str(), "rb");

          for (int iy=0; iy<numy; iy++)
	  {
          for (int ix=0; ix<numx; ix++)
	  { 
		fread(&RR[i], sizeof(char), 1, pFile);
		fread(&II[i], sizeof(char), 1, pFile1);

		if(color_depth == 2)
		{
		  fread(&GG[i], sizeof(char), 1, pFile);
		  fread(&BB[i], sizeof(char), 1, pFile);
		}else{
		  GG[i] = 0;
		  BB[i] = 0;
		}
                 
		RRR[i] = (unsigned int)RR[i];
		GGG[i] = (unsigned int)GG[i];
		BBB[i] = (unsigned int)BB[i];
		III[i] = (unsigned int)II[i];

		x[i] = (float)ix;
		y[i] = (float)iy;
		x[i] = x[i]-xc;
		y[i] = y[i]-yc;

		if(III[i] >= Rmin && RRR[i] < 257)
///////		if(RRR[i] >= Rmin)
		{
		   xy[i] = 1;
// choose x direction to normalize distances, in order to retain the real shape
// in this way x is ALWAYS between -1 and 1, but y is NOT!!!
		   xx[counter] = x[i]/xc;
		   yy[counter] = y[i]/xc;

		   counter++;
		    
		} else {
		   xy[i] = 0;
		}


		i++;
	  }
	  }
	  fclose(pFile);
	  fclose(pFile1);

	  delete [] RR;
	  delete [] GG;
	  delete [] BB;
	  delete [] II;

	  printf("=== NUMBER OF STARS : %ld ===\n", counter);
	  return (counter);

}
