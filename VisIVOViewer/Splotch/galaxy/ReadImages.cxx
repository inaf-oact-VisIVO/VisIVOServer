# include <algorithm>
# include "Galaxy.h"


long ReadImages (paramfile &params, string infile_rgb, string infile_mask, long numx, long numy, float * RRR,
                 float * GGG, float * BBB, float * III, float * xx, float * yy, long nwant)
{

	  FILE * pFile;
	  FILE * pFile1;
	  float toll = 0.0;
          unsigned char * RR;
          unsigned char * GG;
          unsigned char * BB;
          unsigned char * II;
          float * x;
	  float * y;
          float * xM;
	  float * yM;

          long nskip = params.find<long>("Nskip",1);
          if(nskip < 1)exit(100);

// Read Data from raw file

	  if(!(pFile = fopen(infile_rgb.c_str(), "rb")))
	    {
	      printf("   could not open RGB file #%s#\n",infile_rgb.c_str());
	      exit(3);
	    }
	  if(!(pFile1 = fopen(infile_mask.c_str(), "rb")))
	    {
	      printf("   could not open Mask file #%s#\n",infile_mask.c_str());
	      exit(3);
	    }

	  printf("    reading image (%d,%d)\n",numx,numy);

// image size will be read from the FITS file: at the moment it's an input

//       read numx
//       read numy
	 long num = numx*numy;
	 RR = new unsigned char[num];
	 GG = new unsigned char[num];
	 BB = new unsigned char[num];
         II = new unsigned char[num];
	 x = new float [num];
	 y = new float [num];
	 xM = new float [num];
	 yM = new float [num];

	 long counter = 0;
	 float dist = 0.0;
	 float xaux, yaux;
	 float xxx=0.0;
	 float yyy=0.0;

// set the center of the galaxy (at the moment in term of pixels (default center of image)
         float norm = 0.5*(float)numx;
	 float xc = (float)numx/2/norm;
	 float yc = (float)numy/2/norm;

         long i=0;
          for (int iy=0; iy<numy; iy++)
          for (int ix=0; ix<numx; ix++)
	  { 
		fread(&RR[i], sizeof(char), 1, pFile);
		fread(&GG[i], sizeof(char), 1, pFile);
		fread(&BB[i], sizeof(char), 1, pFile);
		fread(&II[i], sizeof(char), 1, pFile1);
                i++;
          }                 
	  fclose(pFile);
	  fclose(pFile1);

	  //  for (long ix=1340; ix<1360; ix++)
	  //    {
	  //      for (long iy=1340; iy<1360; iy++)
	  //	{
	  //	  i = iy*numx + ix;
	  //	  if (II[i] > 0) cout << "1"; else cout << "0";
	  //	}
	  //      cout << endl;
	  //    }


// Now data processing
	  printf("    extracting mask\n",numx,numy);

          i=0;
          for (int iy=0; iy<numy; iy++)
          for (int ix=0; ix<numx; ix++)
          {

		RRR[i] = (float)RR[i]/255.0;
		GGG[i] = (float)GG[i]/255.0;
		BBB[i] = (float)BB[i]/255.0;
		III[i] = (float)II[i]/255.0;
		if(i%nskip == 0 && nskip != 1)III[i] = 0;

		x[i] = (float)ix;
		y[i] = (float)iy;
		x[i] = x[i]/norm-xc;
		y[i] = y[i]/norm-yc;

		if(III[i] > 0.0)
		{
		   xx[counter] = x[i];
		   yy[counter] = y[i];
		   counter++;
		   if(counter > nwant)
		     {
		       printf("Image mask produces more particles than desired (%d>%d), have to stop her ...\n",counter,nwant);
		       exit(3);
		     }
		}
		i++;
	  }

	  delete [] RR;
	  delete [] GG;
	  delete [] BB;
	  delete [] II;
	  delete [] x;
	  delete [] y;

	  printf("    Number of active pixels in mask: %ld\n", counter);

#ifdef WRITE_ASCII
          pFile = fopen("test.txt", "w");
          for(long ii=0; ii<counter; ii++)
          {
             fprintf(pFile, "%f %f 0.0\n", xx[ii],yy[ii]);
          }
          fclose(pFile);
#endif

	  return (counter);

}
