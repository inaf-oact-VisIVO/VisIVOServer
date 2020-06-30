# include "Galaxy.h"
#define MAX(a,b) ((a < b) ?  (b) : (a))

float box_muller(float m, float s);

void CalculateColours (paramfile &params, string ComponentsName,long npart, 
                       float * cred, float * cgreen, float * cblue, float * ciii,
                       float * Red, float * Green, float * Blue, float * III, float * xcoord, 
                       float * ycoord, long nx, long ny)
{

	float xaux, yaux, xcol;
	long ii, jj, iaux;
        float x_rand_max = (float) RAND_MAX;
	float xcolaux;


	float brightness;
	float brightness_fact;
	brightness_fact = params.find<float>("Brightness"+ComponentsName,1.0);
	brightness = 1.0/(float)npart;
	brightness = brightness_fact*pow(brightness, 0.3333333f);
	cout << "--> Coloring " << ComponentsName.c_str() << " with brightness " << brightness << endl;

        float iiimax=-1.0;

	for (long particlei=0; particlei<npart; particlei++)
	{

	   xaux = (0.5*(xcoord[particlei]+1.0)); 
	   yaux = (0.5*(ycoord[particlei]+1.0)); 
	   ii = (int) (xaux*nx);
	   jj = (int) (yaux*ny);

	   if(ii >= nx || ii < 0 || jj >=ny || jj < 0)
	   {
//              xcol = ((float)rand())/x_rand_max;
	      xcolaux = box_muller(0, 0.25);
	      xcol = brightness*fabs(xcolaux);
	      if (xcol > 1.0) xcol = 0.0;

	      
	      cred[particlei]   = xcol;
	      cgreen[particlei] = xcol;
	      cblue[particlei]  = xcol;
	      ciii[particlei]   = xcol;

	   } else {
	      iaux = ii + jj*nx;
	      cred[particlei]   = Red[iaux];
	      cgreen[particlei] = Green[iaux];
	      cblue[particlei]  = Blue[iaux];
	      ciii[particlei]   = brightness*III[iaux];
	   }
	
           iiimax = MAX(iiimax, ciii[particlei]);

	}
        cout << "Maximum Brightness = " << iiimax << endl;
	//for (long particlei=0; particlei<npart; particlei++)ciii[particlei] /= iiimax;
}
