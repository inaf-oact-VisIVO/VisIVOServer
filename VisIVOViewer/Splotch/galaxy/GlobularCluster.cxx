# include "Galaxy.h"

float box_muller(float m, float s);

long GlobularCluster (paramfile &params, string ComponentName, long number_of_points, long npergroup, 
		      float * coordx, float * coordy, float * coordz)
{
  srand(time(NULL));

  float gsigma = params.find<float>("Sigma"+ComponentName,0);
  float sigma[3];
  sigma[0] = params.find<float>("Sigmax"+ComponentName,0);
  sigma[1] = params.find<float>("Sigmay"+ComponentName,0);
  sigma[2] = params.find<float>("Sigmaz"+ComponentName,0);

  printf("    Globular cluster like component with sigma_[x,y,z] = %f,%f,%f and gsigma=%f\n", sigma[0],sigma[1],sigma[2],gsigma);

  long counter=0;

  if(number_of_points > 0)
    {
      long ncluster = number_of_points / npergroup;

      for(long icluster=0;icluster<ncluster;icluster++)
	{
	  float x0 = box_muller(0.0, sigma[0]);
	  float y0 = box_muller(0.0, sigma[1]);
	  float z0 = box_muller(0.0, sigma[2]);

	  for(long istar=0;istar<npergroup;istar++)
	    {
	      coordx[counter] = box_muller(x0, gsigma);
	      coordy[counter] = box_muller(y0, gsigma);
	      coordz[counter] = box_muller(z0, gsigma);
	      counter++;
	    }
	}
    }

  return counter;

}
