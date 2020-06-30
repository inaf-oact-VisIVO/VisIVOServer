# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <string>
using namespace std;


void CalculateDensity (float * hsml, float * rho, float * xcoord, float * ycoord, float * zcoord, long numofpart, float smooth)
{
	float distij;
	long * count = 0;
	float sl = 1.0;
	float toll = 0.2;
	float rhomax = 0.0;


	count = new long [numofpart];

#pragma omp parallel for

	for (long particlei=0; particlei<numofpart; particlei++)
	{
	   count[particlei] = 0;
	   if(particlei%1000 == 0)printf("Calculated %d particles over %d\n", particlei, numofpart);

	   for (long particlej=0; particlej<numofpart; particlej++)
//	   long particlej = 0;
//	   while(count[particlei] <= 100 && particlej <= numofpart)
	   {

		float xxx = xcoord[particlei]-xcoord[particlej];
		float yyy = ycoord[particlei]-ycoord[particlej];
		float zzz = zcoord[particlei]-zcoord[particlej];
		distij = sqrt(xxx*xxx+yyy*yyy+zzz*zzz);
		if(distij <= smooth) count[particlei]++;
//		particlej++;

	   }

	}

	for (long particlei=0; particlei<numofpart; particlei++)
	{
		rho[particlei] = (float)count[particlei];   ///////   /(float)numofpart;
		rhomax = max(rhomax, rho[particlei]);
	}

	for (long particlei=0; particlei<numofpart; particlei++)
	{
		if(rho[particlei] > 0.0)
	        {
		  hsml[particlei] = sl / rho[particlei];
		} else {
		  hsml[particlei] = 0.0;
		}
		if(hsml[particlei] > toll)hsml[particlei] = toll;
	}
}
