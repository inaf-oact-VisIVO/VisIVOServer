#include <stddef.h>
#include <math.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "../libutility/utility.h"

int    icount;
double  mean_gM;
double  mean_gN;


double MONDification(double g_Mr, double g_0, double g_Nr)
{
  return (pow2(g_Mr) - g_Nr * sqrt(pow2(g_0) + pow2(g_Mr)));
}


/*======================================== 
 * simplification as given by Nusser 2001 
 *========================================*/
double get_MONDgNusser(double g_0, double g_N)
{

  /* simplest implementation yet ... no root-finding ! */
  if(g_N < g_0)
    return(sqrt(g_0*g_N));
  else
    return(g_N);
}

/*======================================================
 * proper implementation of mu(x) using Milgrom formula 
 *======================================================*/
#define ACCURACY  1E-7

double get_MONDg(double g_0, double g_N)
{
  double g_low, g_high, g_M;

#ifdef GM_ROOT
  g_low  = g_N;
  g_high = g_N + g_0/g_N;

  if(zbrac((double(*)(double,double,double))(MONDification), 
	   &g_low, &g_high, g_0, g_N))
    {
      g_M = rtbis((double(*)(double,double,double))(MONDification), 
		  g_low, g_high, ACCURACY, g_0, g_N);
    }
  else
    {
      fprintf(io.logfile,"\n zbrac did not work: %g %g %g %g\n\n",
	      g_0, g_N, g_low, g_high);
    }
#else
  g_M = g_N * sqrt( 0.5 + 0.5*sqrt(1+pow2(2*g_0/g_N)) );
#endif

  return(g_M);
}
