#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "utility.h"

/* check astro-ph/0409240 for deconvolution and shot-noise formulae */
/* #define TSCdeconvolve */
/* #define SHOTNOISE */

/*====================================================================
* power_spectrum:   calculate P(k) on domain grid
*====================================================================*/
void power_spectrum(flouble *dens)
{
   int     i, j, k, slen;
   int    *kmod, ks;
   double *dpower, FFTnorm;
   float  *rk, *Pk; 
   int    *jpower;
   int     NGRID, NYQUIST;
   FILE   *fpPk;
   char    Pkfile[MAXSTRING], file_no[MAXSTRING];
   double  shot_noise, rkg, rks, deconv, TSCdeconv, arg;
   
   NGRID    = simu.NGRID_DOM;
   FFTnorm  = (double)(NGRID*NGRID*NGRID);
   
   /* Nyquist frequency in fundamental mode unit */
   NYQUIST = NGRID/2;
   
   /* used with deconvolution of TSC mass assignment function */
   rkg     = NGRID;
   
   /* shot noise */
   shot_noise = pow3(simu.boxsize)/(simu.no_vpart*simu.pmass);

   /* allocate arrays */
   dpower = (double *) calloc(NGRID+1,           sizeof(double));
   jpower = (int *)    calloc(NGRID+1,           sizeof(int));
   kmod   = (int*)     calloc(NGRID,sizeof(long));

   /* this will help extracting information from dens^ */
   for(i=0; i<NGRID; i++)
     {
      if(i <= NYQUIST)
         kmod[i] = i;
      else
         kmod[i] = -(NGRID-i);
     }
   
   for(i=0; i<NYQUIST; i++)
     {
      dpower[i] = (double)0.0;
      jpower[i] = 0;
     }
   
   /* get P(k) */
   for(k=0; k<NGRID; k++)
     {
      for(j=0; j<NGRID; j++)
        {
         for(i=0; i<NGRID; i++)
           {
            
            /* spherical averaging in k-space */
            rks = sqrt( (float)pow2(kmod[i]) +
                        (float)pow2(kmod[j]) +
                        (float)pow2(kmod[k]) );
            ks  = (int) floor(rks);
            
            /* de-convolution of TSC mass assignment function */
#ifdef TSCdeconvolve
            arg        = PI*rks/rkg;
            deconv     = sin(arg)/(arg);
            TSCdeconv  = pow3(deconv);
#else
            TSCdeconv  = 1.0;
#endif
               
            if(ks >= 1 && ks <= NYQUIST)
              {
               /* add |delta|^2 values... */
               dpower[ks] += (double)(pow2(dens[Re(i,j,k,NGRID)]/TSCdeconv)+pow2(dens[Im(i,j,k,NGRID)]/TSCdeconv));               
               jpower[ks] += 1;
              }
           } 
        }
     }
      
   for(ks=1; ks<=NYQUIST; ks++)
     {
      /* convert to a more user friendly notation, i.e. k[h/Mpc] and Pk[(Mpc/h)^3] */
      PkSpectrum.rk[ks]  = ( (double)ks*TWOPI/simu.boxsize + (double)(ks+1)*TWOPI/simu.boxsize ) / 2.;
      PkSpectrum.Pk[ks]  = dpower[ks]/(double)jpower[ks]/FFTnorm * pow3(simu.boxsize)/FFTnorm;
#ifdef SHOTNOISE
      PkSpectrum.Pk[ks] -= shot_noise;
#endif
     }

   /* calc current growth factor */
   PkSpectrum.Dgrowth_now = calc_growth(global.a);
   
   /* store fundamental mode */
   PkSpectrum.Pk_now      = PkSpectrum.Pk[PKMODE];
   
   /* write P(k) to file? */
   if(PkSpectrum.dump_Pk == TRUE)
     {
      strcpy(Pkfile, io.outfile_prefix);
      
#ifdef REDSHIFTNAME
      /* OUTPUT CONVENTION: 3 digits*/
      sprintf(file_no,"z%.3f",fabs(global.z));
#else
      sprintf(file_no,"%05ld",global.no_timestep);
#endif
      strcat(Pkfile, file_no);

      /* add "_Pk" to file name */
      strcat(Pkfile,"_Pk");

      /*  open output file */
      if( (fpPk = fopen(Pkfile,"w")) == NULL)
          {
            fprintf(stderr,"could not open %s\n", Pkfile);
            exit(1);
          }
          
      /* also re-scale back to initial amplitude, i.e. take away linear growth */
      for(ks=1; ks<=NYQUIST; ks++)
         fprintf(fpPk,"%g %g %g\n", 
                 PkSpectrum.rk[ks], PkSpectrum.Pk[ks],
                 PkSpectrum.Pk[ks] * pow2(calc_growth(simu.a_initial)/calc_growth(global.a)));
      
      PkSpectrum.dump_Pk = FALSE;
      
      fclose(fpPk);
     }
   
   /* free temp arrays... */
   free(kmod);
   free(dpower);
   free(jpower);
}

