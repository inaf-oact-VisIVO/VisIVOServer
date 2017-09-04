#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#ifdef WITH_OPENMP
#include <omp.h>
#endif

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "utility.h"
#include "../libio_serial/io_serial.h"
#include "../libmhd/mhd.h"

/* relevant parameters for c2f_slope() */
#define C2F_MIDCENTRE_NODELIMITER
//#define C2F_MIDCENTRE
//#define C2F_REDUCED

/* L1DIM_LENGTH = no. of char's for l1dim in filename (cf. write_filename) */
#define L1DIM_LENGTH 6

char *terminate_amiga = TERMINATE_AMIGA;

/* prepare filename for debugging purposes */
void write_filename(char *f_name, char *prefix, unsigned l1dim)
{
   int  slen, i;                  /* string length                  */
   char file_no[10];              /* file number                    */
   
   /* prepare filename */
   f_name = strcpy(f_name, prefix);
   
   itoa_(l1dim, file_no);
   slen = strlen(file_no);
   if (slen < L1DIM_LENGTH)
     {
      file_no[L1DIM_LENGTH] = '\0';
      
      for(i=L1DIM_LENGTH-1; i >= (L1DIM_LENGTH-slen); i--)
        {
         file_no[i] = file_no[i-(L1DIM_LENGTH-slen)];
        }
      for(i=0; i < (L1DIM_LENGTH-slen); i++)
        {
         file_no[i] = '0';
        }
     }
   f_name = strcat(f_name, file_no);
   f_name = strcat(f_name, ".DAT");
   
}

/*===========================================================================
 * calculate Laplace operator acting on potential at a given node 
 *===========================================================================*/
double Laplace_pot(nptr tsc_nodes[3][3][3], double spacing2)
{
  double Lpot;
  
  Lpot = (  (double)tsc_nodes[1][1][2]->pot + (double)tsc_nodes[1][1][0]->pot
          + (double)tsc_nodes[1][2][1]->pot + (double)tsc_nodes[1][0][1]->pot
          + (double)tsc_nodes[2][1][1]->pot + (double)tsc_nodes[0][1][1]->pot
          - 6. * (double)tsc_nodes[1][1][1]->pot) / spacing2;
  
  return Lpot;
}

/*===========================================================================
 * calculate Laplace operator acting on temp[1] at a given node 
 *===========================================================================*/
double Laplace_temp1(nptr tsc_nodes[3][3][3], double spacing2)
{
  double Ltemp1;
  
  Ltemp1 = (  (double)tsc_nodes[1][1][2]->force.temp[1] + (double)tsc_nodes[1][1][0]->force.temp[1]
            + (double)tsc_nodes[1][2][1]->force.temp[1] + (double)tsc_nodes[1][0][1]->force.temp[1]
            + (double)tsc_nodes[2][1][1]->force.temp[1] + (double)tsc_nodes[0][1][1]->force.temp[1]
            - 6. * (double)tsc_nodes[1][1][1]->force.temp[1]) / spacing2;
  
  return Ltemp1;
}


/* f1mod:  x modulo y for double numbers */
double f1mod(double x, double y)
{
  /* the following part is fully tuned for use with AMIGA !!!!! */
   if(x >= 2.0)
      return(x-2.0);
   else if(x >= 1.0)
      return(x-1.0);
   else
      return(x);   
}


/*==============================================================================
*  get eigenvalues of inertia tensor
*==============================================================================*/
void get_axes(double itensor[3][3], double *axis1, double *axis2, double *axis3)
{
   int           n, i, j, nrot;
   unsigned long idx[4];
   double        a[4][4], d[4], v[4][4], tmp[4];
   
   n = NDIM;
   
   for(i=0; i<4; i++)
      for(j=0; j<4; j++)
         a[i][j] = 0.0;
   
   a[1][1] = itensor[0][0];
   a[2][1] = itensor[1][0];
   a[3][1] = itensor[2][0];
   a[1][2] = itensor[0][1];
   a[2][2] = itensor[1][1];
   a[3][2] = itensor[2][1];
   a[1][3] = itensor[0][2];
   a[2][3] = itensor[1][2];
   a[3][3] = itensor[2][2];
   
   jacobi(a, n, d, v, &nrot);
   
   for(i=1; i<=n; i++)
      tmp[i] = d[i];
   indexx((unsigned long)n, tmp, idx);
   
   *axis1 = d[idx[3]];
   *axis2 = d[idx[2]];
   *axis3 = d[idx[1]];

   itensor[0][0] = v[1][idx[3]];
   itensor[1][0] = v[2][idx[3]];
   itensor[2][0] = v[3][idx[3]];
   
   itensor[0][1] = v[1][idx[2]];
   itensor[1][1] = v[2][idx[2]];
   itensor[2][1] = v[3][idx[2]];
   
   itensor[0][2] = v[1][idx[1]];
   itensor[1][2] = v[2][idx[1]];
   itensor[2][2] = v[3][idx[1]];
}

/*==============================================================================
*  inverts the usage of idx[] returned by indexx sorting routine
*==============================================================================*/
int idx_inv(unsigned long *idx, int numHalos, int i)
{
   int j;
   
   for(j=0; j<numHalos; j++)
     {
      if(idx[j] == i)
         break;
     }
   
   return(j);
}

/*
 ************************************************************
 ************************************************************
 *  Calculates the minima and the maxima
 */
MINMAX MinMax(double x,double xmin,double xmax) {
   
   MINMAX tmpMinMax;		
   
   if (x < xmin)
      tmpMinMax.min = x;
   else
      tmpMinMax.min = xmin;
   
   if (x > xmax)
      tmpMinMax.max = x;
   else
      tmpMinMax.max = xmax;
   
   
   return(tmpMinMax);
   
}
/*
 ************************************************************
 ************************************************************
 *  Calculates the minima and the maxima for the periodic refinements
 */
MINMAX MinMaxBound(double div, double x,double xmin,double xmax) {
   
   MINMAX tmpMinMax;		
   
   
   if ( x < div ) {
      
      if (x > xmax)
         tmpMinMax.max = x;
      else
         tmpMinMax.max = xmax;
      
      tmpMinMax.min = xmin;
      
   } else {
      
      if (x < xmin)
         tmpMinMax.min = x;
      else
         tmpMinMax.min = xmin;
      
      tmpMinMax.max = xmax;
						
   }
   
   return(tmpMinMax);
   
}

/*==============================================================================
 * initialize binning for density profiles
 *==============================================================================*/
void binning_parameter(HALO halo, int *nbins, double *dist_min, double *dist_max)
{
   partptr cur_part;
   double  dX, dY, dZ, Xc, Yc, Zc;
   long    npart, binpart;
   
   /* how many bins should be used for cur_halo */
   *nbins   = (int) (6.2*(log10((double)halo.npart))-3.5);
   *nbins   = MAX(1,*nbins);
   binpart  = (double)(halo.npart-1)/(double)*nbins;
   
   /* halo position in AMIGA units */
   Xc = halo.pos.x;
   Yc = halo.pos.y;
   Zc = halo.pos.z;

   /* first particle */
   cur_part = global.fst_part + halo.ipart[0];
   
   /* minimum distance */
   npart = 0;
#ifndef TRACKER
   while(npart < halo.npart && npart < AHF_MINPART/10.)
#else
   while(npart < halo.npart && npart < TRK_MINPART/10.)
#endif
     {
        npart++;
        cur_part = global.fst_part + halo.ipart[npart];
     }      
   dX = fabs(cur_part->pos[X] - Xc);
   dY = fabs(cur_part->pos[Y] - Yc);
   dZ = fabs(cur_part->pos[Z] - Zc);
   if(dX > 0.5) dX = 1.0-dX;
   if(dY > 0.5) dY = 1.0-dY;
   if(dZ > 0.5) dZ = 1.0-dZ;
   *dist_min  = sqrt(pow2(dX) + pow2(dY) + pow2(dZ));
   
   /* maximum distance */
   cur_part = global.fst_part + halo.ipart[halo.npart-1];
   dX = fabs(cur_part->pos[X] - Xc);
   dY = fabs(cur_part->pos[Y] - Yc);
   dZ = fabs(cur_part->pos[Z] - Zc);
   if(dX > 0.5) dX = 1.0-dX;
   if(dY > 0.5) dY = 1.0-dY;
   if(dZ > 0.5) dZ = 1.0-dZ;
   *dist_max  = sqrt(pow2(dX) + pow2(dY) + pow2(dZ));   
}

/*=============================================================================
 * calculate DRIFT and KICK operators
 *=============================================================================*/
void get_DKoperators(double timestep, double *KICK, double *DRIFT1, double *DRIFT2)
{
  *KICK      = timestep;
  *DRIFT1    = timestep/2;
  *DRIFT2    = timestep/2;
  
  
#ifdef GLASS
  /* this gives repulsive forces */
  *KICK *= -1.0;
#endif
  
#ifdef REVERSE_INTEGRATION
  /* this reverses the time integration */
  *KICK   *= -1.0;
  *DRIFT1 *= -1.0;
  *DRIFT2 *= -1.0;
#endif
  
}

/*===============================================================
* return a timestep based upon user data and NSTEPS parameter
*===============================================================*/
double init_timestep(uparam user)
{
   double no_steps, timestep;
   
   no_steps = NSTEPS;
   
   if(io.header.no_timestep > NSTEPS) no_steps  = 100 * io.header.no_timestep;
   if(user.final_z < 0)               no_steps *= 100;
   
   timestep = (simu.super_t_final-simu.super_t_initial)/(double)(no_steps);
   
   return(timestep);
}

/*========================================================================
 * re-adjust timestep according to various criteria...
 *
 * 1. speeder criterion (i.e. dark matter CFL vriterion)
 * 2. CFL criterion for hydro (if required)
 * 3. cosmological criterion
 *========================================================================*/
double adjust_timestep(double timecounter, double timestep)
{
   double  speeders;
   double  a, da_a;
   double  new_timestep, cfl_timestep;
   
#if (HYDRO_TEST==5 || HYDRO_TEST==6 || HYDRO_TEST==7)
   /* switch off the speeder criterion */
   global.speeding    = FALSE;
   global.no_speeders = 0;
   global.max_dr2     = pow2(MIN(1.5*CELLFRAC_MAX, 1.5*CELLFRAC_MIN));
#endif
   
   /*=========================================
    * ...start with old timestep
    *=========================================*/
   new_timestep = timestep;

   /*=========================================
    * 1. "particles flying too far" criterion
    *=========================================*/
   speeders = (double)global.no_speeders/(double)global.no_part*100.;
   if((global.speeding == TRUE && speeders > SPEEDFRAC) ||
      (global.max_dr2 < pow2(  CELLFRAC_MIN))           ||
      (global.max_dr2 > pow2(2*CELLFRAC_MAX))             )
     {
      new_timestep = timestep * CF_MEAN/sqrt(global.max_dr2);
     }
   
#ifdef HYDRO
   /*=========================================
    * 2. CFL criterion
    *=========================================*/
   a            = calc_super_a(timecounter);
   cfl_timestep = CFL_TUNE * a * global.dom_grid->spacing/global.cfl_speed;
   
   new_timestep = MIN(new_timestep, cfl_timestep);
#endif
   
   /*=========================================
    * 3. cosmological criterion
    *=========================================*/
   da_a = 2*(calc_super_a(timecounter+new_timestep)-calc_super_a(timecounter))/
            (calc_super_a(timecounter+new_timestep)+calc_super_a(timecounter));
   
   /* iteratively reduce dt (not very sophisticated but should work...) */
   while(da_a > CA_CRIT)
     {
      new_timestep *= 0.9;
      da_a          = 2*(calc_super_a(timecounter+new_timestep)-calc_super_a(timecounter))/
         (calc_super_a(timecounter+new_timestep)+calc_super_a(timecounter));
     }

   /*=========================================
    * check if timestep is too big
    *=========================================*/
   if(global.super_t+new_timestep > simu.super_t_final)
      new_timestep = simu.super_t_final-global.super_t;
   
   /*=========================================
    * keep track of new timestep in io.header
    *=========================================*/
   io.header.timestep = new_timestep;
   
   return(new_timestep);
}

/*==============================================================================
 * read AMIGA header 
 *==============================================================================*/
void read_amiga_header(FILE *infile, info_io *io, int *SWAPBYTES)
{
   int       file_sizeof_long, machine_sizeof_long, one;
   int       idummy;
   long      ldummy;
   long long lldummy;
      
#ifdef AMIGA_ONE_FORMAT
   fread(&one, sizeof(int), 1, infile);
   if(one != 1)
     {
      *SWAPBYTES = TRUE;
      fprintf(stderr," => start reading BYTESWAPed AMIGA header ('one' format assuming sizeof(long)=4!) ... ");
     }
   else
     {
      *SWAPBYTES = FALSE;
      fprintf(stderr," => start reading AMIGA header ('one' format assuming sizeof(long)=4!) ... ");
     }

   /* don't bother with this check... */
   machine_sizeof_long = file_sizeof_long = 4;
#else
   /* do we need to do BYTESWAP (simultaneously determine the size of a long...) */
   fread(&(file_sizeof_long), sizeof(int), 1, infile);
   if(file_sizeof_long > 8)
     {
      *SWAPBYTES = TRUE;
      sexchange(&file_sizeof_long, sizeof(int));
      fprintf(stderr," => start reading BYTESWAPed AMIGA header ... ");
     }
   else
     {
      *SWAPBYTES = FALSE;
      fprintf(stderr," => start reading AMIGA header ... ");
     }
   machine_sizeof_long = sizeof(long);
   
   if(machine_sizeof_long != file_sizeof_long)
      fprintf(stderr,"(sizeof(long) mismatch: %d vs. %d!) ",machine_sizeof_long,file_sizeof_long);
#endif
   
   if(*SWAPBYTES == TRUE || machine_sizeof_long != file_sizeof_long)
     {
      /* read in IO header */
      ReadChars(infile,io->header.header,HEADERSTRING);
      ReadInt(infile,&(io->header.multi_mass),*SWAPBYTES);
      ReadInt(infile,&(io->header.double_precision),*SWAPBYTES);
      
      /* read 2x long unsigned */
      if(machine_sizeof_long == 8 && file_sizeof_long == 4)
        {
         ReadInt(infile,&idummy,*SWAPBYTES);
         io->header.no_part    = (long unsigned) idummy;
         ReadInt(infile,&idummy,*SWAPBYTES);
         io->header.no_species = (long unsigned) idummy;
        }
      else if(machine_sizeof_long == 4 && file_sizeof_long == 8)
        {
         ReadLongLong(infile,&lldummy,*SWAPBYTES);
         io->header.no_part    = (long unsigned) lldummy;
         ReadLongLong(infile,&lldummy,*SWAPBYTES);
         io->header.no_species = (long unsigned) lldummy;
        }
      else if((machine_sizeof_long == 8 && file_sizeof_long == 8) ||
              (machine_sizeof_long == 4 && file_sizeof_long == 4))
        {
         ReadLong(infile,&ldummy,*SWAPBYTES);
         io->header.no_part    = (long unsigned) ldummy;
         ReadLong(infile,&ldummy,*SWAPBYTES);
         io->header.no_species = (long unsigned) ldummy;
        }
      else
        {
         fprintf(stderr,"read_amiga_header: machine_sizeof_long=%d vs. file_sizeof_long=%d\n",
                 machine_sizeof_long, file_sizeof_long);
         exit(0);
        }
      
      ReadDouble(infile,&(io->header.no_vpart),*SWAPBYTES);
      ReadDouble(infile,&(io->header.timestep),*SWAPBYTES);
      ReadInt(infile,&(io->header.no_timestep),*SWAPBYTES);
      if(machine_sizeof_long == 4 && file_sizeof_long == 8)
         ReadInt(infile,&idummy,*SWAPBYTES);
      ReadDouble(infile,&(io->header.boxsize),*SWAPBYTES);
      ReadDouble(infile,&(io->header.omega0),*SWAPBYTES);
      ReadDouble(infile,&(io->header.lambda0),*SWAPBYTES);
      ReadDouble(infile,&(io->header.pmass),*SWAPBYTES);
      ReadDouble(infile,&(io->header.cur_reflevel),*SWAPBYTES);
      ReadDouble(infile,&(io->header.cur_frcres),*SWAPBYTES);
      ReadDouble(infile,&(io->header.a_initial),*SWAPBYTES);
      ReadDouble(infile,&(io->header.a_current),*SWAPBYTES);
      ReadDouble(infile,&(io->header.K_initial),*SWAPBYTES);
      ReadDouble(infile,&(io->header.K_current),*SWAPBYTES);
      ReadDouble(infile,&(io->header.U_initial),*SWAPBYTES);
      ReadDouble(infile,&(io->header.U_current),*SWAPBYTES);
      ReadDouble(infile,&(io->header.Eintegral),*SWAPBYTES);
      ReadDouble(infile,&(io->header.Econst),*SWAPBYTES);
      
      ReadDouble(infile,&(io->header.paramNSTEPS),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramNGRID_DOM),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramNth_dom),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramNth_ref),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramE_UPDATE),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramCELLFRAC_MAX),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramCELLFRAC_MIN),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramCA_CRIT),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramMAX_L1DIM),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramDOMSWEEPS),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramREFSWEEPS),*SWAPBYTES);
      
      ReadDouble(infile,&(io->header.paramAHF_MINPART),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramAHF_VTUNE),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramAHF_RISE),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramAHF_SLOPE),*SWAPBYTES);
      ReadDouble(infile,&(io->header.paramAHF_MAXNRISE),*SWAPBYTES);
      
      ReadDouble(infile,&(io->header.min_weight),*SWAPBYTES);
      ReadDouble(infile,&(io->header.max_weight),*SWAPBYTES);
      ReadDouble(infile,&(io->header.t_unit),*SWAPBYTES);
      ReadDouble(infile,&(io->header.B_init),*SWAPBYTES);
      ReadDouble(infile,&(io->header.param_dummy5),*SWAPBYTES);
      ReadDouble(infile,&(io->header.param_dummy6),*SWAPBYTES);
      ReadDouble(infile,&(io->header.param_dummy7),*SWAPBYTES);
      ReadDouble(infile,&(io->header.param_dummy8),*SWAPBYTES);
      
      ReadDouble(infile,&(io->header.version),*SWAPBYTES);
      ReadInt(infile,&(io->header.built),*SWAPBYTES);
      if(machine_sizeof_long == 4 && file_sizeof_long == 8)
         ReadInt(infile,&idummy,*SWAPBYTES);
      
      ReadDouble(infile,&(io->header.omegab),*SWAPBYTES);
      ReadDouble(infile,&(io->header.gamma), *SWAPBYTES);
      ReadDouble(infile,&(io->header.H_frac),*SWAPBYTES);
      ReadDouble(infile,&(io->header.T_init),*SWAPBYTES);

      ReadInt(infile,&(io->header.hydro),*SWAPBYTES);
      if(machine_sizeof_long == 4 && file_sizeof_long == 8)
         ReadInt(infile,&idummy,*SWAPBYTES);
      ReadInt(infile,&(io->header.magneto),*SWAPBYTES);
      if(machine_sizeof_long == 4 && file_sizeof_long == 8)
         ReadInt(infile,&idummy,*SWAPBYTES);

      /* take care of this "long issue" when reading the FILLHEADER stuff */
      if(machine_sizeof_long == 8 && file_sizeof_long == 4)
        {
         /* we read 2x4 bytes fewer than FILLHEADER... */
         ReadChars(infile,io->header.dummy,FILLHEADER +4+4);
        }
      else if(machine_sizeof_long == 4 && file_sizeof_long == 8)
        {
         /* we read 2x4 bytes more than FILLHEADER... */
         ReadChars(infile,io->header.dummy,FILLHEADER -4-4);
        }
      else if((machine_sizeof_long == 8 && file_sizeof_long == 8) ||
              (machine_sizeof_long == 4 && file_sizeof_long == 4))
        {
         ReadChars(infile,io->header.dummy,FILLHEADER);
        }
     } 
   else /* SWAPBYTES */
     {
      /* read header as complete structure */
      if(fread(&(io->header), sizeof(io->header), 1, infile) != 1)
        {
         fprintf(stderr,"\n\ninput: could not read AMIGA IO header\n");
         fflush(stderr);
         fclose(stderr);
         exit(1);
        }
     } /* SWAPBYTES */
  fprintf(stderr,"done <=\n");
}


/*===================================================================================
* get_c2fslope:   calculate the slope f'(x_i) needed with c2f_pot
*
*             f(x)    = f(x_i) + f'(x_i) * (x-x_i)
*
*===================================================================================*/
void get_c2fslope(double func[3][3][3], double slope[NDIM])
{
   double slope_left, slope_right;
   
   /* not that the spacing of the sampling of func is 1 */
   
#ifdef C2F_MIDCENTRE_NODELIMITER
   
   /* mid-centred slope (no delimiter!) */
   slope[X] = (func[1][1][2] - func[1][1][0]) / 2.;
   slope[Y] = (func[1][2][1] - func[1][0][1]) / 2.;
   slope[Z] = (func[2][1][1] - func[0][1][1]) / 2.;
   
#endif /* C2F_MIDCENTRE_NODELIMITER */
   
#ifdef C2F_MIDCENTRE
   
   /* mide-centred slope */
   slope[X] = (func[1][1][2] - func[1][1][0]) / 2.;
   slope[Y] = (func[1][2][1] - func[1][0][1]) / 2.;
   slope[Z] = (func[2][1][1] - func[0][1][1]) / 2.;
   
   /* delimit the slope */
   
   slope_left  = func[1][1][1] - func[1][1][0];
   slope_right = func[1][1][2] - func[1][1][1];
   
   if(slope_left*slope_right < 0.)
      slope[X] = 0.0;
   
   slope_left  = func[1][1][1] - func[1][0][1];
   slope_right = func[1][2][1] - func[1][1][1];
   
   if(slope_left*slope_right < 0.)
      slope[Y] = 0.0;
   
   slope_left  = func[1][1][1] - func[0][1][1];
   slope_right = func[2][1][1] - func[1][1][1];
   
   if(slope_left*slope_right < 0.)
      slope[Z] = 0.0;
   
#endif /* C2F_MIDCENTRE */
   
#ifdef C2F_REDUCED
   
   /* delimit the slope */
   
   slope[X] = 0.0;
   slope[Y] = 0.0;
   slope[Z] = 0.0;
   
   slope_left  = func[1][1][1] - func[1][1][0];
   slope_right = func[1][1][2] - func[1][1][1];
   
   if(slope_left*slope_right > MACHINE_ZERO)
      slope[X] = 2*slope_left*slope_right / (slope_left+slope_right);
   
   slope_left  = func[1][1][1] - func[1][0][1];
   slope_right = func[1][2][1] - func[1][1][1];
   
   if(slope_left*slope_right > MACHINE_ZERO)
      slope[Y] = 2*slope_left*slope_right / (slope_left+slope_right);
   
   slope_left  = func[1][1][1] - func[0][1][1];
   slope_right = func[2][1][1] - func[1][1][1];
   
   if(slope_left*slope_right > MACHINE_ZERO)
      slope[Z] = 2*slope_left*slope_right / (slope_left+slope_right);
   
#endif /* C2F_REDUCED */
   
}

/*===============================================================
 * open AMIGA logfile and start dumping information...
 *===============================================================*/
void init_logfile(uparam user)
{
  int i;
  
  if((io.logfile = fopen(io.logfile_name,"r")) == NULL)
    {
      if((io.logfile = fopen(io.logfile_name,"w")) == NULL) {
        fprintf(stderr,"\n\nstartrun: could not open file %s\n", io.logfile_name);
      exit(1); }
      else {
      fprintf(stderr,"\n...opening new logfile...\n\n"); }
    }
  else
    {
      /* logfile already exists...re-open to append data */
      fclose(io.logfile);
      io.logfile = fopen(io.logfile_name,"a");
      fprintf(stderr,"\n...appending to already existing logfile...\n");
    }
  
  
  fprintf(io.logfile,"================================================================\n");
  fprintf(io.logfile,"\t     AA      M   M    I    GGGG      AA    \n");
  fprintf(io.logfile,"\t   A   A    MM MM    I    G        A   A   \n");
  fprintf(io.logfile,"\t  AAAAA    M M M    I    G GGG    AAAAA   \n");
  fprintf(io.logfile,"\t A   A    M   M    I    G   G    A   A   \n");
  fprintf(io.logfile,"\tA   A    M   M    I     GGG     A   A  (v%3.1f/%d)\n",VERSION,BUILT);
  fprintf(io.logfile,"================================================================\n\n");
  
//  fprintf(io.logfile,"\n");
//  fprintf(io.logfile,"if you feel like stopping AMIGA, please type:\n");
//  fprintf(io.logfile,"=============================================\n");
//  fprintf(io.logfile,"$touch %s\n\n",TERMINATE_AMIGA);
  
  fprintf(io.logfile,"\n");
  fprintf(io.logfile,"user input:\n");
  fprintf(io.logfile,"===========\n");
#ifndef ART	  
  fprintf(io.logfile,"IC file                    = %s\n",io.icfile_name);
#endif
  fprintf(io.logfile,"dump file                  = %s\n",io.dumpfile_name);
  fprintf(io.logfile,"prefix for output files    = %s\n",io.outfile_prefix);
#ifdef LIGHTCONE
  fprintf(io.logfile,"prefix for lightcone files = %s\n",io.lightcone_prefix);
#endif /* LIGHTCONE */
  
  fprintf(io.logfile,"please give number of domain grid cells (1D):     ");
  fprintf(io.logfile,"%d\n",user.NGRID_DOM);
  fprintf(io.logfile,"please give Nth for domain grid:                  ");
  fprintf(io.logfile,"%g\n",user.Nth_dom);
  fprintf(io.logfile,"please give Nth for refinements:                  ");
  fprintf(io.logfile,"%g\n",user.Nth_ref);
  fprintf(io.logfile,"please give final redshift:                       ");
  fprintf(io.logfile,"%g\n",user.final_z);
  
#ifdef LIGHTCONE
  fprintf(io.logfile,"please give lightcone redshift limit:             ");
  fprintf(io.logfile,"%g\n",user.lightcone_z);
  fprintf(io.logfile,"please give lightcone type:                       ");
  fprintf(io.logfile,"%d\n",user.lightcone_type);
  if(lightcone_type==-1) {
    fprintf(io.logfile,"please give patch orientation angles (degrees): ");
    fprintf(io.logfile,"%g %g %g\n",user.sa,user.sb,user.sg);
    fprintf(io.logfile,"please give patch sizes (degrees):              ");
    fprintf(io.logfile,"%g %g\n",user.dcx,user.dcy);
  }
#endif /* LIGHTCONE */
  
#ifdef TIPSY
  fprintf(io.logfile,"please give box size [Mpc/h]:                     ");
  fprintf(io.logfile,"%lf\n", tipsy_boxsize);
  fprintf(io.logfile,"please give omega0:                               ");
  fprintf(io.logfile,"%lf\n", tipsy_omega0);
  fprintf(io.logfile,"please give lambda0:                              ");
  fprintf(io.logfile,"%lf\n", tipsy_lambda0);
  fprintf(io.logfile,"please give initial redshift:                     ");
  fprintf(io.logfile,"%lf\n", tipsy_initalz);
  fprintf(io.logfile,"please give current timestep no:                  ");
  fprintf(io.logfile,"%lf\n", tipsy_currentimeno);
#endif /* TIPSY */
  
  fprintf(io.logfile,"please give dump step:                            ");
  fprintf(io.logfile,"%d\n",user.out_dumps);
  fprintf(io.logfile,"please give total number of output files:         ");
  fprintf(io.logfile,"%d\n",user.no_outputs);
  
#ifndef ISOLATED /* for -DISOLATED the output times are *not* provided by the user! */
  for(i = 0; i < user.no_outputs; i++)
    {
      fprintf(io.logfile,"please give redshift for %5d. output file:          ",i+1);
      fprintf(io.logfile,"%g\n",user.z_out[i]);
    }
#endif
  
#if (defined ART && defined HYDRO)
  fprintf(io.logfile,"please give omegab:                               ");
  fprintf(io.logfile,"%g\n",user.omegab);
  fprintf(io.logfile,"please give gamma:                                ");
  fprintf(io.logfile,"%g\n",user.gamma);
  fprintf(io.logfile,"please give H_frac:                               ");
  fprintf(io.logfile,"%g\n",user.H_frac);
  fprintf(io.logfile,"please give T_init:                               ");
  fprintf(io.logfile,"%g\n",user.T_init);
  fprintf(io.logfile,"please give B_init:                               ");
  fprintf(io.logfile,"%g\n",user.B_init);
#endif
  
#ifdef MOND /* for -DMOND two additional parameters need to be provided */
  fprintf(io.logfile,"please give MOND acceleration g0:                 ");
  fprintf(io.logfile,"%g\n",user.g0);
  fprintf(io.logfile,"please give Hubble parameter H0:                  ");
  fprintf(io.logfile,"%g\n",user.h0);
#endif
  
  fflush(io.logfile);
}


/*=======================================================================
 * we simply copy the user provided parameters over to global. structure
 *=======================================================================*/
void init_global_structure(uparam user)
{
  
  char termfile_name [MAXSTRING];

  fprintf(stderr,"\n => initializing global. structure ... ");

   
  /* global access to the particles */
  global.fst_part    = io.fst_part;
  global.no_part     = io.no_part;
  
  global.fst_gas     = io.fst_gas;
  global.no_gas      = io.no_gas;
  global.offset_gas  = io.offset_gas;
  
  global.fst_star    = io.fst_star;
  global.no_stars    = io.no_stars;
  global.offset_stars= io.offset_stars;
  
  /* the current time variables */
  global.a           = io.header.a_current;
  global.z           = 1./global.a - 1.;
  global.t           = calc_t(global.a);
  global.super_t     = calc_super_t(global.a);
  global.no_timestep = io.header.no_timestep;

  /* we are re-starting a simulation... */
  if(fabs(global.super_t-simu.super_t_initial) > ZERO)
    {
     /* ...but io.header.no_timestep is zero? correct that mistake! */
     if(global.no_timestep == 0)
        global.no_timestep = (int)((global.super_t-simu.super_t_initial) / init_timestep(user)) + 1;
    }
  
  
  /* the termfile_name is a character string defined in param.h */
  strcpy(termfile_name,terminate_amiga);
  global.termfile_name = (char *)malloc((strlen(termfile_name)+1) *sizeof(char));
  strcpy(global.termfile_name, termfile_name);
  
  /* check, if machine is little or big endian before reading IC's */
  global.architecture = test_endian();
  
  /* how much memory per node and particle for this particular run */
  global.bytes_node = sizeof(struct node);
  global.bytes_part = sizeof(struct particle);
  
  /* set terminate flag to FALSE */
  global.terminate = FALSE;
  
  /* restart flag */
  if(io.header.no_timestep == 0)
    global.restart = FALSE;


#ifdef MOND
  global.max_gN          = 0.0;
  global.max_gM          = 0.0;
  global.steps           = 0;
  global.no_ALLevents    = 0;
  global.no_MONDevents   = 0;
#endif
  
  
  fprintf(stderr,"done\n");

}
  
  
/*=======================================================================
 * we simply copy the user provided parameters over to simu. structure
 *=======================================================================*/
void init_simu_structure(uparam user)
{
   
  fprintf(stderr,"\n => initializing simu. structure ... \n");
   
  /* cosmology */
  simu.boxsize = io.header.boxsize;
  simu.omega0  = io.header.omega0;
  simu.lambda0 = io.header.lambda0;
  
  /*----------------------------------------------------------------------------------------
   * the unit stuff
   *
   * simu.pmass
   * simu.t_unit
   * simu.FourPiG
   * simu.mean_dens
   *
   * IMPORTANT: simu.pmass is modified in init_DMhydro() according to the baryon fraction!
   *---------------------------------------------------------------------------------------*/
  simu.pmass   = io.header.pmass;
  simu.t_unit  = io.header.t_unit;
#ifdef NO_EXPANSION   
  simu.FourPiG   = (double) 1.0;         // we explicitly chose the time unit to give FourPiG=1
  simu.mean_dens = (double) 0.0;         // cosmological(!) background density in internal units
#else
  simu.FourPiG   = 1.5*io.header.omega0; // we explicitly chose the time unit to give FourPiG=1.5*omega0
  simu.mean_dens = (double) 1.0;         // cosmological(!) background density in internal units
#endif /* NO_EXPANSION*/
  
  
  
  /* start and end points */
  simu.a_initial       = io.header.a_initial;
  simu.a_final         = 1./(1.+user.final_z);
  simu.z_initial       = 1./simu.a_initial - 1.;
  simu.z_final         = 1./simu.a_final   - 1.;
  
#ifndef CONVERT_TERM
  /*----------------------------------------------------------------
   * create timeline to be able to use calc_x() from cosmology
   * 
   * we need to set cosmology and units prior to create_timeline()!
   *----------------------------------------------------------------*/
   
  /* create_timeline does not work for a_initial=a_final */
  if(fabs(simu.a_initial-simu.a_final) < ZERO)
    simu.a_initial = 0.01; /* redshift z=99 */

  create_timeline(simu.a_initial/10., simu.a_final, &simu.timeline);
#endif
  
  simu.t_initial       = calc_t(simu.a_initial);
  simu.t_final         = calc_t(simu.a_final);
  simu.super_t_initial = calc_super_t(simu.a_initial);
  simu.super_t_final   = calc_super_t(simu.a_final);   
  
  
  /* HYDRO stuff */
#ifdef HYDRO
  simu.gamma           = io.header.gamma;
  simu.omegab          = io.header.omegab;
  simu.omegaDM         = simu.omega0 - simu.omegab;
  simu.f_b             = simu.omegab/simu.omega0;
  simu.H_frac          = io.header.H_frac;
  simu.T_init          = io.header.T_init;
  simu.B_init          = io.header.B_init;
  simu.e_init          = calc_e(io.header.T_init)*pow2(io.header.a_initial/(io.header.boxsize/io.header.t_unit));
#else
  simu.gamma           = 0.0;
  simu.omegab          = 0.0;
  simu.omegaDM         = simu.omega0 - simu.omegab;
  simu.f_b             = simu.omegab/simu.omega0;
  simu.H_frac          = 0.0;
  simu.T_init          = 0.0;
  simu.B_init          = 0.0;
  simu.e_init          = 0.0;
#endif
  
  
    
  /* the particle details */
  simu.no_part    = io.no_part;
  simu.no_gas     = io.no_gas;  
  simu.no_stars   = io.no_stars;
  simu.no_species = io.header.no_species;
  simu.no_vpart   = io.header.no_vpart;
  simu.min_weight = io.header.min_weight;
  simu.max_weight = io.header.max_weight;
  simu.med_weight = io.header.med_weight;
  
  
  /* the grid details */
  simu.NGRID_DOM     = user.NGRID_DOM;
#ifndef FFT
  simu.NGRID_MIN     = MIN_L1DIM;
#else
  simu.NGRID_MIN     = user.NGRID_DOM;
#endif
  simu.NGRID_MAX     = MAX_L1DIM;
  simu.Nth_dom       = user.Nth_dom;
  simu.Nth_ref       = user.Nth_ref;
  simu.SHIFT         = ((double)0.5000000/(double) simu.NGRID_DOM);
#ifdef NP_LIMIT
  simu.np_limit      = TRUE;
#else
  simu.np_limit      = FALSE;
#endif


  /* some useful flags, though they are just used when writing the logfile! */
#ifdef MULTIMASS
#ifdef AHFmmfocus
  simu.mmfocus = 1;
#else
  simu.mmfocus = 0;
#endif
#else
  simu.mmfocus = 0;
#endif  
  
  simu.multi_mass       = io.header.multi_mass;
  simu.double_precision = io.header.double_precision;
  simu.hydro            = io.header.hydro;
  simu.magneto          = io.header.magneto;
  

#ifdef MOND
  simu.g0  = user.g0;
  simu.h0  = user.h0;
  
  /* transfer g0 to internal units */
  simu.g0 *= 3.08568025E12;   /* prepare for km/sec/Mpc */
  simu.g0 /= simu.h0 * simu.boxsize;
#endif
  

#ifdef LIGHTCONE
  simu.z_lightcone = user_data.lightcone_z;
#endif
  
  fprintf(stderr," <= finished initializing simu. structure\n");

}


/*=======================================================================
 * we simply copy the user provided parameters over to io. structure
 *=======================================================================*/
void init_io_structure(uparam user)
{
  char termfile_name[MAXSTRING];
  int  i;
  
  /*====================================================
   * io.
   *====================================================*/
  /* allocate memory for output filenames */
  io.dumpfile_name     = (char *)malloc((strlen(user.outfile_prefix)+1+8) *sizeof(char));
  io.logfile_name      = (char *)malloc((strlen(user.outfile_prefix)+1+3) *sizeof(char));  
  io.icfile_name       = (char *)malloc((strlen(user.icfile_name)+1)      *sizeof(char));
  io.outfile_prefix    = (char *)malloc((strlen(user.outfile_prefix)+1)   *sizeof(char));
#ifdef LIGHTCONE
  io.lightcone_prefix  = (char *)malloc((strlen(user.lightcone_prefix)+1) *sizeof(char));
#endif /* LIGHTCONE */
  
  /* copy filenames over from user structure */
  strcpy(io.icfile_name,      user.icfile_name);
  strcpy(io.outfile_prefix,   user.outfile_prefix);
#ifdef LIGHTCONE
  strcpy(io.lightcone_prefix, user.lightcone_prefix);
#endif /* LIGHTCONE */
  
  /* append "dumpfile" and "log" to outfile_prefix for the dump- and log-file, respectively  */
  strcpy(io.dumpfile_name, user.outfile_prefix);
  strcat(io.dumpfile_name, "dumpfile");
  strcpy(io.logfile_name,  user.outfile_prefix);
  strcat(io.logfile_name,  "log");
  
  /* initialize remainder of io-structure */
  io.out_dumps   = user.out_dumps;
  io.no_outputs  = user.no_outputs;
  io.a_out       = (double *) calloc(io.no_outputs, sizeof(double));
  for(i=0; i<io.no_outputs; i++)
    io.a_out[i] = 1./(1.+user.z_out[i]);
  
  
#ifdef LIGHTCONE
  {
    double lightcone_z;
    int lightcone_type;
    double sa,sb,sg;
    double dcx,dcy;
    
    lightcone_z    = user.lightcone_z;
    lightcone_type = user.lightcone_type;
    sa             = user.sa;
    sb             = user.sb;
    sg             = user.sg;
    dcx            = user.dcx;
    dcy            = user.dcy;
    
    io.conetype = lightcone_type;
    if(lightcone_type==-1) {
      sa*=D2R; sb*=D2R; sg*=D2R; dcx*=(D2R/2.); dcy*=(D2R/2.);
      io.sa=sa;
      io.sb=sb;
      io.sg=sg;
      io.dcx=dcx;
      io.dcy=dcy;
      /* these coefficients are used to check coords in in_patch */
      dcx=sin(dcx);
      io.xcoef=dcx*dcx;
      dcy=sin(dcy);
      io.ycoef=dcy*dcy;
      /* patch geometry: patch pyramid edge vectors */
      /* note: normalization is not necessary       */
      dcx=tan(io.dcx); dcy=tan(io.dcy);
      io.kpatch[0][0]= dcx; io.kpatch[0][1]= dcy; io.kpatch[0][2]=1.;
      io.kpatch[1][0]= dcx; io.kpatch[1][1]= -dcy; io.kpatch[1][2]=1.;
      io.kpatch[2][0]= -dcx; io.kpatch[2][1]= dcy; io.kpatch[2][2]=1.;
      io.kpatch[3][0]= -dcx; io.kpatch[3][1]= -dcy; io.kpatch[3][2]=1.;
      /* normals to patch pyramid faces */
      io.npatch[0][0]=0.; io.npatch[0][1]=cos(io.dcy); 
      io.npatch[0][2]= -sin(io.dcy);
      io.npatch[1][0]=0.; io.npatch[1][1]= -cos(io.dcy); 
      io.npatch[1][2]= -sin(io.dcy);
      io.npatch[2][0]=cos(io.dcx); io.npatch[2][1]=0; 
      io.npatch[2][2]= -sin(io.dcx);
      io.npatch[3][0]= -cos(io.dcx); io.npatch[3][1]=0; 
      io.npatch[3][2]= -sin(io.dcx);
      /* the transform matrix R: */
      /*	sb=sb-PI/2.;	// our beta is measured from the x-y plane */
      sb=PI/2.-sb;		/* first version was wrong */
      /* the Euler alpha-beta-gamma definition for the rotation matrix */
      io.R[0][0]= cos(sa)*cos(sb)*cos(sg)-sin(sa)*sin(sg);
      io.R[0][1]= sin(sa)*cos(sb)*cos(sg)+cos(sa)*sin(sg);
      io.R[0][2]= -sin(sb)*cos(sg);
      io.R[1][0]= -cos(sa)*cos(sb)*sin(sg)-sin(sa)*cos(sg);
      io.R[1][1]= -sin(sa)*cos(sb)*sin(sg)+cos(sa)*cos(sg);
      io.R[1][2]= sin(sb)*sin(sg);
      io.R[2][0]= cos(sa)*sin(sb);
      io.R[2][1]= sin(sa)*sin(sb);
      io.R[2][2]= cos(sb);
    }
  }
#endif /* LIGHTCONE */
  
#if (defined ART && defined HYDRO)
  io.header.omegab = user.omegab;
  io.header.gamma  = user.gamma;
  io.header.H_frac = user.H_frac;
  io.header.T_init = user.T_init;
  io.header.B_init = user.B_init;
#endif
}

/*========================================================================
 * determine...
 *  no_species              no. of different particle species
 *  no_vpart                total mass in box
 *  min_weight              min. particle mass
 *  max_weight              max. particle mass
 *  most_common_weight      most common particle mass
 *
 * -> initialize...
 *  io.header.no_species
 *  io.header.no_vpart
 *  io.header.min_weight
 *  io.header.max_weight
 *  io.header.med_weight
 *  ...and return the most common particle weight
 *========================================================================*/
double init_header_masses()
{
   double  *wspecies, cur_weight;
   long    *nspecies, nmax, imax, ipart;
   int      no_species;
   int      ispecies;
   boolean  old_species;
   double   no_vpart;
   double   min_weight;
   double   max_weight;
   double   K, T, C; // for Kahan summation
   partptr  cur_part;
   double   most_common_weight, med_weight, med_norm;
   double   npart, inv_mean_dens, pmass;
   
   fprintf(stderr,"     init_header_masses():\n");

#ifdef MULTIMASS
   no_species             = 1;
   wspecies               = (double*)calloc(no_species, sizeof(double));
   nspecies               = (long*)  calloc(no_species, sizeof(long));
   
   /* we determine particle weights in units of first particle mass */
   pmass                  = io.fst_part->weight;
   min_weight             = 1.0;
   max_weight             = 1.0;
   no_vpart               = 1.0;
   wspecies[no_species-1] = 1.0;
   nspecies[no_species-1] = 1;
   
   /* for Kahan summation of no_vpart */
   C = 0.0;
   
   for(cur_part=io.fst_part+1; cur_part<io.fst_part+io.no_part; cur_part++)
     {
      /* compare particle weight in units of first particle mass */
      cur_weight  = cur_part->weight/pmass;
             
      /* no_vpart += cur_weight ala Kahan summation */
      K = cur_weight-C;
      T = no_vpart + K;
      C = (T-no_vpart) - K;
      no_vpart = T;

      old_species = FALSE;
      for(ispecies=0; ispecies<no_species; ispecies++)
        {
         if(fabs(cur_weight - wspecies[ispecies])/wspecies[ispecies] < 0.1)
           {
            old_species         = TRUE;
            nspecies[ispecies] += 1;
           }
        }
      
      if(old_species == FALSE)
        {
         no_species++;
          
         if(cur_weight > max_weight) max_weight = cur_weight;
         if(cur_weight < min_weight) min_weight = cur_weight;
         
         wspecies               = (double*) realloc(wspecies, no_species*sizeof(double));
         nspecies               = (long*)   realloc(nspecies, no_species*sizeof(long));

         wspecies[no_species-1] = cur_weight;
         nspecies[no_species-1] = 1;
        }
     }
   
   fprintf(stderr,"     -> the following species have been found...\n");
   med_weight = 0.0;
   med_norm = 0.0;
   for(ispecies=0; ispecies<no_species; ispecies++)
     {
       fprintf(stderr,"       %12d      %16.8g [Msun/h]     %16ld\n",ispecies,wspecies[ispecies]*pmass,nspecies[ispecies]);
       med_weight += nspecies[ispecies]*wspecies[ispecies];
       med_norm   += (double)nspecies[ispecies];
     }
   med_weight /= med_norm;
   fprintf(stderr,"\n");
   
   
   /* determine most common particle weight */
   imax = 0;
   nmax = nspecies[0];
   for(ispecies=1; ispecies<no_species; ispecies++)
     {
      if(nspecies[ispecies] > nmax)
        {
         nmax = nspecies[ispecies];
         imax = ispecies;
        }
     }
   most_common_weight = wspecies[imax];

   free(wspecies); 
   free(nspecies);    
   
#else /* MULTIMASS */
   
#ifdef NO_EXPANSION
   fprintf(stderr,"     YOU ARE USING A NON-COSMOLOGICAL SETTING WITHOUT PROVIDING PARTICLE MASSES\n");
   fprintf(stderr,"      -> not implemented yet...exiting\n");
   exit(0);
#else /* NO_EXPANSION */
   /* use cosmological background density to determine pmass */
   npart                  = (double)io.no_part;
   inv_mean_dens          = pow3(io.header.boxsize)/npart;
   pmass                  = io.header.omega0*rhoc0*inv_mean_dens;
   
   /* these values are in internal units */
   no_species             = 1;
   no_vpart               = npart;
   most_common_weight     = 1.0;
   min_weight             = 1.0;
   max_weight             = 1.0;
   med_weight             = 1.0;
   most_common_weight     = 1.0;

   fprintf(stderr,"     you are using a cosmological setting without providing particle masses\n");
   fprintf(stderr,"      -> will use pmass = %16.8g as particle mass\n\n",pmass);
#endif /* NO_EXPANSION */
   
#endif /* MULTIMASS */
   
   
   
   
   
   /* dump no_vpart, no_species, min_weight, and max_weight */
   fprintf(stderr,"     no_species          = %16d\n",             no_species);
   fprintf(stderr,"     pmass               = %16.8g  Msun/h (%16.8g)\n",               pmass, io.header.pmass);
   fprintf(stderr,"     no_vpart            = %16.8g  Msun/h\n",   no_vpart           * pmass);
   fprintf(stderr,"     most_common_weight  = %16.8g  Msun/h\n\n", most_common_weight * pmass);
   fprintf(stderr,"     min_weight          = %16.8g  Msun/h\n",   min_weight         * pmass);
   fprintf(stderr,"     max_weight          = %16.8g  Msun/h\n\n", max_weight         * pmass);
   fprintf(stderr,"     med_weight          = %16.8g  Msun/h\n\n", med_weight         * pmass);
   
   /* store newly calculated values */
   io.header.no_species = no_species;
   io.header.no_vpart   = no_vpart   * pmass;
   io.header.min_weight = min_weight * pmass;    // these values are in Msun/h
   io.header.max_weight = max_weight * pmass;
  
   // TODO: what is the right choice for med_weight anyways?!
   io.header.med_weight = most_common_weight * pmass; // med_weight * pmass;

   most_common_weight   = most_common_weight * pmass;    // this will be used as the mass unit in AMIGA

  return (most_common_weight);
}

/*========================================================================
 * convert pos[], mom[], and weight to internal units
 *========================================================================*/
void ic_unit_conversion()
{
   partptr       cur_part;
   long          ipart;
   double        x_fac, v_fac, m_fac;
   double        xmin, xmax, ymin, ymax, zmin, zmax;
   double        rho_mean;
   
   fprintf(stderr,"\n=================================================================\n");
   fprintf(stderr,"                    ic_unit_conversion()\n");
   fprintf(stderr,"=================================================================\n");

   /*=====================================================================
    * 1. CHECK POSITIONS RANGE
    *=====================================================================*/
   fprintf(stderr," 1. position range consistency check:\n");
   fprintf(stderr," ------------------------------------\n");
   xmax = -1E40; ymax = -1E40; zmax = -1E40;
   xmin = +1E40; ymin = +1E40; zmin = +1E40;
   
   /* there is no MIN,MAX reduction for OpenMP in C :-( */
   for(cur_part=io.fst_part; cur_part<io.fst_part+io.no_part; cur_part++)
     {
      if(cur_part->pos[X] > xmax) xmax = cur_part->pos[X];
      if(cur_part->pos[Y] > ymax) ymax = cur_part->pos[Y];
      if(cur_part->pos[Z] > zmax) zmax = cur_part->pos[Z];
      
      if(cur_part->pos[X] < xmin) xmin = cur_part->pos[X];
      if(cur_part->pos[Y] < ymin) ymin = cur_part->pos[Y];
      if(cur_part->pos[Z] < zmin) zmin = cur_part->pos[Z];
     }
   
   if(xmin >= 0.0               && ymin >= 0.0               && zmin >= 0.0 && 
      xmax <= io.header.boxsize && ymax <= io.header.boxsize && zmax <= io.header.boxsize)
     {
      fprintf(stderr,"   -> everything's fine!\n");
     }
   
   /* negative coordinates? */
   if(xmin < 0.0 || ymin < 0.0 || zmin < 0.0)
     {
      fprintf(stderr,"   -> negative coordinates:\n");
      fprintf(stderr,"      boxsize = %16.8g\n", io.header.boxsize);
      fprintf(stderr,"      min     = %16.8g    %16.8g    %16.8g\n",xmin,ymin,zmin);
      fprintf(stderr,"      max     = %16.8g    %16.8g    %16.8g\n",xmax,ymax,zmax);

      /* is it possible to simply shift the positions? */
      if(xmax-xmin <= io.header.boxsize || ymax-ymin <= io.header.boxsize || zmax-zmin <= io.header.boxsize)
        {
         fprintf(stderr,"     -> shifting coordinates to [0,%g]\n",io.header.boxsize);
#ifdef WITH_OPENMP
#pragma omp parallel private(cur_part) shared(xmin, xmax, ymin, ymax, zmin, zmax, io)
#pragma omp for schedule(static)
#endif
         for(ipart=0; ipart<io.no_part; ipart++)
           {
            cur_part         = io.fst_part + ipart;
            cur_part->pos[X] = fmod(cur_part->pos[X]+fabs(xmin), io.header.boxsize);
            cur_part->pos[Y] = fmod(cur_part->pos[Y]+fabs(ymin), io.header.boxsize);
            cur_part->pos[Z] = fmod(cur_part->pos[Z]+fabs(zmin), io.header.boxsize);
           }
        }
      else
        {
         fprintf(stderr,"     -> and coordinates exceeding allowed range:\n");
         fprintf(stderr,"        boxsize = %16.8g\n", io.header.boxsize);
         fprintf(stderr,"        min     = %16.8g    %16.8g    %16.8g\n",xmin,ymin,zmin);
         fprintf(stderr,"        max     = %16.8g    %16.8g    %16.8g\n",xmax,ymax,zmax);

         fprintf(stderr,"        not fixed yet...exiting\n");
         exit(0);
        }
     }
   
   /* coordinates exceeding boxsize limit? */
   if(xmax > io.header.boxsize || ymax > io.header.boxsize || zmax > io.header.boxsize)
     {
      fprintf(stderr,"   -> coordinates exceeding allowed range:\n");
      fprintf(stderr,"      boxsize = %16.8g\n", io.header.boxsize);
      fprintf(stderr,"      min     = %16.8g    %16.8g    %16.8g\n",xmin,ymin,zmin);
      fprintf(stderr,"      max     = %16.8g    %16.8g    %16.8g\n",xmax,ymax,zmax);
      
      fprintf(stderr,"      not fixed yet...exiting\n");
      exit(0);
     }
   
   
   
   
   
   /*==================================================================
    * 2. DETERMINE MOST APPROPRIATE MASS UNIT
    *==================================================================*/
   fprintf(stderr,"\n 2. determining most appropriate mass unit:\n");
   fprintf(stderr," ------------------------------------------\n");
   io.header.pmass = init_header_masses();
#ifdef GADGET_PMASS
   io.header.pmass = GADGET_MUNIT;
   fprintf(stderr,"   -> will use GADGET mass unit %g Msun/h as mass unit for AMIGA\n",io.header.pmass);
#endif
#ifdef MANUAL_PMASS
   io.header.pmass = 42570.70209705;   // set the mass unit to whatever you fancy...
#endif
   
   /* AMIGA expects these values to be in internal units! */
   io.header.no_vpart   /= io.header.pmass;
   io.header.min_weight /= io.header.pmass;
   io.header.max_weight /= io.header.pmass;
   io.header.med_weight /= io.header.pmass;
   
   fprintf(stderr,"   -> will use %g Msun/h as mass unit for AMIGA\n",io.header.pmass);

   
   
   /*==================================================================
    * 3. CHOOSE TIME UNIT
    *==================================================================*/
   fprintf(stderr,"\n 3. setting the time unit:\n");
   fprintf(stderr," -------------------------\n");
#ifdef NO_EXPANSION
   rho_mean         = io.header.no_vpart*io.header.pmass/pow3(io.header.boxsize);
   io.header.t_unit = 1./(4.*PI*Grav*rho_mean);
   io.header.t_unit = sqrt(io.header.t_unit);
   
   fprintf(stderr,"   -> non-cosmological setup, will use %16.8g as time unit for AMIGA\n",io.header.t_unit);
#else
   io.header.t_unit = 1./H0;

   fprintf(stderr,"   -> cosmological setup, will use 1/H0 = %g as time unit for AMIGA\n",io.header.t_unit);
#endif
   
   
   
   
   
   /*==================================================================
    *                        UNIT CONVERSION
    *==================================================================*/
   x_fac = 1./io.header.boxsize;
   v_fac = io.header.a_current*io.header.t_unit/io.header.boxsize;
   m_fac = 1./io.header.pmass;
   
#ifdef WITH_OPENMP
#pragma omp parallel private(cur_part) shared(x_fac, v_fac, m_fac, io)
#pragma omp for schedule(static)
#endif
   for(ipart=0; ipart<io.no_part; ipart++)
     {
      cur_part          = io.fst_part + ipart;
      cur_part->pos[X]  = f1mod(cur_part->pos[X]*x_fac + 1.0, 1.0);
      cur_part->pos[Y]  = f1mod(cur_part->pos[Y]*x_fac + 1.0, 1.0);
      cur_part->pos[Z]  = f1mod(cur_part->pos[Z]*x_fac + 1.0, 1.0);
      cur_part->mom[X] *= v_fac;
      cur_part->mom[Y] *= v_fac;
      cur_part->mom[Z] *= v_fac;
#ifdef MULTIMASS
      cur_part->weight *= m_fac;
#endif
     }
   fprintf(stderr,"=================================================================\n");
   fprintf(stderr,"                finished ic_unit_conversion()\n");
   fprintf(stderr,"=================================================================\n");
   
}

/*========================================================================
 * here we simply check whether the DEFINFLAGS in combination with the
 * io.header parameters are actually meaningful...
 *========================================================================*/
void sanity_check()
{
   double pmass;
   
   fprintf(stderr," => performing sanity check of header information:\n");
   
   /* dump AMIGA header (*before* reading particles) */
   fprintf(stderr,"  io.header (as read from input file):\n");
   fprintf(stderr,"  ------------------------------------\n");
   fprintf(stderr,"   %s\n",io.header.header);
   fprintf(stderr,"   header.multi_mass            = %d\n",io.header.multi_mass);
   fprintf(stderr,"   header.double_precision      = %d\n",io.header.double_precision);
   fprintf(stderr,"   header.hydro                 = %d\n",io.header.hydro);
   fprintf(stderr,"   header.magneto               = %d\n",io.header.magneto);
   fprintf(stderr,"   header.no_part               = %ld\n",io.header.no_part);
   fprintf(stderr,"   header.no_species            = %ld\n",io.header.no_species);
   fprintf(stderr,"   header.min_weight            = %g\n",io.header.min_weight);
   fprintf(stderr,"   header.max_weight            = %g\n",io.header.max_weight);
   fprintf(stderr,"   header.med_weight            = %g\n",io.header.med_weight);
   fprintf(stderr,"   header.no_vpart              = %g\n",io.header.no_vpart);
   fprintf(stderr,"   header.timestep              = %g\n",io.header.timestep);
   fprintf(stderr,"   header.no_timestep           = %d\n",io.header.no_timestep);
   fprintf(stderr,"   header.t_unit                = %g\n",io.header.t_unit);
   fprintf(stderr,"   header.pmass                 = %g\n",io.header.pmass);
   fprintf(stderr,"   header.boxsize               = %g\n",io.header.boxsize);
   fprintf(stderr,"   header.omega0                = %g\n",io.header.omega0);
   fprintf(stderr,"   header.omegab                = %g\n",io.header.omegab);
   fprintf(stderr,"   header.lambda0               = %g\n",io.header.lambda0);   
   fprintf(stderr,"   header.gamma                 = %g\n",io.header.gamma);
   fprintf(stderr,"   header.H_frac                = %g\n",io.header.H_frac);
   fprintf(stderr,"   header.T_init                = %g\n",io.header.T_init);
   fprintf(stderr,"   header.cur_reflevel          = %g\n",io.header.cur_reflevel);
   fprintf(stderr,"   header.cur_frcres            = %g\n",io.header.cur_frcres);
   fprintf(stderr,"   header.a_initial             = %g\n",io.header.a_initial);
   fprintf(stderr,"   header.a_current             = %g\n",io.header.a_current);
   fprintf(stderr,"   header.K_initial             = %g\n",io.header.K_initial);
   fprintf(stderr,"   header.K_current             = %g\n",io.header.K_current);
   fprintf(stderr,"   header.U_initial             = %g\n",io.header.U_initial);
   fprintf(stderr,"   header.U_current             = %g\n",io.header.U_current);
   fprintf(stderr,"   header.Eintegral             = %g\n",io.header.Eintegral);
   fprintf(stderr,"   header.Econst                = %g\n",io.header.Econst);
   fprintf(stderr,"   header.paramNSTEPS           = %g\n",io.header.paramNSTEPS);
   fprintf(stderr,"   header.paramNGRID_DOM        = %g\n",io.header.paramNGRID_DOM);   
   fprintf(stderr,"   header.paramNth_dom          = %g\n",io.header.paramNth_dom);
   fprintf(stderr,"   header.paramNth_ref          = %g\n",io.header.paramNth_ref);
   fprintf(stderr,"   header.paramE_UPDATE         = %g\n",io.header.paramE_UPDATE);
   fprintf(stderr,"   header.paramCELLFRAC_MAX     = %g\n",io.header.paramCELLFRAC_MAX);
   fprintf(stderr,"   header.paramCELLFRAC_MIN     = %g\n",io.header.paramCELLFRAC_MIN);
   fprintf(stderr,"   header.paramCA_CRIT          = %g\n",io.header.paramCA_CRIT);
   fprintf(stderr,"   header.paramMAX_L1DIM        = %g\n",io.header.paramMAX_L1DIM);
   fprintf(stderr,"   header.paramDOMSWEEPS        = %g\n",io.header.paramDOMSWEEPS);
   fprintf(stderr,"   header.paramREFSWEEPS        = %g\n",io.header.paramREFSWEEPS);
   fprintf(stderr,"   header.paramAHF_MINPART      = %g\n",io.header.paramAHF_MINPART);
   fprintf(stderr,"   header.paramAHF_VTUNE        = %g\n",io.header.paramAHF_VTUNE);
   fprintf(stderr,"   header.paramAHF_RISE         = %g\n",io.header.paramAHF_RISE);
   fprintf(stderr,"   header.paramAHF_SLOPE        = %g\n",io.header.paramAHF_SLOPE);
   fprintf(stderr,"   header.paramAHF_MAXNRISE     = %g\n",io.header.paramAHF_MAXNRISE);
   
   /*=========
    *  UNITS
    *=========*/
   /* check mass unit */
   if(fabs(io.header.pmass) < ZERO)
     {
      fprintf(stderr," o mass unit not set  ");
      
#ifdef NO_EXPANSION
      fprintf(stderr,"-> non-cosmological setup ... exiting!\n");
      exit(0);
#else
      /* we brute force use the cosmological setting! */
      io.header.pmass = io.header.omega0*rhoc0*pow3(io.header.boxsize)/io.header.no_vpart;
      fprintf(stderr,"-> using cosmological value pmass = %g\n",io.header.pmass);
#endif
     }
   
   /* check time unit */
   if(fabs(io.header.t_unit) < ZERO)
     {
      fprintf(stderr," o time unit not set  ");
      
#ifdef NO_EXPANSION
      rho_mean         = io.header.no_vpart*io.header.pmass/pow3(io.header.boxsize);
      io.header.t_unit = 1./(4.*PI*Grav*rho_mean);
      io.header.t_unit = sqrt(io.header.t_unit);
      fprintf(stderr,"-> using non-cosmological value t_unit = %g\n",io.header.t_unit);
#else
      /* we brute force use the cosmological setting! */
      io.header.t_unit = 1./H0;   
      fprintf(stderr,"-> using cosmological value 1/H0  = %g\n",io.header.t_unit);
#endif
     }
   

   /*====================
    *  io.header masses
    *====================*/
   if(fabs(io.header.no_vpart)   < ZERO || 
      fabs(io.header.no_species) < ZERO || 
      fabs(io.header.min_weight) < ZERO || 
      fabs(io.header.max_weight) < ZERO ||
      fabs(io.header.med_weight) < ZERO)
     {
      fprintf(stderr," o io.header masses not initialized ...\n");
      
      init_header_masses();   
      
       /* init_header_masses() returns the values in physical units Msun/h */
      io.header.no_vpart   /= io.header.pmass;
      io.header.min_weight /= io.header.pmass;
      io.header.max_weight /= io.header.pmass;
      io.header.med_weight /= io.header.pmass;
     }
   

   
   /*=============
    *  MULTIMASS
    *=============*/
#ifndef MULTIMASS
   if(io.header.multi_mass == 1)
     {
      fprintf(stderr,"\n=================================================================\n");      
      fprintf(stderr,"      your are trying to run a multi-mass simulation:\n");      
      fprintf(stderr,"          please recompile AMIGA with -DMULTIMASS\n");
      fprintf(stderr,"=================================================================\n");  
      exit(0);
     }
#endif
#ifdef MULTIMASS
   if(io.header.multi_mass == 0)
     {
      fprintf(stderr,"\n=================================================================\n");      
      fprintf(stderr,"     your are trying to run a sinlge-mass simulation:\n");      
      fprintf(stderr,"        please recompile AMIGA without -DMULTIMASS\n");
      fprintf(stderr,"=================================================================\n");      
      exit(0);
     }
#endif
   
   
   /*============
    *   HYDRO
    *============*/
#ifndef HYDRO
   if(io.header.hydro == 1)
     {
      fprintf(stderr,"\n===============================================================================\n");      
      fprintf(stderr," o your are running a hydro simulation without the hydro-part switched on!\n");      
      fprintf(stderr,"   => is this intentional? ...better exit and let you check!\n");
      fprintf(stderr,"===============================================================================\n");  
      exit(0);
     }  
#endif
#ifdef HYDRO
   if(io.header.hydro == 0)
     {
      if(io.header.omegab > ZERO)
        {
         fprintf(stderr,"\n o you are running a hydro simulation with a pure DM simulation as input\n");
         fprintf(stderr,"     => will use omegab=%8.4g and split a DM particle into DM+gas?\n",io.header.omegab);      
        }
      else
        {
         fprintf(stderr,"\n===============================================================================\n");      
         fprintf(stderr,"   o you are running a hydro simulation with a pure DM simulation as input\n");
         fprintf(stderr,"     => as you also did not provide a credible omegab value AMIGA terminates now!\n");      
         fprintf(stderr,"===============================================================================\n");  
         exit(0);
        }
     }  
   
   io.header.hydro = 1;
#endif
   
   
   /*============
    *    MHD
    *============*/
#ifndef MHD
   if(io.header.magneto == 1)
     {
      fprintf(stderr,"\n===============================================================================\n");      
      fprintf(stderr,"   you are running a MHD simulation without the MHD-part switched on!\n");      
      fprintf(stderr,"   => is this intentional? ...better exit and let you check!\n");
      fprintf(stderr,"===============================================================================\n");  
      exit(0);
     }  
#endif
   
   
   /*============
    *   DOUBLE
    *============*/
#ifdef DOUBLE
   if(io.header.double_precision != 1)
     {
      fprintf(stderr,"\n input file is single precision but simulation will be run in double precision\n");
      fprintf(stderr,"   => will upcast\n");      
     }
   
   io.header.double_precision = 1;
#else
   if(io.header.double_precision == 1)
     {
      fprintf(stderr,"\n input file is double precision but simulation will be run in single precision\n");
      fprintf(stderr,"   => will downcast\n");      
     }
   
   io.header.double_precision = 0;
#endif
   
   fprintf(stderr," <= finished sanity_check()\n");
}

/*========================================================================
 * simply read the filenames and parameters from the input parameter file
 *========================================================================*/
void get_user_data(FILE *startrun_dat, uparamptr user)
{
  int    i;
  double z_out;
  
  fprintf(stderr,"please give name of file with initial conditions: ");
  fscanf(startrun_dat, "%s", user->icfile_name);
  fprintf(stderr,"%s\n",user->icfile_name);
  fprintf(stderr,"please give prefix for output file names:         ");
  fscanf(startrun_dat, "%s", user->outfile_prefix);
  fprintf(stderr,"%s\n",user->outfile_prefix);
#ifdef LIGHTCONE
  fprintf(stderr,"please give prefix for lightcone file names:         ");
  fscanf(startrun_dat, "%s", user->lightcone_prefix);
  fprintf(stderr,"%s\n",user->lightcone_prefix);
#endif /* LIGHTCONE */
  
  fprintf(stderr,"please give number of domain grid cells (1D):     ");
  fscanf(startrun_dat,"%d", &(user->NGRID_DOM));
  fprintf(stderr,"%d\n",user->NGRID_DOM);
  fprintf(stderr,"please give Nth for domain grid:                  ");
  fscanf(startrun_dat,"%lf", &user->Nth_dom);
  fprintf(stderr,"%g\n",user->Nth_dom);
  fprintf(stderr,"please give Nth for refinements:                  ");
  fscanf(startrun_dat,"%lf", &user->Nth_ref);
  fprintf(stderr,"%g\n",user->Nth_ref);
  fprintf(stderr,"please give final redshift:                       ");
  fscanf(startrun_dat,"%lf", &user->final_z);
  fprintf(stderr,"%g\n",user->final_z);
#ifdef LIGHTCONE
  fprintf(stderr,"please give lightcone redshift limit:             ");
  fscanf(startrun_dat,"%lf", &user->lightcone_z);
  fprintf(stderr,"%g\n",user->lightcone_z);
  fprintf(stderr,"please give lightcone type:                       ");
  fscanf(startrun_dat,"%d", &user->lightcone_type);
  fprintf(stderr,"%d\n",user->lightcone_type);
  if(lightcone_type==-1) {
    fprintf(stderr,"please give patch orientation angles (degrees): ");
    fscanf(startrun_dat,"%lf %lf %lf", &user->sa, &user->sb, &user->sg);
    fprintf(stderr,"%g %g %g\n",user->sa,user->sb,user->sg);
    fprintf(stderr,"please give patch sizes (degrees):              ");
    fscanf(startrun_dat,"%lf %lf", &user->dcx, &user->dcy);
    fprintf(stderr,"%g %g\n",user->dcx,user->dcy);
  }
#endif /* LIGHTCONE */
#ifdef TIPSY
  /* NOTE: all the TIPSY relevant parameters are global variabels defined in tdef.h by Justin Read */
  fprintf(stderr,"please give box size [Mpc/h]:                     ");
  fscanf(startrun_dat,"%lf", &tipsy_boxsize);
  fprintf(stderr,"%lf\n", tipsy_boxsize);
  fprintf(stderr,"please give omega0:                               ");
  fscanf(startrun_dat,"%lf", &tipsy_omega0);
  fprintf(stderr,"%lf\n", tipsy_omega0);
  fprintf(stderr,"please give lambda0:                              ");
  fscanf(startrun_dat,"%lf", &tipsy_lambda0);
  fprintf(stderr,"%lf\n", tipsy_lambda0);
  fprintf(stderr,"please give initial redshift:                     ");
  fscanf(startrun_dat,"%lf", &tipsy_initalz);
  fprintf(stderr,"%lf\n", tipsy_initalz);
  fprintf(stderr,"please give current timestep no:                  ");
  fscanf(startrun_dat,"%lf", &tipsy_currentimeno);
  fprintf(stderr,"%lf\n", tipsy_currentimeno);
#endif /* TIPSY */
  
  fprintf(stderr,"please give dump file frequency:                  ");
  fscanf(startrun_dat,"%d",  &user->out_dumps);
  fprintf(stderr,"%d\n",user->out_dumps);
  fprintf(stderr,"please give total number of output files:         ");
  fscanf(startrun_dat,"%d",  &user->no_outputs);
  fprintf(stderr,"%d\n",user->no_outputs);
  
#ifndef ISOLATED /* for -DISOLATED the output times are *not* provided by the user! */
  user->z_out = (double *) calloc(user->no_outputs, sizeof(double));
  for(i = 0; i < user->no_outputs; i++)
    {
      fprintf(stderr,"please give redshift for %5d. output file:      ",i+1);
      fscanf(startrun_dat,"%lf", &z_out);
      fprintf(stderr,"%g\n",z_out);
      user->z_out[i] = z_out;
    }
#endif
  
#if (defined ART && defined HYDRO)
  fprintf(stderr,"please give omegab:                               ");
  fscanf(startrun_dat,"%lf", &user->omegab);
  fprintf(stderr,"%g\n",user->omegab);
  fprintf(stderr,"please give gamma:                                ");
  fscanf(startrun_dat,"%lf", &user->gamma);
  fprintf(stderr,"%g\n",user->gamma);
  fprintf(stderr,"please give H_frac:                               ");
  fscanf(startrun_dat,"%lf", &user->H_frac);
  fprintf(stderr,"%g\n",user->H_frac);
  fprintf(stderr,"please give T_init:                               ");
  fscanf(startrun_dat,"%lf", &user->T_init);
  fprintf(stderr,"%g\n",user->T_init);
  fprintf(stderr,"please give B_init:                               ");
  fscanf(startrun_dat,"%lf", &user->B_init);
  fprintf(stderr,"%g\n",user->B_init);
#endif
  
#ifdef MOND /* for -DMOND two additional parameters need to be provided */
  fprintf(stderr,"please give MOND acceleration g0:                 ");
  fscanf(startrun_dat,"%lf", &user->g0);
  fprintf(stderr,"%g\n",user->g0);
  fprintf(stderr,"please give Hubble parameter H0:                  ");
  fscanf(startrun_dat,"%lf", &user->H0);
  fprintf(stderr,"%g\n",user->H0);
#endif
  
  
}


#if (defined NEWSTARTRUN)
/**
 * cmp_sfckey_part: compares the sfc keys of two particles, used
 * for qsort 
 */
extern int
cmp_sfckey_part(const void *p1, const void *p2) 
{
	if (((partptr)p1)->sfckey < ((partptr)p2)->sfckey)
		return -1;

	if (((partptr)p1)->sfckey > ((partptr)p2)->sfckey)
		return 1;

	return 0;
}
#endif
