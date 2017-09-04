/*==========================================================================
 *
 * This is a re-write from the bottom up of AHF!!
 * In this file I have merged _centres with _halos
 * 
 * 08.02.2009: small adaptations for working with HaloTracker as well
 *
 *	
 *==========================================================================*/


/***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* the important definitions have to be included first */
#include "../common.h"
#include "../param.h"
#include "../tdef.h"

/* ...and now for the actual includes */
#include "ahf_halos.h"
#include "../libutility/utility.h"
#include "../libgrids/grids.h"
#include "../libparticles/particles.h"

#ifdef NEWSTARTRUN
#	include "ahf_halos_sfc.h"
#endif

static long dummy_variable;
void dummy_ahf_halos(void)
{
}

#if (defined AHF || defined TRACKER) 

SRINDEX *spatialRefIndex;
int      totnumIsoRef;       /* Number of spatially isolated refinements */
int     *numIsoRef;          /* Number of spatially isolated refinements */


static int      numHalos;
static HALO    *halos;


static int      removecount = 0;
static int      movecount   = 0;
static int      numDensZero = 0;
static SRINDEX *densZero;

/***************************************************************************/
/* General */
typedef struct {
  double x,y,z,norm;
} XYZnorm;

typedef struct {
  int refLevel, isoRefIndex;
} INDEX;

/*********************  ???? TO GO IN ???? **********************/
/*********************  ???? TO GO IN ???? **********************/
/*********************  ???? TO GO IN ???? **********************/

typedef struct {
  double 	eigenVectors[3][3];
  double	eigenValues[3];
  double	e1,e2;
  double	T;
} INERTIA;

typedef struct {
  XYZ		L;
  XYZ		pos;
  double	dens;
  int		colour;
} GRIDATACOL;

/*********************  ???? TO GO IN ???? **********************/
/*********************  ???? TO GO IN ???? **********************/
/*********************  ???? TO GO IN ???? **********************/

/***************************************************************************/
/* Spatial Refinement */
struct spatREF {
  int		      colour;	      /* The colour of the isolated refinement */
  int		      numNodes;	   /* The number of nodes in the refinement */
  double	      volume;		   /* Volume of the refinement */
  double	      maxDens;	      /* Maximum density of node on the refinement */
  
  XYZ		      boundRefDiv;   /* The divider to calculate the min-max */
  
  XYZ            centre;           /* The actual centre to be used as halo centre           */
  XYZnorm        centreCMpart;     /* The centre of mass of particles on spatialRef[][]     */
  XYZnorm        centreGEOM;       /* The geometric centre of spatialRef[][]                */
  XYZnorm        centreDens;		   /* The density weighted centre of spatialRef[][]         */
  XYZnorm        centrePot;        /* The potential weighted centre of spatialRef[][]       */
  MINMAX         x,y,z;            /* the extent of the refinement		*/
  
  long unsigned	 numParts;      /* Number of particles in the refinement */
  partptr        ll;            /* head of particle linked-list         */
  
  int            numSubStruct;	/* Number of sub halos */
  INDEX	        *subStruct;		/* Indexes of the substructure halos */
  
  int	         numParDom;		/* Number of partent domains */
  INDEX	        *parDom;		   /* Indexes of the parent domains */
  
  INDEX	         daughter;		/* Indexes of the daughter refinemnt */
  
  int            haloIndex;		/* Which halo do I belong to ? */
  
  double	      closeRefDist;	/* Distance to the closest refinement */
  
  /*  GRIDATA	*nodedata;	 The node data ------ not sure if will include?*/
  /*  INERTIA	inertia;	 unity inertia tensor			*/
  /*  INERTIA	winertia;	 The Weighted inertia tensor		*/
  /*  double	dx,dy,dz;	 Length of the REf in the x,y,z dir	*/
  /**/
};
typedef struct spatREF SPATIALREF;

/***************************************************************************/
/* General functions */
int 	nodeCompare (const void *,const void *);
int 	refCompare  (const void *,const void *);
int 	haloCompare (const void *,const void *);

/***************************************************************************/
/* Specific functions */

int	  RefCentre              (gridls*, int, SPATIALREF**);
int	  analyseRef             (int, SPATIALREF**);
int	  spatialRef2halos       (int, SPATIALREF**);
int	  ConstructHalos         (void);
double  get_haloedge           (double*, double*, int, int );
void    return_parts_to_host   (int, int);
void    safe_return_particles  (int);
void    return_particles       (int);
void    prepare_haloes         (void);
void    prepare_particles      (void);
void    gather_parts_from_host (int, int);
void    gather_hostParts       (int);
void    prune_host_halos       (void);



/***************************************************************************/
/* Arrays */
static double *gridl1dim;

/***************************************************************************/
/* Variables */

/***************************************************************************/
int           compare ( struct particle *, struct particle *, XYZ *);
/* Inertia tensor */
static double 	     matrix[NDIM][NDIM]; /* The input inertia tensor 	*/
static double 	     eigenv[NDIM][NDIM]; /* The calculates eigenvectors	*/


#ifdef AHFbinary
void ahf_write_profiles(char *prefix,
                        HALO *halos,
                        unsigned long *idx,
                        int numHalos);
void ahf_write_particles(char *prefix,
                         HALO *halos,
                         unsigned long *idx,
                         int numHalos);
void ahf_write_halos(char *prefix,
                     HALO *halos,
                     unsigned long *idx,
                     int numHalos);
void ahf_write_open_files(FILE **f,
                          FILE **f_info,
                          char *prefix,
                          char *suffix);
#endif

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

/* conversion factors etc. */
#ifdef TRACKER
/* conversion factors declared globally in HaloTracker */
extern double  r_fac, x_fac, v_fac, m_fac, rho_fac, phi_fac, Hubble;
#else
static double  r_fac, x_fac, v_fac, m_fac, rho_fac, phi_fac, Hubble;
#endif
/*================================================
 *                  AHF_HALOS
 *================================================*/
void ahf_halos(gridls *grid_list)
{
  
  long unsigned   *idxtmp,*idx, ipart, ifrac, nhalos;
  double          *fpart;
  int 	          i,j,k,l,slen;
  FILE	         *fout;
  char	          filename[MAXSTRING],fprefix[MAXSTRING], file_no[MAXSTRING];
  partptr          current, cur_part;
  double           r,g,b;
  double           rad, t_relax, age;
  int              ibin, r_conv_i;
  double           a,a2,a3,rho_crit,omega,lambda,ovlim,rho_b,rho_vir;
  int	             numGoodHalos;
  double           Xhost, Yhost, Zhost, Xsat, Ysat, Zsat, dX, dY, dZ, Rhost2, Dsat2;
  long             nsat, *isat;
  SPATIALREF     **spatialRef;
  
  /* pointless trying to find halos when there are no useable grids */
  if(ahf.no_grids <= 0)
    return;
  
  /* conversion factors to get physical units */
  r_fac   = simu.boxsize*global.a;
  x_fac   = simu.boxsize;
  v_fac   = simu.boxsize/simu.t_unit/global.a;
  m_fac   = simu.pmass;
  rho_fac = simu.pmass/pow3(simu.boxsize);
  phi_fac = Grav*simu.pmass/(simu.boxsize*global.a);
  Hubble  = calc_Hubble(global.a);                   /* in km/sec/Mpc */
  
  /* cosmology related stuff */
  a        = global.a;
  a2       = pow2(global.a);
  a3       = pow3(global.a);
  rho_crit = rhoc0 * (simu.lambda0*(1.-1./a2) + simu.omega0*(1./a3-1./a2) + 1./a2);
  omega    =    simu.omega0  / (a+simu.omega0*(1.-a)+simu.lambda0*(a3-a));
  lambda   = a3*simu.lambda0 / (a+simu.omega0*(1.-a)+simu.lambda0*(a3-a));
  
  ovlim    = calc_virial(a);
    
#ifdef LCDM_DVIR
  {
    double omega0, lambda0;
    omega0       = simu.omega0;
    lambda0      = simu.lambda0;
    simu.omega0  = 0.3;
    simu.lambda0 = 0.7;
    ovlim        = calc_virial(a);
    simu.omega0  = omega0;
    simu.lambda0 = lambda0;
  }
#endif
  
  rho_b    =         a3 * omega * rho_crit;  /* comoving(!) background density */
  rho_vir  = ovlim * a3 * omega * rho_crit;  /* comoving(!)   virial   density */
  
  global.ovlim   = ovlim;
  global.rho_b   = rho_b;
  global.rho_vir = rho_vir;  
  
#ifdef VERBOSE
  /***************************************************************************/
  fprintf(stderr,"#################### ahf_halos ####################\n");
  fprintf(stderr,"AHF_MINPART  = %d\n",AHF_MINPART);
  fprintf(stderr,"AHF_VTUNE    = %g\n",AHF_VTUNE);
  fprintf(stderr,"Delta_vir    = %g\n",global.ovlim);
  fprintf(stderr,"Hubble       = %g\n",Hubble);
#endif
#ifdef AHFstep
  fprintf(io.logfile,"#################### ahf_halos ####################\n");
  fprintf(io.logfile,"AHF_MINPART  = %d\n",AHF_MINPART);
  fprintf(io.logfile,"AHF_VTUNE    = %g\n",AHF_VTUNE);
  fprintf(io.logfile,"Delta_vir    = %g\n",global.ovlim);
  fprintf(io.logfile,"Hubble       = %g\n",Hubble);
  fflush(io.logfile);
#endif
  
  
  
  
  /* prepare output filenames */
#	ifdef NEWSTARTRUN
#		ifdef WITH_MPI
  snprintf(fprefix, MAXSTRING, "%s.%04d.",
           global_io.params->outfile_prefix,
           global_mpi.rank);
#		else
#			ifdef AHFrestart
  snprintf(fprefix, MAXSTRING, "%s.%04d.",
           global_io.params->outfile_prefix,
           global_info.rank);
#			else
  snprintf(fprefix, MAXSTRING, "%s.",
           global_io.params->outfile_prefix);
#			endif
#		endif
#		ifdef DEBUG
  io_logging_msg(global_io.log, INT32_C(0),
                 "Using %s as the file prefix for AHF outputs.",
                 fprefix);
#		endif
#	else
  strcpy(fprefix, io.outfile_prefix);
#	endif
  
#ifdef REDSHIFTNAME
  /* OUTPUT CONVENTION: 3 digits*/
  sprintf(file_no,"z%.3f",fabs(global.z));
#else
  sprintf(file_no,"%05ld",global.no_timestep);
#endif
  strcat(fprefix, file_no);
  
  
  /***************************************************************************/
  /* Allocating the memory for the isolated refinements array                */
  /* This information is only a temporary holder and can be set free         */
  /* once the refinements have been joined to form halos                     */
  /***************************************************************************/
  spatialRef=NULL;
  if ((spatialRef = calloc(ahf.no_grids,sizeof(SPATIALREF*)))==NULL)
    {
      fprintf(stderr,"Error in allocating the memory for spatialRef array\n");
      exit(0);
    }
  for(i=0;i<ahf.no_grids;i++)
    {
      spatialRef[i] = NULL;
      if ((spatialRef[i] = calloc(numIsoRef[i],sizeof(SPATIALREF)))==NULL)
        {
          fprintf(stderr,"Error in allocating the memory for spatialRef array\n");
          exit(0);
        }	
    }
  
  /* Initialising the spatialRef array */
  for(i=0;i<ahf.no_grids;i++)
    {
      for(j=0;j<numIsoRef[i];j++)
        {
          
          spatialRef[i][j].numNodes      =  0;
          spatialRef[i][j].maxDens       = -1.0;
          
          spatialRef[i][j].centreDens.x      = 0.0;
          spatialRef[i][j].centreDens.y      = 0.0;
          spatialRef[i][j].centreDens.z      = 0.0;
          spatialRef[i][j].centreDens.norm   = 0.0;
          
          spatialRef[i][j].centrePot.x       = 0.0;
          spatialRef[i][j].centrePot.y       = 0.0;
          spatialRef[i][j].centrePot.z       = 0.0;
          spatialRef[i][j].centrePot.norm    = 0.0;
          
          spatialRef[i][j].centreGEOM.x      = 0.0;
          spatialRef[i][j].centreGEOM.y      = 0.0;
          spatialRef[i][j].centreGEOM.z      = 0.0;
          spatialRef[i][j].centreGEOM.norm   = 0.0;
          
          spatialRef[i][j].centreCMpart.x    = 0.0;
          spatialRef[i][j].centreCMpart.y    = 0.0;
          spatialRef[i][j].centreCMpart.z    = 0.0;
          spatialRef[i][j].centreCMpart.norm = 0.0;
          
          spatialRef[i][j].x.min = 100000.0;
          spatialRef[i][j].x.max = -100000.0;
          spatialRef[i][j].y.min = 100000.0;
          spatialRef[i][j].y.max = -100000.0;
          spatialRef[i][j].z.min = 100000.0;
          spatialRef[i][j].z.max = -100000.0;
          
          spatialRef[i][j].haloIndex            = -1;
          spatialRef[i][j].daughter.refLevel    = -1;
          spatialRef[i][j].daughter.isoRefIndex = -1;
          
          spatialRef[i][j].numParts     = 0;
          spatialRef[i][j].ll           = NULL;
          spatialRef[i][j].numSubStruct = 0;
          spatialRef[i][j].subStruct    = NULL;
          spatialRef[i][j].numParDom    = 0;
          spatialRef[i][j].parDom       = NULL;
          spatialRef[i][j].closeRefDist = -1.0;
          
          spatialRef[i][j].boundRefDiv.x = -1.0;
          spatialRef[i][j].boundRefDiv.y = -1.0;
          spatialRef[i][j].boundRefDiv.z = -1.0;
          
          
        }
    }
  
  /***************************************************************************/
  /* Collect the isolated refinements */
#ifdef VERBOSE
  fprintf(stderr,"\nCollecting the isolated refinements ... ");
  fprintf(io.logfile,"\nCollecting the isolated refinements ... ");
  fflush(io.logfile);
#endif
  if ( RefCentre(grid_list, ahf.no_grids, spatialRef) == 0 )
    {
      fprintf(stderr,"Stuffed up collecting the grid data\n");
      exit(-1);
    }
#ifdef VERBOSE
  fprintf(stderr,"done\n");
  fprintf(io.logfile,"done\n");
  fflush(io.logfile);
#endif
  
  /***************************************************************************/
  /* Analyse the isolated refinements */
#ifdef VERBOSE
  fprintf(stderr,"\nAnalysing the isolated refinements:\n");
  fprintf(io.logfile,"\nAnalysing the isolated refinements:\n");
  fflush(io.logfile);
#endif
  if ( analyseRef(ahf.no_grids, spatialRef) == 0 )
    {
      fprintf(stderr,"Stuffed up analysing the isolated refinements\n");
      exit(-1);
    }
  
  /***************************************************************************/
  /* Connect the isolated refinements */
#ifdef VERBOSE
  fprintf(stderr,"\nConverting the spatialRef[][] tree to a halos[] array:\n");
  fprintf(io.logfile,"\nConverting the spatialRef[][] tree to a halos[] array:\n");
  fflush(io.logfile);
#endif
  if ( spatialRef2halos(ahf.no_grids, spatialRef) == 0 )
    {
      fprintf(stderr,"Stuffed up converting spatialRef[][] to halos[] \n");
      exit(-1);
    }
  
  
  /***************************************************************************/
  /* Freeing the spatialRef array ALEX */
  for(i=0;i<ahf.no_grids;i++)
    free(spatialRef[i]);
  free(spatialRef);
  
  
#ifndef AHFphi_infty
#ifdef AHFstep
  {
    /* the grids are no longer needed => free memory */
    gridls *for_grid;
    int     no_grids;
    
    no_grids = ahf.no_grids    + ahf.min_ref - 1;
    for_grid = global.dom_grid + no_grids;
#ifdef VERBOSE
    fprintf(stderr,"\nFree'ing all grid structures:\n");
    fprintf(io.logfile,"\nFree'ing all grid structures:\n");
#endif
    while(no_grids > 0)
      {
#ifdef VERBOSE
        fprintf(stderr,"  free'ing grid %ld\n",for_grid->l1dim);
        fprintf(io.logfile,"  free'ing grid %ld\n",for_grid->l1dim);
        fflush(io.logfile);
#endif
        free_grid(for_grid, &no_grids);
        for_grid--;
      }
  }
#endif /* AHFstep */
#endif /* AHFphi_infty */
  
#	if (defined WITH_MPI || defined AHFrestart)
  rem_boundary_haloes();
#	endif
  
#ifdef AHFcentrefile
  {
    /***************************************************************************/
    /* Printing the _centres file */
    strcpy(filename,fprefix);
    strcat(filename,".AHF_centres");
    if((fout = fopen(filename,"w")) == NULL)
      {
        fprintf(stderr,"could not open %s\n", filename);
        exit(1);
      }
    for (j=0;j<numHalos;j++)
      {
        fprintf(fout,"P %g %g %g 0 1 0 6\n",
                halos[j].pos.x*simu.boxsize,
                halos[j].pos.y*simu.boxsize,
                halos[j].pos.z*simu.boxsize);
      }
    fclose(fout);
  }
#endif /* AHFcentrefile */
  
  
  /***************************************************************************/
  /* Ordering the halos particles and Removing unbound particles */
#ifdef VERBOSE
  fprintf(stderr,"\nConstructing Halos (%d):\n",numHalos);
  fprintf(stderr,"===================================\n");
  fprintf(io.logfile,"\nConstructing Halos (%d):\n",numHalos);
  fprintf(io.logfile,"===================================\n");
  fflush(io.logfile);
#endif 
#if (!defined NEWSTARTRUN)
  if ( ConstructHalos() == 0 )
    {
      fprintf(stderr,"Stuffed up ordering the halos particles\n");
      exit(-1);
    }
#else
#	pragma omp parallel for \
schedule (dynamic) \
shared (halos, numHalos)
  for (i=0; i<numHalos; i++) {
    ahf_halos_sfc_constructHalo(halos+i);
  }
#endif
  
  /***************************************************************************/
  /* Order halos with respects to their mass */
  idx      = (long unsigned *)  calloc(numHalos,   sizeof(long unsigned));
  idxtmp   = (long unsigned *)  calloc(numHalos+1, sizeof(long unsigned));
  fpart    = (double *)         calloc(numHalos+1, sizeof(double));
  for(i=0; i<numHalos; i++)
    fpart[i+1] = (double)halos[i].npart;
  indexx(numHalos,fpart,idxtmp);
  
  /* indexx sorts ascending, but we want descending */
  for(i=0; i<numHalos; i++)
    idx[numHalos-i-1] = idxtmp[i+1]-1;
  
  free(idxtmp);
  free(fpart);
  
  /***************************************************************************/
  /***************************************************************************/
  /***************************************************************************/
  /* PRINTING RELEVANT FILES */
#ifdef VERBOSE
  fprintf(stderr,"### PRINTING HALO FILES (%d)\n",numHalos);
#endif  
  
#if (  (defined AHFsubstructure || defined AHFgridsubstructure) && !defined WITH_MPI && !defined AHFrestart && !defined AHFbinary   )
  /***************************************************************************/
  /* Printing the _substructure file */
  strcpy(filename,fprefix);
  strcat(filename,".AHF_substructure");
#ifdef VERBOSE
  fprintf(stderr,"%s\n",filename);
#endif
  if((fout = fopen(filename,"w")) == NULL) 
    {
      fprintf(stderr,"could not open %s\n", filename);
      exit(1);
    }
  
  fprintf(fout,"%d\n",numHalos);
  
  for (j=0;j<numHalos;j++)
    {
      
      i = idx[j];
      
#ifdef AHFgridsubstructure
      if ( halos[i].numSubStruct > 0 && halos[i].npart >= AHF_MINPART && (double)halos[i].npart*simu.min_weight/(double)halos[i].M_vir > AHF_MASSMIX )
        {
          /* does *NOT* check if substructure lies within virial radius of host */
          fprintf(fout,"%10d %12d\n",j,halos[i].numSubStruct);
          
          for (k=0;k<halos[i].numSubStruct;k++)
            {
              l = idx_inv(idx, numHalos, halos[i].subStruct[k]);
              fprintf(fout,"%10d ", l);
            }
          
          
          fprintf(fout,"\n");
        }
#else /* AHFgridsubstructure */
      if ( halos[i].npart >= AHF_MINPART && (double)halos[i].npart*simu.min_weight/(double)halos[i].M_vir > AHF_MASSMIX ) 
        {
          /* check if substructure lies within virial radius of host */
          Xhost  = halos[i].pos.x;
          Yhost  = halos[i].pos.y;
          Zhost  = halos[i].pos.z;
          Rhost2 = pow2(halos[i].R_vir);
          
          nsat = 0;
          isat = NULL;
          
          for (k=0;k<numHalos;k++)
            {
              /* make sure to use un-ordered id */
              l = idx[k];
              
              if(l != i)
                {
                  
                  Xsat = halos[l].pos.x;
                  Ysat = halos[l].pos.y;
                  Zsat = halos[l].pos.z;
                  
                  /* distance to host centre */
                  dX = Xsat-Xhost;
                  dY = Ysat-Yhost;
                  dZ = Zsat-Zhost;
                  
                  /* take care of periodic boundaries */
                  if(dX >  0.5) dX -= 1.0;
                  if(dY >  0.5) dY -= 1.0;
                  if(dZ >  0.5) dZ -= 1.0;
                  if(dX < -0.5) dX += 1.0;
                  if(dY < -0.5) dY += 1.0;
                  if(dZ < -0.5) dZ += 1.0;
                  
                  Dsat2 = pow2(dX) + pow2(dY) + pow2(dZ);
                  
                  /* within host's virial radius */
                  if(Dsat2 < Rhost2 && halos[l].npart >= AHF_MINPART && (double)halos[l].npart*simu.min_weight/(double)halos[l].M_vir > AHF_MASSMIX)
                    {
                      nsat++;
                      isat         = (long *) realloc(isat, nsat*sizeof(long));
                      isat[nsat-1] = idx_inv(idx, numHalos, l);
                    }
                }
            }
          
          /* write total number of satellites found */
          fprintf(fout,"%10d %12ld\n",j,nsat);
          
          /* write satellite id's */
          for (k=0;k<nsat;k++)
            {
              l = isat[k];
              fprintf(fout,"%10d ", l);
            }
          free(isat);
          
          fprintf(fout,"\n");
        }
#endif /* AHFgridsubstructure */   
    }
  
  fclose(fout);
#endif /* AHFsubstructure || AHFgridsubstructure */

#ifdef AHFbinary  
  ahf_write_profiles(fprefix, halos, idx, numHalos);
  ahf_write_particles(fprefix, halos, idx, numHalos);
  ahf_write_halos(fprefix, halos, idx, numHalos);
#else
  /***************************************************************************/
  /* Printing the _profiles file */
  strcpy(filename,fprefix);
  strcat(filename,".AHF_profiles");
#ifdef VERBOSE
  fprintf(stderr,"%s\n",filename);
#endif
  if((fout = fopen(filename,"w")) == NULL)
    {
      fprintf(stderr,"could not open %s\n", filename);
      exit(1);
    }
#if (WITH_MPI || AHFrestart)
#	ifdef WITH_MPI
  if(global_mpi.rank == 0)
#	else
    if(global_info.rank == 0)
#	endif
      {
#endif
#ifdef GAS_PARTICLES
        fprintf(fout,"#  r(1)     npart(2)    nvpart(3)    ovdens(4)     dens(5)   vcirc(6)  sigv(7)     Lx(8)        Ly(9)        Lz(10)      a(11)      Eax(12)    Eay(13)     Eaz(14)    b(15)      Ebx(16)     Eby(17)    Ebz(18)    c(19)      Ecx(20)    Ecy(21)    Ecz(22)      Ekin(23)       Epot(24)  nvpart_gas(25)  nvpart_stars(26)");
#else
        fprintf(fout,"#  r(1)     npart(2)    nvpart(3)    ovdens(4)     dens(5)   vcirc(6)  sigv(7)     Lx(8)        Ly(9)        Lz(10)      a(11)      Eax(12)    Eay(13)     Eaz(14)    b(15)      Ebx(16)     Eby(17)    Ebz(18)    c(19)      Ecx(20)    Ecy(21)    Ecz(22)      Ekin(23)       Epot(24)");
#endif
#ifdef AHFphspdens
        fprintf(fout, "\tsigma2_vx_sh(25)\tsigma2_vy_sh(26)");
        fprintf(fout, "\tsigma2_vz_sh(27)\tsigma2_vr_sh(28)");
        fprintf(fout, "\tsigma2_vtheta_sh(29)\tsigma2_vphi_sh(30)");
#ifdef AHFmeanvelocities
        fprintf(fout, "\tmean_vx_sh(31)\tmean_vy_sh(32)");
        fprintf(fout, "\tmean_vz_sh(33)\tmean_vr_sh(34)");
        fprintf(fout, "\tmean_vtheta_sh(35)\tmean_vphi_sh(36)");
        fprintf(fout, "\tmean_vx_sp(37)\tmean_vy_sp(38)");
        fprintf(fout, "\tmean_vz_sp(39)\tmean_vr_sp(40)");
        fprintf(fout, "\tmean_vtheta_sp(41)\tmean_vphi_sp(42)");
#endif
#endif
        fprintf(fout, "\n");
#if (WITH_MPI || AHFrestart)
      }
#endif
  
  for (j=0;j<numHalos;j++)
    {
      i = idx[j];
      
      if(halos[i].npart >= AHF_MINPART && (double)halos[i].npart*simu.min_weight/(double)halos[i].M_vir > AHF_MASSMIX)
        {
          for(r_conv_i=0, ibin=0; ibin<halos[i].prof.nbins; ibin++)
            {
              /* check for converged radius (Power et al. 2003) */
              rad      = halos[i].prof.r[ibin];
              t_relax  = halos[i].prof.npart[ibin]/log(rad/halos[i].spaRes)
              * rad / sqrt(halos[i].prof.v2_circ[ibin]);
              
              /* convert to (Mpc/h) / (km/sec) ~ (Mpc/h) / (kpc/Gyr) */
              t_relax *= r_fac/sqrt(phi_fac);
              
              /* convert to Gyr/h */
              t_relax *= 1E3;
              
              /* age of the Universe in Gyr/h */
              age      = calc_t(global.a) * simu.t_unit*Mpc/Gyr;
              
              /* if not converged, write negative radius into .AHF_profiles */
              if(t_relax < 0.9*age)
                /* The +1 is merely to keep the name consistent:
                 * r_conv_i should give the smallest converged radius, not
                 * the bin before it. Hence to check for converged bin or
                 * not, i<r_conv_i is the thing to do. */
                r_conv_i = ibin+1;
            }
          
          for(ibin=0; ibin<halos[i].prof.nbins; ibin++)
            {
              /* If the radius is not converged, print the radius negative
               * to the .AHF_profiles */
              rad = (ibin<r_conv_i) ?  -halos[i].prof.r[ibin] : halos[i].prof.r[ibin];
              
              fprintf(fout,"%8.4f %10ld %e %12.2f %12.2f %8.2f %8.2f %12.4g %12.4g %12.4g  %10.6f %10.6f %10.6f  %10.6f %10.6f %10.6f  %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f    %e   %e",
                      rad                              * x_fac*1000.,
                      halos[i].prof.npart[ibin],
                      halos[i].prof.nvpart[ibin],
                      halos[i].prof.ovdens[ibin]       * rho_fac/global.rho_b,
                      halos[i].prof.dens[ibin]         * rho_fac/global.rho_b,
                      sqrt(halos[i].prof.v2_circ[ibin] * phi_fac),
                      halos[i].prof.sig_v[ibin]        * v_fac,
                      halos[i].prof.Lx[ibin]           * m_fac*r_fac*v_fac,
                      halos[i].prof.Ly[ibin]           * m_fac*r_fac*v_fac,
                      halos[i].prof.Lz[ibin]           * m_fac*r_fac*v_fac,
                      halos[i].prof.axis1[ibin],
                      halos[i].prof.E1x[ibin],
                      halos[i].prof.E1y[ibin],
                      halos[i].prof.E1z[ibin],
                      halos[i].prof.axis2[ibin],
                      halos[i].prof.E2x[ibin],
                      halos[i].prof.E2y[ibin],
                      halos[i].prof.E2z[ibin],
                      halos[i].prof.axis3[ibin],
                      halos[i].prof.E3x[ibin],
                      halos[i].prof.E3y[ibin],
                      halos[i].prof.E3z[ibin],
                      halos[i].prof.Ekin[ibin]        * m_fac*pow2(v_fac),
                      halos[i].prof.Epot[ibin]        * m_fac*phi_fac
                      );
              
#ifdef AHFphspdens
              fprintf(fout,
                      "\t%g\t%g\t%g\t%g\t%g\t%g",
                      halos[i].prof.sigma2_vx_sh[ibin]*v_fac*v_fac,
                      halos[i].prof.sigma2_vy_sh[ibin]*v_fac*v_fac,
                      halos[i].prof.sigma2_vz_sh[ibin]*v_fac*v_fac,
                      halos[i].prof.sigma2_vr_sh[ibin]*v_fac*v_fac,
                      halos[i].prof.sigma2_vtheta_sh[ibin],
                      halos[i].prof.sigma2_vphi_sh[ibin]);
#ifdef AHFmeanvelocities
              fprintf(fout,
                      "\t%g\t%g\t%g\t%g\t%g\t%g"
                      "\t%g\t%g\t%g\t%g\t%g\t%g",
                      halos[i].prof.mean_vx_sh[ibin]*v_fac,
                      halos[i].prof.mean_vy_sh[ibin]*v_fac,
                      halos[i].prof.mean_vz_sh[ibin]*v_fac,
                      halos[i].prof.mean_vr_sh[ibin]*v_fac,
                      halos[i].prof.mean_vtheta_sh[ibin],
                      halos[i].prof.mean_vphi_sh[ibin],
                      halos[i].prof.mean_vx_sp[ibin]*v_fac,
                      halos[i].prof.mean_vy_sp[ibin]*v_fac,
                      halos[i].prof.mean_vz_sp[ibin]*v_fac,
                      halos[i].prof.mean_vr_sp[ibin]*v_fac,
                      halos[i].prof.mean_vtheta_sp[ibin],
                      halos[i].prof.mean_vphi_sp[ibin]);
#endif
#endif
#ifdef GAS_PARTICLES
              fprintf(fout,"     %e     %e",
                      halos[i].prof.nvpart_gas[ibin],
                      halos[i].prof.nvpart_star[ibin]);
#endif
              fprintf(fout,"\n");
            }
        }
    }
  
  fclose(fout);
  
#ifndef MARE_NOSTRUM
  /***************************************************************************/
  /* Printing the _mass file */
  strcpy(filename,fprefix);
  strcat(filename,".AHF_mass");
#ifdef VERBOSE
  fprintf(stderr,"%s\n",filename);
#endif
  if((fout = fopen(filename,"w")) == NULL) {
    fprintf(stderr,"could not open %s\n", filename);
    exit(1);
  }
  
  nhalos = 1;
  for(j=0; j<numHalos; j++)
    {
      i = idx[j];
      if(halos[i].npart >= AHF_MINPART && (double)halos[i].npart*simu.min_weight/(double)halos[i].M_vir > AHF_MASSMIX)
        {
          fprintf(fout,"%g %ld\n",halos[i].M_vir*m_fac,nhalos);
          nhalos++;
        }
    }
  
  fclose(fout);
#endif /* MARE_NOSTRUM */
  
  /***************************************************************************/
  /* Printing the _halos file */
  strcpy(filename,fprefix);
  strcat(filename,".AHF_halos");
#ifdef VERBOSE
  fprintf(stderr,"%s\n",filename);
#endif
  if((fout = fopen(filename,"w")) == NULL) 
    {
      fprintf(stderr,"could not open %s\n", filename);
      exit(1);
    }
  
#if (WITH_MPI || AHFrestart)
#	ifdef WITH_MPI
  if(global_mpi.rank == 0)
#	else
  if(global_info.rank == 0)
#	endif
    {
#endif
#ifdef GAS_PARTICLES
      fprintf(fout,"#  npart(1)  nvpart(2)      Xc(3)          Yc(4)          Zc(5)      VXc(6)   VYc(7)   VZc(8)    Mvir(9)      Rvir(10) Vmax(11)   Rmax(12) sigV(13) lambda(14)      Lx(15)       Ly(16)       Lz(17)     a(18)     Eax(19)    Eay(20)    Eaz(21)     b(22)     Ebx(23)    Eby(24)    Ebz(25)     c(26)     Ecx(27)    Ecy(28)    Ecz(29) ovdens(30)    Redge(31) nbins(32)     Ekin(33)       Epot(34)   mbp_offset(35)  com_offset(36)   r2(37)      lambdaE(38)     n_gas(39)   M_gas(40)   lambda_gas(41)  Lx_gas(42)  Ly_gas(43)  Lz_gas(44)   a_gas(45)  Eax_gas(46) Eay_gas(47) Eaz_gas(48)   b_gas(49)  Ebx_gas(50) Eby_gas(51) Ebz_gas(52)   c_gas(53)  Ecx_gas(54) Ecy_gas(55) Ecz_gas(56)  Ekin_gas(57)  Epot_gas(58)  lambdaE_gas(59)   phi0(60)      n_star(61)   M_star(62)  lambda_star(63) Lx_star(64) Ly_star(65) Lz_star(66)  a_star(67) Eax_star(68) Eay_star(69) Eaz_star(70)   b_star(71) Ebx_star(72) Eby_star(73) Ebz_star(74)   c_star(75) Ecx_star(76) Ecy_star(77) Ecz_star(78) Ekin_star(79) Epot_star(80) lambdaE_star(81)");
      fprintf(fout,"\n");
      
#else
      fprintf(fout,"  npart_1  nvpart_2      Xc_3          Yc_4          Zc_5      VXc_6   VYc_7   VZc_8    Mvir_9      Rvir_10 Vmax_11   Rmax_12 sigV_13 lambda_14      Lx_15       Ly_16       Lz_17     a_18     Eax_19    Eay_20    Eaz_21     b_22     Ebx_23    Eby_24    Ebz_25     c_26     Ecx_27    Ecy_28    Ecz_29 ovdens_30    Redge_31 nbins_32     Ekin_33       Epot_34   mbp_offset_35  com_offset_36   r2_37      lambdaE_38");
      fprintf(fout,"\n");
#endif  /* GAS_PARTICLES */
#if (WITH_MPI || AHFrestart)
    }
#endif
  
  numGoodHalos = 0;
  for (j=0; j<numHalos; j++)
    {
      i = idx[j];
      
#ifdef AHFmaxhalo
      if(halos[i].npart > AHF_MAXHALO)
        {
#ifdef VERBOSELOG
          fprintf(io.logfile,"ahf_halos: found halo with more particles than %d -> terminating AMIGA!\n",AHF_MAXHALO);
#endif
          global.terminate = TRUE;
        }
#endif
      
      if(halos[i].npart >= AHF_MINPART && (double)halos[i].npart*simu.min_weight/(double)halos[i].M_vir > AHF_MASSMIX)
        {
          numGoodHalos++;
          
#ifdef GAS_PARTICLES
          
          fprintf(fout,"%10ld %e %14.8f %14.8f %14.8f %8.2f %8.2f %8.2f %12.6g %10.2f %8.2f %10.2f %8.2f %10.6f %12.4g %12.4g %12.4g %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %8.2f %12.4g %6d       %12.6g   %12.6g  %10.5f      %10.5f     %10.5f  %12.6f   %10ld   %12.6g   %10.6f   %10.6f  %10.6f  %10.6f   %10.6f   %10.6f  %10.6f  %10.6f  %10.6f   %10.6f  %10.6f  %10.6f  %10.6f   %10.6f  %10.6f  %10.6f  %12.6g   %12.6g   %10.6f    %12.6g   %10ld   %12.6g   %10.6f   %10.6f  %10.6f  %10.6f   %10.6f   %10.6f  %10.6f  %10.6f  %10.6f   %10.6f  %10.6f  %10.6f  %10.6f   %10.6f  %10.6f  %10.6f  %12.6g   %12.6g   %10.6f",
                  halos[i].npart,
                  halos[i].M_vir,
                  halos[i].pos.x*x_fac, halos[i].pos.y*x_fac, halos[i].pos.z*x_fac,
                  halos[i].vel.x*v_fac, halos[i].vel.y*v_fac, halos[i].vel.z*v_fac,
                  halos[i].M_vir*m_fac,
                  halos[i].R_vir*x_fac*1000.,
                  sqrt(halos[i].V2_max*phi_fac),
#ifdef AHFDEBUG
                  halos[i].gatherRad*x_fac*1000.,   /* used for debugging purposes... */
#else
                  halos[i].R_max*x_fac*1000.,
#endif
                  halos[i].velDis*v_fac,
                  halos[i].lambda,
#ifdef AHFabsangmom
                  halos[i].AngMom.x*r_fac*v_fac*m_fac,
                  halos[i].AngMom.y*r_fac*v_fac*m_fac,
                  halos[i].AngMom.z*r_fac*v_fac*m_fac,
#else
                  halos[i].AngMom.x,
                  halos[i].AngMom.y,
                  halos[i].AngMom.z,
#endif
                  halos[i].axis.x,
                  halos[i].E1.x,
                  halos[i].E1.y,
                  halos[i].E1.z,
                  halos[i].axis.y,
                  halos[i].E2.x,
                  halos[i].E2.y,
                  halos[i].E2.z,
                  halos[i].axis.z,
                  halos[i].E3.x,
                  halos[i].E3.y,
                  halos[i].E3.z,
                  halos[i].ovdens,
                  halos[i].R_edge*x_fac*1000.,
                  halos[i].prof.nbins,
                  halos[i].Ekin*m_fac*pow2(v_fac),
                  halos[i].Epot*m_fac*phi_fac,
                  halos[i].mbp_offset*x_fac*1000.,
                  halos[i].com_offset*x_fac*1000.,
                  halos[i].r2*x_fac*1000.,
                  halos[i].lambdaE,
                  halos[i].gas.npart,
                  halos[i].gas.M_vir*m_fac,
                  halos[i].gas.lambda,
                  halos[i].gas.AngMom.x,
                  halos[i].gas.AngMom.y,
                  halos[i].gas.AngMom.z,
                  halos[i].gas.axis.x,
                  halos[i].gas.E1.x,
                  halos[i].gas.E1.y,
                  halos[i].gas.E1.z,                 
                  halos[i].gas.axis.y,
                  halos[i].gas.E2.x,
                  halos[i].gas.E2.y,
                  halos[i].gas.E2.z,                 
                  halos[i].gas.axis.z,
                  halos[i].gas.E3.x,
                  halos[i].gas.E3.y,
                  halos[i].gas.E3.z,
                  halos[i].gas.Ekin*m_fac*pow2(v_fac),
                  halos[i].gas.Epot*m_fac*phi_fac,
                  halos[i].gas.lambdaE,
                  halos[i].Phi0*m_fac*phi_fac,
                  halos[i].star.npart,
                  halos[i].star.M_vir*m_fac,
                  halos[i].star.lambda,
                  halos[i].star.AngMom.x,
                  halos[i].star.AngMom.y,
                  halos[i].star.AngMom.z,
                  halos[i].star.axis.x,
                  halos[i].star.E1.x,
                  halos[i].star.E1.y,
                  halos[i].star.E1.z,                 
                  halos[i].star.axis.y,
                  halos[i].star.E2.x,
                  halos[i].star.E2.y,
                  halos[i].star.E2.z,                 
                  halos[i].star.axis.z,
                  halos[i].star.E3.x,
                  halos[i].star.E3.y,
                  halos[i].star.E3.z,
                  halos[i].star.Ekin*m_fac*pow2(v_fac),
                  halos[i].star.Epot*m_fac*phi_fac,
                  halos[i].star.lambdaE                  
                  );
          fprintf(fout, "\n");
#else /* GAS_PARTICLES */
          
          fprintf(fout,"%10ld %10.0f %14.8f %14.8f %14.8f %8.2f %8.2f %8.2f %12.6g %10.2f %8.2f %10.2f %8.2f %10.6f %12.4g %12.4g %12.4g %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %8.2f %12.4g %6d       %12.6g   %12.6g  %10.5f      %10.5f     %10.5f  %12.6f",
                  halos[i].npart,
                  halos[i].M_vir,
                  halos[i].pos.x*x_fac, halos[i].pos.y*x_fac, halos[i].pos.z*x_fac,
                  halos[i].vel.x*v_fac, halos[i].vel.y*v_fac, halos[i].vel.z*v_fac,
                  halos[i].M_vir*m_fac,
                  halos[i].R_vir*x_fac*1000.,
                  sqrt(halos[i].V2_max*phi_fac),
#ifdef AHFDEBUG
                  halos[i].gatherRad*x_fac*1000.,   /* used for debugging purposes... */
#else
                  halos[i].R_max*x_fac*1000.,
#endif
                  halos[i].velDis*v_fac,
                  halos[i].lambda,
#ifdef AHFabsangmom
                  halos[i].AngMom.x*r_fac*v_fac*m_fac,
                  halos[i].AngMom.y*r_fac*v_fac*m_fac,
                  halos[i].AngMom.z*r_fac*v_fac*m_fac,
#else
                  halos[i].AngMom.x,
                  halos[i].AngMom.y,
                  halos[i].AngMom.z,
#endif
                  halos[i].axis.x,
                  halos[i].E1.x,
                  halos[i].E1.y,
                  halos[i].E1.z,
                  halos[i].axis.y,
                  halos[i].E2.x,
                  halos[i].E2.y,
                  halos[i].E2.z,
                  halos[i].axis.z,
                  halos[i].E3.x,
                  halos[i].E3.y,
                  halos[i].E3.z,
                  halos[i].ovdens,
                  halos[i].R_edge*x_fac*1000.,
                  halos[i].prof.nbins,
                  halos[i].Ekin*m_fac*pow2(v_fac),
                  halos[i].Epot*m_fac*phi_fac,
                  halos[i].mbp_offset*x_fac*1000.,
                  halos[i].com_offset*x_fac*1000.,
                  halos[i].r2*x_fac*1000.,
                  halos[i].lambdaE
                  );
          fprintf(fout, "\n");
#endif /* GAS_PARTICLES */
          
        }
    }
  
  fclose(fout);
  
#ifndef MARE_NOSTRUM /* do not write a .AHF_particles file for MareNostrum subboxes */
  /***************************************************************************/
  /* Printing the _particles file */
  strcpy(filename,fprefix);
  strcat(filename,".AHF_particles");
#ifdef VERBOSE
  fprintf(stderr,"%s\n",filename);
#endif
  if((fout = fopen(filename,"w")) == NULL)
    {
      fprintf(stderr,"could not open %s\n", filename);
      exit(1);
    }
#if (defined PARTICLES_INFO)
  fprintf(fout,"%d\n",numGoodHalos);
#endif
  for (j=0;j<numHalos;j++)
    {      
      i = idx[j];
      
      if(halos[i].npart >= AHF_MINPART && (double)halos[i].npart*simu.min_weight/(double)halos[i].M_vir > AHF_MASSMIX)
        {
          /* start with dumping the actual number of particle IDs to follow... */
          fprintf(fout,"%lu\n", (unsigned long)(halos[i].npart));
          
#ifdef PARTICLES_INFO
          
#if (defined NEWSTARTRUN)
          
          for(ipart=0; ipart<halos[i].npart; ipart++)
            {
              current = global.fst_part + halos[i].ipart[ipart];
#if (defined GAS_PARTICLES)
              fprintf(fout, "%"PRIpartid"\t%d\n",
                      current->id,
                      (current->u >= PGAS) ? (int)PGAS : (int)(-current->u));
#else
              fprintf(fout, "%"PRIpartid"\n", current->id);
#endif
            }
#else /* NEWSTARTRUN */
          
          for(ipart=0; ipart<halos[i].npart; ipart++)
            {
              cur_part = global.fst_part + halos[i].ipart[ipart];
              
              /* gas particle */
              if(cur_part-global.fst_part-global.offset_gas < global.no_gas && cur_part-global.fst_part-global.offset_gas > 0)
                fprintf(fout,"%d  %d  %g\n",cur_part->ID, (int)PGAS, cur_part->weight*m_fac);
              
              /* DM particle */
              else if(cur_part->ID >= gadget.IDmin && cur_part->ID <= gadget.IDmax)
                fprintf(fout,"%d  %d  %g\n",cur_part->ID, (int)(-PDM), cur_part->weight*m_fac);
              
              /* heavy DM particles */
              else if(cur_part->weight > (global.fst_part+global.no_gas+1)->weight)
                fprintf(fout,"%d  %d  %g\n",cur_part->ID, (int)(-PDMbndry), cur_part->weight*m_fac);
              
              /* star particle */
              else
                fprintf(fout,"%d  %d  %g\n",cur_part->ID, (int)(-PSTAR), cur_part->weight*m_fac);
              
              cur_part = cur_part->ll;
            }
          
#endif /* WITH_MPI */
          
#else /* PARTICLES_INFO */
          
#if (defined NEWSTARTRUN)
          
          for(ipart=0; ipart<halos[i].npart; ipart++)
            {
              current = global.fst_part + halos[i].ipart[ipart];
              fprintf(fout,"%"PRIpartid"\n",current->id);
            }
          
#else /* NEWSTARTRUN */
          
#ifdef GADGET_IDS
          /* this will simply dump the GADGET IDs to file ... irrespective of their meaning */
          for(ipart=0; ipart<halos[i].npart; ipart++)
            {
              current = global.fst_part + halos[i].ipart[ipart];
              fprintf(fout,"%ld\n",current->ID);
            }
          
#else /* GADGET_IDS */
          /* here we dump the AMIGA-assigned IDs to be used with amigaExtractHalos.c */
          for(ipart=0; ipart<halos[i].npart; ipart++)
            {
              fprintf(fout,"%ld\n",halos[i].ipart[ipart]);
            }
          
#endif /* GADGET_IDS */
#endif /* WITH_MPI */
#endif /* PARTICLES_INFO */
        }
    }
  
  fclose(fout);
#endif /* MARE_NOSTRUM */
  
#ifdef AHFgeom
  /***************************************************************************/
  /*   FOR TESTING IN STEREO2:
   *
   * -DAHFgeom write one addition file in GEOM format
   *
   *  this file includes all halos as spheres ... as well as ...
   *  -DAHFgeom_SIMUPART: a selection of all simulation particles as points
   *  -DAHFgeom_HALOPART: all particles in every halo as points
   */
  
  /***************************************************************************/
  /* Printing the _halos.geom file */
  strcpy(filename,fprefix);
  strcat(filename,".AHF_halos.geom");
#ifdef VERBOSE
  fprintf(stderr,"%s\n",filename);
#endif
  if((fout = fopen(filename,"w")) == NULL) {
    fprintf(stderr,"could not open %s\n", filename);
    exit(1);
  }
#ifdef AHFgeom_SIMUPART   
  /* 1. write about 500000 simulation particles to AHF_geom file */
#ifdef WITH_MPI
  if (global_info.no_part > 500000)
    ifrac = (long unsigned)((double)global_info.no_part(double)/(500000.));
#else
  if (global.no_part > 500000)
    ifrac = (long unsigned)((double)global.no_part/(double)(500000.));
#endif
  else
    ifrac = (long unsigned)1;
  
#if (defined WITH_MPI || defined AHFrestart)
  for (ipart=0; ipart<global_info.no_part; ipart+=ifrac)
#else
    for (ipart=0; ipart<global.no_part; ipart+=ifrac)
#endif
      {
        cur_part = global.fst_part + ipart;
#ifdef MULTIMASS
        if(cur_part->weight < simu.max_weight)
          fprintf(fout,"p %12.6g %12.6g %12.6g 0 0 1\n",cur_part->pos[X]*x_fac,cur_part->pos[Y]*x_fac,cur_part->pos[Z]*x_fac);
        else
          fprintf(fout,"p %12.6g %12.6g %12.6g 0.3 0 0\n",cur_part->pos[X]*x_fac,cur_part->pos[Y]*x_fac,cur_part->pos[Z]*x_fac);
#else
        fprintf(fout,"p %12.6g %12.6g %12.6g 0 0 1\n",cur_part->pos[X]*x_fac,cur_part->pos[Y]*x_fac,cur_part->pos[Z]*x_fac);
#endif
      }
#endif /* AHFgeom_SIMUPART */
  
  /* 2. write halos */
  for (j=0; j<numHalos; j++)
    {
      i = idx[j];
      
      r = 1.0 - (double)i/(double)numHalos;
      g = (double)i/(double)numHalos;
      b = 0.0;
      
      r = 1.0;
      g = 0.0;
      b = 0.0;
      
      if(halos[i].npart >= AHF_MINPART && (double)halos[i].npart*simu.min_weight/(double)halos[i].M_vir > AHF_MASSMIX)
        {
#ifdef AHFgeom_HALOPART
          current = halos[i].ll;
          while ( current!=NULL ) 
            {
              fprintf(fout,"p %12.6g %12.6g %12.6g  %g %g %g\n",
                      current->pos[X]*x_fac, current->pos[Y]*x_fac, current->pos[Z]*x_fac, r,g,b);            
              current 	= current->ll;
            }
#else
          fprintf(fout,"s %12.6g %12.6g %12.6g %12.6g  %g %g %g\n",
                  halos[i].pos.x*x_fac, halos[i].pos.y*x_fac, halos[i].pos.z*x_fac,halos[i].R_vir*x_fac,r,g,b);
#endif
        }
    }
  
  fclose(fout);
  
#endif /* AHFgeom */ 
  
  
  /* remove halos[], if we are just interested in a snapshot analysis */
  for (i=0;i<numHalos;i++)
    free(halos[i].ipart);
  free(halos);
  
  /* free all sorts of things */
  free(idx);
  free(gridl1dim);
  free(numIsoRef);
  if(numDensZero > 0)
    free(densZero);

#endif /* AHFbinary */
  
#ifdef VERBOSE
  fprintf(io.logfile,"################## ahf_halos finished ##################\n");
  fflush(io.logfile);
#endif
}

/*
 ===============================================================================
 END OF MAIN
 ===============================================================================
 */


int	RefCentre(gridls *grid_list, int num_refgrids, SPATIALREF **spatialRef) 
{
  gridls *cur_grid;
  pqptr   cur_pquad;
  cqptr   cur_cquad, icur_cquad;
  nqptr   cur_nquad, icur_nquad;
  nptr    cur_node;
  nptr    tsc_nodes[3][3][3];
  
  int     i,j,k;
  int	    numNodes, tmpNum;	
  int     refLevel, isoRefIndex;
  
  MINMAX  tmpMinMax;
  long	 x,y,z;
  double	 xx,yy,zz;
  double	 tmpDens, tmpPot;
  
  partptr current, previous, tmpll;
  
  int	    boundRefIndex;
  
  double  boxLen, boxVol;
  double	 px,py,pz,tmpRad,a,b,c,alpha;
  int     colour;
  double	 fl1dim;
  
  double	 xmin,xmax,ymin,ymax,zmin,zmax;
  
  int     iterate;
  
  double  cur_shift;
  
  
  /***************************************************************************/
  /* Calculating the spatial resolution (boxsize) on each refinement level
   */
  
  gridl1dim=NULL;
  if ((gridl1dim = calloc(num_refgrids,sizeof(double)))==NULL) 
    {
      fprintf(stderr,"Error in allocating the memory for halo array\n");
      exit(0);
    }
  cur_grid=global.dom_grid+ahf.min_ref;
  for ( i=0; i<num_refgrids; i++) 
    {
      fl1dim = ((double)(cur_grid->l1dim));
      gridl1dim[i]=fl1dim;
      cur_grid++;
    }
  
  
  /***************************************************************************/
  /* Collecting information about isolated refinements
   * 
   * number of particles
   * ll of particles
   * number of nodes
   * density centre of the refinement
   */
  
  for( 	iterate=0, cur_grid=global.dom_grid+ahf.min_ref; iterate<ahf.no_grids; iterate++, cur_grid++) 
    {
      /* shift of cell centre as compared to edge of box [grid units] */
      cur_shift = 0.5/(double)cur_grid->l1dim;
      
      for(cur_pquad=cur_grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next)
        {
          z = cur_pquad->z;
          for(cur_cquad = cur_pquad->loc;
              cur_cquad < cur_pquad->loc + cur_pquad->length; 
              cur_cquad++, z++)  
            {  
              for(icur_cquad  = cur_cquad; 
                  icur_cquad != NULL; 
                  icur_cquad  = icur_cquad->next)
                {
                  y = icur_cquad->y;
                  for(cur_nquad = icur_cquad->loc;  
                      cur_nquad < icur_cquad->loc + icur_cquad->length; 
                      cur_nquad++, y++) 
                    { 
                      for(icur_nquad  = cur_nquad; 
                          icur_nquad != NULL; 
                          icur_nquad  = icur_nquad->next)
                        {
                          x = icur_nquad->x;
                          for(cur_node = icur_nquad->loc; 
                              cur_node < icur_nquad->loc + icur_nquad->length; 
                              cur_node++, x++)
                            {
                              /* Locate the spatial refinement */
                              refLevel    = spatialRefIndex[cur_node->force.colour].refLevel;
                              isoRefIndex = spatialRefIndex[cur_node->force.colour].isoRefIndex;
                              
                              /* Colour of the isolated refinement */
                              spatialRef[refLevel][isoRefIndex].colour = cur_node->force.colour;
                              
                              /* Count the nodes for the refinement */	
                              spatialRef[refLevel][isoRefIndex].numNodes++;
                              
                              /* Position of the node */
                              xx = (double) x/(double)cur_grid->l1dim + cur_shift;
                              xx = (double)fmod(xx + 1.0, 1.0);
                              yy = (double) y/(double)cur_grid->l1dim + cur_shift;
                              yy = (double)fmod(yy + 1.0, 1.0);
                              zz = (double) z/(double)cur_grid->l1dim + cur_shift;
                              zz = (double)fmod(zz + 1.0, 1.0);
                              
                              /* account for periodic boundary conditions */
                              px = spatialRefIndex[cur_node->force.colour].periodic.x;
                              py = spatialRefIndex[cur_node->force.colour].periodic.y;
                              pz = spatialRefIndex[cur_node->force.colour].periodic.z;
                              
                              if(px == 1 && xx < 0.5) xx += 1.0;
                              if(py == 1 && yy < 0.5) yy += 1.0;
                              if(pz == 1 && zz < 0.5) zz += 1.0;
                              
                              /*-------------------------------------------
                               * Geometrical centre of the refinement
                               *-------------------------------------------*/
                              /* simple spatial centre (no weighing) */
                              spatialRef[refLevel][isoRefIndex].centreGEOM.x    += xx;
                              spatialRef[refLevel][isoRefIndex].centreGEOM.y    += yy;
                              spatialRef[refLevel][isoRefIndex].centreGEOM.z    += zz;
                              spatialRef[refLevel][isoRefIndex].centreGEOM.norm += 1.0;
                              
                              
                              /*-------------------------------------------
                               * Density centre of the refinement
                               *-------------------------------------------*/
                              /* density at node */
                              tmpDens = cur_node->dens + simu.mean_dens;
                              
                              /* double-check density value */
                              if ( tmpDens < 0.0 )
                                {
#ifdef AHFverbose
                                  fprintf(stderr,"RefCentre(): how can we have negative densites?!  l1dim=%ld  x=%g y=%g z=%g  dens=%g  mean_dens=%g\n",
                                          cur_grid->l1dim,(x+0.5)/(double)cur_grid->l1dim,(y+0.5)/(double)cur_grid->l1dim,(z+0.5)/(double)cur_grid->l1dim,cur_node->dens,simu.mean_dens);
                                  tsc_nodes[1][1][1] = cur_node;
                                  get_TSCnodes(cur_grid, cur_pquad, icur_cquad, icur_nquad, tsc_nodes, &z, &y, &x);
                                  if(test_tsc(tsc_nodes) == FALSE)
                                    fprintf(stderr,"              -> boundary node!\n");
#endif
                                  tmpDens = 0.0;
                                }
                              
#ifdef AHFmaxdenscentre
                              /* halo centre = position of max. density node */
                              if ( tmpDens > spatialRef[refLevel][isoRefIndex].maxDens )
                                {       
                                  spatialRef[refLevel][isoRefIndex].centreDens.x    = xx;
                                  spatialRef[refLevel][isoRefIndex].centreDens.y    = yy;
                                  spatialRef[refLevel][isoRefIndex].centreDens.z    = zz;
                                  spatialRef[refLevel][isoRefIndex].centreDens.norm = 1.0;
                                  
                                  /* store maximum density */
                                  spatialRef[refLevel][isoRefIndex].maxDens = tmpDens;
                                } 
                              
#else /* AHFmaxdenscentre */
                              
                              /* halo centre = density weighted centre of isolated refinement */
                              spatialRef[refLevel][isoRefIndex].centreDens.x      += xx*tmpDens;
                              spatialRef[refLevel][isoRefIndex].centreDens.y      += yy*tmpDens;
                              spatialRef[refLevel][isoRefIndex].centreDens.z      += zz*tmpDens;
                              spatialRef[refLevel][isoRefIndex].centreDens.norm   += tmpDens;
                              
                              /* store maximum density */
                              if ( tmpDens > spatialRef[refLevel][isoRefIndex].maxDens )                           
                                spatialRef[refLevel][isoRefIndex].maxDens = tmpDens;
                              
#endif /* AHFmaxdenscentre */
                              
                              
                              /*-------------------------------------------
                               * Potential centre of the refinement
                               *-------------------------------------------*/
                              /* potential at node (switch on -DAHFpotcentre to get meaningful pot-values!) */
                              tmpPot = cur_node->pot;
                              if ( tmpPot > 0.0 || tmpDens < 0.0)
                                tmpPot = 0.0;
                              else
                                tmpPot = fabs(tmpPot);
                              
                              /* halo centre = potential weighted centre of isolated refinement */
                              spatialRef[refLevel][isoRefIndex].centrePot.x      += xx*tmpPot;
                              spatialRef[refLevel][isoRefIndex].centrePot.y      += yy*tmpPot;
                              spatialRef[refLevel][isoRefIndex].centrePot.z      += zz*tmpPot;
                              spatialRef[refLevel][isoRefIndex].centrePot.norm   += tmpPot;                        
                              
                              
                              /* Particles in the refinement (Number and linked list) */
                              /* The current node has particles on it */
                              if ( cur_node->ll != NULL ) 
                                { 
                                  /* It is not the first node in this spatial ref to have particels of it */
                                  if ( spatialRef[refLevel][isoRefIndex].ll != NULL ) 
                                    { 
                                      
                                      tmpll = spatialRef[refLevel][isoRefIndex].ll;
                                      spatialRef[refLevel][isoRefIndex].ll = cur_node->ll;									
                                      
                                      /* Count the number of particles and find the last particle */	
                                      current  = cur_node->ll;
                                      previous = current;
                                      while ( current!=NULL ) 
                                        { 
                                          /*-------------------------------------------
                                           * Centre-Of-Mass of particles on refinement
                                           *-------------------------------------------*/
                                          spatialRef[refLevel][isoRefIndex].centreCMpart.x      += current->pos[X];
                                          spatialRef[refLevel][isoRefIndex].centreCMpart.y      += current->pos[Y];
                                          spatialRef[refLevel][isoRefIndex].centreCMpart.z      += current->pos[Z];
                                          spatialRef[refLevel][isoRefIndex].centreCMpart.norm   += 1.0;
                                          
                                          spatialRef[refLevel][isoRefIndex].numParts++;
                                          previous = current;
                                          current 	= current->ll;
                                        } 
                                      /* Now current is the last partical in this nodes linked list */
                                      previous->ll = tmpll;
                                      tmpll        = NULL;
                                      /* cur_node->ll = NULL; */
                                      
                                    } 
                                  else 
                                    { 
                                      current = cur_node->ll;
                                      while ( current!=NULL ) 
                                        { 
                                          /*-------------------------------------------
                                           * Centre-Of-Mass of particles on refinement
                                           *-------------------------------------------*/
                                          spatialRef[refLevel][isoRefIndex].centreCMpart.x      += current->pos[X];
                                          spatialRef[refLevel][isoRefIndex].centreCMpart.y      += current->pos[Y];
                                          spatialRef[refLevel][isoRefIndex].centreCMpart.z      += current->pos[Z];
                                          spatialRef[refLevel][isoRefIndex].centreCMpart.norm   += 1.0;
                                          
                                          spatialRef[refLevel][isoRefIndex].numParts++;
                                          current 	= current->ll;
                                        } 
                                      
                                      spatialRef[refLevel][isoRefIndex].ll = cur_node->ll;
                                      
                                      /* STU :: Note sure if I should set this equal to NULL?? b/c the next step needs the info */
                                      /* cur_node->ll = NULL; */
                                      
                                    } 
                                } 
                              
                            }
                        }
                    }
                }
            }
        }
    }	
  
  /***************************************************************************/
  /* Finishing the Spatial centre calculation and catching funny refinements
   */
  densZero=NULL;
#ifdef WITH_OPENMP
  /* TODO: this loop may be parallelized? */
#endif
  for ( i=0; i<num_refgrids; i++) 
    {
      for ( j=0; j<numIsoRef[i]; j++) 
        {
          
          /*-------------------------------------------------
           * normalisation of geometrical centre
           *-------------------------------------------------*/
          if(spatialRef[i][j].centreGEOM.norm > 0)
            {
              spatialRef[i][j].centreGEOM.x = f1mod(spatialRef[i][j].centreGEOM.x/spatialRef[i][j].centreGEOM.norm +1.0, 1.0);
              spatialRef[i][j].centreGEOM.y = f1mod(spatialRef[i][j].centreGEOM.y/spatialRef[i][j].centreGEOM.norm +1.0, 1.0);
              spatialRef[i][j].centreGEOM.z = f1mod(spatialRef[i][j].centreGEOM.z/spatialRef[i][j].centreGEOM.norm +1.0, 1.0);
            }
          else
            {
#ifdef AHFverbose
              fprintf(stderr,"RefCentre():  centreGEOM  -> spatialRef[%d][%d] numNodes=%d numParts=%d geom.norm=%g\n",i,j,spatialRef[i][j].numNodes,spatialRef[i][j].numParts,spatialRef[i][j].centreGEOM.norm);
#endif
            }
          
          /*-------------------------------------------------
           * normalisation of density weighted centre
           *-------------------------------------------------*/
          if(spatialRef[i][j].centreDens.norm > 0)
            {
              spatialRef[i][j].centreDens.x = f1mod(spatialRef[i][j].centreDens.x/spatialRef[i][j].centreDens.norm +1.0, 1.0);
              spatialRef[i][j].centreDens.y = f1mod(spatialRef[i][j].centreDens.y/spatialRef[i][j].centreDens.norm +1.0, 1.0);
              spatialRef[i][j].centreDens.z = f1mod(spatialRef[i][j].centreDens.z/spatialRef[i][j].centreDens.norm +1.0, 1.0);
            }
          else
            {
#ifdef AHFverbose
              fprintf(stderr,"RefCentre():  centreDens   -> spatialRef[%d][%d] numNodes=%d numParts=%d dens.norm=%g\n",i,j,spatialRef[i][j].numNodes,spatialRef[i][j].numParts,spatialRef[i][j].centreDens.norm);
#endif
              
              /* assign the only meaningful centre */
              spatialRef[i][j].centreDens.x = spatialRef[i][j].centreGEOM.x;
              spatialRef[i][j].centreDens.y = spatialRef[i][j].centreGEOM.y;
              spatialRef[i][j].centreDens.z = spatialRef[i][j].centreGEOM.z;            
            }
          
          /*-------------------------------------------------
           * normalisation of particles' centre-of-mass
           *-------------------------------------------------*/
          if(spatialRef[i][j].centreCMpart.norm > 0)
            {
              spatialRef[i][j].centreCMpart.x = f1mod(spatialRef[i][j].centreCMpart.x/spatialRef[i][j].centreCMpart.norm +1.0, 1.0);
              spatialRef[i][j].centreCMpart.y = f1mod(spatialRef[i][j].centreCMpart.y/spatialRef[i][j].centreCMpart.norm +1.0, 1.0);
              spatialRef[i][j].centreCMpart.z = f1mod(spatialRef[i][j].centreCMpart.z/spatialRef[i][j].centreCMpart.norm +1.0, 1.0);
            }
          else
            {
#ifdef AHFverbose
              fprintf(stderr,"RefCentre():  centreCMpart -> spatialRef[%d][%d] numNodes=%d numParts=%ld CMpart.norm=%g\n",i,j,spatialRef[i][j].numNodes,spatialRef[i][j].numParts,spatialRef[i][j].centreCMpart.norm);
#endif
              
              /* assign the only meaningful centre */
              spatialRef[i][j].centreCMpart.x = spatialRef[i][j].centreGEOM.x;
              spatialRef[i][j].centreCMpart.y = spatialRef[i][j].centreGEOM.y;
              spatialRef[i][j].centreCMpart.z = spatialRef[i][j].centreGEOM.z;            
            }
          
          
          /*-------------------------------------------------
           * normalisation of potential centre
           *-------------------------------------------------*/
          if(spatialRef[i][j].centrePot.norm > 0)
            {
              spatialRef[i][j].centrePot.x = f1mod(spatialRef[i][j].centrePot.x/spatialRef[i][j].centrePot.norm +1.0, 1.0);
              spatialRef[i][j].centrePot.y = f1mod(spatialRef[i][j].centrePot.y/spatialRef[i][j].centrePot.norm +1.0, 1.0);
              spatialRef[i][j].centrePot.z = f1mod(spatialRef[i][j].centrePot.z/spatialRef[i][j].centrePot.norm +1.0, 1.0);
            }
          else
            {
#if (defined AHFverbose && defined AHFpotcentre)
              fprintf(stderr,"RefCentre():  centrePot    -> spatialRef[%d][%d] numNodes=%d numParts=%d pot.norm=%g\n",i,j,spatialRef[i][j].numNodes,spatialRef[i][j].numParts,spatialRef[i][j].centrePot.norm);
#endif
              
              /* assign the only meaningful centre */
              spatialRef[i][j].centrePot.x = spatialRef[i][j].centreGEOM.x;
              spatialRef[i][j].centrePot.y = spatialRef[i][j].centreGEOM.y;
              spatialRef[i][j].centrePot.z = spatialRef[i][j].centreGEOM.z;            
            }
          
          
          /*-------------------------------------------------
           * catch spatialRef's with unphysical densities
           *-------------------------------------------------*/
          if ( spatialRef[i][j].maxDens <= ZERO ) 
            {
              
              /* store these spatialRef's in a separate list */
              numDensZero++;
              
              if(densZero == NULL)
                {
                  if ((densZero = calloc(numDensZero+1, sizeof(SRINDEX)))==NULL) 
                    {
                      fprintf(stderr,"calloc failed for densZero\n");
                      exit(-1);
                    }
                }
              else
                {
                  if ((densZero = realloc(densZero,(numDensZero+1)*sizeof(SRINDEX)))==NULL) 
                    {
                      fprintf(stderr,"realloc failed for densZero\n");
                      exit(-1);
                    }
                }
              
              densZero[numDensZero-1].refLevel		= i;
              densZero[numDensZero-1].isoRefIndex = j;
              
              /* assign the only meaningful centre */
              spatialRef[i][j].centreDens.x = spatialRef[i][j].centreGEOM.x;
              spatialRef[i][j].centreDens.y = spatialRef[i][j].centreGEOM.y;
              spatialRef[i][j].centreDens.z = spatialRef[i][j].centreGEOM.z;            
            }
          
          /*-------------------------------------------------
           * which centre to use as the halo centre?
           *-------------------------------------------------*/
          /* the default = density weighted centre */
          spatialRef[i][j].centre.x = spatialRef[i][j].centreDens.x;
          spatialRef[i][j].centre.y = spatialRef[i][j].centreDens.y;
          spatialRef[i][j].centre.z = spatialRef[i][j].centreDens.z;
          
#ifdef AHFpotcentre
          spatialRef[i][j].centre.x = spatialRef[i][j].centrePot.x;
          spatialRef[i][j].centre.y = spatialRef[i][j].centrePot.y;
          spatialRef[i][j].centre.z = spatialRef[i][j].centrePot.z;
#endif
          
#ifdef AHFcomcentre
          spatialRef[i][j].centre.x = spatialRef[i][j].centreCMpart.x;
          spatialRef[i][j].centre.y = spatialRef[i][j].centreCMpart.y;
          spatialRef[i][j].centre.z = spatialRef[i][j].centreCMpart.z;
#endif
          
#ifdef AHFgeomcentre
          spatialRef[i][j].centre.x = spatialRef[i][j].centreGEOM.x;
          spatialRef[i][j].centre.y = spatialRef[i][j].centreGEOM.y;
          spatialRef[i][j].centre.z = spatialRef[i][j].centreGEOM.z;
#endif
          
        }
    }
  
  /***************************************************************************/
  /* Gathering information needed to calculate the  boundary of the refinement
   * Calculating the boundRefDiv
   */
  alpha = 1.1;
#ifdef WITH_OPENMP
  /* TODO: this loop may be parallelized? */
#endif
  for ( i=0; i<num_refgrids; i++) 
    {
      boxLen = ((1.0)/((double)(gridl1dim[i])));
      boxVol = boxLen*boxLen*boxLen;
      
      for ( j=0; j<numIsoRef[i]; j++) 
        {
          colour = spatialRef[i][j].colour;
          
          spatialRef[i][j].volume = spatialRef[i][j].numNodes*boxVol;
          
          px = spatialRefIndex[colour].periodic.x;
          py = spatialRefIndex[colour].periodic.y;
          pz = spatialRefIndex[colour].periodic.z;
          
          /* Is this a periodic isolated refinement */
          if ( (px == 1) || (py == 1) || (pz == 1) ) 
            {
              tmpRad = (3.0*spatialRef[i][j].volume)/(4*PI);
              tmpRad = pow(tmpRad,0.333333333);
              tmpRad = tmpRad*alpha;
              
              if ( px == 1 ) 
                {
                  a = tmpRad + spatialRef[i][j].centreDens.x;
                  b = 1.0 - tmpRad + spatialRef[i][j].centreDens.x;
                  
                  if ( b < a ) 
                    {
#ifdef AHFverbose
                      fprintf(stderr,"'boundRefDiv' b < a The refinement is too big!! leave set min/max as 0:1\n");
#endif
                      spatialRef[i][j].boundRefDiv.x = -1.0;
                    }
                  
                  c = (a+b)/2.0;
                  spatialRef[i][j].boundRefDiv.x = fmod(c, 1.0);
                  
                } 
              
              if ( py == 1 ) 
                {
                  a = tmpRad + spatialRef[i][j].centreDens.y;
                  b = 1.0 - tmpRad + spatialRef[i][j].centreDens.y;
                  
                  if ( b < a ) 
                    {
#ifdef AHFverbose
                      fprintf(stderr,"'boundRefDiv' b < a The refinement is too big!! leave set min/max as 0:1\n");
#endif
                      spatialRef[i][j].boundRefDiv.y = -1.0;
                    }
                  
                  c = (a+b)/2.0;
                  spatialRef[i][j].boundRefDiv.y = fmod(c, 1.0);
                  
                } 
              
              if ( pz == 1 )
                {
                  a = tmpRad + spatialRef[i][j].centreDens.z;
                  b = 1.0 - tmpRad + spatialRef[i][j].centreDens.z;
                  
                  if ( b < a ) 
                    {
#ifdef AHFverbose
                      fprintf(stderr,"'boundRefDiv' b < a The refinement is too big!! leave set min/max as 0:1\n");
#endif
                      spatialRef[i][j].boundRefDiv.z = -1.0;
                    }
                  
                  c = (a+b)/2.0;
                  spatialRef[i][j].boundRefDiv.z = fmod(c, 1.0);
                  
                }
              
            } 
        }
    }
  
  
  /***************************************************************************/
  /* Collecting information about isolated refinements
   * 
   * boundary of the refinement
   */
  for( 	iterate=0, cur_grid=global.dom_grid+ahf.min_ref; iterate<ahf.no_grids; iterate++, cur_grid++) 
    {
      /* shift of cell centre as compared to edge of box [grid units] */
      cur_shift = 0.5/(double)cur_grid->l1dim;
      
      for(cur_pquad=cur_grid->pquad; cur_pquad != NULL; cur_pquad=cur_pquad->next)
        {
          z = cur_pquad->z;
          for(cur_cquad = cur_pquad->loc;
              cur_cquad < cur_pquad->loc + cur_pquad->length; 
              cur_cquad++, z++)  
            {  
              for(icur_cquad  = cur_cquad; 
                  icur_cquad != NULL; 
                  icur_cquad  = icur_cquad->next)
                {
                  y = icur_cquad->y;
                  for(cur_nquad = icur_cquad->loc;  
                      cur_nquad < icur_cquad->loc + icur_cquad->length; 
                      cur_nquad++, y++) 
                    { 
                      for(icur_nquad  = cur_nquad; 
                          icur_nquad != NULL; 
                          icur_nquad  = icur_nquad->next)
                        {
                          x = icur_nquad->x;
                          for(cur_node = icur_nquad->loc; 
                              cur_node < icur_nquad->loc + icur_nquad->length; 
                              cur_node++, x++)
                            {
                              
                              /* Locate the spatial refinement */
                              refLevel    = spatialRefIndex[cur_node->force.colour].refLevel;
                              isoRefIndex = spatialRefIndex[cur_node->force.colour].isoRefIndex;
                              
                              /* Position of the node */
                              xx = (double) x/(double)cur_grid->l1dim + cur_shift;
                              xx = (double)fmod(xx + 1.0, 1.0);
                              yy = (double) y/(double)cur_grid->l1dim + cur_shift;
                              yy = (double)fmod(yy + 1.0, 1.0);
                              zz = (double) z/(double)cur_grid->l1dim + cur_shift;
                              zz = (double)fmod(zz + 1.0, 1.0);
                              
                              
                              /*  The extent of the refinement */
                              /* NOTE :: For periodic boundary isolate refinements max < min !!!!  */
                              
                              /**********************************************************/
                              
                              /* It is not a boundary refinement */
                              if ( spatialRef[refLevel][isoRefIndex].boundRefDiv.x  < 0.0 )
                                {
                                  
                                  tmpMinMax = MinMax(xx, spatialRef[refLevel][isoRefIndex].x.min,
                                                     spatialRef[refLevel][isoRefIndex].x.max);
                                  spatialRef[refLevel][isoRefIndex].x.min = tmpMinMax.min;
                                  spatialRef[refLevel][isoRefIndex].x.max = tmpMinMax.max;
                                  
                                } 
                              /* It is a boundary Refinement  */
                              else 
                                { 
                                  
                                  tmpMinMax = MinMaxBound(spatialRef[refLevel][isoRefIndex].boundRefDiv.x, xx, spatialRef[refLevel][isoRefIndex].x.min,
                                                          spatialRef[refLevel][isoRefIndex].x.max);
                                  spatialRef[refLevel][isoRefIndex].x.min = tmpMinMax.min;
                                  spatialRef[refLevel][isoRefIndex].x.max = tmpMinMax.max;
                                  
                                }
                              
                              
                              /**********************************************************/
                              /* It is not a boundary refinement */
                              if ( spatialRef[refLevel][isoRefIndex].boundRefDiv.y  < 0.0 ) 
                                {
                                  
                                  tmpMinMax = MinMax(yy, spatialRef[refLevel][isoRefIndex].y.min,
                                                     spatialRef[refLevel][isoRefIndex].y.max);
                                  spatialRef[refLevel][isoRefIndex].y.min = tmpMinMax.min;
                                  spatialRef[refLevel][isoRefIndex].y.max = tmpMinMax.max;
                                  
                                } 
                              /* It is a boundary Refinement  */
                              else 
                                { 
                                  
                                  tmpMinMax = MinMaxBound(spatialRef[refLevel][isoRefIndex].boundRefDiv.y, yy, spatialRef[refLevel][isoRefIndex].y.min,
                                                          spatialRef[refLevel][isoRefIndex].y.max);
                                  spatialRef[refLevel][isoRefIndex].y.min = tmpMinMax.min;
                                  spatialRef[refLevel][isoRefIndex].y.max = tmpMinMax.max;
                                }
                              
                              /**********************************************************/
                              /* It is not a boundary refinement */
                              if ( spatialRef[refLevel][isoRefIndex].boundRefDiv.z  < 0.0 )
                                {
                                  
                                  tmpMinMax = MinMax(zz, spatialRef[refLevel][isoRefIndex].z.min,
                                                     spatialRef[refLevel][isoRefIndex].z.max);
                                  spatialRef[refLevel][isoRefIndex].z.min = tmpMinMax.min;
                                  spatialRef[refLevel][isoRefIndex].z.max = tmpMinMax.max;
                                  
                                } 
                              /* It is a boundary Refinement  */
                              else 
                                {
                                  
                                  tmpMinMax = MinMaxBound(spatialRef[refLevel][isoRefIndex].boundRefDiv.z, zz, spatialRef[refLevel][isoRefIndex].z.min,
                                                          spatialRef[refLevel][isoRefIndex].z.max);
                                  spatialRef[refLevel][isoRefIndex].z.min = tmpMinMax.min;
                                  spatialRef[refLevel][isoRefIndex].z.max = tmpMinMax.max;
                                }
                              
                            }
                        }
                    }
                }
            }
        }
    }	
  
  /***************************************************************************/
  /* Catching the case where the refinement is periodic but does not extend beyond the end domain node
   * I.e. the periodic refinement is too small!!  refer to the 1000 case 
   * 
   */
  
  for ( i=0; i<num_refgrids; i++) 
    {
      for ( j=0; j<numIsoRef[i]; j++) 
        {
          
          if ( spatialRef[i][j].x.min == 100000.0 ) 
            spatialRef[i][j].x.min = 0.0;
          
          if ( spatialRef[i][j].x.max == -100000.0 )
            spatialRef[i][j].x.max = 1.0;
          
          if ( spatialRef[i][j].y.min == 100000.0 ) 
            spatialRef[i][j].y.min = 0.0;
          
          if ( spatialRef[i][j].y.max == -100000.0 )
            spatialRef[i][j].y.max = 1.0;
          
          if ( spatialRef[i][j].z.min == 100000.0 ) 
            spatialRef[i][j].z.min = 0.0;
          
          if ( spatialRef[i][j].z.max == -100000.0 )
            spatialRef[i][j].z.max = 1.0;
        }
    }
  
  
  /* Freeing the spatialRefIndex array */
  free(spatialRefIndex);
  spatialRefIndex = NULL;
  
  return(TRUE);
  
}


int	analyseRef(int num_refgrids, SPATIALREF **spatialRef)
{
  
  int i,j,k,l,p;
  double	xmin,xmax,ymin,ymax,zmin,zmax;
  double	x,y,z;
  
  int	DETAILSWEEP=0;
  int isoRefIndex;
  
  int isoRefIndexNew, isoRefIndexOLD;
  double dx, dy, dz;
  double dist, tmpMin;
  long   maxParts, maxNodes;
  
  double alpha;	/* caution factor */
  
  gridls		*cur_grid;
  double 	xx,yy,zz;
  double tmpRad;
  
  int EMBEDREF;
  
  int	count1, count2, count3;
  int *tmpArray;
  int tmpCount,totalcount;
  
  int numcomplex;
  
  int *tmpIsoRefIndex;
  int tmpc;
  
  int hirefLevel, hiisoRefIndex;
  int tmpisoRefIndex, tmprefLevel;
  
  int minNodes = 10000000;
  
  /* Isolated refinement under investiagation */
  hirefLevel = 3; hiisoRefIndex = 134;
  
  
  /***************************************************************************/
  /* What is the minimum number of nodes in an isolated refinement? */
  for ( i=0; i<num_refgrids; i++)
    {
      for ( j=0; j<numIsoRef[i]; j++)
        {
          if ( spatialRef[i][j].numNodes < minNodes )
            minNodes = spatialRef[i][j].numNodes;
        }
    }
#ifdef AHFverbose
  fprintf(stderr,"#### Minimum number of nodes = %d\n",minNodes);
#endif
  
  /***************************************************************************/
  /* Searching for embedded refinements */
  for ( i=0; i<num_refgrids-1; i++)
    {
      for ( j=0; j<numIsoRef[i]; j++)
        {
          
          xmin = spatialRef[i][j].x.min;
          xmax = spatialRef[i][j].x.max;
          ymin = spatialRef[i][j].y.min;
          ymax = spatialRef[i][j].y.max;
          zmin = spatialRef[i][j].z.min;
          zmax = spatialRef[i][j].z.max;
          
          /* Is the centre of a finer refinement in this refinement? */
          for ( k=0;k<numIsoRef[i+1];k++)
            {
              
              /* Note we could have used the CofMass centres */
              x = spatialRef[i+1][k].centreDens.x;
              y = spatialRef[i+1][k].centreDens.y;
              z = spatialRef[i+1][k].centreDens.z;
              
              
              /* If yes record the occurance */
              EMBEDREF = 0;
              
              if ( xmin < xmax )
                {
                  
                  if ( (x>xmin) && (x<xmax) )
                    EMBEDREF = 1;		
                  
                }
              else
                {
                  
                  if ( (x>=0) && (x<xmax) )
                    EMBEDREF = 1;		
                  
                  if ( (x>xmin) && (x<=1.0) )
                    EMBEDREF = 1;		
                  
                }	
              if ( EMBEDREF == 1 )
                {
                  
                  if ( ymin < ymax ) 
                    {
                      
                      if ( (y>ymin) && (y<ymax) )
                        EMBEDREF = 2;		
                      
                    } 
                  else 
                    {
                      
                      if ( (y>=0) && (y<ymax) )
                        EMBEDREF = 2;		
                      
                      if ( (y>ymin) && (y<=1.0) )
                        EMBEDREF = 2;		
                      
                    }	
                  
                }
              if ( EMBEDREF == 2 ) 
                {
                  
                  if ( zmin < zmax ) 
                    {
                      
                      if ( (z>zmin) && (z<zmax) )
                        EMBEDREF = 3;		
                      
                    } 
                  else 
                    {
                      
                      if ( (z>=0) && (z<zmax) )
                        EMBEDREF = 3;		
                      
                      if ( (z>zmin) && (z<=1.0) )
                        EMBEDREF = 3;		
                      
                    }	
                  
                }
              if ( EMBEDREF == 3 ) 
                {
                  
                  
                  /***************************************************************************************************
                   * Recording the old substruct info */
                  if ( spatialRef[i][j].numSubStruct > 0 )
                    {
                      
                      tmpArray=NULL;
                      if ((tmpArray = calloc(spatialRef[i][j].numSubStruct,sizeof(int)))==NULL)
                        {
                          fprintf(stderr,"Error in allocating the memory for Recording the old substruct info\n");
                          exit(0);
                        }
                      
                      for ( l=0; l<spatialRef[i][j].numSubStruct; l++)
                        {
                          tmpArray[l] = spatialRef[i][j].subStruct[l].isoRefIndex;
                        }
                    }
                  
                  /* Adding the new substructure */
                  spatialRef[i][j].numSubStruct++;
                  
                  if(spatialRef[i][j].subStruct == NULL)
                    {
                      if ((spatialRef[i][j].subStruct = calloc(spatialRef[i][j].numSubStruct+1, sizeof(INDEX)))==NULL)
                        {
                          fprintf(stderr,"calloc failed in connecting nodes\n");
                          exit(-1);
                        }
                    }
                  else
                    {
                      if ((spatialRef[i][j].subStruct = realloc(spatialRef[i][j].subStruct,(spatialRef[i][j].numSubStruct+1)*sizeof(INDEX)))==NULL)
                        {
                          fprintf(stderr,"realloc failed in connecting nodes\n");
                          exit(-1);
                        }
                    }
                  
                  if ( spatialRef[i][j].numSubStruct-1 > 0 )
                    {
                      for ( l=0; l<spatialRef[i][j].numSubStruct-1; l++)
                        { 
                          spatialRef[i][j].subStruct[l].refLevel    = i+1;
                          spatialRef[i][j].subStruct[l].isoRefIndex = tmpArray[l];
                        }
                      free(tmpArray);
                    }
                  
                  spatialRef[i][j].subStruct[spatialRef[i][j].numSubStruct-1].refLevel    = i+1;
                  spatialRef[i][j].subStruct[spatialRef[i][j].numSubStruct-1].isoRefIndex = k;
                  
                  
                  /***************************************************************************************************
                   * Recording the Parent domain information 
                   * i+1 = finer refinment
                   * k   = looping through the isolated refinements on this level */
                  if ( spatialRef[i+1][k].numParDom > 0 )
                    {
                      
                      tmpArray=NULL;
                      if ((tmpArray = calloc(spatialRef[i+1][k].numParDom,sizeof(int)))==NULL)
                        {
                          fprintf(stderr,"Error in allocating the memory for Recording the Parent domain information\n");
                          exit(0);
                        }
                      
                      for ( l=0; l<spatialRef[i+1][k].numParDom; l++) 
                        tmpArray[l] = spatialRef[i+1][k].parDom[l].isoRefIndex;
                    }
                  
                  /*  Recording the Parent domain */
                  spatialRef[i+1][k].numParDom++;
                  
                  if(spatialRef[i+1][k].parDom == NULL)
                    {
                      if ((spatialRef[i+1][k].parDom = calloc(spatialRef[i+1][k].numParDom+1, sizeof(INDEX)))==NULL)
                        {
                          fprintf(stderr,"calloc failed in connecting nodes\n");
                          exit(-1);
                        }
                    }
                  else
                    {
                      if ((spatialRef[i+1][k].parDom = realloc(spatialRef[i+1][k].parDom,(spatialRef[i+1][k].numParDom+1)*sizeof(INDEX)))==NULL)
                        {
                          fprintf(stderr,"realloc failed in connecting nodes\n");
                          exit(-1);
                        }
                    }
                  
                  if ( spatialRef[i+1][k].numParDom-1 > 0 )
                    {
                      for ( l=0; l<spatialRef[i+1][k].numParDom-1; l++)
                        {
                          spatialRef[i+1][k].parDom[l].refLevel    = i;
                          spatialRef[i+1][k].parDom[l].isoRefIndex = tmpArray[l];
                        }
                      free(tmpArray);
                    }
                  
                  spatialRef[i+1][k].parDom[spatialRef[i+1][k].numParDom-1].refLevel    = i;
                  spatialRef[i+1][k].parDom[spatialRef[i+1][k].numParDom-1].isoRefIndex = j;
                  
                  /* Flaging if we need to do a complex sweep */
                  if ( spatialRef[i+1][k].numParDom > 1 )
                    {
                      DETAILSWEEP = 1;
                    }
                  
                  
                }
            }	
        }
    }
  
  /***************************************************************************/
  /* If the refinements are tricky then the simply algorithm above will not work */
  /* but this one will :) */
  numcomplex = 0;
  if ( DETAILSWEEP == 1 )
    {
#ifdef AHFverbose		
      fprintf(stderr,"analyseRef(): we need to do a complex sweep to find embedded halos\n");
#endif
      
      for ( i=1; i<num_refgrids-1; i++)
        {
          for ( j=0; j<numIsoRef[i]; j++)
            {
              
              if ( spatialRef[i][j].numParDom > 1 )
                {
                  
                  numcomplex++;
                  
                  /***********************************************************/
                  /* parent halo = halo on next coarser level
                   *               with centre closest to current spatialRef */
                  tmpMin   = 10000000000000.0;
                  maxNodes = 0;
                  maxParts = 0;
                  for ( k=0; k<spatialRef[i][j].numParDom; k++)
                    {
                      
                      isoRefIndex    = spatialRef[i][j].parDom[k].isoRefIndex;
                      isoRefIndexNew = -1;
                      
                      /* calculate the distance to the parent */
                      dx = spatialRef[i][j].centreDens.x - spatialRef[i-1][isoRefIndex].centreDens.x;
                      dx = fabs(dx);	
                      
                      dy = spatialRef[i][j].centreDens.y - spatialRef[i-1][isoRefIndex].centreDens.y;
                      dy = fabs(dy);	
                      
                      dz = spatialRef[i][j].centreDens.z - spatialRef[i-1][isoRefIndex].centreDens.z;
                      dz = fabs(dz);	
                      
                      if (dx > 0.5 ) dx = 1.0 - dx;
                      if (dy > 0.5 ) dy = 1.0 - dy;
                      if (dz > 0.5 ) dz = 1.0 - dz;
                      
                      dist = dx*dx + dy*dy + dz*dz;
                      
                      /* isoRefIndexNew :: is the RefIndex of the parent refinement. I.e. the refinement on the coarser level that is the closest */
                      if ( dist < tmpMin )
                        {
                          isoRefIndexNew = isoRefIndex;
                          tmpMin         = dist;
                        }
                    }
                  
#ifdef AHFverbose
                  if(isoRefIndex < 0)
                    fprintf(stderr,"analyseRef(): WARNING ->  no parent   found for spatialRef[%d][%d]  "
                            "x=%g y=%g z=%g   numPart=%lu  numNodes=%d\n",
                            i, j,
                            spatialRef[i][j].centreDens.x,
                            spatialRef[i][j].centreDens.y,
                            spatialRef[i][j].centreDens.z,
                            spatialRef[i][j].numParts,
                            spatialRef[i][j].numNodes);
#endif
                  
                  /*******************************************************/
                  /* Correct the substructure listings */
                  
                  /* Run through the parents that list this substructure */
                  tmpc=0;
                  for ( k=0; k<spatialRef[i][j].numParDom; k++)
                    { 
                      
                      /* Name of parent to check */			
                      isoRefIndex = spatialRef[i][j].parDom[k].isoRefIndex;
                      
                      /* If this refinment is not the real parent then remove this sub-halo from it's substructure listing */
                      /* remove sub-halo I.e. (i,j) from the parent halo I.e. (i-1,isoRefIndex) substructure listing */
                      if ( isoRefIndex != isoRefIndexNew )
                        {
                          
                          /* Create tmp array of the names of the sub halos for this parent*/
                          tmpArray=NULL;
                          if ((tmpArray = calloc(spatialRef[i-1][isoRefIndex].numSubStruct,sizeof(int)))==NULL)
                            {
                              fprintf(stderr,"Error in allocating the memory for remove substrcture array\n");
                              exit(0);
                            }
                          
                          for ( p=0; p<spatialRef[i-1][isoRefIndex].numSubStruct; p++) 
                            tmpArray[p] = spatialRef[i-1][isoRefIndex].subStruct[p].isoRefIndex;
                          
                          /* removing the substructure (i,j) from the parent substructure listing */
                          spatialRef[i-1][isoRefIndex].numSubStruct = spatialRef[i-1][isoRefIndex].numSubStruct - 1;
                          free(spatialRef[i-1][isoRefIndex].subStruct);
                          spatialRef[i-1][isoRefIndex].subStruct = NULL;
                          
                          /* QUICK AND DIRTY FIX ADDED BY AK ON 30/09/2005 */
                          if(spatialRef[i-1][isoRefIndex].numSubStruct == 0)
                            {
#ifdef AHFverbose
                              fprintf(stderr,"   -> NO MORE SUBSTRUCTURE LEFT IN PARENT HALO (%d, %d)\n",i-1,isoRefIndex);
#endif
                              /* however, the loop over p right below requires at least on INDEX structure */
                              if ((spatialRef[i-1][isoRefIndex].subStruct = calloc(1,sizeof(INDEX)))==NULL)
                                {
                                  fprintf(stderr,"calloc failed in substructure\n");
                                  exit(-1);
                                }
                            }
                          else
                            {
                              if ((spatialRef[i-1][isoRefIndex].subStruct = calloc(spatialRef[i-1][isoRefIndex].numSubStruct,sizeof(INDEX)))==NULL)
                                {
                                  fprintf(stderr,"calloc failed in substructure\n");
                                  exit(-1);
                                }
                            }
#ifdef OLD_VERSION
                          if ((spatialRef[i-1][isoRefIndex].subStruct = calloc(spatialRef[i-1][isoRefIndex].numSubStruct,sizeof(INDEX)))==NULL)
                            {
                              fprintf(stderr,"calloc failed in substructure\n");
                              exit(-1);
                            }
#endif                     
                          
                          tmpCount =0;	
                          for ( p=0; p<spatialRef[i-1][isoRefIndex].numSubStruct+1; p++)
                            {
                              /* I.e. include all except our substructure */
                              if ( tmpArray[p] != j )	
                                {	
                                  spatialRef[i-1][isoRefIndex].subStruct[tmpCount].isoRefIndex = tmpArray[p];
                                  spatialRef[i-1][isoRefIndex].subStruct[tmpCount].refLevel = i;
                                  tmpCount++;
                                } 
                            }
                          
                        }
                      
                      free(tmpArray);
                      tmpArray = NULL;
                    }
                  
                  /*******************************************************/
                  /* Finalise the parent halo */
                  spatialRef[i][j].numParDom = 1;
                  
                  if(spatialRef[i][j].parDom == NULL)
                    {
                      if ((spatialRef[i][j].parDom = 
                           calloc(spatialRef[i][j].numParDom+1, sizeof(INDEX)))==NULL) 
                        {
                          fprintf(stderr,"calloc failed in connecting nodes\n");
                          exit(-1);
                        }
                    }
                  else
                    {
                      if ((spatialRef[i][j].parDom = 
                           realloc(spatialRef[i][j].parDom,(spatialRef[i][j].numParDom+1)*sizeof(INDEX)))==NULL) 
                        {
                          fprintf(stderr,"realloc failed in connecting nodes\n");
                          exit(-1);
                        }
                    }
                  
                  /* only store parent if we actually found one */
                  if( isoRefIndexNew > 0 )
                    {
                      spatialRef[i][j].parDom[0].isoRefIndex = isoRefIndexNew;
                      spatialRef[i][j].parDom[0].refLevel    = i-1;
                    }
                  /* otherwise we will try a second time just below... */
                  else
                    {
                      if(spatialRef[i][j].parDom != NULL)
                        free(spatialRef[i][j].parDom);
                      
                      spatialRef[i][j].numParDom = 0;
                      spatialRef[i][j].parDom    = NULL;
                    }
                  
                }
            }
        }
    }
  
  
  
  /*********************************************************************************************************************************************************/
  /*******************************************************************************************/
  /* Double checking
   * Picking up halos that don't have parents */
  for ( i=1; i<num_refgrids; i++)
    {
      for ( j=0; j<numIsoRef[i]; j++)
        {
          
          
          if ( spatialRef[i][j].parDom == NULL )
            { 
              
#ifdef AHFDEBUG		
              fprintf(stderr,"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
              fprintf(stderr,"analyseRef(): no parent for spatialRef[%d][%d].numParDom = %d\n",i,j,spatialRef[i][j].numParDom);
#endif
              
              /* What is the closest potential parent refinement? */
              x = spatialRef[i][j].centreDens.x;
              y = spatialRef[i][j].centreDens.y;
              z = spatialRef[i][j].centreDens.z;
              tmpMin         = 10000000000000.0;
              tmpisoRefIndex = -1.0;
              tmprefLevel    = -1.0;
              
              /* loop over all spatial refinements on next coarser level */
              for ( p=0; p<numIsoRef[i-1]; p++ )
                {
                  /* irrespective of "#ifdef PARDAU_*" we search for the refinement closest in distance */
                  
                  dx = x - spatialRef[i-1][p].centreDens.x;
                  dx = fabs(dx);	
                  
                  dy = y - spatialRef[i-1][p].centreDens.y;
                  dy = fabs(dy);	
                  
                  dz = z - spatialRef[i-1][p].centreDens.z;
                  dz = fabs(dz);	
                  
                  if (dx > 0.5 ) dx = 1.0 - dx;
                  if (dy > 0.5 ) dy = 1.0 - dy;
                  if (dz > 0.5 ) dz = 1.0 - dz;
                  
                  dist = dx*dx + dy*dy + dz*dz;
                  
                  if ( dist < tmpMin ) 
                    {
                      tmpMin         = dist;
                      tmpisoRefIndex = p;
                      tmprefLevel    = i-1;
                    }
                  
                }
              
#ifdef AHFDEBUG		
              fprintf(stderr,"PARENT [%d][%d]\n",tmprefLevel,tmpisoRefIndex);	
              fprintf(stderr,"numSubStruct[%d] daughter[%d][%d]\n",
                      spatialRef[i-1][tmpisoRefIndex].numSubStruct,spatialRef[i-1][tmpisoRefIndex].daughter.isoRefIndex,spatialRef[i-1][tmpisoRefIndex].daughter.refLevel);
#endif
              /* We have found the lost parent refinement 
               * Name the parent */
              spatialRef[i][j].numParDom = 1;
              spatialRef[i][j].parDom    = NULL;
              if ((spatialRef[i][j].parDom = calloc(1,sizeof(INDEX)))==NULL)
                {
                  fprintf(stderr,"Error in allocating the memory for Recording the Parent domain information\n");
                  exit(0);
                }
              spatialRef[i][j].parDom[0].isoRefIndex = tmpisoRefIndex;
              spatialRef[i][j].parDom[0].refLevel    = tmprefLevel;
              
              
#ifdef AHFDEBUG		
              /* Update the parents substructure list */
              fprintf(stderr,"spatialRef[%d][%d].numSubStruct = %d\n",tmprefLevel,tmpisoRefIndex,spatialRef[tmprefLevel][tmpisoRefIndex].numSubStruct);	
#endif
              
              /* Create tmp array of the names of the sub halos for this parent*/
              tmpArray=NULL;
              if ((tmpArray = calloc(spatialRef[tmprefLevel][tmpisoRefIndex].numSubStruct,sizeof(int)))==NULL)
                {
                  fprintf(stderr,"Error in allocating the memory for remove substrcture array\n");
                  exit(0);
                }
              
              for ( p=0; p<spatialRef[tmprefLevel][tmpisoRefIndex].numSubStruct; p++) 
                tmpArray[p] = spatialRef[tmprefLevel][tmpisoRefIndex].subStruct[p].isoRefIndex;
              
              /* adding the substructure to the parent substructure listing */
              spatialRef[tmprefLevel][tmpisoRefIndex].numSubStruct += 1;
              free(spatialRef[tmprefLevel][tmpisoRefIndex].subStruct);
              spatialRef[tmprefLevel][tmpisoRefIndex].subStruct = NULL;
              if ((spatialRef[tmprefLevel][tmpisoRefIndex].subStruct = calloc(spatialRef[tmprefLevel][tmpisoRefIndex].numSubStruct,sizeof(INDEX)))==NULL)
                {
                  fprintf(stderr,"calloc failed in substrcture\n");
                  exit(-1);
                }
              
              for ( p=0; p<spatialRef[tmprefLevel][tmpisoRefIndex].numSubStruct-1; p++) 
                {
                  spatialRef[tmprefLevel][tmpisoRefIndex].subStruct[p].isoRefIndex = tmpArray[p];
                  spatialRef[tmprefLevel][tmpisoRefIndex].subStruct[p].refLevel = i;
                }
              
              free(tmpArray);
              tmpArray = NULL;            
              
              spatialRef[tmprefLevel][tmpisoRefIndex].subStruct[spatialRef[tmprefLevel][tmpisoRefIndex].numSubStruct-1].refLevel = i;
              spatialRef[tmprefLevel][tmpisoRefIndex].subStruct[spatialRef[tmprefLevel][tmpisoRefIndex].numSubStruct-1].isoRefIndex = j;
              
#ifdef AHFDEBUG		
              fprintf(stderr,"spatialRef[%d][%d].numSubStruct = %d\n",tmprefLevel,tmpisoRefIndex,spatialRef[tmprefLevel][tmpisoRefIndex].numSubStruct);	
              for ( p=0; p<spatialRef[tmprefLevel][tmpisoRefIndex].numSubStruct; p++)
                fprintf(stderr,"[%d][%d], ",spatialRef[tmprefLevel][tmpisoRefIndex].subStruct[p].refLevel,spatialRef[tmprefLevel][tmpisoRefIndex].subStruct[p].isoRefIndex);
              
              fprintf(stderr,"\n");
#endif
              
              
              /* Giving this lost refinement the centre from the refinement below
               * We do this because it might not actually be related.
               * If it is related then it will sort itself out - if not the partilces we be removed by unbound */
              spatialRef[i][j].centreDens.x  = spatialRef[tmprefLevel][tmpisoRefIndex].centreDens.x;
              spatialRef[i][j].centreDens.y  = spatialRef[tmprefLevel][tmpisoRefIndex].centreDens.y;
              spatialRef[i][j].centreDens.z  = spatialRef[tmprefLevel][tmpisoRefIndex].centreDens.z;
              spatialRef[i][j].centreGEOM.x = spatialRef[tmprefLevel][tmpisoRefIndex].centreGEOM.x;
              spatialRef[i][j].centreGEOM.y = spatialRef[tmprefLevel][tmpisoRefIndex].centreGEOM.y;
              spatialRef[i][j].centreGEOM.z = spatialRef[tmprefLevel][tmpisoRefIndex].centreGEOM.z;
              
              
            } /* parDom == NULL */
          
        }
    }
  
  /*********************************************************************************************************************************************************/
  /***************************************************************************/
  /* Finding the daughter refinement */
  
  /* just to make the collecting radii a little smaller */
  alpha = 0.9;
  for ( i=0; i<num_refgrids-1; i++)
    {
      for ( j=0; j<numIsoRef[i]; j++)
        {
          
          /**************************************************/
          /* Finding the daughter refinement */
          if ( spatialRef[i][j].numSubStruct > 1 )
            {
              tmpMin         = 10000000000000.0;
              maxNodes       = 0;
              maxParts       = 0;
              isoRefIndexNew = -1;
              
              for ( k=0; k<spatialRef[i][j].numSubStruct; k++)
                {
                  
                  isoRefIndex    = spatialRef[i][j].subStruct[k].isoRefIndex;
                  
#ifdef PARDAU_DISTANCE
                  dx = spatialRef[i][j].centreDens.x - spatialRef[i+1][isoRefIndex].centreDens.x;
                  dx = fabs(dx);	
                  
                  dy = spatialRef[i][j].centreDens.y - spatialRef[i+1][isoRefIndex].centreDens.y;
                  dy = fabs(dy);	
                  
                  dz = spatialRef[i][j].centreDens.z - spatialRef[i+1][isoRefIndex].centreDens.z;
                  dz = fabs(dz);	
                  
                  if (dx > 0.5 ) dx = 1.0 - dx;
                  if (dy > 0.5 ) dy = 1.0 - dy;
                  if (dz > 0.5 ) dz = 1.0 - dz;
                  
                  
                  dist = dx*dx + dy*dy + dz*dz;
                  
                  if ( dist < tmpMin ) 
                    {
                      isoRefIndexNew = isoRefIndex;
                      tmpMin         = dist;
                    }
#endif
#ifdef PARDAU_NODES
                  if(spatialRef[i+1][isoRefIndex].numNodes > maxNodes)
                    {
                      isoRefIndexNew = isoRefIndex;
                      maxNodes       = spatialRef[i+1][isoRefIndex].numNodes;
                    }
#endif
#ifdef PARDAU_PARTS
                  if(spatialRef[i+1][isoRefIndex].numParts > maxParts)
                    {
                      isoRefIndexNew = isoRefIndex;
                      maxParts       = spatialRef[i+1][isoRefIndex].numParts;
                    }
#endif
                }
              
              
              spatialRef[i][j].daughter.isoRefIndex = isoRefIndexNew;
              spatialRef[i][j].daughter.refLevel    = i+1;
              
#ifdef AHFverbose
              if( isoRefIndexNew < 0 )
                fprintf(stderr,"analyseRef(): WARNING ->  no daughter found for spatialRef[%d][%d]  "
                        "x=%g y=%g z=%g   numPart=%lu  numNodes=%d\n",
                        i, j,
                        spatialRef[i][j].centreDens.x,
                        spatialRef[i][j].centreDens.y,
                        spatialRef[i][j].centreDens.z,
                        spatialRef[i][j].numParts,
                        spatialRef[i][j].numNodes);
#endif
              
            }
          else if ( spatialRef[i][j].numSubStruct == 1 )
            {
              spatialRef[i][j].daughter.isoRefIndex = spatialRef[i][j].subStruct[0].isoRefIndex;
              spatialRef[i][j].daughter.refLevel    = i+1;
            }
          else
            {
              /* no daughter at all */
              spatialRef[i][j].daughter.isoRefIndex = -1;
              spatialRef[i][j].daughter.refLevel    = -1;
            }
          
          
          /**************************************************/
          /* Calculating the distance to the closest refinement */
          if ( spatialRef[i][j].numSubStruct > 1 )
            { /* IF */
              
              /* loop over all subhaloes within this host [i][j] */
              for ( k=0; k<spatialRef[i][j].numSubStruct; k++)
                { 
                  
                  isoRefIndex    = spatialRef[i][j].subStruct[k].isoRefIndex;
                  isoRefIndexOLD = isoRefIndex;
                  
                  x = spatialRef[i+1][isoRefIndex].centreDens.x;
                  y = spatialRef[i+1][isoRefIndex].centreDens.y;
                  z = spatialRef[i+1][isoRefIndex].centreDens.z;
                  
                  tmpMin = 10000000000000.0;
                  /* loop over all subhaloes other than k */
                  for ( l=0; l<spatialRef[i][j].numSubStruct; l++) 
                    {
                      
                      if ( k != l ) 
                        {
                          
                          isoRefIndex = spatialRef[i][j].subStruct[l].isoRefIndex;
                          
                          dx = x - spatialRef[i+1][isoRefIndex].centreDens.x;
                          dx = fabs(dx);	
                          
                          dy = y - spatialRef[i+1][isoRefIndex].centreDens.y;
                          dy = fabs(dy);	
                          
                          dz = z - spatialRef[i+1][isoRefIndex].centreDens.z;
                          dz = fabs(dz);	
                          
                          if (dx > 0.5 ) dx = 1.0 - dx;
                          if (dy > 0.5 ) dy = 1.0 - dy;
                          if (dz > 0.5 ) dz = 1.0 - dz;
                          
                          
                          dist = dx*dx + dy*dy + dz*dz;
                          
                          if ( dist < tmpMin )
                            tmpMin = dist;
                          
                        }
                      
                    }
                  
                  spatialRef[i+1][isoRefIndexOLD].closeRefDist = sqrt(tmpMin)/2.;
                  //spatialRef[i+1][isoRefIndexOLD].closeRefDist = sqrt(tmpMin/2.);
                }
            } /* IF */
          
        }
    }
  
  
  /*********************************************************************************************************************************************************/
  /*******************************************************************************************/
  /* Checking the refinement reconstructuion */
  totalcount = 0;
  for ( i=0; i<num_refgrids; i++) 
    {
      
      count1=0; count2=0; count3=0;
      
      /* Count the number of daughters */
      if ( i!= num_refgrids-1 )
        {
          for ( j=0; j<numIsoRef[i]; j++)
            {
              if ( spatialRef[i][j].daughter.isoRefIndex != -1)
                count1++;
            }
        }
      
      /* Count the number of substructure */
      if ( i!= num_refgrids-1 )
        {
          for ( j=0; j<numIsoRef[i]; j++) 
            count2 = count2 + spatialRef[i][j].numSubStruct;			
        }
      
      /* Count the number of parents */
      if ( i!= 0 ) 
        {
          for ( j=0; j<numIsoRef[i]; j++)
            {
              count3 = count3 + spatialRef[i][j].numParDom;
            }
        }
#ifdef VERBOSE			
      fprintf(stderr,"%3d || numSubStruct(%12d) numParDom(%12d) numDaughter(%12d) newHalos(%12d)\n",i,count2,count3,count1,count2-count1);
#endif
      totalcount = totalcount + count2-count1;
    }
#ifdef VERBOSE
  fprintf(stderr,"Number of Dark Matter Halos = %d\n",totalcount+numIsoRef[0]);
#endif
  
  
  /*******************************************************************************************/
  /* Checking for duplication of subhalos */
  
#ifdef VERBOSE
  for ( i=0; i<num_refgrids; i++)
    {
      
      tmpIsoRefIndex=NULL;
      tmpCount=0;
      
      for ( j=0; j<numIsoRef[i]; j++)
        {
          
          for ( k=0; k<spatialRef[i][j].numSubStruct; k++)
            {
              
              for ( l=0; l<tmpCount; l++) 
                {
                  if ( tmpIsoRefIndex[l] == spatialRef[i][j].subStruct[k].isoRefIndex )
                    fprintf(stderr,"#######  [%d][%d]  (while checking for duplication of halos...)\n",
                            spatialRef[i][j].subStruct[k].refLevel,
                            spatialRef[i][j].subStruct[k].isoRefIndex);
                }
              
              if(tmpIsoRefIndex == NULL)
                {
                  if ((tmpIsoRefIndex = calloc(tmpCount+1, sizeof(INDEX)))==NULL)
                    {
                      fprintf(stderr,"calloc failed in substructure %d %ld\n",tmpCount,sizeof(INDEX));
                      exit(-1);
                    }
                }
              else
                {
                  if ((tmpIsoRefIndex = realloc(tmpIsoRefIndex,(tmpCount+1)*sizeof(INDEX)))==NULL)
                    {
                      fprintf(stderr,"realloc failed in substructure %d %ld\n",tmpCount,sizeof(INDEX));
                      exit(-1);
                    }
                }
              tmpIsoRefIndex[tmpCount] = spatialRef[i][j].subStruct[k].isoRefIndex;
              tmpCount++;
            }
          
        }
      if (tmpIsoRefIndex != NULL)
        {
          free(tmpIsoRefIndex);
          tmpIsoRefIndex = NULL;
        }
    }
#endif /* VERBOSE*/
  
  
  /********************************************************************************************************************************************************/
  /********************************************************************************************************************************************************/
  /********************************************************************************************************************************************************/
  /* TESTING 
   
   for ( i=0; i<num_refgrids; i++) {
   for ( j=0; j<numIsoRef[i]; j++) {
   
   if ( i == 0 ) {
   
   fprintf(stderr,"HALO[%d][%d] ",i,j);
   fprintf(stderr,"CENTRE (%g,%g,%g) ", spatialRef[i][j].centreDens.x, spatialRef[i][j].centreDens.y, spatialRef[i][j].centreDens.z);
   fprintf(stderr,"DAUGH[%d][%d]\n",i+1,spatialRef[i][j].daughter.isoRefIndex);
   fprintf(stderr,"PARENT[%d][%d] ",-1, -1);
   fprintf(stderr,"NUM SUBSTRUCT = %d \n",spatialRef[i][j].numSubStruct);
   for ( k=0; k<spatialRef[i][j].numSubStruct; k++) {
   fprintf(stderr,"[%d] ",spatialRef[i][j].subStruct[k].isoRefIndex);
   }
   if ( spatialRef[i][j].numSubStruct != 0 )
   fprintf(stderr,"\n");
   
   } else if ( i == num_refgrids-1 ) {                            
   
   fprintf(stderr,"HALO[%d][%d] ",i,j);
   fprintf(stderr,"CENTRE (%g,%g,%g) ", spatialRef[i][j].centreDens.x, spatialRef[i][j].centreDens.y, spatialRef[i][j].centreDens.z);
   fprintf(stderr,"DAUGH[%d][%d]\n",-1,spatialRef[i][j].daughter.isoRefIndex);
   fprintf(stderr,"PARENT[%d][%d] ",i-1, spatialRef[i][j].parDom[0].isoRefIndex);
   
   } else {
   
   fprintf(stderr,"HALO[%d][%d] ",i,j);
   fprintf(stderr,"CENTRE (%g,%g,%g) ", spatialRef[i][j].centreDens.x, spatialRef[i][j].centreDens.y, spatialRef[i][j].centreDens.z);
   fprintf(stderr,"DAUGH[%d][%d]\n",i+1,spatialRef[i][j].daughter.isoRefIndex);
   fprintf(stderr,"PARENT[%d][%d] ",i-1, spatialRef[i][j].parDom[0].isoRefIndex);
   fprintf(stderr,"NUM SUBSTRUCT = %d \n",spatialRef[i][j].numSubStruct);
   for ( k=0; k<spatialRef[i][j].numSubStruct; k++) {
   fprintf(stderr,"[%d] ",spatialRef[i][j].subStruct[k].isoRefIndex);
   }
   if ( spatialRef[i][j].numSubStruct != 0 )
   fprintf(stderr,"\n");
   
   }
   }
   }
   
   ********************************************************************************************************************************************************/
  /********************************************************************************************************************************************************/
  /********************************************************************************************************************************************************/
  
  return(TRUE);
}



int	spatialRef2halos(int num_refgrids, SPATIALREF **spatialRef)
{
  int      i,j,k,f;
  int      count,numNewHalos;
  int	   isoRefIndex, refLevel;
  int	   SSisoRefIndex, SSrefLevel;
  double   dx,dy,dz,tmpRad;
  int	   primHaloIndex, haloIndex;
  partptr	current, previous, subCurrent;
  int	   countBC,OKcountBC,tmpCount;
  int      countTMP;
  
  int      kcount;
  int     *tmpSubStruct;
  
  int	   deepRef;
  double   oldRad;
  int      jnumpart,inumpart;
  double   maxGathRad, gatherRad2, Rvir2;
  
  int      hostHaloGrids, tmp;
  
  /* re-calculate the number of potential halos */
  count = 0;
  for ( i=0; i<num_refgrids; i++) 
    for ( j=0; j<numIsoRef[i]; j++) 
      if ( spatialRef[i][j].numSubStruct == 0 )
        count++;   
  numHalos = count;
  
#ifdef VERBOSE
  fprintf(stderr," spatialRef2halos():\n");
  fprintf(stderr,"  number of isolated refinements  = %d\n",totnumIsoRef);
  fprintf(stderr,"  first guess for number of halos = %d\n",numHalos);
  fprintf(io.logfile," spatialRef2halos():\n");
  fprintf(io.logfile,"  number of isolated refinements  = %d\n",totnumIsoRef);
  fprintf(io.logfile,"  first guess for number of halos = %d\n",numHalos);
  fflush(io.logfile);
#endif
  
  
  /* allocate memory for halos: 
   *----------------------------
   * rather use totnumIsoRef than numHalos 
   * because there is still this strange bug leading to "count != numHalos"
   */
  halos=NULL;
  if ((halos = calloc(totnumIsoRef+1,sizeof(HALO)))==NULL)
    {
      fprintf(io.logfile,"Error in allocating the memory for halo array\n");
      exit(0);
    }
  
  /* initialize halo properties */
  for ( i=0; i<totnumIsoRef+1; i++)
    {
      halos[i].npart   = 0; 
      halos[i].nll     = 0;
      halos[i].ipart   = NULL;
      halos[i].ll      = NULL;
      
      halos[i].pos.x    = 0;
      halos[i].pos.y    = 0;
      halos[i].pos.z    = 0;
      halos[i].vel.x    = 0;
      halos[i].vel.y    = 0;
      halos[i].vel.z    = 0;
      halos[i].M_vir    = -1.0;
      halos[i].R_vir    = -1.0;
      halos[i].velDis   = 0.0;
      halos[i].lambda   = 0.0;
      halos[i].lambdaE  = 0.0;
      halos[i].R_max    = 0.0;
      halos[i].V2_max   = 0.0;
      halos[i].ovdens   = 0.0;
      halos[i].R_edge   = 0.0;
      halos[i].Ekin     = 0.0;
      halos[i].Epot     = 0.0;
      halos[i].Phi0     = 0.0;
      halos[i].axis.x   = 0.0;
      halos[i].axis.y   = 0.0;
      halos[i].axis.z   = 0.0;
      halos[i].E1.x     = 0.0;
      halos[i].E1.y     = 0.0;
      halos[i].E1.z     = 0.0;
      halos[i].E2.x     = 0.0;
      halos[i].E2.y     = 0.0;
      halos[i].E2.z     = 0.0;
      halos[i].E3.x     = 0.0;
      halos[i].E3.y     = 0.0;
      halos[i].E3.z     = 0.0;      
      halos[i].AngMom.x = 0.0;
      halos[i].AngMom.y = 0.0;
      halos[i].AngMom.z = 0.0;
      halos[i].com_offset = 0.0;
      halos[i].mbp_offset = 0.0;
      halos[i].r2         = 0.0;
      
      
      halos[i].hostHaloGrids    = -1;
      halos[i].hostHalo         = -1;
      
      halos[i].numHostHaloGrids = 0;
      halos[i].hostGrids        = NULL;
      
      halos[i].gatherRad        = 100000000000.0;
      
      halos[i].numSubStruct     =  0;
      halos[i].subStruct        = NULL;
      
      halos[i].spaRes           = 0.0;
      halos[i].refLev           = 0;
      
      halos[i].numNodes = 0;
    }
  
  /**********************************************************************************************************/ 
  /**********************************************************************************************************/ 
  /* connecting the halos */
  /* NOTE :: IF 'halos[i].hostHalo = i' then it is not the substrcture of any other halo */
  count=0; tmpCount=0; countBC=0;
#ifdef VERBOSE
  fprintf(stderr,"  constructing %d (potential) halos from all %d grid levels...\n",numHalos,num_refgrids);
  fprintf(io.logfile,"  constructing %d (potential) halos from all %d grid levels...\n",numHalos,num_refgrids);
  fflush(io.logfile);
#endif
  for ( i=0; i<num_refgrids; i++)
    {
#ifdef VERBOSE
      fprintf(stderr,"      grid level %8d -> %10d isolated refinements\n",i,numIsoRef[i]);
      fprintf(io.logfile,"      grid level %8d -> %10d isolated refinements\n",i,numIsoRef[i]);
      fflush(io.logfile);
#endif
      
      /* treat domain grid separately */
      if ( i == 0 )
        {	
          
          for ( j=0; j<numIsoRef[i]; j++)
            {
              /******************************************************/
              if ( spatialRef[i][j].numSubStruct == 0 )
                {   
                  /* halo centre (there is no finer spatialRef[][] and hence assign centre...) */
                  halos[count].pos.x = spatialRef[i][j].centre.x;
                  halos[count].pos.y = spatialRef[i][j].centre.y;
                  halos[count].pos.z = spatialRef[i][j].centre.z;
                  
                  /* The spatial resolution of this halo */
                  halos[count].spaRes = 1.0/((double)gridl1dim[i]);
                  halos[count].refLev = i;
                  
                  /* particles */
                  halos[count].npart  = spatialRef[i][j].numParts;
                  halos[count].ll     = spatialRef[i][j].ll;
                  spatialRef[i][j].ll = NULL;
                  
                  /* What is my host halo ? */
                  halos[count].hostHalo = -1;
                  
                  /*  Substructure */
                  halos[count].numSubStruct =  0;
                  halos[count].subStruct    = NULL;
                  
                  /* nodes */
                  halos[count].numNodes = spatialRef[i][j].numNodes;
                  
                  count++;
                  tmpCount++; /* Closing of the halo */
                  
                  /******************************************************/
                } 
              else if ( spatialRef[i][j].numSubStruct == 1 )
                {
                  /* do not assign centre as there is a finer refinement... */
                  
                  /* particles */
                  halos[count].npart  = spatialRef[i][j].numParts;
                  halos[count].ll     = spatialRef[i][j].ll;
                  spatialRef[i][j].ll = NULL;
                  
                  /* nodes */
                  halos[count].numNodes = spatialRef[i][j].numNodes;
                  
                  /* What is my host halo ? */
                  halos[count].hostHalo = -1;
                  
                  /*  Substructure (the one substructure on spatialRef[][] is in fact the host!) */
                  halos[count].numSubStruct =  0;
                  halos[count].subStruct    = NULL;
                  
                  /* Telling daughter about host halo */	
                  refLevel    = spatialRef[i][j].daughter.refLevel;
                  isoRefIndex = spatialRef[i][j].daughter.isoRefIndex;
                  
                  if(refLevel != -1 && isoRefIndex != -1)
                    spatialRef[refLevel][isoRefIndex].haloIndex = count;
#ifdef AHFverbose
                  else
                    fprintf(stderr,"      spatialRef2halos():  WARNING -> no daughter for spatialRef[%d][%d]\n",i,j);
#endif
                  
                  count++;
                  
                  /******************************************************/
                }
              else
                {
                  /* What is my host halo ? */
                  halos[count].hostHalo = -1; /* It is a base halo and throw it to the wind */
                  
                  /*count;  If the host halo name is the same as the iterator then it is a base halo */
                  primHaloIndex = count;
                  
                  /* Telling daughter about host halo */	
                  refLevel    = spatialRef[i][j].daughter.refLevel;
                  isoRefIndex = spatialRef[i][j].daughter.isoRefIndex;
                  
                  if(refLevel != -1 && isoRefIndex != -1)
                    spatialRef[refLevel][isoRefIndex].haloIndex = count;
#ifdef AHFverbose
                  else
                    fprintf(stderr,"      spatialRef2halos():  WARNING -> no daughter for spatialRef[%d][%d]\n",i,j);
#endif
                  
                  /* nodes */
                  halos[count].numNodes = spatialRef[i][j].numNodes;
                  
                  /*  Substructure */
                  numNewHalos                 = spatialRef[i][j].numSubStruct - 1;
                  halos[count].numSubStruct   = numNewHalos;
                  halos[count].subStruct      = NULL;
                  if ((halos[count].subStruct = calloc(numNewHalos,sizeof(int)))==NULL)
                    {
                      fprintf(stderr,"Error in allocating the memory for halos[count].subStruct array\n");
                      exit(0);
                    }
                  
                  /* increment halo counter as we now... */
                  count++;
                  
                  /* ...generate the substructure halos of "primHaloIndex" (= the host) */
                  kcount=0;
                  for (k=0;k<numNewHalos+1;k++)
                    {
                      
                      SSrefLevel    = spatialRef[i][j].subStruct[k].refLevel;
                      SSisoRefIndex = spatialRef[i][j].subStruct[k].isoRefIndex;
                      
                      /* Make sure not the daughter */
                      if ( SSisoRefIndex != isoRefIndex ) 
                        {
                          /* What is my host halo ? */
                          halos[count].hostHalo = primHaloIndex;
                          
                          /* Telling the subHalos who they are */
                          spatialRef[SSrefLevel][SSisoRefIndex].haloIndex = count;
                          
                          /* Adding that name to the list of substrcture */
                          halos[primHaloIndex].subStruct[kcount] = count;
                          
                          /* initial guess for virial radius */  
                          halos[count].R_vir = spatialRef[SSrefLevel][SSisoRefIndex].closeRefDist;
                          Rvir2              = pow2(halos[count].R_vir);
                          
                          /* subhalo has no particles yet */
                          halos[count].ll = NULL;
                          
                          count++;   // move to next (sub-)halo
                          kcount++;  // kcount should add up to (spatialRef[i][i].numSubStruct-1)
                        }
                    }
                  
                  /* Assign the remaining particles to the host halo */
                  halos[primHaloIndex].npart = spatialRef[i][j].numParts;
                  halos[primHaloIndex].ll    = spatialRef[i][j].ll;
                  spatialRef[i][j].ll        = NULL;
                  
                }
            } /* for(j<numIsoRef[i]) */
        }
      
      /* i !=0 */
      else
        {
          for ( j=0; j<numIsoRef[i]; j++)
            {
              /* Caution for the boundary condition problem :- now not a problem */
              if ( spatialRef[i][j].numParDom != 0)
                {
                  /******************************************************/
                  if ( spatialRef[i][j].numSubStruct == 0 )
                    {
                      /* haloIndex is the ID of the halo correspondong to the parent refinement level */
                      haloIndex = spatialRef[i][j].haloIndex;
                      
                      if ( haloIndex == -1 )
                        {
                          fprintf(stderr,"  HELP haloIndex HELP :: HALO[%d][%d] \n",i,j);
                          fprintf(stderr,"  HALO-dau[%d][%d] \n", spatialRef[i][j].daughter.refLevel, spatialRef[i][j].daughter.isoRefIndex);
                          fprintf(stderr,"  HALO-par[%d][%d] \n", spatialRef[i][j].parDom[0].isoRefIndex, spatialRef[i][j].parDom[0].refLevel);
                        }
                      
                      /* The spatial resolution of this halo :- to finish ALEX */
                      halos[haloIndex].spaRes = 1.0/((double)gridl1dim[i]);
                      halos[haloIndex].refLev = i;
                      
                      /* halo centre -> change to new centre of this finer refinement! */
                      halos[haloIndex].pos.x = spatialRef[i][j].centre.x;
                      halos[haloIndex].pos.y = spatialRef[i][j].centre.y;
                      halos[haloIndex].pos.z = spatialRef[i][j].centre.z;
                      
                      /* nodes       -> record number of nodes centre is based upon */
                      halos[haloIndex].numNodes = spatialRef[i][j].numNodes;
                      
                      /* particles   -> add additional particles */
                      halos[haloIndex].npart   += spatialRef[i][j].numParts;
                      
                      /* Loop to the last particle in the linked-list */	
                      current = halos[haloIndex].ll;
                      previous = NULL;
                      while ( current!=NULL )
                        {
                          previous = current;
                          current	= current->ll;
                        }
                      
                      /* point last particle to the head of this ref's ll */
                      if ( previous == NULL )
                        halos[haloIndex].ll = spatialRef[i][j].ll;
                      else
                        previous->ll = spatialRef[i][j].ll;
                      
                      spatialRef[i][j].ll = NULL;
                      
                      tmpCount++; /* Closing of the halo */
                      
                      
                      /******************************************************/
                      
                    }
                  else if ( spatialRef[i][j].numSubStruct == 1 )
                    {
                      /* haloIndex is the ID of the halo correspondong to the parent refinement level */
                      haloIndex = spatialRef[i][j].haloIndex;
                      
                      if ( haloIndex == -1 )
                        fprintf(stderr,"  HELP2 haloIndex HELP2\n");
                      
                      /* particles */
                      halos[haloIndex].npart   += spatialRef[i][j].numParts;
                      
                      /* nodes */
                      halos[haloIndex].numNodes = spatialRef[i][j].numNodes;
                      
                      /* Loop to the last particle in the partile ll */	
                      current = halos[haloIndex].ll;
                      previous = NULL;
                      while ( current!=NULL )
                        {
                          previous = current;
                          current	= current->ll;
                        }
                      
                      /* point last particle to the head of this ref's ll */
                      if ( previous == NULL )
                        halos[haloIndex].ll = spatialRef[i][j].ll;
                      else
                        previous->ll = spatialRef[i][j].ll;
                      
                      spatialRef[i][j].ll = NULL;
                      
                      /* Telling daughter about host halo */	
                      refLevel    = spatialRef[i][j].daughter.refLevel;
                      isoRefIndex = spatialRef[i][j].daughter.isoRefIndex;
                      
                      if(refLevel != -1 && isoRefIndex != -1)
                        spatialRef[refLevel][isoRefIndex].haloIndex = haloIndex;
#ifdef AHFverbose
                      else
                        fprintf(stderr,"      spatialRef2halos():  WARNING -> wrong daughter for spatialRef[%d][%d]\n",i,j);
#endif
                      
                      /******************************************************/
                      
                    } 
                  else
                    {
                      /* haloIndex is the ID of the halo correspondong to the parent refinement level */
                      haloIndex     = spatialRef[i][j].haloIndex;
                      primHaloIndex = haloIndex;
                      
                      if ( haloIndex == -1 )
                        fprintf(stderr,"  HELP3 haloIndex HELP3\n");
                      
                      /* nodes */
                      halos[haloIndex].numNodes = spatialRef[i][j].numNodes;
                      
                      /* Telling daughter about host halo */	
                      refLevel    = spatialRef[i][j].daughter.refLevel;
                      isoRefIndex = spatialRef[i][j].daughter.isoRefIndex;
                      
                      if(refLevel != -1 && isoRefIndex != -1)
                        spatialRef[refLevel][isoRefIndex].haloIndex = haloIndex;
#ifdef AHFverbose
                      else
                        fprintf(stderr,"      spatialRef2halos():  WARNING -> wrong daughter for spatialRef[%d][%d]\n",i,j);
#endif
                      
                      
                      /*  Substructure */
                      numNewHalos  = spatialRef[i][j].numSubStruct - 1;
                      kcount       = halos[primHaloIndex].numSubStruct;
                      tmpSubStruct = NULL;
                      if ((tmpSubStruct = calloc(kcount+1,sizeof(int)))==NULL)
                        {
                          fprintf(stderr,"Error in allocating the memory for tmpSubStruct array\n");
                          exit(0);
                        }
                      for ( k=0; k<kcount; k++ )
                        tmpSubStruct[k] = halos[primHaloIndex].subStruct[k];
                      
                      halos[primHaloIndex].numSubStruct = kcount + numNewHalos;
                      
                      if(halos[primHaloIndex].subStruct == NULL)
                        {
                          if ((halos[primHaloIndex].subStruct = calloc(halos[primHaloIndex].numSubStruct+1,sizeof(int)))==NULL)
                            {
                              fprintf(stderr,"Error in allocating the memory for halos[count].subStruct array\n");
                              exit(0);
                            }
                        }
                      else
                        {
                          if ((halos[primHaloIndex].subStruct = realloc(halos[primHaloIndex].subStruct,(halos[primHaloIndex].numSubStruct+1)*sizeof(int)))==NULL)
                            {
                              fprintf(stderr,"Error in reallocating the memory for halos[count].subStruct array\n");
                              exit(0);
                            }
                        }
                      for ( k=0; k<kcount; k++ )
                        halos[primHaloIndex].subStruct[k] = tmpSubStruct[k];
                      
                      free(tmpSubStruct);
                      
                      
                      
                      /* Collecting information for the 'new' substructure */
                      for (k=0;k<spatialRef[i][j].numSubStruct;k++)
                        {
                          
                          SSrefLevel    = spatialRef[i][j].subStruct[k].refLevel;
                          SSisoRefIndex = spatialRef[i][j].subStruct[k].isoRefIndex;
                          
                          if ( SSisoRefIndex!=isoRefIndex )
                            {
                              /* What is my host halo ? */
                              halos[count].hostHalo = primHaloIndex;
                              
                              /* Telling the subHalos who they are */
                              spatialRef[SSrefLevel][SSisoRefIndex].haloIndex = count;
                              
                              /*  Substructure */
                              halos[count].subStruct = NULL;
                              
                              /*  Substructure of primary index */
                              halos[primHaloIndex].subStruct[kcount] = count;
                              
                              /* first guess for subhalo's virial radius */
                              halos[count].R_vir = spatialRef[SSrefLevel][SSisoRefIndex].closeRefDist;
                              Rvir2              = pow2(halos[count].R_vir);
                              
                              /* subhalo has no particles yet */
                              halos[count].ll = NULL;                        
                              
                              count++;
                              kcount++;
                            }
                        }
                      
                      
                      /* Assign the remaining particles to the host halo */
                      halos[haloIndex].npart += spatialRef[i][j].numParts;
                      
                      /* Loop to the last particle in the partile ll */	
                      current = halos[haloIndex].ll;
                      previous = NULL;
                      while ( current!=NULL )
                        {
                          previous = current;
                          current	= current->ll;
                        }
                      
                      /* point last particle to the head of this ref's ll */
                      if ( previous == NULL )
                        halos[haloIndex].ll = spatialRef[i][j].ll;
                      else
                        previous->ll = spatialRef[i][j].ll;
                      
                      spatialRef[i][j].ll = NULL;
                    }
                }
              else 
                {                 
                  OKcountBC=1;
                  for( f=0; f<numDensZero; f++ )
                    {
                      if ( (i==densZero[f].refLevel) && (j==densZero[f].isoRefIndex))
                        OKcountBC=0;
                    }
#ifdef VERBOSE
                  fprintf(stderr,"      s %g %g %g 0.08 0.0 1.0 1.0\n",spatialRef[i][j].centreGEOM.x*simu.boxsize,spatialRef[i][j].centreGEOM.y*simu.boxsize,spatialRef[i][j].centreGEOM.z*simu.boxsize);
                  fprintf(stderr,"      countBC[%d][%d]\n",i,j);
                  fprintf(stderr,"      OKcountBC = %d (0 is good)\n",OKcountBC);
#endif
                  
                  countBC++;
                }
            } /* for(j<numIsoRef[i]) */
        } /* i !=0 */
      
    } /* for(i<num_refgrids) */
  
#ifdef AHFDEBUG
  fprintf(stderr,"done\n");
#endif
  
#ifdef VERBOSE
  if ( countBC > 0 )
    fprintf(stderr,"      countBC = %d => isolated refinements do not have a parent refinement\n",countBC);
  
  fprintf(stderr,"      number of halos: found=%d (max=%d) expected=%d (tmpCount=%d)\n",
          count, totnumIsoRef, numHalos, tmpCount);
  fprintf(stderr,"      numDensZero = %d\n", numDensZero);
  
  fprintf(io.logfile,"      number of halos: found=%d (max=%d) expected=%d (tmpCount=%d)\n",
          count, totnumIsoRef, numHalos, tmpCount);
  fprintf(io.logfile,"      numDensZero = %d\n", numDensZero);      
  fflush(io.logfile);
  
  /* did we find more halos than initially expected? */
  if(count > numHalos) 
    {
      numHalos = count;
      fprintf(stderr,"      => adjusted number of halos: found=%d (max=%d) expected=%d (tmpCount=%d)\n",
              count, totnumIsoRef, numHalos, tmpCount);
      fprintf(io.logfile,"      => adjusted number of halos: found=%d (max=%d) expected=%d (tmpCount=%d)\n",
              count, totnumIsoRef, numHalos, tmpCount);
      fflush(io.logfile);
    }
  
#endif
  
  /**********************************************************************************************************/ 
  /**********************************************************************************************************/ 
  /* Ordering the halos wrt mass 																									*/
  /* The assumption is that sub-halos always have less mass than thier host 			*/
  /* Thus we can move through the list giving particles to the host 
   * without having to re-order the internal particles all the time 							*/ 
  
  /**********************************************************************************************/
  /* Remember old index */
  for ( i=0; i<numHalos; i++ ) 
    halos[i].oldIndex = i;
  
  /**********************************************************************************************/
  /* Actually re-ordering the halos wrt mass */
  qsort(halos,numHalos,sizeof(HALO),haloCompare);
  
  
  /**********************************************************************************************/
  /* Need to use 'old index' to re-order the host halos and the sub-structure listings 
   * I.e. everything that reqires an index internal to the host halos
   * I should have used linked lists!!! */ 
  for ( i=0; i<numHalos; i++ )
    {
      
      /* Re-directing the host halo index */
      
      /* Search through the halos - find the oldIndex that matches hostHalo
       * When found use the search index as the hostHaloGrids index */
      for ( j=0; j<numHalos; j++ )
        {
          if ( halos[i].hostHalo == halos[j].oldIndex )
            {
              halos[i].hostHaloGrids = j;
              halos[i].hostHalo      = -2;
              
              /* Check that hostHaloGrids is not larger that the halo 
               * I.e. the sub-halo is not larger than the host */ 
              if ( i < j )
                {
                  halos[i].hostHaloGrids = -1;
                  /* fprintf(stderr,"halos[%d].hostHaloGrids = %d (i<j)\n",i,j); */
                } 
              
              break;
            }
        }
      
      /* Re-directing the sub-halo index */
      for ( k=0; k<halos[i].numSubStruct; k++ ) 
        {
          
          for ( j=0; j<numHalos; j++ ) 
            {
              if ( halos[i].subStruct[k] == halos[j].oldIndex ) 
                {
                  halos[i].subStruct[k] = j;
                  break;
                }
            }	
        }
    }	/* for ( i=0; i<numHalos; i++ ) */
  
  
  /**********************************************************************************************************/ 
  /**********************************************************************************************************/ 
  /* Calculating the embedded GRID host halos
   * This will catch when you have a premature bridge b/w halos and the grid host halo is prematurely defined */
  /* Loop through the host halos - for each halo start building the 'hostGrids' array */
  for ( i=0; i<numHalos; i++ ) 
    {
      
      count                     = 0;
      halos[i].hostGrids        = NULL;
      halos[i].numHostHaloGrids = 0;
      
      /* This is the halos primary Grid host halo */
      hostHaloGrids = halos[i].hostHaloGrids;
      
      if ( hostHaloGrids != -1 ) 
        {
          
          /* Look at the hostHaloGrids of 'i''s hostHaloGrids */
          tmp = halos[hostHaloGrids].hostHaloGrids;
          
          while ( tmp >= 0 ) 
            {
              
              count++;
              
              
              if(halos[i].hostGrids == NULL)
                {
                  if ((halos[i].hostGrids = calloc(count+1, sizeof(int)))==NULL) 
                    {
                      fprintf(stderr,"calloc failed for hostGrids\n");
                      exit(-1);
                    }
                }
              else
                {
                  if ((halos[i].hostGrids = realloc(halos[i].hostGrids,(count+1)*sizeof(int)))==NULL) 
                    {
                      fprintf(stderr,"realloc failed for hostGrids\n");
                      exit(-1);
                    }
                }
              
              halos[i].hostGrids[count-1] = tmp;
              
              /* What is the new 	hostHaloGrids */
              tmp = halos[tmp].hostHaloGrids;
              
            }
          
          halos[i].numHostHaloGrids = count;
          
        }	/* if ( hostHaloGrids != -1 ) */
    }	/* for ( i=0; i<numHalos; i++ ) */
  
  
  /**********************************************************************************************************/ 
  /**********************************************************************************************************/ 
  /* Calculate the Gathering radius:                                                                        */
  /* the gathering radius is the distance to the closest halo that is more massive than the current halo... */
  
  /* cf. AHF_MAX_GATHER_RAD definition */
  maxGathRad = MIN(AHF_MAX_GATHER_RAD/simu.boxsize, 1./4.);
  
#ifdef WITH_OPENMP
  /* TODO: this loop may be parallelized! */
#endif
  for ( i=numHalos-1; i>=0; i-- )
    {    
      /* root halos do not need to gather anything from larger objects... */
      //if ( halos[i].hostHaloGrids != -1 ) 
        {
          count    = 0;
          
          inumpart = halos[i].npart;
          
          /* loop over all halos that are more massive... */
          for ( j=numHalos-1; j>=0; j-- )
            {
              
              jnumpart = halos[j].npart;
              
              if ( (i!=j) && ( jnumpart > inumpart) )
                {
                  
                  /* Calculate the distance to this more massive halo */
                  dx = fabs(halos[i].pos.x - halos[j].pos.x);
                  dy = fabs(halos[i].pos.y - halos[j].pos.y);
                  dz = fabs(halos[i].pos.z - halos[j].pos.z);
                  if ( dx > 0.5 ) dx = 1.0 - dx;
                  if ( dy > 0.5 ) dy = 1.0 - dy;
                  if ( dz > 0.5 ) dz = 1.0 - dz;
                  
                  tmpRad = dx*dx + dy*dy + dz*dz;
                  
                  /* Is it the smallest distance */
                  if ( tmpRad < halos[i].gatherRad )
                    {
                      halos[i].gatherRad = tmpRad;
                    }
                  
                  count++;
                }
            }
          
          /* Finally calculating the gathering radius */
          halos[i].gatherRad = (sqrt(halos[i].gatherRad))*.5;
          
          if ( count == 0 ) 
            {
#ifdef AHFDEBUG
              fprintf(stderr,"There are no other halos bigger than this halo - should see this just once\n");
#endif
              halos[i].gatherRad = maxGathRad;
            } 
          
        } 
//      else 
//        { /* if ( halos[i].hostHaloGrids != -1 ) */
//          halos[i].gatherRad = maxGathRad;				
//        }
      
      /* restrict the gathering radius */
      halos[i].gatherRad = MAX(halos[i].gatherRad, halos[i].R_vir);
      halos[i].gatherRad = MIN(halos[i].gatherRad, maxGathRad);
      
    } /* for ( i=numHalos-1; i>=0; i-- ) */
  
  
  
  /**********************************************************************************************************/ 
  /**********************************************************************************************************/ 
  /* Calculate the new Host Halo based on the gather radius 
   * Also calculate other hosts */
  
  for ( i=0; i<numHalos; i++ )
    {
      
      /* Now - run through the halos to see if any of them are sub halos of this one :) */
      /* Note:  By definition only halos that are smaller in mass than the host are considered sub-halos */
      gatherRad2 = pow2(halos[i].gatherRad);
      count      = 0;
      for ( j=i; j<numHalos; j++ )
        {
          
          /* Clearly don't need to check root halos b/c they are below the virial density */
          if ( (halos[j].hostHaloGrids != -1) && (i!=j) )
            { /* Halo is not root halos && not this halo */
              
              /* Is this halo 'j' a sub-halo of 'i' */
              /* By sub-halo we mean is it within it's gathering halo */
              dx = halos[i].pos.x - halos[j].pos.x;
              dx = fabs(dx);
              dy = halos[i].pos.y - halos[j].pos.y;
              dy = fabs(dy);
              dz = halos[i].pos.z - halos[j].pos.z;
              dz = fabs(dz);
              if ( dx > 0.5 ) dx = 1.0 - dx;
              if ( dy > 0.5 ) dy = 1.0 - dy;
              if ( dz > 0.5 ) dz = 1.0 - dz;
              
              /* sqdist between the potential host and the potential sub-halo	*/
              tmpRad = dx*dx + dy*dy + dz*dz;
              if ( tmpRad < gatherRad2 )
                {
                  /* fprintf(stderr,"tmpRad[%d](%g) < halos[%d].gatherRad(%g)\n",j, tmpRad, i, halos[i].gatherRad*halos[i].gatherRad ); */
                  
                  /* If it has got this far then 'j' is a sub-halo of 'i' */
                  /* fprintf(stderr,"%d is a sub halo of %d\n",j,i); */
                  /*	However, has this halo been identified as being a sub-halo to another host? */
                  if ( halos[j].hostHalo != -2 )
                    {
                      
                      count++;
                      /* The true host halo of this sub-halo is the closer halo */
                      /* Thus: calculating the distance beteen the halo 'j' and its pontentially 'old' host */
                      dx = halos[halos[j].hostHalo].pos.x - halos[j].pos.x;
                      dx = fabs(dx);
                      dy = halos[halos[j].hostHalo].pos.y - halos[j].pos.y;
                      dy = fabs(dy);
                      dz = halos[halos[j].hostHalo].pos.z - halos[j].pos.z;
                      dz = fabs(dz);
                      if ( dx > 0.5 ) dx = 1.0 - dx;
                      if ( dy > 0.5 ) dy = 1.0 - dy;
                      if ( dz > 0.5 ) dz = 1.0 - dz;
                      
                      oldRad = dx*dx + dy*dy + dz*dz;
                      
                      /* Checking the distance of the two potential hosts */
                      /* Set the host to the halo that is the closest */
                      if ( tmpRad < oldRad ) 
                        halos[j].hostHalo = i;
                      
                    } 
                  else
                    { /* This halo has no other potential host */
                      /* 	mark our halo 'j' as a sub-halo of 'i' */
                      halos[j].hostHalo = i;
                    }
                  
                  /* Checking if the sub-halo has more mass than the host */
                  if ( halos[i].npart < halos[j].npart )
                    {
                      fprintf(stderr,"WARNING! :: The sub-halo has more mass than the host -- investigate\n");
                      fprintf(stderr,"halos[%d].npart[%ld] < halos[%d].npart[%ld]\n",i,halos[i].npart,j,halos[j].npart);
                    }
                  
                } /* if ( tmpRad < halos[i].gatherRad*halos[i].gatherRad ) */
            } /* if ( halos[j].hostHalo != -1 ) */
        } /* for ( j=numHalos-1; j>i; j-- ) */
    } /* for ( i=numHalos-1; i>=0; i-- ) */
  
  
  /* remove excess particlese from all host halos */
  //prune_host_halos();
  
#ifdef AHFverbose
  /**********************************************************************************************************/ 
  /**********************************************************************************************************/ 
  /* Dump preliminary halo properties to logfile */
  fprintf(io.logfile,"\n\n   preliminary halo properties:\n");
  fprintf(io.logfile,"   ============================\n");
  for ( i=numHalos-1; i>=0; i-- )
    {
      fprintf(io.logfile,"   halos[%d].numNodes         = %16d\n",  i,halos[i].numNodes);
      fprintf(io.logfile,"   halos[%d].npart            = %16ld\n", i,halos[i].npart);
      fprintf(io.logfile,"   halos[%d].R_vir            = %16.8g (=closeRefDist)\n", i,halos[i].R_vir*x_fac*1000.);
      fprintf(io.logfile,"   halos[%d].gatherRad        = %16.8g (=closeHostDist)\n",i,halos[i].gatherRad*x_fac*1000.);
      fprintf(io.logfile,"   halos[%d].spaRes           = %16.8g\n",i,halos[i].spaRes);
      fprintf(io.logfile,"   halos[%d].pos.x            = %16.8g\n",i,halos[i].pos.x*x_fac);
      fprintf(io.logfile,"   halos[%d].pos.y            = %16.8g\n",i,halos[i].pos.y*x_fac);
      fprintf(io.logfile,"   halos[%d].pos.z            = %16.8g\n",i,halos[i].pos.z*x_fac);
      fprintf(io.logfile,"   halos[%d].hostHalo         = %16d\n",  i,halos[i].hostHalo);
      fprintf(io.logfile,"   halos[%d].hostHaloGrids    = %16d\n",  i,halos[i].hostHaloGrids);
      fprintf(io.logfile,"   halos[%d].numHostHaloGrids = %16d\n",  i,halos[i].numHostHaloGrids);
      fprintf(io.logfile,"            gatherRad/R_vir   = %16.8g\n",  halos[i].gatherRad/halos[i].R_vir);
      fprintf(io.logfile,"   #########################################################################\n");
    }
  fprintf(io.logfile,"\n\n");
#endif /*AHFverbose */
  
  return(TRUE);
}

/*
 ************************************************************
 ************************************************************
 * Orders the particels within the halos
 */

int compare( struct particle *p, struct particle *q, XYZ *pointer) {
  
  double x1,y1,z1,dist1;
  double x2,y2,z2,dist2;
  double xx,yy,zz;
  double dx1,dy1,dz1;
  double dx2,dy2,dz2;
  
  xx = pointer->x; yy = pointer->y; zz = pointer->z;
  
  x1 = p->pos[0]; y1 = p->pos[1]; z1 = p->pos[2];
  
  dx1 = xx-x1;
  dx1 = fabs(dx1);	
  if (dx1 > 0.5 )
    dx1 = 1.0 - dx1;
  
  dy1 = yy-y1;
  dy1 = fabs(dy1);	
  if (dy1 > 0.5 )
    dy1 = 1.0 - dy1;
  
  dz1 = zz-z1;
  dz1 = fabs(dz1);	
  if (dz1 > 0.5 )
    dz1 = 1.0 - dz1;
  
  x2 = q->pos[0]; y2 = q->pos[1]; z2 = q->pos[2];	
  
  dx2 = xx - x2;
  dx2 = fabs(dx2);	
  if (dx2 > 0.5 )
    dx2 = 1.0 - dx2;
  
  dy2 = yy-y2;
  dy2 = fabs(dy2);	
  if (dy2 > 0.5 )
    dy2 = 1.0 - dy2;
  
  dz2 = zz-z2;
  dz2 = fabs(dz2);	
  if (dz2 > 0.5 )
    dz2 = 1.0 - dz2;
  
  
  dist1 = ( dx1*dx1 + dy1*dy1 + dz1*dz1 )*simu.boxsize*simu.boxsize;
  dist2 = ( dx2*dx2 + dy2*dy2 + dz2*dz2 )*simu.boxsize*simu.boxsize;
  
  
  if ( dist1 < dist2 )
    return(-1);
  else if ( dist1 > dist2 )
    return(1);
  else
    return(0);
  
  
  
}


/****************************************************************************************
 * Prune host halos
 *
 * remove particles that do not lie within gatherRad from host halos
 ****************************************************************************************/ 
void prune_host_halos(void)
{
  long    ihalo, npart;
  partptr cur_part, prev_part;
  double  dx, dy, dz;
  double  Dist2, gatherRad2;
  
  /* loop over all halos */
  for ( ihalo=0; ihalo<numHalos; ihalo++ )
    {
      /* criterion to be considered a host halo */
      if(halos[ihalo].hostHalo         == -1 &&
         halos[ihalo].hostHaloGrids    == -1 &&
         halos[ihalo].numHostHaloGrids == 0)
        {
          /* set the gathering radius */
          gatherRad2 = pow2(halos[ihalo].gatherRad);
          
          /* check whether we called from a) or b) */
          if(halos[ihalo].npart == 0)
            {
#ifdef AHFverbose
              fprintf(stderr,"pruning host halo %12ld nll  =%12d  ->  ",ihalo,halos[ihalo].nll);
              fprintf(io.logfile,"pruning host halo %12ld nll  =%12d  ->  ",ihalo,halos[ihalo].nll);
#endif
              npart = halos[ihalo].nll;
            }
          else
            {
#ifdef AHFverbose
              fprintf(stderr,"pruning host halo %12ld npart=%12ld ->  ",ihalo,halos[ihalo].npart);
              fprintf(io.logfile,"pruning host halo %12ld npart=%12ld ->  ",ihalo,halos[ihalo].npart);
#endif
              npart = halos[ihalo].npart;
            }
          
          /* loop over all particles */
          prev_part = NULL;
          cur_part  = halos[ihalo].ll;
          
          while ( cur_part != NULL )
            {
              /* Calculating the distance between the particle in the host halo and the centre of the current halo */
              dx = halos[ihalo].pos.x - cur_part->pos[X];
              dy = halos[ihalo].pos.y - cur_part->pos[Y];
              dz = halos[ihalo].pos.z - cur_part->pos[Z];
              dx = fabs(dx);
              dy = fabs(dy);
              dz = fabs(dz);
              if ( dx > 0.5 )  dx = 1.0 - dx;
              if ( dy > 0.5 )  dy = 1.0 - dy;
              if ( dz > 0.5 )  dz = 1.0 - dz;
              
              Dist2 = dx*dx + dy*dy + dz*dz;
              
              /* Is this particle outside the gathering radius? */
              if ( Dist2 > gatherRad2 )
                {
                  /* remove particle from host halo */
                  npart--;
                  
                  if(prev_part == NULL)
                    {
                      halos[ihalo].ll = cur_part->ll;
                    }
                  else
                    {
                      prev_part->ll = cur_part->ll;
                    }
                  
                  /* continue through ll */
                  // prev_part remains the same !
                  cur_part = cur_part->ll;
                }
              else
                {
                  /* continue through ll */
                  prev_part = cur_part;
                  cur_part  = cur_part->ll;                  
                }
            } /* while ( current!=NULL ) */
          
          /* check whether we called from a) or b) */
          if(halos[ihalo].npart == 0)
            {
              halos[ihalo].nll = npart;
#ifdef AHFverbose
              fprintf(stderr,"nll  =%12d\n",ihalo,halos[ihalo].nll);
              fprintf(io.logfile,"nll  =%12d\n",ihalo,halos[ihalo].nll);
#endif
            }
          else
            {
              halos[ihalo].npart = npart;
#ifdef AHFverbose
              fprintf(stderr,"npart=%12ld\n",ihalo,halos[ihalo].npart);
              fprintf(io.logfile,"npart=%12ld\n",ihalo,halos[ihalo].npart);
#endif
            }
        }      
    }
}

/****************************************************************************************
 * Prepare halos for construction
 *
 * as currently all particles are soley accessable via the linked list we store
 * the number of particles in the halo in halos[].nll rather than halos[].npart!
 ****************************************************************************************/ 
void prepare_haloes(void)
{
  long ihalo;
  
  for ( ihalo=numHalos-1; ihalo>=0; ihalo-- )
    {
      halos[ihalo].nll   = halos[ihalo].npart;
      halos[ihalo].npart = 0;
    }
}

/****************************************************************************************
 * Prepare halos for construction
 *
 * this routine ensures that each halo obtains all particles that possibly belong to it
 * -> those particles have to be stored in the array halos[].ipart[]
 ****************************************************************************************/ 
void prepare_particles(void)
{
  long ihalo;
  HALO *halo;
  
#ifdef AHFisodensity
  
#ifdef VERBOSE
  fprintf(stderr,"\n  => safe_return_particles() ... ");
  fprintf(io.logfile,"\n  => safe_return_particles() ... \n");
  fflush(io.logfile);
#endif
  
  /*****************************************************************************************
   * return all particles to the hosts' linked-lists but also store them in halos[].ipart[]
   *****************************************************************************************/ 
  for ( ihalo=numHalos-1; ihalo>=0; ihalo-- )
    safe_return_particles(ihalo);
  
#ifdef VERBOSE
  fprintf(stderr,"finished\n");
  fprintf(io.logfile,"  <= safe_return_particles() finished\n\n");
  fflush(io.logfile);
#endif
  
#else /* AHFisodensity */
  
  /****************************************************************************** 
   * 1. We are returning all particles to their respective hosts as
   *    gather_hostParts() below requires to collect >>all possible<< particles
   *
   * NOTE: this process has to be performed in serial!!!!
   ******************************************************************************/ 
#ifdef VERBOSE
  fprintf(stderr,"\n  => return_particles() ... ");
  fprintf(io.logfile,"\n  => return_particles() ... \n");
  fflush(io.logfile);
#endif
  
  for ( ihalo=numHalos-1; ihalo>=0; ihalo-- )
    return_particles(ihalo);
  
#ifdef VERBOSE
  fprintf(stderr,"finished\n");
  fprintf(io.logfile,"  <= return_particles() finished\n\n");
  fflush(io.logfile);
#endif
  
  /******************************************************************************/
  /* 2. Prune host halos (again) prior to gathering particles from them         */	
  /******************************************************************************/ 
#ifdef VERBOSE
  fprintf(stderr,"\n  => prune_host_halos() ... ");
  fprintf(io.logfile,"\n  => prune_host_halos() ... \n");
  fflush(io.logfile);
#endif
  
  //prune_host_halos();
  
#ifdef VERBOSE
  fprintf(stderr,"finished\n");
  fprintf(io.logfile,"  <= prune_host_halos() finished\n\n");
  fflush(io.logfile);
#endif
  
  
  /******************************************************************************/
  /* 3. Gather additional particles from the vicinity (e.g. the hosts)          */	
  /******************************************************************************/ 
  fprintf(stderr,"  => gather_hostParts() ... ");
  fprintf(io.logfile,"  => gather_hostParts() ... \n");
  fflush(io.logfile);
  
#ifdef WITH_OPENMP
#pragma omp parallel private(ihalo) shared(halos, numHalos)
#pragma omp for schedule(dynamic)
#endif
  for ( ihalo=0; ihalo<numHalos; ihalo++ )
#if (defined NEWSTARTRUN && WITH_AHF_HALOS_SFC)
    ahf_halos_sfc_gatherParts(halos+ihalo);
#else
  gather_hostParts(ihalo);
#endif
  
  fprintf(stderr,"finished\n");
  fprintf(io.logfile,"  <= gather_hostParts() finished\n\n");
  fflush(io.logfile);
  
  
  
  /******************************************************************************/
  /* 3. Merge particle linked-list ll with particle array ipart[] */	
  /******************************************************************************/
#if (!defined WITH_AHF_HALOS_SFC)
#ifdef WITH_OPENMP
#pragma omp parallel private(ihalo) shared(halos, numHalos)
#pragma omp for schedule(dynamic)
#endif
  for ( ihalo=0; ihalo<numHalos; ihalo++ )
    {
      halo = &halos[ihalo];
      merge_ll_and_ipart(halo);
    }
#endif /* !WITH_AHF_HALOS_SFC */
#endif /* AHFisodensity */
  
}

int	ConstructHalos (void) 
{
  int  ihalo;
  HALO *halo;
  
  /******************************************************************************/ 
  /* Prepare halos for construction
   * ------------------------------
   * the actual construction of halos[] below requires all
   * prospective halo particles in an array halos[].ipart[]
   *
   * currently the particles are
   *  a) uniquely assigned to only one halo
   *  b) ordered via a linked-list attached to halos[].ll
   ******************************************************************************/ 
  prepare_haloes();
  
  /* we now distinguish between halos[].ll,      halos[].nll  and  
   halos[].ipart[], halos[].npart */
  
  prepare_particles();
  
  /* now each halo has access to all prospective particles via the array ipart[] */
  
  
  /**********************************************************************/ 
  /**********************************************************************/ 
  /* Finally CONSTRUCTING the halos:
   *================================
   * 1. sort particles with respect to distance
   * 2. estimate Rvir and remove everything outside
   * 3. remove unbound particles
   * 4. re-estimate Rvir and remove everything outside
   * 5. calculate halo profiles and properties
   *
   *
   * NOTE: there is no need anymore to perform the loop in reverse order
   *       as each halo has access to all particles that may belong to it
   **********************************************************************/ 
  /**********************************************************************/ 
#ifdef VERBOSE
  fprintf(stderr,"  => halo construction ... ");
  fprintf(io.logfile,"\n  => HALO CONSTRUCTION ... \n");
#endif
#ifdef WITH_OPENMP
#pragma omp parallel private(ihalo, halo) shared(halos, numHalos)
#pragma omp for schedule(dynamic)
#endif
  for ( ihalo=0; ihalo<numHalos; ihalo++ )
    {
#ifdef VERBOSE
      fprintf(io.logfile,"\n  halo %10d   (npart=%16ld   pos = %16.8g %16.8g %16.8g):\n",
              ihalo, 
              halos[ihalo].npart,
              halos[ihalo].pos.x*x_fac, halos[ihalo].pos.y*x_fac, halos[ihalo].pos.z*x_fac);
      fflush(io.logfile);
#endif
      
      halo = &halos[ihalo];
      
      /******************************************************************************/
      /* Sorting the halo particles */	
      sort_halo_particles(halo);
      
#ifndef AHFisodensity
      /******************************************************************************/
      /* estimate R_vir based upon virial overdensity criterion 
       * and remove all particles outside R_vir
       * rem_outsideRvir() resets:
       *    halos[ihalo].npart
       *    halos[ihalo].M_vir
       *    halos[ihalo].R_vir
       *    halos[ihalo].ovdens
       *
       * NOTE: rem_unbound() needs to determine Phi0 as accurately as possible
       *       and hence we should get rid of all extra particles asap...           */	
      rem_outsideRvir(halo, 0);
#endif /* AHFisodensity */
      
      
#ifndef AHFnoremunbound
      /******************************************************************************/
      /* remove unbound particles 
       * rem_outsideRvir() resets:
       *    halos[ihalo].npart
       *    halos[ihalo].M_vir
       *    halos[ihalo].R_vir
       *
       * NOTE: R_vir is being set to the distance of the farthest bound particle
       */	
      rem_unbound(halo);
#else /* AHFnoremunbound */
      /******************************************************************************************************/
      /* We need to fill in the values that are left uncalculated if we do not remove the unbound particles */
      /* rem_nothing() resets:
       *    halos[ihalo].npart
       *    halos[ihalo].M_vir
       *    halos[ihalo].R_vir
       */
      rem_nothing(halo);
#endif /* AHFnoremunbound */
      
      
#ifndef AHFisodensity
      /******************************************************************************/
      /* determine virial overdensity radius again and remove all outliers,
       * this time solely using gravitationally bound particles 
       * rem_outsideRvir() resets:
       *    halos[ihalo].npart
       *    halos[ihalo].M_vir
       *    halos[ihalo].R_vir
       *    halos[ihalo].ovdens
       */	
      rem_outsideRvir(halo, 1);
#endif /* AHFisodensity  */
      
      
      /******************************************************************************/
      /* calculate halo profiles (NOTE: requires an estimate for R_vir!)*/	
#	if (defined WITH_MPI || defined AHFrestart)
      if (halos[ihalo].ignoreme == FALSE)
#	endif
        {
          if ( HaloProfiles(halo) == 0 )
            {
              fprintf(stderr,"Stuffed up calculating profile of halo %d\n",ihalo);
              exit(-1);
            }
#ifdef AHFphspdens
          if ( HaloProfilesPhaseSpace(halo) == 0 )
            {
              fprintf(stderr,"Stuffed up calculating profile of halo %d\n",ihalo);
              exit(-1);
            }
#endif
        }
      
    }
  
#ifdef VERBOSE
  fprintf(stderr,"finished\n\n");
  fprintf(io.logfile,"\n  <= HALO CONSTRUCTION FINISHED\n\n");
#endif
  
  return(TRUE);
}

/*
 ************************************************************
 ************************************************************
 * Merge particles stored via linked-list halos[].ll and in array halos[].ipart[]
 */
void merge_ll_and_ipart(HALO *halo)
{
  partptr       cur_part;
  long unsigned jpart;
  
  /* remember first new ipart[] particle */
  jpart = halo->npart;
  
  /* increment total number of particles */
  halo->npart += halo->nll;
  
  /* reallocate memory for ipart[] array */
  halo->ipart = (long unsigned *) realloc(halo->ipart, halo->npart*sizeof(long unsigned));
  
  /* attach particles from ll to ipart[] */
  cur_part = halo->ll;
  while(cur_part != NULL)
    {
      halo->ipart[jpart] = cur_part - global.fst_part;
      jpart++;
      cur_part = cur_part->ll;
    }
  
  /* even though halos[].ll is not used anymore, reset it to NULL */
  halo->ll  = NULL;
  halo->nll = 0;
}

#ifdef OPT_REALLOC
/*=======================================================================
 * actually gather the particles from the host
 * NOTE: we are gathering from halos[].ll and storing in halos[].ipart[]
 *      >>> this variant is minimizing the usage of realloc() <<<
 *=======================================================================*/
void gather_parts_from_host(int num, int host)
{
  partptr current;
  double  dx,dy,dz,Dist2,gatherRad2;
  long    nhalo, nhost, nbuffer;
  long unsigned *test;
  
  /* set the gathering radius */
  gatherRad2 = pow2(halos[num].gatherRad);
  
  /* we can at max gather all host particles */
  nhalo   = halos[num].npart;
  nhost   = halos[host].nll;
  nbuffer = 1000;
  halos[num].ipart = (long unsigned *) realloc(halos[num].ipart, (halos[num].npart+nbuffer)*sizeof(long unsigned));
  if(halos[num].ipart == NULL)
    {
      fprintf(stderr,"gather_parts_from_host: running out of memory!\n");
      fprintf(stderr,"                        ihalo=%12d npart=%12ld nll=%12ld\n",num,halos[num].npart,halos[num].nll);
      fprintf(stderr,"                        ihost=%12d npart=%12ld nll=%12ld\n",host,halos[host].npart,halos[host].nll);
      fprintf(stderr,"recompile AHFstep with -DREALLOC\n");
      exit(0);
    }
  
  /* Loop over all the particles in the host halo */
  current = halos[host].ll;
  
  while ( current != NULL )
    {
      /* Calculating the distance between the particle in the host halo and the centre of the current halo */
      dx = halos[num].pos.x - current->pos[X];
      dy = halos[num].pos.y - current->pos[Y];
      dz = halos[num].pos.z - current->pos[Z];
      dx = fabs(dx);
      dy = fabs(dy);
      dz = fabs(dz);
      if ( dx > 0.5 )  dx = 1.0 - dx;
      if ( dy > 0.5 )  dy = 1.0 - dy;
      if ( dz > 0.5 )  dz = 1.0 - dz;
      
      Dist2 = dx*dx + dy*dy + dz*dz;
      
      /* Is this particle within the gathering radius? */
      if ( Dist2 <= gatherRad2 )
        {
          /* copy particle to halo's ipart[] array */
          halos[num].npart++;
          
          if(halos[num].npart >= nhalo+nbuffer)
            {
              nhalo += nbuffer;
              
              halos[num].ipart = (long unsigned *) realloc(halos[num].ipart, (halos[num].npart+nbuffer)*sizeof(long unsigned));
              if(halos[num].ipart == NULL)
                {
                  fprintf(stderr,"gather_parts_from_host: running out of memory!\n");
                  fprintf(stderr,"                        ihalo=%12d npart=%12ld nll=%12ld\n",num,halos[num].npart,halos[num].nll);
                  fprintf(stderr,"                        ihost=%12d npart=%12ld nll=%12ld\n",host,halos[host].npart,halos[host].nll);
                  fprintf(stderr,"recompile AHFstep with -DREALLOC\n");
                  exit(0);
                }
            }
          
          halos[num].ipart[halos[num].npart-1] = current - global.fst_part;
        }
      
      /* continue through ll */
      current = current->ll;
    } /* while ( current!=NULL ) */
  
  /* adjust to only store the actual particles */
  halos[num].ipart = (long unsigned *) realloc(halos[num].ipart, halos[num].npart*sizeof(long unsigned));
}
#else /* NO_REALLOC */
/*=======================================================================
 * actually gather the particles from the host
 * NOTE: we are gathering from halos[].ll and storing in halos[].ipart[]
 *=======================================================================*/
void gather_parts_from_host(int num, int host)
{
  partptr current;
  double  dx,dy,dz,Dist2,gatherRad2;
  
  /* set the gathering radius */
  gatherRad2 = pow2(halos[num].gatherRad);
  
  /* Loop over all the particles in the host halo */
  current = halos[host].ll;
  
  while ( current != NULL )
    {
      /* Calculating the distance between the particle in the host halo and the centre of the current halo */
      dx = halos[num].pos.x - current->pos[X];
      dy = halos[num].pos.y - current->pos[Y];
      dz = halos[num].pos.z - current->pos[Z];
      dx = fabs(dx);
      dy = fabs(dy);
      dz = fabs(dz);
      if ( dx > 0.5 )  dx = 1.0 - dx;
      if ( dy > 0.5 )  dy = 1.0 - dy;
      if ( dz > 0.5 )  dz = 1.0 - dz;
      
      Dist2 = dx*dx + dy*dy + dz*dz;
      
      /* Is this particle within the gathering radius? */
      if ( Dist2 <= gatherRad2 )
        {
          /* copy particle to halo's ipart[] array */
          halos[num].npart++;
          halos[num].ipart = (long unsigned *) realloc(halos[num].ipart, halos[num].npart*sizeof(long unsigned));
          halos[num].ipart[halos[num].npart-1] = current - global.fst_part;
        }
      
      /* continue through ll */
      current = current->ll;
      
    } /* while ( current!=NULL ) */
}
#endif /* NO_REALLOC */

/*=======================================================================
 * Gathers particles within the gathering radius from all possible hosts
 *=======================================================================*/
void gather_hostParts(int num)
{
  int i,hostHalo,hostHaloGrids,hostGrids,count;
  
  hostHaloGrids = -1;
  hostHalo 		 = -1;
  hostGrids     = -1; 
  
  /*****************************************************************************************************/	
  /* Gather from hostHalo */
  hostHalo = halos[num].hostHalo;
  
  if ( hostHalo >= 0 )
    {
      /* locate the top-level host for this halo */
      while(halos[hostHalo].hostHalo >= 0)
        hostHalo = halos[hostHalo].hostHalo;
      
#ifdef AHFverbose
      fprintf(io.logfile,"    gather_hostParts: halos[%d]: npart=%12ld nll=%12ld (host         =%12d, nll=%12d, gatherRad=%16.8g kpc/h) -> ",
              num,
              halos[num].npart, halos[num].nll, 
              hostHalo, halos[hostHalo].nll,
              halos[num].gatherRad*x_fac*1000.);
      fflush(io.logfile);
#endif
      /* actually gather particles from hostHalo */
      gather_parts_from_host(num, hostHalo);
      
#ifdef AHFverbose
      fprintf(io.logfile,"npart=%12ld nll=%12ld\n",halos[num].npart,halos[num].nll);
      fflush(io.logfile);
#endif       
    }
  
  
  /*****************************************************************************************************/	
  /* Gather from hostHaloGrids */
  hostHaloGrids = halos[num].hostHaloGrids;
  
  if ( (hostHaloGrids >= 0) && (hostHaloGrids != hostHalo) ) 
    {
#ifdef AHFverbose
      fprintf(io.logfile,"    gather_hostParts: halos[%d]: npart=%12ld nll=%12ld (hostHaloGrids=%12d, nll=%12d, gatherRad=%16.8g kpc/h) -> ",
              num,
              halos[num].npart, halos[num].nll, 
              hostHaloGrids, halos[hostHaloGrids].nll,
              halos[num].gatherRad*x_fac*1000.);
      fflush(io.logfile);
#endif
      /* actually gather particles from hostHaloGrids */
      gather_parts_from_host(num, hostHaloGrids);
      
#ifdef AHFverbose
      fprintf(io.logfile,"npart=%12ld nll=%12ld\n",halos[num].npart,halos[num].nll);
      fflush(io.logfile);
#endif       
    } 
  
  
  /*****************************************************************************************************/	
  /* Gather from other hostHaloGrids 
   * I.e. gathering from the halos host's grids and their hosts etc... */
  
  for ( i=0; i<halos[num].numHostHaloGrids; i++ )
    {
      hostGrids = halos[num].hostGrids[i];
      count=0;
      
      if (hostGrids == hostHalo)
        {
          /* fprintf(stderr,"We have pillaged this host || halos[%d].hostGrids[%d] = %d\n",num,i,hostGrids); */
          continue; /* I.e. we have pillaged this host */
        }
      
      if (hostGrids == hostHaloGrids)
        {
          /* fprintf(stderr,"We have pillaged this host || halos[%d].hostGrids[%d] = %d\n",num,i,hostGrids); */
          continue; /* I.e. we have pillaged this host */
        }
      
      if ( (hostGrids >= 0)  )
        {
#ifdef AHFverbose
          fprintf(io.logfile,"    gather_hostParts: halos[%d]: npart=%12ld nll=%12ld (hostGrids    =%12d, nll=%12d, gatherRad=%16.8g kpc/h) -> ",
                  num,
                  halos[num].npart, halos[num].nll, 
                  hostGrids, halos[hostGrids].nll,
                  halos[num].gatherRad*x_fac*1000.);
          fflush(io.logfile);
#endif
          
          /* actually gather particles from host */
          gather_parts_from_host(num, hostGrids);
          
#ifdef AHFverbose
          fprintf(io.logfile,"npart=%12ld nll=%12ld\n",halos[num].npart,halos[num].nll);
          fflush(io.logfile);
#endif                 
        } 
      
      if ( count == 0 ) 
        break; /* No particles gathered from this grid level - thus you will not gather any from the next level!! */
      
    } /*  ( i=0; i<halos[num].numHostHaloGrids; i++ ) */
  
}

/*=======================================================================
 * Gathers particles within the gathering radius from all possible hosts
 *
 * NOTE: we are gathering from halos[].ipart[] and storing in halos[].ipart[]
 *=======================================================================*/
void gather_hostParts_ipart(int num)
{
  
  int      hostHalo, hostHaloGrids, hostGrids;
  partptr  cur_part;
  long     jpart;
  double   dx, dy, dz, Dist2, gatherRad2;
  
  int      i, count;
  
  hostHaloGrids	= -1;
  hostHalo 		= -1;
  hostGrids 		= -1; 
  
#ifdef AHFverbose
  /*****************************************************************************************************/	
  /* check current halo */
  if ( halos[num].npart == 0)
    fprintf(stderr,"\n    => gather_hostParts(): halo %10d, npart=%12ld, numNodes=%8d",num,halos[num].npart,halos[num].numNodes);
#endif
  
  /*****************************************************************************************************/	
  /* set the gathering radius */
  gatherRad2 = pow2(halos[num].gatherRad);
  
  
  /*****************************************************************************************************/	
  /* 1. gather from hostHalo */
  hostHalo   = halos[num].hostHalo;
  
  if ( hostHalo >= 0 )
    {
      
      /* Loop over all the particles in the host halo */
      for(jpart=0; jpart<halos[hostHalo].npart; jpart++)
        {
          cur_part = global.fst_part + halos[hostHalo].ipart[jpart];
          
          /* Calculating the distance between the particle in the host halo and the centre of the current halo */
          dx = halos[num].pos.x - cur_part->pos[0];
          dx = fabs(dx);
          dy = halos[num].pos.y - cur_part->pos[1];
          dy = fabs(dy);
          dz = halos[num].pos.z - cur_part->pos[2];
          dz = fabs(dz);
          if ( dx > 0.5 )  dx = 1.0 - dx;
          if ( dy > 0.5 )  dy = 1.0 - dy;
          if ( dz > 0.5 )  dz = 1.0 - dz;
          
          Dist2 = dx*dx + dy*dy + dz*dz;
          
          /* Is this particle within the gathering radius? */
          if ( Dist2 < gatherRad2 )
            { 
              /* add particle to halos[num] */
              halos[num].npart++;
              halos[num].ipart                     = (long unsigned *) realloc(halos[num].ipart, halos[num].npart*sizeof(long unsigned));
              halos[num].ipart[halos[num].npart-1] = halos[hostHalo].ipart[jpart];
            }
        }
      
    }
  
  
  /*****************************************************************************************************/	
  /* 2. gather from hostHaloGrids */
  hostHaloGrids = halos[num].hostHaloGrids;
  
  if ( (hostHaloGrids >= 0) && (hostHaloGrids != hostHalo) ) 
    {
      
      /* Loop over all the particles in the host halo */
      for(jpart=0; jpart<halos[hostHaloGrids].npart; jpart++)
        {
          cur_part = global.fst_part + halos[hostHaloGrids].ipart[jpart];
          
          /* Calculating the distance between the particle in the host halo and the centre of the current halo */
          dx = halos[num].pos.x - cur_part->pos[0];
          dx = fabs(dx);
          dy = halos[num].pos.y - cur_part->pos[1];
          dy = fabs(dy);
          dz = halos[num].pos.z - cur_part->pos[2];
          dz = fabs(dz);
          
          if ( dx > 0.5 ) dx = 1.0 - dx;
          if ( dy > 0.5 ) dy = 1.0 - dy;
          if ( dz > 0.5 ) dz = 1.0 - dz;
          
          Dist2 = dx*dx + dy*dy + dz*dz;
          
          /* Is this particle within the gathering radius? */
          if ( Dist2 < gatherRad2 )
            { 
              /* add particle to halos[num] */
              halos[num].npart++;
              halos[num].ipart                     = (long unsigned *) realloc(halos[num].ipart, halos[num].npart*sizeof(long unsigned));
              halos[num].ipart[halos[num].npart-1] = halos[hostHaloGrids].ipart[jpart];
            }
        }
      
    } 
  
  
  /*****************************************************************************************************/	
  /* 3. gather from other hostHaloGrids 
   * I.e. gathering from the halos host's grids and their hosts etc... */
  
  for ( i=0; i<halos[num].numHostHaloGrids; i++ )
    {
      
      hostGrids = halos[num].hostGrids[i];
      count=0;
      
      if (hostGrids == hostHalo)
        {
          /* fprintf(stderr,"We have pillaged this host || halos[%d].hostGrids[%d] = %d\n",num,i,hostGrids); */
          continue; /* I.e. we have pillaged this host */
        }
      
      if (hostGrids == hostHaloGrids)
        {
          /* fprintf(stderr,"We have pillaged this host || halos[%d].hostGrids[%d] = %d\n",num,i,hostGrids); */
          continue; /* I.e. we have pillaged this host */
        }
      
      
      if ( (hostGrids >= 0)  )
        {
          /* Loop over all the particles in the host halo */
          for(jpart=0; jpart<halos[hostGrids].npart; jpart++)
            {
              cur_part = global.fst_part + halos[hostGrids].ipart[jpart];
              
              /* Calculating the distance between the particle in the host halo and the centre of the current halo */
              dx = halos[num].pos.x - cur_part->pos[0];
              dx = fabs(dx);
              dy = halos[num].pos.y - cur_part->pos[1];
              dy = fabs(dy);
              dz = halos[num].pos.z - cur_part->pos[2];
              dz = fabs(dz);
              
              if ( dx > 0.5 ) dx = 1.0 - dx;
              if ( dy > 0.5 ) dy = 1.0 - dy;
              if ( dz > 0.5 ) dz = 1.0 - dz;
              
              Dist2 = dx*dx + dy*dy + dz*dz;
              
              /* Is this particle within the gathering radius? */
              if ( Dist2 < gatherRad2 )
                { 
                  /* if we do not gather from the first level we will not gather from any other level... */
                  count++;
                  
                  /* add particle to halos[num] */
                  halos[num].npart++;
                  halos[num].ipart                     = (long unsigned *) realloc(halos[num].ipart, halos[num].npart*sizeof(long unsigned));
                  halos[num].ipart[halos[num].npart-1] = halos[hostGrids].ipart[jpart];
                }
            }
          
        } 
      
      if ( count == 0 ) 
        break; /* No particles gathered from this grid level - thus you will not gather any from the next level!! */
      
    } /*  ( i=0; i<halos[num].numHostHaloGrids; i++ ) */
  
  return;
}



/*
 ***********************************************************************
 ***********************************************************************
 * fills those values in that normally rem_unbound() would fill in... 
 */
void rem_nothing(HALO *halo)
{
  long unsigned jpart, npart;
  double        Xc, Yc, Zc, Xp, Yp, Zp, dX, dY, dZ, dist2;
  partptr       cur_part;
  double        R_vir, M_vir, weight;
  
  
#ifndef TRACKER
  if(halo->npart < AHF_MINPART)
#else
    if(halo->npart < TRK_MINPART)
#endif
      return;
  
#ifdef VERBOSE
  fprintf(io.logfile,"    rem_nothing:      npart=%12ld -> ",halo->npart);
  fflush(io.logfile);
#endif
  
  
  Xc	   = halo->pos.x;
  Yc	   = halo->pos.y;
  Zc	   = halo->pos.z;
  R_vir	   = -1.0;
  M_vir	   = 0.0;
  npart    = 0;
  
  /* loop over all particles... */
  for(jpart=0; jpart<halo->npart; jpart++)
    {
      /* access particle */
      cur_part = global.fst_part + halo->ipart[jpart];
      
#ifdef MULTIMASS
      weight = (double)cur_part->weight;
#else
      weight = (double)1.0;
#endif
      
      /* cumulative mass */
      M_vir += weight;
      npart++;
      
      /* particle position */
      Xp  = (double)cur_part->pos[X];
      Yp  = (double)cur_part->pos[Y];
      Zp  = (double)cur_part->pos[Z];
      
      /* put particle into halo rest frame */
      dX  = fabs(Xp - Xc);
      dY  = fabs(Yp - Yc);
      dZ  = fabs(Zp - Zc);
      
      /* take care of periodic boundary conditions */
      if(dX >  0.5) dX -= 1.0;
      if(dY >  0.5) dY -= 1.0;
      if(dZ >  0.5) dZ -= 1.0;
      
      /* finally calculate distance squared */
      dist2 = pow2(dX) + pow2(dY) + pow2(dZ);
      
      /* remember distance of current particle (particles are ordered distancewise!) */
      if( dist2 > R_vir ) 
        R_vir  = dist2;
      
    } /* particle for loop */
  
  halo->npart = npart;
  halo->M_vir = M_vir;
  halo->R_vir = sqrt(R_vir);
  halo->Phi0  = 0.0;
  
#ifdef VERBOSE
  fprintf(io.logfile,"%12ld\n",halo->npart);
  fflush(io.logfile);
#endif
}



/*
 ************************************************************
 ************************************************************
 * Removes the unbound particles from the halos 
 */
void rem_unbound(HALO *halo)
{
  
  long unsigned  nremove, npart_old;
  double         Xc, Yc, Zc, VXc, VYc, VZc;
  double         Xp, Yp, Zp, VXp, VYp, VZp;
  double         dX, dY, dZ, dVX, dVY, dVZ;
  double         v2_tune, weight;
  double         Phi, Phi0, Phi_infty, M_r, M_vir, M_vel, vel2, v_esc2, d_prev, R_vir;
  double         I_now, I_mid, I_prev, dr;
  partptr        cur_part, pre_part, host_part, tmp_part;
  gasptr         cur_gas;
  double         x,y,z,dx,dy,dz,distA,distB,dist;
  double         xx,yy,zz;
  int            niter;
  int 	         hostHalo;
  long           jpart;
  long unsigned *bound_ipart, bound_npart;
  double        *mom2;
  long unsigned *idx, no_vbulk;
  
#ifdef TRACKER_DEBRIS
  long unsigned  *ubound_ipart, ubound_npart;
  partptr        upart;
#endif
  
#ifndef TRACKER
  if(halo->npart < AHF_MINPART)
#else
    if(halo->npart < TRK_MINPART)
#endif
      {
#ifdef TRACKER_DEBRIS
        /* reset nupart to 0 for this halo, just in case it was not 0 before */
        halo->nupart=0;
#endif  
        return;
      }
  
#ifdef VERBOSE
  fprintf(io.logfile,"    rem_unbound:      npart=%12ld -> ",halo->npart);
  fflush(io.logfile);
#endif
  
  /* velocity tune parameter */
#ifndef TRACKER
  v2_tune = pow2(AHF_VTUNE);
#else
  v2_tune = pow2(TRK_VTUNE);
#endif
  
  /* how many central particles to use for the initial bulk velocity guess */
#ifndef TRACKER
  no_vbulk = (long unsigned) (AHF_MINPART/2); 
#else
  no_vbulk = (long unsigned) (TRK_MINPART/2); 
#endif
  
  /* remember initial number of particles */
  /*  npart_old = halos[i].npart; */
  
  /* halo centre in AMIGA units */
  Xc = halo->pos.x;
  Yc = halo->pos.y;
  Zc = halo->pos.z;
  
#ifdef TRACKER_DEBRIS
  ubound_npart = 0;
  ubound_ipart = NULL;
#endif
  
  nremove = 4;
  niter   = 0;
  /*----------------------------------------------------------------------------
   * remove particles until all particles are bound or halo mass gets too small
   *----------------------------------------------------------------------------*/
  while( (nremove > 3))
    {      
      /* iteration counter */
      niter++;
      
      /*---------------------------------------------------------------
       * determine Phi0  (the zero point of the potential is infinity)
       *---------------------------------------------------------------*/
      /* reset values */
      I_prev  = 0.0;
      d_prev  = 0.0;
      M_r     = 0.0;
      Phi0    = 0.0;
      
      /************************************************************/
      /* loop over all sorted particles from inside out */
      /* calculate Phi0 */
      for(jpart=0; jpart<halo->npart; jpart++)
        {
          /* access particle */
          cur_part = global.fst_part + halo->ipart[jpart];
          
#ifdef MULTIMASS
          weight = (double)cur_part->weight;
#else
          weight = (double)1.0;
#endif
          
          /* cumulative mass */
          M_r += weight;
          
          /* particle position */
          Xp  = (double)cur_part->pos[X];
          Yp  = (double)cur_part->pos[Y];
          Zp  = (double)cur_part->pos[Z];
          
          /* put particle into halo rest frame */
          dX  = fabs(Xp - Xc);
          dY  = fabs(Yp - Yc);
          dZ  = fabs(Zp - Zc);
          
          /* take care of periodic boundary conditions */
          if(dX >  0.5) dX -= 1.0;
          if(dY >  0.5) dY -= 1.0;
          if(dZ >  0.5) dZ -= 1.0;
          
          /* finally calculate distance (conversion to dist_phys=a*dist via phi_fac!)  */
          dist = sqrt( pow2(dX) + pow2(dY) + pow2(dZ) );
          
          /* accumulate Phi0 */
          if( dist > ZERO )
            {
              
              /* mid-point integration */
              I_now = M_r/pow2(dist);
              I_mid = (I_now + I_prev)/2.;
              dr    =  dist - d_prev;
              
              /* accumulate Phi0 */
              Phi0 += I_mid * dr;
            }
          
          d_prev   = dist;
          I_prev   = I_now;
        } /* particle loop for Phi0 determination */
      
      /* finally calculate Phi0 */
      Phi0 += M_r/dist;
      
      /* remember Phi0 as it will be used when calculating halos[].prof.Epot */
      halo->Phi0 = Phi0;
      
      /*--------------------------------------------------------------------
       * determine Phi             ( v_esc^2 = 2 * |Phi| )
       *--------------------------------------------------------------------*/
      /* reset values */
      nremove = 0;
      I_prev  = 0.0;
      d_prev  = 0.0;
      M_r     = 0.0;
      M_vir   = 0.0;
      Phi     = 0.0;
      
      /* initial bulk velocity of halo */
      if(niter == 1)
        {
          mom2 = (double *)        calloc(no_vbulk + 1, sizeof(double));
          idx  = (long unsigned *) calloc(no_vbulk + 1, sizeof(long unsigned));
          
          for(jpart=0; jpart<no_vbulk; jpart++)
            {
              cur_part       = global.fst_part + halo->ipart[jpart];
              mom2[jpart+1]  = pow2(cur_part->mom[X])+pow2(cur_part->mom[Y])+pow2(cur_part->mom[Z]);
            }
          indexx((long unsigned)(no_vbulk), mom2, idx);
          
          /* use the median of the innermost particle's velocity */
          jpart = idx[(int)(no_vbulk/2)] - 1;
          
          free(mom2);
          free(idx);
        }
      else
        {
          /* use the most central particle's velocity */
          jpart = 0;
        }
      
      cur_part = global.fst_part + halo->ipart[jpart];
#ifdef MULTIMASS
      weight   = (double)cur_part->weight;
#else
      weight   = (double)1.0;
#endif
      M_vel    = weight;
      VXc      = weight*cur_part->mom[X];
      VYc      = weight*cur_part->mom[Y];
      VZc      = weight*cur_part->mom[Z];
      
      /************************************************************/
      /* loop over all sorted particles from inside out */
      bound_npart = 0;
      bound_ipart = calloc(1, sizeof(long unsigned));    // some realloc()'s do not like NULL pointers...
      for(jpart=0; jpart<halo->npart; jpart++)
        {
          /* access particle */
          cur_part = global.fst_part + halo->ipart[jpart];
          
#ifdef MULTIMASS
          weight = (double)cur_part->weight;
#else
          weight = (double)1.0;
#endif
          /* cumulative mass */
          M_r += weight;
          
          /* particle position */
          Xp  = (double)cur_part->pos[X];
          Yp  = (double)cur_part->pos[Y];
          Zp  = (double)cur_part->pos[Z];
          
          /* put particle into halo rest frame :: *no* fabs() this time ! */
          dX  = (Xp - Xc);
          dY  = (Yp - Yc);
          dZ  = (Zp - Zc);
          
          /* take care of periodic boundary conditions */
          if(dX >  0.5) dX -= 1.0;
          if(dY >  0.5) dY -= 1.0;
          if(dZ >  0.5) dZ -= 1.0;
          if(dX < -0.5) dX += 1.0;
          if(dY < -0.5) dY += 1.0;
          if(dZ < -0.5) dZ += 1.0;
          
          /* finally calculate distance (conversion to dist_phys=a*dist via phi_fac!) */
          dist = sqrt( pow2(dX) + pow2(dY) + pow2(dZ) );
          
          /* get potential escape velocity */
          if( dist > ZERO )
            {  	  
              /* mid-point integration */
              I_now = M_r/pow2(dist);
              I_mid = (I_now + I_prev)/2.;
              dr    =  dist - d_prev;
              
              /* accumulate potential */
              Phi += I_mid * dr;
              
              /* get escape velocity */
              v_esc2 = (2*fabs(Phi-Phi0)*phi_fac);            
            }
          
          /* potential="inf" for dist=0 -> set v_esc manually */
          else
            {
              v_esc2 = 1e30;
            }
          
          /* get particle velocity in halo rest frame velocity plus Hubble flow */
          VXp = (double)cur_part->mom[X];
          VYp = (double)cur_part->mom[Y];
          VZp = (double)cur_part->mom[Z];
          
          dVX = (VXp - VXc/M_vel) * v_fac + Hubble * dX * r_fac;
          dVY = (VYp - VYc/M_vel) * v_fac + Hubble * dY * r_fac;
          dVZ = (VZp - VZc/M_vel) * v_fac + Hubble * dZ * r_fac;
          
          /* absolute velocity of current particle */
          vel2  = pow2(dVX) + pow2(dVY) + pow2(dVZ);
          
#ifdef GAS_PARTICLES
#if (defined WITH_MPI || defined AHFrestart)
          /* u = 3/2  kT/m  = 1/2  v_therm^2 ?! */          
          vel2   += (cur_part->u<PGAS ? 0.0 : 2 * cur_part->u/pow2(v_fac));
#else /* WITH_MPI */
          /* cur_part == gas particle => add thermal velocity */
          if(cur_part-global.fst_part-global.offset_gas < global.no_gas && cur_part-global.fst_part-global.offset_gas > 0)
            {
              cur_gas = global.fst_gas + (cur_part-global.fst_part-global.offset_gas);
              
              /* u = 3/2  kT/m  = 1/2  v_therm^2 ?! */
              vel2   += 2*cur_gas->u / pow2(v_fac);            
            }
#endif /* WITH_MPI */
#endif /* GAS_PARTICLES */
          
          /* unbound particle? */
          if( vel2 > v2_tune*v_esc2 )
            {
              /* count number of particles to be removed */
              nremove++;
#ifdef TRACKER_DEBRIS 
              /* store unbound particle */
              ubound_npart++;
              ubound_ipart = (long unsigned *) realloc(ubound_ipart, ubound_npart*sizeof(long unsigned));
              ubound_ipart[ubound_npart-1] = cur_part - global.fst_part;              
              
#endif
            }         
          else 
            {
              /* let particle contribute to bulk/mean velocity of halo */
              M_vel += weight;
              VXc   += weight*cur_part->mom[X];
              VYc   += weight*cur_part->mom[Y];
              VZc   += weight*cur_part->mom[Z];
              
              /* store bound particles temporarily in bound_ipart[] */
              bound_npart++;
              bound_ipart = (long unsigned *) realloc(bound_ipart, bound_npart*sizeof(long unsigned));
              bound_ipart[bound_npart-1] = cur_part - global.fst_part;
              
              /* accumulate virial values (NOTE: for v_esc2 we are accumulating M_r and not M_vir!) */
              M_vir += weight;
              R_vir  = dist;
            }
          
          I_prev   = I_now;
          d_prev   = dist;
          
        } /* particle loop */
      
      /* double-check new number of bound particles */
      if(bound_npart != (halo->npart - nremove))
        fprintf(stderr,"rem_unbound: better check the unbinding procedure! bound_part=%ld vs. halos[num].npart-nremove=%ld\n",
                bound_npart, halo->npart-nremove);
      
      
      /* update number of particles in halo */
      free(halo->ipart);
      halo->ipart = bound_ipart;
      halo->npart = bound_npart;
      halo->M_vir = M_vir;
      halo->R_vir = R_vir;
      
      
      /* there is no point in removing unbound particles once the halo is too small */
#ifndef TRACKER
      if ( ((double)(halo->npart)) < ((double)(AHF_MINPART)) )
#else
        if ( ((double)(halo->npart)) < ((double)(TRK_MINPART)) )
#endif
          {
            /*      fprintf(stderr,"$$BREAK$$ halos[%d].npart(%d) < MINPART(%d)\n",num,halos[num].npart,AHF_MINPART);*/
            break;
          }
      
    } /* while( (nremove > n)  ) */
  
#ifdef TRACKER_DEBRIS
  /* store new list of unbound particles */
  free(halo->iupart);
  halo->iupart = ubound_ipart;
  
  /* update number of unbound particles in halo */
  halo->nupart = ubound_npart;
  fprintf(stderr,"unbound, new npart: %ld %12ld\n\n", halo->nupart, halo->npart);
#endif
  
#ifdef VERBOSE
  fprintf(io.logfile,"%12ld (Rvir=%g kpc/h)\n",halo->npart, halo->R_vir*x_fac*1000.);
  fflush(io.logfile);
  
#endif
}


/*
 ************************************************************
 ************************************************************
 * Calculate integral properties of halos
 */

/* qsort stuff */
int haloCompare(const void *p1, const void *p2) {
  
  if ( ((HALO *)p1)->npart < ((HALO *)p2)->npart )
    return(1);
  else if ( ((HALO *)p1)->npart > ((HALO *)p2)->npart )
    return(-1);
  else
    return(0);
}


/*
 ************************************************************
 ************************************************************
 *  qsort stuff
 */
int refCompare(const void *p1,const  void *p2) {
  
  if ( ((SPATIALREF *)p1)->volume < ((SPATIALREF *)p2)->volume )
    return(1);
  else if ( ((SPATIALREF *)p1)->volume > ((SPATIALREF *)p2)->volume )
    return(-1);
  else
    return(0);
}

void sortd( length, a, ind )
int length;
double a[ ];
int ind[ ];
{
  
  int i, ii, ij, j, m, m1, n2;
  double t;
  
  for ( i = 0 ; i < length; i++ ) ind[ i ] = i;
  m = 1;
  n2 = length / 2;
  m = 2 * m;
  while ( m <= n2 ) m = 2 * m;
  m = m - 1;
three:;
  m1 = m + 1;
  for ( j = m1-1 ; j < length; j++ ) {
    t = a[ ind[ j ] ];
    ij = ind[ j ];
    i = j - m;
    ii = ind[ i ];
  four:;
    if ( t > a[ ii ] ) {
      ind[ i+m ] = ii;
      i = i - m;
      if ( i >= 0 ) {
        ii = ind[ i ];
        goto four;
      }
    }
    ind[ i+m ] = ij;
  }
  m = m / 2;
  if ( m > 0 ) goto three;
  return;
}

/*
 ************************************************************
 ************************************************************
 * Remove particles outside virial radius
 */
void	rem_outsideRvir(HALO *halo, int icall)
{
  
  partptr       cur_part, host_part, first_part;
  double        M_sphere, ovdens_sphere, lR_sphere, R_sphere, V_sphere;
  double        Xc, Yc, Zc, Xp, Yp, Zp, dX, dY, dZ, weight, dummy;
  long unsigned npart_sphere, jpart;
  int           hostHalo, nbins, ibin;
  double        dist_min, dist_max, ldist_min, ldist_max, ldr, *ovdens, *r, R_edge, ovdens_min, cur_dist;
  
#ifdef TRACKER_DEBRIS
  long          nopart;
#endif
  
#ifndef TRACKER
  if(halo->npart < AHF_MINPART)
#else
    if(halo->npart < TRK_MINPART)
#endif
      {
#ifdef TRACKER_DEBRIS
        /* set nopart to 0, to erase possible values from prev. timestep */
        if (icall == 0) halo->nopart = 0;
#endif
        return;
      }  
  
#ifdef VERBOSE
  fprintf(io.logfile,"    rem_outsideRvir:  npart=%12ld -> ",halo->npart);
  fflush(io.logfile);
#endif
  
  /* halo position in AMIGA units */
  Xc = halo->pos.x;
  Yc = halo->pos.y;
  Zc = halo->pos.z;
  
#ifdef AHFprofilerise
  
  /*============================================================================================
   * this part checks for an upturn in the density profile and chops the halo off at that point 
   *
   * NOTE: there are two calls to rem_outsideRvir()! 
   *       - one before rem_unbound() and
   *       - one after  rem_unbound() 
   *
   *       we should only look for the upturn before rem_unbound() as only then
   *       the current halo may or may not be embedded within the background of the host
   *       further, it saves time to not do it...
   *============================================================================================*/
  
  if(icall == 0)
    {
      /* how many bins from where to where for density profile? */
      binning_parameter(*halo, &nbins, &dist_min, &dist_max);
      
#if (definded AHFstep || defined TRACKER)
      fprintf(io.logfile,"(nbins=%d, dist_min=%g dist_max=%g kpc/h) ",nbins,dist_min*x_fac*1000.,dist_max*x_fac*1000.);
      fflush(io.logfile);
#endif
      
      /* logarithmical binning from dist_min to dist_max */
      ldist_min = log10(dist_min);
      ldist_max = log10(dist_max);
      ldr       = (ldist_max-ldist_min)/(double)nbins;
      
      /* create profile arrays */
      ovdens = (double *) calloc(nbins, sizeof(double));
      r      = (double *) calloc(nbins, sizeof(double));
      
      /* first particle */
      jpart    = 0;
      cur_dist = -1.0;
      
      /* calculate binned density profile */
      for(ibin=0; ibin<nbins; ibin++)
        {
          /* get current outer radius using logarithmic radial bins */
          lR_sphere = ldist_min + ((double)ibin+1) * ldr;
          R_sphere  = pow(10., lR_sphere);
          
          /* this heavily assumes that particles are ordered distancewise... */
          while(cur_dist < R_sphere && jpart < halo->npart)
            { 
              /* access particle */
              cur_part = global.fst_part + halo->ipart[jpart];
              
#ifdef MULTIMASS
              weight = (double)cur_part->weight;
#else
              weight = (double)1.0;
#endif
              /* accumulate mass */
              M_sphere += weight;
              
              /* particle position */
              Xp  = (double)cur_part->pos[X];
              Yp  = (double)cur_part->pos[Y];
              Zp  = (double)cur_part->pos[Z];
              
              /* put particle into halo rest frame */
              dX  = (Xp - Xc);
              dY  = (Yp - Yc);
              dZ  = (Zp - Zc);
              
              /* take care of periodic boundary conditions */
              if(dX >  0.5) dX -= 1.0;
              if(dY >  0.5) dY -= 1.0;
              if(dZ >  0.5) dZ -= 1.0;
              if(dX < -0.5) dX += 1.0;
              if(dY < -0.5) dY += 1.0;
              if(dZ < -0.5) dZ += 1.0;
              
              /* distance of current particle (shouldn't sqrt(), but whatever...) */
              cur_dist = sqrt(pow2(dX) + pow2(dY) + pow2(dZ));
              
              /* jump to next particle */
              jpart++;
              
            } /* while(cur_dist < R_sphere && jpart < halos[i].npart) */
          
          /* volume of sphere [0, cur_rad] */
          V_sphere = 4.*PI/3. * pow3(R_sphere);
          
          r[ibin]       = R_sphere;
          ovdens[ibin]  = M_sphere/V_sphere*rho_fac/global.rho_b;        
        } /* ibin */
      
      /* we have ovdens(r) available now... */
      R_edge = get_haloedge(r, ovdens, nbins, 3);        
    } /* icall==0 */
  else
    {
      R_edge = 2*halo->R_vir;
    }
#endif /* AHFprofilerise */
  
  /* (re-)set values... */
  jpart              = 0;
  npart_sphere       = 0;
  M_sphere           = 0.0;
  R_sphere           = -1.0;
  ovdens_sphere      = 2*global.ovlim;
  
  
  while (   jpart < halo->npart && ovdens_sphere >= global.ovlim 
#ifdef AHFprofilerise
         && R_sphere <= R_edge
#endif
         )
    {
      /* access particle */
      cur_part = global.fst_part + halo->ipart[jpart];
      
#ifdef MULTIMASS
      weight = (double) cur_part->weight;
#else
      weight = (double) 1.0;
#endif
      /* mass in current sphere */
      M_sphere += weight;
      npart_sphere++;
      
      /* particle position */
      Xp  = (double)cur_part->pos[X];
      Yp  = (double)cur_part->pos[Y];
      Zp  = (double)cur_part->pos[Z];
      
      /* put particle into halo rest frame */
      dX  = fabs(Xp - Xc);
      dY  = fabs(Yp - Yc);
      dZ  = fabs(Zp - Zc);
      
      /* take care of periodic boundary conditions */
      if(dX >  0.5) dX -= 1.0;
      if(dY >  0.5) dY -= 1.0;
      if(dZ >  0.5) dZ -= 1.0;
      
      /* radius of current sphere */
      R_sphere = sqrt(pow2(dX) + pow2(dY) + pow2(dZ));
      
      /* volume of current sphere */
      V_sphere = 4.*PI/3. * pow3(R_sphere);
      
      /* overdensity in current sphere */
      ovdens_sphere = M_sphere/V_sphere * rho_fac/global.rho_b;
      
      /* move to next particle */
      jpart++;
    }
  
#ifdef TRACKER_DEBRIS
  /* before chopping the outside particles by realloc, store them */
  nopart = halo->npart - npart_sphere;
  if (icall == 0) 
    {
      halo->nopart =  nopart;
      if (halo->iopart) free(halo->iopart);
      halo->iopart = NULL;
      halo->iopart = (long unsigned *) realloc(halo->iopart, halo->nopart*sizeof(long unsigned));
      
      for (jpart=npart_sphere; jpart<halo->npart; jpart++)
        halo->iopart[jpart-npart_sphere] = halo->ipart[jpart];
    }
  else 
    {
      /* just append to already found outliers */
      halo->iopart = (long unsigned *) realloc(halo->iopart, (halo->nopart+nopart)*sizeof(long unsigned));
      
      for (jpart=0; jpart<nopart; jpart++)
        halo->iopart[halo->nopart + jpart] = halo->ipart[npart_sphere+jpart];
      
      halo->nopart += nopart;
    }
  
  fprintf(stderr,"n_outside, n_inside: %ld %ld\n", halo->nopart, npart_sphere);
  
  /*  for (jpart=0; jpart<halo->nopart; jpart++)
   {
   fprintf(stderr,"      outlier ID #%ld: %d\n", jpart+1, (global.fst_part+halo->iopart[jpart])->ID);
   }
   */
#endif
  
  /* update halo properties (a realloc() nicely chops the ipart[] array...) */
  halo->npart  = npart_sphere;
  halo->ipart  = (long unsigned *) realloc(halo->ipart, halo->npart*sizeof(long unsigned));
  halo->M_vir  = M_sphere;
  halo->R_vir  = R_sphere;
  halo->ovdens = ovdens_sphere;
  
  /******************************************************************************************************************/
  /******************************************************************************************************************/
  /******************************************************************************************************************/
  
  
#ifdef VERBOSE
  fprintf(io.logfile,"%12ld (Rvir=%g kpc/h)\n",halo->npart, halo->R_vir*x_fac*1000.);
  fflush(io.logfile);
#endif
  
#ifdef AHFprofilerise
  if(icall == 0)
    {
      free(ovdens);
      free(r);
    }
#endif
}


/*
 ************************************************************
 ************************************************************
 * Calculate profiles of halos
 */
int HaloProfiles(HALO *halo) 
{
  long unsigned ihalo, npart, n_prev, jpart;
  double  Volume, V_prev, dV, cur_rad, lcur_rad, rad_prev, cur_dist, rad_mid;
  double  M_sphere, M_prev, dM, weight,M_vir,r_vir;
  partptr cur_part, most_bound;
  gasptr  cur_gas;
  double  Xc, Yc, Zc, VXc, VYc, VZc;
  double  Xp, Yp, Zp, VXp, VYp, VZp;
  double  dX, dY, dZ, dVX, dVY, dVZ;
  double  R_edge,Lx,Ly,Lz,absAngMom,sig_v;
  double  itensor[3][3];
  double  a11,a22,a33,a12,a13,a23,axis1,axis2,axis3;
  double  pre_dist,I_mid,I_now,I_prev,Phi0,Phi,Epot,Ekin,Tpart,Upart,Epart,Emin;
  double  dist_min, dist_max, ldist_min, ldist_max, dr, ldr;
  int     binpart, nbins, ibin;
  double  dummy;
  double *Vcirc2, *dens_r2, *ovdens;
  double  r2, R_max, V_max, x_max, y_max;
  double  CoM[NDIM];
  double  tmp1, tmp2, tmp3;
  
#ifdef GAS_PARTICLES
  long unsigned ngas;
  double        M_gas;
  double        Lx_gas, Ly_gas, Lz_gas;
  double        Ekin_gas, Epot_gas, v_therm;
  double        a11_gas,a22_gas,a33_gas,a12_gas,a13_gas,a23_gas;
  
  long unsigned nstar;
  double        M_star;
  double        Lx_star, Ly_star, Lz_star;
  double        Ekin_star, Epot_star;
  double        a11_star,a22_star,a33_star,a12_star,a13_star,a23_star;
#endif
  
  /* only consider decent halos 
   * NOTE: npart may have changed since the last if-statement in ConstructHalos() */
#ifndef TRACKER
  if(halo->npart < AHF_MINPART)
#else
    if(halo->npart < TRK_MINPART)
#endif
      return (TRUE);
  
#ifdef VERBOSE
  fprintf(io.logfile,"    HaloProfiles:     npart=%12ld\n",halo->npart);
  fflush(io.logfile);
#endif
  
  /* how many bins from where to where for density profile? */
  binning_parameter(*halo, &nbins, &dist_min, &dist_max);	// halos[i]
  
  /* logarithmical binning from dist_min to dist_max */
  ldist_min = log10(dist_min);
  ldist_max = log10(dist_max);
  ldr       = (ldist_max-ldist_min)/(double)nbins;
  
  /* create profile arrays */
  c_profile(halo, nbins);		// &halos[i]
  Vcirc2  = (double *) calloc(nbins, sizeof(double));
  dens_r2 = (double *) calloc(nbins, sizeof(double));
  ovdens  = (double *) calloc(nbins, sizeof(double));
  
  /* reset values */
  Phi0       = halo->Phi0;
  Phi        = 0.0;
  pre_dist   = 0.0;
  I_prev     = 0.0;
  
  npart      = 0;
  n_prev     = 0;
  M_prev     = 0;
  V_prev     = 0.0;
  rad_prev   = dist_min;
  
  /* reset (cumulative) values based upon particles [0,cur_rad] */
  M_sphere    = 0.0;
  VXc         = 0.0;
  VYc         = 0.0;
  VZc         = 0.0;
  a11         = 0.0;
  a22         = 0.0;
  a33         = 0.0;
  a12         = 0.0;
  a13         = 0.0;
  a23         = 0.0;
  sig_v       = 0.0;
  Lx          = 0.0;
  Ly          = 0.0;
  Lz          = 0.0;
  CoM[X]      = 0.0;
  CoM[Y]      = 0.0;
  CoM[Z]      = 0.0;
  Epot        = 0.0;
  Ekin        = 0.0;
  Emin        = 1e30;
  most_bound  = NULL;
#ifdef GAS_PARTICLES
  ngas        = 0;
  M_gas       = 0.0;
  Lx_gas      = 0.0;
  Ly_gas      = 0.0;
  Lz_gas      = 0.0;
  Ekin_gas    = 0.0;
  Epot_gas    = 0.0;
  a11_gas     = 0.0;
  a22_gas     = 0.0;
  a33_gas     = 0.0;
  a12_gas     = 0.0;
  a13_gas     = 0.0;
  a23_gas     = 0.0;
  
  nstar       = 0;
  M_star      = 0.0;
  Lx_star     = 0.0;
  Ly_star     = 0.0;
  Lz_star     = 0.0;
  Ekin_star   = 0.0;
  Epot_star   = 0.0;
  a11_star    = 0.0;
  a22_star    = 0.0;
  a33_star    = 0.0;
  a12_star    = 0.0;
  a13_star    = 0.0;
  a23_star    = 0.0;
#endif
  
  /* halo position in AMIGA units */
  Xc = halo->pos.x;
  Yc = halo->pos.y;
  Zc = halo->pos.z;
  
  /* first particle */
  jpart    = 0;
  cur_dist = -1.0;
  
  /* loop over all bins */
  for(ibin=0; ibin<nbins; ibin++)
    {
      /* get current outer radius using logarithmic radial bins */
      lcur_rad = ldist_min + ((double)ibin+1) * ldr;
      cur_rad  = pow(10., lcur_rad);
      
      /* reset (differential) values calculated within each shell [cur_rad-1, cur_rad] */
      /*sig_v = 0.0;
       Lx    = 0.0;
       Ly    = 0.0;
       Lz    = 0.0;*/
      
      /* this heavily assumes that particles are ordered distancewise... */
      while(cur_dist < cur_rad && jpart < halo->npart)
        { 
          /* access particle */
          cur_part = global.fst_part + halo->ipart[jpart];
          
          /* increment number of particles counter */
          npart++;
          
          /*----------------------
           * calculate properties
           *----------------------*/
#ifdef MULTIMASS
          weight = (double)cur_part->weight;
#else
          weight = (double)1.0;
#endif
          
          /* accumulate mass */
          M_sphere += weight;
          
          /* particle position */
          Xp  = (double)cur_part->pos[X];
          Yp  = (double)cur_part->pos[Y];
          Zp  = (double)cur_part->pos[Z];
          
          /* put particle into halo rest frame */
          dX  = (Xp - Xc);
          dY  = (Yp - Yc);
          dZ  = (Zp - Zc);
          
          /* take care of periodic boundary conditions */
          if(dX >  0.5) dX -= 1.0;
          if(dY >  0.5) dY -= 1.0;
          if(dZ >  0.5) dZ -= 1.0;
          if(dX < -0.5) dX += 1.0;
          if(dY < -0.5) dY += 1.0;
          if(dZ < -0.5) dZ += 1.0;
          
          /* distance of current particle (shouldn't sqrt(), but whatever...) */
          cur_dist = sqrt(pow2(dX) + pow2(dY) + pow2(dZ));
          
          /* centre of mass with correct boundaries... */
          CoM[X] += weight*(Xc+dX);
          CoM[Y] += weight*(Yc+dY);
          CoM[Z] += weight*(Zc+dZ);
          
#ifdef AHFreducedinertiatensor
          /* reduced inertia tensor of all particles within sphere(!) */
          a11 += weight*dX*dX/pow2(cur_dist);
          a22 += weight*dY*dY/pow2(cur_dist);
          a33 += weight*dZ*dZ/pow2(cur_dist);
          a12 += weight*dX*dY/pow2(cur_dist);
          a13 += weight*dX*dZ/pow2(cur_dist);
          a23 += weight*dY*dZ/pow2(cur_dist);
#else
          /* inertia tensor of all particles within sphere(!) */
          a11 += weight*dX*dX;
          a22 += weight*dY*dY;
          a33 += weight*dZ*dZ;
          a12 += weight*dX*dY;
          a13 += weight*dX*dZ;
          a23 += weight*dY*dZ;
#endif
          /* mean halo velocity of particles in sphere [0,cur_rad] */
          VXc += weight*cur_part->mom[X];
          VYc += weight*cur_part->mom[Y];
          VZc += weight*cur_part->mom[Z];
          
          /* particle velocity in halo rest frame (Hubble correction not needed for L as r x r = 0) */
          VXp = (double)cur_part->mom[X];
          VYp = (double)cur_part->mom[Y];
          VZp = (double)cur_part->mom[Z];
          dVX = (VXp - VXc/M_sphere);
          dVY = (VYp - VYc/M_sphere);
          dVZ = (VZp - VZc/M_sphere);
          
          /* angular momentum of particles within sphere [0, cur_rad] */
          Lx += weight*(dY*dVZ-dZ*dVY);
          Ly += weight*(dZ*dVX-dX*dVZ);
          Lz += weight*(dX*dVY-dY*dVX);
          
          /* add Hubble flow to velocities */
          dVX   += Hubble * dX * r_fac/v_fac;
          dVY   += Hubble * dY * r_fac/v_fac;
          dVZ   += Hubble * dZ * r_fac/v_fac;
          
          /* kinetic energy of particle */
          Tpart  = weight * (pow2(dVX)+pow2(dVY)+pow2(dVZ));
          
          
          
          /* velocity dispersion and kinetic energy of all particles within sphere [0, cur_rad] */
          sig_v += Tpart;
          Ekin  += Tpart;
          
          
          /* get potential energy of particle */
          if(cur_dist > ZERO)
            {
              I_now = M_sphere/pow2(cur_dist);
              I_mid = (I_now + I_prev)/2.;
              Phi  += I_mid * (cur_dist - pre_dist);
              
              Upart = (Phi-Phi0)*weight;
            }
          else
            {
              /* simply use old Phi value... */
              Upart = (Phi-Phi0)*weight;
            }
          
          /* take care of potential energy integration... */
          I_prev   = I_now;
          pre_dist = cur_dist;
          
          /* potential energy of all particles within sphere [0, cur_rad] */
          Epot    += Upart;
          
          /* total energy of current particle */
          Epart = (0.5*Tpart + Upart);
          
          /* remember most bound particle */
          if(Epart < Emin)
            {
              Emin       = Epart;
              most_bound = cur_part;
            }
          
          
#ifdef GAS_PARTICLES
#if (defined NEWSTARTRUN)
          if(isgreaterequal(cur_part->u, PGAS))
#else
          if(cur_part-global.fst_part-global.offset_gas < global.no_gas && cur_part-global.fst_part-global.offset_gas > 0)
#endif
              {
                /* calculate integral properties of gas alone */
                ngas++;
                M_gas   += weight;
                Lx_gas  += weight*(dY*dVZ-dZ*dVY);
                Ly_gas  += weight*(dZ*dVX-dX*dVZ);
                Lz_gas  += weight*(dX*dVY-dY*dVX);
#ifdef AHFreducedinertiatensor
                a11_gas += weight*dX*dX/pow2(cur_dist);
                a22_gas += weight*dY*dY/pow2(cur_dist);
                a33_gas += weight*dZ*dZ/pow2(cur_dist);
                a12_gas += weight*dX*dY/pow2(cur_dist);
                a13_gas += weight*dX*dZ/pow2(cur_dist);
                a23_gas += weight*dY*dZ/pow2(cur_dist);
#else
                a11_gas += weight*dX*dX;
                a22_gas += weight*dY*dY;
                a33_gas += weight*dZ*dZ;
                a12_gas += weight*dX*dY;
                a13_gas += weight*dX*dZ;
                a23_gas += weight*dY*dZ;
#endif
                
#if (!defined NEWSTARTRUN)
                /* current == gas particle => add thermal velocity */
                cur_gas = global.fst_gas + (cur_part-global.fst_part-global.offset_gas);
#endif
                
                /* potential energy of gas alone */
                Epot_gas += Upart;
                
                /* calculate kinetic+internal gas energy */
                Ekin_gas += Tpart;
                
                /* u = 3/2  kT/m  = 1/2  v_therm^2 ?! */
#if (defined NEWSTARTRUN)
                v_therm = weight * 2*cur_part->u/pow2(v_fac);
#else
                v_therm = weight * 2*cur_gas->u /pow2(v_fac);
#endif
                
                /* account for thermal energy */
                sig_v    += v_therm;
                Ekin     += v_therm;
                Ekin_gas += v_therm;
                
              }
          
#if (defined NEWSTARTRUN)
          if(fabs(cur_part->u-PSTAR) < ZERO)
#else
          if(cur_part-global.fst_part-global.offset_stars < global.no_stars && cur_part-global.fst_part-global.offset_stars > 0)
#endif
              {
                /* calculate integral properties of stars alone */
                nstar++;
                M_star  += weight;
                Lx_star += weight*(dY*dVZ-dZ*dVY);
                Ly_star += weight*(dZ*dVX-dX*dVZ);
                Lz_star += weight*(dX*dVY-dY*dVX);
#ifdef AHFreducedinertiatensor
                a11_star += weight*dX*dX/pow2(cur_dist);
                a22_star += weight*dY*dY/pow2(cur_dist);
                a33_star += weight*dZ*dZ/pow2(cur_dist);
                a12_star += weight*dX*dY/pow2(cur_dist);
                a13_star += weight*dX*dZ/pow2(cur_dist);
                a23_star += weight*dY*dZ/pow2(cur_dist);
#else
                a11_star += weight*dX*dX;
                a22_star += weight*dY*dY;
                a33_star += weight*dZ*dZ;
                a12_star += weight*dX*dY;
                a13_star += weight*dX*dZ;
                a23_star += weight*dY*dZ;
#endif
                /* potential energy of gas alone */
                Epot_star += Upart;
                
                /* calculate kinetic+internal gas energy */
                Ekin_star += Tpart;              
              }
          
#endif /* GAS_PARTICLES */
          
          /* jump to next particle */
          jpart++;
          
        } /* while(cur_dist < cur_rad && jpart < halo->npart) */
      
      /* volume of sphere [0, cur_rad] */
      Volume = 4.*PI/3. * pow3(cur_rad);
      
      /* differential values */
      dM = M_sphere - M_prev;  /* mass   in shell [prev_rad, cur_rad] */
      dV = Volume   - V_prev;  /* volume of shell [prev_rad, cur_rad] */ 
      
      /* get eigenavalues of inertia tensor (all part's within sphere!) */
      itensor[0][0] = a11;
      itensor[1][1] = a22;
      itensor[2][2] = a33;
      itensor[0][1] = a12;
      itensor[1][0] = a12;
      itensor[0][2] = a13;
      itensor[2][0] = a13;
      itensor[1][2] = a23;
      itensor[2][1] = a23;
      get_axes(itensor, &axis1, &axis2, &axis3);
      
      rad_mid = (cur_rad+rad_prev)/2.;
      halo->prof.r[ibin]       = cur_rad;
      halo->prof.npart[ibin]   = npart;
      halo->prof.nvpart[ibin]  = M_sphere;
#ifdef GAS_PARTICLES
      halo->prof.nvpart_gas[ibin]  = M_gas;
      halo->prof.nvpart_star[ibin] = M_star;
#endif
      halo->prof.ovdens[ibin]  = M_sphere / Volume;
      halo->prof.dens[ibin]    = dM       / dV;
      halo->prof.v2_circ[ibin] = M_sphere / cur_rad;
      halo->prof.sig_v[ibin]   = sqrt(sig_v/M_sphere);
      halo->prof.Ekin[ibin]    = 0.5*Ekin;
      halo->prof.Epot[ibin]    = 0.5*Epot;
      halo->prof.Lx[ibin]      = Lx;
      halo->prof.Ly[ibin]      = Ly;
      halo->prof.Lz[ibin]      = Lz;
      
      halo->prof.E1x[ibin]     = itensor[0][0];
      halo->prof.E1y[ibin]     = itensor[1][0];
      halo->prof.E1z[ibin]     = itensor[2][0];
      halo->prof.E2x[ibin]     = itensor[0][1];
      halo->prof.E2y[ibin]     = itensor[1][1];
      halo->prof.E2z[ibin]     = itensor[2][1];
      halo->prof.E3x[ibin]     = itensor[0][2];
      halo->prof.E3y[ibin]     = itensor[1][2];
      halo->prof.E3z[ibin]     = itensor[2][2];
#ifdef AHFabsaxes
      halo->prof.axis1[ibin]   = sqrt(axis1/M_sphere);
      halo->prof.axis2[ibin]   = sqrt(axis2/M_sphere);
      halo->prof.axis3[ibin]   = sqrt(axis3/M_sphere);
#else
      halo->prof.axis1[ibin]   = 1.0;
      halo->prof.axis2[ibin]   = sqrt(axis2/axis1);
      halo->prof.axis3[ibin]   = sqrt(axis3/axis1);
#endif
      
      /* store Vcirc and dens_r2 for later determination of Vmax, Rmax, and r2 */
      Vcirc2[ibin]  = halo->prof.v2_circ[ibin];
      dens_r2[ibin] = halo->prof.dens[ibin]*pow2(rad_mid);
      
      /* also store ovdens as get_haloedge() calls find_min() that overrides ovdens[] */
      ovdens[ibin]  = halo->prof.ovdens[ibin];
      
      /* store old values */
      M_prev   = M_sphere;
      V_prev   = Volume;
      n_prev   = npart;
      rad_prev = cur_rad;
      
    } /* ibin */
  
  /* centre of mass */
  CoM[X] /= M_sphere;
  CoM[Y] /= M_sphere;
  CoM[Z] /= M_sphere;
  CoM[X]  = fmod(CoM[X]+1.,1.);
  CoM[Y]  = fmod(CoM[Y]+1.,1.);
  CoM[Z]  = fmod(CoM[Z]+1.,1.);
  
  /* get physical outer radius: either upturn or virial radius */
  R_edge = get_haloedge(halo->prof.r, ovdens, nbins, 1);
  free(ovdens);
  
  /* get scale radius */
#ifdef AHFsplinefit
  find_max_spline (halo->prof.r, dens_r2, nbins, 3, &x_max, &y_max, 10000);    // use spline fit
#else
  find_max        (halo->prof.r, dens_r2, nbins, 3, &x_max, &y_max);
#endif
  r2 = x_max;
  free(dens_r2);
  
  /* get rotation curve properties */
#ifdef AHFsplinefit
  find_max_spline (halo->prof.r, Vcirc2, nbins, 1, &x_max, &y_max, 10000);    // use spline fit
#else
  find_max        (halo->prof.r, Vcirc2, nbins, 1, &x_max, &y_max);
#endif
  R_max = x_max;
  V_max = y_max;
  free(Vcirc2);
  
  /* store final (integral) halo properties */
  halo->R_edge   = R_edge;
  halo->R_max    = R_max;
  halo->V2_max   = V_max;
  halo->r2       = r2;
  
  if(most_bound != NULL)
    {
      dX = most_bound->pos[X]-Xc;
      dY = most_bound->pos[Y]-Yc;
      dZ = most_bound->pos[Z]-Zc;
      if(dX >  0.5) dX -= 1.0;
      if(dY >  0.5) dY -= 1.0;
      if(dZ >  0.5) dZ -= 1.0;
      if(dX < -0.5) dX += 1.0;
      if(dY < -0.5) dY += 1.0;
      if(dZ < -0.5) dZ += 1.0;      
      halo->mbp_offset = sqrt(pow2(dX) + pow2(dY) + pow2(dZ));
    }
  else
    halo->mbp_offset = -1.0;
  
  dX = CoM[X]-Xc;
  dY = CoM[Y]-Yc;
  dZ = CoM[Z]-Zc;
  if(dX >  0.5) dX -= 1.0;
  if(dY >  0.5) dY -= 1.0;
  if(dZ >  0.5) dZ -= 1.0;
  if(dX < -0.5) dX += 1.0;
  if(dY < -0.5) dY += 1.0;
  if(dZ < -0.5) dZ += 1.0;      
  halo->com_offset = sqrt(pow2(dX) + pow2(dY) + pow2(dZ));
  
  halo->M_vir    = M_sphere;
  halo->vel.x    = VXc/M_sphere;
  halo->vel.y    = VYc/M_sphere;
  halo->vel.z    = VZc/M_sphere;
  
  /* we assume that profile stores cumulative values... */
  halo->velDis   = halo->prof.sig_v[nbins-1];
  halo->Ekin     = halo->prof.Ekin[nbins-1];
  halo->Epot     = halo->prof.Epot[nbins-1];
  
  /* standard intertia tensor */
  halo->axis.x   = halo->prof.axis1[nbins-1];
  halo->axis.y   = halo->prof.axis2[nbins-1];
  halo->axis.z   = halo->prof.axis3[nbins-1];
  halo->E1.x	 = halo->prof.E1x[nbins-1];
  halo->E1.y	 = halo->prof.E1y[nbins-1];
  halo->E1.z	 = halo->prof.E1z[nbins-1];
  halo->E2.x	 = halo->prof.E2x[nbins-1];
  halo->E2.y	 = halo->prof.E2y[nbins-1];
  halo->E2.z	 = halo->prof.E2z[nbins-1];
  halo->E3.x	 = halo->prof.E3x[nbins-1];
  halo->E3.y	 = halo->prof.E3y[nbins-1];
  halo->E3.z	 = halo->prof.E3z[nbins-1];
  
  /* angular momentum */
  Lx = halo->prof.Lx[nbins-1];
  Ly = halo->prof.Ly[nbins-1];
  Lz = halo->prof.Lz[nbins-1];
  absAngMom = sqrt(pow2(Lx)+pow2(Ly)+pow2(Lz));
#ifdef AHFabsangmom
  halo->AngMom.x = Lx;
  halo->AngMom.y = Ly;
  halo->AngMom.z = Lz;
#else
  halo->AngMom.x = Lx/absAngMom;
  halo->AngMom.y = Ly/absAngMom;
  halo->AngMom.z = Lz/absAngMom;
#endif
  halo->lambda   = absAngMom/halo->M_vir/sqrt(2. * halo->M_vir * halo->R_vir);
  halo->lambda  *= v_fac*sqrt(r_fac/(Grav*m_fac));
  
  
  /* energy-based lambda: careful with numerics... */
  {
    tmp1 = sqrt(m_fac*halo->M_vir);          /* tmp1 = sqrt(M)         */
    tmp1 = pow5(tmp1);                       /* tmp1 = M^5/2           */
    tmp2 = halo->Ekin*m_fac*pow2(v_fac);     /* tmp2 = Ekin            */
    tmp3 = halo->Epot*m_fac*phi_fac;         /* tmp3 = Epot            */
    tmp2 = fabs(tmp2+tmp3);                  /* tmp2 = E               */
    tmp2 = sqrt(tmp2);                       /* tmp2 = sqrt(|E|)       */
    tmp1 = tmp2/tmp1;                        /* tmp1 = sqrt(|E|)/M^5/2 */
    tmp2 = m_fac*r_fac*v_fac*absAngMom;      /* tmp2 = L               */
    tmp3 = tmp1*tmp2;
    
    halo->lambdaE = tmp3/Grav;
  }
  
#ifdef GAS_PARTICLES
  halo->gas.npart = ngas;
  halo->gas.M_vir = M_gas;
  
  if(ngas > 0)
    {
      /* get spin parameter of gas particles (all gas part's within sphere!) */
      absAngMom           = sqrt(pow2(Lx_gas)+pow2(Ly_gas)+pow2(Lz_gas));
      halo->gas.AngMom.x  = Lx_gas/absAngMom;
      halo->gas.AngMom.y  = Ly_gas/absAngMom;
      halo->gas.AngMom.z  = Lz_gas/absAngMom;
      halo->gas.Ekin      = 0.5*Ekin_gas;
      halo->gas.Epot      = 0.5*Epot_gas;
      halo->gas.lambda    = absAngMom/M_gas/sqrt(2. * halo->M_vir * halo->R_vir);
      halo->gas.lambda   *= v_fac*sqrt(r_fac/(Grav*m_fac));
      /* energy-based lambda: careful with numerics... */
      {
        tmp1 = sqrt(m_fac*halo->M_vir);              /* tmp1 = sqrt(M_tot)             */
        tmp1 = pow3(tmp1);                           /* tmp1 = M_tot^3/2               */
        tmp2 = halo->gas.Ekin*m_fac*pow2(v_fac);     /* tmp2 = Ekin_gas                */
        tmp3 = halo->gas.Epot*m_fac*phi_fac;         /* tmp3 = Epot_gas                */
        tmp2 = fabs(tmp2+tmp3);                      /* tmp2 = E_gas                   */
        tmp2 = sqrt(tmp2);                           /* tmp2 = sqrt(|E_gas|)           */
        tmp1 = tmp2/tmp1;                            /* tmp1 = sqrt(|E_gas|)/M_tot^3/2 */
        tmp2 = m_fac*r_fac*v_fac*absAngMom;          /* tmp2 = L_gas                   */
        tmp2 = tmp2/(m_fac*halo->gas.M_vir);         /* tmp2 = L_gas/M_gas             */
        tmp3 = tmp1*tmp2;
        
        halo->gas.lambdaE = tmp3/Grav;
      }
      
      /* get shape of gas particles (all gas part's within sphere!) */
      itensor[0][0] = a11_gas;
      itensor[1][1] = a22_gas;
      itensor[2][2] = a33_gas;
      itensor[0][1] = a12_gas;
      itensor[1][0] = a12_gas;
      itensor[0][2] = a13_gas;
      itensor[2][0] = a13_gas;
      itensor[1][2] = a23_gas;
      itensor[2][1] = a23_gas;
      get_axes(itensor, &axis1, &axis2, &axis3);
      
      halo->gas.axis.x = 1.0;
      halo->gas.axis.y = sqrt(axis2/axis1);
      halo->gas.axis.z = sqrt(axis3/axis1);
      halo->gas.E1.x   = itensor[0][0];
      halo->gas.E1.y   = itensor[1][0];
      halo->gas.E1.z   = itensor[2][0];
      halo->gas.E2.x   = itensor[0][1];
      halo->gas.E2.y   = itensor[1][1];
      halo->gas.E2.z   = itensor[2][1];
      halo->gas.E3.x   = itensor[0][2];
      halo->gas.E3.y   = itensor[1][2];
      halo->gas.E3.z   = itensor[2][2]; 	
    }
  else
    {
      halo->gas.AngMom.x = 0.0;
      halo->gas.AngMom.y = 0.0;
      halo->gas.AngMom.z = 0.0;
      halo->gas.lambda   = 0.0;
      halo->gas.lambdaE  = 0.0;
      halo->gas.Ekin	 = 0.0;
      halo->gas.Epot	 = 0.0;
      halo->gas.axis.x   = 0.0;
      halo->gas.axis.y   = 0.0;
      halo->gas.axis.z   = 0.0;
      halo->gas.E1.x	 = 0.0;
      halo->gas.E1.y	 = 0.0;
      halo->gas.E1.z	 = 0.0;
      halo->gas.E2.x	 = 0.0;
      halo->gas.E2.y	 = 0.0;
      halo->gas.E2.z	 = 0.0;
      halo->gas.E3.x	 = 0.0;
      halo->gas.E3.y	 = 0.0;
      halo->gas.E3.z	 = 0.0;
    }
  
  halo->star.npart = nstar;
  halo->star.M_vir = M_star;
  
  if(nstar > 0)
    {
      /* get spin parameter of gas particles (all gas part's within sphere!) */
      absAngMom           = sqrt(pow2(Lx_star)+pow2(Ly_star)+pow2(Lz_star));
      halo->star.AngMom.x  = Lx_star/absAngMom;
      halo->star.AngMom.y  = Ly_star/absAngMom;
      halo->star.AngMom.z  = Lz_star/absAngMom;
      halo->star.Ekin      = 0.5*Ekin_star;
      halo->star.Epot      = 0.5*Epot_star;
      halo->star.lambda    = absAngMom/M_star/sqrt(2. * halo->M_vir * halo->R_vir);
      halo->star.lambda   *= v_fac*sqrt(r_fac/(Grav*m_fac));
      /* energy-based lambda: careful with numerics... */
      {
        tmp1 = sqrt(m_fac*halo->M_vir);              /* tmp1 = sqrt(M_tot)               */
        tmp1 = pow3(tmp1);                           /* tmp1 = M_tot^3/2                 */
        tmp2 = halo->star.Ekin*m_fac*pow2(v_fac);    /* tmp2 = Ekin_star                 */
        tmp3 = halo->star.Epot*m_fac*phi_fac;        /* tmp3 = Epot_star                 */
        tmp2 = fabs(tmp2+tmp3);                      /* tmp2 = E_star                    */
        tmp2 = sqrt(tmp2);                           /* tmp2 = sqrt(|E_star|)            */
        tmp1 = tmp2/tmp1;                            /* tmp1 = sqrt(|E_star|)/M_tot^3/2  */
        tmp2 = m_fac*r_fac*v_fac*absAngMom;          /* tmp2 = L_star                    */
        tmp2 = tmp2/(m_fac*halo->star.M_vir);        /* tmp2 = L_star/M_star             */
        tmp3 = tmp1*tmp2;
        
        halo->star.lambdaE = tmp3/Grav;
      }
      
      /* get shape of gas particles (all gas part's within sphere!) */
      itensor[0][0] = a11_star;
      itensor[1][1] = a22_star;
      itensor[2][2] = a33_star;
      itensor[0][1] = a12_star;
      itensor[1][0] = a12_star;
      itensor[0][2] = a13_star;
      itensor[2][0] = a13_star;
      itensor[1][2] = a23_star;
      itensor[2][1] = a23_star;
      get_axes(itensor, &axis1, &axis2, &axis3);
      
      halo->star.axis.x = 1.0;
      halo->star.axis.y = sqrt(axis2/axis1);
      halo->star.axis.z = sqrt(axis3/axis1);
      halo->star.E1.x   = itensor[0][0];
      halo->star.E1.y   = itensor[1][0];
      halo->star.E1.z   = itensor[2][0];
      halo->star.E2.x   = itensor[0][1];
      halo->star.E2.y   = itensor[1][1];
      halo->star.E2.z   = itensor[2][1];
      halo->star.E3.x   = itensor[0][2];
      halo->star.E3.y   = itensor[1][2];
      halo->star.E3.z   = itensor[2][2]; 	
    }
  else
    {
      halo->star.AngMom.x = 0.0;
      halo->star.AngMom.y = 0.0;
      halo->star.AngMom.z = 0.0;
      halo->star.lambda   = 0.0;
      halo->star.lambdaE  = 0.0;
      halo->star.Ekin	 = 0.0;
      halo->star.Epot	 = 0.0;
      halo->star.axis.x   = 0.0;
      halo->star.axis.y   = 0.0;
      halo->star.axis.z   = 0.0;
      halo->star.E1.x	 = 0.0;
      halo->star.E1.y	 = 0.0;
      halo->star.E1.z	 = 0.0;
      halo->star.E2.x	 = 0.0;
      halo->star.E2.y	 = 0.0;
      halo->star.E2.z	 = 0.0;
      halo->star.E3.x	 = 0.0;
      halo->star.E3.y	 = 0.0;
      halo->star.E3.z	 = 0.0;
    }  
#endif
  
  return(TRUE);
}

#ifdef AHFphspdens
int
HaloProfilesPhaseSpace(HALO *halo)
{
  double halo_pos[3];
  double halo_vel[3];
  int nbins, ibin;
  int npart_sp, npart_sh;
  double binpart;
  double lcur_rad, cur_rad, cur_dist;
  double ldist_min, ldist_max, ldr;
  double dist_min, dist_max;
  partptr cur_part, fst_part, lst_part;
  double x, y, z, vx, vy, vz;
  double r, theta, phi, vr, vtheta, vphi;
  double tmp;
  double mean_vx_sp, mean_vy_sp, mean_vz_sp;
  double mean_vr_sp, mean_vtheta_sp, mean_vphi_sp;
  double sum_vx_sp, sum_vy_sp, sum_vz_sp;
  double sum_vr_sp, sum_vtheta_sp, sum_vphi_sp;
  double mean_vx_sh, mean_vy_sh, mean_vz_sh;
  double mean_vr_sh, mean_vtheta_sh, mean_vphi_sh;
  double sigma2_vx_sh, sigma2_vy_sh, sigma2_vz_sh;
  double sigma2_vr_sh, sigma2_vtheta_sh, sigma2_vphi_sh;
  long jpart, fst_jpart, lst_jpart;
  
  /* Don't do it for haloes too small */
#ifndef TRACKER
  if (halo->npart < AHF_MINPART)
#else
    if (halo->npart < TRK_MINPART)
#endif
      return TRUE;
  
  /* Get the halo restframe */
  halo_pos[0] = halo->pos.x;
  halo_pos[1] = halo->pos.y;
  halo_pos[2] = halo->pos.z;
  halo_vel[0] = halo->vel.x;
  halo_vel[1] = halo->vel.y;
  halo_vel[2] = halo->vel.z;
  
  
  /* how many bins from where to where for density profile? */
  binning_parameter(*halo, &nbins, &dist_min, &dist_max);
  
  /* logarithmical binning from dist_min to dist_max */
  ldist_min = log10(dist_min);
  ldist_max = log10(dist_max);
  ldr       = (ldist_max-ldist_min)/(double)nbins;
  
  /* Reset all quantities in the sphere */
  cur_dist = -1.0;
  jpart    = 0;
  npart_sp = 0;
  mean_vx_sp = mean_vy_sp = mean_vz_sp = 0.0;
  mean_vr_sp = mean_vtheta_sp = mean_vphi_sp = 0.0;
  sum_vx_sp = sum_vy_sp = sum_vz_sp = 0.0;
  sum_vr_sp = sum_vtheta_sp = sum_vphi_sp = 0.0;
  
  /* Loop over all bins */
  for (ibin=0; ibin<nbins; ibin++) {
    /* Current outer radius */
    lcur_rad = ldist_min + ((double)ibin+1) * ldr;
    cur_rad = pow(10., lcur_rad);
    
    /* Figure out first and last particle */
    cur_part = global.fst_part + halo->ipart[jpart];
    fst_part = lst_part = cur_part;
    fst_jpart = lst_jpart = jpart;
    npart_sh = 0;
    while (jpart < halo->npart && cur_dist < cur_rad) {
      /* Access particle*/
      cur_part = global.fst_part + halo->ipart[jpart];
      
      /* Get position */
      x = ((double)cur_part->pos[0])-halo_pos[0];
      y = ((double)cur_part->pos[1])-halo_pos[1];
      z = ((double)cur_part->pos[2])-halo_pos[2];
      
      /* Correct for periodic boundary */
      if (x<-0.5) x += 1.0;
      if (y<-0.5) y += 1.0;
      if (z<-0.5) z += 1.0;
      if (x>0.5) x -= 1.0;
      if (y>0.5) y -= 1.0;
      if (z>0.5) z -= 1.0;
      
      /* Get the distance from the center */
      cur_dist = sqrt(x*x + y*y + z*z);
      
      /* Check whether that particle counts */
      if (cur_dist < cur_rad) {
        lst_part = cur_part;
        lst_jpart = jpart;
        jpart++;
        npart_sp++;
        npart_sh++;
      }
    }
    
    /* Loop over all particles in the bin to get mean velocities */
    mean_vx_sh = mean_vy_sh = mean_vz_sh = 0.0;
    mean_vr_sh = mean_vtheta_sh = mean_vphi_sh = 0.0;
    if (npart_sh > 0) {
      for (jpart=fst_jpart;
           jpart<=lst_jpart;
           jpart++) {
        /* Access particle */
        cur_part = global.fst_part + halo->ipart[jpart];
        
        /* Particle position in halo rest frame */
        x = ((double)cur_part->pos[0])-halo_pos[0];
        y = ((double)cur_part->pos[1])-halo_pos[1];
        z = ((double)cur_part->pos[2])-halo_pos[2];
        
        /* Correct for periodic boundary */
        if (x<-0.5) x += 1.0;
        if (y<-0.5) y += 1.0;
        if (z<-0.5) z += 1.0;
        if (x>0.5) x -= 1.0;
        if (y>0.5) y -= 1.0;
        if (z>0.5) z -= 1.0;
        
        /* Used for the Hubble flow correction */
        tmp = Hubble * r_fac/v_fac;
        
        /* Particle velocity in halo rest frame w/ Hubble */
        vx = ((double)cur_part->mom[0]) - halo_vel[0] + x*tmp;
        vy = ((double)cur_part->mom[1]) - halo_vel[1] + y*tmp;
        vz = ((double)cur_part->mom[2]) - halo_vel[2] + z*tmp;
        
        /* Put positions to spherical coordinates */
        r = sqrt(x*x + y*y + z*z);
        theta = acos(z/r);
        phi = atan2(y,x);
        
        /* Put velocities to spherical coordinates */
        vr = 1./r*(x*vx + y*vy + z*vz);
        vtheta = -1.0/sqrt(1-z*z/(r*r)) * (vz/r - z*vr/(r*r));
        vphi = 1./(1 + y*y/(x*x)) * (vy/x - y*vx/(x*x));
        
        /* Keep track of the mean velocities in the shell
         * Here it is only the sum, division by npart takes place
         * after the loop */
        mean_vx_sh += vx;
        mean_vy_sh += vy;
        mean_vz_sh += vz;
        mean_vr_sh += vr;
        mean_vtheta_sh += vtheta;
        mean_vphi_sh += vphi;
        
        /* Keep track of the sum of the velocities in the sphere */
        sum_vx_sp += vx;
        sum_vy_sp += vy;
        sum_vz_sp += vz;
        sum_vr_sp += vr;
        sum_vtheta_sp += vtheta;
        sum_vphi_sp += vphi;
      }
      /* Calculate real mean velocity in the shell */
      tmp = 1./npart_sh;
      mean_vx_sh *= tmp;
      mean_vy_sh *= tmp;
      mean_vz_sh *= tmp;
      mean_vr_sh *= tmp;
      mean_vtheta_sh *= tmp;
      mean_vphi_sh *= tmp;
    } /* if (npart_sh > 0)*/
    
    /* Calculate the real mean velocity within the sphere */
    tmp = 1./npart_sp;
    mean_vx_sp = sum_vx_sp*tmp;
    mean_vy_sp = sum_vy_sp*tmp;
    mean_vz_sp = sum_vz_sp*tmp;
    mean_vr_sp = sum_vr_sp*tmp;
    mean_vtheta_sp = sum_vtheta_sp*tmp;
    mean_vphi_sp = sum_vphi_sp*tmp;
    
    /* Loop again over all particles in the bin to get dispersion */
    sigma2_vx_sh = sigma2_vy_sh = sigma2_vz_sh = 0.0;
    sigma2_vr_sh = sigma2_vtheta_sh = sigma2_vphi_sh = 0.0;
    if (npart_sh > 0) {
      for (jpart=fst_jpart;
           jpart<=lst_jpart;
           jpart++) {
        /* Access particle */
        cur_part = global.fst_part + halo->ipart[jpart];
        
        /* Particle position in halo rest frame */
        x = ((double)cur_part->pos[0])-halo_pos[0];
        y = ((double)cur_part->pos[1])-halo_pos[1];
        z = ((double)cur_part->pos[2])-halo_pos[2];
        
        /* Correct for periodic boundary */
        if (x<-0.5) x += 1.0;
        if (y<-0.5) y += 1.0;
        if (z<-0.5) z += 1.0;
        if (x>0.5) x -= 1.0;
        if (y>0.5) y -= 1.0;
        if (z>0.5) z -= 1.0;
        
        /* Used for the Hubble flow correction */
        tmp = Hubble * r_fac/v_fac;
        
        /* Particle velocity in halo rest frame w/ Hubble */
        vx = ((double)cur_part->mom[0]) - halo_vel[0] + x*tmp;
        vy = ((double)cur_part->mom[1]) - halo_vel[1] + y*tmp;
        vz = ((double)cur_part->mom[2]) - halo_vel[2] + z*tmp;	
        
        /* Put positions to spherical coordinates */
        r = sqrt(x*x + y*y + z*z);
        theta = acos(z/r);
        phi = atan2(y,x);
        
        /* Put velocities to spherical coordinates */
        vr = 1./r*(x*vx + y*vy + z*vz);
        vtheta = -1.0/sqrt(1-z*z/(r*r)) * (vz/r - z*vr/(r*r));
        vphi = 1./(1 + y*y/(x*x)) * (vy/x - y*vx/(x*x));
        
        /* Keep track of shell dispersions*(npart_sh-1) */
        sigma2_vx_sh += pow2(vx - mean_vx_sh);
        sigma2_vy_sh += pow2(vy - mean_vy_sh);
        sigma2_vz_sh += pow2(vz - mean_vz_sh);
        sigma2_vr_sh += pow2(vr - mean_vr_sh);
        sigma2_vtheta_sh += pow2(vtheta - mean_vtheta_sh);
        sigma2_vphi_sh += pow2(vphi - mean_vphi_sh);
      }
      /* Get the real dispersions in the shell */
      tmp = 1./(npart_sh - 1);
      sigma2_vx_sh *= tmp;
      sigma2_vy_sh *= tmp;
      sigma2_vz_sh *= tmp;
      sigma2_vr_sh *= tmp;
      sigma2_vtheta_sh *= tmp;
      sigma2_vphi_sh *= tmp;
    } /* if (npart > 0) */
    
    /* Plug the quantities into the profile array */
    halo->prof.sigma2_vx_sh[ibin] = sigma2_vx_sh;
    halo->prof.sigma2_vy_sh[ibin] = sigma2_vy_sh;
    halo->prof.sigma2_vz_sh[ibin] = sigma2_vz_sh;
    halo->prof.sigma2_vr_sh[ibin] = sigma2_vr_sh;
    halo->prof.sigma2_vtheta_sh[ibin] = sigma2_vtheta_sh;
    halo->prof.sigma2_vphi_sh[ibin] = sigma2_vphi_sh;
#ifdef AHFmeanvelocities
    halo->prof.mean_vx_sh[ibin] = mean_vx_sh;
    halo->prof.mean_vy_sh[ibin] = mean_vy_sh;
    halo->prof.mean_vz_sh[ibin] = mean_vz_sh;
    halo->prof.mean_vr_sh[ibin] = mean_vr_sh;
    halo->prof.mean_vtheta_sh[ibin] = mean_vtheta_sh;
    halo->prof.mean_vphi_sh[ibin] = mean_vphi_sh;
    halo->prof.mean_vx_sp[ibin] = mean_vx_sp;
    halo->prof.mean_vy_sp[ibin] = mean_vy_sp;
    halo->prof.mean_vz_sp[ibin] = mean_vz_sp;
    halo->prof.mean_vr_sp[ibin] = mean_vr_sp;
    halo->prof.mean_vtheta_sp[ibin] = mean_vtheta_sp;
    halo->prof.mean_vphi_sp[ibin] = mean_vphi_sp;
#endif
    
  } /* End of for-loop over the bins */
  
  return TRUE;
}
#endif /* AHFphspdens */


/*==============================================================================
 *  check Rho[] for rise and return position of rise
 *  in addition, check if Rho[] declines like r[]^AHF_SLOPE out to rising point
 *==============================================================================*/
double get_haloedge(double *r, double *Rho, int iRadOut, int ismooth)
{
  int    ir, iRadIn, iRadEnd, nrise;
  int    irm1, irp1, irp2;
  double gradlf, gradrt, R_edge;
  double Rholeft, Rhoright, rleft, rright, Rhomin, rmin;
  
  /* is there a well pronounced minimum in the density profile? */
#ifdef AHFsplinefit
  find_min_spline (r, Rho, iRadOut, ismooth, &rmin, &Rhomin, 10000);
#else
  find_min        (r, Rho, iRadOut, ismooth, &rmin, &Rhomin);
#endif
  
  if(rmin < r[iRadOut-1])
    {
      return(rmin);
    }
  else
    {
      /* loop inwards until above virial overdensity */
      for(ir=iRadOut-2; ir>0; ir--)
        {
          if(Rho[ir] > global.ovlim)
            {
              rleft    = r[ir];
              rright   = r[ir+1];
              Rholeft  = Rho[ir];
              Rhoright = Rho[ir+1];
              
              /* get interpolated virial radius */
              R_edge  = (rleft*(global.ovlim-Rhoright)+rright*(Rholeft-global.ovlim))
              /                  (Rholeft-Rhoright);
              goto found_edge;
            }
        }
      
      /* Rho never below virial overdensity -> search for upturn radius */
      iRadIn = 0;
      
      /* does the profile ever rise indicative of a sub-halo or mis-placed centre */
      nrise = 0;
      for(ir=iRadIn+1; ir<iRadOut; ir++)
        {
          /* profile is rising */
          if(Rho[ir] >= AHF_RISE*Rho[ir-1])
            nrise++;
          else
            nrise = 0;
          
          if(nrise > AHF_MAXNRISE)
            {
              goto found_rise;
            }
        }
      
    found_rise:
      iRadEnd = MIN(ir,iRadOut-1);
      
      /* reduce iRadEnd to the point where profiles declines faster than r^-slope */
      irm1    = MAX(iRadEnd-1,iRadIn);
      irp1    = MIN(iRadEnd+1,iRadOut-1);
      gradlf  = Rho[irm1] * pow(r[irm1], AHF_SLOPE);
      gradrt  = Rho[irp1] * pow(r[irp1], AHF_SLOPE);
      
      while(gradlf < gradrt)
        {
          iRadEnd--;
          
          if(iRadEnd <= iRadIn)
            {
              return(r[iRadIn]);
            }
          
          irm1    = MAX(iRadEnd-1,iRadIn);
          irp1    = MIN(iRadEnd+1,iRadOut-1);
          gradlf  = Rho[irm1] * pow(r[irm1], AHF_SLOPE);
          gradrt  = Rho[irp1] * pow(r[irp1], AHF_SLOPE);
        }
      
      R_edge = r[iRadEnd];
      
    found_edge:
      return (R_edge);
    }
}


/*==============================================================================
 * run through halo-tree locating the most appropriate host halo
 *==============================================================================*/
int get_hostHaloID(int haloID)
{
  int hostHaloID, itmp;
  
  /* identify first host halo */
  hostHaloID = halos[haloID].hostHalo;
  
  if(hostHaloID < 0)
    hostHaloID = halos[haloID].hostHaloGrids;
  
  if ( haloID < hostHaloID )
    {
      fprintf(stderr,"ERROR :: have disordered the halo tree!\n");
      exit(0);
    }
  
  return(hostHaloID);
}

/*====================================================================================
 * this routine actually returns the particles to the host halo
 * NOTE: at this stage we are still dealing with the linked-lists halos[].ll
 *====================================================================================*/
void return_parts_to_host(int ihalo, int hostHalo)
{
  partptr        cur_part, host_part;
  long unsigned jpart;
  
  /* loop over halo particles... */
  cur_part = halos[ihalo].ll;
  
  while( cur_part != NULL )
    {           
      /* insert the particle into the host's linked list */
      host_part              = halos[hostHalo].ll;
      
      halos[hostHalo].ll     = cur_part;
      /* move to next particle in halo (needs to be done *right* here!) */
      cur_part               = cur_part->ll;
      halos[hostHalo].ll->ll = host_part;
      
      /* Increase the number of particles by one for this halo */
      halos[hostHalo].nll++;
      
      /* NOTE:
       * it does not matter where we place the particle within the linked list
       * as there will be a call to sort_halo_particles() before we actually deal
       * with this host halo... */
      
    } /* while */
}

/*==============================================================================
 *  safe particles id's and give then give them back to host:
 *
 *   We need each halo to have access to all particles that potentially belong to it
 *    But at this stage the host halos are devoid of their subhalos' particles!
 *    Therefore, it is mandatory to hand these particles back to the hosts but
 *    without removing them from the subhalo...we are duplicating particles,
 *    but will then be able to perform the halo construction in parallel!
 *==============================================================================*/
void safe_return_particles(int ihalo)
{
  partptr       cur_part, host_part;
  int           hostHalo;  
  long unsigned jpart;
  
  /*-----------------------------------------------------------
   * transfer particles from linked-list over to ipart[] array
   *-----------------------------------------------------------*/
  halos[ihalo].ipart = (long unsigned *) calloc(halos[ihalo].npart, sizeof(long unsigned));
  if(halos[ihalo].ipart == NULL) 
    {
      fprintf(io.logfile,"\n NO MEMORY LEFT FOR ipart[] ARRAY!\n");
      exit(0);
    }
  
  cur_part = halos[ihalo].ll;
  jpart    = 0;
  while( cur_part != NULL )
    {    
      halos[ihalo].ipart[jpart] = cur_part - global.fst_part;
      
      jpart++;
      cur_part = cur_part->ll;
    }
  
  if(jpart != halos[ihalo].npart)
    {
      fprintf(stderr,"\n\n     safe_return_particles() WARNING: something wrong with halo %12d: %12ld vs. %12ld\n",ihalo,jpart,halos[ihalo].npart);
      
      /* try to adjust array (works only if jpart < npart) */
      if(jpart < halos[ihalo].npart)
        {
          fprintf(stderr,"     => adjusting npart and ipart[] array (you are loosing particles!?)!\n");
          halos[ihalo].ipart = (long unsigned *) realloc(halos[ihalo].ipart, jpart);
          halos[ihalo].npart = jpart;
        }
      else
        {
          fprintf(stderr,"     => you better investigate as particles have been written to wrong memory!\n");
        }
    }
  
  
  /*-----------------------------------------------------------
   * return particles to host
   *-----------------------------------------------------------*/
  /* identify host halo */
  hostHalo = get_hostHaloID(ihalo);
  
  /* transfer particles back to host (only if there is a credible host) */
  if (hostHalo >= 0 )
    return_parts_to_host(ihalo, hostHalo);
  
  
  /* do not touch halos[i].npart as we still require to know how many particles are in halos[i].ipart[] */
  halos[ihalo].ll = NULL;
  
  return; 
}

/*====================================================================================
 * return the particles of halo #ihalo to its host #hostHalo...
 * ...and remove them from halo #ihalo
 *
 * NOTE: as we cal return_partices() in a reverse-loop the particles
 *       will end up at the "upper most host halo"!
 *       -> this is important to keep in mind for gather_hostParts()
 *====================================================================================*/
void return_particles(int ihalo)
{
  int hostHalo;  
  
  /*-----------------------------------------------------------
   * return particles to host
   *-----------------------------------------------------------*/
  /* identify host halo */
  hostHalo = get_hostHaloID(ihalo);
  
#ifdef AHFverbose
  fprintf(io.logfile,"    return_particles: halos[%d].nll=%12ld -> ihost=%12d\n",
          ihalo, halos[ihalo].nll, hostHalo);
  fflush(io.logfile);
#endif
  
  /* transfer particles to host's linked-list */
  if (hostHalo >= 0 )
    return_parts_to_host(ihalo, hostHalo);
  
  /* once the particles have been returned remove them from the halo */
  if(hostHalo >= 0)
    {
      halos[ihalo].ll  = NULL;
      halos[ihalo].nll = 0;
    }
}

/*==============================================================================
 * sort the particles belonging to a halo with respects to distance
 *==============================================================================*/
void sort_halo_particles(HALO *halo)
{
  double         Xc, Yc, Zc, Xp, Yp, Zp;
  long unsigned  npart, i;
  double        *dist, dX, dY, dZ, dR;
  long unsigned *idx;
  long unsigned *ipart;
  partptr        cur_part;
  
#ifndef TRACKER
  if(halo->npart < AHF_MINPART)
#else
    if(halo->npart < TRK_MINPART)
#endif
      return;
  
#ifdef VERBOSE
#ifndef TRACKER	/* happens too often for tracker, too much output */
  fprintf(io.logfile,"    sort_particles:   npart=%12ld => ",halo->npart);
  fflush(io.logfile);
#endif
#endif
  
  Xc    = halo->pos.x;
  Yc    = halo->pos.y;
  Zc    = halo->pos.z;
  npart = halo->npart;
  
  ipart    = (long unsigned *)  calloc(npart,   sizeof(long unsigned));
  idx      = (long unsigned *)  calloc(npart+1, sizeof(long unsigned));
  dist     = (double *)         calloc(npart+1, sizeof(double));
  
  /* fill dist[] array with dR^2 */
  for(i=0; i<npart; i++)
    {
      /* access to particle */
      cur_part = global.fst_part + halo->ipart[i];
      
      /* particle position */
      Xp  = (double)cur_part->pos[X];
      Yp  = (double)cur_part->pos[Y];
      Zp  = (double)cur_part->pos[Z];
      
      /* put particle into halo rest frame */
      dX  = (Xp - Xc);
      dY  = (Yp - Yc);
      dZ  = (Zp - Zc);
      
      /* take care of periodic boundary conditions */
      if(dX >  0.5) dX -= 1.0;
      if(dY >  0.5) dY -= 1.0;
      if(dZ >  0.5) dZ -= 1.0;
      if(dX < -0.5) dX += 1.0;
      if(dY < -0.5) dY += 1.0;
      if(dZ < -0.5) dZ += 1.0;
      
      /* distance^2 of current particle */
      dR = (pow2(dX) + pow2(dY) + pow2(dZ));
      
      dist[i+1] = dR;
    }
  
  /* sort particles according to dist[] */
  indexx(npart, dist, idx);
  
  
  /* generate an ordered array ... */
  for(i=0; i<npart; i++)
    ipart[i] = halo->ipart[idx[i+1]-1];
  
  
  /* ... and replace halo.ipart[] */
  for(i=0; i<npart; i++)
    halo->ipart[i] = ipart[i];
  
  free(idx);
  free(dist);
  free(ipart);
  
#ifdef VERBOSE
  Xc    = halo->pos.x;
  Yc    = halo->pos.y;
  Zc    = halo->pos.z;
  
  for (i=0; i<5; i++)
    {
      
      cur_part = global.fst_part + halo->ipart[i];
      Xp  = (double)cur_part->pos[X];
      Yp  = (double)cur_part->pos[Y];
      Zp  = (double)cur_part->pos[Z];
      dX  = (Xp - Xc);
      dY  = (Yp - Yc);
      dZ  = (Zp - Zc);
      if(dX >  0.5) dX -= 1.0;
      if(dY >  0.5) dY -= 1.0;
      if(dZ >  0.5) dZ -= 1.0;
      if(dX < -0.5) dX += 1.0;
      if(dY < -0.5) dY += 1.0;
      if(dZ < -0.5) dZ += 1.0;
      dR = (pow2(dX) + pow2(dY) + pow2(dZ));
      
      fprintf(io.logfile,"dR_min = %16.8g kpc/h\n",sqrt(dR)*x_fac*1000.);
      fflush(io.logfile);
      
      /*  fprintf(stderr,"dR[%ld] = %16.8g kpc/h\n",i,sqrt(dR)*x_fac*1000.); */
    }
  
#endif
  Xc    = halo->pos.x;
  Yc    = halo->pos.y;
  Zc    = halo->pos.z;
  
}

#if (defined AHFrestart || defined WITH_MPI)
void
rem_boundary_haloes(void){
  int i;
  sfc_key_t key;
  sfc_curve_t ctype;
  sfc_key_t minkey;
  sfc_key_t maxkey;
  int level;
  HALO *tmpHalos;
  int k;
  int new_numHalos = 0;
  
  io_logging_msg(global_io.log, INT32_C(0),
                 "\nFiguring out which halo can be ignored.");
  
#	ifdef AHFrestart
  minkey = global_info.minkey;
  maxkey = global_info.maxkey;
  ctype = global_info.ctype;
  level = global_info.level;
#	else
  minkey = global_info.loadbal->fstkey[global_mpi.rank];
  maxkey = global_info.loadbal->lstkey[global_mpi.rank];
  ctype = global_info.loadbal->ctype;
  level = global_info.loadbal->level;
#	endif
  
  /* Loop over all halos and figure out which should be used */
#	pragma omp parallel for \
schedule (dynamic) \
shared(numHalos, halos, level, ctype, \
global_io, minkey, maxkey) \
private(key) \
reduction (+:new_numHalos)
  for(i=0; i<numHalos; i++) {
    /* First get the key*/
    key = sfc_curve_calcKey(ctype,
                            (double)(halos[i].pos.x),
                            (double)(halos[i].pos.y),
                            (double)(halos[i].pos.z),
                            BITS_PER_DIMENSION);
    /* Reduce it */
    key = sfc_curve_contract(level,
                             BITS_PER_DIMENSION,
                             ctype,
                             key);
    
    /* Check if it is with our key range */
    if ( (key>=minkey) && (key<=maxkey) ) {
      halos[i].ignoreme = FALSE;
      new_numHalos++;
    } else {
      halos[i].ignoreme = TRUE;
    }
  }
  io_logging_msg(global_io.log, INT32_C(0),
                 "Done with identifying halos to be ignored.");
  
  /* Now remove them */
  io_logging_msg(global_io.log, INT32_C(0),
                 "Now really getting rid of outsider halos.");
  
  /*That's going to be the new Halo list */
  tmpHalos = (HALO *)malloc(sizeof(HALO)*new_numHalos);
  if (tmpHalos == NULL)
    common_terminate(EXIT_FAILURE);
  
  /* Loop over all halos and only pick the good ones */
  new_numHalos = 0;
  for (k=0; k<numHalos; k++) {
    if (halos[k].ignoreme == FALSE) {
      memcpy((void *)(tmpHalos+new_numHalos),
             (const void *)(halos+k),
             sizeof(HALO));
      new_numHalos++;
    }
  }
  
  /* Status information */
  io_logging_msg(global_io.log, INT32_C(0),
                 "Old number of halos: %i   "
                 "New number of halos: %i",
                 (int)numHalos, (int)new_numHalos);
  
  /* Replace the old halo list with the new one */
  free(halos);
  halos = tmpHalos;
  numHalos = new_numHalos;
  io_logging_msg(global_io.log, INT32_C(0),
                 "Done with getting rid of halos.");
  
  /* Done */
  return;
}
#endif /* (defined AHFrestart || defined WITH_MPI) */


#ifdef AHFbinary

#define LEN_TYPEFIELD 128
#define WRITE_I fwrite(&tmp_int, sizeof(bin_int_t), 1, f)
#define WRITE_F fwrite(&tmp_float, sizeof(bin_float_t), 1, f)

typedef uint32_t bin_int_t;
typedef float bin_float_t;


void ahf_write_open_files(FILE **f,
                          FILE **f_info,
                          char *prefix,
                          char *suffix)
{
	char fname[2048];
	char fname_info[2048];

	/* Generate the filenames */
	sprintf(fname, "%s.%s_bin", prefix, suffix);
	sprintf(fname_info, "%s.%s_info", prefix, suffix);

	/* Open binary output file */
	*f = fopen(fname, "wb");
	if (*f == NULL) {
		fprintf(stderr, "could not open %s with mode `wb'\n", fname);
		exit(EXIT_FAILURE);
	}

	/* Open the info file (in MPI mode, only one does that) */
	*f_info = NULL;
#	if (defined WITH_MPI || defined AHFrestart)
#		ifdef WITH_MPI
	if (global_mpi.rank == 0)
#		else
	if (global_info.rank == 0)
#		endif
#	endif
	{
		*f_info = fopen(fname_info, "w");
		if (*f_info == NULL) {
			fprintf(stderr, "could not open %s with mode `w'\n",
					fname_info);
			exit(EXIT_FAILURE);
		}
	}

#	if (defined VERBOSE)
	fprintf(stderr, "%s\t", fname);
#	endif

	/* Done */
	return;
}

void ahf_write_profiles(char *prefix,
                        HALO *halos,
                        unsigned long *idx,
                        int numHalos)
{
	FILE *f;
	FILE *f_info;
	int i=0, j=0;
	int num_columns;
	uint32_t sizes[2];
	bin_int_t tmp_int;
	bin_float_t tmp_float;
	int r_conv_i, ibin;
	double age, rad, t_relax;
	uint64_t real_num_halos;
	uint64_t total_num_lines;
	uint64_t tmp;
	int8_t typefield[LEN_TYPEFIELD];
	
	/* Store the number of bytes used for integer and float values */
	sizes[0] = sizeof(bin_int_t);
	sizes[1] = sizeof(bin_float_t);

	/* Initialize typefield to float */
	for (i=0; i<LEN_TYPEFIELD; i++)
		typefield[i] = INT8_C(1);

	/* Open the files */
	ahf_write_open_files(&f, &f_info, prefix, "AHF_profiles");

	/* Write info, if required */
	i = 0;
	if (f_info != NULL)
		fprintf(f_info,
		        "r(%i)\nnpart(%i)\nnvpart(%i)\novdens(%i)\ndens(%i)\n"
		        "vcirc(%i)\nsigv(%i)\nLx(%i)\nLy(%i)\nLz(%i)\n"
		        "a(%i)\nEax(%i)\nEay(%i)\nEaz(%i)\n"
		        "b(%i)\nEbx(%i)\nEby(%i)\nEbz(%i)\n"
		        "c(%i)\nEcx(%i)\nEcy(%i)\nEcz(%i)\n"
		        "Ekin(%i)\nEpot(%i)\n",
		        i+1, i+2, i+3, i+4, i+5,
		        i+6, i+7, i+8, i+9, i+10,
		        i+11, i+12, i+13, i+14,
		        i+15, i+16, i+17, i+18,
		        i+19, i+20, i+21, i+22,
		        i+23, i+24);
	i+=24;
	typefield[1] = INT8_C(0);
#	if (defined GAS_PARTICLES)
	if (f_info != NULL)
		fprintf(f_info,
		        "nvpart_gas(%i)\nnvpart_stars(%i)",
		        i+1, i+2);
	i+=2;
#	endif
#	if (defined AHFphspdens)
	if (f_info != NULL)
		fprintf(f_info,
		        "sigma2_vx_sh(%i)\nsigma2_vy_sh(%i)\nsigma2_vz_sh(%i)\n"
		        "sigma2_vr_sh(%i)\nsigma2_vtheta_sh(%i)\nsigma2_vphi_sh(%i)\n",
		        i+1, i+2, i+3,
		        i+4, i+5, i+6);
	i+=6;
#		if (defined AHFmeanvelocities)
	if (f_info != NULL)
		fprintf(f_info,
		        "mean_vx_sh(%i)\nmean_vy_sh(%i)\nmean_vz_sh(%i)\n"
		        "mean_vr_sh(%i)\nmean_vtheta_sh(%i)\nmean_vphi_sh(%i)\n"
		        "mean_vx_sp(%i)\nmean_vy_sp(%i)\nmean_vz_sp(%i)\n"
		        "mean_vr_sp(%i)\nmean_vtheta_sp(%i)\nmean_vphi_sp(%i)\n",
		        i+1, i+2, i+3,
		        i+4, i+5, i+6,
		        i+7, i+8, i+9,
		        i+10, i+11, i+12);
	i+=12;
#		endif
#	endif
	num_columns = i;

	/* Store the sizes */
	fwrite(sizes, sizeof(uint32_t), 2, f);

	/* 
	 * Write a dummy value for the number of halos and total number of
	 * lines 
	 */
	real_num_halos = UINT64_C(0);
	total_num_lines = UINT64_C(0);
	fwrite(&real_num_halos, sizeof(uint64_t), 1, f);
	fwrite(&total_num_lines, sizeof(uint64_t), 1, f);

	/* Store the number of columns */
	fwrite(&num_columns, sizeof(uint32_t), 1, f);

	/* Store the typefield */
	fwrite(typefield, sizeof(int8_t), num_columns, f);

	/* Now loop over all haloes for writing */
	for (j=0; j<numHalos; j++) {
		i = idx[j];

		if (    halos[i].npart >= AHF_MINPART
		     && (double)halos[i].npart*simu.min_weight/(double)halos[i].M_vir
		        >AHF_MASSMIX) {

			/* First identify the converged radial bins */
			for (r_conv_i=0, ibin=0; ibin<halos[i].prof.nbins; ibin++) {
				/* check for converged radius (Power et al. 2003) */
				rad = halos[i].prof.r[ibin];
				t_relax =   halos[i].prof.npart[ibin]/log(rad/halos[i].spaRes)
				           * rad / sqrt(halos[i].prof.v2_circ[ibin]);
              
				/* convert to (Mpc/h) / (km/sec) ~ (Mpc/h) / (kpc/Gyr) */
				t_relax *= r_fac/sqrt(phi_fac);
              
				/* convert to Gyr/h */
				t_relax *= 1E3;
              
				/* age of the Universe in Gyr/h */
				age      = calc_t(global.a) * simu.t_unit*Mpc/Gyr;
              
				/* if not converged, write negative radius into .AHF_profiles */
				if (t_relax < 0.9*age)
					/* The +1 is merely to keep the name consistent:
					 * r_conv_i should give the smallest converged radius, not
					 * the bin before it. Hence to check for converged bin or
					 * not, i<r_conv_i is the thing to do. */
					r_conv_i = ibin+1;
			}

			/* Write the number of profile lines and update total number */
			tmp = (uint64_t)(halos[i].prof.nbins);
			fwrite(&tmp, sizeof(uint64_t), 1, f);
			total_num_lines += tmp;

			/* Now we actually write the stuff to the file */
			for (ibin=0; ibin<halos[i].prof.nbins; ibin++) {
				tmp_float = (bin_float_t)(halos[i].prof.r[ibin]);
				tmp_float = (ibin<r_conv_i) ? - tmp_float : tmp_float;
				tmp_float *= x_fac * 1000.; WRITE_F;
				tmp_int = (bin_int_t)(halos[i].prof.npart[ibin]); WRITE_I;
				tmp_float = (bin_float_t)(halos[i].prof.nvpart[ibin]); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.ovdens[ibin]*rho_fac/global.rho_b); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.dens[ibin]*rho_fac/global.rho_b); WRITE_F;
				tmp_float = (bin_float_t)(sqrt(halos[i].prof.v2_circ[ibin]*phi_fac)); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.sig_v[ibin]*v_fac); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.Lx[ibin]*m_fac*r_fac*v_fac); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.Ly[ibin]*m_fac*r_fac*v_fac); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.Lz[ibin]*m_fac*r_fac*v_fac); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.axis1[ibin]); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.E1x[ibin]); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.E1y[ibin]); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.E1z[ibin]); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.axis2[ibin]); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.E2x[ibin]); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.E2y[ibin]); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.E2z[ibin]); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.axis3[ibin]); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.E3x[ibin]); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.E3y[ibin]); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.E3z[ibin]); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.Ekin[ibin]*m_fac*pow2(v_fac)); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.Epot[ibin]*m_fac*phi_fac); WRITE_F;
#	if (defined GAS_PARTICLES)
				tmp_float = (bin_float_t)(halos[i].prof.nvpart_gas[ibin]); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.nvpart_star[ibin]); WRITE_F;
#	endif
#	if (defined AHFphspdens)
				tmp_float = (bin_float_t)(halos[i].prof.sigma2_vx_sh[ibin]*v_fac*v_fac); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.sigma2_vy_sh[ibin]*v_fac*v_fac); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.sigma2_vz_sh[ibin]*v_fac*v_fac); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.sigma2_vr_sh[ibin]*v_fac*v_fac); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.sigma2_vtheta_sh[ibin]); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.sigma2_vphi_sh[ibin]); WRITE_F;
#		if (defined AHFmeanvelocities)
				tmp_float = (bin_float_t)(halos[i].prof.mean_vx_sh[ibin]*v_fac); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.mean_vy_sh[ibin]*v_fac); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.mean_vz_sh[ibin]*v_fac); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.mean_vr_sh[ibin]*v_fac); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.mean_vtheta_sh[ibin]); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.mean_vphi_sh[ibin]); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.mean_vx_sp[ibin]*v_fac); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.mean_vy_sp[ibin]*v_fac); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.mean_vz_sp[ibin]*v_fac); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.mean_vr_sp[ibin]*v_fac); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.mean_vtheta_sp[ibin]); WRITE_F;
				tmp_float = (bin_float_t)(halos[i].prof.mean_vphi_sp[ibin]); WRITE_F;
#		endif
#	endif
			} /* End of for-loop over profile lines */

			/* And finally increment the halo counter */
			real_num_halos++;
		} /* End of if for suitable halo */
	}

	/* Rewind and put the right numbers in the front (after the sizes) */
	fseek(f, 8L, SEEK_SET);
	fwrite(&real_num_halos, sizeof(uint64_t), 1, f);
	fwrite(&total_num_lines, sizeof(uint64_t), 1, f);
	
	/* Close the files */
	fclose(f);
	if (f_info != NULL)
		fclose(f_info);

#	if (defined VERBOSE)
	/* End the 'Writing file..' statement started when opening the file */
	fprintf(stderr, "done\n");
#	endif

	/* Done */
	return;
}

void ahf_write_particles(char *prefix,
                         HALO *halos,
                         unsigned long *idx,
                         int numHalos)
{
	FILE *f;
	FILE *f_info;
	uint32_t sizes[2];
	bin_int_t tmp_int;
	bin_float_t tmp_float;
	int i=0, j=0, k=0;
	uint32_t num_columns;
	uint64_t real_num_halos;
	uint64_t total_num_particles;
	uint64_t tmp;
	int8_t particle_type;
	double particle_mass;
	partptr cur_part;
	int8_t typefield[LEN_TYPEFIELD];

	/* Store the number of bytes used for integer and float values */
	sizes[0] = sizeof(bin_int_t);
	sizes[1] = sizeof(bin_float_t);

	/* Initialize typefield to int */
	for (i=0; i<LEN_TYPEFIELD; i++)
		typefield[i] = INT8_C(0);

	/* Open the files */
	ahf_write_open_files(&f, &f_info, prefix, "AHF_particles");

	/* Write info, if required */
	i = 0;
	if (f_info != NULL)
		fprintf(f_info, "ID(%i)\n", i+1);
	i += 1;
#	if (defined PARTICLES_INFO)
	if (f_info != NULL)
		fprintf(f_info, "pType(%i)\n", i+1);
	i += 1;
#	endif
	num_columns = i;

	/* Store the sizes */
	fwrite(sizes, sizeof(uint32_t), 2, f);

	/* Write dummy values */
	real_num_halos = UINT64_C(0);
	total_num_particles = UINT64_C(0);
	fwrite(&real_num_halos, sizeof(uint64_t), 1, f);
	fwrite(&total_num_particles, sizeof(uint64_t), 1, f);

	/* Write number of columns */
	fwrite(&num_columns, sizeof(uint32_t), 1, f);

	/* Store the typefield */
	fwrite(typefield, sizeof(int8_t), num_columns, f);

	/* Now loop over all haloes for writing */
	for (j=0; j<numHalos; j++) {
		i = idx[j];

		if (    halos[i].npart >= AHF_MINPART
		     && (double)halos[i].npart*simu.min_weight/(double)halos[i].M_vir
		        >AHF_MASSMIX) {
			/* Write the number of particles in this halo and update total number */
			tmp = (uint64_t)(halos[i].npart);
			fwrite(&tmp, sizeof(uint64_t), 1, f);
			total_num_particles += tmp;

			/* Loop over all particles in this halo */
			for (k=0; k<halos[i].npart; k++) {
				cur_part = global.fst_part + halos[i].ipart[k];
				tmp_int = (bin_int_t)(cur_part->id); WRITE_I;
#	if (defined PARTICLES_INFO)
#		if (!defined GAS_PARTICLES)
				particle_type = 1;
#		else
				particle_type = (cur_part->u >= PGAS) ? (int8_t)PGAS : (int8_t)(-cur_part->u);
#		endif
				fwrite(&particle_type, sizeof(int8_t), 1, f);
#	endif
			} /* End of loop over particles */

			/* Count this halo */
			real_num_halos++;
		} /* end of if selecting 'proper' haloes */
	} /* end of for looping over all haloes*/

	/* Rewind and put the right numbers in the front (after the sizes) */
	fseek(f, 8L, SEEK_SET);
	fwrite(&real_num_halos, sizeof(uint64_t), 1, f);
	fwrite(&total_num_particles, sizeof(uint64_t), 1, f);

	/* Close the files */
	fclose(f);
	if (f_info != NULL)
		fclose(f_info);

#	if (defined VERBOSE)
	/* End the 'Writing file..' statement started when opening the file */
	fprintf(stderr, "done\n");
#	endif

	/* Done */
	return;
}

void ahf_write_halos(char *prefix,
                     HALO *halos,
                     unsigned long *idx,
                     int numHalos)
{
	FILE *f;
	FILE *f_info;
	uint32_t num_columns;
	uint32_t sizes[2];
	int i, j;
	bin_int_t tmp_int;
	bin_float_t tmp_float;
	uint64_t real_num_halos;
	int8_t typefield[LEN_TYPEFIELD];

	/* Store the number of bytes used for integer and float values */
	sizes[0] = sizeof(bin_int_t);
	sizes[1] = sizeof(bin_float_t);

	/* Initialize typefield to float */
	for (i=0; i<LEN_TYPEFIELD; i++)
		typefield[i] = INT8_C(1);

	/* Open the files */
	ahf_write_open_files(&f, &f_info, prefix, "AHF_halos");

	/* Write info, if required */
	i = 0;
	if (f_info != NULL)
		fprintf(f_info,
		        "hid(%i)\nnpart(%i)\nnvpart(%i)\n"
		        "Xc(%i)\nYc(%i)\nZc(%i)\nVXc(%i)\nVYc(%i)\nVZc(%i)\n"
		        "Mvir(%i)\nRvir(%i)\nVmax(%i)\nRmax(%i)\nsigV(%i)\n"
		        "lambda(%i)\nLx(%i)\nLy(%i)\nLz(%i)\n"
		        "a(%i)\nEax(%i)\nEay(%i)\nEaz(%i)\n"
		        "b(%i)\nEbx(%i)\nEby(%i)\nEbz(%i)\n"
		        "c(%i)\nEcx(%i)\nEcy(%i)\nEcz(%i)\n"
		        "ovdens(%i)\nRedge(%i)\nEkin(%i)\nEpot(%i)\n"
		        "mbp_offset(%i)\ncom_offset(%i)\nr2(%i)\nlambdaE(%i)\n",
		        i+1, i+2, i+3,
		        i+4, i+5, i+6, i+7, i+8, i+9,
		        i+10, i+11, i+12, i+13, i+14,
		        i+15, i+16, i+17, i+18,
		        i+19, i+20, i+21, i+22,
		        i+23, i+24, i+25, i+26,
		        i+27, i+28, i+29, i+30,
		        i+31, i+32, i+33, i+34,
		        i+35, i+36, i+37, i+38);
	i+=38;
	typefield[0] = INT8_C(0);
	typefield[1] = INT8_C(0);
#	if (defined GAS_PARTICLES)
	if (f_info != NULL)
		fprintf(f_info,
		        "n_gas(%i)\nM_gas(%i)\nlambda_gas(%i)\n"
		        "Lx_gas(%i)\nLy_gas(%i)\nLz_gas(%i)\n"
		        "a_gas(%i)\nEax_gas(%i)\nEay_gas(%i)\nnEaz_gas(%i)\n"
		        "b_gas(%i)\nEbx_gas(%i)\nEby_gas(%i)\nnEbz_gas(%i)\n"
		        "c_gas(%i)\nEcx_gas(%i)\nEcy_gas(%i)\nnEcz_gas(%i)\n"
		        "Ekin_gas(%i)\nEpot_gas(%i)\nlambdaE_gas(%i)\nphi0(%i)\n"
		        "n_star(%i)\nM_star(%i)\nlambda_star(%i)\n"
		        "Lx_star(%i)\nLy_star(%i)\nLz_star(%i)\n"
		        "a_star(%i)\nEax_star(%i)\nEay_star(%i)\nnEaz_star(%i)\n"
		        "b_star(%i)\nEbx_star(%i)\nEby_star(%i)\nnEbz_star(%i)\n"
		        "c_star(%i)\nEcx_star(%i)\nEcy_star(%i)\nnEcz_star(%i)\n"
		        "Ekin_star(%i)\nEpot_star(%i)\nlambdaE_star(%i)\n",
		        i+1, i+2, i+3,
		        i+4, i+5, i+6,
		        i+7, i+8, i+9, i+10,
		        i+11, i+12, i+13, i+14,
		        i+15, i+16, i+17, i+18,
		        i+19, i+20, i+21, i+22,
		        i+23, i+24, i+25,
		        i+26, i+27, i+28,
		        i+29, i+30, i+31, i+32,
		        i+33, i+34, i+35, i+36,
		        i+37, i+38, i+39, i+40,
		        i+41, i+42, i+43);
	i+=43;
	typefield[38] = INT8_C(0);
	typefield[60] = INT8_C(0);
#	endif
	num_columns = (uint32_t)i;

	/* Store the sizes */
	fwrite(sizes, sizeof(uint32_t), 2, f);

	/* Write a dummy value for the number of halos */
	real_num_halos = UINT64_C(0);
	fwrite(&real_num_halos, sizeof(uint64_t), 1, f);

	/* Store the number of columns used */
	fwrite(&num_columns, sizeof(uint32_t), 1, f);

	/* Store the typefield */
	fwrite(typefield, sizeof(int8_t), num_columns, f);

	/* Now loop over all haloes for writing */
	for (j=0; j<numHalos; j++) {
		i = idx[j];

		if (    halos[i].npart >= AHF_MINPART
		     && (double)halos[i].npart*simu.min_weight/(double)halos[i].M_vir
		        >AHF_MASSMIX) {
			tmp_int = (bin_int_t)real_num_halos; WRITE_I;
			tmp_int = (bin_int_t)(halos[i].npart); WRITE_I;
			tmp_float = (bin_float_t)(halos[i].M_vir); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].pos.x*x_fac); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].pos.y*x_fac); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].pos.z*x_fac); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].vel.x*v_fac); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].vel.y*v_fac); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].vel.z*v_fac); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].M_vir*m_fac); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].R_vir*x_fac*1000.); WRITE_F;
			tmp_float = (bin_float_t)sqrt(halos[i].V2_max*phi_fac); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].R_max*x_fac*1000.); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].velDis*v_fac); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].lambda); WRITE_F;
#	if (defined AHFabsangmom)
			tmp_float = (bin_float_t)(halos[i].AngMom.x*r_fac*v_fac*m_fac); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].AngMom.y*r_fac*v_fac*m_fac); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].AngMom.z*r_fac*v_fac*m_fac); WRITE_F;
#	else
			tmp_float = (bin_float_t)(halos[i].AngMom.x); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].AngMom.y); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].AngMom.z); WRITE_F;
#	endif
			tmp_float = (bin_float_t)(halos[i].axis.x); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].E1.x); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].E1.y); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].E1.z); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].axis.y); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].E2.x); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].E2.y); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].E2.z); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].axis.z); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].E3.x); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].E3.y); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].E3.z); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].ovdens); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].R_edge*x_fac*1000.); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].Ekin*m_fac*pow2(v_fac)); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].Epot*m_fac*pow2(v_fac)); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].mbp_offset*x_fac*1000.); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].com_offset*x_fac*1000.); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].r2*x_fac*1000.); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].lambdaE); WRITE_F;
#	if (defined GAS_PARTICLES)
			tmp_int = (bin_int_t)(halos[i].gas.npart); WRITE_I;
			tmp_float = (bin_float_t)(halos[i].gas.M_vir*m_fac); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].gas.lambda); WRITE_F;
#		if (defined AHFabsangmom)
			tmp_float = (bin_float_t)(halos[i].gas.AngMom.x*r_fac*v_fac*m_fac); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].gas.AngMom.y*r_fac*v_fac*m_fac); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].gas.AngMom.z*r_fac*v_fac*m_fac); WRITE_F;
#		else
			tmp_float = (bin_float_t)(halos[i].gas.AngMom.x); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].gas.AngMom.y); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].gas.AngMom.z); WRITE_F;
#		endif
			tmp_float = (bin_float_t)(halos[i].gas.axis.x); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].gas.E1.x); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].gas.E1.y); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].gas.E1.z); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].gas.axis.y); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].gas.E2.x); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].gas.E2.y); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].gas.E2.z); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].gas.axis.z); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].gas.E3.x); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].gas.E3.y); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].gas.E3.z); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].gas.Ekin*m_fac*pow2(v_fac)); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].gas.Epot*m_fac*phi_fac); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].gas.lambdaE); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].Phi0*m_fac*phi_fac); WRITE_F;
			tmp_int = (bin_int_t)(halos[i].star.npart); WRITE_I;
			tmp_float = (bin_float_t)(halos[i].star.M_vir*m_fac); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].star.lambda); WRITE_F;
#	if (defined AHFabsangmom)
			tmp_float = (bin_float_t)(halos[i].star.AngMom.x*r_fac*v_fac*m_fac); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].star.AngMom.y*r_fac*v_fac*m_fac); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].star.AngMom.z*r_fac*v_fac*m_fac); WRITE_F;
#	else
			tmp_float = (bin_float_t)(halos[i].star.AngMom.x); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].star.AngMom.y); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].star.AngMom.z); WRITE_F;
#	endif
			tmp_float = (bin_float_t)(halos[i].star.axis.x); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].star.E1.x); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].star.E1.y); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].star.E1.z); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].star.axis.y); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].star.E2.x); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].star.E2.y); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].star.E2.z); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].star.axis.z); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].star.E3.x); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].star.E3.y); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].star.E3.z); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].star.Ekin*m_fac*pow2(v_fac)); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].star.Epot*m_fac*phi_fac); WRITE_F;
			tmp_float = (bin_float_t)(halos[i].star.lambdaE); WRITE_F;
#	endif

			real_num_halos++;
		} /* End of if selecting proper haloes */
	} /* End of halo loop */

	/* Rewind and put the right numbers in the front (after the sizes) */
	fseek(f, 8L, SEEK_SET);
	fwrite(&real_num_halos, sizeof(uint64_t), 1, f);

	/* Close the files */
	fclose(f);
	if (f_info != NULL)
		fclose(f_info);

#	if (defined VERBOSE)
	/* End the 'Writing file..' statement started when opening the file */
	fprintf(stderr, "done\n");
#	endif

	/* Done */
	return;
}


#undef BUFFER_TOO_SMALL_ERROR_PROFILES

#endif

#endif /* AHF */
