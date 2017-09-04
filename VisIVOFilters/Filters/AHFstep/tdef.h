#include <stdio.h>
#include <time.h>
#include "define.h"
#include "param.h"
//
//
#ifndef TDEF_INCLUDED
#define TDEF_INCLUDED

#if (defined WITH_MPI || defined NEWSTARTRUN)
#	ifdef HAVE_STDINT_H
#		include <stdint.h>
#	else
#		include <replace_stdint.h>
#	endif
#	include "libsfc/sfc.h"
#endif


/*--------------
 * Useful Types
 *--------------*/

typedef float  fvect[NDIM];      /* "vector" of floats  */
typedef double dvect[NDIM];      /* "vector" of doubles */
typedef char   boolean;          /* C has no boolean data type */
#ifdef DOUBLE
typedef double flouble;
#else
typedef float  flouble;
#endif

/*=========================================================================
 *  particle structure type
 *=========================================================================*/

#if (defined WITH_MPI || defined NEWSTARTRUN)
/** The ID type */
typedef uint64_t part_id_t;
#	define PRIpartid PRIu64
#endif

typedef struct particle *partptr;
typedef struct particle
{
  /* linked-list coord  
   *    :- NOTE -: 
   * for the linkedlist sort to work this element *must* be first in the struct
   */
  partptr ll;


  flouble   pos[NDIM];                	/* position vector             */
  flouble   mom[NDIM];                	/* momentum vector             */
#ifdef MULTIMASS
  flouble   weight;                    /* weight of particle          */
#endif

#if (defined GADGET && !defined NEWSTARTRUN)
#ifdef LGADGET
  long long ID;
#else
  int       ID;
#endif
#endif /* GADGET && !NEWSTARTRUN */

#if (defined WITH_MPI || defined NEWSTARTRUN)
  part_id_t id;
  sfc_key_t sfckey;
#endif
#if (defined GAS_PARTICLES && defined NEWSTARTRUN)
  flouble    u;
#endif 

#if (defined TRACKER && defined PARTICLES_INFO)
//   int      species; 
#endif

#ifdef FORCETEST
  double   forces[NDIM];
  double   pot;
  double   dens;
#endif
  
#ifdef DEVA
  long itype;
  long ID;
#endif
}part;

/*=========================================================================
*  gas particle structure type
*=========================================================================*/
typedef struct gas *gasptr;
typedef struct gas
{
   /* positions and masses are stored along with "normal" particles... */
   flouble u;
}gas;


/*=========================================================================
*  star particle structure type
*=========================================================================*/
typedef struct star *starptr;
typedef struct star
{
   /* positions and masses are stored along with "normal" particles... */
   flouble dummy;
}star;



/*========================================================================
 * Grid structures and associated pointers
 *========================================================================*/

/*---------------------------------------------------
 * The node or gridpoint structure for refinements.
 *---------------------------------------------------*/
typedef struct node *nptr;
typedef struct node
{
  flouble  pot;              	    /* gravitational potential              */
  flouble  dens;                  /* mass density                         */
  partptr  ll;              	    /* head of particle linked-list         */

#ifdef HYDRO
  flouble   u[NHYDRO];             /* vector of hydrodynamic variables              */
  flouble   u_tmp[NHYDRO];         /* immediately update rather than storing fluxes */
  flouble   F[NHYDRO];             /* 3D flux into/out of this cell                 */
#ifdef MHD
  flouble   B[NDIM];               // vector of staggered (!) magnetic field components (see Ziegler 2004)
                                  // B[X] = Bx_(i-1/2,j,k); B[Y] = By_(i,j-1/2,k); B[Z] = Bz_(i,j,k-1/2); mu=1
  flouble   B_tmp[NDIM];          // analogous to u_tmp
  flouble   E[NDIM];              // E-fluxes for integrating the B field according to KNPCT scheme
#endif /* MHD */
#ifdef HYDRO_STORE_FLUX
  flouble   Flux[NDIM][2][NHYDRO]; /* full flux information in all directions       */
#endif
#endif /* HYDRO */

  /*
   * The forces at the node are only required at the end of a timestep.
   * Hence the space allocated to the forces is used to store temporary
   * data at other times.
   */
  union forcestorage
  {
    flouble  forces[NDIM];       /* grid force components                   */
    flouble  temp[NDIM];         /* old potential, source term, whatever... */
    partptr  new_ll;             /* pointer for rebuilding linked list      */
    int      colour;             /* used with -DAHF                         */
  }force;

}node;


/*-------------------------
 * QUAD logical structures
 *-------------------------*/

typedef struct nquad *nqptr;
typedef struct nquad
{
  nptr loc;                       /* pointer to first node in block*/
  long x;                         /* x-coord of first node in block */
  long length;                    /* length of block of nodes */
  nqptr next;               	  /* pointer to next nquad (or null) */
}nquad;

typedef struct cquad *cqptr;
typedef struct cquad
{
  nqptr loc;                      /* pointer to first nquad in block */
  long y;                         /* y-coord of first nquad in block */
  long length;                    /* length of block of nquads */
  cqptr next;               	  /* pointer to next cquad (or null) */
}cquad;

typedef struct pquad *pqptr;
typedef struct pquad
{
  cqptr loc;                	  /* pointer to first cquad in block */
  long z;                         /* z-coord of first cquad in block */
  long length;                    /* length of block of cquads */
  pqptr next;                     /* pointer to next pquad (or null) */
}pquad;

/*---------------------------
 * the grid structure itself
 *---------------------------*/

typedef struct gridlist
{
  pqptr   pquad;                 /* pointer to first pquad on grid           */
  long    no_pquad;              /* how many pquad's are there?              */
  pqptr  *pquad_array;           /* align pquads in memory when using OpenMP */
#ifdef STORE_REFS
  pqptr   old_pquad;             /* pointer to previous time step            */
#endif

  long unsigned    l1dim;        /* virtual grid length                      */
  double  spacing;               /* grid spacing                             */
  double  spacing2;              /* grid spacing squared                     */
  double  critdens;              /* critical density on this grid            */
  double  masstodens;            /* mass to density conversion factor        */
  double  masstopartdens;        /* mass to partdens converion factor        */

  double  old_resid;             /* old residual                             */
  double  cur_resid;             /* current residual                         */
  double  trunc_err;             /* the truncation error                     */
  long    no_sweeps;             /* number of sweeps done on grid            */

  double  timecounter;           /* integration variable                     */
  double  timestep;              /* current timestep for integration variable*/
  int     multistep;             /* keep track of current multistep-phase    */
  
  struct  grid_size
  {
    long unsigned   no_part;     /* number of particles linked to grid       */
    long unsigned   no_nodes;    /* number of nodes actually present         */
  }size;

  struct  grid_time
  {
    time_t potential;             /* deriving the potenital by GS sweeps      */
    time_t density;               /* deriving the proper densities            */
    time_t DK;                    /* drifting and kicking particles           */
    time_t grid;                  /* everything related to grid hierarchy     */
    time_t hydro;                 /* time spent for the hydro-solver          */
  }time;

  struct grid_leavers
  {
    int     no_sendback;
    int     no_keepmoving;
    partptr *send_back;           /* how to treat particles that...           */
    partptr *keep_moving;         /* crossed the boundary.                    */
  }leavers;

  boolean  next;                  /* the grid structure is never destroyed    */

#ifdef WITH_MPI
	/** Number of boundary nodes*/
	uint32_t num_bound;
	/** The boundary nodes */
	nptr bound;
#endif

}gridls;

/*=============================================================================
 * main:   structure carrying useful common information
 *=============================================================================*/
typedef struct info_global
{
  /* access to the domain grid */
  gridls *dom_grid;      /* actual pointer to domain grid strcuture             */
  int     domgrid_no;    /* locate domain grid within grid_list                 */

  /* access to all particles currently in use (NOTE: global.no_part is not necessarily equal to simu.no_part!) */
  partptr       fst_part;
  long unsigned no_part;       /* here we store the actual number of particles        */
                               /* simu.no_part gives the total number of particles in the simulation while 
                                  global.no_part counts the actual particles currently used
                                  (cf. AHFmmfocus!)*/
  gasptr        fst_gas;
  long unsigned no_gas;        /* how many gas particles               */
  long          offset_gas;    /* for those input files that contain gas particles         */
  
  starptr       fst_star;
  long unsigned no_stars;      /* how many star particls               */
  long          offset_stars;  /* for those input files that contain gas particles         */
  
  
  /* access to current timecounter */
  double  a;
  double  t;
  double  z;
  double  super_t;
  int     no_timestep;

  long    fin_l1dim;     /* finest refinement level reached per domain cycle    */
  boolean fst_cycle;
  
  double  max_dr2;       /* maximum change in particle position                 */
  double  max_dp2;       /* maximum change in particle momentum                 */
  double  mean_dr;
  double  mean_dp;
  long    no_kicks;
  long    no_drifts;
  
  double  cfl_speed;

  int     output_count;
  boolean restart;
  
  double  total_time;

  double  bytes_node;
  double  bytes_part;
  
  boolean       speeding;     /* flag that indicates speeding particles           */
  unsigned long no_speeders;  /* number of particles caught speeding              */
  unsigned long no_ssteps;    /* normalisation for speeder particle counter       */

  unsigned long no_leavers;   /* number of particles crossing grid boundaries     */
  unsigned long no_lsteps;    /* normalisation for leaver particle counter        */

  boolean       terminate;    /* flag that terminates code after current step     */
  char         *termfile_name;

  int           architecture; /* little or big endian machine                     */
    
  long unsigned ndummy;
  double        fdummy;

  
  double  ovlim;         /* virial overdensity parameter                             */
  double  rho_b;         /* stores comoving(!) background density                    */
  double  rho_vir;       /* stores comoving(!)   virial   density                    */
  double  max_ovdens;    /* maximum overdensity possible with currrent AMR hierarchy */
  
  int     ioflag;        /* flag that indicates that we wrote an output file         */
  
#ifdef MOND
  /*========================================
   *                MOND
   *========================================*/
  double        max_gN;       /* maximum acceleration Newtonian                   */
  double        mean_gN;
  double        max_gM;       /* maximum acceleration MONDian                     */
  double        mean_gM;
  long unsigned steps;
  long unsigned no_ALLevents;
  long unsigned no_MONDevents;
#endif
} info_global;

/*=============================================================================
 * tline:  structure carrying the functions a(t), omega(t), and lambda(t)
 *=============================================================================*/
typedef struct tline
{
  double age[MAXTIME];
  double hubble[MAXTIME];
  double omega[MAXTIME];
  double lambda[MAXTIME];
  double virial[MAXTIME];
  double growth[MAXTIME];
  double a[MAXTIME];
  double t[MAXTIME];
  double super_t[MAXTIME];
  double rhoc[MAXTIME];
}tline, *tlptr;

/*=============================================================================
 * uparam:  structure holding the information provided by the user
 *=============================================================================*/
typedef struct uparam {
  char icfile_name   [MAXSTRING];
  char dumpfile_name [MAXSTRING];
  char outfile_prefix[MAXSTRING];
#ifdef LIGHTCONE
  char lightcone_prefix[MAXSTRING];
#endif /* LIGHTCONE */

  int     NGRID_DOM;
  
  double  Nth_dom;
  double  Nth_ref;
  
  double  final_z;
  
  int     out_dumps;
  int     no_outputs;
  double *z_out;
  
#ifdef LIGHTCONE
  double lightcone_z;
  int lightcone_type;
  double sa,sb,sg;		/* sky patch orientation angles */
  double dcx,dcy;		/* sky patch sizes              */
#endif /* LIGHTCONE */
  
#ifdef MOND
  double g0;  // cm/s^2
  double h0;  // km/sec/Mpc
#endif
  
#if (defined ART && defined HYDRO)
  /* the hydro-parameters are not stored in the ART binary */
  double omegab;
  double gamma;
  double H_frac;
  double T_init;
  double B_init;
#endif
} uparam, *uparamptr;
  
/*=============================================================================
 * simu:   structure carrying all sorts of information on the simulation
 *=============================================================================*/
typedef struct param_simu
{
   /* cosmological information */
   double        omega0;
   double        lambda0;
   tline         timeline;   

   
   /* timestepping information */
   double        a_initial;
   double        a_final;
   double        z_initial;
   double        z_final;
   double        t_initial;
   double        t_final;
   double        super_t_initial;
   double        super_t_final;
   
   
   /* unit stuff */
   double        mean_dens;
   double        FourPiG;       /* 4piG factor when using -DISOLATED    */
   double        boxsize;       /* the size of the cosmological box     */
   double        pmass;         /* mass of a single particle            */
   double        t_unit;        /* the time unit                        */

   double        l_unit;        /* currently boxsize and pmass are the units...          */
   double        m_unit;        /* ...but at some stage this should be clearly separated */

   
   /* HYDRO information */
   double        gamma;        /* adiabatic index for equation-of-state */
   double        omegab;       /* baryonic matter content               */
   double        omegaDM;      /* dark matter content                   */
   double        f_b;          /* baryon fraction omegab/omega0         */
   double        H_frac;       /* mass fraction of molecular hydrogen   */
   double        T_init;       /* initial gas temperature               */
   double        B_init;
   double        e_init;       /* initial internal energy density       */
   
   
   /* numerical details */
   double        no_vpart;      /* virtual number of particles          */
   long unsigned no_part;       /* real (physical) number of particles  */
   long unsigned no_species;    /* number of different particle species */
   long unsigned no_gas;        /* how many gas particles               */
   long unsigned no_stars;      /* how many star particls               */
   
   
   long unsigned no_halos;      /* number of halos (used by AHF)  */
   
   int           NGRID_DOM;
   int           NGRID_MIN;
   int           NGRID_MAX;
   double        Nth_dom;
   double        Nth_ref;
   double        SHIFT;
   double        min_weight;    /* minimum particle mass in internal units     */ 
   double        max_weight;    /* maximum particle mass in internal units     */
   double        med_weight;
   int           np_limit;

   int           mmfocus;
   int           multi_mass;
   int           double_precision;
   int           hydro;
   int           magneto;
   
#ifdef LIGHTCONE 
   /* the redshift limit of the lightcone */
   double        z_lightcone;
#endif /*LIGHTCONE*/

#ifdef MOND
  double        g0;  // cm/sec^2
  double        h0;  // km/Mpc/sec
#endif
   
} param_simu;

/*=============================================================================
 * files:   structure carrying all information on output/input files
 *=============================================================================*/
typedef struct info_io
{
  /*------------------------------------------------------------------------
   * this first block is filled by copying the appropriate user_data!
   *------------------------------------------------------------------------*/

   
   /*------------------------------------------------------------------------
    * information for managing input and output files
    *   - filenames
    *   - pointer to logfile for global access
    *   - no. of output files and their respective expansion factors
    *------------------------------------------------------------------------*/
   FILE         *logfile;
   
   char         *icfile_name;           /* name of file with IC's             */
   char         *outfile_prefix;        /* prefix for output files            */
   char         *dumpfile_name;         /* name of dump file                  */
   char         *logfile_name;          /* name for log file                  */
   
   int           no_outputs;            /* total number of outputs            */
   int           out_dumps;             /* when to write a dump file          */
   double       *a_out;                 /* output redshifts                   */

  
  
  
  /*-------------------------------------------------------------------------------------------------
   * here we store the particles that are being read in from the input file
   *
   * fst_part is an array of length no_part that holds >>all<< particles
   * (NOTE: a particle structure only carries pos, mom, and weight!)
   *
   * fst_gas is an array of length no_gas that holds >>additional<< information for gas particles
   * (NOTE: the mass/weight is stored in fst_part[offset_gas -> offset_gas+no_gas-1]
   *
   * fst_star is an array of length no_stars that holds >>additional<< information for stars
   * (NOTE: the mass/weight is stored in fst_part[offset_gas+offset_stars -> offset_gas+offset_gas+no_stars-1]
   *------------------------------------------------------------------------------------------------*/
  
  
  /* here we store the particles read in from file (no_part is stored in io.header, too!) */
  partptr       fst_part;      // fst_part gives access to >>all<< particles
  long unsigned no_part;       // = no_DMpart + no_gas + no_stars
  
  gasptr        fst_gas;       // additional information other than mass for gas particles
  long unsigned no_gas;
  long          offset_gas;    // offset in fst_part[] to access gas particles
  
  starptr       fst_star;      // additional information other than mass for star particles
  long unsigned no_stars;
  long          offset_stars;  // offset in fst_part[] to access gas particles
  
  
  /*----------------------------------------------------------------------------
   *                            AMIGA file header
   *
   * here we store information that is relevant for 
   *  - restarting a simulation
   *  - obtaining information about the status of the simulation
   *    at the time the output file was written
   *
   *
   * this header gets written to file and is being constantly extended
   *
   * therefore, the ordering of variables has no physical meaning and
   * is more a reflection of the development of AMIGA than anything else
   *
   * if you require downwards compatibility, please do not change the ordering
   *----------------------------------------------------------------------------*/
  struct io_header_block
    {
       char          header[HEADERSTRING];    // a string at your disposal ;-)
       
       int           multi_mass;
       int           double_precision;        // remember previous settings
       
       long unsigned no_part;                 // this is the total number of particles (=DM+gas+stars)
       long unsigned no_species;
       double        no_vpart;                // information about particles (cf. min/max_weight below!)
       
       double        timestep;
       int           no_timestep;             // stage of simulation and last timestep used
       
       double        boxsize;                 // [length unit]
       double        omega0;
       double        lambda0;                 // information about the cosmological model
       
       double        pmass;                   // [mass unit]
       
       double        cur_reflevel;
       double        cur_frcres;              // information about currently finest refinement level
       
       double        a_initial;               // for historical reasons we store the expansion factor a
       double        a_current;               // rather than supercomoving time in the input/output file
       
       double        K_initial;
       double        K_current;
       double        U_initial;
       double        U_current;
       double        Eintegral;
       double        Econst;                  // needed for layzer_irvine() energy check
       
       double        paramNSTEPS;
       double        paramNGRID_DOM;
       double        paramNth_dom;
       double        paramNth_ref;
       double        paramE_UPDATE;
       double        paramCELLFRAC_MAX;
       double        paramCELLFRAC_MIN;
       double        paramCA_CRIT;
       double        paramMAX_L1DIM;
       double	      paramDOMSWEEPS;
       double        paramREFSWEEPS;          // technical details
       
       double        paramAHF_MINPART;
       double		   paramAHF_VTUNE;
       double        paramAHF_RISE;
       double        paramAHF_SLOPE;
       double        paramAHF_MAXNRISE;       // technical details about AHF
       
       double        min_weight;
       double        max_weight;
       double        t_unit;
       double        B_init;                  // relevant for MHD solver
       double        param_dummy5;
       double        param_dummy6;
       double        param_dummy7;
       double        param_dummy8;            // empty space that once was used (downwards compatibility!)
       
       double        version;
       int           built;                   // the code is under constant development ;-)
       
       double        omegab;
       double        gamma;
       double        H_frac;
       double        T_init;                  // relevant for the HYDRO solver
       
       int           hydro;
       int           magneto;
      
       double        med_weight;
       
       char          dummy[FILLHEADER];       // empty space for future additions (e.g. MHD, ...)
    } header;
  
#ifdef LIGHTCONE 
  /* the lightcone file lives many timesteps, similar to the logfile */
  FILE         *lightcone;
  char		   *lightcone_prefix;
  /* start of the backup list */
  bckptr fst_backup_part;
  /* start of the cubes list */
  cubeptr fst_cube;
  /* a place to find the previous lightcone radius and z */
  double rcone0;
  double z0;
  /* type of lightcone	(see output_lc.c) */
  int conetype;
  /* geometry data for smaller patches */
  double sa,sb,sg;	/* orientation angles                                  */
  double dcx,dcy;	/* patch half-sizes                                    */
  double xcoef,ycoef;	/* used to check if a point is in the patch (sin^2dcx) */
  double R[3][3];	/* coordinate transformation matrix                    */
  double kpatch[4][3];	/* direction vectors for patch pyramid edges           */
  double npatch[4][3];	/* normals to patch pyramid faces                      */
#endif /*LIGHTCONE*/
} info_io;

  





/*======================================
 * definitions used by AHF
 *======================================*/
typedef struct {
 double x,y,z;
} XYZ;

typedef struct {
 double min,max;
} MINMAX;

#ifdef TRACKER
/* needed for halo tracker: */
typedef struct {
   int    id;
   int    species;
   double mass;
} HALOPART;
#endif

/* the HALOPROFILE structure */
typedef struct {
 int             nbins;
 double         *r;
 long unsigned  *npart;
 double         *nvpart;
 double         *ovdens;
 double         *dens;
 double         *v2_circ;
 double         *sig_v;
 double         *Ekin;
 double         *Epot;
 double         *Lx;
 double         *Ly;
 double         *Lz;
 double         *axis1;
 double         *E1x;
 double         *E1y;
 double         *E1z;
 double         *axis2;
 double         *E2x;
 double         *E2y;
 double         *E2z;
 double         *axis3;
 double         *E3x;
 double         *E3y;
 double         *E3z;
#ifdef AHFphspdens
 double *sigma2_vx_sh;
 double *sigma2_vy_sh;
 double *sigma2_vz_sh;
 double *sigma2_vr_sh;
 double *sigma2_vtheta_sh;
 double *sigma2_vphi_sh;
#ifdef AHFmeanvelocities
 double *mean_vx_sh;
 double *mean_vy_sh;
 double *mean_vz_sh;
 double *mean_vr_sh;
 double *mean_vtheta_sh;
 double *mean_vphi_sh;
 double *mean_vx_sp;
 double *mean_vy_sp;
 double *mean_vz_sp;
 double *mean_vr_sp;
 double *mean_vtheta_sp;
 double *mean_vphi_sp;
#endif
#endif
#ifdef GAS_PARTICLES
  double *nvpart_gas;
  double *nvpart_star;
#endif
} HALOPROFILE;


/* the HALO structure */
typedef struct {   
 /*--------------------------------
  * access to the halo's particles
  *--------------------------------*/
 long unsigned  npart;       /* number of particles                              */
 long unsigned *ipart;       /* indizes of halo particles                        */
 long unsigned  nll;         /* temporary storage of no. of parts attached to ll */
 partptr        ll;          /* spatialRef2halos() still requires linked-lists   */
 long unsigned  npart_max;   /* maximum number of particles: used to avoid excessive calls to realloc() */

#ifdef TRACKER
 #if DONT_REORDER
 /* store all particle ids of this halo, for finding them in each snapshot */ 
 /* also need species type and mass for tracking gas and stars properly */
 HALOPART *part;
  #endif

 #ifdef DEBRIS
 long unsigned  nopart;       /* number of outliers                               */
 long unsigned *iopart;       /* indizes of halo particles                        */

 long unsigned  nupart;       /* number of unbound particles                      */
 long unsigned *iupart;       /* indizes of halo particles                        */
 #endif 
#endif
 
 /*----------------------------
  * integral properties
  *----------------------------*/
 XYZ     pos;
 XYZ     vel;
 double  M_vir;
 double  R_vir;
 double  velDis;
 double  lambda;
 double  lambdaE;
 double  R_max;
 double  V2_max;
 double  ovdens;
 double  R_edge;
 double  Ekin;
 double  Epot;
 double  Phi0;
 XYZ     axis;
 XYZ     E1;
 XYZ     E2;
 XYZ     E3;
 XYZ     AngMom;
 
 double  com_offset;   /* offset of centre-of-mass     to halo centre     */
 double  mbp_offset;   /* offset of most-bound-particle to halo centre    */
 double  r2;           /* position in density profile where rho*r^2 peaks */
 
 
#ifdef GAS_PARTICLES
 struct  {
   long unsigned  npart;
   double         M_vir;
   double         lambda;
   double         lambdaE;
   double         Ekin;
   double         Epot;
   XYZ            axis;
   XYZ            E1;
   XYZ            E2;
   XYZ            E3;
   XYZ            AngMom;
 } gas;

  struct  {
    long unsigned  npart;
    double         M_vir;
    double         lambda;
    double         lambdaE;
    double         Ekin;
    double         Epot;
    XYZ            axis;
    XYZ            E1;
    XYZ            E2;
    XYZ            E3;
    XYZ            AngMom;
  } star;
#endif
  
  
 /*----------------------------
  * the halo profile
  *----------------------------*/
 HALOPROFILE prof;
 
 
 /*----------------------------
  * ahf_gridinfo() information
  *----------------------------*/
 int oldIndex;          /* The halos old index                                               */
 
 int hostHaloGrids;     /* Index to the halo that hosts this one based on the grids          */
 int hostHalo;          /* Index to the halo that hosts this one based on the gathering radius' */
 
 int numHostHaloGrids;  /* Additional host grids that might contain particles for this halo  */
 int *hostGrids;        /* Additional host grids that might contain particles for this halo  */
 
 double gatherRad;      /* The radius which we gather all particles to this halo from the two hosts */
 
 /* NOTE :: Substructure is not implemented correctly yet                                    */	
 int  numSubStruct;		/* Number of sub halos                                               */
 int	*subStruct;	  	   /* Indexes of the substructure halos                                 */
 
 double spaRes;         /* The spatial resolution of the halo                                */			
 int    refLev;         /* The deepest refinement level:  0 - is the what we start AHF on    */			
 
 int    numNodes;      /* Number of nodes in this halo                                       */
 
#	if (defined WITH_MPI || defined AHFrestart)
 boolean ignoreme;
#	endif
 
} HALO;
  

typedef struct {
   int x,y,z;
} intXYZ;

/* An index array to save computation */
typedef struct {
   int    refLevel;
   int    isoRefIndex;
   intXYZ periodic;
} SRINDEX;


/*=============================================================================
* in ahf we are going to store some global information about AHF
*=============================================================================*/
typedef struct info_ahf
{
   int       no_grids;    /* Number of grids used in AHF    */
   int       min_ref;     /* coarsest grid to be considered */
   time_t    time;
} info_ahf;


/*======================================================================*/
/* particle coordinates at the previous timestep (the backup structure) */
/*                      (used with -DLIGHTCONE)                         */
/*======================================================================*/
typedef struct oldpart *bckptr;
typedef struct oldpart {
#ifdef DOUBLE
   dvect   pos;
#else
   fvect   pos;                	/* position vector             */
#endif
}oldparticle;


/*==============================================================*/
/* the periodic cubes list for the lightcone                    */
/*              (used with -DLIGHTCONE)                         */
/*==============================================================*/
typedef struct cubes *cubeptr;
typedef struct cubes {
   int I,J,K;
   int used;
} cube;


/*=============================================================================
* energy:   structure carrying all information on energy conservation
*=============================================================================*/
typedef struct info_layzer_irvine
{
  double  K_initial, K_current, Kold;
  double  U_initial, U_current, Uold;
  double  aold;
  double  integral;
  double  econst;
  double  echeck;

  time_t  time;
} info_energy;

/*=============================================================================
 * powerspectrum:   structure carrying all information about on-the-fly P(k)
 *=============================================================================*/
typedef struct info_Pk
{
   double    *rk;
   double    *Pk;

   double    Pk_ini;
   double    Pk_now;
   
   double    Dgrowth_ini;
   double    Dgrowth_now;
   
   int       dump_Pk;
      
   time_t    time;
} info_Pk;

/*=============================================================================
 * info_gadget:   structure carrying all information about GADGET files
 *=============================================================================*/
#define SIZEOFGADGETHEADER 256    /* Size of GADGET header in bytes */
typedef struct info_gadget
{
  int        no_gadget_files;
  int        i_gadget_file;
  int      *(np[6]);
  int        nall1;
  partptr    lst_part;
  gasptr     lst_gas;
  
#ifdef LGADGET
  long long  IDmin;
  long long  IDmax;
#else
  int        IDmin;
  int        IDmax;
#endif
  double     mmin;
  
  gasptr     fst_gas;
  
  starptr    fst_star;
  
  struct io_gadget_header
    {
      int      np[6];
      double   massarr[6];
      double   expansion;
      double   redshift;
      int      flagsfr;
      int      flagfeedback;
      int      nall[6];
      int      flagcooling;
      int      NumFiles;
      double   BoxSize;
      double   Omega0;
      double   OmegaLambda;
      double   HubbleParam;
      char     unused[SIZEOFGADGETHEADER- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];
    } header; 
  
} info_gadget;

  /*==============================================================*/
  /* Variables for reading in the TIPSY format */
  /*===============================================================*/
#ifdef TIPSY
  double tipsy_boxsize,tipsy_omega0,tipsy_lambda0,tipsy_initalz,tipsy_currentimeno;
#define MAXDIM 3
#define forever for(;;)
  
  typedef float Real;
  
  struct gas_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real rho;
    Real temp;
    Real hsmooth;
    Real metals ;
    Real phi ;
  } ;
  
  struct gas_particle *gas_particles;
  
  struct dark_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real eps;
    Real phi ;
  } ;
  
  struct dark_particle *dark_particles;
  
  struct star_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real metals ;
    Real tform ;
    Real eps;
    Real phi ;
  } ;
  
  struct star_particle *star_particles;
  
  struct tipsy_dump {
    double time;
    int    nbodies;
    int    ndim;
    int    nsph;
    int    ndark;
    int    nstar;
    
    /* Jeremy Balin addition */
    char align[ (32 - sizeof(double) - 5 * sizeof(int)) / sizeof(char) ];
    /* total size should be 32 bytes to make alignment okay with
     * 64-bit architectures (ie. alphas) */
    
  } ;
#endif /* TIPSY */
  
  
#endif

