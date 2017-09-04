#ifndef PARAM_INCLUDED
#define PARAM_INCLUDED

#include <float.h>

/*================================================================================
 * switch on/off various features of AMIGA not passed via DEFINEFLAGS in Makefile
 *================================================================================*/
#include "define.h"

/*=============================================================================
 * general parameters
 *=============================================================================*/
#define TERMINATE_AMIGA {"terminateAMIGA"}

#define E_UPDATE     10      /* how often to calculate Layzer-Irvine-Check  10  */

#define PKMODE       2       /* 1=fundamental mode, NGRID_DOM/2=Nyquist mode    */
#define PKGROWLIMIT  5.0     /* PKMODE/Dgrowth >PKGROWLIMIT => AMIGA terminates */

#define MAXTIME      10000   /* interpolation points for timeline (cf. tdef.h)   */

#ifdef MPI_DEBUG
#undef MAXTIME
#define MAXTIME 100
#endif

/** Number of bits used per dimension for the calculation of the Hilbert key */
#define BITS_PER_DIMENSION 21

/** Defines at which level the loadbalancing is done */
//#define LOADBALANCE_DOMAIN_LEVEL 3 /*  8^3 grid */
//#define LOADBALANCE_DOMAIN_LEVEL 4 /*  16^3 grid */
//#define LOADBALANCE_DOMAIN_LEVEL 5 /*  32^3 grid */
#define LOADBALANCE_DOMAIN_LEVEL 6 /* 64^3 grid */
//#define LOADBALANCE_DOMAIN_LEVEL 7 /* 128^3 grid */
//#define LOADBALANCE_DOMAIN_LEVEL 8 /* 256^6 grid */

/** The maximum amount of particles that can be send in on flush */
#define MAX_SEND_PARTICLES 1000000

/** Defines the verbosity level for the io_logging_* functions if
 * NEWSTARTUN is used. Depending on this value some messages might not
 * appear and hence this can be used to reduce the chatter produced by
 * the io_logging_* function. During the starup this will be passed to
 * the logging module (in startrun.c). The lower the number, the less
 * output will be produced. */
#define VERBOSITY 6

/*============================================================================
 * grid parameters
 *============================================================================*/
#ifndef FFT
#define MIN_L1DIM 2               /* minimum domain grid size (1D) */
#endif
#define MAX_L1DIM  16777216       /* maximum domain grid size (1D) */
#define MIN_NNODES 125            /* smallest grid-block 5x5x5     */

/*============================================================================
 * time stepping
 *============================================================================*/
#define NSTEPS       1000    /* number of (initial) steps for time stepping     */
#define CA_CRIT      0.15    /* restricts timestep due to da/a criterion    0.05*/ 
#define CELLFRAC_MAX 0.2     /* how far are particles allowed to move       0.2*/
#define CELLFRAC_MIN 0.05    /* how far should particles move at least      0.05*/
#define CF_MEAN      ((CELLFRAC_MAX+CELLFRAC_MIN)/2)
// the hydro CFL_TUNE parameter is defined below!


/*============================================================================
 * potential solver
 *============================================================================*/
#define DOMSWEEPS  10     /* number of domain and refinement sweeps...   10   */
#define REFSWEEPS  10     /* ...before checking for convergence          10   */
#define W_SOR      1.34   /* successive over-relaxation parameter        1.34 */

#define ETA        0.625  /* parameter for slow convergence              0.625*/  
#define CONVCRIT   0.1    /* convergence criterion                       0.1  */
#define DOMCORRECT 0.25   /* constrained criterion on domain grid        0.25 */

#define SPEEDFRAC  0.01   /* fraction of particles allowed to speed      0.01 */


/*=============================================================================
* AHF related parameters
*=============================================================================*/
#define AHF_MINPART        20     /* minimum number of particles per halo  > 10       */
#define AHF_VTUNE          1.5    /* allow particle velocities VTUNE*v_escape [1,5]   */
#define AHF_MAX_GATHER_RAD 2.0    /* maximal radius for halos in Mpc/h                */
#define AHF_MIN_REF_OFFSET 0      /* offset for first refinement to be used by AHF    */
#define AHF_MASSMIX        0.0    /* how much mass must be in high-res particles      */
#define AHF_MAXHALO        5000   /* most massive halo is heavier -> AMIGA terminates */
#define AHF_RISE           1.00   /* Rho > AHF_RISE*Rho_prev -> rising density        */
#define AHF_SLOPE          0.99   /* outer halo profile at least like r^-AHF_SLOPE    */
#define AHF_MAXNRISE       2      /* try to catch variations in density               */

#define PGAS                0.0   /* identifier for gas particles; has to be exactly 0.0!!!! */
#define PDM                -1.0   /* identifier for dm particles; whatever negative value */
#define PSTAR              -4.0   /* identifier for star particles; whatever negative value */
#define PDMbndry           -5.0   /* identifier for star particles; whatever negative value */

/*=============================================================================
 * Halo tracker related parameters
 *=============================================================================*/
#define TRK_MINPART  5     /* minimum number of particles per halo, usually < AHF_MINPART  */
#define TRK_VTUNE    1.5	/* allow modified v_esc in rem_unbound */

/*=============================================================================
 * GADGET related parameters
 *=============================================================================*/
#define GADGET_MUNIT 1.0e10    /* GADGET mass unit in Msol/h                   */
#ifdef GADGET_LUNIT_KPC
#define GADGET_LUNIT 1.0e-3
#else                          /* GADGET length unit in Mpc/h                  */
#define GADGET_LUNIT 1.0
#endif


#ifndef NEWSTARTRUN
/*=============================================================================
 * TIPSY related parameters
 * NOTE: THESE MUST BE CHANGED TO SUIT YOUR PARTICULAR SIMULATION!!!!!!!!!
 *=============================================================================*/
#define	TIPSY_VUNIT     690.988298942671  /* Tipsy vel unit in km/s [FROM  H0*Lbox/sqrt(8.*!PI/3.)] */
#define	TIPSY_MUNIT     2.2197e15         /* Tipsy mass unit in Msol/h */
#endif /* NEWSTARTRUN */


/*=============================================================================
 * LIGHTCONE related parameters
 *=============================================================================*/
#define D2R           (PI/180.)     /* needed for LIGHTCONE patch trigonometry */

/*============================================================================
 * handling output file names etc.
 *============================================================================*/
#define MAXSTRING    2048   /* used for char statement, i.e. filenames etc.     */

#define AMIGAHEADER  2048   /* maximum size (in bytes) for output file header   */
#define HEADERSTRING 256    /* no. of characters for header string in outfiles  */
#define HEADERSIZE  (HEADERSTRING*sizeof(char)+2*sizeof(long)+6*sizeof(int)+46*sizeof(double))
#define FILLHEADER  (AMIGAHEADER-HEADERSIZE)

/*=============================================================================
 * finally some convenient abreviations...
 *=============================================================================*/
#define NDIM      3          /* DO NOT EVER TOUCH THIS NUMBER!  */
#define CRITMULTI 8.0
#define NP_RATIO  7.5
#ifdef DOUBLE 
#define ZERO         (1E-12)
#define MACHINE_ZERO (1E-20)
#else
#define ZERO         (1e-6)
#define MACHINE_ZERO (1e-10)
#endif
#define MZERO        (1e-10)  /* used when reading GADGET files */

/*=============================================================================
 * hydro stuff
 *=============================================================================*/
#define CFL_TUNE         0.3
#define ETA_DUAL_ENERGY1 0.001
#define ETA_DUAL_ENERGY2 0.3

#define NHYDRO    7          /* number of hydrodynamical variables u = (rho,mvx,mvy,mvz,E,S,e)   */
#define NADVECT   (NHYDRO-1) /* we do not advect edens which is stored >>last<< in u[]           */

#define Udens     0          /* positions in vector u ...                                   */
#define UmomdensX 1
#define UmomdensY 2
#define UmomdensZ 3          /* it is important to consecutiuvely store Umomdens!           */
#define UEdens    4
#define Untrpy    5
#define Uedens    6          /* Uedens >>must<< be the last entry as it is >>not<< advected */

/*=============================================================================
 * some physical constants...
 *=============================================================================*/
#define Gyr       3.1558e16         /* [sec]                   */
#define Mpc       3.08567782e19     /* [km]                    */
#define H0        100.              /* [h*km]/[sec*Mpc]        */
#define rhoc0     2.7755397e11      /* [h^2*Msun]/[Mpc^3]      */
#define Grav      4.3006485e-9      /* [Mpc*km^2]/[Msun*sec^2] */
#define cH0	      2998.0		      /* c/H0 (in h^-1 Mpc)      */
#define kB_per_mp 0.825481286614E-2 /* [(km/sec)^2/K]          */
#define kBoltzman 6.9416792         /* [(km/sec)^2 Msun/K]     */
#define Msun      1.9891e30         /* [kg]                    */

#define bytes2GB  9.313225746154785e-10

/*============================================================================= 
 * no code is complete without defining these numbers ;-)
 *=============================================================================*/
#define PI    3.14159265358979323846264338
#define TWOPI 6.28318530717958647692528677
#define SQRT2 1.41421356237309504880168872

/*=============================================================================
 * The numbers 0, 1 and 2 are used to represent the X, Y and Z coordinates of
 * a three dimensional object (eg vector). To aid readability these numbers are
 * replaced by the names X, Y and Z when used individually.
 *=============================================================================*/
#define X 0     /* x-coord symbol */
#define Y 1     /* y-coord symbol */
#define Z 2     /* z-coord symbol */

/*=============================================================================
 * Boolean parameters
 *=============================================================================*/
#define YES   1
#define NO    0

#define ON    1
#define OFF   0

#define TRUE  1
#define FALSE 0

#define SCDM  0
#define OCDM  1
#define LCDM  2

/*=============================================================================
 * this is written into the logfile just for information
 *=============================================================================*/
#define VERSION 0.0
#define BUILT   271

#endif

