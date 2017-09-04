#ifndef DEFINE_INCLUDED
#define DEFINE_INCLUDED

/*=============================================================================
 * here we switch on/off various features of AMIGA (i.e. DEFINEFLAGS)
 * (#define statements without actually defining a value...)
 *=============================================================================*/

/*-------------------------------------------------- 
 *                  -DSTANDARD
 *--------------------------------------------------*/
#ifdef STANDARD

#define ATS             /* (time-)adaptive time stepping                        */
#define COSMOLOGY       /* write Cosmology.DAT file                             */
#define FFT             /* use FFT on domain grid                               */
#define ENERGY          /* monitor energy conservation                          */
#define REDSHIFTNAME    /* use redshift in output file names                    */
#define PERIODIC        /* use periodic boundary conditions                     */

#endif /* STANDARD */

/*--------------------------------------------------
 *                  -DAHFstep
 *--------------------------------------------------*/
#ifdef AHFstep
#undef  ENERGY       /* switched off to safe start-up time  */
#undef  COSMOLOGY    /* switched off to safe start-up time  */
#ifndef VERBOSE
#define VERBOSE
#endif
#endif

/*--------------------------------------------------
 *                   -DAHF
 *--------------------------------------------------*/
#ifdef AHF

/* the following flags decide how to make Parent-Daughter Assignment in analyseRef() */
//#define PARDAU_DISTANCE  /* use sub-grid with closest distance to follow host */
//#define PARDAU_NODES     /* use sub-grid with most nodes to follow host       */
#define PARDAU_PARTS     /* use sub-grid with most particles to follow host   */


/* the following centre flags determine how to calculate the halo centre... */
/* ...the default (i.e. no flag used!) is a density weighted centre         */
//#define AHFmaxdenscentre /* use cell with maximum density as halo centre                       */
//#define AHFpotcentre     /* use potential weighted centre as halo centre                       */
//#define AHFgeomcentre    /* use geometrical centre as halo centre                              */
#define AHFcomcentre     /* use centre-of-mass of particle on finest refinement as halo centre */


/* miscellanesou flags */
//#define AHFmaxhalo
//#define AHFmmfocus
//#define AHFnoremunbound
//#define MANUAL_DVIR
//#define AHFreducedinertiatensor
//#define AHFcentrefile
//#define AHFabsangmom     /* dump the un-normalized angular momentum to file                */
//#define AHFsplinefit     /* this will use a spline-interpolation to get Rmax, Vmax, and r2 */
//#define AHFprofilerise   /* checks for rising profile in rem_outsideRvir() and chops halo  */

/* the following flags decide whether or not to dump the *.AHF_substructure file    */
#define AHFsubstructure      /* dump substructure information to file             */
#define AHFgridsubstructure  /* ...base that information upon grid hierarchy      */

/* make it standard that the .AHF_particles file contains information about the particle type */
#define PARTICLES_INFO

/* the following flags are mainly used for debugging purposes */
//#define AHFgeom               /* writes .AHF_halos.geom files                      */
//#define AHFgridinfofile       /* writes grid information into a new file           */

#endif


/*--------------------------------------------------
 *                      TRACKER
 *--------------------------------------------------*/
/* rename DEBRIS flag to TRACKER_DEBRIS for more transparency */
#ifdef DEBRIS
 #define TRACKER_DEBRIS
#endif

/* ensure, that flag only is set for TRACKER */
#ifndef TRACKER
 #undef TRACKER_DEBRIS
#endif


/*--------------------------------------------------
 *                   -DHYDRO
 *--------------------------------------------------*/
#ifdef HYDRO
#define DOUBLE
#define MONOATOMIC

/* dual energy formalism */
#define DUAL_ENERGY
#define RYU_CRITERION

/* here you decide what slope limiter to use */
#define SLOPE_vanLeer
//#define SLOPE_MINMOD
//#define SLOPE_SUPERBEE

#define RECONSTRUCT_VELOCITIES       // this is a whole lot more stable for the KNP solver
//#define RECONSTRUCT_PRESSURE
//#define FLAT_RECONSTRUCTION
//#define CHECK_RECONSTRUCTION

/* sets all gas variables to zero */
//#define NO_GAS


/* define the actual hydro test */
/* ---------------------------- */
#ifdef HYDRO_TEST
#undef COSMOLOGY

#if ((HYDRO_TEST>=0 && HYDRO_TEST<4))
#undef  PERIODIC
#undef  DUAL_ENERGY
#define NO_HYDRO_SOURCES
#endif

#if (HYDRO_TEST<3)
#define SHOCK_TUBE
#define SHOCK_RICKER
#undef  RYU_CRITERION
#endif

#if (HYDRO_TEST==3)
#define BLAST_WAVE
#define BLAST_RICKER00
#undef RYU_CRITERION
#endif

#if (HYDRO_TEST==9 || HYDRO_TEST==10 || HYDRO_TEST==11)   // HYDRO_TESTs 9-11: 1D MHD tests with left/right states and without gravity or cosmology
#undef PERIODIC                                           // 9: in x-direction; 10: in y-direction; 11: in z-direction
#define NO_HYDRO_SOURCES
#define SHOCK_TUBE
#define SHOCK_MHD
#undef RYU_CRITERION
#undef DUAL_ENERGY      // to compare better with Ziegler 2004. He does not use dual energy for the 1D tests.
#endif

#if (HYDRO_TEST==9)
#define PERIODIC_Y
#define PERIODIC_Z
#endif

#if (HYDRO_TEST==10)
#define PERIODIC_X
#define PERIODIC_Z
#endif

#if (HYDRO_TEST==11)
#define PERIODIC_X
#define PERIODIC_Y
#endif

#if (HYDRO_TEST==12) //Orszag-Tang vortex MHD test
#define NO_HYDRO_SOURCES
#define ORSZAG_TANG
#undef RYU_CRITERION
#undef DUAL_ENERGY  // TODO: try it with dual energy as well
#endif

#else /* HYDRO_TEST */

#define HYDRO_TEST 10000

#endif /* HYDRO_TEST */

#endif /* HYDRO */

//--------------------------------------------------//
//                     -DMHD                        //
//--------------------------------------------------//
#ifdef MHD

// which constrained transport scheme shall we use?
#define CT_GARDINER_STONE
//#define CT_ZIEGLER1D
//#define CT_ZIEGLER2D
//#define CT_ZIEGLER3D

#endif // MHD

/*--------------------------------------------------
 *                  -DWITH_MPI
 *--------------------------------------------------*/
#ifdef WITH_MPI
#	ifndef NEWSTARTRUN
#		define NEWSTARTRUN
#	endif
#  undef  REF_TEST
#  undef  VERBOSE
#endif


/*--------------------------------------------------
 *                    -DGADGET
 *--------------------------------------------------*/
#ifdef GADGET2
#define GADGET
#endif
#ifdef GADGET
 #ifndef MULTIMASS
  #define MULTIMASS	/* avoid unnecessary compiler warnings */
 #endif
#define GAS_PARTICLES
#endif

/*--------------------------------------------------
 *                    -DTIPSY
 *--------------------------------------------------*/
#ifdef TIPSY
#define  MULTIMASS
#define  GAS_PARTICLES
#endif

/*--------------------------------------------------
 *                    -DDEVA
 *--------------------------------------------------*/
#ifdef DEVA2
#define DEVA
//#define DEVA2_QHULL_FILE "../snapshots/cube00.qhull"
#endif

#ifdef DEVA
#define MULTIMASS
#define GAS_PARTICLES
#endif

/*--------------------------------------------------
 *                -DMARE_NOSTRUM
 *--------------------------------------------------*/
#ifdef MARE_NOSTRUM
#define MULTIMASS
#define GAS_PARTICLES
#endif

/*----------------------------------------------------------------------------
 * -DISOLATED switches off all cosmological expansion and periodic boundaries
 *                     !!!!! WARNING: UNTESTED !!!!!
 *----------------------------------------------------------------------------*/
#ifdef ISOLATED
#undef  REDSHIFTNAME
#undef  COSMOLOGY
#undef  PERIODIC
#endif

/*----------------------------------------------------------------------------
 *  misc definitions
 *----------------------------------------------------------------------------*/
#ifdef NO_GAS
#undef GAS_PARTICLES
#endif

#ifdef VERBOSE2
#define VERBOSE
#endif

#ifdef VERBOSE
#define REF_TEST
#endif

#ifdef CONVERT_TERM
#define CONVERT
#undef AHFmmfocus
#undef COSMOLOGY
#endif

#ifdef OUTDUMPS
#undef REDSHIFTNAME
#endif


#ifdef FORCETEST
#define VERBOSE
#undef COSMOLOGY
#endif

#ifndef TSC
#ifndef CIC   /* forgotten to define mass assignemnt scheme ? => use TSC then... */
#ifndef NGP
#define TSC
#endif
#endif
#endif

/*--------------------------------------------
 * more transparent to read in source-code... 
 *--------------------------------------------*/
#ifndef PM
#define ADAPTIVE
#endif /* PM */

#ifndef CONTINUE
#define TERMINATE
#define TERMINATE2  /* used in leavers.c */
#define VERBOSELOG
#endif /* CONTINUE */

#ifdef PERIODIC
#define PERIODIC_X
#define PERIODIC_Y
#define PERIODIC_Z
#endif /* PERIODIC */

#endif

