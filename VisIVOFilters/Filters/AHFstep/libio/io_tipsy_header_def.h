#ifndef IO_TIPSY_HEADER_DEF_H 
#define IO_TIPSY_HEADER_DEF_H

/**
 * \file io_tipsy_header_def.h
 *
 * Provides the structure definition for the TIPSY header
 * structure. Including useful typedefs.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <inttypes.h>


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              * 
\**********************************************************************/

/** 
 * The size (in bytes) reserved at the beginning of a file for
 * the header 
 */
#define TIPSY_HEADER_SIZE 32

/** 
 * We actually store additional information in the header structure 
 */
#define TIPSY_HEADER_EXTRA 40

/**
 * The header structure itself
 */
struct io_tipsy_header_struct {
  double time;
  int    nbodies;
  int    ndim;
  int    nsph;
  int    ndark;
  int    nstar;
  int    pad;  
  double  TIPSY_OMEGA0;
  double  TIPSY_LAMBDA0;
  double  TIPSY_BOXSIZE;
  double  TIPSY_VUNIT;
  double  TIPSY_MUNIT;
};

/** Convenient typedef */
typedef struct io_tipsy_header_struct io_tipsy_header_struct_t;

/** Convenient typedef */
typedef io_tipsy_header_struct_t *io_tipsy_header_t;


#endif /* IO_TIPSY_HEADER_DEF_H */
