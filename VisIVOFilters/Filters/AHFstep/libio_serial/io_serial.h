#include "../param.h"
#include "../tdef.h"

#ifndef IO_SERIAL_INCLUDED
#define IO_SERIAL_INCLUDED

/* these are the major routines for writing/reading output files in BINARY */
void input            (char *infile_name);
void output           (char *outfile_name, int dumpflag);
void input_grid       (gridls *cur_grid);
void output_grid      (gridls *cur_grid, int dumpflag);
void manage_outputs   (double timecounter, double timestep);

/* wrapper for the actual reading routine */
void read_data        (FILE *icfile);

/* ...the actual reading routines */
void read_amiga       (FILE *icfile);
void read_art         (FILE *icfile);
void read_mlapm       (FILE *icfile);
void read_tipsy       (FILE *icfile);
void read_mare_nostrum(FILE *icfile);
void read_deva        (FILE *icfile);
void read_gadget      (FILE *icfile);
void skim_gadget      (FILE *icfile);


/* routine that controls the io-logging */
#include "io_logging.h"

/* these routines serve mainly debuggin purposes and hence write information in ASCII */
#include "write.h"

/* all those LIGHTCONE routines */
#include "output_lc.h"


#endif
