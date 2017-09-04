#include "../param.h"
#include "../tdef.h"

#ifndef IO_LOGGING_INCLUDED
#define IO_LOGGING_INCLUDED

/* routine that controls the io-logging */
void update_logfile         (double timecounter, double timestep, int no_timestep);
void log_structures  ();
void write_parameterfile    ();

#endif
