#ifndef IO_LOGGING_DEFS_H
#define IO_LOGGING_DEFS_H

/* $Id: io_logging_defs.h,v 1.5 2007/02/11 19:42:42 knolli Exp $ */

/**
 * \file io_logging_defs.h
 *
 * Defines a bunch of things than can be written to the logfile. Only
 * needed in io_logging.c.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/


/**********************************************************************\
 *    Global defines, structure definitions and typedefs              * 
\**********************************************************************/

/**
 * Use that as a seperator of different entries (when reusing a
 * logfile).
 */
#define IO_LOGGING_NEW_ENTRY "\n\n--- NEW ENTRY STARTS HERE ---\n\n"

/**
 * Write that as the last message when stopping to log.
 */
#define IO_LOGGING_STOP "\nLogging stopped.\n\n"

#ifdef WITH_MPI
	/**
	 * The P-AMIGA logo.
	 */
#	define IO_LOGGING_LOGO \
		"\tPPPP           AA      M   M    I    GGGG      AA  \n" \
		"\tP  PP        A   A    MM MM    I    G        A   A \n" \
		"\tPPPP  ===   AAAAA    M M M    I    G GGG    AAAAA  \n" \
		"\tP          A   A    M   M    I    G   G    A   A   \n" \
		"\tP         A   A    M   M    I     GGG     A   A    "
#else
	/**
	 * The AMIGA logo.
	 */
#	define IO_LOGGING_LOGO \
		"\t     AA      M   M    I    GGGG      AA  \n" \
		"\t   A   A    MM MM    I    G        A   A \n" \
		"\t  AAAAA    M M M    I    G GGG    AAAAA  \n" \
		"\t A   A    M   M    I    G   G    A   A   \n" \
		"\tA   A    M   M    I     GGG     A   A    "
#endif /* WITH_MPI */

/**
 * A horizontal bar.
 */
#define IO_LOGGING_BAR \
"========================================================================"

#define IO_LOGGING_BAR2 \
"------------------------------------------------------------------------"

#define IO_LOGGING_PART_PRE "\n\n\n\n\n\n" IO_LOGGING_BAR "\n\n"
#define IO_LOGGING_PART_POST "\n\n" IO_LOGGING_BAR "\n\n"
#define IO_LOGGING_PART_NUM_PRE "\t\tPART "
#define IO_LOGGING_PART_NUM_POST "\n\t"

#define IO_LOGGING_SECTION_PRE "\n\n" IO_LOGGING_BAR "\n"
#define IO_LOGGING_SECTION_POST "\n" IO_LOGGING_BAR "\n\n"
#define IO_LOGGING_SECTION_NUM_PRE "  "
#define IO_LOGGING_SECTION_NUM_POST "  "

#define IO_LOGGING_SUBSECTION_PRE "\n" IO_LOGGING_BAR2 "\n"
#define IO_LOGGING_SUBSECTION_POST "\n" IO_LOGGING_BAR2 "\n\n"
#define IO_LOGGING_SUBSECTION_NUM_PRE "  "
#define IO_LOGGING_SUBSECTION_NUM_POST "  "


#endif /* IO_LOGGING_DEFS_H */
