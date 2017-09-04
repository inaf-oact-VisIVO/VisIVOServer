/* $Id: io_util.c,v 1.6 2007/11/01 09:23:49 knolli Exp $ */

/**
 * \file io_util.c
 *
 * Implements functions to perfom input/output operations.
 */


/**********************************************************************\
 *    Includes                                                        * 
\**********************************************************************/
#include <stdlib.h>
#include <string.h>
#ifdef WITH_MPI
#	include <mpi.h>
#endif

#include "io_util.h"


/**********************************************************************\
 *    Local defines, structure definitions and typedefs               * 
\**********************************************************************/


/**********************************************************************\
 *    Prototypes of local functions                                   * 
\**********************************************************************/


/**********************************************************************\
 *    Implementation of global functions                              * 
\**********************************************************************/
extern void
io_util_sexchange(void *p, size_t s)
{
	int n;
	unsigned char ptmp,*pc;

	pc = (unsigned char *)p;

	for(n=0; n < s/2; n++) {
		ptmp = pc[n];
		pc[n] = pc[s - n - 1];
		pc[s - n - 1] = ptmp;
	}

	return;
}

extern int
io_util_readbytes(FILE *f, void *n, size_t len, io_file_swap_t swap)
{
	if (fread(n, len, 1, f) != 1)
		return 0;

	if (swap == IO_FILE_IS_SWAPPED)
		io_util_sexchange(n, len);

	return 1;
}

extern size_t
io_util_readstring(FILE *f, char *s, size_t n)
{
	int c;
	size_t i;

	s[0] = '\0';
	for (i=0; i<n; i++) {
		c = fgetc(f);
		if (c == EOF)
			return i;
		s[i] = c;
		s[i+1] = '\0';
	}

	return i;
}

extern size_t
io_util_readline(FILE *f, char *s, size_t n)
{
	int c;
	size_t i;

	s[0] = '\0';
	for (i=0; i<n; i++) {
		c = fgetc(f);
		if ((c == '\n') || (c == EOF))
			return i;
		s[i] = c;
		s[i+1] = '\0';
	}
	return i;
}

extern char *
io_util_strdup(const char *str)
{
	char *dummy;

	dummy = (char *)malloc((strlen(str)+1)*sizeof(char));
	if (dummy == NULL)
		return NULL;

	return strcpy(dummy, str);
}

extern char *
io_util_split_pathfname(const char *str, char **path, char **stem)
{
	int64_t len, i;

	if ( (str == NULL) || (path == NULL) || (stem == NULL) )
		return NULL;

	len = strlen(str);

	/* Figure out where to cut */
	if (str[len-1] == '/') {
		*path = io_util_strdup(str);
		*stem = io_util_strdup("");
		return *path;
	}
	i = 1;
	while ( (i <= len-1) && (str[len-1-i] != '/'))
		i++;
	if ( (i == len) && (str[0] != '/') ) {
		*path = io_util_strdup("./");
	} else {
		*path = (char *)malloc(sizeof(char) * (len-i+1));
		if (*path == NULL)
			return NULL;
		strncpy(*path, str, len-i);
		(*path)[len-i] = '\0';
	}

	*stem = (char *)malloc(sizeof(char) * (i+1));
	if ( *stem == NULL) {
		free(*path);
		*path = NULL;
		return NULL;
	}
	strncpy(*stem, str+(len-i), i);
	(*stem)[i] = '\0';

	return *path;
}

extern int32_t
io_util_findfiles(const char *path, const char *stem, char ***fnames)
{
	int32_t numfiles, i;
	int32_t len;
	char *buff;
	FILE *f;

	if ( (path == NULL) || (stem == NULL) || (fnames == NULL) )
		return INT32_C(0);

	len = strlen(path);
	len += strlen(stem);

	buff = (char *)malloc(sizeof(char)*(len+10));
	if (buff == NULL)	
		return INT32_C(0);

	i = 0;
	f = NULL;

	do  {
		if (f != NULL)
			fclose(f);
		sprintf(buff, "%s%s%d", path, stem, i);
		f=fopen(buff, "r");
		i++;
	} while ( (i<999) && (f != NULL) );
	i--;

	*fnames = (char **)malloc(sizeof(char *) * i);
	if (*fnames == NULL) {
		free(buff);
		return INT32_C(0);
	}
	
	for (numfiles=0; numfiles<i; numfiles++) {
		sprintf(buff, "%s%s%d", path, stem, numfiles);
		(*fnames)[numfiles] = io_util_strdup(buff);
		if ( (*fnames)[numfiles] == NULL ) {
			while (numfiles > 0) {
				numfiles--;
				free((*fnames)[numfiles]);
			}
			free(*fnames);
			free (buff);
			return INT32_C(0);
		}
	}

	/* Clean up */
	free(buff);
	
	return numfiles;
}

extern bool *
io_util_getminearr(int32_t numfiles)
{
	int32_t i = 0;
	bool *dummy = NULL;
	int rank = 0;
	int size = 1;

	/* First get the memory */
	dummy = (bool *)malloc(sizeof(bool)*numfiles);
	if (dummy == NULL)
		return NULL;

#ifdef WITH_MPI
	/* Get the MPI dimensions */
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	/* Now figure out the things to read and not to read */
	for (i=0; i<numfiles; i++) {
		if (rank == i%size)
			dummy[i] = true;
		else
			dummy[i] = false;
	}

	return dummy;
}


/**********************************************************************\
 *    Implementation of local functions                               * 
\**********************************************************************/
