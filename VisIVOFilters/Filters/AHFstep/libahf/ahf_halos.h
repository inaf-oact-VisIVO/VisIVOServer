#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../tdef.h"

#if defined AHF || defined TRACKER	 /* function declarations also needed for HaloTracker! */

int	HaloProfiles           (HALO *);
#ifdef AHFphspdens
/** Does the phase-space stuff */
int	HaloProfilesPhaseSpace (HALO *);
#endif
void    rem_nothing            (HALO *);
void    rem_unbound            (HALO *);
void    rem_outsideRvir        (HALO *, int);
void    merge_ll_and_ipart     (HALO *);
void    sort_halo_particles    (HALO *);

#if (defined AHFrestart || defined WITH_MPI)
void
rem_boundary_haloes(void);
#endif

#endif
