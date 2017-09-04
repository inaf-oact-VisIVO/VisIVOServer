#if (defined NEWSTARTRUN)
#ifndef AHF_HALOS_SFC
#define AHF_HALOS_SFC

#include "../tdef.h"

/**
 * \brief  This is the function that will built the haloes from the
 *         given centres and gather radii.
 *
 * \param  halo  The halo to built.
 *
 * \return  Returns nothing, however the external halo structure will be
 *          updated to reflect the changes done in this routine.
 */
void
ahf_halos_sfc_constructHalo(HALO *halo);

/**
 * \brief  This will gather all particles for a given halo which are
 *         inside its  gathering radius.
 *
 * \param  halo  The halo to work on.
 *
 * \return  Returns nothing, however the external halo structure will be
 *          updated to reflect the changes done in this routine.
 */
void
ahf_halos_sfc_gatherParts(HALO *halo);


#endif  /* AHF_HALOS_SFC */
#endif /* defined NEWSTARTRUN */
