/* $Id$ */
#include "piernik.def"

#if ( defined(MULTIGRID) || defined(POISSON_FFT) ) && defined(GRAV)
#define SELF_GRAV
#endif

#if defined(VARIABLE_USER_GP) || defined(SELF_GRAV)
#define VARIABLE_GP
#endif

#if defined (GRAV_PTMASSPURE) || defined (GRAV_PTMASS) || defined (GRAV_PTFLAT) || defined (GRAV_PTMASSSTIFF)
#define GRAV_PTMTYPE
#endif
