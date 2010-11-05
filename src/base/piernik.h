/* $Id$ */
#include "piernik.def"

#if defined(MULTIGRID) || defined(POISSON_FFT)
#define SELF_GRAV
#endif

#ifdef VARIABLE_USER_GP
#define VARIABLE_GP
#endif

#if defined (GRAV_PTMASSPURE) || defined (GRAV_PTMASS) || defined (GRAV_PTFLAT) || defined (GRAV_PTMASSSTIFF)
#define GRAV_PTMTYPE
#endif