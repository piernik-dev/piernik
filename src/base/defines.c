/* $Id$ */
/*
   Here we check for consistency of the piernik.def file.

   This is a .c file to make sure that #error directive stops the make process.
*/

#include "piernik.def"

/*
  Freezing speed

  Exclusive: LOCAL_FR_SPEED, GLOBAL_FR_SPEED

  Used by: src/fluids/ionized/fluxionized.F90, src/fluids/neutral/fluxneutral.F90, src/fluids/dust/fluxdust.F90
*/

#undef FR_SPEED
#undef FR_SPEED2

#ifdef LOCAL_FR_SPEED
#define FR_SPEED
#endif

#ifdef GLOBAL_FR_SPEED
#ifdef FR_SPEED
#define FR_SPEED2
#else
#define FR_SPEED
#endif
#endif

#ifdef FR_SPEED2
#error Both freezing speeds defined
#endif

#ifndef FR_SPEED
#error No freezing speed defined.
#endif

#if defined(PSM) || defined(PLN) || defined(KSG) || defined(KSM) || defined(PGM) || defined(SSY) || defined(SI) || defined(CGS) || defined(WT4)
#error Use run-time parameter constants_set from CONSTANTS namelist instead of { PGM SSY SI CGS WT4 PSM PLN KSG KSM } preprocessor symbols.
#endif

/* basic sanity check for isothermal fluid */

#ifdef ISO_LOCAL
#  ifndef ISO
#     error ISO must be defined with ISO_LOCAL
#  endif
#  ifndef IONIZED
#     error ISO_LOCAL currently works only with ionized fluid
#  endif
#endif

/* at least one of { ionized, neutral, dust } must be defined */

#undef FLUID

#ifdef IONIZED
#define FLUID
#endif

#ifdef DUST
#define FLUID
#endif

#ifdef NEUTRAL
#define FLUID
#endif

#ifndef FLUID
#error None of { IONIZED DUST NEUTRAL } were defined.
#endif

#ifdef GRAV
/*
  Constant gravity field

  Exclusive: GRAV_NULL, GRAV_UNIFORM, GRAV_LINEAR, GRAV_PTMASS, GRAV_PTMASSSTIFF, GRAV_PTMASSPURE, GRAV_PTFLAT, GRAV_USER

  Used by: src/gravity/gravity.F90
*/

#ifdef GRAV_NULL
#define GRVC
#endif

#ifdef GRAV_UNIFORM
#ifdef GRVC
#define GRVC2
#else
#define GRVC
#endif
#endif

#ifdef GRAV_LINEAR
#ifdef GRVC
#define GRVC2
#else
#define GRVC
#endif
#endif

#ifdef GRAV_PTMASS
#ifdef GRVC
#define GRVC2
#else
#define GRVC
#endif
#endif

#ifdef GRAV_PTMASSSTIFF
#ifdef GRVC
#define GRVC2
#else
#define GRVC
#endif
#endif

#ifdef GRAV_PTMASSPURE
#ifdef GRVC
#define GRVC2
#else
#define GRVC
#endif
#endif

#ifdef GRAV_PTFLAT
#ifdef GRVC
#define GRVC2
#else
#define GRVC
#endif
#endif

#ifdef GRAV_USER
#ifdef GRVC
#define GRVC2
#else
#define GRVC
#endif
#endif

#ifdef GRVC2
#error Only one of { GRAV_PTMASSSTIFF GRAV_PTMASSPURE GRAV_PTFLAT GRAV_USER GRAV_NULL GRAV_UNIFORM GRAV_LINEAR GRAV_PTMASS } is allowed.
#endif

#ifndef GRVC
#warning None of { GRAV_PTMASSSTIFF GRAV_PTMASSPURE GRAV_PTFLAT GRAV_USER GRAV_NULL GRAV_UNIFORM GRAV_LINEAR GRAV_PTMASS } were defined. Relying on grav_accel.
/*
  ToDo: add GRAV_ACC symbol then change the #warning above to #error
*/
#endif

#endif /* GRAV */

/*
  Gravity solvers

  Exclusive: MULTIGRID, POISSON_FFT
*/

#if defined(MULTIGRID) && defined(POISSON_FFT)
#error MULTIGRID and POISSON_FFT are not meant to work together.
#endif

/*
  Multigrid solver

  at least one of { GRAV, COSM_RAYS }
*/

#ifdef MULTIGRID
#  if !defined(GRAV) && !defined(COSM_RAYS)
#    warning MULTIGRID defined but none of { GRAV, COSM_RAYS } are used.
#  endif
#endif
