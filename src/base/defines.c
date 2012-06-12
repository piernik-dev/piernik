/* $Id$ */
/*
   This file is not meant to be included anywhere. We just check if all declarations in the piernik.def file are consistent.

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

/*
 * Hydro solvers
 *
 * Exclusive: RTVD, MH
 * Default: RTVD
 */

#undef HYDRO_SOLVER

#ifdef RTVD
#define HYDRO_SOLVER
#endif

#ifdef HLLC
#define HYDRO_SOLVER
#endif

#if defined(RTVD) && defined(HLLC)
#error Choose only one of { RTVD, HLLC }.
#endif

#if !defined(HYDRO_SOLVER)
#error no hydro solver defined, possible choices { RTVD, HLLC }
#endif

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

#if defined HLLC && defined CORIOLIS
#error CORIOLIS for HLLC scheme has not been implemented yet.
#endif

#if defined HLLC && defined SN_SRC
#warning SN_SRC step routines are not called in HLLC scheme - calling via problem_customize_solution is required.
#endif

#ifdef USER_RULES
#include "user_rules.h"
#endif
