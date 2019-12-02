/*
   This file is not meant to be included anywhere. We just check if all declarations in the piernik.def file are consistent.

   This is a .c file to make sure that #error directive stops the make process.
*/

#include "piernik.h"

/*
  Freezing speed

  Exclusive: LOCAL_FR_SPEED, GLOBAL_FR_SPEED

  Used by: src/fluids/ionized/fluxionized.F90, src/fluids/neutral/fluxneutral.F90, src/fluids/dust/fluxdust.F90
*/

#undef FR_SPEED
#undef FR_SPEED2

#ifdef LOCAL_FR_SPEED
#  define FR_SPEED
#endif /* LOCAL_FR_SPEED */

#ifdef GLOBAL_FR_SPEED
#  ifdef FR_SPEED
#    define FR_SPEED2
#  else /* !FR_SPEED */
#    define FR_SPEED
#  endif /* !FR_SPEED */
#endif /* GLOBAL_FR_SPEED */

#ifdef FR_SPEED2
#  error Both freezing speeds defined
#endif /* FR_SPEED2 */

#ifndef FR_SPEED
#  error No freezing speed defined.
#endif /* !FR_SPEED */

#if defined(PSM) || defined(PLN) || defined(KSG) || defined(KSM) || defined(PGM) || defined(SSY) || defined(SI) || defined(CGS) || defined(WT4)
#  error Use run-time parameter constants_set from CONSTANTS namelist instead of { PGM SSY SI CGS WT4 PSM PLN KSG KSM } preprocessor symbols.
#endif /* PSM || PLN || KSG || KSM || PGM || SSY || SI || CGS || WT4 */

/* basic sanity check for isothermal fluid */

#ifdef ISO_LOCAL
#  ifndef ISO
#     error ISO must be defined with ISO_LOCAL
#  endif /* !ISO */
#  ifndef IONIZED
#     error ISO_LOCAL currently works only with ionized fluid
#  endif /* !IONIZED */
#endif /* ISO_LOCAL */

/* at least one of { ionized, neutral, dust } must be defined */

#undef FLUID

#ifdef IONIZED
#  define FLUID
#endif /* IONIZED */

#ifdef DUST
#  define FLUID
#endif /* DUST */

#ifdef NEUTRAL
#  define FLUID
#  ifdef IONIZED
#    error Currently there are no solvers that can manage a mixture of neutral and ionized fluid
#  endif /* IONIZED */
#endif /* NEUTRAL */

#ifndef FLUID
#  error None of { IONIZED DUST NEUTRAL } were defined.
#endif /* !FLUID */

#ifdef NOMAGNETICNORESIST
# warning MAGNETIC is not defined, then RESISTIVE is also cancelled
#endif /* NOMAGNETICNORESIST */

/*
  Multigrid solver

  at least one of { GRAV, COSM_RAYS }
*/

#ifdef MULTIGRID
#  if !defined(SELF_GRAV) && !defined(COSM_RAYS)
#    warning MULTIGRID defined but none of { SELF_GRAV, COSM_RAYS } are used.
#  endif /* !SELF_GRAV && !COSMIC_RAYS */
#endif /* MULTIGRID */

#ifdef USER_RULES
#  include "user_rules.h"
#endif /* USER_RULES */
