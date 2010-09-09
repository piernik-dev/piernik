/*
   Here we check for consistency of the piernik.def file.

   This is a .c file to make sure that #error directive stops the make process.
*/

/*
   Flux limiters

   Exclusive: VANLEER, MONCEN, MINMOD, SUPERBEE

   Used by: src/fluids/fluxes.F90
*/

#include "piernik.def"

#undef FLIM
#undef FLIM2

#ifdef VANLEER
#define FLIM
#endif

#ifdef MONCEN
#ifdef FLIM
#define FLIM2
#else
#define FLIM
#endif
#endif

#ifdef MINMOD
#ifdef FLIM
#define FLIM2
#else
#define FLIM
#endif
#endif

#ifdef SUPERBEE
#ifdef FLIM
#define FLIM2
#else
#define FLIM
#endif
#endif

#ifdef FLIM2
#error Too many flux limiters have been defined. Only one is allowed.
#endif

#ifndef FLIM
#error No flux limiters have been defined.
#endif

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

/*
  Units

  Exclusive: PGM, SSY, SI, CGS, WT4, PSM, PLN, KSG, KSM

  Code generated automatically by:
  ./setup -units | grep "\!" | awk '{print $2}' |\
  awk 'BEGIN {\
    def="UNIT";\
    n=0;\
  } {\
    print "#ifdef", $1;\
    if (n>0) print "#ifdef",def,"\n#define",def"2\n#else";\
    print "#define",def;\
    if (n>0) print "#endif";\
    print "#endif\n";\
    a[n++]=$1;\
  } END {\
    printf "%s %s", "#ifdef",def"2\n#error Only one of { ";\
    for (b in a) printf "%s ", a[b];\
    printf "%s %s %s","} is allowed.\n#endif\n\n#ifndef",def,"\n#error None of { ";\
    for (b in a) printf "%s ", a[b];\
    printf "%s\n","} were defined.\n#endif";\
  }'
*/

#ifdef PSM
#define UNIT
#endif

#ifdef PLN
#ifdef UNIT
#define UNIT2
#else
#define UNIT
#endif
#endif

#ifdef KSG
#ifdef UNIT
#define UNIT2
#else
#define UNIT
#endif
#endif

#ifdef KSM
#ifdef UNIT
#define UNIT2
#else
#define UNIT
#endif
#endif

#ifdef PGM
#ifdef UNIT
#define UNIT2
#else
#define UNIT
#endif
#endif

#ifdef SSY
#ifdef UNIT
#define UNIT2
#else
#define UNIT
#endif
#endif

#ifdef SI
#ifdef UNIT
#define UNIT2
#else
#define UNIT
#endif
#endif

#ifdef CGS
#ifdef UNIT
#define UNIT2
#else
#define UNIT
#endif
#endif

#ifdef WT4
#ifdef UNIT
#define UNIT2
#else
#define UNIT
#endif
#endif

#ifdef UNIT2
#error Only one of { PGM SSY SI CGS WT4 PSM PLN KSG KSM } is allowed
#endif

#ifndef UNIT
#warning None of { PGM SSY SI CGS WT4 PSM PLN KSG KSM } were defined.
#warning Assuming SCALED
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
  ToDo:
    GRAV & co ?
*/
