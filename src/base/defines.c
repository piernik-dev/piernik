/*
   Here we check for consistency of the piernik.def file.

   This is a .c file to make sure that #error directive stops the make process.
*/

/*
   Flux limiters

   Exclusive: VANLEER, MONCEN, MINMOD, SUPERBEE
   Default: VANLEER

   Used by src/fluids/fluxes.F90
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
#error Too many flux limiters have been defined. Only one is  allowed.
#endif

#ifndef FLIM
#error: No flux limiters have been defined
#endif
