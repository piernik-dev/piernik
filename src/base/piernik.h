
#include "piernik.def"

#if defined(MULTIGRID) && defined(GRAV)
#define SELF_GRAV
#endif

#if defined(VARIABLE_USER_GP) || defined(SELF_GRAV)
#define VARIABLE_GP
#endif

#define HDF5
#if defined(I_KNOW_WHAT_I_AM_DOING)
#undef HDF5
#endif

#if !defined(RTVD) && !defined(HLLC)
#define RTVD
/* #  warning no hydro solver defined, possible choices { RTVD, HLLC }, defaulting to RTVD */
#endif

/*
  Disable ISO EOS for benchmaring with Rieman solver because we don't support ISO there yet.

  ToDo: Remove this hack as soon as we add ISO to the Riemann solver.
*/
#ifdef BENCHMARKING_HACK
#  undef ISO
#endif
