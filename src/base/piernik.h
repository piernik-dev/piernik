
#include "piernik.def"

#if defined(MULTIGRID) && defined(GRAV)
#define SELF_GRAV
#endif

#if defined(VARIABLE_USER_GP) || defined(SELF_GRAV)
#define VARIABLE_GP
#endif /* VARIABLE_USER_GP || SELF_GRAV */

#define HDF5
#if defined(I_KNOW_WHAT_I_AM_DOING)
#undef HDF5
#endif /* I_KNOW_WHAT_I_AM_DOING */

#if !defined(RTVD) && !defined(HLLC) && !defined(RIEMANN)
#define RTVD
/* #  warning no hydro solver defined, possible choices { RTVD, HLLC, RIEMANN }, defaulting to RTVD */
#endif /* !RTVD && !HLLC && !RIEMANN */

/*
  Disable ISO EOS for benchmaring with Rieman solver because we don't support ISO there yet.

  ToDo: Remove this hack as soon as we add ISO to the Riemann solver.
*/
#ifdef BENCHMARKING_HACK
#  undef ISO
#endif
