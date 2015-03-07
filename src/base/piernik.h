
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

#if !defined(RTVD) && !defined(HLLC) && !defined(RIEMANN)
#define RTVD
/* #  warning no hydro solver defined, possible choices { RTVD, HLLC, RIEMANN }, defaulting to RTVD */
#endif
