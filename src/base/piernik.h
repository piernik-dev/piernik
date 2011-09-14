#define _PIERNIK_H_VERSION_STR '$Id'

/*

  Define preprocessor symbols that depend on other preprocessor symbols.

 */

#include "piernik.def"

#if ( defined(MULTIGRID) || defined(POISSON_FFT) ) && defined(GRAV)
#define SELF_GRAV
#endif

#if defined(VARIABLE_USER_GP) || defined(SELF_GRAV)
#define VARIABLE_GP
#endif
