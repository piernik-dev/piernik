
/* spelling workaround */

#ifdef IONISED
#  define IONIZED
#endif /* IONISED */

#include "piernik.def"

#ifdef IONIZED
#  ifndef NONMAGNETIC
#    ifndef MAGNETIC
#      define MAGNETIC
#    endif /* !MAGNETIC */
#  endif /* !NONMAGNETIC */
#endif /* IONIZED */

#if !defined(MAGNETIC) && defined(RESISTIVE)
#define NOMAGNETICNORESIST
#undef RESISTIVE
#endif /* !MAGNETIC && RESISTIVE */

#ifdef SELF_GRAV
#  ifndef GRAV
#    define GRAV
#  endif /* !GRAV */
#  ifndef MULTIGRID
#    define MULTIGRID
#  endif /* !MULTIGRID */
#endif /* SELF_GRAV */

#if defined(VARIABLE_USER_GP) || defined(SELF_GRAV)
#define VARIABLE_GP
#endif /* VARIABLE_USER_GP || SELF_GRAV */

#define HDF5
#if defined(I_KNOW_WHAT_I_AM_DOING)
#undef HDF5
#endif /* I_KNOW_WHAT_I_AM_DOING */

