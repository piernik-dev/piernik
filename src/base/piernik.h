
/* spelling workaround */

#ifdef IONISED
#  define IONIZED
#endif /* IONISED */

#include "piernik.def"

#define RNG cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke

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

#ifdef NBODY
#   define NBODY_MULTIGRID
#   define NBODY_1FILE
#endif /* NBODY */

#ifdef NBODY_MULTIGRID
#  ifndef SELF_GRAV
#    define SELF_GRAV
#  endif /* !SELF_GRAV */
#endif /* NBODY_MULTIGRID */

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

#ifdef MPIF08
#  define MPIF mpi_f08
#  define MPIFUN mpi_f08
#else /* !MPIF08 */
#  define MPIF mpi
#  ifdef NO_ALL_MPI_FUNCTIONS_AVAILABLE
/* Ignore the rest of list to avoid import errors for MPICH and the old Fortran interface.
   One may also import something harmless like MPI_OP_NULL at a cost of some more warnings. */
#    define MPIFUN mpi !
#  else /* !NO_ALL_MPI_FUNCTIONS_AVAILABLE */
#    define MPIFUN mpi
#  endif /* !NO_ALL_MPI_FUNCTIONS_AVAILABLE */
#endif /* !MPIF08 */
