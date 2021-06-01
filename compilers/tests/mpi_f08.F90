! Test for avalability of modern MPI interface (Fortran 2008 with type checking)

program F08

#ifndef NO_MPIF08_AVAILABLE

   use mpi_f08

#  ifdef FORBID_F08
#    error The use of mpi_f08 was excluded by FORBID_F08
#  endif

#else

#  ifdef FORBID_F08
#    warning The use of mpi_f08 was excluded by FORBID_F08
#  endif

#  ifdef MPIF08
#    error Both NO_MPIF08_AVAILABLE and MPIF08 are not allowed
#  endif

#endif

end program F08
