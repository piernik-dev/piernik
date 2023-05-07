! Test whether the older MPI Fortran interface is sufficiently modern.
! Older implementations may require some special care, especially since gfortran 10.x

program mpi_new

#ifndef MPIF08

#  ifndef NO_ALL_MPI_FUNCTIONS_AVAILABLE
   use mpi, only: MPI_Send
#  else
#    warning Current MPI library does not provide modern Fortran interface
#  endif

#endif

end program mpi_new
