! Test whether the older MPI Fortran interface is sufficiently modern.
! Older implementations are no longer supported since gfortran 10.x
program mpi_new
   use mpi, only: MPI_Send
end program mpi_new
