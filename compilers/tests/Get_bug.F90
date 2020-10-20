program onesided

   use mpi_f08, only: MPI_COMM_WORLD, MPI_INFO_NULL, MPI_ADDRESS_KIND, MPI_ORDER_FORTRAN, &
        &             MPI_Win, MPI_Datatype, MPI_INTEGER, &
        &             MPI_Init, MPI_Finalize, MPI_Comm_rank, MPI_Comm_size, MPI_Abort, &
        &             MPI_Win_free, MPI_Win_fence, MPI_Get

   implicit none

   integer, parameter :: asize = 10, saoff = 4, sasize = 2, ssoff = 6
   integer, allocatable, dimension(:,:) :: arr, srr
   integer(kind=4) :: err_mpi
   integer(kind=4) :: nproc          !< number of processes
   integer(kind=4) :: proc           !< rank of my process
   type(MPI_Win) :: wina
   type(MPI_Datatype) :: subst
   integer, dimension(2) :: sz, ssz, off
   integer, dimension(4) :: bounds

   call MPI_Init(err_mpi)

   call MPI_Comm_rank(MPI_COMM_WORLD, proc, err_mpi)
   call MPI_Comm_size(MPI_COMM_WORLD, nproc, err_mpi)

   allocate(arr(proc*asize:(proc+1)*asize-1,asize), &
        &   srr(proc*asize:(proc+1)*asize-1,asize))
   arr = 1
   srr = -1

   sz(:) = asize
   ssz(:) = sasize
   off(:) = ssoff
   call MPI_Type_create_subarray(2, sz, ssz, off, MPI_ORDER_FORTRAN, MPI_INTEGER, subst, err_mpi)
   call MPI_Type_commit(subst, err_mpi)

   call MPI_Win_create(arr, size(arr), 8, MPI_INFO_NULL, MPI_COMM_WORLD, wina, err_mpi)
   call MPI_Win_fence(0, wina, err_mpi)
   bounds = [ lbound(srr), ubound(srr) ]
   call MPI_Get(srr, 1, subst, mod(proc + 3, nproc), int(0, kind=MPI_ADDRESS_KIND), 1, subst, wina, err_mpi)
   if (any(bounds  /= [ lbound(srr), ubound(srr) ])) then
      write(*,*)"bounds changed: [", bounds, "] /= [", [ lbound(srr), ubound(srr) ], "]"
      call MPI_Abort(MPI_COMM_WORLD, -2, err_mpi)
   endif
   call MPI_Win_fence(0, wina, err_mpi)
   call MPI_Win_free(wina, err_mpi)

   call MPI_Type_free(subst, err_mpi)

   deallocate(arr, srr)

   call MPI_Finalize(err_mpi)

end program onesided
