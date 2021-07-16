! This program was found to fail on Fedora 32 with MPICH 3.3.2 when the modern mpi_f08 interface is in use.
! The old approach with "use mpi" (after necessary changes such as adding ierror argument and
! removing MPI_Allgatherv from the import statement) works correctly.

program Allgatherv_bug

   use mpi_f08, only: MPI_Init, MPI_Finalize, MPI_Comm_size, MPI_Allgatherv, &
        &             MPI_COMM_WORLD, MPI_IN_PLACE, MPI_DATATYPE_NULL, MPI_INTEGER, MPI_INTEGER_KIND

   implicit none

   integer(kind=MPI_INTEGER_KIND) :: nproc, i, o
   integer(kind=MPI_INTEGER_KIND), allocatable, dimension(:) :: se, cnt, off
   integer :: l, u
   integer(kind=MPI_INTEGER_KIND), parameter :: lu_bug = 7

   call MPI_Init()
   call MPI_Comm_size(MPI_COMM_WORLD, nproc)

   allocate(cnt(0:nproc-1), off(0:nproc-1))
   cnt(:) = 1
   o = 0
   do i = lbound(cnt, 1, kind=MPI_INTEGER_KIND), ubound(cnt, 1, kind=MPI_INTEGER_KIND)
      off(i) = o
      o = o + cnt(i)
   end do
   allocate(se(o))
   se(:) = -1

   l = lbound(se, 1)
   u = ubound(se, 1)
   call MPI_Allgatherv(MPI_IN_PLACE, 0_MPI_INTEGER_KIND, MPI_DATATYPE_NULL, se, cnt, off, MPI_INTEGER, MPI_COMM_WORLD)

   if (l /= lbound(se, 1) .or. u /= ubound(se, 1)) call MPI_Abort(MPI_COMM_WORLD, lu_bug)

   call MPI_Finalize()

end program Allgatherv_bug
