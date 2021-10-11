! This program may help to determine if MPI library works correctly.
!
! Known failure cases indicated by timeout message:
! * In Ubuntu 18.04 with OpenMPI 2.1.x failures may occur randomly for some
!   message counts and sizes indicating message loss. It easily occurs for
!   4096 messages with 512 double precsion FP each.
! * Message count 64k and above may result in extremely slow completion of
!   requests.
!
! OpenMPI 4.0.3 in Ubuntu 20.04 seems to work correctly and is much faster
! for big message count, when compared to its 2.1 release.
!
! This program scans a range of parameters and is too expensive to be called
! from Makefile. With some modifications it is also possible to determine
! communications speed differences between consecutive MPI ranks and detect
! alow links.


program fat_ring

   use mpi, only: MPI_Init, MPI_Finalize, MPI_Comm_size, MPI_Comm_rank, MPI_Wtime, MPI_Barrier, &
        &         MPI_COMM_WORLD, MPI_INTEGER_KIND

   implicit none

   integer(kind=MPI_INTEGER_KIND) :: nproc, proc, ierr
   integer :: num, sz, i, totalloc
   integer, parameter :: maxnum = 14, maxsize = 10, mix = 10000
   integer, parameter :: memlimit = 30  ! 30 means to not allocate more than 2**30 bytes (1 GiB) per process
   type :: arr_T
      real(kind=8), allocatable, dimension(:,:,:) :: s
      real(kind=8), allocatable, dimension(:,:,:) :: r
   end type arr_T
   type(arr_T), allocatable, dimension(:) :: set
   double precision :: t0, t1

   call MPI_Init(ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, proc, ierr)

   do num = max(0, min(8, maxnum - 6)), maxnum
      allocate(set(2**num))
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      if (proc == 0) write(*,*)""
      do sz = 0, maxsize
         call MPI_Barrier(MPI_COMM_WORLD, ierr)
         if (num + 3*sz + 1 + 3 > memlimit) exit
         totalloc = 0
         do i = lbound(set, 1), ubound(set, 1)
            allocate(set(i)%s(2**sz, 2**sz, 2**sz), &
                 &   set(i)%r(2**sz, 2**sz, 2**sz))
            totalloc = totalloc + size(set(i)%s) + size(set(i)%r)
            set(i)%s = proc + mix * i
            set(i)%r = -1
         enddo
         t0 = MPI_Wtime()
         call sendrecv
         t1 = MPI_Wtime()
         call MPI_Barrier(MPI_COMM_WORLD, ierr)
         if (proc == 0) write(*,'(a,i5,a,i8,a,i8,a,i10,a,i5,2(a,f7.3))') &
              "  @ ", proc, " #msgs", size(set)*2, " Size(msg): ", 2**(3*sz), &
              " Allocated/thread= ", totalloc, " doubles (", totalloc / 2**17, &
              " MiB) time(master)= ", t1 - t0, " time(all)= ", MPI_Wtime() - t0
         do i = lbound(set, 1), ubound(set, 1)
            deallocate(set(i)%s, set(i)%r)
         enddo
      enddo
      deallocate(set)
   enddo

   call MPI_Finalize(ierr)

contains

   subroutine sendrecv

      use mpi, only: MPI_DOUBLE_PRECISION, MPI_STATUSES_IGNORE, MPI_STATUS_IGNORE

      implicit none

      integer(kind=MPI_INTEGER_KIND), allocatable, dimension(:) :: req
      double precision, parameter :: timeout = 60.
      double precision :: wt
      integer(kind=MPI_INTEGER_KIND) :: tcnt, cnt_prev
      integer(kind=MPI_INTEGER_KIND), allocatable, dimension(:) :: inds
      logical, parameter :: crash_easily = .false.
      logical :: flag

      allocate(req(2*size(set)))
      do i = lbound(set, 1), ubound(set, 1)
         call MPI_Isend(set(i)%s, size(set(i)%s, kind=MPI_INTEGER_KIND), MPI_DOUBLE_PRECISION, &
              &         mod(        proc + 1, nproc), i, MPI_COMM_WORLD, req(1 + 2*(i - lbound(set,1))), ierr)
         call MPI_Irecv(set(i)%r, size(set(i)%r, kind=MPI_INTEGER_KIND), MPI_DOUBLE_PRECISION, &
              &         mod(nproc + proc - 1, nproc), i, MPI_COMM_WORLD, req(2 + 2*(i - lbound(set,1))), ierr)
      enddo

      wt = MPI_Wtime()
      allocate(inds(size(req)))
      inds(:) = -1
      tcnt = 0
      do while (tcnt < size(req))

         cnt_prev = tcnt
         tcnt = 0
         call MPI_Testsome(size(req), req, tcnt, inds(cnt_prev+1:), MPI_STATUSES_IGNORE, ierr)

         if (tcnt < 0 .or. cnt_prev < 0. .or. ierr /= 0) then
            write(*,*)"    @", proc, " MPI_Testsome failing? ", tcnt, cnt_prev, ierr
            if (crash_easily) call MPI_Abort(MPI_COMM_WORLD, 11, ierr)
         endif

         tcnt = cnt_prev + tcnt

         if (tcnt < size(req) .and. MPI_Wtime() - wt > timeout) then
            write(*, '(a,i5,2(a,i7),a,i9,a,f7.3,a)')"    -@", proc, " : only ", tcnt, " out of ", size(req), &
                 " requests of size ", size(set(lbound(set, 1))%s), " doubles, completed in ", MPI_Wtime() - wt, "s (timeout)"
            if (crash_easily) then
               call MPI_Abort(MPI_COMM_WORLD, 13, ierr)
            else
               do i = 1, size(req)
                  call MPI_Test(req(i), flag, MPI_STATUS_IGNORE, ierr)
                  if (.not. flag) then
                     call MPI_Cancel(req(i), ierr)
                     call MPI_Request_free(req(i), ierr)
                  endif
               enddo
            endif
            exit
         endif

      enddo
      call MPI_Waitall(size(req, kind=MPI_INTEGER_KIND), req, MPI_STATUSES_IGNORE, ierr)
      deallocate(req, inds)

   end subroutine sendrecv

end program fat_ring
