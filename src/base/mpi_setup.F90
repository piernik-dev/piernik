 module mpi_setup

! Written by M.Hanasz, April 2006    -  MPI comunication in "z"
! Modified by M. Hanasz - MPI comunication in "z"   - April 2006 
! Modified by M. Hanasz - MPI comunication in "xyz" - November 2006 

!  use mpi

  implicit none
  include 'mpif.h'
  integer nproc, proc, ierr , rc 
  integer status(MPI_STATUS_SIZE,4)
  integer req(4),err(4)

  
  integer, parameter    :: ndims = 3       ! 3D grid
  integer               :: comm, comm3d
  integer, dimension(3) :: psize, pcoords, coords
  logical               :: periods(3), reorder
  integer    procxl, procxr, procyl, procyr, proczl, proczr, procxyl, procyxl 
  integer    pxleft, pxright,pyleft, pyright,pzleft, pzright

  integer,   parameter             :: buffer_dim=200
  character, dimension(buffer_dim) :: cbuff*32
  integer,   dimension(buffer_dim) :: ibuff
  real,      dimension(buffer_dim) :: rbuff
  
  integer pxsize, pysize, pzsize
  namelist /MPI_BLOCKS/ pxsize, pysize, pzsize
  logical mpi


 contains

!-----------------------------------------------------------------------------

  subroutine mpistart

    implicit none
    
   
    character*4 bnd_xl, bnd_xr, bnd_yl, bnd_yr, bnd_zl, bnd_zr

    integer iproc

      pxsize = 1
      pysize = 1
      pzsize = 1
  

    namelist /BOUNDARIES/ bnd_xl, bnd_xr, bnd_yl, bnd_yr, bnd_zl, bnd_zr
      bnd_xl = 'per'  
      bnd_xr = 'per'  
      bnd_yl = 'per'  
      bnd_yr = 'per'  
      bnd_zl = 'per'  
      bnd_zr = 'per'


    open(1,file='problem.par')
      read(unit=1,nml=MPI_BLOCKS)
      read(unit=1,nml=BOUNDARIES)
    close(1)

    psize(1)   = pxsize     
    psize(2)   = pysize    
    psize(3)   = pzsize    

    call MPI_INIT( ierr )
    call MPI_COMM_RANK(MPI_COMM_WORLD, proc, ierr)
    comm = MPI_COMM_WORLD
    call MPI_COMM_SIZE(comm, nproc, ierr)
    
    if(pxsize*pysize*pzsize .ne. nproc) then 
      if(proc .eq.0)  write(*,*) &
      'nproc =',nproc,' MUST BE EQUAL TO   pxsize*pysize*pzsize =',pxsize*pysize*pzsize
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call MPI_Finalize(ierr)                       
      stop
    endif

    if(pxsize*pysize*pzsize .ne. 1) mpi = .true. 
   

    if(bnd_xl(1:3) .eq. 'per') then
      periods(1) = .true.  ! x periodic
    else
      periods(1) = .false.  
    endif
    if(bnd_yl(1:3) .eq. 'per') then
      periods(2) = .true.  ! y periodic
    else
      periods(2) = .false.  
    endif
    if(bnd_zl(1:3) .eq. 'per') then
      periods(3) = .true.  ! z periodic
    else
      periods(3) = .false.  
    endif

    reorder = .false.     ! allows processes reordered for efficiency

    call MPI_CART_CREATE(comm, ndims, psize, periods, reorder, comm3d, ierr)
    call MPI_CART_COORDS(comm3d, proc, ndims, pcoords, ierr)
!    write(*,*) 'proc=',proc, '    coords=', pcoords
    
! Compute neighbors


    pxleft  = pcoords(1) - 1
    pxright = pcoords(1) + 1
    pyleft  = pcoords(2) - 1
    pyright = pcoords(2) + 1
    pzleft  = pcoords(3) - 1
    pzright = pcoords(3) + 1
        
    coords = (/pxleft,pcoords(2),pcoords(3)/)
    if(pxleft .ne. -1 .or. periods(1)) then
      call MPI_Cart_rank(comm3d,coords,procxl,ierr)
    else
      procxl = MPI_PROC_NULL
    endif

    coords = (/pxright,pcoords(2),pcoords(3)/)
    if(pxright .ne. psize(1) .or. periods(1)) then
      call MPI_Cart_rank(comm3d,coords,procxr,ierr)
    else
      procxr = MPI_PROC_NULL
    endif
!    write(*,*) 'xdir: ',procxl, proc, procxr

    coords = (/pcoords(1),pyleft,pcoords(3)/)
    if(pyleft .ne. -1 .or. periods(2)) then
      call MPI_Cart_rank(comm3d,coords,procyl,ierr)
    else
      procyl = MPI_PROC_NULL
    endif

    coords = (/pcoords(1),pyright,pcoords(3)/)
    if(pyright .ne. psize(2) .or. periods(2)) then
      call MPI_Cart_rank(comm3d,coords,procyr,ierr)
    else
      procyr = MPI_PROC_NULL
    endif
!    write(*,*) 'ydir: ',procyl, proc, procyr

    coords = (/pcoords(1),pcoords(2),pzleft/)
    if(pzleft .ne. -1 .or. periods(3)) then
      call MPI_Cart_rank(comm3d,coords,proczl,ierr)
    else
      proczl = MPI_PROC_NULL
    endif

    coords = (/pcoords(1),pcoords(2),pzright/)
    if(pzright .ne. psize(3) .or. periods(3)) then
      call MPI_Cart_rank(comm3d,coords,proczr,ierr)
    else
      proczr = MPI_PROC_NULL
    endif
!    write(*,*) 'zdir: ',proczl, proc, proczr
!    write(*,*)

    if(bnd_xl(1:3) .eq. 'cor' .and. bnd_yl(1:3) .eq. 'cor' ) then
       if(pcoords(1) .eq. 0 .and. pcoords(2) .gt. 0) then
         coords = (/pcoords(2),pcoords(1),pcoords(3)/)
         call MPI_Cart_rank(comm3d,coords,procxyl,ierr)
       else 
         procxyl = -1
       endif
       if(pcoords(2) .eq. 0 .and. pcoords(1) .gt. 0 ) then
         coords = (/pcoords(2),pcoords(1),pcoords(3)/)
         call MPI_Cart_rank(comm3d,coords,procyxl,ierr) 
       else
         procyxl = -1      
       endif

    endif
    
  end subroutine mpistart

  
!-----------------------------------------------------------------------------


  subroutine mpistop

   implicit none

   call MPI_FINALIZE(rc)

  end subroutine mpistop

!-----------------------------------------------------------------------------


  subroutine mpifind(var, what, loc_arr, loc_proc)
  
    implicit none
    
    
    character what*(*)
    real var, rsend(2), rrecv(2)
    integer loc_arr(3), loc_proc
        
    rsend(1) = var
    rsend(2) = proc
    
    select case (trim(what))
      case('min') 
        CALL MPI_REDUCE(rsend, rrecv, 1, MPI_2DOUBLE_PRECISION, &
                                         MPI_MINLOC, 0, comm, ierr)
      case('max')
        CALL MPI_REDUCE(rsend, rrecv, 1, MPI_2DOUBLE_PRECISION, &
                                         MPI_MAXLOC, 0, comm, ierr)
      case default
        Write(*,*) 'mpifind: actual parameter "', what, '"is not allowed' 
    end select

    if(proc .eq. 0) then
      var = rrecv(1)
      loc_proc = rrecv(2)
    endif
        
    call MPI_BCAST(loc_proc, 1, MPI_INTEGER, 0, comm, ierr)
   
    if(loc_proc .ne. 0) then
      if(proc .eq. loc_proc) then
        CALL MPI_SEND  (loc_arr, 3, MPI_INTEGER,     0, 11, comm, ierr)
      else if(proc .eq. 0) then
        CALL MPI_RECV  (loc_arr, 3, MPI_INTEGER, loc_proc, 11, comm, status, ierr)
      endif
    endif

  end subroutine mpifind

!------------------------------------------------------------------------------------------

end module mpi_setup


