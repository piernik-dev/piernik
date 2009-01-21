! $Id$
#include "piernik.def"

module mag_boundaries

! Written  by M. Hanasz - January/February 2006
! Modified by M. Hanasz - MPI comunication in "z"   - April 2006
! Modified by M. Hanasz - MPI comunication in "xyz" - November 2006
! Modified by M. Hanasz - MPI corner-periodic bcs   - December 2006
! Modified by M. Hanasz - MPI shearing-periodic in "x" - November 2007
  use mpi_setup
  use fluidindex, only : ibx, iby, ibz  
  use start, only  : dimensions
  use grid, only : nb,nx,ny,nz,nxb,nyb,nzb
  use arrays, only : u,b

contains

subroutine bnd_a(A)
   implicit none
   integer :: ireq
   real, dimension(:,:,:,:) :: A

   if(pxsize .gt. 1) then

      CALL MPI_ISEND  (A(1,1,1,1), 1, MAG_YZ_LEFT_DOM,  procxl, 10, comm3d, req(1), ierr)
      CALL MPI_ISEND  (A(1,1,1,1), 1, MAG_YZ_RIGHT_DOM, procxr, 20, comm3d, req(3), ierr)
      CALL MPI_IRECV  (A(1,1,1,1), 1, MAG_YZ_LEFT_BND,  procxl, 20, comm3d, req(2), ierr)
      CALL MPI_IRECV  (A(1,1,1,1), 1, MAG_YZ_RIGHT_BND, procxr, 10, comm3d, req(4), ierr)

      do ireq=1,4
         call MPI_WAIT(req(ireq),status(1,ireq),ierr)
      enddo
   endif

   if(pysize .gt. 1) then
      CALL MPI_ISEND  (A(1,1,1,1), 1, MAG_XZ_LEFT_DOM,  procyl, 30, comm3d, req(1), ierr)
      CALL MPI_ISEND  (A(1,1,1,1), 1, MAG_XZ_RIGHT_DOM, procyr, 40, comm3d, req(3), ierr)
      CALL MPI_IRECV  (A(1,1,1,1), 1, MAG_XZ_LEFT_BND,  procyl, 40, comm3d, req(2), ierr)
      CALL MPI_IRECV  (A(1,1,1,1), 1, MAG_XZ_RIGHT_BND, procyr, 30, comm3d, req(4), ierr)

      do ireq=1,4
         call MPI_WAIT(req(ireq),status(1,ireq),ierr)
      enddo
   endif

   if(pzsize .gt. 1) then
      CALL MPI_ISEND  (A(1,1,1,1), 1, MAG_XY_LEFT_DOM,  proczl, 50, comm3d, req(1), ierr)
      CALL MPI_ISEND  (A(1,1,1,1), 1, MAG_XY_RIGHT_DOM, proczr, 60, comm3d, req(3), ierr)
      CALL MPI_IRECV  (A(1,1,1,1), 1, MAG_XY_LEFT_BND,  proczl, 60, comm3d, req(2), ierr)
      CALL MPI_IRECV  (A(1,1,1,1), 1, MAG_XY_RIGHT_BND, proczr, 50, comm3d, req(4), ierr)

      do ireq=1,4
         call MPI_WAIT(req(ireq),status(1,ireq),ierr)
      enddo
   endif
end subroutine bnd_a

subroutine bnd_b(dim)
#ifdef SHEAR
  use shear, only : eps,delj
#endif /* SHEAR */

  implicit none
  character(len=*) :: dim
  integer i, j
  integer ireq
  real, allocatable :: send_left(:,:,:,:),recv_left(:,:,:,:)
#ifdef SHEAR
  real, allocatable :: send_right(:,:,:,:),recv_right(:,:,:,:)
#endif /* SHEAR */
! MPI block comunication

  select case (dim)
    case ('xdim')
#ifdef SHEAR
      allocate(send_left(3,nb,ny,nz),send_right(3,nb,ny,nz), &
               recv_left(3,nb,ny,nz),recv_right(3,nb,ny,nz))

      send_left (:,:,:,:)  = b(:,nb+1:2*nb,:,:)
      send_right(:,:,:,:)  = b(:,nxb+1:nxb+nb,:,:)

      if(bnd_xl == 'she') then
!
! przesuwamy o calkowita liczbe komorek + periodyczny wb w kierunku y
!
        send_left (:,:,nb+1:nb+nyb,:)         = cshift(send_left (:,:,nb+1:nb+nyb,:),dim=3,shift= delj)
        send_left (:,:,1:nb,:)                = send_left  (:,:,nyb+1:nyb+nb,:)
        send_left (:,:,nb+nyb+1:nyb+2*nb,:)   = send_left  (:,:,nb+1:2*nb,:)
!
! remapujemy  - interpolacja kwadratowa
!
        send_left (:,:,:,:)  = (1.+eps)*(1.-eps) * send_left (:,:,:,:) &
                               -0.5*eps*(1.-eps) * cshift(send_left (:,:,:,:),shift=-1,dim=3) &
                               +0.5*eps*(1.+eps) * cshift(send_left (:,:,:,:),shift=1,dim=3)
      endif ! (bnd_xl == 'she')


      if(bnd_xr == 'she') then
!
! przesuwamy o calkowita liczbe komorek + periodyczny wb w kierunku y
!
        send_right (:,:,nb+1:nb+nyb,:)        = cshift(send_right(:,:,nb+1:nb+nyb,:),dim=3,shift=-delj)
        send_right (:,:,1:nb,:)               = send_right (:,:,nyb+1:nyb+nb,:)
        send_right (:,:,nb+nyb+1:nyb+2*nb,:)  = send_right (:,:,nb+1:2*nb,:)
!
! remapujemy - interpolacja kwadratowa
!
        send_right (:,:,:,:) = (1.+eps)*(1.-eps) * send_right (:,:,:,:) &
                               -0.5*eps*(1.-eps) * cshift(send_right (:,:,:,:),shift=1,dim=3) &
                               +0.5*eps*(1.+eps) * cshift(send_right (:,:,:,:),shift=-1,dim=3)
      endif ! (bnd_xr == 'she')
!
! wysylamy na drugi brzeg
!
        CALL MPI_ISEND   (send_left , 3*ny*nz*nb, MPI_DOUBLE_PRECISION, procxl, 10, comm, req(1), ierr)
        CALL MPI_ISEND   (send_right, 3*ny*nz*nb, MPI_DOUBLE_PRECISION, procxr, 20, comm, req(3), ierr)
        CALL MPI_IRECV   (recv_left , 3*ny*nz*nb, MPI_DOUBLE_PRECISION, procxl, 20, comm, req(2), ierr)
        CALL MPI_IRECV   (recv_right, 3*ny*nz*nb, MPI_DOUBLE_PRECISION, procxr, 10, comm, req(4), ierr)

        do ireq=1,4
           call MPI_WAIT(req(ireq),status(1,ireq),ierr)
        enddo

        b(:,1:nb-1,:,:)               = recv_left(:,1:nb-1,:,:)
        b(:,nxb+nb+1+1:nxb+2*nb,:,:)  = recv_right(:,1+1:nb,:,:)

      deallocate(send_left,send_right,recv_left,recv_right)
!===============================================================================
#else /* SHEAR */

      if(pxsize .gt. 1) then

!        if(procxl .ne. MPI_PROC_NULL) send_left(:,:,:,:)          =  b(:,nb+1:2*nb,:,:)
!        if(procxr .ne. MPI_PROC_NULL) send_right(:,:,:,:)         =  b(:,nxb+1:nxb+nb,:,:)

!        CALL MPI_ISEND   (send_left , 3*ny*nz*nb, MPI_DOUBLE_PRECISION, procxl, 10, comm, req(1), ierr)
!        CALL MPI_ISEND   (send_right, 3*ny*nz*nb, MPI_DOUBLE_PRECISION, procxr, 20, comm, req(3), ierr)
!        CALL MPI_IRECV   (recv_left , 3*ny*nz*nb, MPI_DOUBLE_PRECISION, procxl, 20, comm, req(2), ierr)
!        CALL MPI_IRECV   (recv_right, 3*ny*nz*nb, MPI_DOUBLE_PRECISION, procxr, 10, comm, req(4), ierr)
        CALL MPI_ISEND  (b(1,1,1,1), 1, MAG_YZ_LEFT_DOM,  procxl, 10, comm3d, req(1), ierr)
        CALL MPI_ISEND  (b(1,1,1,1), 1, MAG_YZ_RIGHT_DOM, procxr, 20, comm3d, req(3), ierr)
        CALL MPI_IRECV  (b(1,1,1,1), 1, MAG_YZ_LEFT_BND,  procxl, 20, comm3d, req(2), ierr)
        CALL MPI_IRECV  (b(1,1,1,1), 1, MAG_YZ_RIGHT_BND, procxr, 10, comm3d, req(4), ierr)

        do ireq=1,4
          call MPI_WAIT(req(ireq),status(1,ireq),ierr)
        enddo

!        if(procxl .ne. MPI_PROC_NULL) b(:,1:nb,:,:)               = recv_left(:,:,:,:)
!        if(procxr .ne. MPI_PROC_NULL) b(:,nxb+nb+1:nxb+2*nb,:,:)  = recv_right(:,:,:,:)

      endif
#endif /* SHEAR */


    case ('ydim')
      if(pysize .gt. 1) then
!        allocate(send_left(3,nx,nb,nz),send_right(3,nx,nb,nz), &
!                 recv_left(3,nx,nb,nz),recv_right(3,nx,nb,nz))

!        if(procyl .ne. MPI_PROC_NULL) send_left(:,:,:,:)          =  b(:,:,nb+1:2*nb,:)
!        if(procyr .ne. MPI_PROC_NULL) send_right(:,:,:,:)         =  b(:,:,nyb+1:nyb+nb,:)

!        CALL MPI_ISEND   (send_left , 3*nx*nb*nz, MPI_DOUBLE_PRECISION, procyl, 30, comm, req(1), ierr)
!        CALL MPI_ISEND   (send_right, 3*nx*nb*nz, MPI_DOUBLE_PRECISION, procyr, 40, comm, req(3), ierr)
!        CALL MPI_IRECV   (recv_left , 3*nx*nb*nz, MPI_DOUBLE_PRECISION, procyl, 40, comm, req(2), ierr)
!        CALL MPI_IRECV   (recv_right, 3*nx*nb*nz, MPI_DOUBLE_PRECISION, procyr, 30, comm, req(4), ierr)
        CALL MPI_ISEND  (b(1,1,1,1), 1, MAG_XZ_LEFT_DOM,  procyl, 30, comm3d, req(1), ierr)
        CALL MPI_ISEND  (b(1,1,1,1), 1, MAG_XZ_RIGHT_DOM, procyr, 40, comm3d, req(3), ierr)
        CALL MPI_IRECV  (b(1,1,1,1), 1, MAG_XZ_LEFT_BND,  procyl, 40, comm3d, req(2), ierr)
        CALL MPI_IRECV  (b(1,1,1,1), 1, MAG_XZ_RIGHT_BND, procyr, 30, comm3d, req(4), ierr)

        do ireq=1,4
          call MPI_WAIT(req(ireq),status(1,ireq),ierr)
        enddo

!        if(procyl .ne. MPI_PROC_NULL) b(:,:,1:nb,:)               = recv_left(:,:,:,:)
!        if(procyr .ne. MPI_PROC_NULL) b(:,:,nyb+nb+1:nyb+2*nb,:)  = recv_right(:,:,:,:)

!        deallocate(send_left,send_right,recv_left,recv_right)
      endif

    case ('zdim')
      if(pzsize .gt. 1) then
!        allocate(send_left(3,nx,ny,nb),send_right(3,nx,ny,nb), &
!                 recv_left(3,nx,ny,nb),recv_right(3,nx,ny,nb))

!        if(proczl .ne. MPI_PROC_NULL) send_left(:,:,:,:)          =  b(:,:,:,nb+1:2*nb)
!        if(proczr .ne. MPI_PROC_NULL) send_right(:,:,:,:)         =  b(:,:,:,nzb+1:nzb+nb)

!        CALL MPI_ISEND   (send_left , 3*nx*ny*nb, MPI_DOUBLE_PRECISION, proczl, 50, comm, req(1), ierr)
!        CALL MPI_ISEND   (send_right, 3*nx*ny*nb, MPI_DOUBLE_PRECISION, proczr, 60, comm, req(3), ierr)
!        CALL MPI_IRECV   (recv_left , 3*nx*ny*nb, MPI_DOUBLE_PRECISION, proczl, 60, comm, req(2), ierr)
!        CALL MPI_IRECV   (recv_right, 3*nx*ny*nb, MPI_DOUBLE_PRECISION, proczr, 50, comm, req(4), ierr)
        CALL MPI_ISEND  (b(1,1,1,1), 1, MAG_XY_LEFT_DOM,  proczl, 50, comm3d, req(1), ierr)
        CALL MPI_ISEND  (b(1,1,1,1), 1, MAG_XY_RIGHT_DOM, proczr, 60, comm3d, req(3), ierr)
        CALL MPI_IRECV  (b(1,1,1,1), 1, MAG_XY_LEFT_BND,  proczl, 60, comm3d, req(2), ierr)
        CALL MPI_IRECV  (b(1,1,1,1), 1, MAG_XY_RIGHT_BND, proczr, 50, comm3d, req(4), ierr)


        do ireq=1,4
          call MPI_WAIT(req(ireq),status(1,ireq),ierr)
        enddo

!        if(proczl .ne. MPI_PROC_NULL) b(:,:,:,1:nb)               = recv_left(:,:,:,:)
!        if(proczr .ne. MPI_PROC_NULL) b(:,:,:,nzb+nb+1:nzb+2*nb)  = recv_right(:,:,:,:)

!        deallocate(send_left,send_right,recv_left,recv_right)
      endif
  end select ! (dim)

! MPI + non-MPI corner-periodic boundary condition

  if(bnd_xl .eq. 'cor') then
!   - lower to left
    if(pcoords(1) .eq. 0 .and. pcoords(2) .eq. 0) then
      do i=1,nb
        do j=nb+1,ny
          b(ibx,i,j,:) = -b(iby,j,2*nb+1-i,:)
          b(iby,i,j,:) =  b(ibx,j,2*nb+1-i,:)
          b(ibz,i,j,:) =  b(ibz,j,2*nb+1-i,:)
        enddo
      enddo
    endif

    if(procxyl .gt. 0) then
      allocate(send_left(3,nb,ny,nz), recv_left(3,nx,nb,nz))

      send_left(:,:,:,:) = b(:,nb+1:2*nb,:,:)

      CALL MPI_ISEND   (send_left , 3*nb*ny*nz, MPI_DOUBLE_PRECISION, procxyl, 70, comm, req(1), ierr)
      CALL MPI_IRECV   (recv_left , 3*nx*nb*nz, MPI_DOUBLE_PRECISION, procxyl, 80, comm, req(2), ierr)

      do ireq=1,2
        call MPI_WAIT(req(ireq),status(1,ireq),ierr)
      enddo

      do i=1,nb
        do j=1,ny
          b(ibx,i,j,:) = -recv_left(iby,j,nb+1-i,:)
          b(iby,i,j,:) =  recv_left(ibx,j,nb+1-i,:)
          b(ibz,i,j,:) =  recv_left(ibz,j,nb+1-i,:)
        enddo
      enddo

      deallocate(send_left,recv_left)
    endif
  endif

  if(bnd_yl .eq. 'cor') then
!   - left to lower
    if(pcoords(2) .eq. 0 .and. pcoords(1) .eq. 0 ) then
      do j=1,nb
        do i=nb+1,nx
          b(ibx,i,j,:) =  b(iby,2*nb+1-j,i,:)
          b(iby,i,j,:) = -b(ibx,2*nb+1-j,i,:)
          b(ibz,i,j,:) =  b(ibz,2*nb+1-j,i,:)
        enddo
      enddo
!   - interior to corner
      do j=1,nb
        do i=1,nb
          b(ibx,i,j,:) =  -b(ibx,2*nb+1-i,2*nb+1-j,:)
          b(iby,i,j,:) =  -b(iby,2*nb+1-i,2*nb+1-j,:)
          b(ibz,i,j,:) =   b(ibz,2*nb+1-i,2*nb+1-j,:)
        enddo
      enddo
    endif

    if(procyxl .gt. 0) then
      allocate(send_left(3,nx,nb,nz), recv_left(3,nb,ny,nz))

      send_left(:,:,:,:) = b(:,:,nb+1:2*nb,:)

      CALL MPI_ISEND   (send_left , 3*nx*nb*nz, MPI_DOUBLE_PRECISION, procyxl, 80, comm, req(1), ierr)
      CALL MPI_IRECV   (recv_left , 3*nb*ny*nz, MPI_DOUBLE_PRECISION, procyxl, 70, comm, req(2), ierr)

      do ireq=1,2
        call MPI_WAIT(req(ireq),status(1,ireq),ierr)
      enddo

      do j=1,nb
        do i=1,nx
          b(ibx,i,j,:) =  recv_left(iby,nb+1-j,i,:)
          b(iby,i,j,:) = -recv_left(ibx,nb+1-j,i,:)
          b(ibz,i,j,:) =  recv_left(ibz,nb+1-j,i,:)
        enddo
      enddo

      deallocate(send_left,recv_left)
    endif
  endif


! Non-MPI boundary conditions

  select case (dim)
    case ('xdim')

      select case (bnd_xl(1:3))
        case ('she')
!         Do nothing if 'she'
        case ('mpi')
!         Do nothing if 'mpi'
        case ('per')
          b(:,1:nb,:,:)                  = b(:,nxb+1:nxb+nb,:,:)
        case ('cor')
!         Do nothing if 'cor'
        case ('inf')
!         Do nothing if 'inf'
        case ('ref')
!         Do nothing if 'ref'
        case ('out')
          b(:,1,:,:) = b(:,2,:,:)
!         Do nothing if 'out'
        case default
          write(*,*) 'Boundary condition ',bnd_xl,' not implemented in ',dim
      end select  ! (bnd_xl)

      select case (bnd_xr(1:3))
        case ('she')
!         Do nothing if 'she'
        case ('mpi')
!         Do nothing if 'mpi'
        case ('per')
          b(:,nxb+nb+1:nxb+2*nb,:,:)     = b(:,nb+1:2*nb,:,:)
        case ('cor')
!         Do nothing if 'cor'
        case ('inf')
!         Do nothing if 'inf'
        case ('ref')
!         Do nothing if 'ref'
        case ('out')
!         Do nothing if 'out'
          b(:,nx,:,:) = b(:,nx-1,:,:)
        case default
          write(*,*) 'Boundary condition ',bnd_xr,' not implemented in ',dim
      end select  ! (bnd_xr)


    case ('ydim')

      select case (bnd_yl(1:3))
        case ('mpi')
!         Do nothing if 'mpi'
        case ('per')
          b(:,:,1:nb,:)                  = b(:,:,nyb+1:nyb+nb,:)
        case ('cor')
!         Do nothing if 'cor'
        case ('inf')
!         Do nothing if 'inf'
        case ('ref')
!         Do nothing if 'ref'
        case ('out')
          b(:,:,1,:) = b(:,:,2,:)
!         Do nothing if 'out'
        case default
          write(*,*) 'Boundary condition ',bnd_yl,' not implemented in ',dim
      end select  ! (bnd_yl)

      select case (bnd_yr(1:3))
        case ('mpi')
!         Do nothing if 'mpi'
        case ('per')
          b(:,:,nyb+nb+1:nyb+2*nb,:)     = b(:,:,nb+1:2*nb,:)
        case ('cor')
!         Do nothing if 'cor'
        case ('inf')
!         Do nothing if 'inf'
        case ('ref')
!         Do nothing if 'ref'
        case ('out')
!         Do nothing if 'out'
          b(:,:,ny,:) = b(:,:,ny-1,:)
        case default
          write(*,*) 'Boundary condition ',bnd_yr,' not implemented in ',dim

      end select  ! (bnd_yr)


    case ('zdim')

      select case (bnd_zl(1:3))
        case ('mpi')
!         Do nothing if 'mpi'
        case ('per')
          b(:,:,:,1:nb)                  = b(:,:,:,nzb+1:nzb+nb)
        case ('ref')
!         Do nothing if 'ref'
        case ('out')
          b(:,:,:,1) = b(:,:,:,2)
!         Do nothing if 'out'
        case default
          write(*,*) 'Boundary condition ',bnd_zl,' not implemented in ',dim
      end select  ! (bnd_zl)

      select case (bnd_zr(1:3))
        case ('mpi')
!         Do nothing if 'mpi'
        case ('per')
          b(:,:,:,nzb+nb+1:nzb+2*nb)     = b(:,:,:,nb+1:2*nb)
        case ('ref')
!         Do nothing if 'ref'
        case ('out')
!         Do nothing if 'out'
          b(:,:,:,nz) = b(:,:,:,nz-1)
        case default
          write(*,*) 'Boundary condition ',bnd_zr,' not implemented in ',dim
      end select  ! (bnd_zr)

    end select  ! (dim)

end subroutine bnd_b


!=====================================================================================================

subroutine bnd_emf(var, name, dim)
  use grid, only : nx,ny,nz

  implicit none
  real, dimension(nx,ny,nz) :: var
  real dvarx(ny,nz),  dvary(nx,nz), dvarz(nx,ny)
  character name*4, dim*4
  integer ib


  select case (dim)
    case ('xdim')

      select case (name)
        case ('vxby','vxbz')

          select case (bnd_xl(1:3))
            case ('she')
!             Do nothing if 'she'
            case ('mpi')
!             Do nothing if 'mpi'
            case ('per')
!             Do nothing if 'per'
            case ('cor')
!             Do nothing if 'cor'
            case ('inf')
!             Do nothing if 'inf'
            case ('ref')
              var(nb,:,:)                 = 0.0
              do ib=1,nb-1
                var(nb-ib,:,:)            = -var(nb+ib,:,:)
              enddo
            case ('out')
              dvarx = var(nb+1,:,:)-var(nb,:,:)
              do ib=1,nb-1
                var(ib,:,:)               = var(nb+1,:,:)  - real(nb+1-ib)*dvarx
              enddo
            case default
              write(*,*) 'Boundary condition ', bnd_xl, ' not implemented for ',name,' in ',dim
          end select  ! (bnd_xl)

          select case (bnd_xr(1:3))
            case ('she')
!             Do nothing if 'she'
            case ('mpi')
!             Do nothing if 'mpi'
            case ('per')
!             Do nothing if 'per'
            case ('cor')
!             Do nothing if 'cor'
            case ('inf')
!             Do nothing if 'inf'
            case ('ref')
              var(nb+nxb,:,:)             = 0.0
              do ib=1,nb-1
                var(nb+nxb+ib,:,:)        = -var(nb+nxb-ib,:,:)
              enddo
            case ('out')
              dvarx = var(nb+nxb,:,:)-var(nb+nxb-1,:,:)
              do ib=1,nb-1
                var(nb+nxb+ib,:,:)        = var(nb+nxb,:,:) + real(ib)*dvarx
              enddo
            case default
              write(*,*) 'Boundary condition ', bnd_xr, ' not implemented for ',name,' in ',dim
          end select ! (bnd_xr)

        case ('vybx','vzbx','emfy','emfz')

          select case (bnd_xl(1:3))
            case ('she')
!             Do nothing if 'she'
            case ('mpi')
!             Do nothing if 'mpi'
            case ('per')
!             Do nothing if 'per'
            case ('cor')
!             Do nothing if 'cor'
            case ('inf')
!             Do nothing if 'inf'
            case ('ref')
              var(nb+1,:,:)               = 0.0
              do ib=1,nb
                var(nb+1-ib,:,:)          = -var(nb+1+ib,:,:)
              enddo
            case ('out')
              dvarx = var(nb+2,:,:)-var(nb+1,:,:)
              do ib=1,nb
                var(ib,:,:)               = var(nb+1,:,:)  - real(nb+1-ib)*dvarx
              enddo
            case default
              write(*,*) 'Boundary condition ', bnd_xl, ' not implemented for ',name,' in ',dim
          end select   ! (bnd_xl)

          select case (bnd_xr(1:3))
            case ('she')
!             Do nothing if 'she'
            case ('mpi')
!             Do nothing if 'mpi'
            case ('per')
!             Do nothing if 'per'
            case ('cor')
!             Do nothing if 'cor'
            case ('inf')
!             Do nothing if 'inf'
            case ('ref')
              var(nb+nxb+1,:,:)           = 0.0
              do ib=1,nb-1
                var(nb+nxb+1+ib,:,:)      = -var(nb+nxb+1-ib,:,:)
              enddo
            case ('out')
              dvarx = var(nb+nxb+1,:,:)-var(nb+nxb,:,:)
              do ib=1,nb-1
                var(nb+nxb+1+ib,:,:)      =  var(nb+nxb+1,:,:) + real(ib)*dvarx
              enddo
            case default
              write(*,*) 'Boundary condition ', bnd_xr, ' not implemented for ',name,' in ',dim
          end select  ! (bnd_xr)

        case ('vybz','vzby','emfx')

          select case (bnd_xl(1:3))
            case ('she')
!             Do nothing if 'she'
            case ('mpi')
!             Do nothing if 'mpi'
            case ('per')
!             Do nothing if 'per'
            case ('cor')
!             Do nothing if 'cor'
            case ('inf')
!             Do nothing if 'inf'
            case ('ref')
              do ib=1,nb
                var(nb+1-ib,:,:)          = var(nb+ib,:,:)
              enddo
            case ('out')
              dvarx = var(nb+1,:,:)-var(nb,:,:)
              do ib=1,nb
                var(ib,:,:)               = var(nb+1,:,:)  - real(nb+1-ib)*dvarx
              enddo
            case default
              write(*,*) 'Boundary condition ', bnd_xl, ' not implemented for ',name,' in ',dim
          end select   ! (bnd_xl)

          select case (bnd_xr(1:3))
            case ('she')
!             Do nothing if 'she'
            case ('mpi')
!             Do nothing if 'mpi'
            case ('per')
!             Do nothing if 'per'
            case ('cor')
!             Do nothing if 'cor'
            case ('inf')
!             Do nothing if 'inf'
            case ('ref')
              do ib=1,nb
                var(nb+nxb+ib,:,:)        = var(nb+nxb+1-ib,:,:)
              enddo
            case ('out')
              dvarx = var(nb+nxb+1,:,:)-var(nb+nxb,:,:)
              do ib=1,nb
                var(nb+nxb+ib,:,:)        = var(nb+nxb,:,:) + real(ib)*dvarx
              enddo
            case default
              write(*,*) 'Boundary condition ', bnd_xr, ' not implemented for ',name,' in ',dim
          end select   ! (bnd_xr)

       end select  ! (name)

    case ('ydim')

      select case (name)
        case ('vybz','vybx')

          select case (bnd_yl(1:3))
            case ('mpi')
!             Do nothing if 'mpi'
            case ('per')
!             Do nothing if 'per'
            case ('cor')
!             Do nothing if 'cor'
            case ('inf')
!             Do nothing if 'inf'
            case ('ref')
              var(:,nb,:)                 = 0.0
              do ib=1,nb-1
                var(:,nb-ib,:)            =-var(:,nb+ib,:)
              enddo
            case ('out')
              dvary = var(:,nb+1,:)-var(:,nb,:)
              do ib=1,nb-1
                var(:,ib,:)               = var(:,nb+1,:) - real(nb+1-ib)*dvary
              enddo
            case default
              write(*,*) 'Boundary condition ', bnd_yl, ' not implemented for ',name,' in ',dim
          end select  ! (bnd_yl)

          select case (bnd_yr(1:3))
            case ('mpi')
!             Do nothing if 'mpi'
            case ('per')
!             Do nothing if 'per'
            case ('cor')
!             Do nothing if 'cor'
            case ('inf')
!             Do nothing if 'inf'
            case ('ref')
              var(:,nb+nyb,:)             = 0.0
              do ib=1,nb-1
                var(:,nb+nyb+ib,:)        =-var(:,nb+nyb-ib,:)
              enddo
            case ('out')
              dvary = var(:,nb+nyb,:)-var(:,nb+nyb-1,:)
              do ib=1,nb-1
                var(:,nb+nyb+ib,:)        = var(:,nb+nyb,:) + real(ib)*dvary
              enddo
            case default
              write(*,*) 'Boundary condition ', bnd_yr, ' not implemented for ',name,' in ',dim
          end select ! (bnd_yr)

        case ('vzby','vxby','emfz','emfx')

          select case (bnd_yl(1:3))
            case ('mpi')
!             Do nothing if 'mpi'
            case ('per')
!             Do nothing if 'per'
            case ('cor')
!             Do nothing if 'cor'
            case ('inf')
!             Do nothing if 'inf'
            case ('ref')
              var(:,nb+1,:)               = 0.0
              do ib=1,nb
                var(:,nb+1-ib,:)          = -var(:,nb+1+ib,:)
              enddo
            case ('out')
              dvary = var(:,nb+2,:)-var(:,nb+1,:)
              do ib=1,nb
                var(:,ib,:)               = var(:,nb+1,:) - real(nb+1-ib)*dvary
              enddo
            case default
              write(*,*) 'Boundary condition ', bnd_yl, ' not implemented for ',name,' in ',dim
          end select   ! (bnd_yl)

          select case (bnd_yr(1:3))
            case ('mpi')
!             Do nothing if 'mpi'
            case ('per')
!             Do nothing if 'per'
            case ('cor')
!             Do nothing if 'cor'
            case ('inf')
!             Do nothing if 'inf'
            case ('ref')
              var(:,nb+nyb+1,:)           = 0.0
              do ib=1,nb-1
                var(:,nb+nyb+1+ib,:)      = -var(:,nb+nyb+1-ib,:)
              enddo
            case ('out')
              dvary = var(:,nb+nyb+1,:)-var(:,nb+nyb,:)
              do ib=1,nb-1
                var(:,nb+nyb+1+ib,:)      =  var(:,nb+nyb+1,:) + real(ib)*dvary
              enddo
            case default
              write(*,*) 'Boundary condition ', bnd_yr, ' not implemented for ',name,' in ',dim
          end select  ! (bnd_yr)

        case ('vzbx','vxbz','emfy')

          select case (bnd_yl(1:3))
            case ('mpi')
!             Do nothing if 'mpi'
            case ('per')
!             Do nothing if 'per'
            case ('cor')
!             Do nothing if 'cor'
            case ('inf')
!             Do nothing if 'inf'
            case ('ref')
              do ib=1,nb
                var(:,nb+1-ib,:)          = var(:,nb+ib,:)
              enddo
            case ('out')
              dvary = var(:,nb+1,:)-var(:,nb,:)
              do ib=1,nb
                var(:,ib,:)               = var(:,nb+1,:) - real(nb+1-ib)*dvary
              enddo
            case default
              write(*,*) 'Boundary condition ', bnd_yl, ' not implemented for ',name,' in ',dim
          end select   ! (bnd_yl)

          select case (bnd_yr(1:3))
            case ('mpi')
!             Do nothing if 'mpi'
            case ('per')
!             Do nothing if 'per'
            case ('cor')
!             Do nothing if 'cor'
            case ('inf')
!             Do nothing if 'inf'
            case ('ref')
              do ib=1,nb
                var(:,nb+nyb+ib,:)        = var(:,nb+nyb+1-ib,:)
              enddo
            case ('out')
              dvary = var(:,nb+nyb+1,:)-var(:,nb+nyb,:)
              do ib=1,nb
                var(:,nb+nyb+ib,:)        = var(:,nb+nyb,:) + real(ib)*dvary
              enddo
            case default
              write(*,*) 'Boundary condition ', bnd_yr, ' not implemented for ',name,' in ',dim
          end select   ! (bnd_yr)

       end select  ! (name)


    case ('zdim')

      select case (name)

        case ('vzbx','vzby')

          select case (bnd_zl(1:3))
            case ('mpi')
!             Do nothing if 'mpi'
            case ('per')
!             Do nothing if 'per'
            case ('cor')
!             Do nothing if 'cor'
            case ('inf')
!             Do nothing if 'inf'
            case ('ref')
              var(:,:,nb)                 = 0.0
              do ib=1,nb-1
                var(:,:,nb-ib)            =-var(:,:,nb+ib)
              enddo
            case ('out')
              dvarz = var(:,:,nb+1)-var(:,:,nb)
              do ib=1,nb-1
                var(:,:,ib)               = var(:,:,nb+1) - real(nb+1-ib)*dvarz
              enddo
            case default
              write(*,*) 'Boundary condition ', bnd_zl, ' not implemented for ',name,' in ',dim
          end select  ! (bnd_zl)

         select case (bnd_zr(1:3))
           case ('mpi')
!             Do nothing if 'mpi'
           case ('per')
!             Do nothing if 'per'
            case ('cor')
!             Do nothing if 'cor'
            case ('inf')
!             Do nothing if 'inf'
           case ('ref')
             var(:,:,nb+nzb)             = 0.0
             do ib=1,nb-1
               var(:,:,nb+nzb+ib)        =-var(:,:,nb+nzb-ib)
             enddo
           case ('out')
             dvarz = var(:,:,nb+nzb)-var(:,:,nb+nzb-1)
             do ib=1,nb-1
               var(:,:,nb+nzb+ib)        = var(:,:,nb+nzb) + real(ib)*dvarz
             enddo
            case default
              write(*,*) 'Boundary condition ', bnd_zr, ' not implemented for ',name,' in ',dim
          end select ! (bnd_zr)


        case ('vxbz','vybz','emfy','emfx')

          select case (bnd_zl(1:3))
            case ('mpi')
!             Do nothing if 'mpi'
            case ('per')
!             Do nothing if 'per'
            case ('cor')
!             Do nothing if 'cor'
            case ('inf')
!             Do nothing if 'inf'
            case ('ref')
              var(:,:,nb+1)               = 0.0
              do ib=1,nb
                var(:,:,nb+1-ib)          = -var(:,:,nb+1+ib)
              enddo
            case ('out')
              dvarz = var(:,:,nb+2)-var(:,:,nb+1)
              do ib=1,nb
                var(:,:,ib)               = var(:,:,nb+1) - real(nb+1-ib)*dvarz
              enddo
            case default
              write(*,*) 'Boundary condition ', bnd_zl, ' not implemented for ',name,' in ',dim
          end select   ! (bnd_zl)

          select case (bnd_zr(1:3))
           case ('mpi')
!             Do nothing if 'mpi'
           case ('per')
!             Do nothing if 'per'
            case ('cor')
!             Do nothing if 'cor'
            case ('inf')
!             Do nothing if 'inf'
           case ('ref')
             var(:,:,nb+nzb+1)           = 0.0
             do ib=1,nb-1
               var(:,:,nb+nzb+1+ib)      = -var(:,:,nb+nzb+1-ib)
             enddo
           case ('out')
             dvarz = var(:,:,nb+nzb+1)-var(:,:,nb+nzb)
             do ib=1,nb-1
               var(:,:,nb+nzb+1+ib)      =  var(:,:,nb+nzb+1) + real(ib)*dvarz
             enddo
           case default
             write(*,*) 'Boundary condition ', bnd_zr, ' not implemented for ',name,' in ',dim
          end select   ! (bnd_zr)

        case ('vxby','vybx','emfz')

          select case (bnd_zl(1:3))
            case ('mpi')
!             Do nothing if 'mpi'
            case ('per')
!             Do nothing if 'per'
            case ('cor')
!             Do nothing if 'cor'
            case ('inf')
!             Do nothing if 'inf'
            case ('ref')
              do ib=1,nb
                var(:,:,nb+1-ib)          = var(:,:,nb+ib)
              enddo
            case ('out')
              dvarz = var(:,:,nb+1)-var(:,:,nb)
              do ib=1,nb
                var(:,:,ib)               = var(:,:,nb+1) - real(nb+1-ib)*dvarz
              enddo
            case default
              write(*,*) 'Boundary condition ', bnd_zl, ' not implemented for ',name,' in ',dim
          end select   ! (bnd_zl)

          select case (bnd_zr(1:3))
            case ('mpi')
!             Do nothing if 'mpi'
            case ('per')
!             Do nothing if 'per'
            case ('cor')
!             Do nothing if 'cor'
            case ('inf')
!             Do nothing if 'inf'
            case ('ref')
              do ib=1,nb
                var(:,:,nb+nzb+ib)        = var(:,:,nb+nzb+1-ib)
              enddo
            case ('out')
             dvarz = var(:,:,nb+nzb+1)-var(:,:,nb+nzb)
              do ib=1,nb
                var(:,:,nb+nzb+ib)        = var(:,:,nb+nzb) + real(ib)*dvarz
              enddo
            case default
              write(*,*) 'Boundary condition ', bnd_zr, ' not implemented for ',name,' in ',dim
          end select ! (bnd_zr)

        end select  ! (name)

   end select ! (dim)
end subroutine bnd_emf

  subroutine compute_b_bnd

   call bnd_b('xdim')
   call bnd_b('ydim')
   if(dimensions .eq. '3d') call bnd_b('zdim')

  end subroutine compute_b_bnd

end module mag_boundaries
