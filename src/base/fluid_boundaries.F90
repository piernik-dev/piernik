! $Id$
#include "piernik.def"

module fluid_boundaries

! Written  by M. Hanasz - December 2005 - February 2006
! Modified by M. Hanasz - MPI comunication in "z"      - April 2006
! Modified by M. Hanasz - MPI comunication in "xyz"    - November 2006
! Modified by M. Hanasz - MPI corner-periodic bcs      - December 2006
! Modified by M. Hanasz - MPI shearing-periodic in "x" - November 2007


contains

subroutine bnd_u(dim)
  use mpi_setup
  use start,  only : smalld, smallei, bnd_xl, bnd_xr, bnd_yl, bnd_yr, &
      bnd_zl, bnd_zr, nb, nxd, nyd, nzd, dimensions
#ifdef COSM_RAYS
    use start, only : smallecr
#endif /* COSM_RAYS */



  use arrays, only : x,z,nzb,nyb,nxb,nu,nx,ny,nz, idna, imxa, imya, imza, &
      u, b

#ifndef ISO
  use start,  only : gamma
  use arrays, only : iena
#endif /* ISO */
#ifdef ISO
  use start, only : csi2
#endif /* ISO */
#ifdef GRAV
  use gravity, only : grav_accel
  use start, only : nsub, tune_zeq_bnd
#endif /* GRAV */
#ifdef SHEAR
  use shear, only : eps,delj, unshear_fft_b, unshear_fft, unshear
  use start, only : qshear, omega
#endif /* SHEAR */
#ifndef SPLIT
  use arrays, only : Lu
#endif /* SPLIT */
#ifdef COSM_RAYS
    use arrays, only : iecr
#endif /* COSM_RAYS */


  implicit none
  character(len=*) :: dim
  integer ib

#ifdef GRAV
  integer kb, ksub
  real, dimension(nx,ny) :: db, ekb, eib, csi2b
  real, dimension(nsub+1):: zs, dprofs, gprofs
  real dzs, factor, z1, z2
#endif /* GRAV */
  integer i,j
  integer ireq
  real, allocatable :: send_left(:,:,:,:),recv_left(:,:,:,:)
#ifdef SHEAR_MPI
  real, allocatable :: send_right(:,:,:,:),recv_right(:,:,:,:)
#endif /* SHEAR_MPI */
#ifdef SHEAR_MY
  real, allocatable, dimension(:,:,:) :: temp,tem2
#endif /* SHEAR_MY */
! MPI block comunication

!  write(*,*) '********************************************'
!  write(*,*) '*        You should not enter bnd_u        *'
!  write(*,*) '*Boundary conditions are done on fluxes !!!*'
!  write(*,*) '********************************************'
  select case (dim)
    case ('xdim')
#ifdef SHEAR_MPI
        allocate(send_right(nu,nb,ny,nz), send_left(nu,nb,ny,nz), &
                 recv_left(nu,nb,ny,nz), recv_right(nu,nb,ny,nz) )
        send_left(:,:,:,:)          =  u(:,nb+1:2*nb,:,:)
        send_right(:,:,:,:)         =  u(:,nxb+1:nxb+nb,:,:)
!
! odejmujemy ped_y i energie odpowiadajace niezaburzonej rozniczkowej rotacji na lewym brzegu
!
        if(bnd_xl == 'she') then
          do i=1,nb
            send_left (imya,i,:,:) = send_left(imya,i,:,:) &
                                         +qshear*omega * x(nb+i)     * send_left(idna,i,:,:)
#ifndef ISO
            send_left (iena,i,:,:) = send_left(iena,i,:,:) &
                                    -0.5*(qshear*omega * x(nb+i))**2 * send_left(idna,i,:,:)
!            send_left (iena,i,:,:) = send_left(iena,i,:,:) &
!                                    +(qshear*omega * x(nb+i))*send_left(imya,i,:,:) &
!                                    +0.5*(qshear*omega*x(nb+i))*send_left(idna,i,:,:)
#endif /* ISO */
          enddo
!
! przesuwamy o calkowita liczbe komorek + periodyczny wb w kierunku y
!
          send_left (:,:,nb+1:nb+nyb,:)        = cshift(send_left (:,:,nb+1:nb+nyb,:),dim=3,shift= delj)
          send_left (:,:,1:nb,:)               = send_left (:,:,nyb+1:nyb+nb,:)
          send_left (:,:,nb+nyb+1:nyb+2*nb,:)  = send_left (:,:,nb+1:2*nb,:)
!
! remapujemy  - interpolacja kwadratowa
!
          send_left (:,:,:,:)  = (1.+eps)*(1.-eps) * send_left (:,:,:,:) &
                                 -0.5*eps*(1.-eps) * cshift(send_left(:,:,:,:),shift=-1,dim=3) &
                                 +0.5*eps*(1.+eps) * cshift(send_left(:,:,:,:),shift=1 ,dim=3)
        endif !(bnd_xl == 'she')
!
! odejmujemy ped_y i energie odpowiadajace niezaburzonej rozniczkowej rotacji na prawym brzegu
!
        if(bnd_xr == 'she') then
          do i=1,nb
            send_right(imya,i,:,:) = send_right(imya,i,:,:) &
                                         +qshear*omega * x(nxb+i)     * send_right(idna,i,:,:)
#ifndef ISO
            send_right(iena,i,:,:) = send_right(iena,i,:,:) &
                                    -0.5*(qshear*omega * x(nxb+i))**2 * send_right(idna,i,:,:)
!            send_right(iena,i,:,:) = send_right(iena,i,:,:) &
!                                    +(qshear*omega*x(nxb+i))*send_right(imya,i,:,:) &
!                                    +0.5*(qshear*omega * x(nxb+i))**2 * send_right(idna,i,:,:)

#endif /* ISO */
          enddo
!
! przesuwamy o calkowita liczbe komorek + periodyczny wb w kierunku y
!
          send_right(:,:,nb+1:nb+nyb,:)        = cshift(send_right(:,:,nb+1:nb+nyb,:),dim=3,shift=-delj)
          send_right (:,:,1:nb,:)              = send_right(:,:,nyb+1:nyb+nb,:)
          send_right (:,:,nb+nyb+1:nyb+2*nb,:) = send_right(:,:,nb+1:2*nb,:)
!
! remapujemy  - interpolacja kwadratowa
!
          send_right (:,:,:,:) = (1.+eps)*(1.-eps) * send_right (:,:,:,:) &
                                 -0.5*eps*(1.-eps) * cshift(send_right(:,:,:,:),shift=1 ,dim=3) &
                                 +0.5*eps*(1.+eps) * cshift(send_right(:,:,:,:),shift=-1,dim=3)
        endif !(bnd_xr == 'she')
!
! wysylamy na drugi brzeg
!
        CALL MPI_ISEND   (send_left , nu*ny*nz*nb, MPI_DOUBLE_PRECISION, procxl, 10, comm, req(1), ierr)
        CALL MPI_ISEND   (send_right, nu*ny*nz*nb, MPI_DOUBLE_PRECISION, procxr, 20, comm, req(3), ierr)
        CALL MPI_IRECV   (recv_left , nu*ny*nz*nb, MPI_DOUBLE_PRECISION, procxl, 20, comm, req(2), ierr)
        CALL MPI_IRECV   (recv_right, nu*ny*nz*nb, MPI_DOUBLE_PRECISION, procxr, 10, comm, req(4), ierr)

        do ireq = 1,4
          call MPI_WAIT(req(ireq),status(1,ireq),ierr)
        enddo

!
! dodajemy ped_y i energie odpowiadajace niezaburzonej rozniczkowej rotacji na prawym brzegu
!
        if(bnd_xr == 'she') then
          do i=1,nb
#ifndef ISO
             recv_right (iena,i,:,:) = recv_right (iena,i,:,:) &
                                      +0.5*(qshear*omega * x(nb+nxb+i))**2 * recv_right(idna,i,:,:)
!             recv_right (iena,i,:,:) = recv_right (iena,i,:,:) &
!                                      -(qshear*omega*x(nb+nxb+i))*recv_right(imya,i,:,:) &
!                                      +0.5*(qshear*omega * x(nb+nxb+i))**2 * recv_right(idna,i,:,:)
#endif /* ISO */
             recv_right (imya,i,:,:) = recv_right (imya,i,:,:) &
                                           -qshear*omega * x(nb+nxb+i)     * recv_right(idna,i,:,:)
          enddo

        endif !(bnd_xr == 'she')
!
! dodajemy ped_y i energie odpowiadajace niezaburzonej rozniczkowej rotacji na lewym brzegu
!
        if(bnd_xl == 'she') then

          do i=1,nb
#ifndef ISO
             recv_left(iena,i,:,:) = recv_left(iena,i,:,:) &
                                    +0.5*(qshear*omega * x(i))**2 * recv_left(idna,i,:,:)
!             recv_left(iena,i,:,:) = recv_left(iena,i,:,:) &
!                                    -(qshear*omega*x(i))*recv_left(imya,i,:,:)&
!                                    +0.5*(qshear*omega * x(i))**2 * recv_left(idna,i,:,:)
#endif /* ISO */
             recv_left(imya,i,:,:) = recv_left(imya,i,:,:) &
                                         -qshear*omega * x(i)     * recv_left(idna,i,:,:)
          enddo
        endif !(bnd_xl == 'she')

        u(:,1:nb,:,:)              = recv_left(:,1:nb,:,:)
        u(:,nxb+nb+1:nxb+2*nb,:,:) = recv_right(:,1:nb,:,:)

        u(idna,1:nb,:,:)              = max(u(idna,1:nb,:,:),smalld)
        u(idna,nxb+nb+1:nxb+2*nb,:,:) = max(u(idna,nxb+nb+1:nxb+2*nb,:,:),smalld)
      deallocate(send_left,send_right,recv_left,recv_right)
#else /* SHEAR_MPI */
      if(pxsize .gt. 1) then

!        if(procxl .ne. MPI_PROC_NULL) send_left(:,:,:,:)          =  u(:,nb+1:2*nb,:,:)
!        if(procxr .ne. MPI_PROC_NULL) send_right(:,:,:,:)         =  u(:,nxb+1:nxb+nb,:,:)

!        CALL MPI_ISEND   (send_left , nu*ny*nz*nb, MPI_DOUBLE_PRECISION, procxl, 10, comm, req(1), ierr)
!        CALL MPI_ISEND   (send_right, nu*ny*nz*nb, MPI_DOUBLE_PRECISION, procxr, 20, comm, req(3), ierr)
!        CALL MPI_IRECV   (recv_left , nu*ny*nz*nb, MPI_DOUBLE_PRECISION, procxl, 20, comm, req(2), ierr)
!        CALL MPI_IRECV   (recv_right, nu*ny*nz*nb, MPI_DOUBLE_PRECISION, procxr, 10, comm, req(4), ierr)
        CALL MPI_ISEND   (u(1,1,1,1), 1, MPI_YZ_LEFT_DOM,  procxl, 10, comm3d, req(1), ierr)
        CALL MPI_ISEND   (u(1,1,1,1), 1, MPI_YZ_RIGHT_DOM, procxr, 20, comm3d, req(3), ierr)
        CALL MPI_IRECV   (u(1,1,1,1), 1, MPI_YZ_LEFT_BND,  procxl, 20, comm3d, req(2), ierr)
        CALL MPI_IRECV   (u(1,1,1,1), 1, MPI_YZ_RIGHT_BND, procxr, 10, comm3d, req(4), ierr)

        do ireq=1,4
          call MPI_WAIT(req(ireq),status(1,ireq),ierr)
        enddo

!        if(procxl .ne. MPI_PROC_NULL) u(:,1:nb,:,:)               = recv_left(:,:,:,:)
!        if(procxr .ne. MPI_PROC_NULL) u(:,nxb+nb+1:nxb+2*nb,:,:)  = recv_right(:,:,:,:)

      endif
#endif /* SHEAR_MPI */
    case ('ydim')
      if(pysize .gt. 1) then
!        allocate(send_left(nu,nx,nb,nz),send_right(nu,nx,nb,nz), &
!                 recv_left(nu,nx,nb,nz),recv_right(nu,nx,nb,nz))
!
!        if(procyl .ne. MPI_PROC_NULL) send_left(:,:,:,:)          =  u(:,:,nb+1:2*nb,:)
!        if(procyr .ne. MPI_PROC_NULL) send_right(:,:,:,:)         =  u(:,:,nyb+1:nyb+nb,:)
!
!        CALL MPI_ISEND   (send_left , nu*nx*nb*nz, MPI_DOUBLE_PRECISION, procyl, 30, comm, req(1), ierr)
!        CALL MPI_ISEND   (send_right, nu*nx*nb*nz, MPI_DOUBLE_PRECISION, procyr, 40, comm, req(3), ierr)
!        CALL MPI_IRECV   (recv_left , nu*nx*nb*nz, MPI_DOUBLE_PRECISION, procyl, 40, comm, req(2), ierr)
!        CALL MPI_IRECV   (recv_right, nu*nx*nb*nz, MPI_DOUBLE_PRECISION, procyr, 30, comm, req(4), ierr)

        CALL MPI_ISEND   (u(1,1,1,1), 1, MPI_XZ_LEFT_DOM,  procyl, 30, comm3d, req(1), ierr)
        CALL MPI_ISEND   (u(1,1,1,1), 1, MPI_XZ_RIGHT_DOM, procyr, 40, comm3d, req(3), ierr)
        CALL MPI_IRECV   (u(1,1,1,1), 1, MPI_XZ_LEFT_BND,  procyl, 40, comm3d, req(2), ierr)
        CALL MPI_IRECV   (u(1,1,1,1), 1, MPI_XZ_RIGHT_BND, procyr, 30, comm3d, req(4), ierr)

        do ireq=1,4
          call MPI_WAIT(req(ireq),status(1,ireq),ierr)
        enddo

!        if(procyl .ne. MPI_PROC_NULL) u(:,:,1:nb,:)               = recv_left(:,:,:,:)
!        if(procyr .ne. MPI_PROC_NULL) u(:,:,nyb+nb+1:nyb+2*nb,:)  = recv_right(:,:,:,:)

!        deallocate(send_left,send_right,recv_left,recv_right)
      endif

    case ('zdim')
      if(pzsize .gt. 1) then
!        allocate(send_left(nu,nx,ny,nb),send_right(nu,nx,ny,nb), &
!                 recv_left(nu,nx,ny,nb),recv_right(nu,nx,ny,nb))

!        if(proczl .ne. MPI_PROC_NULL) send_left(:,:,:,:)          =  u(:,:,:,nb+1:2*nb)
!        if(proczr .ne. MPI_PROC_NULL) send_right(:,:,:,:)         =  u(:,:,:,nzb+1:nzb+nb)

!        CALL MPI_ISEND   (send_left , nu*nx*ny*nb, MPI_DOUBLE_PRECISION, proczl, 50, comm, req(1), ierr)
!        CALL MPI_ISEND   (send_right, nu*nx*ny*nb, MPI_DOUBLE_PRECISION, proczr, 60, comm, req(3), ierr)
!        CALL MPI_IRECV   (recv_left , nu*nx*ny*nb, MPI_DOUBLE_PRECISION, proczl, 60, comm, req(2), ierr)
!        CALL MPI_IRECV   (recv_right, nu*nx*ny*nb, MPI_DOUBLE_PRECISION, proczr, 50, comm, req(4), ierr)
        CALL MPI_ISEND   (u(1,1,1,1), 1, MPI_XY_LEFT_DOM,  proczl, 50, comm3d, req(1), ierr)
        CALL MPI_ISEND   (u(1,1,1,1), 1, MPI_XY_RIGHT_DOM, proczr, 60, comm3d, req(3), ierr)
        CALL MPI_IRECV   (u(1,1,1,1), 1, MPI_XY_LEFT_BND,  proczl, 60, comm3d, req(2), ierr)
        CALL MPI_IRECV   (u(1,1,1,1), 1, MPI_XY_RIGHT_BND, proczr, 50, comm3d, req(4), ierr)

        do ireq=1,4
          call MPI_WAIT(req(ireq),status(1,ireq),ierr)
        enddo

!        if(proczl .ne. MPI_PROC_NULL) u(:,:,:,1:nb)               = recv_left(:,:,:,:)
!        if(proczr .ne. MPI_PROC_NULL) u(:,:,:,nzb+nb+1:nzb+2*nb)  = recv_right(:,:,:,:)

!        deallocate(send_left,send_right,recv_left,recv_right)
      endif
  end select ! (dim)

! MPI + non-MPI corner-periodic boundary condition

  if(bnd_xl .eq. 'cor') then
!   - lower to left
    if(pcoords(1) .eq. 0 .and. pcoords(2) .eq. 0) then
      do i=1,nb
        do j=nb+1,ny
          u(idna,i,j,:) =  u(idna,j,2*nb+1-i,:)
          u(imxa,i,j,:) = -u(imya,j,2*nb+1-i,:)
          u(imya,i,j,:) =  u(imxa,j,2*nb+1-i,:)
          u(imza,i,j,:) =  u(imza,j,2*nb+1-i,:)
#ifndef ISO
          u(iena,i,j,:) =  u(iena,j,2*nb+1-i,:)
#endif /* ISO */
#ifdef COSM_RAYS
          u(iecr,i,j,:) =  u(iecr,j,2*nb+1-i,:)
#endif /* COSM_RAYS */
        enddo
      enddo
    endif

    if(procxyl .gt. 0) then
      allocate(send_left(nu,nb,ny,nz), recv_left(nu,nx,nb,nz))

      send_left(:,:,:,:) = u(:,nb+1:2*nb,:,:)

      CALL MPI_ISEND   (send_left , nu*nb*ny*nz, MPI_DOUBLE_PRECISION, procxyl, 70, comm, req(1), ierr)
      CALL MPI_IRECV   (recv_left , nu*nx*nb*nz, MPI_DOUBLE_PRECISION, procxyl, 80, comm, req(2), ierr)

      do ireq=1,2
        call MPI_WAIT(req(ireq),status(1,ireq),ierr)
      enddo

      do i=1,nb
        do j=1,ny
          u(idna,i,j,:) =  recv_left(idna,j,nb+1-i,:)
          u(imxa,i,j,:) = -recv_left(imya,j,nb+1-i,:)
          u(imya,i,j,:) =  recv_left(imxa,j,nb+1-i,:)
          u(imza,i,j,:) =  recv_left(imza,j,nb+1-i,:)
#ifndef ISO
          u(iena,i,j,:) =  recv_left(iena,j,nb+1-i,:)
#endif /* ISO */
#ifdef COSM_RAYS
          u(iecr,i,j,:) =  recv_left(iecr,j,nb+1-i,:)
#endif /* COSM_RAYS */
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
          u(idna,i,j,:) =  u(idna,2*nb+1-j,i,:)
          u(imxa,i,j,:) =  u(imya,2*nb+1-j,i,:)
          u(imya,i,j,:) = -u(imxa,2*nb+1-j,i,:)
          u(imza,i,j,:) =  u(imza,2*nb+1-j,i,:)
#ifndef ISO
          u(iena,i,j,:) =  u(iena,2*nb+1-j,i,:)
#endif /* ISO */
#ifdef COSM_RAYS
          u(iecr,i,j,:) =  u(iecr,2*nb+1-j,i,:)
#endif /* COSM_RAYS */
        enddo
      enddo
!   - interior to corner
      do j=1,nb
        do i=1,nb
          u(idna,i,j,:) =   u(idna,2*nb+1-i,2*nb+1-j,:)
          u(imxa,i,j,:) =  -u(imxa,2*nb+1-i,2*nb+1-j,:)
          u(imya,i,j,:) =  -u(imya,2*nb+1-i,2*nb+1-j,:)
          u(imza,i,j,:) =   u(imza,2*nb+1-i,2*nb+1-j,:)
#ifndef ISO
          u(iena,i,j,:) =   u(iena,2*nb+1-i,2*nb+1-j,:)
#endif /* ISO */
#ifdef COSM_RAYS
          u(iecr,i,j,:) =   u(iecr,2*nb+1-i,2*nb+1-j,:)
#endif /* COSM_RAYS */
        enddo
      enddo
    endif

    if(procyxl .gt. 0) then
      allocate(send_left(nu,nx,nb,nz), recv_left(nu,nb,ny,nz))

      send_left(:,:,:,:) = u(:,:,nb+1:2*nb,:)

      CALL MPI_ISEND   (send_left , nu*nx*nb*nz, MPI_DOUBLE_PRECISION, procyxl, 80, comm, req(1), ierr)
      CALL MPI_IRECV   (recv_left , nu*nb*ny*nz, MPI_DOUBLE_PRECISION, procyxl, 70, comm, req(2), ierr)

      do ireq=1,2
        call MPI_WAIT(req(ireq),status(1,ireq),ierr)
      enddo

      do j=1,nb
        do i=1,nx
          u(idna,i,j,:) =  recv_left(idna,nb+1-j,i,:)
          u(imxa,i,j,:) =  recv_left(imya,nb+1-j,i,:)
          u(imya,i,j,:) = -recv_left(imxa,nb+1-j,i,:)
          u(imza,i,j,:) =  recv_left(imza,nb+1-j,i,:)
#ifndef ISO
          u(iena,i,j,:) =  recv_left(iena,nb+1-j,i,:)
#endif /* ISO */
#ifdef COSM_RAYS
          u(iecr,i,j,:) =  recv_left(iecr,nb+1-j,i,:)
#endif /* COSM_RAYS */
        enddo
      enddo

      deallocate(send_left,recv_left)
    endif
  endif


! Non-MPI boundary conditions

#ifdef SHEAR_MY
   if( (bnd_xl == 'she').and.(bnd_xr == 'she')) then   ! 2d ONLY !!!!!!!
      allocate(temp(nxd,ny,nz))
      allocate(tem2(nxd,ny,nz))
      do i = 1, nu
         if( i == imya ) then
           do j = 1,nx
             u(i,j,:,:) = u(i,j,:,:) + qshear*omega*x(j)*u(1,j,:,:)
           enddo
         endif
#ifndef ISO
         if( i == iena ) then
           do j = 1,nx
             u(i,j,:,:) = u(i,j,:,:) - 0.5*(qshear*omega*x(j))**2 * u(1,j,:,:)
           enddo
         endif
#endif /* ~ISO */
         temp(:,:,:) = unshear(u(i,nb+1:nxd+nb,:,:),x(nb+1:nxd+nb),.true.)
         tem2(:,:,:) = unshear(temp(:,:,:),x(nb+1:nxd+nb),.true.)

         u(i,1:nb,:,:) = tem2(nxd-nb+1:nxd,:,:)
         u(i,nxd+nb+1:nxd+2*nb,:,:) = tem2(1:nb,:,:)
         if( i == imya ) then
           do j = 1,nx
             u(i,j,:,:) = u(i,j,:,:) - qshear*omega*x(j) * u(1,j,:,:)
           enddo
         endif
#ifndef ISO
         if( i == iena ) then
           do j = 1,nx
             u(i,j,:,:) = u(i,j,:,:) + 0.5*(qshear*omega*x(j))**2  * u(1,j,:,:)
           enddo
         endif
#endif /* ~ISO */
      enddo
      deallocate(temp,tem2)
   endif
   u(idna,:,:,:) = max(u(idna,:,:,:),smalld)
#endif /* SHEAR_MY */

  select case (dim)
    case ('xdim')

      select case (bnd_xl)
        case ('she')
!         Do nothing if 'mpi'
        case ('mpi')
!         Do nothing if 'mpi'
        case ('per')
          u(:,1:nb,:,:)                        = u(:,nxb+1:nxb+nb,:,:)
        case ('cor')
!         Do nothing if 'cor'
        case ('ref')
          do ib=1,nb

            u((/idna,imya,imza/),nb+1-ib,:,:)  = u((/idna,imya,imza/),nb+ib,:,:)
            u(imxa,nb+1-ib,:,:)                =-u(imxa,nb+ib,:,:)
#ifndef ISO
            u(iena,nb+1-ib,:,:)                = u(iena,nb+ib,:,:)
#endif /* ISO */
#ifdef COSM_RAYS
            u(iecr,nb+1-ib,:,:)                = u(iecr,nb+ib,:,:)
#endif /* COSM_RAYS */
          enddo
        case ('out')
          do ib=1,nb
            u(:,ib,:,:)                        = u(:,nb+1,:,:)
#ifdef COSM_RAYS
            u(iecr,ib,:,:)                     = smallecr
#endif /* COSM_RAYS */
          enddo
        case ('outd')
          do ib=1,nb

            u((/idna,imya,imza/),ib,:,:)       = u((/idna,imya,imza/),nb+1,:,:)
            u(imxa,ib,:,:)                     = min(u(imxa,nb+1,:,:),0.0)
#ifndef ISO
            u(iena,ib,:,:)                     = u(iena,nb+1,:,:)
#endif /* ISO */
#ifdef COSM_RAYS
            u(iecr,ib,:,:)                     = smallecr
#endif /* COSM_RAYS */
          enddo
        case ('shef')
!         Do nothing if 'mpi'
        case default
          write(*,*) 'Boundary condition ',bnd_xl,' not implemented in ',dim
      end select  ! (bnd_xl)

      select case (bnd_xr)
        case ('shef')
!         Do nothing if 'mpi'
        case ('she')
!         Do nothing if 'mpi'
        case ('mpi')
!         Do nothing if mpi
        case ('per')
          u(:,nxb+nb+1:nxb+2*nb,:,:)            = u(:,nb+1:2*nb,:,:)
        case ('ref')
          do ib=1,nb

            u((/idna,imya,imza/),nb+nxb+ib,:,:) = u((/idna,imya,imza/),nb+nxb+1-ib,:,:)
            u(imxa,nb+nxb+ib,:,:)               =-u(imxa,nb+nxb+1-ib,:,:)
#ifndef ISO
            u(iena,nb+nxb+ib,:,:)               = u(iena,nb+nxb+1-ib,:,:)
#endif /* ISO */
#ifdef COSM_RAYS
            u(iecr,nb+nxb+ib,:,:)               = u(iecr,nb+nxb+1-ib,:,:)
#endif /* COSM_RAYS */
          enddo
        case ('out')
          do ib=1,nb
            u(:,nb+nxb+ib,:,:)                  = u(:,nb+nxb,:,:)
#ifdef COSM_RAYS
            u(iecr,nb+nxb+ib,:,:)               = smallecr
#endif /* COSM_RAYS */
          enddo
        case ('outd')
          do ib=1,nb

            u((/idna,imya,imza/),nb+nxb+ib,:,:) = u((/idna,imya,imza/),nb+nxb,:,:)
            u(imxa,nb+nxb+ib,:,:)               = max(u(imxa,nb+nxb,:,:),0.0)
#ifndef ISO
            u(iena,nb+nxb+ib,:,:)               = u(iena,nb+nxb,:,:)
#endif /* ISO */
#ifdef COSM_RAYS
            u(iecr,nb+nxb+ib,:,:)               = smallecr
#endif /* COSM_RAYS */
          enddo
        case default
          write(*,*) 'Boundary condition ',bnd_xr,' not implemented in ',dim
      end select  ! (bnd_xr)


    case ('ydim')

      select case (bnd_yl)
        case ('mpi')
!         Do nothing if mpi
        case ('per')
          u(:,:,1:nb,:)                         = u(:,:,nyb+1:nyb+nb,:)
        case ('cor')
!         Do nothing if 'cor'
        case ('ref')
          do ib=1,nb

            u((/idna,imxa,imza/),:,nb+1-ib,:)   = u((/idna,imxa,imza/),:,nb+ib,:)
            u(imya,:,nb+1-ib,:)                 =-u(imya,:,nb+ib,:)
#ifndef ISO
            u(iena,:,nb+1-ib,:)                 = u(iena,:,nb+ib,:)
#endif /* ISO */
#ifdef COSM_RAYS
            u(iecr,:,nb+1-ib,:)                 = u(iecr,:,nb+ib,:)
#endif /* COSM_RAYS */
          enddo
        case ('out')
          do ib=1,nb
            u(:,:,ib,:)                         = u(:,:,nb+1,:)
#ifdef COSM_RAYS
            u(iecr,:,ib,:)                      = smallecr
#endif /* COSM_RAYS */
          enddo
        case ('outd')
          do ib=1,nb

            u((/idna,imxa,imza/),:,ib,:)        = u((/idna,imxa,imza/),:,nb+1,:)
            u(imya,:,ib,:)                      = min(u(imya,:,nb+1,:),0.0)
#ifndef ISO
            u(iena,:,ib,:)                      = u(iena,:,nb+1,:)
#endif /* ISO */
#ifdef COSM_RAYS
            u(iecr,:,ib,:)                      = smallecr
#endif /* COSM_RAYS */
          enddo
        case default
          write(*,*) 'Boundary condition ',bnd_yl,' not implemented in ',dim
      end select  ! (bnd_yl)

      select case (bnd_yr)
        case ('mpi')
!         Do nothing if mpi
        case ('per')
          u(:,:,nyb+nb+1:nyb+2*nb,:)            = u(:,:,nb+1:2*nb,:)
        case ('ref')
          do ib=1,nb

            u((/idna,imxa,imza/),:,nb+nyb+ib,:) = u((/idna,imxa,imza/),:,nb+nyb+1-ib,:)
            u(imya,:,nb+nyb+ib,:)               =-u(imya,:,nb+nyb+1-ib,:)
#ifndef ISO
            u(iena,:,nb+nyb+ib,:)               = u(iena,:,nb+nyb+1-ib,:)
#endif /* ISO */
#ifdef COSM_RAYS
            u(iecr,:,nb+nyb+ib,:)               = u(iecr,:,nb+nyb+1-ib,:)
#endif /* COSM_RAYS */
          enddo
        case ('out')
          do ib=1,nb
            u(:,:,nb+nyb+ib,:)                  = u(:,:,nb+nyb,:)
#ifdef COSM_RAYS
            u(iecr,:,nb+nyb+ib,:)               = smallecr
#endif /* COSM_RAYS */
          enddo
        case ('outd')
          do ib=1,nb

            u((/idna,imxa,imza/),:,nb+nyb+ib,:) = u((/idna,imxa,imza/),:,nb+nyb,:)
            u(imya,:,nb+nyb+ib,:)               = max(u(imya,:,nb+nyb,:),0.0)
#ifndef ISO
            u(iena,:,nb+nyb+ib,:)               = u(iena,:,nb+nyb,:)
#endif /* ISO */
#ifdef COSM_RAYS
            u(iecr,:,nb+nyb+ib,:)               = smallecr
#endif /* COSM_RAYS */
          enddo
        case default
          write(*,*) 'Boundary condition ',bnd_yr,' not implemented in ',dim

      end select  ! (bnd_yr)


    case ('zdim')

      select case (bnd_zl)
        case ('mpi')
!         Do nothing if mpi
        case ('per')
          u(:,:,:,1:nb)                         = u(:,:,:,nzb+1:nzb+nb)
        case ('ref')
          do ib=1,nb

            u((/idna,imxa,imya/),:,:,nb+1-ib)   = u((/idna,imxa,imya/),:,:,nb+ib)
            u(imza,:,:,nb+1-ib)                 =-u(imza,:,:,nb+ib)
#ifndef ISO
            u(iena,:,:,nb+1-ib)                 = u(iena,:,:,nb+ib)
#endif /* ISO */
#ifdef COSM_RAYS
            u(iecr,:,:,nb+1-ib)                 = u(iecr,:,:,nb+ib)
#endif /* COSM_RAYS */
          enddo
        case ('out')
          do ib=1,nb
            u(:,:,:,ib)                         = u(:,:,:,nb+1)
#ifdef COSM_RAYS
            u(iecr,:,:,ib)                      = smallecr
#endif /* COSM_RAYS */
          enddo
        case ('outd')
          do ib=1,nb

            u((/idna,imxa,imya/),:,:,ib)        = u((/idna,imxa,imya/),:,:,nb+1)
            u(imza,:,:,ib)                      = min(u(imza,:,:,nb+1),0.0)
#ifndef ISO
            u(iena,:,:,ib)                      = u(iena,:,:,nb+1)
#endif /* ISO */
#ifdef COSM_RAYS
            u(iecr,:,:,ib)                      = smallecr
#endif /* COSM_RAYS */
          enddo
#ifdef GRAV
        case ('outh')
          do ib=1,nb
            kb = nb+2-ib
            db = u(idna,:,:,kb)
            db = max(db,smalld)
#ifdef ISO
            csi2b = csi2
#else /* ISO */
            ekb= 0.5*(u(imxa,:,:,kb)**2+u(imya,:,:,kb)**2+u(imza,:,:,kb)**2)/db
            eib = u(iena,:,:,kb) - ekb
            eib = max(eib,smallei)
            csi2b = (gamma-1)*eib/db
#endif /* ISO */
            z1 = z(kb)
            z2 = z(kb-1)
            dzs = (z2-z1)/real(nsub)

            do ksub=1, nsub+1
              zs(ksub) = z1 + dzs/2 + (ksub-1)*dzs
            enddo

            do j=1,ny
              do i=1,nx

                call grav_accel('zsweep',i,j, zs, nsub, gprofs)
                gprofs=tune_zeq_bnd * gprofs

                dprofs(1) = db(i,j)
                do ksub=1, nsub
                  factor = (1.0 + 0.5*dzs*gprofs(ksub)/csi2b(i,j))  &
                               /(1.0 - 0.5*dzs*gprofs(ksub)/csi2b(i,j))
                  dprofs(ksub+1) = factor * dprofs(ksub)
                enddo

                db(i,j)  = dprofs(nsub+1)
                db(i,j)  = max(db(i,j), smalld)

                u(idna,i,j,kb-1)                =     db(i,j)
                u(imxa:imza,i,j,kb-1)           =     u(imxa:imza,i,j,kb)
! zakomentowac nastepna linie jesli warunek diodowy nie ma byc stosowany razem z hydrostatycznym
!                u(imza,i,j,kb-1)               =     min(u(imza,i,j,kb-1),0.0)
#ifndef ISO
                eib(i,j) = csi2b(i,j)*db(i,j)/(gamma-1)
                eib(i,j) = max(eib(i,j), smallei)
                u(iena,i,j,kb-1)                =     ekb(i,j) + eib(i,j)
#endif /* ISO */
#ifdef COSM_RAYS
                u(iecr,i,j,kb-1)                =     smallecr
#endif /* COSM_RAYS */
              enddo ! i
            enddo ! j
          enddo ! ib
#endif /* GRAV */
        case default
          write(*,*) 'Boundary condition ',bnd_zl,' not implemented in ',dim
      end select  ! (bnd_zl)

      select case (bnd_zr)
        case ('mpi')
!         Do nothing if mpi
        case ('per')
          u(:,:,:,nzb+nb+1:nzb+2*nb)            = u(:,:,:,nb+1:2*nb)
        case ('ref')
          do ib=1,nb

            u((/idna,imxa,imya/),:,:,nb+nzb+ib) = u((/idna,imxa,imya/),:,:,nb+nzb+1-ib)
            u(imza,:,:,nb+nzb+ib)               =-u(imza,:,:,nb+nzb+1-ib)
#ifndef ISO
            u(iena,:,:,nb+nzb+ib)               = u(iena,:,:,nb+nzb+1-ib)
#endif /* ISO */
#ifdef COSM_RAYS
            u(iecr,:,:,nb+nzb+ib)               = u(iecr,:,:,nb+nzb+1-ib)
#endif /* COSM_RAYS */
          enddo
        case ('out')
          do ib=1,nb
            u(:,:,:,nb+nzb+ib)                  = u(:,:,:,nb+nzb)
#ifdef COSM_RAYS
            u(iecr,:,:,nb+nzb+ib)               = smallecr
#endif /* COSM_RAYS */
          enddo
        case ('outd')
          do ib=1,nb

            u((/idna,imxa,imya/),:,:,nb+nzb+ib) = u((/idna,imxa,imya/),:,:,nb+nzb)
            u(imza,:,:,nb+nzb+ib)               = max(u(imza,:,:,nb+nzb),0.0)
#ifndef ISO
            u(iena,:,:,nb+nzb+ib)               = u(iena,:,:,nb+nzb)
#endif /* ISO */
#ifdef COSM_RAYS
            u(iecr,:,:,nb+nzb+ib)               = smallecr
#endif /* COSM_RAYS */
          enddo
#ifdef GRAV
        case ('outh')
          do ib=1,nb
            kb = nb+nzb-1+ib
            db = u(idna,:,:,kb)
            db = max(db,smalld)
#ifdef ISO
            csi2b = csi2
#else /* ISO */
            ekb= 0.5*(u(imxa,:,:,kb)**2+u(imya,:,:,kb)**2+u(imza,:,:,kb)**2)/db
            eib = u(iena,:,:,kb) - ekb
            eib = max(eib,smallei)
            csi2b = (gamma-1)*eib/db
#endif /* ISO */
            z1 = z(kb)
            z2 = z(kb+1)
            dzs = (z2-z1)/real(nsub)

            do ksub=1, nsub+1
              zs(ksub) = z1 + dzs/2 + (ksub-1)*dzs
            enddo

            do j=1,ny
              do i=1,nx

                call grav_accel('zsweep',i,j, zs, nsub, gprofs)
                gprofs=tune_zeq_bnd * gprofs

                dprofs(1) = db(i,j)
                do ksub=1, nsub
                  factor = (1.0 + 0.5*dzs*gprofs(ksub)/csi2b(i,j))  &
                               /(1.0 - 0.5*dzs*gprofs(ksub)/csi2b(i,j))
                  dprofs(ksub+1) = factor * dprofs(ksub)
                enddo

                db(i,j)  = dprofs(nsub+1)
                db(i,j)  = max(db(i,j), smalld)

                u(idna,i,j,kb+1)           =     db(i,j)
                u(imxa:imza,i,j,kb+1)      =     u(imxa:imza,i,j,kb)
! zakomentowac nastepna linie jesli warunek diodowy nie ma byc stosowany razem z hydrostatycznym
!                u(imza,i,j,kb+1)           =     max(u(imza,i,j,kb+1),0.0)
#ifndef ISO
                eib(i,j) = csi2b(i,j)*db(i,j)/(gamma-1)
                eib(i,j) = max(eib(i,j), smallei)
                u(iena,i,j,kb+1)           =     ekb(i,j) + eib(i,j)
#endif /* ISO */
#ifdef COSM_RAYS
                u(iecr,i,j,kb+1)           =     smallecr
#endif /* COSM_RAYS */
              enddo ! i
            enddo ! j
          enddo ! ib
#endif /* GRAV */
        case default
          write(*,*) 'Boundary condition ',bnd_zr,' not implemented in ',dim
      end select  ! (bnd_zr)

    end select  ! (dim)



end subroutine bnd_u

  subroutine compute_u_bnd
   use start,  only : dimensions
#ifndef SPLIT
   use arrays, only : Lu
#endif /* SPLIT */
   implicit none
#ifndef FLX_BND
   call bnd_u('xdim')
   call bnd_u('ydim')
   if(dimensions .eq. '3d') call bnd_u('zdim')
#endif /* ~FLX_BND */

#ifndef SPLIT
   Lu(:,:,:,: ) =  0.0
#endif /* SPLIT */
  end subroutine compute_u_bnd

end module fluid_boundaries
