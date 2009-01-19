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

  use fluidindex, only : nvar, iarr_all_dn,iarr_all_mx,iarr_all_my,iarr_all_mz
#ifndef ISO
  use fluidindex, only : iarr_all_en
  use start,  only : gamma
#endif /* ISO */
  use fluidindex, only : nvar, nfluid
  use arrays, only : x,z,nzb,nyb,nxb,nx,ny,nz, & 
      u, b,  bndxrar, bndyrar

#ifndef ISO
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
  real, dimension(nfluid,nx,ny) :: db, csi2b
#ifndef ISO
  real, dimension(nfluid,nx,ny) :: ekb, eib
  integer ifluid
#endif /* ISO */
  real, dimension(nsub+1):: zs, gprofs
  real, dimension(nfluid,nsub+1) :: dprofs
  real, dimension(nfluid) :: factor
  real dzs,z1,z2
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
        allocate(send_right(nvar,nb,ny,nz), send_left(nvar,nb,ny,nz), &
                 recv_left(nvar,nb,ny,nz), recv_right(nvar,nb,ny,nz) )
        send_left(:,:,:,:)          =  u(:,nb+1:2*nb,:,:)
        send_right(:,:,:,:)         =  u(:,nxb+1:nxb+nb,:,:)
!
! odejmujemy ped_y i energie odpowiadajace niezaburzonej rozniczkowej rotacji na lewym brzegu
!
        if(bnd_xl == 'she') then
          do i=1,nb
            send_left (iarr_all_my,i,:,:) = send_left(iarr_all_my,i,:,:) &
                                         +qshear*omega * x(nb+i)     * send_left(iarr_all_dn,i,:,:)
#ifndef ISO
            send_left (iarr_all_en,i,:,:) = send_left(iarr_all_en,i,:,:) &
                                    -0.5*(qshear*omega * x(nb+i))**2 * send_left(iarr_all_dn,i,:,:)
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
            send_right(iarr_all_my,i,:,:) = send_right(iarr_all_my,i,:,:) &
                                         +qshear*omega * x(nxb+i)     * send_right(iarr_all_dn,i,:,:)
#ifndef ISO
            send_right(iarr_all_en,i,:,:) = send_right(iarr_all_en,i,:,:) &
                                    -0.5*(qshear*omega * x(nxb+i))**2 * send_right(iarr_all_dn,i,:,:)

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
        CALL MPI_ISEND   (send_left , nvar*ny*nz*nb, MPI_DOUBLE_PRECISION, procxl, 10, comm, req(1), ierr)
        CALL MPI_ISEND   (send_right, nvar*ny*nz*nb, MPI_DOUBLE_PRECISION, procxr, 20, comm, req(3), ierr)
        CALL MPI_IRECV   (recv_left , nvar*ny*nz*nb, MPI_DOUBLE_PRECISION, procxl, 20, comm, req(2), ierr)
        CALL MPI_IRECV   (recv_right, nvar*ny*nz*nb, MPI_DOUBLE_PRECISION, procxr, 10, comm, req(4), ierr)

        do ireq = 1,4
          call MPI_WAIT(req(ireq),status(1,ireq),ierr)
        enddo

!
! dodajemy ped_y i energie odpowiadajace niezaburzonej rozniczkowej rotacji na prawym brzegu
!
        if(bnd_xr == 'she') then
          do i=1,nb
#ifndef ISO
             recv_right (iarr_all_en,i,:,:) = recv_right (iarr_all_en,i,:,:) &
                                      +0.5*(qshear*omega * x(nb+nxb+i))**2 * recv_right(iarr_all_dn,i,:,:)
#endif /* ISO */
             recv_right (iarr_all_my,i,:,:) = recv_right (iarr_all_my,i,:,:) &
                                           -qshear*omega * x(nb+nxb+i)     * recv_right(iarr_all_dn,i,:,:)
          enddo

        endif !(bnd_xr == 'she')
!
! dodajemy ped_y i energie odpowiadajace niezaburzonej rozniczkowej rotacji na lewym brzegu
!
        if(bnd_xl == 'she') then

          do i=1,nb
#ifndef ISO
             recv_left(iarr_all_en,i,:,:) = recv_left(iarr_all_en,i,:,:) &
                                    +0.5*(qshear*omega * x(i))**2 * recv_left(iarr_all_dn,i,:,:)
#endif /* ISO */
             recv_left(iarr_all_my,i,:,:) = recv_left(iarr_all_my,i,:,:) &
                                         -qshear*omega * x(i)     * recv_left(iarr_all_dn,i,:,:)
          enddo
        endif !(bnd_xl == 'she')

        u(:,1:nb,:,:)              = recv_left(:,1:nb,:,:)
        u(:,nxb+nb+1:nxb+2*nb,:,:) = recv_right(:,1:nb,:,:)

        u(iarr_all_dn,1:nb,:,:)              = max(u(iarr_all_dn,1:nb,:,:),smalld)
        u(iarr_all_dn,nxb+nb+1:nxb+2*nb,:,:) = max(u(iarr_all_dn,nxb+nb+1:nxb+2*nb,:,:),smalld)
      deallocate(send_left,send_right,recv_left,recv_right)
#else /* SHEAR_MPI */
      if(pxsize .gt. 1) then

        CALL MPI_ISEND   (u(1,1,1,1), 1, MPI_YZ_LEFT_DOM,  procxl, 10, comm3d, req(1), ierr)
        CALL MPI_ISEND   (u(1,1,1,1), 1, MPI_YZ_RIGHT_DOM, procxr, 20, comm3d, req(3), ierr)
        CALL MPI_IRECV   (u(1,1,1,1), 1, MPI_YZ_LEFT_BND,  procxl, 20, comm3d, req(2), ierr)
        CALL MPI_IRECV   (u(1,1,1,1), 1, MPI_YZ_RIGHT_BND, procxr, 10, comm3d, req(4), ierr)

        do ireq=1,4
          call MPI_WAIT(req(ireq),status(1,ireq),ierr)
        enddo

      endif
#endif /* SHEAR_MPI */
    case ('ydim')
      if(pysize .gt. 1) then

        CALL MPI_ISEND   (u(1,1,1,1), 1, MPI_XZ_LEFT_DOM,  procyl, 30, comm3d, req(1), ierr)
        CALL MPI_ISEND   (u(1,1,1,1), 1, MPI_XZ_RIGHT_DOM, procyr, 40, comm3d, req(3), ierr)
        CALL MPI_IRECV   (u(1,1,1,1), 1, MPI_XZ_LEFT_BND,  procyl, 40, comm3d, req(2), ierr)
        CALL MPI_IRECV   (u(1,1,1,1), 1, MPI_XZ_RIGHT_BND, procyr, 30, comm3d, req(4), ierr)

        do ireq=1,4
          call MPI_WAIT(req(ireq),status(1,ireq),ierr)
        enddo

      endif

    case ('zdim')
      if(pzsize .gt. 1) then

        CALL MPI_ISEND   (u(1,1,1,1), 1, MPI_XY_LEFT_DOM,  proczl, 50, comm3d, req(1), ierr)
        CALL MPI_ISEND   (u(1,1,1,1), 1, MPI_XY_RIGHT_DOM, proczr, 60, comm3d, req(3), ierr)
        CALL MPI_IRECV   (u(1,1,1,1), 1, MPI_XY_LEFT_BND,  proczl, 60, comm3d, req(2), ierr)
        CALL MPI_IRECV   (u(1,1,1,1), 1, MPI_XY_RIGHT_BND, proczr, 50, comm3d, req(4), ierr)

        do ireq=1,4
          call MPI_WAIT(req(ireq),status(1,ireq),ierr)
        enddo

      endif
  end select ! (dim)

! MPI + non-MPI corner-periodic boundary condition

  if(bnd_xl .eq. 'cor') then
!   - lower to left
    if(pcoords(1) .eq. 0 .and. pcoords(2) .eq. 0) then
      do i=1,nb
        do j=nb+1,ny
          u(iarr_all_dn,i,j,:) =  u(iarr_all_dn,j,2*nb+1-i,:)
          u(iarr_all_mx,i,j,:) = -u(iarr_all_my,j,2*nb+1-i,:)
          u(iarr_all_my,i,j,:) =  u(iarr_all_mx,j,2*nb+1-i,:)
          u(iarr_all_mz,i,j,:) =  u(iarr_all_mz,j,2*nb+1-i,:)
#ifndef ISO
          u(iarr_all_en,i,j,:) =  u(iarr_all_en,j,2*nb+1-i,:)
#endif /* ISO */
#ifdef COSM_RAYS
          u(iecr,i,j,:) =  u(iecr,j,2*nb+1-i,:)
#endif /* COSM_RAYS */
        enddo
      enddo
    endif

    if(procxyl .gt. 0) then
      allocate(send_left(nvar,nb,ny,nz), recv_left(nvar,nx,nb,nz))

      send_left(:,:,:,:) = u(:,nb+1:2*nb,:,:)

      CALL MPI_ISEND   (send_left , nvar*nb*ny*nz, MPI_DOUBLE_PRECISION, procxyl, 70, comm, req(1), ierr)
      CALL MPI_IRECV   (recv_left , nvar*nx*nb*nz, MPI_DOUBLE_PRECISION, procxyl, 80, comm, req(2), ierr)

      do ireq=1,2
        call MPI_WAIT(req(ireq),status(1,ireq),ierr)
      enddo

      do i=1,nb
        do j=1,ny
          u(iarr_all_dn,i,j,:) =  recv_left(iarr_all_dn,j,nb+1-i,:)
          u(iarr_all_mx,i,j,:) = -recv_left(iarr_all_my,j,nb+1-i,:)
          u(iarr_all_my,i,j,:) =  recv_left(iarr_all_mx,j,nb+1-i,:)
          u(iarr_all_mz,i,j,:) =  recv_left(iarr_all_mz,j,nb+1-i,:)
#ifndef ISO
          u(iarr_all_en,i,j,:) =  recv_left(iarr_all_en,j,nb+1-i,:)
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
          u(iarr_all_dn,i,j,:) =  u(iarr_all_dn,2*nb+1-j,i,:)
          u(iarr_all_mx,i,j,:) =  u(iarr_all_my,2*nb+1-j,i,:)
          u(iarr_all_my,i,j,:) = -u(iarr_all_mx,2*nb+1-j,i,:)
          u(iarr_all_mz,i,j,:) =  u(iarr_all_mz,2*nb+1-j,i,:)
#ifndef ISO
          u(iarr_all_en,i,j,:) =  u(iarr_all_en,2*nb+1-j,i,:)
#endif /* ISO */
#ifdef COSM_RAYS
          u(iecr,i,j,:) =  u(iecr,2*nb+1-j,i,:)
#endif /* COSM_RAYS */
        enddo
      enddo
!   - interior to corner
      do j=1,nb
        do i=1,nb
          u(iarr_all_dn,i,j,:) =   u(iarr_all_dn,2*nb+1-i,2*nb+1-j,:)
          u(iarr_all_mx,i,j,:) =  -u(iarr_all_mx,2*nb+1-i,2*nb+1-j,:)
          u(iarr_all_my,i,j,:) =  -u(iarr_all_my,2*nb+1-i,2*nb+1-j,:)
          u(iarr_all_mz,i,j,:) =   u(iarr_all_mz,2*nb+1-i,2*nb+1-j,:)
#ifndef ISO
          u(iarr_all_en,i,j,:) =   u(iarr_all_en,2*nb+1-i,2*nb+1-j,:)
#endif /* ISO */
#ifdef COSM_RAYS
          u(iecr,i,j,:) =   u(iecr,2*nb+1-i,2*nb+1-j,:)
#endif /* COSM_RAYS */
        enddo
      enddo
    endif

    if(procyxl .gt. 0) then
      allocate(send_left(nvar,nx,nb,nz), recv_left(nvar,nb,ny,nz))

      send_left(:,:,:,:) = u(:,:,nb+1:2*nb,:)

      CALL MPI_ISEND   (send_left , nvar*nx*nb*nz, MPI_DOUBLE_PRECISION, procyxl, 80, comm, req(1), ierr)
      CALL MPI_IRECV   (recv_left , nvar*nb*ny*nz, MPI_DOUBLE_PRECISION, procyxl, 70, comm, req(2), ierr)

      do ireq=1,2
        call MPI_WAIT(req(ireq),status(1,ireq),ierr)
      enddo

      do j=1,nb
        do i=1,nx
          u(iarr_all_dn,i,j,:) =  recv_left(iarr_all_dn,nb+1-j,i,:)
          u(iarr_all_mx,i,j,:) =  recv_left(iarr_all_my,nb+1-j,i,:)
          u(iarr_all_my,i,j,:) = -recv_left(iarr_all_mx,nb+1-j,i,:)
          u(iarr_all_mz,i,j,:) =  recv_left(iarr_all_mz,nb+1-j,i,:)
#ifndef ISO
          u(iarr_all_en,i,j,:) =  recv_left(iarr_all_en,nb+1-j,i,:)
#endif /* ISO */
#ifdef COSM_RAYS
          u(iecr,i,j,:) =  recv_left(iecr,nb+1-j,i,:)
#endif /* COSM_RAYS */
        enddo
      enddo

      deallocate(send_left,recv_left)
    endif
  endif

!wolt - proba--------------------------------------------------

! MPI + non-MPI out-corner-periodic boundary condition

  if(bnd_yr .eq. 'cor') then
!   - upper to onecell outflow
    if(pcoords(2) .eq. psize(2)-1 .and. pcoords(1) .eq. psize(1)-1) then
      do i=1,nx
!        do j=ny-1,ny
          u(iarr_all_dn,i,ny,:) =  u(iarr_all_dn,i,ny-1,:)
          u(iarr_all_mx,i,ny,:) =  u(iarr_all_my,i,ny-1,:)
          u(iarr_all_my,i,ny,:) =  u(iarr_all_mx,i,ny-1,:)
          u(iarr_all_mz,i,ny,:) =  u(iarr_all_mz,i,ny-1,:)
#ifndef ISO
          u(iarr_all_en,i,ny,:) =  u(iarr_all_en,i,ny-1,:)
#endif /* ISO */
#ifdef COSM_RAYS
          u(iecr,i,ny,:) =  u(iecr,i,ny-1,:)
#endif /* COSM_RAYS */
!        enddo
      enddo
    endif

    if(procyxr .gt. 0) then
      allocate(send_left(nvar,nx,nb,nz), recv_left(nvar,nb,ny,nz))

      do j=1,nb
        send_left(:,:,nb+1-j,:)    =  u(:,:,ny+1-min(j,2),:)
        send_left(iarr_all_mx,:,nb+1-j,:) = -u(iarr_all_my,:,nx+1-min(j,2),:)
        send_left(iarr_all_my,:,nb+1-j,:) =  u(iarr_all_mx,:,nx+1-min(j,2),:)
      enddo

      CALL MPI_ISEND   (send_left , nvar*nx*nb*nz, MPI_DOUBLE_PRECISION, procyxr, 170, comm, req(1), ierr)
      CALL MPI_IRECV   (recv_left , nvar*nb*ny*nz, MPI_DOUBLE_PRECISION, procyxr, 180, comm, req(2), ierr)

      do ireq=1,2
        call MPI_WAIT(req(ireq),status(1,ireq),ierr)
      enddo

!      do i=1,nb
!        do j=1,ny
	   u(:,:,ny,:) = u(:,:,ny-1,:)
!          u(iarr_all_dn,i,j,:) =  recv_left(iarr_all_dn,j,nb+1-i,:)
!          u(iarr_all_mx,i,j,:) = -recv_left(iarr_all_my,j,nb+1-i,:)
!          u(iarr_all_my,i,j,:) =  recv_left(iarr_all_mx,j,nb+1-i,:)
!          u(iarr_all_mz,i,j,:) =  recv_left(iarr_all_mz,j,nb+1-i,:)
#ifndef ISO
!          u(iarr_all_en,i,j,:) =  recv_left(iarr_all_en,j,nb+1-i,:)
#endif /* ISO */
#ifdef COSM_RAYS
!          u(iecr,i,j,:) =  recv_left(iecr,j,nb+1-i,:)
#endif /* COSM_RAYS */
!        enddo
!      enddo

      deallocate(send_left,recv_left)
    endif
  endif

  if(bnd_xr .eq. 'cor') then
!   - upper to right
    if(pcoords(1) .eq. psize(1)-1 .and. pcoords(2) .eq. psize(2)-1 ) then
      do i=nx+1-nb,nx
        do j=1,ny+i-nx-1
          u(iarr_all_dn,i,j,:) =  u(iarr_all_dn,max(j,ny-1),i,:)
          u(iarr_all_mx,i,j,:) = -u(iarr_all_my,max(j,ny-1),i,:)
          u(iarr_all_my,i,j,:) =  u(iarr_all_mx,max(j,ny-1),i,:)
          u(iarr_all_mz,i,j,:) =  u(iarr_all_mz,max(j,ny-1),i,:)
#ifndef ISO
          u(iarr_all_en,i,j,:) =  u(iarr_all_en,max(j,ny-1),i,:)
#endif /* ISO */
#ifdef COSM_RAYS
          u(iecr,i,j,:) =  u(iecr,max(j,ny-1),i,:)
#endif /* COSM_RAYS */
        enddo
      enddo
    endif

    if(procxyr .gt. 0) then
      allocate(send_left(nvar,nb,ny,nz), recv_left(nvar,nx,nb,nz))

      send_left(:,:,:,:) = u(:,:,ny-nb+1:ny,:)

      CALL MPI_ISEND   (send_left , nvar*nb*ny*nz, MPI_DOUBLE_PRECISION, procxyr, 180, comm, req(1), ierr)
      CALL MPI_IRECV   (recv_left , nvar*nx*nb*nz, MPI_DOUBLE_PRECISION, procxyr, 170, comm, req(2), ierr)

      do ireq=1,2
        call MPI_WAIT(req(ireq),status(1,ireq),ierr)
      enddo

      do j=1,ny
        do i=1,nb
          u(iarr_all_dn,nx-nb+i,j,:) =  recv_left(iarr_all_dn,j,i,:)
          u(iarr_all_mx,nx-nb+i,j,:) =  recv_left(iarr_all_my,j,i,:)
          u(iarr_all_my,nx-nb+i,j,:) =  recv_left(iarr_all_mx,j,i,:)
          u(iarr_all_mz,nx-nb+i,j,:) =  recv_left(iarr_all_mz,j,i,:)
#ifndef ISO
          u(iarr_all_en,nx-nb+i,j,:) =  recv_left(iarr_all_en,j,i,:)
#endif /* ISO */
#ifdef COSM_RAYS
          u(iecr,nx-nb+i,j,:) =  recv_left(iecr,j,i,:)
#endif /* COSM_RAYS */
        enddo
      enddo

      deallocate(send_left,recv_left)
    endif
  endif

!---------------iflow bnd condit.---------------------------------


  if(bnd_yr .eq. 'inf') then
!   - upper to onecell outflow
    if(pcoords(2) .eq. psize(2)-1) then
      do j=1,nb
        u(:,:,ny-nb+j,:) =  bndyrar(:,:,j,:)
      enddo
    endif

  endif

  if(bnd_xr .eq. 'inf') then
!   - upper to right
    if(pcoords(1) .eq. psize(1)-1) then
      do i=1,nb
        u(:,nx-nb+i,:,:) =  bndxrar(:,i,:,:)
      enddo
    endif

  endif



!===============================================================

! Non-MPI boundary conditions

#ifdef SHEAR_MY
   if( (bnd_xl == 'she').and.(bnd_xr == 'she')) then   ! 2d ONLY !!!!!!!
      allocate(temp(nxd,ny,nz))
      allocate(tem2(nxd,ny,nz))
      do i = 1, nvar
         if( i == iarr_all_my ) then
           do j = 1,nx
             u(i,j,:,:) = u(i,j,:,:) + qshear*omega*x(j)*u(1,j,:,:)
           enddo
         endif
#ifndef ISO
         if( i == iarr_all_en ) then
           do j = 1,nx
             u(i,j,:,:) = u(i,j,:,:) - 0.5*(qshear*omega*x(j))**2 * u(1,j,:,:)
           enddo
         endif
#endif /* ~ISO */
         temp(:,:,:) = unshear(u(i,nb+1:nxd+nb,:,:),x(nb+1:nxd+nb),.true.)
         tem2(:,:,:) = unshear(temp(:,:,:),x(nb+1:nxd+nb),.true.)

         u(i,1:nb,:,:) = tem2(nxd-nb+1:nxd,:,:)
         u(i,nxd+nb+1:nxd+2*nb,:,:) = tem2(1:nb,:,:)
         if( i == iarr_all_my ) then
           do j = 1,nx
             u(i,j,:,:) = u(i,j,:,:) - qshear*omega*x(j) * u(1,j,:,:)
           enddo
         endif
#ifndef ISO
         if( i == iarr_all_en ) then
           do j = 1,nx
             u(i,j,:,:) = u(i,j,:,:) + 0.5*(qshear*omega*x(j))**2  * u(1,j,:,:)
           enddo
         endif
#endif /* ~ISO */
!         do k = 1, nz
!            temp(:,:,k:k) = unshear_fft(u(i,nxd+1:nxd+nb,nb+1:ny-nb,k:k),x(nxd+1:nxd+nb),bdely,.true.)
!            u(i,1:nb,nb+1:ny-nb,k:k) = unshear_fft(temp(:,:,k:k), x(1:nb),bdely)
!            temp(:,:,k:k) = unshear_fft(u(i,nb+1:2*nb,nb+1:ny-nb,k:k),x(nb+1:2*nb),bdely,.true.)
!            u(i,nxd+nb+1:nxd+2*nb,nb+1:ny-nb,k:k) = unshear_fft(temp(:,:,k:k),x(nxd+nb+1:nxd+2*nb),bdely)
!s         enddo
      enddo
      deallocate(temp,tem2)
   endif
   u(iarr_all_dn,:,:,:) = max(u(iarr_all_dn,:,:,:),smalld)
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
        case ('inf')
!         Do nothing if 'inf'
        case ('ref')
          do ib=1,nb

            u((/iarr_all_dn,iarr_all_my,iarr_all_mz/),nb+1-ib,:,:)  = u((/iarr_all_dn,iarr_all_my,iarr_all_mz/),nb+ib,:,:)
            u(iarr_all_mx,nb+1-ib,:,:)                =-u(iarr_all_mx,nb+ib,:,:)
#ifndef ISO
            u(iarr_all_en,nb+1-ib,:,:)                = u(iarr_all_en,nb+ib,:,:)
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

            u((/iarr_all_dn,iarr_all_my,iarr_all_mz/),ib,:,:)       = u((/iarr_all_dn,iarr_all_my,iarr_all_mz/),nb+1,:,:)
            u(iarr_all_mx,ib,:,:)                     = min(u(iarr_all_mx,nb+1,:,:),0.0)
#ifndef ISO
            u(iarr_all_en,ib,:,:)                     = u(iarr_all_en,nb+1,:,:)
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
        case ('cor')
!         Do nothing if 'cor'
        case ('inf')
!         Do nothing if 'inf'
        case ('ref')
          do ib=1,nb

            u((/iarr_all_dn,iarr_all_my,iarr_all_mz/),nb+nxb+ib,:,:) = u((/iarr_all_dn,iarr_all_my,iarr_all_mz/),nb+nxb+1-ib,:,:)
            u(iarr_all_mx,nb+nxb+ib,:,:)               =-u(iarr_all_mx,nb+nxb+1-ib,:,:)
#ifndef ISO
            u(iarr_all_en,nb+nxb+ib,:,:)               = u(iarr_all_en,nb+nxb+1-ib,:,:)
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

            u((/iarr_all_dn,iarr_all_my,iarr_all_mz/),nb+nxb+ib,:,:) = u((/iarr_all_dn,iarr_all_my,iarr_all_mz/),nb+nxb,:,:)
            u(iarr_all_mx,nb+nxb+ib,:,:)               = max(u(iarr_all_mx,nb+nxb,:,:),0.0)
#ifndef ISO
            u(iarr_all_en,nb+nxb+ib,:,:)               = u(iarr_all_en,nb+nxb,:,:)
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
        case ('inf')
!         Do nothing if 'inf'
        case ('ref')
          do ib=1,nb

            u((/iarr_all_dn,iarr_all_mx,iarr_all_mz/),:,nb+1-ib,:)   = u((/iarr_all_dn,iarr_all_mx,iarr_all_mz/),:,nb+ib,:)
            u(iarr_all_my,:,nb+1-ib,:)                 =-u(iarr_all_my,:,nb+ib,:)
#ifndef ISO
            u(iarr_all_en,:,nb+1-ib,:)                 = u(iarr_all_en,:,nb+ib,:)
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

            u((/iarr_all_dn,iarr_all_mx,iarr_all_mz/),:,ib,:)        = u((/iarr_all_dn,iarr_all_mx,iarr_all_mz/),:,nb+1,:)
            u(iarr_all_my,:,ib,:)                      = min(u(iarr_all_my,:,nb+1,:),0.0)
#ifndef ISO
            u(iarr_all_en,:,ib,:)                      = u(iarr_all_en,:,nb+1,:)
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
        case ('cor')
!         Do nothing if 'cor'
        case ('inf')
!         Do nothing if 'inf'
        case ('ref')
          do ib=1,nb

            u((/iarr_all_dn,iarr_all_mx,iarr_all_mz/),:,nb+nyb+ib,:) = u((/iarr_all_dn,iarr_all_mx,iarr_all_mz/),:,nb+nyb+1-ib,:)
            u(iarr_all_my,:,nb+nyb+ib,:)               =-u(iarr_all_my,:,nb+nyb+1-ib,:)
#ifndef ISO
            u(iarr_all_en,:,nb+nyb+ib,:)               = u(iarr_all_en,:,nb+nyb+1-ib,:)
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

            u((/iarr_all_dn,iarr_all_mx,iarr_all_mz/),:,nb+nyb+ib,:) = u((/iarr_all_dn,iarr_all_mx,iarr_all_mz/),:,nb+nyb,:)
            u(iarr_all_my,:,nb+nyb+ib,:)               = max(u(iarr_all_my,:,nb+nyb,:),0.0)
#ifndef ISO
            u(iarr_all_en,:,nb+nyb+ib,:)               = u(iarr_all_en,:,nb+nyb,:)
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

            u((/iarr_all_dn,iarr_all_mx,iarr_all_my/),:,:,nb+1-ib)   = u((/iarr_all_dn,iarr_all_mx,iarr_all_my/),:,:,nb+ib)
            u(iarr_all_mz,:,:,nb+1-ib)                 =-u(iarr_all_mz,:,:,nb+ib)
#ifndef ISO
            u(iarr_all_en,:,:,nb+1-ib)                 = u(iarr_all_en,:,:,nb+ib)
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

            u((/iarr_all_dn,iarr_all_mx,iarr_all_my/),:,:,ib)        = u((/iarr_all_dn,iarr_all_mx,iarr_all_my/),:,:,nb+1)
            u(iarr_all_mz,:,:,ib)                      = min(u(iarr_all_mz,:,:,nb+1),0.0)
#ifndef ISO
            u(iarr_all_en,:,:,ib)                      = u(iarr_all_en,:,:,nb+1)
#endif /* ISO */
#ifdef COSM_RAYS
            u(iecr,:,:,ib)                      = smallecr
#endif /* COSM_RAYS */
          enddo
#ifdef GRAV
        case ('outh')
          do ib=1,nb
            kb = nb+2-ib
            db = u(iarr_all_dn,:,:,kb)
            db = max(db,smalld)
#ifdef ISO
            csi2b = csi2
#else /* ISO */
            ekb= 0.5*(u(iarr_all_mx,:,:,kb)**2+u(iarr_all_my,:,:,kb)**2+u(iarr_all_mz,:,:,kb)**2)/db
            eib = u(iarr_all_en,:,:,kb) - ekb
            eib = max(eib,smallei)
            do ifluid=1,nfluid
               csi2b(ifluid,:,:) = (gamma(ifluid)-1.0)*eib(ifluid,:,:)/db(ifluid,:,:)
            enddo
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

                dprofs(:,1) = db(:,i,j)
                do ksub=1, nsub
                  factor = (1.0 + 0.5*dzs*gprofs(ksub)/csi2b(:,i,j))  &
                               /(1.0 - 0.5*dzs*gprofs(ksub)/csi2b(:,i,j))
                  dprofs(:,ksub+1) = factor * dprofs(:,ksub)
!                 if(i.eq.7 .and. j.eq.7) write(*,999) ksub, zs(ksub), dprofs(ksub)
                enddo

                db(:,i,j)  = dprofs(:,nsub+1)
                db(:,i,j)  = max(db(:,i,j), smalld)

                u(iarr_all_dn,i,j,kb-1)           =     db(:,i,j)
                u(iarr_all_mx,i,j,kb-1)           =     u(iarr_all_mx,i,j,kb)
                u(iarr_all_my,i,j,kb-1)           =     u(iarr_all_my,i,j,kb)
                u(iarr_all_mz,i,j,kb-1)           =     u(iarr_all_mz,i,j,kb)
! zakomentowac nastepna linie jesli warunek diodowy nie ma byc stosowany razem z hydrostatycznym
!                u(iarr_all_mz,i,j,kb-1)               =     min(u(iarr_all_mz,i,j,kb-1),0.0)
#ifndef ISO
                eib(:,i,j) = csi2b(:,i,j)*db(:,i,j)/(gamma-1)
                eib(:,i,j) = max(eib(:,i,j), smallei)
                u(iarr_all_en,i,j,kb-1)                =     ekb(:,i,j) + eib(:,i,j)
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

            u((/iarr_all_dn,iarr_all_mx,iarr_all_my/),:,:,nb+nzb+ib) = u((/iarr_all_dn,iarr_all_mx,iarr_all_my/),:,:,nb+nzb+1-ib)
            u(iarr_all_mz,:,:,nb+nzb+ib)               =-u(iarr_all_mz,:,:,nb+nzb+1-ib)
#ifndef ISO
            u(iarr_all_en,:,:,nb+nzb+ib)               = u(iarr_all_en,:,:,nb+nzb+1-ib)
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

            u((/iarr_all_dn,iarr_all_mx,iarr_all_my/),:,:,nb+nzb+ib) = u((/iarr_all_dn,iarr_all_mx,iarr_all_my/),:,:,nb+nzb)
            u(iarr_all_mz,:,:,nb+nzb+ib)               = max(u(iarr_all_mz,:,:,nb+nzb),0.0)
#ifndef ISO
            u(iarr_all_en,:,:,nb+nzb+ib)               = u(iarr_all_en,:,:,nb+nzb)
#endif /* ISO */
#ifdef COSM_RAYS
            u(iecr,:,:,nb+nzb+ib)               = smallecr
#endif /* COSM_RAYS */
          enddo
#ifdef GRAV
        case ('outh')
          do ib=1,nb
            kb = nb+nzb-1+ib
            db = u(iarr_all_dn,:,:,kb)
            db = max(db,smalld)
#ifdef ISO
            csi2b = csi2
#else /* ISO */
            ekb= 0.5*(u(iarr_all_mx,:,:,kb)**2+u(iarr_all_my,:,:,kb)**2+u(iarr_all_mz,:,:,kb)**2)/db
            eib = u(iarr_all_en,:,:,kb) - ekb
            eib = max(eib,smallei)
            do ifluid=1,nfluid
               csi2b(ifluid,:,:) = (gamma(ifluid)-1.0)*eib(ifluid,:,:)/db(ifluid,:,:)
            enddo
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

                dprofs(:,1) = db(:,i,j)
                do ksub=1, nsub
                  factor = (1.0 + 0.5*dzs*gprofs(ksub)/csi2b(:,i,j))  &
                               /(1.0 - 0.5*dzs*gprofs(ksub)/csi2b(:,i,j))
                  dprofs(:,ksub+1) = factor * dprofs(:,ksub)
!                 if(i.eq.7 .and. j.eq.7) write(*,999) ksub, zs(ksub), dprofs(ksub)
                enddo

                db(:,i,j)  = dprofs(:,nsub+1)
                db(:,i,j)  = max(db(:,i,j), smalld)

                u(iarr_all_dn,i,j,kb+1)      =     db(:,i,j)
                u(iarr_all_mx,i,j,kb+1)      =     u(iarr_all_mx,i,j,kb)
                u(iarr_all_my,i,j,kb+1)      =     u(iarr_all_my,i,j,kb)
                u(iarr_all_mz,i,j,kb+1)      =     u(iarr_all_mz,i,j,kb)
! zakomentowac nastepna linie jesli warunek diodowy nie ma byc stosowany razem z hydrostatycznym
!                u(iarr_all_mz,i,j,kb+1)           =     max(u(iarr_all_mz,i,j,kb+1),0.0)
#ifndef ISO
                eib(:,i,j) = csi2b(:,i,j)*db(:,i,j)/(gamma-1)
                eib(:,i,j) = max(eib(:,i,j), smallei)
                u(iarr_all_en,i,j,kb+1)           =     ekb(:,i,j) + eib(:,i,j)
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
