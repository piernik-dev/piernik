! $Id$
!
! PIERNIK Code Copyright (C) 2006 Michal Hanasz
!
!    This file is part of PIERNIK code.
!
!    PIERNIK is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    PIERNIK is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with PIERNIK.  If not, see <http://www.gnu.org/licenses/>.
!
!    Initial implementation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.h"

module magboundaries
! pulled by MAGNETIC
   implicit none

   private
   public :: bnd_a, bnd_b, bnd_emf, all_mag_boundaries

contains

   subroutine bnd_a(A)

      use constants,  only: MAG, xdim, zdim, LO, HI, BND, BLK
      use dataio_pub, only: die
      use mpisetup,   only: ierr, req, comm3d, procn, status, psize, have_mpi, is_mpi_noncart
      use grid,       only: cg
      use mpi,        only: MPI_COMM_NULL

      implicit none

      real, dimension(:,:,:,:) :: A
      integer                  :: i, itag, jtag

      if (have_mpi .and. is_mpi_noncart) call die("[magboundaries:bnd_a] is_mpi_noncart is not implemented") !procn, psize
      if (comm3d == MPI_COMM_NULL) call die("[magboundaries:bnd_a] comm3d == MPI_COMM_NULL")

      do i = xdim, zdim
         if (psize(i) > 1) then

            jtag = 20*i
            itag = jtag - 10
            call MPI_Isend(A(1,1,1,1), 1, cg%mbc(MAG, i, LO, BLK), procn(i,LO), itag, comm3d, req(1), ierr)
            call MPI_Isend(A(1,1,1,1), 1, cg%mbc(MAG, i, HI, BLK), procn(i,HI), jtag, comm3d, req(3), ierr)
            call MPI_Irecv(A(1,1,1,1), 1, cg%mbc(MAG, i, LO, BND), procn(i,LO), jtag, comm3d, req(2), ierr)
            call MPI_Irecv(A(1,1,1,1), 1, cg%mbc(MAG, i, HI, BND), procn(i,HI), itag, comm3d, req(4), ierr)

            call MPI_Waitall(4,req(:),status(:,:),ierr)
         endif
      enddo

   end subroutine bnd_a

   subroutine bnd_b(dir)

      use arrays,        only: b
      use constants,     only: MAG, xdim, ydim, zdim, LO, HI, BND, BLK, BND_MPI, BND_PER, BND_REF, BND_OUT, BND_OUTD, BND_OUTH, BND_COR, BND_SHE, BND_INF
      use dataio_pub,    only: msg, warn, die
      use fluidindex,    only: ibx, iby, ibz
      use grid,          only: cg
      use mpi,           only: MPI_DOUBLE_PRECISION, MPI_COMM_NULL
      use mpisetup,      only: ierr, req, comm3d, procn, proc, status, psize, procxyl, procyxl, pcoords, comm, master, have_mpi, is_mpi_noncart
#ifdef SHEAR
      use shear,         only: eps,delj
#endif /* SHEAR */

      implicit none

      integer, intent(in) :: dir
      integer           :: i, j, itag, jtag
      real, allocatable :: send_left(:,:,:,:),recv_left(:,:,:,:)
#ifdef SHEAR
      real, allocatable :: send_right(:,:,:,:),recv_right(:,:,:,:)
#endif /* SHEAR */
      logical, save                         :: frun = .true.
      logical, save                         :: bnd_xl_not_provided = .false.
      logical, save                         :: bnd_xr_not_provided = .false.
      logical, save                         :: bnd_yl_not_provided = .false.
      logical, save                         :: bnd_yr_not_provided = .false.
      logical, save                         :: bnd_zl_not_provided = .false.
      logical, save                         :: bnd_zr_not_provided = .false.

! MPI block comunication
      if (comm3d /= MPI_COMM_NULL) then

         if (have_mpi .and. is_mpi_noncart) call die("[magboundaries:bnd_b] is_mpi_noncart is not implemented") !procn, procxyl, procyxl, psize, pcoords
#ifdef SHEAR
         if (dir == xdim) then
            allocate(send_left(3, cg%nb, cg%ny, cg%nz),send_right(3, cg%nb, cg%ny, cg%nz), &
                 &   recv_left(3, cg%nb, cg%ny, cg%nz),recv_right(3, cg%nb, cg%ny, cg%nz))

            send_left (:,:,:,:)  = b(:, cg%is:cg%isb,:,:)
            send_right(:,:,:,:)  = b(:, cg%ieb:cg%ie,:,:)

            if (cg%bnd(xdim, LO) == BND_SHE) then
!
! przesuwamy o calkowita liczbe komorek + periodyczny wb w kierunku y
!
               send_left (:,:,  cg%js:cg%je,:) = cshift(send_left (:,:, cg%js :cg%je, :),dim=3,shift= delj)
               send_left (:,:,      1:cg%nb,:) =        send_left (:,:, cg%jeb:cg%je, :)
               send_left (:,:,cg%je+1:cg%ny,:) =        send_left (:,:, cg%js :cg%jsb,:)
!
! remapujemy  - interpolacja kwadratowa
!
               send_left (:,:,:,:)  = (1.+eps)*(1.-eps) * send_left (:,:,:,:) &
                    &                 -0.5*eps*(1.-eps) * cshift(send_left (:,:,:,:),shift=-1,dim=3) &
                    &                 +0.5*eps*(1.+eps) * cshift(send_left (:,:,:,:),shift= 1,dim=3)
            endif ! (cg%bnd(xdim, LO) == BND_SHE)

            if (cg%bnd(xdim, HI) == BND_SHE) then
!
! przesuwamy o calkowita liczbe komorek + periodyczny wb w kierunku y
!
               send_right (:,:,  cg%js:cg%je,:) = cshift(send_right(:,:, cg%js :cg%je, :),dim=3,shift=-delj)
               send_right (:,:,      1:cg%nb,:) =        send_right(:,:, cg%jeb:cg%je, :)
               send_right (:,:,cg%je+1:cg%ny,:) =        send_right(:,:, cg%js :cg%jsb,:)
!
! remapujemy - interpolacja kwadratowa
!
               send_right (:,:,:,:) = (1.+eps)*(1.-eps) * send_right (:,:,:,:) &
                    &                 -0.5*eps*(1.-eps) * cshift(send_right (:,:,:,:),shift= 1,dim=3) &
                    &                 +0.5*eps*(1.+eps) * cshift(send_right (:,:,:,:),shift=-1,dim=3)
            endif ! (cg%bnd(xdim, HI) == BND_SHE)
!
! wysylamy na drugi brzeg
!
            call MPI_Isend(send_left , 3*cg%ny*cg%nz*cg%nb, MPI_DOUBLE_PRECISION, procn(dir,LO), 10, comm, req(1), ierr)
            call MPI_Isend(send_right, 3*cg%ny*cg%nz*cg%nb, MPI_DOUBLE_PRECISION, procn(dir,HI), 20, comm, req(3), ierr)
            call MPI_Irecv(recv_left , 3*cg%ny*cg%nz*cg%nb, MPI_DOUBLE_PRECISION, procn(dir,LO), 20, comm, req(2), ierr)
            call MPI_Irecv(recv_right, 3*cg%ny*cg%nz*cg%nb, MPI_DOUBLE_PRECISION, procn(dir,HI), 10, comm, req(4), ierr)

            call MPI_Waitall(4,req(:),status(:,:),ierr)

            b(:,        1:cg%nb-1,:,:) = recv_left (:,  1:cg%nb-1,:,:)
            b(:,cg%ie+1+1:cg%nx,  :,:) = recv_right(:,1+1:cg%nb,  :,:)

            if (allocated(send_left))  deallocate(send_left)
            if (allocated(send_right)) deallocate(send_right)
            if (allocated(recv_left))  deallocate(recv_left)
            if (allocated(recv_right)) deallocate(recv_right)

!===============================================================================
         else
#endif /* SHEAR */

            if (psize(dir) > 1) then

               jtag = 20*dir
               itag = jtag - 10
               call MPI_Isend(b(1,1,1,1), 1, cg%mbc(MAG, dir, LO, BLK), procn(dir,LO), itag, comm3d, req(1), ierr)
               call MPI_Isend(b(1,1,1,1), 1, cg%mbc(MAG, dir, HI, BLK), procn(dir,HI), jtag, comm3d, req(3), ierr)
               call MPI_Irecv(b(1,1,1,1), 1, cg%mbc(MAG, dir, LO, BND), procn(dir,LO), jtag, comm3d, req(2), ierr)
               call MPI_Irecv(b(1,1,1,1), 1, cg%mbc(MAG, dir, HI, BND), procn(dir,HI), itag, comm3d, req(4), ierr)

               call MPI_Waitall(4,req(:),status(:,:),ierr)
            endif

#ifdef SHEAR
         endif
#endif /* SHEAR */

! MPI + non-MPI corner-periodic boundary condition

         if (cg%bnd(xdim, LO) == BND_COR) then
!   - lower to left
            if (pcoords(1) == 0 .and. pcoords(2) == 0) then
               do i=1, cg%nb
                  do j=cg%js, cg%ny
                     b(ibx,i,j,:) = -b(iby,j,cg%isb+1-i,:)
                     b(iby,i,j,:) =  b(ibx,j,cg%isb+1-i,:)
                     b(ibz,i,j,:) =  b(ibz,j,cg%isb+1-i,:)
                  enddo
               enddo
            endif

            if (procxyl > 0) then
               allocate(send_left(3, cg%nb, cg%ny, cg%nz), recv_left(3, cg%nx, cg%nb, cg%nz))

               send_left(:,:,:,:) = b(:, cg%is:cg%isb,:,:)

               call MPI_Isend(send_left, 3*cg%nb*cg%ny*cg%nz, MPI_DOUBLE_PRECISION, procxyl, 70, comm, req(1), ierr)
               call MPI_Irecv(recv_left, 3*cg%nx*cg%nb*cg%nz, MPI_DOUBLE_PRECISION, procxyl, 80, comm, req(2), ierr)

               call MPI_Waitall(2,req(:),status(:,:),ierr)

               do i=1, cg%nb
                  do j=1, cg%ny
                     b(ibx,i,j,:) = -recv_left(iby,j, cg%is-i,:)
                     b(iby,i,j,:) =  recv_left(ibx,j, cg%is-i,:)
                     b(ibz,i,j,:) =  recv_left(ibz,j, cg%is-i,:)
                  enddo
               enddo

               if (allocated(send_left)) deallocate(send_left)
               if (allocated(recv_left)) deallocate(recv_left)
            endif
         endif

         if (cg%bnd(ydim, LO) == BND_COR) then
!   - left to lower
            if (pcoords(2) == 0 .and. pcoords(1) == 0 ) then
               do j=1, cg%nb
                  do i=cg%is, cg%nx
                     b(ibx,i,j,:) =  b(iby,cg%isb+1-j,i,:)
                     b(iby,i,j,:) = -b(ibx,cg%isb+1-j,i,:)
                     b(ibz,i,j,:) =  b(ibz,cg%isb+1-j,i,:)
                  enddo
               enddo
!   - interior to corner
               do j=1, cg%nb
                  do i=1, cg%nb
                     b(ibx,i,j,:) =  -b(ibx,cg%isb+1-i,cg%jsb+1-j,:)
                     b(iby,i,j,:) =  -b(iby,cg%isb+1-i,cg%jsb+1-j,:)
                     b(ibz,i,j,:) =   b(ibz,cg%isb+1-i,cg%jsb+1-j,:)
                  enddo
               enddo
            endif

            if (procyxl > 0) then
               allocate(send_left(3, cg%nx, cg%nb, cg%nz), recv_left(3, cg%nb, cg%ny, cg%nz))

               send_left(:,:,:,:) = b(:,:, cg%js:cg%jsb,:)

               call MPI_Isend   (send_left , 3*cg%nx*cg%nb*cg%nz, MPI_DOUBLE_PRECISION, procyxl, 80, comm, req(1), ierr)
               call MPI_Irecv   (recv_left , 3*cg%nb*cg%ny*cg%nz, MPI_DOUBLE_PRECISION, procyxl, 70, comm, req(2), ierr)

               call MPI_Waitall(2,req(:),status(:,:),ierr)

               do j=1, cg%nb
                  do i=1, cg%nx
                     b(ibx,i,j,:) =  recv_left(iby, cg%js-j,i,:)
                     b(iby,i,j,:) = -recv_left(ibx, cg%js-j,i,:)
                     b(ibz,i,j,:) =  recv_left(ibz, cg%js-j,i,:)
                  enddo
               enddo

               if (allocated(send_left)) deallocate(send_left)
               if (allocated(recv_left)) deallocate(recv_left)
            endif
         endif

      endif

! Non-MPI boundary conditions
      if (frun) then
         bnd_xl_not_provided = any( [BND_COR, BND_INF, BND_REF, BND_MPI, BND_SHE] == cg%bnd(xdim, LO))
         bnd_xr_not_provided = any( [BND_COR, BND_INF, BND_REF, BND_MPI, BND_SHE] == cg%bnd(xdim, HI))
         bnd_yl_not_provided = any( [BND_COR, BND_INF, BND_REF, BND_MPI ] == cg%bnd(ydim, LO))
         bnd_yr_not_provided = any( [BND_COR, BND_INF, BND_REF, BND_MPI ] == cg%bnd(ydim, HI))
         bnd_zl_not_provided = any( [BND_INF, BND_REF, BND_MPI ] == cg%bnd(zdim, LO))
         bnd_zr_not_provided = any( [BND_INF, BND_REF, BND_MPI ] == cg%bnd(zdim, HI))
      endif

      if (dir==xdim .and. bnd_xl_not_provided .and. bnd_xr_not_provided) return  ! avoid triple case
      if (dir==ydim .and. bnd_yl_not_provided .and. bnd_yr_not_provided) return  ! avoid triple case
      if (dir==zdim .and. bnd_zl_not_provided .and. bnd_zr_not_provided) return  ! avoid triple case

      select case (dir)
         case (xdim)

            select case (cg%bnd(xdim, LO))
               case (BND_COR, BND_INF, BND_MPI, BND_REF, BND_SHE)
                  ! Do nothing
               case (BND_PER)
                  if (comm3d /= MPI_COMM_NULL) b(:,1:cg%nb,:,:)              = b(:, cg%ieb:cg%ie,:,:)
               case (BND_OUT, BND_OUTD, BND_OUTH)
                  b(:,1,:,:) = b(:,2,:,:)
               case default
                  write(msg,'(2(a,i3))') "[magboundaries:bnd_b]: Boundary condition ",cg%bnd(xdim, LO)," not implemented in ",dir
                  if (master) call warn(msg)
            end select  ! (cg%bnd(xdim, LO))

            select case (cg%bnd(xdim, HI))
               case (BND_COR, BND_INF, BND_MPI, BND_REF, BND_SHE)
                  ! Do nothing
               case (BND_PER)
                  if (comm3d /= MPI_COMM_NULL) b(:, cg%ie+1:cg%nx,:,:) = b(:, cg%is:cg%isb,:,:)
               case (BND_OUT, BND_OUTD, BND_OUTH)
                  b(:, cg%nx,:,:) = b(:, cg%nx-1,:,:)
               case default
                  write(msg,'(2(a,i3))') "[magboundaries:bnd_b]: Boundary condition ",cg%bnd(xdim, HI)," not implemented in ",dir
                  if (master) call warn(msg)
            end select  ! (cg%bnd(xdim, HI))

         case (ydim)

            select case (cg%bnd(ydim, LO))
               case (BND_COR, BND_INF, BND_MPI, BND_REF)
                  ! Do nothing
               case (BND_PER)
                  if (comm3d /= MPI_COMM_NULL) b(:,:,1:cg%nb,:)              = b(:,:, cg%jeb:cg%je,:)
               case (BND_OUT, BND_OUTD, BND_OUTH)
                  b(:,:,1,:) = b(:,:,2,:)
               case default
                  write(msg,'(2(a,i3))') "[magboundaries:bnd_b]: Boundary condition ",cg%bnd(ydim, LO)," not implemented in ",dir
                  if (master) call warn(msg)
            end select  ! (cg%bnd(ydim, LO))

            select case (cg%bnd(ydim, HI))
               case (BND_COR, BND_INF, BND_MPI, BND_REF)
                  ! Do nothing
               case (BND_PER)
                  if (comm3d /= MPI_COMM_NULL) b(:,:, cg%je+1:cg%ny,:) = b(:,:, cg%js:cg%jsb,:)
               case (BND_OUT, BND_OUTD, BND_OUTH)
                  b(:,:, cg%ny,:) = b(:,:, cg%ny-1,:)
               case default
                  write(msg,'(2(a,i3))') "[magboundaries:bnd_b]: Boundary condition ",cg%bnd(ydim, HI)," not implemented in ",dir
                  if (master) call warn(msg)

            end select  ! (cg%bnd(ydim, HI))

         case (zdim)

            select case (cg%bnd(zdim, LO))
               case (BND_MPI, BND_REF)
                  ! Do nothing
               case (BND_PER)
                  if (comm3d /= MPI_COMM_NULL) b(:,:,:,1:cg%nb)              = b(:,:,:, cg%keb:cg%ke)
               case (BND_OUT, BND_OUTD, BND_OUTH)
                  b(:,:,:,1) = b(:,:,:,2)
               case default
                  write(msg,'(2(a,i3))') "[magboundaries:bnd_b]: Boundary condition ",cg%bnd(zdim, LO)," not implemented in ",dir
                  if (master) call warn(msg)
            end select  ! (cg%bnd(zdim, LO))

            select case (cg%bnd(zdim, HI))
               case (BND_MPI, BND_REF)
                  ! Do nothing
               case (BND_PER)
                  if (comm3d /= MPI_COMM_NULL) b(:,:,:, cg%ke+1:cg%nz) = b(:,:,:, cg%ks:cg%ksb)
               case (BND_OUT, BND_OUTD, BND_OUTH)
                  b(:,:,:, cg%nz) = b(:,:,:, cg%nz-1)
               case default
                  write(msg,'(2(a,i3))') "[magboundaries:bnd_b]: Boundary condition ",cg%bnd(zdim, HI)," not implemented in ",dir
                  if (master) call warn(msg)
            end select  ! (cg%bnd(zdim, HI))

      end select  ! (dim)

   end subroutine bnd_b

!=====================================================================================================

   subroutine bnd_emf(var, name, dir)

      use constants,     only: xdim, ydim, zdim, LO, HI, BND_MPI, BND_PER, BND_REF, BND_OUT, BND_OUTD, BND_OUTH, BND_COR, BND_SHE, BND_INF
      use dataio_pub,    only: msg, warn
      use grid,          only: cg
      use mpisetup,      only: master

      implicit none

      real, dimension(:,:,:), intent(inout) :: var
      character(len=*), intent(in)          :: name
      integer, intent(in)                   :: dir
      real, dimension(cg%ny, cg%nz)                :: dvarx
      real, dimension(cg%nx, cg%nz)                :: dvary
      real, dimension(cg%nx, cg%ny)                :: dvarz
      integer                               :: ib
      logical, save                         :: frun = .true.
      logical, save                         :: bnd_xl_not_provided = .false.
      logical, save                         :: bnd_xr_not_provided = .false.
      logical, save                         :: bnd_yl_not_provided = .false.
      logical, save                         :: bnd_yr_not_provided = .false.
      logical, save                         :: bnd_zl_not_provided = .false.
      logical, save                         :: bnd_zr_not_provided = .false.
      integer                               :: ledge, redge, lnbcells, rnbcells, zndiff, rrbase
      real                                  :: bndsign

      if (frun) then
         bnd_xl_not_provided = any( [BND_COR, BND_INF, BND_PER, BND_MPI, BND_SHE] == cg%bnd(xdim, LO))
         bnd_xr_not_provided = any( [BND_COR, BND_INF, BND_PER, BND_MPI, BND_SHE] == cg%bnd(xdim, HI))
         bnd_yl_not_provided = any( [BND_COR, BND_INF, BND_PER, BND_MPI ] == cg%bnd(ydim, LO))
         bnd_yr_not_provided = any( [BND_COR, BND_INF, BND_PER, BND_MPI ] == cg%bnd(ydim, HI))
         bnd_zl_not_provided = any( [BND_INF, BND_PER, BND_MPI ] == cg%bnd(zdim, LO))
         bnd_zr_not_provided = any( [BND_INF, BND_PER, BND_MPI ] == cg%bnd(zdim, HI))
         frun = .false.
      endif

      if (dir==xdim .and. bnd_xl_not_provided .and. bnd_xr_not_provided) return  ! avoid triple case
      if (dir==ydim .and. bnd_yl_not_provided .and. bnd_yr_not_provided) return  ! avoid triple case
      if (dir==zdim .and. bnd_zl_not_provided .and. bnd_zr_not_provided) return  ! avoid triple case

      bndsign = huge(1.0); ledge=huge(1); redge=huge(1); lnbcells=huge(1); rnbcells=huge(1); zndiff=huge(1); rrbase=huge(1)
      ! the code below should not use these values and the compiler should not complain on possible use of uninitialized variables.

      select case (dir)
         case (xdim)

            if (any( [cg%bnd(xdim, LO) == BND_REF,  cg%bnd(xdim, HI) == BND_REF,  cg%bnd(xdim, LO) == BND_OUT,  cg%bnd(xdim, HI) == BND_OUT, &
               &      cg%bnd(xdim, LO) == BND_OUTD, cg%bnd(xdim, HI) == BND_OUTD, cg%bnd(xdim, LO) == BND_OUTH, cg%bnd(xdim, HI) == BND_OUTH] )) then
               select case (name)
                  case ("vxby","vxbz")
                     call compute_bnd_indxs(1, cg%nxb,ledge,redge,lnbcells,rnbcells,bndsign,zndiff,rrbase)
                  case ("vybx","vzbx","emfy","emfz")
                     call compute_bnd_indxs(2, cg%nxb,ledge,redge,lnbcells,rnbcells,bndsign,zndiff,rrbase)
                  case ("vybz","vzby","emfx")
                     call compute_bnd_indxs(3, cg%nxb,ledge,redge,lnbcells,rnbcells,bndsign,zndiff,rrbase)
               end select  ! (name)
            endif

            select case (cg%bnd(xdim, LO))
               case (BND_COR, BND_INF, BND_MPI, BND_PER, BND_SHE)
                  ! Do nothing
               case (BND_REF)
                  if (zndiff == 1) var(ledge,:,:) = 0.0
                  do ib=1,lnbcells
                     var(lnbcells+1-ib,:,:) = bndsign * var(ledge+ib,:,:)
                  enddo
               case (BND_OUT, BND_OUTD, BND_OUTH)
                  ledge = ledge + 1 ; lnbcells = lnbcells + 1
                  dvarx = var(ledge+1,:,:)-var(ledge,:,:)
                  do ib=1,lnbcells
                     var(ib,:,:) = var(cg%is,:,:) - real(cg%is-ib)*dvarx
                  enddo
               case default
                  write(msg,'(a,i3,3a,i3)') "[magboundaries:bnd_emf]: Boundary condition ",cg%bnd(xdim, LO)," not implemented for ",name, " in ", dir
                  if (master) call warn(msg)
            end select  ! (cg%bnd(xdim, LO))

            select case (cg%bnd(xdim, HI))
               case (BND_COR, BND_INF, BND_MPI, BND_PER, BND_SHE)
                  ! Do nothing
               case (BND_REF)
                  if (zndiff == 1) var(redge,:,:) = 0.0
                  do ib=1,rnbcells
                     var(redge+ib,:,:) = bndsign * var(rrbase-ib,:,:)
                  enddo
               case (BND_OUT, BND_OUTD, BND_OUTH)
!                  dvarx = var(rrbase,:,:)-var(rrbase-1,:,:) original
                  dvarx = var(redge,:,:)-var(redge-1,:,:)
                  do ib=1,rnbcells
                     var(redge+ib,:,:) = var(redge,:,:) + real(ib)*dvarx
                  enddo
               case default
                  write(msg,'(a,i3,3a,i3)') "[magboundaries:bnd_emf]: Boundary condition ",cg%bnd(xdim, HI)," not implemented for ",name, " in ", dir
                  if (master) call warn(msg)
            end select ! (cg%bnd(xdim, HI))

         case (ydim)

            if (any( [cg%bnd(ydim, LO) == BND_REF,  cg%bnd(ydim, HI) == BND_REF,  cg%bnd(ydim, LO) == BND_OUT,  cg%bnd(ydim, HI) == BND_OUT, &
               &      cg%bnd(ydim, LO) == BND_OUTD, cg%bnd(ydim, HI) == BND_OUTD, cg%bnd(ydim, LO) == BND_OUTH, cg%bnd(ydim, HI) == BND_OUTH] )) then
               select case (name)
                  case ("vybz","vybx")
                     call compute_bnd_indxs(1, cg%nyb,ledge,redge,lnbcells,rnbcells,bndsign,zndiff,rrbase)
                  case ("vzby","vxby","emfz","emfx")
                     call compute_bnd_indxs(2, cg%nyb,ledge,redge,lnbcells,rnbcells,bndsign,zndiff,rrbase)
                  case ("vzbx","vxbz","emfy")
                     call compute_bnd_indxs(3, cg%nyb,ledge,redge,lnbcells,rnbcells,bndsign,zndiff,rrbase)
               end select  ! (name)
            endif

            select case (cg%bnd(ydim, LO))
               case (BND_COR, BND_INF, BND_MPI, BND_PER)
                  ! Do nothing
               case (BND_REF)
                  if (zndiff == 1) var(:,ledge,:) = 0.0
                  do ib=1,lnbcells
                     var(:,lnbcells+1-ib,:) = bndsign * var(:,ledge+ib,:)
                  enddo
               case (BND_OUT, BND_OUTD, BND_OUTH)
                  ledge = ledge + 1 ; lnbcells = lnbcells + 1
                  dvary = var(:,ledge+1,:)-var(:,ledge,:)
                  do ib=1,lnbcells
                     var(:,ib,:) = var(:, cg%js,:) - real(cg%js-ib)*dvary
                  enddo
               case default
                  write(msg,'(a,i3,3a,i3)') "[magboundaries:bnd_emf]: Boundary condition ",cg%bnd(ydim, LO)," not implemented for ",name, " in ", dir
                  if (master) call warn(msg)
            end select  ! (cg%bnd(ydim, LO))

            select case (cg%bnd(ydim, HI))
               case (BND_COR, BND_INF, BND_MPI, BND_PER)
                  ! Do nothing
               case (BND_REF)
                  if (zndiff == 1) var(:,redge,:) = 0.0
                  do ib=1,rnbcells
                     var(:,redge+ib,:) = bndsign * var(:,rrbase-ib,:)
                  enddo
               case (BND_OUT, BND_OUTD, BND_OUTH)
!                  dvary = var(:,rrbase,:)-var(:,rrbase-1,:) original
                  dvary = var(:,redge,:)-var(:,redge-1,:)
                  do ib=1,rnbcells
                     var(:,redge+ib,:) = var(:,redge,:) + real(ib)*dvary
                  enddo
               case default
                  write(msg,'(a,i3,3a,i3)') "[magboundaries:bnd_emf]: Boundary condition ",cg%bnd(ydim, HI)," not implemented for ",name, " in ", dir
                  if (master) call warn(msg)
            end select ! (cg%bnd(ydim, HI))

         case (zdim)

            if (any( [cg%bnd(zdim, LO) == BND_REF,  cg%bnd(zdim, HI) == BND_REF,  cg%bnd(zdim, LO) == BND_OUT,  cg%bnd(zdim, HI) == BND_OUT, &
               &      cg%bnd(zdim, LO) == BND_OUTD, cg%bnd(zdim, HI) == BND_OUTD, cg%bnd(zdim, LO) == BND_OUTH, cg%bnd(zdim, HI) == BND_OUTH] )) then
               select case (name)
                  case ("vzbx","vzby")
                     call compute_bnd_indxs(1, cg%nzb,ledge,redge,lnbcells,rnbcells,bndsign,zndiff,rrbase)
                  case ("vxbz","vybz","emfy","emfx")
                     call compute_bnd_indxs(2, cg%nzb,ledge,redge,lnbcells,rnbcells,bndsign,zndiff,rrbase)
                  case ("vxby","vybx","emfz")
                     call compute_bnd_indxs(3, cg%nzb,ledge,redge,lnbcells,rnbcells,bndsign,zndiff,rrbase)
               end select  ! (name)
            endif

            select case (cg%bnd(zdim, LO))
               case (BND_INF, BND_MPI, BND_PER)
                  ! Do nothing
               case (BND_REF)
                  if (zndiff == 1) var(:,:,ledge) = 0.0
                  do ib=1,lnbcells
                     var(:,:,lnbcells+1-ib) = bndsign * var(:,:,ledge+ib)
                  enddo
               case (BND_OUT, BND_OUTD, BND_OUTH)
                  ledge = ledge + 1 ; lnbcells = lnbcells + 1
                  dvarz = var(:,:,ledge+1)-var(:,:,ledge)
                  do ib=1,lnbcells
                     var(:,:,ib) = var(:,:, cg%ks) - real(cg%ks-ib)*dvarz
                  enddo
               case default
                  write(msg,'(a,i3,3a,i3)') "[magboundaries:bnd_emf]: Boundary condition ",cg%bnd(zdim, LO)," not implemented for ",name, " in ", dir
                  if (master) call warn(msg)
            end select  ! (cg%bnd(zdim, LO))

            select case (cg%bnd(zdim, HI))
               case (BND_INF, BND_MPI, BND_PER)
                  ! Do nothing
               case (BND_REF)
                  if (zndiff == 1) var(:,:,redge) = 0.0
                  do ib=1,rnbcells
                     var(:,:,redge+ib) = bndsign * var(:,:,rrbase-ib)
                  enddo
               case (BND_OUT, BND_OUTD, BND_OUTH)
!                  dvarz = var(:,:,rrbase)-var(:,:,rrbase-1) original
                  dvarz = var(:,:,redge)-var(:,:,redge-1)
                  do ib=1,rnbcells
                     var(:,:,redge+ib) = var(:,:,redge) + real(ib)*dvarz
                  enddo
               case default
                  write(msg,'(a,i3,3a,i3)') "[magboundaries:bnd_emf]: Boundary condition ",cg%bnd(zdim, HI)," not implemented for ",name, " in ", dir
                  if (master) call warn(msg)
            end select ! (cg%bnd(zdim, HI))

      end select ! (dim)

   end subroutine bnd_emf
!>
!! \brief Routine delivers common boundary cells indexes in cases of reflection or outflow boundary types
!<
   subroutine compute_bnd_indxs(bndcase,ndirb,ledge,redge,lnbcells,rnbcells,bndsign,zndiff,rrbase)

      use grid, only: cg

      implicit none

      integer, intent(in)  :: bndcase     !< 1 - v component compatible with direction; 2 - b component compatible with direction or emf component incompatible with direction; 3 - other cases
      integer, intent(in)  :: ndirb       !< cg%{nxb,nyb,nzb} depanding on the current direction
      integer, intent(out) :: ledge       !< index of the left edge of physical domain for emf
      integer, intent(out) :: redge       !< index of the right edge of physical domain for emf
      integer, intent(out) :: lnbcells    !< number of cells in a loop at left boundary
      integer, intent(out) :: rnbcells    !< number of cells in a loop at right boundary
      real,    intent(out) :: bndsign     !< 1. or -1. to change the sign or not
      integer, intent(out) :: zndiff      !< COMMENT ME
      integer, intent(out) :: rrbase      !< COMMENT ME

      select case (bndcase)
         case (1)
            ledge    = cg%nb
            lnbcells = cg%nb-1
            bndsign  = -1.
         case (2)
            ledge    = cg%nb+1
            lnbcells = cg%nb
            bndsign  = -1.
         case (3)
            ledge    = cg%nb
            lnbcells = cg%nb
            bndsign  = 1.
      end select  ! (name)

      zndiff   = ledge - lnbcells
      rnbcells = cg%nb - zndiff
      redge    = ndirb + ledge
      rrbase   = ndirb + lnbcells + 1  ! = redge + 1 - zndiff

   end subroutine compute_bnd_indxs

   subroutine all_mag_boundaries

      use arrays,    only: b
      use mpisetup,  only: has_dir, comm3d
      use mpi,       only: MPI_COMM_NULL
      use constants, only: xdim, zdim, MAG
      use grid,      only: cg

      implicit none

      integer :: dir
      real, dimension(:,:,:,:), pointer :: pb

      if (comm3d == MPI_COMM_NULL) then
         pb => b
         call cg%internal_boundaries(MAG, pa4d=pb)
      endif

      do dir = xdim, zdim
         if (has_dir(dir)) call bnd_b(dir)
      enddo

   end subroutine all_mag_boundaries

end module magboundaries
