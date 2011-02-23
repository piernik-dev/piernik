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

      use mpisetup,  only: ierr, req, comm3d, procxl, procxr, procyl, procyr, proczl, proczr, status, psize
      use constants, only: xdim, ydim, zdim
      use grid,      only: cg

      implicit none

      real, dimension(:,:,:,:) :: A

      if (psize(xdim) > 1) then

         call MPI_Isend  (A(1,1,1,1), 1, cg%MAG_YZ_LEFT_DOM,  procxl, 10, comm3d, req(1), ierr)
         call MPI_Isend  (A(1,1,1,1), 1, cg%MAG_YZ_RIGHT_DOM, procxr, 20, comm3d, req(3), ierr)
         call MPI_Irecv  (A(1,1,1,1), 1, cg%MAG_YZ_LEFT_BND,  procxl, 20, comm3d, req(2), ierr)
         call MPI_Irecv  (A(1,1,1,1), 1, cg%MAG_YZ_RIGHT_BND, procxr, 10, comm3d, req(4), ierr)

         call MPI_Waitall(4,req(:),status(:,:),ierr)
      endif

      if (psize(ydim) > 1) then
         call MPI_Isend  (A(1,1,1,1), 1, cg%MAG_XZ_LEFT_DOM,  procyl, 30, comm3d, req(1), ierr)
         call MPI_Isend  (A(1,1,1,1), 1, cg%MAG_XZ_RIGHT_DOM, procyr, 40, comm3d, req(3), ierr)
         call MPI_Irecv  (A(1,1,1,1), 1, cg%MAG_XZ_LEFT_BND,  procyl, 40, comm3d, req(2), ierr)
         call MPI_Irecv  (A(1,1,1,1), 1, cg%MAG_XZ_RIGHT_BND, procyr, 30, comm3d, req(4), ierr)

         call MPI_Waitall(4,req(:),status(:,:),ierr)
      endif

      if (psize(zdim) > 1) then
         call MPI_Isend  (A(1,1,1,1), 1, cg%MAG_XY_LEFT_DOM,  proczl, 50, comm3d, req(1), ierr)
         call MPI_Isend  (A(1,1,1,1), 1, cg%MAG_XY_RIGHT_DOM, proczr, 60, comm3d, req(3), ierr)
         call MPI_Irecv  (A(1,1,1,1), 1, cg%MAG_XY_LEFT_BND,  proczl, 60, comm3d, req(2), ierr)
         call MPI_Irecv  (A(1,1,1,1), 1, cg%MAG_XY_RIGHT_BND, proczr, 50, comm3d, req(4), ierr)

         call MPI_Waitall(4,req(:),status(:,:),ierr)
      endif

   end subroutine bnd_a

   subroutine bnd_b(dim)

      use arrays,        only: b
      use dataio_pub,    only: msg, warn
      use fluidindex,    only: ibx, iby, ibz
      use grid,          only: cg
      use mpi,           only: MPI_DOUBLE_PRECISION
      use mpisetup,      only: bnd_xl, bnd_xr, bnd_yl, bnd_yr, bnd_zl, bnd_zr, &
           &                   ierr, req, comm3d, procxl, procxr, procyl, procyr, proczl, proczr, status, &
           &                   psize, procxyl, procyxl, pcoords, comm
      use constants,     only: xdim, ydim, zdim
#ifdef SHEAR
      use shear,         only: eps,delj
#endif /* SHEAR */

      implicit none

      character(len=*)  :: dim
      integer           :: i, j
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

      select case (dim)
         case ("xdim")
#ifdef SHEAR
            allocate(send_left(3, cg%nb, cg%ny, cg%nz),send_right(3, cg%nb, cg%ny, cg%nz), &
                     recv_left(3, cg%nb, cg%ny, cg%nz),recv_right(3, cg%nb, cg%ny, cg%nz))

            send_left (:,:,:,:)  = b(:, cg%is:cg%isb,:,:)
            send_right(:,:,:,:)  = b(:, cg%ieb:cg%ie,:,:)

            if (bnd_xl == "she") then
!
! przesuwamy o calkowita liczbe komorek + periodyczny wb w kierunku y
!
               send_left (:,:, cg%js:cg%je,:)         = cshift(send_left (:,:, cg%js:cg%je,:),dim=3,shift= delj)
               send_left (:,:,1:cg%nb,:)                = send_left  (:,:, cg%jeb:cg%je,:)
               send_left (:,:, cg%je+1:cg%ny,:)   = send_left  (:,:, cg%js:cg%jsb,:)
!
! remapujemy  - interpolacja kwadratowa
!
               send_left (:,:,:,:)  = (1.+eps)*(1.-eps) * send_left (:,:,:,:) &
                                      -0.5*eps*(1.-eps) * cshift(send_left (:,:,:,:),shift=-1,dim=3) &
                                      +0.5*eps*(1.+eps) * cshift(send_left (:,:,:,:),shift=1,dim=3)
            endif ! (bnd_xl == "she")

            if (bnd_xr == "she") then
!
! przesuwamy o calkowita liczbe komorek + periodyczny wb w kierunku y
!
               send_right (:,:, cg%js:cg%je,:)        = cshift(send_right(:,:, cg%js:cg%je,:),dim=3,shift=-delj)
               send_right (:,:,1:cg%nb,:)               = send_right (:,:, cg%jeb:cg%je,:)
               send_right (:,:, cg%je+1:cg%ny,:)  = send_right (:,:, cg%js:cg%jsb,:)
!
! remapujemy - interpolacja kwadratowa
!
               send_right (:,:,:,:) = (1.+eps)*(1.-eps) * send_right (:,:,:,:) &
                                      -0.5*eps*(1.-eps) * cshift(send_right (:,:,:,:),shift=1,dim=3) &
                                      +0.5*eps*(1.+eps) * cshift(send_right (:,:,:,:),shift=-1,dim=3)
            endif ! (bnd_xr == "she")
!
! wysylamy na drugi brzeg
!
            call MPI_Isend   (send_left , 3*cg%ny*cg%nz*cg%nb, MPI_DOUBLE_PRECISION, procxl, 10, comm, req(1), ierr)
            call MPI_Isend   (send_right, 3*cg%ny*cg%nz*cg%nb, MPI_DOUBLE_PRECISION, procxr, 20, comm, req(3), ierr)
            call MPI_Irecv   (recv_left , 3*cg%ny*cg%nz*cg%nb, MPI_DOUBLE_PRECISION, procxl, 20, comm, req(2), ierr)
            call MPI_Irecv   (recv_right, 3*cg%ny*cg%nz*cg%nb, MPI_DOUBLE_PRECISION, procxr, 10, comm, req(4), ierr)

            call MPI_Waitall(4,req(:),status(:,:),ierr)

            b(:,1:cg%nb-1,:,:)               = recv_left(:,1:cg%nb-1,:,:)
            b(:, cg%ie+1+1:cg%nx,:,:)  = recv_right(:,1+1:cg%nb,:,:)

            if (allocated(send_left))  deallocate(send_left)
            if (allocated(send_right)) deallocate(send_right)
            if (allocated(recv_left))  deallocate(recv_left)
            if (allocated(recv_right)) deallocate(recv_right)

!===============================================================================
#else /* !SHEAR */

            if (psize(xdim) > 1) then

               call MPI_Isend  (b(1,1,1,1), 1, cg%MAG_YZ_LEFT_DOM,  procxl, 10, comm3d, req(1), ierr)
               call MPI_Isend  (b(1,1,1,1), 1, cg%MAG_YZ_RIGHT_DOM, procxr, 20, comm3d, req(3), ierr)
               call MPI_Irecv  (b(1,1,1,1), 1, cg%MAG_YZ_LEFT_BND,  procxl, 20, comm3d, req(2), ierr)
               call MPI_Irecv  (b(1,1,1,1), 1, cg%MAG_YZ_RIGHT_BND, procxr, 10, comm3d, req(4), ierr)

               call MPI_Waitall(4,req(:),status(:,:),ierr)

            endif
#endif /* !SHEAR */

         case ("ydim")
            if (psize(ydim) > 1) then

               call MPI_Isend  (b(1,1,1,1), 1, cg%MAG_XZ_LEFT_DOM,  procyl, 30, comm3d, req(1), ierr)
               call MPI_Isend  (b(1,1,1,1), 1, cg%MAG_XZ_RIGHT_DOM, procyr, 40, comm3d, req(3), ierr)
               call MPI_Irecv  (b(1,1,1,1), 1, cg%MAG_XZ_LEFT_BND,  procyl, 40, comm3d, req(2), ierr)
               call MPI_Irecv  (b(1,1,1,1), 1, cg%MAG_XZ_RIGHT_BND, procyr, 30, comm3d, req(4), ierr)

               call MPI_Waitall(4,req(:),status(:,:),ierr)
            endif

         case ("zdim")
            if (psize(zdim) > 1) then
               call MPI_Isend  (b(1,1,1,1), 1, cg%MAG_XY_LEFT_DOM,  proczl, 50, comm3d, req(1), ierr)
               call MPI_Isend  (b(1,1,1,1), 1, cg%MAG_XY_RIGHT_DOM, proczr, 60, comm3d, req(3), ierr)
               call MPI_Irecv  (b(1,1,1,1), 1, cg%MAG_XY_LEFT_BND,  proczl, 60, comm3d, req(2), ierr)
               call MPI_Irecv  (b(1,1,1,1), 1, cg%MAG_XY_RIGHT_BND, proczr, 50, comm3d, req(4), ierr)

               call MPI_Waitall(4,req(:),status(:,:),ierr)
            endif
      end select ! (dim)

! MPI + non-MPI corner-periodic boundary condition

      if (bnd_xl .eq. "cor") then
!   - lower to left
         if (pcoords(1) .eq. 0 .and. pcoords(2) .eq. 0) then
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

            call MPI_Isend   (send_left , 3*cg%nb*cg%ny*cg%nz, MPI_DOUBLE_PRECISION, procxyl, 70, comm, req(1), ierr)
            call MPI_Irecv   (recv_left , 3*cg%nx*cg%nb*cg%nz, MPI_DOUBLE_PRECISION, procxyl, 80, comm, req(2), ierr)

            call MPI_Waitall(2,req(:),status(:,:),ierr)

            do i=1, cg%nb
               do j=1, cg%ny
                  b(ibx,i,j,:) = -recv_left(iby,j, cg%is-i,:)
                  b(iby,i,j,:) =  recv_left(ibx,j, cg%is-i,:)
                  b(ibz,i,j,:) =  recv_left(ibz,j, cg%is-i,:)
               enddo
            enddo

            if (allocated(send_left))  deallocate(send_left)
            if (allocated(recv_left))  deallocate(recv_left)
         endif
      endif

      if (bnd_yl .eq. "cor") then
!   - left to lower
         if (pcoords(2) .eq. 0 .and. pcoords(1) .eq. 0 ) then
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

            if (allocated(send_left))  deallocate(send_left)
            if (allocated(recv_left))  deallocate(recv_left)
         endif
      endif

! Non-MPI boundary conditions
      if (frun) then
         bnd_xl_not_provided = any( [bnd_xl(1:3) == "cor", bnd_xl(1:3) == "inf", bnd_xl(1:3) == "mpi", bnd_xl(1:3) == "ref", bnd_xl(1:3) == "she"] )
         bnd_xr_not_provided = any( [bnd_xr(1:3) == "cor", bnd_xr(1:3) == "inf", bnd_xr(1:3) == "mpi", bnd_xr(1:3) == "ref", bnd_xr(1:3) == "she"] )
         bnd_yl_not_provided = any( [bnd_yl(1:3) == "cor", bnd_yl(1:3) == "inf", bnd_yl(1:3) == "mpi", bnd_yl(1:3) == "ref" ] )
         bnd_yr_not_provided = any( [bnd_yr(1:3) == "cor", bnd_yr(1:3) == "inf", bnd_yr(1:3) == "mpi", bnd_yr(1:3) == "ref" ] )
         bnd_zl_not_provided = any( [bnd_zl(1:3) == "inf", bnd_zl(1:3) == "ref", bnd_zl(1:3) == "mpi" ] )
         bnd_zr_not_provided = any( [bnd_zr(1:3) == "inf", bnd_zr(1:3) == "ref", bnd_zr(1:3) == "mpi" ] )
         frun = .false.
      endif

      if (dim=="xdim" .and. bnd_xl_not_provided .and. bnd_xr_not_provided) return  ! avoid triple case
      if (dim=="ydim" .and. bnd_yl_not_provided .and. bnd_yr_not_provided) return  ! avoid triple case
      if (dim=="zdim" .and. bnd_zl_not_provided .and. bnd_zr_not_provided) return  ! avoid triple case

      select case (dim)
         case ("xdim")

            select case (bnd_xl(1:3))
               case ("cor", "inf", "mpi", "ref", "she")
                  ! Do nothing
               case ("per")
                  b(:,1:cg%nb,:,:)              = b(:, cg%ieb:cg%ie,:,:)
               case ("out")
                  b(:,1,:,:) = b(:,2,:,:)
               case default
                  write(msg,'(4a)') "[magboundaries:bnd_b]: Boundary condition ",bnd_xl," not implemented in ",dim
                  call warn(msg)
            end select  ! (bnd_xl)

            select case (bnd_xr(1:3))
               case ("cor", "inf", "mpi", "ref", "she")
                  ! Do nothing
               case ("per")
                  b(:, cg%ie+1:cg%nx,:,:) = b(:, cg%is:cg%isb,:,:)
               case ("out")
                  b(:, cg%nx,:,:) = b(:, cg%nx-1,:,:)
               case default
                  write(msg,'(4a)') "[magboundaries:bnd_b]: Boundary condition ",bnd_xr," not implemented in ",dim
                  call warn(msg)
            end select  ! (bnd_xr)

         case ("ydim")

            select case (bnd_yl(1:3))
               case ("cor", "inf", "mpi", "ref")
                  ! Do nothing
               case ("per")
                  b(:,:,1:cg%nb,:)              = b(:,:, cg%jeb:cg%je,:)
               case ("out")
                  b(:,:,1,:) = b(:,:,2,:)
               case default
                  write(msg,'(4a)') "[magboundaries:bnd_b]: Boundary condition ",bnd_yl," not implemented in ",dim
                  call warn(msg)
            end select  ! (bnd_yl)

            select case (bnd_yr(1:3))
               case ("cor", "inf", "mpi", "ref")
                  ! Do nothing
               case ("per")
                  b(:,:, cg%je+1:cg%ny,:) = b(:,:, cg%js:cg%jsb,:)
               case ("out")
                  b(:,:, cg%ny,:) = b(:,:, cg%ny-1,:)
               case default
                  write(msg,'(4a)') "[magboundaries:bnd_b]: Boundary condition ",bnd_yr," not implemented in ",dim
                  call warn(msg)

            end select  ! (bnd_yr)

         case ("zdim")

            select case (bnd_zl(1:3))
               case ("mpi", "ref")
                  ! Do nothing
               case ("per")
                  b(:,:,:,1:cg%nb)              = b(:,:,:, cg%keb:cg%ke)
               case ("out")
                  b(:,:,:,1) = b(:,:,:,2)
               case default
                  write(msg,'(4a)') "[magboundaries:bnd_b]: Boundary condition ",bnd_zl," not implemented in ",dim
                  call warn(msg)
            end select  ! (bnd_zl)

            select case (bnd_zr(1:3))
               case ("mpi", "ref")
                  ! Do nothing
               case ("per")
                  b(:,:,:, cg%ke+1:cg%nz) = b(:,:,:, cg%ks:cg%ksb)
               case ("out")
                  b(:,:,:, cg%nz) = b(:,:,:, cg%nz-1)
               case default
                  write(msg,'(4a)') "[magboundaries:bnd_b]: Boundary condition ",bnd_zr," not implemented in ",dim
                  call warn(msg)
            end select  ! (bnd_zr)

      end select  ! (dim)

   end subroutine bnd_b

!=====================================================================================================

   subroutine bnd_emf(var, name, dim)

      use dataio_pub,    only: msg, warn
      use grid,          only: cg
      use mpisetup,      only: bnd_xl, bnd_xr, bnd_yl, bnd_yr, bnd_zl, bnd_zr

      implicit none

      real, dimension(:,:,:), intent(inout) :: var
      character(len=*), intent(in)          :: name, dim
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
         bnd_xl_not_provided = any( [bnd_xl(1:3) == "cor", bnd_xl(1:3) == "inf", bnd_xl(1:3) == "per", bnd_xl(1:3) == "mpi", bnd_xl(1:3) == "she"] )
         bnd_xr_not_provided = any( [bnd_xr(1:3) == "cor", bnd_xr(1:3) == "inf", bnd_xr(1:3) == "per", bnd_xr(1:3) == "mpi", bnd_xr(1:3) == "she"] )
         bnd_yl_not_provided = any( [bnd_yl(1:3) == "cor", bnd_yl(1:3) == "inf", bnd_yl(1:3) == "per", bnd_yl(1:3) == "mpi" ] )
         bnd_yr_not_provided = any( [bnd_yr(1:3) == "cor", bnd_yr(1:3) == "inf", bnd_yr(1:3) == "per", bnd_yr(1:3) == "mpi" ] )
         bnd_zl_not_provided = any( [bnd_zl(1:3) == "inf", bnd_zl(1:3) == "per", bnd_zl(1:3) == "mpi" ] )
         bnd_zr_not_provided = any( [bnd_zr(1:3) == "inf", bnd_zr(1:3) == "per", bnd_zr(1:3) == "mpi" ] )
         frun = .false.
      endif

      if (dim=="xdim" .and. bnd_xl_not_provided .and. bnd_xr_not_provided) return  ! avoid triple case
      if (dim=="ydim" .and. bnd_yl_not_provided .and. bnd_yr_not_provided) return  ! avoid triple case
      if (dim=="zdim" .and. bnd_zl_not_provided .and. bnd_zr_not_provided) return  ! avoid triple case

      bndsign = huge(1.0); ledge=huge(1); redge=huge(1); lnbcells=huge(1); rnbcells=huge(1); zndiff=huge(1); rrbase=huge(1)
      ! the code below should not use these values and the compiler should not complain on possible use of uninitialized variables.

      select case (dim)
         case ("xdim")

            if (any( [bnd_xl(1:3) == "ref", bnd_xr(1:3) == "ref", bnd_xl(1:3) == "out", bnd_xr(1:3) == "out"] )) then
               select case (name)
                  case ("vxby","vxbz")
                     call compute_bnd_indxs(1, cg%nxb,ledge,redge,lnbcells,rnbcells,bndsign,zndiff,rrbase)
                  case ("vybx","vzbx","emfy","emfz")
                     call compute_bnd_indxs(2, cg%nxb,ledge,redge,lnbcells,rnbcells,bndsign,zndiff,rrbase)
                  case ("vybz","vzby","emfx")
                     call compute_bnd_indxs(3, cg%nxb,ledge,redge,lnbcells,rnbcells,bndsign,zndiff,rrbase)
               end select  ! (name)
            endif

            select case (bnd_xl(1:3))
               case ("cor", "inf", "mpi", "per", "she")
                  ! Do nothing
               case ("ref")
                  if (zndiff == 1) var(ledge,:,:) = 0.0
                  do ib=1,lnbcells
                     var(lnbcells+1-ib,:,:) = bndsign * var(ledge+ib,:,:)
                  enddo
               case ("out")
                  ledge = ledge + 1 ; lnbcells = lnbcells + 1
                  dvarx = var(ledge+1,:,:)-var(ledge,:,:)
                  do ib=1,lnbcells
                     var(ib,:,:) = var(cg%is,:,:) - real(cg%is-ib)*dvarx
                  enddo
               case default
                  write(msg,'(6a)') "[magboundaries:bnd_emf]: Boundary condition ",bnd_xl," not implemented for ",name, " in ", dim
                  call warn(msg)
            end select  ! (bnd_xl)

            select case (bnd_xr(1:3))
               case ("cor", "inf", "mpi", "per", "she")
                  ! Do nothing
               case ("ref")
                  if (zndiff == 1) var(redge,:,:) = 0.0
                  do ib=1,rnbcells
                     var(redge+ib,:,:) = bndsign * var(rrbase-ib,:,:)
                  enddo
               case ("out")
!                  dvarx = var(rrbase,:,:)-var(rrbase-1,:,:) original
                  dvarx = var(redge,:,:)-var(redge-1,:,:)
                  do ib=1,rnbcells
                     var(redge+ib,:,:) = var(redge,:,:) + real(ib)*dvarx
                  enddo
               case default
                  write(msg,'(6a)') "[magboundaries:bnd_emf]: Boundary condition ",bnd_xr," not implemented for ",name, " in ", dim
                  call warn(msg)
            end select ! (bnd_xr)

         case ("ydim")

            if (any( [bnd_yl(1:3) == "ref", bnd_yr(1:3) == "ref", bnd_yl(1:3) == "out", bnd_yr(1:3) == "out"] )) then
               select case (name)
                  case ("vybz","vybx")
                     call compute_bnd_indxs(1, cg%nyb,ledge,redge,lnbcells,rnbcells,bndsign,zndiff,rrbase)
                  case ("vzby","vxby","emfz","emfx")
                     call compute_bnd_indxs(2, cg%nyb,ledge,redge,lnbcells,rnbcells,bndsign,zndiff,rrbase)
                  case ("vzbx","vxbz","emfy")
                     call compute_bnd_indxs(3, cg%nyb,ledge,redge,lnbcells,rnbcells,bndsign,zndiff,rrbase)
               end select  ! (name)
            endif

            select case (bnd_yl(1:3))
               case ("cor", "inf", "mpi", "per")
                  ! Do nothing
               case ("ref")
                  if (zndiff == 1) var(:,ledge,:) = 0.0
                  do ib=1,lnbcells
                     var(:,lnbcells+1-ib,:) = bndsign * var(:,ledge+ib,:)
                  enddo
               case ("out")
                  ledge = ledge + 1 ; lnbcells = lnbcells + 1
                  dvary = var(:,ledge+1,:)-var(:,ledge,:)
                  do ib=1,lnbcells
                     var(:,ib,:) = var(:, cg%js,:) - real(cg%js-ib)*dvary
                  enddo
               case default
                  write(msg,'(6a)') "[magboundaries:bnd_emf]: Boundary condition ",bnd_yl," not implemented for ",name, " in ", dim
                  call warn(msg)
            end select  ! (bnd_yl)

            select case (bnd_yr(1:3))
               case ("cor", "inf", "mpi", "per")
                  ! Do nothing
               case ("ref")
                  if (zndiff == 1) var(:,redge,:) = 0.0
                  do ib=1,rnbcells
                     var(:,redge+ib,:) = bndsign * var(:,rrbase-ib,:)
                  enddo
               case ("out")
!                  dvary = var(:,rrbase,:)-var(:,rrbase-1,:) original
                  dvary = var(:,redge,:)-var(:,redge-1,:)
                  do ib=1,rnbcells
                     var(:,redge+ib,:) = var(:,redge,:) + real(ib)*dvary
                  enddo
               case default
                  write(msg,'(6a)') "[magboundaries:bnd_emf]: Boundary condition ",bnd_yr," not implemented for ",name, " in ", dim
                  call warn(msg)
            end select ! (bnd_yr)

         case ("zdim")

            if (any( [bnd_zl(1:3) == "ref", bnd_zr(1:3) == "ref", bnd_zl(1:3) == "out", bnd_zr(1:3) == "out"] )) then
               select case (name)
                  case ("vzbx","vzby")
                     call compute_bnd_indxs(1, cg%nzb,ledge,redge,lnbcells,rnbcells,bndsign,zndiff,rrbase)
                  case ("vxbz","vybz","emfy","emfx")
                     call compute_bnd_indxs(2, cg%nzb,ledge,redge,lnbcells,rnbcells,bndsign,zndiff,rrbase)
                  case ("vxby","vybx","emfz")
                     call compute_bnd_indxs(3, cg%nzb,ledge,redge,lnbcells,rnbcells,bndsign,zndiff,rrbase)
               end select  ! (name)
            endif

            select case (bnd_zl(1:3))
               case ("inf", "mpi", "per")
                  ! Do nothing
               case ("ref")
                  if (zndiff == 1) var(:,:,ledge) = 0.0
                  do ib=1,lnbcells
                     var(:,:,lnbcells+1-ib) = bndsign * var(:,:,ledge+ib)
                  enddo
               case ("out")
                  ledge = ledge + 1 ; lnbcells = lnbcells + 1
                  dvarz = var(:,:,ledge+1)-var(:,:,ledge)
                  do ib=1,lnbcells
                     var(:,:,ib) = var(:,:, cg%ks) - real(cg%ks-ib)*dvarz
                  enddo
               case default
                  write(msg,'(6a)') "[magboundaries:bnd_emf]: Boundary condition ",bnd_zl," not implemented for ",name, " in ", dim
                  call warn(msg)
            end select  ! (bnd_zl)

            select case (bnd_zr(1:3))
               case ("inf", "mpi", "per")
                  ! Do nothing
               case ("ref")
                  if (zndiff == 1) var(:,:,redge) = 0.0
                  do ib=1,rnbcells
                     var(:,:,redge+ib) = bndsign * var(:,:,rrbase-ib)
                  enddo
               case ("out")
!                  dvarz = var(:,:,rrbase)-var(:,:,rrbase-1) original
                  dvarz = var(:,:,redge)-var(:,:,redge-1)
                  do ib=1,rnbcells
                     var(:,:,redge+ib) = var(:,:,redge) + real(ib)*dvarz
                  enddo
               case default
                  write(msg,'(6a)') "[magboundaries:bnd_emf]: Boundary condition ",bnd_zr," not implemented for ",name, " in ", dim
                  call warn(msg)
            end select ! (bnd_zr)

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

      use mpisetup,  only: has_dir
      use constants, only: xdim, ydim, zdim

      implicit none

      if (has_dir(xdim)) call bnd_b("xdim")
      if (has_dir(ydim)) call bnd_b("ydim")
      if (has_dir(zdim)) call bnd_b("zdim")

   end subroutine all_mag_boundaries

end module magboundaries
