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

      use constants,  only: MAG, xdim, zdim, LO, HI, BND, BLK, I_ONE, I_FOUR, I_FIVE, I_TEN
      use dataio_pub, only: die
      use domain,     only: is_mpi_noncart, is_multicg, dom
      use grid,       only: leaves
      use grid_cont,  only: grid_container
      use mpi,        only: MPI_COMM_NULL
      use mpisetup,   only: ierr, req, status, have_mpi
      use types,      only: cdd

      implicit none

      real, dimension(:,:,:,:) :: A  !< vector potential of magnetic field
      integer(kind=4)          :: i, itag, jtag
      type(grid_container), pointer :: cg

      cg => leaves%first%cg
      if (is_multicg) call die("[magboundaries:bnd_a] multiple grid pieces per procesor not implemented yet") !nontrivial MPI_Waitall

      if (have_mpi .and. is_mpi_noncart) call die("[magboundaries:bnd_a] is_mpi_noncart is not implemented") !procn, psize
      if (cdd%comm3d == MPI_COMM_NULL) call die("[magboundaries:bnd_a] cdd%comm3d == MPI_COMM_NULL")

      do i = xdim, zdim
         if (cdd%psize(i) > 1) then

            jtag = I_TEN*i
            itag = jtag - I_FIVE
            call MPI_Isend(A(1,1,1,1), I_ONE, cg%mbc(MAG, i, LO, BLK, dom%nb), cdd%procn(i,LO), itag, cdd%comm3d, req(1), ierr)
            call MPI_Isend(A(1,1,1,1), I_ONE, cg%mbc(MAG, i, HI, BLK, dom%nb), cdd%procn(i,HI), jtag, cdd%comm3d, req(3), ierr)
            call MPI_Irecv(A(1,1,1,1), I_ONE, cg%mbc(MAG, i, LO, BND, dom%nb), cdd%procn(i,LO), jtag, cdd%comm3d, req(2), ierr)
            call MPI_Irecv(A(1,1,1,1), I_ONE, cg%mbc(MAG, i, HI, BND, dom%nb), cdd%procn(i,HI), itag, cdd%comm3d, req(4), ierr)

            call MPI_Waitall(I_FOUR,req(:),status(:,:),ierr)
         endif
      enddo

   end subroutine bnd_a

   subroutine bnd_b(dir, cg)

      use constants,  only: MAG, ndims, xdim, ydim, zdim, LO, HI, BND, BLK, I_ONE, I_TWO, I_FOUR, &
           &                BND_MPI, BND_PER, BND_REF, BND_OUT, BND_OUTD, BND_OUTH, BND_COR, BND_SHE
      use dataio_pub, only: msg, warn, die
      use domain,     only: is_mpi_noncart, is_multicg, dom
      use grid_cont,  only: grid_container
      use mpi,        only: MPI_DOUBLE_PRECISION, MPI_COMM_NULL
      use mpisetup,   only: ierr, req, status, comm, master, have_mpi
      use types,      only: cdd
#ifdef SHEAR
      use constants,  only: half, one
      use shear,      only: eps, delj
#endif /* SHEAR */

      implicit none

      integer(kind=4), intent(in) :: dir
      type(grid_container), pointer, intent(inout) :: cg

      integer(kind=4), parameter :: tag1 = 10, tag2 = 20
      integer(kind=4), parameter :: tag7 = 70, tag8 = 80
      integer           :: i, j
      integer(kind=4) :: itag, jtag
      real, allocatable :: send_left(:,:,:,:),recv_left(:,:,:,:)
#ifdef SHEAR
      real, allocatable :: send_right(:,:,:,:),recv_right(:,:,:,:)
#endif /* SHEAR */
      logical, save                         :: frun = .true.
      logical, dimension(ndims,LO:HI), save :: bnd_not_provided = .false.

! MPI block comunication
      if (cdd%comm3d /= MPI_COMM_NULL) then
         if (is_multicg) call die("[magboundaries:bnd_b] multiple grid pieces per procesor not implemented yet") !nontrivial MPI_Waitall

         if (have_mpi .and. is_mpi_noncart) call die("[magboundaries:bnd_b] is_mpi_noncart is not implemented") !procn, procxyl, procyxl, psize, pcoords
#ifdef SHEAR
         if (dir == xdim) then
            allocate(send_left(3, dom%nb, cg%ny, cg%nz),send_right(3, dom%nb, cg%ny, cg%nz), &
                 &   recv_left(3, dom%nb, cg%ny, cg%nz),recv_right(3, dom%nb, cg%ny, cg%nz))

            send_left (:,:,:,:)  = cg%b(:, cg%is:cg%isb,:,:)
            send_right(:,:,:,:)  = cg%b(:, cg%ieb:cg%ie,:,:)

            if (cg%bnd(xdim, LO) == BND_SHE) then
!
! przesuwamy o calkowita liczbe komorek + periodyczny wb w kierunku y
!
               send_left (:,:,  cg%js:cg%je,:) = cshift(send_left (:,:, cg%js :cg%je, :),dim=3,shift= delj)
               send_left (:,:,      1:dom%nb,:) =        send_left (:,:, cg%jeb:cg%je, :)
               send_left (:,:,cg%je+1:cg%ny,:) =        send_left (:,:, cg%js :cg%jsb,:)
!
! remapujemy  - interpolacja kwadratowa
!
               send_left (:,:,:,:)  = (one+eps)*(one-eps) * send_left (:,:,:,:) &
                    &                 -half*eps*(one-eps) * cshift(send_left (:,:,:,:),shift=-1,dim=3) &
                    &                 +half*eps*(one+eps) * cshift(send_left (:,:,:,:),shift= 1,dim=3)
            endif ! (cg%bnd(xdim, LO) == BND_SHE)

            if (cg%bnd(xdim, HI) == BND_SHE) then
!
! przesuwamy o calkowita liczbe komorek + periodyczny wb w kierunku y
!
               send_right (:,:,  cg%js:cg%je,:) = cshift(send_right(:,:, cg%js :cg%je, :),dim=3,shift=-delj)
               send_right (:,:,      1:dom%nb,:) =        send_right(:,:, cg%jeb:cg%je, :)
               send_right (:,:,cg%je+1:cg%ny,:) =        send_right(:,:, cg%js :cg%jsb,:)
!
! remapujemy - interpolacja kwadratowa
!
               send_right (:,:,:,:) = (one+eps)*(one-eps) * send_right (:,:,:,:) &
                    &                 -half*eps*(one-eps) * cshift(send_right (:,:,:,:),shift= 1,dim=3) &
                    &                 +half*eps*(one+eps) * cshift(send_right (:,:,:,:),shift=-1,dim=3)
            endif ! (cg%bnd(xdim, HI) == BND_SHE)
!
! wysylamy na drugi brzeg
!
            call MPI_Isend(send_left , 3*cg%ny*cg%nz*dom%nb, MPI_DOUBLE_PRECISION, cdd%procn(dir,LO), tag1, comm, req(1), ierr)
            call MPI_Isend(send_right, 3*cg%ny*cg%nz*dom%nb, MPI_DOUBLE_PRECISION, cdd%procn(dir,HI), tag2, comm, req(3), ierr)
            call MPI_Irecv(recv_left , 3*cg%ny*cg%nz*dom%nb, MPI_DOUBLE_PRECISION, cdd%procn(dir,LO), tag2, comm, req(2), ierr)
            call MPI_Irecv(recv_right, 3*cg%ny*cg%nz*dom%nb, MPI_DOUBLE_PRECISION, cdd%procn(dir,HI), tag1, comm, req(4), ierr)

            call MPI_Waitall(I_FOUR,req(:),status(:,:),ierr)

            cg%b(:,        1:dom%nb-1,:,:) = recv_left (:,  1:dom%nb-1,:,:)
            cg%b(:,cg%ie+1+1:cg%n_(xdim),  :,:) = recv_right(:,1+1:dom%nb,  :,:)

            if (allocated(send_left))  deallocate(send_left)
            if (allocated(send_right)) deallocate(send_right)
            if (allocated(recv_left))  deallocate(recv_left)
            if (allocated(recv_right)) deallocate(recv_right)

!===============================================================================
         else
#endif /* SHEAR */

            if (cdd%psize(dir) > 1) then

               jtag = tag2*dir
               itag = jtag - tag1
               call MPI_Isend(cg%b(1,1,1,1), I_ONE, cg%mbc(MAG, dir, LO, BLK, dom%nb), cdd%procn(dir,LO), itag, cdd%comm3d, req(1), ierr)
               call MPI_Isend(cg%b(1,1,1,1), I_ONE, cg%mbc(MAG, dir, HI, BLK, dom%nb), cdd%procn(dir,HI), jtag, cdd%comm3d, req(3), ierr)
               call MPI_Irecv(cg%b(1,1,1,1), I_ONE, cg%mbc(MAG, dir, LO, BND, dom%nb), cdd%procn(dir,LO), jtag, cdd%comm3d, req(2), ierr)
               call MPI_Irecv(cg%b(1,1,1,1), I_ONE, cg%mbc(MAG, dir, HI, BND, dom%nb), cdd%procn(dir,HI), itag, cdd%comm3d, req(4), ierr)

               call MPI_Waitall(I_FOUR,req(:),status(:,:),ierr)
            endif

#ifdef SHEAR
         endif
#endif /* SHEAR */

! MPI + non-MPI corner-periodic boundary condition

         if (cg%bnd(xdim, LO) == BND_COR) then
!   - lower to left
            if (cdd%pcoords(xdim) == 0 .and. cdd%pcoords(ydim) == 0) then
               do i=1, dom%nb
                  do j=cg%js, cg%n_(ydim)
                     cg%b(xdim,i,j,:) = -cg%b(ydim,j,cg%isb+1-i,:)
                     cg%b(ydim,i,j,:) =  cg%b(xdim,j,cg%isb+1-i,:)
                     cg%b(zdim,i,j,:) =  cg%b(zdim,j,cg%isb+1-i,:)
                  enddo
               enddo
            endif

            if (cdd%procxyl > 0) then
               allocate(send_left(3, dom%nb, cg%n_(ydim), cg%n_(zdim)), recv_left(3, cg%n_(xdim), dom%nb, cg%n_(zdim)))

               send_left(:,:,:,:) = cg%b(:, cg%is:cg%isb,:,:)

               call MPI_Isend(send_left, 3*dom%nb*cg%n_(ydim)*cg%n_(zdim), MPI_DOUBLE_PRECISION, cdd%procxyl, tag7, comm, req(1), ierr)
               call MPI_Irecv(recv_left, 3*cg%n_(xdim)*dom%nb*cg%n_(zdim), MPI_DOUBLE_PRECISION, cdd%procxyl, tag8, comm, req(2), ierr)

               call MPI_Waitall(I_TWO,req(:),status(:,:),ierr)

               do i=1, dom%nb
                  do j=1, cg%n_(ydim)
                     cg%b(xdim,i,j,:) = -recv_left(ydim,j, cg%is-i,:)
                     cg%b(ydim,i,j,:) =  recv_left(xdim,j, cg%is-i,:)
                     cg%b(zdim,i,j,:) =  recv_left(zdim,j, cg%is-i,:)
                  enddo
               enddo

               if (allocated(send_left)) deallocate(send_left)
               if (allocated(recv_left)) deallocate(recv_left)
            endif
         endif

         if (cg%bnd(ydim, LO) == BND_COR) then
!   - left to lower
            if (cdd%pcoords(ydim) == 0 .and. cdd%pcoords(xdim) == 0 ) then
               do j=1, dom%nb
                  do i=cg%is, cg%n_(xdim)
                     cg%b(xdim,i,j,:) =  cg%b(ydim,cg%isb+1-j,i,:)
                     cg%b(ydim,i,j,:) = -cg%b(xdim,cg%isb+1-j,i,:)
                     cg%b(zdim,i,j,:) =  cg%b(zdim,cg%isb+1-j,i,:)
                  enddo
               enddo
!   - interior to corner
               do j=1, dom%nb
                  do i=1, dom%nb
                     cg%b(xdim,i,j,:) =  -cg%b(xdim,cg%isb+1-i,cg%jsb+1-j,:)
                     cg%b(ydim,i,j,:) =  -cg%b(ydim,cg%isb+1-i,cg%jsb+1-j,:)
                     cg%b(zdim,i,j,:) =   cg%b(zdim,cg%isb+1-i,cg%jsb+1-j,:)
                  enddo
               enddo
            endif

            if (cdd%procyxl > 0) then
               allocate(send_left(3, cg%n_(xdim), dom%nb, cg%n_(zdim)), recv_left(3, dom%nb, cg%n_(ydim), cg%n_(zdim)))

               send_left(:,:,:,:) = cg%b(:,:, cg%js:cg%jsb,:)

               call MPI_Isend   (send_left , 3*cg%n_(xdim)*dom%nb*cg%n_(zdim), MPI_DOUBLE_PRECISION, cdd%procyxl, tag8, comm, req(1), ierr)
               call MPI_Irecv   (recv_left , 3*dom%nb*cg%n_(ydim)*cg%n_(zdim), MPI_DOUBLE_PRECISION, cdd%procyxl, tag7, comm, req(2), ierr)

               call MPI_Waitall(I_TWO,req(:),status(:,:),ierr)

               do j=1, dom%nb
                  do i=1, cg%n_(xdim)
                     cg%b(xdim,i,j,:) =  recv_left(ydim, cg%js-j,i,:)
                     cg%b(ydim,i,j,:) = -recv_left(xdim, cg%js-j,i,:)
                     cg%b(zdim,i,j,:) =  recv_left(zdim, cg%js-j,i,:)
                  enddo
               enddo

               if (allocated(send_left)) deallocate(send_left)
               if (allocated(recv_left)) deallocate(recv_left)
            endif
         endif

      endif

! Non-MPI boundary conditions
      if (frun) then
         bnd_not_provided(:,         :) = (cg%bnd(:,:) == BND_REF)       .or. (cg%bnd(:,         :) == BND_MPI) !! what about BND_PER?
         bnd_not_provided(xdim:ydim, :) = bnd_not_provided(xdim:ydim, :) .or. (cg%bnd(xdim:ydim, :) == BND_COR)
         bnd_not_provided(xdim,      :) = bnd_not_provided(xdim,      :) .or. (cg%bnd(xdim,      :) == BND_SHE)
      endif

      if (bnd_not_provided(dir,LO) .and. bnd_not_provided(dir,HI)) return  ! avoid triple case

      select case (dir)
         case (xdim)

            select case (cg%bnd(xdim, LO))
               case (BND_COR, BND_MPI, BND_REF, BND_SHE)
                  ! Do nothing
               case (BND_PER)
                  if (cdd%comm3d /= MPI_COMM_NULL) cg%b(:,1:dom%nb,:,:)              = cg%b(:, cg%ieb:cg%ie,:,:)
               case (BND_OUT, BND_OUTD, BND_OUTH)
                  cg%b(:,1,:,:) = cg%b(:,2,:,:)
               case default
                  write(msg,'(2(a,i3))') "[magboundaries:bnd_b]: Boundary condition ",cg%bnd(xdim, LO)," not implemented in ",dir
                  if (master) call warn(msg)
            end select  ! (cg%bnd(xdim, LO))

            select case (cg%bnd(xdim, HI))
               case (BND_COR, BND_MPI, BND_REF, BND_SHE)
                  ! Do nothing
               case (BND_PER)
                  if (cdd%comm3d /= MPI_COMM_NULL) cg%b(:, cg%ie+1:cg%n_(xdim),:,:) = cg%b(:, cg%is:cg%isb,:,:)
               case (BND_OUT, BND_OUTD, BND_OUTH)
                  cg%b(:, cg%n_(xdim),:,:) = cg%b(:, cg%n_(xdim)-1,:,:)
               case default
                  write(msg,'(2(a,i3))') "[magboundaries:bnd_b]: Boundary condition ",cg%bnd(xdim, HI)," not implemented in ",dir
                  if (master) call warn(msg)
            end select  ! (cg%bnd(xdim, HI))

         case (ydim)

            select case (cg%bnd(ydim, LO))
               case (BND_COR, BND_MPI, BND_REF)
                  ! Do nothing
               case (BND_PER)
                  if (cdd%comm3d /= MPI_COMM_NULL) cg%b(:,:,1:dom%nb,:)              = cg%b(:,:, cg%jeb:cg%je,:)
               case (BND_OUT, BND_OUTD, BND_OUTH)
                  cg%b(:,:,1,:) = cg%b(:,:,2,:)
               case default
                  write(msg,'(2(a,i3))') "[magboundaries:bnd_b]: Boundary condition ",cg%bnd(ydim, LO)," not implemented in ",dir
                  if (master) call warn(msg)
            end select  ! (cg%bnd(ydim, LO))

            select case (cg%bnd(ydim, HI))
               case (BND_COR, BND_MPI, BND_REF)
                  ! Do nothing
               case (BND_PER)
                  if (cdd%comm3d /= MPI_COMM_NULL) cg%b(:,:, cg%je+1:cg%n_(ydim),:) = cg%b(:,:, cg%js:cg%jsb,:)
               case (BND_OUT, BND_OUTD, BND_OUTH)
                  cg%b(:,:, cg%n_(ydim),:) = cg%b(:,:, cg%n_(ydim)-1,:)
               case default
                  write(msg,'(2(a,i3))') "[magboundaries:bnd_b]: Boundary condition ",cg%bnd(ydim, HI)," not implemented in ",dir
                  if (master) call warn(msg)

            end select  ! (cg%bnd(ydim, HI))

         case (zdim)

            select case (cg%bnd(zdim, LO))
               case (BND_MPI, BND_REF)
                  ! Do nothing
               case (BND_PER)
                  if (cdd%comm3d /= MPI_COMM_NULL) cg%b(:,:,:,1:dom%nb)              = cg%b(:,:,:, cg%keb:cg%ke)
               case (BND_OUT, BND_OUTD, BND_OUTH)
                  cg%b(:,:,:,1) = cg%b(:,:,:,2)
               case default
                  write(msg,'(2(a,i3))') "[magboundaries:bnd_b]: Boundary condition ",cg%bnd(zdim, LO)," not implemented in ",dir
                  if (master) call warn(msg)
            end select  ! (cg%bnd(zdim, LO))

            select case (cg%bnd(zdim, HI))
               case (BND_MPI, BND_REF)
                  ! Do nothing
               case (BND_PER)
                  if (cdd%comm3d /= MPI_COMM_NULL) cg%b(:,:,:, cg%ke+1:cg%n_(zdim)) = cg%b(:,:,:, cg%ks:cg%ksb)
               case (BND_OUT, BND_OUTD, BND_OUTH)
                  cg%b(:,:,:, cg%n_(zdim)) = cg%b(:,:,:, cg%n_(zdim)-1)
               case default
                  write(msg,'(2(a,i3))') "[magboundaries:bnd_b]: Boundary condition ",cg%bnd(zdim, HI)," not implemented in ",dir
                  if (master) call warn(msg)
            end select  ! (cg%bnd(zdim, HI))

      end select  ! (dim)

   end subroutine bnd_b

!=====================================================================================================

   subroutine bnd_emf(var, emfdir, dir, cg)

      use constants,  only: ndims, xdim, ydim, zdim, LO, HI, I_ZERO, I_ONE, &
                            BND_MPI, BND_PER, BND_REF, BND_OUT, BND_OUTD, BND_OUTH, BND_COR, BND_SHE
      use dataio_pub, only: msg, warn
      use grid_cont,  only: grid_container
      use mpisetup,   only: master

      implicit none

      real, dimension(:,:,:),        intent(inout) :: var
      integer(kind=4),               intent(in)    :: emfdir
      integer(kind=4),               intent(in)    :: dir
      type(grid_container), pointer, intent(inout) :: cg

#ifndef ZERO_BND_EMF
      real, dimension(:,:,:), allocatable          :: dvar
#endif /* !ZERO_BND_EMF */
      real                                         :: bndsign
      logical,                         save        :: frun = .true.
      logical, dimension(ndims,LO:HI), save        :: bnd_not_provided = .false.
      logical                                      :: zndiff
      integer(kind=4), dimension(ndims,LO:HI)      :: l, r
      integer(kind=4), dimension(LO:HI)            :: sbase, edge, nbcells, sidebase
      integer(kind=4)                              :: ssign, side, ib

      if (frun) then
         bnd_not_provided(:, :)         = (cg%bnd(:,:) == BND_PER)       .or. (cg%bnd(:,         :) == BND_MPI)
         bnd_not_provided(xdim:ydim, :) = bnd_not_provided(xdim:ydim, :) .or. (cg%bnd(xdim:ydim, :) == BND_COR)
         bnd_not_provided(xdim, :)      = bnd_not_provided(xdim,      :) .or. (cg%bnd(xdim,      :) == BND_SHE)
         frun = .false.
      endif

      if (bnd_not_provided(dir,LO) .and. bnd_not_provided(dir,HI)) return  ! avoid triple case

      bndsign = huge(1.0); edge=huge(I_ONE); nbcells=huge(I_ONE); zndiff=.false.; sidebase=huge(I_ONE)
      ! the code below should not use these values and the compiler should not complain on possible use of uninitialized variables.

      if (any( [cg%bnd(dir, LO) == BND_REF,  cg%bnd(dir, HI) == BND_REF,  cg%bnd(dir, LO) == BND_OUT,  cg%bnd(dir, HI) == BND_OUT, &
         &      cg%bnd(dir, LO) == BND_OUTD, cg%bnd(dir, HI) == BND_OUTD, cg%bnd(dir, LO) == BND_OUTH, cg%bnd(dir, HI) == BND_OUTH] )) then
         call compute_bnd_indxs(emfdir, cg%n_b(dir),edge,nbcells,sidebase,bndsign,zndiff)
      endif

      l = reshape([lbound(var, kind=4),ubound(var, kind=4)],shape=[ndims,HI]) ; r = l

      do side = LO, HI
         select case (cg%bnd(dir, side))
            case (BND_MPI, BND_PER)
               ! Do nothing
            case (BND_COR)
               if (dir == zdim) then
                  write(msg,'(a,i3,a,i1,a,i3)') "[magboundaries:bnd_emf]: Boundary condition ",cg%bnd(dir, side)," not implemented for ",emfdir, " in ", dir
                  if (master) call warn(msg)
               endif
            case (BND_SHE)
               if (dir /= xdim) then
                  write(msg,'(a,i3,a,i1,a,i3)') "[magboundaries:bnd_emf]: Boundary condition ",cg%bnd(dir, side)," not implemented for ",emfdir, " in ", dir
                  if (master) call warn(msg)
               endif
            case (BND_REF)
               sbase(:)  = [nbcells(LO)+I_ONE, edge(HI)] ; ssign = int(2*side-3, kind=4)
               if (zndiff) then
                  l(dir,:) = edge(side)
                  var(l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = 0.0
               endif
               do ib=1,nbcells(side)
                  l(dir,:) = sbase(side)+ssign*ib
                  r(dir,:) = sidebase(side)-ssign*ib
                  var(l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = bndsign * var(r(xdim,LO):r(xdim,HI),r(ydim,LO):r(ydim,HI),r(zdim,LO):r(zdim,HI))
               enddo
            case (BND_OUT, BND_OUTD, BND_OUTH)
               sbase(:)  = [I_ZERO, edge(HI)]
#ifdef ZERO_BND_EMF
               l(dir,LO) = sbase(side)+1 ; l(dir,HI) = sbase(side)+nbcells(side)
               var(l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = 0.0
#else /* !ZERO_BND_EMF */
               l(dir,:) = 1 ; allocate(dvar(l(xdim,HI) ,l(ydim,HI), l(zdim,HI)))
               edge(side) = edge(side) + HI - side ; nbcells(side) = nbcells(side) + HI - side
!               l(dir,:) = sidebase(side)+HI-side ; r(dir,:) = l(dir,:)-1 original
               l(dir,:) = edge(side)+HI-side ; r(dir,:) = l(dir,:) - I_ONE
               dvar = var(l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) - var(r(xdim,LO):r(xdim,HI),r(ydim,LO):r(ydim,HI),r(zdim,LO):r(zdim,HI))
               r(dir,:) = edge(side)
               do ib=1,nbcells(side)
                  l(dir,:) = sbase(side) + ib
                  var(l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = var(r(xdim,LO):r(xdim,HI),r(ydim,LO):r(ydim,HI),r(zdim,LO):r(zdim,HI)) + real(ib+sbase(side)-edge(side))*dvar
               enddo
               deallocate(dvar)
#endif /* ZERO_BND_EMF */
            case default
               write(msg,'(a,i3,a,i1,a,i3)') "[magboundaries:bnd_emf]: Boundary condition ",cg%bnd(dir, side)," not implemented for ",emfdir, " in ", dir
               if (master) call warn(msg)
         end select
      enddo

   end subroutine bnd_emf
!>
!! \brief Routine delivers common boundary cells indexes in cases of reflection or outflow boundary types
!<
   subroutine compute_bnd_indxs(bndcase, ndirb, edge, nbcells, rrbase, bndsign, zndiff)

      use constants, only: LO, HI, I_ONE, I_TWO
      use domain,    only: dom

      implicit none

      integer(kind=4),                   intent(in)  :: bndcase !> 1 - v component compatible with direction;
                                                                !! 2 - b component compatible with direction or emf component incompatible with direction;
                                                                !< 3 - other cases; BEWARE: magic integers
      integer(kind=4),                   intent(in)  :: ndirb   !< cg%{nxb,nyb,nzb} depanding on the current direction
      integer(kind=4), dimension(LO:HI), intent(out) :: edge    !< index of the left and right edge of physical domain for emf
      integer(kind=4), dimension(LO:HI), intent(out) :: nbcells !< number of cells in a loop at left and right boundaries
      integer(kind=4), dimension(LO:HI), intent(out) :: rrbase  !< COMMENT ME
      real,                              intent(out) :: bndsign !< 1. or -1. to change the sign or not
      logical,                           intent(out) :: zndiff  !< COMMENT ME

      select case (bndcase)
         case (1)
            edge(LO)    = dom%nb
            nbcells(LO) = dom%nb - I_ONE
            bndsign     = -1.
         case (2)
            edge(LO)    = dom%nb + I_ONE
            nbcells(LO) = dom%nb
            bndsign     = -1.
         case (3)
            edge(LO)    = dom%nb
            nbcells(LO) = dom%nb
            bndsign     = 1.
      end select  ! (bndcase)

      zndiff      = (edge(LO) - nbcells(LO) == 1)
!     nbcells(HI) = dom%nb - edge(LO) + nbcells(LO)
      nbcells(HI) = I_TWO * dom%nb - edge(LO)
      edge(HI)    = ndirb + edge(LO)
      rrbase(LO)  = edge(LO)
      rrbase(HI)  = ndirb + nbcells(LO) + I_ONE  ! = edge(HI) + 1 - edge(LO) + nbcells(LO)

   end subroutine compute_bnd_indxs

   subroutine all_mag_boundaries

      use constants,    only: xdim, zdim
      use domain,       only: dom
      use gc_list,      only: cg_list_element, all_cg
      use grid,         only: leaves
      use internal_bnd, only: internal_boundaries_4d
      use mpi,          only: MPI_COMM_NULL
      use types,        only: cdd

      implicit none

      type(cg_list_element), pointer :: cgl
      integer(kind=4) :: dir

      if (cdd%comm3d == MPI_COMM_NULL) then
         do dir = xdim, zdim
            if (dom%has_dir(dir)) call internal_boundaries_4d(all_cg, all_cg%bi, dim=dir) ! should be more selective (modified leaves?)
         enddo
      endif

      cgl => leaves%first
      do while (associated(cgl))
         do dir = xdim, zdim
            if (dom%has_dir(dir)) call bnd_b(dir, cgl%cg)
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine all_mag_boundaries

end module magboundaries
