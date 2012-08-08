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

module fluidboundaries
! pulled by ANY

   implicit none

   private
   public :: bnd_u, all_fluid_boundaries

contains

   subroutine init_fluidboundaries(cg)

      use constants,             only: PIERNIK_INIT_DOMAIN, xdim, zdim, LO, HI, &
           &                           BND_MPI, BND_PER, BND_REF, BND_OUT, BND_OUTD, BND_OUTH, BND_OUTHD, BND_COR, BND_SHE, BND_USER
      use dataio_pub,            only: msg, warn, die, code_progress
      use domain,                only: is_multicg
      use grid_cont,             only: grid_container

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer(kind=4)                           :: dir, side

      if (code_progress < PIERNIK_INIT_DOMAIN) call die("[fluidboundaries:init_fluidboundaries] MPI not initialized.") ! bnd_xl, bnd_xr

      do dir = xdim, zdim
         do side = LO, HI

            select case (cg%bnd(dir, side))
               case (BND_MPI, BND_REF, BND_OUT, BND_OUTD, BND_USER, BND_PER)
                  ! Do nothing
               case (BND_COR)
                  if (dir == zdim) then
                     write(msg,'("[fluid_boundaries:bnd_u]: corner ",i1," boundary condition ",i3," not implemented in ",i1,"-direction")') side, cg%bnd(dir, side), dir
                     call warn(msg)
                  endif
               case (BND_SHE)
                  if (dir /= xdim) then
                     write(msg,'("[fluid_boundaries:bnd_u]: shear ",i1," boundary condition ",i3," not implemented in ",i1,"-direction")') side, cg%bnd(dir, side), dir
                     call warn(msg)
                  endif
               case (BND_OUTH)
                  if (dir == zdim) then
                     if (is_multicg) call die("[fluid_boundaries:bnd_u] hydrostatic:outh_bnd with multiple grid pieces per processor not implemented yet") !nontrivial not really checked
                  else
                     write(msg,'("[fluid_boundaries:bnd_u]: outflow hydrostatic ",i1," boundary condition ",i3," not implemented in ",i1,"-direction")') side, cg%bnd(dir, side), dir
                     call warn(msg)
                  endif
               case (BND_OUTHD)
                  if (dir == zdim) then
                     if (is_multicg) call die("[fluid_boundaries:bnd_u] hydrostatic:outh_bnd with multiple grid pieces per processor not implemented yet") !nontrivial not really checked
                  else
                     write(msg,'("[fluid_boundaries:bnd_u]: outflow hydrostatic ",i1," boundary condition ",i3," not implemented in ",i1,"-direction")') side, cg%bnd(dir, side), dir
                     call warn(msg)
                  endif
               case default
                  write(msg,'("[fluid_boundaries:bnd_u]: unknown ",i1," boundary condition ",i3," not implemented in ",i1,"-direction")') side, cg%bnd(dir, side), dir
                  call warn(msg)
            end select
         enddo
      enddo


   end subroutine init_fluidboundaries

   subroutine bnd_u(dir, cg)

      use constants,             only: FLUID, ndims, xdim, ydim, zdim, LO, HI, BND, BLK, I_ONE, I_TWO, I_FOUR, &
           &                           BND_MPI, BND_PER, BND_REF, BND_OUT, BND_OUTD, BND_COR, BND_SHE, BND_USER, INT4
      use dataio_pub,            only: msg, warn, die
      use cart_comm,             only: cdd
      use domain,                only: dom, is_multicg
      use fluidboundaries_funcs, only: user_fluidbnd
      use fluidindex,            only: flind, iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
      use grid_cont,             only: grid_container
      use mpi,                   only: MPI_DOUBLE_PRECISION, MPI_COMM_NULL
      use mpisetup,              only: mpi_err, req, status, comm
#ifdef COSM_RAYS
      use initcosmicrays,        only: smallecr
      use fluidindex,            only: iarr_all_crs
#endif /* COSM_RAYS */
#ifndef ISO
      use fluidindex,            only: iarr_all_en
#endif /* !ISO */
#ifdef SHEAR_BND
      use shear,                 only: dely, unshear_fft
#ifndef FFTW
      use constants,             only: half
      use global,                only: smalld
      use shear,                 only: qshear, delj, eps, omega
#endif /* !FFTW */
#endif /* SHEAR_BND */
#ifdef GRAV
      use constants,             only: BND_OUTH, BND_OUTHD
      use hydrostatic,           only: outh_bnd
#endif /* GRAV */

      implicit none

      integer(kind=4),               intent(in)    :: dir
      type(grid_container), pointer, intent(inout) :: cg

      integer(kind=4), parameter                   :: tag1 = 10, tag2 = 20
      integer(kind=4), parameter                   :: tag7 = 70, tag8 = 80
      integer(kind=4), dimension(ndims,LO:HI)      :: l, r, seb
      logical, save                                :: frun = .true.
      integer                                      :: i, j
      integer(kind=4)                              :: itag, jtag, side, ssign, ib
      real, allocatable                            :: send_left(:,:,:,:),recv_left(:,:,:,:)
#ifdef SHEAR_BND
      real, allocatable                            :: send_right(:,:,:,:),recv_right(:,:,:,:)
#endif /* SHEAR_BND */

      if (.not. any([xdim, ydim, zdim] == dir)) call die("[fluidboundaries:bnd_u] Invalid direction.")

      if (frun) then
         call init_fluidboundaries(cg)
         frun = .false.
      endif

! MPI block communication
      if (cdd%comm3d /= MPI_COMM_NULL) then
         if (is_multicg) call die("[fluidboundaries:bnd_u] multiple grid pieces per processor not implemented for comm3d")
#ifdef SHEAR_BND
         if (dir == xdim) then
#ifndef FFTW
            allocate(send_right(flind%all, dom%nb,ny,nz), send_left (flind%all, dom%nb,ny,nz), &
                 &    recv_left(flind%all, dom%nb,ny,nz), recv_right(flind%all, dom%nb,ny,nz) )
            send_left(:,:,:,:)          =  cg%u(:, cg%is:cg%isb,:,:)
            send_right(:,:,:,:)         =  cg%u(:, cg%ieb:cg%ie,:,:)
!
! odejmujemy ped_y i energie odpowiadajace niezaburzonej rozniczkowej rotacji na lewym brzegu
!
            if (cg%bnd(xdim, LO) == BND_SHE) then
               do i=1, dom%nb
                  send_left (iarr_all_my,i,:,:) = send_left(iarr_all_my,i,:,:) + qshear*omega * cg%x(dom%nb+i) * send_left(iarr_all_dn,i,:,:)
#ifndef ISO
                  send_left (iarr_all_en,i,:,:) = send_left(iarr_all_en,i,:,:) - half*(qshear*omega * cg%x(dom%nb+i))**2 * send_left(iarr_all_dn,i,:,:)
#endif /* !ISO */
               enddo
!
! przesuwamy o calkowita liczbe komorek + periodyczny wb w kierunku y
!
               if (dom%has_dir(ydim)) then
                  send_left (:,:, cg%js:cg%je,  :) = cshift(send_left (:,:, cg%js:cg%je,:),dim=3,shift= delj)
                  send_left (:,:, 1:dom%nb,      :) = send_left (:,:, cg%jeb:cg%je,:)
                  send_left (:,:, cg%je+1:cg%n_(ydim),:) = send_left (:,:, cg%js:cg%jsb,:)
               endif
!
! remapujemy  - interpolacja kwadratowa
!
               send_left (:,:,:,:)  = (1.+eps)*(1.-eps) * send_left (:,:,:,:) &
                    &                 -half*eps*(1.-eps) * cshift(send_left(:,:,:,:),shift=-1,dim=3) &
                    &                 +half*eps*(1.+eps) * cshift(send_left(:,:,:,:),shift=1 ,dim=3)
            endif !(cg%bnd(xdim, LO) == BND_SHE)
!
! odejmujemy ped_y i energie odpowiadajace niezaburzonej rozniczkowej rotacji na prawym brzegu
!
            if (cg%bnd(xdim, HI) == BND_SHE) then
               do i=1, dom%nb
                  send_right(iarr_all_my,i,:,:) = send_right(iarr_all_my,i,:,:) +qshear*omega * cg%x(cg%nxb+i) * send_right(iarr_all_dn,i,:,:)
#ifndef ISO
                  send_right(iarr_all_en,i,:,:) = send_right(iarr_all_en,i,:,:) -half*(qshear*omega * cg%x(cg%nxb+i))**2 * send_right(iarr_all_dn,i,:,:)
#endif /* !ISO */
               enddo
!
! przesuwamy o calkowita liczbe komorek + periodyczny wb w kierunku y
!
               if (dom%has_dir(ydim)) then
                  send_right (:,:, cg%js:cg%je,  :) = cshift(send_right(:,:, cg%js:cg%je,:),dim=3,shift=-delj)
                  send_right (:,:, 1:dom%nb,      :) = send_right(:,:, cg%jeb:cg%je,:)
                  send_right (:,:, cg%je+1:cg%n_(ydim),:) = send_right(:,:, cg%js:cg%jsb,:)
               endif
!
! remapujemy  - interpolacja kwadratowa
!
               send_right (:,:,:,:) = (1.+eps)*(1.-eps) * send_right (:,:,:,:) &
                    &                 -half*eps*(1.-eps) * cshift(send_right(:,:,:,:),shift=1 ,dim=3) &
                    &                 +half*eps*(1.+eps) * cshift(send_right(:,:,:,:),shift=-1,dim=3)
            endif !(cg%bnd(xdim, HI) == BND_SHE)
!
! wysylamy na drugi brzeg
!
            call MPI_Isend(send_left , flind%all*ny*nz*dom%nb, MPI_DOUBLE_PRECISION, cdd%procn(dir,LO), tag1, comm, req(1), mpi_err)
            call MPI_Isend(send_right, flind%all*ny*nz*dom%nb, MPI_DOUBLE_PRECISION, cdd%procn(dir,HI), tag2, comm, req(3), mpi_err)
            call MPI_Irecv(recv_left , flind%all*ny*nz*dom%nb, MPI_DOUBLE_PRECISION, cdd%procn(dir,LO), tag2, comm, req(2), mpi_err)
            call MPI_Irecv(recv_right, flind%all*ny*nz*dom%nb, MPI_DOUBLE_PRECISION, cdd%procn(dir,HI), tag1, comm, req(4), mpi_err)

            call MPI_Waitall(I_FOUR,req(:),status(:,:),mpi_err)

!
! dodajemy ped_y i energie odpowiadajace niezaburzonej rozniczkowej rotacji na prawym brzegu
!
            if (cg%bnd(xdim, HI) == BND_SHE) then
               do i=1, dom%nb
#ifndef ISO
                  recv_right (iarr_all_en,i,:,:) = recv_right (iarr_all_en,i,:,:) +half*(qshear*omega * cg%x(cg%ie+i))**2 * recv_right(iarr_all_dn,i,:,:)
#endif /* !ISO */
                  recv_right (iarr_all_my,i,:,:) = recv_right (iarr_all_my,i,:,:) -qshear*omega * cg%x(cg%ie+i)     * recv_right(iarr_all_dn,i,:,:)
               enddo
            endif !(cg%bnd(xdim, HI) == BND_SHE)
!
! dodajemy ped_y i energie odpowiadajace niezaburzonej rozniczkowej rotacji na lewym brzegu
!
            if (cg%bnd(xdim, LO) == BND_SHE) then
               do i=1, dom%nb
#ifndef ISO
                  recv_left(iarr_all_en,i,:,:) = recv_left(iarr_all_en,i,:,:) +half*(qshear*omega * cg%x(i))**2 * recv_left(iarr_all_dn,i,:,:)
#endif /* !ISO */
                  recv_left(iarr_all_my,i,:,:) = recv_left(iarr_all_my,i,:,:) -qshear*omega * cg%x(i)     * recv_left(iarr_all_dn,i,:,:)
               enddo
            endif !(cg%bnd(xdim, LO) == BND_SHE)

            cg%u(:, 1:dom%nb,            :,:) = recv_left (:,1:dom%nb,:,:)
            cg%u(:, cg%ie+1:cg%n_(xdim),:,:) = recv_right(:,1:dom%nb,:,:)

            !> \deprecated BEWARE: smalld is called only for the first fluid
            cg%u(iarr_all_dn(1), 1:dom%nb,           :, :) = max(cg%u(iarr_all_dn(1),       1:dom%nb,      :,:), smalld)
            cg%u(iarr_all_dn(1), cg%ie+1:cg%n_(xdim),:,:) = max(cg%u(iarr_all_dn(1), cg%ie+1:cg%n_(xdim),:,:), smalld)
            if (allocated(send_left))  deallocate(send_left)
            if (allocated(send_right)) deallocate(send_right)
            if (allocated(recv_left))  deallocate(recv_left)
            if (allocated(recv_right)) deallocate(recv_right)

#else /* FFTW */

            if ( all(cg%bnd(xdim, LO:HI) == BND_SHE) ) then

               if (allocated(send_right)) deallocate(send_right)
               if (.not.allocated(send_right)) allocate(send_right(flind%all, dom%nb, cg%nyb, cg%n_(zdim)))

               if (allocated(send_left)) deallocate(send_left)
               if (.not.allocated(send_left)) allocate(send_left(flind%all, dom%nb, cg%nyb, cg%n_(zdim)))

               if (allocated(recv_left)) deallocate(recv_left)
               if (.not.allocated(recv_left)) allocate(recv_left(flind%all, dom%nb, cg%nyb, cg%n_(zdim)))

               if (allocated(recv_right)) deallocate(recv_right)
               if (.not.allocated(recv_right)) allocate(recv_right(flind%all, dom%nb, cg%nyb, cg%n_(zdim)))

               do i = lbound(cg%u,1), ubound(cg%u,1)
                  send_left( i,1:dom%nb,:,:) = unshear_fft(cg%u(i, cg%is:cg%isb, cg%js:cg%je,:), cg%x(cg%is:cg%isb),dely,.true.)
                  send_right(i,1:dom%nb,:,:) = unshear_fft(cg%u(i, cg%ieb:cg%ie, cg%js:cg%je,:), cg%x(cg%ieb:cg%ie),dely,.true.)
               enddo

               call MPI_Isend(send_left , flind%all*cg%nyb*cg%n_(zdim)*dom%nb, MPI_DOUBLE_PRECISION, cdd%procn(dir,LO), tag1, comm, req(1), mpi_err)
               call MPI_Isend(send_right, flind%all*cg%nyb*cg%n_(zdim)*dom%nb, MPI_DOUBLE_PRECISION, cdd%procn(dir,HI), tag2, comm, req(3), mpi_err)
               call MPI_Irecv(recv_left , flind%all*cg%nyb*cg%n_(zdim)*dom%nb, MPI_DOUBLE_PRECISION, cdd%procn(dir,LO), tag2, comm, req(2), mpi_err)
               call MPI_Irecv(recv_right, flind%all*cg%nyb*cg%n_(zdim)*dom%nb, MPI_DOUBLE_PRECISION, cdd%procn(dir,HI), tag1, comm, req(4), mpi_err)

               call MPI_Waitall(I_FOUR, req(:),status(:,:),mpi_err)

               do i = lbound(cg%u,1), ubound(cg%u,1)
                  cg%u(i,1:dom%nb,        cg%js:cg%je,:) = unshear_fft(recv_left (i,1:dom%nb,:,:), cg%x(1:dom%nb),dely)
                  cg%u(i, cg%ie+1:cg%n_(xdim), cg%js:cg%je,:) = unshear_fft(recv_right(i,1:dom%nb,:,:), cg%x(cg%ie+1:cg%n_(xdim)),dely)
               enddo

               if (allocated(send_left))  deallocate(send_left)
               if (allocated(send_right)) deallocate(send_right)
               if (allocated(recv_left))  deallocate(recv_left)
               if (allocated(recv_right)) deallocate(recv_right)

            endif
#endif /* FFTW */
         else
#endif /* SHEAR_BND */
            if (cdd%psize(dir) > 1) then

               jtag = tag2 * dir
               itag = jtag - tag1
               call MPI_Isend(cg%u(1,1,1,1), I_ONE, cg%mbc(FLUID, dir, LO, BLK, dom%nb), cdd%procn(dir,LO), itag, cdd%comm3d, req(1), mpi_err)
               call MPI_Isend(cg%u(1,1,1,1), I_ONE, cg%mbc(FLUID, dir, HI, BLK, dom%nb), cdd%procn(dir,HI), jtag, cdd%comm3d, req(2), mpi_err)
               call MPI_Irecv(cg%u(1,1,1,1), I_ONE, cg%mbc(FLUID, dir, LO, BND, dom%nb), cdd%procn(dir,LO), jtag, cdd%comm3d, req(3), mpi_err)
               call MPI_Irecv(cg%u(1,1,1,1), I_ONE, cg%mbc(FLUID, dir, HI, BND, dom%nb), cdd%procn(dir,HI), itag, cdd%comm3d, req(4), mpi_err)

               call MPI_Waitall(I_FOUR, req(:),status(:,:),mpi_err)
            endif
#ifdef SHEAR_BND
         endif
#endif /* SHEAR_BND */

! MPI + non-MPI corner-periodic boundary condition

         if (cg%bnd(xdim, LO) == BND_COR) then
!   - lower to left
            if (cdd%pcoords(xdim) == 0 .and. cdd%pcoords(ydim) == 0) then
               do i=1, dom%nb
                  do j=cg%js, cg%n_(ydim)
                     cg%u(iarr_all_dn,i,j,:) =  cg%u(iarr_all_dn,j,cg%isb+1-i,:)
                     cg%u(iarr_all_mx,i,j,:) = -cg%u(iarr_all_my,j,cg%isb+1-i,:)
                     cg%u(iarr_all_my,i,j,:) =  cg%u(iarr_all_mx,j,cg%isb+1-i,:)
                     cg%u(iarr_all_mz,i,j,:) =  cg%u(iarr_all_mz,j,cg%isb+1-i,:)
#ifndef ISO
                     cg%u(iarr_all_en,i,j,:) =  cg%u(iarr_all_en,j,cg%isb+1-i,:)
#endif /* !ISO */
#ifdef COSM_RAYS
                     cg%u(iarr_all_crs,i,j,:) =  cg%u(iarr_all_crs,j,cg%isb+1-i,:)
#endif /* COSM_RAYS */
                  enddo
               enddo
            endif

            if (cdd%procxyl > 0) then
               allocate(send_left(flind%all, dom%nb, cg%n_(ydim), cg%n_(zdim)), recv_left(flind%all, cg%n_(xdim), dom%nb, cg%n_(zdim)))

               send_left(:,:,:,:) = cg%u(:, cg%is:cg%isb,:,:)

               call MPI_Isend(send_left, flind%all*dom%nb*cg%n_(ydim)*cg%n_(zdim), MPI_DOUBLE_PRECISION, cdd%procxyl, tag7, comm, req(1), mpi_err)
               call MPI_Irecv(recv_left, flind%all*cg%n_(xdim)*dom%nb*cg%n_(zdim), MPI_DOUBLE_PRECISION, cdd%procxyl, tag8, comm, req(2), mpi_err)

               call MPI_Waitall(I_TWO,req(:),status(:,:),mpi_err)

               do i=1, dom%nb
                  do j=1, cg%n_(ydim)
                     cg%u(iarr_all_dn,i,j,:) =  recv_left(iarr_all_dn,j, cg%is-i,:)
                     cg%u(iarr_all_mx,i,j,:) = -recv_left(iarr_all_my,j, cg%is-i,:)
                     cg%u(iarr_all_my,i,j,:) =  recv_left(iarr_all_mx,j, cg%is-i,:)
                     cg%u(iarr_all_mz,i,j,:) =  recv_left(iarr_all_mz,j, cg%is-i,:)
#ifndef ISO
                     cg%u(iarr_all_en,i,j,:) =  recv_left(iarr_all_en,j, cg%is-i,:)
#endif /* !ISO */
#ifdef COSM_RAYS
                     cg%u(iarr_all_crs,i,j,:) =  recv_left(iarr_all_crs,j, cg%is-i,:)
#endif /* COSM_RAYS */
                  enddo
               enddo

               if (allocated(send_left))  deallocate(send_left)
               if (allocated(recv_left))  deallocate(recv_left)
            endif
         endif

         if (cg%bnd(ydim, LO) == BND_COR) then
!   - left to lower
            if (cdd%pcoords(ydim) == 0 .and. cdd%pcoords(xdim) == 0 ) then
               do j=1, dom%nb
                  do i=cg%is, cg%n_(xdim)
                     cg%u(iarr_all_dn,i,j,:) =  cg%u(iarr_all_dn,cg%isb+1-j,i,:)
                     cg%u(iarr_all_mx,i,j,:) =  cg%u(iarr_all_my,cg%isb+1-j,i,:)
                     cg%u(iarr_all_my,i,j,:) = -cg%u(iarr_all_mx,cg%isb+1-j,i,:)
                     cg%u(iarr_all_mz,i,j,:) =  cg%u(iarr_all_mz,cg%isb+1-j,i,:)
#ifndef ISO
                     cg%u(iarr_all_en,i,j,:) =  cg%u(iarr_all_en,cg%isb+1-j,i,:)
#endif /* !ISO */
#ifdef COSM_RAYS
                     cg%u(iarr_all_crs,i,j,:) =  cg%u(iarr_all_crs,cg%isb+1-j,i,:)
#endif /* COSM_RAYS */
                  enddo
               enddo
!   - interior to corner
               do j=1, dom%nb
                  do i=1, dom%nb
                     cg%u(iarr_all_dn,i,j,:) =   cg%u(iarr_all_dn,cg%isb+1-i,cg%jsb+1-j,:)
                     cg%u(iarr_all_mx,i,j,:) =  -cg%u(iarr_all_mx,cg%isb+1-i,cg%jsb+1-j,:)
                     cg%u(iarr_all_my,i,j,:) =  -cg%u(iarr_all_my,cg%isb+1-i,cg%jsb+1-j,:)
                     cg%u(iarr_all_mz,i,j,:) =   cg%u(iarr_all_mz,cg%isb+1-i,cg%jsb+1-j,:)
#ifndef ISO
                     cg%u(iarr_all_en,i,j,:) =   cg%u(iarr_all_en,cg%isb+1-i,cg%jsb+1-j,:)
#endif /* !ISO */
#ifdef COSM_RAYS
                     cg%u(iarr_all_crs,i,j,:) =   cg%u(iarr_all_crs,cg%isb+1-i,cg%jsb+1-j,:)
#endif /* COSM_RAYS */
                  enddo
               enddo
            endif

            if (cdd%procyxl > 0) then
               allocate(send_left(flind%all, cg%n_(xdim), dom%nb, cg%n_(zdim)), recv_left(flind%all, dom%nb, cg%n_(ydim), cg%n_(zdim)))

               send_left(:,:,:,:) = cg%u(:,:, cg%js:cg%jsb,:)

               call MPI_Isend(send_left, flind%all*cg%n_(xdim)*dom%nb*cg%n_(zdim), MPI_DOUBLE_PRECISION, cdd%procyxl, tag8, comm, req(1), mpi_err)
               call MPI_Irecv(recv_left, flind%all*dom%nb*cg%n_(ydim)*cg%n_(zdim), MPI_DOUBLE_PRECISION, cdd%procyxl, tag7, comm, req(2), mpi_err)

               call MPI_Waitall(I_TWO, req(:),status(:,:),mpi_err)

               do j=1, dom%nb
                  do i=1, cg%n_(xdim)
                     cg%u(iarr_all_dn,i,j,:) =  recv_left(iarr_all_dn, cg%js-j,i,:)
                     cg%u(iarr_all_mx,i,j,:) =  recv_left(iarr_all_my, cg%js-j,i,:)
                     cg%u(iarr_all_my,i,j,:) = -recv_left(iarr_all_mx, cg%js-j,i,:)
                     cg%u(iarr_all_mz,i,j,:) =  recv_left(iarr_all_mz, cg%js-j,i,:)
#ifndef ISO
                     cg%u(iarr_all_en,i,j,:) =  recv_left(iarr_all_en, cg%js-j,i,:)
#endif /* !ISO */
#ifdef COSM_RAYS
                     cg%u(iarr_all_crs,i,j,:) =  recv_left(iarr_all_crs, cg%js-j,i,:)
#endif /* COSM_RAYS */
                  enddo
               enddo

               if (allocated(send_left))  deallocate(send_left)
               if (allocated(recv_left))  deallocate(recv_left)
            endif
         endif
      endif

!===============================================================

! Non-MPI boundary conditions
      l = reshape([lbound(cg%u(xdim,:,:,:), kind=4),ubound(cg%u(xdim,:,:,:), kind=4)],shape=[ndims,HI]) ; r = l

      do side = LO, HI

         select case (cg%bnd(dir, side))
         case (BND_MPI, BND_COR, BND_SHE)
            ! Do nothing
         case (BND_USER)
            call user_fluidbnd(dir,side,cg)
         case (BND_PER)
            if (cdd%comm3d /= MPI_COMM_NULL) then
               seb = reshape([[cg%isb, cg%jsb, cg%ksb],[cg%ieb, cg%jeb, cg%keb]],[ndims,HI])
               r(dir, side) = seb(dir,3_INT4-side)
               r(dir,3_INT4-side) = cg%ijkse(dir,3_INT4-side)
               l(dir,:) = [1_INT4, dom%nb] + cg%ijkse(dir,side)*(side-1_INT4)
               cg%u(:,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = cg%u(:,r(xdim,LO):r(xdim,HI),r(ydim,LO):r(ydim,HI),r(zdim,LO):r(zdim,HI))
            endif
         case (BND_REF)
            ssign = 2_INT4*side-3_INT4
            do ib=1_INT4, dom%nb
               l(dir,:) = cg%ijkse(dir,side)+ssign*ib ; r(dir,:) = cg%ijkse(dir,side)+ssign*(1_INT4-ib)
               cg%u(:,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = cg%u(:,r(xdim,LO):r(xdim,HI),r(ydim,LO):r(ydim,HI),r(zdim,LO):r(zdim,HI))
               cg%u(iarr_all_dn+dir, l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = -cg%u(iarr_all_dn+dir,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI))
            enddo
         case (BND_OUT)
            r(dir,:) = cg%ijkse(dir,side)
            ssign = 2_INT4*side-3_INT4
            do ib=1_INT4, dom%nb
               l(dir,:) = cg%ijkse(dir,side)+ssign*ib
               cg%u(:,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = cg%u(:,r(xdim,LO):r(xdim,HI),r(ydim,LO):r(ydim,HI),r(zdim,LO):r(zdim,HI))
#ifdef COSM_RAYS
               cg%u(iarr_all_crs,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = smallecr
#endif /* COSM_RAYS */
            enddo
         case (BND_OUTD)
            r(dir,:) = cg%ijkse(dir,side)
            ssign = 2_INT4*side-3_INT4
            do ib=1_INT4, dom%nb
               l(dir,:) = cg%ijkse(dir,side)+ssign*ib
               cg%u(:,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = cg%u(:,r(xdim,LO):r(xdim,HI),r(ydim,LO):r(ydim,HI),r(zdim,LO):r(zdim,HI))
!> \deprecated BEWARE: use of uninitialized value on first call (a side effect of r1726)
#ifdef COSM_RAYS
               cg%u(iarr_all_crs,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = smallecr
#endif /* COSM_RAYS */
            enddo
            l(dir,:) = [1_INT4, dom%nb] + cg%ijkse(dir,side)*(side-1_INT4)
            if (side == LO) then
               cg%u(iarr_all_dn+dir,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = min(cg%u(iarr_all_dn+dir,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)),0.0)
            else
               cg%u(iarr_all_dn+dir,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)) = max(cg%u(iarr_all_dn+dir,l(xdim,LO):l(xdim,HI),l(ydim,LO):l(ydim,HI),l(zdim,LO):l(zdim,HI)),0.0)
            endif
#ifdef GRAV
         case (BND_OUTH)
            if (dir == zdim) call outh_bnd(side, cg, .false.)
         case (BND_OUTHD)
            if (dir == zdim) call outh_bnd(side, cg, .true.)
#endif /* GRAV */
         case default
            write(msg,'("[fluid_boundaries:bnd_u]: Unrecognized ",i1," boundary condition ",i3," not implemented in ",i1,"-direction")') side, cg%bnd(dir, side), dir
            call warn(msg)
         end select

      enddo

   end subroutine bnd_u

   subroutine all_fluid_boundaries

      use cg_list,        only: cg_list_element
      use cg_list_bnd,    only: leaves
      use cg_list_global, only: all_cg
      use constants,      only: xdim, zdim
      use cart_comm,      only: cdd
      use domain,         only: dom
      use mpi,            only: MPI_COMM_NULL

      implicit none

      type(cg_list_element), pointer :: cgl
      integer(kind=4) :: dir

      if (cdd%comm3d == MPI_COMM_NULL) then
         do dir = xdim, zdim
            if (dom%has_dir(dir)) call all_cg%internal_boundaries_4d(all_cg%fi, dim=dir) ! should be more selective (modified leaves?)
         enddo
      endif

      cgl => leaves%first
      do while (associated(cgl))
         do dir = xdim, zdim
            if (dom%has_dir(dir)) call bnd_u(dir, cgl%cg)
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine all_fluid_boundaries

end module fluidboundaries
