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

      use dataio_pub,            only: msg, warn, die, code_progress
      use constants,             only: PIERNIK_INIT_MPI, xdim, LO, HI, &
           &                           BND_MPI, BND_PER, BND_REF, BND_OUT, BND_OUTD, BND_COR, BND_SHE, BND_INF, BND_USER
      use fluidboundaries_funcs, only: bnd_null, bnd_xl_per, bnd_xl_ref, bnd_xl_out, bnd_xl_outd, bnd_xr_per, bnd_xr_ref, bnd_xr_out, bnd_xr_outd
      use fluidboundaries_pub,   only: user_bnd_xl, user_bnd_xr, func_bnd_xl, func_bnd_xr
      use grid_cont,             only: grid_container
      use mpi,                   only: MPI_COMM_NULL
      use mpisetup,              only: comm3d

      implicit none

      type(grid_container), pointer, intent(in) :: cg

      if (code_progress < PIERNIK_INIT_MPI) call die("[fluidboundaries:init_fluidboundaries] MPI not initialized.") ! bnd_xl, bnd_xr

      select case (cg%bnd(xdim, LO))
         case (BND_COR, BND_MPI, BND_SHE, BND_INF)
            func_bnd_xl => bnd_null
         case (BND_PER)
            if (comm3d == MPI_COMM_NULL) then
               func_bnd_xl => bnd_null
            else
               func_bnd_xl => bnd_xl_per
            endif
         case (BND_USER)
            func_bnd_xl => user_bnd_xl
         case (BND_REF)
            func_bnd_xl => bnd_xl_ref
         case (BND_OUT)
            func_bnd_xl => bnd_xl_out
         case (BND_OUTD)
            func_bnd_xl => bnd_xl_outd
         case default
            func_bnd_xl => bnd_null
            write(msg,'("[fluid_boundaries:init_fluidboundaries]: Left boundary condition ",i3," not implemented in X-direction")') cg%bnd(xdim, LO)
            call warn(msg)
      end select

      select case (cg%bnd(xdim, HI))
         case (BND_COR, BND_MPI, BND_SHE, BND_INF)
            func_bnd_xr => bnd_null
         case (BND_PER)
            if (comm3d == MPI_COMM_NULL) then
               func_bnd_xr => bnd_null
            else
               func_bnd_xr => bnd_xr_per
            endif
         case (BND_USER)
            func_bnd_xr => user_bnd_xr
         case (BND_REF)
            func_bnd_xr => bnd_xr_ref
         case (BND_OUT)
            func_bnd_xr => bnd_xr_out
         case (BND_OUTD)
            func_bnd_xr => bnd_xr_outd
         case default
            func_bnd_xr => bnd_null
            write(msg,'("[fluid_boundaries:init_fluidboundaries]: Right boundary condition ",i3," not implemented in X-direction")') cg%bnd(xdim, HI)
            call warn(msg)
      end select

   end subroutine init_fluidboundaries

   subroutine bnd_u(dir, cg)

      use constants,           only: FLUID, xdim, ydim, zdim, LO, HI, BND, BLK, BND_MPI, BND_PER, BND_REF, BND_OUT, BND_OUTD, BND_OUTH, BND_COR, BND_SHE, BND_INF, BND_USER
      use dataio_pub,          only: msg, warn, die
      use fluidboundaries_pub, only: user_bnd_yl, user_bnd_yr, user_bnd_zl, user_bnd_zr, func_bnd_xl, func_bnd_xr
      use fluidindex,          only: flind, iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
      use grid,                only: cga
      use grid_cont,           only: grid_container
      use mpisetup,            only: ierr, psize, procn, procxyl, procyxl, smalld, pcoords, req, status, comm, comm3d, proc, has_dir
      use mpi,                 only: MPI_DOUBLE_PRECISION, MPI_COMM_NULL
#ifdef COSM_RAYS
      use initcosmicrays,      only: smallecr
      use fluidindex,          only: iarr_all_crs
#endif /* COSM_RAYS */
#ifndef ISO
      use fluidindex,          only: iarr_all_en
#endif /* !ISO */
#ifdef SHEAR_BND
      use shear,               only: qshear, omega, delj, eps, dely, unshear_fft
#endif /* SHEAR_BND */
#ifdef GRAV
      use hydrostatic,         only: outh_bnd
#endif /* GRAV */

      implicit none

      integer, intent(in) :: dir
      type(grid_container), pointer, intent(in) :: cg

      logical, save    :: frun = .true.
      integer :: i, j, ib, itag, jtag
      real, allocatable :: send_left(:,:,:,:),recv_left(:,:,:,:)
#ifdef GRAV
      integer          :: kb
#endif /* GRAV */
#ifdef SHEAR_BND
      real, allocatable :: send_right(:,:,:,:),recv_right(:,:,:,:)
#endif /* SHEAR_BND */

      if (.not. any([xdim, ydim, zdim] == dir)) call die("[fluidboundaries:bnd_u] Invalid direction.")

      if (frun) then
         call init_fluidboundaries(cg)
         frun = .false.
      endif

      if (ubound(cga%cg_all(:), dim=1) > 1) call die("[fluidboundaries:bnd_u] multiple grid pieces per procesor not implemented yet") !nontrivial communication

! MPI block communication
      if (comm3d /= MPI_COMM_NULL) then
#ifdef SHEAR_BND
         if (dir == xdim) then
#ifndef FFTW
            allocate(send_right(flind%all, cg%nb,ny,nz), send_left(flind%all, cg%nb,ny,nz), &
                 &    recv_left(flind%all, cg%nb,ny,nz), recv_right(flind%all, cg%nb,ny,nz) )
            send_left(:,:,:,:)          =  cg%u%arr(:, cg%is:cg%isb,:,:)
            send_right(:,:,:,:)         =  cg%u%arr(:, cg%ieb:cg%ie,:,:)
!
! odejmujemy ped_y i energie odpowiadajace niezaburzonej rozniczkowej rotacji na lewym brzegu
!
            if (cg%bnd(xdim, LO) == BND_SHE) then
               do i=1, cg%nb
                  send_left (iarr_all_my,i,:,:) = send_left(iarr_all_my,i,:,:) + qshear*omega * cg%x(cg%nb+i) * send_left(iarr_all_dn,i,:,:)
#ifndef ISO
                  send_left (iarr_all_en,i,:,:) = send_left(iarr_all_en,i,:,:) - 0.5*(qshear*omega * cg%x(cg%nb+i))**2 * send_left(iarr_all_dn,i,:,:)
#endif /* !ISO */
               enddo
!
! przesuwamy o calkowita liczbe komorek + periodyczny wb w kierunku y
!
               if (has_dir(ydim)) then
                  send_left (:,:, cg%js:cg%je,  :) = cshift(send_left (:,:, cg%js:cg%je,:),dim=3,shift= delj)
                  send_left (:,:, 1:cg%nb,      :) = send_left (:,:, cg%jeb:cg%je,:)
                  send_left (:,:, cg%je+1:cg%ny,:) = send_left (:,:, cg%js:cg%jsb,:)
               endif
!
! remapujemy  - interpolacja kwadratowa
!
               send_left (:,:,:,:)  = (1.+eps)*(1.-eps) * send_left (:,:,:,:) &
                    &                 -0.5*eps*(1.-eps) * cshift(send_left(:,:,:,:),shift=-1,dim=3) &
                    &                 +0.5*eps*(1.+eps) * cshift(send_left(:,:,:,:),shift=1 ,dim=3)
            endif !(cg%bnd(xdim, LO) == BND_SHE)
!
! odejmujemy ped_y i energie odpowiadajace niezaburzonej rozniczkowej rotacji na prawym brzegu
!
            if (cg%bnd(xdim, HI) == BND_SHE) then
               do i=1, cg%nb
                  send_right(iarr_all_my,i,:,:) = send_right(iarr_all_my,i,:,:) +qshear*omega * cg%x(cg%nxb+i) * send_right(iarr_all_dn,i,:,:)
#ifndef ISO
                  send_right(iarr_all_en,i,:,:) = send_right(iarr_all_en,i,:,:) -0.5*(qshear*omega * cg%x(cg%nxb+i))**2 * send_right(iarr_all_dn,i,:,:)
#endif /* !ISO */
               enddo
!
! przesuwamy o calkowita liczbe komorek + periodyczny wb w kierunku y
!
               if (has_dir(ydim)) then
                  send_right (:,:, cg%js:cg%je,  :) = cshift(send_right(:,:, cg%js:cg%je,:),dim=3,shift=-delj)
                  send_right (:,:, 1:cg%nb,      :) = send_right(:,:, cg%jeb:cg%je,:)
                  send_right (:,:, cg%je+1:cg%ny,:) = send_right(:,:, cg%js:cg%jsb,:)
               endif
!
! remapujemy  - interpolacja kwadratowa
!
               send_right (:,:,:,:) = (1.+eps)*(1.-eps) * send_right (:,:,:,:) &
                    &                 -0.5*eps*(1.-eps) * cshift(send_right(:,:,:,:),shift=1 ,dim=3) &
                    &                 +0.5*eps*(1.+eps) * cshift(send_right(:,:,:,:),shift=-1,dim=3)
            endif !(cg%bnd(xdim, HI) == BND_SHE)
!
! wysylamy na drugi brzeg
!
            call MPI_Isend(send_left , flind%all*ny*nz*cg%nb, MPI_DOUBLE_PRECISION, procn(dir,LO), 10, comm, req(1), ierr)
            call MPI_Isend(send_right, flind%all*ny*nz*cg%nb, MPI_DOUBLE_PRECISION, procn(dir,HI), 20, comm, req(3), ierr)
            call MPI_Irecv(recv_left , flind%all*ny*nz*cg%nb, MPI_DOUBLE_PRECISION, procn(dir,LO), 20, comm, req(2), ierr)
            call MPI_Irecv(recv_right, flind%all*ny*nz*cg%nb, MPI_DOUBLE_PRECISION, procn(dir,HI), 10, comm, req(4), ierr)

            call MPI_Waitall(4,req(:),status(:,:),ierr)

!
! dodajemy ped_y i energie odpowiadajace niezaburzonej rozniczkowej rotacji na prawym brzegu
!
            if (cg%bnd(xdim, HI) == BND_SHE) then
               do i=1, cg%nb
#ifndef ISO
                  recv_right (iarr_all_en,i,:,:) = recv_right (iarr_all_en,i,:,:) +0.5*(qshear*omega * cg%x(cg%ie+i))**2 * recv_right(iarr_all_dn,i,:,:)
#endif /* !ISO */
                  recv_right (iarr_all_my,i,:,:) = recv_right (iarr_all_my,i,:,:) -qshear*omega * cg%x(cg%ie+i)     * recv_right(iarr_all_dn,i,:,:)
               enddo
            endif !(cg%bnd(xdim, HI) == BND_SHE)
!
! dodajemy ped_y i energie odpowiadajace niezaburzonej rozniczkowej rotacji na lewym brzegu
!
            if (cg%bnd(xdim, LO) == BND_SHE) then
               do i=1, cg%nb
#ifndef ISO
                  recv_left(iarr_all_en,i,:,:) = recv_left(iarr_all_en,i,:,:) +0.5*(qshear*omega * cg%x(i))**2 * recv_left(iarr_all_dn,i,:,:)
#endif /* !ISO */
                  recv_left(iarr_all_my,i,:,:) = recv_left(iarr_all_my,i,:,:) -qshear*omega * cg%x(i)     * recv_left(iarr_all_dn,i,:,:)
               enddo
            endif !(cg%bnd(xdim, LO) == BND_SHE)

            cg%u%arr(:, 1:cg%nb,      :,:) = recv_left (:,1:cg%nb,:,:)
            cg%u%arr(:, cg%ie+1:cg%nx,:,:) = recv_right(:,1:cg%nb,:,:)

            !> \deprecated BEWARE: smalld is called only for the first fluid
            cg%u%arr(iarr_all_dn(1), 1:cg%nb, :, :)     = max(cg%u%arr(iarr_all_dn(1),       1:cg%nb,:,:), smalld)
            cg%u%arr(iarr_all_dn(1), cg%ie+1:cg%nx,:,:) = max(cg%u%arr(iarr_all_dn(1), cg%ie+1:cg%nx,:,:), smalld)
            if (allocated(send_left))  deallocate(send_left)
            if (allocated(send_right)) deallocate(send_right)
            if (allocated(recv_left))  deallocate(recv_left)
            if (allocated(recv_right)) deallocate(recv_right)

#else /* FFTW */

            if ( all(cg%bnd(xdim, LO:HI) == BND_SHE) ) then

               if (allocated(send_right)) deallocate(send_right)
               if (.not.allocated(send_right)) allocate(send_right(flind%all, cg%nb, cg%nyb, cg%nz))

               if (allocated(send_left)) deallocate(send_left)
               if (.not.allocated(send_left)) allocate(send_left(flind%all, cg%nb, cg%nyb, cg%nz))

               if (allocated(recv_left)) deallocate(recv_left)
               if (.not.allocated(recv_left)) allocate(recv_left(flind%all, cg%nb, cg%nyb, cg%nz))

               if (allocated(recv_right)) deallocate(recv_right)
               if (.not.allocated(recv_right)) allocate(recv_right(flind%all, cg%nb, cg%nyb, cg%nz))

               do i = lbound(cg%u%arr,1), ubound(cg%u%arr,1)
                  send_left( i,1:cg%nb,:,:) = unshear_fft(cg%u%arr(i, cg%is:cg%isb, cg%js:cg%je,:), cg%x(cg%is:cg%isb),dely,.true.)
                  send_right(i,1:cg%nb,:,:) = unshear_fft(cg%u%arr(i, cg%ieb:cg%ie, cg%js:cg%je,:), cg%x(cg%ieb:cg%ie),dely,.true.)
               enddo

               call MPI_Isend(send_left , flind%all*cg%nyb*cg%nz*cg%nb, MPI_DOUBLE_PRECISION, procn(dir,LO), 10, comm, req(1), ierr)
               call MPI_Isend(send_right, flind%all*cg%nyb*cg%nz*cg%nb, MPI_DOUBLE_PRECISION, procn(dir,HI), 20, comm, req(3), ierr)
               call MPI_Irecv(recv_left , flind%all*cg%nyb*cg%nz*cg%nb, MPI_DOUBLE_PRECISION, procn(dir,LO), 20, comm, req(2), ierr)
               call MPI_Irecv(recv_right, flind%all*cg%nyb*cg%nz*cg%nb, MPI_DOUBLE_PRECISION, procn(dir,HI), 10, comm, req(4), ierr)

               call MPI_Waitall(4,req(:),status(:,:),ierr)

               do i = lbound(cg%u%arr,1), ubound(cg%u%arr,1)
                  cg%u%arr(i,1:cg%nb,        cg%js:cg%je,:) = unshear_fft(recv_left (i,1:cg%nb,:,:), cg%x(1:cg%nb),dely)
                  cg%u%arr(i, cg%ie+1:cg%nx, cg%js:cg%je,:) = unshear_fft(recv_right(i,1:cg%nb,:,:), cg%x(cg%ie+1:cg%nx),dely)
               enddo

               if (allocated(send_left))  deallocate(send_left)
               if (allocated(send_right)) deallocate(send_right)
               if (allocated(recv_left))  deallocate(recv_left)
               if (allocated(recv_right)) deallocate(recv_right)

            endif
#endif /* FFTW */
         else
#endif /* SHEAR_BND */
            if (psize(dir) > 1) then

               jtag = 20*dir
               itag = jtag - 10
               call MPI_Isend(cg%u%arr(1,1,1,1), 1, cg%mbc(FLUID, dir, LO, BLK), procn(dir,LO), itag, comm3d, req(1), ierr)
               call MPI_Isend(cg%u%arr(1,1,1,1), 1, cg%mbc(FLUID, dir, HI, BLK), procn(dir,HI), jtag, comm3d, req(2), ierr)
               call MPI_Irecv(cg%u%arr(1,1,1,1), 1, cg%mbc(FLUID, dir, LO, BND), procn(dir,LO), jtag, comm3d, req(3), ierr)
               call MPI_Irecv(cg%u%arr(1,1,1,1), 1, cg%mbc(FLUID, dir, HI, BND), procn(dir,HI), itag, comm3d, req(4), ierr)

               call MPI_Waitall(4,req(:),status(:,:),ierr)
            endif
#ifdef SHEAR_BND
         endif
#endif /* SHEAR_BND */

! MPI + non-MPI corner-periodic boundary condition

         if (cg%bnd(xdim, LO) == BND_COR) then
!   - lower to left
            if (pcoords(xdim) == 0 .and. pcoords(ydim) == 0) then
               do i=1, cg%nb
                  do j=cg%js, cg%ny
                     cg%u%arr(iarr_all_dn,i,j,:) =  cg%u%arr(iarr_all_dn,j,cg%isb+1-i,:)
                     cg%u%arr(iarr_all_mx,i,j,:) = -cg%u%arr(iarr_all_my,j,cg%isb+1-i,:)
                     cg%u%arr(iarr_all_my,i,j,:) =  cg%u%arr(iarr_all_mx,j,cg%isb+1-i,:)
                     cg%u%arr(iarr_all_mz,i,j,:) =  cg%u%arr(iarr_all_mz,j,cg%isb+1-i,:)
#ifndef ISO
                     cg%u%arr(iarr_all_en,i,j,:) =  cg%u%arr(iarr_all_en,j,cg%isb+1-i,:)
#endif /* !ISO */
#ifdef COSM_RAYS
                     cg%u%arr(iarr_all_crs,i,j,:) =  cg%u%arr(iarr_all_crs,j,cg%isb+1-i,:)
#endif /* COSM_RAYS */
                  enddo
               enddo
            endif

            if (procxyl > 0) then
               allocate(send_left(flind%all, cg%nb, cg%ny, cg%nz), recv_left(flind%all, cg%nx, cg%nb, cg%nz))

               send_left(:,:,:,:) = cg%u%arr(:, cg%is:cg%isb,:,:)

               call MPI_Isend(send_left, flind%all*cg%nb*cg%ny*cg%nz, MPI_DOUBLE_PRECISION, procxyl, 70, comm, req(1), ierr)
               call MPI_Irecv(recv_left, flind%all*cg%nx*cg%nb*cg%nz, MPI_DOUBLE_PRECISION, procxyl, 80, comm, req(2), ierr)

               call MPI_Waitall(2,req(:),status(:,:),ierr)

               do i=1, cg%nb
                  do j=1, cg%ny
                     cg%u%arr(iarr_all_dn,i,j,:) =  recv_left(iarr_all_dn,j, cg%is-i,:)
                     cg%u%arr(iarr_all_mx,i,j,:) = -recv_left(iarr_all_my,j, cg%is-i,:)
                     cg%u%arr(iarr_all_my,i,j,:) =  recv_left(iarr_all_mx,j, cg%is-i,:)
                     cg%u%arr(iarr_all_mz,i,j,:) =  recv_left(iarr_all_mz,j, cg%is-i,:)
#ifndef ISO
                     cg%u%arr(iarr_all_en,i,j,:) =  recv_left(iarr_all_en,j, cg%is-i,:)
#endif /* !ISO */
#ifdef COSM_RAYS
                     cg%u%arr(iarr_all_crs,i,j,:) =  recv_left(iarr_all_crs,j, cg%is-i,:)
#endif /* COSM_RAYS */
                  enddo
               enddo

               if (allocated(send_left))  deallocate(send_left)
               if (allocated(recv_left))  deallocate(recv_left)
            endif
         endif

         if (cg%bnd(ydim, LO) == BND_COR) then
!   - left to lower
            if (pcoords(ydim) == 0 .and. pcoords(xdim) == 0 ) then
               do j=1, cg%nb
                  do i=cg%is, cg%nx
                     cg%u%arr(iarr_all_dn,i,j,:) =  cg%u%arr(iarr_all_dn,cg%isb+1-j,i,:)
                     cg%u%arr(iarr_all_mx,i,j,:) =  cg%u%arr(iarr_all_my,cg%isb+1-j,i,:)
                     cg%u%arr(iarr_all_my,i,j,:) = -cg%u%arr(iarr_all_mx,cg%isb+1-j,i,:)
                     cg%u%arr(iarr_all_mz,i,j,:) =  cg%u%arr(iarr_all_mz,cg%isb+1-j,i,:)
#ifndef ISO
                     cg%u%arr(iarr_all_en,i,j,:) =  cg%u%arr(iarr_all_en,cg%isb+1-j,i,:)
#endif /* !ISO */
#ifdef COSM_RAYS
                     cg%u%arr(iarr_all_crs,i,j,:) =  cg%u%arr(iarr_all_crs,cg%isb+1-j,i,:)
#endif /* COSM_RAYS */
                  enddo
               enddo
!   - interior to corner
               do j=1, cg%nb
                  do i=1, cg%nb
                     cg%u%arr(iarr_all_dn,i,j,:) =   cg%u%arr(iarr_all_dn,cg%isb+1-i,cg%jsb+1-j,:)
                     cg%u%arr(iarr_all_mx,i,j,:) =  -cg%u%arr(iarr_all_mx,cg%isb+1-i,cg%jsb+1-j,:)
                     cg%u%arr(iarr_all_my,i,j,:) =  -cg%u%arr(iarr_all_my,cg%isb+1-i,cg%jsb+1-j,:)
                     cg%u%arr(iarr_all_mz,i,j,:) =   cg%u%arr(iarr_all_mz,cg%isb+1-i,cg%jsb+1-j,:)
#ifndef ISO
                     cg%u%arr(iarr_all_en,i,j,:) =   cg%u%arr(iarr_all_en,cg%isb+1-i,cg%jsb+1-j,:)
#endif /* !ISO */
#ifdef COSM_RAYS
                     cg%u%arr(iarr_all_crs,i,j,:) =   cg%u%arr(iarr_all_crs,cg%isb+1-i,cg%jsb+1-j,:)
#endif /* COSM_RAYS */
                  enddo
               enddo
            endif

            if (procyxl > 0) then
               allocate(send_left(flind%all, cg%nx, cg%nb, cg%nz), recv_left(flind%all, cg%nb, cg%ny, cg%nz))

               send_left(:,:,:,:) = cg%u%arr(:,:, cg%js:cg%jsb,:)

               call MPI_Isend(send_left, flind%all*cg%nx*cg%nb*cg%nz, MPI_DOUBLE_PRECISION, procyxl, 80, comm, req(1), ierr)
               call MPI_Irecv(recv_left, flind%all*cg%nb*cg%ny*cg%nz, MPI_DOUBLE_PRECISION, procyxl, 70, comm, req(2), ierr)

               call MPI_Waitall(2,req(:),status(:,:),ierr)

               do j=1, cg%nb
                  do i=1, cg%nx
                     cg%u%arr(iarr_all_dn,i,j,:) =  recv_left(iarr_all_dn, cg%js-j,i,:)
                     cg%u%arr(iarr_all_mx,i,j,:) =  recv_left(iarr_all_my, cg%js-j,i,:)
                     cg%u%arr(iarr_all_my,i,j,:) = -recv_left(iarr_all_mx, cg%js-j,i,:)
                     cg%u%arr(iarr_all_mz,i,j,:) =  recv_left(iarr_all_mz, cg%js-j,i,:)
#ifndef ISO
                     cg%u%arr(iarr_all_en,i,j,:) =  recv_left(iarr_all_en, cg%js-j,i,:)
#endif /* !ISO */
#ifdef COSM_RAYS
                     cg%u%arr(iarr_all_crs,i,j,:) =  recv_left(iarr_all_crs, cg%js-j,i,:)
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

      select case (dir)
      case (xdim)
         call func_bnd_xl(cg)
         call func_bnd_xr(cg)
      case (ydim)

         select case (cg%bnd(ydim, LO))
         case (BND_COR, BND_INF, BND_MPI)
            ! Do nothing
         case (BND_PER)
             if (comm3d /= MPI_COMM_NULL) cg%u%arr(:,:,1:cg%nb,:)                         = cg%u%arr(:,:, cg%jeb:cg%je,:)
         case (BND_USER)
            call user_bnd_yl(cg)
         case (BND_REF)
            do ib=1, cg%nb

               cg%u%arr((/iarr_all_dn,iarr_all_mx,iarr_all_mz/),:, cg%js-ib,:)   = cg%u%arr((/iarr_all_dn,iarr_all_mx,iarr_all_mz/),:, cg%nb+ib,:)
               cg%u%arr(iarr_all_my,:, cg%js-ib,:)                 =-cg%u%arr(iarr_all_my,:, cg%nb+ib,:)
#ifndef ISO
               cg%u%arr(iarr_all_en,:, cg%js-ib,:)                 = cg%u%arr(iarr_all_en,:, cg%nb+ib,:)
#endif /* !ISO */
#ifdef COSM_RAYS
               cg%u%arr(iarr_all_crs,:, cg%js-ib,:)                 = cg%u%arr(iarr_all_crs,:, cg%nb+ib,:)
#endif /* COSM_RAYS */
            enddo
         case (BND_OUT)
            do ib=1, cg%nb
               cg%u%arr(:,:,ib,:)                         = cg%u%arr(:,:, cg%js,:)
#ifdef COSM_RAYS
               cg%u%arr(iarr_all_crs,:,ib,:)                      = smallecr
#endif /* COSM_RAYS */
            enddo
         case (BND_OUTD)
            do ib=1, cg%nb

               cg%u%arr((/iarr_all_dn,iarr_all_mx,iarr_all_mz/),:,ib,:)        = cg%u%arr((/iarr_all_dn,iarr_all_mx,iarr_all_mz/),:, cg%js,:)
               cg%u%arr(iarr_all_my,:,ib,:)                      = min(cg%u%arr(iarr_all_my,:, cg%js,:),0.0)
#ifndef ISO
               cg%u%arr(iarr_all_en,:,ib,:)                      = cg%u%arr(iarr_all_en,:, cg%js,:)
#endif /* !ISO */
#ifdef COSM_RAYS
               cg%u%arr(iarr_all_crs,:,ib,:)                      = smallecr
#endif /* COSM_RAYS */
            enddo
         case default
            write(msg,'("[fluid_boundaries:bnd_u]: Left boundary condition ",i3," not implemented in Y-direction")') cg%bnd(ydim, LO)
            call warn(msg)
         end select  ! (cg%bnd(ydim, LO))

         select case (cg%bnd(ydim, HI))
         case (BND_COR, BND_INF, BND_MPI)
            ! Do nothing
         case (BND_PER)
            if (comm3d /= MPI_COMM_NULL) cg%u%arr(:,:, cg%je+1:cg%ny,:)            = cg%u%arr(:,:, cg%js:cg%jsb,:)
         case (BND_USER)
            call user_bnd_yr(cg)
         case (BND_REF)
            do ib=1, cg%nb

               cg%u%arr((/iarr_all_dn,iarr_all_mx,iarr_all_mz/),:, cg%je+ib,:) = cg%u%arr((/iarr_all_dn,iarr_all_mx,iarr_all_mz/),:, cg%je+1-ib,:)
               cg%u%arr(iarr_all_my,:, cg%je+ib,:)               =-cg%u%arr(iarr_all_my,:, cg%je+1-ib,:)
#ifndef ISO
               cg%u%arr(iarr_all_en,:, cg%je+ib,:)               = cg%u%arr(iarr_all_en,:, cg%je+1-ib,:)
#endif /* !ISO */
#ifdef COSM_RAYS
               cg%u%arr(iarr_all_crs,:, cg%je+ib,:)               = cg%u%arr(iarr_all_crs,:, cg%je+1-ib,:)
#endif /* COSM_RAYS */
            enddo
         case (BND_OUT)
            do ib=1, cg%nb
               cg%u%arr(:,:, cg%je+ib,:)                  = cg%u%arr(:,:, cg%je,:)
#ifdef COSM_RAYS
               cg%u%arr(iarr_all_crs,:, cg%je+ib,:)               = smallecr
#endif /* COSM_RAYS */
            enddo
         case (BND_OUTD)
            do ib=1, cg%nb

               cg%u%arr((/iarr_all_dn,iarr_all_mx,iarr_all_mz/),:, cg%je+ib,:) = cg%u%arr((/iarr_all_dn,iarr_all_mx,iarr_all_mz/),:, cg%je,:)
               cg%u%arr(iarr_all_my,:, cg%je+ib,:)               = max(cg%u%arr(iarr_all_my,:, cg%je,:),0.0)
#ifndef ISO
               cg%u%arr(iarr_all_en,:, cg%je+ib,:)               = cg%u%arr(iarr_all_en,:, cg%je,:)
#endif /* !ISO */
#ifdef COSM_RAYS
               cg%u%arr(iarr_all_crs,:, cg%je+ib,:)               = smallecr
#endif /* COSM_RAYS */
            enddo
         case default
            write(msg,'("[fluid_boundaries:bnd_u]: Right boundary condition ",i3," not implemented in Y-direction")') cg%bnd(ydim, HI)
            call warn(msg)
         end select  ! (cg%bnd(ydim, HI))

      case (zdim)

         select case (cg%bnd(zdim, LO))
         case (BND_MPI)
            ! Do nothing if mpi
         case (BND_USER)
            call user_bnd_zl(cg)
         case (BND_PER)
            if (comm3d /= MPI_COMM_NULL) cg%u%arr(:,:,:,1:cg%nb)                         = cg%u%arr(:,:,:, cg%keb:cg%ke)
         case (BND_REF)
            do ib=1, cg%nb

               cg%u%arr((/iarr_all_dn,iarr_all_mx,iarr_all_my/),:,:, cg%ks-ib)   = cg%u%arr((/iarr_all_dn,iarr_all_mx,iarr_all_my/),:,:, cg%nb+ib)
               cg%u%arr(iarr_all_mz,:,:, cg%ks-ib)                 =-cg%u%arr(iarr_all_mz,:,:, cg%nb+ib)
#ifndef ISO
               cg%u%arr(iarr_all_en,:,:, cg%ks-ib)                 = cg%u%arr(iarr_all_en,:,:, cg%nb+ib)
#endif /* !ISO */
#ifdef COSM_RAYS
               cg%u%arr(iarr_all_crs,:,:, cg%ks-ib)                 = cg%u%arr(iarr_all_crs,:,:, cg%nb+ib)
#endif /* COSM_RAYS */
            enddo
         case (BND_OUT)
            do ib=1, cg%nb
               cg%u%arr(:,:,:,ib)                         = cg%u%arr(:,:,:, cg%ks)
#ifdef COSM_RAYS
               cg%u%arr(iarr_all_crs,:,:,ib)                      = smallecr
#endif /* COSM_RAYS */
            enddo
         case (BND_OUTD)
            do ib=1, cg%nb

               cg%u%arr((/iarr_all_dn,iarr_all_mx,iarr_all_my/),:,:,ib)        = cg%u%arr((/iarr_all_dn,iarr_all_mx,iarr_all_my/),:,:, cg%ks)
!> \deprecated BEWARE: use of uninitialized value on first call (a side effect of r1726)
               cg%u%arr(iarr_all_mz,:,:,ib)                      = min(cg%u%arr(iarr_all_mz,:,:, cg%ks),0.0)
#ifndef ISO
               cg%u%arr(iarr_all_en,:,:,ib)                      = cg%u%arr(iarr_all_en,:,:, cg%ks)
#endif /* !ISO */
#ifdef COSM_RAYS
               cg%u%arr(iarr_all_crs,:,:,ib)                      = smallecr
#endif /* COSM_RAYS */
            enddo
#ifdef GRAV
         case (BND_OUTH)
            do ib=1, cg%nb
               kb = cg%ks-ib
               call outh_bnd(kb+1, kb, "min")
            enddo ! ib
#endif /* GRAV */
         case default
            write(msg,'("[fluid_boundaries:bnd_u]: Left boundary condition ",i3," not implemented in Z-direction")') cg%bnd(zdim, LO)
            call warn(msg)
         end select  ! (cg%bnd(zdim, LO))

         select case (cg%bnd(zdim, HI))
         case (BND_MPI)
            ! Do nothing if mpi
         case (BND_USER)
            call user_bnd_zr(cg)
         case (BND_PER)
            if (comm3d /= MPI_COMM_NULL) cg%u%arr(:,:,:, cg%ke+1:cg%nz)            = cg%u%arr(:,:,:, cg%ks:cg%ksb)
         case (BND_REF)
            do ib=1, cg%nb

               cg%u%arr((/iarr_all_dn,iarr_all_mx,iarr_all_my/),:,:, cg%ke+ib) = cg%u%arr((/iarr_all_dn,iarr_all_mx,iarr_all_my/),:,:, cg%ke+1-ib)
               cg%u%arr(iarr_all_mz,:,:, cg%ke+ib)               =-cg%u%arr(iarr_all_mz,:,:, cg%ke+1-ib)
#ifndef ISO
               cg%u%arr(iarr_all_en,:,:, cg%ke+ib)               = cg%u%arr(iarr_all_en,:,:, cg%ke+1-ib)
#endif /* !ISO */
#ifdef COSM_RAYS
               cg%u%arr(iarr_all_crs,:,:, cg%ke+ib)               = cg%u%arr(iarr_all_crs,:,:, cg%ke+1-ib)
#endif /* COSM_RAYS */
            enddo
         case (BND_OUT)
            do ib=1, cg%nb
               cg%u%arr(:,:,:, cg%ke+ib)                  = cg%u%arr(:,:,:, cg%ke)
#ifdef COSM_RAYS
               cg%u%arr(iarr_all_crs,:,:, cg%ke+ib)               = smallecr
#endif /* COSM_RAYS */
            enddo
         case (BND_OUTD)
            do ib=1, cg%nb

               cg%u%arr((/iarr_all_dn,iarr_all_mx,iarr_all_my/),:,:, cg%ke+ib) = cg%u%arr((/iarr_all_dn,iarr_all_mx,iarr_all_my/),:,:, cg%ke)
!> \deprecated BEWARE: use of uninitialized value on first call (a side effect of r1726)
               cg%u%arr(iarr_all_mz,:,:, cg%ke+ib)               = max(cg%u%arr(iarr_all_mz,:,:, cg%ke),0.0)
#ifndef ISO
               cg%u%arr(iarr_all_en,:,:, cg%ke+ib)               = cg%u%arr(iarr_all_en,:,:, cg%ke)
#endif /* !ISO */
#ifdef COSM_RAYS
               cg%u%arr(iarr_all_crs,:,:, cg%ke+ib)               = smallecr
#endif /* COSM_RAYS */
            enddo
#ifdef GRAV
         case (BND_OUTH)
            do ib=1, cg%nb
               kb = cg%ke-1+ib
               call outh_bnd(kb, kb+1, "max")
            enddo ! ib
#endif /* GRAV */
         case default
            write(msg,'("[fluid_boundaries:bnd_u]: Right boundary condition ",i3," not implemented in Z-direction")') cg%bnd(zdim, HI)
            call warn(msg)
         end select  ! (cg%bnd(zdim, HI))

      end select  ! (dim)

   end subroutine bnd_u

   subroutine all_fluid_boundaries

      use constants,  only: xdim, zdim, FLUID
      use dataio_pub, only: die
      use grid,       only: cga
      use grid_cont,  only: cg_list_element
      use mpi,        only: MPI_COMM_NULL
      use mpisetup,   only: has_dir, comm3d

      implicit none

      type(cg_list_element), pointer :: cgl
      integer :: dir

      if (ubound(cga%cg_all(:), dim=1) > 1) call die("[fluidboundaries:all_fluid_boundaries] multiple grid pieces per procesor not implemented yet") !nontrivial communication

      cgl => cga%cg_leafs%cg_l(1)
      do while (associated(cgl))

         if (comm3d == MPI_COMM_NULL) then
            call cgl%cg%internal_boundaries(FLUID, pa4d=cgl%cg%u%arr)
         endif

         do dir = xdim, zdim
            if (has_dir(dir)) call bnd_u(dir, cgl%cg)
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine all_fluid_boundaries

end module fluidboundaries
