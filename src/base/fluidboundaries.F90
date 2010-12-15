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


      subroutine init_fluidboundaries
         use mpisetup,              only: bnd_xl, bnd_xr
         use dataio_pub,            only: msg, warn
         use fluidboundaries_funcs, only: bnd_null, bnd_xl_per, bnd_xl_ref, bnd_xl_out, bnd_xl_outd, bnd_xr_per, bnd_xr_ref, bnd_xr_out, bnd_xr_outd
         use fluidboundaries_pub,   only: user_bnd_xl, user_bnd_xr, func_bnd_xl, func_bnd_xr
         implicit none

         select case (bnd_xl)
            case ('cor', 'inf', 'mpi', 'she', 'shef')
               func_bnd_xl => bnd_null
            case ('per')
               func_bnd_xl => bnd_xl_per
            case ('user')
               func_bnd_xl => user_bnd_xl
            case ('ref')
               func_bnd_xl => bnd_xl_ref
            case ('out')
               func_bnd_xl => bnd_xl_out
            case ('outd')
               func_bnd_xl => bnd_xl_outd
            case default
               func_bnd_xl => bnd_null
               write(msg,'("[fluid_boundaries:bnd_u]: Boundary condition ",a," not implemented in ",a)') trim(bnd_xl), "xdim"
               call warn(msg)
         end select  ! (bnd_xl)

         select case (bnd_xr)
            case ('cor', 'inf', 'mpi', 'she', 'shef')
               func_bnd_xr => bnd_null
            case ('per')
               func_bnd_xr => bnd_xr_per
            case ('user')
               func_bnd_xr => user_bnd_xr
            case ('ref')
               func_bnd_xr => bnd_xr_ref
            case ('out')
               func_bnd_xr => bnd_xr_out
            case ('outd')
               func_bnd_xr => bnd_xr_outd
            case default
               func_bnd_xr => bnd_null
               write(msg,'("[fluid_boundaries:bnd_u]: Boundary condition ",a," not implemented in ",a)') trim(bnd_xr), "xdim"
               call warn(msg)
         end select  ! (bnd_xr)
      end subroutine init_fluidboundaries

   subroutine bnd_u(dim)
      use arrays,              only: u, b
      use dataio_pub,          only: msg, warn
      use fluidboundaries_pub, only: user_bnd_yl, user_bnd_yr, user_bnd_zl, user_bnd_zr, func_bnd_xl, func_bnd_xr
      use fluidindex,          only: nvar, iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
      use grid,                only: nb, nyb, x, y, z, nzb, nyb, nxb, nx, ny, nz
      use mpisetup,            only: ierr, MPI_XY_RIGHT_DOM, MPI_XY_RIGHT_BND, MPI_XY_LEFT_DOM, MPI_XY_LEFT_BND, &
                                     MPI_XZ_RIGHT_DOM, MPI_XZ_RIGHT_BND, MPI_XZ_LEFT_DOM, MPI_XZ_LEFT_BND, &
                                     MPI_YZ_RIGHT_DOM, MPI_YZ_RIGHT_BND, MPI_YZ_LEFT_DOM, MPI_YZ_LEFT_BND, &
                                     pxsize, pysize, pzsize, proczl, proczr, procyl, procyr, procxl, procxr, &
                                     pcoords, bnd_xr, bnd_xl, bnd_yl, bnd_yr, bnd_zl, bnd_zr, req, status, comm, comm3d, &
                                     procxyl, procyxl, smalld
      use mpi,                 only: MPI_DOUBLE_PRECISION
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

      character(len=*) :: dim
      logical, save    :: frun = .true.
      integer :: i,j, ib
      real, allocatable :: send_left(:,:,:,:),recv_left(:,:,:,:)
#ifdef GRAV
      integer          :: kb
#endif /* GRAV */
#ifdef SHEAR_BND
      real, allocatable :: send_right(:,:,:,:),recv_right(:,:,:,:)
#ifdef FFTW
      real, allocatable, dimension(:,:,:) :: temp
#endif /* FFTW */
#endif /* SHEAR_BND */

      if (frun) then
         call init_fluidboundaries
         frun = .false.
      endif

! MPI block communication
      select case (dim)
      case ('xdim')
#ifdef SHEAR_BND
#ifndef FFTW
         allocate(send_right(nvar%all,nb,ny,nz), send_left(nvar%all,nb,ny,nz), &
                  recv_left(nvar%all,nb,ny,nz), recv_right(nvar%all,nb,ny,nz) )
         send_left(:,:,:,:)          =  u(:,nb+1:2*nb,:,:)
         send_right(:,:,:,:)         =  u(:,nxb+1:nxb+nb,:,:)
!
! odejmujemy ped_y i energie odpowiadajace niezaburzonej rozniczkowej rotacji na lewym brzegu
!
         if (bnd_xl == 'she') then
            do i=1,nb
               send_left (iarr_all_my,i,:,:) = send_left(iarr_all_my,i,:,:) &
                                         +qshear*omega * x(nb+i)     * send_left(iarr_all_dn,i,:,:)
#ifndef ISO
               send_left (iarr_all_en,i,:,:) = send_left(iarr_all_en,i,:,:) &
                                    -0.5*(qshear*omega * x(nb+i))**2 * send_left(iarr_all_dn,i,:,:)
#endif /* !ISO */
            enddo
!
! przesuwamy o calkowita liczbe komorek + periodyczny wb w kierunku y
!
         if (has_dir(ydim)) then
            send_left (:,:,nb+1:nb+nyb,:)        = cshift(send_left (:,:,nb+1:nb+nyb,:),dim=3,shift= delj)
            send_left (:,:,1:nb,:)               = send_left (:,:,nyb+1:nyb+nb,:)
            send_left (:,:,nb+nyb+1:nyb+2*nb,:)  = send_left (:,:,nb+1:2*nb,:)
         endif
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
         if (bnd_xr == 'she') then
            do i=1,nb
               send_right(iarr_all_my,i,:,:) = send_right(iarr_all_my,i,:,:) &
                                         +qshear*omega * x(nxb+i)     * send_right(iarr_all_dn,i,:,:)
#ifndef ISO
               send_right(iarr_all_en,i,:,:) = send_right(iarr_all_en,i,:,:) &
                                    -0.5*(qshear*omega * x(nxb+i))**2 * send_right(iarr_all_dn,i,:,:)

#endif /* !ISO */
            enddo
!
! przesuwamy o calkowita liczbe komorek + periodyczny wb w kierunku y
!
         if (has_dir(ydim)) then
            send_right(:,:,nb+1:nb+nyb,:)        = cshift(send_right(:,:,nb+1:nb+nyb,:),dim=3,shift=-delj)
            send_right (:,:,1:nb,:)              = send_right(:,:,nyb+1:nyb+nb,:)
            send_right (:,:,nb+nyb+1:nyb+2*nb,:) = send_right(:,:,nb+1:2*nb,:)
         endif
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
         CALL MPI_Isend   (send_left , nvar%all*ny*nz*nb, MPI_DOUBLE_PRECISION, procxl, 10, comm, req(1), ierr)
         CALL MPI_Isend   (send_right, nvar%all*ny*nz*nb, MPI_DOUBLE_PRECISION, procxr, 20, comm, req(3), ierr)
         CALL MPI_Irecv   (recv_left , nvar%all*ny*nz*nb, MPI_DOUBLE_PRECISION, procxl, 20, comm, req(2), ierr)
         CALL MPI_Irecv   (recv_right, nvar%all*ny*nz*nb, MPI_DOUBLE_PRECISION, procxr, 10, comm, req(4), ierr)

         call MPI_Waitall(4,req(:),status(:,:),ierr)

!
! dodajemy ped_y i energie odpowiadajace niezaburzonej rozniczkowej rotacji na prawym brzegu
!
         if (bnd_xr == 'she') then
            do i=1,nb
#ifndef ISO
               recv_right (iarr_all_en,i,:,:) = recv_right (iarr_all_en,i,:,:) &
                                      +0.5*(qshear*omega * x(nb+nxb+i))**2 * recv_right(iarr_all_dn,i,:,:)
#endif /* !ISO */
               recv_right (iarr_all_my,i,:,:) = recv_right (iarr_all_my,i,:,:) &
                                           -qshear*omega * x(nb+nxb+i)     * recv_right(iarr_all_dn,i,:,:)
            enddo
         endif !(bnd_xr == 'she')
!
! dodajemy ped_y i energie odpowiadajace niezaburzonej rozniczkowej rotacji na lewym brzegu
!
         if (bnd_xl == 'she') then
            do i=1,nb
#ifndef ISO
               recv_left(iarr_all_en,i,:,:) = recv_left(iarr_all_en,i,:,:) &
                                    +0.5*(qshear*omega * x(i))**2 * recv_left(iarr_all_dn,i,:,:)
#endif /* !ISO */
               recv_left(iarr_all_my,i,:,:) = recv_left(iarr_all_my,i,:,:) &
                                         -qshear*omega * x(i)     * recv_left(iarr_all_dn,i,:,:)
            enddo
         endif !(bnd_xl == 'she')

         u(:,1:nb,:,:)              = recv_left(:,1:nb,:,:)
         u(:,nxb+nb+1:nxb+2*nb,:,:) = recv_right(:,1:nb,:,:)

         !!!! BEWARE: smalld is called only for the first fluid
         u(iarr_all_dn(1),1:nb,:,:)              = max(u(iarr_all_dn(1),1:nb,:,:),smalld)
         u(iarr_all_dn(1),nxb+nb+1:nxb+2*nb,:,:) = max(u(iarr_all_dn(1),nxb+nb+1:nxb+2*nb,:,:),smalld)
         if (allocated(send_left))  deallocate(send_left)
         if (allocated(send_right)) deallocate(send_right)
         if (allocated(recv_left))  deallocate(recv_left)
         if (allocated(recv_right)) deallocate(recv_right)

#else /* FFTW */

         if ( (bnd_xl == 'she').and.(bnd_xr == 'she')) then

         if (allocated(send_right)) deallocate(send_right)
         if (.not.allocated(send_right)) allocate(send_right(nvar%all,nb,nyb,nz))

         if (allocated(send_left)) deallocate(send_left)
         if (.not.allocated(send_left)) allocate(send_left(nvar%all,nb,nyb,nz))

         if (allocated(recv_left)) deallocate(recv_left)
         if (.not.allocated(recv_left)) allocate(recv_left(nvar%all,nb,nyb,nz))

         if (allocated(recv_right)) deallocate(recv_right)
         if (.not.allocated(recv_right)) allocate(recv_right(nvar%all,nb,nyb,nz))

            do i = LBOUND(u,1), UBOUND(u,1)
               send_left(i,1:nb,:,:)   = unshear_fft(u(i,nb+1:2*nb,nb+1:ny-nb,:),x(nb+1:2*nb),dely,.true.)
               send_right(i,1:nb,:,:)  = unshear_fft(u(i,nx-2*nb+1:nx-nb,nb+1:ny-nb,:),x(nx-2*nb+1:nx-nb),dely,.true.)
            enddo

            CALL MPI_Isend   (send_left , nvar%all*nyb*nz*nb, MPI_DOUBLE_PRECISION, procxl, 10, comm, req(1), ierr)
            CALL MPI_Isend   (send_right, nvar%all*nyb*nz*nb, MPI_DOUBLE_PRECISION, procxr, 20, comm, req(3), ierr)
            CALL MPI_Irecv   (recv_left , nvar%all*nyb*nz*nb, MPI_DOUBLE_PRECISION, procxl, 20, comm, req(2), ierr)
            CALL MPI_Irecv   (recv_right, nvar%all*nyb*nz*nb, MPI_DOUBLE_PRECISION, procxr, 10, comm, req(4), ierr)

            call MPI_Waitall(4,req(:),status(:,:),ierr)

            do i = LBOUND(u,1), UBOUND(u,1)
               u(i,1:nb,nb+1:ny-nb,:) = unshear_fft(recv_left(i,1:nb,:,:), x(1:nb),dely)
               u(i,nx-nb+1:nx,nb+1:ny-nb,:) = unshear_fft(recv_right(i,1:nb,:,:),x(nx-nb+1:nx),dely)
            enddo

         if (allocated(send_left))  deallocate(send_left)
         if (allocated(send_right)) deallocate(send_right)
         if (allocated(recv_left))  deallocate(recv_left)
         if (allocated(recv_right)) deallocate(recv_right)

         endif
#endif /* FFTW */
#else /* !SHEAR_BND */
         if (pxsize .gt. 1) then

            CALL MPI_Isend   (u(1,1,1,1), 1, MPI_YZ_LEFT_DOM,  procxl, 10, comm3d, req(1), ierr)
            CALL MPI_Isend   (u(1,1,1,1), 1, MPI_YZ_RIGHT_DOM, procxr, 20, comm3d, req(3), ierr)
            CALL MPI_Irecv   (u(1,1,1,1), 1, MPI_YZ_LEFT_BND,  procxl, 20, comm3d, req(2), ierr)
            CALL MPI_Irecv   (u(1,1,1,1), 1, MPI_YZ_RIGHT_BND, procxr, 10, comm3d, req(4), ierr)

            call MPI_Waitall(4,req(:),status(:,:),ierr)
         endif
#endif /* !SHEAR_BND */
      case ('ydim')
         if (pysize .gt. 1) then

            CALL MPI_Isend   (u(1,1,1,1), 1, MPI_XZ_LEFT_DOM,  procyl, 30, comm3d, req(1), ierr)
            CALL MPI_Isend   (u(1,1,1,1), 1, MPI_XZ_RIGHT_DOM, procyr, 40, comm3d, req(3), ierr)
            CALL MPI_Irecv   (u(1,1,1,1), 1, MPI_XZ_LEFT_BND,  procyl, 40, comm3d, req(2), ierr)
            CALL MPI_Irecv   (u(1,1,1,1), 1, MPI_XZ_RIGHT_BND, procyr, 30, comm3d, req(4), ierr)

            call MPI_Waitall(4,req(:),status(:,:),ierr)
         endif

      case ('zdim')
         if (pzsize .gt. 1) then

            CALL MPI_Isend   (u(1,1,1,1), 1, MPI_XY_LEFT_DOM,  proczl, 50, comm3d, req(1), ierr)
            CALL MPI_Isend   (u(1,1,1,1), 1, MPI_XY_RIGHT_DOM, proczr, 60, comm3d, req(3), ierr)
            CALL MPI_Irecv   (u(1,1,1,1), 1, MPI_XY_LEFT_BND,  proczl, 60, comm3d, req(2), ierr)
            CALL MPI_Irecv   (u(1,1,1,1), 1, MPI_XY_RIGHT_BND, proczr, 50, comm3d, req(4), ierr)

            call MPI_Waitall(4,req(:),status(:,:),ierr)
         endif
      end select ! (dim)

! MPI + non-MPI corner-periodic boundary condition

      if (bnd_xl .eq. 'cor') then
!   - lower to left
         if (pcoords(1) .eq. 0 .and. pcoords(2) .eq. 0) then
            do i=1,nb
               do j=nb+1,ny
                  u(iarr_all_dn,i,j,:) =  u(iarr_all_dn,j,2*nb+1-i,:)
                  u(iarr_all_mx,i,j,:) = -u(iarr_all_my,j,2*nb+1-i,:)
                  u(iarr_all_my,i,j,:) =  u(iarr_all_mx,j,2*nb+1-i,:)
                  u(iarr_all_mz,i,j,:) =  u(iarr_all_mz,j,2*nb+1-i,:)
#ifndef ISO
                  u(iarr_all_en,i,j,:) =  u(iarr_all_en,j,2*nb+1-i,:)
#endif /* !ISO */
#ifdef COSM_RAYS
                  u(iarr_all_crs,i,j,:) =  u(iarr_all_crs,j,2*nb+1-i,:)
#endif /* COSM_RAYS */
               enddo
            enddo
         endif

         if (procxyl .gt. 0) then
            allocate(send_left(nvar%all,nb,ny,nz), recv_left(nvar%all,nx,nb,nz))

            send_left(:,:,:,:) = u(:,nb+1:2*nb,:,:)

            CALL MPI_Isend   (send_left , nvar%all*nb*ny*nz, MPI_DOUBLE_PRECISION, procxyl, 70, comm, req(1), ierr)
            CALL MPI_Irecv   (recv_left , nvar%all*nx*nb*nz, MPI_DOUBLE_PRECISION, procxyl, 80, comm, req(2), ierr)

            call MPI_Waitall(2,req(:),status(:,:),ierr)

            do i=1,nb
               do j=1,ny
                  u(iarr_all_dn,i,j,:) =  recv_left(iarr_all_dn,j,nb+1-i,:)
                  u(iarr_all_mx,i,j,:) = -recv_left(iarr_all_my,j,nb+1-i,:)
                  u(iarr_all_my,i,j,:) =  recv_left(iarr_all_mx,j,nb+1-i,:)
                  u(iarr_all_mz,i,j,:) =  recv_left(iarr_all_mz,j,nb+1-i,:)
#ifndef ISO
                  u(iarr_all_en,i,j,:) =  recv_left(iarr_all_en,j,nb+1-i,:)
#endif /* !ISO */
#ifdef COSM_RAYS
                  u(iarr_all_crs,i,j,:) =  recv_left(iarr_all_crs,j,nb+1-i,:)
#endif /* COSM_RAYS */
               enddo
            enddo

            if (allocated(send_left))  deallocate(send_left)
            if (allocated(recv_left))  deallocate(recv_left)
         endif
      endif

      if (bnd_yl .eq. 'cor') then
!   - left to lower
         if (pcoords(2) .eq. 0 .and. pcoords(1) .eq. 0 ) then
            do j=1,nb
               do i=nb+1,nx
                  u(iarr_all_dn,i,j,:) =  u(iarr_all_dn,2*nb+1-j,i,:)
                  u(iarr_all_mx,i,j,:) =  u(iarr_all_my,2*nb+1-j,i,:)
                  u(iarr_all_my,i,j,:) = -u(iarr_all_mx,2*nb+1-j,i,:)
                  u(iarr_all_mz,i,j,:) =  u(iarr_all_mz,2*nb+1-j,i,:)
#ifndef ISO
                  u(iarr_all_en,i,j,:) =  u(iarr_all_en,2*nb+1-j,i,:)
#endif /* !ISO */
#ifdef COSM_RAYS
                  u(iarr_all_crs,i,j,:) =  u(iarr_all_crs,2*nb+1-j,i,:)
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
#endif /* !ISO */
#ifdef COSM_RAYS
                  u(iarr_all_crs,i,j,:) =   u(iarr_all_crs,2*nb+1-i,2*nb+1-j,:)
#endif /* COSM_RAYS */
               enddo
            enddo
         endif

         if (procyxl .gt. 0) then
            allocate(send_left(nvar%all,nx,nb,nz), recv_left(nvar%all,nb,ny,nz))

            send_left(:,:,:,:) = u(:,:,nb+1:2*nb,:)

            CALL MPI_Isend   (send_left , nvar%all*nx*nb*nz, MPI_DOUBLE_PRECISION, procyxl, 80, comm, req(1), ierr)
            CALL MPI_Irecv   (recv_left , nvar%all*nb*ny*nz, MPI_DOUBLE_PRECISION, procyxl, 70, comm, req(2), ierr)

            call MPI_Waitall(2,req(:),status(:,:),ierr)

            do j=1,nb
               do i=1,nx
                  u(iarr_all_dn,i,j,:) =  recv_left(iarr_all_dn,nb+1-j,i,:)
                  u(iarr_all_mx,i,j,:) =  recv_left(iarr_all_my,nb+1-j,i,:)
                  u(iarr_all_my,i,j,:) = -recv_left(iarr_all_mx,nb+1-j,i,:)
                  u(iarr_all_mz,i,j,:) =  recv_left(iarr_all_mz,nb+1-j,i,:)
#ifndef ISO
                  u(iarr_all_en,i,j,:) =  recv_left(iarr_all_en,nb+1-j,i,:)
#endif /* !ISO */
#ifdef COSM_RAYS
                  u(iarr_all_crs,i,j,:) =  recv_left(iarr_all_crs,nb+1-j,i,:)
#endif /* COSM_RAYS */
               enddo
            enddo

            if (allocated(send_left))  deallocate(send_left)
            if (allocated(recv_left))  deallocate(recv_left)
         endif
      endif

!===============================================================

! Non-MPI boundary conditions

      select case (dim)
      case ('xdim')
         call func_bnd_xl
         call func_bnd_xr
      case ('ydim')

         select case (bnd_yl)
         case ('cor', 'inf', 'mpi')
!         Do nothing
         case ('per')
            u(:,:,1:nb,:)                         = u(:,:,nyb+1:nyb+nb,:)
         case ('user')
            call user_bnd_yl
         case ('ref')
            do ib=1,nb

               u((/iarr_all_dn,iarr_all_mx,iarr_all_mz/),:,nb+1-ib,:)   = u((/iarr_all_dn,iarr_all_mx,iarr_all_mz/),:,nb+ib,:)
               u(iarr_all_my,:,nb+1-ib,:)                 =-u(iarr_all_my,:,nb+ib,:)
#ifndef ISO
               u(iarr_all_en,:,nb+1-ib,:)                 = u(iarr_all_en,:,nb+ib,:)
#endif /* !ISO */
#ifdef COSM_RAYS
               u(iarr_all_crs,:,nb+1-ib,:)                 = u(iarr_all_crs,:,nb+ib,:)
#endif /* COSM_RAYS */
            enddo
         case ('out')
            do ib=1,nb
               u(:,:,ib,:)                         = u(:,:,nb+1,:)
#ifdef COSM_RAYS
               u(iarr_all_crs,:,ib,:)                      = smallecr
#endif /* COSM_RAYS */
            enddo
         case ('outd')
            do ib=1,nb

               u((/iarr_all_dn,iarr_all_mx,iarr_all_mz/),:,ib,:)        = u((/iarr_all_dn,iarr_all_mx,iarr_all_mz/),:,nb+1,:)
               u(iarr_all_my,:,ib,:)                      = min(u(iarr_all_my,:,nb+1,:),0.0)
#ifndef ISO
               u(iarr_all_en,:,ib,:)                      = u(iarr_all_en,:,nb+1,:)
#endif /* !ISO */
#ifdef COSM_RAYS
               u(iarr_all_crs,:,ib,:)                      = smallecr
#endif /* COSM_RAYS */
            enddo
         case default
            write(msg,'("[fluid_boundaries:bnd_u]: Boundary condition ",a," not implemented in ",a)') trim(bnd_yl), trim(dim)
            call warn(msg)
         end select  ! (bnd_yl)

         select case (bnd_yr)
         case ('cor', 'inf', 'mpi')
!         Do nothing
         case ('per')
            u(:,:,nyb+nb+1:nyb+2*nb,:)            = u(:,:,nb+1:2*nb,:)
         case ('user')
            call user_bnd_yr
         case ('ref')
            do ib=1,nb

               u((/iarr_all_dn,iarr_all_mx,iarr_all_mz/),:,nb+nyb+ib,:) = u((/iarr_all_dn,iarr_all_mx,iarr_all_mz/),:,nb+nyb+1-ib,:)
               u(iarr_all_my,:,nb+nyb+ib,:)               =-u(iarr_all_my,:,nb+nyb+1-ib,:)
#ifndef ISO
               u(iarr_all_en,:,nb+nyb+ib,:)               = u(iarr_all_en,:,nb+nyb+1-ib,:)
#endif /* !ISO */
#ifdef COSM_RAYS
               u(iarr_all_crs,:,nb+nyb+ib,:)               = u(iarr_all_crs,:,nb+nyb+1-ib,:)
#endif /* COSM_RAYS */
            enddo
         case ('out')
            do ib=1,nb
               u(:,:,nb+nyb+ib,:)                  = u(:,:,nb+nyb,:)
#ifdef COSM_RAYS
               u(iarr_all_crs,:,nb+nyb+ib,:)               = smallecr
#endif /* COSM_RAYS */
            enddo
         case ('outd')
            do ib=1,nb

               u((/iarr_all_dn,iarr_all_mx,iarr_all_mz/),:,nb+nyb+ib,:) = u((/iarr_all_dn,iarr_all_mx,iarr_all_mz/),:,nb+nyb,:)
               u(iarr_all_my,:,nb+nyb+ib,:)               = max(u(iarr_all_my,:,nb+nyb,:),0.0)
#ifndef ISO
               u(iarr_all_en,:,nb+nyb+ib,:)               = u(iarr_all_en,:,nb+nyb,:)
#endif /* !ISO */
#ifdef COSM_RAYS
               u(iarr_all_crs,:,nb+nyb+ib,:)               = smallecr
#endif /* COSM_RAYS */
            enddo
         case default
            write(msg,'("[fluid_boundaries:bnd_u]: Boundary condition ",a," not implemented in ",a)') trim(bnd_yr), trim(dim)
            call warn(msg)
         end select  ! (bnd_yr)

      case ('zdim')

         select case (bnd_zl)
         case ('mpi')
!         Do nothing if mpi
         case ('user')
            call user_bnd_zl
         case ('per')
            u(:,:,:,1:nb)                         = u(:,:,:,nzb+1:nzb+nb)
         case ('ref')
            do ib=1,nb

               u((/iarr_all_dn,iarr_all_mx,iarr_all_my/),:,:,nb+1-ib)   = u((/iarr_all_dn,iarr_all_mx,iarr_all_my/),:,:,nb+ib)
               u(iarr_all_mz,:,:,nb+1-ib)                 =-u(iarr_all_mz,:,:,nb+ib)
#ifndef ISO
               u(iarr_all_en,:,:,nb+1-ib)                 = u(iarr_all_en,:,:,nb+ib)
#endif /* !ISO */
#ifdef COSM_RAYS
               u(iarr_all_crs,:,:,nb+1-ib)                 = u(iarr_all_crs,:,:,nb+ib)
#endif /* COSM_RAYS */
            enddo
         case ('out')
            do ib=1,nb
               u(:,:,:,ib)                         = u(:,:,:,nb+1)
#ifdef COSM_RAYS
               u(iarr_all_crs,:,:,ib)                      = smallecr
#endif /* COSM_RAYS */
            enddo
         case ('outd')
            do ib=1,nb

               u((/iarr_all_dn,iarr_all_mx,iarr_all_my/),:,:,ib)        = u((/iarr_all_dn,iarr_all_mx,iarr_all_my/),:,:,nb+1)
 ! BEWARE: use of uninitialized value on first call (a side effect of r1726)
               u(iarr_all_mz,:,:,ib)                      = min(u(iarr_all_mz,:,:,nb+1),0.0)
#ifndef ISO
               u(iarr_all_en,:,:,ib)                      = u(iarr_all_en,:,:,nb+1)
#endif /* !ISO */
#ifdef COSM_RAYS
               u(iarr_all_crs,:,:,ib)                      = smallecr
#endif /* COSM_RAYS */
            enddo
#ifdef GRAV
         case ('outh')
            do ib=1,nb
               kb = nb+2-ib
               call outh_bnd(kb, kb-1, "min")
            enddo ! ib
#endif /* GRAV */
         case default
            write(msg,'("[fluid_boundaries:bnd_u]: Boundary condition ",a," not implemented in ",a)') trim(bnd_zl), trim(dim)
            call warn(msg)
         end select  ! (bnd_zl)

         select case (bnd_zr)
         case ('mpi')
!         Do nothing if mpi
         case ('user')
            call user_bnd_zr
         case ('per')
            u(:,:,:,nzb+nb+1:nzb+2*nb)            = u(:,:,:,nb+1:2*nb)
         case ('ref')
            do ib=1,nb

               u((/iarr_all_dn,iarr_all_mx,iarr_all_my/),:,:,nb+nzb+ib) = u((/iarr_all_dn,iarr_all_mx,iarr_all_my/),:,:,nb+nzb+1-ib)
               u(iarr_all_mz,:,:,nb+nzb+ib)               =-u(iarr_all_mz,:,:,nb+nzb+1-ib)
#ifndef ISO
               u(iarr_all_en,:,:,nb+nzb+ib)               = u(iarr_all_en,:,:,nb+nzb+1-ib)
#endif /* !ISO */
#ifdef COSM_RAYS
               u(iarr_all_crs,:,:,nb+nzb+ib)               = u(iarr_all_crs,:,:,nb+nzb+1-ib)
#endif /* COSM_RAYS */
            enddo
         case ('out')
            do ib=1,nb
               u(:,:,:,nb+nzb+ib)                  = u(:,:,:,nb+nzb)
#ifdef COSM_RAYS
               u(iarr_all_crs,:,:,nb+nzb+ib)               = smallecr
#endif /* COSM_RAYS */
            enddo
         case ('outd')
            do ib=1,nb

               u((/iarr_all_dn,iarr_all_mx,iarr_all_my/),:,:,nb+nzb+ib) = u((/iarr_all_dn,iarr_all_mx,iarr_all_my/),:,:,nb+nzb)
 ! BEWARE: use of uninitialized value on first call (a side effect of r1726)
               u(iarr_all_mz,:,:,nb+nzb+ib)               = max(u(iarr_all_mz,:,:,nb+nzb),0.0)
#ifndef ISO
               u(iarr_all_en,:,:,nb+nzb+ib)               = u(iarr_all_en,:,:,nb+nzb)
#endif /* !ISO */
#ifdef COSM_RAYS
               u(iarr_all_crs,:,:,nb+nzb+ib)               = smallecr
#endif /* COSM_RAYS */
            enddo
#ifdef GRAV
         case ('outh')
            do ib=1,nb
               kb = nb+nzb-1+ib
               call outh_bnd(kb, kb+1, "max")
            enddo ! ib
#endif /* GRAV */
         case default
            write(msg,'("[fluid_boundaries:bnd_u]: Boundary condition ",a," not implemented in ",a)') trim(bnd_zr), trim(dim)
            call warn(msg)
         end select  ! (bnd_zr)

      end select  ! (dim)

   end subroutine bnd_u

   subroutine all_fluid_boundaries
      use grid,  only: has_dir, xdim, ydim, zdim
      implicit none
      if (has_dir(xdim)) call bnd_u('xdim')
      if (has_dir(ydim)) call bnd_u('ydim')
      if (has_dir(zdim)) call bnd_u('zdim')

   end subroutine all_fluid_boundaries

end module fluidboundaries
