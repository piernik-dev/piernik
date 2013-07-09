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
#include "macros.h"
!>
!! Implementation of a fast eulerian transport algorithm for differentially rotating disks (Masset 2000)
!!
!! See also:
!!   1. Masset, F. "FARGO: A fast eulerian transport algorithm for differentially rotating disks" (2000) A&A, 141:165-173, arXiv:astro-ph/9910390
!!   2. Kley, W., Bitsch, B., Klahr, H. "Planet migration in three-dimensional radiative discs" (2009) A&A, 506:971-987, arXiv:0908.1863
!!
!<
module fargo
! pulled by ANY
   implicit none
   real,    dimension(:, :, :),  allocatable :: vphi_mean
   real,    dimension(:, :, :),  allocatable :: vphi_cr
   integer, dimension(:, :, :),  allocatable :: nshift

   private
   public :: init_fargo, cleanup_fargo, vphi_mean, vphi_cr, nshift, get_fargo_vels, timestep_fargo, int_shift

contains

   subroutine init_fargo

      use constants,    only: GEO_RPZ
      use dataio_pub,   only: die, warn
      use domain,       only: dom
      use global,       only: use_fargo

      implicit none

      if (.not. use_fargo) return

      if (dom%geometry_type /= GEO_RPZ) call die("[fargo:init_fargo] FARGO works only for cylindrical geometry")

      call warn("[fargo:init_fargo] BEWARE: Fast eulerian transport is an experimental feature")

      ! TODO: Check if we have domain division in ydim

   end subroutine init_fargo

   subroutine cleanup_fargo
      implicit none

      if (allocated(vphi_mean)) deallocate(vphi_mean)
      if (allocated(vphi_cr)) deallocate(vphi_cr)
      if (allocated(nshift)) deallocate(nshift)

   end subroutine cleanup_fargo

   subroutine get_fargo_vels(dt)

      use constants,        only: xdim, ydim, LO, HI, pMAX
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use dataio_pub,       only: die, warn, msg
      use domain,           only: dom
      use grid_cont,        only: grid_container
      use fluidindex,       only: flind
      use fluidtypes,       only: component_fluid
      use mpisetup,         only: master, piernik_MPI_Allreduce

      implicit none

      real, intent(in) :: dt

      type(cg_list_element), pointer    :: cgl
      type(grid_container),  pointer    :: cg
      integer :: ifl, icg, i, max_nshift
      class(component_fluid), pointer :: pfl

      cgl => leaves%first
      cg => cgl%cg
      if (.not. allocated(vphi_mean)) allocate(vphi_mean(cg%lhn(xdim, LO):cg%lhn(xdim, HI), flind%fluids, leaves%cnt))
      if (.not. allocated(vphi_cr)) allocate(vphi_cr(cg%lhn(xdim, LO):cg%lhn(xdim, HI), flind%fluids, leaves%cnt))
      if (.not. allocated(nshift)) allocate(nshift(cg%lhn(xdim, LO):cg%lhn(xdim, HI), flind%fluids, leaves%cnt))

      icg = 1
      do while (associated(cgl))
         cg => cgl%cg
         do i = cg%lhn(xdim, LO), cg%lhn(xdim, HI)
            do ifl = 1, flind%fluids
               pfl => flind%all_fluids(ifl)%fl
               vphi_mean(i, ifl, icg) = sum(cg%u(pfl%imy, i, :, :) / cg%u(pfl%idn, i, :, :))  / size(cg%u(pfl%idn, i, :, :))
            enddo
         enddo

         do ifl = 1, flind%fluids
            pfl => flind%all_fluids(ifl)%fl
            nshift(:, ifl, icg) = nint(vphi_mean(:, ifl, icg) * dt / (cg%x(:) * cg%dl(ydim)))
            vphi_cr(:, ifl, icg) = vphi_mean(:, ifl, icg) - nshift(:, ifl, icg) * (cg%x(:) * cg%dl(ydim)) / dt
         enddo

         cgl => cgl%nxt
         icg = icg + 1
      enddo
      max_nshift = maxval(nshift)
      call piernik_MPI_Allreduce(max_nshift, pMAX)
      if (master .and. max_nshift > dom%nb) then
         write(msg, '(a,I2,a)') "[fargo:get_fargo_vels] FARGO would require ", max_nshift, " ghostcells to work."
         call warn(msg)
         call die("[fargo:get_fargo_vels] FARGO does not support domain division in ydim yet")
      endif

   end subroutine get_fargo_vels

   real function timestep_fargo(cg, dt_min) result(dt)

      use fluidindex,   only: flind
      use grid_cont,    only: grid_container
      use global,       only: cfl
      use constants,    only: ydim
      use fluidtypes,       only: component_fluid

      implicit none
      type(grid_container), pointer, intent(in) :: cg
      real, intent(in) :: dt_min

      real :: dt_shear!, dt_res
      !real :: v_mean, v_cr, nshft, c_fl
      real :: vphi, vphip, dphi
      integer :: i, j, k, ifl
      class(component_fluid), pointer :: pfl


      dt_shear = huge(real(1.0,4))
      do ifl = 1, flind%fluids
         pfl   => flind%all_fluids(ifl)%fl
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  if (cg%leafmap(i, j, k)) then
                     vphi  = cg%u(pfl%imy, i, j, k) / cg%u(pfl%idn, i, j, k) - cg%u(pfl%imy, i-1, j, k) / cg%u(pfl%idn, i-1, j, k)
                     vphi  = max(abs(vphi), 1e-8)
                     vphip = cg%u(pfl%imy, i, j, k) / cg%u(pfl%idn, i, j, k) - cg%u(pfl%imy, i, j-1, K) / cg%u(pfl%idn, i, j-1, k)
                     vphip = max(abs(vphip), 1e-8)
                     dphi = cg%x(i) * cg%dl(ydim)
                     dt_shear = min(dt_shear, 0.5*min(dphi / abs(vphi), dphi / abs(vphip)))
                  endif
               enddo
            enddo
         enddo
      enddo
      dt = cfl * dt_shear

!     dt_res = huge(real(1.0,4))
!     c_fl = 0.0
!     do ifl = 1, flind%fluids
!        pfl   => flind%all_fluids(ifl)%fl
!        do i = cg%is, cg%ie
!           dphi = cg%x(i) * cg%dl(ydim)
!           v_mean = sum(cg%u(pfl%imy, i, :, :) / cg%u(pfl%idn, i, :, :)) / size(cg%u(pfl%idn, i, :, :))
!           nshft = nint(v_mean * dt / dphi)
!           v_cr = v_mean - nshft * dphi / dt
!           c_fl = max(c_fl, abs(v_cr) + pfl%get_cs(i, cg%js, cg%ks, cg%u, cg%b, cg%cs_iso2))
!           dt_res = min(dt_res, dphi / c_fl)
!        enddo
!     enddo
!
!     dt = min(dt, cfl * dt_res)
      if (.false.) print *, dt_min

   end function timestep_fargo

   subroutine int_shift
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, LO, HI
      use fluidindex,       only: flind
      use fluidtypes,       only: component_fluid
      use grid_cont,        only: grid_container

      implicit none

      type(cg_list_element), pointer    :: cgl
      type(grid_container),  pointer    :: cg
      class(component_fluid), pointer   :: pfl
      integer :: icg, ifl, i

      cgl => leaves%first
      icg = 1
      do while (associated(cgl))
         cg => cgl%cg
         do i = cg%lhn(xdim, LO), cg%lhn(xdim, HI)
            do ifl = 1, flind%fluids
               pfl   => flind%all_fluids(ifl)%fl
               cg%u(pfl%beg:pfl%end, i, :, :) = cshift(cg%u(pfl%beg:pfl%end, i, :, :), -nshift(i, ifl, icg), dim=2)
            enddo
         enddo
         cgl => cgl%nxt
         icg = icg + 1
      enddo
   end subroutine int_shift

end module fargo
