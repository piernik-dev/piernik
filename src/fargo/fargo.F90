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
!! \brief Implementation of a fast eulerian transport algorithm for differentially rotating disks (Masset 2000)
!!
!! \details See also:
!!   1. Masset, F. "FARGO: A fast eulerian transport algorithm for differentially rotating disks" (2000) A&A, 141:165-173, arXiv:astro-ph/9910390
!!   2. Kley, W., Bitsch, B., Klahr, H. "Planet migration in three-dimensional radiative discs" (2009) A&A, 506:971-987, arXiv:0908.1863
!!
!<
module fargo
! pulled by ANY
   implicit none

   private
   public :: init_fargo, cleanup_fargo, make_fargosweep, timestep_fargo, fargo_mean_omega

contains

!>
!! \brief FARGO initialization
!!
!! \details Basic sanity checks are performed here
!<
   subroutine init_fargo

      use constants,    only: GEO_RPZ, ydim
      use dataio_pub,   only: die, warn
      use domain,       only: dom, AMR_bsize
      use global,       only: use_fargo
      use mpisetup,     only: master

      implicit none

      if (.not. use_fargo) return
      if (dom%geometry_type /= GEO_RPZ) call die("[fargo:init_fargo] FARGO works only for cylindrical geometry")
      if (master) call warn("[fargo:init_fargo] BEWARE: Fast eulerian transport is an experimental feature")

      if (AMR_bsize(ydim) /= 0 .and. AMR_bsize(ydim) /= dom%n_d(ydim)) &
           call die("[fargo:init_fargo] AMR is not allowed in the azimuthal direction. You must set AMR_bsize(ydim) == dom%n_d(ydim) AND make sure that there will be no azimuthal refinement steps on higher levels")
      !> \todo Implement proper checks (or even refinement fixups) and do not impose restrictions on AMR_bsize(ydim)

   end subroutine init_fargo

!>
!! \brief FARGO finalization
!<
   subroutine cleanup_fargo
      implicit none
   end subroutine cleanup_fargo

!>
!! \brief Perform azimuthal integration steep using FARGO
!!
!! \details As per Masset (2000) we need to split azimuthal integration step
!!    into several phases:
!!    1. source terms and transport with the residual velocity
!!    2. transport with constant residual velocity
!!    3. transport with mean velocity
!<
   subroutine make_fargosweep
      use constants,   only: VEL_RES, VEL_CR, ydim
      use global,      only: dt, skip_sweep
      use sweeps,      only: sweep

      implicit none

      ! TODO we are omitting B and cr update, but FARGO does not work with them yet...
      if (.not.skip_sweep(ydim)) then
         call get_fargo_vels(dt)
         call sweep(ydim, VEL_RES)  ! 1.
         call sweep(ydim, VEL_CR)   ! 2.
         call int_shift             ! 3.
      endif
   end subroutine make_fargosweep

!>
!! \brief Compute mean omega for each radius globally
!<
   subroutine fargo_mean_omega

      use constants,          only: xdim, ydim, zdim, I_ONE, pSUM
      use cg_level_base,      only: base
      use cg_level_connected, only: cg_level_connected_T
      use cg_list,            only: cg_list_element
      use global,             only: use_fargo
      use grid_cont,          only: grid_container
      use fluidindex,         only: flind
      use fluidtypes,         only: component_fluid
      use mpisetup,           only: piernik_MPI_Allreduce

      implicit none

      type(cg_list_element), pointer    :: cgl
      type(grid_container),  pointer    :: cg
      integer :: ifl, i
      class(component_fluid), pointer :: pfl
      type(cg_level_connected_T), pointer :: curl

      curl => base%level
      do while (associated(curl))

         if (.not.allocated(curl%local_omega)) allocate(curl%local_omega(curl%off(xdim):curl%off(xdim)+curl%n_d(xdim)-I_ONE, flind%fluids))
         if (.not.allocated(curl%cell_count)) allocate(curl%cell_count(curl%off(xdim):curl%off(xdim)+curl%n_d(xdim)-I_ONE))
         curl%local_omega(:,:) = 0.0
         curl%cell_count(:) = 0

         if (use_fargo) then
            cgl => curl%first
            do while (associated(cgl))
               cg => cgl%cg
               do i = cg%is, cg%ie
                  do ifl = 1, flind%fluids
                     pfl => flind%all_fluids(ifl)%fl
                     curl%local_omega(i, ifl) = curl%local_omega(i, ifl) + sum(cg%u(pfl%imy, i, cg%js:cg%je, cg%ks:cg%ke) / cg%u(pfl%idn, i, cg%js:cg%je, cg%ks:cg%ke) / cg%x(i))
                  enddo
                  curl%cell_count(i) = curl%cell_count(i) + product(cg%n_b(ydim:zdim))
               enddo
               cgl => cgl%nxt
            enddo
         endif

         if (use_fargo) then
            call piernik_MPI_Allreduce(curl%local_omega, pSUM)
            call piernik_MPI_Allreduce(curl%cell_count, pSUM)
            do i = lbound(curl%cell_count, dim=1), ubound(curl%cell_count, dim=1)
               curl%local_omega(i, :) = curl%local_omega(i, :) / curl%cell_count(i)
            enddo
         endif

         curl => curl%finer

      enddo

   end subroutine fargo_mean_omega

!>
!! \brief Compute shift, constant residual and residual azimuthal velocity
!!
!! \details This routines fills FARGO auxiliary arrays with valid data
!<
   subroutine get_fargo_vels(dt)

      use constants,          only: xdim, ydim, I_ONE
      use cg_level_base,      only: base
      use cg_level_connected, only: cg_level_connected_T
      use domain,             only: dom
      use fluidindex,         only: flind

      implicit none

      real, intent(in) :: dt
      integer :: ifl
      type(cg_level_connected_T), pointer :: curl
      real :: dl_ydim

      call fargo_mean_omega

      curl => base%level
      do while (associated(curl))

         if (.not. allocated(curl%omega_mean)) allocate(curl%omega_mean(curl%off(xdim):curl%off(xdim)+curl%n_d(xdim)-I_ONE, flind%fluids))
         if (.not. allocated(curl%omega_cr))   allocate(curl%omega_cr  (curl%off(xdim):curl%off(xdim)+curl%n_d(xdim)-I_ONE, flind%fluids))
         if (.not. allocated(curl%nshift))     allocate(curl%nshift    (curl%off(xdim):curl%off(xdim)+curl%n_d(xdim)-I_ONE, flind%fluids))

         dl_ydim = dom%L_(ydim)/curl%n_d(ydim)

         do ifl = 1, flind%fluids
            curl%omega_mean(:, ifl) = curl%local_omega(:, ifl)
            curl%nshift(:, ifl) = nint(curl%omega_mean(:, ifl) * dt / dl_ydim)
            curl%omega_cr(:, ifl) = curl%omega_mean(:, ifl) - curl%nshift(:, ifl) * dl_ydim / dt
         enddo

         if (allocated(curl%local_omega)) deallocate(curl%local_omega)

         curl => curl%finer
      enddo

   end subroutine get_fargo_vels

!>
!! \brief Compute FARGO time constraint
!!
!! \details "An additional time step limitation is given by the requirement that
!!    the shift should not disconnect two neighbouring grid cells in radial and
!!    in the vertical direction" -- Kley et al. 2009
!<
   real function timestep_fargo(cg, dt_min) result(dt)

      use fluidindex,   only: flind
      use grid_cont,    only: grid_container
      use global,       only: cfl
      use constants,    only: ydim, big, small, half
      use fluidtypes,   only: component_fluid

      implicit none
      type(grid_container), pointer, intent(in) :: cg
      real, intent(in) :: dt_min

      real :: dt_shear!, dt_res
      !real :: v_mean, v_cr, nshft, c_fl
      real :: omega!, omegap, dphi
      integer :: i, j, k, ifl
      class(component_fluid), pointer :: pfl


      dt_shear = big
      do ifl = lbound(flind%all_fluids, 1), ubound(flind%all_fluids, 1)
         pfl   => flind%all_fluids(ifl)%fl
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  if (cg%leafmap(i, j, k)) then
                     omega  = cg%u(pfl%imy, i, j, k) / cg%u(pfl%idn, i, j, k) / cg%x(i) - cg%u(pfl%imy, i-1, j, k) / cg%u(pfl%idn, i-1, j, k) / cg%x(i-1)
                     omega  = max(abs(omega), small)
                     !omegap = (cg%u(pfl%imy, i, j, k) / cg%u(pfl%idn, i, j, k) - cg%u(pfl%imy, i, j-1, K) / cg%u(pfl%idn, i, j-1, k)) / cg%x(i)
                     !omegap = max(abs(omegap), small)
                     !dphi = cg%dl(ydim)
                     !dt_shear = min(dt_shear, half*min(dphi / abs(omega), dphi / abs(omegap)))
                     dt_shear = min(dt_shear, half * cg%dl(ydim) / omega)
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

!>
!! \brief Perform integer part of azimuthal transport step
!!
!! \details Routine that shifts data for an integer number of cells.
!! \todo reduce MPI communication
!<
   subroutine int_shift

      use cg_level_base,      only: base
      use cg_level_connected, only: cg_level_connected_T
      use cg_list,            only: cg_list_element
      use constants,          only: ydim
      use domain,             only: dom
      use fluidindex,         only: flind
      use fluidtypes,         only: component_fluid
      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
#ifdef SELF_GRAV
      use constants,          only: sgp_n
      use named_array_list,   only: qna
#endif /* SELF_GRAV */

      implicit none

      type(cg_list_element), pointer    :: cgl
      type(grid_container),  pointer    :: cg
      class(component_fluid), pointer   :: pfl
      integer :: ifl, i, max_nshift, iter
      type(cg_level_connected_T), pointer :: curl

      curl => base%level
      do while (associated(curl))
         max_nshift = maxval(curl%nshift)

         do iter = 1, max(ceiling(float(max_nshift) / float(int(dom%nb))), 1)
            cgl => curl%first
            do while (associated(cgl))
               cg => cgl%cg
               do i = cg%is, cg%ie
                  if (all(curl%nshift(i, :) == 0)) cycle
                  do ifl = 1, flind%fluids
                     pfl   => flind%all_fluids(ifl)%fl
                     cg%u(pfl%beg:pfl%end, i, :, :) = cshift(cg%u(pfl%beg:pfl%end, i, :, :), -min(curl%nshift(i, ifl), int(dom%nb)), dim=2)
                  enddo
#ifdef SELF_GRAV
                  cg%sgp(i, :, :) = cshift(cg%sgp(i, :, :), -min(curl%nshift(i, 1), int(dom%nb)), dim=1) ! TODO: what about ifl?
#endif /* SELF_GRAV */
               enddo
               curl%nshift(:, :) = max(curl%nshift(:, :) - dom%nb, 0)
               cgl => cgl%nxt
            enddo
            call curl%arr4d_boundaries(wna%fi, dir=ydim, nocorners=.true.)
            call curl%bnd_u(ydim)
#ifdef SELF_GRAV
            call curl%arr3d_boundaries(qna%ind(sgp_n))!, dir=ydim) !, nocorners=.true.)
#endif /* SELF_GRAV */
         enddo

         curl => curl%finer
      enddo

   end subroutine int_shift

end module fargo
