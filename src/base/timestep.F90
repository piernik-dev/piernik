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

!>
!! \brief This module gathers all applicable timestep limits and computes next timestep.
!<

module timestep

   implicit none

   private
   public :: time_step, cfl_manager

   real :: c_all_old
   procedure(), pointer :: cfl_manager => init_time_step

contains

!>
!! \brief Initialization routine
!!
!! \details This routine sets cfl_manager according to global::cflcontrol parameter.
!<

   subroutine init_time_step

      use constants,  only: PIERNIK_INIT_GLOBAL
      use dataio_pub, only: msg, die, warn, code_progress
      use global,     only: cflcontrol

      implicit none

      if (code_progress < PIERNIK_INIT_GLOBAL) call die("[grid:init_grid] globals not initialized.")

      select case (cflcontrol)
         case ('warn')
            cfl_manager => cfl_warn
         case ('auto', 'adaptive')
            cfl_manager => cfl_auto
         case ('none', '')
            if (associated(cfl_manager)) nullify(cfl_manager)
         case default
            write(msg, '(3a)')"[timestep:init_time_step] Unknown cfl_manager '",trim(cflcontrol),"'. Assuming 'none'."
            call warn(msg)
      end select
      if (.not.associated(cfl_manager)) call die("[timestep:init_time_step] cfl_manager was not associated.")
      call cfl_manager

   end subroutine init_time_step

!>
!! \brief Timestep calculation
!!
!! \details This routine calls various routines associated with different modules.
!! These routines return limit for timestep due to various physical and numerical conditions.
!! At the end the timestep is checked against remaining simulation time, minimum, and maximum allowed values etc.
!<
   subroutine time_step(dt, flind)

      use constants,            only: one, two, zero, half, I_ONE
      use dataio,               only: write_crashed
      use dataio_pub,           only: tend, msg, warn
      use fluidtypes,           only: var_numbers
      use global,               only: t, dt_old, dt_max_grow, dt_initial, dt_min, nstep
      use grid,                 only: leaves
      use cg_list,              only: cg_list_element
      use grid_cont,            only: grid_container
      use mpi,                  only: MPI_DOUBLE_PRECISION, MPI_MIN, MPI_MAX, MPI_IN_PLACE
      use mpisetup,             only: comm, mpi_err, master
      use timestep_pub,         only: c_all
      use timestepinteractions, only: timestep_interactions
#ifdef COSM_RAYS
      use timestepcosmicrays,   only: timestep_crs, dt_crs
#endif /* COSM_RAYS */
#ifdef RESISTIVE
      use resistivity,          only: dt_resist, timestep_resist
#endif /* RESISTIVE */
#ifdef DEBUG
      use piernikdebug,         only: has_const_dt, constant_dt
      use dataio_pub,           only: printinfo
#endif /* DEBUG */

      implicit none

      real, intent(inout) :: dt !< the timestep
      type(var_numbers), intent(in) :: flind

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      real :: c_, dt_
      integer :: ifl

! Timestep computation

      dt_old = dt

      c_all = zero
      dt = huge(one)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         !> \todo make the timestep_* routines members of fluidtypes::component_fluid
         do ifl = lbound(flind%all_fluids, dim=1), ubound(flind%all_fluids, dim=1)
            call timestep_fluid(cg, flind%all_fluids(ifl)%fl, dt_, c_)
            dt    = min(dt, dt_)
            c_all = max(c_all, c_)
         enddo

#ifdef COSM_RAYS
         call timestep_crs(cg)
         dt=min(dt,dt_crs)
#endif /* COSM_RAYS */

#ifdef RESISTIVE
         call timestep_resist(cg)
         dt = min(dt,dt_resist)
#endif /* RESISTIVE */

#ifndef BALSARA
         dt = min(dt,timestep_interactions(cg))
#endif /* BALSARA */

         cgl => cgl%nxt
      enddo

      call MPI_Allreduce(MPI_IN_PLACE, dt,    I_ONE, MPI_DOUBLE_PRECISION, MPI_MIN, comm, mpi_err)
      call MPI_Allreduce(MPI_IN_PLACE, c_all, I_ONE, MPI_DOUBLE_PRECISION, MPI_MAX, comm, mpi_err)

      ! finally apply some sanity factors
      if (nstep <=1) then
         if (dt_initial > zero) dt = min(dt, dt_initial)
      else
         if (dt_old > zero) dt = min(dt, dt_old*dt_max_grow)
      endif

      if (associated(cfl_manager)) call cfl_manager
      c_all_old = c_all

      if (dt < dt_min) then ! something nasty had happened
         if (master) then
            write(msg,'(2(a,es12.4))')"[timestep:time_step] dt = ",dt,", less than allowed minimum = ",dt_min
            call warn(msg)
         endif
         call write_crashed("[timestep:time_step] dt < dt_min")
      endif

      dt  = min(dt, (half*(tend-t)) + (two*epsilon(one)*((tend-t))))
#ifdef DEBUG
      ! We still need all above for c_all
      if (has_const_dt) then
         dt    = constant_dt
         write(msg,*) "[timestep:time_step]: (constant_dt) c_all = ", c_all
         call printinfo(msg)
      endif
#endif /* DEBUG */

   end subroutine time_step

!------------------------------------------------------------------------------------------
!>
!! \brief This routine detects sudden timestep changes due to strong velocity changes and interprets them as possible CFL criterion violations
!<

   subroutine cfl_warn

      use constants,  only: I_ONE
      use dataio_pub, only: msg, warn
      use global,     only: cfl, cfl_max, cfl_violated
      use mpi,        only: MPI_LOGICAL
      use mpisetup,   only: comm, mpi_err, master, FIRST
      use timestep_pub, only: c_all

      implicit none

      real :: stepcfl

      stepcfl = cfl
      if (c_all_old > 0.) stepcfl = c_all/c_all_old*cfl

      if (master) then
         msg = ''
         cfl_violated = .false.
         if (stepcfl > cfl_max) then
            write(msg,'(a,g10.3)') "[timestep:cfl_warn] Possible violation of CFL: ",stepcfl
            cfl_violated = .true.
         else if (stepcfl < 2*cfl - cfl_max) then
            write(msg,'(2(a,g10.3))') "[timestep:cfl_warn] Low CFL: ", stepcfl, " << ", cfl
         endif
         if (len_trim(msg) > 0) call warn(msg)
      endif

      call MPI_Bcast(cfl_violated, I_ONE, MPI_LOGICAL, FIRST, comm, mpi_err)

   end subroutine cfl_warn

!------------------------------------------------------------------------------------------
!>
!! \brief This routine detects sudden timestep changes due to strong velocity changes and interprets them as possible CFL criterion violations
!! Timestep changes are used to estimate safe value of CFL for the next timestep (EXPERIMENTAL)
!<

   subroutine cfl_auto

      use constants,  only: one, half, zero
      use dataio_pub, only: msg, warn
      use global,     only: cfl, cfl_max, dt, dt_old
      use mpisetup,   only: master
      use timestep_pub, only: c_all

      implicit none

      real, save :: stepcfl=zero, cfl_c=one
      real       :: stepcfl_old

      stepcfl_old = stepcfl
      stepcfl = cfl
      if (c_all_old > zero) then
         stepcfl = c_all/c_all_old*cfl
      else
         stepcfl_old = cfl
      endif

      if (stepcfl > zero .and. dt_old > zero) then
         cfl_c = min(one, half * (cfl_c + min(one, stepcfl_old/stepcfl *dt/dt_old)))
         dt = dt * cfl_c
      else
         cfl_c = one
      endif

      if (master) then
         msg = ''
         if (stepcfl > cfl_max) then
            write(msg,'(a,g10.3)') "[timestep:cfl_auto] Possible violation of CFL: ",stepcfl
         else if (stepcfl < 2*cfl - cfl_max) then
            write(msg,'(2(a,g10.3))') "[timestep:cfl_auto] Low CFL: ", stepcfl, " << ", cfl
         endif
         if (len_trim(msg) > 0) call warn(msg)
      endif

   end subroutine cfl_auto

!------------------------------------------------------------------------------------------
!>
!! \brief %Timestep computation for the fluid (FIXME)
!!
!! %Timestep for the fluid is set as the minimum %timestep for all of the MPI blocks times the Courant number.
!! To compute the %timestep in each MPI block, the fastest speed at which information travels in each direction is computed as
!! \f{equation}
!! c_x=\max\limits_{i,j,k}{\left(v_x^{i,j,k}+c_f^{i,j,k}\right)},
!! \f}
!! where \f$v_x^{i,j,k}\f$ is the maximum speed in \f$x\f$ direction for the cell \f$(i,j,k)\f$ and \f$c_f^{i,j,k}\f$ is the speed of sound for
!! ionized fluid computed as \f$c_f^{i,j,k}=\sqrt{\left|\frac{2p_{mag}+\gamma p}{\rho^{i,j,k}}\right|}\f$, where \f$p\f$ stands for pressure,
!! \f$p_{mag}\f$ is pressure of magnetic field, \f$\gamma\f$ is the adiabatic index of the ionized fluid and \f$\rho^{i,j,k}\f$ is fluid density in the cell
!! \f$(i,j,k)\f$. For directions \f$y, z\f$ the computations are made in similar way.
!!
!! %Timestep for each MPI block is then computed as
!! \f{equation}
!! dt=\min{\left(\left|\frac{dx}{c_x}\right|,\left|\frac{dy}{c_y}\right|,\left|\frac{dz}{c_z}\right|\right)},
!! \f}
!! where \f$dx\f$, \f$dy\f$ and \f$dz\f$ are the cell lengths in each direction.
!!
!! Information about the computed %timesteps is exchanged between MPI blocks in order to choose the minimum %timestep for the fluid.
!! The final %timestep is multiplied by the Courant number specified in parameters of each task.
!<

   subroutine timestep_fluid(cg, fl, dt, c_fl)

      use constants,  only: big, xdim, ydim, zdim, ndims, GEO_RPZ, LO, ndims
      use domain,     only: dom
      use fluidtypes, only: component_fluid
      use global,     only: cfl
      use grid_cont,  only: grid_container

      implicit none

      type(grid_container), pointer, intent(in) :: cg !< current grid container
      real, intent(out)                         :: dt !< resulting timestep
      real, intent(out)                         :: c_fl !< maximum speed at which information travels in the fluid

      real, dimension(ndims) :: c !< maximum velocity in all directions
      real, dimension(ndims) :: v !< maximum velocity of fluid in all directions

      ! locals
      class(component_fluid), pointer, intent(in) :: fl
      integer :: i, j, k
      real, dimension(ndims) :: dt_proc !< timestep for the current cg
      integer :: d

      c(:) = 0.0
      c_fl = 0.0

      do k = cg%ks, cg%ke
         do j = cg%js, cg%je
            do i = cg%is, cg%ie
               if (cg%u(fl%idn,i,j,k) > 0.0) then
                  v(:) = abs(cg%u(fl%imx:fl%imz, i, j, k) / cg%u(fl%idn, i, j, k))
               else
                  v(:) = 0.0
               endif

               c(:) = max( c(:), v(:) + fl%get_cs(cg, i, j, k) )
               c_fl = max(c_fl, maxval(c(:)))
            enddo
         enddo
      enddo

      do d = xdim, zdim
         if (dom%has_dir(d) .and. c(d) /= 0.) then
            dt_proc(d) = cg%dl(d)/c(d)
            if (dom%geometry_type == GEO_RPZ .and. d == ydim) dt_proc(d) = dt_proc(d) * dom%edge(xdim, LO)
         else
            dt_proc(d) = big
         endif
      enddo

      dt = cfl * minval(dt_proc)
      call fl%set_c(c_fl)

   end subroutine timestep_fluid

end module timestep
