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

!>  \brief Integrators for particles

module particle_solvers
! pulled by NBODY

   implicit none

   private
   public :: init_psolver, psolver, hermit_4ord, leapfrog_2ord, update_particle_kinetic_energy

   procedure(particle_solver_P), pointer :: psolver => NULL()

   abstract interface
      subroutine particle_solver_P(forward)
         logical, optional, intent(in) :: forward
      end subroutine particle_solver_P
   end interface

contains

!> \brief Initialize psolver
   subroutine init_psolver

      use dataio_pub,   only: msg, die
      use particle_pub, only: default_ti, time_integrator

      implicit none

      psolver => null()
      select case (trim(time_integrator))
#ifndef _CRAYFTN
         case ('hermit4')
            psolver => hermit_4ord
         case ('leapfrog2')
            psolver => leapfrog_2ord
         case (default_ti) ! be quiet
#endif /* !_CRAYFTN */
         case default
            write(msg, '(3a)')"[particle_solvers:init_psolver] Unknown integrator '",trim(time_integrator),"'"
            call die(msg)
      end select

   end subroutine init_psolver

   !>
   !!  from:   Moving Stars Around (Piet Hut and Jun Makino, 2003).
   !!  http://www.ids.ias.edu/~piet/act/comp/algorithms/starter/index.html
   !<

   subroutine hermit_4ord(forward)

      use cg_leaves,  only: leaves
      use constants,  only: ndims, xdim, zdim
      use dataio_pub, only: die
      use domain,     only: is_refined, is_multicg
      use global,     only: t, dt
      use grid_cont,  only: grid_container
      use particle_utils, only: count_all_particles
      use particle_types, only: particle

      implicit none

      logical, optional, intent(in) :: forward

      real, parameter :: dt_param = 0.03        ! control parameter to determine time step size
      real, parameter :: dt_dia = 1             ! time interval between diagnostics output
      real, parameter :: dt_out = 0.01          ! time interval between output of snapshots

      real, dimension(:),    allocatable :: mass
      real, dimension(:, :), allocatable :: pos, vel, acc, jerk
      type(grid_container),  pointer     :: cg
      type(particle), pointer            :: pset

      real    :: epot, coll_time, t_dia, t_out, t_end, einit, dth, th
      integer :: nsteps, n, ndim, lun_out, lun_err, i, p

      if (is_refined) call die("[particle_solvers:hermit_4ord] AMR not implemented for hermit 4ord solver yet")
      if (is_multicg) call die("[particle_solvers:hermit_4ord] multi_cg not implemented for hermit 4ord solver yet")
      cg => leaves%first%cg

      open(newunit=lun_out, file='nbody_out.log', status='unknown',  position='append')
      open(newunit=lun_err, file='nbody_err.log', status='unknown',  position='append')

      n = count_all_particles()
      th = t - dt
      allocate(mass(n), pos(n, ndims), vel(n, ndims), acc(n, ndims), jerk(n, ndims))

      !Maybe not ideal with a list here?
      pset => cg%pset%first
      p=1
      do while (associated(pset))
         mass(p) = pset%pdata%mass
         do ndim = xdim, zdim
            pos(p, ndim) = pset%pdata%pos(ndim)
            vel(p, ndim) = pset%pdata%vel(ndim)
         enddo
         pset => pset%nxt
         p = p + 1
      enddo

      write(lun_err, *) "Starting a Hermite integration for a ", n, "-body system,"
      write(lun_err, *) "  from time t = ", th, " with time step control parameter dt_param = ", dt_param, "  until time ", th + dt, " ,"
      write(lun_err, *) "  with diagnostics output interval dt_dia = ", dt_dia, ","
      write(lun_err, *) "  and snapshot output interval dt_out = ", dt_out, "."

      call get_acc_jerk_pot_coll(mass, pos, vel, acc, jerk, n, epot, coll_time)

      nsteps = 0

      call write_diagnostics(mass, pos, vel, acc, jerk, n, th, epot, nsteps, einit, .True.)

      t_dia = th + dt_dia  ! next time for diagnostics output
      t_out = th + dt_out  ! next time for snapshot output
      t_end = th + dt      ! final time, to finish the integration

      do
         do while (th < t_dia .and. th < t_out .and. th < t_end)
            dth = dt_param * coll_time
            call evolve_step(mass, pos, vel, acc, jerk, n, th, dth, epot, coll_time)
            nsteps = nsteps + 1
         enddo
         if (th >= t_dia) then
            call write_diagnostics(mass, pos, vel, acc, jerk, n, th, epot, nsteps, einit, .False.)
            t_dia = t_dia + dt_dia
         endif
         if (th >= t_out) then
            write(lun_out, *) "#", th
            do i = 1, n
               write(lun_out, '(7(E13.6,1X))') mass(i), pos(i,:), vel(i,:)
            enddo
            t_out = t_out + dt_out
         endif
         if (th >= t_end) exit
      enddo

      pset => cg%pset%first
      p=1
      do while (associated(pset))
         do ndim = xdim, zdim
            pset%pdata%pos(ndim) = pos(p, ndim)
            pset%pdata%vel(ndim) = vel(p, ndim)
         enddo
         pset => pset%nxt
         p = p + 1
      enddo

      deallocate (mass, pos, vel, acc, jerk)
      close(lun_out)
      close(lun_err)
      return
      if (forward) return ! suppress compiler warnings

   contains

      subroutine write_diagnostics(mass, pos, vel, acc, jerk, n, t, epot, nsteps, einit, init_flag)

         implicit none

         integer,                   intent(in)    :: n, nsteps
         real, dimension(n),        intent(in)    :: mass
         real, dimension(n, ndims), intent(in)    :: pos, vel, acc, jerk
         real,                      intent(in)    :: t, epot
         real,                      intent(inout) :: einit
         logical,                   intent(in)    :: init_flag

         real :: ekin, etot

         ekin = sum(0.5 * mass(:) * sum(vel(:,:)**2, dim=2))
         etot = ekin + epot

         if (init_flag) einit = etot  ! at first pass, pass the initial energy back to the calling function

         write(lun_err, *) "at time t = ", t, " , after ", nsteps, " steps :"
         write(lun_err, *) " E_kin = ", ekin, " , E_pot = ", epot, " , E_tot = ", etot
         write(lun_err, *) "                absolute energy error: E_tot - E_init = ", etot - einit
         write(lun_err, *) "                relative energy error: (E_tot - E_init) / E_init = ", (etot - einit) / einit

         if (.false.) write(lun_err, *) pos, acc, jerk ! suppress compiler warnings

      end subroutine write_diagnostics

   end subroutine hermit_4ord

   subroutine evolve_step(mass, pos, vel, acc, jerk, n, t, dt, epot, coll_time)

      use constants, only: ndims

      implicit none

      integer,                  intent(in)    :: n
      real, dimension(n),       intent(in)    :: mass
      real, dimension(n,ndims), intent(out)   :: pos, vel
      real, dimension(n,ndims), intent(inout) :: acc, jerk
      real,                     intent(in)    :: dt
      real,                     intent(inout) :: t, epot, coll_time

      real, dimension(n,ndims)                :: old_pos, old_vel, old_acc, old_jerk

      old_pos = pos
      old_vel = vel
      old_acc = acc
      old_jerk = jerk

      call predict_step(               pos, vel, acc, jerk, n, dt)
      call get_acc_jerk_pot_coll(mass, pos, vel, acc, jerk, n, epot, coll_time)
      call correct_step(               pos, vel, acc, jerk, old_pos, old_vel, old_acc, old_jerk, n, dt);

      t = t + dt

   contains

      subroutine predict_step(pos, vel, acc, jerk, n, dt)

         use constants, only: half, onesth

         implicit none

         integer,                   intent(in)  :: n
         real, dimension(n, ndims), intent(out) :: pos, vel
         real, dimension(n, ndims), intent(in)  :: acc, jerk
         real,                      intent(in)  :: dt

         real                                   :: hdt, hdt2

         hdt  = half   * dt**2
         hdt2 = onesth * dt**3

         pos(:,:) = pos(:,:) + vel(:,:)*dt + acc(:,:)*hdt + jerk(:,:)*hdt2
         vel(:,:) = vel(:,:) + acc(:,:)*dt + jerk(:,:)*hdt

      end subroutine predict_step

      subroutine correct_step(pos, vel, acc, jerk, old_pos, old_vel, old_acc, old_jerk, n, dt)

         use constants, only: half, onet

         implicit none

         integer,                   intent(in)  :: n
         real, dimension(n, ndims), intent(out) :: pos, vel
         real, dimension(n, ndims), intent(in)  :: acc, jerk, old_pos, old_vel, old_acc, old_jerk
         real,                      intent(in)  :: dt

         real                                   :: hdt, hdt2

         hdt  = half * dt
         hdt2 = onet * hdt**2

         vel(:,:) = old_vel(:,:) + (old_acc(:,:) + acc(:,:))*hdt + (old_jerk(:,:) - jerk(:,:)) * hdt2
         pos(:,:) = old_pos(:,:) + (old_vel(:,:) + vel(:,:))*hdt + (old_acc(:,:)  - acc(:,:) ) * hdt2

      end subroutine correct_step

   end subroutine evolve_step

   subroutine get_acc_jerk_pot_coll(mass, pos, vel, acc, jerk, n, epot, coll_time)

      use constants, only: ndims

      implicit none

      integer,                  intent(in)  :: n
      real, dimension(n),       intent(in)  :: mass
      real, dimension(n,ndims), intent(in)  :: pos
      real, dimension(n,ndims), intent(in)  :: vel
      real, dimension(n,ndims), intent(out) :: acc
      real, dimension(n,ndims), intent(out) :: jerk
      real,                     intent(out) :: epot
      real,                     intent(out) :: coll_time

      real, dimension(ndims)                :: rji, vji, da, dj
      integer                               :: i, j
      real                                  :: coll_time_q  ! collision time to 4th power

      real :: r   ! | rji |
      real :: r2  ! | rji |^2
      real :: r3  ! | rji |^3
      real :: v2  ! | vji |^2
      real :: rv_r2 ! ( rij . vij ) / | rji |^2
      real :: da2

      acc(:,:)  = 0.0
      jerk(:,:) = 0.0
      epot      = 0.0

      coll_time_q = huge(1.0)  ! collision time to 4th power

      do i = 1, n
         do j = i+1, n
            rji(:) = pos(j, :) - pos(i, :)
            vji(:) = vel(j, :) - vel(i, :)

            r2 = sum(rji**2)
            v2 = sum(vji**2)
            rv_r2 = sum(rji*vji) / r2

            r = sqrt(r2)
            r3 = r * r2

            ! add the {i,j} contribution to the total potential energy for the
            ! system
            epot = epot - mass(i) * mass(j) / r

            da(:) = rji(:) / r3
            dj(:) = (vji(:) - 3.0 * rv_r2 * rji(:)) / r3

            acc(i,:) = acc(i,:) + mass(j) * da(:)
            acc(j,:) = acc(j,:) - mass(i) * da(:)
            jerk(i,:) = jerk(i,:) + mass(j) * dj(:)
            jerk(j,:) = jerk(j,:) - mass(i) * dj(:)

            if (v2 > 0.0) &
               coll_time_q = min(coll_time_q, r2**2 / v2**2) ! first collision time estimate, based on unaccelerated linear motion

            da2 = sum(da**2) * (mass(i) + mass(j))**2  ! square of the pair-wise acceleration between particles i and j
            coll_time_q = min(coll_time_q, r2 / da2)   ! second collision time estimate, based on free fall
         enddo
      enddo
      coll_time = sqrt(sqrt(coll_time_q))
      return

   end subroutine get_acc_jerk_pot_coll

   subroutine leapfrog_2ord(forward)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: half, two
      use global,           only: dt
      use particle_diag,    only: particle_diagnostics
      use particle_gravity, only: update_particle_gravpot_and_acc
      use particle_pub,     only: twodtscheme

      implicit none

      logical, optional, intent(in)  :: forward
      real                           :: dt_kick   !< timestep for kicks
      type(cg_list_element), pointer :: cgl

      if (twodtscheme) then
         if (forward) then
            call kick(dt)                        !1. kick
         else
            call drift(two*dt)                   !2. drift
            call update_particle_gravpot_and_acc
            call kick(dt)                        !3. kick
            call update_particle_kinetic_energy
            call particle_diagnostics(.true.)
         endif
      else
         dt_kick = dt * half
         call kick(dt_kick)                   !1. kick
         call drift(dt)                       !2. drift
         call update_particle_gravpot_and_acc
         call kick(dt_kick)                   !3. kick
         call update_particle_kinetic_energy
         call particle_diagnostics(.true.)
      endif

      ! Since the nbdn field is updated in update_particle_gravpot_and_acc (prior to call to source_terms_grav), it contains slightly outdated info.
      ! But call to map_particles is relatively expensive, so in case of an urgent need of up-to-date nbdn it is advised to hook an extra
      ! call to map_particles somewhere in the I/O routines to avoid the overhead in each step.

   contains

      subroutine kick(kdt)

         use cg_cost_data,   only: I_PARTICLE
         use constants,      only: PPP_PART
         use particle_types, only: particle
         use ppp,            only: ppp_main

         implicit none

         real,        intent(in) :: kdt
         type(particle), pointer :: pset
         character(len=*), parameter :: k_label = "part_kick"

         call ppp_main%start(k_label, PPP_PART)
         cgl => leaves%first
         do while (associated(cgl))
            call cgl%cg%costs%start
            pset => cgl%cg%pset%first
            do while (associated(pset))
               pset%pdata%vel = pset%pdata%vel + pset%pdata%acc * kdt
               pset => pset%nxt
            enddo
            call cgl%cg%costs%stop(I_PARTICLE)
            cgl => cgl%nxt
         enddo
         call ppp_main%stop(k_label, PPP_PART)

      end subroutine kick

      subroutine drift(ddt)

         use cg_cost_data,   only: I_PARTICLE
         use constants,      only: PPP_PART
         use particle_utils, only: part_leave_cg, is_part_in_cg, detach_particle
         use particle_types, only: particle
         use ppp,            only: ppp_main

         implicit none

         real,            intent(in) :: ddt
         type(particle), pointer     :: pset
         character(len=*), parameter :: d_label = "part_drift"

         call ppp_main%start(d_label, PPP_PART)
         cgl => leaves%first
         do while (associated(cgl))
            call cgl%cg%costs%start
            pset => cgl%cg%pset%first
            do while (associated(pset))
               if (pset%pdata%phy .and. pset%pdata%fin) then
                  pset%pdata%pos = pset%pdata%pos + pset%pdata%vel * ddt
                  call pset%pdata%is_outside()
                  call is_part_in_cg(cgl%cg, pset%pdata%pos, .not.pset%pdata%outside, pset%pdata%in, pset%pdata%phy, pset%pdata%out, pset%pdata%fin)
                  pset => pset%nxt
               else !Remove ghosts
                  call detach_particle(cgl%cg, pset)
               endif
            enddo
            call cgl%cg%costs%stop(I_PARTICLE)
            cgl => cgl%nxt
         enddo
         call ppp_main%stop(d_label, PPP_PART)

         call part_leave_cg()

      end subroutine drift

   end subroutine leapfrog_2ord

   subroutine update_particle_kinetic_energy

      use cg_cost_data,   only: I_PARTICLE
      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: half
      use particle_types, only: particle

      implicit none

      type(cg_list_element), pointer :: cgl
      type(particle), pointer        :: pset
      real    :: v2 !< particle velocity squared

      cgl => leaves%first
      do while (associated(cgl))
         call cgl%cg%costs%start
         pset => cgl%cg%pset%first
         do while (associated(pset))
            v2 = sum(pset%pdata%vel(:)**2)
            !energy       = 1/2  *      m         *  v**2 +     Ep(x,y,z)
            pset%pdata%energy = half * pset%pdata%mass * v2 + pset%pdata%energy
            pset => pset%nxt
         enddo
         call cgl%cg%costs%stop(I_PARTICLE)
         cgl => cgl%nxt
      enddo

   end subroutine update_particle_kinetic_energy

end module particle_solvers
