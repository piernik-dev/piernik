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

module particle_integrators
! pulled by GRAV

   use particle_types, only: particle_solver_T

   implicit none

   private
   public :: hermit4
#ifdef NBODY
   public :: leapfrog2
#endif /* NBODY */

   type, extends(particle_solver_T) :: hermit4_T
   contains
      procedure, nopass :: evolve => hermit_4ord
   end type hermit4_T

#ifdef NBODY
   type, extends(particle_solver_T) :: leapfrog2_T
   contains
      procedure, nopass :: evolve => leapfrog2ord
   end type leapfrog2_T
#endif /* NBODY */

   type(hermit4_T), target   :: hermit4
#ifdef NBODY
   type(leapfrog2_T), target :: leapfrog2
#endif /* NBODY */

contains

   !>
   !!  from:   Moving Stars Around (Piet Hut and Jun Makino, 2003).
   !!  http://www.ids.ias.edu/~piet/act/comp/algorithms/starter/index.html
   !<

   subroutine hermit_4ord(pset, t_glob, dt_tot, forward)

      use constants, only: ndims, xdim, zdim
      use particle_types, only: particle_set

      implicit none

      class(particle_set), intent(inout) :: pset  !< particle list
      real,                intent(in)    :: t_glob, dt_tot
      logical, optional,   intent(in)    :: forward

      real, parameter :: dt_param = 0.03        ! control parameter to determine time step size
      real, parameter :: dt_dia = 1             ! time interval between diagnostics output
      real, parameter :: dt_out = 0.01          ! time interval between output of snapshots

      real, dimension(:),    allocatable :: mass
      real, dimension(:, :), allocatable :: pos, vel, acc, jerk

      real    :: epot, coll_time
      real    :: t_dia, t_out, t_end, einit, dt, t
      integer :: nsteps, n, ndim, lun_out, lun_err, i

      open(newunit=lun_out, file='nbody_out.log', status='unknown',  position='append')
      open(newunit=lun_err, file='nbody_err.log', status='unknown',  position='append')

      n = size(pset%p, dim=1)
      t = t_glob
      allocate(mass(n), pos(n, ndims), vel(n, ndims), acc(n, ndims), jerk(n, ndims))

      mass(:) = pset%p(:)%mass
      do ndim = xdim, zdim
         pos(:, ndim) = pset%p(:)%pos(ndim)
         vel(:, ndim) = pset%p(:)%vel(ndim)
      enddo

      write(lun_err, *) "Starting a Hermite integration for a ", n, "-body system,"
      write(lun_err, *) "  from time t = ", t, " with time step control parameter dt_param = ", dt_param, "  until time ", t + dt_tot, " ,"
      write(lun_err, *) "  with diagnostics output interval dt_dia = ", dt_dia, ","
      write(lun_err, *) "  and snapshot output interval dt_out = ", dt_out, "."

      call get_acc_jerk_pot_coll(mass, pos, vel, acc, jerk, n, epot, coll_time)

      nsteps = 0

      call write_diagnostics(mass, pos, vel, acc, jerk, n, t, epot, nsteps, einit, .True.)

      t_dia = t + dt_dia  ! next time for diagnostics output
      t_out = t + dt_out  ! next time for snapshot output
      t_end = t + dt_tot  ! final time, to finish the integration

      do
         do while (t < t_dia .and. t < t_out .and. t < t_end)
            dt = dt_param * coll_time
            call evolve_step(mass, pos, vel, acc, jerk, n, t, dt, epot, coll_time)
            nsteps = nsteps + 1
         enddo
         if (t >= t_dia) then
            call write_diagnostics(mass, pos, vel, acc, jerk, n, t, epot, nsteps, einit, .False.)
            t_dia = t_dia + dt_dia
         endif
         if (t >= t_out) then
            write(lun_out, *) "#", t
            do i = 1, n
               write(lun_out, '(7(E13.6,1X))') mass(i), pos(i,:), vel(i,:)
            enddo
            t_out = t_out + dt_out
         endif
         if (t >= t_end) exit
      enddo

      do ndim = xdim, zdim
         pset%p(:)%pos(ndim) = pos(:, ndim)
         pset%p(:)%vel(ndim) = vel(:, ndim)
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
            write(lun_err, *) " E_kin = ", ekin, " , E_pot = ", epot, &
               & " , E_tot = ", etot
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

#ifdef NBODY
   subroutine leapfrog2ord(pset, t_glob, dt_tot, forward)

      use constants,        only: half, two
      use particle_gravity, only: update_particle_gravpot_and_acc
      use particle_types,   only: particle_set, twodtscheme

      implicit none

      class(particle_set), intent(inout) :: pset                 !< particle list
      real,                intent(in)    :: t_glob               !< initial time of simulation
      real,                intent(in)    :: dt_tot               !< timestep of simulation
      logical, optional,   intent(in)    :: forward

      real                               :: dt_kick              !< timestep for kicks
      real                               :: total_energy         !< total energy of set of the particles
      integer                            :: n                    !< number of particles

      n = size(pset%p, dim=1)
      if (twodtscheme) then
         dt_kick = dt_tot
      else
         dt_kick = dt_tot * half
      endif

      if (forward .or. .not.twodtscheme) then
         call kick                           (n,     dt_kick) !1. kick
      endif
      if (.not.forward .or. .not.twodtscheme) then
         call drift                          (n, two*dt_kick) !2. drift

         call update_particle_gravpot_and_acc

         call kick                           (n,     dt_kick) !3. kick

         call update_particle_kinetic_energy (n, total_energy)

         call leapfrog2_diagnostics          (n, total_energy)
      endif

      contains

         subroutine kick(n, kdt)

            implicit none

            integer, intent(in) :: n
            real,    intent(in) :: kdt
            integer             :: i

            do i = 1, n
               pset%p(i)%vel = pset%p(i)%vel + pset%p(i)%acc * kdt
            enddo

         end subroutine kick

         subroutine drift(n, ddt)

            implicit none

            integer, intent(in) :: n
            real,    intent(in) :: ddt
            integer             :: i

            do i = 1, n
               pset%p(i)%pos = pset%p(i)%pos + pset%p(i)%vel * ddt
            enddo

         end subroutine drift

         subroutine update_particle_kinetic_energy(n, total_energy)

            use constants, only: half, zero, xdim, zdim

            implicit none

            integer, intent(in)  :: n            !< number of particles
            real,    intent(out) :: total_energy !< total energy of set of particles
            integer(kind=4)      :: cdim
            integer              :: p
            real                 :: v2           !< particle velocity squared

            total_energy = zero

            do p = 1, n
               v2 = zero
               do cdim = xdim, zdim
                  v2 = v2 + pset%p(p)%vel(cdim)**2
               enddo
               !energy       = 1/2  *      m         *  v**2 +     Ep(x,y,z)
               pset%p(p)%energy = half * pset%p(p)%mass * v2 + pset%p(p)%energy
               total_energy = total_energy + pset%p(p)%energy
            enddo

         end subroutine update_particle_kinetic_energy

         subroutine leapfrog2_diagnostics(n, total_energy)

            use constants,        only: ndims
            use dataio_pub,       only: msg, printinfo
            !use func,             only: operator(.equals.)
            use particle_gravity, only: get_acc_model

            implicit none

            integer,    intent(in) :: n
            real,       intent(in) :: total_energy
            real, dimension(ndims) :: acc2
            real                   :: ang_momentum         !< angular momentum of set of the particles
            integer                :: i, lun_out
            !real                   :: d_energy             !< error of energy of set of the particles in succeeding timesteps
            !real                   :: d_ang_momentum       !< error of angular momentum in succeeding timensteps
            !real, save             :: initial_energy       !< total initial energy of set of the particles
            !real, save             :: init_ang_momentum    !< initial angular momentum of set of the particles
            !logical, save          :: first_run_lf = .true.

            call get_ang_momentum_2(n, ang_momentum)
            write(msg,'(a,2f8.5)') '[particle_integrators:leapfrog2_diagnostics] Energia, ANG_MOMENTUM----------: ', total_energy, ang_momentum
            call printinfo(msg)

            !if (first_run_lf) then
            !   initial_energy = total_energy
            !   d_energy       = 0.0
            !   init_ang_momentum = ang_momentum
            !   d_ang_momentum    = 0.0
            !else
            !   d_energy = log(abs((total_energy - initial_energy)/initial_energy))
            !   if (init_ang_momentum .equals. zero) then
            !      d_ang_momentum = ang_momentum
            !   else
            !      d_ang_momentum = (ang_momentum - init_ang_momentum)/init_ang_momentum
            !      if (d_ang_momentum .equals. zero) then
            !         d_ang_momentum = 0.0
            !      else
            !         d_ang_momentum = log(abs(d_ang_momentum))
            !      endif
            !   endif
            !endif

            open(newunit=lun_out, file='leapfrog_out.log', status='unknown',  position='append')

            do i = 1, n
               call get_acc_model(i, 0.0, acc2)
               write(lun_out, '(a,I3.3,1X,19(E13.6,1X))') 'particle', i, t_glob+dt_kick, dt_kick, pset%p(i)%mass, pset%p(i)%pos, pset%p(i)%vel, pset%p(i)%acc, acc2(:)!, pset%p(i)%energy, total_energy, initial_energy, d_energy, ang_momentum, init_ang_momentum, d_ang_momentum
            enddo

            close(lun_out)

         end subroutine leapfrog2_diagnostics

         subroutine get_ang_momentum_2(n, ang_momentum)

            use constants, only: xdim, ydim, zdim

            implicit none

            integer, intent(in)  :: n
            real,    intent(out) :: ang_momentum
            integer              :: i
            real                 :: L1, L2, L3

            ang_momentum = 0.0

            do i = 1, n, 1
               L1 = pset%p(i)%pos(ydim) * pset%p(i)%vel(zdim) - pset%p(i)%pos(zdim) * pset%p(i)%vel(ydim)
               L2 = pset%p(i)%pos(zdim) * pset%p(i)%vel(xdim) - pset%p(i)%pos(xdim) * pset%p(i)%vel(zdim)
               L3 = pset%p(i)%pos(xdim) * pset%p(i)%vel(ydim) - pset%p(i)%pos(ydim) * pset%p(i)%vel(xdim)
               ang_momentum = ang_momentum + pset%p(i)%mass * sqrt(L1**2 + L2**2 + L3**2)
            enddo

         end subroutine get_ang_momentum_2

   end subroutine leapfrog2ord

   subroutine get_acc_pot(mass, pos, acc, n, epot)

      use constants, only: ndims

      implicit none

      integer,                  intent(in)  :: n
      real, dimension(n),       intent(in)  :: mass
      real, dimension(n,ndims), intent(in)  :: pos
      real, dimension(n,ndims), intent(out) :: acc
      real,                     intent(out) :: epot

      integer                               :: i, j
      real, dimension(ndims)                :: rji, da
      real                                  :: r   ! | rji |
      real                                  :: r2  ! | rji |^2
      real                                  :: r3  ! | rji |^3

      acc(:,:) = 0.0
      epot = 0.0

      do i = 1, n
         do j = i+1, n
            rji(:) = pos(j, :) - pos(i, :)

            r2 = sum(rji**2)
            r = sqrt(r2)
            r3 = r * r2

            ! add the {i,j} contribution to the total potential energy for the system

            epot = epot - mass(i) * mass(j)/r

            da(:) = rji(:) / r3


            acc(i,:) = acc(i,:) + mass(j) * da(:)
            acc(j,:) = acc(j,:) - mass(i) * da(:)
         enddo
      enddo

   end subroutine get_acc_pot
#endif /* NBODY */

end module particle_integrators
