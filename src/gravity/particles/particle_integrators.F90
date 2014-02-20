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
   public :: hermit4, leapfrog2

   type, extends(particle_solver_T) :: hermit4_T
   contains
      procedure, nopass :: evolve => hermit_4ord
   end type hermit4_T
   
   type, extends(particle_solver_T) :: leapfrog2_T
   contains
      procedure, nopass :: evolve => leapfrog2ord
   end type leapfrog2_T

   type(hermit4_T), target :: hermit4
   type(leapfrog2_T), target :: leapfrog2

contains

   !>
   !!  from:   Moving Stars Around (Piet Hut and Jun Makino, 2003).
   !!  http://www.ids.ias.edu/~piet/act/comp/algorithms/starter/index.html
   !<

   subroutine hermit_4ord(pset, t_glob, dt_tot)
      use constants, only: ndims, xdim, zdim
      use particle_types, only: particle_set
      implicit none
      class(particle_set), intent(inout) :: pset  !< particle list
      real, intent(in) :: t_glob, dt_tot

      real, parameter :: dt_param = 0.03        ! control parameter to determine time step size
      real, parameter :: dt_dia = 1             ! time interval between diagnostics output
      real, parameter :: dt_out = 0.01          ! time interval between output of snapshots

      real, dimension(:), allocatable :: mass
      real, dimension(:, :), allocatable :: pos, vel, acc, jerk

      real :: epot, coll_time
      real :: t_dia, t_out, t_end, einit, dt, t
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

      contains

         subroutine write_diagnostics(mass, pos, vel, acc, jerk, n, t, epot, nsteps, einit, init_flag)
            implicit none
            integer, intent(in) :: n, nsteps
            real, intent(in) :: t, epot
            real, intent(inout) :: einit
            logical, intent(in) :: init_flag
            real, dimension(n), intent(in) :: mass
            real, dimension(n, ndims), intent(in) :: pos, vel, acc, jerk

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


   subroutine leapfrog2ord(pset, t_glob, dt_tot)
      use constants, only: ndims, xdim, zdim
      use particle_types, only: particle_set
      implicit none
      class(particle_set), intent(inout) :: pset  !< particle list
      real, intent(in) :: t_glob, dt_tot
      
      real, dimension(:), allocatable :: mass
      real, dimension(:, :), allocatable :: pos, vel, acc, vel_h, acc2
      
      real :: t_dia, t_out, t_end, einit, dt, t, dth, eta, eps, a
      integer :: nsteps, n, ndim, lun_out, lun_err, i
      
      open(newunit=lun_out, file='nbody_out.log', status='unknown',  position='append')
      
      n = size(pset%p, dim=1)
      t = t_glob
      allocate(mass(n), pos(n, ndims), vel(n, ndims), acc(n, ndims), vel_h(n, ndims), acc2(n,ndims))
      
      mass(:) = pset%p(:)%mass
      do ndim = xdim, zdim
         pos(:, ndim) = pset%p(:)%pos(ndim)
         vel(:, ndim) = pset%p(:)%vel(ndim)
      enddo
      
      
      t = 0.0
      eta = 35.0
      eps = 1.0e-5
      t_end=t+dt_tot
      
      call get_acc(mass, pos, vel, acc, n)
      
      call get_acc_mod(acc, n, a)
      
      
      !a tu trzeba policzyć przyspieszenia
      dt = sqrt(2.0*eta*eps/a)	
      dth = dt/2.0
      
      !main loop
      do while (t<t_end)
         !1.kick(dth)
         vel_h=vel
         call kick(vel_h, acc, dth, n) !velocity
         !2.drift(dt)
         call drift(pos,vel_h,dt,n) !position
         !3.acceleration + |a|
         call get_acc(mass, pos, vel, acc, n)
         call get_acc_mod(acc, n, a)
         !4.kick(dth)
         call kick(vel, acc, dth, n)   !velocity
         !5.t=t+dt
         t = t + dt
         !6.dt=sqrt(2.0*eta*eps/a)		!dt[n+1]
         dt	= sqrt(2.0*eta*eps/a)
         !7.dth=dt/2.0
			dth =	0.5*dt
         
      end do
      
      do ndim = xdim, zdim
         pset%p(:)%pos(ndim) = pos(:, ndim)
         pset%p(:)%vel(ndim) = vel(:, ndim)
      enddo

      deallocate (mass, pos, vel, acc)
      close(lun_out)
      
      contains
         
         !Kick
         subroutine kick(vel, acc, t, n)
            use constants, only: ndims
            implicit none
            real, intent(in) :: t
            integer, intent(in) :: n
            real, dimension(n, ndims), intent(in) :: acc
            real, dimension(n, ndims), intent(out) :: vel
		    
            vel = vel + acc*t
            
         end subroutine kick

         !Drift
         subroutine drift(pos, vel, t, n)
            use constants, only: ndims
            implicit none
            real, intent(in) :: t
            integer, intent(in) :: n
            real, dimension(n, ndims), intent(out) :: pos
            real, dimension(n, ndims), intent(in) :: vel
           
            pos = pos + vel*t
         end subroutine drift
      
   end subroutine leapfrog2ord


   subroutine get_acc(mass, pos, vel, acc, n)
      use constants, only: ndims
      implicit none
      integer, intent(in) :: n
      real, dimension(n), intent(in) :: mass
      real, dimension(n,ndims), intent(in) :: pos
      real, dimension(n,ndims), intent(in) :: vel
      real, dimension(n,ndims), intent(out) :: acc
      
      integer :: i, j
      real, dimension(ndims) :: rji, vji, da
      
      real :: r   ! | rji |
      real :: r2  ! | rji |^2
      real :: r3  ! | rji |^3
      
      acc(:,:) = 0.0
      
      do i = 1, n
         do j = i+1, ndims
            rji(:) = pos(j, :) - pos(i, :)
            vji(:) = vel(j, :) - vel(i, :)

            r2 = sum(rji**2)
            
            r = sqrt(r2)
            r3 = r * r2

            ! add the {i,j} contribution to the total potential energy for the
            ! system
            !epot = epot - mass(i) * mass(j) / r

            da(:) = rji(:) / r3
            

            acc(i,:) = acc(i,:) + mass(j) * da(:)
            acc(j,:) = acc(j,:) - mass(i) * da(:)
         enddo
      enddo
   end subroutine get_acc

   subroutine get_acc_mod(acc, n, a)
      use constants, only: ndims
      implicit none
      integer, intent(in) :: n
      integer  :: i, j
      real, dimension(n, ndims), intent(in) :: acc
      real, dimension(n) :: acc2
      real, intent(out)  :: a
      
      acc2 = 0.0
      
      do i = 1, n
         do j = 1, ndims
            acc2(i) = acc2(i) + acc(i,j)**2
         enddo
      enddo
      a = sqrt(maxval(acc2))
         
   end subroutine get_acc_mod

   subroutine evolve_step(mass, pos, vel, acc, jerk, n, t, dt, epot, coll_time)
      use constants, only: ndims
      implicit none
      integer, intent(in) :: n
      real, dimension(n), intent(in) :: mass
      real, dimension(n, ndims), intent(out) :: vel, pos
      real, intent(in) :: dt
      real, intent(inout) :: t, epot, coll_time
      real, dimension(n, ndims), intent(inout) :: acc, jerk

      real, dimension(n, ndims) :: old_pos, old_vel, old_acc, old_jerk

      old_pos = pos
      old_vel = vel
      old_acc = acc
      old_jerk = jerk

      call predict_step(pos, vel, acc, jerk, n, dt)
      call get_acc_jerk_pot_coll(mass, pos, vel, acc, jerk, n, epot, coll_time)
      call correct_step(pos, vel, acc, jerk, old_pos, old_vel, old_acc, old_jerk, n, dt);

      t = t + dt

      contains

         subroutine predict_step(pos, vel, acc, jerk, n, dt)
            implicit none
            integer, intent(in) :: n
            real, intent(in) :: dt
            real, dimension(n, ndims), intent(in) :: acc, jerk
            real, dimension(n, ndims), intent(out) :: vel, pos

            real :: hdt, hdt2
            real, parameter :: sixth = 1./6.

            hdt = 0.5*dt**2
            hdt2 = dt**3 * sixth

            pos(:,:) = pos(:,:) + vel(:,:)*dt + acc(:,:)*hdt + jerk(:,:)*hdt2
            vel(:,:) = vel(:,:) + acc(:,:)*dt + jerk(:,:)*hdt
         end subroutine predict_step

         subroutine correct_step(pos, vel, acc, jerk, old_pos, old_vel, old_acc, old_jerk, n, dt)
            implicit none
            integer, intent(in) :: n
            real, intent(in) :: dt
            real, dimension(n, ndims), intent(in) :: acc, jerk, old_pos, old_vel, old_acc, old_jerk
            real, dimension(n, ndims), intent(out) :: vel, pos

            real :: hdt, hdt2
            real, parameter :: third = 1./3.

            hdt = 0.5*dt
            hdt2 = hdt**2 * third

            vel(:,:) = old_vel(:,:) + (old_acc(:,:) + acc(:,:))*hdt + (old_jerk(:,:) - jerk(:,:)) * hdt2
            pos(:,:) = old_pos(:,:) + (old_vel(:,:) + vel(:,:))*hdt + (old_acc(:,:)  - acc(:,:) ) * hdt2
         end subroutine correct_step

   end subroutine evolve_step

   subroutine get_acc_jerk_pot_coll(mass, pos, vel, acc, jerk, n, epot, coll_time)
      use constants, only: ndims
      implicit none
      integer, intent(in) :: n
      real, dimension(n), intent(in) :: mass
      real, dimension(n,ndims), intent(in) :: pos
      real, dimension(n,ndims), intent(in) :: vel
      real, dimension(n,ndims), intent(out) :: acc
      real, dimension(n,ndims), intent(out) :: jerk
      real, intent(out) :: epot
      real, intent(out) :: coll_time

      real, dimension(ndims) :: rji, vji, da, dj
      integer :: i, j
      real :: coll_time_q  ! collision time to 4th power

      real :: r   ! | rji |
      real :: r2  ! | rji |^2
      real :: r3  ! | rji |^3
      real :: v2  ! | vji |^2
      real :: rv_r2 ! ( rij . vij ) / | rji |^2
      real :: da2

      acc(:,:) = 0.0
      jerk(:,:) = 0.0
      epot = 0.0

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

end module particle_integrators
