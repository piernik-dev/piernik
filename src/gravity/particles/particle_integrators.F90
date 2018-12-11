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
   public :: leapfrog2, lf_c, timestep_nbody, dt_nbody, is_setacc_cic, is_setacc_int
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
   real                      :: lf_c               !< timestep should depends of grid and velocities of particles (used to extrapolation of the gravitational potential)
   real                      :: dt_nbody           !< timestep depends on particles
   logical                   :: is_setacc_cic, is_setacc_int
#endif /* NBODY */

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
      real,                intent(in)    :: t_glob, dt_tot

#ifdef NBODY
      real, parameter :: dt_param = 0.0001      ! control parameter to determine time step size
#else /* !NBODY */
      real, parameter :: dt_param = 0.03        ! control parameter to determine time step size
#endif /* !NBODY */
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
#ifdef NBODY
            dt = dt_param              !constant timestep
#else /* !NBODY */
            !variable timestep
            dt = dt_param * coll_time
#endif /* !NBODY */
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
#ifdef NBODY
         print *, "Hermit dt=", dt
#endif /* NBODY */
      enddo
#ifdef NBODY
      print *, "Hermit nsteps=", nsteps
#endif /* NBODY */

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

#ifdef NBODY
   subroutine leapfrog2ord(pset, t_glob, dt_tot)

      use cg_list,        only: cg_list_element
      use constants,      only: ndims
      use dataio_pub,     only: die, msg, printinfo
      use domain,         only: is_refined, is_multicg
      use func,           only: operator(.equals.)
      use global,         only: dt_old
      use grid_cont,      only: grid_container
      use particle_types, only: particle_set

      implicit none

      class(particle_set), intent(inout) :: pset                 !< particle list
      real,                intent(in)    :: t_glob               !< initial time of simulation
      real,                intent(in)    :: dt_tot               !< timestep of simulation

      real, dimension(:), allocatable    :: mass                 !< 1D array of mass of the particles
      real                               :: dt_tot_h             !< half of timestep, dt_tot_h = 0.5*dt_tot
      real                               :: total_energy         !< total energy of set of the particles
      !real, save                         :: initial_energy       !< total initial energy of set of the particles
      !real                               :: d_energy             !< error of energy of set of the particles in succeeding timesteps
      real                               :: ang_momentum         !< angular momentum of set of the particles
      !real, save                         :: init_ang_momentum    !< initial angular momentum of set of the particles
      !real                               :: d_ang_momentum       !< error of angular momentum in succeeding timensteps
      integer                            :: i
      integer                            :: n                    !< number of particles
      logical                            :: external_pot         !< if .true. gravitational potential will be deleted and replaced by external potential of point mass
      logical, save                      :: first_run_lf = .true.
      integer, save                      :: counter
      integer                            :: lun_out
      real, dimension(:,:), allocatable  :: acc2

      open(newunit=lun_out, file='leapfrog_out.log', status='unknown',  position='append')

      n = size(pset%p, dim=1)

      allocate(mass(n))
      allocate(acc2(n, ndims))

      if (t_glob < 0.0) i=1 ! supress compiler warnings

      mass(:) = pset%p(:)%mass

      !cgl => leaves%first
      if (is_refined) call die("[particle_integrators:leapfrog2ord] AMR not implemented for particles yet")
      if (is_multicg) call die("[particle_integrators:leapfrog2ord] multi_cg not implemented for particles yet")

      !do while (associated(cgl))
         !cg => cgl%cg
         !cgl => cgl%nxt
      !enddo

      !obliczenie zewnętrznego potencjalu na siatce
      !external_pot = .true.
      external_pot = .false.

      call get_energy(pset, total_energy, n)
      !write(*,*) "Energia----------: ", total_energy

      call get_ang_momentum_2(pset, n, ang_momentum)
      !write(*,*) "ANG_MOMENTUM-----: ", ang_momentum
      write(msg,'(a,2f8.5)') '[particle_integrators:leapfrog2ord] Energia, ANG_MOMENTUM----------: ', total_energy, ang_momentum
      call printinfo(msg)
      !stop

      !if (first_run_lf) then
      !   initial_energy = total_energy
      !  d_energy       = 0.0
      !   init_ang_momentum = ang_momentum
      !   d_ang_momentum    = 0.0
      !   write(*,*) "Pierwszy"
      !else
      !   d_energy = log(abs((total_energy - initial_energy)/initial_energy))
      !   if (init_ang_momentum.equals. zero) then
      !      d_ang_momentum = ang_momentum
      !   else
      !      d_ang_momentum = (ang_momentum - init_ang_momentum)/init_ang_momentum
      !      if (d_ang_momentum.equals. zero) then
      !         d_ang_momentum = 0.0
      !      else
      !         d_ang_momentum = log(abs(d_ang_momentum))
      !      endif
      !   endif
      !endif

      write(msg,'(a,3f8.5)') '[particle_integrators:leapfrog2ord] pset%p(1)%pos = ', pset%p(1)%pos
      call printinfo(msg)
      write(msg,'(a,3f8.5)') '[particle_integrators:leapfrog2ord] pset%p(1)%vel = ', pset%p(1)%vel
      call printinfo(msg)
      !write(*,*) "1:", pset%p(1)%vel
      !write(*,*) "2:", pset%p(2)%vel

      counter = 1

      do i = 1, n
      !write(*,*) "n=",n," i=",i
         call get_acc_model(pset, acc2, 0.0, n)
         write(lun_out, '(I3,1X,19(E13.6,1X))') i, t_glob+dt_tot, dt_old, mass(i), pset%p(i)%pos, pset%p(i)%vel, pset%p(i)%acc, acc2(i,:)!, pset%p(i)%energy, total_energy, initial_energy, d_energy, ang_momentum, init_ang_momentum, d_ang_momentum
         !write(lun_out, '(I3,1X,19(E13.6,1X))') i, t_glob+dt_tot, dt_old, pset%p(i)%pos, pset%p(i)%vel, pset%p(i)%acc, acc2(i,:)!, pset%p(i)%energy, total_energy, initial_energy, d_energy, ang_momentum, init_ang_momentum, d_ang_momentum
      enddo

      !call save_particles(n, lf_t, mass, pset, counter)
      !write(*,*) "[p_i]:-----------------------------"
      !write(*,*) "[p_i]:dt_tot= ", dt_tot
      !write(*,*) "[p_i]:dt_old= ", dt_old
      write(msg,'(a,2f8.5)') '[particle_integrators:leapfrog2ord] [p_i]: dt_tot, dt_old = ', dt_tot, dt_old
      call printinfo(msg)
      dt_tot_h = 0.5*dt_tot

      if(first_run_lf) then
         first_run_lf = .false.
      else
         !3.kick(dt_old)
         call kick(pset, 0.5*dt_old, n)
      endif
      !1. Kick (dt_tot_h)
      call kick(pset, dt_tot_h, n)

      !2.drift(lf_dt)
      call drift(pset, dt_tot, n)
      !stop

      !call save_particles(n, lf_t, mass, pset, counter)

      close(lun_out)
      !write(*,*) "[p_i]:-----------------------------"

      contains

         !Kick
         subroutine kick(pset, t, n)

            use particle_types, only: particle_set

            implicit none

            class(particle_set), intent(inout) :: pset  !< particle list
            real,                intent(in)    :: t
            integer,             intent(in)    :: n
            integer                            :: i

            do i = 1, n
               !pset%p(i)%vel = pset%p(i)%vel + pset%p(i)%acc * t
               pset%p(i)%vel = pset%p(i)%vel + 0.0 * t
            enddo

         end subroutine kick

         !Drift
         subroutine drift(pset, t, n)

            use particle_types, only: particle_set

            implicit none

            class(particle_set), intent(inout) :: pset  !< particle list
            real,                intent(in)    :: t
            integer,             intent(in)    :: n
            integer                            :: i

            do i=1, n
               pset%p(i)%pos = pset%p(i)%pos + pset%p(i)%vel * t
            enddo

         end subroutine drift

         subroutine get_energy(pset, total_energy, n)

            use constants,      only: ndims, half, zero
            use particle_types, only: particle_set

            implicit none

            class(particle_set), intent(inout) :: pset         !< particle list
            real,                intent(out)   :: total_energy !< total energy of set of particles
            integer,             intent(in)    :: n            !< number of particles
            integer                            :: cdim, p
            real, dimension(n)                 :: v            !< kwadraty predkosci czastek

            v = zero
            total_energy = zero

            do p = 1, n
               do cdim = 1, ndims
                  v(p) = v(p)+ pset%p(p)%vel(cdim)**2
               enddo
               !energy       = 1/2  *      m         *  v**2 +     Ep(x,y,z)
               pset%p(p)%energy = half * pset%p(p)%mass *  v(p) + pset%p(p)%energy
               total_energy = total_energy + pset%p(p)%energy
            enddo

         end subroutine get_energy

         subroutine save_particles(n, lf_t, mass, pset, counter)

            use particle_types, only: particle_set

            implicit none

            class(particle_set), intent(in)    :: pset  !< particle list
            integer,             intent(in)    :: n
            real,                intent(in)    :: lf_t
            real, dimension(n),  intent(in)    :: mass
            integer,             intent(inout) :: counter
            integer                            :: i, data_file = 757
            character(len=17)                  :: filename
            character(len=3)                   :: counter_char

            if(counter<10) then
               write(counter_char, '(I1)') counter
               write(filename,'(A9,A3,A1,A4)') 'particles','_00',counter_char,".dat"
            endif
            if ((counter >=10) .and.(counter <100)) then
               write(counter_char, '(I2)') counter
               write(filename,'(A9,A2,A2,A4)') 'particles','_0',counter_char,".dat"
            endif
            if (counter >=100) then
               write(counter_char, '(I3)') counter
               write(filename,'(A9,A1,A3,A4)') 'particles','_',counter_char,".dat"
            endif

            open(unit = data_file, file=filename)
               write(data_file, *) "#t =", lf_t
               do i = 1, n
                  write(data_file,*) i, mass(i), pset%p(i)%pos, pset%p(i)%vel
               enddo
            close(data_file)

            counter = counter + 1

         end subroutine save_particles

         subroutine get_acc_model(pset, acc2, eps, n)

            use constants, only: ndims, xdim, zdim
            use grid_cont, only: grid_container

            implicit none

            class(particle_set),      intent(in)  :: pset  !< particle list
            integer,                  intent(in)  :: n
            real,                     intent(in)  :: eps
            real, dimension(n,ndims), intent(out) :: acc2
            integer                               :: p
            integer(kind=4)                       :: dir

            do p = 1, n
               do dir = xdim, zdim
               acc2(p, dir) = -der_xyz(pset%p(p)%pos, 1.0e-8, eps, dir)
               enddo
            enddo

         end subroutine get_acc_model

      function der_xyz(pos, d, eps, dir)

         use constants, only: idm, ndims

         implicit none

         real(kind=8),dimension(1,ndims), intent(in) :: pos
         real(kind=8),                    intent(in) :: d, eps
         integer(kind=4),                 intent(in) :: dir
         real(kind=8)                                :: der_xyz

         der_xyz = ( phi_pm_part(pos(1,:)+real(idm(dir,:))*d, eps, 10.0) - phi_pm_part(pos(1,:)-real(idm(dir,:))*d, eps, 10.0) ) / (2.0*d)

      end function der_xyz

      subroutine get_ang_momentum_2(pset, n, ang_momentum)

         use constants,      only: xdim, ydim, zdim
         use particle_types, only: particle_set

         implicit none

         class(particle_set), intent(in)  :: pset
         integer,             intent(in)  :: n
         real,                intent(out) :: ang_momentum
         integer                          :: i
         real                             :: L1, L2, L3

         ang_momentum = 0.0

         do i = 1, n, 1
            L1 = pset%p(i)%pos(ydim) * pset%p(i)%vel(zdim) - pset%p(i)%pos(zdim) * pset%p(i)%vel(ydim)
            L2 = pset%p(i)%pos(zdim) * pset%p(i)%vel(xdim) - pset%p(i)%pos(xdim) * pset%p(i)%vel(zdim)
            L3 = pset%p(i)%pos(xdim) * pset%p(i)%vel(ydim) - pset%p(i)%pos(ydim) * pset%p(i)%vel(xdim)
            ang_momentum = ang_momentum + pset%p(i)%mass * sqrt(L1**2 + L2**2 + L3**2)
         enddo

      end subroutine get_ang_momentum_2

   end subroutine leapfrog2ord

   function phi_pm_part(pos, eps, mass)

      use constants, only: ndims, xdim, ydim, zdim
      use units,     only: newtong

      implicit none

      real, dimension(ndims), intent(in) :: pos
      real,                   intent(in) :: eps, mass
      real                               :: r, phi_pm_part

      r = sqrt(pos(xdim)**2 + pos(ydim)**2 + pos(zdim)**2 + eps**2)
      phi_pm_part = -newtong*mass / r

   end function phi_pm_part

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

            ! add the {i,j} contribution to the total potential energy for the
            ! system

            epot = epot - mass(i) * mass(j)/r

            da(:) = rji(:) / r3


            acc(i,:) = acc(i,:) + mass(j) * da(:)
            acc(j,:) = acc(j,:) - mass(i) * da(:)
         enddo
      enddo

   end subroutine get_acc_pot
#endif /* NBODY */

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
   subroutine timestep_nbody(dt_nbody, pset)

      use constants,      only: ndims, zero
      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use grid_cont,      only: grid_container
      use particle_types, only: particle_set
#ifdef VERBOSE
      use dataio_pub,     only: printinfo
#endif /* VERBOSE */

      implicit none

      class(particle_set),   intent(inout) :: pset
      type(grid_container),  pointer       :: cg
      type(cg_list_element), pointer       :: cgl

      integer                              :: n_part
      real                                 :: max_acc

      real,    dimension(:,:), allocatable :: dist
      integer, dimension(:,:), allocatable :: cells
      real                                 :: dt_nbody
      integer, save                        :: kroki = 0

#ifdef VERBOSE
      call printinfo('[particle_integrators:timestep_nbody] Commencing timestep_nbody')
#endif /* VERBOSE */

      n_part = size(pset%p, dim=1)
      !write(msg,'(a,i6)') '[particle_integrators:timestep_nbody] Number of particles: ', n_part
      !call printinfo(msg)
      allocate(cells(n_part, ndims), dist(n_part, ndims))

      cgl => leaves%first
      cg  => cgl%cg

      !call pot_from_part(cg, zero, n_part, pset)
      call pot2grid(cg, zero)

      call save_pot_pset(cg, pset)
      kroki = kroki + 1
      !write(*,*) "++++++++KROKI+++++++: ", kroki

      call find_cells(n_part, cg, cells, dist)

      call potential2(n_part, cg, cells, dist, pset)    !szukanie energii potencjalnej w punktach-polozeniach czastek

      if (is_setacc_int) then
         call get_acc_int(n_part, cg, cells, dist, pset)
      elseif (is_setacc_cic) then
         call get_acc_cic(n_part, cg, cells, pset)
      endif

      call get_acc_max(n_part, pset, max_acc)

      call get_var_timestep_c(dt_nbody, max_acc, lf_c, pset, cg)
#ifdef VERBOSE
      call printinfo('[particle_integrators:timestep_nbody] Finish timestep_nbody')
#endif /* VERBOSE */

   contains

      subroutine pot2grid(cg, eps2)

         use constants, only: xdim, ydim, zdim, CENTER
         use grid_cont, only: grid_container

         implicit none

         type(grid_container), pointer, intent(inout) :: cg
         real,                          intent(in)    :: eps2
         integer                                      :: i, j, k

         do i = lbound(cg%gpot, dim=1), ubound(cg%gpot, dim=1)
            do j = lbound(cg%gpot, dim=2), ubound(cg%gpot, dim=2)
               do k = lbound(cg%gpot, dim=3), ubound(cg%gpot, dim=3)
                  cg%gpot(i,j,k) = phi_pm_part([cg%coord(CENTER,xdim)%r(i), &
                                                cg%coord(CENTER,ydim)%r(j), &
                                                cg%coord(CENTER,zdim)%r(k)], eps2, 1.0)
               enddo
            enddo
         enddo

      end subroutine pot2grid

      subroutine pot_from_part(cg, eps2, n_part, pset)

         use constants,      only: xdim, ydim, zdim, CENTER
         use grid_cont,      only: grid_container
         use particle_types, only: particle_set

         implicit none

         class(particle_set)                          :: pset  !< particle list
         type(grid_container), pointer, intent(inout) :: cg
         real,                          intent(in)    :: eps2
         integer,                       intent(in)    :: n_part
         integer                                      :: i, j, k, p

         cg%gpot = 0.0

         open(unit=999, file='pset.dat')
            do p = 1, n_part
               do i = lbound(cg%gpot, dim=1), ubound(cg%gpot, dim=1)
                  do j = lbound(cg%gpot, dim=2), ubound(cg%gpot, dim=2)
                     do k = lbound(cg%gpot, dim=3), ubound(cg%gpot, dim=3)
                        cg%gpot(i,j,k) = cg%gpot(i,j,k) + phi_pm_part([cg%coord(CENTER,xdim)%r(i) - pset%p(p)%pos(xdim), &
                                                                       cg%coord(CENTER,ydim)%r(j) - pset%p(p)%pos(ydim), &
                                                                       cg%coord(CENTER,zdim)%r(k) - pset%p(p)%pos(zdim)], eps2, pset%p(p)%mass)
                     enddo
                  enddo
               enddo
               write(999,*) p, pset%p(p)%pos
            enddo
         close(999)

      end subroutine pot_from_part

      subroutine find_cells(n_part, cg, cells, dist)

         use constants, only: ndims, xdim, CENTER, LO, HI
         use grid_cont, only: grid_container

         implicit none

         integer,                          intent(in)  :: n_part
         type(grid_container),             intent(in)  :: cg
         integer, dimension(n_part,ndims), intent(out) :: cells
         real,    dimension(n_part,ndims), intent(out) :: dist
         integer                                       :: i, cdim

         do i = 1, n_part
            do cdim = xdim, ndims
               if ((pset%p(i)%pos(cdim) >= cg%ijkse(cdim, LO)) .or. (pset%p(i)%pos(cdim) <= cg%ijkse(cdim, HI))) then
                  pset%p(i)%outside = .false.
               else
                  pset%p(i)%outside = .true.
               endif

               cells(i, cdim) = int( 0.5 + (pset%p(i)%pos(cdim) - cg%coord(CENTER,cdim)%r(0)) * cg%idl(cdim) )

               dist(i, cdim)  = pset%p(i)%pos(cdim) - ( cg%coord(CENTER, cdim)%r(0) + cells(i,cdim) * cg%dl(cdim) )
            enddo
         enddo

      end subroutine find_cells

      subroutine potential2(n_part, cg, cells, dist, pset)

         use constants,        only: gpot_n, ndims, half, xdim, ydim, zdim
         use grid_cont,        only: grid_container
         use named_array_list, only: qna
         use particle_func,    only: df_d_o2, d2f_d2_o2, d2f_dd_o2

         implicit none

         integer,                          intent(in)    :: n_part
         type(grid_container), pointer,    intent(in)    :: cg
         integer, dimension(n_part,ndims), intent(in)    :: cells
         real,    dimension(n_part,ndims), intent(in)    :: dist
         class(particle_set),              intent(inout) :: pset  !< particle list
         integer                                         :: p
         integer(kind=4)                                 :: ig
         real, dimension(n_part)                         :: dpot, d2pot

         ig = qna%ind(gpot_n)
         do p = 1, n_part
            dpot(p) = df_d_o2([cells(p, :)], cg, ig, xdim) * dist(p, xdim) + &
                      df_d_o2([cells(p, :)], cg, ig, ydim) * dist(p, ydim) + &
                      df_d_o2([cells(p, :)], cg, ig, zdim) * dist(p, zdim)

            d2pot(p) = d2f_d2_o2([cells(p, :)], cg, ig, xdim) * dist(p, xdim)**2 + &
                       d2f_d2_o2([cells(p, :)], cg, ig, ydim) * dist(p, ydim)**2 + &
                       d2f_d2_o2([cells(p, :)], cg, ig, zdim) * dist(p, zdim)**2 + &
                       2.0*d2f_dd_o2([cells(p, :)], cg, ig, xdim, ydim) * dist(p, xdim)*dist(p, ydim) + &
                       2.0*d2f_dd_o2([cells(p, :)], cg, ig, xdim, zdim) * dist(p, xdim)*dist(p, zdim)
            pset%p(p)%energy = cg%q(ig)%point(cells(p,:)) + dpot(p) + half * d2pot(p)
         enddo

      end subroutine potential2

      subroutine get_acc_int(n_part, cg, cells, dist, pset)

         use constants,        only: gpot_n, ndims, xdim, ydim, zdim
         use grid_cont,        only: grid_container
         use named_array_list, only: qna
         use particle_func,    only: df_d_p, d2f_d2_p, d2f_dd_p
         use particle_types,   only: particle_set

         implicit none

         integer,                          intent(in)    :: n_part
         type(grid_container), pointer,    intent(in)    :: cg
         integer, dimension(n_part,ndims), intent(in)    :: cells
         real, dimension(n_part, ndims),   intent(in)    :: dist
         class(particle_set),              intent(inout) :: pset
         integer                                         :: i
         integer(kind=4)                                 :: ig

         ig = qna%ind(gpot_n)
         do i = 1, n_part
            if ( (pset%p(i)%outside) .eqv. .false.) then
               pset%p(i)%acc(xdim) = - (df_d_p([cells(i, :)], cg, ig, xdim) + &
                                      d2f_d2_p([cells(i, :)], cg, ig, xdim)       * dist(i, xdim) + &
                                      d2f_dd_p([cells(i, :)], cg, ig, xdim, ydim) * dist(i, ydim) + &
                                      d2f_dd_p([cells(i, :)], cg, ig, xdim, zdim) * dist(i, zdim))

               pset%p(i)%acc(ydim) = -( df_d_p([cells(i, :)], cg, ig, ydim) + &
                                      d2f_d2_p([cells(i, :)], cg, ig, ydim)       * dist(i, ydim) + &
                                      d2f_dd_p([cells(i, :)], cg, ig, xdim, ydim) * dist(i, xdim) + &
                                      d2f_dd_p([cells(i, :)], cg, ig, ydim, zdim) * dist(i, zdim))

               pset%p(i)%acc(zdim) = -( df_d_p([cells(i, :)], cg, ig, zdim) + &
                                      d2f_d2_p([cells(i, :)], cg, ig, zdim)       * dist(i, zdim) + &
                                      d2f_dd_p([cells(i, :)], cg, ig, xdim, zdim) * dist(i, xdim) + &
                                      d2f_dd_p([cells(i, :)], cg, ig, ydim, zdim) * dist(i, ydim))
            !else
            !   call !funkcja liczaca pochodne z potencjalu policzonego z rozwiniecia multipolowego
            endif
         enddo

      end subroutine get_acc_int

      subroutine get_acc_max(n_part, pset, max_acc)

         use constants,      only: xdim, ndims
         use particle_types, only: particle_set

         implicit none

         integer,             intent(in)    :: n_part
         class(particle_set), intent(inout) :: pset
         real,                intent(out)   :: max_acc
         integer                            :: i, cdim
         real, dimension(n_part)            :: acc

         acc = 0.0

         do i = 1, n_part
            do cdim = xdim, ndims
               acc(i) = acc(i) + pset%p(i)%acc(cdim)**2
            enddo
         enddo

         max_acc = sqrt(maxval(acc))

      end subroutine get_acc_max

      subroutine get_acc_cic(n_part, cg, cells, pset)

         use constants,      only: ndims, CENTER, xdim, ydim, zdim, half, zero
         use grid_cont,      only: grid_container
         use particle_types, only: particle_set

         implicit none

         integer,                          intent(in)    :: n_part
         type(grid_container), pointer,    intent(in)    :: cg
         integer, dimension(n_part,ndims), intent(in)    :: cells
         class(particle_set),              intent(inout) :: pset

         integer                                         :: i, j, k, c, cdim
         integer                                         :: p
         integer(kind=8), dimension(n_part,ndims)        :: cic_cells
         real,            dimension(n_part,ndims)        :: dxyz
         real(kind=8),    dimension(n_part,8)            :: wijk, fx, fy, fz

         !write(*,*) "[get_acc_cic]: particles = ", n_part
         do i = 1, n_part
            pset%p(i)%acc = zero
            if ((pset%p(i)%outside) .eqv. .false.) then
               do cdim = xdim, ndims
                  if (pset%p(i)%pos(cdim) < cg%coord(CENTER, cdim)%r(cells(i,cdim))) then
                     cic_cells(i, cdim) = cells(i, cdim) - 1
                  else
                     cic_cells(i, cdim) = cells(i, cdim)
                  endif
                  dxyz(i, cdim) = abs(pset%p(i)%pos(cdim) - cg%coord(CENTER, cdim)%r(cic_cells(i,cdim)))

               enddo

               wijk(i, 1) = (cg%dx - dxyz(i, xdim))*(cg%dy - dxyz(i, ydim))*(cg%dz - dxyz(i, zdim)) !a(i  ,j  ,k  )
               wijk(i, 2) = (cg%dx - dxyz(i, xdim))*(cg%dy - dxyz(i, ydim))*         dxyz(i, zdim)  !a(i+1,j  ,k  )
               wijk(i, 3) = (cg%dx - dxyz(i, xdim))*         dxyz(i, ydim) *(cg%dz - dxyz(i, zdim)) !a(i  ,j+1,k  )
               wijk(i, 4) = (cg%dx - dxyz(i, xdim))*         dxyz(i, ydim) *         dxyz(i, zdim)  !a(i  ,j  ,k+1)
               wijk(i, 5) =          dxyz(i, xdim) *(cg%dy - dxyz(i, ydim))*(cg%dz - dxyz(i, zdim)) !a(i+1,j+1,k  )
               wijk(i, 6) =          dxyz(i, xdim) *(cg%dy - dxyz(i, ydim))*         dxyz(i, zdim)  !a(i  ,j+1,k+1)
               wijk(i, 7) =          dxyz(i, xdim) *         dxyz(i, ydim) *(cg%dz - dxyz(i, zdim)) !a(i+1,j  ,k+1)
               wijk(i, 8) =          dxyz(i, xdim) *         dxyz(i, ydim) *         dxyz(i, zdim)  !a(i+1,j+1,k+1)
            !else multipole expansion for particles outside domain
            endif

         enddo

         wijk = wijk/cg%dvol

         do p = 1, n_part
            c = 1
            do i = 0, 1
               do j = 0, 1
                  do k = 0, 1
                     fx(p, c) = -(cg%gpot(cic_cells(p, xdim)+1+i, cic_cells(p, ydim)  +j, cic_cells(p, zdim)  +k) - cg%gpot(cic_cells(p, xdim)-1+i, cic_cells(p, ydim)  +j, cic_cells(p, zdim)  +k))
                     fy(p, c) = -(cg%gpot(cic_cells(p, xdim)  +i, cic_cells(p, ydim)+1+j, cic_cells(p, zdim)  +k) - cg%gpot(cic_cells(p, xdim)  +i, cic_cells(p, ydim)-1+j, cic_cells(p, zdim)  +k))
                     fz(p, c) = -(cg%gpot(cic_cells(p, xdim)  +i, cic_cells(p, ydim)  +j, cic_cells(p, zdim)+1+k) - cg%gpot(cic_cells(p, xdim)  +i, cic_cells(p, ydim)  +j, cic_cells(p, zdim)-1+k))
                     c = c + 1
                  enddo
               enddo
            enddo
         enddo

         fx = half*fx*cg%idx
         fy = half*fy*cg%idy
         fz = half*fz*cg%idz

         do p = 1, n_part
            do c = 1, 8
               pset%p(p)%acc(xdim) = pset%p(p)%acc(xdim) + wijk(p, c) * fx(p, c)
               pset%p(p)%acc(ydim) = pset%p(p)%acc(ydim) + wijk(p, c) * fy(p, c)
               pset%p(p)%acc(zdim) = pset%p(p)%acc(zdim) + wijk(p, c) * fz(p, c)
            enddo
            !write(*,*) "------", p, pset%p(p)%acc(xdim), pset%p(p)%acc(ydim), pset%p(p)%acc(zdim)
         enddo
         !czastka nr statyczna:
         !pset%p(1)%acc = 0.0

      end subroutine get_acc_cic

      subroutine get_var_timestep_c(dt_nbody, max_acc, lf_c, pset, cg)

         use constants,      only: ndims, xdim, zdim, big, one
         use func,           only: operator(.notequals.)
         use grid_cont,      only: grid_container
         use particle_types, only: particle_set

         implicit none

         type(grid_container), pointer, intent(in)  :: cg
         class(particle_set),           intent(in)  :: pset  !< particle list
         real,                          intent(in)  :: max_acc, lf_c
         real,                          intent(out) :: dt_nbody
         real                                       :: eta, eps, factor
         real, dimension(ndims)                     :: maxv, minv, max_v
         integer                                    :: cdim

         eta = 1.0
         eps = 1.0e-4
         factor = big

         if(max_acc.notequals.0.0) then
            dt_nbody = sqrt(2.0*eta*eps/max_acc)

            do cdim = xdim, zdim
               maxv(cdim)  = abs(maxval(pset%p(:)%vel(cdim)))
               minv(cdim)  = abs(minval(pset%p(:)%vel(cdim)))
               max_v(cdim) = max(maxv(cdim), minv(cdim))
            enddo


            if (any(max_v*dt_nbody > cg%dl)) then

               if (any(max_v.notequals.0.0)) then
                  do cdim = xdim, zdim
                     if ((max_v(cdim).notequals.0.0)) then
                        factor = min(cg%dl(cdim)/max_v(cdim), factor)
                     endif
                  enddo
               endif
            else
               factor = one
            endif

            dt_nbody  = lf_c * factor * dt_nbody
         else
            dt_nbody = zero
         endif
         !write(*,*) "[get_var_timestep_c]: dt_nbody =", dt_nbody

      end subroutine get_var_timestep_c

      subroutine save_pot_pset(cg, pset)

         use constants,      only: xdim, ydim, zdim, CENTER
         use dataio_pub,     only: printinfo, printio
         use grid_cont,      only: grid_container
         use particle_types, only: particle_set

         implicit none

         class(particle_set),           intent(in) :: pset  !< particle list
         type(grid_container), pointer, intent(in) :: cg

         integer                                   :: numer1, numer2, i, j, p
         logical                                   :: finish, save_potential

         save_potential = .true.
         finish         = .false.
         numer1 = 64
         numer2 = 20

         open(unit=90, file='pset.dat', action='write', position='append')
            do p=1, ubound(pset%p, dim=1)
               write(90,*) p, pset%p(p)%pos
            enddo
         close(90)
         if (save_potential) then
            call printio('[particle_integrators:save_pot_pset] Writing potential to files: potencjal1.dat, potencjal2.dat')
            open(unit=88, file='potencjal1.dat')
            open(unit=89, file='potencjal2.dat')
               do i = lbound(cg%gpot,dim=1), ubound(cg%gpot,dim=1)
                  do j = lbound(cg%gpot,dim=2), ubound(cg%gpot,dim=2)
                     !do k=lbound(cg%gpot,dim=3),ubound(cg%gpot,dim=3)
                     write(88,*) i, j, numer1, cg%coord(CENTER, xdim)%r(i), cg%coord(CENTER, ydim)%r(j), cg%coord(CENTER, zdim)%r(numer1), cg%gpot(i,j,numer1)
                     write(89,*) i, j, numer2, cg%coord(CENTER, xdim)%r(i), cg%coord(CENTER, ydim)%r(j), cg%coord(CENTER, zdim)%r(numer2), cg%gpot(i,j,numer2)
                     !enddo
                  enddo
                  write(88,*)
                  write(89,*)
               enddo
            close(88)
            close(89)

            if (finish) then
               call printinfo('[particle_integrators:save_pot_pset] Condition to end - finishing.')
               stop
            endif
         endif

      end subroutine save_pot_pset

   end subroutine timestep_nbody
#endif /* NBODY */

end module particle_integrators
