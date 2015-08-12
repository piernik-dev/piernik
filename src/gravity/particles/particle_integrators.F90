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
   use particle_types, only: particle_solver_T, particle_set
   use constants,      only: cbuff_len

#ifdef NBODY
   use domain,         only: is_refined, is_multicg
#endif /* NBODY */

   implicit none
   !------------------

   !------------------
   private
   public :: hermit4, leapfrog2, acc_interp_method, lf_c, get_timestep_nbody, dt_nbody


   type, extends(particle_solver_T) :: hermit4_T
   contains
      procedure, nopass :: evolve => hermit_4ord
   end type hermit4_T
   
   type, extends(particle_solver_T) :: leapfrog2_T
   contains
      procedure, nopass :: evolve => leapfrog2ord
   end type leapfrog2_T

   type(hermit4_T),   target :: hermit4
   type(leapfrog2_T), target :: leapfrog2


#ifdef NBODY
   character(len=cbuff_len)  :: acc_interp_method  !< acceleration interpolation method
   real                       :: lf_c               !< timestep should depends of grid and velocities of particles (used to extrapolation of the gravitational potential)
   real                       :: dt_nbody           !< timestep depends on particles

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
      
      real, intent(in) :: t_glob, dt_tot

      real, parameter :: dt_param = 0.0001      ! control parameter to determine time step size
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
            !dt = dt_param * coll_time  !variable timestep
            dt = dt_param              !constant timestep
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
         print *, "Hermit dt=", dt
      enddo
      print *, "Hermit nsteps=", nsteps

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
      !use constants,      only: xdim, zdim!, half, zero
      use particle_types, only: particle_set
      use domain,         only: is_refined, is_multicg
      use cg_list,        only: cg_list_element
      use grid_cont,      only: grid_container
      use dataio_pub,     only: die, printinfo
      use global,         only: dt_old


      implicit none 
      class(particle_set), intent(inout) :: pset                 !< particle list

      real, intent(in)                   :: t_glob               !< initial time of simulation
      real, intent(in)                   :: dt_tot               !< timestep of simulation
      real, dimension(:), allocatable    :: mass                 !< 1D array of mass of the particles
      real                               :: dt_tot_h             !< half of timestep, dt_tot_h = 0.5*dt_tot
      real                               :: energy               !< total energy of set of particles
      real                               :: init_energy          !< total energy of set of particles at t_glob
      real                               :: d_energy = 0.0       !< error of energy of set of particles in succeeding timesteps, at t=t_glob=0.0
      real                               :: ang_momentum         !< angular momentum of set of particles
      real                               :: init_ang_mom         !< angular momentum of set of particles at t_glob
      real                               :: d_ang_momentum = 0.0 !< error of angular momentum in succeeding timensteps, at t=t_glob=0.0

      integer                            :: i
      integer                            :: n                    !< number of particles
      logical                            :: external_pot         !< if .true. gravitational potential will be deleted and replaced by external potential of point mass
      logical, save                      :: first_run_lf = .true.
      integer, save                      :: counter
      integer                            :: lun_out


      if (first_run_lf) then
         select case (acc_interp_method)
            case('cic', 'CIC')
               call printinfo("[particle_integrators:leapfrog2ord] Acceleration interpolation method: CIC") 
            case('lagrange', 'Lagrange')
               call printinfo("[particle_integrators:leapfrog2ord] Acceleration interpolation method: Lagrange polynomials")
         end select

      endif


      open(newunit=lun_out, file='leapfrog_out.log', status='unknown',  position='append')

      n = size(pset%p, dim=1)


      allocate(mass(n))

      if (t_glob < 0.0) i=1 ! supress compiler warnings
      

      mass(:) = pset%p(:)%mass


!      cgl => leaves%first
      if (is_refined) call die("[particle_integrators:leapfrog2ord] AMR not implemented for particles yet")
      if (is_multicg) call die("[particle_integrators:leapfrog2ord] multi_cg not implemented for particles yet")

      !do while (associated(cgl))
!            cg => cgl%cg
            !cgl => cgl%nxt
      !enddo


      !obliczenie zewnętrznego potencjalu na siatce
      !external_pot = .true.
      external_pot = .false.






      !call get_ang_momentum_2(pset, n, ang_momentum)
      !init_ang_mom = ang_momentum



      !call get_energy(pset, cg, cells, dist, n, energy)
      !init_energy = energy


      counter = 1


      do i = 1, n
         write(lun_out, '(I3,1X,16(E13.6,1X))') i, t_glob+dt_tot, dt_old, mass(i), pset%p(i)%pos, pset%p(i)%vel, pset%p(i)%acc, energy, d_energy, ang_momentum, d_ang_momentum
      enddo

      !call save_particles(n, lf_t, mass, pset, counter)
      write(*,*) "[p_i]:-----------------------------"
      write(*,*) "[p_i]:dt_tot= ", dt_tot
      write(*,*) "[p_i]:dt_old= ", dt_old
      dt_tot_h = 0.5*dt_tot

      if(first_run_lf) then
         first_run_lf = .false.
      else
         !3.kick(dt_old)
         call kick(pset, 0.5*dt_old, n)
         call energy(pset, n)
      endif
      !1. Kick (dt_tot_h)
      call kick(pset, dt_tot_h, n)

      !2.drift(lf_dt)
      call drift(pset, dt_tot, n)

      !call save_particles(n, lf_t, mass, pset, counter)

      close(lun_out)
      write(*,*) "[p_i]:-----------------------------"


      contains

         !Kick
         subroutine kick(pset, t, n)
            use particle_types, only: particle_set
            implicit none
            class(particle_set), intent(inout)     :: pset  !< particle list
            real, intent(in)                       :: t
            integer, intent(in)                    :: n
            integer                                :: i


            do i = 1, n
               pset%p(i)%vel = pset%p(i)%vel + pset%p(i)%acc * t
            enddo

         end subroutine kick

         !Drift
         subroutine drift(pset, t, n)
            use particle_types, only: particle_set
            implicit none
            class(particle_set), intent(inout) :: pset  !< particle list
            real, intent(in)                   :: t
            integer                            :: i
            integer, intent(in)                :: n

            do i=1, n
               pset%p(i)%pos = pset%p(i)%pos + pset%p(i)%vel * t
            enddo

         end subroutine drift


         subroutine energy(pset, n)
            use constants,    only: ndims, half, zero
            use particle_types, only: particle_set
            implicit none
            class(particle_set), intent(inout) :: pset  !< particle list
            integer, intent(in)                :: n      !< number of particles
            integer                              :: cdim
            real, dimension(n)                  :: v     !< kwadraty predkosci czastek 

            v = zero
            do p = 1, n
               do cdim = 1, ndims
                  v(p) = pset%p(p)%vel(cdim)**2
               enddo
               !energia      = 1/2  *      m         *  v**2 +     Ep(x,y,z)
               pset%p(p)%pot = half * pset%p(p)%mass *  v(p) + pset%p(p)%pot
            enddo
         end subroutine energy





         subroutine save_particles(n, lf_t, mass, pset, counter)
            use particle_types, only: particle_set
            implicit none
               class(particle_set), intent(in) :: pset  !< particle list
               integer, intent (in)            :: n
               integer, intent (inout)         :: counter
               real, dimension(n), intent(in)  :: mass
               integer                         :: i, data_file=757
               real, intent(in)                :: lf_t
               character(len=17)               :: filename
               character(len=3)                :: counter_char

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
                  do i=1,n
                     write(data_file,*) i, mass(i), pset%p(i)%pos, pset%p(i)%vel
                  enddo
               close(data_file)

               counter = counter + 1

         end subroutine save_particles

         subroutine potential2(pset, cg, cells, dist, n)
            use constants,    only: ndims, half, xdim, ydim, zdim
            use grid_cont,    only: grid_container
               implicit none
               type(grid_container), pointer, intent(in)     :: cg
               class(particle_set), intent(inout)            :: pset  !< particle list
               integer, intent(in)                           :: n
               integer, dimension(n, ndims), intent(in)      :: cells
               real(kind=8), dimension(n, ndims), intent(in) :: dist
               integer                                       :: p
               real,dimension(n)                             :: dpot, d2pot
               
               do p = 1, n
                  dpot(p) = df_dx_o2([cells(p, :)], cg) * dist(p, xdim) + &
                         df_dy_o2([cells(p, :)], cg) * dist(p, ydim) + &
                         df_dz_o2([cells(p, :)], cg) * dist(p, zdim)
                  
                  d2pot(p) = d2f_dx2_o2([cells(p, :)], cg) * dist(p, xdim)**2 + &
                          d2f_dy2_o2([cells(p, :)], cg) * dist(p, ydim)**2 + &
                          d2f_dz2_o2([cells(p, :)], cg) * dist(p, zdim)**2 + &
                          2.0*d2f_dxdy_o2([cells(p, :)], cg) * dist(p, xdim)*dist(p, ydim) + &
                          2.0*d2f_dxdz_o2([cells(p, :)], cg) * dist(p, xdim)*dist(p, zdim)
                  pset%p(p)%pot = cg%gpot(cells(p, xdim), cells(p, ydim), cells(p, zdim)) + dpot(p) + half * d2pot(p)
               enddo
         end subroutine potential2

!         subroutine potential(pset, cg, cells, dist, n)!poprawic te funcje, bo teraz nie dziala prawidlowo :/
!            use constants,    only: ndims, half, xdim, ydim, zdim
!            use grid_cont,    only: grid_container
!            implicit none
!            type(grid_container), pointer, intent(in)     :: cg
!            class(particle_set), intent(inout)            :: pset  !< particle list
!            integer, intent(in)                           :: n
!            integer, dimension(n, ndims), intent(in)      :: cells
!            real(kind=8), dimension(n, ndims), intent(in) :: dist
!            integer                                       :: i
!            integer (kind=8)                              :: p, q, r
!            real,dimension(n)                             :: dpot, d2pot


 !           do i = 1, n

!               dpot(i) = df_dx_o2([cells(i, :)], cg) * dist(i, xdim) + &
!                      df_dy_o2([cells(i, :)], cg) * dist(i, ydim) + &
!                      df_dz_o2([cells(i, :)], cg) * dist(i, zdim)
!               
!               d2pot(i) = d2f_dx2_o2([cells(i, :)], cg) * dist(i, xdim)**2 + &
!                       d2f_dy2_o2([cells(i, :)], cg) * dist(i, ydim)**2 + &
!                       d2f_dz2_o2([cells(i, :)], cg) * dist(i, zdim)**2 + &
!                       2.0*d2f_dxdy_o2([cells(i, :)], cg) * dist(i, xdim)*dist(i, ydim) + &
!                       2.0*d2f_dxdz_o2([cells(i, :)], cg) * dist(i, xdim)*dist(i, zdim)
!            enddo

!            do i = 1, n, 1
!               p = cells(i, xdim)
!               q = cells(i, ydim)
!               r = cells(i, zdim)
!               pset%p(i)%pot = cg%gpot(p, q, r) + dpot(i) + half * d2pot(i)
!            enddo

!         end subroutine potential


!         subroutine get_energy(pset, cg, cells, dist, n, energy)
!            use constants,    only: ndims
!            use grid_cont,    only: grid_container
!            implicit none
!               type(grid_container), pointer, intent(in)     :: cg
!               class(particle_set), intent(inout)            :: pset  !< particle list
!               integer                                       :: i, j
!               integer, intent(in)                           :: n
!               integer, dimension(n, ndims), intent(in)      :: cells
!               real(kind=8), dimension(n, ndims), intent(in) :: dist
!               real, intent(out)                             :: energy
!               real                                          :: velocity = 0.0

!               call potential(pset, cg, cells, dist, n)

!               energy = 0.0

!               do i=1, n
!                  do j=1, ndims
!                     velocity = velocity + pset%p(i)%vel(j)**2
!                  enddo
!
!                  energy = energy + 0.5*velocity + pset%p(i)%pot
!                  velocity = 0.0
!               enddo
!         end subroutine get_energy


!         subroutine get_acc_model(pset, acc2, eps, n)
!            use constants, only: ndims, xdim, ydim, zdim
!            use grid_cont,  only: grid_container
!            implicit none
!               class(particle_set), intent(in)        :: pset  !< particle list
!               integer, intent(in)                    :: n
!               real, intent(in)                       :: eps
!               real, dimension(n, ndims), intent(out) :: acc2

!               do i=1,n
!                  acc2(i, xdim) = -der_x(pset%p(i)%pos, 1.0e-8, eps)
!                  acc2(i, ydim) = -der_y(pset%p(i)%pos, 1.0e-8, eps)
!                  acc2(i, zdim) = -der_z(pset%p(i)%pos, 1.0e-8, eps)
!               enddo

!         end subroutine get_acc_model


!         function der_x(pos, d, eps)
!         implicit none
!            real(kind=8)                :: x, y, z, der_x, d, eps
!            real(kind=8),dimension(1,3) :: pos
!            x = pos(1,1)
!            y = pos(1,2)
!            z = pos(1,3)
!            der_x = ( phi_pm(x+d, y, z, eps) - phi_pm(x-d, y, z, eps) ) / (2.0*d)
!         end function der_x

         !Pochodna wzgledem y
!         function der_y(pos, d, eps)
!            implicit none
!               real(kind=8) :: x, y, z, der_y, d, eps
!               real(kind=8),dimension(1,3) :: pos

!               x = pos(1,1)
!               y = pos(1,2)
!               z = pos(1,3)
!               der_y = ( phi_pm(x, y+d, z, eps) - phi_pm(x, y-d, z, eps) ) / (2.0*d)
!         end function der_y

         !Pochodna wzgledem z
!         function der_z(pos, d, eps)
!            implicit none
!               real(kind=8)                :: x, y, z, der_z, d, eps
!               real(kind=8),dimension(1,3) :: pos

!               x = pos(1,1)
!               y = pos(1,2)
!               z = pos(1,3)
!               der_z = ( phi_pm(x, y, z+d, eps) - phi_pm(x, y, z-d, eps) ) / (2.0*d)
!         end function der_z

         
         subroutine get_ang_momentum_2(pset, n, ang_momentum)
            use constants, only : xdim, ydim, zdim
            use particle_types, only: particle_set
            implicit none
               class(particle_set), intent(in) :: pset
               integer                         :: i
               integer, intent(in)             :: n
               real, intent(out)               :: ang_momentum
               real                            :: L1, L2, L3

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
      integer, intent(in)                    :: n
      real, dimension(n), intent(in)         :: mass
      real, dimension(n, ndims), intent(in)  :: pos
      real, dimension(n, ndims), intent(out) :: acc
      real, intent(out)                      :: epot
      
      integer                                :: i, j
      real, dimension(ndims)                 :: rji, da
      
      real :: r   ! | rji |
      real :: r2  ! | rji |^2
      real :: r3  ! | rji |^3
      
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

#ifdef NBODY
   subroutine get_timestep_nbody(dt_nbody, pset)
      use constants,      only: ndims, zero
      use particle_types, only: particle_set
      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use grid_cont,      only: grid_container

      implicit none

      interface
         function dxi(cell, cg)
            use constants, only: ndims, xdim, ydim, zdim
            use grid_cont, only: grid_container 
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, dimension(ndims), intent(in)     :: cell
               real                                      :: dxi
         end function dxi

         function d2dxi2(cell, cg)
            use constants, only: ndims, xdim, ydim, zdim
            use grid_cont, only: grid_container 
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, dimension(ndims), intent(in)     :: cell
               real                                      :: d2dxi2
         end function d2dxi2

         function d2dxixj(cell, cg)
            use constants, only: ndims, xdim, ydim, zdim
            use grid_cont, only: grid_container 
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, dimension(ndims), intent(in)     :: cell
               real                                      :: d2dxixj
         end function d2dxixj

      end interface

      class(particle_set),  intent(inout)  :: pset
      type(grid_container),  pointer       :: cg
      type(cg_list_element), pointer       :: cgl

      integer                              :: order               !< order of Lagrange polynomials (if acc_interp_method = 'lagrange')
      real                                 :: eta, eps
      integer                              :: n_part
      real                                 :: max_acc

      real,    dimension(:,:), allocatable :: dist
      integer, dimension(:,:), allocatable :: cells
      real                                 :: dt_nbody
      logical                              :: save_potential, finish




      procedure(dxi), pointer     :: df_dx_p    => NULL(), &
                                     df_dy_p    => NULL(), &
                                     df_dz_p    => NULL()
      procedure(d2dxi2), pointer  :: d2f_dx2_p  => NULL(), &
                                     d2f_dy2_p  => NULL(), &
                                     d2f_dz2_p  => NULL() 
      procedure(d2dxixj), pointer :: d2f_dxdy_p => NULL(), &
                                     d2f_dxdz_p => NULL(), &
                                     d2f_dydz_p => NULL()

      write(*,*) "Poszukiwanie kroku nbody"
      eta = 1.0
      eps = 1.0e-4
      !write(*,*) "Przed n_part"
      n_part = size(pset%p, dim=1)
      write(*,*) "Number of particles: ", n_part
      allocate(cells(n_part, ndims), dist(n_part, ndims))

      !write(*,*) "Przed cg"
      cgl => leaves%first
      cg  => cgl%cg

      !save_potential = .true.
      save_potential = .false.
      !finish         = .true.
      finish         = .false.

      !call pot2grid(cg, zero)

      call save_pot(save_potential, finish, cg)

      if (acc_interp_method == 'lagrange') then
         order = 4
         call check_ord(order, df_dx_p, d2f_dx2_p, df_dy_p, d2f_dy2_p, & 
                  df_dz_p, d2f_dz2_p, d2f_dxdy_p, d2f_dxdz_p, d2f_dydz_p)
      endif

      call find_cells(cells, dist, cg, n_part)

      call potential2(pset, cg, cells, dist, n)    !szukanie energii potencjalnej w punktach-polozeniach czastek
      
      select case (acc_interp_method)
         case('lagrange', 'Lagrange', 'polynomials')
            call get_acc_int(cells, dist, pset, cg, n_part, &
                                 df_dx_p, d2f_dx2_p, df_dy_p, d2f_dy2_p,& 
                                 df_dz_p, d2f_dz2_p, d2f_dxdy_p, d2f_dxdz_p, d2f_dydz_p)
         case('cic', 'CIC')
            call get_acc_cic(pset, cg, cells, n_part)
      end select

      call get_acc_max(pset, n_part, max_acc)
      !write(*,*) "[get_timestep_nbody]: max_acc=", max_acc
      !write(*,*) "[get_timestep_nbody]:  eta   =", eta
      !write(*,*) "[get_timestep_nbody]:  eps   =", eps
      !write(*,*) "[get_timestep_nbody]:  lf_c  =", lf_c

      call get_var_timestep_c(dt_nbody, eta, eps, max_acc, lf_c, pset, cg)
      write(*,*) "Koniec szukania kroku nbody"


   contains

      function phi_pm(x, y, z, eps)
         use units,    only: newtong
         implicit none
            real, intent(in) :: x, y, z, eps
            real             :: r, phi_pm, G,M, mu
               G = 1.0
               M = 1.0
               mu = newtong*M
               !write(*,*) "[phi_pm: newtong=]", newtong
               !write(*,*) "[phi_pm: mu     =]", mu
               r = sqrt(x**2 + y**2 + z**2 + eps**2)
               !stop

               phi_pm = -mu / r
      end function phi_pm


      subroutine pot2grid(cg, eps2)
         use constants, only: xdim, ydim, zdim, CENTER
         use grid_cont, only: grid_container
         implicit none
            type(grid_container), pointer, intent(inout) :: cg
            integer                                      :: i, j, k
            real, intent(in)                             :: eps2

               do i = lbound(cg%gpot, dim=1), ubound(cg%gpot, dim=1)
                  do j = lbound(cg%gpot, dim=2), ubound(cg%gpot, dim=2)
                     do k = lbound(cg%gpot, dim=3), ubound(cg%gpot, dim=3)
                        cg%gpot(i,j,k) = phi_pm(cg%coord(CENTER,xdim)%r(i), &
                                                cg%coord(CENTER,ydim)%r(j), &
                                                cg%coord(CENTER,zdim)%r(k),eps2)
                     enddo
                  enddo
               enddo

      end subroutine pot2grid


      subroutine check_ord(order, df_dx_p, d2f_dx2_p, df_dy_p, d2f_dy2_p,& 
                  df_dz_p, d2f_dz2_p, d2f_dxdy_p, d2f_dxdz_p, d2f_dydz_p)
         implicit none
            integer,intent(in) :: order
            procedure(dxi), pointer, intent(inout)     :: df_dx_p, df_dy_p, df_dz_p
            procedure(d2dxi2), pointer, intent(inout)  :: d2f_dx2_p, d2f_dy2_p, d2f_dz2_p
            procedure(d2dxixj), pointer, intent(inout) :: d2f_dxdy_p, d2f_dxdz_p, d2f_dydz_p

               if (order == 2) then
                  df_dx_p => df_dx_o2
                  df_dy_p => df_dy_o2
                  df_dz_p => df_dz_o2
                  d2f_dx2_p => d2f_dx2_o2
                  d2f_dy2_p => d2f_dy2_o2
                  d2f_dz2_p => d2f_dz2_o2
                  d2f_dxdy_p => d2f_dxdy_o2
                  d2f_dxdz_p => d2f_dxdz_o2
                  d2f_dydz_p => d2f_dydz_o2
               else
                  df_dx_p => df_dx_o4
                  df_dy_p => df_dy_o4
                  df_dz_p => df_dz_o4
                  d2f_dx2_p => d2f_dx2_o4
                  d2f_dy2_p => d2f_dy2_o4
                  d2f_dz2_p => d2f_dz2_o4
                  d2f_dxdy_p => d2f_dxdy_o4
                  d2f_dxdz_p => d2f_dxdz_o4
                  d2f_dydz_p => d2f_dydz_o4
               endif
      end subroutine check_ord


      subroutine find_cells(cells, dist, cg, n_part)

         use constants,      only: ndims, xdim, CENTER, LO, HI
         use grid_cont,      only: grid_container
         !use particle_types,   only: particle_set
         implicit none
            !class(particle_set)                          :: pset  !< particle list
            type(grid_container)                          :: cg
            integer                                       :: i, cdim
            integer, intent(in)                           :: n_part
            integer,dimension(n_part, ndims), intent(out) :: cells
            real,dimension(n_part, ndims), intent(out)    :: dist

            !write(*,*) "Finding cells"
            !write(*,*) "[find_cells]: Particles =", n_part 
            do i = 1, n_part
               do cdim = xdim, ndims
                  if ((pset%p(i)%pos(cdim) >= cg%ijkse(cdim, LO)) .or. (pset%p(i)%pos(cdim) <= cg%ijkse(cdim, HI))) then
                     pset%p(i)%outside = .false.
                  else
                     pset%p(i)%outside = .true.
                  endif
               
                  cells(i, cdim) = int( 0.5 + (pset%p(i)%pos(cdim) - cg%coord(CENTER,cdim)%r(0)) / cg%dl(cdim) )

                  dist(i, cdim)  = pset%p(i)%pos(cdim) - ( cg%coord(CENTER, cdim)%r(0) + cells(i,cdim) * cg%dl(cdim) )
               enddo

            enddo

      end subroutine find_cells


      function df_dx_o2(cell, cg)
         use constants, only: ndims, xdim, ydim, zdim
         use grid_cont, only: grid_container
         implicit none
            type(grid_container), pointer, intent(in) :: cg
            integer, dimension(ndims), intent(in)     :: cell
            real,target                               :: df_dx_o2

               !o(R^2)
               df_dx_o2 = ( cg%gpot(cell(xdim)+1, cell(ydim), cell(zdim)) - &
                           cg%gpot(cell(xdim)-1, cell(ydim), cell(zdim)) ) / (2.0*cg%dx)

      end function df_dx_o2


      function df_dy_o2(cell, cg)
         use constants, only: ndims, xdim, ydim, zdim
         use grid_cont, only: grid_container 
         implicit none
            type(grid_container), pointer, intent(in) :: cg
            integer, dimension(ndims), intent(in)     :: cell
            real,target                               :: df_dy_o2

               !o(R^2)
               df_dy_o2 = ( cg%gpot(cell(xdim), cell(ydim)+1, cell(zdim)) - &
                           cg%gpot(cell(xdim), cell(ydim)-1, cell(zdim)) ) / (2.0*cg%dy)

      end function df_dy_o2


      function df_dz_o2(cell, cg)
         use constants, only: ndims, xdim, ydim, zdim
         use grid_cont, only: grid_container 
         implicit none
            type(grid_container), pointer, intent(in) :: cg
            integer, dimension(ndims), intent(in)     :: cell
            real,target                               :: df_dz_o2

               !o(R^2)
               df_dz_o2 = ( cg%gpot(cell(xdim), cell(ydim), cell(zdim)+1) - &
                           cg%gpot(cell(xdim), cell(ydim), cell(zdim)-1) ) / (2.0*cg%dz)

      end function df_dz_o2


      function d2f_dx2_o2(cell, cg)
         use constants, only: ndims, xdim, ydim, zdim
         use grid_cont, only: grid_container 
         implicit none
            type(grid_container), pointer, intent(in) :: cg
            integer, dimension(ndims), intent(in)     :: cell
            real,target                               :: d2f_dx2_o2

               !o(R^2)
               d2f_dx2_o2 = ( cg%gpot(cell(xdim)+1, cell(ydim), cell(zdim)) - &
                           2.0*cg%gpot(cell(xdim), cell(ydim), cell(zdim)) + &
                           cg%gpot(cell(xdim)-1, cell(ydim), cell(zdim)) ) / (cg%dx**2)

      end function d2f_dx2_o2


      function d2f_dy2_o2(cell, cg)
         use constants, only: ndims, xdim, ydim, zdim
         use grid_cont, only: grid_container 
         implicit none
            type(grid_container), pointer, intent(in) :: cg
            integer, dimension(ndims), intent(in)     :: cell
            real,target                               :: d2f_dy2_o2

               !o(R^2)
               d2f_dy2_o2 = ( cg%gpot(cell(xdim), cell(ydim)+1, cell(zdim)) - &
                           2.0*cg%gpot(cell(xdim), cell(ydim), cell(zdim)) + &
                           cg%gpot(cell(xdim), cell(ydim)-1, cell(zdim)) ) / (cg%dy**2)

      end function d2f_dy2_o2


      function d2f_dz2_o2(cell, cg)
         use constants, only: ndims, xdim, ydim, zdim
         use grid_cont, only: grid_container 
         implicit none
            type(grid_container), pointer, intent(in) :: cg
            integer, dimension(ndims), intent(in)     :: cell
            real,target                               :: d2f_dz2_o2

               !o(R^2)
               d2f_dz2_o2 = ( cg%gpot(cell(xdim), cell(ydim), cell(zdim)+1) - &
                           2.0*cg%gpot(cell(xdim), cell(ydim), cell(zdim)) + &
                           cg%gpot(cell(xdim), cell(ydim), cell(zdim)-1) ) / (cg%dz**2)

      end function d2f_dz2_o2


      function d2f_dxdy_o2(cell, cg)
         use constants, only: ndims, xdim, ydim, zdim
         use grid_cont, only: grid_container 
         implicit none
            type(grid_container), pointer, intent(in) :: cg
            integer, dimension(ndims), intent(in)     :: cell
            real,target                               :: d2f_dxdy_o2

               !o(R^2)
               d2f_dxdy_o2 = ( cg%gpot(cell(xdim)+1, cell(ydim)+1, cell(zdim)) - &
                           cg%gpot(cell(xdim)+1, cell(ydim)-1, cell(zdim)) - &
                           cg%gpot(cell(xdim)-1, cell(ydim)+1, cell(zdim)) + &
                           cg%gpot(cell(xdim)-1, cell(ydim)-1, cell(zdim)) ) / (4.0*cg%dx*cg%dy)

      end function d2f_dxdy_o2


      function d2f_dxdz_o2(cell, cg)
         use constants, only: ndims, xdim, ydim, zdim
         use grid_cont, only: grid_container 
         implicit none
            type(grid_container), pointer, intent(in) :: cg
            integer, dimension(ndims), intent(in)     :: cell
            real,target                               :: d2f_dxdz_o2

               !o(R^2)
               d2f_dxdz_o2 = ( cg%gpot(cell(xdim)+1, cell(ydim), cell(zdim)+1) - &
                           cg%gpot(cell(xdim)+1, cell(ydim), cell(zdim)-1) - &
                           cg%gpot(cell(xdim)-1, cell(ydim), cell(zdim)+1) + &
                           cg%gpot(cell(xdim)-1, cell(ydim), cell(zdim)-1) ) / (4.0*cg%dx*cg%dz)

      end function d2f_dxdz_o2


      function d2f_dydz_o2(cell, cg)
         use constants, only: ndims, xdim, ydim, zdim
         use grid_cont, only: grid_container 
         implicit none
            type(grid_container), pointer, intent(in) :: cg
            integer, dimension(ndims), intent(in)     :: cell
            real,target                               :: d2f_dydz_o2

               !o(R^2)
               d2f_dydz_o2 = ( cg%gpot(cell(xdim), cell(ydim)+1, cell(zdim)+1) - &
                           cg%gpot(cell(xdim), cell(ydim)+1, cell(zdim)-1) - &
                           cg%gpot(cell(xdim), cell(ydim)-1, cell(zdim)+1) + &
                           cg%gpot(cell(xdim), cell(ydim)-1, cell(zdim)-1) ) / (4.0*cg%dy*cg%dz)

      end function d2f_dydz_o2


      function df_dx_o4(cell, cg)
         use constants, only: ndims, xdim, ydim, zdim
         use grid_cont, only: grid_container 
         implicit none
            type(grid_container), pointer, intent(in) :: cg
            integer, dimension(ndims), intent(in)     :: cell
            real,target                               :: df_dx_o4

               !o(R^4)
               df_dx_o4 = 2.0 * (cg%gpot(cell(xdim)+1, cell(ydim), cell(zdim)) - &
                        cg%gpot(cell(xdim)-1, cell(ydim), cell(zdim)) ) / (3.0*cg%dx) - &
                        (cg%gpot(cell(xdim)+2, cell(ydim), cell(zdim)) - &
                        cg%gpot(cell(xdim)-2, cell(ydim), cell(zdim)) ) / (12.0*cg%dx)

      end function df_dx_o4


      function df_dy_o4(cell, cg)
         use constants, only: ndims, xdim, ydim, zdim
         use grid_cont, only: grid_container 
         implicit none
            type(grid_container), pointer, intent(in) :: cg
            integer, dimension(ndims), intent(in)     :: cell
            real,target                               :: df_dy_o4

               !o(R^4)
               df_dy_o4 = 2.0 * ( cg%gpot(cell(xdim), cell(ydim)+1, cell(zdim)) - &
                        cg%gpot(cell(xdim), cell(ydim)-1, cell(zdim)) ) / (3.0*cg%dy) - &
                        (cg%gpot(cell(xdim), cell(ydim)+2, cell(zdim)) - &
                        cg%gpot(cell(xdim), cell(ydim)-2, cell(zdim)) ) / (12.0*cg%dy)

      end function df_dy_o4


      function df_dz_o4(cell, cg)
         use constants, only: ndims, xdim, ydim, zdim
         use grid_cont, only: grid_container 
         implicit none
            type(grid_container), pointer, intent(in) :: cg
            integer, dimension(ndims), intent(in)     :: cell
            real,target                               :: df_dz_o4

               !o(R^4)
               df_dz_o4 = 2.0* (cg%gpot(cell(xdim), cell(ydim), cell(zdim)+1) - &
                        cg%gpot(cell(xdim), cell(ydim), cell(zdim)-1) ) / (3.0*cg%dz) - &
                        ( cg%gpot(cell(xdim), cell(ydim), cell(zdim)+2) - &
                        cg%gpot(cell(xdim), cell(ydim), cell(zdim)-2) ) / (12.0*cg%dz)

      end function df_dz_o4


      function d2f_dx2_o4(cell, cg)
         use constants, only: ndims, xdim, ydim, zdim
         use grid_cont, only: grid_container 
         implicit none
            type(grid_container), pointer, intent(in) :: cg
            integer, dimension(ndims), intent(in)     :: cell
            real,target                               :: d2f_dx2_o4

               !o(R^4)
               d2f_dx2_o4 = 4.0 * ( cg%gpot(cell(xdim)+1, cell(ydim), cell(zdim)) + &
                           cg%gpot(cell(xdim)-1, cell(ydim), cell(zdim)) - &
                           2.0 * cg%gpot(cell(xdim), cell(ydim), cell(zdim)) ) / (3.0*cg%dx**2) - &
                           (cg%gpot(cell(xdim)+2, cell(ydim), cell(zdim)) + &
                           cg%gpot(cell(xdim)-2, cell(ydim), cell(zdim)) - &
                           2.0 * cg%gpot(cell(xdim), cell(ydim), cell(zdim)) ) / (12.0*cg%dx**2)

      end function d2f_dx2_o4


      function d2f_dy2_o4(cell, cg)
         use constants, only: ndims, xdim, ydim, zdim
         use grid_cont, only: grid_container 
         implicit none
            type(grid_container), pointer, intent(in) :: cg
            integer, dimension(ndims), intent(in)     :: cell
            real,target                               :: d2f_dy2_o4

               !o(R^4)
               d2f_dy2_o4 = 4.0*( cg%gpot(cell(xdim), cell(ydim)+1, cell(zdim)) + &
                           cg%gpot(cell(xdim), cell(ydim)-1, cell(zdim)) - &
                           2.0*cg%gpot(cell(xdim), cell(ydim), cell(zdim)) ) / (3.0*cg%dy**2) - &
                           (cg%gpot(cell(xdim), cell(ydim)+2, cell(zdim)) + &
                           cg%gpot(cell(xdim), cell(ydim)-2, cell(zdim)) - &
                           2.0*cg%gpot(cell(xdim), cell(ydim), cell(zdim)) ) / (12.0*cg%dy**2)

      end function d2f_dy2_o4


      function d2f_dz2_o4(cell, cg)
         use constants, only: ndims, xdim, ydim, zdim
         use grid_cont, only: grid_container 
         implicit none
            type(grid_container), pointer, intent(in) :: cg
            integer, dimension(ndims), intent(in)     :: cell
            real,target                               :: d2f_dz2_o4

               !o(R^4)
               d2f_dz2_o4 = 4.0*( cg%gpot(cell(xdim), cell(ydim), cell(zdim)+1) + &
                           cg%gpot(cell(xdim), cell(ydim), cell(zdim)-1) - &
                           2.0*cg%gpot(cell(xdim), cell(ydim), cell(zdim)) ) / (3.0*cg%dz**2) - &
                           ( cg%gpot(cell(xdim), cell(ydim), cell(zdim)+2) + &
                           cg%gpot(cell(xdim), cell(ydim), cell(zdim)-2) - &
                           2.0*cg%gpot(cell(xdim), cell(ydim), cell(zdim)) ) / (12.0*cg%dz**2)

      end function d2f_dz2_o4


      function d2f_dxdy_o4(cell, cg)
         use constants, only: ndims, xdim, ydim, zdim
         use grid_cont, only: grid_container 
         implicit none
            type(grid_container), pointer, intent(in) :: cg
            integer, dimension(ndims), intent(in)     :: cell
            real,target                               :: d2f_dxdy_o4

               !o(R^4)
               d2f_dxdy_o4 = ( cg%gpot(cell(xdim)+1, cell(ydim)+1, cell(zdim)) + &
                           cg%gpot(cell(xdim)-1, cell(ydim)-1, cell(zdim)) - &
                           cg%gpot(cell(xdim)+1, cell(ydim)-1, cell(zdim)) - &
                           cg%gpot(cell(xdim)-1, cell(ydim)+1, cell(zdim)) ) / (3.0*cg%dx*cg%dy) - &
                           ( cg%gpot(cell(xdim)+2, cell(ydim)+2, cell(zdim)) + &
                           cg%gpot(cell(xdim)-2, cell(ydim)-2, cell(zdim)) - &
                           cg%gpot(cell(xdim)+2, cell(ydim)-2, cell(zdim)) - &
                           cg%gpot(cell(xdim)-2, cell(ydim)+2, cell(zdim)) ) / (48.0*cg%dx*cg%dy)
    
      end function d2f_dxdy_o4


      function d2f_dxdz_o4(cell, cg)
         use constants, only: ndims, xdim, ydim, zdim
         use grid_cont, only: grid_container 
         implicit none
            type(grid_container), pointer, intent(in) :: cg
            integer, dimension(ndims), intent(in)     :: cell
            real,target                               :: d2f_dxdz_o4

               !o(R^4)
               d2f_dxdz_o4 = ( cg%gpot(cell(xdim)+1, cell(ydim), cell(zdim)+1) + &
                           cg%gpot(cell(xdim)-1, cell(ydim), cell(zdim)-1) - &
                           cg%gpot(cell(xdim)+1, cell(ydim), cell(zdim)-1) - &
                           cg%gpot(cell(xdim)-1, cell(ydim), cell(zdim)+1) ) / (3.0*cg%dx*cg%dz) - &
                           ( cg%gpot(cell(xdim)+2, cell(ydim), cell(zdim)+2) + &
                           cg%gpot(cell(xdim)-2, cell(ydim), cell(zdim)-2) - &
                           cg%gpot(cell(xdim)+2, cell(ydim), cell(zdim)-2) - &
                           cg%gpot(cell(xdim)-2, cell(ydim), cell(zdim)+2) ) / (48.0*cg%dx*cg%dz)

      end function d2f_dxdz_o4


      function d2f_dydz_o4(cell, cg)
         use constants, only: ndims, xdim, ydim, zdim
         use grid_cont, only: grid_container 
         implicit none
            type(grid_container), pointer, intent(in) :: cg
            integer, dimension(ndims), intent(in)     :: cell
            real,target                               :: d2f_dydz_o4

               !o(R^4)
               d2f_dydz_o4 = ( cg%gpot(cell(xdim), cell(ydim)+1, cell(zdim)+1) + &
                           cg%gpot(cell(xdim), cell(ydim)-1, cell(zdim)-1) - &
                           cg%gpot(cell(xdim), cell(ydim)+1, cell(zdim)-1) - &
                           cg%gpot(cell(xdim), cell(ydim)-1, cell(zdim)+1) ) / (3.0*cg%dy*cg%dz) - &
                           ( cg%gpot(cell(xdim), cell(ydim)+2, cell(zdim)+2) + &
                           cg%gpot(cell(xdim), cell(ydim)-2, cell(zdim)-2) - &
                           cg%gpot(cell(xdim), cell(ydim)+2, cell(zdim)-2) - &
                           cg%gpot(cell(xdim), cell(ydim)-2, cell(zdim)+2) ) / (48.0*cg%dy*cg%dz)

      end function d2f_dydz_o4


      subroutine get_acc_int(cells, dist, pset, cg, n_part, &
                              df_dx_p, d2f_dx2_p, df_dy_p, d2f_dy2_p,& 
                              df_dz_p, d2f_dz2_p, d2f_dxdy_p, d2f_dxdz_p, d2f_dydz_p)
         use constants,      only: ndims, xdim, ydim, zdim
         use grid_cont,      only: grid_container
         use particle_types, only: particle_set
         implicit none
            class(particle_set), intent(inout)            :: pset
            type(grid_container), pointer, intent(in)     :: cg
            integer, intent(in)                           :: n_part
            integer, dimension(n_part, ndims), intent(in) :: cells
            real, dimension(n_part, ndims), intent(in)    :: dist
            integer                                       :: i

               
            procedure(dxi), pointer, intent(in)     :: df_dx_p, df_dy_p, df_dz_p
            procedure(d2dxi2), pointer, intent(in)  :: d2f_dx2_p, d2f_dy2_p, d2f_dz2_p
            procedure(d2dxixj), pointer, intent(in) :: d2f_dxdy_p, d2f_dxdz_p, d2f_dydz_p


               do i=1, n_part
                  if ((pset%p(i)%outside) .eqv. .false.) then
                     pset%p(i)%acc(xdim) = - (df_dx_p([cells(i, :)], cg) + &
                                    d2f_dx2_p([cells(i, :)], cg)  * dist(i, xdim) + &
                                    d2f_dxdy_p([cells(i, :)], cg) * dist(i, ydim) + &
                                    d2f_dxdz_p([cells(i, :)], cg) * dist(i, zdim))

                     pset%p(i)%acc(ydim) = -( df_dy_p([cells(i, :)], cg) + &
                                    d2f_dy2_p([cells(i, :)], cg)  * dist(i, ydim) + &
                                    d2f_dxdy_p([cells(i, :)], cg) * dist(i, xdim) + &
                                    d2f_dydz_p([cells(i, :)], cg) * dist(i, zdim))

                     pset%p(i)%acc(zdim) = -( df_dz_p([cells(i, :)], cg) + &
                                    d2f_dz2_p([cells(i, :)], cg)  * dist(i, zdim) + &
                                    d2f_dxdz_p([cells(i, :)], cg) * dist(i, xdim) + &
                                    d2f_dydz_p([cells(i, :)], cg) * dist(i, ydim))
                  !else
                  !   call !funkcja liczaca pochodne z potencjalu policzonego z rozwiniecia multipolowego
                  endif
               enddo
               !pset%p(1)%acc = 0.0

               !stara wersja wykorzystujaca tablice
               !acc(:, xdim) = -( df_dx_p(cells, cg, n) + &
               !            !   d2f_dx2_p(cells, cg, n)  * dist(:, xdim) + &
               !            !   d2f_dxdy_p(cells, cg, n) * dist(:, ydim) + &
               !            !   d2f_dxdz_p(cells, cg, n) * dist(:, zdim))

               !acc(:, ydim) = -( df_dy_p(cells, cg, n) + &
               !            !   d2f_dy2_p(cells, cg, n)  * dist(:, ydim) + &
               !            !   d2f_dxdy_p(cells, cg, n) * dist(:, xdim) + &
               !            !   d2f_dydz_p(cells, cg, n) * dist(:, zdim))

               !acc(:, zdim) = -( df_dz_p(cells, cg, n) + &
               !            !   d2f_dz2_p(cells, cg, n)  * dist(:, zdim) + &
               !            !   d2f_dxdz_p(cells, cg, n) * dist(:, xdim) + &
               !            !   d2f_dydz_p(cells, cg, n) * dist(:, ydim))

      end subroutine get_acc_int

      subroutine get_acc_max(pset, n_part, max_acc)
         use constants,      only: xdim, ndims
         use particle_types, only: particle_set
         implicit none
            class(particle_set), intent(inout) :: pset
            integer, intent(in)                :: n_part
            integer                            :: i, cdim
            real, dimension(n_part)            :: acc
            real, intent(out)                  :: max_acc

               acc = 0.0

               do i = 1, n_part
                  do cdim = xdim, ndims
                        acc(i) = acc(i) + pset%p(i)%acc(cdim)**2
                  enddo
               enddo

               max_acc = sqrt(maxval(acc))

      end subroutine get_acc_max


      subroutine get_acc_cic(pset, cg, cells, n_part)
         use constants,      only: ndims, CENTER, xdim, ydim, zdim, half, zero
         use grid_cont,      only: grid_container
         use particle_types, only: particle_set
         implicit none
            type(grid_container), pointer, intent(in)     :: cg
            class(particle_set), intent(inout)            :: pset

            integer, intent(in)                           :: n_part
            integer                                       :: i, j, k, c, cdim
            integer                                       :: p
            integer, dimension(n_part, ndims), intent(in) :: cells
            integer(kind=8), dimension(n_part, ndims)     :: cic_cells
            real, dimension(n_part, ndims)                :: dxyz
            real(kind=8), dimension(n_part, 8)            :: wijk, fx, fy, fz


               write(*,*) "[get_acc_cic]: particles = ", n_part
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
                     !write(*,*) cells(i,:)
                     !write(*,*) cic_cells(i,:)
                     
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
               !write(*,*) "wijk"
               !write(*,*) wijk(1,:)

!stare---------------------------------------------------------------------------
               do p = 1, n_part
               !write(*,*) "1:Czastka ", p
                  c = 1
                  do i = 0, 1
                     do j = 0, 1
                        do k = 0, 1
                        !write (*,*) i, j, k, c
                        !write(*,'(A6, 3(I4), A3, 3(I4), A1)') "x:   [", i+1, j, k, "] [", i-1, j, k, "]"
                        !write(*,'(A6, 3(I4), A3, 3(I4), A1)') "y:   [", i, j+1, k, "] [", i, j-1, k, "]" 
                        !write(*,'(A6, 3(I4), A3, 3(I4), A1)') "z:   [", i, j, k+1, "] [", i, j, k-1, "]" 
                        
                           fx(p, c) = -(cg%gpot(cic_cells(p, xdim)+1+i, cic_cells(p, ydim)  +j, cic_cells(p, zdim)  +k) - cg%gpot(cic_cells(p, xdim)-1+i, cic_cells(p, ydim)  +j, cic_cells(p, zdim)  +k))
                           fy(p, c) = -(cg%gpot(cic_cells(p, xdim)  +i, cic_cells(p, ydim)+1+j, cic_cells(p, zdim)  +k) - cg%gpot(cic_cells(p, xdim)  +i, cic_cells(p, ydim)-1+j, cic_cells(p, zdim)  +k))
                           fz(p, c) = -(cg%gpot(cic_cells(p, xdim)  +i, cic_cells(p, ydim)  +j, cic_cells(p, zdim)+1+k) - cg%gpot(cic_cells(p, xdim)  +i, cic_cells(p, ydim)  +j, cic_cells(p, zdim)-1+k))
                           c = c + 1
                        enddo
                     enddo
                  enddo
               enddo
!---------------------------------------------------------------------------------

               fx = half*fx*cg%idx
               fy = half*fy*cg%idy
               fz = half*fz*cg%idz
               !write(*,*) "fx ", fx
               !write(*,*) "fy ", fy
               !write(*,*) "fz ", fz

               !write(*,*) "acc(xdim)", pset%p(1)%acc(xdim)
               !write(*,*) "acc(ydim)", pset%p(1)%acc(ydim)
               !write(*,*) "acc(zdim)", pset%p(1)%acc(zdim)
               
               do p = 1, n_part
                  !write(*,*) "2:Czastka ", p
                  do c = 1, 8
                  !write(*,*) pset%p(p)%acc(xdim), wijk(p, c), fx(p, c), pset%p(p)%acc(xdim) + wijk(p, c) * fx(p, c)
                     pset%p(p)%acc(xdim) = pset%p(p)%acc(xdim) + wijk(p, c) * fx(p, c)
                     
                     pset%p(p)%acc(ydim) = pset%p(p)%acc(ydim) + wijk(p, c) * fy(p, c)
                     pset%p(p)%acc(zdim) = pset%p(p)%acc(zdim) + wijk(p, c) * fz(p, c)
                  enddo
                  write(*,*) "------", p, pset%p(p)%acc(xdim), pset%p(p)%acc(ydim), pset%p(p)%acc(zdim)
                  !stop
                  !pset%p(p)%acc(xdim) = wijk(p,1) * fx(p,1) + wijk(p, 8)*fx(p,8) + &
                  !                  wijk(p,2) * fx(p,2) + wijk(p, 7)*fx(p,7) + &
                  !                  wijk(p,3) * fx(p,3) + wijk(p, 6)*fx(p,6) + &
                  !                  wijk(p,4) * fx(p,4) + wijk(p, 5)*fx(p,5)
               enddo
               !write(*,*) "acc(xdim)", pset%p(1)%acc(xdim)
               !write(*,*) "acc(ydim)", pset%p(1)%acc(ydim)
               !write(*,*) "acc(zdim)", pset%p(1)%acc(zdim)
               !stop

      end subroutine get_acc_cic



      subroutine get_var_timestep_c(dt_nbody, eta, eps, max_acc, lf_c, pset, cg)
         use constants,      only: ndims, xdim, zdim, big, one
         use particle_types, only: particle_set
         use grid_cont,      only: grid_container
         use func,           only: operator(.notequals.)
         implicit none
            type(grid_container), pointer, intent(in) :: cg
            class(particle_set), intent(in)           :: pset  !< particle list

            real, intent(in)       :: eta, eps, max_acc, lf_c
            real, intent(out)      :: dt_nbody
            real                   :: factor
            real, dimension(ndims) :: maxv, minv, max_v
            integer                :: cdim

               factor = big
               
               if(max_acc.notequals.0.0) then
               dt_nbody = sqrt(2.0*eta*eps/max_acc)
               !write(*,*) "[get_var_timestep_c]: dt_nbody =", dt_nbody

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

               !write(*,*) "[get_var_timestep_c]:  factor  =", factor
               dt_nbody  = lf_c * factor * dt_nbody
               else
               dt_nbody = zero
               endif
               write(*,*) "[get_var_timestep_c]: dt_nbody =", dt_nbody

      end subroutine get_var_timestep_c


      subroutine save_pot(save_potential, finish, cg)
         use constants, only: xdim, ydim, zdim, CENTER
         use grid_cont, only: grid_container
            implicit none
            type(grid_container), pointer, intent(in) :: cg
            logical :: save_potential, finish
            integer :: i, j, k

            if(save_potential) then
                  write(*,*) "Zapis potencjalu do pliku"
                  open(unit=88, file='potencjal.dat')
                     do i=lbound(cg%gpot,dim=1),ubound(cg%gpot,dim=1)
                        do j=lbound(cg%gpot,dim=2),ubound(cg%gpot,dim=2)
                           do k=lbound(cg%gpot,dim=3),ubound(cg%gpot,dim=3)
                              write(88,*) i, j, k, cg%coord(CENTER, xdim)%r(i), &
                                       cg%coord(CENTER, ydim)%r(i), cg%coord(CENTER, zdim)%r(i), &
                                       cg%gpot(i,j,k)
                           enddo
                        enddo
                     write(88,*)
                  enddo
               close(88)

               if(finish) then
                  write(*,*) "Warunek zakonczenia-zatrzymano"
                  stop
               endif
            endif
      end subroutine save_pot

   end subroutine get_timestep_nbody

#endif /*NBODY */

end module particle_integrators
