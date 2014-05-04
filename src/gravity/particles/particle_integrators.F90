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
   !------------------

   !------------------
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

      real, parameter :: dt_param = 0.0001        ! control parameter to determine time step size
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
      use constants, only: ndims, xdim, zdim
      use particle_types, only: particle_set
      use domain, only: dom
      implicit none 
      class(particle_set), intent(inout) :: pset  !< particle list
           
      
      interface
         function df_dxi(cells, potential, n_cell, delta_xi, n_particles)
            use constants, only: ndims
            integer, intent(in) :: n_particles
            integer,dimension(n_particles, ndims),intent(in) :: cells
            integer, dimension(ndims), intent(in) :: n_cell
            real,dimension(n_cell(1), n_cell(2), n_cell(3)),intent(in) :: potential
            real,intent(in) :: delta_xi
            real,dimension(n_particles) :: df_dxi
         end function df_dxi

         function d2f_dxi_2(cells, potential, n_cell, delta_xi, n_particles)
            use constants, only: ndims
            integer, intent(in) :: n_particles
            integer,dimension(n_particles, ndims),intent(in) :: cells
            integer, dimension(ndims), intent(in) :: n_cell
            real,dimension(n_cell(1), n_cell(2), n_cell(3)),intent(in) :: potential
            real,intent(in) :: delta_xi
            real,dimension(n_particles) :: d2f_dxi_2
         end function d2f_dxi_2

         function d2f_dxi_dxj(cells, potential, n_cell, delta_xi, delta_xj, n_particles)
            use constants, only: ndims
            integer, intent(in) :: n_particles
            integer,dimension(n_particles, ndims),intent(in) :: cells
            integer, dimension(ndims), intent(in) :: n_cell
            real,dimension(n_cell(1), n_cell(2), n_cell(3)),intent(in) :: potential
            real,intent(in) :: delta_xi, delta_xj
            real,dimension(n_particles) :: d2f_dxi_dxj
         end function d2f_dxi_dxj
      
      end interface
   
      
      real, intent(in) :: t_glob, dt_tot
      integer, dimension(:,:),allocatable :: cells
      real, dimension(3) :: l_borders, r_borders
      real, dimension(:), allocatable :: mass, mins, maxs,delta_cells
      real, dimension(:, :), allocatable :: pos, vel, acc, vel_h, d_particles
      real, dimension(:, :, :), allocatable :: pot
      integer, dimension(:), allocatable :: n_cell
      real :: t_dia, t_out, t_end, einit, dt, t, dth, eta, eps, a, epot
      integer :: nsteps, n, ndim, lun_out, lun_err, i, j, k, nx, ny, nz, order
      real :: eps2, xmin, xmax, ymin, ymax, zmin, zmax, n_orbit, tend, dx, dy, dz, ax, ay, az, axx, ayy, azz, energy, init_energy, d_energy, ang_momentum, init_ang_mom, d_ang_momentum, zero

      procedure(df_dxi),pointer :: df_dx_p, df_dy_p, df_dz_p
      procedure(d2f_dxi_2),pointer :: d2f_dx2_p, d2f_dy2_p, d2f_dz2_p
      procedure(d2f_dxi_dxj),pointer :: d2f_dxdy_p, d2f_dxdz_p, d2f_dydz_p
      
      
      open(newunit=lun_out, file='leapfrog_out.log', status='unknown',  position='append')
      
      n = size(pset%p, dim=1)
      
      allocate(mass(n), pos(n, ndims), vel(n, ndims), acc(n, ndims), vel_h(n, ndims), cells(n, ndims), mins(ndims), maxs(ndims), delta_cells(ndims), n_cell(ndims), d_particles(n,ndims))
      
      
      mass(:) = pset%p(:)%mass

      do ndim = xdim, zdim
         pos(:, ndim) = pset%p(:)%pos(ndim)
         vel(:, ndim) = pset%p(:)%vel(ndim)
      enddo

      
      t = t_glob
      t_end = t + dt_tot
      print *, "Leafrog: t_end= ", t_end


      !mins(:) = dom%edge(:,1)
      mins=-3.0
      !maxs(:) = dom%edge(:,2)
      maxs=3.0
      !n_cell(:) = dom%n_d
      n_cell(1)=300
      n_cell(2)=300
      n_cell(3)=300


      zero = 0.0
      order = 4


      eta = 1.0
      eps = 1.0e-6
      eps2 = 0.00



      !alokacja potencjalu
      allocate( pot(n_cell(1), n_cell(2), n_cell(3)) )
      write(*,*) "Zaalokowano potencjal"

      !obliczenie potencjalu na siatce
      call pot_grid(pot, mins, maxs, n_cell, delta_cells, eps2)

      init_ang_mom = get_ang_momentum(pos, vel, mass, n)

      call cell_nr(pos, cells, mins, delta_cells, n)

      call check_ord(order, df_dx_p, d2f_dx2_p, df_dy_p, d2f_dy2_p,& 
                        df_dz_p, d2f_dz2_p, d2f_dxdy_p, d2f_dxdz_p, d2f_dydz_p)


      !initial acceleration
      call get_acc(cells, pos, acc, pot, n_cell, mins, delta_cells, n)


      !timestep
      !dt = sqrt(2.0*eta*eps/a)            !variable
      dt = 0.001                            !constant
      dth = dt/2.0
      

      nsteps = 0


      !main loop
      do while (t<t_end)
         !1.kick(dth)
         vel_h(:,:) = vel
         call kick(vel_h, acc, dth, n) !velocity
         !2.drift(dt)
         call drift(pos, vel_h, dt, n) !position
         !3.acceleration + |a|
         call cell_nr(pos, cells, mins, delta_cells, n)
         call get_acc(cells, pos, acc, pot, n_cell, mins, delta_cells, n)
         !call get_acc_mod(acc, n, a)
         !4.kick(dth)
         vel(:,:) = vel_h
         call kick(vel, acc, dth, n)   !velocity
         !5.t
         t = t + dt
         !6.dt		!dt[n+1]
         !dt	= sqrt(2.0*eta*eps/a)
         !7.dth
         !dth =	0.5*dt
         nsteps = nsteps + 1
         d_ang_momentum = log(abs((get_ang_momentum(pos, vel, mass, n) - init_ang_mom)/init_ang_mom))

         do i = 1, n
            write(lun_out, '(I3,1X,8(E13.6,1X))') i, mass(i), pos(i,:), vel(i,:), d_ang_momentum
         enddo

      end do
      
      print *, "Leapfrog: nsteps=", nsteps
      
      do ndim = xdim, zdim
         pset%p(:)%pos(ndim) = pos(:, ndim)
         pset%p(:)%vel(ndim) = vel(:, ndim)
      enddo

      deallocate (mass, pos, vel, acc, vel_h, cells, mins, maxs, delta_cells, n_cell, d_particles)
      close(lun_out)
      
      !return
      contains
         
         !Kick
         subroutine kick(vel, acc, t, n)
            use constants, only: ndims
            implicit none
            real, intent(in) :: t
            integer, intent(in) :: n
            real, dimension(n, ndims), intent(in) :: acc
            real, dimension(n, ndims), intent(inout) :: vel
		    
            vel = vel + acc*t
            
         end subroutine kick

         !Drift
         subroutine drift(pos, vel, t, n)
            use constants, only: ndims
            implicit none
            real, intent(in) :: t
            integer, intent(in) :: n
            real, dimension(n, ndims), intent(inout) :: pos
            real, dimension(n, ndims), intent(in) :: vel
           
            pos = pos + vel*t
         end subroutine drift



         function phi_pm(x, y, z, eps)
            implicit none
               real(kind=8) :: x, y, z, r, phi_pm, eps, G,M
                  G = 1.0
                  M = 1.0
                  r = sqrt(x**2 + y**2 + z**2)
                  phi_pm = -G*M / (sqrt(r**2 + eps**2))
         end function phi_pm


         subroutine pot_grid(pot, mins, maxs, n_cell, delta_cells, eps2)
            use constants, only: ndims
            implicit none
               integer :: i, j, k
               integer, dimension(ndims), intent(in) :: n_cell
               real,dimension(ndims),intent(in) :: mins, maxs
               real,dimension(ndims),intent(out) :: delta_cells
               real, intent(in) :: eps2
               real,dimension(n_cell(1), n_cell(2), n_cell(3)),intent(out) :: pot
               !open(unit=77,file='potencjal.dat')
               delta_cells = (maxs - mins) / n_cell

                  do i = 1, n_cell(1), 1
                     do j = 1, n_cell(2), 1
                        do k = 1, n_cell(3), 1
                           pot(i, j, k) = phi_pm(mins(1) + i*delta_cells(1), mins(2) + j*delta_cells(2), mins(3) + k*delta_cells(3), eps2)
                           !write(77,*) i,j,k,pot(i,j,k)
                        enddo
                     enddo
                  enddo
               !close(77)
         end subroutine pot_grid

         subroutine check_ord(order, df_dx_p, d2f_dx2_p, df_dy_p, d2f_dy2_p,& 
                  df_dz_p, d2f_dz2_p, d2f_dxdy_p, d2f_dxdz_p, d2f_dydz_p)
            implicit none
               integer,intent(in) :: order
               procedure(df_dxi),pointer :: df_dx_p, df_dy_p, df_dz_p
               procedure(d2f_dxi_2),pointer :: d2f_dx2_p, d2f_dy2_p, d2f_dz2_p
               procedure(d2f_dxi_dxj),pointer :: d2f_dxdy_p, d2f_dxdz_p, d2f_dydz_p
                  if (order==2) then
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

         subroutine get_acc(cells, pos, acc, pot, n_cell, mins, delta_cells, n)
            use constants, only : ndims
            implicit none
               integer:: i
               real::  dx_cell, dy_cell, dz_cell
               integer, intent(in) :: n
               integer,dimension(n, ndims), intent(in) :: cells
               real,dimension(n, ndims),intent(in) :: pos
               real,dimension(n, ndims),intent(out) :: acc
               integer,dimension(ndims),intent(in) :: n_cell
               real, dimension(ndims),intent(in) :: mins, delta_cells
               real,dimension(n_cell(1), n_cell(2), n_cell(3)),intent(in) :: pot
               real,dimension(n, ndims) :: delta
               real,dimension(n) :: ax, ay, az

                  do i = 1, ndims, 1
                     delta(:,i) = - ( cells(:,i) * delta_cells(i) + mins(i) - pos(:,i) )
                  enddo

                  dx_cell = delta_cells(1)
                  dy_cell = delta_cells(2)
                  dz_cell = delta_cells(3)

                  ax = -( df_dx_p(cells, pot, n_cell, dx_cell, n) + &
                     d2f_dx2_p(cells, pot, n_cell, dx_cell, n) * delta(:,1) + &
                     d2f_dxdy_p(cells, pot, n_cell, dx_cell, dy_cell, n) * delta(:,2) + &
                     d2f_dxdz_p(cells, pot, n_cell, dx_cell, dz_cell, n) * delta(:,3))
                  ay = -( df_dy_p(cells, pot, n_cell, dy_cell, n) + &
                     d2f_dy2_p(cells, pot, n_cell, dy_cell, n) * delta(:,2) + &
                     d2f_dxdy_p(cells, pot, n_cell, dx_cell, dy_cell, n) * delta(:,1) + &
                     d2f_dydz_p(cells, pot, n_cell, dy_cell, dz_cell, n) * delta(:,3))
                  az = -( df_dz_p(cells, pot, n_cell, dz_cell, n) + &
                     d2f_dz2_p(cells, pot, n_cell, dz_cell, n) * delta(:,3) + &
                     d2f_dxdz_p(cells, pot, n_cell, dx_cell, dz_cell, n) * delta(:,1) +&
                     d2f_dydz_p(cells, pot, n_cell, dy_cell, dz_cell, n) * delta(:,2))
                  acc(:,1) = ax
                  acc(:,2) = ay
                  acc(:,3) = az
         end subroutine get_acc

         subroutine cell_nr(pos, cells, mins, delta_cells, n)
            use constants, only: ndims
            implicit none
               integer :: i
               integer, intent(in) :: n
               integer,dimension(n, ndims),intent(out) :: cells
               real,dimension(n, ndims),intent(in) :: pos
               real, dimension(ndims), intent(in)  :: mins, delta_cells

               do i=1, ndims
                  cells(:,i) = int( (pos(:,i) - mins(i) - 0.5*delta_cells(i)) / delta_cells(i) ) + 1
               enddo
         end subroutine cell_nr

         function df_dx_o2(cells, pot, n_cell, dx_cell, n)
            use constants, only: ndims
            implicit none
               integer, intent(in) :: n
               integer,dimension(n, ndims),intent(in) :: cells
               integer, dimension(ndims), intent(in):: n_cell
               integer :: i, p, q, r
               real,dimension(n_cell(1), n_cell(2), n_cell(3)),intent(in) :: pot
               real,intent(in) :: dx_cell
               real,dimension(n),target :: df_dx_o2

               do i=1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^2)
                  df_dx_o2(i) = ( pot(p+1, q, r) - pot(p-1, q, r) ) / (2.0*dx_cell)
               enddo
         end function df_dx_o2

         function df_dx_o4(cells, pot, n_cell, dx_cell, n)
            use constants, only: ndims
            implicit none
               integer, intent(in) :: n
               integer,dimension(n, ndims),intent(in) :: cells
               integer, dimension(ndims), intent(in):: n_cell
               integer :: i, p, q, r
               real,dimension(n_cell(1), n_cell(2), n_cell(3)),intent(in) :: pot
               real,intent(in) :: dx_cell
               real,dimension(n),target :: df_dx_o4

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^4)
                  df_dx_o4(i) = ( 2.0* (pot(p+1, q, r) - pot(p-1, q, r) ) ) / (3.0*dx_cell) - &
                              ( pot(p+2, q, r) - pot(p-2, q, r) ) / (12.0*dx_cell)
               enddo
         end function df_dx_o4

         function df_dy_o2(cells, pot, n_cell, dy_cell, n)
            use constants, only: ndims
            implicit none
               integer, intent(in) :: n
               integer,dimension(n, ndims),intent(in) :: cells
               integer, dimension(ndims), intent(in):: n_cell
               integer :: i, p, q, r
               real,dimension(n_cell(1), n_cell(2), n_cell(3)),intent(in) :: pot
               real,intent(in) :: dy_cell
               real,dimension(n),target ::df_dy_o2

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^2)
                  df_dy_o2(i) = (pot(p, q+1, r) - pot(p, q-1, r) ) / (2.0*dy_cell)
               enddo
         end function df_dy_o2


         function df_dy_o4(cells, pot, n_cell, dy_cell, n)
            use constants, only: ndims
            implicit none
               integer, intent(in) :: n
               integer,dimension(n, ndims),intent(in) :: cells
               integer, dimension(ndims), intent(in):: n_cell
               integer :: i, p, q, r
               real,dimension(n_cell(1), n_cell(2), n_cell(3)),intent(in) :: pot
               real,intent(in) :: dy_cell
               real, dimension(n), target :: df_dy_o4

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^4)
                  df_dy_o4(i) = ( 2.0 * ( pot(p, q+1, r) - pot(p, q-1, r) ) ) / (3.0*dy_cell) - &
                        ( pot(p, q+2, r) - pot(p, q-2, r) ) / (12.0*dy_cell)
               enddo
         end function df_dy_o4

         function df_dz_o2(cells, pot, n_cell, dz_cell, n)
            use constants, only: ndims
            implicit none
               integer, intent(in) :: n
               integer,dimension(n, ndims),intent(in) :: cells
               integer, dimension(ndims), intent(in):: n_cell
               integer :: i, p, q, r
               real,dimension(n_cell(1), n_cell(2), n_cell(3)),intent(in) :: pot
               real,intent(in) :: dz_cell
               real, dimension(n), target :: df_dz_o2

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^2)
                  df_dz_o2(i) = ( pot(p, q, r+1) - pot(p, q, r-1) ) / (2.0*dz_cell)
               enddo
         end function df_dz_o2


         function df_dz_o4(cells, pot, n_cell, dz_cell, n)
            use constants, only: ndims
            implicit none
               integer, intent(in) :: n
               integer,dimension(n, ndims),intent(in) :: cells
               integer, dimension(ndims), intent(in):: n_cell
               integer :: i, p, q, r
               real,dimension(n_cell(1), n_cell(2), n_cell(3)),intent(in) :: pot
               real,intent(in) :: dz_cell
               real, dimension(n), target :: df_dz_o4

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^4)
                  df_dz_o4(i) = ( 2.0* (pot(p, q, r+1) - pot(p, q, r-1) ) ) / (3.0*dz_cell) - &
                           ( pot(p, q, r+2) - pot(p, q, r-2) ) / (12.0*dz_cell)
               enddo
         end function df_dz_o4 

         function d2f_dx2_o2(cells, pot, n_cell, dx_cell, n)
            use constants, only: ndims
            implicit none
               integer, intent(in) :: n
               integer,dimension(n, ndims),intent(in) :: cells
               integer, dimension(ndims), intent(in):: n_cell
               integer :: i, p, q, r
               real,dimension(n_cell(1), n_cell(2), n_cell(3)),intent(in) :: pot
               real,intent(in) :: dx_cell
               real,dimension(n),target :: d2f_dx2_o2

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^2)
                  d2f_dx2_o2(i) = (pot(p+1, q, r) - 2.0*pot(p, q, r) + pot(p-1, q, r) ) / (dx_cell**2)
               enddo
         end function d2f_dx2_o2

         function d2f_dx2_o4(cells, pot, n_cell, dx_cell, n)
            use constants, only: ndims
            implicit none
               integer, intent(in) :: n
               integer,dimension(n, ndims),intent(in) :: cells
               integer, dimension(ndims), intent(in) :: n_cell
               integer :: i, p, q, r
               real,dimension(n_cell(1), n_cell(2), n_cell(3)),intent(in) :: pot
               real,intent(in) :: dx_cell
               real,dimension(n),target :: d2f_dx2_o4

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^4)
                  d2f_dx2_o4(i) = 4.0 * ( pot(p+1, q, r) + pot(p-1, q, r) - &
                           2.0 * pot(p, q, r) ) / (3.0*dx_cell**2) - &
                           ( pot(p+2, q, r) + pot(p-2, q, r) - 2.0 * pot(p, q, r) ) / (12.0*dx_cell**2)
               enddo
         end function d2f_dx2_o4

         function d2f_dy2_o2(cells, pot, n_cell, dy_cell, n)
            use constants, only: ndims
            implicit none
               integer, intent(in) :: n
               integer,dimension(n, ndims),intent(in) :: cells
               integer, dimension(ndims), intent(in):: n_cell
               integer :: i, p, q, r
               real,dimension(n_cell(1), n_cell(2), n_cell(3)),intent(in) :: pot
               real,intent(in) :: dy_cell
               real,dimension(n),target :: d2f_dy2_o2

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^2)
                  d2f_dy2_o2(i) = ( pot(p, q+1, r) - 2.0*pot(p, q, r) + pot(p, q-1, r) ) / (dy_cell**2)
               enddo
         end function d2f_dy2_o2

         function d2f_dy2_o4(cells, pot, n_cell, dy_cell, n)
            use constants, only: ndims
            implicit none
               integer, intent(in) :: n
               integer,dimension(n, ndims),intent(in) :: cells
               integer, dimension(ndims), intent(in) :: n_cell
               integer :: i, p, q, r
               real,dimension(n_cell(1), n_cell(2), n_cell(3)),intent(in) :: pot
               real,intent(in) :: dy_cell
               real,dimension(n),target :: d2f_dy2_o4

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^4)
                  d2f_dy2_o4(i) = 4.0*( pot(p, q+1, r) + pot(p, q-1, r) - &
                        2.0*pot(p, q, r) ) / (3.0*dy_cell**2) - &
                        ( pot(p, q+2, r) + pot(p, q-2, r) - 2.0*pot(p, q, r) ) / (12.0*dy_cell**2)
               enddo
         end function d2f_dy2_o4

         function d2f_dz2_o2(cells, pot, n_cell, dz_cell, n)
            use constants, only: ndims
            implicit none
               integer, intent(in) :: n
               integer,dimension(n, ndims),intent(in) :: cells
               integer, dimension(ndims), intent(in) :: n_cell
               integer :: i, p, q, r
               real,dimension(n_cell(1), n_cell(2), n_cell(3)),intent(in) :: pot
               real,intent(in) :: dz_cell
               real,dimension(n),target :: d2f_dz2_o2

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^2)
                  d2f_dz2_o2(i) = ( pot(p, q, r+1) - 2.0*pot(p, q, r) + pot(p, q, r-1) ) / (dz_cell**2)
               enddo
         end function d2f_dz2_o2

         function d2f_dz2_o4(cells, pot, n_cell, dz_cell, n)
            use constants, only: ndims
            implicit none
               integer, intent(in) :: n
               integer,dimension(n, ndims),intent(in) :: cells
               integer, dimension(ndims), intent(in) :: n_cell
               integer :: i, p, q, r
               real,dimension(n_cell(1), n_cell(2), n_cell(3)),intent(in) :: pot
               real,intent(in) :: dz_cell
               real,dimension(n),target :: d2f_dz2_o4

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^4)
                  d2f_dz2_o4(i) = 4.0*( pot(p, q, r+1) + pot(p, q, r-1) - &
                           2.0*pot(p, q, r) ) / (3.0*dz_cell**2) - &
                           ( pot(p, q, r+2) + pot(p, q, r-2) - 2.0*pot(p, q, r) ) / (12.0*dz_cell**2)
               enddo
         end function d2f_dz2_o4

         function d2f_dxdy_o2(cells, pot, n_cell, dx_cell, dy_cell, n)
            use constants, only: ndims
            implicit none
               integer, intent(in) :: n
               integer,dimension(n, ndims),intent(in) :: cells
               integer, dimension(ndims), intent(in) :: n_cell
               integer :: i, p, q, r
               real,dimension(n_cell(1), n_cell(2), n_cell(3)),intent(in) :: pot
               real,intent(in) :: dx_cell, dy_cell
               real,dimension(n),target :: d2f_dxdy_o2

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^2)
                  d2f_dxdy_o2(i) = ( pot(p+1, q+1, r) - pot(p+1, q-1, r) - &
                              pot(p-1, q+1, r) + pot(p-1, q-1, r) ) / (4.0*dx_cell*dy_cell)
               enddo
         end function d2f_dxdy_o2


         function d2f_dxdy_o4(cells, pot, n_cell, dx_cell, dy_cell, n)
            use constants, only: ndims
            implicit none
               integer, intent(in) :: n
               integer,dimension(n, ndims),intent(in) :: cells
               integer, dimension(ndims), intent(in):: n_cell
               integer :: i, p, q, r
               real,dimension(n_cell(1), n_cell(2), n_cell(3)),intent(in) :: pot
               real,intent(in) :: dx_cell, dy_cell
               real,dimension(n),target :: d2f_dxdy_o4

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^4)
                  d2f_dxdy_o4(i) = ( pot(p+1, q+1, r) + pot(p-1, q-1, r) - pot(p+1, q-1, r) - &
                              pot(p-1, q+1, r) ) / (3.0*dx_cell*dy_cell) - &
                              ( pot(p+2, q+2, r) + pot(p-2, q-2, r) - pot(p+2, q-2, r) - &
                              pot(p-2, q+2, r) ) / (48.0*dx_cell*dy_cell)
               enddo
         end function d2f_dxdy_o4

         function d2f_dxdz_o2(cells, pot, n_cell, dx_cell, dz_cell, n)
            use constants, only: ndims
            implicit none
               integer, intent(in) :: n
               integer,dimension(n, ndims),intent(in) :: cells
               integer, dimension(ndims), intent(in) :: n_cell
               integer :: i, p, q, r
               real,dimension(n_cell(1), n_cell(2), n_cell(3)),intent(in) :: pot
               real,intent(in) :: dx_cell, dz_cell
               real,dimension(n),target :: d2f_dxdz_o2

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^2)
                  d2f_dxdz_o2(i) = ( pot(p+1, q, r+1) - pot(p+1, q, r-1) - &
                                 pot(p-1, q, r+1) + pot(p-1, q, r-1) ) / (4.0*dx_cell*dz_cell)
               enddo
         end function d2f_dxdz_o2

         function d2f_dxdz_o4(cells, pot, n_cell, dx_cell, dz_cell, n)
            use constants, only: ndims
            implicit none
               integer, intent(in) :: n
               integer,dimension(n, ndims),intent(in) :: cells
               integer, dimension(ndims), intent(in) :: n_cell
               integer :: i, p, q, r
               real,dimension(n_cell(1), n_cell(2), n_cell(3)),intent(in) :: pot
               real,intent(in) :: dx_cell, dz_cell
               real,dimension(n),target :: d2f_dxdz_o4

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^4)
                  d2f_dxdz_o4(i) = ( pot(p+1, q, r+1) + pot(p-1, q, r-1) - pot(p+1, q, r-1) - &
                              pot(p-1, q, r+1) ) / (3.0*dx_cell*dz_cell) - &
                              ( pot(p+2, q, r+2) + pot(p-2, q, r-2) - pot(p+2, q, r-2) - &
                              pot(p-2, q, r+2) ) / (48.0*dx_cell*dz_cell)
               enddo
         end function d2f_dxdz_o4
      
      
         function d2f_dydz_o2(cells, pot, n_cell, dy_cell, dz_cell, n)
            use constants, only: ndims
            implicit none
               integer, intent(in) :: n
               integer,dimension(n, ndims),intent(in) :: cells
               integer, dimension(ndims), intent(in) :: n_cell
               integer :: i, p, q, r
               real,dimension(n_cell(1), n_cell(2), n_cell(3)),intent(in) :: pot
               real,intent(in) :: dy_cell, dz_cell
               real,dimension(n),target :: d2f_dydz_o2

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^2)
                  d2f_dydz_o2(i) = ( pot(p, q+1, r+1) - pot(p, q+1, r-1) - &
                                 pot(p, q-1, r+1) + pot(p, q-1, r-1) ) / (4.0*dy_cell*dz_cell)
               enddo
         end function d2f_dydz_o2

         function d2f_dydz_o4(cells, pot, n_cell, dy_cell, dz_cell, n)
            use constants, only: ndims
            implicit none
               integer, intent(in) :: n
               integer,dimension(n, ndims),intent(in) :: cells
               integer, dimension(ndims), intent(in) :: n_cell
               integer :: i, p, q, r
               real,dimension(n_cell(1), n_cell(2), n_cell(3)),intent(in) :: pot
               real,intent(in) :: dy_cell, dz_cell
               real,dimension(n),target :: d2f_dydz_o4

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^4)
                  d2f_dydz_o4(i) = ( pot(p, q+1, r+1) + pot(p, q-1, r-1) - pot(p, q+1, r-1) - &
                              pot(p, q-1, r+1) ) / (3.0*dy_cell*dz_cell) - &
                              ( pot(p, q+2, r+2) + pot(p, q-2, r-2) - pot(p, q+2, r-2) - &
                              pot(p, q-2, r+2) ) / (48.0*dy_cell*dz_cell)
               enddo
         end function d2f_dydz_o4

         function get_ang_momentum(pos, vel, mass, n)
            use constants, only : ndims
            implicit none
               integer :: i, j
               integer, intent(in) :: n
               real, dimension(n, ndims), intent(in) :: pos, vel
               real, dimension(n) :: mass
               real :: ang_mom = 0.0, get_ang_momentum, r2, p2, r, p, rp

               do i = 1, n, 1
                  r2 = 0.0
                  p2 = 0.0
                  rp = 0.0
                  do j = 1, ndims
                     r2 = r2 + pos(i, j)**2
                     p2 = p2 + mass(i)**2 * vel(i, j)**2
                     rp = rp + pos(i, j) * mass(i) * vel(i, j)
                  enddo

                  ang_mom = ang_mom + sqrt( r2 * p2 * ( 1.0 - ( rp / sqrt(r2*p2) )**2 ) )
               enddo
               get_ang_momentum = ang_mom
         end function get_ang_momentum

  end subroutine leapfrog2ord      
      
   subroutine get_acc_pot(mass, pos, acc, n, epot)
      use constants, only: ndims
      implicit none
      integer, intent(in) :: n
      real, dimension(n), intent(in) :: mass
      real, dimension(n, ndims), intent(in) :: pos
      real, dimension(n, ndims), intent(out) :: acc
      real  :: eps
      real, intent(out) :: epot
      
      integer :: i, j
      real, dimension(ndims) :: rji, vji, da
      
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
