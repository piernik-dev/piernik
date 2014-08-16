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
      use constants, only: ndims
      use particle_types, only: particle_set
      use domain, only: dom
      use cg_leaves,    only: leaves
      use cg_list,      only: cg_list_element
      use grid_cont,    only: grid_container


      implicit none 
      class(particle_set), intent(inout) :: pset  !< particle list
      type(grid_container), pointer :: cg
      type(cg_list_element),  pointer  :: cgl

      interface
         function df_dxi(neighb, cg, n_particles)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
            type(grid_container), pointer, intent(in) :: cg
            integer, intent(in) :: n_particles
            integer(kind=8),dimension(n_particles, ndims),intent(in) :: neighb
            real,dimension(n_particles) :: df_dxi
         end function df_dxi

         function d2f_dxi_2(neighb, cg, n_particles)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
            type(grid_container), pointer, intent(in) :: cg
            integer, intent(in) :: n_particles
            integer(kind=8),dimension(n_particles, ndims),intent(in) :: neighb
            real,dimension(n_particles) :: d2f_dxi_2
         end function d2f_dxi_2

         function d2f_dxi_dxj(neighb, cg, n_particles)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
            type(grid_container), pointer, intent(in) :: cg
            integer, intent(in) :: n_particles
            integer(kind=8),dimension(n_particles, ndims),intent(in) :: neighb
            real,dimension(n_particles) :: d2f_dxi_dxj
         end function d2f_dxi_dxj
      
      end interface
      

      integer(kind=8), dimension(:,:), allocatable :: neighb
      real(kind=8), dimension(:,:), allocatable :: dist, acc2, acc3
      
      real, intent(in) :: t_glob, dt_tot
      integer, dimension(:,:),allocatable :: cells
      real, dimension(:), allocatable :: mass, mins, maxs,delta_cells
      real, dimension(:, :), allocatable :: pos, vel, acc, vel_h, d_particles
      integer, dimension(:), allocatable :: n_cell
      real :: t_end, dt, t, dth, eta, eps, a, eps2, energy, init_energy, d_energy = 0.0, ang_momentum, init_ang_mom, d_ang_momentum = 0.0, zero
      integer :: nsteps, n, lun_out, i, order


      procedure(df_dxi),pointer :: df_dx_p, df_dy_p, df_dz_p
      procedure(d2f_dxi_2),pointer :: d2f_dx2_p, d2f_dy2_p, d2f_dz2_p
      procedure(d2f_dxi_dxj),pointer :: d2f_dxdy_p, d2f_dxdz_p, d2f_dydz_p
      
      
      open(newunit=lun_out, file='leapfrog_out.log', status='unknown',  position='append')
      
      n = size(pset%p, dim=1)
      !
      allocate(neighb(n, ndims), dist(n, ndims), acc2(n, ndims))

      

      allocate(mass(n), pos(n, ndims), vel(n, ndims), acc(n, ndims), acc3(n, ndims), vel_h(n, ndims), cells(n, ndims), mins(ndims), maxs(ndims), delta_cells(ndims), n_cell(ndims), d_particles(n,ndims))

      

      mass(:) = pset%p(:)%mass

      !do ndim = xdim, zdim
      !   pos(:, ndim) = pset%p(:)%pos(ndim)
      !   vel(:, ndim) = pset%p(:)%vel(ndim)
      !enddo

      
      t = t_glob
      t_end = t + dt_tot
      !print *, "Leafrog: t_end= ", t_end


      mins(:) = dom%edge(:,1)

      maxs(:) = dom%edge(:,2)

      n_cell(:) = dom%n_d


      zero = 0.0
      order = 4


      eta = 1.0 !1.0
      eps = 1.0e-4
      eps2 = 0.00


      cgl => leaves%first
      do while (associated(cgl))
            cg => cgl%cg
           
            cgl => cgl%nxt
      enddo
      



      !obliczenie potencjalu na siatce
      call pot_grid2(cg, mins, delta_cells, eps)
      write(*,*) "Obliczono potencjal"


      
     ! open(unit=88, file='potencjal.dat')
     ! do i=lbound(cg%gpot,dim=1),ubound(cg%gpot,dim=1)
     !    do j=lbound(cg%gpot,dim=2),ubound(cg%gpot,dim=2)
     !       do k=lbound(cg%gpot,dim=3),ubound(cg%gpot,dim=3)
     !          write(88,*) i,j,k,cg%gpot(i,j,k)
     !       enddo
     !    enddo
     !    write(88,*)
     ! enddo
     ! close(88)
      call get_ang_momentum_2(pset, n, ang_momentum)
      init_ang_mom = ang_momentum

      call cell_nr2(pset, neighb, dist, mins, cg, n)




      call check_ord(order, df_dx_p, d2f_dx2_p, df_dy_p, d2f_dy2_p,& 
                        df_dz_p, d2f_dz2_p, d2f_dxdy_p, d2f_dxdz_p, d2f_dydz_p)


      !initial acceleration
      call get_acc2(neighb, dist, acc, cg, n)
      write(*,*) "Obliczono przyspieszenie"
      
      
      call get_energy(pset, cg, neighb, dist, n, energy)
      init_energy = energy
      
      
      call get_acc_num2(pset, acc2, eps, n)
      !write(*,*) "Znaleziono potencjal modelowy"
      
      call grav_pot2acc_cic2(pset, cg, neighb, acc3, n)
     


      call get_acc_mod(acc, n, a)
      
      
      !timestep
      dt = sqrt(2.0*eta*eps/a) 
      !write(*,*) "dt"           !variable
      !dt = 0.01                            !constant
      dth = dt/2.0
      

      nsteps = 0


      
      !main loop
      do while (t < t_end)


         if (t + dt > t_end) then
            dt = t_end - t
            dth = 0.5 * dt
         endif


         do i=1, n
            write(lun_out, '(I3,1X,13(E13.6,1X))') n, t, dt, mass(i), pset%p(i)%pos, acc(i,:), acc2(i,:), acc3(i,:), energy, d_energy, ang_momentum, d_ang_momentum
         enddo

         !------------------------------------------------------------------
         !acc=0.0
         !call kick(vel, acc, dth, n) !velocity vel_h
         !call drift(pos, vel, dt, n) !position vel_h
         !call get_acc_num(pos, acc2, eps, n)
         !call cell_nr(pos, cells, mins, delta_cells, n)
         !call get_acc(cells, pos, acc, pot, n_cell, mins, delta_cells, n)
         !call get_acc_mod(acc, n, a)
         !call kick(vel, [zero,zero,zero], dth, n)   !velocity
         !------------------------------------------------------------------

         !1.kick(dth)
         !acc(:,:) = 0.0                                                 !odkomentowac dla ruchu po prostej
         call kick2(pset, acc, dth, n)

         !2.drift(dt)         
         call drift2(pset, dt, n)


         call cell_nr2(pset, neighb, dist, mins, cg, n)                 !finding cells

         !3.acceleration + |a|
         call get_acc2(neighb, dist, acc, cg, n)                        !Lagrange polynomials acceleration
         call get_acc_num2(pset, acc2, eps, n)                          !centered finite differencing acceleration (if gravitational potential is known explicite)
         call grav_pot2acc_cic2(pset, cg, neighb, acc3, n)              !CIC acceleration
         call get_acc_mod(acc, n, a)                                    !max(|a_i|)
         
         !4.kick(dth)
         !call kick2(pset, acc, dth, n)                                  !zakomentowac dla ruchu po prostej
         call kick2(pset, acc, dth, n)
         !call kick2(pset, [zero,zero,zero], dth, n)                     !odkomentowac dla ruchu po prostej

         
         call get_energy(pset, cg, neighb, dist, n, energy)
         d_energy = log(abs((energy - init_energy)/init_energy))

         call get_ang_momentum_2(pset, n, ang_momentum)
         d_ang_momentum = log(abs((ang_momentum - init_ang_mom)/init_ang_mom))

         
         !5.t
         t = t + dt
         !6.dt		!dt[n+1]
         dt = sqrt(2.0*eta*eps/a)
         !7.dth
         dth = 0.5*dt
         nsteps = nsteps + 1

      end do
      
      write(*,*) "Leapfrog: nsteps=", nsteps
      


      deallocate (acc, acc2, acc3, neighb, dist, mins, maxs, delta_cells, n_cell)
      close(lun_out)
      

      contains
         
         !Kick
         subroutine kick2(pset, acc, t, n)
            use constants, only: ndims
            use particle_types, only: particle_set
            implicit none
            class(particle_set), intent(inout) :: pset  !< particle list
            real, intent(in) :: t
            integer, intent(in) :: n
            integer :: i
            real, dimension(n, ndims), intent(in) :: acc

            do i=1,n
               pset%p(i)%vel = pset%p(i)%vel + acc(i,:)*t
            enddo
         end subroutine kick2

         !Drift
         subroutine drift2(pset, t, n)
            use particle_types, only: particle_set
            implicit none
            class(particle_set), intent(inout) :: pset  !< particle list
            real, intent(in) :: t
            integer :: i
            integer, intent(in) :: n

            do i=1,n
               pset%p(i)%pos = pset%p(i)%pos + pset%p(i)%vel * t
            enddo
         end subroutine drift2


         function phi_pm(x, y, z, eps)
            implicit none
               real(kind=8) :: x, y, z, r, phi_pm, eps, G,M
                  G = 1.0
                  M = 1.0
                  r = sqrt(x**2 + y**2 + z**2 + eps**2)

                  phi_pm = -G*M / r
         end function phi_pm


         subroutine pot_grid2(cg, mins, delta_cells, eps)
            use constants, only: ndims
            use grid_cont, only: grid_container
            implicit none
               type(grid_container), pointer, intent(inout) :: cg
               integer :: i, j, k
               real,dimension(ndims),intent(in) :: mins
               real,dimension(ndims),intent(out) :: delta_cells
               real, intent(in) :: eps

               !open(unit=77,file='potencjal.dat')


                  delta_cells(1) = cg%dx  !potem trzeba to usunac
                  delta_cells(2) = cg%dy
                  delta_cells(3) = cg%dz


                  do i = lbound(cg%gpot, dim=1), ubound(cg%gpot, dim=1)
                     do j = lbound(cg%gpot, dim=2), ubound(cg%gpot, dim=2)
                        do k = lbound(cg%gpot, dim=3), ubound(cg%gpot, dim=3)
                           cg%gpot(i, j, k) = phi_pm(mins(1) + (i-0.5)*cg%dx, &
                                                     mins(2) + (j-0.5)*cg%dy, &
                                                     mins(3) + (k-0.5)*cg%dz, eps)
                           !write(77,*) i,j,k,cg%gpot(i,j,k)
                        enddo
                     enddo
                  enddo

               !close(77)
         end subroutine pot_grid2


         subroutine check_ord(order, df_dx_p, d2f_dx2_p, df_dy_p, d2f_dy2_p,& 
                  df_dz_p, d2f_dz2_p, d2f_dxdy_p, d2f_dxdz_p, d2f_dydz_p)
            implicit none
               integer,intent(in) :: order
               procedure(df_dxi),pointer :: df_dx_p, df_dy_p, df_dz_p
               procedure(d2f_dxi_2),pointer :: d2f_dx2_p, d2f_dy2_p, d2f_dz2_p
               procedure(d2f_dxi_dxj),pointer :: d2f_dxdy_p, d2f_dxdz_p, d2f_dydz_p
                  if (order==2) then
                     df_dx_p => df_dx_o2_2
                     df_dy_p => df_dy_o2_2
                     df_dz_p => df_dz_o2_2
                     d2f_dx2_p => d2f_dx2_o2_2
                     d2f_dy2_p => d2f_dy2_o2_2
                     d2f_dz2_p => d2f_dz2_o2_2
                     d2f_dxdy_p => d2f_dxdy_o2_2
                     d2f_dxdz_p => d2f_dxdz_o2_2
                     d2f_dydz_p => d2f_dydz_o2_2
                  else
                     df_dx_p => df_dx_o4_2
                     df_dy_p => df_dy_o4_2
                     df_dz_p => df_dz_o4_2
                     d2f_dx2_p => d2f_dx2_o4_2
                     d2f_dy2_p => d2f_dy2_o4_2
                     d2f_dz2_p => d2f_dz2_o4_2
                     d2f_dxdy_p => d2f_dxdy_o4_2
                     d2f_dxdz_p => d2f_dxdz_o4_2
                     d2f_dydz_p => d2f_dydz_o4_2
                  endif
         end subroutine check_ord


         subroutine grav_pot2acc_cic2(pset, cg, neighb, acc3, n)
            use constants, only: ndims, CENTER, xdim, ydim, zdim
            use grid_cont,        only: grid_container
            use particle_types, only: particle_set

            implicit none
            

            type(grid_container), pointer, intent(in) :: cg
            class(particle_set), intent(in) :: pset  !< particle list
            
            integer, intent(in) :: n
            integer :: i, j, k, c, cdim
            integer(kind=8) :: p
            
            integer(kind=8), dimension(n, ndims), intent(in) :: neighb
            integer(kind=8), dimension(n, ndims) :: neighb2
            real, dimension(n, ndims) :: dxyz
            real :: d3
            real, dimension(n, ndims), intent(out) :: acc3
            
            real(kind=8),dimension(n, 8) :: aijk, fx, fy, fz

            acc3 = 0.0
            d3 = cg%dx*cg%dy*cg%dz
            
            do i = 1, n
               do cdim = 1, ndims
                  if (pset%p(i)%pos(cdim) < cg%coord(CENTER, cdim)%r(neighb(i,cdim)-1)) then
                     neighb2(i,cdim) = neighb(i,cdim)-1
                  else
                     neighb2(i,cdim) = neighb(i,cdim)
                  endif
                  dxyz(i, cdim) = pset%p(i)%pos(cdim) - cg%coord(CENTER, cdim)%r(neighb2(i,cdim)-1)
               enddo
               aijk(i, 1) = (cg%dx - dxyz(i, xdim))*(cg%dy - dxyz(i, ydim))*(cg%dz - dxyz(i, zdim))
               aijk(i, 2) = (cg%dx - dxyz(i, xdim))*(cg%dy - dxyz(i, ydim))*         dxyz(i, zdim)
               aijk(i, 3) = (cg%dx - dxyz(i, xdim))*         dxyz(i, ydim) *(cg%dz - dxyz(i, zdim))
               aijk(i, 4) = (cg%dx - dxyz(i, xdim))*         dxyz(i, ydim) *         dxyz(i, zdim)
               aijk(i, 5) =          dxyz(i, xdim) *(cg%dy - dxyz(i, ydim))*(cg%dz - dxyz(i, zdim))
               aijk(i, 6) =          dxyz(i, xdim) *(cg%dy - dxyz(i, ydim))*         dxyz(i, zdim)
               aijk(i, 7) =          dxyz(i, xdim) *         dxyz(i, ydim) *(cg%dz - dxyz(i, zdim))
               aijk(i, 8) =          dxyz(i, xdim) *         dxyz(i, ydim) *         dxyz(i, zdim)
            enddo
            
            aijk = aijk/d3
            
            
            c = 1
            do p = 1, n
               do i = 0, 1
                  do j = 0, 1
                     do k = 0, 1
                        fx(p, c) = -(cg%gpot(neighb2(p, xdim)+1+i, neighb2(p, ydim)  +j, neighb2(p, zdim)  +k) - cg%gpot(neighb2(p, xdim)-1+i, neighb2(p, ydim)  +j, neighb2(p, zdim)  +k))
                        fy(p, c) = -(cg%gpot(neighb2(p, xdim)  +i, neighb2(p, ydim)+1+j, neighb2(p, zdim)  +k) - cg%gpot(neighb2(p, xdim)  +i, neighb2(p, ydim)-1+j, neighb2(p, zdim)  +k))
                        fz(p, c) = -(cg%gpot(neighb2(p, xdim)  +i, neighb2(p, ydim)  +j, neighb2(p, zdim)+1+k) - cg%gpot(neighb2(p, xdim)  +i, neighb2(p, ydim)  +j, neighb2(p, zdim)-1+k))
                        c = c + 1
                     enddo
                  enddo
               enddo
            enddo
            fx = 0.5*fx*cg%idx
            fy = 0.5*fy*cg%idy
            fz = 0.5*fz*cg%idz
            
            do p = 1, n
               do c = 1, 8
                  acc3(p, xdim) = acc3(p, xdim) + aijk(p, c)*fx(p, c)
                  acc3(p, ydim) = acc3(p, ydim) + aijk(p, c)*fy(p, c)
                  acc3(p, zdim) = acc3(p, zdim) + aijk(p, c)*fz(p, c)
               enddo
            enddo

         end subroutine grav_pot2acc_cic2


         subroutine potential(pset, cg, neighb, dist, n)
         use constants,    only: ndims
         use grid_cont,    only: grid_container
         implicit none
            type(grid_container), pointer, intent(in) :: cg
            class(particle_set), intent(inout) :: pset  !< particle list
            integer, intent(in) :: n
            integer(kind=8), dimension(n, ndims), intent(in) :: neighb
            real(kind=8), dimension(n, ndims), intent(in) :: dist
            integer :: i
            integer (kind=8) :: p, q, r
            real,dimension(n):: pot_x,pot_y,pot_z


            pot_x = df_dx_p(neighb, cg, n) * dist(:,1) + &
                     0.5*d2f_dx2_p(neighb, cg, n) * dist(:,1)**2
            pot_y = df_dy_p(neighb, cg, n) * dist(:,2) + &
                     0.5*d2f_dy2_p(neighb, cg, n) * dist(:,2)**2
            pot_z = df_dz_p(neighb, cg, n) * dist(:,3) + &
                     0.5*d2f_dy2_p(neighb, cg, n) * dist(:,3)**2

            do i = 1, n, 1
               p = neighb(i, 1)
               q = neighb(i, 2)
               r = neighb(i, 3)
               pset%p(i)%pot = cg%gpot(p,q,r) +sqrt(pot_x(i)**2 + pot_y(i)**2 + pot_z(i)**2)
            enddo
         end subroutine potential

         subroutine get_energy(pset, cg, neighb, dist, n, energy)
            use constants,    only: ndims
            use grid_cont,    only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               class(particle_set), intent(inout) :: pset  !< particle list
               integer :: i, j
               integer, intent(in) :: n
               integer(kind=8), dimension(n, ndims), intent(in) :: neighb
               real(kind=8), dimension(n, ndims), intent(in) :: dist
               real, intent(out) :: energy
               real :: velocity = 0.0

               call potential(pset, cg, neighb, dist, n)

               energy = 0.0

               do i=1, n
                  do j=1, ndims
                     velocity = velocity + pset%p(i)%vel(j)**2
                  enddo

                  energy = energy + 0.5*velocity + pset%p(i)%pot
                  velocity = 0.0
               enddo
         end subroutine get_energy


         subroutine get_acc2(neighb, dist, acc, cg, n)
            use constants, only : ndims, xdim, ydim, zdim
            use cg_list,      only: cg_list_element
            use grid_cont,    only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims), intent(in) :: neighb
               real(kind=8), dimension(n, ndims), intent(in) :: dist
               real(kind=8),dimension(n, ndims),intent(out) :: acc



                  acc(:, xdim) = -( df_dx_p(neighb, cg, n) + &
                              d2f_dx2_p(neighb, cg, n) * dist(:,1) + &
                              d2f_dxdy_p(neighb, cg, n) * dist(:,2) + &
                              d2f_dxdz_p(neighb, cg, n) * dist(:,3))

                  acc(:, ydim) = -( df_dy_p(neighb, cg, n) + &
                              d2f_dy2_p(neighb, cg, n) * dist(:,2) + &
                              d2f_dxdy_p(neighb, cg, n) * dist(:,1) + &
                              d2f_dydz_p(neighb, cg, n) * dist(:,3))

                  acc(:, zdim) = -( df_dz_p(neighb, cg, n) + &
                     d2f_dz2_p(neighb, cg, n) * dist(:,3) + &
                     d2f_dxdz_p(neighb, cg, n) * dist(:,1) +&
                     d2f_dydz_p(neighb, cg, n) * dist(:,2))

         end subroutine get_acc2



         

         
         subroutine get_acc_num2(pset, acc2, eps, n)
            use constants, only: ndims, xdim, ydim, zdim
            use grid_cont,  only: grid_container
            implicit none
               class(particle_set), intent(in) :: pset  !< particle list
               integer, intent(in) :: n
               real, intent(in) :: eps
               real, dimension(n, ndims), intent(out) :: acc2

               do i=1,n
                  acc2(i, xdim) = -der_x(pset%p(i)%pos, 1.0e-8, eps)
                  acc2(i, ydim) = -der_y(pset%p(i)%pos, 1.0e-8, eps)
                  acc2(i, zdim) = -der_z(pset%p(i)%pos, 1.0e-8, eps)
               enddo

         end subroutine get_acc_num2


         function der_x(pos, d, eps)
         implicit none
            real(kind=8) :: x, y, z, der_x, d, eps
            real(kind=8),dimension(1,3) :: pos
            x = pos(1,1)
            y = pos(1,2)
            z = pos(1,3)
            der_x = ( phi_pm(x+d, y, z, eps) - phi_pm(x-d, y, z, eps) ) / (2.0*d)
         end function der_x

         !Pochodna wzgledem y
         function der_y(pos, d, eps)
            implicit none
               real(kind=8) :: x, y, z, der_y, d, eps
               real(kind=8),dimension(1,3) :: pos

               x = pos(1,1)
               y = pos(1,2)
               z = pos(1,3)
               der_y = ( phi_pm(x, y+d, z, eps) - phi_pm(x, y-d, z, eps) ) / (2.0*d)
         end function der_y

         !Pochodna wzgledem z
         function der_z(pos, d, eps)
            implicit none
               real(kind=8) :: x, y, z, der_z, d, eps
               real(kind=8),dimension(1,3) :: pos

               x = pos(1,1)
               y = pos(1,2)
               z = pos(1,3)
               der_z = ( phi_pm(x, y, z+d, eps) - phi_pm(x, y, z-d, eps) ) / (2.0*d)
         end function der_z

         
         subroutine cell_nr2(pset, neighb, dist, mins, cg, n)
            use constants, only: ndims, xdim, ydim, zdim
            use grid_cont,  only: grid_container
            implicit none
               class(particle_set), intent(in) :: pset  !< particle list
               type(grid_container), pointer, intent(in) :: cg
               integer :: i, cdim
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(out) :: neighb
               real(kind=8),dimension(n, ndims),intent(out) :: dist
               real, dimension(ndims), intent(in)  :: mins
               real, dimension(ndims) :: delta_cells

               delta_cells(xdim) = cg%dx  !jak to zastapic przez samo cg?
               delta_cells(ydim) = cg%dy
               delta_cells(zdim) = cg%dz


               do i=1, n
                  do cdim = 1, ndims
                     neighb(i,cdim) = int( (pset%p(i)%pos(cdim) - mins(cdim)) / delta_cells(cdim) ) + 1
                  enddo
               enddo

               do i = 1, n
                     dist(i,:) =  pset%p(i)%pos - (neighb(i,:)-0.5) * delta_cells - mins
               enddo

               open(unit=777, file='dist.dat', status='unknown',  position='append')
                  do i=1,n
                     write(777,*) i, neighb(i,:), dist(i,:)
                  enddo
               close(777)
         end subroutine cell_nr2


         function df_dx_o2_2(neighb, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: neighb
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: df_dx_o2_2

               do i=1, n, 1
                  p = neighb(i, 1)
                  q = neighb(i, 2)
                  r = neighb(i, 3)

                  !o(R^2)
                  df_dx_o2_2(i) = ( cg%gpot(p+1, q, r) - cg%gpot(p-1, q, r) ) / (2.0*cg%dx)
               enddo
         end function df_dx_o2_2


         function df_dy_o2_2(neighb, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container 
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: neighb
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: df_dy_o2_2

               do i = 1, n, 1
                  p = neighb(i, 1)
                  q = neighb(i, 2)
                  r = neighb(i, 3)

                  !o(R^2)
                  df_dy_o2_2(i) = (cg%gpot(p, q+1, r) - cg%gpot(p, q-1, r) ) / (2.0*cg%dy)
               enddo
         end function df_dy_o2_2
         

      function df_dz_o2_2(neighb, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container 
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: neighb
               integer :: i
               integer(kind=8) :: p, q, r
               real, dimension(n), target :: df_dz_o2_2

               do i = 1, n, 1
                  p = neighb(i, 1)
                  q = neighb(i, 2)
                  r = neighb(i, 3)

                  !o(R^2)
                  df_dz_o2_2(i) = ( cg%gpot(p, q, r+1) - cg%gpot(p, q, r-1) ) / (2.0*cg%dz)
               enddo
         end function df_dz_o2_2
         


         function d2f_dx2_o2_2(neighb, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: neighb
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: d2f_dx2_o2_2

               do i = 1, n, 1
                  p = neighb(i, 1)
                  q = neighb(i, 2)
                  r = neighb(i, 3)

                  !o(R^2)
                  d2f_dx2_o2_2(i) = (cg%gpot(p+1, q, r) - 2.0*cg%gpot(p, q, r) + cg%gpot(p-1, q, r) ) / (cg%dx**2)
               enddo
         end function d2f_dx2_o2_2


         function d2f_dy2_o2_2(neighb, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: neighb
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: d2f_dy2_o2_2

               do i = 1, n, 1
                  p = neighb(i, 1)
                  q = neighb(i, 2)
                  r = neighb(i, 3)

                  !o(R^2)
                  d2f_dy2_o2_2(i) = ( cg%gpot(p, q+1, r) - 2.0*cg%gpot(p, q, r) + cg%gpot(p, q-1, r) ) / (cg%dy**2)
               enddo
         end function d2f_dy2_o2_2


         function d2f_dz2_o2_2(neighb, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: neighb
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: d2f_dz2_o2_2

               do i = 1, n, 1
                  p = neighb(i, 1)
                  q = neighb(i, 2)
                  r = neighb(i, 3)

                  !o(R^2)
                  d2f_dz2_o2_2(i) = ( cg%gpot(p, q, r+1) - 2.0*cg%gpot(p, q, r) + cg%gpot(p, q, r-1) ) / (cg%dz**2)
               enddo
         end function d2f_dz2_o2_2

         function d2f_dxdy_o2_2(neighb, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: neighb
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: d2f_dxdy_o2_2

               do i = 1, n, 1
                  p = neighb(i, 1)
                  q = neighb(i, 2)
                  r = neighb(i, 3)

                  !o(R^2)
                  d2f_dxdy_o2_2(i) = ( cg%gpot(p+1, q+1, r) - cg%gpot(p+1, q-1, r) - &
                              cg%gpot(p-1, q+1, r) + cg%gpot(p-1, q-1, r) ) / (4.0*cg%dx*cg%dy)
               enddo
         end function d2f_dxdy_o2_2


         function d2f_dxdz_o2_2(neighb, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: neighb
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: d2f_dxdz_o2_2

               do i = 1, n, 1
                  p = neighb(i, 1)
                  q = neighb(i, 2)
                  r = neighb(i, 3)

                  !o(R^2)
                  d2f_dxdz_o2_2(i) = ( cg%gpot(p+1, q, r+1) - cg%gpot(p+1, q, r-1) - &
                                 cg%gpot(p-1, q, r+1) + cg%gpot(p-1, q, r-1) ) / (4.0*cg%dx*cg%dz)
               enddo
         end function d2f_dxdz_o2_2


         function d2f_dydz_o2_2(neighb, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: neighb
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: d2f_dydz_o2_2

               do i = 1, n, 1
                  p = neighb(i, 1)
                  q = neighb(i, 2)
                  r = neighb(i, 3)

                  !o(R^2)
                  d2f_dydz_o2_2(i) = ( cg%gpot(p, q+1, r+1) - cg%gpot(p, q+1, r-1) - &
                                 cg%gpot(p, q-1, r+1) + cg%gpot(p, q-1, r-1) ) / (4.0*cg%dy*cg%dz)
               enddo
         end function d2f_dydz_o2_2


         function df_dx_o4_2(neighb, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: neighb
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: df_dx_o4_2

               do i = 1, n, 1
                  p = neighb(i, 1)
                  q = neighb(i, 2)
                  r = neighb(i, 3)

                  !o(R^4)
                  df_dx_o4_2(i) = ( 2.0* (cg%gpot(p+1, q, r) - cg%gpot(p-1, q, r) ) ) / (3.0*cg%dx) - &
                              ( cg%gpot(p+2, q, r) - cg%gpot(p-2, q, r) ) / (12.0*cg%dx)
               enddo
         end function df_dx_o4_2


         function df_dy_o4_2(neighb, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: neighb
               integer :: i
               integer(kind=8) :: p, q, r
               real, dimension(n), target :: df_dy_o4_2

               do i = 1, n, 1
                  p = neighb(i, 1)
                  q = neighb(i, 2)
                  r = neighb(i, 3)

                  !o(R^4)
                  df_dy_o4_2(i) = ( 2.0 * ( cg%gpot(p, q+1, r) - cg%gpot(p, q-1, r) ) ) / (3.0*cg%dy) - &
                        ( cg%gpot(p, q+2, r) - cg%gpot(p, q-2, r) ) / (12.0*cg%dy)
               enddo
         end function df_dy_o4_2


         function df_dz_o4_2(neighb, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: neighb
               integer :: i
               integer(kind=8) :: p, q, r
               real, dimension(n), target :: df_dz_o4_2

               do i = 1, n, 1
                  p = neighb(i, 1)
                  q = neighb(i, 2)
                  r = neighb(i, 3)

                  !o(R^4)
                  df_dz_o4_2(i) = ( 2.0* (cg%gpot(p, q, r+1) - cg%gpot(p, q, r-1) ) ) / (3.0*cg%dz) - &
                           ( cg%gpot(p, q, r+2) - cg%gpot(p, q, r-2) ) / (12.0*cg%dz)
               enddo
         end function df_dz_o4_2


         function d2f_dx2_o4_2(neighb, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: neighb
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: d2f_dx2_o4_2

               do i = 1, n, 1
                  p = neighb(i, 1)
                  q = neighb(i, 2)
                  r = neighb(i, 3)

                  !o(R^4)
                  d2f_dx2_o4_2(i) = 4.0 * ( cg%gpot(p+1, q, r) + cg%gpot(p-1, q, r) - &
                           2.0 * cg%gpot(p, q, r) ) / (3.0*cg%dx**2) - &
                           ( cg%gpot(p+2, q, r) + cg%gpot(p-2, q, r) - 2.0 * cg%gpot(p, q, r) ) / (12.0*cg%dx**2)
               enddo
         end function d2f_dx2_o4_2


         function d2f_dy2_o4_2(neighb, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: neighb
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: d2f_dy2_o4_2

               do i = 1, n, 1
                  p = neighb(i, 1)
                  q = neighb(i, 2)
                  r = neighb(i, 3)

                  !o(R^4)
                  d2f_dy2_o4_2(i) = 4.0*( cg%gpot(p, q+1, r) + cg%gpot(p, q-1, r) - &
                        2.0*cg%gpot(p, q, r) ) / (3.0*cg%dy**2) - &
                        ( cg%gpot(p, q+2, r) + cg%gpot(p, q-2, r) - 2.0*cg%gpot(p, q, r) ) / (12.0*cg%dy**2)
               enddo
         end function d2f_dy2_o4_2


         function d2f_dz2_o4_2(neighb, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: neighb
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: d2f_dz2_o4_2

               do i = 1, n, 1
                  p = neighb(i, 1)
                  q = neighb(i, 2)
                  r = neighb(i, 3)

                  !o(R^4)
                  d2f_dz2_o4_2(i) = 4.0*( cg%gpot(p, q, r+1) + cg%gpot(p, q, r-1) - &
                           2.0*cg%gpot(p, q, r) ) / (3.0*cg%dz**2) - &
                           ( cg%gpot(p, q, r+2) + cg%gpot(p, q, r-2) - 2.0*cg%gpot(p, q, r) ) / (12.0*cg%dz**2)
               enddo
         end function d2f_dz2_o4_2


         function d2f_dxdy_o4_2(neighb, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: neighb
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: d2f_dxdy_o4_2

               do i = 1, n, 1
                  p = neighb(i, 1)
                  q = neighb(i, 2)
                  r = neighb(i, 3)

                  !o(R^4)
                  d2f_dxdy_o4_2(i) = ( cg%gpot(p+1, q+1, r) + cg%gpot(p-1, q-1, r) - cg%gpot(p+1, q-1, r) - &
                              cg%gpot(p-1, q+1, r) ) / (3.0*cg%dx*cg%dy) - &
                              ( cg%gpot(p+2, q+2, r) + cg%gpot(p-2, q-2, r) - cg%gpot(p+2, q-2, r) - &
                              cg%gpot(p-2, q+2, r) ) / (48.0*cg%dx*cg%dy)
               enddo
         end function d2f_dxdy_o4_2


         function d2f_dxdz_o4_2(neighb, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: neighb
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: d2f_dxdz_o4_2

               do i = 1, n, 1
                  p = neighb(i, 1)
                  q = neighb(i, 2)
                  r = neighb(i, 3)

                  !o(R^4)
                  d2f_dxdz_o4_2(i) = ( cg%gpot(p+1, q, r+1) + cg%gpot(p-1, q, r-1) - cg%gpot(p+1, q, r-1) - &
                              cg%gpot(p-1, q, r+1) ) / (3.0*cg%dx*cg%dz) - &
                              ( cg%gpot(p+2, q, r+2) + cg%gpot(p-2, q, r-2) - cg%gpot(p+2, q, r-2) - &
                              cg%gpot(p-2, q, r+2) ) / (48.0*cg%dx*cg%dz)
               enddo
         end function d2f_dxdz_o4_2
         
         function d2f_dydz_o4_2(neighb, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: neighb
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: d2f_dydz_o4_2

               do i = 1, n, 1
                  p = neighb(i, 1)
                  q = neighb(i, 2)
                  r = neighb(i, 3)

                  !o(R^4)
                  d2f_dydz_o4_2(i) = ( cg%gpot(p, q+1, r+1) + cg%gpot(p, q-1, r-1) - cg%gpot(p, q+1, r-1) - &
                              cg%gpot(p, q-1, r+1) ) / (3.0*cg%dy*cg%dz) - &
                              ( cg%gpot(p, q+2, r+2) + cg%gpot(p, q-2, r-2) - cg%gpot(p, q+2, r-2) - &
                              cg%gpot(p, q-2, r+2) ) / (48.0*cg%dy*cg%dz)
               enddo
         end function d2f_dydz_o4_2



         
         subroutine get_ang_momentum_2(pset, n, ang_momentum)
            use constants, only : xdim, ydim, zdim
            use particle_types, only: particle_set
            implicit none
               class(particle_set), intent(in) :: pset
               integer :: i
               integer, intent(in) :: n
               real, intent(out) :: ang_momentum
               real :: L1,L2,L3

               ang_momentum = 0.0
               
               do i = 1, n, 1
                  L1 = pset%p(i)%pos(ydim) * pset%p(i)%vel(zdim) - pset%p(i)%pos(zdim) * pset%p(i)%vel(ydim)
                  L2 = pset%p(i)%pos(zdim) * pset%p(i)%vel(xdim) - pset%p(i)%pos(xdim) * pset%p(i)%vel(zdim)
                  L3 = pset%p(i)%pos(xdim) * pset%p(i)%vel(ydim) - pset%p(i)%pos(ydim) * pset%p(i)%vel(xdim)
                  ang_momentum = ang_momentum + pset%p(i)%mass * sqrt(L1**2 + L2**2 + L3**2)
               enddo

         end subroutine get_ang_momentum_2


         subroutine get_acc_mod(acc, n, a)
            use constants, only: ndims
            implicit none
            integer, intent(in) :: n
            integer  :: cdim
            real, dimension(n, ndims), intent(in) :: acc
            real, dimension(n) :: ac
            real, intent(out)  :: a

            ac = 0.0

            do cdim = 1, ndims
                  ac(:) = ac(:) + acc(:, cdim)**2
            enddo

            a = sqrt(maxval(ac))

         end subroutine get_acc_mod

  end subroutine leapfrog2ord      
      
   subroutine get_acc_pot(mass, pos, acc, n, epot)
      use constants, only: ndims
      implicit none
      integer, intent(in) :: n
      real, dimension(n), intent(in) :: mass
      real, dimension(n, ndims), intent(in) :: pos
      real, dimension(n, ndims), intent(out) :: acc
      real, intent(out) :: epot
      
      integer :: i, j
      real, dimension(ndims) :: rji, da
      
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

end module particle_integrators
