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
      use constants, only: ndims, CENTER, xdim, ydim, zdim, zero
      use particle_types, only: particle_set
      use domain, only: dom
      use cg_leaves,    only: leaves
      use cg_list,      only: cg_list_element
      use grid_cont,    only: grid_container
#ifdef SELF_GRAV
      use cg_list_dataop,   only: ind_val
      use constants,        only: gp_n, gpot_n, hgpot_n, one, half, sgp_n, sgpm_n
      use named_array_list, only: qna
#endif /* SELF_GRAV */

      implicit none 
      class(particle_set), intent(inout) :: pset  !< particle list
      type(grid_container), pointer :: cg
      type(cg_list_element),  pointer  :: cgl

      interface
         function df_dxi(cells, cg, n_particles)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
            type(grid_container), pointer, intent(in) :: cg
            integer, intent(in) :: n_particles
            integer(kind=8),dimension(n_particles, ndims),intent(in) :: cells
            real,dimension(n_particles) :: df_dxi
         end function df_dxi

         function d2f_dxi_2(cells, cg, n_particles)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
            type(grid_container), pointer, intent(in) :: cg
            integer, intent(in) :: n_particles
            integer(kind=8),dimension(n_particles, ndims),intent(in) :: cells
            real,dimension(n_particles) :: d2f_dxi_2
         end function d2f_dxi_2

         function d2f_dxi_dxj(cells, cg, n_particles)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
            type(grid_container), pointer, intent(in) :: cg
            integer, intent(in) :: n_particles
            integer(kind=8),dimension(n_particles, ndims),intent(in) :: cells
            real,dimension(n_particles) :: d2f_dxi_dxj
         end function d2f_dxi_dxj
      
      end interface

      integer(kind=8), dimension(:,:), allocatable :: cells
      real, dimension(:,:), allocatable :: dist, acc, acc2, acc3
      
      real, intent(in) :: t_glob, dt_tot
      real, dimension(:), allocatable :: mass, mins, maxs

      real :: t_end, dt, t, dth, t_out, eta, eps, a, eps2, energy, init_energy, d_energy = 0.0, ang_momentum, init_ang_mom, d_ang_momentum = 0.0!, zero
      integer :: nsteps, n, lun_out, i, order, j, k, cdim
      real, parameter :: dt_out = 0.05          ! time interval between output of snapshots
      logical :: save_potential,finish,external_pot,var_timestep

      procedure(df_dxi),pointer :: df_dx_p, df_dy_p, df_dz_p
      procedure(d2f_dxi_2),pointer :: d2f_dx2_p, d2f_dy2_p, d2f_dz2_p
      procedure(d2f_dxi_dxj),pointer :: d2f_dxdy_p, d2f_dxdz_p, d2f_dydz_p
      
      
      open(newunit=lun_out, file='leapfrog_out.log', status='unknown',  position='append')
      
      n = size(pset%p, dim=1)


      allocate(mass(n), acc(n, ndims), acc2(n, ndims), acc3(n, ndims), cells(n, ndims), dist(n, ndims), mins(ndims), maxs(ndims))



      mass(:) = pset%p(:)%mass

      !do ndim = xdim, zdim
      !   pos(:, ndim) = pset%p(:)%pos(ndim)
      !   vel(:, ndim) = pset%p(:)%vel(ndim)
      !enddo


      t = t_glob
      t_end = t + dt_tot
      t_out = t_glob + dt_out



      mins(:) = dom%edge(:,1)

      maxs(:) = dom%edge(:,2)


      order = 2


      eta = 5.0 !1.0
      eps = 1.0e-4
      eps2 = 0.00


      cgl => leaves%first
      do while (associated(cgl))
            cg => cgl%cg
            cgl => cgl%nxt
      enddo


      !obliczenie zewnętrznego potencjalu na siatce
      external_pot = .true.
      !external_pot = .false.

      if (external_pot) then
         call pot2grid(cg, mins, eps2)
         write(*,*) "Obliczono potencjal zewnetrzny"
      endif



      !save_potential = .true.
      save_potential = .false.
      !finish         = .true.
      finish         = .false.

      if(save_potential) then
         open(unit=88, file='potencjal.dat')

         do i=lbound(cg%gpot,dim=1),ubound(cg%gpot,dim=1)
            do j=lbound(cg%gpot,dim=2),ubound(cg%gpot,dim=2)
               do k=lbound(cg%gpot,dim=3),ubound(cg%gpot,dim=3)
                  write(88,*) i, j, k, cg%coord(CENTER, xdim)%r(i), &
                              cg%coord(CENTER, xdim)%r(i), cg%coord(CENTER, xdim)%r(i), &
                              cg%gpot(i,j,k)
               enddo
            enddo
            write(88,*)
         enddo
         close(88)

         if(finish) stop
      endif


      call get_ang_momentum_2(pset, n, ang_momentum)
      init_ang_mom = ang_momentum


      call find_cells(pset, cells, dist, mins, cg, n)


      call check_ord(order, df_dx_p, d2f_dx2_p, df_dy_p, d2f_dy2_p,& 
                        df_dz_p, d2f_dz2_p, d2f_dxdy_p, d2f_dxdz_p, d2f_dydz_p)




      !przyspieszenie interpolowane
      call get_acc_int(cells, dist, acc, cg, n)



      call get_energy(pset, cg, cells, dist, n, energy)
      init_energy = energy

      if (external_pot) then
         call get_acc_model(pset, acc2, eps, n)
      else
         acc2(:,:) = 0.0
      endif
      
      call get_acc_cic_o4(pset, cg, cells, acc3, n)


      call get_acc_max(acc3, n, a)
      write(*,*) "get_acc_max, a=", a


      !var_timestep=.true.
      var_timestep=.false.
      !timestep
      if (var_timestep) then
         dt = sqrt(2.0*eta*eps/a)               !variable
         write(*,*) "leapfrog: dt = ",dt
      else
         dt = 0.1
         !dt = 0.01                            !constant
      endif
      dth = dt/2.0


      nsteps = 0



      !main loop
      do while (t < t_end)


         if (t + dt > t_end) then
            dt = t_end - t
            dth = 0.5 * dt
         endif

         !if (t >=t_out) then
            do i=1, n
               !write(lun_out,*) "# t=", t
               write(lun_out, '(I3,1X,19(E13.6,1X))') i, t, dt, mass(i), pset%p(i)%pos, acc(i,:), acc2(i,:), acc3(i,:), energy, d_energy, ang_momentum, d_ang_momentum
            enddo
            !t_out = t_out + dt_out
         !endif


         !1.kick(dth)
         call kick(pset, acc3, dth, n)
         !call kick(pset, [zero,zero,zero], dth, n)

         !2.drift(dt)         
         call drift(pset, dt, n)


         !ekstrapolacja potencjalu w przypadku samograwitacji
         !czegos tu brakuje :/
!#ifdef SELF_GRAV
!         call leaves%q_copy(qna%ind(sgp_n), qna%ind(sgpm_n))
!         call leaves%q_lin_comb([ ind_val(qna%ind(gp_n), 1.), ind_val(qna%ind(sgp_n), one+dt),  ind_val(qna%ind(sgpm_n), -dt) ], qna%ind(gpot_n))
!         call leaves%q_lin_comb([ ind_val(qna%ind(gp_n), 1.), ind_val(qna%ind(sgp_n), one+dth), ind_val(qna%ind(sgpm_n), -dth)], qna%ind(hgpot_n))
!         write(*,*) "Self_GRAV!"
!#endif /* SELF_GRAV */


         call find_cells(pset, cells, dist, mins, cg, n)                !finding cells

         !3.acceleration + |a|
         call get_acc_int(cells, dist, acc, cg, n)                      !Lagrange polynomials acceleration


         if (external_pot) then
            call get_acc_model(pset, acc2, eps, n)                       !centered finite differencing acceleration (if gravitational potential is known explicite)
         endif


         call get_acc_cic_o4(pset, cg, cells, acc3, n)                     !CIC acceleration
         call get_acc_max(acc3, n, a)                                    !max(|a_i|)


         !4.kick(dth)
         !call kick(pset, acc, dth, n)                                  !zakomentowac dla ruchu po prostej
         call kick(pset, acc3, dth, n)
         !call kick(pset, [zero,zero,zero], dth, n)                     !odkomentowac dla ruchu po prostej


         call get_energy(pset, cg, cells, dist, n, energy)
         d_energy = log(abs((energy - init_energy)/init_energy))


         call get_ang_momentum_2(pset, n, ang_momentum)

         !write(*,*) ang_momentum, init_ang_mom
         !d_ang_momentum = log(abs((ang_momentum - init_ang_mom)/init_ang_mom))

         !5.t
         t = t + dt
         !6.dt		!dt[n+1]
         if (var_timestep) then
            dt = sqrt(2.0*eta*eps/a)
            dth = 0.5*dt
         endif

         nsteps = nsteps + 1

      end do

      write(*,*) "Leapfrog: nsteps=", nsteps



      deallocate (acc, acc2, acc3, cells, dist, mins, maxs)
      close(lun_out)


      contains

         !Kick
         subroutine kick(pset, acc, t, n)
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
         end subroutine kick

         !Drift
         subroutine drift(pset, t, n)
            use particle_types, only: particle_set
            implicit none
            class(particle_set), intent(inout) :: pset  !< particle list
            real, intent(in) :: t
            integer :: i
            integer, intent(in) :: n

            do i=1,n
               pset%p(i)%pos = pset%p(i)%pos + pset%p(i)%vel * t
            enddo
         end subroutine drift


         function phi_pm(x, y, z, eps)
            implicit none
               real,intent(in) :: x, y, z, eps
               real::r, phi_pm, G,M 
                  G = 1.0
                  M = 1.0
                  r = sqrt(x**2 + y**2 + z**2 + eps**2)

                  phi_pm = -G*M / r
         end function phi_pm


         subroutine pot2grid(cg, mins, eps2)
            use constants, only: ndims,xdim,ydim,CENTER
            use grid_cont, only: grid_container
            implicit none
               type(grid_container), pointer, intent(inout) :: cg
               integer :: i, j, k
               real,dimension(ndims),intent(in) :: mins
               real, intent(in) :: eps2

               open(unit=77,file='potencjal.dat')


                  do i = lbound(cg%gpot, dim=1), ubound(cg%gpot, dim=1)
                     do j = lbound(cg%gpot, dim=2), ubound(cg%gpot, dim=2)
                        do k = lbound(cg%gpot, dim=3), ubound(cg%gpot, dim=3)
                           !cg%gpot(i, j, k) = phi_pm(mins(1) + (i-0.5)*cg%dx, &
                           !                          mins(2) + (j-0.5)*cg%dy, &
                           !                          mins(3) + (k-0.5)*cg%dz, eps)
                           cg%gpot(i,j,k) = phi_pm(cg%coord(CENTER,xdim)%r(i),&
                                                   cg%coord(CENTER,ydim)%r(j),&
                                                   cg%coord(CENTER,zdim)%r(k),eps2)
                           write(77,*) cg%coord(CENTER,xdim)%r(i),cg%coord(CENTER,ydim)%r(j),&
                                       cg%coord(CENTER,zdim)%r(k),cg%gpot(i,j,k)
                        enddo
                     enddo
                  enddo

               close(77)
         end subroutine pot2grid


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


         subroutine get_acc_cic(pset, cg, cells, acc3, n)
            use constants, only: ndims, CENTER, xdim, ydim, zdim, half
            use grid_cont,        only: grid_container
            use particle_types, only: particle_set

            implicit none
            type(grid_container), pointer, intent(in) :: cg
            class(particle_set), intent(in) :: pset  !< particle list

            integer, intent(in) :: n
            integer :: i, j, k, c, cdim
            integer(kind=8) :: p
            integer(kind=8), dimension(n, ndims), intent(in) :: cells
            integer(kind=8), dimension(n, ndims) :: cells2
            real, dimension(n, ndims) :: dxyz
            real, dimension(n, ndims), intent(out) :: acc3
            real(kind=8), dimension(n, 8) :: aijk, fx, fy, fz

            acc3 = 0.0


            do i = 1, n
               do cdim = xdim, ndims
                  if (pset%p(i)%pos(cdim) < cg%coord(CENTER, cdim)%r(cells(i,cdim))) then
                     cells2(i,cdim) = cells(i,cdim)-1
                  else
                     cells2(i,cdim) = cells(i,cdim)!+1
                  endif
                  dxyz(i, cdim) = abs(pset%p(i)%pos(cdim) - cg%coord(CENTER, cdim)%r(cells2(i,cdim)))

               enddo
               aijk(i, 1) = (cg%dx - dxyz(i, xdim))*(cg%dy - dxyz(i, ydim))*(cg%dz - dxyz(i, zdim)) !a(i  ,j  ,k  )
               aijk(i, 2) = (cg%dx - dxyz(i, xdim))*(cg%dy - dxyz(i, ydim))*         dxyz(i, zdim)  !a(i+1,j  ,k  )
               aijk(i, 3) = (cg%dx - dxyz(i, xdim))*         dxyz(i, ydim) *(cg%dz - dxyz(i, zdim)) !a(i  ,j+1,k  )
               aijk(i, 4) = (cg%dx - dxyz(i, xdim))*         dxyz(i, ydim) *         dxyz(i, zdim)  !a(i  ,j  ,k+1)
               aijk(i, 5) =          dxyz(i, xdim) *(cg%dy - dxyz(i, ydim))*(cg%dz - dxyz(i, zdim)) !a(i+1,j+1,k  )
               aijk(i, 6) =          dxyz(i, xdim) *(cg%dy - dxyz(i, ydim))*         dxyz(i, zdim)  !a(i  ,j+1,k+1)
               aijk(i, 7) =          dxyz(i, xdim) *         dxyz(i, ydim) *(cg%dz - dxyz(i, zdim)) !a(i+1,j  ,k+1)
               aijk(i, 8) =          dxyz(i, xdim) *         dxyz(i, ydim) *         dxyz(i, zdim)  !a(i+1,j+1,k+1)
            enddo

            aijk = aijk/cg%dvol




            do p = 1, n
               c = 1
               do i = 0, 1
                  do j = 0, 1
                     do k = 0, 1
                        fx(p, c) = -(cg%gpot(cells2(p, xdim)+1+i, cells2(p, ydim)  +j, cells2(p, zdim)  +k) - cg%gpot(cells2(p, xdim)-1+i, cells2(p, ydim)  +j, cells2(p, zdim)  +k))
                        fy(p, c) = -(cg%gpot(cells2(p, xdim)  +i, cells2(p, ydim)+1+j, cells2(p, zdim)  +k) - cg%gpot(cells2(p, xdim)  +i, cells2(p, ydim)-1+j, cells2(p, zdim)  +k))
                        fz(p, c) = -(cg%gpot(cells2(p, xdim)  +i, cells2(p, ydim)  +j, cells2(p, zdim)+1+k) - cg%gpot(cells2(p, xdim)  +i, cells2(p, ydim)  +j, cells2(p, zdim)-1+k))
                        c = c + 1
                     enddo
                  enddo
               enddo
            enddo

            fx = half*fx*cg%idx
            fy = half*fy*cg%idy
            fz = half*fz*cg%idz

            do p = 1, n
               do c = 1, 8
                  acc3(p, xdim) = acc3(p, xdim) + aijk(p, c)*fx(p, c)
                  acc3(p, ydim) = acc3(p, ydim) + aijk(p, c)*fy(p, c)
                  acc3(p, zdim) = acc3(p, zdim) + aijk(p, c)*fz(p, c)
               enddo
            enddo

         end subroutine get_acc_cic


         subroutine get_acc_cic_o4(pset, cg, cells, acc3, n)
            use constants, only: ndims, CENTER, xdim, ydim, zdim, half
            use grid_cont,        only: grid_container
            use particle_types, only: particle_set

            implicit none
            type(grid_container), pointer, intent(in) :: cg
            class(particle_set), intent(in) :: pset  !< particle list

            integer, intent(in) :: n
            integer :: i, j, k, c, cdim
            integer(kind=8) :: p
            integer(kind=8), dimension(n, ndims), intent(in) :: cells
            integer(kind=8), dimension(n, ndims) :: cells2
            real, dimension(n, ndims) :: dxyz
            real, dimension(n, ndims), intent(out) :: acc3
            real(kind=8), dimension(n, 8) :: aijk, fx, fy, fz

            acc3 = 0.0


            do i = 1, n
               do cdim = xdim, ndims
                  if (pset%p(i)%pos(cdim) < cg%coord(CENTER, cdim)%r(cells(i,cdim))) then
                     cells2(i,cdim) = cells(i,cdim)-1
                  else
                     cells2(i,cdim) = cells(i,cdim)!+1
                  endif
                  dxyz(i, cdim) = abs(pset%p(i)%pos(cdim) - cg%coord(CENTER, cdim)%r(cells2(i,cdim)))

               enddo
               aijk(i, 1) = (cg%dx - dxyz(i, xdim))*(cg%dy - dxyz(i, ydim))*(cg%dz - dxyz(i, zdim)) !a(i  ,j  ,k  )
               aijk(i, 2) = (cg%dx - dxyz(i, xdim))*(cg%dy - dxyz(i, ydim))*         dxyz(i, zdim)  !a(i+1,j  ,k  )
               aijk(i, 3) = (cg%dx - dxyz(i, xdim))*         dxyz(i, ydim) *(cg%dz - dxyz(i, zdim)) !a(i  ,j+1,k  )
               aijk(i, 4) = (cg%dx - dxyz(i, xdim))*         dxyz(i, ydim) *         dxyz(i, zdim)  !a(i  ,j  ,k+1)
               aijk(i, 5) =          dxyz(i, xdim) *(cg%dy - dxyz(i, ydim))*(cg%dz - dxyz(i, zdim)) !a(i+1,j+1,k  )
               aijk(i, 6) =          dxyz(i, xdim) *(cg%dy - dxyz(i, ydim))*         dxyz(i, zdim)  !a(i  ,j+1,k+1)
               aijk(i, 7) =          dxyz(i, xdim) *         dxyz(i, ydim) *(cg%dz - dxyz(i, zdim)) !a(i+1,j  ,k+1)
               aijk(i, 8) =          dxyz(i, xdim) *         dxyz(i, ydim) *         dxyz(i, zdim)  !a(i+1,j+1,k+1)
            enddo

            aijk = aijk/cg%dvol




            do p = 1, n
               c = 1
               do i = 0, 1
                  do j = 0, 1
                     do k = 0, 1
                        fx(p, c) = -( (2.0/3.0)*(cg%gpot(cells2(p, xdim)+1+i, cells2(p, ydim)  +j, cells2(p, zdim)  +k) - cg%gpot(cells2(p, xdim)-1+i, cells2(p, ydim)  +j, cells2(p, zdim)  +k))*cg%idx - &
                                    (1.0/12.0) * (cg%gpot(cells2(p, xdim)+2+i, cells2(p, ydim)  +j, cells2(p, zdim)  +k) - cg%gpot(cells2(p, xdim)-2+i, cells2(p, ydim)  +j, cells2(p, zdim)  +k))*cg%idx)
                        fy(p, c) = -( (2.0/3.0)*(cg%gpot(cells2(p, xdim)  +i, cells2(p, ydim)+1+j, cells2(p, zdim)  +k) - cg%gpot(cells2(p, xdim)  +i, cells2(p, ydim)-1+j, cells2(p, zdim)  +k))*cg%idy - &
                                    (1.0/12.0) * (cg%gpot(cells2(p, xdim)  +i, cells2(p, ydim)+2+j, cells2(p, zdim)  +k) - cg%gpot(cells2(p, xdim)  +i, cells2(p, ydim)-2+j, cells2(p, zdim)  +k))*cg%idy)
                        fz(p, c) = -( (2.0/3.0)*(cg%gpot(cells2(p, xdim)  +i, cells2(p, ydim)  +j, cells2(p, zdim)+1+k) - cg%gpot(cells2(p, xdim)  +i, cells2(p, ydim)  +j, cells2(p, zdim)-1+k))*cg%idz - &
                                    (1.0/12.0) * (cg%gpot(cells2(p, xdim)  +i, cells2(p, ydim)  +j, cells2(p, zdim)+2+k) - cg%gpot(cells2(p, xdim)  +i, cells2(p, ydim)  +j, cells2(p, zdim)-2+k))*cg%idz)
                        c = c + 1
                     enddo
                  enddo
               enddo
            enddo

            !fx = half*fx*cg%idx
            !fy = half*fy*cg%idy
            !fz = half*fz*cg%idz

            do p = 1, n
               do c = 1, 8
                  acc3(p, xdim) = acc3(p, xdim) + aijk(p, c)*fx(p, c)
                  acc3(p, ydim) = acc3(p, ydim) + aijk(p, c)*fy(p, c)
                  acc3(p, zdim) = acc3(p, zdim) + aijk(p, c)*fz(p, c)
               enddo
            enddo

         end subroutine get_acc_cic_o4

         subroutine potential(pset, cg, cells, dist, n)
            use constants,    only: ndims, half, xdim, ydim, zdim
            use grid_cont,    only: grid_container
            implicit none
            type(grid_container), pointer, intent(in) :: cg
            class(particle_set), intent(inout) :: pset  !< particle list
            integer, intent(in) :: n
            integer(kind=8), dimension(n, ndims), intent(in) :: cells
            real(kind=8), dimension(n, ndims), intent(in) :: dist
            integer :: i
            integer (kind=8) :: p, q, r
            real,dimension(n):: pot_x,pot_y,pot_z,dpot,d2pot


            !pot_x = df_dx_p(cells, cg, n) * dist(:,1) + &
            !         0.5*d2f_dx2_p(cells, cg, n) * dist(:,1)**2
            !pot_y = df_dy_p(cells, cg, n) * dist(:,2) + &
            !         0.5*d2f_dy2_p(cells, cg, n) * dist(:,2)**2
            !pot_z = df_dz_p(cells, cg, n) * dist(:,3) + &
            !         0.5*d2f_dy2_p(cells, cg, n) * dist(:,3)**2
                     
            dpot = df_dx_p(cells, cg, n) * dist(:, xdim) + &
                   df_dy_p(cells, cg, n) * dist(:, ydim) + &
                   df_dz_p(cells, cg, n) * dist(:, zdim)
            
            d2pot = d2f_dx2_p(cells, cg, n) * dist(:, xdim)**2 + &
                    d2f_dy2_p(cells, cg, n) * dist(:, ydim)**2 + &
                    d2f_dz2_p(cells, cg, n) * dist(:, zdim)**2 + &
                    2.0*d2f_dxdy_p(cells, cg, n) * dist(:, xdim)*dist(:, ydim) + &
                    2.0*d2f_dxdz_p(cells, cg, n) * dist(:, xdim)*dist(:, zdim)
            
            do i = 1, n, 1
               p = cells(i, xdim)
               q = cells(i, ydim)
               r = cells(i, zdim)
               !pset%p(i)%pot = cg%gpot(p,q,r) +sqrt(pot_x(i)**2 + pot_y(i)**2 + pot_z(i)**2)
               pset%p(i)%pot = cg%gpot(p,q,r) + dpot(i) + half*d2pot(i)
            enddo

           ! pot=cg%gpot(p,q,r) + dpot+half*d2pot
         end subroutine potential

         subroutine get_energy(pset, cg, cells, dist, n, energy)
            use constants,    only: ndims
            use grid_cont,    only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               class(particle_set), intent(inout) :: pset  !< particle list
               integer :: i, j
               integer, intent(in) :: n
               integer(kind=8), dimension(n, ndims), intent(in) :: cells
               real(kind=8), dimension(n, ndims), intent(in) :: dist
               real, intent(out) :: energy
               real :: velocity = 0.0

               call potential(pset, cg, cells, dist, n)

               energy = 0.0

               do i=1, n
                  do j=1, ndims
                     velocity = velocity + pset%p(i)%vel(j)**2
                  enddo

                  energy = energy + 0.5*velocity + pset%p(i)%pot
                  velocity = 0.0
               enddo
         end subroutine get_energy


         subroutine get_acc_int(cells, dist, acc, cg, n)
            use constants, only : ndims, xdim, ydim, zdim
            use cg_list,      only: cg_list_element
            use grid_cont,    only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims), intent(in) :: cells
               real(kind=8), dimension(n, ndims), intent(in) :: dist
               real(kind=8),dimension(n, ndims),intent(out) :: acc




                  acc(:, xdim) = -( df_dx_p(cells, cg, n) + &
                              d2f_dx2_p(cells, cg, n) * dist(:,1) + &
                              d2f_dxdy_p(cells, cg, n) * dist(:,2) + &
                              d2f_dxdz_p(cells, cg, n) * dist(:,3))

                  acc(:, ydim) = -( df_dy_p(cells, cg, n) + &
                              d2f_dy2_p(cells, cg, n) * dist(:,2) + &
                              d2f_dxdy_p(cells, cg, n) * dist(:,1) + &
                              d2f_dydz_p(cells, cg, n) * dist(:,3))

                  acc(:, zdim) = -( df_dz_p(cells, cg, n) + &
                              d2f_dz2_p(cells, cg, n) * dist(:,3) + &
                              d2f_dxdz_p(cells, cg, n) * dist(:,1) +&
                              d2f_dydz_p(cells, cg, n) * dist(:,2))
                  !write(*,*) "acc_int: ",acc(:,xdim), acc(:,ydim), acc(:,zdim)

         end subroutine get_acc_int


         subroutine get_acc_model(pset, acc2, eps, n)
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

         end subroutine get_acc_model


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

         
         subroutine find_cells(pset, cells, dist, mins, cg, n)
            use constants, only: ndims, xdim, ydim, zdim, CENTER
            use grid_cont,  only: grid_container
            implicit none
               class(particle_set), intent(in) :: pset  !< particle list
               type(grid_container), pointer, intent(in) :: cg
               integer :: i, cdim
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(out) :: cells
               real(kind=8),dimension(n, ndims),intent(out) :: dist
               real, dimension(ndims), intent(in)  :: mins



               do i=1, n
                  do cdim = xdim, ndims
                     cells(i,cdim) = int( 0.5 + (pset%p(i)%pos(cdim) - cg%coord(CENTER,cdim)%r(0)) / cg%dl(cdim) )! + 1
                  enddo
               enddo

               do i=1,n
                  do cdim = xdim, ndims
                     dist(i, cdim) = pset%p(i)%pos(cdim) - ( cg%coord(CENTER, cdim)%r(0) + cells(i,cdim) * cg%dl(cdim) )
                  enddo
               enddo

               do i=1,n
                  do cdim = xdim, ndims
                     if (pset%p(i)%pos(cdim)<mins(cdim)) then
                        write(*,*) "przekroczono"
                        stop
                     endif
                  enddo
               enddo

               open(unit=777, file='dist.dat', status='unknown',  position='append')
                  do i=1,n
                     write(777,*) i, cells(i,:), dist(i,:)
                  enddo
               close(777)
         end subroutine find_cells


         function df_dx_o2(cells, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: cells
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: df_dx_o2

               do i=1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)
                  !write(*,*) i, p,q,r, (cg%gpot(p+1, q, r) - cg%gpot(p-1, q, r) ) / (2.0*cg%dx)

                  !o(R^2)
                  df_dx_o2(i) = ( cg%gpot(p+1, q, r) - cg%gpot(p-1, q, r) ) / (2.0*cg%dx)
               enddo
         end function df_dx_o2
         
         function df_dx_o2_1(cells, cg)!, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               !integer, intent(in) :: n
               integer,dimension(ndims),intent(in) :: cells
               integer :: i
               integer(kind=8) :: p, q, r
               real :: df_dx_o2_1

               
                  p = cells(1)
                  q = cells(2)
                  r = cells(3)

                  !o(R^2)
                  df_dx_o2_1 = ( cg%gpot(p+1, q, r) - cg%gpot(p-1, q, r) ) / (2.0*cg%dx)
           
         end function df_dx_o2_1


         function df_dy_o2(cells, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container 
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: cells
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: df_dy_o2

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^2)
                  df_dy_o2(i) = (cg%gpot(p, q+1, r) - cg%gpot(p, q-1, r) ) / (2.0*cg%dy)
               enddo
         end function df_dy_o2
         
         function df_dy_o2_1(cells, cg)!, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container 
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               !integer, intent(in) :: n
               integer,dimension(ndims),intent(in) :: cells
               integer :: i
               integer(kind=8) :: p, q, r
               real:: df_dy_o2_1

               
                  p = cells(1)
                  q = cells(2)
                  r = cells(3)

                  !o(R^2)
                  df_dy_o2_1 = (cg%gpot(p, q+1, r) - cg%gpot(p, q-1, r) ) / (2.0*cg%dy)
               
         end function df_dy_o2_1
         

      function df_dz_o2(cells, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container 
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: cells
               integer :: i
               integer(kind=8) :: p, q, r
               real, dimension(n), target :: df_dz_o2

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^2)
                  df_dz_o2(i) = ( cg%gpot(p, q, r+1) - cg%gpot(p, q, r-1) ) / (2.0*cg%dz)
               enddo
         end function df_dz_o2
         
         function df_dz_o2_1(cells, cg)!, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container 
            implicit none
               type(grid_container), pointer, intent(in) :: cg
              ! integer, intent(in) :: n
               integer,dimension(ndims),intent(in) :: cells
               integer :: i
               integer(kind=8) :: p, q, r
               real :: df_dz_o2_1

               
                  p = cells(1)
                  q = cells(2)
                  r = cells(3)

                  !o(R^2)
                  df_dz_o2_1 = ( cg%gpot(p, q, r+1) - cg%gpot(p, q, r-1) ) / (2.0*cg%dz)
            
         end function df_dz_o2_1
         


         function d2f_dx2_o2(cells, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: cells
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: d2f_dx2_o2

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^2)
                  d2f_dx2_o2(i) = (cg%gpot(p+1, q, r) - 2.0*cg%gpot(p, q, r) + cg%gpot(p-1, q, r) ) / (cg%dx**2)
               enddo
         end function d2f_dx2_o2


         function d2f_dy2_o2(cells, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: cells
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: d2f_dy2_o2

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^2)
                  d2f_dy2_o2(i) = ( cg%gpot(p, q+1, r) - 2.0*cg%gpot(p, q, r) + cg%gpot(p, q-1, r) ) / (cg%dy**2)
               enddo
         end function d2f_dy2_o2


         function d2f_dz2_o2(cells, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: cells
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: d2f_dz2_o2

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^2)
                  d2f_dz2_o2(i) = ( cg%gpot(p, q, r+1) - 2.0*cg%gpot(p, q, r) + cg%gpot(p, q, r-1) ) / (cg%dz**2)
               enddo
         end function d2f_dz2_o2

         function d2f_dxdy_o2(cells, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: cells
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: d2f_dxdy_o2

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^2)
                  d2f_dxdy_o2(i) = ( cg%gpot(p+1, q+1, r) - cg%gpot(p+1, q-1, r) - &
                              cg%gpot(p-1, q+1, r) + cg%gpot(p-1, q-1, r) ) / (4.0*cg%dx*cg%dy)
               enddo
         end function d2f_dxdy_o2


         function d2f_dxdz_o2(cells, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: cells
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: d2f_dxdz_o2

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^2)
                  d2f_dxdz_o2(i) = ( cg%gpot(p+1, q, r+1) - cg%gpot(p+1, q, r-1) - &
                                 cg%gpot(p-1, q, r+1) + cg%gpot(p-1, q, r-1) ) / (4.0*cg%dx*cg%dz)
               enddo
         end function d2f_dxdz_o2


         function d2f_dydz_o2(cells, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: cells
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: d2f_dydz_o2

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^2)
                  d2f_dydz_o2(i) = ( cg%gpot(p, q+1, r+1) - cg%gpot(p, q+1, r-1) - &
                                 cg%gpot(p, q-1, r+1) + cg%gpot(p, q-1, r-1) ) / (4.0*cg%dy*cg%dz)
               enddo
         end function d2f_dydz_o2


         function df_dx_o4(cells, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: cells
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: df_dx_o4

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^4)
                  df_dx_o4(i) = 2.0* (cg%gpot(p+1, q, r) - cg%gpot(p-1, q, r) ) / (3.0*cg%dx) - &
                              ( cg%gpot(p+2, q, r) - cg%gpot(p-2, q, r) ) / (12.0*cg%dx)
               enddo
         end function df_dx_o4


         function df_dy_o4(cells, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: cells
               integer :: i
               integer(kind=8) :: p, q, r
               real, dimension(n), target :: df_dy_o4

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^4)
                  df_dy_o4(i) = 2.0 * ( cg%gpot(p, q+1, r) - cg%gpot(p, q-1, r) ) / (3.0*cg%dy) - &
                        ( cg%gpot(p, q+2, r) - cg%gpot(p, q-2, r) ) / (12.0*cg%dy)
               enddo
         end function df_dy_o4


         function df_dz_o4(cells, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: cells
               integer :: i
               integer(kind=8) :: p, q, r
               real, dimension(n), target :: df_dz_o4

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^4)
                  df_dz_o4(i) = 2.0* (cg%gpot(p, q, r+1) - cg%gpot(p, q, r-1) ) / (3.0*cg%dz) - &
                           ( cg%gpot(p, q, r+2) - cg%gpot(p, q, r-2) ) / (12.0*cg%dz)
               enddo
         end function df_dz_o4


         function d2f_dx2_o4(cells, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: cells
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: d2f_dx2_o4

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^4)
                  d2f_dx2_o4(i) = 4.0 * ( cg%gpot(p+1, q, r) + cg%gpot(p-1, q, r) - &
                           2.0 * cg%gpot(p, q, r) ) / (3.0*cg%dx**2) - &
                           ( cg%gpot(p+2, q, r) + cg%gpot(p-2, q, r) - 2.0 * cg%gpot(p, q, r) ) / (12.0*cg%dx**2)
               enddo
         end function d2f_dx2_o4


         function d2f_dy2_o4(cells, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: cells
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: d2f_dy2_o4

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^4)
                  d2f_dy2_o4(i) = 4.0*( cg%gpot(p, q+1, r) + cg%gpot(p, q-1, r) - &
                        2.0*cg%gpot(p, q, r) ) / (3.0*cg%dy**2) - &
                        ( cg%gpot(p, q+2, r) + cg%gpot(p, q-2, r) - 2.0*cg%gpot(p, q, r) ) / (12.0*cg%dy**2)
               enddo
         end function d2f_dy2_o4


         function d2f_dz2_o4(cells, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: cells
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: d2f_dz2_o4

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^4)
                  d2f_dz2_o4(i) = 4.0*( cg%gpot(p, q, r+1) + cg%gpot(p, q, r-1) - &
                           2.0*cg%gpot(p, q, r) ) / (3.0*cg%dz**2) - &
                           ( cg%gpot(p, q, r+2) + cg%gpot(p, q, r-2) - 2.0*cg%gpot(p, q, r) ) / (12.0*cg%dz**2)
               enddo
         end function d2f_dz2_o4


         function d2f_dxdy_o4(cells, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: cells
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: d2f_dxdy_o4

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^4)
                  d2f_dxdy_o4(i) = ( cg%gpot(p+1, q+1, r) + cg%gpot(p-1, q-1, r) - cg%gpot(p+1, q-1, r) - &
                              cg%gpot(p-1, q+1, r) ) / (3.0*cg%dx*cg%dy) - &
                              ( cg%gpot(p+2, q+2, r) + cg%gpot(p-2, q-2, r) - cg%gpot(p+2, q-2, r) - &
                              cg%gpot(p-2, q+2, r) ) / (48.0*cg%dx*cg%dy)
               enddo
         end function d2f_dxdy_o4


         function d2f_dxdz_o4(cells, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: cells
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: d2f_dxdz_o4

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^4)
                  d2f_dxdz_o4(i) = ( cg%gpot(p+1, q, r+1) + cg%gpot(p-1, q, r-1) - cg%gpot(p+1, q, r-1) - &
                              cg%gpot(p-1, q, r+1) ) / (3.0*cg%dx*cg%dz) - &
                              ( cg%gpot(p+2, q, r+2) + cg%gpot(p-2, q, r-2) - cg%gpot(p+2, q, r-2) - &
                              cg%gpot(p-2, q, r+2) ) / (48.0*cg%dx*cg%dz)
               enddo
         end function d2f_dxdz_o4
         
         function d2f_dydz_o4(cells, cg, n)
            use constants, only: ndims
            use grid_cont,  only: grid_container
            implicit none
               type(grid_container), pointer, intent(in) :: cg
               integer, intent(in) :: n
               integer(kind=8),dimension(n, ndims),intent(in) :: cells
               integer :: i
               integer(kind=8) :: p, q, r
               real,dimension(n),target :: d2f_dydz_o4

               do i = 1, n, 1
                  p = cells(i, 1)
                  q = cells(i, 2)
                  r = cells(i, 3)

                  !o(R^4)
                  d2f_dydz_o4(i) = ( cg%gpot(p, q+1, r+1) + cg%gpot(p, q-1, r-1) - cg%gpot(p, q+1, r-1) - &
                              cg%gpot(p, q-1, r+1) ) / (3.0*cg%dy*cg%dz) - &
                              ( cg%gpot(p, q+2, r+2) + cg%gpot(p, q-2, r-2) - cg%gpot(p, q+2, r-2) - &
                              cg%gpot(p, q-2, r+2) ) / (48.0*cg%dy*cg%dz)
               enddo
         end function d2f_dydz_o4



         
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


         subroutine get_acc_max(acc, n, a)
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

         end subroutine get_acc_max

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
