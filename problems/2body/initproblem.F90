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
#define RNG is:ie, js:je, ks:ke
module initproblem

! Initial condition for Sedov-Taylor explosion
! Written by: M. Hanasz, March 2006

   implicit none

   private
   public  :: read_problem_par, problem_initial_conditions, problem_pointers

   integer(kind=4) :: n_sn
   real            :: d0, p0, bx0, by0, bz0, Eexpl, x0, y0, z0, r0, dt_sn, r, t_sn
   

   namelist /PROBLEM_CONTROL/ d0, p0, bx0, by0, bz0, Eexpl, x0, y0, z0, r0, n_sn, dt_sn

contains
!-----------------------------------------------------------------------------
   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers
!-----------------------------------------------------------------------------
   subroutine read_problem_par

      implicit none

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine problem_initial_conditions

      use cg_leaves,    only: leaves
      use cg_list,      only: cg_list_element
      use constants,    only: xdim, ydim, zdim, LO, HI
      use dataio_pub,   only: printinfo
      use fluidindex,   only: flind
      use particle_pub, only: pset

      implicit none

      integer                          :: i, j, k, p, n_particles, particles
      real                              :: x, y, z, dtheta, e
      real,dimension(3)                :: pos_init, vel_init
      real, parameter :: pi2=6.283185307
      logical, save                   :: first_run = .true.
      type(cg_list_element),  pointer :: cgl
      n_particles=1
      
      do p = lbound(flind%all_fluids, dim=1), ubound(flind%all_fluids, dim=1)
         cgl => leaves%first
         do while (associated(cgl))
            associate(cg => cgl%cg)
               do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
                  do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
                     do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                        associate( fl => flind%all_fluids(p)%fl )
                           cg%u(fl%idn,i,j,k) = 1.0
                           cg%u(fl%imx,i,j,k) = 0.0
                           cg%u(fl%imy,i,j,k) = 0.0
                           cg%u(fl%imz,i,j,k) = 0.0
                        end associate
                     enddo
                  enddo
               enddo
            end associate
            cgl => cgl%nxt
         enddo
      enddo


      e = 0.9

      particles = 1
      write(*,*) "Particles: ", particles
      dtheta = pi2/particles
      write(*,*) "dtheta: ", dtheta

      pos_init(1) = 1.0
      pos_init(2) = 0.0
      pos_init(3) = 0.0
      
      vel_init = velocities(pos_init, e)

      if (first_run) then
         do i = 1, particles, 1
            
            call pset%add(1.0, pos_init, vel_init)
            pos_init = positions(dtheta, pos_init)
            vel_init = rotate(dtheta, vel_init)
            
         enddo
         !call printinfo('To see results type: gnuplot -p -e ''plot "nbody_out.log" u 2:3'' ')
         first_run = .false.
      endif

      contains

      function positions(dtheta, pos_init)
         implicit none
            real, dimension(3) :: positions, pos_init
            real :: dtheta
               positions = rotate(dtheta, pos_init)
      end function positions


      function velocities(pos_init, e)
         implicit none
            real, dimension(3) :: pos_init, velocities
            real :: e, r, vx, vy
            real, parameter:: mu=1.0
            
            r = sqrt(pos_init(1)**2 + pos_init(2)**2 + pos_init(3)**2)
            vx = 0.0
            vy = sqrt( (2.0 * mu * (1 - e**2) ) / (r * (r**2 + (1 - e**2) ) ) )
            velocities(1) = vx
            velocities(2) = vy
            velocities(3) = 0.0
      end function velocities
      
      function rotate (theta, vector)
         implicit none
            real, dimension(3) :: vector, rotate
            real :: theta
               rotate(1) = vector(1)*cos(theta) - vector(2)*sin(theta)
               rotate(2) = vector(1)*sin(theta) + vector(2)*cos(theta)
               rotate(3) = vector(3)
      end function rotate
   end subroutine problem_initial_conditions
!-----------------------------------------------------------------------------
end module initproblem
