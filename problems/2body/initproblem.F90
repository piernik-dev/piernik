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

   real            :: d0, r

   namelist /PROBLEM_CONTROL/ d0

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

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: xdim, ydim, zdim, LO, HI
      use dataio_pub,     only: printinfo
      use fluidindex,     only: flind
      use particle_types, only: pset
#ifdef NBODY
      use particle_pub,   only: ht_integrator
      !use particles_io_hdf5
#endif /* NBODY */

      implicit none

      integer                         :: i, j, k, p
#ifdef NBODY
      integer                         :: n_particles        !< number of particles
      real                            :: e                  !< orbit eccentricity
      logical,save                    :: first_run = .true.
      character(len=2)                :: plane
#endif /* NBODY */
      type(cg_list_element),  pointer :: cgl

      do p = lbound(flind%all_fluids, dim=1), ubound(flind%all_fluids, dim=1)
         cgl => leaves%first
         do while (associated(cgl))
            associate(cg => cgl%cg)
               do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
                  do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
                     do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                        associate( fl => flind%all_fluids(p)%fl )
#ifdef NBODY
                           cg%u(fl%idn,i,j,k) = 1.0e-6
#else /* !NBODY */
                           cg%u(fl%idn,i,j,k) = 1.0
#endif /* !NBODY */
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

#ifndef NBODY
      call pset%add(1.0, [ 0.9700436,  -0.24308753,  0.0], [ 0.466203685,  0.43236573, 0.0])
      call pset%add(1.0, [-0.9700436,   0.24308753,  0.0], [ 0.466203685,  0.43236573, 0.0])
      call pset%add(1.0, [ 0.0,         0.0,         0.0], [-0.932407370, -0.86473146, 0.0])
      call printinfo('To see results type: gnuplot -p -e ''plot "nbody_out.log" u 2:3'' ')

   end subroutine problem_initial_conditions
#endif /* !NBODY */

#ifdef NBODY
      if (ht_integrator) then
         if(first_run) then
            call pset%add(1.0, [ 0.9700436, -0.24308753, 0.0], [ 0.466203685, 0.43236573, 0.0], [0.0, 0.0, 0.0],0.0)
            call pset%add(1.0, [-0.9700436, 0.24308753, 0.0], [ 0.466203685, 0.43236573, 0.0], [0.0, 0.0, 0.0],0.0)
            call pset%add(1.0, [ 0.0, 0.0, 0.0], [-0.932407370, -0.86473146, 0.0], [0.0, 0.0, 0.0], 0.0 )
         endif
         first_run=.false.
      else
         e = 0.0
         n_particles = 1
         plane = 'XY'
         call orbits(n_particles, e, first_run, plane)
         !call relax_time(n_particles, first_run)
         !call read_buildgal
         !call twobodies(n_particles, e, first_run, plane)
      endif

      contains

!< \brief rotate (x,y,z) vector by an angle theta
      function positions(dtheta, pos_init, plane)

         implicit none

         real, dimension(3) :: positions, pos_init
         real               :: dtheta
         character(len=2)   :: plane

         positions = rotate(dtheta, pos_init, plane)

      end function positions

!<
!! \brief compute velocity of particle
!!
!! \details compute velocity of particle with position pos_init and eccentricity e <0,1)
!!
!! \warning it works properly only in XY plane 
!>
      function velocities(pos_init, e)

         use constants,  only: zero, one, dpi
         use dataio_pub, only: die
         use func,       only: operator(.equals.)
         use units,      only: newtong

         implicit none

         real, dimension(3) :: pos_init, velocities
         real               :: a        !< semi-major axis of initial elliptical orbit of particle
         real               :: r        !< lenght of radius vector
         real               :: e
         real               :: mu
         real, parameter    :: M = 10.0
         real               :: lenght  !usunac

         mu = newtong*M

         if( (e < zero) .or. (e >= one) ) then
            call die("[initproblem:velocities] Invalid eccentricity")
         else
            r = sqrt(pos_init(1)**2 + pos_init(2)**2 + pos_init(3)**2)

            if (e.equals.0.0) then
               velocities(2) = sqrt(mu/r)
               write(*,*) "Orbita kolowa"
            else
               a = r/(1.0 + e)
               velocities(2) = sqrt(mu*(2.0/r - 1.0/a))
               write(*,'(A11,F4.2,A3,F5.3,A3,F5.3)') "#Elipsa: e=", e, " a=",a, " b=", a*sqrt(1.0 - e**2)

               lenght = dpi*sqrt((a**3)/mu)  !usunac
               write(*,*) "lenght=", lenght
            endif
         endif
         velocities(1) = 0.0
         velocities(3) = 0.0

      end function velocities

!< 
!! \brief rotate vector over one of the axes by an angle theta
!!
!! \todo add to selection of axis (next variable)
!>
      function rotate (theta, vector, plane)

         implicit none

         real, dimension(3) :: vector, rotate
         real               :: theta
         character(len=2)   :: plane

         select case (plane)
            case('XY', 'YX', 'xy', 'yx')
               rotate(1) = vector(1)*cos(theta) - vector(2)*sin(theta)
               rotate(2) = vector(1)*sin(theta) + vector(2)*cos(theta)
               rotate(3) = vector(3)
            case('XZ', 'ZX', 'xz', 'zx')
               rotate(1) = vector(1)*cos(theta) - vector(3)*sin(theta)
               rotate(2) = vector(2)
               rotate(3) = vector(1)*sin(theta) + vector(3)*cos(theta)
            case('YZ', 'ZY', 'yz', 'zy')
               rotate(1) = vector(1)
               rotate(2) = vector(2)*cos(theta) - vector(3)*sin(theta)
               rotate(3) = vector(2)*sin(theta) + vector(3)*cos(theta)
         end select

      end function rotate

      function change_plane(vector, plane)

         implicit none

         real, dimension(3) :: vector, change_plane
         character(len=2)   :: plane

         select case (plane)
            case('XZ', 'ZX', 'xz', 'zx')
               change_plane(1) = vector(1)
               change_plane(2) = vector(3)
               change_plane(3) = vector(2)
            case('YZ', 'ZY', 'yz', 'zy')
               change_plane(3) = vector(2)
               change_plane(2) = vector(1)
               change_plane(1) = vector(3)
         end select

      end function change_plane

      subroutine twobodies(n_particles, e, first_run, plane)

         use dataio_pub,     only: msg, printinfo
         use particle_types, only: pset

         implicit none

         integer,          intent(in)    :: n_particles
         real,             intent(in)    :: e
         logical,          intent(inout) :: first_run
         character(len=2), intent(in)    :: plane
         real, dimension(3)              :: init_pos_body_one, init_pos_body_two, init_vel_body_one, init_vel_body_two
         real                            :: m1, m2

         write(msg,'(a,i6)') '[initproblem:twobodies] Number of particles: ', n_particles
         call printinfo(msg)

         m1 = 10.0
         init_pos_body_one = 0.0
         init_vel_body_one = 0.0

         m2 = 1.0
         init_pos_body_two = [2.0, 0.0, 0.0]
         init_vel_body_two = vel_2bodies(m1, init_pos_body_one, init_pos_body_two, e)

         !init_vel_body_two = 0.0
         write(*,*) m1, " @ ", init_pos_body_one, ", with ", init_vel_body_one 
         write(*,*) m2, " @ ", init_pos_body_two, ", with ", init_vel_body_two 

         if(first_run) then
            call pset%add(m1, init_pos_body_one, init_vel_body_one, [0.0, 0.0, 0.0], 0.0) !dominujace cialo
            call pset%add(m2, init_pos_body_two, init_vel_body_two, [0.0, 0.0, 0.0], 0.0)
            !call pset%add(10.0, init_pos_body_one, init_vel_body_one, [0.0, 0.0, 0.0], 0.0) !dominujace cialo
            first_run=.false.
            write(*,*) "[2b]: Obliczono pozycje czastek "
         endif

      end subroutine twobodies

      function vel_2bodies(mass, init_pos_body_one, init_pos_body_two, e)

         use constants,  only: zero, one, dpi
         use dataio_pub, only: die
         use func,       only: operator(.equals.)
         use units,      only: newtong

         implicit none

         real, dimension(3) :: init_pos_body_one, init_pos_body_two, vel_2bodies
         real               :: a        !< semi-major axis of initial elliptical orbit of particle
         real               :: r        !< lenght of radius vector
         real               :: e
         real               :: mu
         real               :: mass
         real               :: lenght  !usunac

         mu = newtong * mass
         write(*,*) "NEWTONG=", newtong, "mu=", mu

         if( (e < zero) .or. (e >= one) ) then
            call die("[initproblem:velocities] Invalid eccentricity")
         else
            r = sqrt((init_pos_body_one(1) - init_pos_body_two(1))**2 + &
                     (init_pos_body_one(2) - init_pos_body_two(2))**2 + &
                     (init_pos_body_one(3) - init_pos_body_two(3))**2)

            if (e.equals.0.0) then
               vel_2bodies(2) = sqrt(mu/r)
               write(*,*) "Orbita kolowa"
            else
               a = r/(1.0 + e)
               vel_2bodies(2) = sqrt(mu*(2.0/r - 1.0/a))
               write(*,*) "predkosc 1: ", vel_2bodies(2)
               vel_2bodies(2) = sqrt((mu-mu*e)/r)
               write(*,*) "predkosc 2: ", vel_2bodies(2)

               write(*,'(A11,F4.2,A3,F5.3,A3,F5.3)') "#Elipsa: e=", e, " a=",a, " b=", a*sqrt(1.0 - e**2)
               lenght = dpi*sqrt((a**3)/mu)  !usunac
               write(*,*) "lenght=", lenght
            endif
         endif
         vel_2bodies(1) = 0.0
         vel_2bodies(3) = 0.0

      end function vel_2bodies

      subroutine orbits(n_particles, e, first_run, plane)

         !use constants,        only: dpi
         use dataio_pub,       only: msg, printinfo
         use gravity,          only: sum_potential
         use particle_gravity, only: update_particle_gravpot_and_acc
         use particle_types,   only: pset

         implicit none

         integer,          intent(in)    :: n_particles
         real,             intent(in)    :: e
         logical,          intent(inout) :: first_run
         character(len=2), intent(in)    :: plane
         real, dimension(3)              :: pos_init, vel_init
         !real                            :: dtheta

         write(msg,'(a,i6)') '[initproblem:orbits] Number of particles that will be added: ', n_particles
         call printinfo(msg)

         !dtheta = dpi/n_particles

         pos_init(1) = 2.0
         pos_init(2) = 1.0
         pos_init(3) = 0.0

         !vel_init = velocities(pos_init, e)
         vel_init = [-0.5, 0.0, 0.0]

         if(first_run) then
            !call pset%add(1.0, [ 0.9700436, -0.24308753, 0.0], [ 0.466203685, 0.43236573, 0.0], [0.0, 0.0, 0.0], 0.0)
            !call pset%add(1.0, [-0.9700436, 0.24308753, 0.0], [ 0.466203685, 0.43236573, 0.0], [0.0, 0.0, 0.0], 0.0)
            !call pset%add(1.0, [0.0, 0.0, 0.0], [-0.932407370, -0.86473146, 0.0], [0.0, 0.0, 0.0], 0.0 )
            do i = 1, n_particles, 1
               call pset%add(1.0, pos_init, vel_init, [0.0, 0.0, 0.0], 0.0 ) !orbita eliptyczna
               !call pset%add(1.0, [4.0, 2.0, 0.0],[-0.5, 0.0, 0.0], [0.0, 0.0, 0.0], 0.0)
               !call pset%add(1.0, [3.0, 2.0, 0.0],[0.0, -1.0, 0.0],  [0.0, 0.0, 0.0], 0.0)

               !pos_init = positions(dtheta, pos_init, plane)
               !vel_init = rotate(dtheta, vel_init, plane)
            enddo

            !call pset%add(10.0, [0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],0.0)
            !call pset%add(1.0, [0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],0.0) ! to "dziala"

            first_run = .false.

            write(msg,'(a,i6)') '[initproblem:orbits] Number of particles added to the domain: ', size(pset%p, dim=1)
            call printinfo(msg)
            write(msg,'(a,3f5.2)') '[initproblem:orbits] Initial position of the particle: ', pos_init
            call printinfo(msg)
            write(msg,'(a,3f5.2)') '[initproblem:orbits] Initial velocity of the particle: ', vel_init
            call printinfo(msg)

            call update_particle_gravpot_and_acc
            call sum_potential

         endif
      end subroutine orbits

!< \brief create a set of particles at random positions inside a sphere
      subroutine relax_time(n_particles, first_run)

         use dataio_pub,        only: msg, printinfo
         use particle_types,    only: pset
#ifdef HDF5
         use particles_io_hdf5, only: write_hdf5, read_hdf5
#endif /* HDF5 */

         implicit none

         integer, intent(in)            :: n_particles
         logical, intent(inout)         :: first_run

         integer                        :: i, j
         integer, parameter             :: seed = 86437
         real, dimension(n_particles,3) :: pos_init
         real, dimension(3,2)           :: domain
         real                           :: factor, r_dom
         real, parameter                :: onesixth = 1.0/6.0

         domain(1,1) = -5.0
         domain(2,1) = -5.0
         domain(3,1) = -5.0
         domain(1,2) = 5.0
         domain(2,2) = 5.0
         domain(3,2) = 5.0

         write(msg,'(a,i6)') '[initproblem:relax_time] Number of particles: ', n_particles
         call printinfo(msg)

         call srand(seed)
         r_dom = onesixth*sqrt(domain(1,2)**2 + domain(2,2)**2 + domain(3,2)**2)

         if(first_run) then
            do i = 1, n_particles
               r = r_dom
               do while ((r>=r_dom))
                  do j = 1, 3
                     factor = rand(0)
                     pos_init(i, j) = sign(rand(0)*domain(j, 2),factor-0.5)
                  enddo
                  r = sqrt(pos_init(i,1)**2 + pos_init(i,2)**2 + pos_init(i,3)**2)
               enddo
               call pset%add(1.0, pos_init(i,:), [0.0,0.0,0.0],[0.0, 0.0, 0.0],0.0 )
            enddo
            first_run = .false.
            write(msg,'(a,i6)') '[initproblem:relax_time] Particles position computed'
            call printinfo(msg)
#ifdef HDF5
            call write_hdf5(pos_init, n_particles)
#endif /* HDF5 */
         endif

      end subroutine relax_time

   end subroutine problem_initial_conditions

   subroutine read_buildgal

      use particle_types, only: pset

      implicit none

      integer                           :: i, j, nbodies, dims=3
      integer                           :: galfile = 1
      character(len=6)                  :: galname="SPIRAL"
      real, dimension(:,:), allocatable :: pos, vel
      real, dimension(:),   allocatable :: mass

      open(unit=galfile,file=galname,action="read",status="old")
         read(galfile,*) nbodies
         write(*,*) "nbodies=",nbodies

         allocate(mass(nbodies),pos(nbodies,3),vel(nbodies,3))

         read(galfile,*) (mass(i),i=1,nbodies), &
         ((pos(i,j),j=1,dims),i=1,nbodies),&
         ((vel(i,j),j=1,dims),i=1,nbodies)

      close(galfile)
      
      open(unit=2,file='galtest.dat')
         do i = 1, nbodies
            if (modulo(i, 1000) .eq. 0) then
               write(*,*) i
            endif
            call pset%add(mass(i), pos(i,:), vel(i,:),[0.0, 0.0, 0.0], 0.0)
         enddo
      close(2)

   end subroutine read_buildgal
#endif /* NBODY */
!-----------------------------------------------------------------------------
end module initproblem
