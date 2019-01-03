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
module initproblem

! Initial condition for Sedov-Taylor explosion
! Written by: M. Hanasz, March 2006

   use constants, only: cbuff_len

   implicit none

   private
   public  :: read_problem_par, problem_initial_conditions, problem_pointers

   character(len=cbuff_len) :: topic_2body
   real                     :: fdens               !< fluid density
   real                     :: e                   !< orbit eccentricity
   real                     :: mass1               !< (higher) mass of the primary particle
   real                     :: mass2               !< (lower) mass of secondary particles

   namelist /PROBLEM_CONTROL/ topic_2body, fdens, e, mass1, mass2

contains
!-----------------------------------------------------------------------------
   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers
!-----------------------------------------------------------------------------
   subroutine read_problem_par

      use dataio_pub, only: nh      ! QA_WARN required for diff_nml
      use mpisetup,   only: cbuff, rbuff, master, slave, piernik_MPI_Bcast

      implicit none

      ! namelist default parameter values
      topic_2body = 'default'
      fdens       = 1.0e-6
      e           = 0.0
      mass1       = 10.0
      mass2       = 1.0

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=PROBLEM_CONTROL, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL")
         read(nh%cmdl_nml,nml=PROBLEM_CONTROL, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         call nh%compare_namelist()

         cbuff(1) = topic_2body
         rbuff(1) = fdens
         rbuff(2) = e
         rbuff(3) = mass1
         rbuff(4) = mass2

      endif

      call piernik_MPI_Bcast(cbuff, cbuff_len)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         topic_2body = cbuff(1)
         fdens       = rbuff(1)
         e           = rbuff(2)
         mass1       = rbuff(3)
         mass2       = rbuff(4)

      endif

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine problem_initial_conditions

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use dataio_pub, only: die, msg, printinfo
      use fluidindex, only: flind

      implicit none

      integer                         :: p
      type(cg_list_element),  pointer :: cgl

      cgl => leaves%first
      do while (associated(cgl))
         associate(cg => cgl%cg)
            do p = 1, flind%fluids
               associate(fl => flind%all_fluids(p)%fl)
                  cg%u(fl%idn,:,:,:) = fdens
                  cg%u(fl%imx,:,:,:) = 0.0
                  cg%u(fl%imy,:,:,:) = 0.0
                  cg%u(fl%imz,:,:,:) = 0.0
               end associate
            enddo
         end associate
         cgl => cgl%nxt
      enddo

      select case (trim(topic_2body))
         case ('orbits')
            call orbits
         case ('relaxtime')
            call relax_time
         case ('buildgal')
            call read_buildgal
         case ('twobodies')
            call twobodies
         case default
            write(msg, '(3a)')"[initproblem:problem_initial_conditions] Unknown topic_2body '",trim(topic_2body),"'"
            call die(msg)
      end select

   end subroutine problem_initial_conditions

   subroutine twobodies

      use constants,      only: ndims
      use particle_types, only: pset

      implicit none

      real, dimension(ndims,2) :: init_pos_body, init_vel_body
      real, dimension(2)       :: m
      integer                  :: p

      m = [mass1, mass2]
      init_pos_body = 0.0
      init_vel_body = 0.0

      init_pos_body(:,2) = [2.0, 0.0, 0.0]
      init_vel_body(:,2) = vel_2bodies(m(1), init_pos_body)

      do p = 1, 2
         write(*,*) m(p), " @ ", init_pos_body(:,p), ", with ", init_vel_body(:,p)
         call pset%add(m(p), init_pos_body(:,p), init_vel_body(:,p), [0.0, 0.0, 0.0], 0.0)
      enddo

   end subroutine twobodies

   function vel_2bodies(mass, init_pos_body)

      use constants,  only: ndims, zero, one, two, ydim
      use dataio_pub, only: die
      use func,       only: operator(.equals.)
      use units,      only: newtong

      implicit none

      real,                     intent(in) :: mass
      real, dimension(ndims,2), intent(in) :: init_pos_body
      real, dimension(ndims)               :: vel_2bodies
      real                                 :: a        !< semi-major axis of initial elliptical orbit of particle
      real                                 :: r        !< lenght of radius vector
      real                                 :: mu

      vel_2bodies = zero
      mu = newtong * mass

      if( (e < zero) .or. (e >= one) ) call die("[initproblem:vel_2bodies] Invalid eccentricity")

      r = sqrt(sum((init_pos_body(:,1) - init_pos_body(:,2))**2))

      if (e .equals. zero) then
         vel_2bodies(2) = sqrt(mu/r)
      else
         a = r/(one + e)
         vel_2bodies(ydim) = sqrt(mu*(two/r - one/a))
         write(*,*) "velocity 1: ", vel_2bodies(2)
         vel_2bodies(ydim) = sqrt((mu-mu*e)/r)
         write(*,*) "velocity 2: ", vel_2bodies(2)

         write(*,'(A11,F4.2,A3,F5.3,A3,F5.3)') "#Elipsa: e=", e, " a=",a, " b=", a*sqrt(one - e**2)
      endif

   end function vel_2bodies

   subroutine orbits

      use constants,        only: ndims !, dpi, zdim
      use dataio_pub,       only: msg, printinfo
      use gravity,          only: sum_potential
      use particle_gravity, only: update_particle_gravpot_and_acc
      use particle_pub,     only: npart
      use particle_types,   only: pset

      implicit none

      integer                :: p
      real, dimension(ndims) :: pos_init, vel_init

      pos_init = [2.0, 1.0, 0.0]

      !vel_init = velocities(pos_init)
      vel_init = [-0.5, 0.0, 0.0]

      do p = 1, npart, 1
         call pset%add(mass2, pos_init, vel_init, [0.0, 0.0, 0.0], 0.0 ) !elliptical orbit
         !call pset%add(mass2, [4.0, 2.0, 0.0],[-0.5, 0.0, 0.0], [0.0, 0.0, 0.0], 0.0)
         !call pset%add(mass2, [3.0, 2.0, 0.0],[0.0, -1.0, 0.0],  [0.0, 0.0, 0.0], 0.0)

         !pos_init = rotate(pos_init, dpi/npart, zdim)
         !vel_init = rotate(vel_init, dpi/npart, zdim)
      enddo

      !call pset%add(mass1, [0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],0.0)
      !call pset%add(mass2, [0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],0.0) ! it "works"

      write(msg,'(a,i6)') '[initproblem:orbits] Number of particles added to the domain: ', size(pset%p, dim=1)
      call printinfo(msg)
      write(msg,'(a,3f5.2)') '[initproblem:orbits] Initial position of the particle: ', pos_init
      call printinfo(msg)
      write(msg,'(a,3f5.2)') '[initproblem:orbits] Initial velocity of the particle: ', vel_init
      call printinfo(msg)

      call update_particle_gravpot_and_acc
      call sum_potential

   end subroutine orbits

!<
!! \brief compute velocity of particle
!! \details compute velocity of particle with position pos_init and eccentricity e <0,1)
!! \warning it works properly only in XY plane
!>
   function velocities(pos_init)

      use constants,  only: dpi, one, two, ydim, zero
      use dataio_pub, only: die
      use func,       only: operator(.equals.)
      use units,      only: newtong

      implicit none

      real, dimension(3) :: pos_init, velocities
      real               :: a        !< semi-major axis of initial elliptical orbit of particle
      real               :: r        !< lenght of radius vector
      real               :: mu
      real, parameter    :: M = 10.0
      real               :: lenght  !usunac

      mu = newtong*M
      velocities(:) = zero

      if( (e < zero) .or. (e >= one) ) then
         call die("[initproblem:velocities] Invalid eccentricity")
      else
         r = sqrt(sum(pos_init(:)**2))

         if (e .equals. zero) then
            velocities(ydim) = sqrt(mu/r)
            write(*,*) "Circular orbit"
         else
            a = r/(one + e)
            velocities(ydim) = sqrt(mu*(two/r - one/a))
            write(*,'(A11,F4.2,A3,F5.3,A3,F5.3)') "#Ellipse: e=", e, " a=",a, " b=", a*sqrt(one - e**2)

            lenght = dpi*sqrt((a**3)/mu)  !usunac
            write(*,*) "lenght=", lenght
         endif
      endif

   end function velocities

!<
!! \brief rotate vector over one of the axes by an angle theta
!! \todo add to selection of axis (next variable)
!>
   function rotate (vector, theta, dir)

      use constants, only: ndims, xdim, ydim

      implicit none

      real, dimension(ndims), intent(in) :: vector
      real,                   intent(in) :: theta
      integer(kind=4),        intent(in) :: dir
      real, dimension(ndims)             :: rotate
      integer(kind=4)                    :: dir1, dir2

      rotate(dir) = vector(dir)
      dir1 = mod(dir+xdim,ndims)
      dir2 = mod(dir+ydim,ndims)
      rotate(dir1) = vector(dir1)*cos(theta) - vector(dir2)*sin(theta)
      rotate(dir2) = vector(dir1)*sin(theta) + vector(dir2)*cos(theta)

   end function rotate

!< \brief create a set of particles at random positions inside a sphere
   subroutine relax_time

      use constants,         only: ndims, onesth, one, two
      use domain,            only: dom
      use particle_pub,      only: npart
      use particle_types,    only: pset
#ifdef VERBOSE
      use dataio_pub,        only: printinfo
#endif /* VERBOSE */
#ifdef HDF5
      use particles_io_hdf5, only: write_hdf5
#endif /* HDF5 */

      implicit none

      integer                      :: i
      real, dimension(npart,ndims) :: pos_init
      real, dimension(ndims)       :: nrand
      real                         :: r_dom
      logical                      :: outsphere
#ifndef RANDOMIZE
      integer                      :: seed = 86437

      call random_seed(seed)
#endif /* !RANDOMIZE */

      r_dom = onesth*sqrt(sum(dom%L_(:)**2))

      do i = 1, npart
         outsphere = .true.
         do while (outsphere)
            call random_number(nrand)
            pos_init(i,:) = dom%C_ + (two*nrand-one)*r_dom
            outsphere = sqrt(sum(pos_init(i,:)**2)) >= r_dom
         enddo
         call pset%add(mass2, pos_init(i,:), [0.0,0.0,0.0], [0.0, 0.0, 0.0], 0.0)
      enddo
#ifdef VERBOSE
      call printinfo('[initproblem:relax_time] Particles position computed')
#endif /* VERBOSE */
#ifdef HDF5
      call write_hdf5(pos_init, npart)
#endif /* HDF5 */

   end subroutine relax_time

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

end module initproblem
