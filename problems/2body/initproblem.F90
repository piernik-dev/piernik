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

   namelist /PROBLEM_CONTROL/ topic_2body

contains
!-----------------------------------------------------------------------------
   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers
!-----------------------------------------------------------------------------
   subroutine read_problem_par

      use dataio_pub, only: nh      ! QA_WARN required for diff_nml
      use mpisetup,   only: cbuff, master, slave, piernik_MPI_Bcast

      implicit none

      ! namelist default parameter values
      topic_2body = 'default'

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

      endif

      call piernik_MPI_Bcast(cbuff, cbuff_len)

      if (slave) then

         topic_2body = cbuff(1)

      endif

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
      use dataio_pub,     only: die, msg
      use particle_pub,   only: ht_integrator
      !use particles_io_hdf5
#endif /* NBODY */

      implicit none

      integer                         :: i, j, k, p
#ifdef NBODY
      real                            :: e                  !< orbit eccentricity
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
#endif /* !NBODY */

#ifdef NBODY
      if (ht_integrator) then
         call pset%add(1.0, [ 0.9700436,  -0.24308753,  0.0], [ 0.466203685,  0.43236573, 0.0], [0.0, 0.0, 0.0], 0.0)
         call pset%add(1.0, [-0.9700436,   0.24308753,  0.0], [ 0.466203685,  0.43236573, 0.0], [0.0, 0.0, 0.0], 0.0)
         call pset%add(1.0, [ 0.0,         0.0,         0.0], [-0.932407370, -0.86473146, 0.0], [0.0, 0.0, 0.0], 0.0)
      else
         e = 0.0
         plane = 'XY'

         select case (trim(topic_2body))
            case ('orbits')
               call orbits(e, plane)
            case ('relaxtime')
               call relax_time
            case ('buildgal')
               call read_buildgal
            case ('twobodies')
               call twobodies(e)
            case default
               write(msg, '(3a)')"[initproblem:problem_initial_conditions] Unknown topic_2body '",trim(topic_2body),"'"
               call die(msg)
         end select
      endif
#endif /* NBODY */

   end subroutine problem_initial_conditions

#ifdef NBODY
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
!! \details compute velocity of particle with position pos_init and eccentricity e <0,1)
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

   subroutine twobodies(e)

      use particle_types, only: pset

      implicit none

      real,             intent(in)    :: e
      real, dimension(3)              :: init_pos_body_one, init_pos_body_two, init_vel_body_one, init_vel_body_two
      real                            :: m1, m2

      m1 = 10.0
      init_pos_body_one = 0.0
      init_vel_body_one = 0.0

      m2 = 1.0
      init_pos_body_two = [2.0, 0.0, 0.0]
      init_vel_body_two = vel_2bodies(m1, init_pos_body_one, init_pos_body_two, e)

      !init_vel_body_two = 0.0
      write(*,*) m1, " @ ", init_pos_body_one, ", with ", init_vel_body_one
      write(*,*) m2, " @ ", init_pos_body_two, ", with ", init_vel_body_two

      call pset%add(m1, init_pos_body_one, init_vel_body_one, [0.0, 0.0, 0.0], 0.0) !dominujace cialo
      call pset%add(m2, init_pos_body_two, init_vel_body_two, [0.0, 0.0, 0.0], 0.0)
      !call pset%add(10.0, init_pos_body_one, init_vel_body_one, [0.0, 0.0, 0.0], 0.0) !dominujace cialo
      write(*,*) "[2b]: Obliczono pozycje czastek "

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

   subroutine orbits(e, plane)

      use constants,        only: ndims !, dpi
      use dataio_pub,       only: msg, printinfo
      use gravity,          only: sum_potential
      use particle_gravity, only: update_particle_gravpot_and_acc
      use particle_pub,     only: npart
      use particle_types,   only: pset

      implicit none

      real,             intent(in) :: e
      character(len=2), intent(in) :: plane
      real, dimension(ndims)       :: pos_init, vel_init
      !real                         :: dtheta
      integer                      :: p

      write(msg,'(a,i6)') '[initproblem:orbits] Number of particles that will be added: ', npart
      call printinfo(msg)

      !dtheta = dpi/npart

      pos_init = [2.0, 1.0, 0.0]

      !vel_init = velocities(pos_init, e)
      vel_init = [-0.5, 0.0, 0.0]

      !call pset%add(1.0, [ 0.9700436, -0.24308753, 0.0], [ 0.466203685, 0.43236573, 0.0], [0.0, 0.0, 0.0], 0.0)
      !call pset%add(1.0, [-0.9700436, 0.24308753, 0.0], [ 0.466203685, 0.43236573, 0.0], [0.0, 0.0, 0.0], 0.0)
      !call pset%add(1.0, [0.0, 0.0, 0.0], [-0.932407370, -0.86473146, 0.0], [0.0, 0.0, 0.0], 0.0 )
      do p = 1, npart, 1
         call pset%add(1.0, pos_init, vel_init, [0.0, 0.0, 0.0], 0.0 ) !orbita eliptyczna
         !call pset%add(1.0, [4.0, 2.0, 0.0],[-0.5, 0.0, 0.0], [0.0, 0.0, 0.0], 0.0)
         !call pset%add(1.0, [3.0, 2.0, 0.0],[0.0, -1.0, 0.0],  [0.0, 0.0, 0.0], 0.0)

         !pos_init = positions(dtheta, pos_init, plane)
         !vel_init = rotate(dtheta, vel_init, plane)
      enddo

      !call pset%add(10.0, [0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],0.0)
      !call pset%add(1.0, [0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],0.0) ! to "dziala"

      write(msg,'(a,i6)') '[initproblem:orbits] Number of particles added to the domain: ', size(pset%p, dim=1)
      call printinfo(msg)
      write(msg,'(a,3f5.2)') '[initproblem:orbits] Initial position of the particle: ', pos_init
      call printinfo(msg)
      write(msg,'(a,3f5.2)') '[initproblem:orbits] Initial velocity of the particle: ', vel_init
      call printinfo(msg)

      call update_particle_gravpot_and_acc
      call sum_potential

   end subroutine orbits

!< \brief create a set of particles at random positions inside a sphere
   subroutine relax_time

      use constants,         only: onesth
      use dataio_pub,        only: msg, printinfo
      use particle_pub,      only: npart
      use particle_types,    only: pset
#ifdef HDF5
      use particles_io_hdf5, only: write_hdf5, read_hdf5
#endif /* HDF5 */

      implicit none

      integer                  :: i, j
      integer, parameter       :: seed = 86437
      real, dimension(npart,3) :: pos_init
      real, dimension(3,2)     :: domain
      real                     :: factor, r_dom, r

      domain(:,1) = -5.0
      domain(:,2) = 5.0

      write(msg,'(a,i6)') '[initproblem:relax_time] Number of particles: ', npart
      call printinfo(msg)

      call srand(seed)
      r_dom = onesth*sqrt(domain(1,2)**2 + domain(2,2)**2 + domain(3,2)**2)

      do i = 1, npart
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
      write(msg,'(a,i6)') '[initproblem:relax_time] Particles position computed'
      call printinfo(msg)
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
#endif /* NBODY */

end module initproblem
