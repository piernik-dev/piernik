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

   use constants, only: cbuff_len, cwdlen

   implicit none

   private
   public  :: read_problem_par, problem_initial_conditions, problem_initial_nbody, problem_pointers

   character(len=cbuff_len) :: topic_2body
   character(len=cwdlen)    :: bgfile              !< buildgal file name
   real                     :: fdens               !< fluid density
   real                     :: e                   !< orbit eccentricity
   real                     :: mass1               !< (higher) mass of the primary particle
   real                     :: mass2               !< (lower) mass of secondary particles

   namelist /PROBLEM_CONTROL/ topic_2body, fdens, e, mass1, mass2, bgfile

contains
!-----------------------------------------------------------------------------
   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers
!-----------------------------------------------------------------------------
   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use dataio_pub, only: nh, restarted_sim
      use mpisetup,   only: cbuff, rbuff, master, slave

      implicit none

      ! namelist default parameter values
      topic_2body = 'default'
      bgfile      = 'SPIRAL'
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
      call piernik_MPI_Bcast(bgfile, cwdlen)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         topic_2body = cbuff(1)

         fdens       = rbuff(1)
         e           = rbuff(2)
         mass1       = rbuff(3)
         mass2       = rbuff(4)

      endif

      ! Perhaps a more general approach should be implemeted once we decide to drop most of precompiler conditionals on NBODY.
      ! Right now it is used in just one setup.
      ! restarted_sim is set up early in the initpiernik, in init_dataio_parameters, much earlier than init_dataio
      if (.not. restarted_sim) call problem_initial_nbody

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine problem_initial_conditions

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use fluidindex,       only: flind

      implicit none

      integer                        :: p
      type(cg_list_element), pointer :: cgl

      cgl => leaves%first
      do while (associated(cgl))
         associate(cg => cgl%cg)
            do p = 1, flind%fluids
               associate(fl => flind%all_fluids(p)%fl)
                  cg%u(fl%idn,RNG) = fdens
                  cg%u(fl%imx,RNG) = 0.0
                  cg%u(fl%imy,RNG) = 0.0
                  cg%u(fl%imz,RNG) = 0.0
               end associate
            enddo
         end associate
         cgl => cgl%nxt
      enddo

   end subroutine problem_initial_conditions

   subroutine problem_initial_nbody

      use dataio_pub, only: die, msg

      implicit none

      select case (trim(topic_2body))
         case ('twobodies')
            call twobodies
         case ('buildgal')
            call read_buildgal
         case default
            write(msg, '(3a)')"[initproblem:problem_initial_conditions] Unknown topic_2body '",trim(topic_2body),"'"
            call die(msg)
      end select

   end subroutine problem_initial_nbody

   subroutine twobodies

      use constants,      only: ndims
      use dataio_pub,     only: msg, printinfo
      use particle_utils, only: add_part_in_proper_cg

      implicit none

      real, dimension(ndims,2) :: init_pos_body, init_vel_body
      real, dimension(2)       :: m
      integer(kind=4)          :: p

      m = [mass1, mass2]
      init_pos_body = 0.0
      init_vel_body = 0.0

      init_pos_body(:,2) = [2.0, 0.0, 0.0]
      init_vel_body(:,2) = vel_2bodies(m(1), init_pos_body(:,1)-init_pos_body(:,2))

      do p = 1, 2
         write(msg,'(f8.5,a,3f8.5,a,3f8.5)') m(p), " @ ", init_pos_body(:,p), ", with ", init_vel_body(:,p) ; call printinfo(msg)
         call add_part_in_proper_cg(p, m(p), init_pos_body(:,p), init_vel_body(:,p), [0.0, 0.0, 0.0], 0.0)
      enddo

   end subroutine twobodies

!<
!! \brief compute velocity of particle
!! \details compute velocity of particle with position pos_init and eccentricity e <0,1)
!! \warning it works properly only in XY plane
!>
   function vel_2bodies(mass, rel_pos)

      use constants,  only: ndims, one, ydim, zero
      use dataio_pub, only: die, msg, printinfo
      use func,       only: operator(.equals.)
      use units,      only: newtong

      implicit none

      real,                   intent(in) :: mass
      real, dimension(ndims), intent(in) :: rel_pos
      real, dimension(ndims)             :: vel_2bodies
      real                               :: a        !< semi-major axis of initial elliptical orbit of particle
      real                               :: r        !< length of radius vector
      real                               :: mu

      vel_2bodies = zero
      mu = newtong * mass

      if ( (e < zero) .or. (e >= one) ) call die("[initproblem:vel_2bodies] Invalid eccentricity")

      r = sqrt(sum(rel_pos**2))

      if (e .equals. zero) then
         vel_2bodies(ydim) = sqrt(mu/r)
         call printinfo("[initproblem:vel_2bodies] Circular orbit")
      else
         a = r/(one + e)
         vel_2bodies(ydim) = sqrt((one-e)*mu/r)

         write(msg,'(A11,F4.2,A3,F5.3,A3,F5.3)') "#Ellipse: e=", e, " a=",a, " b=", a*sqrt(one - e**2) ; call printinfo(msg)
      endif

   end function vel_2bodies

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

! the routine read_buildgal should perhaps go somewhere to particles or IO (together with other readers for different particle containers)

   subroutine read_buildgal

      use constants,      only: ndims, I_ONE
      use dataio_pub,     only: msg, printio
      use particle_utils, only: add_part_in_proper_cg
      use mpisetup,       only: master

      implicit none

      integer                           :: j
      integer(kind=4)                   :: i, nbodies
      integer, parameter                :: galfile = 1
      real, dimension(:,:), allocatable :: pos, vel
      real, dimension(:),   allocatable :: mass
      real, dimension(ndims)            :: posi, veli

      open(unit=galfile, file=bgfile, action='read', status='old')
      read(galfile,*) nbodies
         if (master) then
            write(msg,'(3a,i8,a)') 'Reading ', trim(bgfile), ' file with ', nbodies, ' particles'
            call printio(msg)
         endif

         allocate(mass(nbodies), pos(nbodies,ndims), vel(nbodies,ndims))

         read(galfile,*) (mass(i), i = 1, nbodies), ((pos(i,j), j = 1, ndims), i = 1, nbodies), ((vel(i,j), j = 1, ndims), i = 1, nbodies)

      close(galfile)

      i = 0
      do j = 1, nbodies
         i = i + I_ONE
         if (i > nbodies) exit
#ifdef VERBOSE
         if (modulo(i, 10000) == 0) then
            write(msg,'(i8,a)') i, ' particles read' ; call printio(msg)
         endif
#endif /* VERBOSE */
         posi = pos(i,:)
         veli = vel(i,:)
         call add_part_in_proper_cg(i, mass(i), posi, veli, [0.0, 0.0, 0.0], 0.0)
      enddo
      ! ToDo: check whether all particles were added exactly once
      deallocate(mass,pos,vel)

   end subroutine read_buildgal

end module initproblem
