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

! A blob made of particles thrown into space near a point-mass.
! This setup is meant to stress-test and give hints on optimization of particle implementation in Piernik.

   use constants, only: ndims, INT4

   implicit none

   private
   public  :: read_problem_par, problem_initial_conditions, problem_pointers

   real                   :: blob_mass       !< Total mass of all particles created initially
   integer(kind=INT4)     :: blob_particles  !< Number of particles that the blob_mass is evenly split to
   real, dimension(ndims) :: blob_center     !< The center of the swarm of particles
   real                   :: blob_size       !< The size of the swarm of particles
   real, dimension(ndims) :: blob_velocity   !< The velocity of the particles

   ! This is an early implementation, so no fancy shaping and customizing for keplerian velocity, etc. yet.

   namelist /PROBLEM_CONTROL/ blob_mass, blob_particles, blob_center, blob_size, blob_velocity

contains

   !< \brief No special problem pointers to be set here

   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

   !< \brief Read the parameters

   subroutine read_problem_par

      use dataio_pub, only: nh, restarted_sim, warn
      use mpisetup,   only: ibuff, rbuff, master, slave, piernik_MPI_Bcast

      implicit none

      ! namelist default parameter values
      blob_mass      = 1e-20
      blob_particles = 1
      blob_center    = [ 1., 0., 0. ]
      blob_size      = 0.1
      blob_velocity  = [ 0., 1., 0. ]

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

         ibuff(1) = blob_particles

         rbuff(1) = blob_mass
         rbuff(2) = blob_size
         rbuff(3:3+ndims-1) = blob_center
         rbuff(3+ndims:3+2*ndims-1) = blob_velocity

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         blob_particles = ibuff(1)

         blob_mass      = rbuff(1)
         blob_size      = rbuff(2)
         blob_center    = rbuff(3:3+ndims-1)
         blob_velocity  = rbuff(3+ndims:3+2*ndims-1)

      endif

      if (blob_particles < 0) then
         blob_particles = 0
         call warn("[initproblem:read_problem_par] blob_particles has to be >= 0 (fixed)")
      endif

      ! Perhaps a more general approach should be implemeted once we decide to drop most of precompiler conditionals on NBODY.
      ! Right now it is used in just one setup.
      ! restarted_sim is set up early in the initpiernik, in init_dataio_parameters, much earlier than init_dataio
      if (.not. restarted_sim) call problem_create_nbodies

   end subroutine read_problem_par

   !> \brief Set up the ambient gas.

   subroutine problem_initial_conditions

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use fluidindex, only: flind
      use func,       only: ekin
      use global,     only: smalld, smallei

      implicit none

      integer                        :: p
      type(cg_list_element), pointer :: cgl

      cgl => leaves%first
      do while (associated(cgl))
         associate(cg => cgl%cg)
            do p = 1, flind%fluids
               associate(fl => flind%all_fluids(p)%fl)
                  cg%u(fl%idn, RNG) = smalld
                  cg%u(fl%imx, RNG) = 0.0
                  cg%u(fl%imy, RNG) = 0.0
                  cg%u(fl%imz, RNG) = 0.0
                  cg%u(fl%ien, RNG) = smallei + ekin(cg%u(fl%imx,RNG), cg%u(fl%imy,RNG), cg%u(fl%imz,RNG), cg%u(fl%idn,RNG))
                  ! Tricky: smallei will control the initial temperature and thus may limit the timestep
               end associate
            enddo
         end associate
         cgl => cgl%nxt
      enddo

   end subroutine problem_initial_conditions

   !> \brief Set up the particles

   subroutine problem_create_nbodies

      use constants,      only: INT4, ndims, half, I_ONE
      use particle_utils, only: add_part_in_proper_cg

      implicit none

      real, dimension(ndims) :: pos
      real :: d
      integer(kind=INT4) :: p

      d = real(blob_particles)**(1./3.)

      do p = I_ONE, blob_particles
         ! The prescription below gives a particle exactly at blob_center for blob_particles == 1, 3**3, 5**3, ...
         pos = blob_center + blob_size * ([        (p - half)        / blob_particles, &
              &                             modulo((p - half) * d    / blob_particles, 1.), &
              &                             modulo((p - half) * d**2 / blob_particles, 1.) ] - half)
         call add_part_in_proper_cg(p, blob_mass/blob_particles, pos, blob_velocity, [0.0, 0.0, 0.0], 0.0)
      enddo

   end subroutine problem_create_nbodies

end module initproblem
