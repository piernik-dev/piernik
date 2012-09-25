! $Id$
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
#include "macros.h"

!>
!! \brief Main particle module
!!
!! This module contains all particle related things that are needed by the rest
!! of the code. No other module should directly import anything from particle_*
!<

module particle_pub
! pulled by GRAV
   use particle_types,  only: particle_set, particle_solver_T

   implicit none
   private
   public :: pset, psolver, init_particles, cleanup_particles

   type(particle_set) :: pset !< default particle list
   class(particle_solver_T), pointer :: psolver

contains

!> \brief Read namelist and initialize pset with 0 particles

   subroutine init_particles

      use constants,             only: cbuff_len
      use dataio_pub,            only: nh  ! QA_WARN required for diff_nml
      use dataio_pub,            only: msg, die
      use mpisetup,              only: master, slave, cbuff, piernik_mpi_bcast
      use particle_integrators,  only: hermit4

      implicit none
      character(len=cbuff_len) :: time_integrator
      character(len=cbuff_len), parameter :: default_ti = "none"

      namelist /PARTICLES/ time_integrator

      time_integrator = default_ti

      if (master) then
         diff_nml(PARTICLES)

         cbuff(1) = time_integrator
      endif

      call piernik_MPI_Bcast(cbuff, cbuff_len)

      if (slave) then
         time_integrator = cbuff(1)
      endif

      psolver => null()
      select case (trim(time_integrator))
         case ('hermit4')
            psolver => hermit4
         case (default_ti) ! be quiet
         case default
            write(msg, '(3a)')"[particle_pub:init_particles] Unknown integrator '",trim(time_integrator),"'"
            call die(msg)
      end select

      call pset%init()

   end subroutine init_particles

!> \brief Deallocate pset

   subroutine cleanup_particles

     implicit none

     call pset%cleanup

   end subroutine cleanup_particles

end module particle_pub
