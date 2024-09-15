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

!>
!! \brief Main particle module
!!
!! This module contains all particle related things that are needed by the rest
!! of the code. No other module should directly import anything from particle_*
!<

module particle_pub
! pulled by NBODY

   use constants, only: cbuff_len

   implicit none

   private
   public :: init_particles, cleanup_particles
   public :: eps, lf_c, ignore_dt_fluid, is_setacc_cic, is_setacc_int, is_setacc_tsc, mask_gpot1b, dump_diagnose, r_soft, twodtscheme, default_ti, time_integrator

   real    :: lf_c               !< timestep should depends of grid and velocities of particles (used to extrapolation of the gravitational potential)
   logical :: ignore_dt_fluid    !< option to test only nbody part of the code, never true by default
   logical :: dump_diagnose        !< dump diagnose for each particle to a seperate log file
   logical :: twodtscheme
   logical :: is_setacc_cic, is_setacc_int, is_setacc_tsc, mask_gpot1b
   real    :: eps, r_soft
   character(len=cbuff_len)            :: time_integrator
   character(len=cbuff_len), parameter :: default_ti = "none"

contains

!> \brief Read namelist of particle parameters
   subroutine init_particles

      use bcast,         only: piernik_MPI_Bcast
      use constants,     only: I_NGP, I_CIC, I_TSC, V_INFO
      use dataio_pub,    only: msg, die, printinfo, nh
      use mpisetup,      only: master, slave, cbuff, ibuff, lbuff, rbuff
      use particle_func, only: check_ord
      use particle_maps, only: set_map

      implicit none

      character(len=cbuff_len)            :: interpolation_scheme
      character(len=cbuff_len), parameter :: default_is = "ngp"

      character(len=cbuff_len) :: acc_interp_method  !< acceleration interpolation method
      integer                  :: order              !< order of Lagrange polynomials (if acc_interp_method = 'lagrange')

      namelist /PARTICLES/ time_integrator, interpolation_scheme, acc_interp_method, lf_c, eps, mask_gpot1b, ignore_dt_fluid, dump_diagnose

      time_integrator      = default_ti
      interpolation_scheme = default_is

      acc_interp_method    = 'cic'
      lf_c                 = 1.0
      eps                  = 0.1
      r_soft               = 0.0
      twodtscheme          = .false.
      ignore_dt_fluid      = .false.
      dump_diagnose        = .false.
#ifdef NBODY_GRIDDIRECT
      mask_gpot1b          = .true.
#else /* !NBODY_GRIDDIRECT */
      mask_gpot1b          = .false.
#endif /* !NBODY_GRIDDIRECT */

      if (master) then
         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=PARTICLES)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=PARTICLES, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PARTICLES")
         read(nh%cmdl_nml,nml=PARTICLES, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PARTICLES", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=PARTICLES)
         close(nh%lun)
         call nh%compare_namelist()

         cbuff(1) = time_integrator
         cbuff(2) = interpolation_scheme
         cbuff(3) = acc_interp_method

         rbuff(1) = lf_c
         rbuff(2) = eps
         rbuff(3) = r_soft

         lbuff(1) = mask_gpot1b
         lbuff(2) = twodtscheme
         lbuff(3) = ignore_dt_fluid
         lbuff(4) = dump_diagnose
      endif

      call piernik_MPI_Bcast(cbuff, cbuff_len)
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then
         time_integrator      = cbuff(1)
         interpolation_scheme = cbuff(2)
         acc_interp_method    = cbuff(3)

         lf_c                 = rbuff(1)
         eps                  = rbuff(2)
         r_soft               = rbuff(3)

         mask_gpot1b          = lbuff(1)
         twodtscheme          = lbuff(2)
         ignore_dt_fluid      = lbuff(3)
         dump_diagnose        = lbuff(4)
      endif

      select case (trim(interpolation_scheme))
         case ('ngp', 'NGP', 'neareast grid point')
            call set_map(I_NGP)
         case ('cic', 'CIC', 'cloud in cell')
            call set_map(I_CIC)
         case ('tsc', 'TSC', 'triangular shaped cloud')
            call set_map(I_TSC)
         case default
            write(msg, '(3a)')"[particle_pub:init_particles] Unknown interpolation scheme '",trim(interpolation_scheme),"'"
            call die(msg)
      end select

      is_setacc_int = .false.
      is_setacc_cic = .false.
      is_setacc_tsc = .false.
      select case (acc_interp_method)
         case ('lagrange', 'Lagrange', 'polynomials')
            is_setacc_int = .true.
            order = 4
            call check_ord(order)
            msg = "[particle_pub:init_particles] Acceleration interpolation method: Lagrange polynomials"
         case ('cic', 'CIC')
            is_setacc_cic = .true.
            msg = "[particle_pub:init_particles] Acceleration interpolation method: CIC"
         case ('tsc', 'TSC')
            is_setacc_tsc = .true.
            msg = "[particle_pub:init_particles] Acceleration interpolation method: TSC"
      end select
      if (master) call printinfo(trim(msg), V_INFO)

   end subroutine init_particles

!> \brief Deallocate pset

   subroutine cleanup_particles

      implicit none

      ! no need to do anything here

   end subroutine cleanup_particles

end module particle_pub
