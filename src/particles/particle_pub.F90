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
! pulled by GRAV
   use particle_types, only: particle_solver_T, pset

   implicit none
   private
   public :: psolver, init_particles, cleanup_particles
#ifdef NBODY
   public :: npart, lf_c, ignore_dt_fluid

   integer                           :: npart              !< number of particles
   real                              :: lf_c               !< timestep should depends of grid and velocities of particles (used to extrapolation of the gravitational potential)
   logical                           :: ignore_dt_fluid    !< option to test only nbody part of the code, never true by default
#endif /* NBODY */
   class(particle_solver_T), pointer :: psolver

contains

!> \brief Read namelist and initialize pset with 0 particles

   subroutine init_particles

      use constants,             only: cbuff_len, I_NGP, I_CIC, I_TSC
      use dataio_pub,            only: nh  ! QA_WARN required for diff_nml
      use dataio_pub,            only: msg, die
      use mpisetup,              only: master, slave, cbuff, piernik_mpi_bcast
      use particle_integrators,  only: hermit4
#ifdef NBODY
      use dataio_pub,            only: printinfo
      use mpisetup,              only: ibuff, lbuff, rbuff
      use particle_func,         only: check_ord
      use particle_gravity,      only: is_setacc_cic, is_setacc_int, mask_gpot1b, is_setacc_int, is_setacc_tsc
      use particle_integrators,  only: leapfrog2
      use particle_types,        only: twodtscheme, dump_diagnose
#endif /* NBODY */

      implicit none
      character(len=cbuff_len) :: time_integrator
      character(len=cbuff_len) :: interpolation_scheme
      character(len=cbuff_len), parameter :: default_ti = "none"
      character(len=cbuff_len), parameter :: default_is = "ngp"

#ifdef NBODY
      character(len=cbuff_len) :: acc_interp_method  !< acceleration interpolation method
      integer                  :: order              !< order of Lagrange polynomials (if acc_interp_method = 'lagrange')

      namelist /PARTICLES/ npart, time_integrator, interpolation_scheme, acc_interp_method, lf_c, mask_gpot1b, ignore_dt_fluid, dump_diagnose
#else /* !NBODY */
      namelist /PARTICLES/ time_integrator, interpolation_scheme
#endif /* !NBODY */

      time_integrator      = default_ti
      interpolation_scheme = default_is
#ifdef NBODY
      npart                = 0
      acc_interp_method    = 'cic'
      lf_c                 = 1.0
      twodtscheme          = .false.
      ignore_dt_fluid      = .false.
      dump_diagnose        = .false.
#ifdef NBODY_GRIDDIRECT
      mask_gpot1b          = .true.
#else /* !NBODY_GRIDDIRECT */
      mask_gpot1b          = .false.
#endif /* !NBODY_GRIDDIRECT */
#endif /* NBODY */

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
#ifdef NBODY
         cbuff(3) = acc_interp_method
         ibuff(1) = npart
         rbuff(1) = lf_c
         lbuff(1) = mask_gpot1b
         lbuff(2) = twodtscheme
         lbuff(3) = ignore_dt_fluid
         lbuff(4) = dump_diagnose
#endif /* NBODY */
      endif

      call piernik_MPI_Bcast(cbuff, cbuff_len)
#ifdef NBODY
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(rbuff)
#endif /* NBODY */

      if (slave) then
         time_integrator = cbuff(1)
         interpolation_scheme = cbuff(2)
#ifdef NBODY
         acc_interp_method    = cbuff(3)
         npart                = ibuff(1)
         lf_c                 = rbuff(1)
         mask_gpot1b          = lbuff(1)
         twodtscheme          = lbuff(2)
         ignore_dt_fluid      = lbuff(3)
         dump_diagnose        = lbuff(4)
#endif /* NBODY */
      endif

      psolver => null()
      select case (trim(time_integrator))
#ifndef _CRAYFTN
         case ('hermit4')
            psolver => hermit4
#ifdef NBODY
         case ('leapfrog2')
            psolver => leapfrog2
#endif /* NBODY */
         case (default_ti) ! be quiet
#endif /* !_CRAYFTN */
         case default
            write(msg, '(3a)')"[particle_pub:init_particles] Unknown integrator '",trim(time_integrator),"'"
            call die(msg)
      end select

      call pset%init()

      select case (trim(interpolation_scheme))
         case ('ngp', 'NGP', 'neareast grid point')
            call pset%set_map(I_NGP)
         case ('cic', 'CIC', 'cloud in cell')
            call pset%set_map(I_CIC)
         case ('tsc', 'TSC', 'triangular shaped cloud')
            call pset%set_map(I_TSC)
         case default
            write(msg, '(3a)')"[particle_pub:init_particles] Unknown interpolation scheme '",trim(interpolation_scheme),"'"
            call die(msg)
      end select

#ifdef NBODY
      is_setacc_int = .false.
      is_setacc_cic = .false.
      is_setacc_tsc = .false.
      select case (acc_interp_method)
         case ('lagrange', 'Lagrange', 'polynomials')
            is_setacc_int = .true.
            order = 4
            call check_ord(order)
            call printinfo("[particle_pub:init_particles] Acceleration interpolation method: Lagrange polynomials")
         case ('cic', 'CIC')
            is_setacc_cic = .true.
            call printinfo("[particle_pub:init_particles] Acceleration interpolation method: CIC")
         case ('tsc', 'TSC')
            is_setacc_tsc = .true.
            call printinfo("[particle_pub:init_particles] Acceleration interpolation method: TSC")
      end select
#endif /* NBODY */

   end subroutine init_particles

!> \brief Deallocate pset

   subroutine cleanup_particles

     implicit none

     call pset%cleanup

   end subroutine cleanup_particles

end module particle_pub