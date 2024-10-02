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
!! \brief Implementation of the non-inertial forces for rotating grid
!<

module non_inertial

   implicit none

   private
   public  :: init_non_inertial, non_inertial_force, get_omega

   real    :: omega !< angular frequency around the z axis

contains

!>
!! \brief perform basic checks
!<

   subroutine init_non_inertial

      use bcast,      only: piernik_MPI_Bcast
      use constants,  only: PIERNIK_INIT_GRID, GEO_XYZ
      use dataio_pub, only: die, code_progress, nh
      use domain,     only: dom
      use mpisetup,   only: rbuff, master, slave
#ifdef VERBOSE
      use dataio_pub, only: printinfo
#endif /* VERBOSE */

      implicit none

      namelist /NONINERTIAL/ omega

#ifdef VERBOSE
      if (master) call printinfo("[non_inertial:init_non_inertial] Commencing module initialization")
#endif /* VERBOSE */

      omega = 0.0

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=NONINERTIAL)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=''
         read(unit=nh%lun, nml=NONINERTIAL, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "NONINERTIAL")
         read(nh%cmdl_nml,nml=NONINERTIAL, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "NONINERTIAL", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=NONINERTIAL)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1) = omega

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         omega    = rbuff(1)

      endif

      if (code_progress < PIERNIK_INIT_GRID) call die("[non_inertial:init_non_inertial] grid not initialized.") ! this is a weak check, the real dependency is init_geometry at the moment

      if (dom%geometry_type /= GEO_XYZ) call die("[non_inertial:init_non_inertial] Only cartesian geometry is implemented")

   end subroutine init_non_inertial

!>
!! \brief Compute the non-inertial acceleration for a given row of cells.
!!
!! The equations for the non-inertial acceleration in the x and y directions are given by:
!! \f{equation}
!! acc_x = 2 * \Omega * v_y + \Omega^2 * x
!! acc_y =  -2 * \Omega * v_x + \Omega^2 * y
!! \f}
!!
!! \details This is a low-order estimate of the accelerations, because this routine uses density and velocity fields
!! from the beginning of the time step. This is a simple approach, but ignores any changes due to other accelerations during the time step.
!<

   function non_inertial_force(sweep, u, cg) result(rotacc)

      use fluidindex, only: flind, iarr_all_dn, iarr_all_my
      use constants,  only: xdim, ydim !, zdim
      use grid_cont,  only: grid_container

      implicit none

      integer(kind=4), intent(in)               :: sweep  !< string of characters that points out the current sweep direction
      real, dimension(:,:), intent(in)          :: u      !< current fluid state vector
      type(grid_container), pointer, intent(in) :: cg     !< current grid piece
      integer                                   :: ifl
      real, dimension(size(u,1), flind%fluids)  :: rotacc !< an array for non-inertial accelerations

      ! non-inertial (Coriolis and centrifugal) forces for corotating coords
      do ifl = 1, flind%fluids
         select case (sweep)
            case (xdim)
               rotacc(:, ifl) = +2.0 * omega * u(:, iarr_all_my(ifl))/u(:, iarr_all_dn(ifl)) + omega**2 * cg%x
            case (ydim)
               rotacc(:, ifl) = -2.0 * omega * u(:, iarr_all_my(ifl))/u(:, iarr_all_dn(ifl)) + omega**2 * cg%y
   !         case (zdim) !no z-component of non-inertial forces
            case default
               rotacc(:, ifl) = 0.0
         end select
      enddo

   end function non_inertial_force

!>
!! \brief return angular rotation parameter
!<


   real function get_omega()

      implicit none

      get_omega = omega

   end function get_omega

end module non_inertial
