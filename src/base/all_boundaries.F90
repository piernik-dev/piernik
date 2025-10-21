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
!! \brief Global update of boundary routines
!!
!! \details moved to a separate file due to concerns about circular dependencies
!<
module all_boundaries

   implicit none

   private
   public :: all_bnd, all_bnd_vital_q, all_fluid_boundaries
#ifdef MAGNETIC
   public :: all_mag_boundaries
#endif /* MAGNETIC */

contains

!>
!! Subroutine calling all type boundaries after initialization of new run or restart reading
!! \todo make sure that all_fluid_boundaries and all_mag_boundaries can handle BND_USER boundaries right now, or do the boundaries later
!<

   subroutine all_bnd

      implicit none

      call all_fluid_boundaries
#ifdef MAGNETIC
      call all_mag_boundaries
#endif /* MAGNETIC */

   end subroutine all_bnd

   subroutine all_bnd_vital_q

      use cg_leaves,        only: leaves
      use named_array_list, only: qna
      use ppp,              only: ppp_main

      implicit none

      integer(kind=4) :: iq
      character(len=*), parameter :: abq_label = "all_boundaries_vital_q"

      call ppp_main%start(abq_label)

      do iq = lbound(qna%lst(:), dim=1, kind=4), ubound(qna%lst(:), dim=1, kind=4)
         if (qna%lst(iq)%vital) call leaves%leaf_arr3d_boundaries(iq)
      enddo

      call ppp_main%stop(abq_label)

   end subroutine all_bnd_vital_q

   subroutine all_fluid_boundaries(dir, nocorners, istep)

      use cg_leaves,          only: leaves
!      use cg_level_finest,    only: finest
      use constants,          only: xdim, zdim, uh_n, first_stage
      use global,             only: integration_order
      use domain,             only: dom
      use named_array_list,   only: wna
      use ppp,                only: ppp_main

      implicit none

      integer(kind=4), optional, intent(in) :: dir       !< select only this direction
      logical,         optional, intent(in) :: nocorners !< .when .true. then don't care about proper edge and corner update
      integer,         optional, intent(in) :: istep

      integer(kind=4) :: d, ind
      character(len=*), parameter :: abf_label = "all_fluid_boundaries"

      if (present(dir)) then
         if (.not. dom%has_dir(dir)) return
      endif

      call ppp_main%start(abf_label)

!      call finest%level%restrict_to_base

      ! should be more selective (modified leaves?)
      ind = wna%fi
      if (present(istep)) then
         if (istep == first_stage(integration_order)) ind = wna%ind(uh_n)
      endif
      call leaves%leaf_arr4d_boundaries(ind, dir=dir, nocorners=nocorners)

      if (present(dir)) then
         call leaves%bnd_u(dir)
      else
         do d = xdim, zdim
            if (dom%has_dir(d)) call leaves%bnd_u(d)
         enddo
      endif

      call ppp_main%stop(abf_label)

   end subroutine all_fluid_boundaries

#ifdef MAGNETIC
   subroutine all_mag_boundaries(istep)

      use cg_leaves,        only: leaves
!!$      use cg_list_global,   only: all_cg
      use constants,        only: xdim, zdim, psi_n, BND_INVALID, PPP_MAG, psih_n, magh_n, first_stage
      use domain,           only: dom
      use global,           only: psi_bnd, integration_order
      use named_array_list, only: wna, qna
      use ppp,              only: ppp_main

      implicit none

      integer, optional, intent(in) :: istep

      integer(kind=4) :: dir, ind
      character(len=*), parameter :: abm_label = "all_mag_boundaries"

      call ppp_main%start(abm_label, PPP_MAG)


      do dir = xdim, zdim
         if (dom%has_dir(dir)) call leaves%bnd_b(dir)
      enddo

      ind = wna%bi
      if (present(istep)) then
         if (istep == first_stage(integration_order)) ind = wna%ind(magh_n)
      endif
      call leaves%leaf_arr4d_boundaries(ind)

      if (qna%exists(psi_n)) then  ! assumed that qna%exists(psih_n) too
         ind = qna%ind(psi_n)
         if (present(istep)) then
            if (istep == first_stage(integration_order)) ind = qna%ind(psih_n)
         endif

         call leaves%leaf_arr3d_boundaries(ind)
         if (psi_bnd == BND_INVALID) then
            call leaves%external_boundaries(ind)
         else
            call leaves%external_boundaries(ind, bnd_type=psi_bnd)
         endif
      endif

      call ppp_main%stop(abm_label, PPP_MAG)

   end subroutine all_mag_boundaries
#endif /* MAGNETIC */

end module all_boundaries
