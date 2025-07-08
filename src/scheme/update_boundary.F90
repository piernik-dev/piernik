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
!! The job this module is to update boundaries. This was previously inside the module sweeps . Moved here so that this can be called inside 
!! solver after each time step on the entire 
!<

module update_boundary

! pulled by ANY

   implicit none

   private
   public :: update_boundaries

contains

!>
!! \brief Call all boundaries, try to avoid unnecessary parts.
!!
!! For some reasons dir=cdim affect mcrwind tests if sweeps_mgu
!! \todo Find out why. Is it related to position of magnetic field components?
!!
!! \todo Once it gets simplified enough merge it back to sweep.
!<

   subroutine update_boundaries(cdim, istep)

      use all_boundaries, only: all_fluid_boundaries
!      use cg_leaves,      only: leaves
      use constants,      only: first_stage, DIVB_HDC,xdim,zdim
      use domain,         only: dom
      use global,         only: sweeps_mgu, integration_order, divB_0_method
#ifdef MAGNETIC
      use all_boundaries, only: all_mag_boundaries
#endif /* MAGNETIC */

      implicit none

      integer(kind=4),optional, intent(in) :: cdim
      integer,                  intent(in) :: istep


      integer                              :: ub_i


      if (.not. present(cdim) .or. cdim==-1) then
         do ub_i=xdim,zdim 
            if (dom%has_dir(ub_i)) then
               if (sweeps_mgu) then
                  if (istep == first_stage(integration_order)) then
                     call all_fluid_boundaries(nocorners = .true., dir = ub_i)
                  else
                     call all_fluid_boundaries(nocorners = .true.)
                  endif
               else
                  ! nocorners and dir = cdim can be used safely only when ord_fluid_prolong == 0 .and. cc_mag
                  ! essential speedups here are possible but it requires c/f boundary prolongation that does not require corners

                  ! if (istep == first_stage(integration_order)) then
                  !    call all_fluid_boundaries(nocorners = .true.)
                  ! else
                     call all_fluid_boundaries !(nocorners = .true., dir = cdim)
                  ! endif
               endif
            endif
         end do
      else
         if (dom%has_dir(cdim)) then
            if (sweeps_mgu) then
               if (istep == first_stage(integration_order)) then
                  call all_fluid_boundaries(nocorners = .true., dir = cdim)
               else
                  call all_fluid_boundaries(nocorners = .true.)
               endif
            else
               ! nocorners and dir = cdim can be used safely only when ord_fluid_prolong == 0 .and. cc_mag
               ! essential speedups here are possible but it requires c/f boundary prolongation that does not require corners

               ! if (istep == first_stage(integration_order)) then
               !    call all_fluid_boundaries(nocorners = .true.)
               ! else
                  call all_fluid_boundaries !(nocorners = .true., dir = cdim)
               ! endif
            endif
         endif
      endif 

      if (divB_0_method == DIVB_HDC) then
#ifdef MAGNETIC
         call all_mag_boundaries ! ToDo: take care of psi boundaries
#endif /* MAGNETIC */
      endif

   end subroutine update_boundaries

end module update_boundary