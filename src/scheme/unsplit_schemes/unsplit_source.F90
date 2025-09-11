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

module unsplit_source

! pulled by ANY

   implicit none

   private
   public  :: apply_source

contains

! This routine has to conform to the interface defined in sweeps::sweep

   subroutine apply_source(cg, istep)

      use grid_cont,          only: grid_container
      use named_array_list,   only: wna
      use constants,          only: pdims, ORTHO1, ORTHO2, LO, HI, uh_n, rk_coef, first_stage, xdim, zdim, last_stage
      use global,             only: dt, integration_order
      use domain,             only: dom
      use fluidindex,         only: flind, iarr_all_dn, iarr_all_mx, iarr_all_swp
      use sources,            only: internal_sources, care_for_positives
      use diagnostics,        only: my_allocate, my_deallocate
#ifdef MAGNETIC
      use constants,          only: magh_n
      use fluidindex,         only: iarr_mag_swp
#endif /* MAGNETIC */

        implicit none

      type(grid_container), pointer,     intent(in) :: cg
      integer,                           intent(in) :: istep

      integer(kind=4)                                             :: uhi, ddim
      integer                                                     :: i1, i2
      real, dimension(:,:),allocatable                            :: u
      real, dimension(:,:), pointer                               :: pu
      real, allocatable, target                                   :: vx(:,:)
      real, dimension(:,:),allocatable                            :: u1
#ifdef MAGNETIC
      real, dimension(:,:), pointer                               :: pb
      real, dimension(:,:),allocatable                            :: b
      integer                                                     :: bhi
      bhi = wna%ind(magh_n)
#else /* !MAGNETIC */
      real, dimension(1, 1)                                       :: b_ugly ! ugly

      b_ugly = 0.0
#endif /* !MAGNETIC */

      uhi = wna%ind(uh_n)
      do ddim = xdim, zdim
         if (.not. dom%has_dir(ddim)) cycle
         call my_allocate(u, [cg%n_(ddim), size(cg%u, 1, kind=4)])
         call my_allocate(vx, [size(u, 1, kind=4), flind%fluids])
         call my_allocate(u1, [size(u, 1, kind=4),size(u, 2, kind=4)])
#ifdef MAGNETIC
         call my_allocate(b, [cg%n_(ddim), size(cg%b, 1, kind=4)])
#endif /* MAGNETIC */
         do i2 = cg%ijkse(pdims(ddim, ORTHO2), LO), cg%ijkse(pdims(ddim, ORTHO2), HI)
            do i1 = cg%ijkse(pdims(ddim, ORTHO1), LO), cg%ijkse(pdims(ddim, ORTHO1), HI)
               pu => cg%w(wna%fi)%get_sweep(ddim, i1, i2)
#ifdef MAGNETIC
               pb => cg%w(wna%bi)%get_sweep(ddim, i1, i2)
#endif /* MAGNETIC */
               if (istep == first_stage(integration_order) .or. integration_order < 2) then
                  pu => cg%w(uhi)%get_sweep(ddim, i1, i2)
#ifdef MAGNETIC
                  pb => cg%w(bhi)%get_sweep(ddim, i1, i2)
#endif /* MAGNETIC */
                    endif
                  u(:, iarr_all_swp(ddim,:)) = transpose(pu(:,:))
#ifdef MAGNETIC
                  b(:, iarr_mag_swp(ddim,:)) = transpose(pb(:,:))
#endif /* MAGNETIC */
                  u1 = u
                  vx = u(:, iarr_all_mx) / u(:, iarr_all_dn) ! this may also be useful for gravitational acceleration
#ifdef MAGNETIC
                  call internal_sources(size(u, 1, kind=4), u, u1, b, cg, istep, ddim, i1, i2, rk_coef(istep) * dt, vx)

                  if (istep == last_stage(integration_order)) then
                     call care_for_positives(size(u, 1, kind=4), u1, b, cg, ddim, i1, i2)
                  endif
#else
                  call internal_sources(size(u, 1, kind=4), u, u1, b_ugly, cg, istep, ddim, i1, i2, rk_coef(istep) * dt, vx)

                  if (istep == last_stage(integration_order)) then
                     call care_for_positives(size(u, 1, kind=4), u1, b_ugly, cg, ddim, i1, i2)
                  endif

#endif /* MAGNETIC */
                  pu(:,:) = transpose(u1(:, iarr_all_swp(ddim,:)))

               enddo

            enddo

            call my_deallocate(vx); call my_deallocate(u1); call my_deallocate(u)
#ifdef MAGNETIC
            call my_deallocate(b)
#endif /* MAGNETIC */
      enddo

   end subroutine apply_source

end module unsplit_source
