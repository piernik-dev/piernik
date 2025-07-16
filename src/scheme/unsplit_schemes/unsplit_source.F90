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

    subroutine apply_source(cg,istep)
        use grid_cont,          only: grid_container
        use named_array_list,   only: wna, qna
        use constants,          only: pdims, ORTHO1, ORTHO2, LO, HI, magh_n, uh_n, rk_coef, cs_i2_n, first_stage, last_stage, xdim, ydim, zdim
        use global,             only: dt, integration_order, nstep
        use domain,             only: dom
        use fluidindex,         only: flind, iarr_all_dn, iarr_all_mx, iarr_all_swp, iarr_mag_swp
        use sources,            only: internal_sources, care_for_positives

        implicit none

        type(grid_container), pointer,     intent(in) :: cg
        integer,                           intent(in) :: istep

        integer                                                     :: ddim, i1, i2, uhi, bhi
        real, dimension(:,:),allocatable                            :: u
        real, dimension(:,:), pointer                               :: pu,pb
        real, allocatable, target                                   :: vx(:,:)
        real, dimension(1, 1)                                       :: b_ugly ! ugly
        real, dimension(:,:),allocatable                            :: b
        real, dimension(:,:),allocatable                            :: u1

        uhi = wna%ind(uh_n)
        bhi = wna%ind(magh_n)

        do ddim=xdim,zdim

            if (.not. allocated(u)) then
                allocate(u(cg%n_(ddim), size(cg%u,1)))
            else
                deallocate(u)
                allocate(u(cg%n_(ddim), size(cg%u,1)))

            endif
            if (.not. allocated(b)) then
                allocate(b(cg%n_(ddim), size(cg%b,1)))
            else
                deallocate(b)
                allocate(b(cg%n_(ddim), size(cg%b,1)))

            endif
            if (.not. allocated(vx)) then
                allocate(vx(size(u,1), flind%fluids))
            else
                deallocate(vx)
                allocate(vx(size(u,1), flind%fluids))
            endif
            if (.not. allocated(u1)) then
                allocate(u1(size(u, 1),size(u, 2)))
            else
                deallocate(u1)
                allocate(u1(size(u, 1),size(u, 2)))
            endif
            if (.not. dom%has_dir(ddim)) cycle
            do i2 = cg%ijkse(pdims(ddim, ORTHO2), LO), cg%ijkse(pdims(ddim, ORTHO2), HI)
                do i1 = cg%ijkse(pdims(ddim, ORTHO1), LO), cg%ijkse(pdims(ddim, ORTHO1), HI)  

                    pu => cg%w(uhi)%get_sweep(ddim,i1,i2)
                    pb => cg%w(bhi)%get_sweep(ddim,i1,i2)
                    if (istep == first_stage(integration_order) .or. integration_order < 2 ) then
                            pu => cg%w(wna%fi)%get_sweep(ddim,i1,i2)
                            pb => cg%w(wna%bi)%get_sweep(ddim,i1,i2)
                    endif
                    

                    u(:, iarr_all_swp(ddim,:)) = transpose(pu(:,:))
                    b(:, iarr_mag_swp(ddim,:)) = transpose(pb(:,:))

                    u1 = u

                    vx = u(:, iarr_all_mx) / u(:, iarr_all_dn) ! this may also be useful for gravitational acceleration
#ifdef MAGNETIC
                    call internal_sources(size(u, 1, kind=4), u, u1, b, cg, istep, ddim, i1, i2, rk_coef(istep) * dt, vx)

                    call care_for_positives(size(u, 1, kind=4), u1, b, cg, ddim, i1, i2)
#else
                    call internal_sources(size(u, 1, kind=4), u, u1, b_ugly, cg, istep, ddim, i1, i2, rk_coef(istep) * dt, vx)

                    call care_for_positives(size(u, 1, kind=4), u1, b_ugly, cg, ddim, i1, i2)
#endif /* MAGNETIC */
                    pu(:,:) = transpose(u1(:, iarr_all_swp(ddim,:)))
                end do        
            end do

            deallocate(vx,u1,u,b)
            
        end do
    end subroutine apply_source
end module unsplit_source 