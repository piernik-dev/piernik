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
!    HLLD Riemann solver for ideal magnetohydrodynamics
!    Varadarajan Parthasarathy, CAMK, Warszawa. 2015.
!
!--------------------------------------------------------------------------------------------------------------

#include "piernik.def"

module fluidupdate
! pulled by RIEMANN

  use hlld,  only: fluxes, riemann_hlld
  implicit none
  private
  public :: fluid_update, sweep_dsplit, rk2
  

contains

  subroutine fluid_update

    use cg_list,      only: cg_list_element
    use cg_leaves,    only: leaves
    use constants,    only: xdim, zdim
    use domain,       only: dom
    use global,       only: dt, dtm, t
    use user_hooks,   only: problem_customize_solution
    use dataio_pub,    only: halfstep

    implicit none

    logical, save                   :: first_run = .true.
    type(cg_list_element), pointer  :: cgl
    integer(kind=4)                 :: ddim
    !integer(kind=4)                 :: cdim

    halfstep = .false.
    if (first_run) then
       dtm = 0.0

    else

       dtm = dt

    endif

    t = t + dt

    cgl => leaves%first

    do while (associated(cgl))
       do ddim = xdim, zdim
          !if (dom%has_dir(ddim)) call sweep_dsplit(cgl%cg,dt,ddim, cdim)
          if (dom%has_dir(ddim)) call sweep_dsplit(cgl%cg,dt,ddim)
       enddo

       if (associated(problem_customize_solution)) call problem_customize_solution(.true.)
       cgl => cgl%nxt
       enddo

     t = t + dt
     dtm = dt
     halfstep = .true.

     cgl => leaves%first

     do while (associated(cgl))
        do ddim = zdim, xdim, -1
           !if (dom%has_dir(ddim)) call sweep_dsplit(cgl%cg,dt,ddim,cdim)
           if (dom%has_dir(ddim)) call sweep_dsplit(cgl%cg,dt,ddim)
        enddo
        if (associated(problem_customize_solution)) call problem_customize_solution(.false.)
        cgl => cgl%nxt
     enddo

       if (first_run) first_run = .false.

  end subroutine fluid_update

!-------------------------------------------------------------------------------------------------------------------

  !subroutine sweep_dsplit(cg, dt, ddim, cdim)
  subroutine sweep_dsplit(cg, dt, ddim)

    use constants,        only: pdims, xdim, zdim, ORTHO1, ORTHO2, LO, HI
    use all_boundaries,   only: all_fluid_boundaries
    use fluidindex,       only: iarr_all_swp
    use grid_cont,        only: grid_container
    use named_array_list, only: wna

    implicit none

    type(grid_container), pointer, intent(in) :: cg
    real,                          intent(in) :: dt
    integer(kind=4),               intent(in) :: ddim
    !integer(kind=4),               intent(in) :: cdim

    integer                                   :: n
    real, dimension(size(cg%u,1),cg%n_(ddim)) :: u1d
    !real, dimension(xdim:zdim, cg%n_(ddim)), pointer   :: b1d
    real, dimension(:,:), pointer             :: b1d
    real, dimension(xdim:zdim, cg%n_(ddim))   :: bb1d
    real, dimension(:,:), pointer             :: pu
    real, dimension(:), pointer               :: cs2
    integer                                   :: i1, i2
  
    

    do i2 = cg%lhn(pdims(ddim,ORTHO2),LO), cg%lhn(pdims(ddim,ORTHO2),HI)
       do i1 = cg%lhn(pdims(ddim, ORTHO1), LO), cg%lhn(pdims(ddim, ORTHO1), HI)
          pu => cg%w(wna%fi)%get_sweep(ddim,i1,i2)

          u1d(iarr_all_swp(ddim,:),:) = pu(:,:)

          !call rk2(n,u1d,b1d,bb1d,cs2,cdim,dt/cg%dl(ddim))
          call rk2(n,u1d,b1d,bb1d,cs2,dt/cg%dl(ddim))

          pu(:,:) = u1d(iarr_all_swp(ddim,:),:)

       enddo

    enddo

    call all_fluid_boundaries

  end subroutine sweep_dsplit

!-----------------------------------------------------------------------------------------------------------------------

  !subroutine rk2(n,u,b,bb,cs2,cdim,dtodx)
  subroutine rk2(n,u,b,bb,cs2,dtodx)

    use constants,   only: half, xdim, ydim, zdim
    use fluidindex,  only: flind, iarr_mag_swp
    use fluidtypes,  only: component_fluid


    implicit none

    !real,                 intent(in) :: dx, dt
    integer,                        intent(in)   :: n
    real, dimension(:,:),           intent(inout)  :: u
    real, dimension(:),   pointer,  intent(in)   :: cs2
    real, dimension(:,:),  intent(in)            :: bb
    real, dimension(:,:), pointer, intent(in)    :: b
    !integer(kind=4),       intent(in)            :: cdim
    real,                    intent(in)          :: dtodx

    class(component_fluid), pointer              :: fl
    !real, dimension(:,:),                        :: u_predict
    real, dimension(n, flind%all)                :: u_predict
    real, dimension(size(u,1),size(u,2)), target :: flx
    real, dimension(:,:), pointer                :: p_flx, p_b
    integer                                      :: nx, p
    integer(kind=4)                              :: ddim
    !real                                         :: dt, dx
    integer(kind=4)                              :: ibx, iby, ibz
    

    nx    = size(u,2)

    ibx = iarr_mag_swp(ddim,xdim)
    iby = iarr_mag_swp(ddim,ydim)
    ibz = iarr_mag_swp(ddim,zdim)

 
    !flx = fluxes(u0,b,bb,cs2,ddim)
    flx  = fluxes(u,b,bb,cs2,ddim)

    u_predict  =  u + dtodx*flx(:,:)

    do p = 1, flind%fluids
       fl    => flind%all_fluids(p)%fl
       p_flx => flx(fl%beg:fl%end,:)
       p_b   => b(ibx:ibz,:)
       call riemann_hlld(nx, p_flx, p_b, fl%gam, ddim)
    end do

    !u  =  half*(u0 + u_predict + dtodx*flx(:,:))
    u  =  half*(u + u_predict + dtodx*flx(:,:))

    
  end subroutine rk2

!---------------------------------------------------------------------------------------------------------------------------------------------------------------------
end module fluidupdate
