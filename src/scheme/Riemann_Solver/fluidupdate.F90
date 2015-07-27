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

  implicit none
  private
  public :: fluid_update

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

          if (dom%has_dir(ddim)) call sweep(cgl%cg,dt,ddim)

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

           if (dom%has_dir(ddim)) call sweep(cgl%cg,dt,ddim)

        enddo

        if (associated(problem_customize_solution)) call problem_customize_solution(.false.)

        cgl => cgl%nxt

     enddo

       if (first_run) first_run = .false.

  end subroutine fluid_update

!-------------------------------------------------------------------------------------------------------------------

  function utoq(u,b) result(q)

    use constants,  only: half, two
    use fluidindex, only: flind
    use fluidtypes, only: component_fluid
    use func,       only: ekin

    implicit none

    real, dimension(:,:), intent(in)      :: u,b

    real, dimension(size(u,1),size(u,2))  :: q

    integer :: p

    class(component_fluid), pointer       :: fl

    do p = 1, flind%fluids
       fl => flind%all_fluids(p)%fl

       q(fl%idn,:)  =  u(fl%idn,:)
       q(fl%imx,:)  =  u(fl%idn,:)/u(fl%imx,:)
       q(fl%imy,:)  =  u(fl%idn,:)/u(fl%imy,:)
       q(fl%imz,:)  =  u(fl%idn,:)/u(fl%imz,:)

       if (fl%has_energy) then

          q(fl%ien,:) = fl%gam_1*(u(fl%ien,:) - ekin(u(fl%imx,:), u(fl%imy,:), u(fl%imz,:), u(fl%idn,:)))

          if (fl%is_magnetized) q(fl%ien,:) = q(fl%ien,:) + (two-fl%gam)*half*sum(b(:,:)**2,dim=1)

       endif

    enddo

  end function utoq

!-------------------------------------------------------------------------------------------------------------------------

  subroutine sweep(cg,dt,ddim)

    use constants,        only: pdims, xdim, zdim, cs_i2_n, ORTHO1, ORTHO2, LO, HI
    use all_boundaries,   only: all_fluid_boundaries
    use fluidindex,       only: iarr_all_swp
    use grid_cont,        only: grid_container
    use named_array_list, only: qna, wna

    implicit none

    type(grid_container), pointer, intent(in) :: cg
    real,                          intent(in) :: dt
    integer(kind=4),               intent(in) :: ddim

    ! Sweep later ..

  end subroutine sweep

!-----------------------------------------------------------------------------------------------------------------------

  subroutine rk2()

    ! later
    
  end subroutine rk2

!---------------------------------------------------------------------------------------------------------------------------------------------------------------------
end module fluidupdate
