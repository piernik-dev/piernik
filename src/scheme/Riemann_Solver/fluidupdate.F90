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
!    Dr. Artur Gawryszczak, CAMK, Warszawa.
!
!    Energy fix up routines for CT and its related comments are not used in the current version.
!    The algorithm is simply present for experimental purposes.
!--------------------------------------------------------------------------------------------------------------

#include "piernik.h"

module fluidupdate
! pulled by RIEMANN

  implicit none
  private
  public :: fluid_update

contains

!>
!! \brief Advance the solution by two timesteps using directional splitting
!!
!! Spaghetti warning: copied from rtvd_split
!<

  subroutine fluid_update

    use constants,    only: DIVB_HDC
    use dataio_pub,   only: halfstep
    use fluxlimiters, only: set_limiters
    use global,       only: dt, dtm, t, limiter, limiter_b, divB_0_method
    use mass_defect,  only: update_magic_mass
    use hdc,          only: update_chspeed

    implicit none

    logical, save   :: first_run = .true.

!!!call repeat_fluidstep

    if (divB_0_method == DIVB_HDC) call update_chspeed
    if (first_run) call set_limiters(limiter, limiter_b)

    halfstep = .false.
    t = t + dt
    call make_3sweeps(.true.) ! X -> Y -> Z

! Sources should be hooked to problem_customize_solution with forward argument

    halfstep = .true.
    t = t + dt
    dtm = dt
    call make_3sweeps(.false.) ! Z -> Y -> X
    call update_magic_mass

    if (first_run) first_run = .false.

  end subroutine fluid_update

!-------------------------------------------------------------------------------------------------------------------

  subroutine make_3sweeps(forward)

    use constants,      only: xdim, zdim, I_ONE, DIVB_HDC
    use global,         only: divB_0_method
    use hdc,            only: glmdamping, eglm
    use user_hooks,     only: problem_customize_solution
#ifdef GRAV
      use global,              only: t, dt
      use gravity,             only: source_terms_grav
      use particle_pub,        only: pset, psolver
#endif /* GRAV */
#if defined(COSM_RAYS) && defined(MULTIGRID)
    use all_boundaries,      only: all_fluid_boundaries
    use initcosmicrays,      only: use_split
    use multigrid_diffusion, only: multigrid_solve_diff
#endif /* COSM_RAYS && MULTIGRID */
#ifdef SHEAR
      use shear,               only: shear_3sweeps
#endif /* SHEAR */

    implicit none

    logical, intent(in) :: forward  !< If .true. then do X->Y->Z sweeps, if .false. then reverse that order

    integer(kind=4)                 :: ddim

#ifdef SHEAR
      call shear_3sweeps
#endif /* SHEAR */

#ifdef GRAV
      call source_terms_grav
#endif /* GRAV */

#if defined(COSM_RAYS) && defined(MULTIGRID)
    if (.not. use_split) then
       call multigrid_solve_diff
       call all_fluid_boundaries
    endif
#endif /* COSM_RAYS && MULTIGRID */

    if (forward) then
       do ddim = xdim, zdim, 1
          call make_sweep(ddim, forward)
       enddo
    else
       do ddim = zdim, xdim, -I_ONE
          call make_sweep(ddim, forward)
       enddo
    endif
#ifdef GRAV
      if (associated(psolver)) call pset%evolve(psolver, t-dt, dt)
#endif /* GRAV */
   if (associated(problem_customize_solution)) call problem_customize_solution(forward)

    call eglm
    if (divB_0_method == DIVB_HDC) call glmdamping

  end subroutine make_3sweeps

!-------------------------------------------------------------------------------------------------------------------

  subroutine make_sweep(dir, forward)

    use domain,           only: dom
    use global,           only: force_cc_mag
    use sweeps,           only: sweep
#ifdef COSM_RAYS
    use crdiffusion,      only: cr_diff
    use initcosmicrays,   only: use_split
#endif /* COSM_RAYS */
#ifdef MAGNETIC
    use constants,        only: DIVB_CT
    use ct,               only: magfield
    use global,           only: divB_0_method
#endif /* MAGNETIC */

    implicit none

    integer(kind=4), intent(in) :: dir      !< direction, one of xdim, ydim, zdim
    logical,         intent(in) :: forward  !< if .false. then reverse operation order in the sweep

    ! ToDo: check if changes of execution order here (block loop, direction loop, boundary update can change
    ! cost or allow for reduction of required guardcells
    if (.not. force_cc_mag) call bfc2bcc

    if (dom%has_dir(dir)) then
       if (.not. forward) then
#ifdef COSM_RAYS
          if (use_split) call cr_diff(dir)
#endif /* COSM_RAYS */
#ifdef MAGNETIC
          if (divB_0_method == DIVB_CT) call magfield(dir)
#endif /* MAGNETIC */
       endif

       call sweep(dir)

       if (forward) then
#ifdef MAGNETIC
          if (divB_0_method == DIVB_CT) call magfield(dir)
#endif /* MAGNETIC */
#ifdef COSM_RAYS
          if (use_split) call cr_diff(dir)
#endif /* COSM_RAYS */
       endif

    endif

  end subroutine make_sweep

!-------------------------------------------------------------------------------------------------------------------

  subroutine bfc2bcc

     use cg_leaves,        only: leaves
     use cg_list,          only: cg_list_element
     use constants,        only: xdim, ydim, zdim, LO, HI, half
     use dataio_pub,       only: die
     use domain,           only: dom
     use global,           only: force_cc_mag
     use grid_cont,        only: grid_container
     use named_array_list, only: wna

     implicit none

     type(cg_list_element), pointer :: cgl
     type(grid_container),  pointer :: cg

     if (force_cc_mag) call die("[fluidupdate:bfc2bcc] no  point in converting cell-centered magnetic field to cell centers like it was face-centered")

     cgl => leaves%first
     do while (associated(cgl))
        cg => cgl%cg

        cg%w(wna%bcci)%arr(:,:,:,:) = half * cg%b(:, :, :, :)

        cg%w(wna%bcci)%arr(xdim, cg%lhn(xdim, LO):cg%lhn(xdim, HI)-1, :, :) = &
             cg%w(wna%bcci)%arr(xdim, cg%lhn(xdim, LO):cg%lhn(xdim, HI)-1, :, :) + &
             half * cg%b(xdim, cg%lhn(xdim, LO)+dom%D_x:cg%lhn(xdim, HI)-1+dom%D_x, :, :)

        cg%w(wna%bcci)%arr(ydim, :,cg%lhn(ydim, LO):cg%lhn(ydim, HI)-1, :) = &
             cg%w(wna%bcci)%arr(ydim, :, cg%lhn(ydim, LO):cg%lhn(ydim, HI)-1, :) + &
             half * cg%b(ydim, :, cg%lhn(ydim, LO)+dom%D_y:cg%lhn(ydim, HI)-1+dom%D_y, :)

        cg%w(wna%bcci)%arr(zdim, :, :, cg%lhn(zdim, LO):cg%lhn(zdim, HI)-1) = &
             cg%w(wna%bcci)%arr(zdim, :, :, cg%lhn(zdim, LO):cg%lhn(zdim, HI)-1) + &
             half * cg%b(zdim, :, :, cg%lhn(zdim, LO)+dom%D_z:cg%lhn(zdim, HI)-1+dom%D_z)

        cgl => cgl%nxt
     enddo

  end subroutine bfc2bcc

end module fluidupdate
