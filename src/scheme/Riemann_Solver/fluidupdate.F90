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

#include "piernik.def"

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

    use constants,    only: GEO_XYZ, DIVB_HDC, I_ZERO, I_ONE
    use dataio_pub,   only: halfstep, die
    use domain,       only: dom, is_refined
    use fluidindex,   only: flind
    use fluxlimiters, only: set_limiters
    use global,       only: dt, dtm, t, limiter, limiter_b, divB_0_method
    use mass_defect,  only: update_magic_mass
    use hdc,          only: update_chspeed

    implicit none

    integer(kind=4) :: nmag, i
    logical, save   :: first_run = .true.

    ! is_multicg should be safe
    if (is_refined) call die("[fluid_update] This Rieman solver is not compatible with mesh refinements yet!")
    if (dom%geometry_type /= GEO_XYZ) call die("[fluid_update] Non-cartesian geometry is not implemented yet in this Riemann solver.")
#ifdef ISO
#  error Isothermal EOS is not implemented yet in this Riemann solver.
#endif /* ISO */
    nmag = I_ZERO
    do i = 1, flind%fluids
       if (flind%all_fluids(i)%fl%is_magnetized) nmag = nmag + I_ONE
    enddo
    if (nmag > 1) call die("[fluidupdate:fluid_update] At most one magnetized fluid is implemented")

!!!call repeat_fluidstep

    if (divB_0_method == DIVB_HDC) call update_chspeed
    halfstep = .false.
    if (first_run) then
       dtm = 0.0
       call set_limiters(limiter, limiter_b)
    else
       dtm = dt
    endif
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

#ifdef MAGNETIC
   subroutine magfield(dir)

      use ct,          only: advectb
      use constants,   only: ndims, I_ONE
      use dataio_pub,  only: die
      use global,      only: force_cc_mag
#ifdef RESISTIVE
      use resistivity, only: diffuseb
#endif /* RESISTIVE */

       implicit none

       integer(kind=4), intent(in) :: dir

       integer(kind=4)             :: bdir, dstep

       if (force_cc_mag) call die("[fluidupdate:magfield] forcing cell-centered magnetic field is not allowed for constrained transport")

       do dstep = 0, 1
          bdir  = I_ONE + mod(dir+dstep,ndims)
          call advectb(bdir, dir)
#ifdef RESISTIVE
         call diffuseb(bdir, dir)
#endif /* RESISTIVE */
         call mag_add(dir, bdir)
      enddo

   end subroutine magfield

!-------------------------------------------------------------------------------------------------------------------

 subroutine mag_add(dim1, dim2)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use dataio_pub,       only: die
      use global,           only: force_cc_mag
      use grid_cont,        only: grid_container
      use all_boundaries,   only: all_mag_boundaries
      use user_hooks,       only: custom_emf_bnd
#ifdef RESISTIVE
     use constants,        only: wcu_n
     use dataio_pub,       only: die
     use domain,           only: is_multicg
     use named_array_list, only: qna
#endif /* RESISTIVE */

      implicit none

      integer(kind=4), intent(in)    :: dim1, dim2

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
#ifdef RESISTIVE
      real, dimension(:,:,:), pointer :: wcu
#endif /* RESISTIVE */

      if (force_cc_mag) call die("[fluidupdate:mag_add] forcing cell-centered magnetic field is not allowed for constrained transport")

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
#ifdef RESISTIVE
! DIFFUSION FULL STEP
         wcu => cg%q(qna%ind(wcu_n))%arr
         if (is_multicg) call die("[fluidupdate:mag_add] multiple grid pieces per processor not implemented yet") ! not tested custom_emf_bnd
         if (associated(custom_emf_bnd)) call custom_emf_bnd(wcu)
         cg%b(dim2,:,:,:) = cg%b(dim2,:,:,:) -              wcu*cg%idl(dim1)
         cg%b(dim2,:,:,:) = cg%b(dim2,:,:,:) + pshift(wcu,dim1)*cg%idl(dim1)
         cg%b(dim1,:,:,:) = cg%b(dim1,:,:,:) +              wcu*cg%idl(dim2)
         cg%b(dim1,:,:,:) = cg%b(dim1,:,:,:) - pshift(wcu,dim2)*cg%idl(dim2)
#endif /* RESISTIVE */
! ADVECTION FULL STEP
         if (associated(custom_emf_bnd)) call custom_emf_bnd(cg%wa)
         cg%b(dim2,:,:,:) = cg%b(dim2,:,:,:) - cg%wa*cg%idl(dim1)
         cg%wa = mshift(cg%wa,dim1)
         cg%b(dim2,:,:,:) = cg%b(dim2,:,:,:) + cg%wa*cg%idl(dim1)
         cg%b(dim1,:,:,:) = cg%b(dim1,:,:,:) - cg%wa*cg%idl(dim2)
         cg%wa = pshift(cg%wa,dim2)
         cg%b(dim1,:,:,:) = cg%b(dim1,:,:,:) + cg%wa*cg%idl(dim2)
         cgl => cgl%nxt
      enddo

      call all_mag_boundaries

   end subroutine mag_add

!-------------------------------------------------------------------------------------------------------------------

!>
!! \brief Function pshift makes one-cell, forward circular shift of 3D array in any direction
!! \param tab input array
!! \param d direction of the shift, where 1,2,3 corresponds to \a x,\a y,\a z respectively
!! \return real, dimension(size(tab,1),size(tab,2),size(tab,3))
!!
!! The function was written in order to significantly improve
!! the performance at the cost of the flexibility of original \p CSHIFT.
!<
   function pshift(tab, d)

      use dataio_pub,    only: warn

      implicit none

      real, dimension(:,:,:), intent(inout) :: tab
      integer(kind=4),        intent(in)    :: d

      integer :: ll
      real, dimension(size(tab,1),size(tab,2),size(tab,3)) :: pshift

      ll = size(tab,d)

      if (ll==1) then
         pshift = tab
         return
      endif

      if (d==1) then
         pshift(1:ll-1,:,:) = tab(2:ll,:,:); pshift(ll,:,:) = tab(1,:,:)
      else if (d==2) then
         pshift(:,1:ll-1,:) = tab(:,2:ll,:); pshift(:,ll,:) = tab(:,1,:)
      else if (d==3) then
         pshift(:,:,1:ll-1) = tab(:,:,2:ll); pshift(:,:,ll) = tab(:,:,1)
      else
         call warn('[fluidupdate:pshift]: Dim ill defined in pshift!')
      endif

      return
   end function pshift

!>
!! \brief Function mshift makes one-cell, backward circular shift of 3D array in any direction
!! \param tab input array
!! \param d direction of the shift, where 1,2,3 corresponds to \a x,\a y,\a z respectively
!! \return real, dimension(size(tab,1),size(tab,2),size(tab,3))
!!
!! The function was written in order to significantly improve
!! the performance at the cost of the flexibility of original \p CSHIFT.
!<
   function mshift(tab,d)

      use dataio_pub,    only: warn

      implicit none

      real, dimension(:,:,:), intent(inout) :: tab
      integer(kind=4),        intent(in)    :: d

      integer :: ll
      real, dimension(size(tab,1) , size(tab,2) , size(tab,3)) :: mshift

      ll = size(tab,d)

      if (ll==1) then
         mshift = tab
         return
      endif

      if (d==1) then
         mshift(2:ll,:,:) = tab(1:ll-1,:,:); mshift(1,:,:) = tab(ll,:,:)
      else if (d==2) then
         mshift(:,2:ll,:) = tab(:,1:ll-1,:); mshift(:,1,:) = tab(:,ll,:)
      else if (d==3) then
         mshift(:,:,2:ll) = tab(:,:,1:ll-1); mshift(:,:,1) = tab(:,:,ll)
      else
         call warn('[fluidupdate:mshift]: Dim ill defined in mshift!')
      endif

      return
   end function mshift

#endif /* MAGNETIC */

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
