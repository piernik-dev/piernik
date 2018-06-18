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

    use constants,    only: GEO_XYZ, DIVB_HDC
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
#ifdef GRAV
#  error Graviy is not implemented yet in this Riemann solver.
#endif /* GRAV */
#ifdef ISO
#  error Isothermal EOS is not implemented yet in this Riemann solver.
#endif /* ISO */
    nmag = 0
    do i = 1, flind%fluids
       if (flind%all_fluids(i)%fl%is_magnetized) nmag = nmag + 1
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
    use hdc,            only: glmdamping!,eglm
    use user_hooks,     only: problem_customize_solution
#if defined(COSM_RAYS) && defined(MULTIGRID)
    use all_boundaries,      only: all_fluid_boundaries
    use initcosmicrays,      only: use_split
    use multigrid_diffusion, only: multigrid_solve_diff
#endif /* COSM_RAYS && MULTIGRID */

    implicit none

    logical, intent(in) :: forward  !< If .true. then do X->Y->Z sweeps, if .false. then reverse that order

    integer(kind=4)                 :: ddim

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
    if (associated(problem_customize_solution)) call problem_customize_solution(forward)

    if (divB_0_method == DIVB_HDC) then

       call glmdamping

    endif
  end subroutine make_3sweeps

!-------------------------------------------------------------------------------------------------------------------

  subroutine make_sweep(dir, forward)

    use all_boundaries,   only: all_bnd
    use cg_leaves,        only: leaves
    use cg_list,          only: cg_list_element
    use constants,        only: psi_n
    use domain,           only: dom
    use global,           only: skip_sweep, dt, force_cc_mag
    use named_array_list, only: qna
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

    type(cg_list_element), pointer  :: cgl

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

       cgl => leaves%first
       do while (associated(cgl))
          ! Warning: 2.5D MHD may need all directional calls anyway
          if (.not. skip_sweep(dir)) call sweep_dsplit(cgl%cg,dt,dir)
          cgl => cgl%nxt
       enddo
       if (qna%exists(psi_n)) call leaves%leaf_arr3d_boundaries(qna%ind(psi_n))

       call all_bnd
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

  subroutine sweep_dsplit(cg, dt, ddim)

    use constants,        only: pdims, xdim, zdim, ORTHO1, ORTHO2, LO, HI, psi_n, INVALID
    use fluidindex,       only: iarr_all_swp, iarr_mag_swp
    use global,           only: force_cc_mag
    use grid_cont,        only: grid_container
    use named_array_list, only: wna, qna
#ifdef COSM_RAYS
    use crhelpers,        only: div_v, set_div_v1d
    use fluidindex,       only: flind
#endif /* COSM_RAYS */

    implicit none

    type(grid_container), pointer, intent(in) :: cg
    real,                          intent(in) :: dt
    integer(kind=4),               intent(in) :: ddim

    real, dimension(size(cg%u,1), cg%n_(ddim)) :: u1d
    real, dimension(xdim:zdim, cg%n_(ddim))    :: b_cc1d
    real, dimension(1, cg%n_(ddim))            :: psi_d ! artificial rank-2 to conform to flux limiter interface
    real, dimension(:,:), pointer              :: pu, pb
    integer                                    :: i1, i2
    integer                                    :: bi
    real, dimension(:), pointer                :: ppsi
    integer                                    :: psii
    real, dimension(:), pointer                :: div_v1d => null()

    if (force_cc_mag) then
       bi = wna%bi
    else
       bi = wna%bcci
    endif

    psii = INVALID
    if (qna%exists(psi_n)) psii = qna%ind(psi_n)
    psi_d = 0.
    nullify(ppsi)

#ifdef COSM_RAYS
    call div_v(flind%ion%pos, cg)
#endif /* COSM_RAYS */

    do i2 = cg%lhn(pdims(ddim, ORTHO2), LO), cg%lhn(pdims(ddim,ORTHO2), HI)
       do i1 = cg%lhn(pdims(ddim, ORTHO1), LO), cg%lhn(pdims(ddim, ORTHO1), HI)
          pu => cg%w(wna%fi)%get_sweep(ddim,i1,i2)
          u1d(iarr_all_swp(ddim,:),:) = pu(:,:)
          pb => cg%w(bi)%get_sweep(ddim,i1,i2)
          b_cc1d(iarr_mag_swp(ddim,:),:) = pb(:,:)
#ifdef COSM_RAYS
          call set_div_v1d(div_v1d, ddim, i1, i2, cg)
#endif /* COSM_RAYS */
          if (psii /= INVALID) then
             ppsi => cg%q(psii)%get_sweep(ddim,i1,i2)
             psi_d(1, :) = ppsi(:)
             call solve(u1d, b_cc1d, dt/cg%dl(ddim), psi_d, div_v1d)
             ppsi(:) = psi_d(1,:)
          else
             call solve(u1d, b_cc1d, dt/cg%dl(ddim), psi_d, div_v1d)
          endif
          pu(:,:) = u1d(iarr_all_swp(ddim,:),:)
          pb(:,:) = b_cc1d(iarr_mag_swp(ddim,:),:) ! ToDo figure out how to manage CT energy fixup without extra storage
       enddo
    enddo

  end subroutine sweep_dsplit

!---------------------------------------------------------------------------------------------------------------------

  subroutine solve(u, b_cc, dtodx, psi, div_v1d)

     use constants,  only: half
     use dataio_pub, only: die
     use global,     only: h_solver
     use interpolations, only: interpol

     implicit none

     real, dimension(:,:), intent(inout) :: u
     real, dimension(:,:), intent(inout) :: b_cc
     real,                 intent(in)    :: dtodx
     real, dimension(:,:), intent(inout) :: psi
     real, dimension(:), pointer, intent(in) :: div_v1d

     real, dimension(size(b_cc,1),size(b_cc,2)), target :: b_cc_l, b_cc_r, mag_cc
     real, dimension(size(b_cc,1),size(b_cc,2))         :: bclflx, bcrflx, db1, db2, db3
     real, dimension(size(u,1),size(u,2)), target       :: flx, ql, qr
     real, dimension(size(u,1),size(u,2))               :: flx_l, flx_r
     real, dimension(size(u,1),size(u,2))               :: du1, du2, du3

     real, dimension(size(psi,1),size(psi,2))           :: psilflx, psirflx, dpsi1, dpsi2, dpsi3
     real, dimension(size(psi,1),size(psi,2)), target   :: psi_l, psi_r
     real, dimension(size(psi,1),size(psi,2)),target    :: psi_cc


     integer                                            :: nx

     nx  = size(u,2)
     if (size(b_cc,2) /= nx) call die("[fluidupdate:rk2] size b_cc and u mismatch")
     mag_cc = huge(1.)

     ! Only muscl and rk2 schemes should be considered for production use.
     ! Other schemes are left here for educational purposes, just to show how to construct alternative approaches.
     select case (h_solver)
     case ("rk2")
        call interpol(u,b_cc,psi,ql,qr,b_cc_l,b_cc_r,psi_l,psi_r)
        call riemann_wrap                   ! Now we advance the left and right states by a timestep.
        call du_db(du1, db1,dpsi1)
        call interpol(u+half*du1,b_cc+half*db1,psi+half*dpsi1,ql,qr,b_cc_l,b_cc_r,psi_l,psi_r)
        call riemann_wrap                   ! second call for Riemann problem uses states evolved to half timestep
        call update
     case ("muscl")
        call interpol(u,b_cc,psi,ql,qr,b_cc_l,b_cc_r,psi_l,psi_r)
        call musclflx(nx, ql, b_cc_l, psi_l, flx_l, bclflx, psilflx)
        call musclflx(nx, qr, b_cc_r, psi_r, flx_r, bcrflx, psirflx)
        call ulr_fluxes_qlr
        call riemann_wrap
        call update
     case default
        call die("[fluidupdate:sweep_dsplit] No recognized solver")
     end select

   contains

        ! some shortcuts

     subroutine du_db(du, db, dpsi)

       use constants,  only: DIVB_HDC
       use global,     only: divB_0_method

       implicit none

       real, dimension(size(u,1),size(u,2)),       intent(out) :: du
       real, dimension(size(b_cc,1),size(b_cc,2)), intent(out) :: db
       real, dimension(size(psi,1),size(psi,2)),   intent(out) :: dpsi

       du(:,2:nx) = dtodx*(flx(:,1:nx-1) - flx(:,2:nx))
       du(:,1) = du(:,2)

       db(:,2:nx) = dtodx*(mag_cc(:,1:nx-1) - mag_cc(:,2:nx))
       db(:,1) = db(:,2)

       if (divB_0_method == DIVB_HDC) then
          dpsi(:,2:nx) = dtodx*(psi_cc(:,1:nx-1) - psi_cc(:,2:nx))
          dpsi(:,1) = dpsi(:,2)
       else
          dpsi = 0.
       endif

     end subroutine du_db

     subroutine ulr_fluxes_qlr

        use constants, only: DIVB_HDC
        use global,    only: divB_0_method

        implicit none

        call addflux(ql, flx_l, half * dtodx, nx)
        call addflux(qr, flx_r, half * dtodx, nx)
        call addflux(b_cc_l, bclflx, half * dtodx, nx)
        call addflux(b_cc_r, bcrflx, half * dtodx, nx)

        if (divB_0_method == DIVB_HDC) then
           call addflux(psi_l, psilflx, half * dtodx, nx)
           call addflux(psi_r, psirflx, half * dtodx, nx)
        endif

     end subroutine ulr_fluxes_qlr

     subroutine addflux(a, fa, dt, nx)

        implicit none

        real, dimension(:,:), intent(inout) :: a
        real, dimension(:,:), intent(in)    :: fa
        real                                :: dt
        integer                             :: nx

        a(:, 2:nx) = a(:, 2:nx) + dt * (fa(:, 1:nx-1) - fa(:, 2:nx))
        a(:, 1) = a(:, 2)

     end subroutine addflux

     subroutine musclflx(n, q, b_cc, psi, qf, b_ccf, psif)

        use constants,  only: half, xdim, ydim, zdim, DIVB_HDC, zero
        use fluidindex, only: flind
        use fluidtypes, only: component_fluid
        use global,     only: divB_0_method
        use hdc,        only: chspeed

        implicit none

        integer,              intent(in)  :: n
        real, dimension(:,:), intent(in)  :: q
        real, dimension(:,:), intent(in)  :: b_cc
        real, dimension(:,:), intent(in)  :: psi
        real, dimension(:,:), intent(out) :: qf
        real, dimension(:,:), intent(out) :: b_ccf, psif

        class(component_fluid), pointer   :: fl
        real                              :: en
        integer                           :: ip, i

        qf    = zero
        b_ccf = zero
        psif  = zero

        do ip = 1, flind%fluids

           fl => flind%all_fluids(ip)%fl

           do i = 1, n

              qf(fl%idn,i) = q(fl%idn,i)*q(fl%imx,i)
              if (fl%has_energy) then
                 if (fl%is_magnetized) then
                    qf(fl%imx,i) = q(fl%idn,i)*q(fl%imx,i)*q(fl%imx,i) + (q(fl%ien,i) + half*sum(b_cc(xdim:zdim,i)**2,dim=1)) - b_cc(xdim,i)**2
                 else
                    qf(fl%imx,i) = q(fl%idn,i)*q(fl%imx,i)*q(fl%imx,i) + q(fl%ien,i)
                 endif
              else
                 qf(fl%imx,i) = q(fl%idn,i)*q(fl%imx,i)*q(fl%imx,i)
              endif
              if (fl%is_magnetized) then
                 qf(fl%imy,i) = q(fl%idn,i)*q(fl%imx,i)*q(fl%imy,i) - b_cc(xdim,i)*b_cc(ydim,i)
                 qf(fl%imz,i) = q(fl%idn,i)*q(fl%imx,i)*q(fl%imz,i) - b_cc(xdim,i)*b_cc(zdim,i)
                 b_ccf(ydim,i)  = b_cc(ydim,i)*q(fl%imx,i) - b_cc(xdim,i)*q(fl%imy,i)
                 b_ccf(zdim,i)  = b_cc(zdim,i)*q(fl%imx,i) - b_cc(xdim,i)*q(fl%imz,i)
                 if (divB_0_method .eq. DIVB_HDC) then
                    b_ccf(xdim,i) = psi(1,i)
                    psif(1,i)   = (chspeed**2)*b_cc(xdim,i)
                 endif
              else
                 qf(fl%imy,i) = q(fl%idn,i)*q(fl%imx,i)*q(fl%imy,i)
                 qf(fl%imz,i) = q(fl%idn,i)*q(fl%imx,i)*q(fl%imz,i)
              endif
              if (fl%has_energy) then
                 if (fl%is_magnetized) then
                    en = (q(fl%ien,i)/(fl%gam_1)) + half*q(fl%idn,i)*sum(q(fl%imx:fl%imz,i)**2) + half*sum(b_cc(xdim:zdim,i)**2)
                    qf(fl%ien,i) = (en + (q(fl%ien,i) + half*sum(b_cc(xdim:zdim,i)**2,dim=1)))*q(fl%imx,i) - b_cc(xdim,i)*dot_product(q(fl%imx:fl%imz,i),b_cc(xdim:zdim,i))
                 else
                    en = (q(fl%ien,i)/(fl%gam_1)) + half*q(fl%idn,i)*sum(q(fl%imx:fl%imz,i)**2)
                    qf(fl%ien,i) = (en + (q(fl%ien,i)))*q(fl%imx,i)
                 endif
              endif

           enddo

        enddo

     end subroutine musclflx

     subroutine riemann_wrap()

       use constants,  only: xdim, zdim, DIVB_HDC
       use fluidindex, only: flind
       use fluidtypes, only: component_fluid
       use global,     only: divB_0_method
       use hlld,       only: riemann_hlld

       implicit none

       integer :: i
       class(component_fluid), pointer :: fl
       real, dimension(size(b_cc,1),size(b_cc,2)), target :: b0, bf0
       real, dimension(:,:), pointer :: p_flx, p_bcc, p_bccl, p_bccr, p_ql, p_qr
       real, dimension(size(psi,1),size(psi,2)), target ::  p0, pf0
       real, dimension(:,:), pointer :: p_psif, p_psi_l, p_psi_r

       do i = 1, flind%fluids
          fl    => flind%all_fluids(i)%fl
          p_flx => flx(fl%beg:fl%end,:)
          p_ql  => ql(fl%beg:fl%end,:)
          p_qr  => qr(fl%beg:fl%end,:)
          if (fl%is_magnetized) then
             p_bccl => b_cc_l(xdim:zdim,:)
             p_bccr => b_cc_r(xdim:zdim,:)
             p_bcc  => mag_cc(xdim:zdim,:)
             if (divB_0_method == DIVB_HDC) then
                p_psi_l => psi_l(:,:)
                p_psi_r => psi_r(:,:)
                p_psif  => psi_cc(:,:)
             else  ! CT
                p0 = 0.
                p_psi_l => p0
                p_psi_r => p0
                p_psif  => pf0
            endif
          else ! ignore all magnetic field
             b0 = 0.
             p_bccl => b0
             p_bccr => b0
             p_bcc  => bf0
             p0 = 0.
             p_psi_l => p0
             p_psi_r => p0
             p_psif  => pf0
          endif

          call riemann_hlld(nx, p_flx, p_ql, p_qr, p_bcc, p_bccl, p_bccr, p_psi_l, p_psi_r, p_psif, fl%gam) ! whole mag_cc is not needed now for simple schemes but rk2 and rk4 still rely on it
       enddo
     end subroutine riemann_wrap

     subroutine update(weights)

       use constants,        only: xdim, ydim, zdim, DIVB_HDC
       use fluidindex,       only: flind
       use global,           only: divB_0_method

#ifdef COSM_RAYS
       use fluidindex,       only: iarr_all_dn, iarr_all_mx, iarr_all_en
       use global,           only: dt
       use initcosmicrays,   only: iarr_crs, smallecr
       use sourcecosmicrays, only: src_gpcr
#ifdef COSM_RAYS_SOURCES
       use initcosmicrays,   only: iarr_crn
       use sourcecosmicrays, only: src_crn
#endif /* COSM_RAYS_SOURCES */
#endif /* COSM_RAYS */

       implicit none

       real, optional, dimension(:), intent(in) :: weights

       real, dimension(:), allocatable :: w
       integer :: iend  !< last component of any fluid (i.e. exclude CR or tracers here)

#ifdef COSM_RAYS
       real, dimension(size(u,2),size(u,1))           :: u1
       real, dimension(nx, flind%fluids), target      :: vx
       real, dimension(nx)                            :: grad_pcr
       real, dimension(nx, flind%crs%all)             :: decr
#ifdef COSM_RAYS_SOURCES
       real, dimension(nx, flind%crn%all)             :: srccrn
#endif /* COSM_RAYS_SOURCES */
#endif /* COSM_RAYS */

       iend = flind%all_fluids(flind%fluids)%fl%end

       if (present(weights)) then
          allocate(w(size(weights)))
          w = weights/sum(weights)
       else
          allocate(w(1))
          w(1) = 1.
       endif

       u(:iend,2:nx) = u(:iend,2:nx) + w(1) * dtodx * (flx(:iend,1:nx-1) - flx(:iend,2:nx))
       if (size(w)>=2) u(:iend,2:nx) = u(:iend,2:nx) + w(2) * du1(:iend,2:nx)
       if (size(w)>=3) u(:iend,2:nx) = u(:iend,2:nx) + w(3) * du2(:iend,2:nx)
       if (size(w)>=4) u(:iend,2:nx) = u(:iend,2:nx) + w(4) * du3(:iend,2:nx)
       u(:iend,1) = u(:iend,2)
       u(:iend,nx) = u(:iend,nx-1)

       b_cc(ydim:zdim,2:nx) = b_cc(ydim:zdim,2:nx) + w(1) * dtodx * (mag_cc(ydim:zdim,1:nx-1) - mag_cc(ydim:zdim,2:nx))
       if (size(w)>=2)  b_cc(ydim:zdim,2:nx) = b_cc(ydim:zdim,2:nx) + w(2) * db1(ydim:zdim,2:nx)
       if (size(w)>=3)  b_cc(ydim:zdim,2:nx) = b_cc(ydim:zdim,2:nx) + w(3) * db2(ydim:zdim,2:nx)
       if (size(w)>=4)  b_cc(ydim:zdim,2:nx) = b_cc(ydim:zdim,2:nx) + w(4) * db3(ydim:zdim,2:nx)

       if (divB_0_method == DIVB_HDC) then

          b_cc(xdim, 2:nx) = b_cc(xdim, 2:nx) + w(1) * dtodx * (mag_cc(xdim, 1:nx-1) - mag_cc(xdim, 2:nx))
          if (size(w)>=2)  b_cc(xdim,2:nx) = b_cc(xdim,2:nx) + w(2) * db1(xdim,2:nx)
          if (size(w)>=3)  b_cc(xdim,2:nx) = b_cc(xdim,2:nx) + w(3) * db2(xdim,2:nx)
          if (size(w)>=4)  b_cc(xdim,2:nx) = b_cc(xdim,2:nx) + w(4) * db3(xdim,2:nx)

          psi(1, 2:nx) = psi(1, 2:nx) + w(1) * dtodx * (psi_cc(1, 1:nx-1) - psi_cc(1, 2:nx))
          if (size(w)>=2)  psi_cc(1,2:nx) = psi_cc(1,2:nx) + w(2) * dpsi1(1,2:nx)
          if (size(w)>=3)  psi_cc(1,2:nx) = psi_cc(1,2:nx) + w(3) * dpsi2(1,2:nx)
          if (size(w)>=4)  psi_cc(1,2:nx) = psi_cc(1,2:nx) + w(4) * dpsi3(1,2:nx)
          psi(:,1) = psi(:,2)
          psi(:,nx) = psi(:, nx-1)

          !damping
          !psi = psi*exp(-glm_alpha*chspeed*dtodx)

       endif

       b_cc(:, 1)  = b_cc(:, 2)
       b_cc(:, nx) = b_cc(:, nx-1)

       ! This is lowest order implementation of CR
       ! It agrees with the implementation in RTVD in the limit of small CR energy amounts
       ! ToDo: integrate it into h_solver schemes
#if defined COSM_RAYS && defined IONIZED

       if (nx > 1) then
          ! transposition for compatibility with RTVD-based routines
          u1 = transpose(u)

          vx = u1(:, iarr_all_mx) / u1(:, iarr_all_dn) ! this may also be useful for gravitational acceleration
          ! Replace dt/dtodx by dx == cg%dl(ddim)
          call src_gpcr(u1, nx, dt/dtodx, div_v1d, decr, grad_pcr)
          u1(:, iarr_crs(:)) = u1(:, iarr_crs(:)) + decr(:,:) * dt
          u1(:, iarr_crs(:)) = max(smallecr, u1(:, iarr_crs(:)))
          u1(:, iarr_all_mx(flind%ion%pos)) = u1(:, iarr_all_mx(flind%ion%pos)) + grad_pcr * dt
#ifndef ISO
          u1(:, iarr_all_en(flind%ion%pos)) = u1(:, iarr_all_en(flind%ion%pos)) + vx(:, flind%ion%pos) * grad_pcr * dt
#endif /* !ISO */
       endif
#ifdef COSM_RAYS_SOURCES
       call src_crn(u1, nx, srccrn, dt) ! n safe
       u1(:, iarr_crn) = u1(:, iarr_crn) + srccrn(:,:)*dt
#endif /* COSM_RAYS_SOURCES */

       u = transpose(u1)

#else
       if (.false.) div_v1d = div_v1d + 0.  ! suppress compiler warnings
#endif /* COSM_RAYS && IONIZED */

       deallocate(w)

     end subroutine update

   end subroutine solve

 end module fluidupdate
