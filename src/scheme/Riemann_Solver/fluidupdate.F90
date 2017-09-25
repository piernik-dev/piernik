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

  subroutine fluid_update

    use all_boundaries, only: all_bnd
    use cg_leaves,      only: leaves
    use cg_list,        only: cg_list_element
    use constants,      only: xdim, zdim, GEO_XYZ
    use dataio_pub,     only: halfstep, die, warn
    use domain,         only: dom, is_refined
    use fluxlimiters,   only: set_limiters
    use global,         only: skip_sweep, dt, dtm, t, limiter, limiter_b
    use mpisetup,       only: master
    use user_hooks,     only: problem_customize_solution

    implicit none

    logical, save                   :: first_run = .true.
    type(cg_list_element), pointer  :: cgl
    integer(kind=4)                 :: ddim

    ! is_multicg should be safe
    if (is_refined) call die("[fluid_update] This Rieman solver is not compatible with mesh refinements yet!")
    if (dom%geometry_type /= GEO_XYZ) call die("[fluid_update] Non-cartesian geometry is not implemented yet in this Riemann solver.")
#ifdef GRAV
#  error Graviy is not implemented yet in this Riemann solver.
#endif /* GRAV */
#ifdef ISO
#  error Isothermal EOS is not implemented yet in this Riemann solver.
#endif /* ISO */
#ifdef COSM_RAYS
#  error   Isothermal EOS is notCosmic rays are not implemented yet in this Riemann solver.
#endif /* COSM_RAYS */

    halfstep = .false.
    if (first_run) then
       dtm = 0.0
       call set_limiters(limiter, limiter_b)
       if (master) call warn("[fluid_update] This is an experimental implementation of the Riemann solver. Expect unexpected.")
    else
       dtm = dt
    endif
    t = t + dt

    ! ToDo: check if changes of execution order here (block loop, direction loop, boundary update can change
    ! cost or allow for reduction of required guardcells
    call bfc2bcc
    do ddim = xdim, zdim, 1
      if (dom%has_dir(ddim)) then

         cgl => leaves%first
         do while (associated(cgl))

            ! Warning: 2.5D MHD may need all directional calls anyway
            if (.not. skip_sweep(ddim)) call sweep_dsplit(cgl%cg,dt,ddim)
            if (associated(problem_customize_solution)) call problem_customize_solution(.true.)
            cgl => cgl%nxt
         enddo
         call all_bnd
#ifdef MAGNETIC
         call magfield(ddim)
#endif /* MAGNETIC */
      endif
    enddo

    t = t + dt
    dtm = dt
    halfstep = .true.

    call bfc2bcc
    do ddim = zdim, xdim, -1
       if (dom%has_dir(ddim)) then
#ifdef MAGNETIC
          call magfield(ddim)
#endif /* MAGNETIC */
          cgl => leaves%first
          do while (associated(cgl))

             if (.not. skip_sweep(ddim)) call sweep_dsplit(cgl%cg,dt,ddim)
             if (associated(problem_customize_solution)) call problem_customize_solution(.false.)
             cgl => cgl%nxt
          enddo
          call all_bnd
       endif
    enddo

    if (first_run) first_run = .false.

  end subroutine fluid_update

!-------------------------------------------------------------------------------------------------------------------

#ifdef MAGNETIC
   subroutine magfield(dir)

      use ct,     only: advectb, ctb
      use constants,   only: ndims, I_ONE
#ifdef RESISTIVE
      use resistivity, only: diffuseb
#endif /* RESISTIVE */

       implicit none

       integer(kind=4), intent(in) :: dir

       integer(kind=4)             :: bdir, dstep

       do dstep = 0, 1
          bdir  = I_ONE + mod(dir+dstep,ndims)
          call advectb(bdir, dir)
#ifdef RESISTIVE
         call diffuseb(bdir, dir)
#endif /* RESISTIVE */
         call mag_add(dir, bdir)
      enddo

      ! call energy_fixup use cfl_resist = 0.3

   end subroutine magfield

!-------------------------------------------------------------------------------------------------------------------

 subroutine mag_add(dim1, dim2)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
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
     use domain,           only: dom
     use grid_cont,        only: grid_container
     use named_array_list, only: wna

     implicit none

     type(cg_list_element), pointer :: cgl
     type(grid_container),  pointer :: cg

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

  subroutine energy_fixup

     use all_boundaries,   only: all_bnd
     use cg_leaves,        only: leaves
     use cg_list,          only: cg_list_element
     use constants,        only: xdim, ydim, zdim, LO, HI, half
     use domain,           only: dom
     use fluidindex,       only: flind
     use func,             only: emag
     use grid_cont,        only: grid_container
     use named_array_list, only: wna

     implicit none

     type(cg_list_element), pointer :: cgl
     type(grid_container),  pointer :: cg

     integer :: i, j, k

     cgl => leaves%first
     do while (associated(cgl))
        cg => cgl%cg

        do k = cg%lh1(zdim, LO), cg%lh1(zdim, HI)
           do j = cg%lh1(ydim, LO), cg%lh1(ydim, HI)-1
              do i = cg%lh1(xdim, LO), cg%lh1(xdim, HI)-1
                 cg%u(flind%ion%ien,i,j,k) = cg%u(flind%ion%ien,i,j,k) + ( &
                      emag(half*(cg%b(xdim, i, j, k) + cg%b(xdim, i+dom%D_x, j, k)), &
                      &    half*(cg%b(ydim, i, j, k) + cg%b(ydim, i, j+dom%D_y, k)), &
                      &    half*(cg%b(zdim, i, j, k) + cg%b(zdim, i, j, k+dom%D_z)) ) - &
                      emag(cg%w(wna%bcci)%arr(xdim, i, j, k), &
                      &    cg%w(wna%bcci)%arr(ydim, i, j, k), &
                      &    cg%w(wna%bcci)%arr(zdim, i, j, k) ) ) /4.
                 ! 1/8. to 1/6. seems to give best survivability of the otvortex but it is still far from being good
              enddo
           enddo
        enddo
        cgl => cgl%nxt
     enddo

     call all_bnd ! overkill

  end subroutine energy_fixup

  subroutine sweep_dsplit(cg, dt, ddim)

    use constants,        only: pdims, xdim, zdim, ORTHO1, ORTHO2, LO, HI
    use fluidindex,       only: iarr_all_swp, iarr_mag_swp
    use grid_cont,        only: grid_container
    use named_array_list, only: wna

    implicit none

    type(grid_container), pointer, intent(in) :: cg
    real,                          intent(in) :: dt
    integer(kind=4),               intent(in) :: ddim

    real, dimension(size(cg%u,1), cg%n_(ddim)) :: u1d
    real, dimension(xdim:zdim, cg%n_(ddim))   :: b_cc1d
    real, dimension(:,:), pointer             :: pu, pb
    integer                                   :: i1, i2

    do i2 = cg%lhn(pdims(ddim, ORTHO2), LO), cg%lhn(pdims(ddim,ORTHO2), HI)
       do i1 = cg%lhn(pdims(ddim, ORTHO1), LO), cg%lhn(pdims(ddim, ORTHO1), HI)
          pu => cg%w(wna%fi)%get_sweep(ddim,i1,i2)
          u1d(iarr_all_swp(ddim,:),:) = pu(:,:)
          pb => cg%w(wna%bcci)%get_sweep(ddim,i1,i2)
          b_cc1d(iarr_mag_swp(ddim,:),:) = pb(:,:)
          call solve(u1d,b_cc1d, dt/cg%dl(ddim))
          pu(:,:) = u1d(iarr_all_swp(ddim,:),:)
          pb(:,:) = b_cc1d(iarr_mag_swp(ddim,:),:) ! ToDo figure out how to manage CT energy fixup without extra storage
       enddo
    enddo

  end subroutine sweep_dsplit

!---------------------------------------------------------------------------------------------------------------------

  subroutine solve(u, b_cc, dtodx)

     use constants,  only: half
     use dataio_pub, only: die
     use global,     only: h_solver

     implicit none

     real, dimension(:,:), intent(inout) :: u
     real, dimension(:,:), intent(inout) :: b_cc
     real,                 intent(in)    :: dtodx

     real, dimension(size(b_cc,1),size(b_cc,2)), target :: b_cc_l, b_cc_r, mag_cc
     real, dimension(size(b_cc,1),size(b_cc,2))         :: b_ccl, b_ccr, db1, db2, db3
     real, dimension(size(u,1),size(u,2)), target       :: flx, ql, qr
     real, dimension(size(u,1),size(u,2))               :: ul, ur, du1, du2, du3
     integer                                            :: nx

     nx  = size(u,2)
     if (size(b_cc,2) /= nx) call die("[fluidupdate:rk2] size b_cc and u mismatch")
     mag_cc = huge(1.)

     ! Only muscl and rk2 schemes should be considered for production use.
     ! Other schemes are left here for educational purposes, just to show how to construct alternative approaches.
     select case (h_solver)
        case ("muscl")
           call slope
           call ulr_fluxes_qlr
           call riemann_wrap
           call update
        case ("rk2")
           call slope
           call ulr_to_qlr                     ! Just the slope is used to feed 1st call to Riemann solver
           call riemann_wrap                   ! Now we advance the left and right states by half timestep.
           call du_db(du1, db1)                ! The slope is already calculated and can be reused
           call ulr_to_qlr(half*du1, half*db1)
           call riemann_wrap                   ! second call for Riemann problem uses states evolved to half timestep
           call update
        case ("rk2_s")                         ! RK2 with alternative approach to calculating slopes for 2nd step
           call slope
           call ulr_to_qlr
           call riemann_wrap
           call du_db(du1, db1)
           call slope(half*du1, half*db1)
           call ulr_to_qlr
           call riemann_wrap
           call update
        case ("rk2_muscl")
           call slope
           call ulr_fluxes_qlr
           call riemann_wrap                   ! MUSCL-Hancock is used to feed 1st call to Riemann solver
           call du_db(du1, db1)                ! Now we can calculate state for half-timestep and recalculate slopes
           call ulr_to_qlr(half*du1, half*db1)
           call riemann_wrap                   ! second call for Riemann problem needs just the slope from states evolved to half timestep
           call update
        case ("rk2_muscl_s")                   ! MUSCL-RK2 with alternative approach to calculating slopes for 2nd step
           call slope
           call ulr_fluxes_qlr
           call riemann_wrap
           call du_db(du1, db1)
           call slope(half*du1, half*db1)
           call ulr_to_qlr
           call riemann_wrap
           call update
        case ("euler")                         ! Gives quite sharp advection, especially for CFL=0.5, but it is unstable and produces a lot of noise. Do not use, except for educational purposes.
           call slope
           call ulr_to_qlr
           call riemann_wrap
           call update
        case ("heun")
           call slope
           call ulr_to_qlr
           call riemann_wrap
           call du_db(du1, db1)
           call ulr_to_qlr(du1, db1)
           call riemann_wrap
           call update([1., 1.])
        case ("rk4")
           call slope
           call ulr_to_qlr
           call riemann_wrap
           call du_db(du1, db1)
           call ulr_to_qlr(half*du1, half*db1)
           call riemann_wrap
           call du_db(du2, db2)
           call ulr_to_qlr(half*du2, half*db2)
           call riemann_wrap
           call du_db(du3, db3)
           call ulr_to_qlr(du3, db3)
           call riemann_wrap
           call update([1., 1., 2., 2.])
        case ("rk4_s")                         ! RK4 with alternative approach to calculating slopes for 2nd-4th step
           call slope
           call ulr_to_qlr
           call riemann_wrap
           call du_db(du1, db1)
           call slope(half*du1, half*db1)
           call ulr_to_qlr
           call riemann_wrap
           call du_db(du2, db2)
           call slope(half*du2, half*db2)
           call ulr_to_qlr
           call riemann_wrap
           call du_db(du3, db3)
           call slope(du3, db3)
           call ulr_to_qlr
           call riemann_wrap
           call update([1., 1., 2., 2.])
        case default
           call die("[fluidupdate:sweep_dsplit] No recognized solver")
     end select

     contains

        ! some shortcuts

        subroutine slope(uu, bb)

           use constants,  only: half
           use dataio_pub, only: die
           use fluxlimiters, only: flimiter, blimiter

           implicit none

           real, optional, dimension(size(u,1),size(u,2)),       intent(in) :: uu
           real, optional, dimension(size(b_cc,1),size(b_cc,2)), intent(in) :: bb

           real, dimension(size(u,1),size(u,2))       :: du
           real, dimension(size(b_cc,1),size(b_cc,2)) :: db

           if (present(uu) .neqv. present(bb)) call die("[fluidupdate:solve:slope] either none or both optional arguments must be present")

           if (present(uu)) then
              du  = flimiter(u + uu)
              ul  = u + uu - half*du
              ur  = u + uu + half*du
           else
              du  = flimiter(u)
              ul  = u - half*du
              ur  = u + half*du
           endif

           if (present(bb)) then
              db  = blimiter(b_cc + bb)
              b_ccl = b_cc + bb - half*db
              b_ccr = b_cc + bb + half*db
           else
              db  = blimiter(b_cc)
              b_ccl = b_cc - half*db
              b_ccr = b_cc + half*db
           endif

        end subroutine slope

        subroutine ulr_to_qlr(du, db)

           use dataio_pub, only: die

           implicit none

           real, optional, dimension(size(u,1),size(u,2)),       intent(in) :: du
           real, optional, dimension(size(b_cc,1),size(b_cc,2)), intent(in) :: db

           real, dimension(size(u,1),size(u,2))               :: u_l, u_r

           if (present(du) .neqv. present(db)) call die("[fluidupdate:solve:ulr_to_qlr] either mone or both optional arguments must be present")

           if (present(du)) then
              u_l = ur + du
              u_r(:,1:nx-1) = ul(:,2:nx) + du(:,2:nx)
           else
              u_l = ur
              u_r(:,1:nx-1) = ul(:,2:nx)
           endif
           u_r(:,nx) = u_r(:,nx-1)

           if (present(db)) then
              b_cc_l = b_ccr + db
              b_cc_r(:,1:nx-1) = b_ccl(:,2:nx) + db(:,2:nx)
           else
              b_cc_l = b_ccr
              b_cc_r(:,1:nx-1) = b_ccl(:,2:nx)
           endif
           b_cc_r(:,nx) = b_cc_r(:,nx-1)

           ql = utoq(u_l,b_cc_l)
           qr = utoq(u_r,b_cc_r)

        end subroutine ulr_to_qlr

        subroutine ulr_fluxes_qlr

           use hlld,       only: fluxes

           implicit none

           real, dimension(size(u,1)+size(b_cc,1),size(u,2)), target :: flx
           real, dimension(size(u,1),size(u,2))                      :: u_l, u_r

           flx  = fluxes(ul,b_ccl) - fluxes(ur,b_ccr)

           u_l = ur + half*dtodx*flx(:size(u,1),:)
           u_r(:,1:nx-1) = ul(:,2:nx) + half*dtodx*flx(:size(u,1),2:nx)
           u_r(:,nx) = u_r(:,nx-1)

           b_cc_l(:,2:nx) = b_ccr(:,2:nx) + half*dtodx*flx(size(u,1)+1:,2:nx)
           b_cc_l(:,1) = b_cc_l(:,2)
           b_cc_r(:,1:nx-1) = b_ccl(:,2:nx) + half*dtodx*flx(size(u,1)+1:,2:nx)
           b_cc_r(:,nx) = b_cc_r(:,nx-1)

           ql = utoq(u_l,b_cc_l)
           qr = utoq(u_r,b_cc_r)

        end subroutine ulr_fluxes_qlr

        subroutine du_db(du, db)

           implicit none

           real, dimension(size(u,1),size(u,2)),       intent(out) :: du
           real, dimension(size(b_cc,1),size(b_cc,2)), intent(out) :: db

           du(:,2:nx) = dtodx*(flx(:,1:nx-1) - flx(:,2:nx))
           du(:,1) = du(:,2)

           db(:,2:nx) = dtodx*(mag_cc(:,1:nx-1) - mag_cc(:,2:nx))
           db(:,1) = db(:,2)

        end subroutine du_db

        subroutine riemann_wrap()

           use constants,  only: xdim, zdim
           use fluidindex, only: flind
           use fluidtypes, only: component_fluid
           use hlld,       only: riemann_hlld

           implicit none

           integer :: i
           class(component_fluid), pointer :: fl
           real, dimension(size(b_cc,1),size(b_cc,2)), target :: b0
           real, dimension(:,:), pointer :: p_flx, p_bcc, p_bccl, p_bccr, p_ql, p_qr

           do i = 1, flind%fluids
              fl    => flind%all_fluids(i)%fl
              p_flx => flx(fl%beg:fl%end,:)
              p_ql  => ql(fl%beg:fl%end,:)
              p_qr  => qr(fl%beg:fl%end,:)
              p_bcc => mag_cc(xdim:zdim,:)
              if (fl%is_magnetized) then
                 p_bccl => b_cc_l(xdim:zdim,:)
                 p_bccr => b_cc_r(xdim:zdim,:)
              else ! ignore all magnetic field
                 b0 = 0.
                 p_bccl => b0
                 p_bccr => b0
              endif
              call riemann_hlld(nx, p_flx, p_ql, p_qr, p_bcc, p_bccl, p_bccr, fl%gam) ! whole mag_cc is not needed now for simple schemes but rk2 and rk4 still rely on it
           enddo

        end subroutine riemann_wrap

        subroutine update(weights)

           implicit none

           real, optional, dimension(:), intent(in) :: weights

           real, dimension(:), allocatable :: w

           if (present(weights)) then
              allocate(w(size(weights)))
              w = weights/sum(weights)
           else
              allocate(w(1))
              w(1) = 1.
           endif

           u(:,2:nx) = u(:,2:nx) + w(1) * dtodx * (flx(:,1:nx-1) - flx(:,2:nx))
           if (size(w)>=2) u(:,2:nx) = u(:,2:nx) + w(2) * du1(:,2:nx)
           if (size(w)>=3) u(:,2:nx) = u(:,2:nx) + w(3) * du2(:,2:nx)
           if (size(w)>=4) u(:,2:nx) = u(:,2:nx) + w(4) * du3(:,2:nx)
           u(:,1) = u(:,2)
           u(:,nx) = u(:,nx-1)

           b_cc(:,2:nx) = b_cc(:,2:nx) + w(1) * dtodx * (mag_cc(:,1:nx-1) - mag_cc(:,2:nx))
           if (size(w)>=2)  b_cc(:,2:nx) = b_cc(:,2:nx) + w(2) * db1(:,2:nx)
           if (size(w)>=3)  b_cc(:,2:nx) = b_cc(:,2:nx) + w(3) * db2(:,2:nx)
           if (size(w)>=4)  b_cc(:,2:nx) = b_cc(:,2:nx) + w(4) * db3(:,2:nx)
           b_cc(:,1) = b_cc(:,2)
           b_cc(:,nx) = b_cc(:,nx-1)

           deallocate(w)

        end subroutine update

  end subroutine solve

!--------------------------------------------------------------------------------------------------------------

   function utoq(u,b_cc) result(q)

     use constants,  only: half, xdim, zdim
     use fluidindex, only: flind
     use fluidtypes, only: component_fluid
     use func,       only: ekin

     implicit none

     real, dimension(:,:),   intent(in)    :: u , b_cc

     real, dimension(size(u,1),size(u,2))  :: q
     integer  :: p

     class(component_fluid), pointer       :: fl

     do p = 1, flind%fluids
        fl => flind%all_fluids(p)%fl

        q(fl%idn,:) =  u(fl%idn,:)
        q(fl%imx,:) =  u(fl%imx,:)/u(fl%idn,:)
        q(fl%imy,:) =  u(fl%imy,:)/u(fl%idn,:)
        q(fl%imz,:) =  u(fl%imz,:)/u(fl%idn,:)
        ! J.CoPhy 208 (2005),Pg 317, Eq. 2. Gas pressure: p = (gamma-1)*(e-half*rho*v^2-half*B^2) and Total pressure: p_T = p + half*B^2. (1) and (2) are markers for HD and MHD.
        if (fl%has_energy) then
            q(fl%ien,:) =  fl%gam_1*(u(fl%ien,:) - ekin(u(fl%imx,:), u(fl%imy,:), u(fl%imz,:), u(fl%idn,:))) ! Primitive variable for gas pressure (p) without magnetic fields. (1)
            if (fl%is_magnetized) then
               q(fl%ien,:) =  q(fl%ien,:) - half*fl%gam_1*sum(b_cc(xdim:zdim,:)**2, dim=1) ! Primitive variable for gas pressure (p) with magnetic fields. The requirement of total pressure is dealt in the fluxes and hlld routines. (2)
            endif
        endif

     enddo

   end function utoq

end module fluidupdate
