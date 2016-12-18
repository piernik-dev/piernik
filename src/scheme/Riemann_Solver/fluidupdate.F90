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
!--------------------------------------------------------------------------------------------------------------

#include "piernik.def"

module fluidupdate
! pulled by RIEMANN

  implicit none
  private
  public :: fluid_update

contains

  subroutine fluid_update

    use cg_list,        only: cg_list_element
    use cg_leaves,      only: leaves
    use constants,      only: xdim, zdim
    use domain,         only: dom
    use all_boundaries, only: all_bnd
    use global,         only: skip_sweep, dt, dtm, t
    use user_hooks,     only: problem_customize_solution
    use dataio_pub,     only: halfstep

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

    ! ToDo: check if changes of execution order here (block loop, direction loop, boundary update can change
    ! cost or allow for reduction of required guardcells
    cgl => leaves%first
    do while (associated(cgl))

       do ddim = xdim, zdim, 1
          ! Warning: 2.5D MHD may need all directional calls anyway
          if (.not. skip_sweep(ddim) .and. dom%has_dir(ddim)) call sweep_dsplit(cgl%cg,dt,ddim)
       enddo
       if (associated(problem_customize_solution)) call problem_customize_solution(.true.)
       cgl => cgl%nxt
    enddo
    call all_bnd

    t = t + dt
    dtm = dt
    halfstep = .true.

    cgl => leaves%first
    do while (associated(cgl))

       do ddim = zdim, xdim, -1
          if (.not. skip_sweep(ddim) .and. dom%has_dir(ddim)) call sweep_dsplit(cgl%cg,dt,ddim)
       enddo
       if (associated(problem_customize_solution)) call problem_customize_solution(.false.)
       cgl => cgl%nxt
    enddo
    call all_bnd

    if (first_run) first_run = .false.

  end subroutine fluid_update

!-------------------------------------------------------------------------------------------------------------------

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
          pb => cg%w(wna%bi)%get_sweep(ddim,i1,i2)
          u1d(iarr_all_swp(ddim,:),:) = pu(:,:)
          b_cc1d(iarr_mag_swp(ddim,:),:) = pb(:,:)
          call solve(u1d,b_cc1d, dt/cg%dl(ddim))
          pu(:,:) = u1d(iarr_all_swp(ddim,:),:)
          pb(:,:) = b_cc1d(iarr_mag_swp(ddim,:),:)
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

     select case (h_solver)
        case ("muscl")
           call slope
           call ulr_fluxes_qlr
           call riemann_wrap
           call update
        case ("rk2")
           call slope

           call ulr_to_qlr

           ! Just the slope is used to feed 1st call to Riemann solver
           call riemann_wrap

           ! Now we advance the left and right states by half timestep.
           ! The slope is already calculated and can be reused
           call du_db(du1, db1)

           call ulr_to_qlr(half*du1, half*db1)

           ! second call for Riemann problem uses states evolved to half timestep
           call riemann_wrap

           call update

        case ("rk2_muscl")
           call rk2_muscl(u, b_cc, dtodx)
        case ("euler") ! Gives quite sharp advection, especially for CFL=0.5, but it is unstable and produces a lot of noise. Do not use, except for educational purposes.
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

           u(:,2:nx) = u(:,2:nx) + half*(du1(:,2:nx) + dtodx*(flx(:,1:nx-1) - flx(:,2:nx)))
           u(:,1) = u(:,2)
           u(:,nx) = u(:,nx-1)

           b_cc(:,2:nx) = b_cc(:,2:nx) + half*(db1(:,2:nx) + dtodx*(mag_cc(:,1:nx-1) - mag_cc(:,2:nx)))
           b_cc(:,1) = b_cc(:,2)
           b_cc(:,nx) = b_cc(:,nx-1)

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

           u(:,2:nx) = u(:,2:nx) + (du1(:,2:nx) + 2*du2(:,2:nx) + 2*du3(:,2:nx) + dtodx*(flx(:,1:nx-1) - flx(:,2:nx)))/6.
           u(:,1) = u(:,2)
           u(:,nx) = u(:,nx-1)

           b_cc(:,2:nx) = b_cc(:,2:nx) + (db1(:,2:nx) + 2*db2(:,2:nx) + 2*db3(:,2:nx) + dtodx*(mag_cc(:,1:nx-1) - mag_cc(:,2:nx)))/6.
           b_cc(:,1) = b_cc(:,2)
           b_cc(:,nx) = b_cc(:,nx-1)

        case default
           call die("[fluidupdate:sweep_dsplit] No recognized solver")
     end select

     contains

        ! some shortcuts

        subroutine slope()

           use constants,  only: half

           implicit none

           real, dimension(size(u,1),size(u,2))       :: du
           real, dimension(size(b_cc,1),size(b_cc,2)) :: db

           du  = calculate_slope_vanleer(u)
           ul  = u - half*du
           ur  = u + half*du

           db  = calculate_slope_vanleer(b_cc)
           b_ccl = b_cc - half*db
           b_ccr = b_cc + half*db

        end subroutine slope

        subroutine ulr_to_qlr(du, db)

           implicit none

           real, optional, dimension(size(u,1),size(u,2)),       intent(in) :: du
           real, optional, dimension(size(b_cc,1),size(b_cc,2)), intent(in) :: db

           real, dimension(size(u,1),size(u,2))               :: u_l, u_r

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
              call riemann_hlld(nx, p_flx, p_ql, p_qr, p_bcc, p_bccl, p_bccr, fl%gam)
           enddo

        end subroutine riemann_wrap

        subroutine update

           implicit none

           u(:,2:nx) = u(:,2:nx) + dtodx*(flx(:,1:nx-1) - flx(:,2:nx))
           u(:,1) = u(:,2)
           u(:,nx) = u(:,nx-1)

           b_cc(:,2:nx) = b_cc(:,2:nx) + dtodx*(mag_cc(:,1:nx-1) - mag_cc(:,2:nx))
           b_cc(:,1) = b_cc(:,2)
           b_cc(:,nx) = b_cc(:,nx-1)

        end subroutine update

  end subroutine solve

!---------------------------------------------------------------------------------------------------------------------

  function calculate_slope_vanleer(u) result(dq)

      implicit none

      real, dimension(:,:), intent(in)     :: u

      real, dimension(size(u,1),size(u,2)) :: dlft, drgt, dcen, dq
      integer :: n


      n = size(u,2)

      dlft(:,2:n)   = (u(:,2:n) - u(:,1:n-1)) ; dlft(:,1) = dlft(:,2)    ! (14.38)
      drgt(:,1:n-1) = dlft(:,2:n) ;             drgt(:,n) = drgt(:,n-1)

      dcen = dlft*drgt

      where (dcen>0.0)
         dq = 2.0*dcen / (dlft+drgt)       ! (14.54)
      elsewhere
         dq = 0.0
      endwhere

   end function calculate_slope_vanleer

!-----------------------------------------------------------------------------------------------------------------------

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

  !---------------------------------------------------------------------------------------------------------------------

  subroutine rk2_muscl(u,b_cc, dtodx)

    use constants,   only: half, xdim, zdim
    use fluidindex,  only: flind
    use fluidtypes,  only: component_fluid
    use hlld,        only: fluxes, riemann_hlld

    implicit none

    real, dimension(:,:),           intent(inout)   :: u
    real, dimension(:,:),           intent(inout)   :: b_cc
    real,                           intent(in)      :: dtodx

    class(component_fluid), pointer                           :: fl
    real, dimension(size(b_cc,1),size(u,2)), target           :: b_ccl, b_ccr, mag_cc
    real, dimension(size(b_cc,1),size(u,2)), target           :: db
    real, dimension(size(u,1),size(u,2))                      :: u_l, u_r
    real, dimension(size(u,1)+size(b_cc,1),size(u,2)), target :: flx
    real, dimension(size(u,1),size(u,2)), target              :: ql, qr, du, ul, ur
    real, dimension(:,:), pointer                             :: p_flx, p_bcc, p_bccl, p_bccr, p_ql, p_qr
    integer                                                   :: nx, i
    real, dimension(size(u,1),size(u,2))                      :: uhalf

    nx  = size(u,2)

    du  = calculate_slope_vanleer(u)
    ul  = u - half*du
    ur  = u + half*du

    db  = calculate_slope_vanleer(b_cc)
    b_ccl = b_cc - half*db
    b_ccr = b_cc + half*db

    flx  = fluxes(ul,b_ccl) - fluxes(ur,b_ccr)

    u_l = ur + half*dtodx*flx(:size(u,1),:)
    u_r(:,1:nx-1) = ul(:,2:nx) + half*dtodx*flx(:size(u,1),2:nx)
    u_r(:,nx) = u_r(:,nx-1)

    ql = utoq(u_l,b_ccl)
    qr = utoq(u_r,b_ccr)

    ! MUSCL-Hancock is used to feed 1st call to Riemann solver
    do i = 1, flind%fluids
       fl    => flind%all_fluids(i)%fl
       p_flx => flx(fl%beg:fl%end,:)
       p_ql  => ql(fl%beg:fl%end,:)
       p_qr  => qr(fl%beg:fl%end,:)
       p_bcc => mag_cc(xdim:zdim,:)
       p_bccl => b_ccl(xdim:zdim,:)
       p_bccr => b_ccr(xdim:zdim,:)
       call riemann_hlld(nx, p_flx, p_ql, p_qr, p_bcc, p_bccl, p_bccr, fl%gam)
    enddo

    ! Now we can calculate state for half-timestep and recalculate slopes
    uhalf(:,2:nx) = u(:,2:nx) + half*dtodx*(flx(:size(u,1),1:nx-1) - flx(:size(u,1),2:nx))
    uhalf(:,1) = uhalf(:,2) ; uhalf(:,nx) = uhalf(:,nx-1)

    du  = calculate_slope_vanleer(uhalf)
    ul = uhalf - half*du
    ur = uhalf + half*du

    db  = calculate_slope_vanleer(b_cc)
    b_ccl = b_cc - half*db
    b_ccr = b_cc + half*db

    u_l = ur
    u_r(:,1:nx-1) = ul(:,2:nx)
    u_r(:,nx) = u_r(:,nx-1)

    ql = utoq(u_l,b_ccl)
    qr = utoq(u_r,b_ccr)

    ! second call for Riemann problem needs just the slope from states evolved to half timestep
    do i = 1, flind%fluids
       fl    => flind%all_fluids(i)%fl
       p_flx => flx(fl%beg:fl%end,:)
       p_ql  => ql(fl%beg:fl%end,:)
       p_qr  => qr(fl%beg:fl%end,:)
       p_bcc => mag_cc(xdim:zdim,:)
       p_bccl => b_ccl(xdim:zdim,:)
       p_bccr => b_ccr(xdim:zdim,:)
       call riemann_hlld(nx, p_flx, p_ql, p_qr, p_bcc, p_bccl, p_bccr, fl%gam)
    enddo

    u(:,2:nx) = u(:,2:nx) + dtodx*(flx(:size(u,1),1:nx-1) - flx(:size(u,1),2:nx))
    u(:,1) = u(:,2)
    u(:,nx) = u(:,nx-1)

  end subroutine rk2_muscl

!----------------------------------------------------------------------------------------------------------------------------------

end module fluidupdate
