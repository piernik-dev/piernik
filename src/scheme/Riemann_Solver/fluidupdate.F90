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
  public :: fluid_update, sweep_dsplit, rk2, utoq, calculate_slope_vanleer, euler, muscl, muscl_check, euler_check, rk2_check

contains

  subroutine fluid_update

    use cg_list,      only: cg_list_element
    use cg_leaves,    only: leaves
    use constants,    only: xdim, zdim
    use domain,       only: dom
    use global,       only: dt, dtm, t
    use user_hooks,   only: problem_customize_solution
    use dataio_pub,   only: halfstep

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

       do ddim = xdim, zdim, 1
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
          if (dom%has_dir(ddim)) call sweep_dsplit(cgl%cg,dt,ddim)
       enddo
       if (associated(problem_customize_solution)) call problem_customize_solution(.false.)
       cgl => cgl%nxt
    enddo


    if (first_run) first_run = .false.

  end subroutine fluid_update

!-------------------------------------------------------------------------------------------------------------------

  subroutine sweep_dsplit(cg, dt, ddim)

    use constants,        only: pdims, xdim, zdim, ORTHO1, ORTHO2, LO, HI
    use dataio_pub,       only: die
    use all_boundaries,   only: all_fluid_boundaries
    use fluidindex,       only: iarr_all_swp
    use global,           only: h_solver
    use grid_cont,        only: grid_container
    use named_array_list, only: wna
    use dataio_pub,       only: warn

    implicit none

    interface
       subroutine solver_1d(f_data, b_data, dt_dx)
          implicit none
          real, dimension(:,:), intent(inout) :: f_data
          real, dimension(:,:), intent(inout) :: b_data
          real,                 intent(in)    :: dt_dx
       end subroutine solver_1d
    end interface

    type(grid_container), pointer, intent(in) :: cg
    real,                          intent(in) :: dt
    integer(kind=4),               intent(in) :: ddim

    real, dimension(size(cg%u,1), cg%n_(ddim)) :: u1d
    real, dimension(xdim:zdim, cg%n_(ddim))   :: b_cc1d
    real, dimension(:,:), pointer             :: pu
    integer                                   :: i1, i2
    logical, save                             :: firstcall = .true.

    procedure(solver_1d), pointer, save :: solve => NULL()

    b_cc1d = 0.
    if (firstcall) then
       call warn("[fluidupdate:sweep] magnetic field unimplemented yet. Forcing to be 0")
       select case (h_solver)
          case ("muscl")
             solve => muscl
          case ("euler_check")
             solve => euler_check
          case ("rk2_check")
             solve => rk2_check
          case ("rk2_flx")
             solve => rk2_flx
          case ("rk2_uh")
             solve => rk2_uh
          case ("euler")
             solve => euler
          case ("rk2")
             solve => rk2
           case ("rk4")
             solve => rk4
         case default
             call die("[fluidupdate:sweep_dsplit] No recognized solver")
       end select
    endif
    firstcall = .false.

    do i2 = cg%lhn(pdims(ddim, ORTHO2), LO), cg%lhn(pdims(ddim,ORTHO2), HI)
       do i1 = cg%lhn(pdims(ddim, ORTHO1), LO), cg%lhn(pdims(ddim, ORTHO1), HI)
          pu => cg%w(wna%fi)%get_sweep(ddim,i1,i2)
          u1d(iarr_all_swp(ddim,:),:) = pu(:,:)
          call solve(u1d,b_cc1d, dt/cg%dl(ddim))
          pu(:,:) = u1d(iarr_all_swp(ddim,:),:)
       enddo
    enddo

    call all_fluid_boundaries

  end subroutine sweep_dsplit

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
        if (fl%has_energy) then
           q(fl%ien,:) =  fl%gam_1*(u(fl%ien,:) - ekin(u(fl%imx,:), u(fl%imy,:), u(fl%imz,:), u(fl%idn,:)) - half*sum(b_cc(xdim:zdim,:)**2)) + half*sum(b_cc(xdim:zdim,:)**2)
        endif


     enddo

   end function utoq

  !-------------------------------------------------------------------------------------------------

  subroutine rk2(u,b_cc, dtodx)

    use constants,   only: half, xdim, zdim
    use fluidindex,  only: flind
    use fluidtypes,  only: component_fluid
    use hlld,        only: fluxes, riemann_hlld

    implicit none

    real, dimension(:,:),           intent(inout)   :: u
    real, dimension(:,:),           intent(inout)   :: b_cc
    real,                           intent(in)      :: dtodx

    class(component_fluid), pointer                 :: fl
    real, dimension(xdim:zdim,size(u,2)), target    :: b_ccl, b_ccr, mag_cc
    real, dimension(xdim:zdim,size(u,2)), target    :: db
    real, dimension(size(u,1),size(u,2)), target    :: flx, ql, qr, du, ul, ur, u_l, u_r
    real, dimension(:,:), pointer                   :: p_flx, p_bcc, p_bccl, p_bccr, p_ql, p_qr
    integer                                         :: nx
    real, dimension(size(u,1),size(u,2))            :: du1

    nx  = size(u,2)

    call u_slope_q(u)
    call riemann_wrap
    du1(:,2:nx) = dtodx*(flx(:,1:nx-1) - flx(:,2:nx))
    du1(:,1) = du1(:,2) ; du1(:,nx) = du1(:,nx-1)

    call u_slope_q(u+du1)
    call riemann_wrap
    u(:,2:nx) = u(:,2:nx) + half*(du1(:,2:nx) + dtodx*(flx(:,1:nx-1) - flx(:,2:nx)))
    u(:,1) = u(:,2) ; u(:,nx) = u(:,nx-1)

    contains

       subroutine u_slope_q(uu)

          implicit none

          real, dimension(:,:), intent(in) :: uu

          du  = calculate_slope_vanleer(uu)
          ul  = uu - half*du
          ur  = uu + half*du

          db  = calculate_slope_vanleer(b_cc)
          b_ccl = b_cc - half*db
          b_ccr = b_cc + half*db

          mag_cc = b_cc

          u_l = ur
          u_r(:,1:nx-1) = ul(:,2:nx)
          u_r(:,nx) = u_r(:,nx-1)

          ql = utoq(u_l,b_ccl)
          qr = utoq(u_r,b_ccr)

       end subroutine u_slope_q

       subroutine riemann_wrap()

          implicit none

          integer :: i

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

       end subroutine riemann_wrap

  end subroutine rk2

  !-------------------------------------------------------------------------------------------------

  subroutine rk4(u,b_cc, dtodx)

    use constants,   only: half, xdim, zdim
    use fluidindex,  only: flind
    use fluidtypes,  only: component_fluid
    use hlld,        only: fluxes, riemann_hlld

    implicit none

    real, dimension(:,:),           intent(inout)   :: u
    real, dimension(:,:),           intent(inout)   :: b_cc
    real,                           intent(in)      :: dtodx

    class(component_fluid), pointer                 :: fl
    real, dimension(xdim:zdim,size(u,2)), target    :: b_ccl, b_ccr, mag_cc
    real, dimension(xdim:zdim,size(u,2)), target    :: db
    real, dimension(size(u,1),size(u,2)), target    :: flx, ql, qr, du, ul, ur, u_l, u_r
    real, dimension(:,:), pointer                   :: p_flx, p_bcc, p_bccl, p_bccr, p_ql, p_qr
    integer                                         :: nx
    real, dimension(size(u,1),size(u,2))            :: du1, du2, du3

    nx  = size(u,2)

    call u_slope_q(u)
    call riemann_wrap
    du1(:,2:nx) = dtodx*(flx(:,1:nx-1) - flx(:,2:nx))
    du1(:,1) = du1(:,2) ; du1(:,nx) = du1(:,nx-1)

    call u_slope_q(u+half*du1)
    call riemann_wrap
    du2(:,2:nx) = dtodx*(flx(:,1:nx-1) - flx(:,2:nx))
    du2(:,1) = du2(:,2) ; du2(:,nx) = du2(:,nx-1)

    call u_slope_q(u+half*du2)
    call riemann_wrap
    du3(:,2:nx) = dtodx*(flx(:,1:nx-1) - flx(:,2:nx))
    du3(:,1) = du3(:,2) ; du3(:,nx) = du3(:,nx-1)

    call u_slope_q(u+du3)
    call riemann_wrap
    u(:,2:nx) = u(:,2:nx) + (du1(:,2:nx) + 2*du2(:,2:nx) + 2*du3(:,2:nx) + dtodx*(flx(:,1:nx-1) - flx(:,2:nx)))/6.
    u(:,1) = u(:,2) ; u(:,nx) = u(:,nx-1)

    contains

       subroutine u_slope_q(uu)

          implicit none

          real, dimension(:,:), intent(in) :: uu

          du  = calculate_slope_vanleer(uu)
          ul  = uu - half*du
          ur  = uu + half*du

          db  = calculate_slope_vanleer(b_cc)
          b_ccl = b_cc - half*db
          b_ccr = b_cc + half*db

          mag_cc = b_cc

          u_l = ur
          u_r(:,1:nx-1) = ul(:,2:nx)
          u_r(:,nx) = u_r(:,nx-1)

          ql = utoq(u_l,b_ccl)
          qr = utoq(u_r,b_ccr)

       end subroutine u_slope_q

       subroutine riemann_wrap()

          implicit none

          integer :: i

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

       end subroutine riemann_wrap

  end subroutine rk4

  !-------------------------------------------------------------------------------------------------
  ! Gives very sharp advection, especially for CFL=0.5, seems to produce a lot of noise

  subroutine euler(u,b_cc, dtodx)

    use constants,   only: half, xdim, zdim
    use fluidindex,  only: flind
    use fluidtypes,  only: component_fluid
    use hlld,        only: fluxes, riemann_hlld

    implicit none

    real, dimension(:,:),           intent(inout)   :: u
    real, dimension(:,:),           intent(inout)   :: b_cc
    real,                           intent(in)      :: dtodx

    class(component_fluid), pointer                 :: fl
    real, dimension(xdim:zdim,size(u,2)), target    :: b_ccl, b_ccr, mag_cc
    real, dimension(xdim:zdim,size(u,2)), target    :: db
    real, dimension(size(u,1),size(u,2)), target    :: flx, ql, qr, du, ul, ur, u_l, u_r
    real, dimension(:,:), pointer                   :: p_flx, p_bcc, p_bccl, p_bccr, p_ql, p_qr
    integer                                         :: nx, i !, ii

    nx  = size(u,2)

    du  = calculate_slope_vanleer(u)
    ul  = u - half*du
    ur  = u + half*du

    db  = calculate_slope_vanleer(b_cc)
    b_ccl = b_cc - half*db
    b_ccr = b_cc + half*db

    mag_cc = b_cc

    u_l = ur
    u_r(:,1:nx-1) = ul(:,2:nx)
    u_r(:,nx) = u_r(:,nx-1)

    ql = utoq(u_l,b_ccl)
    qr = utoq(u_r,b_ccr)

    ! Just the slope is used to feed Riemann solver
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

    u(:,2:nx) = u(:,2:nx) + dtodx*(flx(:,1:nx-1) - flx(:,2:nx))
    u(:,1) = u(:,2) ; u(:,nx) = u(:,nx-1)

  end subroutine euler

  !---------------------------------------------------------------------------------------------------------------------------------------------------------

  subroutine muscl(u,b_cc, dtodx)

    use constants,   only: half, xdim, zdim
    use fluidindex,  only: flind
    use fluidtypes,  only: component_fluid
    use hlld,        only: fluxes, riemann_hlld

    implicit none

    real, dimension(:,:),           intent(inout)   :: u
    real, dimension(:,:),           intent(inout)   :: b_cc
    real,                           intent(in)      :: dtodx

    class(component_fluid), pointer                 :: fl
    real, dimension(xdim:zdim,size(u,2)), target    :: b_ccl, b_ccr, mag_cc
    real, dimension(xdim:zdim,size(u,2)), target    :: db
    real, dimension(size(u,1),size(u,2))            :: u_l, u_r
    real, dimension(size(u,1),size(u,2)), target    :: flx, ql, qr, du, ul, ur !, u_l, u_r
    real, dimension(:,:), pointer                   :: p_flx, p_bcc, p_bccl, p_bccr, p_ql, p_qr
    integer                                         :: nx, i

    nx  = size(u,2)

    du  = calculate_slope_vanleer(u)
    ul  = u - half*du
    ur  = u + half*du

    db  = calculate_slope_vanleer(b_cc)
    b_ccl = b_cc - half*db
    b_ccr = b_cc + half*db

    mag_cc = b_cc

    flx  = fluxes(ul,b_ccl) - fluxes(ur,b_ccr)

    u_l = ur + half*dtodx*flx
    u_r(:,1:nx-1) = ul(:,2:nx) + half*dtodx*flx(:,2:nx) ; u_r(:,nx) = u_r(:,nx-1)

    ql = utoq(u_l,b_ccl)
    qr = utoq(u_r,b_ccr)

    do i = 1, flind%fluids
       fl    => flind%all_fluids(i)%fl
       p_flx => flx(fl%beg:fl%end,:)  ! flx is intent (out) and ignored at end of call
       p_ql  => ql(fl%beg:fl%end,:)
       p_qr  => qr(fl%beg:fl%end,:)
       p_bcc => mag_cc(xdim:zdim,:)
       p_bccl => b_ccl(xdim:zdim,:)
       p_bccr => b_ccr(xdim:zdim,:)
       call riemann_hlld(nx, p_flx, p_ql, p_qr, p_bcc, p_bccl, p_bccr, fl%gam)
    enddo

    u(:,2:nx) = u(:,2:nx) + dtodx*(flx(:,1:nx-1) - flx(:,2:nx))
    u(:,1) = u(:,2) ; u(:,nx) = u(:,nx-1)

  end subroutine muscl

  ! --------------------------------------------------------------------------------------------------------------------------------------------------------

   subroutine euler_check(u,b_cc, dtodx)

    use constants,   only: half, xdim, zdim
    use fluidindex,  only: flind
    use fluidtypes,  only: component_fluid
    use hlld,        only: fluxes, riemann_hlld

    implicit none

    real, dimension(:,:),           intent(inout)   :: u
    real, dimension(:,:),           intent(inout)   :: b_cc
    real,                           intent(in)      :: dtodx

    class(component_fluid), pointer                 :: fl
    real, dimension(xdim:zdim,size(u,2)), target    :: b_ccl, b_ccr, mag_cc
    real, dimension(xdim:zdim,size(u,2)), target    :: db
    !real, dimension(size(u,1),size(u,2))            :: u_l, u_r
    real, dimension(size(u,1),size(u,2)), target    :: flx, ql, qr, du, ul, ur, u_l, u_r
    real, dimension(:,:), pointer                   :: p_flx, p_bcc, p_bccl, p_bccr, p_ql, p_qr
    integer                                         :: nx, i!, ii
    real, dimension(size(u,1),size(u,2)) :: dflx

    nx  = size(u,2)

    du  = calculate_slope_vanleer(u)
    ul  = u - half*du
    ur  = u + half*du

    db  = calculate_slope_vanleer(b_cc)
    b_ccl = b_cc - half*db
    b_ccr = b_cc + half*db

    mag_cc = b_cc

    flx  = fluxes(ul,b_ccl) - fluxes(ur,b_ccr)

    u_l = ul !+ half*dtodx*flx
    u_r(:,1:nx-1) = ul(:,2:nx) !+ half*dtodx*flx(:,2:nx) ;
    u_r(:,nx) = u_r(:,nx-1)

    ql = utoq(u_l,b_ccl)
    qr = utoq(u_r,b_ccr)

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

    !write(*,*) "p_flx", p_flx(1,:)
    !write(*,*) "u2", u(:,2:nx)
    !write(*,*) "dtodx", dtodx
    !write(*,*) "flx1", flx(:,1:nx-1)
    !write(*,*) "flx2", flx(:,2:nx) write(*,*) "flxdiff", flx(:,1:nx-1) - flx(:,2:nx)
    !u(:,2:nx) = u(:,2:nx) + dtodx*(flx(:,1:nx-1) - flx(:,2:nx))
    !write(*,*) "flx(1,2:nx)    ", flx(1,2:nx)
    !write(*,*) "cshift(flx,1,2)", cshift(flx(1,:),1,1)
    !write(*,*) "flx", flx
    !write(*,*) "p_flx", p_flx
    dflx =  cshift(flx,-1,2) - flx
    !u = u + dtodx*(flx - cshift(flx,1,2))
    u = u + dtodx*dflx
    !dflx = flx - cshift(flx,1,2)
    !dflx =  cshift(flx,-1,2) - flx
    write(*,*) "dflx", dflx(1,:)
    !write(*,*) "u_before", u(:,:)

    !u(:,:) = u(:,:) + dtodx*flx(:,:)
    !u(:,2:nx) = u(:,2:nx) + dtodx*(flx(:,1:nx-1) - flx(:,2:nx))
    !write(*,*) "flx12diff", flx(:,1:nx-1) - flx(:,2:nx)
    !write(*,*) "u2nx_64", size(u(:,2:nx))
    !u(:,1) = u(:,2) ; u(:,nx) = u(:,nx-1)
    !write(*,*) "u_after", u(:,:)
    !write(*,*) "flx", flx(:,:)
    !write(*,*) "u", u(:,:)
    !write(*,*) "flx", flx(:,:)

    !do ii = lbound(flx, 1), ubound(flx, 1)
    !   write(*,*) "flx1", flx - cshift(flx,1,2)
    !enddo

    !do ii = lbound(flx, 2), ubound(flx, 2)
    !   write(*,*) "flx2", flx - cshift(flx,1,2)
    !enddo



    !u = u + dtodx*flx



  end subroutine euler_check

  !--------------------------------------------------------------------------------------------------------------------------------

  subroutine muscl_check(u,b_cc, dtodx)

    use constants,   only: half, xdim, zdim
    use fluidindex,  only: flind
    use fluidtypes,  only: component_fluid
    use hlld,        only: fluxes, riemann_hlld

    implicit none

    real, dimension(:,:),           intent(inout)   :: u
    real, dimension(:,:),           intent(inout)   :: b_cc
    real,                           intent(in)      :: dtodx

    class(component_fluid), pointer                 :: fl
    real, dimension(xdim:zdim,size(u,2)), target    :: b_ccl, b_ccr, mag_cc
    real, dimension(xdim:zdim,size(u,2)), target    :: db
    real, dimension(size(u,1),size(u,2))            :: u_l, u_r
    real, dimension(size(u,1),size(u,2)), target    :: flx, ql, qr, du, ul, ur !, u_l, u_r
    real, dimension(:,:), pointer                   :: p_flx, p_bcc, p_bccl, p_bccr, p_ql, p_qr
    integer                                         :: nx, i!, ii
     real, dimension(size(u,1),size(u,2)) :: dflx

    nx  = size(u,2)

    du  = calculate_slope_vanleer(u)
    ul  = u - half*du
    ur  = u + half*du

    db  = calculate_slope_vanleer(b_cc)
    b_ccl = b_cc - half*db
    b_ccr = b_cc + half*db

    mag_cc = b_cc

    flx  = fluxes(ul,b_ccl) - fluxes(ur,b_ccr)

    u_l = ur + half*dtodx*flx
    u_r(:,1:nx-1) = ul(:,2:nx) + half*dtodx*flx(:,2:nx) ; u_r(:,nx) = u_r(:,nx-1)

    ql = utoq(u_l,b_ccl)
    qr = utoq(u_r,b_ccr)

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

    !write(*,*) "p_flx", p_flx(1,:)

    u(:,2:nx) = u(:,2:nx) + dtodx*(flx(:,1:nx-1) - flx(:,2:nx))
    dflx(:,2:nx) = flx(:,1:nx-1) - flx(:,2:nx)
     write(*,*) "dflx", dflx(1,:)

    !write(*,*) "flx_diff", flx(:,1:nx-1) - flx(:,2:nx)

    !do ii = lbound(flx,1), ubound(flx,1)
    !   write(*,*) "flx1", flx(:,1:nx-1) - flx(:,2:nx)
    !enddo

    !do ii = lbound(flx,2), ubound(flx,2)
    !   write(*,*) "flx2", flx(:,1:nx-1) - flx(:,2:nx)
    !enddo

    u(:,1) = u(:,2) ; u(:,nx) = u(:,nx-1)

  end subroutine muscl_check

  !-------------------------------------------------------------------------------------------------

  subroutine rk2_check(u,b_cc, dtodx)

    use constants,   only: half, xdim, zdim
    use fluidindex,  only: flind
    use fluidtypes,  only: component_fluid
    use hlld,        only: fluxes, riemann_hlld

    implicit none

    real, dimension(:,:),           intent(inout)   :: u
    real, dimension(:,:),           intent(inout)   :: b_cc
    real,                           intent(in)      :: dtodx

    class(component_fluid), pointer                 :: fl
    real, dimension(xdim:zdim,size(u,2)), target    :: b_ccl, b_ccr, mag_cc
    real, dimension(xdim:zdim,size(u,2)), target    :: db
    real, dimension(size(u,1),size(u,2)), target    :: flx, ql, qr, du, ul, ur, u_l, u_r
    real, dimension(:,:), pointer                   :: p_flx, p_bcc, p_bccl, p_bccr, p_ql, p_qr
    integer                                         :: nx, i !, ii
    real, dimension(size(u,1),size(u,2))            :: uhalf
    real                                            ::a

    nx  = size(u,2)
    a = half

    du  = calculate_slope_vanleer(u)
    ul  = u - half*du
    ur  = u + half*du

    db  = calculate_slope_vanleer(b_cc)
    b_ccl = b_cc - half*db
    b_ccr = b_cc + half*db

    mag_cc = b_cc

    flx  = fluxes(ul,b_ccl) - fluxes(ur,b_ccr)

    u_l = ur + half*dtodx*flx
    u_r(:,1:nx-1) = ul(:,2:nx) + half*dtodx*flx(:,2:nx) ; u_r(:,nx) = u_r(:,nx-1)

    ql = utoq(u_l,b_ccl)
    qr = utoq(u_r,b_ccr)

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

    uhalf(:,2:nx) = u(:,2:nx) + a*dtodx*(flx(:,1:nx-1) - flx(:,2:nx))
    uhalf(:,1) = uhalf(:,2) ; uhalf(:,nx) = uhalf(:,nx-1)


    du  = calculate_slope_vanleer(uhalf)
    ul = uhalf - half*du
    ur = uhalf + half*du

    db  = calculate_slope_vanleer(b_cc)
    b_ccl = b_cc - half*db
    b_ccr = b_cc + half*db

    mag_cc = b_cc

    flx  = fluxes(ul,b_ccl) - fluxes(ur,b_ccr)

    u_l = ur + half*dtodx*flx
    u_r(:,1:nx-1) = ul(:,2:nx) + half*dtodx*flx(:,2:nx) ; u_r(:,nx) = u_r(:,nx-1)

    ql = utoq(u_l,b_ccl)
    qr = utoq(u_r,b_ccr)

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

    u(:,2:nx) = u(:,2:nx) + dtodx*(flx(:,1:nx-1) - flx(:,2:nx))
    u(:,1) = u(:,2) ; u(:,nx) = u(:,nx-1)

  end subroutine rk2_check

  !-------------------------------------------------------------------------------------------------

  subroutine rk2_flx(u,b_cc, dtodx)

    use constants,   only: half, xdim, zdim
    use fluidindex,  only: flind
    use fluidtypes,  only: component_fluid
    use hlld,        only: fluxes, riemann_hlld

    implicit none

    real, dimension(:,:),           intent(inout)   :: u
    real, dimension(:,:),           intent(inout)   :: b_cc
    real,                           intent(in)      :: dtodx

    class(component_fluid), pointer                 :: fl
    real, dimension(xdim:zdim,size(u,2)), target    :: b_ccl, b_ccr, mag_cc
    real, dimension(xdim:zdim,size(u,2)), target    :: db
    real, dimension(size(u,1),size(u,2)), target    :: flx, ql, qr, du, ul, ur, u_l, u_r
    real, dimension(:,:), pointer                   :: p_flx, p_bcc, p_bccl, p_bccr, p_ql, p_qr
    integer                                         :: nx, i !, ii

    nx  = size(u,2)

    du  = calculate_slope_vanleer(u)
    ul  = u - half*du
    ur  = u + half*du

    db  = calculate_slope_vanleer(b_cc)
    b_ccl = b_cc - half*db
    b_ccr = b_cc + half*db

    mag_cc = b_cc

    u_l = ur
    u_r(:,1:nx-1) = ul(:,2:nx); u_r(:,nx) = u_r(:,nx-1)

    ql = utoq(u_l,b_ccl)
    qr = utoq(u_r,b_ccr)

    ! Just the slope is used to feed 1st call to Riemann solver
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

    ! Now we advance the left and right states by half timestep.
    ! The slope is already calculated and can be reused
    u_l(:,2:nx) = ur(:,2:nx) + half*dtodx*(flx(:,1:nx-1) - flx(:,2:nx))
    u_l(:,1) = u_l(:,2)

    u_r(:,1:nx-1) = ul(:,2:nx) + half*dtodx*(flx(:,1:nx-1) - flx(:,2:nx))
    u_r(:,nx) = u_r(:,nx-1)

    ql = utoq(u_l,b_ccl)
    qr = utoq(u_r,b_ccr)

    ! second call for Riemann problem uses states evolved to half timestep
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

    u(:,2:nx) = u(:,2:nx) + dtodx*(flx(:,1:nx-1) - flx(:,2:nx))
    u(:,1) = u(:,2) ; u(:,nx) = u(:,nx-1)

  end subroutine rk2_flx

  !-------------------------------------------------------------------------------------------------

  subroutine rk2_uh(u,b_cc, dtodx)

    use constants,   only: half, xdim, zdim
    use fluidindex,  only: flind
    use fluidtypes,  only: component_fluid
    use hlld,        only: fluxes, riemann_hlld

    implicit none

    real, dimension(:,:),           intent(inout)   :: u
    real, dimension(:,:),           intent(inout)   :: b_cc
    real,                           intent(in)      :: dtodx

    class(component_fluid), pointer                 :: fl
    real, dimension(xdim:zdim,size(u,2)), target    :: b_ccl, b_ccr, mag_cc
    real, dimension(xdim:zdim,size(u,2)), target    :: db
    real, dimension(size(u,1),size(u,2)), target    :: flx, ql, qr, du, ul, ur, u_l, u_r
    real, dimension(:,:), pointer                   :: p_flx, p_bcc, p_bccl, p_bccr, p_ql, p_qr
    integer                                         :: nx, i !, ii
    real, dimension(size(u,1),size(u,2))            :: uhalf

    nx  = size(u,2)

    du  = calculate_slope_vanleer(u)
    ul  = u - half*du
    ur  = u + half*du

    db  = calculate_slope_vanleer(b_cc)
    b_ccl = b_cc - half*db
    b_ccr = b_cc + half*db

    mag_cc = b_cc

    flx  = fluxes(ul,b_ccl) - fluxes(ur,b_ccr)

    u_l = ur + half*dtodx*flx
    u_r(:,1:nx-1) = ul(:,2:nx) + half*dtodx*flx(:,2:nx) ; u_r(:,nx) = u_r(:,nx-1)

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
    uhalf(:,2:nx) = u(:,2:nx) + half*dtodx*(flx(:,1:nx-1) - flx(:,2:nx))
    uhalf(:,1) = uhalf(:,2) ; uhalf(:,nx) = uhalf(:,nx-1)

    du  = calculate_slope_vanleer(uhalf)
    ul = uhalf - half*du
    ur = uhalf + half*du

    db  = calculate_slope_vanleer(b_cc)
    b_ccl = b_cc - half*db
    b_ccr = b_cc + half*db

    mag_cc = b_cc

    flx  = fluxes(ul,b_ccl) - fluxes(ur,b_ccr)

    u_l = ur
    u_r(:,1:nx-1) = ul(:,2:nx) ; u_r(:,nx) = u_r(:,nx-1)

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

    u(:,2:nx) = u(:,2:nx) + dtodx*(flx(:,1:nx-1) - flx(:,2:nx))
    u(:,1) = u(:,2) ; u(:,nx) = u(:,nx-1)

  end subroutine rk2_uh

end module fluidupdate
