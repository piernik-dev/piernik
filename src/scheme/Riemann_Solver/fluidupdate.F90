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
  public :: fluid_update, sweep_dsplit, rk2, utoq, calculate_slope_vanleer
  

contains

  subroutine fluid_update

    use constants,    only: xdim, zdim
    use domain,       only: dom
    use global,       only: dt, dtm, t
    use user_hooks,   only: problem_customize_solution
    use dataio_pub,   only: halfstep

    implicit none

    logical, save                   :: first_run = .true.
    integer(kind=4)                 :: ddim
    

    halfstep = .false.
    if (first_run) then
       dtm = 0.0
    else
       dtm = dt
    endif
    t = t + dt

    do ddim = xdim, zdim
       if (dom%has_dir(ddim)) call sweep_dsplit(dt,ddim)
    enddo
    if (associated(problem_customize_solution)) call problem_customize_solution(.true.)

    t = t + dt
    dtm = dt
    halfstep = .true.

    do ddim = zdim, xdim, -1
       if (dom%has_dir(ddim)) call sweep_dsplit(dt,ddim)
    enddo
    if (associated(problem_customize_solution)) call problem_customize_solution(.false.)

    if (first_run) first_run = .false.

  end subroutine fluid_update

!-------------------------------------------------------------------------------------------------------------------

  subroutine sweep_dsplit(dt, ddim)

    use cg_list,          only: cg_list_element
    use constants,        only: pdims, xdim, zdim, ORTHO1, ORTHO2, LO, HI
    use all_boundaries,   only: all_fluid_boundaries
    use fluidindex,       only: iarr_all_swp
    use grid_cont,        only: grid_container
    use cg_leaves,        only: leaves
    use named_array_list, only: wna
    use dataio_pub,       only: die

    implicit none

    real,                          intent(in) :: dt
    integer(kind=4),               intent(in) :: ddim
    real, allocatable, dimension(:,:)         :: u1d
    real, allocatable, dimension(:,:)         :: b_cc1d
    real, dimension(:,:), pointer             :: pu
    integer                                   :: i1, i2
    type(cg_list_element), pointer            :: cgl

    cgl => leaves%first
    do while (associated(cgl))

       if (allocated(u1d) .or. allocated(b_cc1d)) call die("[fluidupdate:sweep_dsplit] sweep vectors not clean")
       !OPT: can check if the size is already right and avoid reallocation on AMR domains
       allocate(u1d(size(cgl%cg%u,1),cgl%cg%n_(ddim)), b_cc1d(xdim:zdim, cgl%cg%n_(ddim)))

       b_cc1d = 0.
       
       do i2 = cgl%cg%lhn(pdims(ddim,ORTHO2),LO), cgl%cg%lhn(pdims(ddim,ORTHO2),HI)
          do i1 = cgl%cg%lhn(pdims(ddim, ORTHO1), LO), cgl%cg%lhn(pdims(ddim, ORTHO1), HI)
             pu => cgl%cg%w(wna%fi)%get_sweep(ddim,i1,i2)
             !b1d => cgl%cg%w(wna%bi)%get_sweep(ddim,i1,i2)
             ! WARNING: what position of the b components do we expect? If cell-centered then interpolation is required above
 !write(*,*)"sweep",ddim, i1,i2,dt
             u1d(iarr_all_swp(ddim,:),:) = pu(:,:)

             !call rk2(u1d,b_cc1d,ddim, dt/cgl%cg%dl(ddim))
             call euler(u1d,b_cc1d,ddim, dt/cgl%cg%dl(ddim))

             pu(:,:) = u1d(iarr_all_swp(ddim,:),:)
             
          enddo
       enddo

       deallocate(u1d, b_cc1d)

       cgl => cgl%nxt
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
         dq = 2.0*dcen / (dlft+drgt)       ! (14.54) ?
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

     real, dimension(:,:),   intent(in)    :: u, b_cc
     real, dimension(size(u,1),size(u,2))  :: q
     integer  :: p

     class(component_fluid), pointer       :: fl

     do p = 1, flind%fluids
        fl => flind%all_fluids(p)%fl

        q(fl%idn,:) =  u(fl%idn,:)
        q(fl%imx,:) =  u(fl%imx,:)/u(fl%idn,:)
        q(fl%imy,:) =  u(fl%imy,:)/u(fl%idn,:)
        q(fl%imz,:) =  u(fl%imz,:)/u(fl%idn,:)

        q(fl%ien,:) =  fl%gam_1*(u(fl%ien,:) - ekin(u(fl%imx,:), u(fl%imy,:), u(fl%imz,:), u(fl%idn,:)) - half*sum(b_cc(xdim:zdim,:))**2) + half*sum(b_cc(xdim:zdim,:))**2
       
     enddo
     
   end function utoq
   

!-----------------------------------------------------------------------------------------------------------------------

  subroutine rk2(u,b_cc,ddim, dtodx)

    use constants,   only: half, xdim, zdim
    use fluidindex,  only: flind
    use fluidtypes,  only: component_fluid


    implicit none

    real, dimension(:,:),           intent(inout)   :: u
    real, dimension(:,:),           intent(inout)   :: b_cc
    integer(kind=4),                intent(in)      :: ddim
    real,                           intent(in)      :: dtodx

    class(component_fluid), pointer                 :: fl
    real, dimension(xdim:zdim,size(u,2)), target    :: b_ccl, b_ccr, mag_cc
    real, dimension(xdim:zdim,size(u,2)), target    :: db
    real, dimension(size(u,1),size(u,2))            :: u_predict
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

    flx  = fluxes(ul,b_ccl,ddim) - fluxes(ur,b_ccr,ddim)
 
    ql = utoq(ul,b_ccl)
    qr = utoq(ur,b_ccr)

    u_predict  =  u + dtodx*flx(:,:)
    
    do i = 1, flind%fluids
       fl    => flind%all_fluids(i)%fl
       p_ql  => ql(fl%beg:fl%end,:)
       p_qr  => qr(fl%beg:fl%end,:)
       p_flx => flx(fl%beg:fl%end,:)
       p_bcc => mag_cc(xdim:zdim,:)
       p_bccl => b_ccl(xdim:zdim,:)
       p_bccr => b_ccr(xdim:zdim,:)
       call riemann_hlld(nx, p_flx, p_ql, p_qr, mag_cc, p_bccl, p_bccr, fl%gam)
    end do

    !u  =  half*(u + u_predict + dtodx*flx(:,:))
    u  =  half*(u + u_predict + dtodx*(fluxes(p_ql,p_bccl,ddim)-fluxes(p_qr,p_bccr,ddim)))
    

    
  end subroutine rk2

  !---------------------------------------------------------------------------------------------------------------------------------------------------------------------

  subroutine euler(u,b_cc,ddim, dtodx)

    use constants,   only: half, xdim, zdim
    !use fluidindex,  only: flind
    use fluidtypes,  only: component_fluid


    implicit none

    real, dimension(:,:),           intent(inout)   :: u
    real, dimension(:,:),           intent(inout)   :: b_cc
    integer(kind=4),                intent(in)      :: ddim
    real,                           intent(in)      :: dtodx

    !class(component_fluid), pointer                 :: fl
    real, dimension(xdim:zdim,size(u,2)), target    :: b_ccl, b_ccr, mag_cc
    real, dimension(xdim:zdim,size(u,2)), target    :: db
    real, dimension(size(u,1),size(u,2)), target    :: ql, qr, du, ul, ur, flx
    integer                                         :: nx


    nx  = size(u,2)

    du  = calculate_slope_vanleer(u)
    ul  = u - half*du
    ur  = u + half*du
    !write(*,*) "ul", ul
    !write(*,*) "ur", ur
    
    db  = calculate_slope_vanleer(b_cc)
    b_ccl = b_cc - half*db
    b_ccr = b_cc + half*db

    mag_cc = b_cc
 
    ql = utoq(ul,b_ccl)
    qr = utoq(ur,b_ccr)

    flx = fluxes(ul,b_ccl,ddim) - fluxes(ur,b_ccr,ddim)
    !write(*,*) "flx", flx(:,:)
    
    u  =  u + dtodx*flx(:,:)
    write(*,*) "u", u

  end subroutine euler

    !---------------------------------------------------------------------------------------------------------------------------------------------------------
    
end module fluidupdate
