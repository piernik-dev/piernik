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
!  References:
!
!  A multi-state HLLD approximate Riemann solver for ideal magnetohydrodynamics.
!  Takahiro Miyoshi, Kanya Kusano
!  Journal of Computational Physics 208 (2005) 315-344
!
!  ->Solve one dimensional Riemann problem using adiabatic HLLD scheme
!
!  Varadarajan Parthasarathy, CAMK, Warszawa. 2015.
!
!---------------------------------------------------------------------------------------------------------------------------


#include "piernik.def"
!>
!!  \brief This module implements HLLD Riemann solver following the work of Miyoshi & Kusano (2005)
!<
module hlld
! pulled by RIEMANN

  implicit none

  private
  public :: riemann_hlld

contains

  function fluxes(u,b,bb,cs2) result(f)

    use constants,  only: half, xdim, ydim, zdim
    use fluidindex, only: flind
    use fluidtypes, only: component_fluid
    use func,       only: ekin

    implicit none

    real, dimension(:,:),      intent(in)    :: u
    real, dimension(:,:),      intent(in)    :: b
    real, dimension(:,:),      intent(out)   :: bb
    real, dimension(:), pointer, intent(in)  :: cs2

    real, dimension(size(u,1), size(u,2))    :: f
    real, dimension(size(u,2))               :: vx, vy, vz, p
    integer :: ip
    class(component_fluid),    pointer       :: fl

    do ip = 1, flind%fluids

       fl => flind%all_fluids(ip)%fl

       vx  =  u(fl%imx,:)/u(fl%idn,:)
       vy  =  u(fl%imy,:)/u(fl%idn,:)
       vz  =  u(fl%imz,:)/u(fl%idn,:)

       if(fl%has_energy) then
          p = fl%gam_1*(u(fl%ien,:) - ekin(u(fl%imx,:), u(fl%imy,:), u(fl%imz,:), u(fl%idn,:)) - half*sum(b**2,dim=2))
       else
          if(associated(cs2)) then
             p = cs2*u(fl%idn,:)
          else
             p = 0.
          endif
       endif
       
       f(fl%idn,:)  =  u(fl%imx,:)
       f(fl%imx,:)  =  u(fl%imx,:)*vx + p - b(xdim,:)*b(xdim,:)
       f(fl%imy,:)  =  u(fl%imy,:)*vx - b(xdim,:)*b(ydim,:)
       f(fl%imz,:)  =  u(fl%imz,:)*vx - b(xdim,:)*b(zdim,:)
       bb(:,ydim)    =  b(ydim,:)*vx - b(xdim,:)*vy
       bb(:,zdim)    =  b(zdim,:)*vx - b(xdim,:)*vz
       if(fl%has_energy) then
          f(fl%ien,:) = (u(fl%ien,:) + p(:))*vx(:) - b(xdim,:)*(b(xdim,:)*vx(:) + b(ydim,:)*vy(:) + b(zdim,:)*vz(:))
       endif
    enddo

    return

  end function fluxes
  

!-------------------------------------------------------------------------------------------------------------------------------------------------

  subroutine riemann_hlld(n, gamma, uleft, uright, b, cdim, f)

    ! external procedures
    
    use constants,  only: half, xdim, ydim, zdim, idn, imx, imy, imz, ien
    use fluidindex, only: iarr_mag_swp, flind
    !use func,       only: emag, ekin
    !use grid_cont,  only: grid_container
    !use fluxes,     only: all_fluxes, flimiter
    
    ! arguments
    
    implicit none

    integer,                       intent(in)    :: n      
    real,                          intent(in)    :: gamma  
    real, dimension(:,:), pointer, intent(in)    :: uleft, uright
    real, dimension(:,:), pointer, intent(out)   :: f
    real, dimension(:,:),          intent(in)    :: b
    integer(kind=4)                              :: ibx, iby, ibz
    integer(kind=4),               intent(in)    :: cdim

    !
    ! Local variables
    !
    
    ! left state

    real, dimension(n)                           :: c_fastl  !< fast magneto-acoustic wave-left
    real, dimension(n)                           :: SL       !< speed of left moving wave
    real, dimension(n)                           :: SL_star  !< speed of left moving Alfven wave
    real, dimension(n)                           :: gampr_l  !< gamma*pressure
    real, dimension(n)                           :: ul       !< conserved variable-left
    real, dimension(n,flind%all)                 :: fl       !< flux left
    real, dimension(n)                           :: SLSM     !< SL - SM
    real, dimension(n)                           :: SLVL     !< SL - vx_left
    real, dimension(n)                           :: SMVL     !< SM - vx_left
    real, dimension(n)                           :: DNLSLVL  !< rho_left*SLVL
    
    ! right state

    real, dimension(n)                           :: c_fastr  !< fast magneto-acoustic wave-right
    real, dimension(n)                           :: SR       !< speed of right moving wave
    real, dimension(n)                           :: SR_star  !< speed of right moving Alfven wave
    real, dimension(n)                           :: gampr_r  !< gamma*pressure
    real, dimension(n)                           :: ur       !< conserved variable-right
    real, dimension(n,flind%all)                 :: fr       !< flux right
    real, dimension(n)                           :: SRSM     !< SR - SM
    real, dimension(n)                           :: SRVR     !< SR - vx_right
    real, dimension(n)                           :: SMVR     !< SM - vx_right
    real, dimension(n)                           :: DNRSRVR   !< rho_right*SRVR

    
    integer                                      :: i
    real, parameter                              :: four  = 4.0
    !real, dimension(n)                           :: c_fast        !< fast magneto-acoustic wave
    real, dimension(n)                           :: SM            !< entropy wave Eq. (38)
    real, dimension(n)                           :: SM_nr, SM_dr  !< numerator and denominator of Eq. (38)
    real, dimension(n)                           :: SRSL          !< SR - SL
    
    
    ibx = iarr_mag_swp(cdim,xdim)
    iby = iarr_mag_swp(cdim,ydim)
    ibz = iarr_mag_swp(cdim,zdim)

    ! solver

    do i = 1,n

       ul  =  uleft(imx,i)
       ur  =  uright(imx,i)
    
       gampr_l  =  gamma*uleft(ien,i)
       gampr_r  =  gamma*uright(ien,i)

       c_fastl  =   (gampr_l+(uleft(ibx,i)**2+uleft(iby,i)**2+uleft(ibz,i)**2))  &
            + sqrt((gampr_l+(uleft(ibx,i)**2+uleft(iby,i)**2+uleft(ibz,i)**2))**2-(four*gampr_l*uleft(ibx,i)**2))

       c_fastl  =  sqrt(half*c_fastl/uleft(idn,i))

       c_fastr  =   (gampr_r+(uright(ibx,i)**2+uright(iby,i)**2+uright(ibz,i)**2))  &
            + sqrt((gampr_r+(uright(ibx,i)**2+uright(iby,i)**2+uright(ibz,i)**2))**2-(four*gampr_r*uright(ibx,i)**2))

       c_fastr  =  sqrt(half*c_fastr/uright(idn,i))

       ! Eq. (67)
    
       SL  =  min(ul,ur) - max(c_fastl,c_fastr)
       SR  =  max(ul,ur) + max(c_fastl,c_fastr)

       ! Speed of contact discontinuity Eq. (38)

       SM_nr  =  (SR*ur - SL*ul) - (fr(imx,i) - fl(imx,i))
       SM_dr  =  (SR*uright(idn,i) - SL*uleft(idn,i)) - (fr(idn,i) - fl(idn,i))
       SM     =  SM_nr/SM_dr

       ! Speed differences

       SLSM  =  SL - SM
       SRSM  =  SR - SM

       SLVL  =  SL - (uleft(imx,i)/uleft(idn,i))
       SRVR  =  SR - (uright(imx,i)/uright(idn,i))

       SMVL  =  SM - (uleft(imx,i)/uleft(idn,i))
       SMVR  =  SM - (uright(imx,i)/uright(idn,i))

       SRSL  =  SR - SL
    
       ! Co-efficients

       DNLSLVL  =  uleft(idn,i)*SLVL
       DNRSRVR  =  uright(idn,i)*SRVR
       
       ! 

       
    end do
    

    
    

  end subroutine riemann_hlld

end module hlld
