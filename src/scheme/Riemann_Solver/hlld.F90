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
  public :: riemann_hlld, fluxes

contains

  subroutine fluxes(n, q, u, f, gamma, b, cdim)
 
    use constants,   only: half,one,idn, imx, imy, imz, ien, xdim, ydim, zdim
    use fluidindex, only: iarr_mag_swp
    use func,       only: emag, ekin
    
    implicit none
     
    integer,                       intent(in)    :: n  
    real,                          intent(in)    :: gamma  
    real, dimension(:,:),          intent(in)    :: q, u
    real, dimension(:,:),          intent(out)   :: f
    real, dimension(:,:),          intent(in)    :: b         
    integer(kind=4)                              :: ibx, iby, ibz
    integer(kind=4),               intent(in)    :: cdim

    ! Local variables

    integer :: i
    real    :: vb, pr, prtot

    ibx = iarr_mag_swp(cdim,xdim)
    iby = iarr_mag_swp(cdim,ydim)
    ibz = iarr_mag_swp(cdim,zdim)

    do i = 1, n

       vb   = sum((q(imx:imz,i)/q(idn,i))*q(ibx:ibz,i))
       pr   = (gamma-one)*(q(ien,i) - ekin(q(imx,i), q(imy,i), q(imz,i), q(idn,i)) - emag(q(ibx,i),q(iby,i),q(ibz,i)))
       prtot = pr + emag(q(ibx,i),q(iby,i),q(ibz,i))

   ! MHD fluxes
       
       f(idn,i) = u(imx,i)
       f(imx,i) = u(imx,i)*(q(imx,i)/q(idn,i)) - q(ibx,i)*q(ibx,i) - prtot
       f(imy,i) = u(imy,i)*(q(imx,i)/q(idn,i)) - q(ibx,i)*q(iby,i)
       f(imz,i) = u(imz,i)*(q(imx,i)/q(idn,i)) - q(ibx,i)*q(ibz,i)
       f(ien,i) = (u(ien,i)+prtot)*(q(imx,i)/q(idn,i)) - q(ibx,i)*vb

       f(iby,i) = q(iby,i)*(q(imx,i)/q(idn,i)) - q(ibx,i)*(q(imy,i)/q(idn,i))
       f(ibz,i) = q(ibz,i)*(q(imx,i)/q(idn,i)) - q(ibz,i)*(q(imz,i)/q(idn,i))

    end do

  end subroutine fluxes

!-------------------------------------------------------------------------------------------------------------------------------------------------

  subroutine riemann_hlld(n, gamma, qleft, qright, b, cdim, f)

    use constants,  only: one, two, zero, xdim, ydim, zdim, idn, imx, imy, imz, ien
    !use domain,     only: dom
    use fluidindex, only: iarr_mag_swp
    use func,       only: emag, ekin
    !use fluidtypes  only: component_fluid
    use grid_cont,  only: grid_container
    use fluxes,     only: all_fluxes, flimiter
    

    implicit none
    
    real, parameter  ::  four  = 4.0
    integer,                       intent(in)    :: n      
    real,                          intent(in)    :: gamma  
    real, dimension(:,:), pointer, intent(in)    :: qleft, qright
    real, dimension(:,:),          intent(out)   :: f
    real, dimension(:,:),          intent(in)    :: b
    integer(kind=4)                              :: ibx, iby, ibz
    integer(kind=4),               intent(in)    :: cdim
    
   
    integer                                      :: i          
    !real, dimension(:,:)                         :: fleft
    !real, dimension(:,:)                         :: fright
    real, dimension(n)                           :: SL, SR
    real, dimension(n)                           :: denl,prl,ul, magprl
    real, dimension(n)                           :: denr,prr,ur, magprr
    real, dimension(n)                           :: cfastl
    real, dimension(n)                           :: cfastr

    ibx = iarr_mag_swp(cdim,xdim)
    iby = iarr_mag_swp(cdim,ydim)
    ibz = iarr_mag_swp(cdim,zdim)

    ! Copy normal component of magnetic field to left and right states

    qleft(ibx,:)  = b(ibx,:)  ! Check Artur
    qright(ibx,:) = b(ibx,:)  ! Check Artur

    ! Left variables

    denl  = qleft(idn,:)
    ul    = qleft(imx,:)/qleft(idn,:)
    magprl = emag(qleft(ibx,:), qleft(iby,:), qleft(ibz,:))
    prl   = (gamma-one)*(qleft(ien,:) - ekin(qleft(imx,:), qleft(imy,:), qleft(imz,:), denl) - magprl)
    cfastl = (gamma*prl+two*magprl)+sqrt((gamma*prl+two*magprl)*(gamma*prl+two*magprl) - (four*gamma*prl*qleft(ibx,:)*qleft(ibx,:)))
    cfastl = sqrt(cfastl/two*denl)
       
    ! Right variables
    
    denr  = qright(idn,:)
    ur    = qright(imx,:)/qright(idn,:)
    magprr =  emag(qright(ibx,:), qright(iby,:), qright(ibz,:))
    prr   = (gamma-one)*(qright(ien,:) - ekin(qright(imx,:), qright(imy,:), qright(imz,:), denr) - magprr)
    cfastr = (gamma*prr+two*magprr)+sqrt((gamma*prr+two*magprr)*(gamma*prr+two*magprr) - (four*gamma*prr*qright(ibx,:)*qright(ibx,:)))
    cfastr = sqrt(cfastr/two*denr)

    ! Compute wave speed
    
    SL = min(ul,ur) - max(cfastl,cfastr)
    SR = max(ul,ur) + max(cfastl,cfastr)

    ! HLLD flux
    
    where(SL >= zero) 

      ! f(:,:) = fleft(:,:) 

    elsewhere (SR <= zero) 

       !f(:,:) = fright(:,:)

    elsewhere ! SL < 0 < SR

      ! print*,'finish it'

    endwhere
    

  end subroutine riemann_hlld

end module hlld
