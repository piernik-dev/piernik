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
!    Hyperbolic divergence cleaning
!    A. Dedner et al., Journal of Computational Physics 175, 645-673 (2002)
!    Dr. Varadarajan Parthasarathy, Dr. Artur Gawryszczak, CAMK (Warszawa)
!    November 2017
!
#include "piernik.h"

module hdc
! pulled by MAGNETIC && RIEMANN
  implicit none

  private
  public :: chspeed

contains

  function chspeed() ! Temporary fix

    use fluidindex, only: flind
    use fluidtypes, only: phys_prop

    type(phys_prop),          pointer            :: sn
    integer(kind=4)                              :: ifl
    real                                         :: chspeed

    do ifl = lbound(flind%all_fluids, 1, kind=4), ubound(flind%all_fluids, 1, kind=4)
       sn => flind%all_fluids(ifl)%fl%snap
    end do

    chspeed = sn%cs_max%val

    return
    
  end function chspeed
   
end module hdc

! function chspeed(u,b_cc,gamma) ! Calculates Eq 15 from Mignone, Tzeferacos, JCoPh 229 (2010) 2117-2138

!     use constants,  only: half, zero, xdim, ydim, zdim
!     use fluidindex, only: flind
!     use fluidtypes, only: component_fluid
!     use func,       only: ekin
    
!     implicit none
    
!     real, dimension(:,:), intent(in)             :: u
!     real, dimension(:,:), intent(in)             :: b_cc
!     real                , intent(in)             :: gamma
    
    
    

!     real, dimension(size(u,2))                           :: vx, vy, vz, pr
!     integer                                              :: ip, boff
!     real, dimension(size(u,2))                           :: v_amp, cf
!     class(component_fluid), pointer                      :: fl
!     real, parameter                                      :: four = 4.0
!     real, dimension(size(u,2))                           :: chspeed

!     boff = size(u, 1) ! assume xdim == 1

!     chspeed = zero

!     do ip = 1, flind%fluids

!        fl => flind%all_fluids(ip)%fl

!        vx  =  u(fl%imx,:)/u(fl%idn,:)
!        vy  =  u(fl%imy,:)/u(fl%idn,:)
!        vz  =  u(fl%imz,:)/u(fl%idn,:)

!        if (fl%has_energy) then
!           ! Gas pressure without magnetic fields. Pg 317, Eq. 2. (1) and (2) are markers for HD and MHD
!           pr = fl%gam_1*(u(fl%ien,:) - ekin(u(fl%imx,:), u(fl%imy,:), u(fl%imz,:), u(fl%idn,:))) ! (1)
!           if (fl%is_magnetized) then
!              ! Gas pressure with magnetic fields. Pg 317, Eq. 2.
!              pr = pr - half*fl%gam_1*sum(b_cc(xdim:zdim,:)**2, dim=1) ! (2)
!           endif
!        else
!           ! Dust
!           pr = 0.
!        endif

!        v_amp = sqrt( sum(vx(:)*vx(:) + vy(:)*vy(:) + vz(:)*vz(:) ) )
!        cf    = sqrt(half*((gamma*pr(:) +sum(b_cc(xdim:zdim,:)**2,dim=1) ) + sqrt((gamma*pr(:) + sum(b_cc(xdim:zdim,:)**2,dim=1))**2 - four*gamma*pr(:)*b_cc(xdim,:)**2))/u(fl%idn,:))
!        chspeed = max(chspeed,v_amp+cf)

!     end do

!     return
    
!   end function chspeed
  
