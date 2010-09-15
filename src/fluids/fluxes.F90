! $Id$
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
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.def"
!>
!! \brief (MH/JD) (doxy comments ready) Module that collects all flux components from each fluid
!!
!!The %fluxes for all fluids are combined in the same order as conservarive variables in the array \a u(:,:,:,:)
!!\f{equation}
!!\vec{F}(\vec{u},\vec{B})
!!  = \big(\vec{F}^i(\vec{u}^i,\vec{B}), \vec{F}^n(\vec{u}^n),
!!  \vec{F}^d(\vec{u}^d)\big)^T,
!!\f}
!!where the elementary flux vectors like \f$\vec{F}^i(\vec{u}^i,\vec{B})$,
!!$\vec{F}^n(\vec{u}^n)\f$,  \f$\vec{F}^d(\vec{u}^d)\f$ are %fluxes computed
!!independently for each fluid. In multidimensional computations the %fluxes
!!\f$\vec{G}(\vec{u},\vec{B})\f$  and \f$\vec{H}(\vec{u},\vec{B})\f$,
!!corresponding to the transport of conservative quantities in \f$y\f$ and
!!\f$z\f$--directions, are constructed in a similar way.
!!
!! \warning This module should not be changed by user.
!<
module fluxes

#ifdef IONIZED
  use initionized,    only : iarr_ion
  use fluxionized,    only : flux_ion
#endif /* IONIZED */
#ifdef NEUTRAL
  use initneutral,    only : iarr_neu
  use fluxneutral,    only : flux_neu
#endif /* NEUTRAL */
#ifdef DUST
  use initdust,       only : iarr_dst
  use fluxdust,       only : flux_dst
#endif /* DUST */
#ifdef COSM_RAYS
  use initcosmicrays, only : iarr_crs
  use fluxcosmicrays, only : flux_crs
#endif /* COSM_RAYS */

contains

!>
!! \brief Subroutine which changes flux and cfr from mhdflux regarding specified fluids.
!! \param flux flux
!! \param cfr freezing speed
!! \param uu currently used fluid table
!! \param bb magnetic field x,y,z-components table
!! \param n number of cells in the current sweep
!<
subroutine all_fluxes(n, flux, cfr, uu, bb, cs_iso2)

    use fluidindex, only : nvar, nmag

    implicit none

    integer,                      intent(in)  :: n
    real, dimension(nvar%all,n),  intent(out) :: flux, cfr, uu
    real, dimension(nmag,n),      intent(in)  :: bb
    real, dimension(n), optional, intent(in)  :: cs_iso2

    real, dimension(n)              :: vion

#ifdef IONIZED
    real, dimension(nvar%ion%all,n) :: fluxion,cfrion,uuion
#else
    integer :: dummy
#endif /* IONIZED */

#ifdef NEUTRAL
    real, dimension(nvar%neu%all,n) :: fluxneu,cfrneu,uuneu
#endif /* NEUTRAL */

#ifdef DUST
    real, dimension(nvar%dst%all,n) :: fluxdst,cfrdst,uudst
#endif /* DUST */

#ifdef COSM_RAYS
    real, dimension(nvar%crs%all,n) :: fluxcrs,uucrs
#endif /* COSM_RAYS */

   vion(:) = 0.0

#ifdef IONIZED
   uuion(:,:)=uu(iarr_ion,:)

   call flux_ion(fluxion,cfrion,vion,uuion,bb,n,cs_iso2)

   flux(iarr_ion,:) = fluxion
   cfr(iarr_ion,:)  = cfrion
   uu(iarr_ion,:)   = uuion
#else
   if (.false.) dummy = size(bb)*size(cs_iso2) ! suppress compiler warnings on unused arguments
#endif /* IONIZED */

#ifdef NEUTRAL
   uuneu(:,:)=uu(iarr_neu,:)

   call flux_neu(fluxneu,cfrneu,uuneu,n)

   flux(iarr_neu,:) = fluxneu
   cfr(iarr_neu,:)  = cfrneu
   uu(iarr_neu,:)   = uuneu
#endif /* NEUTRAL */

#ifdef DUST
   uudst=uu(iarr_dst,:)
   call flux_dst(fluxdst,cfrdst,uudst,n)
   flux(iarr_dst,:) = fluxdst
   cfr(iarr_dst,:)  = cfrdst
   uu(iarr_dst,:)   = uudst
#endif /* DUST */

#ifdef COSM_RAYS
   uucrs=uu(iarr_crs,:)
   call flux_crs(fluxcrs,vion,uucrs,n)
   flux(iarr_crs,:) = fluxcrs
   cfr(iarr_crs,:)  = spread(cfrion(1,:),1,nvar%crs%all)
   uu(iarr_crs,:)   = uucrs
#endif /* COSM_RAYS */


end subroutine all_fluxes

!==========================================================================================

!>
!! \brief This subroutine applies flux limiter.
!!
!! Flux limiter is a function used when interpolation of %fluxes onto cell boundaries is made to avoid the spurious oscillations.
!!
!! You can choose between van Leer's, monotonized central, minmod or superbee flux limiters. The chosen flux limiter has to be defined in
!! file piernik.def.
!!
!! BEWARE: There is no default choice. If one forgets to define a flux limiter, none will be applied and error will be reported.
!!
!! The van Leer flux limiter can be noted as:
!! \f{equation}
!! \Delta \vec{F}_{i+1/2}^{(2)L} =
!! \left\{\begin{array}{lll}
!! \frac{2\Delta \vec{F}_{i+1/2}^{L-} \Delta \vec{F}_{i+1/2}^{L+}}{\Delta \vec{F}_{i+1/2}^{L-} + \Delta \vec{F}_{i+1/2}^{L+}}
!! &\textrm{ if } &\Delta \vec{F}_{i+1/2}^{L-} \Delta \vec{F}_{i+1/2}^{L+} > 0, \\
!! 0 & \textrm{ if } &\Delta \vec{F}_{i+1/2}^{L-} \Delta \vec{F}_{i+1/2}^{L+} < 0,
!! \end{array}\right.
!! \f}
!! \f{equation}
!! \Delta \vec{F}_{i+1/2}^{(2)R} =
!! \left\{\begin{array}{lll}
!! \frac{2\Delta \vec{F}_{i+1/2}^{R-} \Delta \vec{F}_{i+1/2}^{R+}}{\Delta \vec{F}_{i+1/2}^{R-} + \Delta \vec{F}_{i+1/2}^{R+}}
!!  &\textrm{ if } &\Delta \vec{F}_{i+1/2}^{R-} \Delta \vec{F}_{i+1/2}^{R+} > 0, \\
!! 0 & \textrm{ if } &\Delta \vec{F}_{i-1/2}^{R-} \Delta \vec{F}_{i+1/2}^{R+} < 0,
!! \end{array}\right.
!! \f}
!! where
!! \f{eqnarray}
!! \Delta \vec{F}_{i+1/2}^{L-} = \frac{1}{2} (\vec{F}_{i+1}^L-\vec{F}_{i}^L), &\qquad &
!! \Delta \vec{F}_{i+1/2}^{L+} = \frac{1}{2} (\vec{F}_{i+2}^L-\vec{F}_{i+1}^L), \\
!! \Delta \vec{F}_{i+1/2}^{R-} = \frac{1}{2} (\vec{F}_{i}^R-\vec{F}_{i-1}^R), &\qquad &
!! \Delta \vec{F}_{i+1/2}^{R+} = \frac{1}{2} (\vec{F}_{i+1}^R-\vec{F}_{i}^R),
!! \f}
!! are %fluxes of left- and right-moving waves interpolated to cell boundaries.
!<
  subroutine flimiter(f,a,b,m,n)
    implicit none
    integer, intent(in)  :: m  !< number of conservative variables
    integer, intent(in)  :: n  !< array size
    real, dimension(m,n) :: f  !< second order flux correction for left- or right- moving waves
    real, dimension(m,n) :: a  !< second order correction of left- or right- moving waves flux on the left cell boundary
    real, dimension(m,n) :: b  !< second order correction of left- or right- moving waves flux on the right cell boundary
    real, dimension(m,n) :: c  !< a*b
#ifdef VANLEER
      c = a*b
      where (c .gt. 0.0)
        f = f+2.0*c/(a+b)
      endwhere
#endif /* VANLEER */
#ifdef MONCEN
        f = f+(sign(1.0,a)+sign(1.0,b))*min(2.*abs(a),2.*abs(b),0.5*abs(a+b))*0.5
#endif /* MONCEN */
#ifdef MINMOD
        f = f+(sign(1.0,a)+sign(1.0,b))*min(abs(a),abs(b))*0.5
#endif /* MINMOD */
#ifdef SUPERBEE
      where (abs(a) .gt. abs(b))
        f = f+(sign(1.0,a)+sign(1.0,b))*min(abs(a), abs(2.0*b))*0.5
      elsewhere
        f = f+(sign(1.0,a)+sign(1.0,b))*min(abs(2.0*a), abs(b))*0.5
      endwhere
#endif /* SUPERBEE */

    return
  end subroutine flimiter

end module fluxes
