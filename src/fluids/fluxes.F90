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
!    Initial implementation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.h"
!>
!! \brief (MH/JD) (doxy comments ready) Module that collects all flux components from each fluid
!!
!!The %fluxes for all fluids are combined in the same order as conservative variables in the array \a u(:,:,:,:)
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
! pulled by ANY
   implicit none
   private
   public  :: all_fluxes, flimiter, set_limiter

   logical, save                              :: fluxes_initialized = .false.

   interface
      subroutine limiter(f,a,b)
         implicit none
         real, dimension(:,:), intent(in)      :: a
         real, dimension(:,:), intent(in)      :: b
         real, dimension(:,:), intent(inout)   :: f
      end subroutine limiter
   end interface

   procedure(limiter), pointer :: flimiter

contains

   subroutine init_fluxes
      implicit none

      fluxes_initialized = .true.
   end subroutine init_fluxes
!>
!! \brief Subroutine which changes flux and cfr from mhdflux regarding specified fluids.
!! \param flux flux
!! \param cfr freezing speed
!! \param uu currently used fluid table
!! \param bb magnetic field x,y,z-components table
!! \param n number of cells in the current sweep
!<
   subroutine all_fluxes(n, flux, cfr, uu, bb, cs_iso2)
      use types,          only: component_fluid
#ifdef IONIZED
      use fluxionized,    only: flux_ion
#endif /* IONIZED */
#ifdef NEUTRAL
      use fluxneutral,    only: flux_neu
#endif /* NEUTRAL */
#ifdef DUST
      use fluxdust,       only: flux_dst
#endif /* DUST */
#ifdef COSM_RAYS
      use fluxcosmicrays, only: flux_crs
#endif /* COSM_RAYS */
      use fluidindex, only: nvar, nmag

      implicit none

      integer,                      intent(in)  :: n
      real, dimension(nvar%all,n),  intent(out), target :: flux, cfr, uu
      real, dimension(nmag,n),      intent(in),  target :: bb
      real, dimension(n), optional, intent(in),  target :: cs_iso2

      real, dimension(nvar%fluids,n), target            :: vx

      real, dimension(:,:), pointer                     :: pflux, pcfr, puu, pbb
      real, dimension(:), pointer                       :: pcs2, pvx
      type(component_fluid), pointer                    :: pfl

#ifndef IONIZED
      integer :: dummy
#endif /* !IONIZED */

#ifdef IONIZED
      pfl   => nvar%ion
      puu   =>   uu(pfl%beg:pfl%end,:)
      pcfr  =>  cfr(pfl%beg:pfl%end,:)
      pflux => flux(pfl%beg:pfl%end,:)
      pvx   =>   vx(pfl%pos,:)
      pbb   =>   bb(:,:)
      if (present(cs_iso2)) then
         pcs2  => cs_iso2(:)
      else
         pcs2  => null()
      endif

      call flux_ion(pflux, pcfr, puu, n, pvx, pbb, pcs2)
#else /* !IONIZED */
      if (.false.) dummy = size(bb)*size(cs_iso2) ! suppress compiler warnings on unused arguments
#endif /* !IONIZED */

#ifdef NEUTRAL
      pfl   => nvar%neu
      puu   =>   uu(pfl%beg:pfl%end,:)
      pcfr  =>  cfr(pfl%beg:pfl%end,:)
      pflux => flux(pfl%beg:pfl%end,:)
      pvx   =>  vx(pfl%pos,:)
      pbb   => null()
      pcs2  => null()

      call flux_neu(pflux, pcfr, puu, n, pvx, pbb, pcs2)
#endif /* NEUTRAL */

#ifdef DUST
      pfl   => nvar%dst
      puu   =>   uu(pfl%beg:pfl%end,:)
      pcfr  =>  cfr(pfl%beg:pfl%end,:)
      pflux => flux(pfl%beg:pfl%end,:)
      pvx   =>  vx(pfl%pos,:)
      pbb   => null()
      pcs2  => null()

      call flux_dst(pflux, pcfr, puu, n, pvx, pbb, pcs2)
#endif /* DUST */

#ifdef COSM_RAYS
      puu   => uu(nvar%crs%beg:nvar%crs%end,:)
      pflux => flux(nvar%crs%beg:nvar%crs%end,:)
      pvx   => vx(nvar%ion%pos,:)

      call flux_crs(pflux,pvx,puu,n)

      cfr(nvar%crs%beg:nvar%crs%end,:)  = spread(cfr(nvar%ion%iarr(1),:),1,nvar%crs%all)
#endif /* COSM_RAYS */

   end subroutine all_fluxes

!==========================================================================================
!/*
!>
!! \brief This subroutine applies flux limiter.
!!
!! Flux limiter is a function used when interpolation of %fluxes onto cell boundaries is made to avoid the spurious oscillations.
!!
!! You can choose between van Leer's, monotonized central, minmod or superbee flux limiters. The chosen flux limiter has to be defined in
!! file piernik.def. Default is "vanleer"
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
!*/
   subroutine set_limiter(lname)
      use dataio_pub, only: msg, die
#ifdef VERBOSE
      use dataio_pub, only: printinfo
#endif /* VERBOSE */
      implicit none
      character(len=*), intent(in) :: lname
      if (associated(flimiter)) call die("[fluxes:set_limiter] flimiter already associated")
      select case (lname)
         case ('vanleer', 'VANLEER')
            flimiter => vanleer_limiter
         case ('minmod', 'MINMOD')
            flimiter => minmod_limiter
         case ('moncen', 'MONCEN')
            flimiter => moncen_limiter
         case ('superbee', 'SUPERBEE')
            flimiter => superbee_limiter
         case default
            write(msg,'(2a)') "[fluxes:set_limiter] unknown limiter ", lname
            call die(msg)
      end select
#ifdef VERBOSE
      write(msg,'(2a)') "[fluxes:set_limiter] limiter set to ", lname
      call printinfo(msg)
#endif /* VERBOSE */
   end subroutine set_limiter

   subroutine vanleer_limiter(f,a,b)
      implicit none
      real, dimension(:,:), intent(in)      :: a !< second order correction of left- or right- moving waves flux on the left cell boundary
      real, dimension(:,:), intent(in)      :: b !< second order correction of left- or right- moving waves flux on the right cell boundary
      real, dimension(:,:), intent(inout)   :: f !< second order flux correction for left- or right- moving waves
      ! locals
      real, dimension(size(a,1), size(a,2)) :: c !< a*b

      c = a*b                                                                    ! ToDO: OPTIMIZE ME
      where (c > 0.0)
         f = f+2.0*c/(a+b)
      endwhere
      return
   end subroutine vanleer_limiter

   subroutine moncen_limiter(f,a,b)
      implicit none
      real, dimension(:,:), intent(in)      :: a
      real, dimension(:,:), intent(in)      :: b
      real, dimension(:,:), intent(inout)   :: f

      f = f+(sign(1.0,a)+sign(1.0,b))*min(2.*abs(a),2.*abs(b),0.5*abs(a+b))*0.5  ! ToDo: OPTIMIZE ME
      return
   end subroutine moncen_limiter

   subroutine minmod_limiter(f,a,b)
      implicit none
      real, dimension(:,:), intent(in)      :: a
      real, dimension(:,:), intent(in)      :: b
      real, dimension(:,:), intent(inout)   :: f

      f = f+(sign(1.0,a)+sign(1.0,b))*min(abs(a),abs(b))*0.5                     ! ToDo: OPTIMIZE ME
      return
   end subroutine minmod_limiter

   subroutine superbee_limiter(f,a,b)
      implicit none
      real, dimension(:,:), intent(in)      :: a
      real, dimension(:,:), intent(in)      :: b
      real, dimension(:,:), intent(inout)   :: f

      where (abs(a) > abs(b))                                                    ! ToDo: OPTIMIZE ME
         f = f+(sign(1.0,a)+sign(1.0,b))*min(abs(a), abs(2.0*b))*0.5
      elsewhere
         f = f+(sign(1.0,a)+sign(1.0,b))*min(abs(2.0*a), abs(b))*0.5
      endwhere
      return
   end subroutine superbee_limiter
end module fluxes
