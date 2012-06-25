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
!! \brief Module that collects all flux components from each fluid
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
   public  :: all_fluxes, flimiter, set_limiter, init_fluxes, cleanup_fluxes

   logical, save                              :: fluxes_initialized = .false.

   interface
      subroutine limiter(f,a,b)
         implicit none
         real, dimension(:,:), intent(in)      :: a
         real, dimension(:,:), intent(in)      :: b
         real, dimension(:,:), intent(inout)   :: f
      end subroutine limiter
   end interface

   interface
      subroutine flux_interface(flux,cfr,uu,n,vx,ps,bb,cs_iso2)
         implicit none
         integer(kind=4), intent(in)                  :: n         !< number of cells in the current sweep
         real, dimension(:,:), intent(inout), pointer :: flux      !< flux of fluid
         real, dimension(:,:), intent(in),    pointer :: uu        !< part of u for fluid
         real, dimension(:,:), intent(inout), pointer :: cfr       !< freezing speed for fluid
         real, dimension(:,:), intent(in),    pointer :: bb        !< magnetic field x,y,z-components table
         real, dimension(:),   intent(inout), pointer :: vx        !< velocity of fluid for current sweep
         real, dimension(:),   intent(inout), pointer :: ps        !< pressure of fluid for current sweep
         real, dimension(:),   intent(in),    pointer :: cs_iso2   !< isothermal sound speed squared
      end subroutine flux_interface
   end interface

   type :: flux_func
      procedure(flux_interface), pointer, nopass :: flux_func
   end type flux_func

   type(flux_func), dimension(:), allocatable :: flist
   procedure(limiter), pointer :: flimiter

contains

   subroutine init_fluxes
      use constants,   only: NEU, DST, ION
      use fluxdust,    only: flux_dst
      use fluidindex,  only: flind
      use fluxionized, only: flux_ion
      use fluxneutral, only: flux_neu
      implicit none
      integer :: i

      allocate(flist(flind%fluids))

      do i = 1, flind%fluids
         select case (flind%all_fluids(i)%fl%tag)
            case (NEU)
               flist(flind%all_fluids(i)%fl%pos)%flux_func => flux_neu
            case (DST)
               flist(flind%all_fluids(i)%fl%pos)%flux_func => flux_dst
            case (ION)
               flist(flind%all_fluids(i)%fl%pos)%flux_func => flux_ion
         end select
      enddo

      fluxes_initialized = .true.
   end subroutine init_fluxes

   subroutine cleanup_fluxes

      implicit none

      if (allocated(flist))  deallocate(flist)

   end subroutine cleanup_fluxes

!>
!! \brief Subroutine which changes flux and cfr from mhdflux regarding specified fluids.
!! \param flux flux
!! \param cfr freezing speed
!! \param uu currently used fluid table
!! \param bb magnetic field x,y,z-components table
!! \param n number of cells in the current sweep
!! \param cs_iso2 isothermal sound speed squared
!! \param pp
!<

   subroutine all_fluxes(n, flux, cfr, uu, bb, pp, vx, cs_iso2)
      use fluidtypes,     only: component_fluid
#ifdef COSM_RAYS
      use fluxcosmicrays, only: flux_crs
#endif /* COSM_RAYS */
#ifdef TRACER
      use fluxtracer,     only: flux_tracer
      use inittracer,     only: trace_fluid
#endif /* TRACER */
      use fluidindex,     only: flind, nmag

      implicit none

      integer(kind=4), intent(in)                           :: n        !< size of input arrays
      real, dimension(flind%all,n),    target,  intent(out) :: flux     !< array storing all fluxes
      real, dimension(flind%all,n),    target,  intent(out) :: cfr      !< array storing all freezing speeds
      real, dimension(flind%all,n),    target,  intent(out) :: uu       !< array with current fluid state
      real, dimension(nmag,n),         target,  intent(in)  :: bb       !< array with current magnetic field state
      real, dimension(flind%fluids,n), target,  intent(out) :: vx       !< array storing velocity in current sweep direction (reused later)
      real, dimension(flind%fluids,n), target,  intent(out) :: pp       !< array storing pressure in current sweeo (reused later)
      real, dimension(:),              pointer, intent(in)  :: cs_iso2  !< array with current sound speed squared

      real, dimension(:,:),            pointer              :: pflux, pcfr, puu, pbb
      real, dimension(:),              pointer              :: pvx, ppp
#ifdef TRACER
      real, dimension(:),              pointer              :: pu1d, pfl1d
#endif /* TRACER */
      class(component_fluid),          pointer              :: pfl
      integer                                               :: p
!>
!! \todo pbb may need more careful treatment
!!  currently it doesn't matter for dst and neu and
!!  is properly set for ion
!<
      pbb   =>   bb(:,:)

      do p = 1, flind%fluids
         pfl   => flind%all_fluids(p)%fl
         puu   =>   uu(pfl%beg:pfl%end,:)
         pcfr  =>  cfr(pfl%beg:pfl%end,:)
         pflux => flux(pfl%beg:pfl%end,:)
         pvx   =>   vx(pfl%pos,:)
         ppp   =>   pp(pfl%pos,:)

         call flist(pfl%pos)%flux_func(pflux, pcfr, puu, n, pvx, ppp, pbb, cs_iso2)
      enddo

#ifdef COSM_RAYS
      puu   => uu(flind%crs%beg:flind%crs%end,:)
      pflux => flux(flind%crs%beg:flind%crs%end,:)
      pvx   => vx(flind%ion%pos,:)

      call flux_crs(pflux,pvx,puu,n)

      cfr(flind%crs%beg:flind%crs%end,:)  = spread(cfr(flind%ion%iarr(1),:),1,flind%crs%all)
#endif /* COSM_RAYS */

#ifdef TRACER
      do p = 1, size(trace_fluid)
         pu1d  => uu(flind%trc%beg + p - 1, :)
         pfl1d => flux(flind%trc%beg + p - 1, :)
         pvx   => vx(flind%all_fluids(trace_fluid(p))%fl%pos, :)
         call flux_tracer(pfl1d, pu1d, pvx)
         cfr(flind%trc%beg + p - 1,:)  = cfr(flind%all_fluids(trace_fluid(p))%fl%iarr(1),:)
      enddo
#endif /* TRACER */

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

      ! OPT: We don't really need to know c(:,:) everywhere, we even don't need a(:,:)*b(:,:) to be a separate array
      ! \todo benchmark an alternative code, and leave a note about results
      ! where ((a(:,:)>0 .and. b(:,:)>0) .or. (a(:,:)<0 .and b(:,:)<0))
      !   f = f+2.0*a*b/(a+b)
      ! endwhere
      c = a*b                                                                    !> \todo OPTIMIZE ME
      where (c > 0.0)
         f = f+2.0*c/(a+b)
      endwhere
      return
   end subroutine vanleer_limiter

   subroutine moncen_limiter(f,a,b)
      use constants, only: one, two, half
      implicit none
      real, dimension(:,:), intent(in)      :: a
      real, dimension(:,:), intent(in)      :: b
      real, dimension(:,:), intent(inout)   :: f

      f = f+(sign(one,a)+sign(one,b))*min(two*abs(a),two*abs(b),half*abs(a+b))*half  !> \todo OPTIMIZE ME
      return
   end subroutine moncen_limiter

   subroutine minmod_limiter(f,a,b)
      use constants, only: one, half
      implicit none
      real, dimension(:,:), intent(in)      :: a
      real, dimension(:,:), intent(in)      :: b
      real, dimension(:,:), intent(inout)   :: f

      f = f+(sign(one,a)+sign(one,b))*min(abs(a),abs(b))*half                    !> \todo OPTIMIZE ME
      return
   end subroutine minmod_limiter

   subroutine superbee_limiter(f,a,b)
      use constants, only: one, two, half
      implicit none
      real, dimension(:,:), intent(in)      :: a
      real, dimension(:,:), intent(in)      :: b
      real, dimension(:,:), intent(inout)   :: f

      where (abs(a) > abs(b))                                                    !> \todo OPTIMIZE ME
         f = f+(sign(one,a)+sign(one,b))*min(abs(a), abs(two*b))*half
      elsewhere
         f = f+(sign(one,a)+sign(one,b))*min(abs(two*a), abs(b))*half
      endwhere
      return
   end subroutine superbee_limiter
end module fluxes
