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
!! \brief Routines used to interpolate grid fields onto particle positions and vice versa
!<

module particle_maps
! pulled by ANY
   implicit none

   private
   public :: set_map, map_particles, quadratic_spline

   procedure(map_scheme), pointer :: map_particles => NULL()

   abstract interface
      subroutine map_scheme(iv, factor)
         implicit none
         integer(kind=4), intent(in) :: iv     !< index in cg%q array, where we want the particles to be projected
         real,            intent(in) :: factor !< typically fpiG
      end subroutine map_scheme
   end interface

contains

!>
!! \brief Set function that projects particles onto grid
!<

   subroutine set_map(ischeme)

      use constants,  only: I_NGP, I_CIC, I_TSC
      use dataio_pub, only: die

      implicit none

      integer(kind=4), intent(in) :: ischeme !< index of interpolation scheme

      select case (ischeme)
         case (I_NGP)
            map_particles => map_ngp
         case (I_CIC)
            map_particles => map_cic
         case (I_TSC)
            map_particles => map_tsc
         case default
            call die("[particle_interpolation:set_map] Interpolation scheme selector's logic in particle_pub:init_particles is broken. Go fix it!")
      end select

   end subroutine set_map

!>
!! \brief Project the particles onto density map
!!
!! \details With the help of multigrid self-gravity solver the gravitational potential of the particle set can be found
!!
!! \warning Current implementation of the multipole solver isn't aware of particles, co if any of them exist outside of the domain, their potential will be ignored.
!! \todo Fix the multipole solver
!!
!! \todo Add an option for less compact mapping that nullifies self-forces
!!
!! \warning Particles outside periodic domain are ignored
!<

   subroutine map_ngp(iv, factor)

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: xdim, ydim, zdim, ndims, LO, HI, GEO_XYZ
      use dataio_pub, only: die
      use domain,     only: dom

      implicit none

      integer(kind=4), intent(in)       :: iv     !< index in cg%q array, where we want the particles to be projected
      real,            intent(in)       :: factor !< typically fpiG

      type(cg_list_element), pointer    :: cgl
      integer                           :: p
      integer(kind=8), dimension(ndims) :: ijkp

      cgl => leaves%first
      do while (associated(cgl))

         ijkp(:) = cgl%cg%ijkse(:,LO)
         do p = lbound(cgl%cg%pset%p, dim=1), ubound(cgl%cg%pset%p, dim=1)
            associate( field => cgl%cg%q(iv)%arr, part => cgl%cg%pset%p(p), cg => cgl%cg )

            if (dom%geometry_type /= GEO_XYZ) call die("[particle_maps:map_ngp] Unsupported geometry")
            where (dom%has_dir(:)) ijkp(:) = floor((part%pos(:) - cg%fbnd(:, LO)) * cg%idl(:)) + cg%ijkse(:, LO)
            if (all(ijkp >= cg%ijkse(:,LO)) .and. all(ijkp <= cg%ijkse(:,HI))) &
                 field(ijkp(xdim), ijkp(ydim), ijkp(zdim)) = field(ijkp(xdim), ijkp(ydim), ijkp(zdim)) + factor * part%mass / cg%dvol

            end associate
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine map_ngp

   subroutine map_cic(iv, factor)

      use cg_leaves, only: leaves
      use cg_list,   only: cg_list_element
      use constants, only: xdim, ydim, zdim, ndims, LO, HI, CENTER
      use domain,    only: dom

      implicit none

      integer(kind=4), intent(in)              :: iv     !< index in cg%q array, where we want the particles to be projected
      real,            intent(in)              :: factor !< typically fpiG

      type(cg_list_element), pointer           :: cgl
      integer                                  :: p, cdim
      integer(kind=8)                          :: cn, i, j, k
      integer(kind=8), dimension(ndims, LO:HI) :: ijkp
      integer(kind=8), dimension(ndims)        :: cur_ind
      real                                     :: weight

      cgl => leaves%first
      do while (associated(cgl))

         do p = lbound(cgl%cg%pset%p, dim=1), ubound(cgl%cg%pset%p, dim=1)
            associate( field => cgl%cg%q(iv)%arr, part => cgl%cg%pset%p(p), cg => cgl%cg )

               if (part%outside) cycle

               do cdim = xdim, zdim
                  if (dom%has_dir(cdim)) then
                     cn = nint((part%pos(cdim) - cg%coord(CENTER, cdim)%r(1))*cg%idl(cdim)) + 1
                     if (cg%coord(CENTER, cdim)%r(cn) > part%pos(cdim)) then
                        ijkp(cdim, LO) = cn - 1
                     else
                        ijkp(cdim, LO) = cn
                     endif
                     ijkp(cdim, HI) = ijkp(cdim, LO) + 1
                  else
                     ijkp(cdim, :) = cg%ijkse(cdim, :)
                  endif
               enddo
               do i = ijkp(xdim, LO), ijkp(xdim, HI)
                  do j = ijkp(ydim, LO), ijkp(ydim, HI)
                     do k = ijkp(zdim, LO), ijkp(zdim, HI)
                        weight = factor * part%mass / cg%dvol
                        cur_ind(:) = [i, j, k]
                        do cdim = xdim, zdim
                           if (dom%has_dir(cdim)) &
                              weight = weight*( 1.0 - abs(part%pos(cdim) - cg%coord(CENTER, cdim)%r(cur_ind(cdim)))*cg%idl(cdim) )
                        enddo
                        field(i,j,k) = field(i,j,k) +  weight
                      enddo
                  enddo
               enddo

            end associate
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine map_cic

   subroutine map_tsc(iv, factor)

      use cg_leaves, only: leaves
      use cg_list,   only: cg_list_element
      use constants, only: xdim, ydim, zdim, ndims, LO, HI, IM, I0, IP, CENTER, half
      use domain,    only: dom

      implicit none

      integer(kind=4), intent(in)    :: iv     !< index in cg%q array, where we want the particles to be projected
      real,            intent(in)    :: factor !< typically fpiG

      type(cg_list_element), pointer           :: cgl
      integer                                  :: p, cdim
      integer(kind=8)                          :: i, j, k
      integer(kind=8), dimension(ndims, IM:IP) :: ijkp
      integer(kind=8), dimension(ndims)        :: cur_ind
      real                                     :: weight, delta_x, weight_tmp

      cgl => leaves%first
      do while (associated(cgl))

         do p = lbound(cgl%cg%pset%p, dim=1), ubound(cgl%cg%pset%p, dim=1)
            associate( field => cgl%cg%q(iv)%arr, part => cgl%cg%pset%p(p), cg => cgl%cg )

               if (part%outside) cycle

               do cdim = xdim, zdim
                  if (dom%has_dir(cdim)) then
                     ijkp(cdim, I0) = nint((part%pos(cdim) - cg%coord(CENTER, cdim)%r(1))*cg%idl(cdim)) + 1   !!! BEWARE hardcoded magic
                     ijkp(cdim, IM) = max(ijkp(cdim, I0) - 1, int(cg%lhn(cdim, LO), kind=8))
                     ijkp(cdim, IP) = min(ijkp(cdim, I0) + 1, int(cg%lhn(cdim, HI), kind=8))
                  else
                     ijkp(cdim, IM) = cg%ijkse(cdim, LO)
                     ijkp(cdim, I0) = cg%ijkse(cdim, LO)
                     ijkp(cdim, IP) = cg%ijkse(cdim, HI)
                  endif
               enddo

               do i = ijkp(xdim, IM), ijkp(xdim, IP)
                  do j = ijkp(ydim, IM), ijkp(ydim, IP)
                     do k = ijkp(zdim, IM), ijkp(zdim, IP)

                        cur_ind(:) = [i, j, k]
                        weight = 1.0

                        do cdim = xdim, zdim
                           if (.not.dom%has_dir(cdim)) cycle
                           delta_x = ( part%pos(cdim) - cg%coord(CENTER, cdim)%r(cur_ind(cdim)) ) * cg%idl(cdim)
                           if (cur_ind(cdim) /= ijkp(cdim, 0)) then   !!! BEWARE hardcoded magic
                              weight_tmp = 1.125 - 1.5 * abs(delta_x) + half * delta_x**2
                           else
                              weight_tmp = 0.75 - delta_x**2
                           endif
                           weight = weight_tmp * weight
                        enddo

                        field(i, j, k) = field(i, j, k) + factor * (part%mass / cg%dvol) * weight
                      enddo
                  enddo
               enddo
            end associate
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine map_tsc

!>
!! \brief Quadratic spline interpolation
!<
    function quadratic_spline(f, dp, ijkp) result(gp)

      use constants, only: xdim, ydim, zdim, IM, I0, IP, ndims
      use domain,    only: dom

      implicit none

      real, dimension(:,:,:,:), pointer, intent(in) :: f     !< field
      real, dimension(ndims),            intent(in) :: dp    !< point position in [0,1] range
      integer, dimension(ndims),         intent(in) :: ijkp  !< nearest grid point where dp resides
      real, dimension(size(f,1))                    :: gp    !< interpolated value

      ! locals
      real, dimension(ndims, IM:IP) :: fac

!     This should be moved to particle specific modules, to avoid depending on  cg
!     do cdim = xdim, zdim
!        dp(cdim) = (xp(cdim) - cg%coord(cdim)%r(ijkp(cdim))) * cg%idl(ijkp(cdim))
!     enddo

      !  Interpolation coefficients
      fac(:, :) = 0.0
      where (dom%has_dir)
         fac(:, IM) = 0.5 * (0.5 - dp(:))**2
         fac(:, I0) = 0.75 - dp(:)**2
         fac(:, IP) = 0.5 * (0.5 + dp(:))**2
      elsewhere
         fac(:, IM) = 0.0
         fac(:, I0) = 1.0
         fac(:, IP) = 0.0
      endwhere

      ! Some computations can be avoided by spliting following statement into pieces
      ! and using ifs
      gp(:) = fac(xdim, I0)*fac(ydim, I0)*fac(zdim, I0)*f(ijkp(xdim),ijkp(ydim),ijkp(zdim),:) + &
              fac(xdim, I0)*fac(ydim, I0)*( f(ijkp(xdim)             ,ijkp(ydim)             ,ijkp(zdim)+dom%D_(zdim),:)*fac(zdim, IP) + &
                                            f(ijkp(xdim)             ,ijkp(ydim)             ,ijkp(zdim)-dom%D_(zdim),:)*fac(zdim, IM) ) + &
              fac(xdim, I0)*fac(zdim, I0)*( f(ijkp(xdim)             ,ijkp(ydim)+dom%D_(ydim),ijkp(zdim)             ,:)*fac(ydim, IP) + &
                                            f(ijkp(xdim)             ,ijkp(ydim)-dom%D_(ydim),ijkp(zdim)             ,:)*fac(ydim, IM) ) + &
              fac(ydim, I0)*fac(zdim, I0)*( f(ijkp(xdim)+dom%D_(xdim),ijkp(ydim)             ,ijkp(zdim)             ,:)*fac(xdim, IP) + &
                                            f(ijkp(xdim)-dom%D_(xdim),ijkp(ydim)             ,ijkp(zdim)             ,:)*fac(xdim, IM) ) + &
              fac(xdim, IP)*fac(ydim, IP)*( f(ijkp(xdim)+dom%D_(xdim),ijkp(ydim)+dom%D_(ydim),ijkp(zdim)+dom%D_(zdim),:)*fac(zdim, IP) + &
                                            f(ijkp(xdim)+dom%D_(xdim),ijkp(ydim)+dom%D_(ydim),ijkp(zdim)-dom%D_(zdim),:)*fac(zdim, IM) ) + &
              fac(xdim, IP)*fac(ydim, IM)*( f(ijkp(xdim)+dom%D_(xdim),ijkp(ydim)-dom%D_(ydim),ijkp(zdim)+dom%D_(zdim),:)*fac(zdim, IP) + &
                                            f(ijkp(xdim)+dom%D_(xdim),ijkp(ydim)-dom%D_(ydim),ijkp(zdim)-dom%D_(zdim),:)*fac(zdim, IM) ) + &
              fac(xdim, IM)*fac(ydim, IP)*( f(ijkp(xdim)-dom%D_(xdim),ijkp(ydim)+dom%D_(ydim),ijkp(zdim)+dom%D_(zdim),:)*fac(zdim, IP) + &
                                            f(ijkp(xdim)-dom%D_(xdim),ijkp(ydim)+dom%D_(ydim),ijkp(zdim)-dom%D_(zdim),:)*fac(zdim, IM) ) + &
              fac(xdim, IM)*fac(ydim, IM)*( f(ijkp(xdim)-dom%D_(xdim),ijkp(ydim)-dom%D_(ydim),ijkp(zdim)+dom%D_(zdim),:)*fac(zdim, IP) + &
                                            f(ijkp(xdim)-dom%D_(xdim),ijkp(ydim)-dom%D_(ydim),ijkp(zdim)-dom%D_(zdim),:)*fac(zdim, IM) ) + &
              fac(xdim, I0)*fac(ydim, IP)*( f(ijkp(xdim)             ,ijkp(ydim)+dom%D_(ydim),ijkp(zdim)+dom%D_(zdim),:)*fac(zdim, IP) + &
                                            f(ijkp(xdim)             ,ijkp(ydim)+dom%D_(ydim),ijkp(zdim)-dom%D_(zdim),:)*fac(zdim, IM) ) + &
              fac(xdim, I0)*fac(ydim, IM)*( f(ijkp(xdim)             ,ijkp(ydim)-dom%D_(ydim),ijkp(zdim)+dom%D_(zdim),:)*fac(zdim, IP) + &
                                            f(ijkp(xdim)             ,ijkp(ydim)-dom%D_(ydim),ijkp(zdim)-dom%D_(zdim),:)*fac(zdim, IM) ) + &
              fac(ydim, I0)*fac(zdim, IP)*( f(ijkp(xdim)+dom%D_(xdim),ijkp(ydim)             ,ijkp(zdim)+dom%D_(zdim),:)*fac(xdim, IP) + &
                                            f(ijkp(xdim)-dom%D_(xdim),ijkp(ydim)             ,ijkp(zdim)+dom%D_(zdim),:)*fac(xdim, IM) ) + &
              fac(ydim, I0)*fac(zdim, IM)*( f(ijkp(xdim)+dom%D_(xdim),ijkp(ydim)             ,ijkp(zdim)-dom%D_(zdim),:)*fac(xdim, IP) + &
                                            f(ijkp(xdim)-dom%D_(xdim),ijkp(ydim)             ,ijkp(zdim)-dom%D_(zdim),:)*fac(xdim, IM) ) + &
              fac(zdim, I0)*fac(xdim, IP)*( f(ijkp(xdim)+dom%D_(xdim),ijkp(ydim)+dom%D_(ydim),ijkp(zdim)             ,:)*fac(ydim, IP) + &
                                            f(ijkp(xdim)+dom%D_(xdim),ijkp(ydim)-dom%D_(ydim),ijkp(zdim)             ,:)*fac(ydim, IM) ) + &
              fac(zdim, I0)*fac(xdim, IM)*( f(ijkp(xdim)-dom%D_(xdim),ijkp(ydim)+dom%D_(ydim),ijkp(zdim)             ,:)*fac(ydim, IP) + &
                                            f(ijkp(xdim)-dom%D_(xdim),ijkp(ydim)-dom%D_(ydim),ijkp(zdim)             ,:)*fac(ydim, IM) )
    end function quadratic_spline

end module particle_maps
