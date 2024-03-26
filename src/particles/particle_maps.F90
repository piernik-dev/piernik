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
! pulled by NBODY
   implicit none

   private
   public :: set_map, map_particles, map_ngp, map_cic, map_tsc, quadratic_spline

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
            call die("[particle_maps:set_map] Interpolation scheme selector's logic in particle_pub:init_particles is broken. Go fix it!")
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

      use cg_cost_data,   only: I_PARTICLE
      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: xdim, ydim, zdim, ndims, LO, GEO_XYZ, PPP_PART
      use dataio_pub,     only: die
      use domain,         only: dom
      use particle_func,  only: ijk_of_particle, in_range
      use particle_types, only: particle
      use ppp,            only: ppp_main

      implicit none

      integer(kind=4), intent(in)       :: iv     !< index in cg%q array, where we want the particles to be projected
      real,            intent(in)       :: factor !< typically fpiG

      type(cg_list_element), pointer    :: cgl
      type(particle), pointer           :: pset
      integer(kind=4), dimension(ndims) :: ijkp
      character(len=*), parameter       :: map_label = "map_ngp"

      call ppp_main%start(map_label, PPP_PART)

      cgl => leaves%first
      do while (associated(cgl))
         call cgl%cg%costs%start

         ijkp(:) = cgl%cg%ijkse(:,LO)
         pset => cgl%cg%pset%first
         do while (associated(pset))
            associate( field => cgl%cg%q(iv)%arr, part => pset%pdata, cg => cgl%cg )

            if (dom%geometry_type /= GEO_XYZ) call die("[particle_maps:map_ngp] Unsupported geometry")
            where (dom%has_dir(:)) ijkp(:) = ijk_of_particle(part%pos, dom%edge(:,LO), cg%idl)
            if (in_range(ijkp, cg%ijkse)) field(ijkp(xdim), ijkp(ydim), ijkp(zdim)) = field(ijkp(xdim), ijkp(ydim), ijkp(zdim)) + factor * part%mass / cg%dvol

            end associate
            pset => pset%nxt
         enddo

         call cgl%cg%costs%stop(I_PARTICLE)
         cgl => cgl%nxt
      enddo

      call ppp_main%stop(map_label, PPP_PART)

   end subroutine map_ngp

   subroutine map_cic(iv, factor)

      use cg_cost_data,   only: I_PARTICLE
      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: half, xdim, ydim, zdim, ndims, LO, HI, CENTER, PPP_PART
      use domain,         only: dom
      use particle_types, only: particle
      use ppp,            only: ppp_main

      implicit none

      integer(kind=4), intent(in)      :: iv     !< index in cg%q array, where we want the particles to be projected
      real,            intent(in)      :: factor !< typically fpiG

      type(cg_list_element), pointer   :: cgl
      type(particle), pointer          :: pset
      integer(kind=4)                  :: cdim
      integer                          :: cn, i, j, k
      integer, dimension(ndims, LO:HI) :: ijkp
      integer, dimension(ndims)        :: cur_ind
      real                             :: weight
      character(len=*), parameter      :: map_label = "map_cic"

      call ppp_main%start(map_label, PPP_PART)

      cgl => leaves%first
      do while (associated(cgl))
         call cgl%cg%costs%start

         pset => cgl%cg%pset%first
         do while (associated(pset))
            associate( field => cgl%cg%q(iv)%arr, part => pset%pdata, cg => cgl%cg )

               if (part%outside) then
                  pset => pset%nxt
                  cycle
               endif

               do cdim = xdim, zdim
                  if (dom%has_dir(cdim)) then
                     cn = nint((part%pos(cdim) - dom%edge(cdim, LO)) * cg%idl(cdim) + half) + 1
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
            pset => pset%nxt
         enddo

         call cgl%cg%costs%stop(I_PARTICLE)
         cgl => cgl%nxt
      enddo

      call ppp_main%stop(map_label, PPP_PART)

   end subroutine map_cic

   subroutine map_tsc(iv, factor)

      use cg_cost_data,   only: I_PARTICLE
      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: xdim, ydim, zdim, ndims, LO, HI, IM, I0, IP, half, I_ONE, PPP_PART
      use dataio_pub,     only: die
      use domain,         only: dom
      use particle_func,  only: ijk_of_particle, l_neighb_part, r_neighb_part
      use particle_types, only: particle
      use ppp,            only: ppp_main

      implicit none

      integer(kind=4),              intent(in) :: iv     !< index in cg%q array, where we want the particles to be projected
      real,                         intent(in) :: factor !< typically fpiG

      type(cg_list_element), pointer           :: cgl
      type(particle), pointer                  :: pset
      integer                                  :: i, j, k
      integer(kind=4), dimension(ndims, IM:IP) :: ijkp
      real                                     :: weight_x, weight_xy, delta, fac
      character(len=*), parameter              :: map_label = "map_tsc"
      logical                                  :: full_span

      ! This routine is performance-critical, so every optimization matters
      if (dom%eff_dim /= ndims) call die("[particle_maps:map_tsc] Only 3D version is implemented")

      call ppp_main%start(map_label, PPP_PART)

      cgl => leaves%first
      do while (associated(cgl))
         call cgl%cg%costs%start

         ijkp(:,IM) = cgl%cg%ijkse(:,LO)
         ijkp(:,I0) = cgl%cg%ijkse(:,LO)
         ijkp(:,IP) = cgl%cg%ijkse(:,HI)
         pset => cgl%cg%pset%first
         do while (associated(pset))
            associate( field => cgl%cg%q(iv)%arr, part => pset%pdata, cg => cgl%cg )
               if (part%outside) then
                  pset => pset%nxt
                  cycle
               endif

               ! Nearly the same code as in particle_gravity::update_particle_acc_tsc
               ! The same performance constraints apply here

               ijkp(:,I0) = ijk_of_particle(part%pos, dom%edge(:,LO), cg%idl)
               ijkp(:,IM) = l_neighb_part(ijkp(:,I0), cg%lhn(:,LO))
               ijkp(:,IP) = r_neighb_part(ijkp(:,I0), cg%lhn(:,HI))

               full_span = (ijkp(zdim, IM) == ijkp(zdim, I0) - I_ONE) .and. (ijkp(zdim, IP) == ijkp(zdim, I0) + I_ONE)

               do i = ijkp(xdim, IM), ijkp(xdim, IP)

                  delta = (part%pos(xdim) - cg%x(i)) * cg%idl(xdim)
                  weight_x = merge(0.75 - delta**2, 1.125 - 1.5 * abs(delta) + half * delta**2, i == ijkp(xdim, I0))  !!! BEWARE hardcoded magic

                  do j = ijkp(ydim, IM), ijkp(ydim, IP)

                     delta = (part%pos(ydim) - cg%y(j)) * cg%idl(ydim)
                     weight_xy = weight_x * merge(0.75 - delta**2, 1.125 - 1.5 * abs(delta) + half * delta**2, j == ijkp(ydim, I0))  !!! BEWARE hardcoded magic

                     fac = factor * (part%mass / cg%dvol) * weight_xy
                     if (full_span) then  ! aggressive manual unrolling
                        k = ijkp(zdim, I0)
                        delta = (part%pos(zdim) - cg%z(k)) * cg%idl(zdim) ! delta for neighbors differs by -/+ 1 so no need to recalculate it for every k
                        field(i, j, k-1) = field(i, j, k-1) + fac * (1.625 - 1.5 * abs(delta + 1.) + delta + half * delta**2)
                        field(i, j, k  ) = field(i, j, k  ) + fac * (0.75 - delta**2)
                        field(i, j, k+1) = field(i, j, k+1) + fac * (1.625 - 1.5 * abs(delta - 1.) - delta + half * delta**2)
                     else
                        do k = ijkp(zdim, IM), ijkp(zdim, IP)
                           delta = (part%pos(zdim) - cg%z(k)) * cg%idl(zdim)
                           field(i, j, k) = field(i, j, k) + fac * merge(0.75 - delta**2, 1.125 - 1.5 * abs(delta) + half * delta**2, k == ijkp(zdim, I0))  !!! BEWARE hardcoded magic
                        enddo
                     endif
                  enddo
               enddo
            end associate
            pset => pset%nxt
         enddo

         call cgl%cg%costs%stop(I_PARTICLE)
         cgl => cgl%nxt
      enddo

      call leaves%leaf_arr3d_boundaries(iv)

      call ppp_main%stop(map_label, PPP_PART)

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

      ! Some computations can be avoided by splitting following statement into pieces
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
