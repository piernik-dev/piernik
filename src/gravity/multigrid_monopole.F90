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
!! \brief Monopole solver for isolated boundaries
!<

module mg_monopole
! pulled by MULTIGRID && SELF_GRAV

   ! needed for global vars in this module
   use constants, only: ndims, xdim

   implicit none

   private
   public :: isolated_monopole, find_img_CoM

   integer, parameter           :: imass = xdim - 1    !< index for mass in CoM(:)
   real, dimension(imass:ndims) :: CoM                 !< Total mass and center of mass coordinates

contains

!>
!! \brief Set boundary potential from monopole source. Fill cg%mg%bnd_[xyz] arrays in leaves with expected values of the gravitational potential at external face of computational domain.
!! \details This is a simplified approach that can be used for tests and as a fast replacement for the
!! multipole boundary solver for nearly spherically symmetric source distributions.
!! The isolated_monopole subroutine ignores the radial profile of the monopole
!<

   subroutine isolated_monopole

      use cg_list,      only: cg_list_element
      use cg_leaves,    only: leaves
      use constants,    only: xdim, ydim, zdim, LO, HI, GEO_XYZ !, GEO_RPZ
      use dataio_pub,   only: die
      use domain,       only: dom
      use grid_cont,    only: grid_container
      use units,        only: newtong
#ifdef NBODY
      use dataio_pub,     only: msg, warn
      use mpisetup,       only: proc
      use particle_types, only: particle
#endif /* NBODY */

      implicit none

      integer :: i, j, k, lh
      real    :: r2
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
#ifdef NBODY
      type(particle), pointer    :: pset
#endif /* NBODY */

      if (dom%geometry_type /= GEO_XYZ) call die("[multigrid_monopole:isolated_monopole] non-cartesian geometry not implemented yet")

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         do lh = LO, HI
            if (cg%ext_bnd(xdim, lh)) then
               do j = cg%js, cg%je
                  do k = cg%ks, cg%ke
                     r2 = (cg%y(j)-CoM(ydim))**2 + (cg%z(k) - CoM(zdim))**2
                     cg%mg%bnd_x(j, k, lh) = - newtong * CoM(imass) / sqrt(r2 + (cg%fbnd(xdim, lh)-CoM(xdim))**2)
                  enddo
               enddo
            endif
            if (cg%ext_bnd(ydim, lh)) then
               do i = cg%is, cg%ie
                  do k = cg%ks, cg%ke
                     r2 = (cg%x(i)-CoM(xdim))**2 + (cg%z(k) - CoM(zdim))**2
                     cg%mg%bnd_y(i, k, lh) = - newtong * CoM(imass) / sqrt(r2 + (cg%fbnd(ydim, lh)-CoM(ydim))**2)
                  enddo
               enddo
            endif
            if (cg%ext_bnd(zdim, lh)) then
               do i = cg%is, cg%ie
                  do j = cg%js, cg%je
                     r2 = (cg%x(i)-CoM(xdim))**2 + (cg%y(j) - CoM(ydim))**2
                     cg%mg%bnd_z(i, j, lh) = - newtong * CoM(imass) / sqrt(r2 + (cg%fbnd(zdim, lh)-CoM(zdim))**2)
                  enddo
               enddo
            endif
         enddo
         cgl => cgl%nxt

#ifdef NBODY
         pset => cg%pset%first
         do while (associated(pset))
            if (pset%pdata%outside) then
               write(msg, '(a,i8,a,i5,a)')"[multigrid_monopole:isolated_monopole] Particle #", i, " on process ", proc, "ignored"
               call warn(msg)
            endif
            pset => pset%nxt
         enddo
#endif /* NBODY */

      enddo

   end subroutine isolated_monopole

!>
!! \brief Find total mass and its center
!!
!! \details This routine does the summation only on external boundaries
!<

   subroutine find_img_CoM

      use allreduce,    only: piernik_MPI_Allreduce
      use cg_leaves,    only: leaves
      use cg_list,      only: cg_list_element
      use constants,    only: ndims, xdim, ydim, zdim, LO, HI, GEO_XYZ, pSUM, zero, V_DEBUG !, GEO_RPZ
      use dataio_pub,   only: die, msg, printinfo
      use domain,       only: dom
      use func,         only: operator(.notequals.)
      use grid_cont,    only: grid_container
      use mpisetup,     only: master
      use units,        only: fpiG
#ifdef NBODY
      use particle_types, only: particle
#endif /* NBODY */

      implicit none

      real, dimension(imass:ndims)   :: lsum, dsum
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer  :: cg
      integer                        :: lh, d
#ifdef NBODY
      type(particle), pointer    :: pset
#endif /* NBODY */

      if (dom%geometry_type /= GEO_XYZ) call die("[multigrid_monopole:find_img_CoM] non-cartesian geometry not implemented yet")

      lsum(:) = 0.

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         do lh = LO, HI
            if (cg%ext_bnd(xdim, lh)) then
               d = cg%ijkse(xdim, lh)
               dsum(imass)     =        sum( cg%mg%bnd_x(cg%js:cg%je, cg%ks:cg%ke, lh), mask=cg%leafmap(d, :, :) )
               dsum(xdim:zdim) = [ dsum(imass) * cg%fbnd(xdim, lh), &
                    &              sum( sum( cg%mg%bnd_x(cg%js:cg%je, cg%ks:cg%ke, lh), mask=cg%leafmap(d, :, :), dim=2) * cg%y(cg%js:cg%je) ), &
                    &              sum( sum( cg%mg%bnd_x(cg%js:cg%je, cg%ks:cg%ke, lh), mask=cg%leafmap(d, :, :), dim=1) * cg%z(cg%ks:cg%ke) ) ]
               lsum(:)         = lsum(:) + dsum(:) * cg%dyz
            endif
            if (cg%ext_bnd(ydim, lh)) then
               d = cg%ijkse(ydim, lh)
               dsum(imass)     =        sum( cg%mg%bnd_y(cg%is:cg%ie, cg%ks:cg%ke, lh), mask=cg%leafmap(:, d, :) )
               dsum(xdim:zdim) = [ sum( sum( cg%mg%bnd_y(cg%is:cg%ie, cg%ks:cg%ke, lh), mask=cg%leafmap(:, d, :), dim=2) * cg%x(cg%is:cg%ie) ), &
                    &              dsum(imass) * cg%fbnd(ydim, lh), &
                    &              sum( sum( cg%mg%bnd_y(cg%is:cg%ie, cg%ks:cg%ke, lh), mask=cg%leafmap(:, d, :), dim=1) * cg%z(cg%ks:cg%ke) ) ]
               lsum(:)         = lsum(:) + dsum(:) * cg%dxz
            endif
            if (cg%ext_bnd(zdim, lh)) then
               d = cg%ijkse(zdim, lh)
               dsum(imass)     =        sum( cg%mg%bnd_z(cg%is:cg%ie, cg%js:cg%je, lh), mask=cg%leafmap(:, :, d) )
               dsum(xdim:zdim) = [ sum( sum( cg%mg%bnd_z(cg%is:cg%ie, cg%js:cg%je, lh), mask=cg%leafmap(:, :, d), dim=2) * cg%x(cg%is:cg%ie) ), &
                    &              sum( sum( cg%mg%bnd_z(cg%is:cg%ie, cg%js:cg%je, lh), mask=cg%leafmap(:, :, d), dim=1) * cg%y(cg%js:cg%je) ), &
                    &              dsum(imass) * cg%fbnd(zdim, lh) ]
               lsum(:)         = lsum(:) + dsum(:) * cg%dxy
            endif
         enddo

         ! Add only those particles, which are placed outside the domain. Particles inside the domain were already mapped on the grid.
         !> \warning Do we need to use the fppiG factor here?
#ifdef NBODY
         pset => cg%pset%first
         do while (associated(pset))
            if (pset%pdata%outside) lsum(:) = lsum(:) + [ pset%pdata%mass, pset%pdata%mass * pset%pdata%pos(:) ]
            pset => pset%nxt
         enddo
#endif /* NBODY */

         cgl => cgl%nxt
      enddo

      CoM(imass:ndims) = lsum(imass:ndims)
      call piernik_MPI_Allreduce(CoM(imass:ndims), pSUM)

      if (CoM(imass).notequals.zero) then
         CoM(xdim:zdim) = CoM(xdim:zdim) / CoM(imass)
      else
         call die("[multigrid_monopole:find_img_CoM] Total mass == 0")
      endif
      if (master) then
         write(msg, '(a,g14.6,a,3g14.6,a)')"[multigrid_monopole:find_img_CoM] Total mass = ", CoM(imass)/fpiG," at (",CoM(xdim:zdim),")"
         call printinfo(msg, V_DEBUG)
      endif

   end subroutine find_img_CoM

end module mg_monopole
