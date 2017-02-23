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

!> \brief Module containing the grid container type and its associated methods related to operations on fluid data

module grid_cont_dataop

   use grid_cont_basic, only: grid_container_basic

   implicit none

   private
   public :: grid_container_dataop, get_cs

   !> \brief Everything required for autonomous computation of a single sweep on a portion of the domain on a single process
   type, extends(grid_container_basic) :: grid_container_dataop
   contains

      procedure       :: get_cs                               !< calculate the speed of sound at given position
      procedure       :: set_constant_b_field                 !< set constant magnetic field on whole block
      procedure       :: emag_point                           !< return energy asociated with magnetic field at specified point
      procedure       :: emag_range                           !< return energy asociated with magnetic field at specified range
      generic, public :: emag => emag_point, emag_range

   end type grid_container_dataop

contains

!< \brief set constant magnetic field on whole block

   subroutine set_constant_b_field(this, b)

      use constants, only: xdim, ydim, zdim, I_ONE

      implicit none

      class(grid_container_dataop), intent(inout) :: this !< object invoking type-bound procedure
      real, dimension(xdim:zdim),   intent(in)    :: b    !< the value of the magnetic field vector in whole block

      integer :: d

      if (associated(this%b)) then
         do d = xdim, zdim
            this%b(d, this%is:this%ie, this%js:this%je, this%ks:this%ke) = b(d)
         enddo
      endif

      if (associated(this%bf)) then
         this%bf(xdim)%arr(this%is:this%ie+I_ONE, this%js:this%je,       this%ks:this%ke      ) = b(xdim)
         this%bf(ydim)%arr(this%is:this%ie,       this%js:this%je+I_ONE, this%ks:this%ke      ) = b(ydim)
         this%bf(zdim)%arr(this%is:this%ie,       this%js:this%je,       this%ks:this%ke+I_ONE) = b(zdim)
      endif

   end subroutine set_constant_b_field

!< \brief return energy asociated with magnetic field at specified point

   function emag_point(this, ijk) result(e_mag)

      use constants,  only: xdim, ydim, zdim, I_ONE
      use dataio_pub, only: die
      use func,       only: emag

      implicit none

      class(grid_container_dataop), intent(in)  :: this
      integer, dimension(:),        intent(in)  :: ijk

      real :: e_mag

      if (associated(this%b)) then
         e_mag = emag(this%b(xdim, ijk(xdim), ijk(ydim), ijk(zdim)), &
              &       this%b(ydim, ijk(xdim), ijk(ydim), ijk(zdim)), &
              &       this%b(zdim, ijk(xdim), ijk(ydim), ijk(zdim)))
      else if (associated(this%bf)) then
         e_mag = emag((this%bf(xdim)%arr(ijk(xdim), ijk(ydim), ijk(zdim)) + this%bf(xdim)%arr(ijk(xdim)+I_ONE, ijk(ydim), ijk(zdim)))/2., &
              &       (this%bf(ydim)%arr(ijk(xdim), ijk(ydim), ijk(zdim)) + this%bf(ydim)%arr(ijk(xdim), ijk(ydim)+I_ONE, ijk(zdim)))/2., &
              &       (this%bf(zdim)%arr(ijk(xdim), ijk(ydim), ijk(zdim)) + this%bf(zdim)%arr(ijk(xdim), ijk(ydim), ijk(zdim)+I_ONE))/2.)
      else
         call die("[grid_container:emag_point] no magnetic field declared here")
         e_mag = 0.
      endif

   end function emag_point

!< \brief return energy asociated with magnetic field at specified range

   function emag_range(this, ijk) result(e_mag)

      use constants,  only: xdim, ydim, zdim, LO, HI, I_ONE
      use dataio_pub, only: die
      use func,       only: emag

      implicit none

      class(grid_container_dataop), intent(in)  :: this
      integer, dimension(:,:),      intent(in)  :: ijk

      real, dimension(ijk(xdim, LO):ijk(xdim, HI), ijk(ydim, LO):ijk(ydim, HI), ijk(zdim, LO):ijk(zdim, HI)) :: e_mag

      if (associated(this%b)) then
         e_mag = emag(this%b(xdim, ijk(xdim, LO):ijk(xdim, HI), ijk(ydim, LO):ijk(ydim, HI), ijk(zdim, LO):ijk(zdim, HI)), &
              &       this%b(ydim, ijk(xdim, LO):ijk(xdim, HI), ijk(ydim, LO):ijk(ydim, HI), ijk(zdim, LO):ijk(zdim, HI)), &
              &       this%b(zdim, ijk(xdim, LO):ijk(xdim, HI), ijk(ydim, LO):ijk(ydim, HI), ijk(zdim, LO):ijk(zdim, HI)))
      else if (associated(this%bf)) then
         e_mag = emag((this%bf(xdim)%arr(ijk(xdim, LO)      :ijk(xdim, HI),       ijk(ydim, LO):ijk(ydim, HI), ijk(zdim, LO):ijk(zdim, HI)) + &
              &        this%bf(xdim)%arr(ijk(xdim, LO)+I_ONE:ijk(xdim, HI)+I_ONE, ijk(ydim, LO):ijk(ydim, HI), ijk(zdim, LO):ijk(zdim, HI)))/2., &
              &       (this%bf(ydim)%arr(ijk(xdim, LO):ijk(xdim, HI), ijk(ydim, LO)      :ijk(ydim, HI),       ijk(zdim, LO):ijk(zdim, HI)) + &
              &        this%bf(ydim)%arr(ijk(xdim, LO):ijk(xdim, HI), ijk(ydim, LO)+I_ONE:ijk(ydim, HI)+I_ONE, ijk(zdim, LO):ijk(zdim, HI)))/2., &
              &       (this%bf(zdim)%arr(ijk(xdim, LO):ijk(xdim, HI), ijk(ydim, LO):ijk(ydim, HI), ijk(zdim, LO)      :ijk(zdim, HI)      ) + &
              &        this%bf(zdim)%arr(ijk(xdim, LO):ijk(xdim, HI), ijk(ydim, LO):ijk(ydim, HI), ijk(zdim, LO)+I_ONE:ijk(zdim, HI)+I_ONE))/2.)
      else
         call die("[grid_container:emag_range] no magnetic field declared here")
         e_mag = 0.
      endif

   end function emag_range

!< \brief calculate the speed of sound at given position

   real function get_cs(this, fl, i, j, k) ! ion_cs

      use constants,  only: zero, two, ION, NEU, DST
      use dataio_pub, only: die
      use fluidtypes, only: component_fluid
#ifndef ISO
      use func,       only: ekin
#endif /* !ISO */
#ifdef MAGNETIC
      use constants,  only: xdim, ydim, zdim, half
      use domain,     only: dom
      use func,       only: emag
#endif /* !MAGNETIC */

      implicit none

      class(grid_container_dataop), intent(in) :: this
      class(component_fluid),       intent(in) :: fl
      integer,                      intent(in) :: i, j, k

#ifdef MAGNETIC
      real :: bx, by, bz
#endif /* MAGNETIC */
      real :: pmag, p, ps

      select case (fl%tag)
         case (DST)
            get_cs = zero
         case (ION)
#ifdef MAGNETIC
            bx = half*(this%b(xdim,i,j,k) + this%b(xdim, i+dom%D_x, j,         k        ))
            by = half*(this%b(ydim,i,j,k) + this%b(ydim, i,         j+dom%D_y, k        ))
            bz = half*(this%b(zdim,i,j,k) + this%b(zdim, i,         j,         k+dom%D_z))

            pmag = emag(bx, by, bz)
#else /* !MAGNETIC */
            ! all_mag_boundaries has not been called so we cannot trust b(xdim, ie+dom%D_x:), b(ydim,:je+dom%D_y and b(zdim,:,:, ke+dom%D_z
            pmag = zero
#endif /* !MAGNETIC */
#ifdef ISO
            p  = this%cs_iso2(i, j, k) * this%u(fl%idn, i, j, k)
            ps = p + pmag
            get_cs = sqrt(abs((two * pmag + p) / this%u(fl%idn, i, j, k)))
#else /* !ISO */
            ps = (this%u(fl%ien, i, j, k) - &
                 &   ekin(this%u(fl%imx, i, j, k), this%u(fl%imy, i, j, k), this%u(fl%imz, i, j, k), this%u(fl%idn, i, j, k)) &
                 & ) * fl%gam_1 + (two - fl%gam) * pmag
            p  = ps - pmag
            get_cs = sqrt(abs((two * pmag + fl%gam * p) / this%u(fl%idn, i, j, k)))
#endif /* !ISO */

         case (NEU)
#ifdef ISO
            get_cs = sqrt(this%cs_iso2(i, j, k))
#else /* !ISO */
            p = (this%u(fl%ien, i, j, k) - &
                 &   ekin(this%u(fl%imx, i, j, k), this%u(fl%imy, i, j, k), this%u(fl%imz, i, j, k), this%u(fl%idn, i, j, k)) &
                 & ) * fl%gam_1
            get_cs = sqrt(abs((fl%gam * p) / this%u(fl%idn, i, j, k)))
#endif /* !ISO */
         case default
            call die("[grid_container_dataop:get_cs] unknown fluid (!ISO)")
            get_cs = zero
      end select

   end function get_cs

end module grid_cont_dataop
