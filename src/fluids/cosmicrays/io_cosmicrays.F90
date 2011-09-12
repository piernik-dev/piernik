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
#include "macros.h"

module io_cosmicrays
! pulled by COSM_RAYS
   implicit none

   private
#ifdef NEW_HDF5
   public :: cr_add_hdf5
#endif /* NEW_HDF5 */

contains

!>
!! \deprecated These routines will probably be completely replaced by a different approach that allows for consistent use of multiple grids per thread
!<

#ifdef NEW_HDF5
   subroutine cr_add_hdf5(flind_crs)

      use dataio_pub, only: die
      use grid,       only: all_cg
      use grid_cont,  only: grid_container
      use initcosmicrays, only: iarr_crs
      use list_hdf5,  only: add_lhdf5, lhdf5_info

      implicit none

      integer, intent(in)              :: flind_crs
      type(lhdf5_info) :: item
      integer          :: i
      type(grid_container), pointer :: cg

      cg => all_cg%first%cg
      if (all_cg%cnt > 1) call die("[io_cosmicrays:cr_add_hdf5] multiple grid pieces per procesor not implemented yet") !nontrivial add_lhdf5(item)

      item%p    => get_cr
      if (.not.allocated(item%ivec)) allocate(item%ivec(10))
      if (.not.allocated(item%rvec)) allocate(item%rvec(0))

      do i = 1, flind_crs

         write(item%key,'(A,I1)')  "cr",i
         item%ivec  = [cg%n_b(:), cg%is, cg%ie, cg%js, cg%je, cg%ks, cg%ke, iarr_crs(i)]
         call add_lhdf5(item)

      enddo
   end subroutine cr_add_hdf5

   subroutine get_cr(ivec,rvec,outtab)

      use dataio_pub, only: die
      use grid,       only: all_cg
      use grid_cont,  only: grid_container
      use dataio_pub, only: die

      implicit none

      integer, dimension(:), intent(in)  :: ivec
      real,    dimension(:), intent(in)  :: rvec
      real, dimension(:,:,:), allocatable, intent(out) :: outtab
      type(grid_container), pointer :: cg

      cg => all_cg%first%cg
      if (all_cg%cnt > 1) call die("[io_cosmicrays:cr_add_hdf5] multiple grid pieces per procesor not implemented yet") !nontrivial

      if (allocated(outtab)) call die("[io_cosmicrays:get_cr]: outtab already allocated")
      allocate(outtab(ivec(1),ivec(2),ivec(3)))
      outtab(:,:,:) = cg%u%arr(ivec(10),ivec(4):ivec(5),ivec(6):ivec(7),ivec(8):ivec(9))

      return

      if (.false.) write(0,*) rvec

   end subroutine get_cr

#endif /* NEW_HDF5 */

end module io_cosmicrays
