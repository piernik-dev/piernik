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
!! \brief Update internal boundaries
!!
!! \detailed This module contains subroutines that are responsible for preparing internal boundary cells for cg%q(:) and cg%w(:) arrays.
!! For simplicity, let's assume, that all boundaries, that rely on MPI communication are internal.
!! This implies that periodic, corner shear and fine-coarse boundaries are also "internal"
!!
!! Note that this routine may not properly update some layers of guardcells when number of guardcell layers exceedes number of active cells.
!! Appropriate checks should be made in divide_domain routine.
!!
!! \todo integrate here as much stuff from fluidboundaries, magboundaries, etc.  as possible.
!!
!<

module internal_bnd
! pulled by ANY

   implicit none

   private
   public :: internal_boundaries_3d, internal_boundaries_4d, internal_boundaries

contains

!> \brief A wrapper that calls internal_boundaries for 3D arrays (cg%q(:))

   subroutine internal_boundaries_3d(ind, nb, dim)

      use constants,  only: ARR

      implicit none

      integer(kind=4),           intent(in) :: ind  !> index of cg%q(:) 3d array
      integer,         optional, intent(in) :: nb   !> number of grid cells to exchange (not implemented for comm3d)
      integer(kind=4), optional, intent(in) :: dim  !> do the internal boundaries only in the specified dimension

      call internal_boundaries(ind, .true., ARR, nb, dim)

   end subroutine internal_boundaries_3d

!> \brief A wrapper that calls internal_boundaries for 4D arrays (cg%u, cg%b, cg%w(:))

   subroutine internal_boundaries_4d(type, nb, dim)

      use constants,  only: FLUID, MAG, CR, INT4, fluid_n, mag_n, wcr_n
      use dataio_pub, only: die
      use grid,       only: all_cg

      implicit none

      integer(kind=4),           intent(in) :: type !> FLUID, MAG, CR \todo put all of them into cg%w(:)
      integer,         optional, intent(in) :: nb   !> number of grid cells to exchange (not implemented for comm3d)
      integer(kind=4), optional, intent(in) :: dim  !> do the internal boundaries only in the specified dimension

      integer(kind=4) :: ind

      select case (type)
         case (FLUID)
            ind = all_cg%first%cg%get_na_ind_4d(fluid_n)
         case (MAG)
            ind = all_cg%first%cg%get_na_ind_4d(mag_n)
         case (CR)
            ind = all_cg%first%cg%get_na_ind_4d(wcr_n)
         case default
            call die("[internal_bnd:internal_boundaries_4d] What?")
            ind = 0_INT4 ! suppress compiler warnings
      end select

      call internal_boundaries(ind, .false., type, nb, dim)

   end subroutine internal_boundaries_4d

!>
!! \brief This routine exchanges guardcells for BND_MPI and BND_PER boundaries on rank-3 and rank-4 arrays
!! \details This routine should not be called directly. Appropriate wrappers for rank-3 and rank-4 arrays are provided above.
!! The corners should be properly updated if this%[io]_bnd(:, ind) was set up appropriately and this routine is called separately for each dimension.
!!
!! \todo Check how much performance is lost due to using MPI calls even for local copies. Decide whether it is worth to convert local MPI calls to direct memory copies.
!<

   subroutine internal_boundaries(ind, tgt3d, type, nb, dim)

      use constants,  only: FLUID, MAG, CR, ARR, xdim, zdim, I_ONE, I_TWO
      use dataio_pub, only: die, warn
      use domain,     only: cdd, dom
      use gc_list,    only: cg_list_element
      use grid,       only: all_cg
      use grid_cont,  only: grid_container
      use mpi,        only: MPI_COMM_NULL
      use mpisetup,   only: comm, ierr, req, status

      implicit none

      integer(kind=4),           intent(in) :: ind    !> index of cg%q(:) 3d array or cg%w(:) 4d array
      logical,                   intent(in) :: tgt3d  !> .true. for ARR
      integer(kind=4),           intent(in) :: type   !> FLUID, MAG, CR, ARR, second index in [io]_bnd arrays
      integer,         optional, intent(in) :: nb     !> number of grid cells to exchange (not implemented for comm3d)
      integer(kind=4), optional, intent(in) :: dim    !> do the internal boundaries only in the specified dimension

      integer                               :: g, d, n
      integer(kind=4)                       :: nr     !> index of first free slot in req and status arrays
      logical, dimension(xdim:zdim)         :: dmask
      type(grid_container),     pointer     :: cg
      type(cg_list_element),    pointer     :: cgl
      real, dimension(:,:,:),   pointer     :: pa3d
      real, dimension(:,:,:,:), pointer     :: pa4d

      if (cdd%comm3d /= MPI_COMM_NULL) then
         call warn("[internal_bnd:internal_boundaries] comm3d is implemented somewhere else.")
         return
         ! ToDo: move comm3d variants here
      endif

      if (tgt3d .and. type /= ARR) call die("[internal_bnd:internal_boundaries] tgt3d .and. type /= ARR")

      dmask(:) = dom%has_dir(:)
      if (present(dim)) then
         dmask(:) = .false.
         dmask(dim) = dom%has_dir(dim)
      endif

      n = dom%nb
      if (present(nb)) then
         n = nb
         if (n<=0 .or. n>dom%nb) call die("[internal_bnd:internal_boundaries] wrong number of guardcell layers")
      endif

      nr = 0
      cgl => all_cg%first
      if (tgt3d) then
         if (ind > ubound(cgl%cg%q(:), dim=1) .or. ind < lbound(cgl%cg%q(:), dim=1)) call die("[internal_bnd:internal_boundaries] wrong 3d index")
      else
         if (ind > ubound(cgl%cg%w(:), dim=1) .or. ind < lbound(cgl%cg%w(:), dim=1)) call die("[internal_bnd:internal_boundaries] wrong 4d index")
      endif
      do while (associated(cgl))
         cg => cgl%cg

         if (tgt3d) then
            if (cg%q(ind)%name /= all_cg%first%cg%q(ind)%name) call die("[internal_bnd:internal_boundaries] 3d array name mismatch")
         else
            if (cg%w(ind)%name /= all_cg%first%cg%w(ind)%name) call die("[internal_bnd:internal_boundaries] 4d array name mismatch")
         endif

         do d = xdim, zdim
            if (dmask(d)) then
               if (allocated(cg%i_bnd(d, n)%seg)) then
                  if (.not. allocated(cg%o_bnd(d, n)%seg)) call die("[internal_bnd:internal_boundaries] cg%i_bnd without cg%o_bnd")
                  if (ubound(cg%i_bnd(d, n)%seg(:), dim=1) /= ubound(cg%o_bnd(d, n)%seg(:), dim=1)) &
                       call die("[internal_bnd:internal_boundaries] cg%i_bnd differs in number of entries from cg%o_bnd")
                  do g = 1, ubound(cg%i_bnd(d, n)%seg(:), dim=1)
                     if (tgt3d) then
                        pa3d => cg%q(ind)%arr
                        call MPI_Irecv(pa3d, I_ONE, cg%q_i_mbc(d, n)%mbc(g), cg%i_bnd(d, n)%seg(g)%proc, cg%i_bnd(d, n)%seg(g)%tag, comm, req(nr+I_ONE), ierr)
                        call MPI_Isend(pa3d, I_ONE, cg%q_o_mbc(d, n)%mbc(g), cg%o_bnd(d, n)%seg(g)%proc, cg%o_bnd(d, n)%seg(g)%tag, comm, req(nr+I_TWO), ierr)
                     else
                        pa4d => cg%w(ind)%arr
                        call MPI_Irecv(pa4d, I_ONE, cg%w(ind)%w_i_mbc(d, n)%mbc(g), cg%i_bnd(d, n)%seg(g)%proc, cg%i_bnd(d, n)%seg(g)%tag, comm, req(nr+I_TWO), ierr)
                        call MPI_Isend(pa4d, I_ONE, cg%w(ind)%w_o_mbc(d, n)%mbc(g), cg%o_bnd(d, n)%seg(g)%proc, cg%o_bnd(d, n)%seg(g)%tag, comm, req(nr+I_ONE), ierr)
                     endif
                     nr = nr + I_TWO
                  enddo
               else
                  if (allocated(cg%o_bnd(d, n)%seg)) call die("[grid_container:internal_boundaries] cg%o_bnd without cg%i_bnd")
               endif
            endif
         enddo

         if (nr >  ubound(req(:), dim=1)) call die("[grid_container:internal_boundaries] nr > size(req) at exit")
         cgl => cgl%nxt
      enddo
      call MPI_Waitall(nr, req(:nr), status(:,:nr), ierr)

   end subroutine internal_boundaries

end module internal_bnd
