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

module internal_bnd
! pulled by ANY

   implicit none

   private
   public :: arr3d_boundaries, internal_boundaries_3d, internal_boundaries_4d

contains

!> \brief A wrapper that calls internal_boundaries for 3D arrays (cg%q(:))

   subroutine internal_boundaries_3d(ind, nb, dim)

      use constants,  only: ARR

      implicit none

      integer(kind=4), intent(in) :: ind  !> index of cg%q(:) 3d array
      integer, optional, intent(in) :: nb !> number of grid cells to exchange (not implemented for comm3d)
      integer(kind=4), optional, intent(in) :: dim  !> do the internal boundaries only in the specified dimension

      call internal_boundaries(ind, .true., ARR, nb, dim)

   end subroutine internal_boundaries_3d

!> \brief A wrapper that calls internal_boundaries for 4D arrays (cg%u, cg%b, cg%w(:))

   subroutine internal_boundaries_4d(type, nb, dim)

      use constants,  only: FLUID, MAG, CR, INT4, fluid_n, mag_n, wcr_n
      use dataio_pub, only: die
      use grid,       only: all_cg

      implicit none

      integer(kind=4), intent(in) :: type !> FLUID, MAG, CR \todo put all of them into cg%w(:)
      integer, optional, intent(in) :: nb !> number of grid cells to exchange (not implemented for comm3d)
      integer(kind=4), optional, intent(in) :: dim  !> do the internal boundaries only in the specified dimension

      integer :: ind

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
      use domain,     only: has_dir, cdd, dom
      use gc_list,    only: cg_list_element
      use grid,       only: all_cg
      use grid_cont,  only: grid_container
      use mpi,        only: MPI_COMM_NULL
      use mpisetup,   only: comm, ierr, req, status

      implicit none

      integer(kind=4), intent(in) :: ind  !> index of cg%q(:) 3d array or cg%w(:) 4d array
      logical, intent(in)         :: tgt3d !> .true. for ARR
      integer(kind=4), intent(in) :: type !> FLUID, MAG, CR, ARR, second index in [io]_bnd arrays
      integer, optional, intent(in) :: nb !> number of grid cells to exchange (not implemented for comm3d)
      integer(kind=4), optional, intent(in) :: dim  !> do the internal boundaries only in the specified dimension

      integer :: g, d, n
      integer(kind=4) :: nr    !> index of first free slot in req and status arrays
      logical, dimension(xdim:zdim) :: dmask
      type(grid_container), pointer :: cg
      type(cg_list_element), pointer :: cgl
      real, pointer, dimension(:,:,:)   :: pa3d
      real, pointer, dimension(:,:,:,:) :: pa4d

      if (cdd%comm3d /= MPI_COMM_NULL) then
         call warn("[internal_bnd:internal_boundaries] comm3d is implemented somewhere else.")
         return
         ! ToDo: move comm3d variants here
      endif

      if (tgt3d .and. type /= ARR) call die("[internal_bnd:internal_boundaries] tgt3d .and. type /= ARR")

      dmask(:) = has_dir(:)
      if (present(dim)) then
         dmask(:) = .false.
         dmask(dim) = has_dir(dim)
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
               if (allocated(cg%i_bnd(d, type, n)%seg)) then
                  if (.not. allocated(cg%o_bnd(d, type, n)%seg)) call die("[internal_bnd:internal_boundaries] cg%i_bnd without cg%o_bnd")
                  if (ubound(cg%i_bnd(d, type, n)%seg(:), dim=1) /= ubound(cg%o_bnd(d, type, n)%seg(:), dim=1)) &
                       call die("[internal_bnd:internal_boundaries] cg%i_bnd differs in number of entries from cg%o_bnd")
                  do g = 1, ubound(cg%i_bnd(d, type, n)%seg(:), dim=1)
                     if (tgt3d) then
                        pa3d => cg%q(ind)%arr
                        call MPI_Irecv(pa3d, I_ONE, cg%i_bnd(d, type, n)%seg(g)%mbc, cg%i_bnd(d, type, n)%seg(g)%proc, cg%i_bnd(d, type, n)%seg(g)%tag, comm, req(nr+I_ONE), ierr)
                        call MPI_Isend(pa3d, I_ONE, cg%o_bnd(d, type, n)%seg(g)%mbc, cg%o_bnd(d, type, n)%seg(g)%proc, cg%o_bnd(d, type, n)%seg(g)%tag, comm, req(nr+I_TWO), ierr)
                     else
                        pa4d => cg%w(ind)%arr
                        call MPI_Irecv(pa4d, I_ONE, cg%i_bnd(d, type, n)%seg(g)%mbc, cg%i_bnd(d, type, n)%seg(g)%proc, cg%i_bnd(d, type, n)%seg(g)%tag, comm, req(nr+I_TWO), ierr)
                        call MPI_Isend(pa4d, I_ONE, cg%o_bnd(d, type, n)%seg(g)%mbc, cg%o_bnd(d, type, n)%seg(g)%proc, cg%o_bnd(d, type, n)%seg(g)%tag, comm, req(nr+I_ONE), ierr)
                     endif
                     nr = nr + I_TWO
                  enddo
               else
                  if (allocated(cg%o_bnd(d, type, n)%seg)) call die("[grid_container:internal_boundaries] cg%o_bnd without cg%i_bnd")
               endif
            endif
         enddo

         if (nr >  ubound(req(:), dim=1)) call die("[grid_container:internal_boundaries] nr > size(req) at exit")
         cgl => cgl%nxt
      enddo
      call MPI_Waitall(nr, req(:nr), status(:,:nr), ierr)

   end subroutine internal_boundaries

!-----------------------------------------------------------------------------
!
! This routine sets up all guardcells (internal and external) for given rank-3 arrays.
!

   subroutine arr3d_boundaries(ind, nb, area_type, dname)

      use constants,  only: ARR, xdim, ydim, zdim, LO, HI, BND, BLK, BND_PER, BND_MPI, BND_SHE, BND_COR, AT_NO_B, I_ONE
      use dataio_pub, only: die, msg
      use domain,     only: has_dir, cdd
      use gc_list,    only: cg_list_element
      use grid,       only: all_cg
      use grid_cont,  only: grid_container
      use mpi,        only: MPI_REQUEST_NULL, MPI_IN_PLACE, MPI_LOGICAL, MPI_LOR, MPI_COMM_NULL
      use mpisetup,   only: ierr, comm, proc, req, status

      implicit none

      integer(kind=4), intent(in) :: ind  !> Negative value: index of cg%q(:) 3d array
      integer, optional, intent(in) :: nb !> number of grid cells to exchange (not implemented for comm3d)
      integer(kind=4), intent(in), optional          :: area_type
      character(len=*), intent(in), optional         :: dname

      integer :: i, d, n
      integer(kind=4) :: lh
      logical :: dodie, do_permpi
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg
      real, dimension(:,:,:), pointer :: pa3d

      dodie = .false.

      !> \todo fill corners with big_float ?

      if (cdd%comm3d == MPI_COMM_NULL) then

         do_permpi = .true.
         if (present(area_type)) then
            if (area_type /= AT_NO_B) do_permpi = .false.
         endif

         if (do_permpi) call internal_boundaries_3d(ind, nb=nb)

      endif

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         n = cg%nb
         if (present(nb)) then
            n = nb
            if (n<=0 .or. n>cg%nb) call die("[internal_bnd:arr3d_boundaries] wrong number of guardcell layers")
         endif

         req(:) = MPI_REQUEST_NULL

         if (ind > ubound(cgl%cg%q(:), dim=1) .or. ind < lbound(cgl%cg%q(:), dim=1)) call die("[internal_bnd:arr3d_boundaries] wrong 3d index")
         pa3d =>cg%q(ind)%arr

         do d = xdim, zdim
            if (has_dir(d)) then
               do lh = LO, HI

                  select case (cg%bnd(d, lh))
                     case (BND_PER)
                        if (cdd%comm3d /= MPI_COMM_NULL) then
                           if (present(area_type)) then
                              if (area_type /= AT_NO_B) cycle
                           endif
                           do i = 1, ceiling(n/real(cg%n_b(d))) ! Repeating is important for domains that are narrower than their guardcells (e.g. cg%n_b(d) = 2)
                              select case (2*d+lh)
                                 case (2*xdim+LO)
                                    pa3d(1:cg%nb, :, :) = pa3d(cg%ieb:cg%ie, :, :) ! local copy is cheap (and don't occur so often in large runs) so don't boyher with the value of n
                                 case (2*ydim+LO)
                                    pa3d(:, 1:cg%nb, :) = pa3d(:, cg%jeb:cg%je, :)
                                 case (2*zdim+LO)
                                    pa3d(:, :, 1:cg%nb) = pa3d(:, :, cg%keb:cg%ke)
                                 case (2*xdim+HI)
                                    pa3d(cg%ie+1:cg%n_(xdim), :, :) = pa3d(cg%is:cg%isb, :, :)
                                 case (2*ydim+HI)
                                    pa3d(:, cg%je+1:cg%n_(ydim), :) = pa3d(:, cg%js:cg%jsb, :)
                                 case (2*zdim+HI)
                                    pa3d(:, :, cg%ke+1:cg%n_(zdim)) = pa3d(:, :, cg%ks:cg%ksb)
                              end select
                           enddo
                        endif
                     case (BND_MPI)
                        if (cdd%comm3d /= MPI_COMM_NULL) then
                           if (cdd%psize(d) > 1) then
                              call MPI_Isend(pa3d(1, 1, 1), I_ONE, cg%mbc(ARR, d, lh, BLK, n), cdd%procn(d, lh), int(2*d+(LO+HI-lh), kind=4), cdd%comm3d, req(4*(d-xdim)+1+2*(lh-LO)), ierr)
                              call MPI_Irecv(pa3d(1, 1, 1), I_ONE, cg%mbc(ARR, d, lh, BND, n), cdd%procn(d, lh), int(2*d+       lh,  kind=4), cdd%comm3d, req(4*(d-xdim)+2+2*(lh-LO)), ierr)
                           else
                              call die("[grid:arr3d_boundaries] bnd_[xyz][lr] == 'mpi' && cdd%psize([xyz]dim) <= 1")
                           endif
                        endif
                     case (BND_SHE) !> \todo move appropriate code from poissonsolver::poisson_solve or do nothing. or die until someone really needs SHEAR.
                        write(msg,*) "[grid:arr3d_boundaries] 'she' not implemented for ",dname
                        dodie = .true.
                     case (BND_COR)
                        if (present(area_type)) then
                           if (area_type /= AT_NO_B) cycle
                        endif
                        write(msg,*) "[grid:arr3d_boundaries] 'cor' not implemented for ", dname
                        dodie = .true.
                     case default ! Set gradient == 0 on the external boundaries
                        if (present(area_type)) then
                           if (area_type /= AT_NO_B) cycle
                        endif
                        do i = 1, cg%nb
                           select case (2*d+lh)
                              case (2*xdim+LO)
                                 pa3d(i, :, :) = pa3d(cg%is, :, :)
                              case (2*ydim+LO)
                                 pa3d(:, i, :) = pa3d(:, cg%js, :)
                              case (2*zdim+LO)
                                 pa3d(:, :, i) = pa3d(:, :, cg%ks)
                              case (2*xdim+HI)
                                 pa3d(cg%ie+i, :, :) = pa3d(cg%ie, :, :)
                              case (2*ydim+HI)
                                 pa3d(:, cg%je+i, :) = pa3d(:, cg%je, :)
                              case (2*zdim+HI)
                                 pa3d(:, :, cg%ke+i) = pa3d(:, :, cg%ke)
                           end select
                        enddo
                  end select

               enddo
            endif
            !> \warning outside xdim-zdim loop MPI_Waitall may change the operations order and as a result may leave mpi-corners uninitiallized
            if (cdd%comm3d /= MPI_COMM_NULL) call MPI_Waitall(size(req(:)), req(:), status(:,:), ierr)
         enddo

         cgl => cgl%nxt
      enddo

      call MPI_Allreduce(MPI_IN_PLACE, dodie, I_ONE, MPI_LOGICAL, MPI_LOR, comm, ierr)
      if (dodie) call die(msg)

   end subroutine arr3d_boundaries

end module internal_bnd
