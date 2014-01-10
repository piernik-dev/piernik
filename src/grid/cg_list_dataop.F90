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

!> \brief This module contains grid container list and some methods to manipulate data contained in cg

module cg_list_dataop

   use cg_list, only: cg_list_T

   implicit none

   private
   public :: cg_list_dataop_T, ind_val, dirty_label, expanded_domain

   !> \brief Arbitrary list of grid containers
   type, extends(cg_list_T) :: cg_list_dataop_T
    contains

      ! Misc
      procedure :: get_extremum                      !< Find minimum or maximum value over a s list
      procedure :: set_dirty                         !< Pollute selected array with an insane value dirtyH.
      procedure :: check_dirty                       !< Check for detectable traces of set_dirty calls.
      procedure :: check_all_dirty                   !< Check all registered arrays for detectable traces of set_dirty calls.
      procedure :: check_for_dirt                    !< Check all named arrays for constants:big_float
      procedure :: corners2wa                        !< Cut out a cross in the middle of the grid block and move the data to show what is stored in corners and partially in edges and faces as well

      ! Arithmetic on the fields
      procedure :: set_q_value                       !< reset given field to the value
      procedure :: q_copy                            !< copy a given field to another
      procedure :: qw_copy                           !< copy a given rank-3 field to a component of rank-4 field
      procedure :: wq_copy                           !< copy a component of rank-4 field to a given rank-3 field
      procedure :: q_add                             !< add a field to another
      procedure :: q_add_val                         !< add a value to a field
      procedure :: q_lin_comb                        !< assign linear combination of q fields
      procedure :: subtract_average                  !< subtract average value from the list
      procedure :: norm_sq                           !< calculate L2 norm
      procedure :: mul_by_r                          !< multiply data by the r (first) coordinate
      procedure :: div_by_r                          !< divide data by the r (first) coordinate
      procedure :: scalar_product                    !< return scalar product of the given fields

      ! Multigrid
      generic,   public  :: reset_boundaries => zero_boundaries, dirty_mg_boundaries
      procedure, private :: zero_boundaries                   !< Clear boundary values
      procedure, private :: dirty_mg_boundaries               !< Set boundary values
!> \todo merge lists

   end type cg_list_dataop_T

   !> \brief Index - value pairs for calling arithmetic on the grids with q_lin_comb
   type ind_val
      integer :: ind  !< index in cg%q
      real    :: val  !< value for multiplication
   end type ind_val

   integer, parameter    :: dl_len = 64 !< length of label buffer
   character(len=dl_len) :: dirty_label !< buffer for label for check_dirty subroutine

   type(cg_list_dataop_T):: expanded_domain !< grid pieces that were created in the area, where computational domain was recently expanded

contains

!>
!! \brief Find minimum or maximum value over a specified list of grid containers
!!
!! \details It should be possible to find an extremum over a given level or leaf blocks or something
!<
   subroutine get_extremum(this, ind, minmax, prop, dir)

      use cg_list,     only: cg_list_element
      use constants,   only: MINL, MAXL, I_ONE, ndims, xdim, ydim, zdim, big_float, LO
      use dataio_pub,  only: msg, warn, die
      use domain,      only: dom
      use grid_cont,   only: grid_container
      use mpi,         only: MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_STATUS_IGNORE, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_MAXLOC, MPI_IN_PLACE
      use mpisetup,    only: comm, mpi_err, master, proc, FIRST
!      use named_array_list, only: qna
      use types,       only: value

      implicit none

      class(cg_list_dataop_T),   intent(in)    :: this    !< object invoking type-bound procedure
      integer(kind=4),           intent(in)    :: ind     !< Index in cg%q(:)
      integer(kind=4),           intent(in)    :: minmax  !< minimum or maximum ?
      type(value),               intent(out)   :: prop    !< precise location of the extremum to be found
      integer(kind=4), optional, intent(in)    :: dir     !< order the cell size in dir direction


      type(grid_container),   pointer          :: cg_x
      type(cg_list_element),  pointer          :: cgl
      real, dimension(:,:,:), pointer          :: tab
      integer,                       parameter :: tag1 = 11, tag2 = tag1 + 1, tag3 = tag2 + 1
      integer, dimension(MINL:MAXL), parameter :: op = [ MPI_MINLOC, MPI_MAXLOC ]
      enum, bind(C)
         enumerator :: I_V, I_P !< value and proc
      end enum
      real, dimension(I_V:I_P)  :: v_red

      if (.not. present(dir)) prop%assoc = - big_float ! means uninitialized

      prop%loc(:) = 0
      select case (minmax)
         case (MINL)
            prop%val = huge(1.)
         case (MAXL)
            prop%val = -huge(1.)
         case default
            prop%val = 0.0   !! set dummy value
            write(msg,*) "[cg_list_dataop:get_extremum]: I don't know what to do with minmax = ", minmax
            call warn(msg)
      end select

      nullify(cg_x)
      if (associated(this%first)) then
         cgl => this%first
         if (ind > ubound(cgl%cg%q(:), dim=1) .or. ind < lbound(cgl%cg%q(:), dim=1)) call die("[cg_list_dataop:get_extremum] Wrong index")
         do while (associated(cgl))

            tab => cgl%cg%q(ind)%span(cgl%cg%ijkse)
            select case (minmax)
               case (MINL)
                  if (minval(tab) < prop%val) then
                     prop%val = minval(tab)
                     prop%loc = minloc(tab) + cgl%cg%ijkse(:, LO) - 1
                     ! it isn't too much intiutive that minloc and maxloc return values as the array indexing start from 1, but here tab does so anyway
                     cg_x => cgl%cg
                  endif
               case (MAXL)
                  if (maxval(tab) > prop%val) then
                     prop%val = maxval(tab)
                     prop%loc = maxloc(tab) + cgl%cg%ijkse(:, LO) - 1
                     cg_x => cgl%cg
                  endif
            end select
            cgl => cgl%nxt
         enddo
         v_red(I_V) = prop%val;
      else
         v_red(I_V) = -prop%val   !! No cgl means nothing to do there
      endif

      v_red(I_P) = real(proc)

      call MPI_Allreduce(MPI_IN_PLACE, v_red, I_ONE, MPI_2DOUBLE_PRECISION, op(minmax), comm, mpi_err)

      prop%val = v_red(I_V)
      prop%proc = int(v_red(I_P))

      if (proc == prop%proc) then
         where (.not. dom%has_dir(:)) prop%coords(:) = 0.
         if (associated(cg_x)) then
            if (dom%has_dir(xdim)) prop%coords(xdim) = cg_x%x(prop%loc(xdim))
            if (dom%has_dir(ydim)) prop%coords(ydim) = cg_x%y(prop%loc(ydim))
            if (dom%has_dir(zdim)) prop%coords(zdim) = cg_x%z(prop%loc(zdim))
            if (present(dir))      prop%assoc        = cg_x%dl(dir)
!         else
!            write(msg,'(a,a)') "[cg_list_dataop:get_extremum] cg_x not associated for q array name ", qna%lst(ind)%name
!            call die(msg)
         endif
      endif

      if (prop%proc /= 0) then
         if (proc == prop%proc) then ! slave
            call MPI_Send (prop%loc,    ndims, MPI_INTEGER,          FIRST, tag1, comm, mpi_err)
            call MPI_Send (prop%coords, ndims, MPI_DOUBLE_PRECISION, FIRST, tag2, comm, mpi_err)
            if (present(dir)) call MPI_Send (prop%assoc, I_ONE, MPI_DOUBLE_PRECISION, FIRST, tag3, comm, mpi_err)
         endif
         if (master) then
            call MPI_Recv (prop%loc,    ndims, MPI_INTEGER,          prop%proc, tag1, comm, MPI_STATUS_IGNORE, mpi_err)
            call MPI_Recv (prop%coords, ndims, MPI_DOUBLE_PRECISION, prop%proc, tag2, comm, MPI_STATUS_IGNORE, mpi_err)
            if (present(dir)) call MPI_Recv (prop%assoc, I_ONE, MPI_DOUBLE_PRECISION, prop%proc, tag3, comm, MPI_STATUS_IGNORE, mpi_err)
         endif
      endif

   end subroutine get_extremum

!> \brief reset given field to the value (usually 0. or dirty)

   subroutine set_q_value(this, ind, val)

      use cg_list, only: cg_list_element

      implicit none

      class(cg_list_dataop_T), intent(in) :: this    !< object invoking type-bound procedure
      integer(kind=4),         intent(in) :: ind     !< Index in cg%q(:)
      real,                    intent(in) :: val     !< value to put

      type(cg_list_element), pointer      :: cgl

      cgl => this%first
      do while (associated(cgl))
         cgl%cg%q(ind)%arr(:, :, :) = val
         cgl => cgl%nxt
      enddo

   end subroutine set_q_value

!>
!! \brief copy a given field to another
!!
!! \todo OPT: this is quite often used routine. Find when the guardcells can be omitted and save some memory copying.
!<

   subroutine q_copy(this, i_from, i_to)

      use cg_list, only: cg_list_element

      implicit none

      class(cg_list_dataop_T), intent(in) :: this    !< object invoking type-bound procedure
      integer(kind=4),         intent(in) :: i_from  !< Index of source in cg%q(:)
      integer(kind=4),         intent(in) :: i_to    !< Index of destination in cg%q(:)

      type(cg_list_element), pointer      :: cgl

      cgl => this%first
      do while (associated(cgl))
         cgl%cg%q(i_to)%arr(:, :, :) = cgl%cg%q(i_from)%arr(:, :, :)
         cgl => cgl%nxt
      enddo

   end subroutine q_copy

!> \brief copy a given rank-3 field to a component of a rank-4 field (cg%w(w_to)%arr(w_ind,:,:,:) = cg%q(q_from)%arr(:, :, :))

   subroutine qw_copy(this, q_from, w_to, w_ind)

      use cg_list, only: cg_list_element

      implicit none

      class(cg_list_dataop_T), intent(in) :: this    !< object invoking type-bound procedure
      integer(kind=4),         intent(in) :: q_from  !< Index of source in cg%q(:)
      integer(kind=4),         intent(in) :: w_to    !< Index of destination in cg%w(:)
      integer(kind=4),         intent(in) :: w_ind   !< First index of destination in cg%w(w_to)%arr(:,:,:,:)

      type(cg_list_element), pointer      :: cgl

      cgl => this%first
      do while (associated(cgl))
         cgl%cg%w(w_to)%arr(w_ind, :, :, :) = cgl%cg%q(q_from)%arr(:, :, :)
         cgl => cgl%nxt
      enddo

   end subroutine qw_copy

!> \brief copy a component of rank-4 field to a given rank-3 field (cg%q(q_to)%arr(:, :, :) = cg%w(w_from)%arr(w_ind,:,:,:))

   subroutine wq_copy(this, w_from, w_ind, q_to)

      use cg_list, only: cg_list_element

      implicit none

      class(cg_list_dataop_T), intent(in) :: this    !< object invoking type-bound procedure
      integer(kind=4),         intent(in) :: w_from  !< Index of source in cg%w(:)
      integer(kind=4),         intent(in) :: w_ind   !< First index of source in cg%w(w_from)%arr(:,:,:,:)
      integer(kind=4),         intent(in) :: q_to    !< Index of destination in cg%q(:)

      type(cg_list_element), pointer      :: cgl

      cgl => this%first
      do while (associated(cgl))
         cgl%cg%q(q_to)%arr(:, :, :) = cgl%cg%w(w_from)%arr(w_ind, :, :, :)
         cgl => cgl%nxt
      enddo

   end subroutine wq_copy

!> \brief Add a field to another (e.g. apply a correction)

   subroutine q_add(this, i_add, i_to)

      use cg_list, only: cg_list_element

      implicit none

      class(cg_list_dataop_T), intent(in) :: this    !< object invoking type-bound procedure
      integer(kind=4),         intent(in) :: i_add   !< Index of field to be added in cg%q(:)
      integer(kind=4),         intent(in) :: i_to    !< Index of field to be modified in cg%q(:)


      type(cg_list_element), pointer      :: cgl

      cgl => this%first
      do while (associated(cgl))
         cgl%cg%q(i_to)%arr(:, :, :) = cgl%cg%q(i_to)%arr(:, :, :) + cgl%cg%q(i_add)%arr(:, :, :)
         cgl => cgl%nxt
      enddo

   end subroutine q_add

!> \brief Add a value to a field (e.g. correct for average value)

   subroutine q_add_val(this, i_add, val)

      use cg_list, only: cg_list_element

      implicit none

      class(cg_list_dataop_T), intent(in) :: this    !< object invoking type-bound procedure
      integer(kind=4),         intent(in) :: i_add   !< Index of field to be modified in cg%q(:)
      real,                    intent(in) :: val     !< Value to be added


      type(cg_list_element), pointer      :: cgl

      cgl => this%first
      do while (associated(cgl))
         cgl%cg%q(i_add)%arr(:, :, :) = cgl%cg%q(i_add)%arr(:, :, :) + val
         cgl => cgl%nxt
      enddo

   end subroutine q_add_val

!>
!! \brief Assign linear combination of q fields on the whole list
!!
!! \details On the whole list cg%q(ind) is assigned to sum of iv%val * cg%q(iv%ind)
!!
!! \todo add an option to select region of cg
!<

   subroutine q_lin_comb(this, iv, ind)

      use cg_list,    only: cg_list_element
      use dataio_pub, only: die, warn

      implicit none

      class(cg_list_dataop_T),     intent(in) :: this    !< object invoking type-bound procedure
      type(ind_val), dimension(:), intent(in) :: iv      !< list of (coefficient, index) pairs
      integer(kind=4),             intent(in) :: ind     !< Index in cg%q(:)

      integer                                 :: i
      type(ind_val), dimension(size(iv))      :: iv_safe !< sanitized copy of iv
      logical                                 :: swapped
      type(cg_list_element), pointer          :: cgl

      if (size(iv) <= 0) then
         call warn("[cg_list_dataop::q_lin_comb] Nothing to do")
         return
      endif

      iv_safe(:) = iv(:)
      ! if own field (ind) is involved then move it to the first position to avoid side effects and allow in-place operation
      swapped =.false.
      do i = lbound(iv, dim=1)+1, ubound(iv, dim=1)
         if (ind == iv(i)%ind) then
            if (swapped) call die("[cg_list_dataop::q_lin_comb] Cannot use own field twice due to side effects")
            iv_safe(lbound(iv, dim=1)) = iv(i)
            iv_safe(i) = iv(lbound(iv, dim=1))
            swapped = .true.
         endif
      enddo

      cgl => this%first
      do while (associated(cgl))
         cgl%cg%q(ind)%arr(:, :, :) = iv_safe(lbound(iv_safe, dim=1))%val * cgl%cg%q(iv_safe(lbound(iv_safe, dim=1))%ind)%arr(:, :, :)
         do i = lbound(iv_safe, dim=1)+1, ubound(iv_safe, dim=1)
            cgl%cg%q(ind)%arr(:, :, :) = cgl%cg%q(ind)%arr(:, :, :) + iv_safe(i)%val * cgl%cg%q(iv_safe(i)%ind)%arr(:, :, :)
         enddo
         cgl => cgl%nxt
      enddo

    end subroutine q_lin_comb

!>
!! \brief Compute the average value over a list and subtract it
!!
!! \details Typically it is used on a list of leaves or on a single level
!<

   subroutine subtract_average(this, iv)

      use cg_list,          only: cg_list_element
      use constants,        only: GEO_XYZ, GEO_RPZ, pSUM
      use dataio_pub,       only: die
      use domain,           only: dom
      use grid_cont,        only: grid_container
      use mpisetup,         only: piernik_MPI_Allreduce
#ifdef DEBUG
      use dataio_pub,       only: msg, printinfo
      use mpisetup,         only: master
      use named_array_list, only: qna
#endif /* DEBUG */

      implicit none

      class(cg_list_dataop_T), intent(in) :: this !< list for which we want to subtract its average from
      integer(kind=4),         intent(in) :: iv   !< index of variable in cg%q(:) which we want to have zero average

      real                                :: avg, vol
      integer                             :: i
      type(cg_list_element), pointer      :: cgl
      type(grid_container),  pointer      :: cg

      avg = 0.
      vol = 0.
      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg
         select case (dom%geometry_type)
            case (GEO_XYZ)
               avg = avg + sum(cg%q(iv)%span(cg%ijkse), mask=cg%leafmap) * cg%dvol
               vol = vol + count(cg%leafmap) * cg%dvol
            case (GEO_RPZ)
               do i = cg%is, cg%ie
                  avg = avg + sum(cg%q(iv)%arr(i, cg%js:cg%je, cg%ks:cg%ke), mask=cg%leafmap(i, :, :)) * cg%dvol * cg%x(i)
                  vol = vol + count(cg%leafmap(i, :, :)) * cg%dvol * cg%x(i)
               enddo
            case default
               call die("[cg_list_dataop:subtract_average] Unsupported geometry.")
         end select
         cgl => cgl%nxt
      enddo
      call piernik_MPI_Allreduce(avg, pSUM)
      call piernik_MPI_Allreduce(vol, pSUM) !! \todo calculate this in some init routine
      avg = avg / vol

      call this%q_add_val(iv, -avg)

#ifdef DEBUG
      if (master) then
         write(msg, '(2a,2(a,g15.7))')"[cg_list_dataop:subtract_average] Average of '",qna%lst(iv)%name,"' over volume ",vol," is ",avg
         call printinfo(msg)
      endif
#endif /* DEBUG */

   end subroutine subtract_average

!>
!! \brief Calculate L2 norm
!!
!! \todo modify the code for reusing in subtract_average?
!!
!! \todo move to cg_leaves and use select type instead of the nomask parameter
!<

   real function norm_sq(this, iv, nomask) result (norm)

      use cg_list,    only: cg_list_element
      use constants,  only: GEO_XYZ, GEO_RPZ, pSUM
      use dataio_pub, only: die
      use domain,     only: dom
      use grid_cont,  only: grid_container
      use mpisetup,   only: piernik_MPI_Allreduce

      implicit none

      class(cg_list_dataop_T), intent(in) :: this   !< list for which we want to calculate the L2 norm
      integer(kind=4),         intent(in) :: iv     !< index of variable in cg%q(:) for which we want to find the norm
      logical, optional,       intent(in) :: nomask !< Treat the list as a complete level and do not use leafmask

      integer                             :: i
      type(cg_list_element), pointer      :: cgl
      type(grid_container),  pointer      :: cg
      logical                             :: mask

      mask = .true.
      if (present(nomask)) then
         if (nomask) mask = .false.
      endif

      norm = 0.
      cgl => this%first
      do while (associated(cgl))
         cg => cgl%cg
         select case (dom%geometry_type)
            case (GEO_XYZ)
               if (mask) then
                  norm = norm + sum(cg%q(iv)%span(cg%ijkse)**2, mask=cg%leafmap) * cg%dvol
               else
                  norm = norm + sum(cg%q(iv)%span(cg%ijkse)**2) * cg%dvol
               endif
            case (GEO_RPZ)
               do i = cg%is, cg%ie
                  if (mask) then
                     norm = norm + sum(cg%q(iv)%arr(i, cg%js:cg%je, cg%ks:cg%ke)**2, mask=cg%leafmap(i, :, :)) * cg%dvol * cg%x(i)
                  else
                     norm = norm + sum(cg%q(iv)%arr(i, cg%js:cg%je, cg%ks:cg%ke)**2) * cg%dvol * cg%x(i)
                  endif
               enddo
            case default
               call die("[cg_list_dataop:norm_sq] Unsupported geometry.")
         end select
         cgl => cgl%nxt
      enddo
      call piernik_MPI_Allreduce(norm, pSUM)
      norm = sqrt(norm)

   end function norm_sq

!< \brief Multiply two vectors and return the scalar result

   real function scalar_product(this, var1, var2)

      use cg_list,    only: cg_list_element
      use constants,  only: pSUM
      use mpisetup,   only: piernik_MPI_Allreduce

      implicit none

      class(cg_list_dataop_T), intent(in) :: this   !< list for which we want to calculate the scalar product
      integer(kind=4),         intent(in) :: var1   !< first vector
      integer(kind=4),         intent(in) :: var2   !< second vector

      type(cg_list_element), pointer      :: cgl

      scalar_product = 0.

      cgl => this%first
      do while (associated(cgl))
         scalar_product = scalar_product + sum(cgl%cg%q(var1)%span(cgl%cg%ijkse)*cgl%cg%q(var2)%span(cgl%cg%ijkse), mask=cgl%cg%leafmap)
         ! do we need any weighting by volume of elements?
         cgl => cgl%nxt
      enddo
      call piernik_MPI_Allreduce(scalar_product, pSUM)

   end function scalar_product

!< \brief multiply data by the r (first) coordinate (useful to convert phi-momentum into angular momentum)

   subroutine mul_by_r(this, ind)

      use cg_list,   only: cg_list_element

      implicit none

      class(cg_list_dataop_T), intent(in) :: this  !< object invoking type-bound procedure
      integer(kind=4),         intent(in) :: ind   !< index of variable in cg%q(:) which we want to pollute

      type(cg_list_element), pointer :: cgl
      integer :: i

      cgl => this%first
      do while (associated(cgl))
         do i = lbound(cgl%cg%q(ind)%arr, dim=1), ubound(cgl%cg%q(ind)%arr, dim=1)
            cgl%cg%q(ind)%arr(i, :, :) = cgl%cg%q(ind)%arr(i, :, :) * cgl%cg%x(i)
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine mul_by_r

!> \brief divide data by the r (first) coordinate (useful to convert angular momentum into phi-momentum)

   subroutine div_by_r(this, ind)

      use cg_list,    only: cg_list_element
      use constants,  only: zero
      use dataio_pub, only: die
      use func,       only: operator(.equals.)

      implicit none

      class(cg_list_dataop_T), intent(in) :: this  !< object invoking type-bound procedure
      integer(kind=4),         intent(in) :: ind   !< index of variable in cg%q(:) which we want to pollute

      type(cg_list_element), pointer :: cgl
      integer :: i

      cgl => this%first
      do while (associated(cgl))
         do i = lbound(cgl%cg%q(ind)%arr, dim=1), ubound(cgl%cg%q(ind)%arr, dim=1)
            if (cgl%cg%x(i).equals.zero) call die("[cg_list_dataop:div_by_r] cannot divide by radius == 0")
            cgl%cg%q(ind)%arr(i, :, :) = cgl%cg%q(ind)%arr(i, :, :) / cgl%cg%x(i)
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine div_by_r

!> \brief Clear boundary values

   subroutine zero_boundaries(this)

      implicit none

      class(cg_list_dataop_T), intent(inout) :: this  !< list for which clear the boundary values (typically a single level)

      call this%dirty_mg_boundaries(0.)

   end subroutine zero_boundaries

!> \brief Mark boundary values with given value

   subroutine dirty_mg_boundaries(this, value)

      use cg_list, only: cg_list_element

      implicit none

      class(cg_list_dataop_T), intent(inout) :: this   !< list for which clear the boundary values (typically a single level)
      real,                    intent(in)    :: value  !< value to pollute

      type(cg_list_element), pointer         :: cgl

      cgl => this%first
      do while (associated(cgl))
         cgl%cg%mg%bnd_x(:,:,:) = value
         cgl%cg%mg%bnd_y(:,:,:) = value
         cgl%cg%mg%bnd_z(:,:,:) = value
         cgl => cgl%nxt
      enddo

   end subroutine dirty_mg_boundaries

!>
!! \brief This routine pollutes selected array with an insane value dirtyH.
!!
!! \details If anything in the multigrid works by accident, through compiler-dependent initialization or unintentional relying on outdated values,
!! the insane value should pollute the solution in an easily visible way.
!<

   subroutine set_dirty(this, iv)

      use constants, only: dirtyH
      use global,    only: dirty_debug

      implicit none

      class(cg_list_dataop_T), intent(inout) :: this !< list for which we want to apply pollution
      integer(kind=4),         intent(in)    :: iv   !< index of variable in cg%q(:) which we want to pollute

      if (.not. dirty_debug) return

      call this%set_q_value(iv, dirtyH)

   end subroutine set_dirty

!> \brief Check all registered arrays for detectable traces of set_dirty calls.

   subroutine check_all_dirty(this, label, expand, warn_only)

      use named_array_list, only: qna, wna

      implicit none

      class(cg_list_dataop_T),   intent(inout) :: this      !< level which we are checking
      character(len=*),          intent(in)    :: label     !< label to indicate the origin of call
      integer(kind=4), optional, intent(in)    :: expand    !< also check guardcells
      logical, optional,         intent(in)    :: warn_only !< do not die when dirty value has been spotted

      integer(kind=4) :: iv, iw

      do iv = lbound(qna%lst, dim=1, kind=4), ubound(qna%lst, dim=1, kind=4)
         call this%check_dirty(iv, label, expand=expand, warn_only=warn_only)
      enddo

      do iw = lbound(wna%lst, dim=1, kind=4), ubound(wna%lst, dim=1, kind=4)
         do iv = 1, wna%lst(iw)%dim4
            call this%check_dirty(iw, label, expand=expand, subfield=iv, warn_only=warn_only)
         enddo
      enddo

   end subroutine check_all_dirty

!> \brief This routine checks for detectable traces of set_dirty calls.

   subroutine check_dirty(this, iv, label, expand, subfield, warn_only)

      use cg_list,          only: cg_list_element
      use constants,        only: dirtyL, pSUM
      use dataio_pub,       only: warn, msg, die
      use domain,           only: dom
      use global,           only: dirty_debug, show_n_dirtys, no_dirty_checks
      use mpisetup,         only: proc, master, piernik_MPI_Allreduce
      use named_array_list, only: qna, wna

      implicit none

      class(cg_list_dataop_T),   intent(inout) :: this      !< level which we are checking
      integer(kind=4),           intent(in)    :: iv        !< index of variable in cg%q(:) which we want to check
      character(len=*),          intent(in)    :: label     !< label to indicate the origin of call
      integer(kind=4), optional, intent(in)    :: expand    !< also check guardcells
      integer(kind=4), optional, intent(in)    :: subfield  !< when present use it to check cg%w array
      logical, optional,         intent(in)    :: warn_only !< do not die when dirty value has been spotted

      integer                                  :: i, j, k, ng, cnt
      type(cg_list_element), pointer           :: cgl

      if (.not. dirty_debug .or. no_dirty_checks) return

      if (present(subfield)) then
         if (iv < lbound(wna%lst, dim=1) .or. iv > ubound(wna%lst, dim=1)) call die("[cg_list_dataop:check_dirty] Invalid w-variable index.")
         if (subfield < 1 .or. subfield > wna%lst(iv)%dim4) call die("[cg_list_dataop:check_dirty] Invalid w-variable component.")
      else
         if (iv < lbound(qna%lst, dim=1) .or. iv > ubound(qna%lst, dim=1)) call die("[cg_list_dataop:check_dirty] Invalid q-variable index.")
      endif

      ng = 0
      if (present(expand)) ng = min(dom%nb, expand)

      cnt = 0
      cgl => this%first
      do while (associated(cgl))
         do k = cgl%cg%ks-ng*dom%D_z, cgl%cg%ke+ng*dom%D_z
            do j = cgl%cg%js-ng*dom%D_y, cgl%cg%je+ng*dom%D_y
               do i = cgl%cg%is-ng*dom%D_x, cgl%cg%ie+ng*dom%D_x
                  ! if (count([i<cgl%cg%is .or. i>cgl%cg%ie, j<cgl%cg%js .or. j>cgl%cg%je, k<cgl%cg%ks .or. k>cgl%cg%ke]) <=1) then ! excludes corners
                  if (present(subfield)) then
                     if (associated(cgl%cg%w(iv)%arr)) then
                        if (abs(cgl%cg%w(iv)%arr(subfield, i, j, k)) > dirtyL) then
                           if (cnt <= show_n_dirtys) then
                              if (cnt < show_n_dirtys) then
                                 write(msg, '(3a,i4,a,i3,a,i5,3a,4i6,a,g20.12)') &
                                      &                   "[cg_list_dataop:check_dirty] ", trim(label), "@", proc, " lvl^", cgl%cg%level_id, " cg#", cgl%cg%grid_id, &
                                      &                   " '", trim(wna%lst(iv)%name), "'(", subfield, i, j, k, ") = ", cgl%cg%w(iv)%arr(subfield, i, j, k)
                                 call warn(msg)
                              endif
                           endif
                           cnt = cnt + 1
                        endif
                     endif
                  else
                     if (associated(cgl%cg%q(iv)%arr)) then
                        if (abs(cgl%cg%q(iv)%arr(i, j, k)) > dirtyL) then
                           if (cnt <= show_n_dirtys) then
                              if (cnt < show_n_dirtys) then
                                 write(msg, '(3a,i4,a,i3,a,i5,3a,3i6,a,g20.12)') &
                                      &                   "[cg_list_dataop:check_dirty] ", trim(label), "@", proc, " lvl^", cgl%cg%level_id, " cg#", cgl%cg%grid_id, &
                                      &                   " '", trim(qna%lst(iv)%name), "'(", i, j, k, ") = ", cgl%cg%q(iv)%arr(i, j, k)
                                 call warn(msg)
                              endif
                           endif
                           cnt = cnt + 1
                        endif
                     endif
                  endif
                  ! endif
               enddo
            enddo
         enddo
         cgl => cgl%nxt
      enddo

      call piernik_MPI_Allreduce(cnt, pSUM)
      if (cnt /= 0 .and. master) then
         write(msg,'(3a,i8,a)')"[cg_list_dataop:check_dirty] @'",trim(label),"' Found ", cnt, " dirty value(s) in '"
         if (present(subfield)) then
            write(msg(len_trim(msg)+1:), '(2a,i3,a)') trim(wna%lst(iv)%name),"(",subfield,")'"
         else
            write(msg(len_trim(msg)+1:), '(2a)') trim(qna%lst(iv)%name),"'"
         endif
         if (present(warn_only)) then
            if (warn_only) then
               call warn(msg)
            else
               call die(msg)
            endif
         else
            call die(msg)
         endif
      endif

   end subroutine check_dirty

!> \brief Check values of all named arrays for big_float

   subroutine check_for_dirt(this)

      use cg_list,          only: cg_list_element
      use constants,        only: big_float
      use dataio_pub,       only: warn, msg
      use named_array_list, only: qna, wna

      implicit none

      class(cg_list_dataop_T), intent(in) :: this          !< object invoking type-bound procedure

      integer                             :: i
      type(cg_list_element), pointer      :: cgl

      cgl => this%first
      do while (associated(cgl))
         do i = lbound(qna%lst(:), dim=1), ubound(qna%lst(:), dim=1)
            if (cgl%cg%q(i)%check()) then
               write(msg,'(3a,I12,a)') "[cg_list_dataop:check_for_dirt] Array ", trim(qna%lst(i)%name), " has ", &
                  & count(cgl%cg%q(i)%arr >= big_float), " wrong values."
               call warn(msg)
            endif
         enddo
         do i = lbound(wna%lst(:), dim=1), ubound(wna%lst(:), dim=1)
            if (cgl%cg%w(i)%check()) then
               write(msg,'(3a,I12,a)') "[cg_list_dataop:check_for_dirt] Array ", trim(wna%lst(i)%name), " has ", &
                  & count(cgl%cg%w(i)%arr >= big_float), " wrong values."
               call warn(msg)
            endif
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine check_for_dirt

!>
!! \brief Cut out a cross in the middle of the grid block and move the data to shov what is stored in corners and partially in edges and faces as well.
!! Put the result in cg%wa. Useful for debugging boundaries
!<

   subroutine corners2wa(this, ind)

      use cg_list,   only: cg_list_element
      use constants, only: xdim, ydim, zdim, LO, HI
      use domain,    only: dom

      implicit none

      class(cg_list_dataop_T), intent(in) :: this  !< object invoking type-bound procedure
      integer(kind=4),         intent(in) :: ind   !< index of variable in cg%q(:) which we want to pollute

      type(cg_list_element), pointer :: cgl

      cgl => this%first
      do while (associated(cgl))

         associate ( cg => cgl%cg )
         cg%prolong_xyz = cg%q(ind)%arr

         if (dom%has_dir(xdim)) then
            cg%wa(cg%ijkse(xdim, LO):cg%ijkse(xdim, LO)+cg%n_b(xdim)/2-1, :, :) = cg%prolong_xyz(cg%lhn(xdim, LO):cg%lhn(xdim, LO)+cg%n_b(xdim)/2-1, :, :)
            cg%wa(cg%ijkse(xdim, HI)-cg%n_b(xdim)/2+1:cg%ijkse(xdim, HI), :, :) = cg%prolong_xyz(cg%lhn(xdim, HI)-cg%n_b(xdim)/2+1:cg%lhn(xdim, HI), :, :)
            cg%prolong_xyz = cg%wa
         endif

         if (dom%has_dir(ydim)) then
            cg%wa(:, cg%ijkse(ydim, LO):cg%ijkse(ydim, LO)+cg%n_b(ydim)/2-1, :) = cg%prolong_xyz(:, cg%lhn(ydim, LO):cg%lhn(ydim, LO)+cg%n_b(ydim)/2-1, :)
            cg%wa(:, cg%ijkse(ydim, HI)-cg%n_b(ydim)/2+1:cg%ijkse(ydim, HI), :) = cg%prolong_xyz(:, cg%lhn(ydim, HI)-cg%n_b(ydim)/2+1:cg%lhn(ydim, HI), :)
            cg%prolong_xyz = cg%wa
         endif

         if (dom%has_dir(zdim)) then
            cg%wa(:, :, cg%ijkse(zdim, LO):cg%ijkse(zdim, LO)+cg%n_b(zdim)/2-1) = cg%prolong_xyz(:, :, cg%lhn(zdim, LO):cg%lhn(zdim, LO)+cg%n_b(zdim)/2-1)
            cg%wa(:, :, cg%ijkse(zdim, HI)-cg%n_b(zdim)/2+1:cg%ijkse(zdim, HI)) = cg%prolong_xyz(:, :, cg%lhn(zdim, HI)-cg%n_b(zdim)/2+1:cg%lhn(zdim, HI))
         endif
         end associate

         cgl => cgl%nxt
      enddo

   end subroutine corners2wa

end module cg_list_dataop
