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

!> \brief This module contains list of all grid containers and related methods

module cg_list_global

   use cg_list_bnd, only: cg_list_bnd_T

   implicit none

   private
   public :: all_cg

   !>
   !! \brief A list of grid containers that are supposed to have the same variables registered
   !!
   !! \details The main purpose of this type is to provide a type for the set of all grid containers with methods and properties
   !! that should not be available for any arbitrarily composed subset of grid containers. Typically there will be only one variable
   !! of this type available in the code: all_cg.
   !!
   !! It should be possible to use this type for more fancy, multi-domain  grid configurations such as:
   !! - Yin-Yang grid (covering sphere with two domains shaped as parts of the sphere to avoid polar singularities),
   !! - Cylindrical grid with Cartesian core covering singularity at the axis.
   !! - Simulations with mixed dimensionality (e.g. 2d grid for dust particles and 3d grid for gas) should probably also use separate cg_list
   !! for their data (and additional routine for coupling the two grid sets).
   !<
   type, extends(cg_list_bnd_T) :: cg_list_global_T

      integer(kind=4) :: ord_prolong_nb                !< Maximum number of boundary cells required for prolongation

    contains
      procedure         :: reg_var         !< Add a variable (cg%q or cg%w) to all grid containers
      procedure         :: register_fluids !< Register all crucial fields, which we cannot live without
      procedure         :: check_na        !< Check if all named arrays are consistently registered
   end type cg_list_global_T

   type(cg_list_global_T) :: all_cg   !< all grid containers; \todo restore protected

contains

!>
!! \brief Use this routine to add a variable (cg%q or cg%w) to all grid containers.
!!
!! \details Register a rank-3 array of given name in each grid container (in cg%q) and decide what to do with it on restart.
!! When dim4 is present then create a rank-4 array instead.(in cg%w)
!<

   subroutine reg_var(this, name, vital, restart_mode, ord_prolong, dim4, position, multigrid)

      use constants,   only: INVALID, VAR_CENTER, AT_NO_B, AT_IGNORE, I_ZERO, I_ONE, I_TWO, O_INJ, O_LIN, O_I2, O_D2, O_I3, O_I4, O_D3, O_D4
      use dataio_pub,  only: die, warn, msg
      use domain,      only: dom
      use cg_list,     only: cg_list_element
      use named_array, only: qna, wna, na_var

      implicit none

      class(cg_list_global_T),                     intent(inout) :: this          !< object invoking type-bound procedure
      character(len=*),                        intent(in)    :: name          !< Name of the variable to be registered
      logical,                       optional, intent(in)    :: vital         !< .false. for arrays that don't need to be prolonged or restricted automatically
      integer(kind=4),               optional, intent(in)    :: restart_mode  !< Write to the restart if not AT_IGNORE. Several write modes can be supported.
      integer(kind=4),               optional, intent(in)    :: ord_prolong   !< Prolongation order for the variable
      integer(kind=4),               optional, intent(in)    :: dim4          !< If present then register the variable in the cg%w array.
      integer(kind=4), dimension(:), optional, intent(in)    :: position      !< If present then use this value instead of VAR_CENTER
      logical,                       optional, intent(in)    :: multigrid     !< If present and .true. then allocate cg%q(:)%arr and cg%w(:)%arr also below base level

      type(cg_list_element), pointer :: cgl
      logical :: mg, vit
      integer :: nvar
      integer(kind=4) :: op, d4, rm
      integer(kind=4), allocatable, dimension(:) :: pos

      vit = .false.
      if (present(vital)) vit = vital

      rm = AT_IGNORE
      if (present(restart_mode)) rm = restart_mode

      op = O_INJ
      if (present(ord_prolong)) op = ord_prolong

      mg = .false.
      if (present(multigrid)) mg = multigrid

      if (present(dim4)) then
         if (mg) call die("[cg_list_global:reg_var] there are no rank-4 multigrid arrays yet")
         d4 = dim4
         nvar = dim4
      else
         d4 = int(INVALID, kind=4)
         nvar = 1
      endif

      if (allocated(pos)) call die("[cg_list_global:reg_var] pos(:) already allocated")
      allocate(pos(nvar))
      pos(:) = VAR_CENTER
      if (present(position)) then
         if (any(size(position) == [1, nvar])) then
            pos = position  !> \deprecated BEWARE: lhs reallocation
         else
            write(msg,'(2(a,i3))')"[cg_list_global:reg_var] position should be an array of 1 or ",nvar," values. Got ",size(position)
            call die(msg)
         endif
      endif
      if (any(pos(:) /= VAR_CENTER) .and. rm == AT_NO_B) then
         write(msg,'(3a)')"[cg_list_global:reg_var] no boundaries for restart with non cel-centered variable '",name,"' may result in loss of information in the restart files."
         call warn(msg)
      endif

      if (present(dim4)) then
         call wna%add2lst(na_var(name, vit, rm, op, pos, d4, mg))
      else
         call qna%add2lst(na_var(name, vit, rm, op, pos, d4, mg))
      endif

      select case (op)
         case (O_INJ)
            this%ord_prolong_nb = max(this%ord_prolong_nb, I_ZERO)
         case (O_LIN, O_I2, O_D2)
            this%ord_prolong_nb = max(this%ord_prolong_nb, I_ONE)
         case (O_I3, O_I4, O_D3, O_D4)
            this%ord_prolong_nb = max(this%ord_prolong_nb, I_TWO)
         case default
            call die("[cg_list_global:reg_var] Unknown prolongation order")
      end select
      if (this%ord_prolong_nb > dom%nb) call die("[cg_list_global:reg_var] Insufficient number of guardcells for requested prolongation stencil")

      cgl => this%first
      do while (associated(cgl))
         if (present(dim4)) then
            call cgl%cg%add_na_4d(dim4)
         else
            call cgl%cg%add_na(mg)
         endif
         cgl => cgl%nxt
      enddo

      deallocate(pos)

   end subroutine reg_var

!> \brief Register all crucial fields, which we cannot live without

   subroutine register_fluids(this, nfluids)

      use constants,  only: wa_n, fluid_n, uh_n, mag_n, u0_n, b0_n, ndims, AT_NO_B, AT_OUT_B, VAR_XFACE, VAR_YFACE, VAR_ZFACE, PIERNIK_INIT_FLUIDS
      use dataio_pub, only: die, code_progress
      use global,     only: repeat_step
#ifdef ISO
      use constants,  only: cs_i2_n
#endif /* ISO */

      implicit none

      class(cg_list_global_T), intent(inout) :: this          !< object invoking type-bound procedure
      integer(kind=4),     intent(in)    :: nfluids       !< number of components in the main array of fluids (should be flind%all)

      integer(kind=4), parameter, dimension(ndims) :: xyz_face = [ VAR_XFACE, VAR_YFACE, VAR_ZFACE ]

      if (code_progress < PIERNIK_INIT_FLUIDS) call die("[cg_list_global:register_fluids] Fluids are not yet initialized")

      call this%reg_var(wa_n,                                                           multigrid=.true.)  !! Auxiliary array. Multigrid required only for CR diffusion
      call this%reg_var(fluid_n, vital = .true., restart_mode = AT_NO_B,  dim4 = nfluids)                  !! Main array of all fluids' components, "u"
      call this%reg_var(uh_n,                                             dim4 = nfluids)                  !! Main array of all fluids' components (for t += dt/2)
      call this%reg_var(mag_n,   vital = .true., restart_mode = AT_OUT_B, dim4 = ndims, position=xyz_face) !! Main array of magnetic field's components, "b"
      if (repeat_step) then
         call this%reg_var(u0_n,                                          dim4 = nfluids)                  !! Copy of main array of all fluids' components
         call this%reg_var(b0_n,                                          dim4 = ndims, position=xyz_face) !! Copy of main array of magnetic field's components
      endif
#ifdef ISO
      call all_cg%reg_var(cs_i2_n, vital = .true., restart_mode = AT_NO_B)
#endif /* ISO */

   end subroutine register_fluids

!> \brief Check if all named arrays are consistently registered

   subroutine check_na(this)

      use constants,   only: INVALID, base_level_id
      use dataio_pub,  only: msg, die
      use cg_list,     only: cg_list_element
      use named_array, only: qna, wna

      implicit none

      class(cg_list_global_T), intent(in) :: this          !< object invoking type-bound procedure

      integer :: i
      type(cg_list_element), pointer :: cgl
      logical :: bad

      cgl => this%first
      do while (associated(cgl))
         if (allocated(qna%lst) .neqv. allocated(cgl%cg%q)) then
            write(msg,'(2(a,l2))')"[cg_list_global:check_na] allocated(qna%lst) .neqv. allocated(cgl%cg%q):",allocated(qna%lst)," .neqv. ",allocated(cgl%cg%q)
            call die(msg)
         else if (allocated(qna%lst)) then
            if (size(qna%lst(:)) /= size(cgl%cg%q)) then
               write(msg,'(2(a,i5))')"[cg_list_global:check_na] size(qna) /= size(cgl%cg%q)",size(qna%lst(:))," /= ",size(cgl%cg%q)
               call die(msg)
            else
               do i = lbound(qna%lst(:), dim=1), ubound(qna%lst(:), dim=1)
                  if (qna%lst(i)%dim4 /= INVALID) then
                     write(msg,'(3a,i10)')"[cg_list_global:check_na] qna%lst(",i,"_, named '",qna%lst(i)%name,"' has dim4 set to ",qna%lst(i)%dim4
                     call die(msg)
                  endif
                  if (associated(cgl%cg%q(i)%arr) .and. cgl%cg%level_id < base_level_id .and. .not. qna%lst(i)%multigrid) then
                     write(msg,'(a,i3,3a)')"[cg_list_global:check_na] non-multigrid cgl%cg%q(",i,"), named '",qna%lst(i)%name,"' allocated on coarse level"
                     call die(msg)
                  endif
               enddo
            endif
         endif
         if (allocated(wna%lst) .neqv. allocated(cgl%cg%w)) then
            write(msg,'(2(a,l2))')"[cg_list_global:check_na] allocated(wna%lst) .neqv. allocated(cgl%cg%w)",allocated(wna%lst)," .neqv. ",allocated(cgl%cg%w)
            call die(msg)
         else if (allocated(wna%lst)) then
            if (size(wna%lst(:)) /= size(cgl%cg%w)) then
               write(msg,'(2(a,i5))')"[cg_list_global:check_na] size(wna) /= size(cgl%cg%w)",size(wna%lst(:))," /= ",size(cgl%cg%w)
               call die(msg)
            else
               do i = lbound(wna%lst(:), dim=1), ubound(wna%lst(:), dim=1)
                  bad = .false.
                  if (associated(cgl%cg%w(i)%arr)) bad = wna%lst(i)%dim4 /= size(cgl%cg%w(i)%arr, dim=1) .and. cgl%cg%level_id >= base_level_id
                  if (wna%lst(i)%dim4 <= 0 .or. bad) then
                     write(msg,'(a,i3,2a,2(a,i7))')"[cg_list_global:check_na] wna%lst(",i,"_ named '",wna%lst(i)%name,"' has inconsistent dim4: ",&
                          &         wna%lst(i)%dim4," /= ",size(cgl%cg%w(i)%arr, dim=1)
                     call die(msg)
                  endif
                  if (associated(cgl%cg%w(i)%arr) .and. cgl%cg%level_id < base_level_id) then
                     write(msg,'(a,i3,3a)')"[cg_list_global:check_na] cgl%cg%w(",i,"), named '",wna%lst(i)%name,"' allocated on coarse level"
                     call die(msg)
                  endif
               enddo
            endif
         endif
         cgl => cgl%nxt
      enddo

   end subroutine check_na

end module cg_list_global
