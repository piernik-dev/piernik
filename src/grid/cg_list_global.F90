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

#if defined(__INTEL_COMPILER)
   !! \deprecated remove this clause as soon as Intel Compiler gets required
   !! features and/or bug fixes
   use cg_list,     only: cg_list_t   ! QA_WARN ICE is the alternative :)
#endif /*__INTEL_COMPILER) */
   use cg_list_bnd, only: cg_list_bnd_t
   use constants,   only: dsetnamelen

   implicit none

   private
   public :: all_cg, all_cg_n

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
   type, extends(cg_list_bnd_t) :: cg_list_global_t
      integer(kind=4) :: ord_prolong_nb  !< Maximum number of boundary cells required for prolongation
   contains
      procedure :: init             !< Initialize
      procedure :: reg_var          !< Add a variable (cg%q or cg%w) to all grid containers
      procedure :: register_fluids  !< Register all crucial fields, which we cannot live without
      procedure :: check_na         !< Check if all named arrays are consistently registered
      procedure :: delete_all       !< Delete the grid container from all lists
      procedure :: mark_orphans     !< Find grid pieces that do not belong to any list except for all_cg
   end type cg_list_global_t

   type(cg_list_global_t)                :: all_cg   !< all grid containers; \todo restore protected
   character(len=dsetnamelen), parameter :: all_cg_n = "all_cg" !< name of the all_cg list

contains

!> \brief Initialize

   subroutine init(this)

      use constants,        only: I_ZERO
      use list_of_cg_lists, only: all_lists

      implicit none

      class(cg_list_global_t), intent(inout) :: this           !< object invoking type-bound procedure

      call all_lists%register(this, all_cg_n)
      this%ord_prolong_nb = I_ZERO
      call this%reset_costs

   end subroutine init

!> \brief destroy the global list, all grid containers and all lists

   subroutine delete_all(this)

      use dataio_pub,       only: die
      use grid_cont,        only: grid_container
      use list_of_cg_lists, only: all_lists

      implicit none

      class(cg_list_global_t), intent(inout) :: this           !< object invoking type-bound procedure

      type(grid_container),  pointer         :: cg

      !> \todo implement what is said in the description

      do while (associated(this%first))
         if (associated(this%last%cg)) then
            cg => this%last%cg
            call all_lists%forget(cg)
         else
            call die("[cg_list_global:delete_from_all] Attempted to remove an empty element")
         endif
      enddo

   end subroutine delete_all

!>
!! \brief Use this routine to add a variable (cg%q or cg%w) to all grid containers.
!!
!! \details Register a rank-3 array of given name in each grid container (in cg%q) and decide what to do with it on restart.
!! When dim4 is present then create a rank-4 array instead.(in cg%w)
!<

   subroutine reg_var(this, name, vital, restart_mode, ord_prolong, dim4, position, multigrid)

      use cg_list,          only: cg_list_element
      use constants,        only: INVALID, VAR_CENTER, AT_NO_B, AT_IGNORE, I_ZERO, I_ONE, I_TWO, I_THREE, O_INJ, O_LIN, O_I2, O_D2, O_I3, O_I4, O_I5, O_I6, O_D3, O_D4, O_D5, O_D6
      use dataio_pub,       only: die, warn, msg
      use domain,           only: dom
      use memory_usage,     only: check_mem_usage
      use named_array_list, only: qna, wna, na_var, na_var_4d

      implicit none

      class(cg_list_global_t),                 intent(inout) :: this          !< object invoking type-bound procedure
      character(len=*),                        intent(in)    :: name          !< Name of the variable to be registered
      logical,                       optional, intent(in)    :: vital         !< .false. for arrays that don't need to be prolonged or restricted automatically
      integer(kind=4),               optional, intent(in)    :: restart_mode  !< Write to the restart if >= AT_IGNORE. Several write modes can be supported.
      integer(kind=4),               optional, intent(in)    :: ord_prolong   !< Prolongation order for the variable
      integer(kind=4),               optional, intent(in)    :: dim4          !< If present then register the variable in the cg%w array.
      integer(kind=4), dimension(:), optional, intent(in)    :: position      !< If present then use this value instead of VAR_CENTER
      logical,                       optional, intent(in)    :: multigrid     !< If present and .true. then allocate cg%q(:)%arr and cg%w(:)%arr also below base level

      type(cg_list_element), pointer             :: cgl
      logical                                    :: mg, vit
      integer                                    :: nvar
      integer(kind=4)                            :: op, d4, rm
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
            pos = position
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
         call wna%add2lst(na_var_4d(name, vit, rm, op, mg, dim4=d4))
      else
         call qna%add2lst(na_var(name, vit, rm, op, mg))
      endif

      select case (op)
         case (O_INJ)
            this%ord_prolong_nb = max(this%ord_prolong_nb, I_ZERO)
         case (O_LIN, O_I2, O_D2)
            this%ord_prolong_nb = max(this%ord_prolong_nb, I_ONE)
         case (O_I3, O_I4, O_D3, O_D4)
            this%ord_prolong_nb = max(this%ord_prolong_nb, I_TWO)
         case (O_I5, O_I6, O_D5, O_D6)
            this%ord_prolong_nb = max(this%ord_prolong_nb, I_THREE)
         case default
            call die("[cg_list_global:reg_var] Unknown prolongation order")
      end select
      if (I_TWO*this%ord_prolong_nb > dom%nb) call die("[cg_list_global:reg_var] Insufficient number of guardcells for requested prolongation stencil. Expected crash in cg_level_connected::prolong_bnd_from_coarser")
      ! I-TWO because our refinement factor is 2. and we want to fill all layers of fine guardcells
      ! Technically it is possible to maintain high order prolongation and thin layer of guardcells,
      ! but fine boundaries that coincide with coarse boundaries cannot be fully reconstructed
      ! unless we communicate the missing part from another block (do multi-parent prolongation).

      cgl => this%first
      do while (associated(cgl))
         if (present(dim4)) then
            call cgl%cg%add_na_4d(d4)  ! Strange: passing dim4 here resulted in an access to already freed memory. Possibly a gfortran bug.
         else
            call cgl%cg%add_na(mg)
         endif
         cgl => cgl%nxt
      enddo
      call check_mem_usage

      deallocate(pos)

   end subroutine reg_var

!> \brief Register all crucial fields, which we cannot live without

   subroutine register_fluids(this)

      use constants,  only: wa_n, fluid_n, uh_n, AT_NO_B, PIERNIK_INIT_FLUIDS
      use dataio_pub, only: die, code_progress
      use fluidindex, only: flind
      use global,     only: ord_fluid_prolong
#ifdef ISO
      use constants,  only: cs_i2_n
#endif /* ISO */
#ifdef MAGNETIC
      use constants,  only: mag_n, magh_n, ndims, AT_OUT_B, VAR_XFACE, VAR_YFACE, VAR_ZFACE, VAR_CENTER, psi_n, psih_n
      use global,     only: cc_mag, ord_mag_prolong
#endif /* MAGNETIC */

      implicit none

      class(cg_list_global_t), intent(inout)          :: this          !< object invoking type-bound procedure

#ifdef MAGNETIC
      integer(kind=4), dimension(ndims), parameter :: xyz_face = [ VAR_XFACE, VAR_YFACE, VAR_ZFACE ]
      integer(kind=4), dimension(ndims), parameter :: xyz_center = [ VAR_CENTER, VAR_CENTER, VAR_CENTER ]
      integer(kind=4), dimension(ndims) :: pia

      pia = merge(xyz_center, xyz_face, cc_mag)
#endif /* MAGNETIC */

      if (code_progress < PIERNIK_INIT_FLUIDS) call die("[cg_list_global:register_fluids] Fluids are not yet initialized")

      call this%reg_var(wa_n, multigrid=.true.)  !! Auxiliary array. Multigrid required only for CR diffusion
      call this%reg_var(fluid_n, vital = .true., restart_mode = AT_NO_B,  dim4 = flind%all, ord_prolong = ord_fluid_prolong) !! Main array of all fluids' components, "u"
      call this%reg_var(uh_n,                                             dim4 = flind%all, ord_prolong = ord_fluid_prolong) !! Main array of all fluids' components (for t += dt/2)
      call set_fluid_names
#ifdef COSM_RAYS
      call set_cr_names
#endif /* COSM_RAYS */
#ifdef CRESP
      call set_cresp_names
#endif /* CRESP */

#ifdef MAGNETIC
      call this%reg_var(mag_n,  vital = .true.,  dim4 = ndims, ord_prolong = ord_mag_prolong, restart_mode = AT_OUT_B, position=pia)  !! Main array of magnetic field's components, "b"
      call this%reg_var(magh_n, vital = .false., dim4 = ndims) !! Array for copy of magnetic field's components, "b" used in half-timestep in RK2
      call set_magnetic_names

      if (cc_mag) then
         call this%reg_var(psi_n,  vital = .true., ord_prolong = ord_mag_prolong, restart_mode = AT_OUT_B)  !! an array for div B cleaning
         call this%reg_var(psih_n, vital = .false.)  !! its copy for use in RK2
      endif
#endif /* MAGNETIC */

#ifdef ISO
      call all_cg%reg_var(cs_i2_n, vital = .true., restart_mode = AT_NO_B)
#endif /* ISO */

   contains

      subroutine set_fluid_names

         use fluidindex,       only: flind
         use fluids_pub,       only: has_dst, has_ion, has_neu
         use named_array_list, only: wna, na_var_4d

         implicit none

         select type (lst => wna%lst)
            type is (na_var_4d)
               if (has_ion) then
                  call lst(wna%fi)%set_compname(flind%ion%idn, "deni")
                  call lst(wna%fi)%set_compname(flind%ion%imx, "momxi")
                  call lst(wna%fi)%set_compname(flind%ion%imy, "momyi")
                  call lst(wna%fi)%set_compname(flind%ion%imz, "momzi")
                  if (flind%ion%has_energy) call lst(wna%fi)%set_compname(flind%ion%ien, "enei")
               endif

               if (has_neu) then
                  call lst(wna%fi)%set_compname(flind%neu%idn, "denn")
                  call lst(wna%fi)%set_compname(flind%neu%imx, "momxn")
                  call lst(wna%fi)%set_compname(flind%neu%imy, "momyn")
                  call lst(wna%fi)%set_compname(flind%neu%imz, "momzn")
                  if (flind%neu%has_energy) call lst(wna%fi)%set_compname(flind%neu%ien, "enen")
               endif

               if (has_dst) then
                  call lst(wna%fi)%set_compname(flind%dst%idn, "dend")
                  call lst(wna%fi)%set_compname(flind%dst%imx, "momxd")
                  call lst(wna%fi)%set_compname(flind%dst%imy, "momyd")
                  call lst(wna%fi)%set_compname(flind%dst%imz, "momzd")
               endif
         end select

      end subroutine set_fluid_names

#ifdef COSM_RAYS
      subroutine set_cr_names

         use constants,        only: dsetnamelen, I_ONE
         use cr_data,          only: cr_names, cr_spectral
         use named_array_list, only: wna, na_var_4d

         implicit none

         integer(kind=4) :: i, k
         character(len=dsetnamelen) :: var

         select type (lst => wna%lst)
            type is (na_var_4d)
               k = flind%crn%beg
               do i = I_ONE, size(cr_names, kind=4)
                  if (.not. cr_spectral(i)) then
                     if (len_trim(cr_names(i)) > 0) then
                        write(var, '(2a)') "cr_", trim(cr_names(i))
                     else
                        write(var, '(a,i2.2)') "cr", i
                     endif
                     call lst(wna%fi)%set_compname(k, var)
                     k = k + I_ONE
                  endif
               enddo
            class default
               call die("[cg_list_global:set_cr_names] Unknown list type")
         end select

      end subroutine set_cr_names

#endif /* COSM_RAYS */

#ifdef CRESP
      subroutine set_cresp_names

         use constants,        only: dsetnamelen
         ! use cr_data,          only: cr_names
         use named_array_list, only: wna, na_var_4d

         implicit none

         integer(kind=4) :: i
         character(len=dsetnamelen) :: var

         ! The "e-" part of the name is used for the CR energy density should be cr_names(1) currently
         ! After merge of Antoine's branch the Isotope names should go there
         select type (lst => wna%lst)
            type is (na_var_4d)
               do i = flind%cre%nbeg, flind%cre%nend
                  write(var, '(a,i2.2)') "cr_e-n", i - flind%cre%nbeg + 1
                  call lst(wna%fi)%set_compname(i, var)
               enddo
               do i = flind%cre%ebeg, flind%cre%eend
                  write(var, '(a,i2.2)') "cr_e-e", i - flind%cre%ebeg + 1
                  call lst(wna%fi)%set_compname(i, var)
               enddo
            class default
               call die("[cg_list_global:set_cresp_names] Unknown list type")
         end select

      end subroutine set_cresp_names
#endif /* CRESP */

#ifdef MAGNETIC
      subroutine set_magnetic_names

         use constants,        only: xdim, ydim, zdim
         use named_array_list, only: wna, na_var_4d

         implicit none

         select type (lst => wna%lst)
            type is (na_var_4d)

               call lst(wna%bi)%set_compname(xdim, "magx")
               call lst(wna%bi)%set_compname(ydim, "magy")
               call lst(wna%bi)%set_compname(zdim, "magz")
         end select

      end subroutine set_magnetic_names
#endif /* MAGNETIC */

   end subroutine register_fluids

!> \brief Check if all named arrays are consistently registered

   subroutine check_na(this)

      use constants,        only: base_level_id
      use dataio_pub,       only: msg, die
      use cg_list,          only: cg_list_element
      use named_array_list, only: qna, wna

      implicit none

      class(cg_list_global_t), intent(in) :: this          !< object invoking type-bound procedure

      integer(kind=4)                     :: i
      type(cg_list_element), pointer      :: cgl
      logical                             :: bad

      cgl => this%first
      do while (associated(cgl))
         if (associated(cgl%cg)) then
            if (allocated(qna%lst) .neqv. allocated(cgl%cg%q)) then
               write(msg,'(2(a,l2))')"[cg_list_global:check_na] allocated(qna%lst) .neqv. allocated(cgl%cg%q):",allocated(qna%lst)," .neqv. ",allocated(cgl%cg%q)
               call die(msg)
            else if (allocated(qna%lst)) then
               if (size(qna%lst(:)) /= size(cgl%cg%q)) then
                  write(msg,'(2(a,i5))')"[cg_list_global:check_na] size(qna) /= size(cgl%cg%q)",size(qna%lst(:))," /= ",size(cgl%cg%q)
                  call die(msg)
               else
                  do i = lbound(qna%lst(:), dim=1, kind=4), ubound(qna%lst(:), dim=1, kind=4)
                     if (associated(cgl%cg%q(i)%arr) .and. cgl%cg%l%id < base_level_id .and. .not. qna%lst(i)%multigrid) then
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
                  do i = lbound(wna%lst(:), dim=1, kind=4), ubound(wna%lst(:), dim=1, kind=4)
                     bad = .false.
                     if (associated(cgl%cg%w(i)%arr)) bad = wna%get_dim4(i) /= size(cgl%cg%w(i)%arr, dim=1) .and. cgl%cg%l%id >= base_level_id
                     if (wna%get_dim4(i) <= 0 .or. bad) then
                        write(msg,'(a,i3,2a,2(a,i7))')"[cg_list_global:check_na] wna%lst(",i,"_ named '",wna%lst(i)%name,"' has inconsistent dim4: ",&
                             &         wna%get_dim4(i)," /= ",size(cgl%cg%w(i)%arr, dim=1)
                        call die(msg)
                     endif
                     if (associated(cgl%cg%w(i)%arr) .and. cgl%cg%l%id < base_level_id) then
                        write(msg,'(a,i3,3a)')"[cg_list_global:check_na] cgl%cg%w(",i,"), named '",wna%lst(i)%name,"' allocated on coarse level"
                        call die(msg)
                     endif
                  enddo
               endif
            endif
         endif
         cgl => cgl%nxt
      enddo

   end subroutine check_na

!>
!! \brief Find grid pieces that do not belong to any list except for all_cg
!!
!! \todo Convert warn() to die()
!<

   subroutine mark_orphans(this)

      use cg_list,          only: cg_list_element
      use constants,        only: INVALID
      use dataio_pub,       only: warn, msg
      use list_of_cg_lists, only: all_lists

      implicit none

      class(cg_list_global_t), intent(in) :: this          !< object invoking type-bound procedure

      type(cg_list_element), pointer :: cgl
      integer :: i
      integer, parameter :: VERY_INVALID = 2*INVALID
      integer, save :: na = 0, nm = 0, nf = 0, no = 0

      ! scan all lists except for all_cg for cg and set their membership to a bogus value
      do i = lbound(all_lists%entries(:), dim=1), ubound(all_lists%entries(:), dim=1)
         cgl => all_lists%entries(i)%lp%first
         if (all_lists%entries(i)%lp%label /= all_cg_n) then
            do while (associated(cgl))
               cgl%cg%membership = VERY_INVALID
               cgl => cgl%nxt
            enddo
         endif
      enddo

      ! mark all cg's with INVALID. If some aren't listed on all_cg then they should remain with %membership set to VERY_INVALID
      cgl => this%first
      do while (associated(cgl))
         if (associated(cgl%cg)) then
            cgl%cg%membership = INVALID
         else
            na = na + 1
         endif
         cgl => cgl%nxt
      enddo

      ! scan all lists except for all_cg
      do i = lbound(all_lists%entries(:), dim=1), ubound(all_lists%entries(:), dim=1)
         if (all_lists%entries(i)%lp%label /= all_cg_n) then
            cgl => all_lists%entries(i)%lp%first
            do while (associated(cgl))
               if (cgl%cg%membership == VERY_INVALID) then
                  write(msg, '(a,i7,a,i3,a)')"[cg_list_global:mark_orphans] Grid #",cgl%cg%grid_id, " at level ",cgl%cg%l%id," is hidden."
                  call warn(msg)
                  cgl%cg%membership = 0
                  nm = nm + 1
               endif
               if (cgl%cg%membership == INVALID) cgl%cg%membership = 0
               cgl%cg%membership = cgl%cg%membership + 1
               cgl => cgl%nxt
            enddo
         endif
      enddo

      ! Now search for not associated grid containers
      cgl => this%first
      do while (associated(cgl))
         if (associated(cgl%cg)) then
            if (cgl%cg%membership < 1) then
               call all_lists%forget(cgl%cg)
               nf = nf + 1
            endif
         endif
         cgl => cgl%nxt
      enddo

      cgl => this%first
      do while (associated(cgl))
         if (associated(cgl%cg)) then
            if (cgl%cg%membership < 1) then
               write(msg, '(a,i7,a,i3,a)')"[cg_list_global:mark_orphans] Grid #",cgl%cg%grid_id, " at level ",cgl%cg%l%id," is orphaned."
               call warn(msg)
               no = no + 1
            endif
         endif
         cgl => cgl%nxt
      enddo

      if (any([na, nm, nf, no] /= 0)) then
         write(msg, '(4(a,i6))')"[cg_list_global:mark_orphans] na = ", na, ", nm =", nm, ", nf = ", nf, ", no = ", no
         call warn(msg)
      endif

   end subroutine mark_orphans

end module cg_list_global
