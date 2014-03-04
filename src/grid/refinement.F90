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

!> \brief This module contains variables and initialization routines related to refinement

module refinement

   use constants,          only: ndims, LO, HI, cbuff_len
   use refinement_flag,    only: level_min, level_max
   use refinement_filters, only: ref_crit

   implicit none

   private
   public :: n_updAMR, oop_thr, refine_points, auto_refine_derefine, cleanup_refinement, &
        &    refine_boxes, init_refinement, emergency_fix, set_n_updAMR, strict_SFC_ordering, prefer_n_bruteforce, &
        &    refines2list, user_ref2list

   integer(kind=4), protected :: n_updAMR            !< how often to update the refinement structure
   logical,         protected :: strict_SFC_ordering !< Enforce strict SFC ordering to allow optimized neighbour search
   real,            protected :: oop_thr             !< Maximum allowed ratio of Out-of-Place grid pieces (according to current ordering scheme)
   logical,         protected :: prefer_n_bruteforce !< if .false. then try DFC algorithms for neighbor searches

   ! some refinement primitives
   integer, parameter :: nshapes = 10

   !> \brief Refinement point
   type :: ref_point
      integer(kind=4)        :: level  !> desired level of refinement
      real, dimension(ndims) :: coords !> coordinates, where to refine
   end type ref_point
   type(ref_point), dimension(nshapes), protected :: refine_points

   !> \brief Refinement box
   type :: ref_box
      integer(kind=4)               :: level  !> desired level of refinement
      real, dimension(ndims, LO:HI) :: coords !> coordinates, where to refine
   end type ref_box
   type(ref_box), dimension(nshapes), protected :: refine_boxes

   !> \brief Parameters of automagic refinement
   type :: ref_auto_param
      character(len=cbuff_len) :: rvar  !< name of the refinement variable
      character(len=cbuff_len) :: rname !< name of the refinement routine
      real :: ref_thr                   !< refinement threshold
      real :: deref_thr                 !< derefinement threshold
      real :: aux                       !< auxiliary parameter (can be smoother or filter strength)
   end type ref_auto_param
   integer, parameter :: n_ref_auto_param = 10
   type(ref_auto_param), dimension(10), protected :: refine_vars
   type(ref_crit), dimension(:), allocatable :: ref_crit_list

   character(len=cbuff_len), parameter :: inactive_name = "none"

   logical :: emergency_fix !< set to .true. if you want to call update_refinement ASAP

   namelist /AMR/ level_min, level_max, n_updAMR, strict_SFC_ordering, &
        &         prefer_n_bruteforce, oop_thr, refine_points, refine_boxes, refine_vars

contains

!> \brief Initialization of parameters of refinement mechanics

   subroutine init_refinement

      use constants,  only: base_level_id, PIERNIK_INIT_DOMAIN, xdim, ydim, zdim, I_ONE, LO, HI, cbuff_len
      use dataio_pub, only: nh      ! QA_WARN required for diff_nml
      use dataio_pub, only: die, code_progress, warn
      use domain,     only: AMR_bsize, dom
      use mpisetup,   only: cbuff, ibuff, lbuff, rbuff, master, slave, piernik_MPI_Bcast

      implicit none

      integer :: d
      logical :: allow_AMR

      if (code_progress < PIERNIK_INIT_DOMAIN) call die("[refinement:init_refinement] Domain not initialized.")

      level_min = base_level_id
      level_max = level_min
      n_updAMR  = huge(I_ONE)
      strict_SFC_ordering = .false.
      allow_AMR = .true.
      prefer_n_bruteforce = .false.
      oop_thr = 0.1
      do d = xdim, zdim
         if (dom%has_dir(d))  then
            if (AMR_bsize(d) < dom%nb) then
               if (allow_AMR .and. master) call warn("[refinement:init_refinement] Refinements disabled (AMR_bsize too small)")
               allow_AMR = .false.
            else
               if (mod(dom%n_d(d), AMR_bsize(d)) /= 0) then
                  if (allow_AMR .and. master) call warn("[refinement:init_refinement] Refinements disabled (domain not divisible by AMR_bsize)")
                  allow_AMR = .false.
               endif
            endif
         endif
      enddo
      refine_points(:) = ref_point(base_level_id-1, [ 0., 0., 0.] )
      refine_boxes (:) = ref_box  (base_level_id-1, reshape([ 0., 0., 0., 0., 0., 0.], [ndims, HI-LO+I_ONE] ) )
      refine_vars  (:) = ref_auto_param (inactive_name, inactive_name, 0., 0., 0.)

      if (1 + 9*nshapes +3*n_ref_auto_param > ubound(rbuff, dim=1)) call die("[refinement:init_refinement] increase rbuff size") ! should be detected at compile time but it is only a warning
      if (2*n_ref_auto_param > ubound(cbuff, dim=1)) call die("[refinement:init_refinement] increase cbuff size")
      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=AMR)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=AMR, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "AMR")
         read(nh%cmdl_nml,nml=AMR, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "AMR", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=AMR)
         close(nh%lun)
         call nh%compare_namelist()

         ! sanitizing
         if (allow_AMR) then
            level_min = max(level_min, base_level_id)
            level_max = max(level_max, level_min)
         else
            level_min = base_level_id
            level_max = base_level_id
            n_updAMR  = huge(I_ONE)
         endif
         where (.not. dom%has_dir(:)) AMR_bsize(:) = huge(I_ONE)

         cbuff(1                 :  n_ref_auto_param) = refine_vars(:)%rvar
         cbuff(1+n_ref_auto_param:2*n_ref_auto_param) = refine_vars(:)%rname

         ibuff(1) = level_min
         ibuff(2) = level_max
         ibuff(3) = n_updAMR
         ibuff(4        :3+  nshapes) = refine_points(:)%level
         ibuff(4+nshapes:3+2*nshapes) = refine_boxes (:)%level

         lbuff(1) = allow_AMR
         lbuff(2) = strict_SFC_ordering
         lbuff(3) = prefer_n_bruteforce

         rbuff(1) = oop_thr
         rbuff(2          :1+  nshapes) = refine_points(:)%coords(xdim)
         rbuff(2+  nshapes:1+2*nshapes) = refine_points(:)%coords(ydim)
         rbuff(2+2*nshapes:1+3*nshapes) = refine_points(:)%coords(zdim)
         rbuff(2+3*nshapes:1+4*nshapes) = refine_boxes (:)%coords(xdim, LO)
         rbuff(2+4*nshapes:1+5*nshapes) = refine_boxes (:)%coords(xdim, HI)
         rbuff(2+5*nshapes:1+6*nshapes) = refine_boxes (:)%coords(ydim, LO)
         rbuff(2+6*nshapes:1+7*nshapes) = refine_boxes (:)%coords(ydim, HI)
         rbuff(2+7*nshapes:1+8*nshapes) = refine_boxes (:)%coords(zdim, LO)
         rbuff(2+8*nshapes:1+9*nshapes) = refine_boxes (:)%coords(zdim, HI)
         rbuff(2+9*nshapes                   :1+9*nshapes+  n_ref_auto_param) = refine_vars(:)%ref_thr
         rbuff(2+9*nshapes+  n_ref_auto_param:1+9*nshapes+2*n_ref_auto_param) = refine_vars(:)%deref_thr
         rbuff(2+9*nshapes+2*n_ref_auto_param:1+9*nshapes+3*n_ref_auto_param) = refine_vars(:)%aux

      endif

      call piernik_MPI_Bcast(cbuff, cbuff_len)
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         refine_vars(:)%rvar  = cbuff(1                 :  n_ref_auto_param)
         refine_vars(:)%rname = cbuff(1+n_ref_auto_param:2*n_ref_auto_param)

         level_min = ibuff(1)
         level_max = ibuff(2)
         n_updAMR  = ibuff(3)
         refine_points(:)%level = ibuff(4        :3+  nshapes)
         refine_boxes (:)%level = ibuff(4+nshapes:3+2*nshapes)

         allow_AMR           = lbuff(1)
         strict_SFC_ordering = lbuff(2)
         prefer_n_bruteforce = lbuff(3)

         oop_thr = rbuff(1)
         refine_points(:)%coords(xdim)     = rbuff(2          :1+  nshapes)
         refine_points(:)%coords(ydim)     = rbuff(2+  nshapes:1+2*nshapes)
         refine_points(:)%coords(zdim)     = rbuff(2+2*nshapes:1+3*nshapes)
         refine_boxes (:)%coords(xdim, LO) = rbuff(2+3*nshapes:1+4*nshapes)
         refine_boxes (:)%coords(xdim, HI) = rbuff(2+4*nshapes:1+5*nshapes)
         refine_boxes (:)%coords(ydim, LO) = rbuff(2+5*nshapes:1+6*nshapes)
         refine_boxes (:)%coords(ydim, HI) = rbuff(2+6*nshapes:1+7*nshapes)
         refine_boxes (:)%coords(zdim, LO) = rbuff(2+7*nshapes:1+8*nshapes)
         refine_boxes (:)%coords(zdim, HI) = rbuff(2+8*nshapes:1+9*nshapes)
         refine_vars  (:)%ref_thr          = rbuff(2+9*nshapes                   :1+9*nshapes+  n_ref_auto_param)
         refine_vars  (:)%deref_thr        = rbuff(2+9*nshapes+  n_ref_auto_param:1+9*nshapes+2*n_ref_auto_param)
         refine_vars  (:)%aux              = rbuff(2+9*nshapes+2*n_ref_auto_param:1+9*nshapes+3*n_ref_auto_param)

      endif

      if (.not. allow_AMR) AMR_bsize=0

      ! Such large refinements may require additional work in I/O routines, visualization, computing MPI tags and so on.
      if (level_max > 40) call warn("[refinement:init_refinement] BEWARE: At such large refinements, integer overflows may happen under certain conditions.")

      emergency_fix = .false.

   end subroutine init_refinement

!> \brief free the memory

   subroutine cleanup_refinement

     implicit none

     if (allocated(ref_crit_list)) deallocate(ref_crit_list)

   end subroutine cleanup_refinement

!> \brief convert refinement criteria parameters to a list

   subroutine refines2list

      use constants,  only: INVALID
      use dataio_pub, only: die, msg

      implicit none

      integer :: i, c, iv
      integer, dimension(:), allocatable :: ic

      do i = 1, n_ref_auto_param
         call identify_field(refine_vars(i)%rvar, iv, ic)
         if (iv /= INVALID .and. trim(refine_vars(i)%rname) /= trim(inactive_name)) then
            if (.not. allocated(ic)) call die("[refinement:refines2list] .not. allocated(ic))")
            if (size(ic) <= 0) then
               write(msg,'(3a)')"[refinement:refines2list] iv /= INVALID and  unrecognized field name '",trim(refine_vars(i)%rname),"'"
               call die(msg)
            else
               do c = lbound(ic, dim=1), ubound(ic, dim=1)
                  call user_ref2list(iv, ic(c), refine_vars(i)%ref_thr, refine_vars(i)%deref_thr, refine_vars(i)%aux, refine_vars(i)%rname)
               enddo
            endif
         endif
         if (allocated(ic)) deallocate(ic)
      enddo

   end subroutine refines2list

!> \brief Add a user-defined criteria to the list

   subroutine user_ref2list(iv, ic, ref_thr, deref_thr, aux, rname)

      use constants,          only: INVALID
      use dataio_pub,         only: warn, die
      use refinement_filters, only: refine_on_gradient, refine_on_relative_gradient

      implicit none

      integer,          intent(in) :: iv        !< field index in cg%q or cg%w array
      integer,          intent(in) :: ic        !< component index of 4D array or INVALID for 3D arrays
      real,             intent(in) :: ref_thr   !< refinement threshold
      real,             intent(in) :: deref_thr !< derefinement threshold
      real,             intent(in) :: aux       !< auxiliary parameter
      character(len=*), intent(in) :: rname     !< name of the refinement routine

      if (iv == INVALID) then
         call warn("[refinement:user_ref2list] invalid field. Ignored.")
         return
      endif
      if (.not. allocated(ref_crit_list)) allocate(ref_crit_list(0))
      ref_crit_list = [ ref_crit_list, ref_crit(iv, ic, ref_thr, deref_thr, aux, null()) ]
      select case (trim(rname))
         case ("grad")
            ref_crit_list(ubound(ref_crit_list, dim=1))%refine => refine_on_gradient
         case ("relgrad")
            ref_crit_list(ubound(ref_crit_list, dim=1))%refine => refine_on_relative_gradient

!> \todo Implement Richardson extrapolation method, as described in M. Berger papers

!>
!! \todo Implement Loechner criteria
!! Original paper: https://www.researchgate.net/publication/222452974_An_adaptive_finite_element_scheme_for_transient_problems_in_CFD
!! Cartesian grid implementation: http://flash.uchicago.edu/~jbgallag/2012/flash4_ug/node14.html#SECTION05163100000000000000 (note that some indices in the left part of denominator are slightly messed up)
!<
         case (trim(inactive_name)) ! do nothing
         case default
            call die("[refinement:user_ref2list] unknown refinement detection routine")
      end select
      !> \todo try to detect doubled criteria

   end subroutine user_ref2list

!> \brief Change the protected parameter n_updAMR

   subroutine set_n_updAMR(n)

      use dataio_pub, only: die
      use mpisetup,   only: piernik_MPI_Bcast

      implicit none

      integer(kind=4), intent(in) :: n

      integer(kind=4) :: nn

      nn = n
      call piernik_MPI_Bcast(nn)
      if (nn /= n) call die("[refinement:set_n_updAMR] n_updAMR goes out of sync")
      n_updAMR = nn

   end subroutine set_n_updAMR

!>
!! \brief Apply automatic refinement citeria
!!
!! \details The leaves argument should normally be cg_leaves::leaves. W cannot use it directly here because of circular dependencies
!! Note that if you pass only subset of leaves (i.e. single level or somehow filtered cg list), the (de)refinement will be marked only on that list.
!! The rest of the domain will stay unaffected or be corrected for refinement defects.
!<

   subroutine auto_refine_derefine(leaves)

      use cg_list,          only: cg_list_element, cg_list_T
      use constants,        only: INVALID
      use dataio_pub,       only: die
      use named_array_list, only: qna, wna

      implicit none

      class(cg_list_T), intent(inout) :: leaves

      integer :: i
      logical :: var3d
      type(cg_list_element), pointer :: cgl
      real, dimension(:,:,:), pointer :: p3d

      if (.not. allocated(ref_crit_list)) return
      do i = lbound(ref_crit_list, dim=1), ubound(ref_crit_list, dim=1)
         var3d = (ref_crit_list(i)%ic == INVALID)

         if (var3d) then
            if (ref_crit_list(i)%iv<lbound(qna%lst, dim=1) .or. ref_crit_list(i)%iv>ubound(qna%lst, dim=1)) &
                 call die("[refinement:auto_refine_derefine] 3D index out of range")
         else
            if (ref_crit_list(i)%iv<lbound(wna%lst, dim=1) .or. ref_crit_list(i)%iv>ubound(wna%lst, dim=1)) &
                 call die("[refinement:auto_refine_derefine] 4D index out of range")
            if (ref_crit_list(i)%ic <= 0 .or. ref_crit_list(i)%ic > wna%lst(ref_crit_list(i)%iv)%dim4) &
                 call die("[refinement:auto_refine_derefine] component out of range")
         endif

         cgl => leaves%first
         do while (associated(cgl))
            if (any(cgl%cg%leafmap)) then
               if (var3d) then
                  p3d => cgl%cg%q(ref_crit_list(i)%iv)%arr
               else
                  associate (a=>cgl%cg%w(ref_crit_list(i)%iv)%arr)
                     p3d(lbound(a, dim=2):, lbound(a, dim=3):, lbound(a, dim=4):) => cgl%cg%w(ref_crit_list(i)%iv)%arr(ref_crit_list(i)%ic, :, :, :)
                  end associate
               endif
               call ref_crit_list(i)%refine(cgl%cg, p3d)
            endif
            cgl => cgl%nxt
         enddo
      enddo

   end subroutine auto_refine_derefine

!> \brief Identify field name and return indices to cg%q or cg%w arrays

   subroutine identify_field(vname, iv, ic)

      use constants,        only: INVALID, cbuff_len
      use dataio_pub,       only: msg, warn
      use fluidindex,       only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz, iarr_all_en
      use named_array_list, only: qna, wna

      implicit none

      character(len=cbuff_len),           intent(in)  :: vname !< string specifying the field on
      integer,                            intent(out) :: iv    !< field index in cg%q or cg%w array
      integer, dimension(:), allocatable, intent(out) :: ic    !< component index array (cg%w(iv)%arr(ic,:,:,:)) or INVALID for 3D arrays

      iv = INVALID

      if (trim(vname) == trim(inactive_name)) return ! ignore this

      if (qna%exists(trim(vname))) then
         iv = qna%ind(trim(vname))
         if (iv /= INVALID) then
            allocate(ic(1))
            ic = INVALID
            return ! this is a 3d array name
         endif
      endif

      if (trim(vname) == "dens") then
         allocate(ic(lbound(iarr_all_dn, dim=1):ubound(iarr_all_dn, dim=1)))
         iv = wna%fi
         ic = iarr_all_dn
         return
      else if (trim(vname) == "velx") then
         allocate(ic(lbound(iarr_all_mx, dim=1):ubound(iarr_all_mx, dim=1)))
         iv = wna%fi
         ic = iarr_all_mx
         return
      else if (trim(vname) == "vely") then
         allocate(ic(lbound(iarr_all_my, dim=1):ubound(iarr_all_my, dim=1)))
         iv = wna%fi
         ic = iarr_all_my
         return
      else if (trim(vname) == "velz") then
         allocate(ic(lbound(iarr_all_mz, dim=1):ubound(iarr_all_mz, dim=1)))
         iv = wna%fi
         ic = iarr_all_mz
         return
      else if (trim(vname) == "ener") then
         allocate(ic(lbound(iarr_all_en, dim=1):ubound(iarr_all_en, dim=1)))
         iv = wna%fi
         ic = iarr_all_en
         return
      endif
      !> \todo identify here all {den,vl[xyz],ene}{d,n,i}
      !> \todo introduce possibility to operate on pressure or other indirect fields

      write(msg,'(3a)')"[refinement:identify_field] Unidentified refinement variable: '",trim(vname),"'"
      call warn(msg)

      allocate(ic(0))

   end subroutine identify_field

end module refinement
