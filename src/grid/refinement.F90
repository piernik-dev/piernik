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

   use constants,       only: ndims, LO, HI, cbuff_len
   use refinement_flag, only: level_min, level_max

   implicit none

   private
   public :: n_updAMR, allow_face_rstep, allow_corner_rstep, oop_thr, refine_points, auto_refine_derefine, &
        &    refine_boxes, init_refinement, emergency_fix, set_n_updAMR, strict_SFC_ordering, prefer_n_bruteforce, &
        &    refines2list, user_ref2list

   integer(kind=4), protected :: n_updAMR            !< how often to update the refinement structure
   logical,         protected :: allow_face_rstep    !< Allows >1 refinement step across faces (do not use it for any physical problems)
   logical,         protected :: allow_corner_rstep  !< Allows >1 refinement step across edges and corners (do not use it for any physical problems)
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

   type :: ref_crit
      integer              :: iv !< field index in cg%q or cg%w array
      integer              :: ic !< component index (cg%w(iv)%arr(ic,:,:,:)) or INVALID for 3D arrays
      real                 :: ref_thr                   !< refinement threshold
      real                 :: deref_thr                 !< derefinement threshold
      real                 :: aux                       !< auxiliary parameter (can be smoother or filter strength)
      procedure(refine_crit), pass, pointer :: refine
   end type ref_crit
   type(ref_crit), dimension(:), allocatable :: ref_crit_list

   interface
      subroutine refine_crit(this, cg, p3d)

         use grid_cont, only: grid_container
         import ref_crit

         implicit none

         class(ref_crit),                 intent(in)    :: this !< this contains refinement parameters
         type(grid_container), pointer,   intent(inout) :: cg   !< current grid piece
         real, dimension(:,:,:), pointer, intent(in)    :: p3d  !< pointer to array to be examined for (de)refinement
      end subroutine refine_crit
   end interface

   character(len=cbuff_len), parameter :: inactive_name = "none"

   logical :: emergency_fix !< set to .true. if you want to call update_refinement ASAP

   namelist /AMR/ level_min, level_max, n_updAMR, allow_face_rstep, allow_corner_rstep, strict_SFC_ordering, &
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
      allow_face_rstep    = .false.
      allow_corner_rstep  = .false.
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
         where (.not. dom%has_dir(:)) AMR_bsize(:) = huge(1)

         cbuff(1                 :  n_ref_auto_param) = refine_vars(:)%rvar
         cbuff(1+n_ref_auto_param:2*n_ref_auto_param) = refine_vars(:)%rname

         ibuff(1) = level_min
         ibuff(2) = level_max
         ibuff(3) = n_updAMR
         ibuff(4        :3+  nshapes) = refine_points(:)%level
         ibuff(4+nshapes:3+2*nshapes) = refine_boxes (:)%level

         lbuff(1) = allow_face_rstep
         lbuff(2) = allow_corner_rstep
         lbuff(3) = allow_AMR
         lbuff(4) = strict_SFC_ordering
         lbuff(5) = prefer_n_bruteforce

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

         allow_face_rstep    = lbuff(1)
         allow_corner_rstep  = lbuff(2)
         allow_AMR           = lbuff(3)
         strict_SFC_ordering = lbuff(4)
         prefer_n_bruteforce = lbuff(5)

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

!> convert refinement criteria parameters to a list

   subroutine refines2list

      use constants,  only: INVALID
      use dataio_pub, only: die

      implicit none

      integer :: i, iv, ic

      do i = 1, n_ref_auto_param
         call identify_field(refine_vars(i)%rvar, iv, ic)
         if (iv /= INVALID .and. trim(refine_vars(i)%rname) /= trim(inactive_name)) &
              call user_ref2list(iv, ic, refine_vars(i)%ref_thr, refine_vars(i)%deref_thr, refine_vars(i)%aux, refine_vars(i)%rname)
      enddo

   end subroutine refines2list

!> \brief Add a user-defined criteria to the list

   subroutine user_ref2list(iv, ic, ref_thr, deref_thr, aux, rname)

      use constants,  only: INVALID
      use dataio_pub, only: warn, die

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
!> \todo Implement ||grad u|| normalized by average(|u|)

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
      use named_array_list, only: qna

      implicit none

      character(len=cbuff_len), intent(in)  :: vname !< string specifying the field on
      integer,                  intent(out) :: iv    !< field index in cg%q or cg%w array
      integer,                  intent(out) :: ic    !< component index (cg%w(iv)%arr(ic,:,:,:)) or INVALID for 3D arrays
      iv = INVALID
      ic = INVALID

      if (trim(vname) == trim(inactive_name)) return ! ignore this

      iv = qna%ind(trim(vname))
      if (iv /= INVALID) return ! this is a 3d array name

      !> \todo identify here all {den,vl[xyz],ene}{d,n,i}
      !> \todo interpret dens, vel[xyz], ener as multiple components to be checked

      write(msg,'(3a)')"[refinement:identify_field] Unidentified refinement variable: '",trim(vname),"'"

   end subroutine identify_field

!>
!! \brief Refine/derefine based on ||grad u||
!! This is sensitive to gradients, but the thresholds must be rescaled, when you change units of the problem.
!<

   subroutine refine_on_gradient(this, cg, p3d)

      use domain,    only: dom
      use grid_cont, only: grid_container

      implicit none

      class(ref_crit),                 intent(in)    :: this !< this contains refinement parameters
      type(grid_container), pointer,   intent(inout) :: cg   !< current grid piece
      real, dimension(:,:,:), pointer, intent(in)    :: p3d  !< pointer to array to be examined for (de)refinement needs (should contain at least one layer of updated guardcells)

      integer :: i, j, k
      real :: r, max_r

      max_r = -huge(1.)

      !> \todo implement how far we should look for (de)refinements

      do k = cg%ks, cg%ke
         do j = cg%js, cg%je
            do i = cg%is, cg%ie
               r = grad2(i, j, k)
               max_r = max(max_r, r)
               cg%refinemap(i, j, k) = cg%refinemap(i, j, k) .or. (r >= this%ref_thr**2)
               ! we can avoid calculating square root here
            enddo
         enddo
      enddo

      ! check additional 1 perimeter of cells for derefinement
      max_r = max(max_r, &
           &      maxval(grad2([cg%is-2*dom%D_x, cg%is-dom%D_x, cg%ie+dom%D_x, cg%ie+2*dom%D_x], &
           &                   [(j, j=cg%js-2*dom%D_y, cg%je+2*dom%D_y)], &
           &                   [(k, k=cg%ks-2*dom%D_z, cg%ke+2*dom%D_z)])), &
           &      maxval(grad2([(i, i=cg%is, cg%ie)], &
           &                   [cg%js-2*dom%D_y, cg%js-dom%D_y, cg%je+dom%D_y, cg%je+2*dom%D_y], &
           &                   [(k, k=cg%ks-2*dom%D_z, cg%ke+2*dom%D_z)])), &
           &      maxval(grad2([(i, i=cg%is, cg%ie)], &
           &                   [(j, j=cg%js, cg%je)], &
           &                   [cg%ks-2*dom%D_z, cg%ks-dom%D_z, cg%ke+dom%D_z, cg%ke+2*dom%D_z])) )

      cg%refine_flags%derefine = cg%refine_flags%derefine .or. (max_r < this%deref_thr**2)

   contains

      elemental real function grad2(i, j, k)

         use domain, only: dom

         implicit none

         integer,                         intent(in) :: i, j, k !< indices
!         real, dimension(:,:,:), pointer, intent(in) :: p3d     !< pointer to array

         grad2 = (p3d(i+dom%D_x, j, k) - p3d(i-dom%D_x, j, k))**2 + &
              &  (p3d(i, j+dom%D_y, k) - p3d(i, j-dom%D_y, k))**2 + &
              &  (p3d(i, j, k+dom%D_z) - p3d(i, j, k-dom%D_z))**2

      end function grad2

   end subroutine refine_on_gradient

end module refinement
