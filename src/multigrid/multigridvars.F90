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

!!$ ============================================================================
!>
!! \brief Variables and data structures required by multigrid routines.
!!
!! \details These variables are not meant to be accessed (or worse: altered) from the outside of multigrid routines (Fortran has no friends :-( )
!<

module multigridvars
! pulled by MULTIGRID

   use constants, only: xdim, zdim, LO, HI
   use grid_cont, only: grid_container

   implicit none

   public ! QA_WARN no secrets are kept here
   private :: xdim, zdim, LO, HI, grid_container ! QA_WARN prevent re-exporting

   ! multigrid constants
   enum, bind(C)
      enumerator :: source = 1                                        !< Index of the density field
      enumerator :: solution                                          !< Index of the iterated solution (potential) fields
      enumerator :: defect                                            !< Index of the defect field (effectively the density not accounted in current solution)
      enumerator :: correction                                        !< Index of the correction to the potential to be applied at the end of V-cycle
   end enum

   ! these constants should be moved to constants module

   ! namelist parameters
   integer            :: ord_prolong                                  !< Prolongation operator order; allowed values are -4 -2, 0 (default), 2 and 4; -2 is often fast
   integer            :: ord_prolong_face_norm                        !< Face prolongation operator order in the direction normal to the face; allowed values are 0, 1  and 2
   integer            :: ord_prolong_face_par                         !< Face prolongation operator order in the directions parallel to the face; allowed values are -2 .. 2
   logical            :: stdout                                       !< print verbose messages to stdout
   logical            :: verbose_vcycle                               !< Print one line of log per V-cycle, summary otherwise

   ! dimensions
   integer(kind=4) :: ngridvars                                       !< number of variables required for implementation of multigrid

   ! boundaries
   enum, bind(C)                                                      !< constants for enumerating multigrid boundary types
      enumerator :: bnd_periodic                                      !< periodic
      enumerator :: bnd_dirichlet                                     !< 0-value boundary type (uniform Dirichlet)
      enumerator :: bnd_isolated                                      !< isolated boundary type
      enumerator :: bnd_neumann                                       !< 0-gradient boundary type (uniform Neumann)
      enumerator :: bnd_givenval                                      !< given value boundary type (general Dirichlet)
      enumerator :: bnd_invalid = bnd_periodic - 1                    !< invalid
   end enum
   logical, dimension(xdim:zdim, LO:HI) :: is_external                !< .true. for non-"mpi" local domain boundaries
   enum, bind(C)
      enumerator :: extbnd_donothing                                  !< Do not touch external boundaries
      enumerator :: extbnd_zero                                       !< Fill external boundaries with zeroes
      enumerator :: extbnd_extrapolate                                !< Perform extrapolation in external boundaries
      enumerator :: extbnd_mirror                                     !< Zero-gradient, mirroring external boundaries
      enumerator :: extbnd_antimirror = - extbnd_mirror               !< mirroring external boundaries with opposite sign
   end enum

   ! miscellaneous
   real                    :: ts                                      !< time for runtime profiling
   real                    :: tot_ts                                  !< total multigrid time
   logical                 :: is_mg_uneven                            !< .true. when domain shapes differ across procesors, even on the coarsest grids
   logical                 :: single_base                             !< .true. when the whole base level is located on a single cpu
   logical                 :: need_general_pf                         !< .false. only for most regular domain decomposition

   integer, parameter :: prefix_len = 3                               !< length of prefix for distinguishing V-cycles in the log
   type :: vcycle_stats
      real, allocatable, dimension(:) :: factor                       !< norm reduction factor
      real, allocatable, dimension(:) :: time                         !< time spent
      integer                         :: count                        !< number of executed V-cycles
      real                            :: norm_rhs                     !< norm of the source
      real                            :: norm_final                   !< norm of the defect relative to the source
      character(len=prefix_len)       :: cprefix                      !< prefix for distinguishing V-cycles in the log (e.g inner or outer potential, CR component)
   end type vcycle_stats

   type, extends(grid_container) :: plvl                              !< single level container

     type(plvl), pointer :: finer, coarser                           !< pointers to level+1 and level-1 (null() if such level does not eist)

    contains

      procedure :: restrict_level
      procedure :: prolong_level0

   end type plvl

   type(plvl), dimension(:), allocatable, target :: lvl               !< a stack of multigrid arrays
   type(plvl), pointer                           :: base              !< pointer to coarsest level
   type(plvl), pointer                           :: roof              !< pointer to finest level

contains

!!$ ============================================================================
!>
!! \brief Simplest restriction (averaging).
!! \todo implement high order restriction and test its influence on V-cycle convergence rate
!<

   subroutine restrict_level(this, iv)

      use constants,  only: xdim, ydim, zdim, LO, HI, LONG, I_ONE
      use dataio_pub, only: msg, warn, die
      use domain,     only: dom
      use mpisetup,   only: proc, comm, ierr, req, status
      use mpi,        only: MPI_DOUBLE_PRECISION

      implicit none

      class(plvl), intent(inout), target  :: this  !< object invoking type-bound procedure
      integer(kind=4), intent(in)         :: iv    !< variable to be restricted

      integer(kind=4), parameter :: tag1 = 1
      class(plvl), pointer :: coarse
      integer :: g, g1, d
      integer(kind=8), dimension(:,:), pointer :: fse, cse ! shortcuts for fine segment and coarse segment
      integer(kind=8) :: i, j, k, ic, jc, kc
      integer(kind=8), dimension(xdim:zdim) :: off1
      real :: norm
      integer(kind=4) :: nr

      if (iv < lbound(this%mg%var(:,:,:,:), dim=4) .or. iv > ubound(this%mg%var(:,:,:,:), dim=4)) call die("[multigridvars:restrict_level] Invalid variable index.")

      coarse => this%coarser
      if (.not. associated(coarse)) then
         write(msg,'(a,i3)')"[multigridvars:restrict_level] no coarse level here: ", this%mg%level
         call warn(msg) ! can't restrict base level
      else
         !OPT find a way to reduce this to areas with nonlocal incoming restriction
         coarse%mg%var(:,:,:, iv) = 0. ! disables check_dirty
      endif

      nr = 0
      if (allocated(coarse%mg%f_tgt%seg)) then
         do g = 1, ubound(coarse%mg%f_tgt%seg(:), dim=1)
            if (coarse%mg%f_tgt%seg(g)%proc /= proc) then
               nr = nr + I_ONE
               call MPI_Irecv(coarse%mg%f_tgt%seg(g)%buf(1, 1, 1), size(coarse%mg%f_tgt%seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, coarse%mg%f_tgt%seg(g)%proc, tag1, comm, req(nr), ierr)
            endif
         enddo
      endif

      do g = 1, ubound(this%mg%c_tgt%seg(:), dim=1)

         fse => this%mg%c_tgt%seg(g)%se

         off1(:) = mod(this%off(:), 2_LONG)
         norm = 1./2**dom%eff_dim
         if (this%mg%c_tgt%seg(g)%proc == proc) then

            nullify(cse)
            do g1 = 1, ubound(coarse%mg%f_tgt%seg(:), dim=1) !> \todo should be set up in init_multigrid
               if (coarse%mg%f_tgt%seg(g1)%proc == proc) then
                  if (.not. associated(cse)) then
                     cse => coarse%mg%f_tgt%seg(g1)%se
                  else
                     call die("[multigridvars:restrict_level] multiple local coarse grid targets")
                  endif
               endif
            enddo
            if (.not. associated(cse)) call die("[multigridvars:restrict_level] cannot find local coarse grid")

            do d = xdim, zdim ! debug
               if (dom%has_dir(d)) then
                  if ((fse(d, LO)-off1(d)-this%nb-1)/2+cse(d, LO) < cse(d, LO)) call die("mv:rl <cse")
                  if ((fse(d, HI)-off1(d)-this%nb-1)/2+cse(d, LO) > cse(d, HI)) call die("mv:rl >cse")
               endif
            enddo
            ! OPT: completely unoptimized,
            ! note that e.g. 10 cells on fine grid may contribute to 5 or 6 cells on coarse grid, depending on offset and both cases require a bit different code
            ! \todo use array sections to restrict the fully covered interior of coarse segment and finish the boundaries, where necessary
            ! OPT: old code operations: Ir:Dr:Dw:Dr_m:Dw_m = 13:8:1:1:0.1, this code 130:20:8:1.4:0.2
            ! OPT: computation of ic consumes ~30% Ir
            ! OPT: convert the loop over i (at least) to array section operation
            ! previous implementation of restriction (sum of array sections without loops) was probably faster and did produce slightly different results
            ! (with typical relative difference at level 1e-15 due to different numerical roundoffs)
            ! It is possible that future optimisations will affect the results and bit-to-bit comparisions (with h5diff) will fail
            do k = fse(zdim, LO), fse(zdim, HI)
               kc = (k-off1(zdim)-this%ks)/2+cse(zdim, LO)
               do j = fse(ydim, LO), fse(ydim, HI)
                  jc = (j-off1(ydim)-this%js)/2+cse(ydim, LO)
                  do i = fse(xdim, LO), fse(xdim, HI)
                     ic = (i-off1(xdim)-this%is)/2+cse(xdim, LO)
                     coarse%mg%var(ic, jc, kc, iv) = coarse%mg%var(ic, jc, kc, iv) + this%mg%var(i, j, k, iv) * norm
                  enddo
               enddo
            enddo
            !\todo add geometrical terms to improve convergence on cylindrical grids
         else
            ! OPT: see the notes above
            this%mg%c_tgt%seg(g)%buf(:, :, :) = 0.
            off1(:) = mod(this%off(:)+fse(:, LO) - this%ijkse(:, LO), 2_LONG)
            do k = fse(zdim, LO), fse(zdim, HI)
               kc = (k-fse(zdim, LO)+off1(zdim))/2 + 1
               do j = fse(ydim, LO), fse(ydim, HI)
                  jc = (j-fse(ydim, LO)+off1(ydim))/2 + 1
                  do i = fse(xdim, LO), fse(xdim, HI)
                     ic = (i-fse(xdim, LO)+off1(xdim))/2 + 1
                     this%mg%c_tgt%seg(g)%buf(ic, jc, kc) = this%mg%c_tgt%seg(g)%buf(ic, jc, kc) + this%mg%var(i, j, k, iv) * norm
                  enddo
               enddo
            enddo
            nr = nr + I_ONE
            call MPI_Isend(this%mg%c_tgt%seg(g)%buf(1, 1, 1), size(this%mg%c_tgt%seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, this%mg%c_tgt%seg(g)%proc, tag1, comm, req(nr), ierr)
         endif
      enddo

      if (nr>0) call MPI_Waitall(nr, req(:nr), status(:,:nr), ierr)

      if (allocated(coarse%mg%f_tgt%seg)) then
         do g = 1, ubound(coarse%mg%f_tgt%seg(:), dim=1)
            if (coarse%mg%f_tgt%seg(g)%proc /= proc) then
               cse => coarse%mg%f_tgt%seg(g)%se
               coarse%mg%var     (cse(xdim, LO):cse(xdim, HI), cse(ydim, LO):cse(ydim, HI), cse(zdim, LO):cse(zdim, HI), iv) = &
                    coarse%mg%var(cse(xdim, LO):cse(xdim, HI), cse(ydim, LO):cse(ydim, HI), cse(zdim, LO):cse(zdim, HI), iv) + coarse%mg%f_tgt%seg(g)%buf(:, :, :)
            endif
         enddo
      endif

   end subroutine restrict_level

!!$ ============================================================================
!>
!! \brief 0th order prolongation : injection
!<

   subroutine prolong_level0(this, iv)

      use constants,  only: xdim, ydim, zdim, LO, HI, LONG, I_ONE
      use dataio_pub, only: msg, warn, die
      use domain,     only: dom
      use mpisetup,   only: proc, comm, ierr, req, status
      use mpi,        only: MPI_DOUBLE_PRECISION

      implicit none

      class(plvl), intent(inout), target  :: this  !< object invoking type-bound procedure
      integer(kind=4), intent(in)         :: iv    !< variable to be prolonged

      integer(kind=4), parameter :: tag1 = 1
      class(plvl), pointer :: fine
      integer :: g, g1, d
      integer(kind=8), dimension(:,:), pointer :: fse, cse ! shortcuts for fine segment and coarse segment
      integer(kind=8) :: i, j, k, ic, jc, kc
      integer(kind=8), dimension(xdim:zdim) :: off1
      integer(kind=4) :: nr

      if (iv < lbound(this%mg%var(:,:,:,:), dim=4) .or. iv > ubound(this%mg%var(:,:,:,:), dim=4)) call die("[multigridvars:prolong_level0] Invalid variable index.")

      fine => this%finer
      if (.not. associated(fine)) then
         write(msg,'(a,i3)')"[multigridvars:restrict_level] no fine level here: ", this%mg%level
         call warn(msg) ! can't prolong finest level
      else
         ! OPT: try to remove or limit this  ~20% Ir, ~50% Dw_m
         fine%mg%var(:,:,:, iv) = 0. ! disables check_dirty
      endif

      nr = 0
      if (allocated(fine%mg%c_tgt%seg)) then
         do g = 1, ubound(fine%mg%c_tgt%seg(:), dim=1)
            if (fine%mg%c_tgt%seg(g)%proc /= proc) then
               nr = nr + I_ONE
               call MPI_Irecv(fine%mg%c_tgt%seg(g)%buf(1, 1, 1), size(fine%mg%c_tgt%seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, fine%mg%c_tgt%seg(g)%proc, tag1, comm, req(nr), ierr)
            endif
         enddo
      endif

      off1(:) = mod(fine%off(:), 2_LONG)
      do g = 1, ubound(this%mg%f_tgt%seg(:), dim=1)

         cse => this%mg%f_tgt%seg(g)%se

         if (this%mg%f_tgt%seg(g)%proc == proc) then

            nullify(fse)
            do g1 = 1, ubound(fine%mg%c_tgt%seg(:), dim=1) !> \todo should be set up in init_multigrid
               if (fine%mg%c_tgt%seg(g1)%proc == proc) then
                  if (.not. associated(fse)) then
                     fse => fine%mg%c_tgt%seg(g1)%se
                  else
                     call die("[multigridvars:prolong_level0] multiple local fine grid targets")
                  endif
               endif
            enddo
            if (.not. associated(fse)) call die("[multigridvars:prolong_level0] cannot find local fine grid")

            ! No guardcells required here

            ! Possible optimization candidate: reduce L1 and L2 cache misses on both read and write (RBGS only, secondary importance)
            do d = xdim, zdim ! debug
               if (dom%has_dir(d)) then
                  if ((fse(d, LO)-off1(d)-this%nb-1)/2+cse(d, LO) < cse(d, LO)) call die("mv:rl <cse")
                  if ((fse(d, HI)-off1(d)-this%nb-1)/2+cse(d, LO) > cse(d, HI)) call die("mv:rl >cse")
               endif
            enddo
            ! OPT: completely unoptimized
            ! OPT: old code operations: Ir:Dr:Dw:Dr_m:Dw_m = 37:9:11:1.5:1.5, new code = 72:10:9:0.1:1.7
            ! OPT: computation of ic consumes ~25% Ir
            do k = fse(zdim, LO), fse(zdim, HI)
               kc = (k-off1(zdim)-this%ks)/2+cse(zdim, LO)
               do j = fse(ydim, LO), fse(ydim, HI)
                  jc = (j-off1(ydim)-this%js)/2+cse(ydim, LO)
                  do i = fse(xdim, LO), fse(xdim, HI)
                     ic = (i-off1(xdim)-this%is)/2+cse(xdim, LO)
                     fine%mg%var(i, j, k, iv) = this%mg%var(ic, jc, kc, iv)
                  enddo
               enddo
            enddo
         else
            nr = nr + I_ONE
            this%mg%f_tgt%seg(g)%buf(:, :, :) = this%mg%var(cse(xdim, LO):cse(xdim, HI), cse(ydim, LO):cse(ydim, HI), cse(zdim, LO):cse(zdim, HI), iv)
            call MPI_Isend(this%mg%f_tgt%seg(g)%buf(1, 1, 1), size(this%mg%f_tgt%seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, this%mg%f_tgt%seg(g)%proc, tag1, comm, req(nr), ierr)
         endif

      enddo

      if (nr>0) call MPI_Waitall(nr, req(:nr), status(:,:nr), ierr)

      if (allocated(fine%mg%c_tgt%seg)) then
         do g = 1, ubound(fine%mg%c_tgt%seg(:), dim=1)
            if (fine%mg%c_tgt%seg(g)%proc /= proc) then
               fse => fine%mg%c_tgt%seg(g)%se
               off1(:) = mod(fine%off(:)+fse(:, LO) - this%ijkse(:, LO), 2_LONG)
               do k = fse(zdim, LO), fse(zdim, HI)
                  kc = (k-fse(zdim, LO)+off1(zdim))/2 + 1
                  do j = fse(ydim, LO), fse(ydim, HI)
                     jc = (j-fse(ydim, LO)+off1(ydim))/2 + 1
                     do i = fse(xdim, LO), fse(xdim, HI)
                        ic = (i-fse(xdim, LO)+off1(xdim))/2 + 1
                        fine%mg%var(i, j, k, iv) = fine%mg%c_tgt%seg(g)%buf(ic, jc, kc)
                     enddo
                  enddo
               enddo
            endif
         enddo
      endif

   end subroutine prolong_level0

end module multigridvars
