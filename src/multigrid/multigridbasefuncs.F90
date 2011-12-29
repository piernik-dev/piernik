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
!! \brief This module contains various routines (interpolation, boundaries and some global reduction)
!! that are useful for all flavours of multigrid solvers.
!<

module multigridbasefuncs
! pulled by MULTIGRID
   implicit none

   private

   public :: prolong_level, restrict_all, norm_sq, subtract_average, prolong_faces, zero_boundaries

contains

!!$ ============================================================================
!>
!! \brief Clear boundary values
!<

   subroutine zero_boundaries(curl)

      use gc_list,     only: cg_list_element
      use cg_list_lev, only: cg_list_level

      implicit none

      type(cg_list_level), pointer, intent(in) :: curl  !< level for which clear the boundary values

      type(cg_list_element), pointer :: cgl

      cgl => curl%first
      do while (associated(cgl))
         cgl%cg%mg%bnd_x(:,:,:) = 0.
         cgl%cg%mg%bnd_y(:,:,:) = 0.
         cgl%cg%mg%bnd_z(:,:,:) = 0.
         cgl => cgl%nxt
      enddo

   end subroutine zero_boundaries

!!$ ============================================================================
!>
!! \brief Multigrid elementary operators: prolongation, restriction, norm etc.
!<

   subroutine prolong_level(coarse, iv)

      use dataio_pub,            only: die
      use cg_list_lev,           only: cg_list_level
      use multigridhelpers,      only: dirty_debug, check_dirty, dirtyH
      use multigridvars,         only: roof, ord_prolong, is_mg_uneven
      use multigridexperimental, only: prolong_level_hord

      implicit none

      type(cg_list_level), pointer, intent(in) :: coarse !< level to prolong from
      integer(kind=4), intent(in) :: iv    !< variable to be prolonged

      type(cg_list_level), pointer :: fine

      if (associated(coarse, roof)) return ! can't prolong finest level

      if (.not. associated(coarse)) call die("[multigridbasefuncs:prolong_level] coarse == null()")
      fine   => coarse%finer
      if (.not. associated(fine)) call die("[multigridbasefuncs:prolong_level] fine == null()")

      !! \warning: need set_dirty_level routine
      if (dirty_debug) fine%first%cg%q(iv)%arr(:, :, :) = dirtyH

      call check_dirty(coarse, iv, "prolong-")

      if (ord_prolong == 0 .or. is_mg_uneven) then
         call coarse%prolong_level0(iv)
      else
         call prolong_level_hord(coarse, iv) ! experimental part
      endif

      call check_dirty(fine, iv, "prolong+")

   end subroutine prolong_level


!!$ ============================================================================
!>
!! \brief Restriction operators
!<

   subroutine restrict_all(iv)

      use cg_list_lev,      only: cg_list_level
      use multigridhelpers, only: check_dirty
      use multigridvars,    only: roof, base

      implicit none

      integer(kind=4), intent(in)      :: iv    !< variable to be restricted

      type(cg_list_level), pointer :: curl

      call check_dirty(roof, iv, "restrict_all-")

      curl => roof
      do while (associated(curl) .and. .not. associated(curl, base))
         call curl%restrict_level(iv)
         curl => curl%coarser
      enddo

      call check_dirty(base, iv, "restrict_all+")

   end subroutine restrict_all

!!$ ============================================================================
!>
!! \brief Calculate L2 norm
!<

   subroutine norm_sq(iv, norm)

      use constants,  only: GEO_XYZ, GEO_RPZ, I_ONE
      use dataio_pub, only: die
      use domain,     only: dom
      use gc_list,    only: cg_list_element
      use grid,       only: leaves
      use mpi,        only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_IN_PLACE
      use mpisetup,   only: comm, ierr

      implicit none

      integer(kind=4), intent(in)  :: iv   !< index of variable in cg%q(:) for which we want to find the norm
      real,            intent(out) :: norm !< the calculated norm

      integer :: i
      type(cg_list_element), pointer :: cgl

      norm = 0.
      cgl => leaves%first
      do while (associated(cgl))
         select case (dom%geometry_type)
            case (GEO_XYZ)
               norm = norm + sum(cgl%cg%q(iv)%arr(cgl%cg%is:cgl%cg%ie, cgl%cg%js:cgl%cg%je, cgl%cg%ks:cgl%cg%ke)**2) * cgl%cg%dvol
            case (GEO_RPZ)
               do i = cgl%cg%is, cgl%cg%ie
                  norm = norm + sum(cgl%cg%q(iv)%arr(i, cgl%cg%js:cgl%cg%je, cgl%cg%ks:cgl%cg%ke)**2) * cgl%cg%dvol * cgl%cg%x(i)
               enddo
            case default
               call die("[multigridbasefuncs:norm_sq] Unsupported geometry.")
         end select
         cgl => cgl%nxt
      enddo
      call MPI_Allreduce(MPI_IN_PLACE, norm, I_ONE, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
      norm = sqrt(norm)

   end subroutine norm_sq

!!$ ============================================================================
!>
!! \brief Compute the global average value and subtract it from the whole domain
!<

   subroutine subtract_average(curl, iv)

      use constants,   only: GEO_XYZ, GEO_RPZ, I_ONE
      use dataio_pub,  only: die
      use domain,      only: dom
      use gc_list,     only: cg_list_element
      use cg_list_lev, only: cg_list_level
      use mpi,         only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_IN_PLACE
      use mpisetup,    only: comm, ierr

      implicit none

      type(cg_list_level), pointer, intent(in) :: curl !< level for which we want to subtract its average from
      integer(kind=4),              intent(in) :: iv   !< index of variable in cg%q(:) which we want to have zero average

      real                :: avg, vol
      integer             :: i
      type(cg_list_element), pointer :: cgl

      avg = 0.
      vol = 0.
      cgl => curl%first
      do while (associated(cgl))
         select case (dom%geometry_type)
            case (GEO_XYZ)
               avg = avg + sum(cgl%cg%q(iv)%arr(cgl%cg%is:cgl%cg%ie, cgl%cg%js:cgl%cg%je, cgl%cg%ks:cgl%cg%ke)) * cgl%cg%dvol
            case (GEO_RPZ)
               do i = cgl%cg%is, cgl%cg%ie
                  avg = avg + sum(cgl%cg%q(iv)%arr(i, cgl%cg%js:cgl%cg%je, cgl%cg%ks:cgl%cg%ke)) * cgl%cg%dvol * cgl%cg%x(i)
               enddo
            case default
               call die("[multigridbasefuncs:subtract_average] Unsupported geometry.")
         end select
         vol = vol + cgl%cg%vol
         cgl => cgl%nxt
      enddo
      call MPI_Allreduce(MPI_IN_PLACE, avg, I_ONE, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, vol, I_ONE, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr) !! \todo calculate this in some init routine
      avg = avg / vol

      cgl => curl%first
      do while (associated(cgl))
         cgl%cg%q(iv)%arr(:, :, :) = cgl%cg%q(iv)%arr(:, :, :) - avg
         cgl => cgl%nxt
      enddo

   end subroutine subtract_average

!!$ ============================================================================
!>
!! \brief Prolong solution data at level (lev-1) to faces at level lev
!!
!! \details It looks that parallel prolongation order 0 is the best from range -2 .. 2. Quite good is also order 2.(tests were performed with the maclaurin problem at r4051).
!! There were no experiments with higher order prolongation, but it seems there is not much room for any improvements.
!! The tests were performed on uniform grid, where the interpolation affects only convergence factors.
!! On a refinement step the interpolation type affects also sulution by influencing the way how fine and coarse grids are coupled.
!!
!! \todo move to cg_list_level or multigrid_gravity
!!
!! OPT: completely unoptimized
!<

   subroutine prolong_faces(fine, soln)

      use constants,         only: xdim, ydim, zdim, LO, HI, LONG, I_ONE, half
      use dataio_pub,        only: die, warn
      use domain,            only: dom, is_multicg
      use cg_list_lev,       only: cg_list_level
      use grid_cont,         only: pr_segment
      use mpi,               only: MPI_DOUBLE_PRECISION
      use mpisetup,          only: proc, comm, ierr, req, status, master
      use multigridhelpers,  only: check_dirty
      use multigridmpifuncs, only: mpi_multigrid_bnd
      use multigridvars,     only: ord_prolong_face_norm, ord_prolong_face_par, base, extbnd_antimirror, is_external, need_general_pf

      implicit none

      type(cg_list_level), pointer, intent(in) :: fine !< level for which approximate the solution
      integer(kind=4), intent(in) :: soln !< index of solution in cg%q(:) ! \todo change the name

      integer                       :: i, j, k, d, lh, g, g1, gc, ib, jb, ibh, jbh, l
      type(cg_list_level), pointer           :: coarse
      integer, parameter            :: s_wdth  = 3           ! interpolation stencil width
      integer(kind=4), parameter    :: s_rng = (s_wdth-1)/2  ! stencil range around 0
      real, parameter, dimension(s_wdth) :: p0  = [ 0.,       1.,     0.     ] ! injection
      real, parameter, dimension(s_wdth) :: p1  = [ 0.,       3./4.,  1./4.  ] ! 1D linear prolongation stencil
      real, parameter, dimension(s_wdth) :: p2i = [ -1./8.,   1.,     1./8.  ] ! 1D integral cubic prolongation stencil
      real, parameter, dimension(s_wdth) :: p2d = [ -3./32., 30./32., 5./32. ] ! 1D direct cubic prolongation stencil
      real, dimension(-s_rng:s_rng)                    :: p
      real, dimension(-s_rng:s_rng, -s_rng:s_rng, 2, 2):: pp   ! 2D prolongation stencil
      real                          :: pp_norm

      integer, dimension(xdim:zdim) :: ii
      integer, dimension(xdim:zdim) :: off1
      integer(kind=4) :: nr
      integer(kind=8), dimension(:,:), pointer :: cse ! shortcut for coarse segment
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: se
      type c_bnd
         real, dimension(:, :), pointer :: bnd
      end type c_bnd
      type(c_bnd), dimension(xdim:zdim, LO:HI) :: p_bnd
      real :: opfn1, opfn3
      integer(kind=4) :: b_rng
      integer, dimension(xdim:zdim), parameter :: d1 = [ ydim, xdim, xdim ] , d2 = [ zdim, zdim, ydim ]
      type(pr_segment), pointer :: pseg

      ! coefficients for intepolation perpendicular to the face can be expressed as
      ! 1/2*| 1 0 ... ] + a2*| 1 -1 0 ... ] + a4*| 1 -2 1 0 ... ] + a6*| 1 -3 3 -1 0 ... ] + a8*| 1 -4 6 -4 1 0 ... ] + ...
      ! 1./2.   * ( 0   0   1   1   0 0 ) ! average, direct linear (a2 =  0.,      ...)
      ! 1./16.  * ( 0  -1   9   9  -1 0 ) ! direct cubic           (a2 =  1./16.,  a4 = 0.,      ...)
      ! 1./256. * ( 3 -25 150 150 -25 3 ) ! direct quintic         (a2 = 19./256., a4 = 3./256., a6 = 0., ...)

      ! factor 2./8. (a2 = 1./8.) was found experimentally
      opfn1 = 1. + 2./8. ! an additional factor equal to 0.5 is applied later as normalization of  pp(:,:,:,:)
      opfn3 = -2./8.
      ! the optimal value of a2 seems to depend at least on ord_prolong_face_par
      ! for ord_prolong_face_par =  0 it is ~1.7 (best convergence, reduction form 8/9 cycles to 7 cycles in maclaurin test for norm_tol=1e-10 and nsmool=8)
      ! for ord_prolong_face_par =  2 it is ~2,2 (slightly worse convergence, 8 cycles instead of 9/10 cycles)
      ! for ord_prolong_face_par = -2 it is ~2.5 (convergence in 8 cycles instead of 10 cycles)
      ! for ord_prolong_face_par =  1 it is ~3.6 (worst convergence, 11 -> 9 cycles)

      if (associated(fine, base)) then
         call warn("[multigridbasefuncs:prolong_faces] Cannot prolong anything to base level")
         return
      endif
      if (is_multicg) call die("[multigridbasefuncs:prolong_faces] multicg not implemented")

      if (.not. associated(fine)) call die("[multigridbasefuncs:prolong_faces] fine == null()")
      coarse => fine%coarser
      if (.not. associated(coarse)) call die("[multigridbasefuncs:prolong_faces] coarse == null()")

      if (need_general_pf) then

         if (ord_prolong_face_par /= 0) then
            ord_prolong_face_par = 0
            if (master) call warn("[multigridbasefuncs:prolong_faces] only injection is supported for the current domain decomposition type.")
            ! The tests made with comm3d suggests that paralel interpolation only degrades convergence rate.
            ! Implement ord_prolong_face_par /= 0 if and only if it improves the coupling between fine and coarse solutions
         endif

         do lh = LO, HI ! \todo convert cg_list_level%mg%bnd_[xyz] to an array and make obsolete the following pointer assignments
            p_bnd(xdim, lh)%bnd => fine%first%cg%mg%bnd_x(:, :, lh)
            p_bnd(ydim, lh)%bnd => fine%first%cg%mg%bnd_y(:, :, lh)
            p_bnd(zdim, lh)%bnd => fine%first%cg%mg%bnd_z(:, :, lh)
            do d = xdim, zdim
               p_bnd(d, lh)%bnd(:, :) = 0. ! \todo mark dirty somehow, eg use a "magic" small value, like 1.234567890123456789e-123
            enddo
         enddo

         nr = 0

         ! Gather coarse data for internal boundaries from others
         do d = xdim, zdim
            do lh = LO, HI
               if (allocated(fine%first%cg%mg%pfc_tgt(d, lh)%seg)) then
                  do g = 1, ubound(fine%first%cg%mg%pfc_tgt(d, lh)%seg(:), dim=1)
                     if (fine%first%cg%mg%pfc_tgt(d, lh)%seg(g)%proc /= proc) then
                        nr = nr + I_ONE
                        call MPI_Irecv(fine%first%cg%mg%pfc_tgt(d, lh)%seg(g)%buf(1, 1, 1), size(fine%first%cg%mg%pfc_tgt(d, lh)%seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, &
                             &         fine%first%cg%mg%pfc_tgt(d, lh)%seg(g)%proc, HI*d+lh, comm, req(nr), ierr)
                     endif
                  enddo
               endif
            enddo
         enddo

         ! Send own coarse data to others
         off1(:) = int(mod(fine%first%cg%off(:), 2_LONG), kind=4)
         do d = xdim, zdim
            do lh = LO, HI
               if (allocated(coarse%first%cg%mg%pff_tgt(d, lh)%seg)) then
                  do g = 1, ubound(coarse%first%cg%mg%pff_tgt(d, lh)%seg(:), dim=1)

                     pseg => coarse%first%cg%mg%pff_tgt(d, lh)%seg(g)
                     cse => pseg%se

                     if (pseg%proc /= proc) then
                        nr = nr + I_ONE
                        se(:,:) = cse(:,:)
                        pseg%buf(:, :, :) = 0. ! this can be avoided by extracting first assignment from the loop
                        do l = 1, ubound(pseg%f_lay(:), dim=1)
                           se(d,:) = pseg%f_lay(l)%layer
                           pseg%buf(:, :, :) = pseg%buf(:, :, :) + pseg%f_lay(l)%coeff * &
                                coarse%first%cg%q(soln)%arr(se(xdim, LO):se(xdim, HI), se(ydim, LO):se(ydim, HI), se(zdim, LO):se(zdim, HI))
                        enddo
                        call MPI_Isend(pseg%buf(1, 1, 1), size(pseg%buf(:, :, :)), MPI_DOUBLE_PRECISION, pseg%proc, HI*d+lh, comm, req(nr), ierr)
                     endif
                  enddo
               endif
            enddo
         enddo

         ! Make a copy of local coarse data to buffer (\todo: can be merged with interpolation step)
         do d = xdim, zdim
            do lh = LO, HI
               if (allocated(fine%first%cg%mg%pfc_tgt(d, lh)%seg)) then
                  do g = 1, ubound(fine%first%cg%mg%pfc_tgt(d, lh)%seg(:), dim=1)
                     if (fine%first%cg%mg%pfc_tgt(d, lh)%seg(g)%proc == proc) then
                        nullify(cse)
                        gc = 0
                        do g1 = 1, ubound(coarse%first%cg%mg%pff_tgt(d, lh)%seg(:), dim=1) !> \todo should be set up in mpi_multigrid_prep_grav
                           if (coarse%first%cg%mg%pff_tgt(d, lh)%seg(g1)%proc == proc) then
                              if (.not. associated(cse)) then
                                 cse => coarse%first%cg%mg%pff_tgt(d, lh)%seg(g1)%se
                                 gc = g1
                              else
                                 call die("[multigridbasefuncs:prolong_faces] multiple local coarse grid targets")
                              endif
                           endif
                        enddo
                        if (.not. associated(cse)) call die("[multigridbasefuncs:prolong_faces] missing coarse grid target")

                        se(:,:) = cse(:,:)
                        pseg => fine%first%cg%mg%pfc_tgt(d, lh)%seg(g)
                        pseg%buf(:,:,:) = 0. ! this can be avoided by extracting first assignment from the loop
                        do l = 1, ubound(coarse%first%cg%mg%pff_tgt(d, lh)%seg(gc)%f_lay(:), dim=1)
                           se(d,:) = coarse%first%cg%mg%pff_tgt(d, lh)%seg(gc)%f_lay(l)%layer
                           pseg%buf(:,:,:) =  pseg%buf(:,:,:) + coarse%first%cg%mg%pff_tgt(d, lh)%seg(gc)%f_lay(l)%coeff * &
                                coarse%first%cg%q(soln)%arr(se(xdim, LO):se(xdim, HI), se(ydim, LO):se(ydim, HI), se(zdim, LO):se(zdim, HI))
                        enddo
                     endif
                  enddo
               endif
            enddo
         enddo

         if (nr>0) call MPI_Waitall(nr, req(:nr), status(:,:nr), ierr)

         ! Interpolate content of buffers to boundary face-layes fine%first%cg%mg%bnd_[xyz](:, :, :)
         do d = xdim, zdim
            if (dom%has_dir(d)) then
               do lh = LO, HI
                  ! make external homogenous Dirichlet boundaries, where applicable, interpolate (inject) otherwise.
                  ! \todo for homogenous neumann boundary lay2 should be .false. and pc_bnd(1, :, :) should contain layer of the coarse data next to the external boundary
                  if (is_external(d, lh)) then
                     p_bnd(d, lh)%bnd(:,:) = 0. ! BEWARE homogenous Dirichlet by default
                  else
                     p_bnd(d, lh)%bnd(:,:) = 0
                     if (allocated(fine%first%cg%mg%pfc_tgt(d, lh)%seg)) then
                        do g = 1, ubound(fine%first%cg%mg%pfc_tgt(d, lh)%seg(:), dim=1)

                           ii(d) = 1
                           do i = 1, ubound(fine%first%cg%mg%pfc_tgt(d, lh)%seg(g)%buf, dim=d1(d))
                              ii(d1(d)) = i
                              do j = 1, ubound(fine%first%cg%mg%pfc_tgt(d, lh)%seg(g)%buf, dim=d2(d))
                                 ii(d2(d)) = j

! ? + off (0 or 1)
                                 ib = 2*i - 1 + int(fine%first%cg%mg%pfc_tgt(d, lh)%seg(g)%se(d1(d), LO), kind=4) - int(mod(fine%first%cg%off(d1(d)), 2_LONG), 4)
                                 jb = 2*j - 1 + int(fine%first%cg%mg%pfc_tgt(d, lh)%seg(g)%se(d2(d), LO), kind=4) - int(mod(fine%first%cg%off(d2(d)), 2_LONG), 4)

                                 ibh = ib
                                 if (dom%has_dir(d1(d)) .and. ubound(p_bnd(d, lh)%bnd(:,:), dim=1) > ib) ibh = ib + 1
                                 if (dom%has_dir(d1(d)) .and. lbound(p_bnd(d, lh)%bnd(:,:), dim=1) > ib) ib  = ib + 1
                                 jbh = jb
                                 if (dom%has_dir(d2(d)) .and. ubound(p_bnd(d, lh)%bnd(:,:), dim=2) > jb) jbh = jb + 1
                                 if (dom%has_dir(d2(d)) .and. lbound(p_bnd(d, lh)%bnd(:,:), dim=2) > jb) jb  = jb + 1

                                 p_bnd(d, lh)%bnd(ib:ibh, jb:jbh) = p_bnd(d, lh)%bnd(ib:ibh, jb:jbh) + fine%first%cg%mg%pfc_tgt(d, lh)%seg(g)%buf(ii(xdim), ii(ydim), ii(zdim))
                              enddo
                           enddo
                        enddo
                     endif
                  endif
               enddo
            endif
         enddo

      else
         select case (ord_prolong_face_par)
            case (0)
               p(:) = p0(:)
            case (1,-1)
               p(:) = p1(:)
            case (2)
               p(:) = p2i(:)
            case (-2)
               p(:) = p2d(:)
            case default
               p(:) = p0(:)
         end select

         do i = -s_rng, s_rng
            pp(i,:,1,1) = half*p( i)*p(:)       ! half because of face averaging
            pp(i,:,1,2) = half*p( i)*p(1:-1:-1) ! or use matmul()
            pp(i,:,2,1) = half*p(-i)*p(:)
            pp(i,:,2,2) = half*p(-i)*p(1:-1:-1)
         enddo

         if (ord_prolong_face_norm > 1) ord_prolong_face_norm = 1
         b_rng = s_rng
         if (ord_prolong_face_norm > 0) b_rng = max(b_rng, int(ord_prolong_face_norm+1, kind=4))
         call mpi_multigrid_bnd(coarse, soln, b_rng, extbnd_antimirror, corners=(ord_prolong_face_par/=0)) !> \deprecated BEWARE for higher prolongation order more guardcell are required
         call check_dirty(coarse, soln, "prolong_faces", s_rng)

         if (ord_prolong_face_norm > 0) then
            if (dom%has_dir(xdim)) then
               pp_norm = 2.*sum(pp(-dom%D_y:dom%D_y, -dom%D_z:dom%D_z, 1, 1)) ! normalization is required for ord_prolong_face_par == 1 and -2
               do j = coarse%first%cg%js, coarse%first%cg%je
                  do k = coarse%first%cg%ks, coarse%first%cg%ke
                     fine%first%cg%mg%bnd_x(-fine%first%cg%js+2*j,        -fine%first%cg%ks+2*k,        LO)=sum(pp(-dom%D_y:dom%D_y, -dom%D_z:dom%D_z, 1, 1) * ( &
                          opfn1*(coarse%first%cg%q(soln)%arr(coarse%first%cg%is,  j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(coarse%first%cg%is-1,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z)) + &
                          opfn3*(coarse%first%cg%q(soln)%arr(coarse%first%cg%is+1,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(coarse%first%cg%is-2,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z)))) / pp_norm
                     fine%first%cg%mg%bnd_x(-fine%first%cg%js+2*j+dom%D_y,-fine%first%cg%ks+2*k,        LO)=sum(pp(-dom%D_y:dom%D_y, -dom%D_z:dom%D_z, 2, 1) * ( &
                          opfn1*(coarse%first%cg%q(soln)%arr(coarse%first%cg%is,  j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(coarse%first%cg%is-1,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z)) + &
                          opfn3*(coarse%first%cg%q(soln)%arr(coarse%first%cg%is+1,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(coarse%first%cg%is-2,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z)))) / pp_norm
                     fine%first%cg%mg%bnd_x(-fine%first%cg%js+2*j,        -fine%first%cg%ks+2*k+dom%D_z,LO)=sum(pp(-dom%D_y:dom%D_y, -dom%D_z:dom%D_z, 1, 2) * ( &
                          opfn1*(coarse%first%cg%q(soln)%arr(coarse%first%cg%is  ,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(coarse%first%cg%is-1,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z)) + &
                          opfn3*(coarse%first%cg%q(soln)%arr(coarse%first%cg%is+1,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(coarse%first%cg%is-2,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z)))) / pp_norm
                     fine%first%cg%mg%bnd_x(-fine%first%cg%js+2*j+dom%D_y,-fine%first%cg%ks+2*k+dom%D_z,LO)=sum(pp(-dom%D_y:dom%D_y, -dom%D_z:dom%D_z, 2, 2) * ( &
                          opfn1*(coarse%first%cg%q(soln)%arr(coarse%first%cg%is,  j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(coarse%first%cg%is-1,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z)) + &
                          opfn3*(coarse%first%cg%q(soln)%arr(coarse%first%cg%is+1,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(coarse%first%cg%is-2,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z)))) / pp_norm
                     fine%first%cg%mg%bnd_x(-fine%first%cg%js+2*j,        -fine%first%cg%ks+2*k,        HI)=sum(pp(-dom%D_y:dom%D_y, -dom%D_z:dom%D_z, 1, 1) * ( &
                          opfn1*(coarse%first%cg%q(soln)%arr(coarse%first%cg%ie,  j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(coarse%first%cg%ie+1,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z)) + &
                          opfn3*(coarse%first%cg%q(soln)%arr(coarse%first%cg%ie-1,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(coarse%first%cg%ie+2,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z)))) / pp_norm
                     fine%first%cg%mg%bnd_x(-fine%first%cg%js+2*j+dom%D_y,-fine%first%cg%ks+2*k,        HI)=sum(pp(-dom%D_y:dom%D_y, -dom%D_z:dom%D_z, 2, 1) * ( &
                          opfn1*(coarse%first%cg%q(soln)%arr(coarse%first%cg%ie,  j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(coarse%first%cg%ie+1,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z)) + &
                          opfn3*(coarse%first%cg%q(soln)%arr(coarse%first%cg%ie-1,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(coarse%first%cg%ie+2,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z)))) / pp_norm
                     fine%first%cg%mg%bnd_x(-fine%first%cg%js+2*j,        -fine%first%cg%ks+2*k+dom%D_z,HI)=sum(pp(-dom%D_y:dom%D_y, -dom%D_z:dom%D_z, 1, 2) * ( &
                          opfn1*(coarse%first%cg%q(soln)%arr(coarse%first%cg%ie,  j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(coarse%first%cg%ie+1,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z)) + &
                          opfn3*(coarse%first%cg%q(soln)%arr(coarse%first%cg%ie-1,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(coarse%first%cg%ie+2,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z)))) / pp_norm
                     fine%first%cg%mg%bnd_x(-fine%first%cg%js+2*j+dom%D_y,-fine%first%cg%ks+2*k+dom%D_z,HI)=sum(pp(-dom%D_y:dom%D_y, -dom%D_z:dom%D_z, 2, 2) * ( &
                          opfn1*(coarse%first%cg%q(soln)%arr(coarse%first%cg%ie,  j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(coarse%first%cg%ie+1,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z)) + &
                          opfn3*(coarse%first%cg%q(soln)%arr(coarse%first%cg%ie-1,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(coarse%first%cg%ie+2,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z)))) / pp_norm
                  enddo
               enddo
            endif

            if (dom%has_dir(ydim)) then
               pp_norm = 2.*sum(pp(-dom%D_x:dom%D_x, -dom%D_z:dom%D_z, 1, 1))
               do i = coarse%first%cg%is, coarse%first%cg%ie
                  do k = coarse%first%cg%ks, coarse%first%cg%ke
                     fine%first%cg%mg%bnd_y(-fine%first%cg%is+2*i,        -fine%first%cg%ks+2*k,        LO)=sum(pp(-dom%D_x:dom%D_x, -dom%D_z:dom%D_z, 1, 1) * ( &
                          opfn1*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%js,  k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%js-1,k-dom%D_z:k+dom%D_z)) + &
                          opfn3*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%js+1,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%js-2,k-dom%D_z:k+dom%D_z)))) / pp_norm
                     fine%first%cg%mg%bnd_y(-fine%first%cg%is+2*i+dom%D_x,-fine%first%cg%ks+2*k,        LO)=sum(pp(-dom%D_x:dom%D_x, -dom%D_z:dom%D_z, 2, 1) * ( &
                          opfn1*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%js,  k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%js-1,k-dom%D_z:k+dom%D_z)) + &
                          opfn3*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%js+1,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%js-2,k-dom%D_z:k+dom%D_z)))) / pp_norm
                     fine%first%cg%mg%bnd_y(-fine%first%cg%is+2*i,        -fine%first%cg%ks+2*k+dom%D_z,LO)=sum(pp(-dom%D_x:dom%D_x, -dom%D_z:dom%D_z, 1, 2) * ( &
                          opfn1*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%js,  k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%js-1,k-dom%D_z:k+dom%D_z)) + &
                          opfn3*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%js+1,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%js-2,k-dom%D_z:k+dom%D_z)))) / pp_norm
                     fine%first%cg%mg%bnd_y(-fine%first%cg%is+2*i+dom%D_x,-fine%first%cg%ks+2*k+dom%D_z,LO)=sum(pp(-dom%D_x:dom%D_x, -dom%D_z:dom%D_z, 2, 2) * ( &
                          opfn1*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%js,  k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%js-1,k-dom%D_z:k+dom%D_z)) + &
                          opfn3*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%js+1,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%js-2,k-dom%D_z:k+dom%D_z)))) / pp_norm
                     fine%first%cg%mg%bnd_y(-fine%first%cg%is+2*i,        -fine%first%cg%ks+2*k,        HI)=sum(pp(-dom%D_x:dom%D_x, -dom%D_z:dom%D_z, 1, 1) * ( &
                          opfn1*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%je,  k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%je+1,k-dom%D_z:k+dom%D_z)) + &
                          opfn3*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%je-1,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%je+2,k-dom%D_z:k+dom%D_z)))) / pp_norm
                     fine%first%cg%mg%bnd_y(-fine%first%cg%is+2*i+dom%D_x,-fine%first%cg%ks+2*k,        HI)=sum(pp(-dom%D_x:dom%D_x, -dom%D_z:dom%D_z, 2, 1) * ( &
                          opfn1*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%je,  k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%je+1,k-dom%D_z:k+dom%D_z)) + &
                          opfn3*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%je-1,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%je+2,k-dom%D_z:k+dom%D_z)))) / pp_norm
                     fine%first%cg%mg%bnd_y(-fine%first%cg%is+2*i,        -fine%first%cg%ks+2*k+dom%D_z,HI)=sum(pp(-dom%D_x:dom%D_x, -dom%D_z:dom%D_z, 1, 2) * ( &
                          opfn1*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%je,  k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%je+1,k-dom%D_z:k+dom%D_z)) + &
                          opfn3*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%je-1,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%je+2,k-dom%D_z:k+dom%D_z)))) / pp_norm
                     fine%first%cg%mg%bnd_y(-fine%first%cg%is+2*i+dom%D_x,-fine%first%cg%ks+2*k+dom%D_z,HI)=sum(pp(-dom%D_x:dom%D_x, -dom%D_z:dom%D_z, 2, 2) * ( &
                          opfn1*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%je,  k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%je+1,k-dom%D_z:k+dom%D_z)) + &
                          opfn3*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%je-1,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%je+2,k-dom%D_z:k+dom%D_z)))) / pp_norm
                  enddo
               enddo
            endif

            if (dom%has_dir(zdim)) then
               pp_norm = 2.*sum(pp(-dom%D_x:dom%D_x, -dom%D_y:dom%D_y, 1, 1))
               do i = coarse%first%cg%is, coarse%first%cg%ie
                  do j = coarse%first%cg%js, coarse%first%cg%je
                     fine%first%cg%mg%bnd_z(-fine%first%cg%is+2*i,        -fine%first%cg%js+2*j,        LO)=sum(pp(-dom%D_x:dom%D_x, -dom%D_y:dom%D_y, 1, 1) * ( &
                          opfn1*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ks) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ks-1)) + &
                          opfn3*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ks+1) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ks-2)))) / pp_norm
                     fine%first%cg%mg%bnd_z(-fine%first%cg%is+2*i+dom%D_x,-fine%first%cg%js+2*j,        LO)=sum(pp(-dom%D_x:dom%D_x, -dom%D_y:dom%D_y, 2, 1) * ( &
                          opfn1*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ks) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ks-1)) + &
                          opfn3*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ks+1) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ks-2)))) / pp_norm
                     fine%first%cg%mg%bnd_z(-fine%first%cg%is+2*i,        -fine%first%cg%js+2*j+dom%D_y,LO)=sum(pp(-dom%D_x:dom%D_x, -dom%D_y:dom%D_y, 1, 2) * ( &
                          opfn1*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ks) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ks-1)) + &
                          opfn3*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ks+1) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ks-2)))) / pp_norm
                     fine%first%cg%mg%bnd_z(-fine%first%cg%is+2*i+dom%D_x,-fine%first%cg%js+2*j+dom%D_y,LO)=sum(pp(-dom%D_x:dom%D_x, -dom%D_y:dom%D_y, 2, 2) * ( &
                          opfn1*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ks) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ks-1)) + &
                          opfn3*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ks+1) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ks-2)))) / pp_norm
                     fine%first%cg%mg%bnd_z(-fine%first%cg%is+2*i,        -fine%first%cg%js+2*j,        HI)=sum(pp(-dom%D_x:dom%D_x, -dom%D_y:dom%D_y, 1, 1) * ( &
                          opfn1*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ke) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ke+1)) + &
                          opfn3*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ke-1) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ke+2)))) / pp_norm
                     fine%first%cg%mg%bnd_z(-fine%first%cg%is+2*i+dom%D_x,-fine%first%cg%js+2*j,        HI)=sum(pp(-dom%D_x:dom%D_x, -dom%D_y:dom%D_y, 2, 1) * ( &
                          opfn1*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ke) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ke+1)) + &
                          opfn3*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ke-1) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ke+2)))) / pp_norm
                     fine%first%cg%mg%bnd_z(-fine%first%cg%is+2*i,        -fine%first%cg%js+2*j+dom%D_y,HI)=sum(pp(-dom%D_x:dom%D_x, -dom%D_y:dom%D_y, 1, 2) * ( &
                          opfn1*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ke) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ke+1)) + &
                          opfn3*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ke-1) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ke+2)))) / pp_norm
                     fine%first%cg%mg%bnd_z(-fine%first%cg%is+2*i+dom%D_x,-fine%first%cg%js+2*j+dom%D_y,HI)=sum(pp(-dom%D_x:dom%D_x, -dom%D_y:dom%D_y, 2, 2) * ( &
                          opfn1*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ke) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ke+1)) + &
                          opfn3*(coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ke-1) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ke+2)))) / pp_norm
                  enddo
               enddo
            endif

         else

            if (dom%has_dir(xdim)) then
               pp_norm = 2.*sum(pp(-dom%D_y:dom%D_y, -dom%D_z:dom%D_z, 1, 1)) ! normalization is required for ord_prolong_face_par == 1 and -2
               do j = coarse%first%cg%js, coarse%first%cg%je
                  do k = coarse%first%cg%ks, coarse%first%cg%ke
                     fine%first%cg%mg%bnd_x(-fine%first%cg%js+2*j,        -fine%first%cg%ks+2*k,        LO)=sum(pp(-dom%D_y:dom%D_y, -dom%D_z:dom%D_z, 1, 1) * &
                        & (coarse%first%cg%q(soln)%arr(coarse%first%cg%is,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(coarse%first%cg%is-1,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z))) / pp_norm
                     fine%first%cg%mg%bnd_x(-fine%first%cg%js+2*j+dom%D_y,-fine%first%cg%ks+2*k,        LO)=sum(pp(-dom%D_y:dom%D_y, -dom%D_z:dom%D_z, 2, 1) * &
                        & (coarse%first%cg%q(soln)%arr(coarse%first%cg%is,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(coarse%first%cg%is-1,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z))) / pp_norm
                     fine%first%cg%mg%bnd_x(-fine%first%cg%js+2*j,        -fine%first%cg%ks+2*k+dom%D_z,LO)=sum(pp(-dom%D_y:dom%D_y, -dom%D_z:dom%D_z, 1, 2) * &
                        & (coarse%first%cg%q(soln)%arr(coarse%first%cg%is,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(coarse%first%cg%is-1,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z))) / pp_norm
                     fine%first%cg%mg%bnd_x(-fine%first%cg%js+2*j+dom%D_y,-fine%first%cg%ks+2*k+dom%D_z,LO)=sum(pp(-dom%D_y:dom%D_y, -dom%D_z:dom%D_z, 2, 2) * &
                        & (coarse%first%cg%q(soln)%arr(coarse%first%cg%is,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(coarse%first%cg%is-1,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z))) / pp_norm
                     fine%first%cg%mg%bnd_x(-fine%first%cg%js+2*j,        -fine%first%cg%ks+2*k,        HI)=sum(pp(-dom%D_y:dom%D_y, -dom%D_z:dom%D_z, 1, 1) * &
                        & (coarse%first%cg%q(soln)%arr(coarse%first%cg%ie,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(coarse%first%cg%ie+1,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z))) / pp_norm
                     fine%first%cg%mg%bnd_x(-fine%first%cg%js+2*j+dom%D_y,-fine%first%cg%ks+2*k,        HI)=sum(pp(-dom%D_y:dom%D_y, -dom%D_z:dom%D_z, 2, 1) * &
                        & (coarse%first%cg%q(soln)%arr(coarse%first%cg%ie,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(coarse%first%cg%ie+1,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z))) / pp_norm
                     fine%first%cg%mg%bnd_x(-fine%first%cg%js+2*j,        -fine%first%cg%ks+2*k+dom%D_z,HI)=sum(pp(-dom%D_y:dom%D_y, -dom%D_z:dom%D_z, 1, 2) * &
                        & (coarse%first%cg%q(soln)%arr(coarse%first%cg%ie,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(coarse%first%cg%ie+1,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z))) / pp_norm
                     fine%first%cg%mg%bnd_x(-fine%first%cg%js+2*j+dom%D_y,-fine%first%cg%ks+2*k+dom%D_z,HI)=sum(pp(-dom%D_y:dom%D_y, -dom%D_z:dom%D_z, 2, 2) * &
                        & (coarse%first%cg%q(soln)%arr(coarse%first%cg%ie,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(coarse%first%cg%ie+1,j-dom%D_y:j+dom%D_y,k-dom%D_z:k+dom%D_z))) / pp_norm
                  enddo
               enddo
            endif

            if (dom%has_dir(ydim)) then
               pp_norm = 2.*sum(pp(-dom%D_x:dom%D_x, -dom%D_z:dom%D_z, 1, 1))
               do i = coarse%first%cg%is, coarse%first%cg%ie
                  do k = coarse%first%cg%ks, coarse%first%cg%ke
                     fine%first%cg%mg%bnd_y(-fine%first%cg%is+2*i,        -fine%first%cg%ks+2*k,        LO)=sum(pp(-dom%D_x:dom%D_x, -dom%D_z:dom%D_z, 1, 1) * &
                        & (coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%js,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%js-1,k-dom%D_z:k+dom%D_z))) / pp_norm
                     fine%first%cg%mg%bnd_y(-fine%first%cg%is+2*i+dom%D_x,-fine%first%cg%ks+2*k,        LO)=sum(pp(-dom%D_x:dom%D_x, -dom%D_z:dom%D_z, 2, 1) * &
                        & (coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%js,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%js-1,k-dom%D_z:k+dom%D_z))) / pp_norm
                     fine%first%cg%mg%bnd_y(-fine%first%cg%is+2*i,        -fine%first%cg%ks+2*k+dom%D_z,LO)=sum(pp(-dom%D_x:dom%D_x, -dom%D_z:dom%D_z, 1, 2) * &
                        & (coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%js,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%js-1,k-dom%D_z:k+dom%D_z))) / pp_norm
                     fine%first%cg%mg%bnd_y(-fine%first%cg%is+2*i+dom%D_x,-fine%first%cg%ks+2*k+dom%D_z,LO)=sum(pp(-dom%D_x:dom%D_x, -dom%D_z:dom%D_z, 2, 2) * &
                        & (coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%js,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%js-1,k-dom%D_z:k+dom%D_z))) / pp_norm
                     fine%first%cg%mg%bnd_y(-fine%first%cg%is+2*i,        -fine%first%cg%ks+2*k,        HI)=sum(pp(-dom%D_x:dom%D_x, -dom%D_z:dom%D_z, 1, 1) * &
                        & (coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%je,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%je+1,k-dom%D_z:k+dom%D_z))) / pp_norm
                     fine%first%cg%mg%bnd_y(-fine%first%cg%is+2*i+dom%D_x,-fine%first%cg%ks+2*k,        HI)=sum(pp(-dom%D_x:dom%D_x, -dom%D_z:dom%D_z, 2, 1) * &
                        & (coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%je,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%je+1,k-dom%D_z:k+dom%D_z))) / pp_norm
                     fine%first%cg%mg%bnd_y(-fine%first%cg%is+2*i,        -fine%first%cg%ks+2*k+dom%D_z,HI)=sum(pp(-dom%D_x:dom%D_x, -dom%D_z:dom%D_z, 1, 2) * &
                        & (coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%je,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%je+1,k-dom%D_z:k+dom%D_z))) / pp_norm
                     fine%first%cg%mg%bnd_y(-fine%first%cg%is+2*i+dom%D_x,-fine%first%cg%ks+2*k+dom%D_z,HI)=sum(pp(-dom%D_x:dom%D_x, -dom%D_z:dom%D_z, 2, 2) * &
                        & (coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%je,k-dom%D_z:k+dom%D_z) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,coarse%first%cg%je+1,k-dom%D_z:k+dom%D_z))) / pp_norm
                  enddo
               enddo
            endif

            if (dom%has_dir(zdim)) then
               pp_norm = 2.*sum(pp(-dom%D_x:dom%D_x, -dom%D_y:dom%D_y, 1, 1))
               do i = coarse%first%cg%is, coarse%first%cg%ie
                  do j = coarse%first%cg%js, coarse%first%cg%je
                     fine%first%cg%mg%bnd_z(-fine%first%cg%is+2*i,        -fine%first%cg%js+2*j,        LO)=sum(pp(-dom%D_x:dom%D_x, -dom%D_y:dom%D_y, 1, 1) * &
                        & (coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ks) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ks-1))) / pp_norm
                     fine%first%cg%mg%bnd_z(-fine%first%cg%is+2*i+dom%D_x,-fine%first%cg%js+2*j,        LO)=sum(pp(-dom%D_x:dom%D_x, -dom%D_y:dom%D_y, 2, 1) * &
                        & (coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ks) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ks-1))) / pp_norm
                     fine%first%cg%mg%bnd_z(-fine%first%cg%is+2*i,        -fine%first%cg%js+2*j+dom%D_y,LO)=sum(pp(-dom%D_x:dom%D_x, -dom%D_y:dom%D_y, 1, 2) * &
                        & (coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ks) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ks-1))) / pp_norm
                     fine%first%cg%mg%bnd_z(-fine%first%cg%is+2*i+dom%D_x,-fine%first%cg%js+2*j+dom%D_y,LO)=sum(pp(-dom%D_x:dom%D_x, -dom%D_y:dom%D_y, 2, 2) * &
                        & (coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ks) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ks-1))) / pp_norm
                     fine%first%cg%mg%bnd_z(-fine%first%cg%is+2*i,        -fine%first%cg%js+2*j,        HI)=sum(pp(-dom%D_x:dom%D_x, -dom%D_y:dom%D_y, 1, 1) * &
                        & (coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ke) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ke+1))) / pp_norm
                     fine%first%cg%mg%bnd_z(-fine%first%cg%is+2*i+dom%D_x,-fine%first%cg%js+2*j,        HI)=sum(pp(-dom%D_x:dom%D_x, -dom%D_y:dom%D_y, 2, 1) * &
                        & (coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ke) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ke+1))) / pp_norm
                     fine%first%cg%mg%bnd_z(-fine%first%cg%is+2*i,        -fine%first%cg%js+2*j+dom%D_y,HI)=sum(pp(-dom%D_x:dom%D_x, -dom%D_y:dom%D_y, 1, 2) * &
                        & (coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ke) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ke+1))) / pp_norm
                     fine%first%cg%mg%bnd_z(-fine%first%cg%is+2*i+dom%D_x,-fine%first%cg%js+2*j+dom%D_y,HI)=sum(pp(-dom%D_x:dom%D_x, -dom%D_y:dom%D_y, 2, 2) * &
                        & (coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ke) + coarse%first%cg%q(soln)%arr(i-dom%D_x:i+dom%D_x,j-dom%D_y:j+dom%D_y,coarse%first%cg%ke+1))) / pp_norm
                  enddo
               enddo
            endif
         endif
      endif

   end subroutine prolong_faces

end module multigridbasefuncs
