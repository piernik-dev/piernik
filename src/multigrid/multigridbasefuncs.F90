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

   public :: prolong_level, prolong_faces, zero_boundaries

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

      use cg_list_lev,           only: cg_list_level
      use constants,             only: O_INJ
      use dataio_pub,            only: die
      use multigridhelpers,      only: dirty_debug, check_dirty, dirtyH
      use multigridvars,         only: roof, ord_prolong, is_mg_uneven

      implicit none

      type(cg_list_level), pointer, intent(in) :: coarse !< level to prolong from
      integer,                      intent(in) :: iv     !< variable to be prolonged

      type(cg_list_level), pointer :: fine

      if (associated(coarse, roof)) return ! can't prolong finest level

      if (.not. associated(coarse)) call die("[multigridbasefuncs:prolong_level] coarse == null()")
      fine   => coarse%finer
      if (.not. associated(fine)) call die("[multigridbasefuncs:prolong_level] fine == null()")

      !! \warning: need set_dirty_level routine
      if (dirty_debug) fine%first%cg%q(iv)%arr(:, :, :) = dirtyH

      call check_dirty(coarse, iv, "prolong-")

      if (ord_prolong == O_INJ .or. is_mg_uneven) then
         call coarse%prolong0_q_1var(iv)
      else
         call prolong_level_hord(coarse, iv) ! experimental part
      endif

      call check_dirty(fine, iv, "prolong+")

   end subroutine prolong_level

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

      use constants,         only: xdim, ydim, zdim, LO, HI, LONG, I_ONE, half, O_INJ, O_LIN, O_D2, O_I2
      use dataio_pub,        only: die, warn
      use domain,            only: dom, is_multicg
      use cg_list_lev,       only: cg_list_level
      use grid_cont,         only: pr_segment
      use mpi,               only: MPI_DOUBLE_PRECISION
      use mpisetup,          only: comm, ierr, req, status, master
      use multigridhelpers,  only: check_dirty
      use multigridmpifuncs, only: mpi_multigrid_bnd
      use multigridvars,     only: ord_prolong_face_norm, ord_prolong_face_par, base, extbnd_antimirror, need_general_pf

      implicit none

      type(cg_list_level), pointer, intent(in) :: fine !< level for which approximate the solution
      integer,                      intent(in) :: soln !< index of solution in cg%q(:) ! \todo change the name

      integer                       :: i, j, k, d, lh, g, ib, jb, ibh, jbh, l
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
      if (is_multicg) call die("[multigridbasefuncs:prolong_faces] multicg not implemented") ! fine%first%cg%

      if (.not. associated(fine)) call die("[multigridbasefuncs:prolong_faces] fine == null()")
      coarse => fine%coarser
      if (.not. associated(coarse)) call die("[multigridbasefuncs:prolong_faces] coarse == null()")

      if (need_general_pf) then

         if (ord_prolong_face_par /= O_INJ) then
            ord_prolong_face_par = O_INJ
            if (master) call warn("[multigridbasefuncs:prolong_faces] only injection is supported for the current domain decomposition type.")
            ! The tests made with comm3d suggests that paralel interpolation only degrades convergence rate.
            ! Implement ord_prolong_face_par /= O_INJ if and only if it improves the coupling between fine and coarse solutions
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
                  do g = lbound(fine%first%cg%mg%pfc_tgt(d, lh)%seg(:), dim=1), ubound(fine%first%cg%mg%pfc_tgt(d, lh)%seg(:), dim=1)
                     nr = nr + I_ONE
                     call MPI_Irecv(fine%first%cg%mg%pfc_tgt(d, lh)%seg(g)%buf(1, 1, 1), size(fine%first%cg%mg%pfc_tgt(d, lh)%seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, &
                          &         fine%first%cg%mg%pfc_tgt(d, lh)%seg(g)%proc, HI*d+lh, comm, req(nr), ierr)
                  enddo
               endif
            enddo
         enddo

         ! Send own coarse data to others
         off1(:) = int(mod(fine%first%cg%off(:), 2_LONG), kind=4)
         do d = xdim, zdim
            do lh = LO, HI
               if (associated(coarse%first)) then !!!
                  if (allocated(coarse%first%cg%mg%pff_tgt(d, lh)%seg)) then
                     do g = lbound(coarse%first%cg%mg%pff_tgt(d, lh)%seg(:), dim=1), ubound(coarse%first%cg%mg%pff_tgt(d, lh)%seg(:), dim=1)

                        pseg => coarse%first%cg%mg%pff_tgt(d, lh)%seg(g)
                        cse => pseg%se

                        nr = nr + I_ONE
                        se(:,:) = cse(:,:)
                        pseg%buf(:, :, :) = 0. ! this can be avoided by extracting first assignment from the loop
                        do l = lbound(pseg%f_lay(:), dim=1), ubound(pseg%f_lay(:), dim=1)
                           se(d,:) = pseg%f_lay(l)%layer
                           pseg%buf(:, :, :) = pseg%buf(:, :, :) + pseg%f_lay(l)%coeff * coarse%first%cg%q(soln)%span(se)
                        enddo
                        call MPI_Isend(pseg%buf(1, 1, 1), size(pseg%buf(:, :, :)), MPI_DOUBLE_PRECISION, pseg%proc, HI*d+lh, comm, req(nr), ierr)
                     enddo
                  endif
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
                  if (fine%first%cg%ext_bnd(d, lh)) then
                     p_bnd(d, lh)%bnd(:,:) = 0. ! BEWARE homogenous Dirichlet by default
                  else
                     p_bnd(d, lh)%bnd(:,:) = 0
                     if (allocated(fine%first%cg%mg%pfc_tgt(d, lh)%seg)) then
                        do g = lbound(fine%first%cg%mg%pfc_tgt(d, lh)%seg(:), dim=1), ubound(fine%first%cg%mg%pfc_tgt(d, lh)%seg(:), dim=1)

                           ii(d) = 1
                           do i = lbound(fine%first%cg%mg%pfc_tgt(d, lh)%seg(g)%buf, dim=d1(d)), ubound(fine%first%cg%mg%pfc_tgt(d, lh)%seg(g)%buf, dim=d1(d))
                              ii(d1(d)) = i
                              do j = lbound(fine%first%cg%mg%pfc_tgt(d, lh)%seg(g)%buf, dim=d2(d)), ubound(fine%first%cg%mg%pfc_tgt(d, lh)%seg(g)%buf, dim=d2(d))
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
            case (O_INJ)
               p(:) = p0(:)
            case (O_LIN)
               p(:) = p1(:)
            case (O_I2)
               p(:) = p2i(:)
            case (O_D2)
               p(:) = p2d(:)
            case default
               call die("[multigridbasefuncs:prolong_faces] invalid ord_prolong_face_par")
         end select

         do i = -s_rng, s_rng
            pp(i,:,1,1) = half*p( i)*p(:)       ! half because of face averaging
            pp(i,:,1,2) = half*p( i)*p(1:-1:-1) ! or use matmul()
            pp(i,:,2,1) = half*p(-i)*p(:)
            pp(i,:,2,2) = half*p(-i)*p(1:-1:-1)
         enddo

         if (ord_prolong_face_norm > O_LIN) ord_prolong_face_norm = O_LIN
         if (ord_prolong_face_norm < O_INJ) ord_prolong_face_norm = O_INJ
         b_rng = s_rng
         if (ord_prolong_face_norm > O_INJ) b_rng = max(b_rng, int(ord_prolong_face_norm+1, kind=4))
         call mpi_multigrid_bnd(coarse, soln, b_rng, extbnd_antimirror, corners=(ord_prolong_face_par/=0)) !> \deprecated BEWARE for higher prolongation order more guardcell are required
         call check_dirty(coarse, soln, "prolong_faces", s_rng)

         if (ord_prolong_face_norm /= O_INJ) then
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

!>
!! \brief high-order order prolongation interpolation
!!
!! \details
!! <table border="1" cellpadding="4" cellspacing="0">
!!   <tr><td> Cell-face prolongation stencils for fast convergence on uniform grid </td>
!!       <td> -1./12. </td><td> 7./12. </td><td> 7./12. </td><td> -1./12. </td><td> integral cubic </td></tr>
!!   <tr><td> Slightly slower convergence, less wide stencil  </td>
!!       <td>         </td><td> 1./2.  </td><td> 1./2.  </td><td>         </td><td> average; integral and direct linear </td></tr>
!! </table>
!!\n Prolongation of cell faces from cell centers are required for FFT local solver, red-black Gauss-Seidel relaxation don't use it.
!!
!!\n Cell-centered prolongation stencils, for odd fine cells, for even fine cells reverse the order.
!! <table border="1" cellpadding="4" cellspacing="0">
!!   <tr><td> 35./2048. </td><td> -252./2048. </td><td> 1890./2048. </td><td> 420./2048. </td><td> -45./2048. </td><td> direct quartic </td></tr>
!!   <tr><td>           </td><td>   -7./128.  </td><td>  105./128.  </td><td>  35./128.  </td><td>  -5./128.  </td><td> direct cubic </td></tr>
!!   <tr><td>           </td><td>   -3./32.   </td><td>   30./32.   </td><td>   5./32.   </td><td>            </td><td> direct quadratic </td></tr>
!!   <tr><td>           </td><td>             </td><td>    1.       </td><td>            </td><td>            </td><td> injection (0-th order), direct and integral approach </td></tr>
!!   <tr><td>           </td><td>             </td><td>    3./4.    </td><td>   1./4.    </td><td>            </td><td> linear, direct and integral approach </td></tr>
!!   <tr><td>           </td><td>   -1./8.    </td><td>    1.       </td><td>   1./8.    </td><td>            </td><td> integral quadratic </td></tr>
!!   <tr><td>           </td><td>   -5./64.   </td><td>   55./64.   </td><td>  17./64.   </td><td>  -3./64.   </td><td> integral cubic </td></tr>
!!   <tr><td>   3./128. </td><td>  -11./64.   </td><td>    1.       </td><td>  11./64.   </td><td>  -3./128.  </td><td> integral quartic </td></tr>
!! </table>
!!
!!\n General rule is that the second big positive coefficient should be assigned to closer neighbor of the coarse parent cell.
!!\n Thus a single coarse contributes to fine cells in the following way:
!! <table border="1" cellpadding="4" cellspacing="0">
!!   <tr><td> fine level   </td>
!!       <td> -3./128. </td><td> 3./128. </td><td> -11./64. </td><td>  11./64. </td><td> 1. </td><td> 1. </td>
!!       <td> 11./64. </td><td> -11./64. </td><td> 3./128.  </td><td> -3./128. </td><td> integral quartic coefficients </td></tr>
!!   <tr><td> coarse level </td>
!!       <td colspan="2">                </td><td colspan="2">                 </td><td colspan="2"> 1.  </td>
!!       <td colspan="2">                </td><td colspan="2">                 </td><td>                               </td></tr>
!! </table>
!!\n
!!\n The term "n-th order integral interpolation" here means that the prolonged values satisfy the following condition:
!!\n Integral over a cell of a n-th order polynomial fit to the nearest 5 points in each dimension on coarse level
!! is equal to the sum of similar integrals over fine cells covering the coarse cell.
!!\n
!!\n The term "n-th order direct interpolation" here means that the prolonged values are n-th order polynomial fit
!! to the nearest 5 points in each dimension on coarse level evaluated for fine cell centers.
!!\n
!!\n It seems that for 3D Cartesian grid with isolated boundaries and relaxation implemented in approximate_solution
!! direct quadratic and cubic interpolations give best norm reduction factors per V-cycle (maclaurin problem).
!!  For other boundary types, FFT implementation of approximate_solution, specific source distribution or
!!  some other multigrid scheme may give faster convergence rate.
!!\n
!!\n Estimated prolongation costs for integral quartic stencil:
!!\n  "gather" approach: loop over fine cells, each one collects weighted values from 5**3 coarse cells (125*n_fine multiplications
!!\n  "scatter" approach: loop over coarse cells, each one contributes weighted values to 10**3 fine cells (1000*n_coarse multiplications, roughly equal to cost of gather)
!!\n  "directionally split" approach: do the prolongation (either gather or scatter type) first in x direction (10*n_coarse multiplications -> 2*n_coarse intermediate cells
!!                                  result), then in y direction (10*2*n_coarse multiplications -> 4*n_coarse intermediate cells result), then in z direction
!!                                  (10*4*n_coarse multiplications -> 8*n_coarse = n_fine cells result). Looks like 70*n_coarse multiplications.
!!                                  Will require two additional arrays for intermediate results.
!!\n  "FFT" approach: do the convolution in Fourier space. Unfortunately it is not periodic box, so we cannot use power of 2 FFT sizes. No idea how fast or slow this can be.
!!\n
!!\n For AMR or nested grid low-order prolongation schemes (injection and linear interpolation at least) are known to produce artifacts
!! on fine-coarse boundaries. For uniform grid the simplest operators are probably the fastest and give best V-cycle convergence rates.
!!
!! \deprecated These routines assume simplest domain decomposition where fine grids cover exactly the same area as coarse grids
!!
!! \todo This will finally belong to cg_list_level type
!<

   subroutine prolong_level_hord(coarse, iv)

      use cg_list_lev,       only: cg_list_level
      use constants,         only: ndims, O_INJ, O_LIN, O_D2, O_D3, O_D4, O_I2, O_I3, O_I4
      use dataio_pub,        only: die, warn, msg
      use domain,            only: dom, is_multicg
      use grid_cont,         only: grid_container
      use mpisetup,          only: master
      use multigridmpifuncs, only: mpi_multigrid_bnd
      use multigridvars,     only: ord_prolong, extbnd_antimirror

      implicit none

      type(cg_list_level), pointer, intent(in) :: coarse !< level to prolong from
      integer,                      intent(in) :: iv     !< variable to be prolonged

      logical, save :: firstcall = .true.
      real :: P_2, P_1, P0, P1, P2
      type(grid_container), pointer :: cg_f, cg_c

      if (firstcall) then
         if (master) then
            write(msg,'(a,i3,a)')"[multigridbasefuncs:prolong_level_hord] prolongation order ",ord_prolong," is experimental"
            call warn(msg)
         endif
         firstcall = .false.
      endif
      if (abs(ord_prolong) > 2*dom%nb) call die("[multigridbasefuncs:prolong_level_hord] not enough guardcells for given prolongation operator order")

      call mpi_multigrid_bnd(coarse, iv, int((abs(ord_prolong)+1)/2, kind=4), extbnd_antimirror) ! exchange guardcells with corners

      if (dom%eff_dim<ndims) call die("[multigridbasefuncs:prolong_level_hord] 1D and 2D not finished")
      if (is_multicg) call die("[multigridbasefuncs:prolong_level_hord] multicg not implemented")

      select case (ord_prolong)
         case (O_D4)
            P_2 = 35./2048.; P_1 = -252./2048.; P0 = 1890./2048.; P1 = 420./2048.; P2 = -45./2048.
         case (O_D3)
            P_2 = 0.;        P_1 = -7./128.;    P0 = 105./128.;   P1 = 35./128.;   P2 = -5./128.
         case (O_D2)
            P_2 = 0.;        P_1 = -3./32.;     P0 = 30./32.;     P1 = 5./32.;     P2 = 0.
         case (O_LIN)
            P_2 = 0.;        P_1 = 0.;          P0 = 3./4.;       P1 = 1./4.;      P2 = 0.
         case (O_INJ)
            P_2 = 0.;        P_1 = 0.;          P0 = 1.;          P1 = 0.;         P2 = 0.
         case (O_I2)
            P_2 = 0.;        P_1 = -1./8.;      P0 = 1.;          P1 = 1./8.;      P2 = 0.
         case (O_I3)
            P_2 = 0.;        P_1 = -5./64.;     P0 = 55./64;      P1 = 17./64.;    P2 = -3./64.
         case (O_I4)
            P_2 = 3./128.;   P_1 = -11./64.;    P0 = 1.;          P1 = 11./64.;    P2 = -3./128.
            ! case 0 is handled through cg_list_level%prolong0_q_1var
         case default
            call die("[multigridbasefuncs:prolong_level_hord] Unsupported order")
            return
      end select

      ! this design doesn't allow multiple blocks or inter-process fine-coarse communication
      if (.not. associated(coarse)) call die("[multigridbasefuncs:prolong_level_hord] coarse == null()")
      if (.not. associated(coarse%finer)) call die("[multigridbasefuncs:prolong_level_hord] fine == null()")
      cg_c => coarse%first%cg
      cg_f => coarse%finer%first%cg

      ! convolve with the prolongation operator
      ! the two cases here are for optimization (see also call mpi_multigrid_bnd above)
      if (P_2 == 0. .and. P2 == 0.) then
         cg_f%prolong_x(         cg_f%is  :cg_f%ie-1:2, cg_c%js-1:cg_c%je+1, cg_c%ks-1:cg_c%ke+1) = &  ! x-odd cells
              + P1 * cg_c%q(iv)%arr(cg_c%is-1:cg_c%ie-1,   cg_c%js-1:cg_c%je+1, cg_c%ks-1:cg_c%ke+1)   &
              + P0 * cg_c%q(iv)%arr(cg_c%is  :cg_c%ie,     cg_c%js-1:cg_c%je+1, cg_c%ks-1:cg_c%ke+1)   &
              + P_1* cg_c%q(iv)%arr(cg_c%is+1:cg_c%ie+1,   cg_c%js-1:cg_c%je+1, cg_c%ks-1:cg_c%ke+1)
         cg_f%prolong_x(         cg_f%is+1:cg_f%ie:2,   cg_c%js-1:cg_c%je+1, cg_c%ks-1:cg_c%ke+1) = &  ! x-even cells
              + P_1* cg_c%q(iv)%arr(cg_c%is-1:cg_c%ie-1,   cg_c%js-1:cg_c%je+1, cg_c%ks-1:cg_c%ke+1)   &
              + P0 * cg_c%q(iv)%arr(cg_c%is  :cg_c%ie,     cg_c%js-1:cg_c%je+1, cg_c%ks-1:cg_c%ke+1)   &
              + P1 * cg_c%q(iv)%arr(cg_c%is+1:cg_c%ie+1,   cg_c%js-1:cg_c%je+1, cg_c%ks-1:cg_c%ke+1)

         cg_f%prolong_xy(           cg_f%is:cg_f%ie, cg_f%js  :cg_f%je-1:2, cg_c%ks-1:cg_c%ke+1) = &    ! y-odd cells
              + P1 * cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js-1:cg_c%je-1,   cg_c%ks-1:cg_c%ke+1)   &
              + P0 * cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js  :cg_c%je,     cg_c%ks-1:cg_c%ke+1)   &
              + P_1* cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js+1:cg_c%je+1,   cg_c%ks-1:cg_c%ke+1)
         cg_f%prolong_xy(           cg_f%is:cg_f%ie, cg_f%js+1:cg_f%je:2,   cg_c%ks-1:cg_c%ke+1) = &    ! y-even cells
              + P_1* cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js-1:cg_c%je-1,   cg_c%ks-1:cg_c%ke+1)   &
              + P0 * cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js  :cg_c%je,     cg_c%ks-1:cg_c%ke+1)   &
              + P1 * cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js+1:cg_c%je+1,   cg_c%ks-1:cg_c%ke+1)

         cg_f%q(iv)%arr(                cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_f%ks  :cg_f%ke-1:2) = &   ! z-odd cells
              + P1 * cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks-1:cg_c%ke-1  )   &
              + P0 * cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks  :cg_c%ke    )   &
              + P_1* cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks+1:cg_c%ke+1  )
         cg_f%q(iv)%arr(                cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_f%ks+1:cg_f%ke:2  ) = &   ! z-even cells
              + P_1* cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks-1:cg_c%ke-1  )   &
              + P0 * cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks  :cg_c%ke    )   &
              + P1 * cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks+1:cg_c%ke+1  )
      else
         cg_f%prolong_x(         cg_f%is  :cg_f%ie-1:2, cg_c%js-2:cg_c%je+2, cg_c%ks-2:cg_c%ke+2) = &  ! x-odd cells
              + P2 * cg_c%q(iv)%arr(cg_c%is-2:cg_c%ie-2,   cg_c%js-2:cg_c%je+2, cg_c%ks-2:cg_c%ke+2)   &
              + P1 * cg_c%q(iv)%arr(cg_c%is-1:cg_c%ie-1,   cg_c%js-2:cg_c%je+2, cg_c%ks-2:cg_c%ke+2)   &
              + P0 * cg_c%q(iv)%arr(cg_c%is  :cg_c%ie,     cg_c%js-2:cg_c%je+2, cg_c%ks-2:cg_c%ke+2)   &
              + P_1* cg_c%q(iv)%arr(cg_c%is+1:cg_c%ie+1,   cg_c%js-2:cg_c%je+2, cg_c%ks-2:cg_c%ke+2)   &
              + P_2* cg_c%q(iv)%arr(cg_c%is+2:cg_c%ie+2,   cg_c%js-2:cg_c%je+2, cg_c%ks-2:cg_c%ke+2)
         cg_f%prolong_x(         cg_f%is+1:cg_f%ie:2,   cg_c%js-2:cg_c%je+2, cg_c%ks-2:cg_c%ke+2) = &  ! x-even cells
              + P_2* cg_c%q(iv)%arr(cg_c%is-2:cg_c%ie-2,   cg_c%js-2:cg_c%je+2, cg_c%ks-2:cg_c%ke+2)   &
              + P_1* cg_c%q(iv)%arr(cg_c%is-1:cg_c%ie-1,   cg_c%js-2:cg_c%je+2, cg_c%ks-2:cg_c%ke+2)   &
              + P0 * cg_c%q(iv)%arr(cg_c%is  :cg_c%ie,     cg_c%js-2:cg_c%je+2, cg_c%ks-2:cg_c%ke+2)   &
              + P1 * cg_c%q(iv)%arr(cg_c%is+1:cg_c%ie+1,   cg_c%js-2:cg_c%je+2, cg_c%ks-2:cg_c%ke+2)   &
              + P2 * cg_c%q(iv)%arr(cg_c%is+2:cg_c%ie+2,   cg_c%js-2:cg_c%je+2, cg_c%ks-2:cg_c%ke+2)

         cg_f%prolong_xy(           cg_f%is:cg_f%ie, cg_f%js  :cg_f%je-1:2, cg_c%ks-2:cg_c%ke+2) = &    ! y-odd cells
              + P2 * cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js-2:cg_c%je-2,   cg_c%ks-2:cg_c%ke+2)   &
              + P1 * cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js-1:cg_c%je-1,   cg_c%ks-2:cg_c%ke+2)   &
              + P0 * cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js  :cg_c%je,     cg_c%ks-2:cg_c%ke+2)   &
              + P_1* cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js+1:cg_c%je+1,   cg_c%ks-2:cg_c%ke+2)   &
              + P_2* cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js+2:cg_c%je+2,   cg_c%ks-2:cg_c%ke+2)
         cg_f%prolong_xy(           cg_f%is:cg_f%ie, cg_f%js+1:cg_f%je:2,   cg_c%ks-2:cg_c%ke+2) = &    ! y-even cells
              + P_2* cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js-2:cg_c%je-2,   cg_c%ks-2:cg_c%ke+2)   &
              + P_1* cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js-1:cg_c%je-1,   cg_c%ks-2:cg_c%ke+2)   &
              + P0 * cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js  :cg_c%je,     cg_c%ks-2:cg_c%ke+2)   &
              + P1 * cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js+1:cg_c%je+1,   cg_c%ks-2:cg_c%ke+2)   &
              + P2 * cg_f%prolong_x(cg_f%is:cg_f%ie, cg_c%js+2:cg_c%je+2,   cg_c%ks-2:cg_c%ke+2)

         cg_f%q(iv)%arr(                cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_f%ks  :cg_f%ke-1:2) = &   ! z-odd cells
              + P2 * cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks-2:cg_c%ke-2  )   &
              + P1 * cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks-1:cg_c%ke-1  )   &
              + P0 * cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks  :cg_c%ke    )   &
              + P_1* cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks+1:cg_c%ke+1  )   &
              + P_2* cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks+2:cg_c%ke+2  )
         cg_f%q(iv)%arr(                cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_f%ks+1:cg_f%ke:2  ) = &   ! z-even cells
              + P_2* cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks-2:cg_c%ke-2  )   &
              + P_1* cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks-1:cg_c%ke-1  )   &
              + P0 * cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks  :cg_c%ke    )   &
              + P1 * cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks+1:cg_c%ke+1  )   &
              + P2 * cg_f%prolong_xy(cg_f%is:cg_f%ie, cg_f%js:cg_f%je, cg_c%ks+2:cg_c%ke+2  )
      endif

   end subroutine prolong_level_hord

end module multigridbasefuncs
