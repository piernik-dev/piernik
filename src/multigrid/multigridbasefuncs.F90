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

      use multigridvars, only: plvl

      implicit none

      type(plvl), pointer, intent(in) :: curl  !< level for which clear the boundary values

      curl%bnd_x(:,:,:) = 0.
      curl%bnd_y(:,:,:) = 0.
      curl%bnd_z(:,:,:) = 0.

   end subroutine zero_boundaries

!!$ ============================================================================
!>
!! \brief Multigrid elementary operators: prolongation, restriction, norm etc.
!<

   subroutine prolong_level(coarse, iv)

      use dataio_pub,            only: die
      use multigridhelpers,      only: dirty_debug, check_dirty, dirtyH
      use multigridvars,         only: plvl, roof, ord_prolong, ngridvars, is_mg_uneven
      use multigridexperimental, only: prolong_level_hord

      implicit none

      type(plvl), pointer, intent(in) :: coarse !< level to prolong from
      integer(kind=4), intent(in) :: iv    !< variable to be prolonged

      type(plvl), pointer :: fine

      if (associated(coarse, roof)) return ! can't prolong finest level
      if (iv < 1 .or. iv > ngridvars) call die("[multigridbasefuncs:prolong_level] Invalid variable index.")

      if (.not. associated(coarse)) call die("[multigridbasefuncs:prolong_level] coarse == null()")
      fine   => coarse%finer
      if (.not. associated(fine)) call die("[multigridbasefuncs:prolong_level] fine == null()")

      if (dirty_debug) fine%mgvar(:, :, :, iv) = dirtyH

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

      use dataio_pub,         only: die
      use multigridhelpers,   only: check_dirty
      use multigridvars,      only: roof, base, ngridvars, plvl

      implicit none

      integer(kind=4), intent(in)      :: iv    !< variable to be restricted

      type(plvl), pointer :: curl

      if (iv < 1 .or. iv > ngridvars) call die("[multigridbasefuncs:restrict_all] Invalid variable index.")

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

      use constants,     only: GEO_XYZ, GEO_RPZ, I_ONE
      use dataio_pub,    only: die
      use domain,        only: geometry_type
      use mpi,           only: MPI_DOUBLE_PRECISION, MPI_SUM
      use mpisetup,      only: comm, ierr
      use multigridvars, only: ngridvars, roof

      implicit none

      integer(kind=4), intent(in)  :: iv   !< index of variable in lvl()%mgvar for which we want to find the norm
      real,    intent(out) :: norm !< the calculated norm

      real                 :: lsum
      integer              :: i

      if (iv <= 0 .or. iv > ngridvars) call die("[multigridbasefuncs:norm_sq] Invalid variable index")

      select case (geometry_type)
         case (GEO_XYZ)
            lsum = sum(roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, iv)**2) * roof%dvol
         case (GEO_RPZ)
            lsum = 0.
            do i = roof%is, roof%ie
               lsum = lsum + sum(roof%mgvar(i, roof%js:roof%je, roof%ks:roof%ke, iv)**2) * roof%dvol * roof%x(i)
            enddo
         case default
            call die("[multigridbasefuncs:norm_sq] Unsupported geometry.")
      end select
      call MPI_Allreduce(lsum, norm, I_ONE, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
      norm = sqrt(norm)

   end subroutine norm_sq

!!$ ============================================================================
!>
!! \brief Compute the global average value and subtract it from the whole domain
!<

   subroutine subtract_average(curl, iv)

      use constants,     only: GEO_XYZ, GEO_RPZ, I_ONE
      use dataio_pub,    only: die
      use domain,        only: geometry_type
      use mpi,           only: MPI_DOUBLE_PRECISION, MPI_SUM
      use mpisetup,      only: comm, ierr
      use multigridvars, only: plvl, ngridvars

      implicit none

      type(plvl), pointer, intent(in) :: curl !< level for which we want to subtract its average from
      integer(kind=4), intent(in) :: iv   !< index of variable in lvl()%mgvar which we want to have zero average

      real                :: lsum, avg, vol
      integer             :: i

      if (iv < 1 .or. iv > ngridvars) call die("[multigridbasefuncs:subtract_average] Invalid variable index.")

      select case (geometry_type)
         case (GEO_XYZ)
            lsum = sum(curl%mgvar(curl%is:curl%ie, curl%js:curl%je, curl%ks:curl%ke, iv)) * curl%dvol
         case (GEO_RPZ)
            lsum = 0.
            do i = curl%is, curl%ie
               lsum = lsum + sum(curl%mgvar(i, curl%js:curl%je, curl%ks:curl%ke, iv)) * curl%dvol * curl%x(i)
            enddo
         case default
            call die("[multigridbasefuncs:subtract_average] Unsupported geometry.")
      end select
      call MPI_Allreduce(lsum, avg, I_ONE, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
      call MPI_Allreduce(curl%vol, vol, I_ONE, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr) ! This probably can be calculated in init_multigrid
      avg = avg / vol

      curl%mgvar(:, :, :, iv) = curl%mgvar(:, :, :, iv) - avg

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
!! \todo move to plvl or multigrid_gravity
!!
!! OPT: completely unoptimized
!<

   subroutine prolong_faces(fine, soln)

      use constants,          only: xdim, ydim, zdim, LO, HI, LONG, I_ONE, half
      use dataio_pub,         only: die, warn
      use domain,             only: has_dir, D_x, D_y, D_z
      use mpi,                only: MPI_DOUBLE_PRECISION
      use mpisetup,           only: proc, comm, ierr, req, status, master
      use multigridhelpers,   only: check_dirty
      use multigridmpifuncs,  only: mpi_multigrid_bnd
      use multigridvars,      only: plvl, ord_prolong_face_norm, ord_prolong_face_par, base, extbnd_antimirror, is_external, need_general_pf, pr_segment

      implicit none

      type(plvl), pointer, intent(in) :: fine !< level for which approximate the solution
      integer(kind=4), intent(in) :: soln !< index of solution in lvl()%mgvar ! \todo change the name

      integer                       :: i, j, k, d, lh, g, g1, gc, ib, jb, ibh, jbh, l
      type(plvl), pointer           :: coarse
      integer, parameter            :: s_wdth  = 3           ! interpolation stencil width
      integer, parameter            :: s_rng = (s_wdth-1)/2  ! stencil range around 0
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
      integer :: b_rng
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

         do lh = LO, HI ! \todo convert plvl%bnd_[xyz] to an array and make obsolete the following pointer assignments
            p_bnd(xdim, lh)%bnd => fine%bnd_x(:, :, lh)
            p_bnd(ydim, lh)%bnd => fine%bnd_y(:, :, lh)
            p_bnd(zdim, lh)%bnd => fine%bnd_z(:, :, lh)
            do d = xdim, zdim
               p_bnd(d, lh)%bnd(:, :) = 0. ! \todo mark dirty somehow, eg use a "magic" small value, like 1.234567890123456789e-123
            enddo
         enddo

         nr = 0

         ! Gather coarse data for internal boundaries from others
         do d = xdim, zdim
            do lh = LO, HI
               if (allocated(fine%pfc_tgt(d, lh)%seg)) then
                  do g = 1, ubound(fine%pfc_tgt(d, lh)%seg(:), dim=1)
                     if (fine%pfc_tgt(d, lh)%seg(g)%proc /= proc) then
                        nr = nr + I_ONE
                        call MPI_Irecv(fine%pfc_tgt(d, lh)%seg(g)%buf(1, 1, 1), size(fine%pfc_tgt(d, lh)%seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, &
                             &         fine%pfc_tgt(d, lh)%seg(g)%proc, HI*d+lh, comm, req(nr), ierr)
                     endif
                  enddo
               endif
            enddo
         enddo

         ! Send own coarse data to others
         off1(:) = int(mod(fine%off(:), 2_LONG), kind=4)
         do d = xdim, zdim
            do lh = LO, HI
               if (allocated(coarse%pff_tgt(d, lh)%seg)) then
                  do g = 1, ubound(coarse%pff_tgt(d, lh)%seg(:), dim=1)

                     pseg => coarse%pff_tgt(d, lh)%seg(g)
                     cse => pseg%se

                     if (pseg%proc /= proc) then
                        nr = nr + I_ONE
                        se(:,:) = cse(:,:)
                        pseg%buf(:, :, :) = 0. ! this can be avoided by extracting first assignment from the loop
                        do l = 1, ubound(pseg%f_lay(:), dim=1)
                           se(d,:) = pseg%f_lay(l)%layer
                           pseg%buf(:, :, :) = pseg%buf(:, :, :) + pseg%f_lay(l)%coeff * &
                                coarse%mgvar(se(xdim, LO):se(xdim, HI), se(ydim, LO):se(ydim, HI), se(zdim, LO):se(zdim, HI), soln)
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
               if (allocated(fine%pfc_tgt(d, lh)%seg)) then
                  do g = 1, ubound(fine%pfc_tgt(d, lh)%seg(:), dim=1)
                     if (fine%pfc_tgt(d, lh)%seg(g)%proc == proc) then
                        nullify(cse)
                        gc = 0
                        do g1 = 1, ubound(coarse%pff_tgt(d, lh)%seg(:), dim=1) !> \todo should be set up in mpi_multigrid_prep_grav
                           if (coarse%pff_tgt(d, lh)%seg(g1)%proc == proc) then
                              if (.not. associated(cse)) then
                                 cse => coarse%pff_tgt(d, lh)%seg(g1)%se
                                 gc = g1
                              else
                                 call die("[multigridbasefuncs:prolong_faces] multiple local coarse grid targets")
                              endif
                           endif
                        enddo
                        if (.not. associated(cse)) call die("[multigridbasefuncs:prolong_faces] missing coarse grid target")

                        se(:,:) = cse(:,:)
                        pseg => fine%pfc_tgt(d, lh)%seg(g)
                        pseg%buf(:,:,:) = 0. ! this can be avoided by extracting first assignment from the loop
                        do l = 1, ubound(coarse%pff_tgt(d, lh)%seg(gc)%f_lay(:), dim=1)
                           se(d,:) = coarse%pff_tgt(d, lh)%seg(gc)%f_lay(l)%layer
                           pseg%buf(:,:,:) =  pseg%buf(:,:,:) + coarse%pff_tgt(d, lh)%seg(gc)%f_lay(l)%coeff * &
                                coarse%mgvar(se(xdim, LO):se(xdim, HI), se(ydim, LO):se(ydim, HI), se(zdim, LO):se(zdim, HI), soln)
                        enddo
                     endif
                  enddo
               endif
            enddo
         enddo

         if (nr>0) call MPI_Waitall(nr, req(:nr), status(:,:nr), ierr)

         ! Interpolate content of buffers to boundary face-layes fine%bnd_[xyz](:, :, :)
         do d = xdim, zdim
            if (has_dir(d)) then
               do lh = LO, HI
                  ! make external homogenous Dirichlet boundaries, where applicable, interpolate (inject) otherwise.
                  ! \todo for homogenous neumann boundary lay2 should be .false. and pc_bnd(1, :, :) should contain layer of the coarse data next to the external boundary
                  if (is_external(d, lh)) then
                     p_bnd(d, lh)%bnd(:,:) = 0. ! BEWARE homogenous Dirichlet by default
                  else
                     p_bnd(d, lh)%bnd(:,:) = 0
                     if (allocated(fine%pfc_tgt(d, lh)%seg)) then
                        do g = 1, ubound(fine%pfc_tgt(d, lh)%seg(:), dim=1)

                           ii(d) = 1
                           do i = 1, ubound(fine%pfc_tgt(d, lh)%seg(g)%buf, dim=d1(d))
                              ii(d1(d)) = i
                              do j = 1, ubound(fine%pfc_tgt(d, lh)%seg(g)%buf, dim=d2(d))
                                 ii(d2(d)) = j

! ? + off (0 or 1)
                                 ib = 2*i - 1 + int(fine%pfc_tgt(d, lh)%seg(g)%se(d1(d), LO), kind=4) - int(mod(fine%off(d1(d)), 2_LONG), 4)
                                 jb = 2*j - 1 + int(fine%pfc_tgt(d, lh)%seg(g)%se(d2(d), LO), kind=4) - int(mod(fine%off(d2(d)), 2_LONG), 4)

                                 ibh = ib
                                 if (has_dir(d1(d)) .and. ubound(p_bnd(d, lh)%bnd(:,:), dim=1) > ib) ibh = ib + 1
                                 if (has_dir(d1(d)) .and. lbound(p_bnd(d, lh)%bnd(:,:), dim=1) > ib) ib  = ib + 1
                                 jbh = jb
                                 if (has_dir(d2(d)) .and. ubound(p_bnd(d, lh)%bnd(:,:), dim=2) > jb) jbh = jb + 1
                                 if (has_dir(d2(d)) .and. lbound(p_bnd(d, lh)%bnd(:,:), dim=2) > jb) jb  = jb + 1

                                 p_bnd(d, lh)%bnd(ib:ibh, jb:jbh) = p_bnd(d, lh)%bnd(ib:ibh, jb:jbh) + fine%pfc_tgt(d, lh)%seg(g)%buf(ii(xdim), ii(ydim), ii(zdim))
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
         if (ord_prolong_face_norm > 0) b_rng = max(b_rng, ord_prolong_face_norm+1)
         call mpi_multigrid_bnd(coarse, soln, b_rng, extbnd_antimirror, corners=(ord_prolong_face_par/=0)) !> \deprecated BEWARE for higher prolongation order more guardcell are required
         call check_dirty(coarse, soln, "prolong_faces", s_rng)

         if (ord_prolong_face_norm > 0) then
            if (has_dir(xdim)) then
               pp_norm = 2.*sum(pp(-D_y:D_y, -D_z:D_z, 1, 1)) ! normalization is required for ord_prolong_face_par == 1 and -2
               do j = coarse%js, coarse%je
                  do k = coarse%ks, coarse%ke
                     fine%bnd_x(-fine%js+2*j,    -fine%ks+2*k,    LO)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 1) * ( &
                          opfn1*(coarse%mgvar(coarse%is,  j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-1,j-D_y:j+D_y,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(coarse%is+1,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-2,j-D_y:j+D_y,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_x(-fine%js+2*j+D_y,-fine%ks+2*k,    LO)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 1) * ( &
                          opfn1*(coarse%mgvar(coarse%is,  j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-1,j-D_y:j+D_y,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(coarse%is+1,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-2,j-D_y:j+D_y,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_x(-fine%js+2*j,    -fine%ks+2*k+D_z,LO)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 2) * ( &
                          opfn1*(coarse%mgvar(coarse%is  ,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-1,j-D_y:j+D_y,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(coarse%is+1,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-2,j-D_y:j+D_y,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_x(-fine%js+2*j+D_y,-fine%ks+2*k+D_z,LO)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 2) * ( &
                          opfn1*(coarse%mgvar(coarse%is,  j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-1,j-D_y:j+D_y,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(coarse%is+1,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-2,j-D_y:j+D_y,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_x(-fine%js+2*j,    -fine%ks+2*k,    HI)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 1) * ( &
                          opfn1*(coarse%mgvar(coarse%ie,  j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+1,j-D_y:j+D_y,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(coarse%ie-1,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+2,j-D_y:j+D_y,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_x(-fine%js+2*j+D_y,-fine%ks+2*k,    HI)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 1) * ( &
                          opfn1*(coarse%mgvar(coarse%ie,  j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+1,j-D_y:j+D_y,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(coarse%ie-1,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+2,j-D_y:j+D_y,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_x(-fine%js+2*j,    -fine%ks+2*k+D_z,HI)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 2) * ( &
                          opfn1*(coarse%mgvar(coarse%ie,  j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+1,j-D_y:j+D_y,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(coarse%ie-1,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+2,j-D_y:j+D_y,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_x(-fine%js+2*j+D_y,-fine%ks+2*k+D_z,HI)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 2) * ( &
                          opfn1*(coarse%mgvar(coarse%ie,  j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+1,j-D_y:j+D_y,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(coarse%ie-1,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+2,j-D_y:j+D_y,k-D_z:k+D_z,soln)))) / pp_norm
                  enddo
               enddo
            endif

            if (has_dir(ydim)) then
               pp_norm = 2.*sum(pp(-D_x:D_x, -D_z:D_z, 1, 1))
               do i = coarse%is, coarse%ie
                  do k = coarse%ks, coarse%ke
                     fine%bnd_y(-fine%is+2*i,    -fine%ks+2*k,    LO)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 1) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,coarse%js,  k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-1,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,coarse%js+1,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-2,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_y(-fine%is+2*i+D_x,-fine%ks+2*k,    LO)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 1) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,coarse%js,  k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-1,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,coarse%js+1,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-2,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_y(-fine%is+2*i,    -fine%ks+2*k+D_z,LO)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 2) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,coarse%js,  k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-1,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,coarse%js+1,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-2,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_y(-fine%is+2*i+D_x,-fine%ks+2*k+D_z,LO)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 2) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,coarse%js,  k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-1,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,coarse%js+1,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-2,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_y(-fine%is+2*i,    -fine%ks+2*k,    HI)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 1) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,coarse%je,  k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+1,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,coarse%je-1,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+2,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_y(-fine%is+2*i+D_x,-fine%ks+2*k,    HI)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 1) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,coarse%je,  k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+1,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,coarse%je-1,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+2,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_y(-fine%is+2*i,    -fine%ks+2*k+D_z,HI)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 2) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,coarse%je,  k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+1,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,coarse%je-1,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+2,k-D_z:k+D_z,soln)))) / pp_norm
                     fine%bnd_y(-fine%is+2*i+D_x,-fine%ks+2*k+D_z,HI)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 2) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,coarse%je,  k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+1,k-D_z:k+D_z,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,coarse%je-1,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+2,k-D_z:k+D_z,soln)))) / pp_norm
                  enddo
               enddo
            endif

            if (has_dir(zdim)) then
               pp_norm = 2.*sum(pp(-D_x:D_x, -D_y:D_y, 1, 1))
               do i = coarse%is, coarse%ie
                  do j = coarse%js, coarse%je
                     fine%bnd_z(-fine%is+2*i,    -fine%js+2*j,    LO)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 1) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks,  soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-1,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks+1,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-2,soln)))) / pp_norm
                     fine%bnd_z(-fine%is+2*i+D_x,-fine%js+2*j,    LO)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 1) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks,  soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-1,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks+1,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-2,soln)))) / pp_norm
                     fine%bnd_z(-fine%is+2*i,    -fine%js+2*j+D_y,LO)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 2) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks,  soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-1,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks+1,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-2,soln)))) / pp_norm
                     fine%bnd_z(-fine%is+2*i+D_x,-fine%js+2*j+D_y,LO)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 2) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks,  soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-1,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks+1,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-2,soln)))) / pp_norm
                     fine%bnd_z(-fine%is+2*i,    -fine%js+2*j,    HI)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 1) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke,  soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+1,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke-1,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+2,soln)))) / pp_norm
                     fine%bnd_z(-fine%is+2*i+D_x,-fine%js+2*j,    HI)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 1) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke,  soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+1,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke-1,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+2,soln)))) / pp_norm
                     fine%bnd_z(-fine%is+2*i,    -fine%js+2*j+D_y,HI)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 2) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke,  soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+1,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke-1,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+2,soln)))) / pp_norm
                     fine%bnd_z(-fine%is+2*i+D_x,-fine%js+2*j+D_y,HI)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 2) * ( &
                          opfn1*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke,  soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+1,soln)) + &
                          opfn3*(coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke-1,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+2,soln)))) / pp_norm
                  enddo
               enddo
            endif

         else

            if (has_dir(xdim)) then
               pp_norm = 2.*sum(pp(-D_y:D_y, -D_z:D_z, 1, 1)) ! normalization is required for ord_prolong_face_par == 1 and -2
               do j = coarse%js, coarse%je
                  do k = coarse%ks, coarse%ke
                     fine%bnd_x(-fine%js+2*j,    -fine%ks+2*k,    LO)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 1) * (coarse%mgvar(coarse%is,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
                     fine%bnd_x(-fine%js+2*j+D_y,-fine%ks+2*k,    LO)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 1) * (coarse%mgvar(coarse%is,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
                     fine%bnd_x(-fine%js+2*j,    -fine%ks+2*k+D_z,LO)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 2) * (coarse%mgvar(coarse%is,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
                     fine%bnd_x(-fine%js+2*j+D_y,-fine%ks+2*k+D_z,LO)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 2) * (coarse%mgvar(coarse%is,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
                     fine%bnd_x(-fine%js+2*j,    -fine%ks+2*k,    HI)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 1) * (coarse%mgvar(coarse%ie,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
                     fine%bnd_x(-fine%js+2*j+D_y,-fine%ks+2*k,    HI)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 1) * (coarse%mgvar(coarse%ie,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
                     fine%bnd_x(-fine%js+2*j,    -fine%ks+2*k+D_z,HI)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 2) * (coarse%mgvar(coarse%ie,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
                     fine%bnd_x(-fine%js+2*j+D_y,-fine%ks+2*k+D_z,HI)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 2) * (coarse%mgvar(coarse%ie,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
                  enddo
               enddo
            endif

            if (has_dir(ydim)) then
               pp_norm = 2.*sum(pp(-D_x:D_x, -D_z:D_z, 1, 1))
               do i = coarse%is, coarse%ie
                  do k = coarse%ks, coarse%ke
                     fine%bnd_y(-fine%is+2*i,    -fine%ks+2*k,    LO)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 1) * (coarse%mgvar(i-D_x:i+D_x,coarse%js,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-1,k-D_z:k+D_z,soln))) / pp_norm
                     fine%bnd_y(-fine%is+2*i+D_x,-fine%ks+2*k,    LO)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 1) * (coarse%mgvar(i-D_x:i+D_x,coarse%js,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-1,k-D_z:k+D_z,soln))) / pp_norm
                     fine%bnd_y(-fine%is+2*i,    -fine%ks+2*k+D_z,LO)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 2) * (coarse%mgvar(i-D_x:i+D_x,coarse%js,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-1,k-D_z:k+D_z,soln))) / pp_norm
                     fine%bnd_y(-fine%is+2*i+D_x,-fine%ks+2*k+D_z,LO)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 2) * (coarse%mgvar(i-D_x:i+D_x,coarse%js,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-1,k-D_z:k+D_z,soln))) / pp_norm
                     fine%bnd_y(-fine%is+2*i,    -fine%ks+2*k,    HI)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 1) * (coarse%mgvar(i-D_x:i+D_x,coarse%je,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+1,k-D_z:k+D_z,soln))) / pp_norm
                     fine%bnd_y(-fine%is+2*i+D_x,-fine%ks+2*k,    HI)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 1) * (coarse%mgvar(i-D_x:i+D_x,coarse%je,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+1,k-D_z:k+D_z,soln))) / pp_norm
                     fine%bnd_y(-fine%is+2*i,    -fine%ks+2*k+D_z,HI)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 2) * (coarse%mgvar(i-D_x:i+D_x,coarse%je,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+1,k-D_z:k+D_z,soln))) / pp_norm
                     fine%bnd_y(-fine%is+2*i+D_x,-fine%ks+2*k+D_z,HI)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 2) * (coarse%mgvar(i-D_x:i+D_x,coarse%je,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+1,k-D_z:k+D_z,soln))) / pp_norm
                  enddo
               enddo
            endif

            if (has_dir(zdim)) then
               pp_norm = 2.*sum(pp(-D_x:D_x, -D_y:D_y, 1, 1))
               do i = coarse%is, coarse%ie
                  do j = coarse%js, coarse%je
                     fine%bnd_z(-fine%is+2*i,    -fine%js+2*j,    LO)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 1) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-1,soln))) / pp_norm
                     fine%bnd_z(-fine%is+2*i+D_x,-fine%js+2*j,    LO)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 1) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-1,soln))) / pp_norm
                     fine%bnd_z(-fine%is+2*i,    -fine%js+2*j+D_y,LO)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 2) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-1,soln))) / pp_norm
                     fine%bnd_z(-fine%is+2*i+D_x,-fine%js+2*j+D_y,LO)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 2) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-1,soln))) / pp_norm
                     fine%bnd_z(-fine%is+2*i,    -fine%js+2*j,    HI)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 1) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+1,soln))) / pp_norm
                     fine%bnd_z(-fine%is+2*i+D_x,-fine%js+2*j,    HI)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 1) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+1,soln))) / pp_norm
                     fine%bnd_z(-fine%is+2*i,    -fine%js+2*j+D_y,HI)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 2) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+1,soln))) / pp_norm
                     fine%bnd_z(-fine%is+2*i+D_x,-fine%js+2*j+D_y,HI)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 2) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+1,soln))) / pp_norm
                  enddo
               enddo
            endif
         endif
      endif

   end subroutine prolong_faces

end module multigridbasefuncs
