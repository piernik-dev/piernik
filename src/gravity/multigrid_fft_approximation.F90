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
#include "macros.h"

!>
!! \brief Multigrid Poisson solver: Approximate local solver based on FFT
!!
!! \details Local solver based on FFT may reduce the amount of required communication because it relies only on proided boundary values
!! (i.e. does not need the whole grid to be prolonged)
!!
!! \warning These routines need to be reimplemented in order to work properly in most general grid decompositions.
!<

module multigrid_fftapprox
! pulled by MULTIGRID && GRAV
#if defined(__INTEL_COMPILER)
      !! \deprecated remove this clause as soon as Intel Compiler gets required
      !! features and/or bug fixes
   use cg_list_bnd,   only: cg_list_bnd_T   ! QA_WARN intel
#endif /* __INTEL_COMPILER */

   implicit none

   private
   public :: approximate_solution_fft, mpi_multigrid_prep_grav, fft_convolve

   include "fftw3.f"
   ! constants from fftw3.f
   !   integer, parameter :: FFTW_MEASURE=0, FFTW_PATIENT=32, FFTW_ESTIMATE=64
   !   integer, parameter :: FFTW_RODFT01=8, FFTW_RODFT10=9

contains

!>
!! \brief Set up communication for face prolongation
!!
!! \todo implement also prolongation of coarsened multipoles
!!
!! \todo move this to cg_level_connected_T
!<

   subroutine mpi_multigrid_prep_grav

      use constants,     only: xdim, ydim, zdim, ndims, LO, HI, LONG, zero, one, half, O_INJ, O_LIN, O_I2, INT4
      use dataio_pub,    only: warn, die
      use domain,        only: dom
      use cg_list,       only: cg_list_element
      use cg_level_connected, only: cg_level_connected_T, coarsest
      use grid_cont,     only: pr_segment, grid_container, is_overlap
      use mpisetup,      only: proc, nproc, FIRST, LAST, procmask, inflate_req
      use multigridvars, only: ord_prolong_face_norm, need_general_pf
#ifdef DEBUG
      use constants,     only: two
      use piernikdebug,  only: aux_R
#endif /* DEBUG */

      implicit none

      integer :: d, g, j, lh, hl, l, nl
      integer(kind=8), dimension(xdim:zdim) :: ijks, per
      logical, dimension(xdim:zdim) :: dmask
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: coarsened, b_layer
      type(pr_segment), pointer :: seg
      type(cg_level_connected_T), pointer   :: curl        !< current level (a pointer sliding along the linked list)
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg            !< current grid container
      logical :: is_internal_fine
      integer, parameter :: max_opfn = 2
      real, dimension(0:max_opfn) :: opfn_c_ff, opfn_c_cf

      if (.not. need_general_pf) return

      call inflate_req(size([LO, HI]) * 2 * nproc * ndims)

      if (ord_prolong_face_norm > max_opfn) ord_prolong_face_norm = max_opfn
      select case (ord_prolong_face_norm)
         case (O_INJ)
            opfn_c_ff(:) = [ half, zero, zero ]
            opfn_c_cf(:) = [ one,  zero, zero ]
         case (O_LIN)
#ifdef DEBUG
            ! the maximum convergence is at aux_R(1) = 0.25 +/- 0.05
            opfn_c_ff(:) = half * [ one + aux_R(1), -aux_R(1), zero ]
            ! the maximum convergence is at aux_R(2) = -0.05 +/- 0.05 and 0
            opfn_c_cf(:) = [ one + two*aux_R(2), -aux_R(2), zero ]
#else /* !DEBUG */
            opfn_c_ff(:) = [  5., -one, zero ] / 8.  ! adjusted experimentally
            opfn_c_cf(:) = [ 18.,  one, zero ] / 20. ! adjusted experimentally
#endif /* !DEBUG */
         case (O_I2)
#ifdef DEBUG
            ! the maximum convergence is at aux_R(1) = 0.35 +/- 0.1 and aux_R(3) = 0.20 +- 0.05
            opfn_c_ff(:) = half * [ one+aux_R(1), -two*aux_R(1)+aux_R(3), aux_R(1)-aux_R(3) ]
            ! the maximum convergence is at aux_R(2) = -0.01 +/- 0.02 and aux_R(4) = -0.125 +- 0.025
            opfn_c_cf(:) = [ one + two*aux_R(2), -two*aux_R(2)+aux_R(4), aux_R(2)-aux_R(4) ]
#else /* !DEBUG */
            opfn_c_ff(:) = [ 27., -10., 3.  ] / 40. ! adjusted experimentally
            opfn_c_cf(:) = [  8., -one, one ] / 8.
#endif /* !DEBUG */
         case default
            call die("mg:mmpg opfn_c_[cf]f(:)")
      end select

      curl => coarsest
      do while (associated(curl))

         if (ubound(curl%pse(proc)%c(:), dim=1) > 1) call die("[multigrid_fft_approximation:mpi_multigrid_prep_grav] Multiple blocks per process not implemented yet")

         per(:) = 0
         where (dom%periodic(:)) per(:) = curl%n_d(:)

         do d = xdim, zdim
            dmask(:) = .false.
            dmask(d) = .true.
            do lh = LO, HI
               if (dom%has_dir(d)) then
                  hl = LO+HI-lh

                  ! find coarse target for receiving data to be prolonged
                  if (associated(curl%coarser)) then
                     cgl => curl%first
                     do while (associated(cgl))
                        cg => cgl%cg
                        if (.not. cg%ext_bnd(d, lh)) then
                           ijks(:) = cg%ijkse(:, LO) - cg%off(:)  ! add this to convert absolute cell coordinates to local indices. (+nb - off(:))
                           procmask(:) = 0
                           ! two layers of cells are required for even locations (+ two layers per each interpolation order)
                           ! one layer of cells is required for odd locations (the local domain face is exactly at the centers of coarse cells)
                           coarsened(:, :) = cg%my_se(:, :)/2
                           coarsened(d, hl) = coarsened(d, lh)
                           select case (lh)
                              case (LO)
                                 if (mod(cg%off(d),    2_LONG) == 0) coarsened(d, :) = coarsened(d, :) + [ -1_INT4-ord_prolong_face_norm,        ord_prolong_face_norm ]
                              case (HI)
                                 if (mod(cg%h_cor1(d), 2_LONG) == 0) coarsened(d, :) = coarsened(d, :) + [        -ord_prolong_face_norm, 1_INT4+ord_prolong_face_norm ]
                           end select

                           do j = FIRST, LAST
                              if (ubound(curl%coarser%pse(j)%c, dim=1) > 0) then
                                 !> \warning: it is not anymore guaranteed that curl%coarser%pse(j)%c(1)%se(:,:) exists
                                 ! Typically it does on coarser level (counterintuitive but true), but not necessarily on the finer one
                                 if (is_overlap(coarsened(:,:), curl%coarser%pse(j)%c(1)%se(:,:), per)) procmask(j) = 1
                              endif
                           enddo
                           allocate(cg%mg%pfc_tgt(d, lh)%seg(count(procmask(:) /= 0)))

                           g = 0
                           do j = FIRST, LAST
                              if (procmask(j) /= 0) then
                                 g = g + 1
                                 if (.not. allocated(cg%mg%pfc_tgt(d, lh)%seg) .or. g>ubound(cg%mg%pfc_tgt(d, lh)%seg, dim=1)) call die("mg:mmpg pfc_tgt g>")
                                 seg => cg%mg%pfc_tgt(d, lh)%seg(g)
                                 if (allocated(seg%buf)) then
                                    call warn("mg:mmpg o seg%buf a a")
                                    deallocate(seg%buf)
                                 endif
                                 seg%proc = j
                                 ! find cross-section of own face segment with refined coarse segment
                                 b_layer(:, :) = cg%my_se(:, :)
                                 b_layer(d, hl) = b_layer(d, lh)
                                 !b_layer(d, lh) = b_layer(d, lh) + 2*lh-LO-HI ! extend to two layers of buffer

                                 where (.not. dmask(:)) ! find extents perpendicular to d
                                    seg%se(:, LO) = max(b_layer(:, LO), curl%coarser%pse(j)%c(1)%se(:, LO)*2  )
                                    seg%se(:, HI) = min(b_layer(:, HI), curl%coarser%pse(j)%c(1)%se(:, HI)*2+1)
                                 endwhere
                                 seg%se(d, :) = b_layer(d, :)
                                 allocate(seg%buf(seg%se(xdim, HI)/2-seg%se(xdim, LO)/2 + 1, &
                                      &           seg%se(ydim, HI)/2-seg%se(ydim, LO)/2 + 1, &
                                      &           seg%se(zdim, HI)/2-seg%se(zdim, LO)/2 + 1))

                                 seg%se(:, LO) = seg%se(:, LO) - cg%off(:) !+ ijks(:)
                                 seg%se(:, HI) = seg%se(:, HI) - cg%off(:) !+ ijks(:)
                              endif
                           enddo
                        endif

                        cgl => cgl%nxt
                     enddo
                  endif

                  ! find fine target(s) for sending the data to be prolonged
                  if (associated(curl%finer)) then
                     cgl => curl%first
                     do while (associated(cgl))
                        cg => cgl%cg
                        ijks(:) = cg%ijkse(:, LO) - cg%off(:)  ! add this to convert absolute cell coordinates to local indices. (+nb - off(:))
                        procmask(:) = 0
                        do j = FIRST, LAST
                           is_internal_fine = dom%periodic(d)
                           if (ubound(curl%finer%pse(j)%c, dim=1) > 0) then
                              coarsened(:, :) = curl%finer%pse(j)%c(1)%se(:,:)/2
                              coarsened(d, hl) = coarsened(d, lh)
                              select case (lh)
                                 case (LO)
                                    if (mod(curl%finer%pse(j)%c(1)%se(d, LO),     2_LONG) == 0) &
                                         coarsened(d, :) = coarsened(d, :) + [ -1_INT4-ord_prolong_face_norm,        ord_prolong_face_norm ]
                                    is_internal_fine = is_internal_fine .or. (curl%finer%pse(j)%c(1)%se(d, lh) /= 0)
                                 case (HI)
                                    if (mod(curl%finer%pse(j)%c(1)%se(d, HI) + 1, 2_LONG) == 0) &
                                         coarsened(d, :) = coarsened(d, :) + [        -ord_prolong_face_norm, 1_INT4+ord_prolong_face_norm ]
                                    is_internal_fine = is_internal_fine .or. (curl%finer%pse(j)%c(1)%se(d, lh) + 1 < curl%finer%n_d(d))
                              end select
                           else
                              coarsened(:,:) = -2*dom%nb
                           endif
                           if (is_internal_fine) then
                              if (is_overlap(coarsened(:, :), cg%my_se(:, :), per)) procmask(j) = 1
                           endif
                        enddo
                        allocate(cg%mg%pff_tgt(d, lh)%seg(count(procmask(:) /= 0)))

                        g = 0
                        do j = FIRST, LAST
                           if (procmask(j) /= 0) then
                              g = g + 1
                              if (.not. allocated(cg%mg%pff_tgt(d, lh)%seg) .or. g>ubound(cg%mg%pff_tgt(d, lh)%seg, dim=1)) call die("mg:mmpg pff_tgt g>")
                              seg => cg%mg%pff_tgt(d, lh)%seg(g)
                              if (allocated(seg%buf)) then
                                 call warn("mg:mmpg o seg%buf a a")
                                 deallocate(seg%buf)
                              endif
                              seg%proc = j

                              ! find cross-section of own segment with coarsened fine face segment
                              if (ubound(curl%finer%pse(j)%c, dim=1) > 0) then
                                 coarsened(:, :) = curl%finer%pse(j)%c(1)%se(:,:)
                                 coarsened(d, hl) = coarsened(d, lh)
                                 coarsened(:, :) = coarsened(:, :)/2
                                 where (.not. dmask(:))
                                    seg%se(:, LO) = max(cg%my_se(:, LO), coarsened(:, LO))
                                    seg%se(:, HI) = min(cg%my_se(:, HI), coarsened(:, HI))
                                 endwhere
                                 seg%se(d, :) = coarsened(d, :)
                                 allocate(seg%buf(seg%se(xdim, HI)-seg%se(xdim, LO) + 1, &
                                      &           seg%se(ydim, HI)-seg%se(ydim, LO) + 1, &
                                      &           seg%se(zdim, HI)-seg%se(zdim, LO) + 1))

                                 coarsened(:, :) = curl%finer%pse(j)%c(1)%se(:,:)
                                 coarsened(d, hl) = coarsened(d, lh)
                                 coarsened(d, lh) = coarsened(d, lh) + 2*lh-LO-HI ! extend to two layers of buffer
                                 coarsened(:, :) = coarsened(:, :)/2
                                 coarsened(d, :) = coarsened(d, :) + [ -ord_prolong_face_norm, ord_prolong_face_norm ]

                                 seg%se(:, LO) = max(cg%my_se(:, LO), coarsened(:, LO)) + ijks(:)
                                 seg%se(:, HI) = min(cg%my_se(:, HI), coarsened(:, HI)) + ijks(:)

                                 coarsened(d, :) = coarsened(d, :) - [ -ord_prolong_face_norm, ord_prolong_face_norm ] ! revert broadening
                                 allocate(seg%f_lay(seg%se(d, HI) - seg%se(d, LO) + 1))
                                 do l = 1, size(seg%f_lay(:))
                                    seg%f_lay(l)%layer = l + int(seg%se(d, LO), kind=4) - 1
                                    nl = int(minval(abs(seg%f_lay(l)%layer - ijks(d) - coarsened(d, :))), kind=4)
                                    if (mod(curl%finer%pse(j)%c(1)%se(d, lh) + lh - LO, 2_LONG) == 0) then ! fine face at coarse face
                                       seg%f_lay(l)%coeff = opfn_c_ff(nl)
                                    else                                                              ! fine face at coarse center
                                       seg%f_lay(l)%coeff = opfn_c_cf(nl)
                                    endif
                                 enddo
                              endif

                           endif
                        enddo

                        cgl => cgl%nxt
                     enddo
                  endif

               endif
            enddo
         enddo

         curl => curl%finer
      enddo

   end subroutine mpi_multigrid_prep_grav

!>
!! \brief FFT given-boundary Poisson solver applied to local domain. Should require less communication than RBGS implementation.
!!
!! \todo test a configuration with wider area being subjected to FFT (sizes would no longer be 2**n) to avoid the need of relaxation
!<

   subroutine approximate_solution_fft(curl, src, soln)

      use cg_level_connected, only: cg_level_connected_T, coarsest
      use constants,     only: LO, HI, ndims, xdim, ydim, zdim, GEO_XYZ, half, I_ONE, idm2, BND_NEGREF, fft_none, fft_dst, dirtyL
      use dataio_pub,    only: die, warn
      use domain,        only: dom
      use cg_list,       only: cg_list_element
      use global,        only: dirty_debug
      use grid_cont,     only: grid_container
      use multigridvars, only: single_base, fft_full_relax, multidim_code_3D, nsmool, nsmoof
      use named_array,   only: p3

      implicit none

      type(cg_level_connected_T), pointer, intent(in) :: curl !< pointer to a level for which we approximate the solution
      integer,                      intent(in) :: src  !< index of source in cg%q(:)
      integer,                      intent(in) :: soln !< index of solution in cg%q(:)

      integer :: nf, n, nsmoo
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      integer, dimension(ndims,LO:HI) :: pdn, D_2
      integer                         :: dir

      if (associated(curl%first)) then
         if (associated(curl%first%nxt)) call warn("[multigrid_fft_approximation:approximate_solution_fft] This routine will likely fail due to numerous incompatibilities with newest grid features")
      endif

      if (curl%fft_type == fft_none) call die("[multigrid_fft_approximation:approximate_solution_fft] unknown FFT type")

      if (dom%geometry_type /= GEO_XYZ) call die("[multigrid_fft_approximation:approximate_solution_fft] FFT is not allowed in non-cartesian coordinates.")

      do nf = 1, nsmoof
         cgl => curl%first
         do while (associated(cgl))
            cg => cgl%cg ; p3 => cg%q(src)%span(cg%ijkse)
            cg%mg%src(:, :, :) = p3
            cgl => cgl%nxt
         enddo

         if (curl%fft_type == fft_dst) then !correct boundaries on non-periodic local domain
            if (nf == 1 .and. .not. associated(curl, coarsest)) then
               call make_face_boundaries(curl, soln)
            else
               call curl%arr3d_boundaries(soln, nb = I_ONE, bnd_type = BND_NEGREF)
               cgl => curl%first
               do while (associated(cgl))
                  cg => cgl%cg
                  if (dom%has_dir(xdim)) then
                     cg%mg%bnd_x(:, :, LO) = half* sum (cg%q(soln)%arr(cg%is-1:cg%is, cg%js:cg%je, cg%ks:cg%ke), 1)
                     cg%mg%bnd_x(:, :, HI) = half* sum (cg%q(soln)%arr(cg%ie:cg%ie+1, cg%js:cg%je, cg%ks:cg%ke), 1)
                  endif
                  if (dom%has_dir(ydim)) then
                     cg%mg%bnd_y(:, :, LO) = half* sum (cg%q(soln)%arr(cg%is:cg%ie, cg%js-1:cg%js, cg%ks:cg%ke), 2)
                     cg%mg%bnd_y(:, :, HI) = half* sum (cg%q(soln)%arr(cg%is:cg%ie, cg%je:cg%je+1, cg%ks:cg%ke), 2)
                  endif
                  if (dom%has_dir(zdim)) then
                     cg%mg%bnd_z(:, :, LO) = half* sum (cg%q(soln)%arr(cg%is:cg%ie, cg%js:cg%je, cg%ks-1:cg%ks), 3)
                     cg%mg%bnd_z(:, :, HI) = half* sum (cg%q(soln)%arr(cg%is:cg%ie, cg%js:cg%je, cg%ke:cg%ke+1), 3)
                  endif
                  cgl => cgl%nxt
               enddo
            endif

            cgl => curl%first
            do while (associated(cgl))
               cg => cgl%cg
               if (dirty_debug) then
                  if (dom%has_dir(xdim) .and. any(abs(curl%first%cg%mg%bnd_x(:, :, :)) > dirtyL)) call warn("approximate_solution_fft dirty bnd_x")
                  if (dom%has_dir(ydim) .and. any(abs(curl%first%cg%mg%bnd_y(:, :, :)) > dirtyL)) call warn("approximate_solution_fft dirty bnd_y")
                  if (dom%has_dir(zdim) .and. any(abs(curl%first%cg%mg%bnd_z(:, :, :)) > dirtyL)) call warn("approximate_solution_fft dirty bnd_z")
               endif

               if (dom%has_dir(xdim)) then
                  cg%mg%src(1,      :, :) = cg%mg%src(1,      :, :) - cg%mg%bnd_x(:, :, LO) * 2. * cg%idx2
                  cg%mg%src(cg%nxb, :, :) = cg%mg%src(cg%nxb, :, :) - cg%mg%bnd_x(:, :, HI) * 2. * cg%idx2
               endif
               if (dom%has_dir(ydim)) then
                  cg%mg%src(:, 1,      :) = cg%mg%src(:, 1,      :) - cg%mg%bnd_y(:, :, LO) * 2. * cg%idy2
                  cg%mg%src(:, cg%nyb, :) = cg%mg%src(:, cg%nyb, :) - cg%mg%bnd_y(:, :, HI) * 2. * cg%idy2
               endif
               if (dom%has_dir(zdim)) then
                  cg%mg%src(:, :, 1     ) = cg%mg%src(:, :, 1     ) - cg%mg%bnd_z(:, :, LO) * 2. * cg%idz2
                  cg%mg%src(:, :, cg%nzb) = cg%mg%src(:, :, cg%nzb) - cg%mg%bnd_z(:, :, HI) * 2. * cg%idz2
               endif
               cgl => cgl%nxt
            enddo
         endif

         call fft_convolve(curl)

         cgl => curl%first
         do while (associated(cgl))
            cg => cgl%cg
            p3 => cg%q(soln)%span(cg%ijkse) ; p3 = cg%mg%src(:, :, :)
            cgl => cgl%nxt
         enddo

         call curl%check_dirty(soln, "approx_soln fft+")

         !> \deprecated BEWARE use dom%has_dir() here in a way that does not degrade performance

         if (associated(curl, coarsest) .and. single_base) then !> do not relax if it is a whole domain (BEWARE: simplified check)
            nsmoo = 0
         else
            nsmoo = nsmool
         endif

         !relax the boundaries
         do n = 1, nsmoo
            call curl%arr3d_boundaries(soln, nb = I_ONE, bnd_type = BND_NEGREF)
            ! Possible optimization: This is a quite costly part of the local FFT solver
            cgl => curl%first
            do while (associated(cgl))
               cg => cgl%cg
               if (fft_full_relax) then
                  p3 => cg%q(soln)%span(cg%ijkse)
                  if (dom%eff_dim == ndims .and. .not. multidim_code_3D) then
                     p3 = cg%mg%rx * (cg%q(soln)%span(cg%ijkse-idm2(xdim,:,:)) + cg%q(soln)%span(cg%ijkse+idm2(xdim,:,:))) + &
                          cg%mg%ry * (cg%q(soln)%span(cg%ijkse-idm2(ydim,:,:)) + cg%q(soln)%span(cg%ijkse+idm2(ydim,:,:))) + &
                          cg%mg%rz * (cg%q(soln)%span(cg%ijkse-idm2(zdim,:,:)) + cg%q(soln)%span(cg%ijkse+idm2(zdim,:,:))) - &
                          cg%mg%r  *  cg%q(src)%span (cg%ijkse)
                  else

                     call die("[multigrid_fft_approximation:approximate_solution_fft] fft_full_relax is allowed only for 3D at the moment")

                     ! An additional array (cg%mg%src would be good enough) is required here to assemble partial results or use red-black passes
                     p3 = - cg%mg%r * cg%q(src)%span(cg%ijkse)
                     if (dom%has_dir(xdim)) p3 = p3 + cg%mg%rx * (cg%q(soln)%span(cg%ijkse-idm2(xdim,:,:)) + cg%q(soln)%span(cg%ijkse+idm2(xdim,:,:)))
                     if (dom%has_dir(ydim)) p3 = p3 + cg%mg%ry * (cg%q(soln)%span(cg%ijkse-idm2(ydim,:,:)) + cg%q(soln)%span(cg%ijkse+idm2(ydim,:,:)))
                     if (dom%has_dir(zdim)) p3 = p3 + cg%mg%rz * (cg%q(soln)%span(cg%ijkse-idm2(zdim,:,:)) + cg%q(soln)%span(cg%ijkse+idm2(zdim,:,:)))
                  endif
               else
                  ! relax only two layers of cells (1 is  significantly worse, 3 does not improve much)
                  ! edges are relaxed twice, corners are relaxed three times which seems to be good

                  D_2 = spread(dom%D_,2,2)
                  do dir = xdim, zdim
                     if (dom%has_dir(dir)) then
                        pdn = cg%ijkse ; pdn(dir,HI) = pdn(dir,HI) + dom%D_(dir) ! -X/-Y/-Z
                        p3 => cg%q(soln)%span(pdn)
                        p3 = cg%mg%rx * (cg%q(soln)%span(pdn-idm2(xdim,:,:)*D_2) + cg%q(soln)%span(pdn+idm2(xdim,:,:)*D_2)) + &
                             cg%mg%ry * (cg%q(soln)%span(pdn-idm2(ydim,:,:)*D_2) + cg%q(soln)%span(pdn+idm2(ydim,:,:)*D_2)) + &
                             cg%mg%rz * (cg%q(soln)%span(pdn-idm2(zdim,:,:)*D_2) + cg%q(soln)%span(pdn+idm2(zdim,:,:)*D_2)) - cg%mg%r  *  cg%q(src)%span(pdn)

                        pdn = cg%ijkse ; pdn(dir,HI) = pdn(dir,LO) - dom%D_(dir) ! +X/+Y/+Z
                        p3 => cg%q(soln)%span(pdn)
                        p3 = cg%mg%rx * (cg%q(soln)%span(pdn-idm2(xdim,:,:)*D_2) + cg%q(soln)%span(pdn+idm2(xdim,:,:)*D_2)) + &
                             cg%mg%ry * (cg%q(soln)%span(pdn-idm2(ydim,:,:)*D_2) + cg%q(soln)%span(pdn+idm2(ydim,:,:)*D_2)) + &
                             cg%mg%rz * (cg%q(soln)%span(pdn-idm2(zdim,:,:)*D_2) + cg%q(soln)%span(pdn+idm2(zdim,:,:)*D_2)) - cg%mg%r  *  cg%q(src)%span(pdn)
                     endif
                  enddo

               endif

               cgl => cgl%nxt
            enddo
         enddo

         call curl%check_dirty(soln, "approx_soln relax+")

      enddo

   end subroutine approximate_solution_fft

!> \brief This routine prepares boundary values for local-FFT solver

   subroutine make_face_boundaries(curl, soln)

      use cg_level_connected, only: cg_level_connected_T, coarsest
      use dataio_pub,    only: warn
      use mpisetup,      only: nproc
      use multigridvars, only: single_base, bnd_periodic, bnd_givenval, grav_bnd

      implicit none

      type(cg_level_connected_T), pointer, intent(in) :: curl !< pointer to a level for which we approximate the solution
      integer,                      intent(in) :: soln !< index of solution in cg%q(:)

      if (grav_bnd == bnd_periodic .and. (nproc == 1 .or. (associated(curl, coarsest) .and. single_base) ) ) then
         call curl%reset_boundaries
      else
         if (.not. associated(curl, coarsest)) then
            call prolong_faces(curl, soln)
         else
            if (grav_bnd /= bnd_givenval) call curl%reset_boundaries
            call warn("m:mfb WTF?")
         endif
      endif

   end subroutine make_face_boundaries

!> \brief Do the FFT convolution

   subroutine fft_convolve(curl)

      use constants,   only: fft_rcr, fft_dst
      use dataio_pub,  only: die
      use cg_list,     only: cg_list_element
      use cg_level_connected, only: cg_level_connected_T
      use grid_cont,   only: grid_container

      implicit none

      type(cg_level_connected_T), pointer, intent(in) :: curl !< pointer to a level at which make the convolution

      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      cgl => curl%first
      do while (associated(cgl))
         cg => cgl%cg

         ! do the convolution in Fourier space; cg%mg%src(:,:,:) -> cg%mg%fft{r}(:,:,:)
         call dfftw_execute(cg%mg%planf)

         select case (curl%fft_type)
            case (fft_rcr)
               cg%mg%fft  = cg%mg%fft  * cg%mg%Green3D
            case (fft_dst)
               cg%mg%fftr = cg%mg%fftr * cg%mg%Green3D
            case default
               call die("[multigrid_gravity:fft_convolve] Unknown FFT type.")
         end select

         call dfftw_execute(cg%mg%plani) ! cg%mg%fft{r}(:,:,:) -> cg%mg%src(:,:,:)
         cgl => cgl%nxt
      enddo

   end subroutine fft_convolve

!>
!! \brief Prolong solution data at level (lev-1) to faces at level lev
!!
!! \details It looks that parallel prolongation order 0 is the best from range -2 .. 2. Quite good is also order 2.(tests were performed with the maclaurin problem at r4051).
!! There were no experiments with higher order prolongation, but it seems there is not much room for any improvements.
!! The tests were performed on uniform grid, where the interpolation affects only convergence factors.
!! On a refinement step the interpolation type affects also sulution by influencing the way how fine and coarse grids are coupled.
!!
!! \todo move to cg_level_connected_T or multigrid_gravity
!!
!! OPT: completely unoptimized
!<

   subroutine prolong_faces(fine, soln)

      use cg_level_connected, only: cg_level_connected_T, coarsest
      use constants,     only: xdim, ydim, zdim, LO, HI, LONG, I_ONE, half, O_INJ, O_LIN, O_D2, O_I2, BND_NEGREF
      use dataio_pub,    only: die, warn
      use domain,        only: dom, is_multicg
      use grid_cont,     only: pr_segment
      use mpi,           only: MPI_DOUBLE_PRECISION
      use mpisetup,      only: comm, mpi_err, req, status, master, inflate_req
      use multigridvars, only: ord_prolong_face_norm, ord_prolong_face_par, need_general_pf

      implicit none

      type(cg_level_connected_T), pointer, intent(in) :: fine !< level for which approximate the solution
      integer,                      intent(in) :: soln !< index of solution in cg%q(:) ! \todo change the name

      integer                       :: i, j, k, d, lh, g, ib, jb, ibh, jbh, l
      type(cg_level_connected_T), pointer           :: coarse
      integer, parameter            :: s_wdth  = 3           ! interpolation stencil width
      integer(kind=4), parameter    :: s_rng = (s_wdth-1)/2  ! stencil range around 0
      real, parameter, dimension(s_wdth) :: p0  = [ 0.,       1.,     0.     ] ! injection
      real, parameter, dimension(s_wdth) :: p1  = [ 0.,       3./4.,  1./4.  ] ! 1D linear prolongation stencil
      real, parameter, dimension(s_wdth) :: p2i = [ -1./8.,   1.,     1./8.  ] ! 1D integral cubic prolongation stencil
      real, parameter, dimension(s_wdth) :: p2d = [ -3./32., 30./32., 5./32. ] ! 1D direct cubic prolongation stencil
      real, dimension(-s_rng:s_rng)                    :: p
      real, dimension(-s_rng:s_rng, -s_rng:s_rng, 2, 2):: pp   ! 2D prolongation stencil
      real                          :: ipp_norm

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

      if (associated(fine, coarsest)) then
         call warn("[multigrid_fft_approximation:prolong_faces] Cannot prolong anything to coarsest level")
         return
      endif
      if (is_multicg) call die("[multigrid_fft_approximation:prolong_faces] multicg not implemented") ! fine%first%cg%

      if (.not. associated(fine)) call die("[multigrid_fft_approximation:prolong_faces] fine == null()")
      coarse => fine%coarser
      if (.not. associated(coarse)) call die("[multigrid_fft_approximation:prolong_faces] coarse == null()")

      if (need_general_pf) then

         if (ord_prolong_face_par /= O_INJ) then
            ord_prolong_face_par = O_INJ
            if (master) call warn("[multigrid_fft_approximation:prolong_faces] only injection is supported for the current domain decomposition type.")
            ! The tests made with comm3d suggests that paralel interpolation only degrades convergence rate.
            ! Implement ord_prolong_face_par /= O_INJ if and only if it improves the coupling between fine and coarse solutions
         endif

         do lh = LO, HI ! \todo convert cg_level_connected_T%mg%bnd_[xyz] to an array and make obsolete the following pointer assignments
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
                     if (nr > size(req, dim=1)) call inflate_req
                     call MPI_Irecv(fine%first%cg%mg%pfc_tgt(d, lh)%seg(g)%buf(1, 1, 1), size(fine%first%cg%mg%pfc_tgt(d, lh)%seg(g)%buf(:, :, :)), MPI_DOUBLE_PRECISION, &
                          &         fine%first%cg%mg%pfc_tgt(d, lh)%seg(g)%proc, HI*d+lh, comm, req(nr), mpi_err)
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
                        if (nr > size(req, dim=1)) call inflate_req
                        se(:,:) = cse(:,:)
                        pseg%buf(:, :, :) = 0. ! this can be avoided by extracting first assignment from the loop
                        do l = lbound(pseg%f_lay(:), dim=1), ubound(pseg%f_lay(:), dim=1)
                           se(d,:) = pseg%f_lay(l)%layer
                           pseg%buf(:, :, :) = pseg%buf(:, :, :) + pseg%f_lay(l)%coeff * coarse%first%cg%q(soln)%span(se)
                        enddo
                        call MPI_Isend(pseg%buf(1, 1, 1), size(pseg%buf(:, :, :)), MPI_DOUBLE_PRECISION, pseg%proc, HI*d+lh, comm, req(nr), mpi_err)
                     enddo
                  endif
               endif
            enddo
         enddo

         if (nr>0) call MPI_Waitall(nr, req(:nr), status(:,:nr), mpi_err)

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
               call die("[multigrid_fft_approximation:prolong_faces] invalid ord_prolong_face_par")
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
         call coarse%arr3d_boundaries(soln, nb = b_rng, bnd_type = BND_NEGREF, corners = (ord_prolong_face_par/=0)) !> \deprecated BEWARE for higher prolongation order more guardcell are required
         call coarse%check_dirty(soln, "prolong_faces", s_rng)
         associate ( &
            ffcg => fine%first%cg, &
            cfcg => coarse%first%cg, &
            coarse_soln => coarse%first%cg%q(soln)%arr, &
            bnd_x => fine%first%cg%mg%bnd_x, &
            bnd_y => fine%first%cg%mg%bnd_y, &
            bnd_z => fine%first%cg%mg%bnd_z, &
            D_x => dom%D_x, &
            D_y => dom%D_y, &
            D_z => dom%D_z &
         )
         if (ord_prolong_face_norm /= O_INJ) then
            if (dom%has_dir(xdim)) then
               !pp_norm = 2.*sum(pp(-D_y:D_y, -D_z:D_z, 1, 1)) ! normalization is required for ord_prolong_face_par == 1 and -2
               ipp_norm = 0.5 / sum(pp(-D_y:D_y, -D_z:D_z, 1, 1))
               do j = cfcg%js, cfcg%je
                  do k = cfcg%ks, cfcg%ke
                     bnd_x(-ffcg%js+2*j,        -ffcg%ks+2*k,        LO)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 1) * ( &
                          opfn1*(coarse_soln(cfcg%is,  j-D_y:j+D_y,k-D_z:k+D_z) + coarse_soln(cfcg%is-1,j-D_y:j+D_y,k-D_z:k+D_z)) + &
                          opfn3*(coarse_soln(cfcg%is+1,j-D_y:j+D_y,k-D_z:k+D_z) + coarse_soln(cfcg%is-2,j-D_y:j+D_y,k-D_z:k+D_z)))) * ipp_norm
                     bnd_x(-ffcg%js+2*j+D_y,-ffcg%ks+2*k,        LO)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 1) * ( &
                          opfn1*(coarse_soln(cfcg%is,  j-D_y:j+D_y,k-D_z:k+D_z) + coarse_soln(cfcg%is-1,j-D_y:j+D_y,k-D_z:k+D_z)) + &
                          opfn3*(coarse_soln(cfcg%is+1,j-D_y:j+D_y,k-D_z:k+D_z) + coarse_soln(cfcg%is-2,j-D_y:j+D_y,k-D_z:k+D_z)))) * ipp_norm
                     bnd_x(-ffcg%js+2*j,        -ffcg%ks+2*k+D_z,LO)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 2) * ( &
                          opfn1*(coarse_soln(cfcg%is  ,j-D_y:j+D_y,k-D_z:k+D_z) + coarse_soln(cfcg%is-1,j-D_y:j+D_y,k-D_z:k+D_z)) + &
                          opfn3*(coarse_soln(cfcg%is+1,j-D_y:j+D_y,k-D_z:k+D_z) + coarse_soln(cfcg%is-2,j-D_y:j+D_y,k-D_z:k+D_z)))) * ipp_norm
                     bnd_x(-ffcg%js+2*j+D_y,-ffcg%ks+2*k+D_z,LO)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 2) * ( &
                          opfn1*(coarse_soln(cfcg%is,  j-D_y:j+D_y,k-D_z:k+D_z) + coarse_soln(cfcg%is-1,j-D_y:j+D_y,k-D_z:k+D_z)) + &
                          opfn3*(coarse_soln(cfcg%is+1,j-D_y:j+D_y,k-D_z:k+D_z) + coarse_soln(cfcg%is-2,j-D_y:j+D_y,k-D_z:k+D_z)))) * ipp_norm
                     bnd_x(-ffcg%js+2*j,        -ffcg%ks+2*k,        HI)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 1) * ( &
                          opfn1*(coarse_soln(cfcg%ie,  j-D_y:j+D_y,k-D_z:k+D_z) + coarse_soln(cfcg%ie+1,j-D_y:j+D_y,k-D_z:k+D_z)) + &
                          opfn3*(coarse_soln(cfcg%ie-1,j-D_y:j+D_y,k-D_z:k+D_z) + coarse_soln(cfcg%ie+2,j-D_y:j+D_y,k-D_z:k+D_z)))) * ipp_norm
                     bnd_x(-ffcg%js+2*j+D_y,-ffcg%ks+2*k,        HI)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 1) * ( &
                          opfn1*(coarse_soln(cfcg%ie,  j-D_y:j+D_y,k-D_z:k+D_z) + coarse_soln(cfcg%ie+1,j-D_y:j+D_y,k-D_z:k+D_z)) + &
                          opfn3*(coarse_soln(cfcg%ie-1,j-D_y:j+D_y,k-D_z:k+D_z) + coarse_soln(cfcg%ie+2,j-D_y:j+D_y,k-D_z:k+D_z)))) * ipp_norm
                     bnd_x(-ffcg%js+2*j,        -ffcg%ks+2*k+D_z,HI)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 2) * ( &
                          opfn1*(coarse_soln(cfcg%ie,  j-D_y:j+D_y,k-D_z:k+D_z) + coarse_soln(cfcg%ie+1,j-D_y:j+D_y,k-D_z:k+D_z)) + &
                          opfn3*(coarse_soln(cfcg%ie-1,j-D_y:j+D_y,k-D_z:k+D_z) + coarse_soln(cfcg%ie+2,j-D_y:j+D_y,k-D_z:k+D_z)))) * ipp_norm
                     bnd_x(-ffcg%js+2*j+D_y,-ffcg%ks+2*k+D_z,HI)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 2) * ( &
                          opfn1*(coarse_soln(cfcg%ie,  j-D_y:j+D_y,k-D_z:k+D_z) + coarse_soln(cfcg%ie+1,j-D_y:j+D_y,k-D_z:k+D_z)) + &
                          opfn3*(coarse_soln(cfcg%ie-1,j-D_y:j+D_y,k-D_z:k+D_z) + coarse_soln(cfcg%ie+2,j-D_y:j+D_y,k-D_z:k+D_z)))) * ipp_norm
                  enddo
               enddo
            endif

            if (dom%has_dir(ydim)) then
               !pp_norm = 2.*sum(pp(-D_x:D_x, -D_z:D_z, 1, 1))
               ipp_norm = 0.5 / sum(pp(-D_x:D_x, -D_z:D_z, 1, 1))
               do i = cfcg%is, cfcg%ie
                  do k = cfcg%ks, cfcg%ke
                     bnd_y(-ffcg%is+2*i,        -ffcg%ks+2*k,        LO)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 1) * ( &
                          opfn1*(coarse_soln(i-D_x:i+D_x,cfcg%js,  k-D_z:k+D_z) + coarse_soln(i-D_x:i+D_x,cfcg%js-1,k-D_z:k+D_z)) + &
                          opfn3*(coarse_soln(i-D_x:i+D_x,cfcg%js+1,k-D_z:k+D_z) + coarse_soln(i-D_x:i+D_x,cfcg%js-2,k-D_z:k+D_z)))) * ipp_norm
                     bnd_y(-ffcg%is+2*i+D_x,-ffcg%ks+2*k,        LO)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 1) * ( &
                          opfn1*(coarse_soln(i-D_x:i+D_x,cfcg%js,  k-D_z:k+D_z) + coarse_soln(i-D_x:i+D_x,cfcg%js-1,k-D_z:k+D_z)) + &
                          opfn3*(coarse_soln(i-D_x:i+D_x,cfcg%js+1,k-D_z:k+D_z) + coarse_soln(i-D_x:i+D_x,cfcg%js-2,k-D_z:k+D_z)))) * ipp_norm
                     bnd_y(-ffcg%is+2*i,        -ffcg%ks+2*k+D_z,LO)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 2) * ( &
                          opfn1*(coarse_soln(i-D_x:i+D_x,cfcg%js,  k-D_z:k+D_z) + coarse_soln(i-D_x:i+D_x,cfcg%js-1,k-D_z:k+D_z)) + &
                          opfn3*(coarse_soln(i-D_x:i+D_x,cfcg%js+1,k-D_z:k+D_z) + coarse_soln(i-D_x:i+D_x,cfcg%js-2,k-D_z:k+D_z)))) * ipp_norm
                     bnd_y(-ffcg%is+2*i+D_x,-ffcg%ks+2*k+D_z,LO)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 2) * ( &
                          opfn1*(coarse_soln(i-D_x:i+D_x,cfcg%js,  k-D_z:k+D_z) + coarse_soln(i-D_x:i+D_x,cfcg%js-1,k-D_z:k+D_z)) + &
                          opfn3*(coarse_soln(i-D_x:i+D_x,cfcg%js+1,k-D_z:k+D_z) + coarse_soln(i-D_x:i+D_x,cfcg%js-2,k-D_z:k+D_z)))) * ipp_norm
                     bnd_y(-ffcg%is+2*i,        -ffcg%ks+2*k,        HI)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 1) * ( &
                          opfn1*(coarse_soln(i-D_x:i+D_x,cfcg%je,  k-D_z:k+D_z) + coarse_soln(i-D_x:i+D_x,cfcg%je+1,k-D_z:k+D_z)) + &
                          opfn3*(coarse_soln(i-D_x:i+D_x,cfcg%je-1,k-D_z:k+D_z) + coarse_soln(i-D_x:i+D_x,cfcg%je+2,k-D_z:k+D_z)))) * ipp_norm
                     bnd_y(-ffcg%is+2*i+D_x,-ffcg%ks+2*k,        HI)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 1) * ( &
                          opfn1*(coarse_soln(i-D_x:i+D_x,cfcg%je,  k-D_z:k+D_z) + coarse_soln(i-D_x:i+D_x,cfcg%je+1,k-D_z:k+D_z)) + &
                          opfn3*(coarse_soln(i-D_x:i+D_x,cfcg%je-1,k-D_z:k+D_z) + coarse_soln(i-D_x:i+D_x,cfcg%je+2,k-D_z:k+D_z)))) * ipp_norm
                     bnd_y(-ffcg%is+2*i,        -ffcg%ks+2*k+D_z,HI)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 2) * ( &
                          opfn1*(coarse_soln(i-D_x:i+D_x,cfcg%je,  k-D_z:k+D_z) + coarse_soln(i-D_x:i+D_x,cfcg%je+1,k-D_z:k+D_z)) + &
                          opfn3*(coarse_soln(i-D_x:i+D_x,cfcg%je-1,k-D_z:k+D_z) + coarse_soln(i-D_x:i+D_x,cfcg%je+2,k-D_z:k+D_z)))) * ipp_norm
                     bnd_y(-ffcg%is+2*i+D_x,-ffcg%ks+2*k+D_z,HI)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 2) * ( &
                          opfn1*(coarse_soln(i-D_x:i+D_x,cfcg%je,  k-D_z:k+D_z) + coarse_soln(i-D_x:i+D_x,cfcg%je+1,k-D_z:k+D_z)) + &
                          opfn3*(coarse_soln(i-D_x:i+D_x,cfcg%je-1,k-D_z:k+D_z) + coarse_soln(i-D_x:i+D_x,cfcg%je+2,k-D_z:k+D_z)))) * ipp_norm
                  enddo
               enddo
            endif

            if (dom%has_dir(zdim)) then
               !pp_norm = 2.*sum(pp(-D_x:D_x, -D_y:D_y, 1, 1))
               ipp_norm = 0.5 / sum(pp(-D_x:D_x, -D_y:D_y, 1, 1))
               do i = cfcg%is, cfcg%ie
                  do j = cfcg%js, cfcg%je
                     bnd_z(-ffcg%is+2*i,        -ffcg%js+2*j,        LO)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 1) * ( &
                          opfn1*(coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ks) + coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ks-1)) + &
                          opfn3*(coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ks+1) + coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ks-2)))) * ipp_norm
                     bnd_z(-ffcg%is+2*i+D_x,-ffcg%js+2*j,        LO)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 1) * ( &
                          opfn1*(coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ks) + coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ks-1)) + &
                          opfn3*(coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ks+1) + coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ks-2)))) * ipp_norm
                     bnd_z(-ffcg%is+2*i,        -ffcg%js+2*j+D_y,LO)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 2) * ( &
                          opfn1*(coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ks) + coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ks-1)) + &
                          opfn3*(coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ks+1) + coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ks-2)))) * ipp_norm
                     bnd_z(-ffcg%is+2*i+D_x,-ffcg%js+2*j+D_y,LO)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 2) * ( &
                          opfn1*(coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ks) + coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ks-1)) + &
                          opfn3*(coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ks+1) + coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ks-2)))) * ipp_norm
                     bnd_z(-ffcg%is+2*i,        -ffcg%js+2*j,        HI)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 1) * ( &
                          opfn1*(coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ke) + coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ke+1)) + &
                          opfn3*(coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ke-1) + coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ke+2)))) * ipp_norm
                     bnd_z(-ffcg%is+2*i+D_x,-ffcg%js+2*j,        HI)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 1) * ( &
                          opfn1*(coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ke) + coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ke+1)) + &
                          opfn3*(coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ke-1) + coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ke+2)))) * ipp_norm
                     bnd_z(-ffcg%is+2*i,        -ffcg%js+2*j+D_y,HI)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 2) * ( &
                          opfn1*(coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ke) + coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ke+1)) + &
                          opfn3*(coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ke-1) + coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ke+2)))) * ipp_norm
                     bnd_z(-ffcg%is+2*i+D_x,-ffcg%js+2*j+D_y,HI)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 2) * ( &
                          opfn1*(coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ke) + coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ke+1)) + &
                          opfn3*(coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ke-1) + coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ke+2)))) * ipp_norm
                  enddo
               enddo
            endif

         else

            if (dom%has_dir(xdim)) then
               !pp_norm = 2.*sum(pp(-D_y:D_y, -D_z:D_z, 1, 1)) ! normalization is required for ord_prolong_face_par == 1 and -2
               ipp_norm = 0.5 / sum(pp(-D_y:D_y, -D_z:D_z, 1, 1))
               do j = cfcg%js, cfcg%je
                  do k = cfcg%ks, cfcg%ke
                     bnd_x(-ffcg%js+2*j,        -ffcg%ks+2*k,        LO)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 1) * &
                        & (coarse_soln(cfcg%is,j-D_y:j+D_y,k-D_z:k+D_z) + coarse_soln(cfcg%is-1,j-D_y:j+D_y,k-D_z:k+D_z))) * ipp_norm
                     bnd_x(-ffcg%js+2*j+D_y,-ffcg%ks+2*k,        LO)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 1) * &
                        & (coarse_soln(cfcg%is,j-D_y:j+D_y,k-D_z:k+D_z) + coarse_soln(cfcg%is-1,j-D_y:j+D_y,k-D_z:k+D_z))) * ipp_norm
                     bnd_x(-ffcg%js+2*j,        -ffcg%ks+2*k+D_z,LO)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 2) * &
                        & (coarse_soln(cfcg%is,j-D_y:j+D_y,k-D_z:k+D_z) + coarse_soln(cfcg%is-1,j-D_y:j+D_y,k-D_z:k+D_z))) * ipp_norm
                     bnd_x(-ffcg%js+2*j+D_y,-ffcg%ks+2*k+D_z,LO)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 2) * &
                        & (coarse_soln(cfcg%is,j-D_y:j+D_y,k-D_z:k+D_z) + coarse_soln(cfcg%is-1,j-D_y:j+D_y,k-D_z:k+D_z))) * ipp_norm
                     bnd_x(-ffcg%js+2*j,        -ffcg%ks+2*k,        HI)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 1) * &
                        & (coarse_soln(cfcg%ie,j-D_y:j+D_y,k-D_z:k+D_z) + coarse_soln(cfcg%ie+1,j-D_y:j+D_y,k-D_z:k+D_z))) * ipp_norm
                     bnd_x(-ffcg%js+2*j+D_y,-ffcg%ks+2*k,        HI)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 1) * &
                        & (coarse_soln(cfcg%ie,j-D_y:j+D_y,k-D_z:k+D_z) + coarse_soln(cfcg%ie+1,j-D_y:j+D_y,k-D_z:k+D_z))) * ipp_norm
                     bnd_x(-ffcg%js+2*j,        -ffcg%ks+2*k+D_z,HI)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 2) * &
                        & (coarse_soln(cfcg%ie,j-D_y:j+D_y,k-D_z:k+D_z) + coarse_soln(cfcg%ie+1,j-D_y:j+D_y,k-D_z:k+D_z))) * ipp_norm
                     bnd_x(-ffcg%js+2*j+D_y,-ffcg%ks+2*k+D_z,HI)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 2) * &
                        & (coarse_soln(cfcg%ie,j-D_y:j+D_y,k-D_z:k+D_z) + coarse_soln(cfcg%ie+1,j-D_y:j+D_y,k-D_z:k+D_z))) * ipp_norm
                  enddo
               enddo
            endif

            if (dom%has_dir(ydim)) then
               !pp_norm = 2.*sum(pp(-D_x:D_x, -D_z:D_z, 1, 1))
               ipp_norm = 0.5 / sum(pp(-D_x:D_x, -D_z:D_z, 1, 1))
               do i = cfcg%is, cfcg%ie
                  do k = cfcg%ks, cfcg%ke
                     bnd_y(-ffcg%is+2*i,        -ffcg%ks+2*k,        LO)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 1) * &
                        & (coarse_soln(i-D_x:i+D_x,cfcg%js,k-D_z:k+D_z) + coarse_soln(i-D_x:i+D_x,cfcg%js-1,k-D_z:k+D_z))) * ipp_norm
                     bnd_y(-ffcg%is+2*i+D_x,-ffcg%ks+2*k,        LO)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 1) * &
                        & (coarse_soln(i-D_x:i+D_x,cfcg%js,k-D_z:k+D_z) + coarse_soln(i-D_x:i+D_x,cfcg%js-1,k-D_z:k+D_z))) * ipp_norm
                     bnd_y(-ffcg%is+2*i,        -ffcg%ks+2*k+D_z,LO)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 2) * &
                        & (coarse_soln(i-D_x:i+D_x,cfcg%js,k-D_z:k+D_z) + coarse_soln(i-D_x:i+D_x,cfcg%js-1,k-D_z:k+D_z))) * ipp_norm
                     bnd_y(-ffcg%is+2*i+D_x,-ffcg%ks+2*k+D_z,LO)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 2) * &
                        & (coarse_soln(i-D_x:i+D_x,cfcg%js,k-D_z:k+D_z) + coarse_soln(i-D_x:i+D_x,cfcg%js-1,k-D_z:k+D_z))) * ipp_norm
                     bnd_y(-ffcg%is+2*i,        -ffcg%ks+2*k,        HI)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 1) * &
                        & (coarse_soln(i-D_x:i+D_x,cfcg%je,k-D_z:k+D_z) + coarse_soln(i-D_x:i+D_x,cfcg%je+1,k-D_z:k+D_z))) * ipp_norm
                     bnd_y(-ffcg%is+2*i+D_x,-ffcg%ks+2*k,        HI)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 1) * &
                        & (coarse_soln(i-D_x:i+D_x,cfcg%je,k-D_z:k+D_z) + coarse_soln(i-D_x:i+D_x,cfcg%je+1,k-D_z:k+D_z))) * ipp_norm
                     bnd_y(-ffcg%is+2*i,        -ffcg%ks+2*k+D_z,HI)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 2) * &
                        & (coarse_soln(i-D_x:i+D_x,cfcg%je,k-D_z:k+D_z) + coarse_soln(i-D_x:i+D_x,cfcg%je+1,k-D_z:k+D_z))) * ipp_norm
                     bnd_y(-ffcg%is+2*i+D_x,-ffcg%ks+2*k+D_z,HI)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 2) * &
                        & (coarse_soln(i-D_x:i+D_x,cfcg%je,k-D_z:k+D_z) + coarse_soln(i-D_x:i+D_x,cfcg%je+1,k-D_z:k+D_z))) * ipp_norm
                  enddo
               enddo
            endif

            if (dom%has_dir(zdim)) then
               !pp_norm = 2.*sum(pp(-D_x:D_x, -D_y:D_y, 1, 1))
               ipp_norm = 0.5 / sum(pp(-D_x:D_x, -D_y:D_y, 1, 1))
               do i = cfcg%is, cfcg%ie
                  do j = cfcg%js, cfcg%je
                     bnd_z(-ffcg%is+2*i,        -ffcg%js+2*j,        LO)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 1) * &
                        & (coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ks) + coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ks-1))) * ipp_norm
                     bnd_z(-ffcg%is+2*i+D_x,-ffcg%js+2*j,        LO)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 1) * &
                        & (coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ks) + coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ks-1))) * ipp_norm
                     bnd_z(-ffcg%is+2*i,        -ffcg%js+2*j+D_y,LO)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 2) * &
                        & (coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ks) + coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ks-1))) * ipp_norm
                     bnd_z(-ffcg%is+2*i+D_x,-ffcg%js+2*j+D_y,LO)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 2) * &
                        & (coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ks) + coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ks-1))) * ipp_norm
                     bnd_z(-ffcg%is+2*i,        -ffcg%js+2*j,        HI)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 1) * &
                        & (coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ke) + coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ke+1))) * ipp_norm
                     bnd_z(-ffcg%is+2*i+D_x,-ffcg%js+2*j,        HI)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 1) * &
                        & (coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ke) + coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ke+1))) * ipp_norm
                     bnd_z(-ffcg%is+2*i,        -ffcg%js+2*j+D_y,HI)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 2) * &
                        & (coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ke) + coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ke+1))) * ipp_norm
                     bnd_z(-ffcg%is+2*i+D_x,-ffcg%js+2*j+D_y,HI)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 2) * &
                        & (coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ke) + coarse_soln(i-D_x:i+D_x,j-D_y:j+D_y,cfcg%ke+1))) * ipp_norm
                  enddo
               enddo
            endif
         endif
         end associate
      endif

   end subroutine prolong_faces

end module multigrid_fftapprox
