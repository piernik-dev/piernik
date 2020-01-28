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

!>
!!\brief This module serves purpose to update single cells using CRESP routines over the whole grid
!<

module cresp_grid
! pulled by COSM_RAY_ELECTRONS

   implicit none

   private
   public :: cresp_update_grid, cresp_init_grid, cfl_cresp_violation, cresp_clean_grid

   logical :: cfl_cresp_violation

   contains

!----------------------------------------------------------------------------------------------------

   subroutine cresp_init_grid

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use cg_list_global,   only: all_cg
      use cresp_crspectrum, only: cresp_allocate_all, cresp_init_state, p_rch_init
      use cresp_NR_method,  only: cresp_initialize_guess_grids
      use dataio,           only: vars
      use dataio_pub,       only: printinfo
      use func,             only: emag
      use grid_cont,        only: grid_container
      use initcosmicrays,   only: iarr_cre_n, iarr_cre_e, ncre
      use initcrspectrum,   only: norm_init_spectrum, dfpq, check_if_dump_fpq
      use mpisetup,         only: master
      use named_array_list, only: wna

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      call cresp_initialize_guess_grids
      call cresp_allocate_all

      call check_if_dump_fpq(vars)

      if (dfpq%any_dump) then
         if (dfpq%dump_f) call all_cg%reg_var(dfpq%f_nam, dim4=ncre+1)
         if (dfpq%dump_p) call all_cg%reg_var(dfpq%p_nam, dim4=2)
         if (dfpq%dump_q) call all_cg%reg_var(dfpq%q_nam, dim4=ncre)
      endif

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         cg%u(iarr_cre_n,:,:,:)  = 0.0
         cg%u(iarr_cre_e,:,:,:)  = 0.0

         if (dfpq%any_dump) then
            if (dfpq%dump_f) cg%w(wna%ind(dfpq%f_nam))%arr = 0.0
            if (dfpq%dump_p) cg%w(wna%ind(dfpq%p_nam))%arr = 0.0
            if (dfpq%dump_q) cg%w(wna%ind(dfpq%q_nam))%arr = 0.0
         endif

         cgl => cgl%nxt
      enddo

      call p_rch_init               !< sets the right pointer for p_rch function, based on used Taylor expansion coefficient

      call cresp_init_state(norm_init_spectrum%n, norm_init_spectrum%e)   !< initialize spectrum here, f_init should be 1.0

      if (master) call printinfo(" [cresp_grid:cresp_init_grid] CRESP initialized")

   end subroutine cresp_init_grid

   subroutine cresp_update_grid

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, onet
      use cresp_crspectrum, only: cresp_update_cell
      use crhelpers,        only: divv_i
      use func,             only: emag
      use global,           only: dt
      use grid_cont,        only: grid_container
      use initcosmicrays,   only: iarr_cre_e, iarr_cre_n
      use initcrspectrum,   only: spec_mod_trms, synch_active, adiab_active, cresp, crel, dfpq, fsynchr, u_b_max
      use named_array_list, only: wna
#ifdef DEBUG
      use cresp_crspectrum, only: cresp_detect_negative_content
#endif /* DEBUG */

      implicit none

      integer                        :: i, j, k
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer  :: cg
      type(spec_mod_trms)            :: sptab

      cgl => leaves%first
      cfl_cresp_violation = .false.

      do while (associated(cgl))
         cg => cgl%cg
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  sptab%ud = 0.0 ; sptab%ub = 0.0 ; sptab%ucmb = 0.0
                  cresp%n = cg%u(iarr_cre_n, i, j, k)
                  cresp%e = cg%u(iarr_cre_e, i, j, k)
                  if (synch_active) sptab%ub = min(emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k)) * fsynchr, u_b_max)    !< WARNING assusmes that b is in mGs
                  if (adiab_active) sptab%ud = cg%q(divv_i)%point([i,j,k]) * onet
#ifdef CRESP_VERBOSED
                  print *, 'Output of cosmic ray electrons module for grid cell with coordinates i,j,k:', i, j, k
#endif /* CRESP_VERBOSED */
                  call cresp_update_cell(2 * dt, cresp%n, cresp%e, sptab, cfl_cresp_violation)
#ifdef DEBUG
                  call cresp_detect_negative_content([i, j, k])
#endif /* DEBUG */
                  if (cfl_cresp_violation) return ! nothing to do here!
                  cg%u(iarr_cre_n, i, j, k) = cresp%n
                  cg%u(iarr_cre_e, i, j, k) = cresp%e
                  if (dfpq%any_dump) then
                     if (dfpq%dump_f) cg%w(wna%ind(dfpq%f_nam))%arr(:, i, j, k) = crel%f
                     if (dfpq%dump_p) cg%w(wna%ind(dfpq%p_nam))%arr(:, i, j, k) = crel%p(crel%i_cut)
                     if (dfpq%dump_q) cg%w(wna%ind(dfpq%q_nam))%arr(:, i, j, k) = crel%q
                  endif
               enddo
            enddo
         enddo
         cgl=>cgl%nxt
      enddo

   end subroutine cresp_update_grid

!----------------------------------------------------------------------------------------------------

   subroutine cresp_clean_grid

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use cresp_crspectrum, only: detect_clean_spectrum
      use grid_cont,        only: grid_container
      use initcosmicrays,   only: iarr_cre_e, iarr_cre_n
      use initcrspectrum,   only: cresp, nullify_empty_bins

      implicit none

      integer                        :: i, j, k
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      if (.not.nullify_empty_bins) return

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  cresp%n = cg%u(iarr_cre_n, i, j, k)
                  cresp%e = cg%U(iarr_cre_e, i, j, k)

                  call detect_clean_spectrum

                  cg%u(iarr_cre_n, i, j, k) = cresp%n
                  cg%u(iarr_cre_e, i, j, k) = cresp%e
               enddo
            enddo
         enddo
         cgl=>cgl%nxt
      enddo

   end subroutine cresp_clean_grid

end module cresp_grid
