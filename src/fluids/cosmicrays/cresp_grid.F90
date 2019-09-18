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

   use global,         only: dt, t
   use initcosmicrays, only: iarr_cre_e, iarr_cre_n

   implicit none

   private
   public :: cresp_update_grid, cresp_init_grid, cfl_cresp_violation, cresp_clean_grid, fsynchr

   real(kind=8) :: fsynchr
   logical      :: cfl_cresp_violation, register_p, register_q, register_f

   contains

!----------------------------------------------------------------------------------------------------

   subroutine cresp_init_grid

      use cg_leaves,          only: leaves
      use cg_list,            only: cg_list_element
      use cg_list_global,     only: all_cg
      use cresp_crspectrum,   only: cresp_allocate_all, e_threshold_lo, e_threshold_up, fail_count_interpol, &
                              &     fail_count_NR_2dim, fail_count_comp_q, cresp_init_state, p_rch_init
      use cresp_NR_method,    only: cresp_initialize_guess_grids
      use dataio,             only: vars
      use dataio_pub,         only: warn, printinfo, msg
      use grid_cont,          only: grid_container
      use initcosmicrays,     only: iarr_cre_n, iarr_cre_e, ncre
      use initcrspectrum,     only: e_small, e_small_approx_p_lo, e_small_approx_p_up, norm_init_spectrum, f_init, &
                                    hdf_save_fpq, nam_cresp_f, nam_cresp_p, nam_cresp_q, check_if_dump_fpq, dump_f, dump_p, dump_q
      use mpisetup,           only: master
      use named_array_list,   only: wna
      use units,              only: clight, me, sigma_T

      implicit none

      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

      register_f = .false.
      register_p = .false.
      register_q = .false.

      call cresp_initialize_guess_grids
      call cresp_allocate_all

      fail_count_interpol = 0
      fail_count_NR_2dim  = 0
      fail_count_comp_q   = 0

      e_threshold_lo = e_small * e_small_approx_p_lo
      e_threshold_up = e_small * e_small_approx_p_up

      fsynchr =  (4. / 3. ) * sigma_T / (me * clight)
      write (msg, *) "[cresp_grid:cresp_init_grid] 4/3 * sigma_T / ( me * c ) = ", fsynchr
      if (master) call printinfo(msg)

      call check_if_dump_fpq(vars)

      if (hdf_save_fpq) then
         if (dump_f) call all_cg%reg_var(nam_cresp_f, dim4=ncre+1)
         if (dump_p) call all_cg%reg_var(nam_cresp_p, dim4=2)
         if (dump_q) call all_cg%reg_var(nam_cresp_q, dim4=ncre)
      endif

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         cg%u(iarr_cre_n,:,:,:)  = 0.0
         cg%u(iarr_cre_e,:,:,:)  = 0.0

         if (hdf_save_fpq) then
            if (dump_f) cg%w(wna%ind(nam_cresp_f))%arr = 0.0
            if (dump_p) cg%w(wna%ind(nam_cresp_p))%arr = 0.0
            if (dump_q) cg%w(wna%ind(nam_cresp_q))%arr = 0.0
         endif

         cgl => cgl%nxt
      enddo

      call p_rch_init               !< sets the right pointer for p_rch function, based on used Taylor expansion coefficient

      call cresp_init_state(norm_init_spectrum%n, norm_init_spectrum%e, f_init)   !< initialize spectrum here, f_init should be 1.0

      if (master) call printinfo(" [cresp_grid:cresp_init_grid] CRESP initialized")

   end subroutine cresp_init_grid

   subroutine cresp_update_grid

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, onet
      use cresp_crspectrum, only: cresp_update_cell
      use crhelpers,        only: divv_n
      use func,             only: emag
      use grid_cont,        only: grid_container
      use initcrspectrum,   only: spec_mod_trms, synch_active, adiab_active, cresp, hdf_save_fpq, crel, nam_cresp_f, nam_cresp_p, nam_cresp_q, dump_f, dump_p, dump_q
      use named_array,      only: p4
      use named_array_list, only: qna, wna
#ifdef DEBUG
      use cresp_crspectrum, only: cresp_detect_negative_content
#endif /* DEBUG */

      implicit none

      integer                        :: i, j, k
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer  :: cg
      type(spec_mod_trms)            :: sptab

      i = 0; j = 0;  k = 0
      cgl => leaves%first
      cfl_cresp_violation = .false.

      do while (associated(cgl))
         cg => cgl%cg
         p4 => cg%w(wna%fi)%arr
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  sptab%ud = 0.0 ; sptab%ub = 0.0 ; sptab%ucmb = 0.0
                  cresp%n    = p4(iarr_cre_n, i, j, k)
                  cresp%e    = p4(iarr_cre_e, i, j, k)
                  if (synch_active) sptab%ub = emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k)) * fsynchr    !< WARNING assusmes that b is in mGs
                  if (adiab_active) sptab%ud = cg%q(qna%ind(divv_n))%point([i,j,k]) * onet
#ifdef CRESP_VERBOSED
                  print *, 'Output of cosmic ray electrons module for grid cell with coordinates i,j,k:', i, j, k
#endif /* CRESP_VERBOSED */
                  call cresp_update_cell(2 * dt, cresp%n, cresp%e, sptab, cfl_cresp_violation)
#ifdef DEBUG
                  call cresp_detect_negative_content([i, j, k])
#endif /* DEBUG */
                  if ( cfl_cresp_violation ) return ! nothing to do here!
                  p4(iarr_cre_n, i, j, k) = cresp%n
                  p4(iarr_cre_e, i, j, k) = cresp%e
                  if (hdf_save_fpq) then
                     if (dump_f) cg%w(wna%ind(nam_cresp_f))%arr(:, i, j, k) = crel%f
                     if (dump_p) cg%w(wna%ind(nam_cresp_p))%arr(:, i, j, k) = [crel%p(crel%i_lo), crel%p(crel%i_up)]
                     if (dump_q) cg%w(wna%ind(nam_cresp_q))%arr(:, i, j, k) = crel%q
                  endif
               enddo
            enddo
         enddo
         cg%u(iarr_cre_n, :,:,:) = p4(iarr_cre_n, :,:,:)
         cg%u(iarr_cre_e, :,:,:) = p4(iarr_cre_e, :,:,:)
         cgl=>cgl%nxt
      enddo

   end subroutine cresp_update_grid

!----------------------------------------------------------------------------------------------------

   subroutine cresp_clean_grid

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use cresp_crspectrum, only: detect_clean_spectrum
      use grid_cont,        only: grid_container
      use initcrspectrum,   only: cresp, nullify_empty_bins
      use named_array,      only: p4
      use named_array_list, only: wna

      implicit none

      integer                         :: i, j, k
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      logical                         :: empty_cell

      if (nullify_empty_bins) then ! else nothing is done here
         i = 0; j = 0;  k = 0
         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg
            p4 => cg%w(wna%fi)%arr
            do k = cg%ks, cg%ke
               do j = cg%js, cg%je
                  do i = cg%is, cg%ie
                     cresp%n    = p4(iarr_cre_n, i, j, k)
                     cresp%e    = p4(iarr_cre_e, i, j, k)
                     empty_cell = .true.

                     call detect_clean_spectrum(cresp%n, cresp%e, empty_cell)

                     p4(iarr_cre_n, i, j, k) = cresp%n
                     p4(iarr_cre_e, i, j, k) = cresp%e
                  enddo
               enddo
            enddo
            cg%u(iarr_cre_n, :,:,:) = p4(iarr_cre_n, :,:,:)
            cg%u(iarr_cre_e, :,:,:) = p4(iarr_cre_e, :,:,:)
            cgl=>cgl%nxt
         enddo
      endif

   end subroutine cresp_clean_grid

end module cresp_grid
