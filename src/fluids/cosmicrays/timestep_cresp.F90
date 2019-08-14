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
!! \brief Computation of %timestep for energy & number density spectrum evolution in momentum space via CRESP algorithm
!<

module timestep_cresp
! pulled by COSM_RAY_ELECTRONS

   implicit none

   private
   public :: dt_cre, cresp_timestep, dt_cre_min_ub, dt_cre_min_ud, dt_cre_K

   real(kind=8)    :: dt_cre, dt_cre_min_ub, dt_cre_min_ud, dt_cre_K
   integer(kind=4) :: i_up_max_prev, i_up_max

contains

   function assume_p_up(cell_i_up)

      use initcosmicrays, only: ncre
      use initcrspectrum, only: p_fix, p_mid_fix

      implicit none

      integer(kind=4), intent(in) :: cell_i_up
      real(kind=8)                :: assume_p_up

      assume_p_up = p_fix(ncre-1)
      if (cell_i_up .eq. ncre) then
         assume_p_up = p_mid_fix(ncre) ! for i = 0 & ncre p_fix(i) = 0.0
      else
         assume_p_up = p_fix(cell_i_up)
      endif

   end function assume_p_up
!----------------------------------------------------------------------------------------------------

   function evaluate_i_up(e_cell, n_cell) ! obtain i_up index from energy densities in cell

      use constants,      only: zero
      use initcosmicrays, only: ncre
      use initcrspectrum, only: e_small

      implicit none

      real(kind=8), dimension(:), intent(in) :: e_cell, n_cell
      integer                                :: evaluate_i_up, i

      do i = ncre, 1, -1  ! we start counting from ncre since upper cutoff is rather expected at higher index numbers. Might change it though.
         if (e_cell(i) .gt. e_small .and. n_cell(i) .gt. zero) then ! better compare to zero or to eps?
            evaluate_i_up = i
            return ! if cell is empty, evaluate_i_up returns 0, which is handled by cresp_timestep
         endif
      enddo
      evaluate_i_up = 0

   end function evaluate_i_up

!----------------------------------------------------------------------------------------------------

   subroutine cresp_timestep

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, half, zero, big
      use cresp_grid,       only: fsynchr
      use cresp_crspectrum, only: cresp_find_prepare_spectrum
      use crhelpers,        only: div_v, divv_n
      use fluidindex,       only: flind
      use func,             only: emag
      use grid_cont,        only: grid_container
      use initcosmicrays,   only: K_cre_paral, K_cre_perp, cfl_cr, iarr_cre_e, iarr_cre_n
      use initcrspectrum,   only: spec_mod_trms, synch_active, adiab_active, use_cresp, cresp
      use named_array_list, only: qna

      implicit none

      integer(kind=4)                :: i, j, k, i_up_max_tmp
      type(grid_container),  pointer :: cg
      type(cg_list_element), pointer :: cgl
      real(kind=8)                   :: dt_cre_tmp, K_cre_max_sum, abs_max_ud
      type(spec_mod_trms)            :: sptab
      logical                        :: empty_cell

      dt_cre        = big
      dt_cre_tmp    = big
      dt_cre_min_ub = big
      dt_cre_min_ud = big

      if (.not. use_cresp) return

      abs_max_ud   = zero
      i_up_max     = 1
      i_up_max_tmp = 1

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         if (adiab_active) call div_v(flind%ion%pos, cg)
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  sptab%ud = 0.0 ; sptab%ub = 0.0 ; sptab%ucmb = 0.0 ; empty_cell = .false.
                  if (synch_active) then
                     sptab%ub = emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k)) * fsynchr
                     cresp%n = cg%u(iarr_cre_n, i, j, k)
                     cresp%e = cg%u(iarr_cre_e, i, j, k)
                     call cresp_find_prepare_spectrum(cresp%n, cresp%e, empty_cell, i_up_max_tmp) ! needed for synchrotron timestep

                     if (empty_cell) cycle         ! then nothing to do in this iteration

                     call cresp_timestep_synchrotron(dt_cre_tmp, sptab%ub, i_up_max_tmp)
                     dt_cre = min(dt_cre, dt_cre_tmp)
                     i_up_max = max(i_up_max, i_up_max_tmp)
                  endif
                  if (adiab_active) abs_max_ud = max(abs_max_ud, abs(cg%q(qna%ind(divv_n))%point([i,j,k])))
               enddo
            enddo
         enddo
         cgl=>cgl%nxt
      enddo

      if (adiab_active) call cresp_timestep_adiabatic(dt_cre_tmp, abs_max_ud)

      dt_cre = min(dt_cre, dt_cre_tmp)
      if (i_up_max_prev .ne. i_up_max) then ! dt_cre_K saved, computed again only if in the whole domain highest i_up changes.
         i_up_max_prev = i_up_max
         K_cre_max_sum = K_cre_paral(i_up_max) + K_cre_perp(i_up_max) ! assumes the same K for energy and number density
         if ( K_cre_max_sum <= 0) then                                ! K_cre dependent on momentum - maximal for highest bin number
            dt_cre_K = big
         else                                                         ! We use cfl_cr here (CFL number for diffusive CR transport)
            dt_cre_K = cfl_cr * half / K_cre_max_sum                  ! cfl_cre used only for spectrum evolution
            if (cg%dxmn < sqrt(big)/dt_cre_K) dt_cre_K = dt_cre_K * cg%dxmn**2
         endif
      endif

      dt_cre = min(dt_cre, dt_cre_K)
      dt_cre = half * dt_cre                  ! dt comes in to cresp_crspectrum with factor * 2

   end subroutine cresp_timestep

!----------------------------------------------------------------------------------------------------

   subroutine cresp_timestep_adiabatic(dt_cre_ud, u_d_abs)

      use constants,       only: logten, three
      use initcrspectrum,  only: w, eps, cfl_cre

      implicit none

      real(kind=8), intent(out) :: dt_cre_ud
      real(kind=8), intent(in)  :: u_d_abs    ! assumes that u_d > 0 always

      if (u_d_abs .gt. eps) then
         dt_cre_ud = cfl_cre * three * logten * w / u_d_abs
         dt_cre_min_ud = min(dt_cre_ud, dt_cre_min_ud)      ! remember to max dt_cre_min_ud at beginning of the search!
      endif

   end subroutine cresp_timestep_adiabatic

!----------------------------------------------------------------------------------------------------

   subroutine cresp_timestep_synchrotron(dt_cre_ub, u_b, i_up_cell)

      use constants,      only: zero
      use initcrspectrum, only: w, cfl_cre

      implicit none

      real(kind=8),    intent(out) :: dt_cre_ub
      real(kind=8),    intent(in)  :: u_b
      integer(kind=4), intent(in)  :: i_up_cell

 ! Synchrotron cooling timestep (is dependant only on p_up, highest value of p):
      if (u_b .gt. zero) then
         dt_cre_ub = cfl_cre * w / (assume_p_up(i_up_cell) * u_b)
         dt_cre_min_ub = min(dt_cre_ub, dt_cre_min_ub)    ! remember to max dt_cre_min_ub at the beginning of the search
      endif

   end subroutine cresp_timestep_synchrotron

end module timestep_cresp
