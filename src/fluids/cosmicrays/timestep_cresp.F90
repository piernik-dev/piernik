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

   use constants,      only: one, zero
   use initcrspectrum, only: ncre, cfl_cre

   implicit none

   private
   public :: dt_cre, cresp_timestep, dt_cre_min_ub, dt_cre_min_ud

   real(kind=8) :: dt_cre, dt_cre_min_ub, dt_cre_min_ud

 contains

   function assume_p_up(cell_i_up)

      use initcrspectrum, only: p_fix, p_mid_fix, ncre

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
      use initcrspectrum, only: ncre, e_small

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
! Subroutine consistent with rules depicted in crspectrum.pdf
!----------------------------------------------------------------------------------------------------

   subroutine cresp_timestep(dt_comp, sptab, n_cell, e_cell, i_up_cell)

      use constants,        only: zero, I_ZERO
      use cresp_crspectrum, only: cresp_find_prepare_spectrum
      use initcrspectrum,   only: spec_mod_trms, ncre, w, eps

      implicit none

      real(kind=8),               intent(out) :: dt_comp
      type(spec_mod_trms)                     :: sptab
      real(kind=8), dimension(:), intent(in)  :: n_cell, e_cell
      integer(kind=4),            intent(out) :: i_up_cell
      real(kind=8)                            :: dt_cre_ud, dt_cre_ub, p_u
      real(kind=8), dimension(ncre)           :: n_inout, e_inout
      logical                                 :: empty_cell

      empty_cell = .true.
      n_inout = n_cell
      e_inout = e_cell
      i_up_cell = I_ZERO
      dt_cre_ud = huge(one)
      dt_cre_ub = huge(one)
      call cresp_find_prepare_spectrum(n_inout, e_inout, empty_cell, i_up_cell)

! cell is assumed empty if evaluate_i_up over whole ncre range returns 0 -> nothing to do here
      if (i_up_cell .gt. 0) then
! Adiabatic cooling timestep:

         if ( abs(sptab%ud) .gt. eps) then
            dt_cre_ud = cfl_cre * w / sptab%ud
            dt_cre_ud = abs(dt_cre_ud)
            dt_cre_min_ud = min(dt_cre_ud, dt_cre_min_ud)
         endif

! Synchrotron cooling timestep (is dependant only on p_up, highest value of p):
         if (sptab%ub .gt. zero) then
!           i_up_cell = evaluate_i_up(e_cell, n_cell)
            p_u = assume_p_up(i_up_cell) ! TODO: fix problems with negative p_u
            dt_cre_ub = cfl_cre * w / (p_u * sptab%ub)
            dt_cre_min_ub = min(dt_cre_ub, dt_cre_min_ub)
         endif
      endif

! Here shortest among calculated timesteps is chosen.
      dt_comp = min(dt_cre_ud, dt_cre_ub)

! Should dt_cre_ud or dt_cre_ub be greater than current one, next timestep shall be independently increased by piernik in timestep
      dt_cre = min(dt_cre, dt_comp) ! Assures that minimal timestep among computed for current cell and previous ones us chosen

   end subroutine cresp_timestep

end module timestep_cresp
