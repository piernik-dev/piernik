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
   public :: dt_cre, cresp_timestep, dt_cre_synch, dt_cre_adiab, dt_cre_K

   real :: dt_cre, dt_cre_synch, dt_cre_adiab, dt_cre_K

contains

   subroutine cresp_timestep

      use all_boundaries,   only: all_fluid_boundaries
      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim, half, zero, big, pMIN
      use cresp_crspectrum, only: cresp_find_prepare_spectrum
      use crhelpers,        only: div_v, divv_n
      use fluidindex,       only: flind
      use func,             only: emag
      use grid_cont,        only: grid_container
      use initcosmicrays,   only: K_cre_paral, K_cre_perp, cfl_cr, iarr_cre_e, iarr_cre_n
      use initcrspectrum,   only: spec_mod_trms, synch_active, adiab_active, use_cresp, cresp, fsynchr
      use mpisetup,         only: piernik_MPI_Allreduce
      use named_array_list, only: qna

      implicit none

      integer(kind=4)                :: i, j, k, i_up_max_tmp, i_up_max
      type(grid_container),  pointer :: cg
      type(cg_list_element), pointer :: cgl
      type(spec_mod_trms)            :: sptab
      real                           :: K_cre_max_sum, abs_max_ud
      logical                        :: empty_cell

      dt_cre       = big
      dt_cre_K     = big
      dt_cre_synch = big
      dt_cre_adiab = big

      if (.not. use_cresp) return

      abs_max_ud   = zero
      i_up_max     = 1
      i_up_max_tmp = 1
      if (adiab_active) call all_fluid_boundaries()

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         if (adiab_active) then
            call div_v(flind%ion%pos, cg)
            abs_max_ud = max(abs_max_ud, maxval(abs(cg%q(qna%ind(divv_n))%span(cg%ijkse))))
         endif

         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  sptab%ud = zero ; sptab%ub = zero ; sptab%ucmb = zero ; empty_cell = .false.
                  sptab%ub = emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k)) * fsynchr
                  cresp%n = cg%u(iarr_cre_n, i, j, k)
                  cresp%e = cg%u(iarr_cre_e, i, j, k)
                  call cresp_find_prepare_spectrum(cresp%n, cresp%e, empty_cell, i_up_max_tmp) ! needed for synchrotron timestep
                  i_up_max = max(i_up_max, i_up_max_tmp)

                  if (.not. empty_cell .and. synch_active) call cresp_timestep_synchrotron(sptab%ub, i_up_max_tmp)
               enddo
            enddo
         enddo

         cgl=>cgl%nxt
      enddo

      if (adiab_active) call cresp_timestep_adiabatic(abs_max_ud)

      K_cre_max_sum = K_cre_paral(i_up_max) + K_cre_perp(i_up_max) ! assumes the same K for energy and number density
      if (K_cre_max_sum > zero) then                               ! K_cre dependent on momentum - maximal for highest bin number
         dt_cre_K = cfl_cr * half / K_cre_max_sum                  ! We use cfl_cr here (CFL number for diffusive CR transport), cfl_cre used only for spectrum evolution
         if (cg%dxmn < sqrt(big)/dt_cre_K) dt_cre_K = dt_cre_K * cg%dxmn**2
      endif

      call piernik_MPI_Allreduce(dt_cre_adiab, pMIN)
      call piernik_MPI_Allreduce(dt_cre_synch, pMIN)
      call piernik_MPI_Allreduce(dt_cre_K,     pMIN)
      dt_cre = half * min(dt_cre_adiab, dt_cre_synch, dt_cre_K)       ! dt comes in to cresp_crspectrum with factor * 2

   end subroutine cresp_timestep

!----------------------------------------------------------------------------------------------------

   subroutine cresp_timestep_adiabatic(u_d_abs)

      use constants,      only: logten, three
      use initcrspectrum, only: w, eps, cfl_cre

      implicit none

      real, intent(in) :: u_d_abs    ! assumes that u_d > 0 always

      if (u_d_abs > eps) dt_cre_adiab = cfl_cre * three * logten * w / u_d_abs

   end subroutine cresp_timestep_adiabatic

!----------------------------------------------------------------------------------------------------

   subroutine cresp_timestep_synchrotron(u_b, i_up_cell)

      use constants,      only: zero
      use initcrspectrum, only: w, cfl_cre

      implicit none

      real,            intent(in) :: u_b
      integer(kind=4), intent(in) :: i_up_cell
      real                        :: dt_cre_ub

 ! Synchrotron cooling timestep (is dependant only on p_up, highest value of p):
      if (u_b > zero) then
         dt_cre_ub = cfl_cre * w / (assume_p_up(i_up_cell) * u_b)
         dt_cre_synch = min(dt_cre_ub, dt_cre_synch)    ! remember to max dt_cre_synch at the beginning of the search
      endif

   end subroutine cresp_timestep_synchrotron

!----------------------------------------------------------------------------------------------------

   real function assume_p_up(cell_i_up)

      use initcosmicrays, only: ncre
      use initcrspectrum, only: p_fix, p_mid_fix

      implicit none

      integer(kind=4), intent(in) :: cell_i_up

      if (cell_i_up == ncre) then
         assume_p_up = p_mid_fix(ncre) ! for i = 0 & ncre p_fix(i) = 0.0
      else
         assume_p_up = p_fix(cell_i_up)
      endif

   end function assume_p_up

end module timestep_cresp
