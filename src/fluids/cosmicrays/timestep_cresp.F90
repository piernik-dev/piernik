module timestep_cresp
! pulled by COSM_RAY_ELECTRONS

   use initcrspectrum, only: ncre, cfl_cre
   use constants,      only: one, zero

   implicit none

   private
   public :: dt_cre, cresp_timestep, dt_cre_min_ub, dt_cre_min_ud

   real(kind=8) :: dt_cre, dt_cre_min_ub, dt_cre_min_ud

!----------------------------------------------------------------------------------------------------
!
 contains
!
!----------------------------------------------------------------------------------------------------

   function assume_p_up(cell_i_up)
      use initcrspectrum, only: p_fix, p_mid_fix, ncre
      implicit none

      integer(kind=4), intent(in)            :: cell_i_up
      real(kind=8)                           :: assume_p_up

      assume_p_up = p_fix(ncre-1)
      if (cell_i_up .eq. ncre) then
         assume_p_up = p_mid_fix(ncre) ! for i = 0 & ncre p_fix(i) = 0.0
         return
      else
         assume_p_up = p_fix(cell_i_up)
         return
      endif

   end function assume_p_up
!----------------------------------------------------------------------------------------------------

   function evaluate_i_up(e_cell, n_cell) ! obtain i_up index from energy densities in cell

      use initcrspectrum, only: ncre, e_small
      use constants, only: zero

      implicit none

      real(kind=8), dimension(:), intent(in) :: e_cell, n_cell
      integer :: evaluate_i_up, i

      do i=ncre, 1, -1  ! we start counting from ncre since upper cutoff is rather expected at higher index numbers. Might change it though.
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

      use initcrspectrum,   only: spec_mod_trms, ncre, w, eps
      use cresp_crspectrum, only: cresp_find_prepare_spectrum
      use constants, only: zero, I_ZERO

      implicit none

      real(kind=8), intent(out)  :: dt_comp
      real(kind=8)               :: dt_cre_ud, dt_cre_ub
      real(kind=8) :: p_u
      real(kind=8), dimension(:), intent(in) :: n_cell, e_cell
      real(kind=8), dimension(ncre)          :: n_inout, e_inout
      integer(kind=4), intent(inout) :: i_up_cell
      type(spec_mod_trms) sptab
      logical      :: empty_cell

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
