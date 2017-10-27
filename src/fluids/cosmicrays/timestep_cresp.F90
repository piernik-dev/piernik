module timestep_cresp
! pulled by COSM_RAY_ELECTRONS

 use initcrspectrum, only: ncre, cfl_cre
 use constants,      only: one, zero

      implicit none
    private
    public :: dt_cre, cresp_timestep

    real(kind=8) :: dt_cre

!----------------------------------------------------------------------------------------------------
!
 contains
!
!----------------------------------------------------------------------------------------------------

  function approximate_p_up(n_cell, e_cell, cell_i_up)
   use initcrspectrum, only: p_fix
   use cresp_NR_method, only: intpol_pf_from_NR_grids
   use units, only: clight
   implicit none
    integer(kind=4), intent(in)            :: cell_i_up
    real(kind=8), dimension(:), intent(in) :: n_cell, e_cell
    real(kind=8) :: approximate_p_up, alpha_bnd, n_bnd
    real(kind=8), dimension(1:2) :: pf_ratio
    character(len=2) :: which_bound="up"
    logical :: interpolation_successful, intpol_fail
        intpol_fail = .true.
        if (cell_i_up .ne. 1) then
            alpha_bnd = (e_cell(cell_i_up)/(n_cell(cell_i_up)*clight*p_fix(cell_i_up-1)))
            n_bnd     = n_cell(cell_i_up)
            pf_ratio  = intpol_pf_from_NR_grids(which_bound, alpha_bnd, n_bnd, interpolation_successful, intpol_fail) ! we use just an interpolated ratio
            if ( intpol_fail ) then
                approximate_p_up = p_fix(cell_i_up) ! if interpolation fails, upper p_fix boundary is provided - might cause problems with p_up moving beyond p_fix
                return
            else
                approximate_p_up = pf_ratio(1) * p_fix(cell_i_up-1)
                return
            endif
        else
            approximate_p_up = p_fix(cell_i_up) ! approximate_p_up = pf_ratio(1) * p_fix(cell_i_up-1)
        endif
  end function approximate_p_up
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
                return
            endif  ! no need for other conditions - if there IS a bin that has literally no energy, the algorithm will most likely crash.
        enddo
  end function evaluate_i_up

!----------------------------------------------------------------------------------------------------
! Subroutine consistent with rules depicted in crspectrum.pdf
!----------------------------------------------------------------------------------------------------

  subroutine cresp_timestep(dt_comp, sptab, n_cell, e_cell, i_up_cell)
   use initcrspectrum,   only: spec_mod_trms, ncre, w, e_small
   use constants, only: zero, I_ZERO
   implicit none
    real(kind=8), intent(out)  :: dt_comp
    real(kind=8)               :: dt_cre_ud, dt_cre_ub
    real(kind=8) :: p_u
    real(kind=8), dimension(:), intent(in) :: n_cell, e_cell
    integer(kind=4), intent(inout) :: i_up_cell
    type(spec_mod_trms) sptab
        i_up_cell = I_ZERO
        dt_cre_ud = huge(one)
        dt_cre_ub = huge(one)
        if (maxval(e_cell) .gt. e_small) then ! any timestep evaluation makes sense only if there's any information to be migrated between bins
! Adiabatic cooling timestep:
            if (abs(sptab%ud) .gt. zero) then
                dt_cre_ud = cfl_cre * w / sptab%ud
                dt_cre_ud = abs(dt_cre_ud)
            endif
! Synchrotron cooling timestep (is dependant only on p_up, highest value of p):
            if (sptab%ub .gt. zero) then
                i_up_cell = evaluate_i_up(e_cell, n_cell)
                p_u = approximate_p_up(n_cell, e_cell, i_up_cell)
                dt_cre_ub = cfl_cre * w / (p_u * sptab%ub)
            endif
        endif
! Here shortest among calculated timesteps is chosen.
        dt_comp = min(dt_cre_ud, dt_cre_ub)

! Should dt_cre_ud or dt_cre_ub be greater than current one, next timestep shall be independently increased by piernik in timestep
        dt_cre = min(dt_cre, dt_comp) ! Assures that minimal timestep among computed for current cell and previous ones us chosen
 end subroutine cresp_timestep

end module timestep_cresp
