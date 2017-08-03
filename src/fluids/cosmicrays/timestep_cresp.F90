module timestep_cresp
! pulled by COSM_RAY_ELECTRONS

 use initcrspectrum, only: ncre, cfl_cre
 use constants,      only: one, zero

      implicit none
!       private         p_tmp
 
      real(kind=8) :: dt_cre, dt_tmp
!       real(kind=8), allocatable, dimension(:) :: p_tmp   ! do not use p variable in this module anywhere outside
!       real(kind=8), allocatable, dimension(:) :: dts_new

!-------------------      
!      
 contains
!
!-------------------
! if we only saved p_up value, this might save a LOT of time.

  function approximate_p_up(n_cell, e_cell)
   use initcrspectrum, only: p_fix
   use cresp_NR_method, only: intpol_pf_from_NR_grids, alpha_tab_up, n_tab_up
   use units, only: clight
   implicit none
    integer :: cell_i_up
    real(kind=8), dimension(:) :: n_cell, e_cell
    real(kind=8) :: approximate_p_up
    real(kind=8), dimension(1:2) :: pf_ratio
        cell_i_up = evaluate_i_up(e_cell)
!         alpha_inp = (e_cell(cell_i_up)/(n_inp(cell_i_up)*clight*p_fix(cell_i_up-1)))
!         n_inp     = n_cell(cell_i_up)
        pf_ratio  = intpol_pf_from_NR_grids("up",(e_cell(cell_i_up)/(n_cell(cell_i_up)*clight*p_fix(cell_i_up-1))), &
                                                        n_cell(cell_i_up), alpha_tab_up, n_tab_up) ! we use just an interpolated ratio, who knows if it'll work
        approximate_p_up = pf_ratio(1) * p_fix(cell_i_up-1)
    end function approximate_p_up
!----------------------------------------------------------------------------------------------------
  function evaluate_i_up(e_cell) ! obain i_up index from energy densities in cell
  use initcrspectrum, only: ncre
  use constants, only: zero
  implicit none
    real(kind=8), dimension(:) :: e_cell
    integer :: evaluate_i_up, i
        do i=ncre, 1, -1  ! we start counting from ncre since upper cutoff is rather expected at higher index numbers. Might change it though.
            if (e_cell(i) .gt. zero) then ! better compare to zero or to eps?
                evaluate_i_up = i
                return
            endif  ! no need for other conditions - if there IS a bin that has literally no energy, the algorithm will most likely crash.
        enddo
  end function evaluate_i_up
  
!----------------------------------------------------------------------------------------------------
! Subroutine consistent with rules depicted in crspectrum.pdf
!----------------------------------------------------------------------------------------------------

  subroutine cresp_timestep(dt_comp, sptab, n_cell, e_cell)
   use initcrspectrum,   only: spec_mod_trms, ncre, w
   use cresp_arrays_handling, only: allocate_with_index
   use constants, only: zero
   implicit none
    real(kind=8)               :: dt_comp
    real(kind=8)               :: dt_cre_ud, dt_cre_ub
    real(kind=8) :: p_u
    real(kind=8), dimension(:), intent(in) :: n_cell, e_cell
    type(spec_mod_trms) sptab
        dt_cre_ud = huge(one)
        dt_cre_ub = huge(one)

        if (maxval(e_cell) .gt. zero) then ! any timestep evaluation makes sense only if there's any information to be migrated between bins
! Adiabatic cooling timestep:
            if (sptab%ud .ne. zero) then
                dt_cre_ud = cfl_cre * w / sptab%ud
                dt_cre_ud = abs(dt_cre_ud)
            endif
!    dt_cre_ud = minval(dts_new)   ! it was already multiplied by cfl_cre
!    dts_new = huge(one)
   
! Synchrotron cooling timestep (is dependant only on p_up, highest value of p):
            p_u = approximate_p_up(n_cell, e_cell)
            if (sptab%ub .gt. zero) then
                dt_cre_ub = cfl_cre * w / (p_u * sptab%ub)
            endif
!    dt_cre_ub = minval(dts_new)
#ifdef VERBOSE
            print *, '[@timestep_cresp:] Computed timesteps:'
            print *, 'dt_cre_ud = ', dt_cre_ud
            print *, 'dt_cre_ub = ', dt_cre_ub
#endif /* VERBOSE */
        endif
! Here shortest among calculated timesteps is chosen.
        dt_comp = min(dt_cre_ud, dt_cre_ub)

! Should dt_cre_ud or dt_cre_ub be greater than current one, next timestep shall be independently increased by piernik in timestep
        dt_cre = min(dt_cre, dt_comp) ! Gives minimal timestep among computed for current cell and previous ones
   
!     if (allocated(p_tmp))   deallocate(p_tmp)
!     if (allocated(dts_new)) deallocate(dts_new)
    
 end subroutine cresp_timestep
  
end module timestep_cresp
