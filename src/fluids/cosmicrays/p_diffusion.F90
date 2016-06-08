! In this module resides algorithm used to update cutoff momenta of cosmic ray electron spectrum.
 module cresp_p_diffusion
! pulled by COSM_RAY_ELECTRONS
!    logical, dimension(:)             :: mask_flux_pos, mask_flux_neg   ! arrays of masks determining direction of momenta diffusion
   use constants, only: I_ONE
   implicit none
! --------------
contains
!   
!   
!   
!   

  subroutine cresp_p_cut_diff
  use initcosmicrays, only: iarr_cre, iarr_cre_pl, iarr_cre_pu
  use cg_leaves,        only: leaves
  use cg_list,          only: cg_list_element
  use fluidtypes,       only: var_numbers
  use grid_cont,        only: grid_container
  implicit none
  integer       :: i,j,k!, id, jd, kd
  integer       :: i_lo, i_up
   type(grid_container), pointer   :: cg
   type(cg_list_element), pointer  :: cgl
 
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         do k = cg%ks, cg%ke
           do j = cg%js, cg%je
              do i = cg%is, cg%ie
                i_lo = i_lo_detect(cg%u(iarr_cre_pl,i,j,k))
                i_up = i_up_detect(cg%u(iarr_cre_pu,i,j,k))
! #ifdef VERBOSE                
!                 print *,'@p_diffusion(', i,j,k, '): i_lo = ', i_lo , ', i_up = ,', i_up 
! #endif /* VERBOSE */
              enddo
           enddo
         enddo
         
      cgl => cgl%nxt  
      enddo
  
  end subroutine cresp_p_cut_diff
  
!------------------------------------------  
  function i_lo_detect(p_lo)
  use initcrspectrum, only: ncre, p_fix
  implicit none
   real(kind=8), intent(in) :: p_lo
   integer   :: i
   integer(kind=4)   :: i_lo_detect
  
  i = 1
  do while (i .le. ncre)
      i_lo_detect = i
!       print *, 'p_lo', i, p_lo, p_fix(i)
      if ( p_lo .lt. p_fix(i)) exit ! checked and works properly
      i = i + I_ONE

  enddo
  
  
!   return p_up_ind
  end function i_lo_detect
 
!------------------------------------------
  function i_up_detect(p_up)
  use initcrspectrum, only: ncre, p_fix
  use constants, only: I_ONE
  implicit none 
   real(kind=8), intent(in) :: p_up
   integer   :: i
   integer(kind=4)   :: i_up_detect
  
  i = ncre - I_ONE
  do while (i .ge. 0)
    i_up_detect = i+I_ONE
!     print *, 'p_up ', i, p_up, p_fix(i)
    if (p_up .gt.  p_fix(i)) exit   ! checked and works properly - if p_fix(i) is lt p_up, then e(max) must be in i-th + 1 bin.
    i = i - I_ONE
  enddo
  
!   return p_up_ind
  end function i_up_detect

end module cresp_p_diffusion
