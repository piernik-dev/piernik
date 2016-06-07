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
  use initcosmicrays, only: iarr_cre
  use cg_leaves,        only: leaves
  use cg_list,          only: cg_list_element
  use fluidtypes,       only: var_numbers
  use grid_cont,        only: grid_container
  implicit none
  integer       :: i,j,k!, id, jd, kd
  integer       :: i_lo, i_up
   type(grid_container), pointer   :: cg
   type(cg_list_element), pointer  :: cgl
   type(var_numbers)    :: flind
  
!       wcri = wna%ind(wcr_n)
      
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
!          wcr => cg%w(wcri)%arr
           
         do k = cg%ks, cg%ke
           do j = cg%js, cg%je
              do i = cg%is, cg%ie
              
                i_lo = i_lo_detect(cg%u(flind%cre%plo,i,j,k))
                i_up = i_up_detect(cg%u(flind%cre%plo,i,j,k))
! #ifdef VERBOSE                
                print *,'@p_diffusion(', i,j,k, '): i_lo = ', i_lo , ', i_up = ,', i_up
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
   integer   :: m
   integer(kind=4)   :: i_lo_detect
  
  m = 0
  do while (m .le. ncre)
      if ( p_lo .lt. p_fix(m)) exit
      i_lo_detect = m
      m = m + I_ONE
      print *,m
  enddo
  
  
!   return p_up_ind
  end function i_lo_detect
 
!------------------------------------------
  function i_up_detect(p_up)
  use initcrspectrum, only: ncre, p_fix
  use constants, only: I_ONE
  implicit none 
   real(kind=8), intent(in) :: p_up
   integer   :: n
   integer(kind=4)   :: i_up_detect
  
  n = ncre
  do while (n .ge. 0)
    i_up_detect = n
    n = n - I_ONE
    if (p_up .lt.  p_fix(n)) exit
  enddo
  
!   return p_up_ind
  end function i_up_detect

end module cresp_p_diffusion
