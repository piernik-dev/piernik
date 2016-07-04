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

  subroutine cresp_p_cut_diff(credim)
  use initcosmicrays, only: iarr_cre, iarr_cre_pl, iarr_cre_pu
  use initcrspectrum, only: ncre
  use cg_leaves,        only: leaves
  use cg_list,          only: cg_list_element
  use fluidtypes,       only: var_numbers
  use grid_cont,        only: grid_container
  use domain,  only: dom
  use constants, only: I_ZERO, I_ONE, xdim, ydim, zdim, ndims
  implicit none
  integer       :: i,j,k, il, jl, kl, ih, jh, kh, ild, jld, kld
   type(grid_container), pointer   :: cg
   type(cg_list_element), pointer  :: cgl
   integer(kind=4),allocatable, dimension(:,:,:)  :: i_lo_dom !dimension(dom%n_d(1),dom%n_d(2),dom%n_d(3))
   integer(kind=4),allocatable, dimension(:,:,:)  :: i_up_dom !dimension(dom%n_d(1),dom%n_d(2),dom%n_d(3))
   integer(kind=4),intent(in) :: credim
   integer, dimension(ndims)            :: idm, ndm, hdm, ldm
   
   print *, 'dim@cresp_p_cut_diff = ', credim
   
   allocate(i_lo_dom(0:dom%n_d(1)-I_ONE,0:dom%n_d(2)-I_ONE,0:dom%n_d(3)-I_ONE))
   allocate(i_up_dom(0:dom%n_d(1)-I_ONE,0:dom%n_d(2)-I_ONE,0:dom%n_d(3)-I_ONE))

   i_lo_dom(:,:,:) = I_ZERO
   i_up_dom(:,:,:) = I_ZERO

      cgl => leaves%first
! Determinig cutoff momenta to show which cree and cren are going to be used
      do while (associated(cgl))
         cg => cgl%cg
         do k = cg%ks, cg%ke
           do j = cg%js, cg%je
              do i = cg%is, cg%ie
                i_lo_dom(i,j,k) = i_lo_detect(cg%u(iarr_cre_pl,i,j,k))
                i_up_dom(i,j,k) = i_up_detect(cg%u(iarr_cre_pu,i,j,k))
! #ifdef VERBOSE                
!                 print *,'@p_diffusion(', i,j,k, '): i_lo = ', cg%u(iarr_cre_pl,i,j,k), i_lo_dom(i,j,k) , ', i_up = ,', cg%u(iarr_cre_pu,i,j,k), i_up_dom(i,j,k)
! #endif /* VERBOSE */
              enddo
           enddo
         enddo
      cgl => cgl%nxt  
      enddo
      
      
!       ldm        = cg%ijkse(:,LO) ;      ldm(credim) = cg%lhn(credim,LO) + dom%D_(credim)      ! ldm =           1 + D_
!       hdm        = cg%ijkse(:,HI) ;      hdm(credim) = cg%lhn(credim,HI)                      ! hdm = cg%n_ + idm - D_
      do k = ldm(zdim), hdm(zdim)       ; kl = k-1 ; kh = k+1 ; kld = k-idm(zdim)
         do j = ldm(ydim), hdm(ydim)    ; jl = j-1 ; jh = j+1 ; jld = j-idm(ydim)
            do i = ldm(xdim), hdm(xdim) ; il = i-1 ; ih = i+1 ; ild = i-idm(xdim)
!            
            enddo
         enddo
      enddo
!       
  if(allocated(i_lo_dom)) deallocate(i_lo_dom)
  if(allocated(i_up_dom)) deallocate(i_up_dom)
  
  end subroutine cresp_p_cut_diff
  
!------------------------------------------  
  function i_lo_detect(p_lo)
  use initcrspectrum, only: ncre, p_fix, w
  implicit none
   real(kind=8), intent(in) :: p_lo
   integer(kind=4)   :: i_lo_detect
! We are intrested in bin number in this case, so (*) is increased by one.
   
    i_lo_detect = int(floor(log10(p_lo/p_fix(1))/w)) + 1
    i_lo_detect = max(0, i_lo_detect)
    i_lo_detect = min(i_lo_detect, ncre - 1) + 1  ! (*)
!     print *,'iul = ', i_lo_detect, int(floor(log10(p_lo/p_fix(1))/w))

  end function i_lo_detect
 
!------------------------------------------
  function i_up_detect(p_up)
  use initcrspectrum, only: ncre, p_fix, w
!   use constants, only: I_ONE
  implicit none 
   real(kind=8), intent(in) :: p_up
   integer(kind=4)   :: i_up_detect

    i_up_detect = int(floor(log10(p_up/p_fix(1))/w)) + 2
    i_up_detect = max(1,i_up_detect)
    i_up_detect = min(i_up_detect,ncre)
!       print *, 'iud = ', i_up_detect
  
  end function i_up_detect

end module cresp_p_diffusion
