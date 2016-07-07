! In this module resides algorithm used to update cutoff momenta of cosmic ray electron spectrum.
 module cresp_p_diffusion
! pulled by COSM_RAY_ELECTRONS
!    logical, dimension(:)             :: mask_flux_pos, mask_flux_neg   ! arrays of masks determining direction of momenta diffusion
   use constants, only: I_ONE, ndims
   use fluidindex,       only: flind
   use initcrspectrum,   only: flux_cre  ! WARNING - flux_cre does not 'remember' diffusive fluxes in all directions, only the current one
   implicit none
!    logical         :: p_cut_computed
! --------------
contains
!   
!   
!   
!   

  subroutine cresp_p_cut_diff(credim)
  use initcosmicrays, only: iarr_cre_pl, iarr_cre_pu
  use initcrspectrum, only: ncre
  use cg_leaves,        only: leaves
  use cg_list,          only: cg_list_element
  use fluidtypes,       only: var_numbers
  use grid_cont,        only: grid_container
  use domain,  only: dom
  use constants, only: I_ZERO, I_ONE, xdim, ydim, zdim
  implicit none
  integer       :: i,j,k, il, jl, kl, ih, jh, kh!, ild, jld, kld
   type(grid_container), pointer   :: cg
   type(cg_list_element), pointer  :: cgl
   integer(kind=4),allocatable, dimension(:,:,:)  :: i_lo_dom !dimension(dom%n_d(1),dom%n_d(2),dom%n_d(3))
   integer(kind=4),allocatable, dimension(:,:,:)  :: i_up_dom !dimension(dom%n_d(1),dom%n_d(2),dom%n_d(3))
   integer(kind=4),intent(in) :: credim
   integer, dimension(ndims)            :: hdm, ldm!, ndm, idm 
   logical, dimension(ndims)            :: current_direction
   real, dimension(0:dom%n_d(1)-I_ONE,0:dom%n_d(2)-I_ONE,0:dom%n_d(3)-I_ONE) :: p_lo_new, p_up_new
   
!    print *, 'dim@cresp_p_cut_diff = ', credim
    
      if (.not.dom%has_dir(credim)) return
      
   current_direction = dom%has_dir .and. ( [ xdim,ydim,zdim ] /= credim )
   
   allocate(i_lo_dom(0:dom%n_d(1)-I_ONE,0:dom%n_d(2)-I_ONE,0:dom%n_d(3)-I_ONE))
   allocate(i_up_dom(0:dom%n_d(1)-I_ONE,0:dom%n_d(2)-I_ONE,0:dom%n_d(3)-I_ONE))

   i_lo_dom(:,:,:) = I_ZERO
   i_up_dom(:,:,:) = I_ZERO


      cgl => leaves%first
! Determinig cutoff momenta to show which cree and cren are going to be used
      do while (associated(cgl))
         cg => cgl%cg
!          if (p_cut_computed .eqv. .false.) then
         
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

!       print *, one_posnegzero(0.000250)
!       print *, one_posnegzero(-0.000005)
!       print *, one_posnegzero(0.0)

! WARNING - wrong indexing and whole array flux_cre is equal zero
      do k = ldm(zdim), hdm(zdim)       ; kl = k-1 ; kh = k+1
         do j = ldm(ydim), hdm(ydim)    ; jl = j-1 ; jh = j+1
            do i = ldm(xdim), hdm(xdim) ; il = i-1 ; ih = i+1
! p_lo update - value of p_lo will migrate accordingly to the sign of flux (1,-1 or 0), computed in one_posnegzero

              if (current_direction(xdim)) then
                p_lo_new(i,j,k) = min(cg%u(iarr_cre_pl,i + one_posnegzero(flux_cre(i_lo_dom(i,j,k),i,j,k)),j,k) &
                  , cg%u(iarr_cre_pl,i,j,k) )
              endif
              if (current_direction(ydim)) then
                p_lo_new(i,j,k) = min(cg%u(iarr_cre_pl,i ,j + one_posnegzero(flux_cre(i_lo_dom(i,j,k),i,j,k)), k) &
                  , cg%u(iarr_cre_pl,i,j,k) )
              endif
              if (current_direction(zdim)) then
                p_lo_new(i,j,k) = min(cg%u(iarr_cre_pl,i ,j, k + one_posnegzero(flux_cre(i_lo_dom(i,j,k),i,j,k))) &
                  , cg%u(iarr_cre_pl,i,j,k) )
              endif

! p_up update - value of p_up will migrate accordingly to the sign of flux (1,-1 or 0), computed in one_posnegzero           

              if (current_direction(xdim)) then
                 p_up_new(i,j,k) = min(cg%u(iarr_cre_pl,i + one_posnegzero(flux_cre(i_up_dom(i,j,k),i,j,k)),j,k) &
                  , cg%u(iarr_cre_pl,i,j,k) )
              endif
              if (current_direction(ydim)) then
                p_up_new(i,j,k) = min(cg%u(iarr_cre_pl,i ,j + one_posnegzero(flux_cre(i_up_dom(i,j,k),i,j,k)), k) &
                  , cg%u(iarr_cre_pl,i,j,k) )
              endif
              if (current_direction(zdim)) then
                p_up_new(i,j,k) = min(cg%u(iarr_cre_pl,i ,j, k + one_posnegzero(flux_cre(i_up_dom(i,j,k),i,j,k))) &
                  , cg%u(iarr_cre_pl,i,j,k) )
              endif
              
             enddo
         enddo
      enddo
      
      
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
  
!------------------------------------------  
  function one_posnegzero(value)
! Function computes and returns 0 / 1 / -1 depending on the sign of total cre energy, 
! minimalizing usage of if conditions.
   use constants,    only: one
   implicit none
   real(kind=8)      :: value
   integer(kind = 4) :: one_posnegzero
    
    one_posnegzero = int(sign(one,value))
    if (value .eq. 0.0) one_posnegzero = 0
  
   return
  end function one_posnegzero

end module cresp_p_diffusion
