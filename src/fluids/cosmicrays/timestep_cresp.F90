module timestep_cresp
!pulled by COSM_RAY_ELECTRONS

 use initcrspectrum, only: ncre, cfl_cre
 use constants,      only: one, zero

      implicit none
      
!       public          dt_cre      
           
      private         p_tmp
 
      real(kind=8) :: dt_cre, dt_tmp
      real(kind=8), allocatable, dimension(:) :: p_tmp   ! do not use p variable in this module anywhere outside
      real(kind=8), allocatable, dimension(:) :: dts_new

!-------------------      
!      
contains 
!
!-------------------

!----------------------------------------------------------
! Subroutine consistent with rules depicted in crspectrum.pdf

 subroutine cresp_timestep_new(dt_comp, p_l, p_u, sptab)
!    use initcosmicrays,   only: p_fix, ncre
   use initcrspectrum,   only: spec_mod_trms, p_fix, ncre
   use cresp_arrays_handling, only: allocate_with_index
   implicit none
   
   real(kind=8)               :: dt_comp
   real(kind=8)               :: dt_cre_ud, dt_cre_ub
   real(kind=8) :: p_l, p_u
   type(spec_mod_trms) sptab

  call allocate_with_index(p_tmp,0,ncre)
  call allocate_with_index(dts_new,1,ncre)
   
   dt_cre_ud = huge(one)
   dt_cre_ub = huge(one)
   dts_new = huge(one)
   
! Adiabatic cooling timestep:
   if (sptab%ud .ne. zero) then
     where (p_fix(1:ncre) .gt. zero .and. p_fix(0:ncre-1) .gt. zero)
       dts_new = cfl_cre * log10(p_fix(1:ncre) / p_fix(0:ncre-1)) / sptab%ud
       dts_new = abs(dts_new)
     end where
   endif
   dt_cre_ud = minval(dts_new)   ! it was already multiplied by cfl_cre
   
   dts_new = huge(one)
   
! Synchrotron cooling timestep (is dependant only on p_up, highest value of p):
   if (sptab%ub .gt. zero) then
     where (p_fix(1:ncre).gt.zero .and. p_fix(0:ncre-1) .gt. zero)   ! We only need p_(i+1/2) to compute this timestep.
       dts_new = cfl_cre * log10(p_fix(1:ncre)/p_fix(0:ncre-1)) / (p_u * sptab%ub) 
     end where
   endif
   dt_cre_ub = minval(dts_new)
#ifdef VERBOSE
   print *, '[@timestep_cresp:] Computed timesteps:'
   print *, 'dt_cre_ud = ', dt_cre_ud, 'dt_cre_ub = ', dt_cre_ub
#endif /* VERBOSE */

! Here shortest among calculated timesteps is chosen.
   dt_comp = min(dt_cre_ud, dt_cre_ub)
   
! Should dt_cre_ud or dt_cre_ub be greater than current one, next timestep shall be independently increased by piernik in timestep
   dt_cre = min(dt_cre, dt_comp) ! Gives minimal timestep among computed for current cell and previous ones
   
    if (allocated(p_tmp))   deallocate(p_tmp)
    if (allocated(dts_new)) deallocate(dts_new)
    
 end subroutine cresp_timestep_new

!----------------------------------------------------------
 
   subroutine cresp_timestep(dt_comp, u_b)    ! cresp_timestep_old
    use initcrspectrum,    only: cfl_cre
    implicit none
    real(kind=8)                  :: dt_comp
    real(kind=8)                  :: dts_min
    real(kind=8)                  :: u_b
      
      dt_comp = huge(one)
      dts_new = huge(one) ! whole dts_new array
      where (abs(b_losses(p_tmp(1:ncre), u_b)) .ne. zero)
        dts_new =  (p_tmp(1:ncre)-p_tmp(0:ncre-1))/abs(b_losses(p_tmp(1:ncre), u_b))
        where ((p_tmp(1:ncre)-p_tmp(0:ncre-1)).eq.zero)   !!!
          dts_new = huge(one) 
        end where
      end where
      
      dts_min = cfl_cre*minval(dts_new)   ! min of array
      
      if (dt_comp.ge.dts_min) then   ! gives minimal timestep in bin space
         dt_comp = dts_min
      endif
      
      dts_new = zero
      
   end subroutine cresp_timestep

!----------------------------------------------------------
   
  function b_losses(p, u_b)
    implicit none
    real(kind=8), intent(in)                :: u_b
    real(kind=8), dimension(:), intent(in)  :: p
    real(kind=8), dimension(size(p)) :: b_losses
   
    b_losses = u_b*p**2  !!! b_sync_ic = 8.94e-25*(u_b+u_cmb)*gamma_l**2 ! erg/cm

  end function
  
end module timestep_cresp
