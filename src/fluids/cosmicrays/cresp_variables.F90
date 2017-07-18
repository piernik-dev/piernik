module cresp_variables ! & constants
! pulled by COSM_RAY_ELECTRONS

  implicit none
!   integer          , parameter :: order = 1

  real(kind=8)     , parameter :: div_v= 0.0
  real(kind=8)     , parameter :: omega_d = 0.0
  real(kind=8)     , parameter :: clight = 1.0

  ! these will most probably be in types and will be modified by the driver (piernik)
! 
!   integer                      :: taylor_coeff_2nd, taylor_coeff_3rd

!   real(kind=8), parameter      :: cnst_c  = 1.0e0 ! speed of light
!   real(kind=8)    :: p_lo_d, p_up_d, u_d_d, u_b_d   ! additional variables that are passed down to timestep by grid, but since grid is absent in driver module, they must be compensated for this way.
  !real(kind=8), parameter      :: cnst_me = 1.0d0 ! mass of electron
  
  
  ! used in driver and crspectrum module

!   real(kind=8), dimension(:),allocatable   :: n, e, r

end module cresp_variables
