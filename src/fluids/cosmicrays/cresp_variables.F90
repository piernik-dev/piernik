module cresp_variables ! & constants
! pulled by COSM_RAY_ELECTRONS

   use constants, only: fpi, three

   implicit none

   public                                                ! QA_WARN no secrets are kept here

   real, parameter :: div_v        = 0.0
   real, parameter :: omega_d      = 0.0
   real, parameter :: clight_cresp = 1.0

   real, parameter :: fpcc         = fpi * clight_cresp
   real, parameter :: fpcc2        = fpi * clight_cresp**2
   real, parameter :: fp3cc        = fpi / three * clight_cresp

! these will most probably be in types and will be modified by the driver (piernik)
!
!   integer              :: taylor_coeff_2nd, taylor_coeff_3rd

!   real, parameter      :: cnst_c  = 1.0e0 ! speed of light
!   real                 :: p_lo_d, p_up_d, u_d_d, u_b_d   ! additional variables that are passed down to timestep by grid, but since grid is absent in driver module, they must be compensated for this way.
!   real, parameter      :: cnst_me = 1.0d0 ! mass of electron


end module cresp_variables
