module cresp_variables ! & constants
! pulled by COSM_RAY_ELECTRONS

  implicit none
  integer          , parameter :: order = 3
  integer          , parameter :: ncre = 5

  type bin_old
    integer                           :: i_lo
    integer                           :: i_up
    real(kind=8), dimension(0:ncre)   :: p
    real(kind=8), dimension(0:ncre)   :: f
    real(kind=8), dimension(1:ncre)   :: q
  end type bin_old

  type(bin_old) crel

  type cresp_vector
    real(kind=8), dimension(1:2*ncre+2) :: cresp_ind
    real(kind=8)   :: uB        ! magnetic energy density
    real(kind=8)   :: uD        ! adiabatic coeff (div_v)
  end type cresp_vector
  
type (cresp_vector) x
!   real(kind=8), dimension(1:2*ncre+1) :: x
  
  real(kind=8)     , parameter :: t_max = 6.0e0 !10d0
  real(kind=8)     , parameter :: f_init  = 1.0e0
  real(kind=8)     , parameter :: q_init  = 5.0
  real(kind=8)                 :: p_min_fix = 1.0e4
  real(kind=8)                 :: p_max_fix = 1.0e5
  real(kind=8)                 :: p_lo = 1.05e4 !p_min_fix  !/10. !* 20.0
  real(kind=8)                 :: p_up = 1.05e5 !p_max_fix ! / 2.0    !9.9d0
  real(kind=8)     , parameter :: dt_ini = 0.1
  real(kind=8)     , parameter :: q_big = 10e0
  real(kind=8)     , parameter :: cfl_cr  = 0.1e0 ! cfl factor for CR

  ! these will most probably be in types and will be modificated by the driver (piernik)



  integer                      :: c2nd, c3rd
!   real(kind=8)                 :: p_lo_next, p_up_next  ! momemntum for spectrum cut-offs
  real(kind=8)                 :: n_tot, n_tot0, e_tot, e_tot0
!   real(kind=8),dimension(1:ncre):: ndt, edt

  ! physical constants
  real(kind=8), parameter      :: cnst_pi = 3.14159265358979311599796346854419e0
  real(kind=8), parameter      :: cnst_c  = 1.0e0 ! speed of light
  !real(kind=8), parameter      :: cnst_me = 1.0d0 ! mass of electron
  
  
  ! used in driver and crspectrum module

!   real(kind=8), dimension(:),allocatable   :: n, e, r

end module cresp_variables
