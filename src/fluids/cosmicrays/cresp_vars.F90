module vars ! & constants
! pulled by NO_FLAG
  implicit none

  integer          , parameter :: order = 3
  integer          , parameter :: ncre = 5

  real(kind=8)     , parameter :: t_max = 6.0d0 !10d0
  real(kind=8)     , parameter :: f_init  = 1.0d0
  real(kind=8)     , parameter :: q_init  = 5.0
  real(kind=8)                 :: p_min_fix = 1.0d4
  real(kind=8)                 :: p_max_fix = 1.0d5
  real(kind=8)                 :: p_lo = 1.05d4 !p_min_fix  !/10. !* 20.0
  real(kind=8)                 :: p_up = 1.05d5 !p_max_fix ! / 2.0    !9.9d0
  real(kind=8)                 :: w
  real(kind=8)     , parameter :: dt_ini = 0.1
  real(kind=8)     , parameter :: q_big = 10d0
  real(kind=8)     , parameter      :: cfl_cr  = 0.1d0 ! cfl factor for CR

  ! these will most probably be in types and will be modificated by the driver (piernik)
  real(kind=8)                 :: u_d0 = -2.5d-1 ! 5.0d-1
  real(kind=8)                 :: u_d
  real(kind=8)                 :: u_b = 0.d0 ! 1d-7 !0d0 !5d-7!
  real(kind=8)                 :: div_v = 0d0
  real(kind=8)                 :: omega_d = 0.5d0 !0.1d0    ! frequency of div(v) oscilations



  integer                      :: c2nd, c3rd
!   real(kind=8)                 :: p_lo_next, p_up_next  ! momemntum for spectrum cut-offs
  real(kind=8)                 :: n_tot, n_tot0, e_tot, e_tot0
!   real(kind=8),dimension(1:ncre):: ndt, edt

  ! physical constants
  real(kind=8), parameter      :: cnst_pi = 3.14159265358979311599796346854419d0
  real(kind=8), parameter      :: cnst_c  = 1.0d0 ! speed of light
  !real(kind=8), parameter      :: cnst_me = 1.0d0 ! mass of electron
  
  
  ! used in driver and crspectrum module

!   real(kind=8), dimension(:),allocatable   :: n, e, r

end module vars

