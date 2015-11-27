module cresp_variables ! & constants
! pulled by COSM_RAY_ELECTRONS
  
  
  integer          , parameter :: cresp_taylor_order = 3 ! TODO - include somewhere in piernik
  integer          , parameter :: ncre_passed = 5
  
  integer                      :: ind_e_beg, ind_e_end, ind_n_beg, ind_n_end, ind_p_lo, ind_p_up
  
  type bin_old
    integer                           :: i_lo
    integer                           :: i_up
    real(kind=8), allocatable, dimension(:)   :: p  !, dimension(0:ncre_passed)   :: p
    
    real(kind=8), allocatable, dimension(:)   :: f  ! dimension(0:ncre_passed)   :: f
    
    real(kind=8), allocatable, dimension(:)   :: q  !dimension(1:ncre_passed)   :: q
    
  end type bin_old
! 
!   type(bin_old) crel
!      allocate(crel%p(0:5))
!      allocate(crel%f(0:5))
!      allocate(crel%q(1:5))
!      
     
  type cresp_vector_template
    real(kind=8), allocatable, dimension(:) :: cresp_ind
    real(kind=8)   :: uB        ! magnetic energy density
    real(kind=8)   :: uD        ! adiabatic coeff (div_v)
  end type cresp_vector_template
  
! type (cresp_vector_template) cresp_vector

  real(kind=8)     , parameter :: q_big = 10e0 ! must be consulted, whether we move it to cresp_crspectrum or make it user defined
!   real(kind=8)     , parameter :: cfl_cr  = 0.1e0 ! cfl factor for CR ! it is not used in the module, just the driver (dt calculation)

  ! these will most probably be in types and will be modificated by the driver (piernik)



  integer                      :: taylor_coeff_2nd, taylor_coeff_3rd

  real(kind=8), parameter      :: cnst_c  = 1.0e0 ! speed of light
  !real(kind=8), parameter      :: cnst_me = 1.0d0 ! mass of electron
  
  
  ! used in driver and crspectrum module

! !   real(kind=8), dimension(:),allocatable   :: n, e, r
contains

!   subroutine init_cresp_types
  

end module cresp_variables
