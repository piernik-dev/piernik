module initcrspectrum  
! pulled by COSM_RAY_ELECTRONS

! contains routines reading namelist in problem.par file dedicated to cosmic ray electron spectrum and initializes types used.

! present in namelist COSMIC_RAY_SPECTRUM
 integer(kind=4)    :: ncre !       = 5
 real(kind=8)       :: p_min_fix ! = 1.0e2
 real(kind=8)       :: p_max_fix ! = 1.0e5
 real(kind=8)       :: f_init    ! = 1.0
 real(kind=8)       :: q_init    ! = 5.0
 real(kind=8)       :: q_big     ! maximal amplitude of q
 real(kind=8)       :: p_lo_init ! = 1.2e1
 real(kind=8)       :: p_up_init ! = 1.3e3
 real(kind=8)       :: cfl_cre   ! = 1.0
 real(kind=8)       :: cre_eff
 real(kind=8)       :: K_cre_e_paral ! = 0
 real(kind=8)       :: K_cre_e_perp ! = 0
 real(kind=8)       :: K_cre_n_paral ! = 0
 real(kind=8)       :: K_cre_n_perp ! = 0
 real(kind=8)       :: K_pow_index  ! = 0 
 integer(kind=4)    :: expan_order
 
 real(kind=8),allocatable, dimension(:) :: p_fix
 real(kind=8)       :: w
!-----------
! Types used in module: 
   type bin_old
    integer                           :: i_lo
    integer                           :: i_up
    real(kind=8),allocatable,dimension(:)   :: p
    real(kind=8),allocatable,dimension(:)   :: f
    real(kind=8),allocatable,dimension(:)   :: q
  end type bin_old

  type(bin_old) crel

  type spec_mod_trms
   real(kind=8)    :: ub
   real(kind=8)    :: ud
   real(kind=8)    :: ucmb
  end type spec_mod_trms
      
  integer(kind=4)  :: taylor_coeff_2nd, taylor_coeff_3rd
  
  real, allocatable :: fdif_cre(:,:,:,:)    ! array of diffusion fluxes for p_cut diffusion
  real, allocatable :: fadv_cre(:,:,:,:)    ! array of advection fluxes for p_cut diffusion
  
  integer,allocatable, dimension(:) :: cresp_edges
!----------------------------------------------  
!  
contains
!
!----------------------------------------------

  subroutine init_cresp
  use constants,             only: I_ZERO, I_ONE, zero, one !, ten
  use cresp_arrays_handling, only: allocate_with_index
  use domain,  only: dom

  implicit none
   
   integer                           :: vecdim  ! vector dimension
!    real(kind=8)                    :: w
   integer                           :: i       ! enumerator
 
  namelist /COSMIC_RAY_SPECTRUM/ cfl_cre, p_lo_init, p_up_init, f_init, q_init, q_big, ncre, &
           &                         p_min_fix, p_max_fix, cre_eff, K_cre_e_paral, K_cre_e_perp, &
           &                         K_cre_n_paral, K_cre_n_perp, K_pow_index, expan_order
           
  open(unit=101, file="problem.par", status="unknown")
  read(unit=101, nml=COSMIC_RAY_SPECTRUM)
  close(unit=101)
  
  if (ncre .ne. I_ZERO)  then
   print *,'[@init_cresp] Initial CRESP parameters read:'
   print *,'ncre        = ', ncre
   print *,'p_min_fix   = ', p_min_fix
   print *,'p_max_fix   = ', p_max_fix
   print *,'p_lo_init   = ', p_lo_init
   print *,'p_up_init   = ', p_up_init
   print *,'f_init      = ', f_init
   print *,'q_init      = ', q_init
   print *,'q_big       = ', q_big
   print *,'cfl_cre     = ', cfl_cre
   print *,'expan_order = ', expan_order
   print *,'K_cre_e_paral=', K_cre_e_paral
   print *,'K_cre_e_perp=', K_cre_e_perp
   print *,'K_cre_n_paral=', K_cre_n_paral
   print *,'K_cre_n_perp=', K_cre_n_perp
   
! arrays initialization and stuff

  vecdim = ncre
  call allocate_with_index(p_fix,0,vecdim)
  call allocate_with_index(cresp_edges ,0,vecdim)
  
  cresp_edges = (/ (i,i=0,ncre) /)
  p_fix = zero 
    w  = (log10(p_max_fix/p_min_fix))/dble(ncre-2)
    p_fix(1:ncre-1)  =  p_min_fix*10.0**(w*dble(cresp_edges(1:ncre-1)-1))
    p_fix(0)    = zero
    p_fix(ncre) = zero
  print *, 'p_fix = ', p_fix
   
 else
  print *,'ncre   = ', ncre, '; cr-electrons NOT initnialized. If COSM_RAY_ELECTRONS flag is on, please check your parameters.'
 endif

  taylor_coeff_2nd = (mod(2,expan_order) / 2 + mod(3,expan_order))       ! coefficient which is always equal to 1 when order = 2 or = 3 and 0 if order = 1
  taylor_coeff_3rd = (expan_order - 1)*(expan_order- 2) / 2              ! coefficient which is equal to 1 only when order = 3
 
 call init_cresp_types
!  print *,'[@init_cresp:] type crel%[p,f,q,i_lo i_up] :'
!  print *, crel%p
!  print *, crel%f
!  print *, crel%q
!  print *, crel%i_lo, crel%i_up
 
  allocate(fdif_cre(ncre,0:dom%n_d(1)-I_ONE,0:dom%n_d(2)-I_ONE,0:dom%n_d(3)-I_ONE))
  allocate(fadv_cre(ncre,0:dom%n_d(1)-I_ONE,0:dom%n_d(2)-I_ONE,0:dom%n_d(3)-I_ONE))
 
 end subroutine init_cresp
 
!----------------------------

 subroutine init_cresp_types
   use cresp_arrays_handling,   only: allocate_with_index
   use constants,               only: zero
   implicit none

   call allocate_with_index(crel%p,0,ncre)
   call allocate_with_index(crel%f,0,ncre)
   call allocate_with_index(crel%q,1,ncre)
   
   crel%p = zero
   crel%q = zero
   crel%f = zero
   crel%i_lo = zero
   crel%i_up = zero
   
 end subroutine init_cresp_types
 
 function compute_K(K0, alpha, alpha_n, length)
 implicit none
 real :: K0, alpha, alpha_n
 integer :: length, i
 real, dimension(length) :: compute_K
 
 do i = 1, length+1
    compute_K(i)   = K0*(p_min_fix*10.0**(w*dble(cresp_edges(i-1)-1))/p_max_fix) ** (alpha - alpha_n)
 enddo 
 print *, '@compute_K = ', compute_K

 end function compute_K
end module initcrspectrum