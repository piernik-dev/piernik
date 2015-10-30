 module cresp_types
! pulled by COSM_RAY_ELECTRONS
 use cresp_variables, only: ncre
 
 implicit none
 
! -------------
! 
!contains
!
!--------------
! 
!  type crsgridcell
!   real(kind=8), dimension(1:ncre)   :: n, e
! !   real(kind=8), dimension(0:ncre)   :: p,  f !,p_fix0
!   real(kind=8)                      :: p_lo, p_up
! end type crsgridcell
! 
! type (crsgridcell) u             ! cosmic ray electrons

type bin_old
  integer                           :: i_lo
  integer                           :: i_up
  real(kind=8), dimension(0:ncre)   :: p
  real(kind=8), dimension(0:ncre)   :: f
  real(kind=8), dimension(1:ncre)   :: q
end type bin_old

type(bin_old)crel

! type cresp_vector
!  real(kind=8), dimension(1:2*ncre+2) :: cresp_ind
! end type cresp_vector
! 
! type (cresp_vector) x

end module cresp_types