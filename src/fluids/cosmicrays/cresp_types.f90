 module cresp_types
! pulled by COSM_RAY_ELECTRONS
 
 use vars!, only: nbin
 
 implicit none
 
! -------------
! 
!contains
!
!--------------

 type crsgridcell
  real(kind=8), dimension(1:nbin)   :: q
  real(kind=8), dimension(0:nbin)   :: p,  f !,p_fix0
  integer                           :: i_lo, i_up
end type crsgridcell

type (crsgridcell) crel             ! cosmic ray electrons

! type tosave
!   integer                           :: crel%i_lo
!   integer                           :: crel%i_up
!   real(kind=8), dimension(0:nbin)   :: crel%p
!   real(kind=8), dimension(0:nbin)   :: crel%f
!   real(kidn=8), dimension(1:nbin)   :: crel%q
! end type tosave
end module cresp_types