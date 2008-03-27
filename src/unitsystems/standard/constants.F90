#include "mhd.def"

module constants     ! module containg numerical and physical constants !!!

  real, parameter :: one  = 1.0
  real, parameter :: half = 0.5
  real, parameter :: onet = 1./3.
  real, parameter :: twot = 2./3.
  real, parameter :: oneq = 1./4.
  real, parameter :: thrq = 3./4.
  real, parameter :: small= 1.e-50
  real, parameter :: big  = 1.e+50
  
  real, parameter :: pi = 3.141592653589793238
  real, parameter :: dpi = 2.*pi
  real, parameter :: fpi = 4.*pi
  real, parameter :: e  = 2.718281828459045235
  real, parameter :: pc = 3.086e18
  real, parameter :: cmkm = 1.0e5 ! km -> cm 
  real, parameter :: sekmyr = 3.154e13 ! Myr -> s 
  real, parameter :: cmps2 = 3.23e8
  real, parameter :: kpc = 1000.0
  real, parameter :: r_gc_sun = 8500 ! pc 
  real, parameter :: vsun =     220.0  

  real, parameter :: sek =      1.0e-6/365.2652/24.0/3600.0
  real, parameter :: cm =       1.0/3.0856e18
  
  
  real, parameter :: G_cgs = 6.6725985E-8  ! (cgs)
  real, parameter :: k_B = 1.38066e-16 ! erg/K 
  real, parameter :: m_H = 1.66053e-24 ! g 

  real, parameter :: G_one = 1.0  ! (scaled units)
!krix
  real, parameter :: hydro_mass = m_H * cmkm**2
  real, parameter :: chcf = hydro_mass / sekmyr
!krix


  real, parameter :: fpiG  = fpi*G_one
  real, parameter :: newtong = G_one

end module constants
