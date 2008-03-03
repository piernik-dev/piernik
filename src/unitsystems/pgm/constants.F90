#include "mhd.def"

module constants     ! module containg numerical and physical constants !!!
		     !	actual units system is pc, G_one, Myr	!

! units systems:
!
! cm, g, s
!
! g ---> mean particle mass pm:= 1.36 * mprot = 1.36 * 1.672614e-24 g = 2.2748e-24 g
!
! pc, Msun, myr: pc = 3.0856e18 cm, Msun = 1.989e33 g, myr = 3.154e13 s
! density:    		Msun/pc3 = 	6.7677922e-23 g/cm3
!					29.763 pm/cm3
! column density: 	Msun/pc2 = 	2.0891e-4
!					0.91838e20 pm/cm3
! velocity: 		pc/myr = 	0.97775e5 cm/s = 0.97775 km/s
! energy: 		Msun*pc2/myr2 = 1.9015e43 erg
!
! G_one, pc, myr -> [m] = 1.9713302e21 g
!-----------------> [m]/pc3 = 6.7076688e-35 g/cm3
!-----------------> [energy] = 1.8817e31 erg

  real, parameter :: one  = 1.0
  real, parameter :: half = 0.5
  real, parameter :: onet = 1./3.
  real, parameter :: twot = 2./3.
  real, parameter :: oneq = 1./4.
  real, parameter :: thrq = 3./4.
  real, parameter :: small= 1.e-50
  real, parameter :: big  = 1.e+50
!  real, dimension(2,3)  :: cn
!  data                     cn  /one,thrq,onet, &
!                                one,oneq,twot/

  real, parameter :: pi = 3.141592653589793238
  real, parameter :: dpi = 2.*pi
  real, parameter :: fpi = 4.*pi
  real, parameter :: e  = 2.718281828459045235

! uses: length --> pc, newtong --> 1.0, time --> myr, temperature --> kelvin
! length units:
	real, parameter :: cm = 	1.0/3.0856e18
	real, parameter :: metr = 	1.0e2/3.0856e18
	real, parameter :: km =		1.0e5/3.0856e18
	real, parameter :: au =		1.49597892e13/3.0856e18
	real, parameter :: pc =		1.0
	real, parameter :: kpc =	1000.0*pc
	real, parameter :: ly =		9.4605e17/3.0856e18
! time units:
	real, parameter :: sek =	1.0e-6/365.2652/24.0/3600.0
	real, parameter :: minute =	1.0e-6/365.2652/24.0/60.0
	real, parameter :: hour =	1.0e-6/365.2652/24.0
	real, parameter :: day =	1.0e-6/365.2652
	real, parameter :: year =	1.0e-6
	real, parameter :: myr =	1.0
! mass units: (unusually here)
	real, parameter :: G_one =	1.0 !this is not a mass unit, nevertheless it is usefull here
	real, parameter :: gram =	6.6732e-8*cm**3/G_one/sek**2
	real, parameter :: kg =		1.0e3
	real, parameter :: Msun =	1.989e33*gram
	real, parameter :: gmu =	2.32e7*Msun
	real, parameter :: me = 	9.109558e-28*gram
	real, parameter :: mp = 	1.672614e-24*gram
	real, parameter :: mH = 	1.673559e-24*gram
	real, parameter :: amu =	1.660531e-24*gram
! force units:
	real, parameter :: newton =	kg*metr/sek**2
	real, parameter :: dyna =	gram*cm/sek**2
! energy units:
	real, parameter :: joul =	kg*metr**2/sek**2
	real, parameter :: erg =	gram*cm**2/sek**2
	real, parameter :: eV   =       1.6022e-12*erg	
! density units:
	real, parameter :: ppcm3 =	1.36 * mp / cm**3
	real, parameter :: ppcm2 =	1.36 * mp / cm**2
! temperature units:
	real, parameter :: kelvin =	1.0
! physical constants:
	real, parameter :: kboltz =	1.380622e-16*erg/kelvin
	real, parameter :: gasRconst =	8.31434e7*erg/kelvin !/mol
	real, parameter :: clight =	2.997924562e10*cm/sek
	real, parameter :: newtong =	G_one
	real, parameter :: planck =	6.626196e-27*erg*sek
	real, parameter :: r_gc_sun =	8.5*kpc
	real, parameter :: vsun =	220.0*km/sek
	real, parameter :: sunradius =	6.9598e10*cm
	real, parameter :: Lsun =	3.826e33*erg/sek
	real, parameter :: Mearth =	5.977e27*gram
	real, parameter :: earthradius=	6378.17*km
! utilities (to remove in the future):
	real, parameter :: k_B = 1.38066e-16*erg/kelvin ! erg/K ---> kboltz
	real, parameter :: m_H = 1.66053e-24*gram ! g     ---> amu
	real, parameter :: hydro_mass = mH
	real, parameter :: chcf = hydro_mass / myr
	real, parameter :: fpiG = fpi*newtong
	real, parameter :: sekmyr = myr/sek
	real, parameter :: cmps2 = cm/sek**2

end module constants
