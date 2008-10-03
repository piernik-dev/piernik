! $Id$
#include "piernik.def"

module constants     ! module containg numerical and physical constants !!!
		     !	   actual units system is CGS		!


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

! uses: length --> cm, mass --> gram, time --> sek, temperature --> kelvin
! length units:
	real, parameter :: cm = 	1.0				! cgs
	real, parameter :: metr = 	1.0e2*cm			! cgs
	real, parameter :: km =		1.0e5*cm
	real, parameter :: au =		1.49597892e13*cm
	real, parameter :: pc =		3.0856e18*cm			! cgs
	real, parameter :: kpc =	1000.0*pc
	real, parameter :: ly =		9.4605e17*cm
! time units:
	real, parameter :: sek =	1.0				! cgs
	real, parameter :: minute =	60.0*sek
	real, parameter :: hour =	3600.0*sek
	real, parameter :: day =	86400.0*sek
	real, parameter :: year =	365.2652*day			! cgs
	real, parameter :: myr =	1.0e6*year			! cgs
! mass units:
	real, parameter :: gram =	1.0				! cgs
	real, parameter :: kg =		1.0e3*gram			! cgs
	real, parameter :: Msun =	1.989e33			! cgs
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
	real, parameter :: newtong =	6.6732e-8*cm**3/gram/sek**2
	real, parameter :: planck =	6.626196e-27*erg*sek
	real, parameter :: r_gc_sun =	8.5*kpc
	real, parameter :: vsun =	220.0*km/sek
	real, parameter :: sunradius =	6.9598e10*cm
	real, parameter :: Lsun =	3.826e33*erg/sek
	real, parameter :: Mearth =	5.977e27*gram
	real, parameter :: earthradius=	6378.17*km
! utilities (to remove in the future):
	real, parameter :: k_B = 1.38066e-16*erg/kelvin ! ---> kboltz
	real, parameter :: m_H = 1.66053e-24*gram       ! ---> amu
	real, parameter :: hydro_mass = mH
	real, parameter :: chcf = hydro_mass / myr
	real, parameter :: fpiG = fpi*newtong
	real, parameter :: sekmyr = myr/sek
	real, parameter :: cmps2 = cm/sek**2

end module constants
