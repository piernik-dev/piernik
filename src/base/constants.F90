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

#ifdef PSM
! uses: length --> pc, mass --> Msun, time --> myr, temperature --> kelvin
! length units:
	real, parameter :: cm = 	1.0/3.0856e18			! psm
	real, parameter :: metr = 	1.0e2*cm			! psm
	real, parameter :: pc =		1.0				! psm
! time units:
	real, parameter :: sek =	1.0e-6/365.2652/24.0/3600.0	! psm
	real, parameter :: year =	1.0e-6				! psm
	real, parameter :: myr =	1.0				! psm
! mass units:
	real, parameter :: gram =	1.0/1.989e33			! psm
	real, parameter :: kg =		1.0e3*gram			! psm
	real, parameter :: Msun =	1.0				! psm

#elif defined (PGM)
! uses: length --> pc, newtong --> 1.0, time --> myr, temperature --> kelvin
! length units:
	real, parameter :: cm = 	1.0/3.0856e18			! pgm
	real, parameter :: metr = 	1.0e2*cm			! pgm
	real, parameter :: pc =		1.0				! pgm
! time units:
	real, parameter :: sek =	1.0e-6/365.2652/24.0/3600.0	! pgm
	real, parameter :: year =	1.0e-6				! pgm
	real, parameter :: myr =	1.0				! pgm
! mass units:
	real, parameter :: G_one =	1.0 !this is not a mass unit, nevertheless it is usefull here
	real, parameter :: gram =	6.6732e-8*cm**3/G_one/sek**2	! pgm
	real, parameter :: kg =		1.0e3*gram			! pgm
	real, parameter :: Msun =	1.989e33*gram			! pgm

#elif defined (DMY)
! uses: length --> cm^16, mass --> Msun, time --> year, temperature --> kelvin
! length units:
	real, parameter :: cm = 	1.0e-16				! dmy
	real, parameter :: metr = 	1.0e2*cm			! dmy
	real, parameter :: pc =		3.0856e18*cm			! dmy
! time units:
	real, parameter :: sek =	1.0/365.2652/86400.0		! dmy
	real, parameter :: year =	1.0				! dmy
	real, parameter :: myr =	1.0e6*year			! dmy
! mass units:
	real, parameter :: gram =	1.0/1.989e33			! dmy
	real, parameter :: kg =		1.0e3*gram			! dmy
	real, parameter :: Msun =	1.0				! dmy

#elif defined(SI)
! uses: length --> metr, mass --> kg, time --> sek, temperature --> kelvin
! length units:
	real, parameter :: cm = 	1.0e-2				! SI
	real, parameter :: metr = 	1.0				! SI
	real, parameter :: pc =		3.0856e18*cm			! SI
! time units:
	real, parameter :: sek =	1.0				! SI
	real, parameter :: year =	365.2652*day			! SI
	real, parameter :: myr =	1.0e6*year			! SI
! mass units:
	real, parameter :: gram =	1.0e-3				! SI
	real, parameter :: kg =		1.0				! SI
	real, parameter :: Msun =	1.989e30			! SI

#elif defined(CGS)
! uses: length --> cm, mass --> gram, time --> sek, temperature --> kelvin
! length units:
	real, parameter :: cm = 	1.0				! cgs
	real, parameter :: metr = 	1.0e2*cm			! cgs
	real, parameter :: pc =		3.0856e18*cm			! cgs
! time units:
	real, parameter :: sek =	1.0				! cgs
	real, parameter :: year =	365.2652*day			! cgs
	real, parameter :: myr =	1.0e6*year			! cgs
! mass units:
	real, parameter :: gram =	1.0				! cgs
	real, parameter :: kg =		1.0e3*gram			! cgs
	real, parameter :: Msun =	1.989e33			! cgs

#endif /* PSM, PGM, DMY, SI, CGS */

#ifndef STANDARD
! length units:
	real, parameter :: km =		1.0e5*cm
	real, parameter :: au =		1.49597892e13*cm
	real, parameter :: kpc =	1000.0*pc
	real, parameter :: ly =		9.4605e17*cm
! time units:
	real, parameter :: minute =	60.0*sek
	real, parameter :: hour =	3600.0*sek
	real, parameter :: day =	86400.0*sek
! mass units:
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
#ifdef PGM
	real, parameter :: newtong =	G_one				! pgm
#else /* PGM */
	real, parameter :: newtong =	6.6732e-8*cm**3/gram/sek**2
#endif /* PGM */
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

#else /* STANDARD */
  real, parameter :: sek =      1.0e-6/365.2652/24.0/3600.0
  real, parameter :: cm =       1.0/3.0856e18
  
  
  real, parameter :: G_cgs = 6.6725985E-8  ! (cgs)
  real, parameter :: k_B = 1.38066e-16 ! erg/K 
  real, parameter :: m_H = 1.66053e-24 ! g 

  real, parameter :: G_one = 1.0  ! (scaled units)


  real, parameter :: fpiG  = fpi*G_one
  real, parameter :: newtong = G_one
  real, parameter :: cmps2 = 3.23e8
  real, parameter :: pc = 3.086e18
  real, parameter :: cmkm = 1.0e5      ! km -> cm 
  real, parameter :: sekmyr = 3.154e13 ! Myr -> s 
  real, parameter :: kpc = 1000.0
  real, parameter :: r_gc_sun = 8500   ! pc 

  real, parameter :: vsun =     220.0  
  real, parameter :: hydro_mass = m_H * cmkm**2
  real, parameter :: chcf = hydro_mass / sekmyr
#endif /* STANDARD */

end module constants
