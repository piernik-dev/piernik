! $Id$
!
! PIERNIK Code Copyright (C) 2006 Michal Hanasz
!
!    This file is part of PIERNIK code.
!
!    PIERNIK is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    PIERNIK is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with PIERNIK.  If not, see <http://www.gnu.org/licenses/>.
!
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.def"

!>
!! \brief R [DW] Module containing numerical and physical %constants
!! \details Module constants contains numerical and physical %constatns for several units systems.
!! To use one system a proper precompiler directive should be defined.
!! Available units systems:
!!
!! @b PSM (Parsec - Solar_mass - Megayear) - good for global galactic simulations
!! @n length --> pc,     mass --> Msun,        time --> myr,        miu0 --> 4*pi,    temperature --> kelvin
!!
!! @b PLN (PLaNetary) - good for planetary nebulae
!! @n length --> AU,     mass --> Mjup,        time --> yr,         miu0 --> 4*pi,    temperature --> kelvin
!!
!! @b KSG (Kiloparsec - Solar_mass - Gigayear) - good for galactic and intergalactic simulations
!! @n length --> kpc,    mass --> 10^6*Msun,   time --> Gyr,        miu0 --> 4*pi,    temperature --> kelvin
!!
!! @b PGM (Parsec - Gravity=1 - Megayear) - modification of PSM system, where %gravity constant is one
!! @n length --> pc,     newtong --> 1.0,      time --> myr,        miu0 --> 4*pi,    temperature --> kelvin
!!
!! @b SSY (10^Sixteenth cm - Solar_mass - Year) - good for circumstellar simulations (planetaries etc.)
!! @n length --> 10^16 cm,  mass --> Msun,     time --> year,       miu0 --> 4*pi,    temperature --> kelvin
!!
!! @b SI - (franc. Système International d'Unités)
!! @n length --> metr,   mass --> kg,          time --> sek,        miu0 --> 4*pi,    temperature --> kelvin
!!
!! @b CGS - (Centimetre - Gram - Second)
!! @n length --> cm,     mass --> gram,        time --> sek,        miu0 --> 4*pi,    temperature --> kelvin
!!
!! @b SCALED - a suit of %constants without physical units. This is automatically set while neither of the former systems is chosen.
!<
module constants

   real, parameter :: one        = 1.0                   !< one
   real, parameter :: half       = 0.5                   !< a half
   real, parameter :: onet       = 1./3.                 !< one third
   real, parameter :: twot       = 2./3.                 !< two thirds
   real, parameter :: oneq       = 1./4.                 !< one fourth
   real, parameter :: thrq       = 3./4.                 !< three fourths
   real, parameter :: small      = 1.e-29                !< a constant used as the lower limit number
   real, parameter :: big        = 1.e+29                !< a constant used as the upper limit number

   real, parameter :: pi         = 3.141592653589793238  !< Pi (Archimedes' constant)
   real, parameter :: dpi        = 2.*pi                 !< doubled Pi
   real, parameter :: fpi        = 4.*pi                 !< four Pi
   real, parameter :: e          = 2.718281828459045235  !< Napier's constant (base of Natural logarithm)

#ifdef PSM
! PSM  uses: length --> pc,     mass --> Msun,        time --> myr,        miu0 --> 4*pi,    temperature --> kelvin
! length units:
   real, parameter :: cm         = 1.0/3.0856e18         !< centimetre, length unit
   real, parameter :: metr       = 1.0e2*cm              !< metre, length unit
   real, parameter :: pc         = 1.0                   !< parsec, length unit
! time units:
   real, parameter :: sek        = 1.0e-6/365.2652/24.0/3600.0       !< second, time unit
   real, parameter :: year       = 1.0e-6                !< year, time unit
   real, parameter :: myr        = 1.0                   !< megayear, time unit
! mass units:
   real, parameter :: gram       = 1.0/1.989e33          !< gram, mass unit
   real, parameter :: kg         = 1.0e3*gram            !< kilogram, mass unit
   real, parameter :: Msun       = 1.0                   !< mass of Sun

#elif defined (PLN)
! PLN  uses: length --> AU,     mass --> Mjup,        time --> yr,         miu0 --> 4*pi,    temperature --> kelvin
! length units:
   real, parameter :: cm         = 1.0/1.49598e13        !< centimetre, length unit
   real, parameter :: metr       = 1.0e2*cm              !< metre, length unit
   real, parameter :: AU         = 1.0                   !< astronomical unit
   real, parameter :: pc         = 206264.806*AU         !< parsec, length unit
! time units:
   real, parameter :: sek        = 1.0/365.2652/24.0/3600.0       !< second, time unit
   real, parameter :: year       = 1.0                   !< year, time unit
   real, parameter :: myr        = 1.0e6                 !< megayear, time unit
! mass units:
   real, parameter :: gram       = 1.0/1.8986e30         !< gram, mass unit
   real, parameter :: kg         = 1.0e3*gram            !< kilogram, mass unit
   real, parameter :: Msun       = 1047.4                !< mass of Sun


#elif defined (KSG)
! KSG  uses: length --> kpc,    mass --> 10^6*Msun,   time --> Gyr,        miu0 --> 4*pi,    temperature --> kelvin
! length units:
   real, parameter :: cm         = 1.0/3.0856e21         !< centimetre, length unit
   real, parameter :: metr       = 1.0e2*cm              !< metre, length unit
   real, parameter :: pc         = 1.0e-3                !< parsec, length unit
! time units:
   real, parameter :: sek        = 1.0e-9/365.2652/24.0/3600.0       !< second, time unit
   real, parameter :: year       = 1.0e-9                !< year, time unit
   real, parameter :: myr        = 1.0e-3                !< megayear, time unit
! mass units:
   real, parameter :: gram       = 1.0e-6/1.989e33       !< gram, mass unit
   real, parameter :: kg         = 1.0e3*gram            !< kilogram, mass unit
   real, parameter :: Msun       = 1.0e-6                !< mass of Sun

#elif defined (PGM)
! PGM  uses: length --> pc,     newtong --> 1.0,      time --> myr,        miu0 --> 4*pi,    temperature --> kelvin
! length units:
   real, parameter :: cm         = 1.0/3.0856e18            !< centimetre, length unit
   real, parameter :: metr       = 1.0e2*cm                 !< metre, length unit
   real, parameter :: pc         = 1.0                      !< parsec, length unit
! time units:
   real, parameter :: sek        = 1.0e-6/365.2652/24.0/3600.0    !< second, time unit
   real, parameter :: year       = 1.0e-6                   !< year, time unit
   real, parameter :: myr        = 1.0                      !< megayear, time unit
! mass units:
   real, parameter :: G_one      = 1.0 !this is not a mass unit, nevertheless it is usefull to be set here
   real, parameter :: gram       = 6.6732e-8*cm**3/G_one/sek**2      !< gram, mass unit
   real, parameter :: kg         = 1.0e3*gram               !< kilogram, mass unit
   real, parameter :: Msun       = 1.989e33*gram            !< mass of Sun

#elif defined (SSY)
! SSY  uses: length --> 10^16 cm,  mass --> Msun,     time --> year,       miu0 --> 4*pi,    temperature --> kelvin
! length units:
   real, parameter :: cm         = 1.0e-16                  !< centimetre, length unit
   real, parameter :: metr       = 1.0e2*cm                 !< metre, length unit
   real, parameter :: pc         = 3.0856e18*cm             !< parsec, length unit
! time units:
   real, parameter :: sek        = 1.0/365.2652/86400.0     !< second, time unit
   real, parameter :: year       = 1.0                      !< year, time unit
   real, parameter :: myr        = 1.0e6*year               !< megayear, time unit
! mass units:
   real, parameter :: gram       = 1.0/1.989e33             !< gram, mass unit
   real, parameter :: kg         = 1.0e3*gram               !< kilogram, mass unit
   real, parameter :: Msun       = 1.0                      !< mass of Sun

#elif defined(SI)
! SI   uses: length --> metr,   mass --> kg,          time --> sek,        miu0 --> 4*pi,    temperature --> kelvin
! length units:
   real, parameter :: cm         = 1.0e-2                   !< centimetre, length unit
   real, parameter :: metr       = 1.0                      !< metre, length unit
   real, parameter :: pc         = 3.0856e18*cm             !< parsec, length unit
! time units:
   real, parameter :: sek        = 1.0                      !< second, time unit
   real, parameter :: year       = 365.2652*24.*3600.*sek   !< year, time unit
   real, parameter :: myr        = 1.0e6*year               !< megayear, time unit
! mass units:
   real, parameter :: gram       = 1.0e-3                   !< gram, mass unit
   real, parameter :: kg         = 1.0                      !< kilogram, mass unit
   real, parameter :: Msun       = 1.989e30                 !< mass of Sun

#elif defined(CGS)
! CGS  uses: length --> cm,     mass --> gram,        time --> sek,        miu0 --> 4*pi,    temperature --> kelvin
! length units:
   real, parameter :: cm         = 1.0                      !< centimetre, length unit
   real, parameter :: metr       = 1.0e2*cm                 !< metre, length unit
   real, parameter :: pc         = 3.0856e18*cm             !< parsec, length unit
! time units:
   real, parameter :: sek        = 1.0                      !< second, time unit
   real, parameter :: year       = 365.2652*24.*3600.*sek   !< year, time unit
   real, parameter :: myr        = 1.0e6*year               !< megayear, time unit
! mass units:
   real, parameter :: gram       = 1.0                      !< gram, mass unit
   real, parameter :: kg         = 1.0e3*gram               !< kilogram, mass unit
   real, parameter :: Msun       = 1.989e33                 !< mass of Sun
#else
#define SCALED
#endif /* PSM, KSG, PGM, SSY, SI, CGS */

#ifndef SCALED
! length units:
   real, parameter :: km         = 1.0e5*cm                 !< kilometer, length unit
#ifndef PLN
   real, parameter :: au         = 1.49597892e13*cm         !< astronomical unit (length unit)
#endif /* PLN */
   real, parameter :: kpc        = 1000.0*pc                !< kiloparsec, length unit
   real, parameter :: lyr        = 9.4605e17*cm             !< light year, length unit
! time units:
   real, parameter :: minute     = 60.0*sek                 !< minute, time unit
   real, parameter :: hour       = 3600.0*sek               !< hour, time unit
   real, parameter :: day        = 86400.0*sek              !< day, time unit
! mass units:
   real, parameter :: gmu        = 2.32e7*Msun              !< galactic mass unit
   real, parameter :: me         = 9.109558e-28*gram        !< electron mass
   real, parameter :: mp         = 1.672614e-24*gram        !< proton mass
   real, parameter :: mH         = 1.673559e-24*gram        !< hydrogen atom mass
   real, parameter :: amu        = 1.660531e-24*gram        !< atomic mass unit
! force units:
   real, parameter :: newton     = kg*metr/sek**2           !< 1N (SI force unit)
   real, parameter :: dyna       = gram*cm/sek**2           !< 1 dyna (cgs force unit)
! energy units:
   real, parameter :: joul       = kg*metr**2/sek**2        !< 1J (SI energy unit)
   real, parameter :: erg        = gram*cm**2/sek**2        !< 1 erg (cgs energy unit)
   real, parameter :: eV         = 1.6022e-12*erg           !< 1 eV
! density units:
   real, parameter :: ppcm3      = 1.36 * mp / cm**3        !< spatial density unit
   real, parameter :: ppcm2      = 1.36 * mp / cm**2        !< column density unit
! temperature units:
   real, parameter :: kelvin     = 1.0                      !< kelvin, temperature unit
! physical constants:
   real, parameter :: kboltz     = 1.3806504e-16*erg/kelvin !< Boltzmann constant
   real, parameter :: gasRconst  = 8.314472e7*erg/kelvin    !< gas constant R =  8.314472e7*erg/kelvin/mol
   real, parameter :: clight     = 2.997924562e10*cm/sek    !< speed of light in vacuum

   real, parameter :: Gs         = sqrt(4.*pi*gram/cm)/sek  !< 1 Gs (cgs magnetic induction unit)
   real, parameter :: mGs        = Gs*1.e-6                 !< 1 microgauss
   real, parameter :: Tesla      = 1.e4*Gs                  !< 1 T (SI magnetic induction unit)
#ifdef PGM
   real, parameter :: newtong    = G_one                    !< Newtonian constant of gravitation (equal to G_one while PGM defined)
#else /* PGM */
   real, parameter :: newtong    = 6.6732e-8*cm**3/gram/sek**2 !< Newtonian constant of gravitation
#endif /* PGM */
   real, parameter :: fpiG       = fpi*newtong              !< four Pi times Newtonian constant of gravitation (commonly used in self-gravity routines)
   real, parameter :: planck     = 6.626196e-27*erg*sek     !< Planck constant
   real, parameter :: r_gc_sun   = 8.5*kpc                  !< Sun distance from the Galaxy Center
   real, parameter :: vsun       = 220.0*km/sek             !< velocity value of Sun in the Galaxy
   real, parameter :: sunradius  = 6.9598e10*cm             !< radius of Sun
   real, parameter :: Lsun       = 3.826e33*erg/sek         !< luminosity of Sun
   real, parameter :: Mearth     = 5.977e27*gram            !< mass of Earth
   real, parameter :: earthradius= 6378.17*km               !< radius of Earth

#else /* SCALED */
!>
!! \todo to check validity of declaration in SCALED units following constants: sek, cm, pc, kpc (e.g. to find coincidence with another unit system)
!<
   real, parameter :: sek        = 1.0e-6/365.2652/24.0/3600.0     !< second expressed in units of megayears
   real, parameter :: cm         = 1.0/3.0856e18            !< centimetre expressed in units of parsecs

   real, parameter :: kboltz     = 1.3806504e-16            !< Boltzmann constant (in erg/kelvin while SCALED)
   real, parameter :: gasRconst  = 1.0                      !< gas constant (scaled)
   real, parameter :: amu        = 1.660531e-24             !< atomic mass unit (in grams while SCALED)
   real, parameter :: mH         = 1.0                      !< hydrogen atom mass (scaled, used to compute gas temperature)

   real, parameter :: G_one      = 1.0                      !< Newtonian constant of gravitation in scaled units

   real, parameter :: fpiG       = fpi*G_one                !< four Pi times Newtonian constant of gravitation
   real, parameter :: newtong    = G_one                    !< Newtonian constant of gravitation
   real, parameter :: pc         = 3.086e18                 !< parsec expressed in units of centimetres
   real, parameter :: kpc        = 1000.0                   !< kiloparsec expressed in units of parsec
   real, parameter :: r_gc_sun   = 8500                     !< Sun distance from the Galaxy Center expressed in parsecs

#endif /* SCALED */

end module constants
