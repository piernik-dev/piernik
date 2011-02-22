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
!    Initial implementation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.h"
#include "macros.h"
!>
!! \brief (DW) [R] Module containing physical %constants
!! \details Module units contains physical %constants for several units systems.
!! To use one system a proper value of units_set should be set.
!! Available units systems defined by units_set value:
!! @n
!! @n @b PSM (Parsec - Solar_mass - Megayear) - good for global galactic simulations
!! @n length --> pc,     mass --> Msun,        time --> myr,        miu0 --> 4*pi,    temperature --> kelvin
!! @n
!! @n @b PLN (PLaNetary) - good for planetary nebulae
!! @n length --> AU,     mass --> Mjup,        time --> yr,         miu0 --> 4*pi,    temperature --> kelvin
!! @n
!! @n @b KSG (Kiloparsec - Solar_mass - Gigayear) - good for galactic and intergalactic simulations
!! @n length --> kpc,    mass --> 10^6*Msun,   time --> Gyr,        miu0 --> 4*pi,    temperature --> kelvin
!! @n
!! @n @b KSM (Kiloparsec - Solar_mass - Megayear)
!! @n length --> kpc,    mass --> Msun,        time --> myr,        miu0 --> 4*pi,    temperature --> kelvin
!! @n
!! @n @b PGM (Parsec - Gravity=1 - Megayear) - modification of PSM system, where %gravity constant is one
!! @n length --> pc,     newtong --> 1.0,      time --> myr,        miu0 --> 4*pi,    temperature --> kelvin
!! @n
!! @n @b SSY (10^Sixteenth cm - Solar_mass - Year) - good for circumstellar simulations (planetaries etc.)
!! @n length --> 10^16 cm,  mass --> Msun,     time --> year,       miu0 --> 4*pi,    temperature --> kelvin
!! @n
!! @n @b SI - (franc. Système International d'Unités) '
!! @n length --> metr,   mass --> kg,          time --> sek,        miu0 --> 4*pi,    temperature --> kelvin
!! @n
!! @n @b CGS - (Centimetre - Gram - Second)
!! @n length --> cm,     mass --> gram,        time --> sek,        miu0 --> 4*pi,    temperature --> kelvin
!! @n
!! @n @b WT4 - unit system for Wengen Test #4
!! @n length --> 6.25AU, mass --> 0.1 M_sun,   time --> 2.5**3.5  pi years (=> G ≈ 1.)
!! @n
!! @n @b USER - units system defined by user
!! @n following variables from UNITS namelist should be specified: miu0, kelvin, cm, gram, sek
!! @n
!! @n @b SCALED - a suit of %units without physical units. This is automatically set while neither of the former systems is chosen.
!<
module units

   use dataio_pub, only: cbuff_len

   implicit none

   public                                                ! QA_WARN no secrets are kept here
   private :: au_cm, pc_au, pc_cm, msun_g, mjup_g, day_s, yr_day, yr_s, newton_cgs, kB_cgs  ! QA_WARN don't use those vars outside units!

   character(len=cbuff_len) :: units_set                !< type of units set
   character(len=cbuff_len) :: s_len_u                   !< name of length unit
   character(len=cbuff_len) :: s_time_u                  !< name of time unit
   character(len=cbuff_len) :: s_mass_u                  !< name of mass unit

   real(kind=8), parameter :: au_cm       =  1.49597870691d13   ! Astonomical unit [cm]
   real(kind=8), parameter :: pc_au       =  206264.806248712d0 ! Parsec [AU] 1 pc/1 AU = 1./atan(pi/180. * 1./3600.)
   real(kind=8), parameter :: pc_cm       =  pc_au*au_cm        ! Parsec [cm]
   real(kind=8), parameter :: Msun_g      =  1.98892e33         ! Solar mass [g]
   real(kind=8), parameter :: Mjup_g      =  1.8986e30          ! Jovian mass [g]
   real(kind=8), parameter :: day_s       =  24.d0*3600.d0      ! Earth's day [s]
   real(kind=8), parameter :: yr_day      =  365.256363051d0    ! sideral year [days]
   real(kind=8), parameter :: yr_s        =  yr_day*day_s       ! Year [s]
   real(kind=8), parameter :: newton_cgs  =  6.67428e-8         ! Gravitational constant [ cm^3 / g s^2 ]
   real(kind=8), parameter :: kB_cgs      =  1.3806504e-16      ! Boltzmann constant [ g cm^2 / K s^2 ]

   real, protected :: cm      !< centimetre, length unit
   real, protected :: gram    !< gram, mass unit
   real, protected :: sek     !< second, time unit
   real, protected :: miu0    !< permeability
   real, protected :: kelvin  !< kelvin, temperature unit

! length units:
   real, protected :: metr                                  !< metre, length unit
   real, protected :: km                                    !< kilometer, length unit
   real, protected :: au                                    !< astronomical unit (length unit)
   real, protected :: pc                                    !< parsec, length unit
   real, protected :: kpc                                   !< kiloparsec, length unit
   real, protected :: lyr                                   !< light year, length unit
! time units:
   real, protected :: minute                                !< minute, time unit
   real, protected :: hour                                  !< hour, time unit
   real, protected :: day                                   !< day, time unit
   real, protected :: year                                  !< year, time unit
   real, protected :: myr                                   !< megayear, time unit
! mass units:
   real, protected :: kg                                    !< kilogram, mass unit
   real, protected :: me                                    !< electron mass
   real, protected :: mp                                    !< proton mass
   real, protected :: mH                                    !< hydrogen atom mass
   real, protected :: amu                                   !< atomic mass unit
   real, protected :: Msun                                  !< Solar mass, mass unit
   real, protected :: gmu                                   !< galactic mass unit
! force units:
   real, protected :: newton                                !< 1N (SI force unit)
   real, protected :: dyna                                  !< 1 dyna (cgs force unit)
! energy units:
   real, protected :: joul                                  !< 1J (SI energy unit)
   real, protected :: erg                                   !< 1 erg (cgs energy unit)
   real, protected :: eV                                    !< 1 eV
! density units:
   real, protected :: ppcm3                                 !< spatial density unit
   real, protected :: ppcm2                                 !< column density unit
! physical constants:
   real, protected :: kboltz                                !< boltzmann constant
   real, protected :: gasRconst                             !< gas constant R =  8.314472e7*erg/kelvin/mol
   real, protected :: N_A                                   !< Avogadro constant
   real, protected :: clight                                !< speed of light in vacuum
   real, protected :: Gs                                    !< 1 Gs (cgs magnetic induction unit)
   real, protected :: mGs                                   !< 1 microgauss
   real, protected :: Tesla                                 !< 1 T (SI magnetic induction unit)
   real, protected :: newtong                               !< newtonian constant of gravitation
   real, protected :: fpiG                                  !< four Pi times Newtonian constant of gravitation (commonly used in self-gravity routines)
   real, protected :: planck                                !< Planck constant
   real, protected :: r_gc_sun                              !< Sun distance from the Galaxy Center
   real, protected :: vsun                                  !< velocity value of Sun in the Galaxy
   real, protected :: sunradius                             !< radius of Sun
   real, protected :: Lsun                                  !< luminosity of Sun
   real, protected :: Mearth                                !< mass of Earth
   real, protected :: earthradius                           !< radius of Earth

contains
!>
!! \brief Routine initializing units module
!!
!! \details
!! @b UNITS
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>units_set </td><td>'scaled'</td><td>string of characters</td><td>type of units set     </td></tr>
!! <tr><td>miu0      </td><td>4*pi    </td><td>real                </td><td>\copydoc units::miu0  </td></tr>
!! <tr><td>kelvin    </td><td>1       </td><td>real                </td><td>\copydoc units::kelvin</td></tr>
!! <tr><td>cm        </td><td>1       </td><td>real                </td><td>\copydoc units::cm    </td></tr>
!! <tr><td>gram      </td><td>1       </td><td>real                </td><td>\copydoc units::gram  </td></tr>
!! <tr><td>sek       </td><td>1       </td><td>real                </td><td>\copydoc units::sek   </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_units

      use constants,  only: one, pi, fpi, small
      use mpisetup,   only: cbuff, rbuff, buffer_dim, comm, ierr, master, slave
      use mpi,        only: MPI_CHARACTER, MPI_DOUBLE_PRECISION
      use dataio_pub, only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml  ! QA_WARN required for diff_nml
      use dataio_pub, only: warn, printinfo, msg, die, code_progress, PIERNIK_INIT_MPI

      implicit none

      logical, save            :: scale_me = .false.
      logical                  :: to_stdout

      namelist /UNITS/ units_set, miu0, kelvin, cm, gram, sek

      if (code_progress < PIERNIK_INIT_MPI) call die("[units:init_units] MPI not initialized.")

      units_set='scaled'

      miu0   = fpi
      kelvin = one
      cm     = small
      gram   = small
      sek    = small

      if (master) then

         diff_nml(UNITS)

         cbuff(1) = units_set

         rbuff(1) = miu0
         rbuff(2) = kelvin
         rbuff(3) = cm
         rbuff(4) = gram
         rbuff(5) = sek

      endif

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_Bcast(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (slave) then

         units_set = cbuff(1)

         miu0   = rbuff(1)
         kelvin = rbuff(2)
         cm     = rbuff(3)
         gram   = rbuff(4)
         sek    = rbuff(5)

      endif

      to_stdout = .false.
#ifdef VERBOSE
      to_stdout = .true.
#endif  /* VERBOSE */
      s_len_u  = ' undefined'; s_time_u = s_len_u; s_mass_u = s_len_u

!>
!! \deprecated BEWARE: miu0 and kelvin may be overwritten by values from problem.par even though we choose units_set value one of the following
!! nevertheless, they are not used so far (r3612)
!<
      select case (trim(units_set))
         case ("PSM", "psm")
            ! PSM  uses: length --> pc,     mass --> Msun,        time --> myr,        miu0 --> 4*pi,    temperature --> kelvin
            cm         = 1.0/pc_cm            !< centimetre, length unit
            sek        = 1.0/(1.0e6*yr_s)     !< second, time unit
            gram       = 1.0/msun_g           !< gram, mass unit
            s_len_u  = ' [pc]'
            s_time_u = ' [Myr]'
            s_mass_u = ' [M_sun]'

         case ("PLN", "pln")
            ! PLN  uses: length --> AU,     mass --> Mjup,        time --> yr,         miu0 --> 4*pi,    temperature --> kelvin
            cm         = 1.0/au_cm            !< centimetre, length unit
            sek        = 1.0/yr_s             !< second, time unit
            gram       = 1.0/mjup_g           !< gram, mass unit
            s_len_u  = ' [AU]'
            s_time_u = ' [yr]'
            s_mass_u = ' [M_jup]'

         case ("KSG", "ksg")
            ! KSG  uses: length --> kpc,    mass --> 10^6*Msun,   time --> Gyr,        miu0 --> 4*pi,    temperature --> kelvin
            cm         = 1.0/(1.0e3*pc_cm)    !< centimetre, length unit
            sek        = 1.0/(1.0e9*yr_s)     !< second, time unit
            gram       = 1.0/(1.0e6*msun_g)   !< gram, mass unit
            s_len_u  = ' [kpc]'
            s_time_u = ' [Gyr]'
            s_mass_u = ' [10^6 M_sun]'

         case ("KSM", "ksm")
            ! KSM  uses: length --> kpc,    mass --> Msun,        time --> myr,        miu0 --> 4*pi,    temperature --> kelvin
            cm =         1.0/(1.0e3*pc_cm)    !< centimetre, length unit
            sek =        1.0/(1.0e6*yr_s)     !< second, time unit
            gram =       1.0/msun_g           !< gram, mass unit
            s_len_u  = ' [kpc]'
            s_time_u = ' [Myr]'
            s_mass_u = ' [M_sun]'

         case ("PGM", "pgm")
            ! PGM  uses: length --> pc,     newtong --> 1.0,      time --> myr,        miu0 --> 4*pi,    temperature --> kelvin
            cm         = 1.0/pc_cm            !< centimetre, length unit
            sek        = 1.0/(1.0e6*yr_s)     !< second, time unit
            gram       = newton_cgs*cm**3/1.0/sek**2      !< gram, mass unit  G = 1.0
            s_len_u  = ' [pc]'
            s_time_u = ' [Myr]'
            s_mass_u = ' [-> G=1]'

         case ("SSY", "ssy")
            ! SSY  uses: length --> 10^16 cm,  mass --> Msun,     time --> year,       miu0 --> 4*pi,    temperature --> kelvin
            cm         = 1.0/1.0e16           !< centimetre, length unit
            sek        = 1.0/yr_s             !< second, time unit
            gram       = 1.0/msun_g           !< gram, mass unit
            s_len_u  = ' [10^16 cm]'
            s_time_u = ' [yr]'
            s_mass_u = ' [M_sun]'

         case ("SI", "si")
            ! SI   uses: length --> metr,   mass --> kg,          time --> sek,        miu0 --> 4*pi,    temperature --> kelvin
            cm         = 1.0/1.0e2            !< centimetre, length unit
            sek        = 1.0                  !< second, time unit
            gram       = 1.0/1.0e3            !< gram, mass unit
            s_len_u  = ' [m]'
            s_time_u = ' [s]'
            s_mass_u = ' [kg]'

         case ("CGS", "cgs")
            ! CGS  uses: length --> cm,     mass --> gram,        time --> sek,        miu0 --> 4*pi,    temperature --> kelvin
            cm         = 1.0                  !< centimetre, length unit
            sek        = 1.0                  !< second, time unit
            gram       = 1.0                  !< gram, mass unit
            s_len_u  = ' [cm]'
            s_time_u = ' [s]'
            s_mass_u = ' [g]'

         case ("WT4", "wt4")
            ! WT4  uses: length --> 6.25AU, mass --> 0.1 M_sun,   time --> 2.5**3.5 /pi years (=> G \approx 1. in Wengen Test #4),
            cm          = 1./(6.25*au_cm)     !< centimetre, length unit
            ! It's really weird that use of 2.5**3.5 here can cause Internal Compiler Error at multigridmultipole.F90:827
            sek         = 1./(24.7052942200655/pi * yr_s) !< year, time unit; 24.7052942200655 = 2.5**3.5
            gram        = 1/(0.1*msun_g)      !< gram, mass unit
            s_len_u  = ' [6.25 AU]'
            s_time_u = ' [2.5**3.5 /pi years]'
            s_mass_u = ' [0.1 M_sun]'

         case ("USER", "user")
            if (master) call warn("[units:init_units] PIERNIK will use 'cm', 'sek', 'gram' defined in problem.par")
            if (any([cm == small, sek == small, gram == small])) &
               call die("[units:init_units] units_set=='user', yet one of {'cm','sek','gram'} is not set in problem.par") ! Don't believe in coincidence
            to_stdout = .true.               ! Force output in case someone is not aware what he/she is doing
            s_len_u   = ' [user unit]'
            s_time_u  = s_len_u; s_mass_u  = s_len_u

         case default
            if (master) call warn("[units:init_units] you haven't chosen units set. That means physical vars taken from 'units' are worthless or equal 1")
            cm   = small
            gram = small
            sek  = small

            scale_me = .true.

      end select

      if (master .and. .not. scale_me) then
         write(msg,'(a,es14.7,a)') '[units:init_units] cm   = ', cm,   trim(s_len_u)
         call printinfo(msg, to_stdout)
         write(msg,'(a,es14.7,a)') '[units:init_units] sek  = ', sek,  trim(s_time_u)
         call printinfo(msg, to_stdout)
         write(msg,'(a,es14.7,a)') '[units:init_units] gram = ', gram, trim(s_mass_u)
         call printinfo(msg, to_stdout)
      endif

! length units:
      metr       = 1.0e2*cm                 !< metre, length unit
      km         = 1.0e5*cm                 !< kilometer, length unit
      au         = au_cm*cm                 !< astronomical unit (length unit)
      pc         = pc_cm*cm                 !< parsec, length unit
      kpc        = 1000.0*pc                !< kiloparsec, length unit
      lyr        = 9.4605e17*cm             !< light year, length unit
! time units:
      minute     = 60.0*sek                 !< minute, time unit
      hour       = 3600.0*sek               !< hour, time unit
      day        = day_s*sek                !< day, time unit
      year       = yr_s*sek                 !< year, time unit
      myr        = 1.0e6*year               !< megayear, time unit
! mass units:
      kg         = 1.0e3*gram               !< kilogram, mass unit
      me         = 9.109558e-28*gram        !< electron mass
      mp         = 1.672614e-24*gram        !< proton mass
      mH         = 1.673559e-24*gram        !< hydrogen atom mass
      amu        = 1.660531e-24*gram        !< atomic mass unit
      Msun       = msun_g*gram              !< Solar mass, mass unit
      gmu        = 2.32e7*Msun              !< galactic mass unit
! force units:
      newton     = kg*metr/sek**2           !< 1N (SI force unit)
      dyna       = gram*cm/sek**2           !< 1 dyna (cgs force unit)
! energy units:
      joul       = kg*metr**2/sek**2        !< 1J (SI energy unit)
      erg        = gram*cm**2/sek**2        !< 1 erg (cgs energy unit)
      eV         = 1.6022e-12*erg           !< 1 eV
! density units:
      ppcm3      = 1.36 * mp / cm**3        !< spatial density unit
      ppcm2      = 1.36 * mp / cm**2        !< column density unit
! temperature units:
      kelvin     = 1.0                      !< kelvin, temperature unit
! physical constants:
      kboltz     = kB_cgs*erg/kelvin        !< boltzmann constant
      gasRconst  = 8.314472e7*erg/kelvin    !< gas constant R =  8.314472e7*erg/kelvin/mol = k_B * N_A
      N_A        = gasRconst / kboltz       !< Avogadro constant
      clight     = 2.997924562e10*cm/sek    !< speed of light in vacuum
      Gs         = sqrt(miu0*gram/cm)/sek   !< 1 Gs (cgs magnetic induction unit)
      mGs        = Gs*1.e-6                 !< 1 microgauss
      Tesla      = 1.e4*Gs                  !< 1 T (SI magnetic induction unit)
      newtong    = newton_cgs*cm**3/gram/sek**2 !< newtonian constant of gravitation
      fpiG       = fpi*newtong              !< four Pi times Newtonian constant of gravitation (commonly used in self-gravity routines)
      planck     = 6.626196e-27*erg*sek     !< Planck constant
      r_gc_sun   = 8.5*kpc                  !< Sun distance from the Galaxy Center
      vsun       = 220.0*km/sek             !< velocity value of Sun in the Galaxy
      sunradius  = 6.9598e10*cm             !< radius of Sun
      Lsun       = 3.826e33*erg/sek         !< luminosity of Sun
      Mearth     = 5.977e27*gram            !< mass of Earth
      earthradius= 6378.17*km               !< radius of Earth

      ! Following physicical constants are used in various modules.
      ! They need to have some sane values.
      if (scale_me) then
         kboltz    = 1.0  ! dataio
         gasRconst = 1.0  ! dataio
         mH        = 1.0  ! dataio
         fpiG      = fpi  ! multigrid_gravity, poissonsolver
         newtong   = 1.0  ! multigridmultipole, gravity, poissonsolver
      endif

   end subroutine init_units

end module units
