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

!>
!! \brief Calculates radiative energy loss.
!! \details PURPOSE:  This routine finds new energy in the source step by solving the implicit eqn:
!!      de / dt = - d**2 * COOL(e) + HEAT
!! where COOL is an empirical cooling function of e only, and HEAT is an empirical heating function.
!<

module thermal
! pulled by THERM

   use constants, only: cbuff_len

   implicit none

   private
   public ::  maxdeint, init_thermal, thermal_active, cfl_coolheat, src_thermal_exec

   character(len=cbuff_len) :: cool_model, heat_model
   logical                  :: thermal_active
   real                     :: alpha_cool, L0_cool, G0_heat, G1_heat, G2_heat, cfl_coolheat
   real                     :: Teql         !> temperature of cooling / heating equilibrium

contains

!--------------------------------------------------------------------------

   subroutine init_thermal

      use constants,      only: PIERNIK_INIT_MPI
      use dataio_pub,     only: code_progress, die, nh, printinfo
      use mpisetup,       only: cbuff, lbuff, rbuff, master, slave, piernik_MPI_Bcast
      use units,          only: cm, erg, sek, mH

      implicit none

      real :: G0, G1, G2 !> standard heating model coefficients in cgs units
      real :: Lambda0    !> power law cooling model coefficient in cgs units
      real :: x_ion      !> ionization degree

      namelist /THERMAL/ thermal_active, cool_model, heat_model, Lambda0, alpha_cool, Teql, G0, G1, G2, x_ion, cfl_coolheat

      if (code_progress < PIERNIK_INIT_MPI) call die("[thermal:init_thermal] mpi not initialized.")

#ifdef VERBOSE
      if (master) call printinfo("[thermal:init_thermal] Commencing thermal module initialization")
#endif /* VERBOSE */

      thermal_active = .True.
      cool_model     = 'power_law'
      heat_model     = 'G012'
      alpha_cool     = 1.0
      Teql           = 1000.0
      Lambda0        = 1.0e-25
      G0             = 1.0e-25
      G1             = 1.0e-25
      G2             = 1.0e-27
      x_ion          = 1.0
      cfl_coolheat   = 0.1

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=THERMAL)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=''
         read(unit=nh%lun, nml=THERMAL, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "THERMAL")
         read(nh%cmdl_nml,nml=THERMAL, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "THERMAL", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=THERMAL)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1) = Lambda0
         rbuff(2) = alpha_cool
         rbuff(3) = Teql
         rbuff(4) = G0
         rbuff(5) = G1
         rbuff(6) = G2
         rbuff(7) = x_ion
         rbuff(8) = cfl_coolheat

         lbuff(1) = thermal_active

         cbuff(1) = cool_model
         cbuff(2) = heat_model

      endif

      call piernik_MPI_Bcast(cbuff, cbuff_len)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         cool_model     = cbuff(1)
         heat_model     = cbuff(2)

         thermal_active = lbuff(1)

         Lambda0        = rbuff(1)
         alpha_cool     = rbuff(2)
         Teql           = rbuff(3)
         G0             = rbuff(4)
         G1             = rbuff(5)
         G2             = rbuff(6)
         x_ion          = rbuff(7)
         cfl_coolheat   = rbuff(8)

      endif

      G0_heat = G0      * erg / sek * cm**3 / mH**2 * x_ion**2
      G1_heat = G1      * erg / sek         / mH    * x_ion
      G2_heat = G2      * erg / sek / cm**3
      L0_cool = Lambda0 * erg / sek * cm**3 / mH**2 * x_ion**2

   end subroutine init_thermal

!>
!! \brief Computation of cooling and heating source terms
!<
   subroutine src_thermal_exec(uu, nn, bb, usrc)

      use constants,  only: xdim, ydim, zdim
      use domain,     only: dom
      use fluidindex, only: flind, nmag
      use fluidtypes, only: component_fluid
      use func,       only: emag, ekin

      implicit none

      integer(kind=4),                intent(in)  :: nn                 !< array size
      real, dimension(nn, flind%all), intent(in)  :: uu                 !< vector of conservative variables
      real, dimension(nn, nmag),      intent(in)  :: bb                 !< local copy of magnetic field
      real, dimension(nn, flind%all), intent(out) :: usrc               !< u array update component for sources
!locals

      real, dimension(nn)                         :: eint_src, kin_ener, int_ener, mag_ener
      class(component_fluid), pointer             :: pfl
      integer                                     :: ifl

      usrc = 0.0
      if (.not.thermal_active) return

      do ifl = 1, flind%fluids
         pfl => flind%all_fluids(ifl)%fl
         if (pfl%has_energy) then
            kin_ener = ekin(uu(:, pfl%imx), uu(:, pfl%imy), uu(:, pfl%imz), uu(:, pfl%idn))
            if (pfl%is_magnetized) then
               mag_ener = emag(bb(:, xdim), bb(:, ydim), bb(:, zdim))
               int_ener = uu(:, pfl%ien) - kin_ener - mag_ener
            else
               int_ener = uu(:, pfl%ien) - kin_ener
            endif
            call cool_heat(pfl%gam, nn, uu(:,pfl%idn), int_ener, eint_src)
            usrc(:, pfl%ien) = usrc(:, pfl%ien) + 1./dom%eff_dim * eint_src
         endif
      enddo

   end subroutine src_thermal_exec


   subroutine cool_heat(gamma, n, dens, eint, esrc)

      use units,      only: kboltz, mH

      implicit none

      integer(kind=4),    intent(in)  :: n
      real,               intent(in)  :: gamma
      real, dimension(n), intent(in)  :: dens, eint
      real, dimension(n), intent(out) :: esrc
      real, dimension(n)              :: cfunc, hfunc, temp

      temp(:)  = 0.0
      cfunc(:) = 0.0
      hfunc(:) = 0.0
      esrc(:)  = 0.0

      temp = (gamma-1) * mH / kboltz * eint / dens

      call cool(n, temp, cfunc)

      call heat(n, dens, hfunc)
      esrc =  dens**2*cfunc + hfunc
!      esrc = MIN(esrc, esrc_upper_lim * eint)
!      esrc = MAX(esrc, esrc_lower_lim * eint)

   end subroutine cool_heat

   subroutine maxdeint(cg, max_deint)

      use constants,        only: xdim, ydim, zdim
      use fluidindex,       only: flind
      use fluidtypes,       only: component_fluid
      use func,             only: emag, ekin
      use grid_cont,        only: grid_container

      implicit none

      type(grid_container), pointer, intent(in)              :: cg
      real,                          intent(out)             :: max_deint
      integer                                                :: i, j, ifl
      class(component_fluid), pointer                        :: pfl
      real, dimension(cg%ks:cg%ke)                           :: int_ener, kin_ener, mag_ener, eint_src
      real, dimension(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) :: deint

      do ifl = 1, flind%fluids
         pfl => flind%all_fluids(ifl)%fl
         do j = cg%js, cg%je
            do i = cg%is, cg%ie
               if (pfl%has_energy) then
                  kin_ener = ekin(cg%u(pfl%imx,i,j,cg%ks:cg%ke), cg%u(pfl%imy,i,j,cg%ks:cg%ke), cg%u(pfl%imz,i,j,cg%ks:cg%ke), cg%u(pfl%idn,i,j,cg%ks:cg%ke))
                  if (pfl%is_magnetized) then
                     mag_ener = emag(cg%b(xdim,i,j,cg%ks:cg%ke), cg%b(ydim,i,j,cg%ks:cg%ke), cg%b(zdim,i,j,cg%ks:cg%ke))
                     int_ener = cg%u(pfl%ien,i,j,cg%ks:cg%ke) - kin_ener - mag_ener
                  else
                     int_ener = cg%u(pfl%ien,i,j,cg%ks:cg%ke) - kin_ener
                  endif
               endif
               call cool_heat(pfl%gam, cg%nzb, cg%u(pfl%idn,i,j,cg%ks:cg%ke), int_ener, eint_src)
               deint(i,j,:) = eint_src/int_ener
            enddo
         enddo
      enddo
      max_deint = maxval(abs(deint))

   end subroutine maxdeint

   subroutine cool(n, temp, coolf)

      use dataio_pub, only: msg, warn
      use mpisetup,   only: master

      implicit none

      integer(kind=4),    intent(in)  :: n
      real, dimension(n), intent(in)  :: temp
      real, dimension(n), intent(out) :: coolf


      select case (cool_model)
         case ('power_law')
            coolf = -L0_cool * (temp/Teql)**(alpha_cool)
         case ('null')
            return
        case default
          write(msg,'(3a)') 'Cool model: ',cool_model,' not implemented'
          if (master) call warn(msg)
      end select

   end subroutine cool

   subroutine heat(n, dens, heatf)

      use dataio_pub, only: msg, warn
      use mpisetup,   only: master

      implicit none

      integer(kind=4),    intent(in)  :: n
      real, dimension(n), intent(in)  :: dens
      real, dimension(n), intent(out) :: heatf

      select case (heat_model)
        case ('G012')
          heatf =  G0_heat * dens**2 + G1_heat * dens + G2_heat
        case ('null')
          return
        case default
           write(msg,'(3a)') 'Heat model: ',heat_model,' not implemented'
           if (master) call warn(msg)
      end select

   end subroutine heat

!--------------------------------------------------------------------------

end module thermal
