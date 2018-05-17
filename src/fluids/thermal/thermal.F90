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
   public ::  maxdeint, cool_heat, init_thermal, cool_model, heat_model, A_heat, A_cool, alpha_cool, thermal_active, cfl_coolheat, beta_cool

   character(len=cbuff_len) :: cool_model, heat_model
   logical                  :: thermal_active
   real                     :: A_cool, A_heat, alpha_cool, beta_cool, cfl_coolheat

contains

!--------------------------------------------------------------------------

   subroutine init_thermal

      use constants,      only: PIERNIK_INIT_MPI
      use dataio_pub,     only: code_progress, die, nh, printinfo
      use mpisetup,       only: cbuff, lbuff, rbuff, master, slave, piernik_MPI_Bcast

      implicit none

      namelist /THERMAL/ thermal_active, cool_model, heat_model, A_cool, A_heat, alpha_cool, beta_cool, cfl_coolheat

      if (code_progress < PIERNIK_INIT_MPI) call die("[thermal:init_thermal] mpi not initialized.")

#ifdef VERBOSE
      if (master) call printinfo("[thermal:init_thermal] Commencing thermal module initialization")
#endif /* VERBOSE */

      thermal_active = .True.
      cool_model     = 'power_law'
      heat_model     = 'beta_coef'
      A_cool         = 1.0
      A_heat         = 1.0
      alpha_cool     = 1.0
      beta_cool      = 1.0
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

         rbuff(1) = A_cool
         rbuff(2) = alpha_cool
         rbuff(3) = A_heat
         rbuff(4) = beta_cool
         rbuff(5) = cfl_coolheat

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

         A_cool         = rbuff(1)
         alpha_cool     = rbuff(2)
         A_heat         = rbuff(3)
         beta_cool      = rbuff(4)
         cfl_coolheat   = rbuff(5)

      endif

   end subroutine init_thermal

   subroutine cool_heat(gamma, n, dens, eint, esrc)

      use units,      only: kboltz, mH

      implicit none

      integer,            intent(in)  :: n
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
         do j=cg%js,cg%je
            do i=cg%is,cg%ie
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
      max_deint = MAXVAL(ABS(deint))

   end subroutine maxdeint

   subroutine cool(n, temp, coolf)

      use dataio_pub, only: msg, warn
      use mpisetup,   only: master

      implicit none

      integer,            intent(in)  :: n
      real, dimension(n), intent(in)  :: temp
      real, dimension(n), intent(out) :: coolf


      select case (cool_model)
         case ('power_law')
            coolf = -A_cool * temp**(alpha_cool)
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

      integer,            intent(in)  :: n
      real, dimension(n), intent(in)  :: dens
      real, dimension(n), intent(out) :: heatf

      select case (heat_model)
        case ('beta_coef')
          heatf =  A_heat * dens**beta_cool  !> \todo trzeba dodac czytanie parametrow grzania z problem.par
        case ('null')
          return
        case default
           write(msg,'(3a)') 'Heat model: ',heat_model,' not implemented'
           if (master) call warn(msg)
      end select

   end subroutine heat

!--------------------------------------------------------------------------

end module thermal
