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

   use constants, only: cbuff_len, INVALID

   implicit none

   private
   public ::  maxdeint, init_thermal, thermal_active, cfl_coolheat, src_thermal_exec, EIS

   character(len=cbuff_len) :: cool_model, heat_model
   logical                  :: thermal_active
   real                     :: alpha_cool, L0_cool, G0_heat, G1_heat, G2_heat, cfl_coolheat
   real                     :: Teql         !> temperature of cooling / heating equilibrium
   integer(kind=4), protected    :: itemp = INVALID, ieis = INVALID

contains

!--------------------------------------------------------------------------

   subroutine init_thermal

      use cg_list_global,   only: all_cg
      use constants,        only: PIERNIK_INIT_MPI
      use dataio_pub,       only: code_progress, die, nh, printinfo
      use mpisetup,         only: cbuff, lbuff, rbuff, master, slave, piernik_MPI_Bcast
      use named_array_list, only: qna
      use units,            only: cm, erg, sek, mH

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

      call all_cg%reg_var('Temperature')          ! Make it cleaner
      itemp = qna%ind('Temperature')
      call all_cg%reg_var('Temp_EIS')          ! Make it cleaner
      ieis = qna%ind('Temp_EIS')

   end subroutine init_thermal

!>
!! \brief Computation of cooling and heating source terms
!<
   subroutine src_thermal_exec(uu, nn, cg, i1, i2, sweep, bb, usrc)

      use constants,  only: xdim, ydim, zdim
      use domain,     only: dom
      use fluidindex, only: flind, nmag
      use fluidtypes, only: component_fluid
      use func,       only: emag, ekin
      use grid_cont,  only: grid_container

      implicit none

      type(grid_container), pointer,  intent(in)  :: cg                 !< current grid piece
      integer(kind=4),                intent(in)  :: nn, i1, i2                 !< array size
      integer(kind=4),    intent(in)  :: sweep              !< direction (x, y or z) we are doing calculations for  
      real, dimension(nn, flind%all), intent(in)  :: uu                 !< vector of conservative variables
      real, dimension(nn, nmag),      intent(in)  :: bb                 !< local copy of magnetic field
      real, dimension(nn, flind%all), intent(out) :: usrc               !< u array update component for sources
!locals

      real, dimension(nn)                         :: eint_src, kin_ener, int_ener, mag_ener
      class(component_fluid), pointer             :: pfl
      integer                                     :: ifl
      logical                                     :: store

      usrc = 0.0
      store=.true.
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
            call cool_heat(pfl%gam, nn, uu(:,pfl%idn), cg, i1, i2, sweep, int_ener, eint_src, store)
            usrc(:, pfl%ien) = usrc(:, pfl%ien) + 1./dom%eff_dim * eint_src
         endif
      enddo

   end subroutine src_thermal_exec


   subroutine cool_heat(gamma, n, dens, cg, i1, i2, sweep, eint, esrc, store)

     
      use domain,     only: dom
      use global,     only: nstep, dt
      use grid_cont,  only: grid_container
      use units,      only: kboltz, mH

      implicit none

      type(grid_container), pointer,  intent(in)  :: cg                 !< current grid piece
      integer(kind=4),    intent(in)  :: n, i1, i2
      integer(kind=4),    intent(in)  :: sweep              !< direction (x, y or z) we are doing calculations for  
      
      real,               intent(in)  :: gamma
      real, dimension(n), intent(in)  :: dens, eint
      logical,            intent(in)  :: store
      real, dimension(n), intent(out) :: esrc
      real, dimension(n)              :: cfunc, hfunc, temp, eint2
      real, dimension(:),  pointer    :: ta, t_eis

      temp(:)  = 0.0
      cfunc(:) = 0.0
      hfunc(:) = 0.0
      esrc(:)  = 0.0
      ta => cg%q(itemp)%get_sweep(sweep, i1, i2)
      t_eis => cg%q(ieis)%get_sweep(sweep, i1, i2)
      if (nstep .eq. 0) then
         temp = (gamma-1) * mH / kboltz * eint / dens
      else
         temp=t_eis
      endif
      !temp = (gamma-1) * mH / kboltz * eint / dens
      call cool(n, temp, cfunc)
      !print *, 'cool', dens**2*cfunc
      call heat(n, dens, hfunc)
      esrc =  dens**2*cfunc + hfunc
      if (store) then
         eint2 = eint + 1./dom%eff_dim * esrc *dt
         ta = (gamma-1) * mH / kboltz * eint2 / dens
         print *, 'esrc', esrc
      endif

!      esrc = MIN(esrc, esrc_upper_lim * eint)
!      esrc = MAX(esrc, esrc_lower_lim * eint)

    end subroutine cool_heat
   

   subroutine maxdeint(cg, max_deint, min_tcool)

      use constants,        only: xdim, ydim, zdim
      use fluidindex,       only: flind
      use fluidtypes,       only: component_fluid
      use func,             only: emag, ekin
      use grid_cont,        only: grid_container
      use units,            only: kboltz, mH

      implicit none

      type(grid_container), pointer, intent(in)              :: cg
      real,                          intent(out)             :: max_deint, min_tcool
      integer                                                :: i, j, ifl
      class(component_fluid), pointer                        :: pfl
      real, dimension(cg%ks:cg%ke)                           :: int_ener, kin_ener, mag_ener, eint_src, temp, tcool
      real, dimension(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) :: deint
      logical                                                :: store

      store=.false.
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
               call cool_heat(pfl%gam, cg%nzb, cg%u(pfl%idn,i,j,cg%ks:cg%ke), cg, i, j, 3, int_ener, eint_src, store)
               temp = (pfl%gam-1) * mH / kboltz * int_ener / cg%u(pfl%idn,i,j,cg%ks:cg%ke)
               tcool = - kboltz * temp / ((pfl%gam-1) * mH *  L0_cool * (temp/Teql)**(alpha_cool))
               deint(i,j,:) = eint_src/int_ener
               !if (maxval(abs(deint(i,j,:))) .gt. 1000000) then
               !   print *, i,j, 'OOOPS', maxval(eint_src), minval(int_ener)
               !endif
            enddo
         enddo
      enddo
      max_deint = maxval(abs(deint))
      min_tcool = minval(abs(tcool))

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
            !print *, 'L0_cool', L0_cool, alpha_cool
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
           !print *, 'heat', G0_heat, G1_heat, G2_heat, heatf
        case ('null')
          return
        case default
           write(msg,'(3a)') 'Heat model: ',heat_model,' not implemented'
           if (master) call warn(msg)
      end select

    end subroutine heat


     subroutine EIS(dt)

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use fluidindex, only: flind
      use grid_cont,  only: grid_container
      use named_array_list, only: wna
!      use units,      only: kboltz, mH

      implicit none

      type(cg_list_element), pointer  :: cgl
      type(grid_container),  pointer  :: cg
      
!      integer(kind=4),    intent(in)  :: sweep              !< direction (x, y or z) we are doing calculations for                     
!      real, dimension(n), intent(out) :: esrc
      !      real, dimension(n)              :: cfunc, hfunc, temp, T, Tn, eint2
      real, intent(in)                :: dt
      real, dimension(:, :, :),  pointer    :: t_eis, ta, T
      real                            :: gamma
      real, dimension(:,:,:), pointer :: dens



      gamma = flind%all_fluids(1)%fl%gam
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         dens =>  cg%w(wna%fi)%span(flind%all_fluids(1)%fl%idn,cg%ijkse)!lh_out)
         ta => cg%q(itemp)%span(cg%ijkse)!lh_out)
         !print *, 'TEMP1', ta
         t_eis => cg%q(ieis)%span(cg%ijkse)!lh_out)
         T => cg%q(ieis)%span(cg%ijkse)!lh_out)
         print *, 'temp0', ta(1,1,1), t_eis(1,1,1)
         !print *, 'DENS', dens
         call temp_EIS(dt, shape(ta), gamma, ta, dens, T)
         !print *, 'boltzman', kboltz, mH
         !print *, 'YOP'
         !eint2 = eint + 1./dom%eff_dim * esrc *dt
         !Tn = (gamma-1) * mH / kboltz * eint2 / dens
         t_eis = T
      
         print *, 'TEMP_EIS', T(1,1,1)
         cgl => cgl%nxt
      enddo


    end subroutine EIS

        
   subroutine temp_EIS(dt, n, gamma, temp, dens, Tnew)

      use dataio_pub, only: msg, warn
      use func,       only: operator(.equals.)
      use mpisetup,   only: master
      use units,      only: kboltz, mH

      implicit none

      integer(kind=4), dimension(3),     intent(in)  :: n
      real,               intent(in)  :: gamma, dt
      real, dimension(n(1), n(2), n(3)), intent(in)  :: temp, dens
      real, dimension(n(1), n(2), n(3))              :: TEF, Teql1, tcool
      real, dimension(n(1), n(2), n(3)), intent(out) :: Tnew

      Teql1(:,:,:)=Teql
      select case (cool_model)
      case ('power_law')
            if (alpha_cool .equals. 1.0) then
               TEF = log(Teql1/temp)
               !print *, 'TEF',  TEF + (gamma-1) * dens / kboltz * dt * coolf_pl(Teql1) * mH
               Tnew = invTEF_pl_a1( TEF + (gamma-1) * dens / kboltz * dt * coolf_pl(Teql1) * mH / Teql1)
               !print *, 'Tnew', Tnew
            else
               TEF= (1/(1-alpha_cool)) * (1 - (Teql/temp)**(alpha_cool-1) )
               tcool = kboltz * temp / ((gamma-1) * mH * coolf_pl(Teql1))
               print *, 'invTEF', tcool(1,1,1), dt!, gamma, (1-(1-alpha_cool)*dt/tcool)**(1/(1-alpha_cool))!TEF + (gamma-1) * dens / kboltz * dt * coolf_pl(Teql1) * mH / Teql1, dt, invTEF_pl(0.0*TEF), invTEF_pl(-1.0+TEF)
               !Tnew = invTEF_pl( TEF + (gamma-1) * dens / kboltz * dt * coolf_pl(Teql1) * mH / Teql1)
               Tnew = temp * (1 - (1-alpha_cool)*dt / tcool) **(1/(1-alpha_cool))
            endif
            !coolf = coolf_pl(Tnew)
            !coolf = -L0_cool * (temp/Teql)**(alpha_cool)
         case ('null')
            return
         case default
            write(msg,'(3a)') 'Cool model: ',cool_model,' not implemented'
            if (master) call warn(msg)  
       end select
         

     contains
       
            function coolf_pl(T) result(coolf)
              real, dimension(n(1), n(2), n(3)), intent(in) :: T
              real, dimension(n(1), n(2), n(3))             :: coolf
              coolf = L0_cool * (T/Teql1)**(alpha_cool)
            end function coolf_pl
            
            function invTEF_pl_a1(Y) result(T)
              real, dimension(n(1), n(2), n(3)), intent(in) :: Y
              real, dimension(n(1), n(2), n(3))             :: T
              T = Teql * exp(-Y)
            end function invTEF_pl_a1
            
            function invTEF_pl(Y) result(T)
              real, dimension(n(1), n(2), n(3)), intent(in) :: Y
              real, dimension(n(1), n(2), n(3))             :: T
              T = Teql * (1 - (1-alpha_cool)*Y) ** (1/(1-alpha_cool))
            end function invTEF_pl

   end subroutine temp_EIS

!--------------------------------------------------------------------------

end module thermal
