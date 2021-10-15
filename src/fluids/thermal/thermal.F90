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
   public ::  init_thermal, thermal_active, cfl_coolheat, thermal_sources, itemp, fit_cooling_curve

   character(len=cbuff_len)        :: cool_model, cool_curve, heat_model, scheme, cool_file
   logical                         :: thermal_active
   real                            :: alpha_cool, L0_cool, G0_heat, G1_heat, G2_heat, cfl_coolheat
   real                            :: Teq          !> cooling parameter
   real                            :: Teql         !> temperature of cooling / heating equilibrium
   integer(kind=4), protected      :: itemp = INVALID
   real                            :: x_ion        !> ionization degree
   integer                         :: isochoric    !> 1 for isochoric, 2 for isobaric
   real                            :: d_isochoric  ! constant density used in isochoric case
   real, dimension(:), allocatable :: Tref, alpha, lambda0
   integer                         :: nfuncs

contains

   subroutine init_thermal

      use cg_list_global,   only: all_cg
      use constants,        only: PIERNIK_INIT_MPI
      use dataio_pub,       only: code_progress, die, nh, printinfo
      use mpisetup,         only: cbuff, lbuff, rbuff,ibuff,  master, slave, piernik_MPI_Bcast
      use named_array_list, only: qna
      use units,            only: cm, erg, sek, mH

      implicit none

      real :: G0, G1, G2  !> standard heating model coefficients in cgs units
      real :: Lambda_0    !> power law cooling model coefficient in cgs units

      namelist /THERMAL/ thermal_active, heat_model, Lambda_0, alpha_cool, Teq, G0, G1, G2, x_ion, cfl_coolheat, isochoric, scheme, cool_model, cool_curve, cool_file, d_isochoric

      if (code_progress < PIERNIK_INIT_MPI) call die("[thermal:init_thermal] mpi not initialized.")

#ifdef VERBOSE
      if (master) call printinfo("[thermal:init_thermal] Commencing thermal module initialization")
#endif /* VERBOSE */

      thermal_active = .True.
      cool_model     = 'power_law'
      cool_curve     = 'power_law'
      cool_file      = ''
      heat_model     = 'G012'
      scheme         = 'EIS'
      alpha_cool     = 1.0
      Teq            = 1000.0
      Lambda_0       = 1.0e-25
      G0             = 1.0e-25
      G1             = 1.0e-25
      G2             = 1.0e-27
      x_ion          = 1.0
      isochoric      = 1
      d_isochoric    = 1.0
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

         rbuff(1) = Lambda_0
         rbuff(2) = alpha_cool
         rbuff(3) = Teq
         rbuff(4) = G0
         rbuff(5) = G1
         rbuff(6) = G2
         rbuff(7) = x_ion
         rbuff(8) = cfl_coolheat
         rbuff(9) = d_isochoric

         lbuff(1) = thermal_active

         cbuff(1) = cool_model
         cbuff(2) = cool_curve
         cbuff(3) = cool_file
         cbuff(4) = heat_model
         cbuff(5) = scheme

         ibuff(1) = isochoric

      endif

      call piernik_MPI_Bcast(cbuff, cbuff_len)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(ibuff)

      if (slave) then

         cool_model     = cbuff(1)
         cool_curve     = cbuff(2)
         cool_file      = cbuff(3)
         heat_model     = cbuff(4)
         scheme         = cbuff(5)

         thermal_active = lbuff(1)

         Lambda_0       = rbuff(1)
         alpha_cool     = rbuff(2)
         Teq            = rbuff(3)
         G0             = rbuff(4)
         G1             = rbuff(5)
         G2             = rbuff(6)
         x_ion          = rbuff(7)
         cfl_coolheat   = rbuff(8)
         d_isochoric    = rbuff(9)

         isochoric      = ibuff(1)

      endif

      G0_heat = G0       * erg / sek * cm**3 / mH**2 * x_ion**2
      G1_heat = G1       * erg / sek         / mH    * x_ion
      G2_heat = G2       * erg / sek / cm**3
      L0_cool = Lambda_0 * erg / sek * cm**3 / mH**2 * x_ion**2

      call all_cg%reg_var('Temperature')          ! Make it cleaner
      itemp = qna%ind('Temperature')

   end subroutine init_thermal

   subroutine fit_cooling_curve

      use dataio_pub,  only: msg, warn, die, printinfo
      use func,        only: operator(.equals.)
      use mpisetup,    only: master
      use units,       only: cm, erg, sek, mH

      implicit none

      integer                                 :: nbins
      integer, parameter                      :: coolfile = 1
      real(kind=8), dimension(:), allocatable :: logT, lambda, cool, heat
      integer                                 :: i
      real                                    :: T, d1

      if (cool_model .ne. 'piecewise_power_law') return
      if (master) call printinfo('[thermal] Cooling & heating handled with a single cool - heat curve fitted with a piecewise power law function')

      if (cool_curve .eq. 'tabulated') then
         open(unit=coolfile, file=cool_file, action='read', status='old')
         read(coolfile,*) nbins
         if (master) then
            write(msg,'(3a,i8,a)') '[thermal] Reading ', trim(cool_file), ' file with ', nbins, ' points'
            call printinfo(msg)
         endif
         allocate(logT(nbins), lambda(nbins), cool(nbins), heat(nbins))
         do i = 1, nbins
            read(coolfile,*) logT(i), cool(i)
         enddo
      else
         nbins = 700
         allocate(logT(nbins), lambda(nbins), cool(nbins), heat(nbins))
         do i = 1, nbins
            logT(i)   = 1 + i * 7.0 / nbins  ! linear spacing between 10 and 10**8 K
         enddo
      endif

      d1 = 0.0
      Teql = 0.0
      do i = 1, nbins
         T = 10**logT(i)

         select case (cool_curve)
            case ('power_law')
               cool(i)   = L0_cool * (T / Teq )**alpha_cool
            case ('Heintz')
               if ((master) .and. (i .eq. 1)) call printinfo('[thermal] Heintz cooling function used. Cooling power law parameters not used.')
               cool(i) = (7.3 * 10.0**(-21) * exp(-118400/(T+1500)) + 7.9 * 10.0**(-27) * exp(-92/T) ) * erg / sek * cm**3 / mH**2 * x_ion**2
            case ('tabulated')
               cool(i) = cool(i) * erg / sek * cm**3 / mH**2 * x_ion**2
            case default
               call die('[init_thermal] Cooling curve function not implemented')
         end select

         if (isochoric .eq. 1) then
            d1 = d_isochoric
            heat(i) = G0_heat * d1**2 + G1_heat * d1 + G2_heat
         else if (isochoric .eq. 2) then
            write(msg,'(3a)') 'isobaric case is not working well around equilibrium temperature'
            if ((master) .and. (i .eq. 1)) call warn(msg)
            d1 = Teq / 10**logT(i)
            heat(i) = G0_heat * d1**2 + G1_heat * d1 + G2_heat
         endif

         lambda(i) = cool(i) - heat(i)/d1**2
         if (i .gt. 1) then
            if ((lambda(i-1) * lambda(i) .le. 0.0) .and. (Teql .equals. 0.0)) Teql = T
            if ((T .gt. Teql) .and. (Teql .gt. 0.0) .and. (lambda(i-1) * lambda(i) .lt. 0.0)) call die('[init_thermal] More than 1 Teql')
         endif
      enddo
      if (master) then
         write(msg, '(a,f10.2)') '[thermal] Equilibrium Temperature = ', Teql
         call printinfo(msg)
      endif

      nfuncs = 1000
      call fit_proc(nbins, logT, lambda)   ! Find nfuncs
      if (master) then
         write(msg, '(a,i4)') '[thermal] fit nfuncs = ', nfuncs
         call printinfo(msg)
      endif
      deallocate(Tref, alpha,lambda0)
      call fit_proc(nbins, logT, lambda)       ! Perform fit

      deallocate(logT, lambda, cool, heat)

   end subroutine fit_cooling_curve

   subroutine fit_proc(nbins, logT, lambda)

      use constants,  only: big
      use dataio_pub, only: die
      use func,       only: operator(.equals.)

      implicit none

      integer,                intent(in)    :: nbins
      real, dimension(nbins), intent(inout) :: logT, lambda
      integer                               :: i, j, k
      real                                  :: a, b, r, rlim
      real, dimension(nbins)                :: fit, loglambda
      logical                               :: eq_point

      rlim = 10.0**(-6)
      do i = 1, nbins
         if (lambda(i) .equals. 0.0) then
            loglambda(i) = -big
         else
            loglambda(i) = log10(abs(lambda(i)))
         endif
      enddo
      if (allocated(Tref)) deallocate(Tref,alpha,lambda0)
      allocate(Tref(nfuncs), alpha(nfuncs), lambda0(nfuncs))

      i = 1
      k = 0
      do j = 2, nbins
         if (logT(j) .equals. log10(Teql)) then                                       ! log(lambda) goes to -inf at T=Teql
            cycle
         else if ((logT(j-1) .le. log10(Teql)) .and. logT(j) .gt. log10(Teql)) then   ! Look for the point right after Teql
            if (isochoric .eq. 1) then                                                      ! Linear fit of lambda between the point right before Teql, and 0
               a =  - lambda(i)/ (log10(Teql) - logT(i))
               b = 0.0
               k = k + 1
               Tref(k) = 10**logT(i)
               alpha(k) = b
               lambda0(k) = a/log(10.0)/Teql
               lambda(j+1) = a * logT(j+1)                                            ! Reassign lambda after Teql to match the linear function
               loglambda(j+1) = log10(abs(a * logT(j+1)))
            endif
            i = j
         else
            a = (loglambda(j) - loglambda(i)) / (logT(j) - logT(i))
            b = loglambda(j) - a*logT(j)
            fit(:) = a*logT + b
            r = sum( abs((loglambda(i:j) - fit(i:j))/fit(i:j)) ) / (j-i+1)
            eq_point = .false.
            if (j .lt. nbins) then
               if (logT(j+1) .ge. log10(Teql) .and. logT(j) .lt. log10(Teql)) eq_point = .true.
            endif
            if ((r .gt. rlim) .or. (eq_point)) then
               k = k + 1
               if (k .gt. nfuncs) call die('[init_thermal]: too many piecewise functions')
               Tref(k) = 10**logT(i)
               if (i .gt. 1) then
                  if ((isochoric .eq. 2) .and. ((logT(i-1) .le. log10(Teql)) .and. logT(i) .gt. log10(Teql))) Tref(k) = Teql
               endif
               alpha(k) = a
               lambda0(k) = lambda(j)/abs(lambda(j)) * 10**(b+a*logT(i))
               i = j
            endif
         endif
      enddo

      nfuncs = k

   end subroutine fit_proc

   subroutine thermal_sources(dt)

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: xdim, ydim, zdim
      use dataio_pub,       only: msg, warn
      use fluidindex,       only: flind
      use fluidtypes,       only: component_fluid
      use func,             only: ekin, emag, operator(.equals.)
      use grid_cont,        only: grid_container
      use mpisetup,         only: master
      use named_array_list, only: wna
      use units,            only: kboltz, mH

      implicit none

      real,                    intent(in) :: dt
      real, dimension(:, :, :), pointer   :: ta, dens, ener
      real, dimension(:,:,:), allocatable :: int_ener, kinmag_ener
      real                                :: dt_cool, t1, tcool, cfunc, hfunc, esrc, diff, kbgmh, ikbgmh, Tnew
      integer                             :: ifl, i, x, y, z
      integer, dimension(3)               :: n

      type(cg_list_element),  pointer     :: cgl
      type(grid_container),   pointer     :: cg
      class(component_fluid), pointer     :: pfl

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do ifl = 1, flind%fluids
            pfl => flind%all_fluids(ifl)%fl
            kbgmh  = kboltz / (pfl%gam_1 * mH)
            ikbgmh = pfl%gam_1 * mH / kboltz

            dens => cg%w(wna%fi)%span(pfl%idn,cg%lhn)
            ener => cg%w(wna%fi)%span(pfl%ien,cg%lhn)
            ta   => cg%q(itemp)%span(cg%lhn)

            n = shape(ta)
            allocate(int_ener(n(1),n(2),n(3)), kinmag_ener(n(1),n(2),n(3)))
            int_ener = 0.0
            if (pfl%has_energy) then
               kinmag_ener = ekin(cg%w(wna%fi)%span(pfl%imx,cg%lhn), cg%w(wna%fi)%span(pfl%imy,cg%lhn), cg%w(wna%fi)%span(pfl%imz,cg%lhn), dens)
               if (pfl%is_magnetized) then
                  kinmag_ener = kinmag_ener + emag(cg%w(wna%bi)%span(xdim,cg%lhn), cg%w(wna%bi)%span(ydim,cg%lhn), cg%w(wna%bi)%span(zdim,cg%lhn))
               endif
               int_ener = ener - kinmag_ener
            endif

            select case (scheme)

               case ('Explicit')
                  if (cool_model .eq. 'piecewise_power_law') then
                     write(msg,'(3a)') 'Warning: Make sure you are not using the heating in the cooling curve.'
                     if (master) call warn(msg)
                  endif
                  write(msg,'(3a)') 'Warning: substepping with a different timestep for every cell in the Explicit scheme leads to perturbations. Take a very small cfl_coolheat (~10^-6) or use a constant timestep.'
                  if (master) call warn(msg)
                  do x = 1, n(1)
                     do y = 1, n(2)
                        do z = 1, n(3)
                           call cool(ta(x,y,z), cfunc)
                           call heat(dens(x,y,z), hfunc)
                           esrc = dens(x,y,z)**2*cfunc + hfunc
                           dt_cool = min(dt, cfl_coolheat*abs(1./(esrc/int_ener(x,y,z))))
                           t1 = 0.0
                           do while (t1 .lt. dt)
                              call cool(ta(x,y,z), cfunc)
                              esrc = dens(x,y,z)**2*cfunc + hfunc
                              int_ener(x,y,z) = int_ener(x,y,z) + esrc * dt_cool
                              ener(x,y,z) = ener(x,y,z) + esrc * dt_cool
                              ta(x,y,z) = ikbgmh * int_ener(x,y,z) / dens(x,y,z)
                              t1 = t1 + dt_cool
                              if (t1 + dt_cool .gt. dt) dt_cool = dt - t1
                           enddo
                        enddo
                     enddo
                  enddo

               case ('EIS')
                  do x = 1, n(1)
                     do y = 1, n(2)
                        do z = 1, n(3)
                           do i = 1, nfuncs
                              if (i .eq. nfuncs) then
                                 if (ta(x,y,z) .ge. Tref(i))   tcool = kbgmh * ta(x,y,z) / (dens(x,y,z) * abs(lambda0(i)) * (ta(x,y,z)/Tref(i))**alpha(i))
                              else if (i .eq. 1) then
                                 if (ta(x,y,z) .le. Tref(i+1)) tcool = kbgmh * ta(x,y,z) / (dens(x,y,z) * abs(lambda0(i)) * (ta(x,y,z)/Tref(i))**alpha(i))
                              else if ((ta(x,y,z) .ge. Tref(i)) .and. (ta(x,y,z) .le. Tref(i+1))) then
                                 if (alpha(i) .equals. 0.0) then
                                    diff = MAX(abs(ta(x,y,z) - Teql), 0.000001)
                                    tcool = kbgmh * ta(x,y,z) / (abs(lambda0(i)) * diff * dens(x,y,z))
                                 else
                                    tcool = kbgmh * ta(x,y,z) / (dens(x,y,z) * abs(lambda0(i)) * (ta(x,y,z)/Tref(i))**alpha(i))
                                 endif
                              endif
                           enddo
                           dt_cool = min(dt, tcool/10.0)
                           t1 = 0.0
                           do while (t1 .lt. dt)
                              ta(x,y,z) = int_ener(x,y,z) * ikbgmh / dens(x,y,z)
                              call temp_EIS(tcool, dt_cool, pfl%gam, kbgmh, ta(x,y,z), dens(x,y,z), Tnew)
                              int_ener(x,y,z) = dens(x,y,z) * kbgmh * Tnew
                              ener(x,y,z) = kinmag_ener(x,y,z) + int_ener(x,y,z)

                              t1 = t1 + dt_cool
                              if (t1 + dt_cool .gt. dt) dt_cool = dt - t1
                           enddo
                        enddo
                     enddo
                  enddo

               case ('EE')
                  if (cool_model .eq. 'piecewise_power_law') then
                     write(msg,'(3a)') 'Warning: Make sure you are not using the heating in both the explicit and EI schemes.'
                     if (master) call warn(msg)
                  endif
                  do x = 1, n(1)
                     do y = 1, n(2)
                        do z = 1, n(3)
                           tcool = kbgmh * ta(x,y,z) / (dens(x,y,z) * abs(L0_cool) * (ta(x,y,z)/Teq)**alpha_cool)
                           dt_cool = min(dt, tcool/100.0)
                           t1 = 0.0
                           do while (t1 .lt. dt)
                              ta(x,y,z) = int_ener(x,y,z) * ikbgmh / dens(x,y,z)
                              call temp_EIS(tcool, dt_cool, pfl%gam, kbgmh, ta(x,y,z), dens(x,y,z), Tnew)
                              int_ener(x,y,z) = dens(x,y,z) * kbgmh * Tnew
                              ener(x,y,z) = kinmag_ener(x,y,z) + int_ener(x,y,z)

                              call heat(dens(x,y,z), hfunc)
                              int_ener(x,y,z) = int_ener(x,y,z) + hfunc * dt_cool
                              ener(x,y,z)     = ener(x,y,z)     + hfunc * dt_cool

                              t1 = t1 + dt_cool
                              if (t1 + dt_cool .gt. dt) dt_cool = dt - t1
                           enddo
                        enddo
                     enddo
                  enddo

               case default
                  write(msg,'(3a)') 'scheme: ',scheme,' not implemented'
                  if (master) call warn(msg)

            end select

            deallocate(int_ener, kinmag_ener)
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine thermal_sources

   subroutine cool(temp, coolf)

      use dataio_pub, only: msg, warn
      use mpisetup,   only: master

      implicit none

      real, intent(in)  :: temp
      real, intent(out) :: coolf
      integer           :: i

      select case (cool_model)
         case ('power_law')
            coolf = -L0_cool * (temp/Teq)**(alpha_cool)
         case ('piecewise_power_law')
            coolf = 0.0
            do i = 1, nfuncs
               if (i .eq. nfuncs) then
                  if (temp .ge. Tref(i)) then
                     coolf = - lambda0(i) * (temp/Tref(i))**alpha(i)
                  endif
               else if (i .eq. 1) then
                  if (temp .le. Tref(i+1)) then
                     coolf = - lambda0(i) * (temp/Tref(i))**alpha(i)
                  endif
               else
                  if ((temp .ge. Tref(i)) .and. (temp .le. Tref(i+1))) then
                     coolf = - lambda0(i) * (temp/Tref(i))**alpha(i)
                  endif
               endif
            enddo
         case ('null')
            return
         case default
            write(msg,'(3a)') 'Cool model: ',cool_model,' not implemented'
            if (master) call warn(msg)
      end select

   end subroutine cool

   subroutine heat(dens, heatf)

      use dataio_pub, only: msg, warn
      use mpisetup,   only: master

      implicit none

      real, intent(in)  :: dens
      real, intent(out) :: heatf

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

   subroutine temp_EIS(tcool, dt, gamma, kbgmh, temp, dens, Tnew)

      use dataio_pub, only: msg, warn
      use func,       only: operator(.equals.)
      use mpisetup,   only: master

      implicit none

      real, intent(in)        :: tcool, dt, gamma, temp, dens, kbgmh
      real, intent(out)       :: Tnew
      real                    :: lambda1, T1, alpha0, Y0, tcool2, TN, iso2, diff
      integer                 :: i, sign
      real, dimension(nfuncs) :: Y

      select case (cool_model)

         case ('power_law')
            if (alpha_cool .equals. 1.0) then
               if (isochoric == 1) then
                  Tnew = temp * exp(-dt/tcool)                 !isochoric
               else
                  Tnew = temp * (1 - 1/gamma * dt/tcool)      !isobar
               endif
               !Tnew = Tnew + dt * (gamma-1) * mH * G0_heat / kboltz  !heating
            else
               if (isochoric == 1) then
                  Tnew = temp * (1 - (1-alpha_cool) * dt / tcool) **(1.0/(1-alpha_cool))             !isochoric
               else
                  Tnew = temp * (1 - (isochoric-alpha_cool)*dt / tcool / gamma) **(1/(isochoric-alpha_cool))    !isobar
               endif
               !Tnew = Tnew + dt * (gamma-1) * mH * G0_heat / kboltz   ! isochoric heating
               !Tnew = Tnew * sqrt(1 + 2 * dt * (gamma-1) * mH * G0_heat *dens / kboltz / Tnew / gamma)  ! isobar heating
            endif

         case ('piecewise_power_law')
            iso2 = 1.0
            if (isochoric .eq. 2) then
               iso2 = 1.0/gamma
            endif
            TN = 10**8
            Y = 0.0
            T1 = 0.0
            sign = 0.0
            lambda1 = 0.0
            diff = 0.0
            Y0 = 0.0
            Y(nfuncs) = - 1 / (isochoric-alpha(nfuncs)) * (TN/Tref(nfuncs))**(alpha(nfuncs)-isochoric) * (1 - (Tref(nfuncs)/TN)**(alpha(nfuncs)-isochoric))
            do i = nfuncs-1, 1, -1
               if (alpha(i) .equals. 0.0) then
                  Y(i) = Y(i+1) - lambda0(nfuncs)/lambda0(i) * (TN/Tref(nfuncs))**alpha(nfuncs) / TN * log((Teql - Tref(i)) / (Tref(i+1)-Teql))
               else
                  Y(i) = Y(i+1) - 1 / (isochoric-alpha(i)) * lambda0(nfuncs)/lambda0(i) * (TN/Tref(nfuncs))**alpha(nfuncs) * (Tref(i)/TN)**isochoric * (1 - (Tref(i)/Tref(i+1))**(alpha(i)-isochoric))
               endif
            enddo
            do i = 1, nfuncs
               if (i .eq. nfuncs) then
                  if ((temp .ge. Tref(i))) then
                     Y0 = Y(i) + 1/(isochoric-alpha(i)) * lambda0(nfuncs)/lambda0(i) * (TN/Tref(nfuncs))**alpha(nfuncs) * (Tref(i)/TN)**isochoric * (1 - (Tref(i)/temp)**(alpha(i)-isochoric))
                     T1 = Tref(i)
                     alpha0 = alpha(i)
                     lambda1 = lambda0(i)
                  endif
               else if (i .eq. 1) then
                  if (temp .le. Tref(i+1)) then
                     Y0 = Y(i) + 1/(isochoric-alpha(i)) * lambda0(nfuncs)/lambda0(i) * (TN/Tref(nfuncs))**alpha(nfuncs) * (Tref(i)/TN)**isochoric * (1 - (Tref(i)/temp)**(alpha(i)-isochoric))
                     T1 = Tref(i)
                     lambda1 = lambda0(i)
                     alpha0 = alpha(i)
                  endif
               else
                  if ((temp .ge. Tref(i)) .and. (temp .le. Tref(i+1))) then
                     if (alpha(i) .equals. 0.0) then
                        if (temp .gt. Teql) then
                           sign = 1
                        else
                           sign = -1
                        endif
                        diff = MAX(abs(temp-Teql), 0.000001)
                        Y0 = Y(i) + lambda0(nfuncs)/lambda0(i) * (TN/Tref(nfuncs))**alpha(nfuncs) / TN * log((abs(Teql - Tref(i)) / diff))
                     else
                        Y0 = Y(i) + 1/(isochoric-alpha(i)) * lambda0(nfuncs)/lambda0(i) * (TN/Tref(nfuncs))**alpha(nfuncs) * (Tref(i)/TN)**isochoric * (1 - (Tref(i)/temp)**(alpha(i)-isochoric))
                     endif
                     T1 = Tref(i)
                     alpha0 = alpha(i)
                     lambda1 = lambda0(i)
                  endif
               endif
            enddo
            if (alpha0 .equals. 0.0) then
               tcool2 = kbgmh * temp / (lambda1 * diff * dens)
               tcool2 = min(tcool2, 1.0*10**6)
               Y0 = Y0 + (temp/TN) * lambda0(nfuncs)/lambda1 * (TN/Tref(nfuncs))**alpha(nfuncs) / diff * dt/tcool2
            else
               tcool2 = kbgmh * temp / (lambda1 * (temp/T1)**alpha0 * dens)
               Y0 = Y0 + (temp/TN)**isochoric * lambda0(nfuncs)/lambda1 * (TN/Tref(nfuncs))**alpha(nfuncs) * (T1/temp)**alpha0 * dt/tcool2 * iso2
            endif
            do i = 1, nfuncs
               if (i .eq. nfuncs) then
                  if ((temp .ge. Tref(i))) then
                     Tnew = Tref(i) * (1 - (isochoric-alpha(i)) * lambda0(i)/lambda0(nfuncs) * (Tref(nfuncs)/TN)**alpha(nfuncs) * (TN/Tref(i))**isochoric * (Y0 - Y(i)) )**(1/(isochoric-alpha(i)))
                  endif
               else if (i .eq. 1) then
                  if (temp .le. Tref(i+1)) then
                     Tnew = Tref(i) * (1 - (isochoric-alpha(i)) * lambda0(i)/lambda0(nfuncs) * (Tref(nfuncs)/TN)**alpha(nfuncs) * (TN/Tref(i))**isochoric * (Y0 - Y(i)) )**(1/(isochoric-alpha(i)))
                  endif
               else
                  if ((temp .ge. Tref(i)) .and. (temp .le. Tref(i+1))) then
                     if (alpha0 .equals. 0.0) then
                        Tnew = Teql + sign * (Teql-Tref(i)) * exp(-TN * (Tref(nfuncs)/TN)**alpha(nfuncs) * lambda0(i)/lambda0(nfuncs) * (Y0 - Y(i)))
                     else
                        Tnew = Tref(i) * (1 - (isochoric-alpha(i)) * lambda0(i)/lambda0(nfuncs) * (Tref(nfuncs)/TN)**alpha(nfuncs) * (TN/Tref(i))**isochoric * (Y0 - Y(i)) )**(1/(isochoric-alpha(i)))
                     endif
                  endif
               endif
            enddo
            if (Tnew .lt. 100.0) then
               Tnew = 100.0                        ! To improve
            endif

         case ('null')
            return

         case default
            write(msg,'(3a)') 'Cool model: ',cool_model,' not implemented'
            if (master) call warn(msg)

       end select

   end subroutine temp_EIS

end module thermal
