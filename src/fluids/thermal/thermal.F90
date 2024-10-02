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
   public ::  init_thermal, thermal_active, cfl_coolheat, thermal_sources, itemp, fit_cooling_curve, cleanup_thermal, calc_tcool

   character(len=cbuff_len)        :: cool_model, cool_curve, heat_model, scheme, cool_file
   logical                         :: thermal_active
   real                            :: alpha_cool, L0_cool, G0_heat, G1_heat, G2_heat, cfl_coolheat
   real                            :: Teq          !> cooling parameter
   real                            :: Teql         !> temperature of cooling / heating equilibrium
   integer(kind=4), protected      :: itemp = INVALID
   real                            :: x_ion        !> ionization degree
   integer(kind=4)                 :: isochoric    !> 1 for isochoric, 2 for isobaric
   real                            :: d_isochoric  ! constant density used in isochoric case
   real                            :: TN, ltntrna
   real, dimension(:), allocatable :: Tref, alpha, lambda0, Y
   integer                         :: nfuncs

contains

   subroutine init_thermal

      use bcast,            only: piernik_MPI_Bcast
      use cg_list_global,   only: all_cg
      use constants,        only: PIERNIK_INIT_MPI
      use dataio_pub,       only: code_progress, die, nh, printinfo, warn
      use mpisetup,         only: cbuff, lbuff, rbuff,ibuff,  master, slave
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

      thermal_active = .true.
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

      call all_cg%reg_var('Temperature')          ! Make it cleaner
      itemp = qna%ind('Temperature')
      if (.not. thermal_active) return

      G0_heat = G0       * erg / sek * cm**3 / mH**2 * x_ion**2
      G1_heat = G1       * erg / sek         / mH    * x_ion
      G2_heat = G2       * erg / sek / cm**3
      L0_cool = Lambda_0 * erg / sek * cm**3 / mH**2 * x_ion**2

      call fit_cooling_curve()

      if (scheme == 'Explicit') call warn('[thermal:init_thermal][scheme: Explicit] Warning: substepping with a different timestep for every cell in the Explicit scheme leads to perturbations. Take a very small cfl_coolheat (~10^-6) or use a constant timestep.')

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

      if (cool_model /= 'piecewise_power_law') return
      if (master) then
         call printinfo('[thermal:fit_cooling_curve] Cooling & heating handled with a single cool - heat curve fitted with a piecewise power law function')
         if (scheme == 'Explicit') call warn('[thermal:fit_cooling_curve][scheme: Explicit] Warning: Make sure you are not using the heating in the cooling curve.')
         if (scheme == 'EE')       call warn('[thermal:fit_cooling_curve][scheme: EE] Warning: Make sure you are not using the heating in both the explicit and EI schemes.')
      endif

      if (cool_curve == 'tabulated') then
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
               if ((master) .and. (i == 1)) call printinfo('[thermal:fit_cooling_curve] Heintz cooling function used. Cooling power law parameters not used.')
               cool(i) = (7.3 * 10.0**(-21) * exp(-118400/(T+1500)) + 7.9 * 10.0**(-27) * exp(-92/T) ) * erg / sek * cm**3 / mH**2 * x_ion**2
            case ('tabulated')
               cool(i) = cool(i) * erg / sek * cm**3 / mH**2 * x_ion**2
            case default
               call die('[init_thermal] Cooling curve function not implemented')
         end select

         if (isochoric == 1) then
            d1 = d_isochoric
            heat(i) = G0_heat * d1**2 + G1_heat * d1 + G2_heat
         else if (isochoric == 2) then
            write(msg,'(3a)') 'isobaric case is not working well around equilibrium temperature'
            if ((master) .and. (i == 1)) call warn(msg)
            d1 = Teq / 10**logT(i)
            heat(i) = G0_heat * d1**2 + G1_heat * d1 + G2_heat
         endif

         lambda(i) = cool(i) - heat(i)/d1**2
         if (i > 1) then
            if ((lambda(i-1) * lambda(i) <= 0.0) .and. (Teql .equals. 0.0)) Teql = T
            if ((T > Teql) .and. (Teql > 0.0) .and. (lambda(i-1) * lambda(i) < 0.0)) call die('[thermal:fit_cooling_curve] More than 1 Teql')
         endif
      enddo
      if (master) then
         write(msg, '(a,f10.2)') '[thermal] Equilibrium Temperature = ', Teql
         call printinfo(msg)
      endif

      call fit_proc(nbins, logT, lambda)   ! Find nfuncs and Perform fit
      deallocate(logT, lambda, cool, heat)

   end subroutine fit_cooling_curve

   subroutine fit_proc(nbins, logT, lambda)

      use constants,  only: big
      use dataio_pub, only: msg, printinfo !,die
      use func,       only: operator(.equals.)
      use mpisetup,   only: master

      implicit none

      integer,                intent(in)    :: nbins
      real, dimension(nbins), intent(inout) :: logT, lambda
      integer                               :: i, j, k, iter
      real                                  :: a, b, r, rlim
      real, dimension(nbins)                :: fit, loglambda
      logical                               :: eq_point, set_nfuncs, fill_array

      rlim = 10.0**(-6)
      do i = 1, nbins
         if (lambda(i) .equals. 0.0) then
            loglambda(i) = -big
         else
            loglambda(i) = log10(abs(lambda(i)))
         endif
      enddo

      do iter = 1, 2
         set_nfuncs = (iter == 1)
         fill_array = (iter == 2)
         if (fill_array) allocate(Tref(nfuncs), alpha(nfuncs), lambda0(nfuncs))

         i = 1
         k = 0
         do j = 2, nbins
            if (logT(j) .equals. log10(Teql)) then                                       ! log(lambda) goes to -inf at T=Teql
               cycle
            else if ((logT(j-1) <= log10(Teql)) .and. logT(j) > log10(Teql)) then   ! Look for the point right after Teql
               if (isochoric == 1) then                                                      ! Linear fit of lambda between the point right before Teql, and 0
                  a =  - lambda(i)/ (log10(Teql) - logT(i))
                  b = 0.0
                  k = k + 1
                  if (fill_array) then
                     Tref(k) = 10**logT(i)
                     alpha(k) = b
                     lambda0(k) = a/log(10.0)/Teql
                  endif
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
               if (j < nbins) then
                  if (logT(j+1) >= log10(Teql) .and. logT(j) < log10(Teql)) eq_point = .true.
               endif
               if ((r > rlim) .or. (eq_point)) then
                  k = k + 1
                  !if (k > nfuncs) call die('[init_thermal]: too many piecewise functions')
                  if (fill_array) then
                     Tref(k) = 10**logT(i)
                     if (i > 1) then
                        if ((isochoric == 2) .and. ((logT(i-1) <= log10(Teql)) .and. logT(i) > log10(Teql))) Tref(k) = Teql
                     endif
                     alpha(k) = a
                     lambda0(k) = lambda(j)/abs(lambda(j)) * 10**(b+a*logT(i))
                  endif
                  i = j
               endif
            endif
         enddo

         if (set_nfuncs) then
            nfuncs = k
            if (master) then
               write(msg, '(a,i4)') '[thermal] fit nfuncs = ', nfuncs
               call printinfo(msg)
            endif
         endif

         if (fill_array) then
            TN = 10**8
            ltntrna = lambda0(nfuncs) * (TN/Tref(nfuncs))**alpha(nfuncs)

            if (.false.) then ! this might be used in some future implementation
               allocate(Y(nfuncs))
               Y = 0.0
               Y(nfuncs) = - 1 / (isochoric-alpha(nfuncs)) * ((TN/Tref(nfuncs))**(alpha(nfuncs)-isochoric) - 1)
               do i = nfuncs-1, 1, -1
                  if (alpha(i) .equals. 0.0) then
                     Y(i) = Y(i+1) - ltntrna / lambda0(i) / TN * log((Teql - Tref(i)) / (Tref(i+1)-Teql))
                  else
                     Y(i) = Y(i+1) - ltntrna / lambda0(i) / (isochoric-alpha(i)) * (Tref(i)/TN)**isochoric * (1 - (Tref(i)/Tref(i+1))**(alpha(i)-isochoric))
                  endif
               enddo
            endif
         endif
      enddo

   end subroutine fit_proc

   subroutine thermal_sources(dt)

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: xdim, ydim, zdim
      use dataio_pub, only: msg, warn
      use fluidindex, only: flind
      use fluidtypes, only: component_fluid
      use func,       only: ekin, emag
      use grid_cont,  only: grid_container
      use mpisetup,   only: master
      use units,      only: kboltz, mH

      implicit none

      real,                    intent(in) :: dt
      real, dimension(:, :, :), pointer   :: ta, dens, ener
      real, dimension(:,:,:), allocatable :: kinmag_ener
      real                                :: dt_cool, t1, tcool, cfunc, hfunc, esrc, kbgmh, ikbgmh, Tnew, int_ener
      integer                             :: ifl, i, j, k
      integer, dimension(3)               :: n

      type(cg_list_element),  pointer     :: cgl
      type(grid_container),   pointer     :: cg
      class(component_fluid), pointer     :: pfl

      if (.not. thermal_active) return
      hfunc = huge(1.)  ! suppress spurious compiler warning triggered by -Wmaybe-uninitialized

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do ifl = 1, flind%fluids
            pfl => flind%all_fluids(ifl)%fl
            if (.not. pfl%has_energy) cycle
            kbgmh  = kboltz / (pfl%gam_1 * mH)
            ikbgmh = pfl%gam_1 * mH / kboltz

            dens => cg%u(pfl%idn,:,:,:)
            ener => cg%u(pfl%ien,:,:,:)
            ta   => cg%q(itemp)%arr(:,:,:)

            n = shape(ta)
            allocate(kinmag_ener(n(xdim),n(ydim),n(zdim)))
            kinmag_ener = ekin(cg%u(pfl%imx,:,:,:), cg%u(pfl%imy,:,:,:), cg%u(pfl%imz,:,:,:), dens)
            if (pfl%is_magnetized) kinmag_ener = kinmag_ener + emag(cg%b(xdim,:,:,:), cg%b(ydim,:,:,:), cg%b(zdim,:,:,:))

            select case (scheme)

               case ('Explicit')
                  do i = 1, n(xdim)
                     do j = 1, n(ydim)
                        do k = 1, n(zdim)
                           int_ener = ener(i,j,k) - kinmag_ener(i,j,k)
                           call cool(ta(i,j,k), cfunc)
                           call heat(dens(i,j,k), hfunc)
                           esrc = dens(i,j,k)**2 * cfunc + hfunc
                           dt_cool = min(dt, cfl_coolheat*abs(1./(esrc/int_ener)))
                           t1 = 0.0
                           do while (t1 < dt)
                              call cool(ta(i,j,k), cfunc)
                              esrc = dens(i,j,k)**2 * cfunc + hfunc
                              int_ener = int_ener + esrc * dt_cool
                              ener(i,j,k) = ener(i,j,k) + esrc * dt_cool
                              ta(i,j,k) = ikbgmh * int_ener / dens(i,j,k)
                              t1 = t1 + dt_cool
                              if (t1 + dt_cool > dt) dt_cool = dt - t1
                           enddo
                        enddo
                     enddo
                  enddo

               case ('EIS')
                  do i = 1, n(xdim)
                     do j = 1, n(ydim)
                        do k = 1, n(zdim)
                           int_ener = ener(i,j,k) - kinmag_ener(i,j,k)
                           ta(i,j,k) = int_ener * ikbgmh / dens(i,j,k)
                           if (ta(i,j,k) .lt. 10.0) ta(i,j,k) = 10.0
                           if (ta(i,j,k) .gt. 10.0**8) ta(i,j,k) = 10.0**8
                           call calc_tcool(ta(i,j,k), dens(i,j,k), kbgmh, tcool)
                           !tcool = 1.0
                           dt_cool = min(dt, tcool/10.0)
                           t1 = 0.0
                           do while (t1 < dt)
                              call temp_EIS(tcool, dt_cool, igamma(pfl%gam), kbgmh, ta(i,j,k), dens(i,j,k), Tnew)
                              !Tnew = ta(i,j,k)
                              int_ener    = dens(i,j,k) * kbgmh * Tnew
                              !if (int_ener .gt. 10.0) print *, 'wow!', i, j, k, int_ener, Tnew, ta(i,j,k), dens(i,j,k), tcool
                              ener(i,j,k) = kinmag_ener(i,j,k) + int_ener
                              ta(i,j,k) = Tnew
                              t1 = t1 + dt_cool
                              if (t1 + dt_cool > dt) dt_cool = dt - t1
                           enddo
                        enddo
                     enddo
                  enddo

               case ('EE')
                  do i = 1, n(xdim)
                     do j = 1, n(ydim)
                        do k = 1, n(zdim)
                           int_ener = ener(i,j,k) - kinmag_ener(i,j,k)
                           tcool    = kbgmh * ta(i,j,k) / (dens(i,j,k) * abs(L0_cool) * (ta(i,j,k)/Teq)**alpha_cool)
                           dt_cool  = min(dt, tcool/100.0)
                           t1 = 0.0
                           do while (t1 < dt)
                              ta(i,j,k) = int_ener * ikbgmh / dens(i,j,k)
                              call temp_EIS(tcool, dt_cool, igamma(pfl%gam), kbgmh, ta(i,j,k), dens(i,j,k), Tnew)
                              int_ener    = dens(i,j,k) * kbgmh * Tnew
                              ener(i,j,k) = kinmag_ener(i,j,k) + int_ener

                              call heat(dens(i,j,k), hfunc)
                              int_ener    = int_ener    + hfunc * dt_cool
                              ener(i,j,k) = ener(i,j,k) + hfunc * dt_cool

                              t1 = t1 + dt_cool
                              if (t1 + dt_cool > dt) dt_cool = dt - t1
                           enddo
                        enddo
                     enddo
                  enddo

               case default
                  write(msg,'(3a)') 'scheme: ',scheme,' not implemented'
                  if (master) call warn(msg)

            end select

            deallocate(kinmag_ener)
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine thermal_sources

   real function igamma(gam) result(fiso)

      implicit none

      real, intent(in) :: gam

      fiso = 1.0
      if (isochoric == 2) fiso = 1.0/gam

   end function igamma

   subroutine find_temp_bin(temp, ii)

      implicit none

      real,    intent(in)  :: temp
      integer, intent(out) :: ii
      integer              :: i

      ii = 0
      if (temp >= Tref(nfuncs)) then
         ii = nfuncs
      else if (temp < Tref(2)) then
         ii = 1
      else
         do i = 2, nfuncs - 1
            if ((temp >= Tref(i)) .and. (temp < Tref(i+1))) ii = i
         enddo
      endif

    end subroutine find_temp_bin

    subroutine calc_tcool(temp, dens, kbgmh, tcool)

      use func,        only: operator(.equals.)

      implicit none

      real,    intent(in)  :: temp, dens, kbgmh
      real,    intent(out) :: tcool
      integer              :: ii
      real                 :: alpha1, Tref1, lambda1, diff

      if (cool_model == 'piecewise_power_law') then
         call find_temp_bin(temp, ii)
         alpha1  = alpha(ii)
         Tref1   = Tref(ii)
         lambda1 = lambda0(ii)
      else
         alpha1  = alpha_cool
         Tref1   = Teq
         lambda1 = L0_cool
      endif

      if (alpha1 .equals. 0.0) then
         diff = max(abs(temp - Teql), 0.000001)
         tcool = kbgmh * temp / (dens * abs(lambda1) * diff)
      else
         tcool = kbgmh * temp / (dens * abs(lambda1) * (temp/Tref1)**alpha1)
      endif

    end subroutine calc_tcool

   subroutine cool(temp, coolf)

      use dataio_pub, only: msg, warn
      use mpisetup,   only: master

      implicit none

      real, intent(in)  :: temp
      real, intent(out) :: coolf
      integer           :: ii

      select case (cool_model)
         case ('power_law')
            coolf = -L0_cool * (temp/Teq)**alpha_cool
         case ('piecewise_power_law')
            coolf = 0.0
            call find_temp_bin(temp, ii)
            coolf = - lambda0(ii) * (temp/Tref(ii))**alpha(ii)
         case ('null')
            coolf = 0.0
         case default
            write(msg,'(3a)') 'Cool model: ',cool_model,' not implemented'
            if (master) call warn(msg)
            coolf = huge(1.)  ! this may crash the code or at least disturb the output to catch attention
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
            heatf = 0.0
         case default
            write(msg,'(3a)') 'Heat model: ',heat_model,' not implemented'
            if (master) call warn(msg)
            heatf = huge(1.)  ! this may crash the code or at least disturb the output to catch attention
      end select

   end subroutine heat

   subroutine temp_EIS(tcool, dt, fiso, kbgmh, temp, dens, Tnew)

      use dataio_pub, only: msg, warn
      use func,       only: operator(.equals.)
      use mpisetup,   only: master

      implicit none

      real, intent(in)  :: tcool, dt, fiso, temp, dens, kbgmh
      real, intent(out) :: Tnew
      real              :: lambda1, T1, alpha0, Y0f, tcool2, diff
      integer           :: ii

      select case (cool_model)

         case ('power_law')
            if (alpha_cool .equals. 1.0) then
               if (isochoric == 1) then
                  Tnew = temp * exp(-dt/tcool)             !isochoric
               else
                  Tnew = temp * (1 - fiso * dt/tcool)      !isobar
               endif
               !Tnew = Tnew + dt * (gamma-1) * mH * G0_heat / kboltz  !heating
            else
               Tnew = temp * (1 - (isochoric-alpha_cool) * dt / tcool * fiso)**(1.0/(isochoric-alpha_cool))
               !Tnew = Tnew + dt * (gamma-1) * mH * G0_heat / kboltz   ! isochoric heating
               !Tnew = Tnew * sqrt(1 + 2 * dt * (gamma-1) * mH * G0_heat *dens / kboltz / Tnew / gamma)  ! isobar heating
            endif

         case ('piecewise_power_law')
            call find_temp_bin(temp, ii)
            T1 = Tref(ii)
            alpha0 = alpha(ii)
            lambda1 = lambda0(ii)
            if (alpha0 .equals. 0.0) then
               diff = max(abs(temp-Teql), 0.000001)
               !Y0 = Y(ii) + ltntrna / lambda1 / TN * log((abs(Teql - T1) / diff))
               Y0f = log((abs(Teql - T1) / diff)) / lambda1
               tcool2 = kbgmh * temp / (lambda1 * diff * dens)
               tcool2 = min(tcool2, 1.0e6 * TN / ltntrna)
               Y0f = Y0f + temp / lambda1 / diff * dt/tcool2
               !Tnew = Teql - sign(1.0, Teql - temp) * (Teql-T1) * exp(-TN * lambda1 / ltntrna * (Y0 - Y(ii)))
               Tnew = Teql - sign(1.0, Teql - temp) * (Teql-T1) * exp(-lambda1 * Y0f)
            else
               Tnew = temp * (1 - (isochoric-alpha0) * sign(1.0,lambda1)* fiso * dt / tcool)**(1.0/(isochoric-alpha0))
            endif
         case ('null')
            return

         case default
            write(msg,'(3a)') 'Cool model: ',cool_model,' not implemented'
            if (master) call warn(msg)

       end select

       if (Tnew < 10.0) Tnew = 10.0                        ! To improve

   end subroutine temp_EIS

   subroutine cleanup_thermal

      implicit none

      if (allocated(Tref))    deallocate(Tref)
      if (allocated(alpha))   deallocate(alpha)
      if (allocated(lambda0)) deallocate(lambda0)
      if (allocated(Y))       deallocate(Y)

   end subroutine cleanup_thermal

end module thermal
