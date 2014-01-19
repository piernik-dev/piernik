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

module initproblem

   use constants, only: ndims

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real                   :: clump_mass, clump_K, clump_r, epsC, epsM, dtrig, dmax
   real, dimension(ndims) :: clump_pos, clump_vel, clump_pos_xyz
   logical                :: crashNotConv, exp_speedup, verbose
   integer(kind=4)        :: maxitC, maxitM
   integer, parameter     :: REL_CALC = 1, REL_SET = REL_CALC + 1

   namelist /PROBLEM_CONTROL/  clump_mass, clump_pos, clump_vel, clump_K, clump_r, epsC, epsM, maxitC, maxitM, crashNotConv, exp_speedup, verbose, dtrig

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use dataio_user, only: user_vars_hdf5, user_attrs_wr, user_attrs_rd
      use user_hooks,  only: late_initial_conditions, problem_domain_update

      implicit none

      late_initial_conditions => sg_late_init
      problem_domain_update => sg_dist_to_edge
      user_attrs_wr => sg_attrs_wr
      user_attrs_rd => sg_attrs_rd
      user_vars_hdf5 => sg_vars

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use constants,  only: GEO_XYZ, GEO_RPZ, xdim, ydim, zdim
      use dataio_pub, only: nh   ! QA_WARN required for diff_nml
      use dataio_pub, only: die, warn
      use domain,     only: dom
      use func,       only: operator(.notequals.)
      use mpisetup,   only: rbuff, ibuff, lbuff, master, slave, piernik_MPI_Bcast

      implicit none

      ! namelist default parameter values
      clump_mass   = 1.0e10                !< Mass of the clump
      clump_pos(:) = dom%C_(:)             !< position of the clump
      clump_vel(:) = 0.                    !< uniform velocity field
      clump_K      = 1.                    !< polytropic constant K for p = K rho**gamma formula
      clump_r      = 0.                    !< initial radius of the clump
      epsC         = 1.e-5                 !< tolerance limit for energy level change
      epsM         = 1.e-10                !< tolerance limit for clump mass change
      dtrig        = -0.001                !< fraction of dmax that causes extension of computational domain
      maxitC       = 100                   !< iteration limit for energy level
      maxitM       = 100                   !< iteration limit for clump mass
      crashNotConv = .true.                !< Crash if unable to converge initial conditions
      exp_speedup  = .false.               !< Use exponential fit to speed up convergence
      verbose      = .false.               !< Turn on some extra messages
      !\todo add rotation

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=PROBLEM_CONTROL, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL")
         read(nh%cmdl_nml,nml=PROBLEM_CONTROL, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1) = clump_mass
         rbuff(2:4) = clump_vel
         rbuff(5) = clump_K
         rbuff(6) = clump_r
         rbuff(7) = epsC
         rbuff(8) = epsM
         rbuff(9) = dtrig
         rbuff(10:12) = clump_pos

         ibuff(1) = maxitC
         ibuff(2) = maxitM

         lbuff(1) = crashNotConv
         lbuff(2) = exp_speedup
         lbuff(3) = verbose

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         clump_mass   = rbuff(1)
         clump_vel    = rbuff(2:4)
         clump_K      = rbuff(5)
         clump_r      = rbuff(6)
         epsC         = rbuff(7)
         epsM         = rbuff(8)
         dtrig        = rbuff(9)
         clump_pos    = rbuff(10:12)

         maxitC       = ibuff(1)
         maxitM       = ibuff(2)

         crashNotConv = lbuff(1)
         exp_speedup  = lbuff(2)
         verbose      = lbuff(3)

      endif

      if (clump_mass <= 0.) call die("[initproblem:read_problem_par] Negative mass of the clump.")
      if (clump_K <= 0.)    call die("[initproblem:read_problem_par] Negative polytropic constant.")
#ifdef ISO
      call die("[initproblem:read_problem_par] Isothermal EOS not supported.")
#endif /* ISO */

      clump_r = max(clump_r, maxval(dom%L_(:)/dom%n_d(:), mask=dom%has_dir(:)))

      dmax = 0.

      select case (dom%geometry_type)
         case (GEO_XYZ)
            clump_pos_xyz = clump_pos
         case (GEO_RPZ)
            clump_pos_xyz = [ clump_pos(xdim)*cos(clump_pos(ydim)), clump_pos(xdim)*sin(clump_pos(ydim)), clump_pos(zdim) ]
            if (any(clump_vel(xdim:ydim) .notequals. 0.)) then
               clump_vel(xdim:ydim) = 0.
               if (master) call warn("[initproblem:read_problem_par] clump_vel(xdim:ydim) /= 0 not implemented for GEO_RPZ. Forced to 0.")
            endif
         case default
            call die("[initproblem:read_problem_par] Unsupported geometry")
      end select

   end subroutine read_problem_par

!-----------------------------------------------------------------------------
!  Iterate density field to get it self-consistent with potential, rotation etc.
!  Based on: Hachisu, ApJS, 61, 479
!
! BEWARE: hardcoded numbers
!
   subroutine problem_initial_conditions

      use cg_leaves,         only: leaves
      use cg_list,           only: cg_list_element
      use constants,         only: pi, xdim, ydim, zdim, pSUM, pMIN, pMAX, GEO_XYZ, GEO_RPZ
      use dataio_pub,        only: msg, die, warn, printinfo
      use domain,            only: dom
      use fluidindex,        only: flind
      use fluidtypes,        only: component_fluid
      use func,              only: ekin, emag, operator(.equals.), operator(.notequals.)
      use global,            only: smalld, smallei, t
      use mpisetup,          only: master, piernik_MPI_Allreduce
      use multigrid_gravity, only: multigrid_solve_grav
      use units,             only: newtong

      implicit none

      class(component_fluid), pointer :: fl
      real, parameter                 :: virial_tol = 0.02
      integer, parameter              :: LOW=1, HIGH=LOW+1, TRY=3, NLIM=3
      integer                         :: i, j, k, tt,tmax, iC, iC_cg, iM, il, ih, jl, jh, kl, kh
      logical                         :: doneC, doneM
      real, dimension(LOW:HIGH)       :: Cint, totME
      character(len=HIGH)             :: ind
      real, dimension(TRY)            :: Cint_try, totME_try
      character(len=TRY)              :: i_try
      real                            :: Cint_old, Clim, Clim_old
      real, dimension(NLIM)           :: Clast
      integer, parameter              :: cmt_len=9 ! length of the " Exp warp" string
      character(len=cmt_len)          :: Ccomment
      real                            :: t_save
      real                            :: Msph
      type(cg_list_element),  pointer :: cgl

      fl => flind%ion

      t_save = t
      Cint_old = huge(1.)

      iC = 0
      totME(1) = clump_mass / (4./3. * pi * clump_r**3)

      Msph = 0.
      cgl => leaves%first
      do while (associated(cgl))
         associate (cg => cgl%cg)

         iC_cg = 0
         cg%b(:,    :, :, :) = 0.
         cg%u(fl%idn, :, :, :) = smalld
         cg%u(fl%ien, :, :, :) = smallei

         ! Initialize density with uniform sphere
         il = cg%ie+1
         ih = cg%is-1
         do i = cg%is, cg%ie
            if (abs(cg%x(i) - clump_pos(xdim)) <= clump_r) then
               il = min(i, il)
               ih = max(i, ih)
            endif
         enddo
         jl = cg%je+1
         jh = cg%js-1
         do j = cg%js, cg%je
            select case (dom%geometry_type)
               case (GEO_XYZ)
                  if (abs(cg%y(j) - clump_pos(ydim)) <= clump_r) then
                     jl = min(j, jl)
                     jh = max(j, jh)
                  endif
               case (GEO_RPZ)
                  if (cg%x(cg%is-1)*abs(cg%y(j) - clump_pos(ydim)) <= clump_r) then
                     jl = min(j, jl)
                     jh = max(j, jh)
                  endif
               case default
                  call die("[initproblem:problem_initial_conditions] Unsupported geometry")
            end select
         enddo
         kl = cg%ke+1
         kh = cg%ks-1
         do k = cg%ks, cg%ke
            if (abs(cg%z(k) - clump_pos(zdim)) <= clump_r) then
               kl = min(k, kl)
               kh = max(k, kh)
            endif
         enddo

         do k = kl, kh
            do j = jl, jh
               do i = il, ih
                  select case (dom%geometry_type)
                     case (GEO_XYZ)
                        if (sum(([cg%x(i), cg%y(j), cg%z(k)] -clump_pos(:))**2) < clump_r**2) then
                           cg%u(fl%idn, i, j, k) = totME(1)
                           iC_cg = iC_cg + 1
                        endif
                     case (GEO_RPZ)
                        if (sum(([cg%x(i)*cos(cg%y(j)), cg%x(i)*sin(cg%y(j)), cg%z(k)] - clump_pos_xyz(:))**2) < clump_r**2) then
                           cg%u(fl%idn, i, j, k) = totME(1)
                           iC_cg = iC_cg + 1
                        endif
                  end select
               enddo
            enddo
         enddo

         iC = iC + iC_cg
         Msph = Msph + iC_cg * totME(1) * cg%dvol

         end associate
         cgl => cgl%nxt
      enddo

      call piernik_MPI_Allreduce (iC,   pSUM)
      call piernik_MPI_Allreduce (Msph, pSUM)
      if (master .and. verbose) then
         write(msg,'(a,es13.7,a,i7,a)')"[initproblem:problem_initial_conditions] Starting with uniform sphere with M = ", Msph, " (", iC, " cells)"
         call printinfo(msg, .true.)
      endif

      ! Find C - the level of enthalpy at which density vanishes
      iC = 1
      doneC = .false.
      Clast(:) = 0. ; Clim = 0. ; Clim_old = 0.
      Ccomment = ''

      do while (.not. doneC)

         t = iC * sqrt(tiny(1.0)) ! trick to allow solution extrapolation in multigrid_solve_grav

         call multigrid_solve_grav([fl%idn])
         Cint = [ huge(1.), -huge(1.) ]
         cgl => leaves%first
         do while (associated(cgl))
            associate (cg => cgl%cg)

            if (exp_speedup .and. (Clim_old .notequals. 0.)) then ! extrapolate potential assuming exponential convergence (extremely risky)
               if (abs(1. - Clim/Clim_old) < min(sqrt(epsC), 100.*epsC, 0.01)) then
                  cg%sgp = (cg%sgp*cg%hgpot - cg%gpot**2)/(cg%sgp + cg%hgpot - 2.*cg%gpot)
                  Ccomment = ' Exp warp'
                  Clast(:) = 0. ; Clim = 0.
               else
                  Ccomment = ''
               endif
            endif
            cg%hgpot = cg%gpot
            cg%gpot  = cg%sgp

            Cint = [ min(Cint(LOW),  minval(cg%sgp(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), mask=cg%leafmap)), &
                 &   max(Cint(HIGH), maxval(cg%sgp(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), mask=cg%leafmap)) ] ! rotation will modify this
            end associate
            cgl => cgl%nxt
         enddo

         call piernik_MPI_Allreduce (Cint(LOW),  pMIN)
         call piernik_MPI_Allreduce (Cint(HIGH), pMAX)

         call totalMEnthalpic(Cint(LOW),  totME(LOW),  REL_CALC)
         call totalMEnthalpic(Cint(HIGH), totME(HIGH), REL_CALC)
         ind = '-+'

         if (iC > 1) then ! try previous C
            tmax = HIGH
            i_try(1:tmax) = '12'
            do tt = 1, tmax ! replicated code
               call totalMEnthalpic(Cint_try(tt), totME_try(tt), REL_CALC)
               if (totME_try(tt) > clump_mass .and. totME_try(tt) < totME(HIGH)) then
                  Cint(HIGH)     = Cint_try(tt)
                  totME(HIGH)    = totME_try(tt)
                  ind(HIGH:HIGH) = i_try(tt:tt)
               else if (totME_try(tt) < clump_mass .and. totME_try(tt) > totME(LOW)) then
                  Cint(LOW)    = Cint_try(tt)
                  totME(LOW)   = totME_try(tt)
                  ind(LOW:LOW) = i_try(tt:tt)
               endif
            enddo
            if (Cint(LOW) > Cint(HIGH)) then
               Cint_try(LOW:HIGH) = Cint(LOW:HIGH)
               Cint(LOW:HIGH) = Cint_try(HIGH:LOW:-1)
               totME_try(LOW:HIGH) = totME(LOW:HIGH)
               totME(LOW:HIGH) = totME_try(HIGH:LOW:-1)
            endif
         endif

         if (master .and. verbose) then
            write(msg,'(2(a,i4),2(a,2es15.7),2a)')"[initproblem:problem_initial_conditions] iter = ",iC,"/",0," dM= ",totME-clump_mass, " C= ", Cint, " ind = ",ind
            call printinfo(msg, .true.)
         endif

         ! Find what C corresponds to mass M in current potential well
         iM = 1
         doneM = .false.

         do while (.not. doneM)

            if (clump_mass > totME(LOW) .and. clump_mass < totME(HIGH)) then
               ind = '<>'
               tmax = TRY
               Cint_try(LOW)  = Cint(LOW) + (Cint(HIGH) - Cint(LOW))*(clump_mass - totME(LOW))/(totME(HIGH) - totME(LOW) ) ! secant
               Cint_try(HIGH) = (Cint(LOW) + Cint(HIGH))/2.                                                                ! bisection
               Cint_try(TRY)  = 2*Cint_try(1) - Cint(LOW)                                                                  ! 2*overshoot secant
               i_try = 'sbS'
               do tt = 1, tmax
                  call totalMEnthalpic(Cint_try(tt), totME_try(tt), REL_CALC)
                  if (totME_try(tt) >= clump_mass .and. totME_try(tt) < totME(HIGH)) then
                     Cint(HIGH)     = Cint_try(tt)
                     totME(HIGH)    = totME_try(tt)
                     ind(HIGH:HIGH) = i_try(tt:tt)
                  else if (totME_try(tt) <= clump_mass .and. totME_try(tt) > totME(LOW)) then
                     Cint(LOW)    = Cint_try(tt)
                     totME(LOW)   = totME_try(tt)
                     ind(LOW:LOW) = i_try(tt:tt)
                  endif
               enddo
            else ! the interval does not contain target mass M
               tmax = HIGH
               if (clump_mass > totME(HIGH)) then     ! amoeba crawls up
                  Cint_try(LOW)  = Cint(HIGH)
                  Cint_try(HIGH) = 3*Cint(HIGH) - 2*Cint(LOW)
                  i_try(1:tmax)  = '+u'
               else if (clump_mass < totME(LOW)) then ! amoeba crawls down
                  Cint_try(HIGH) = Cint(LOW)
                  Cint_try(LOW)  = 3*Cint(LOW) - 2*Cint(HIGH)
                  i_try(1:tmax)  = 'd-'
               endif
               do tt = LOW, HIGH
                  Cint(tt) = Cint_try(tt)
                  call totalMEnthalpic(Cint(tt), totME(tt), REL_CALC)
                  ind(tt:tt) = i_try(tt:tt)
               enddo
            endif

            if (master .and. verbose) then
               write(msg,'(2(a,i4),2(a,2es15.7),2a)')"[initproblem:problem_initial_conditions] iter = ",iC,"/",iM," dM= ",totME-clump_mass, " C= ", Cint, " ind = ",ind
               call printinfo(msg, .true.)
            endif
            tt = LOW
            if (abs(1. - totME(LOW)/clump_mass) > abs(1. - totME(HIGH)/clump_mass)) tt = HIGH
            if (abs(1. - totME(tt)/clump_mass) < epsM) doneM = .true.

            iM = iM + 1
            if (iM > maxitM .and. .not. doneM) then
               if (crashNotConv) then
                  call die("[initproblem:problem_initial_conditions] M-iterations not converged.")
               else
                  if (master) call warn("[initproblem:problem_initial_conditions] M-iterations not converged. Continue anyway.")
                  doneM = .true.
               endif
            endif

         enddo
         call totalMEnthalpic(Cint(tt), totME(tt), REL_SET)
         call virialCheck(huge(1.0))

         Clast(1:NLIM-1) = Clast(2:NLIM)
         Clast(NLIM) = Cint(tt)
         if (any(Clast(:) .equals. 0.)) then
            if (master) then
               write(msg,'(a,i4,2(a,es15.7),a)')"[initproblem:problem_initial_conditions] iter = ",iC,"     M=",totME(tt), " C=", Cint(tt), Ccomment
               call printinfo(msg, .true.)
            endif
         else
            if (Clim .notequals. 0.) Clim_old = Clim
            ! exponential estimate: \lim C \simeq \frac{C_{t} C_{t-2} - C_{t-1}^2}{C_{t} - 2 C_{t-1} + C{t-2}}
            Clim = (Clast(NLIM)*Clast(NLIM-2) - Clast(NLIM-1)**2)/(Clast(NLIM) - 2.*Clast(NLIM-1) + Clast(NLIM-2))
            if (master) then
               write(msg, '(a,i4,4(a,es15.7))')"[initproblem:problem_initial_conditions] iter = ",iC,"     M=",totME(tt), " C=", Cint(tt), " Clim=", Clim, " Clim-C=",Clim-Cint(tt)
               call printinfo(msg, .true.)
            endif
         endif

         if (abs(1. - Cint(tt)/Cint_old) < epsC) doneC = .true.
         Cint_old = Cint(tt)
         Cint_try(LOW:HIGH) = Cint ! try them as first guesses in next iteration

         iC = iC + 1
         if (iC > maxitC .and. .not. doneC) then
            if (crashNotConv) then
               call die("[initproblem:problem_initial_conditions] C-iterations not converged.")
            else
               if (master) call warn("[initproblem:problem_initial_conditions] C-iterations not converged. Continue anyway.")
               doneC = .true.
            endif
         endif

      enddo

      call virialCheck(virial_tol)

      ! final touch
      t = t_save ! restore initial time
      call multigrid_solve_grav([fl%idn])
      cgl => leaves%first
      do while (associated(cgl))
         associate (cg => cgl%cg)
         cg%gpot = cg%sgp

         where (cg%u(fl%idn, :, :, :) < smalld) cg%u(fl%idn, :, :, :) = smalld
         cg%u(fl%imx, :, :, :) = clump_vel(xdim) * cg%u(fl%idn,:,:,:)
         cg%u(fl%imy, :, :, :) = clump_vel(ydim) * cg%u(fl%idn,:,:,:)
         cg%u(fl%imz, :, :, :) = clump_vel(zdim) * cg%u(fl%idn,:,:,:)
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  cg%u(fl%ien,i,j,k) = max(smallei, presrho(cg%u(fl%idn, i, j, k)) / fl%gam_1                              + &
                       &              ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k)) + &
                       &              emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k)))
               enddo
            enddo
         enddo

         dmax = max(dmax, maxval(cg%u(fl%idn, cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke), mask=cg%leafmap))

         end associate
         cgl => cgl%nxt
      enddo
      call piernik_MPI_Allreduce(dmax, pMAX)

      if (master) then
         write(msg, '(a,g13.7)')"[initproblem:problem_initial_conditions] Relaxation finished. Largest orbital period: ",2.*pi*sqrt( (minval(dom%L_(:))/2.)**3/(newtong * clump_mass) )
         call printinfo(msg, .true.)
         write(msg, '(a,g13.7)')"[initproblem:problem_initial_conditions] Peak density ", dmax
         call printinfo(msg)
      endif

   end subroutine problem_initial_conditions

!-------------------------------------------------------------------------------
! check the value of | 2T + W + 3P | / | W |. Should be small

   subroutine virialCheck(tol)

      use cg_leaves,     only: leaves
      use cg_list,       only: cg_list_element
      use constants,     only: pSUM, GEO_XYZ, GEO_RPZ
      use dataio_pub,    only: msg, die, warn, printinfo
      use domain,        only: dom
      use fluidindex,    only: flind
      use grid_cont,     only: grid_container
      use mpisetup,      only: master, piernik_MPI_Allreduce
      use multigridvars, only: grav_bnd, bnd_isolated

      implicit none

      real, intent(in)                :: tol

      integer                         :: i, j, k
      integer, parameter              :: nTWP = 3
      real, dimension(nTWP)           :: TWP, TWPcg
      real                            :: vc
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

      TWP(:) = 0.

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         TWPcg(:) = 0.
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  if (cg%leafmap(i, j, k)) then
                     select case (dom%geometry_type)
                        case (GEO_XYZ)
                           ! TWPcg(1) = TWPcg(1) + cg%u(flind%ion%idn, i, j, k) * 0.                !T, will be /= 0. for rotating clump
                           TWPcg(2) = TWPcg(2) + cg%u(flind%ion%idn, i, j, k) * cg%sgp(i, j, k) * 0.5 !W
                           TWPcg(3) = TWPcg(3) + presrho(cg%u(flind%ion%idn, i, j, k))             !P
                        case (GEO_RPZ)
                           ! TWPcg(1) = TWPcg(1) + cg%u(flind%ion%idn, i, j, k) * 0. *cg%x(i)                 !T, will be /= 0. for rotating clump
                           TWPcg(2) = TWPcg(2) + cg%u(flind%ion%idn, i, j, k) * cg%sgp(i, j, k) * 0.5 *cg%x(i) !W
                           TWPcg(3) = TWPcg(3) + presrho(cg%u(flind%ion%idn, i, j, k)) *cg%x(i)              !P
                        case default
                           call die("[initproblem:virialCheck] Unsupported geometry")
                     end select
                  endif
               enddo
            enddo
         enddo
         TWP = TWP + TWPcg * cg%dvol

         cgl => cgl%nxt
      enddo

      call piernik_MPI_Allreduce(TWP, pSUM)

      vc = abs(2.*TWP(1) + TWP(2) + 3*TWP(3))/abs(TWP(2))
      if (master .and. (verbose .or. tol < 1.0)) then
         write(msg,'(a,es15.7,a,3es15.7,a)')"[initproblem:virialCheck] VC=",vc, " TWP=(",TWP(:),")"
         call printinfo(msg, .true.)
      endif

      if (vc > tol .and. grav_bnd == bnd_isolated) then
         if (master) then
            if (3*abs(TWP(3)) < abs(TWP(2))) then
               call warn("[initproblem:virialCheck] Virial imbalance occured because the clump is not resolved.")
            else
               call warn("[initproblem:virialCheck] Virial imbalance occured because the clump overfills the domain.")
            endif
         endif
         if (crashNotConv) call die("[initproblem:virialCheck] Virial defect too high.")
      endif

   end subroutine virialCheck

!-------------------------------------------------------------------------------
! Try two values of integral constant C and return corresponding masses

   subroutine totalMEnthalpic(C, totME, mode)

      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use constants,   only: pSUM, GEO_XYZ, GEO_RPZ
      use dataio_pub,  only: die
      use domain,      only: dom
      use fluidindex,  only: flind
      use grid_cont,   only: grid_container
      use mpisetup,    only: piernik_MPI_Allreduce

      implicit none

      real,    intent(in)            :: C
      real,    intent(out)           :: totME
      integer, intent(in)            :: mode

      integer                        :: i, j, k
      real                           :: rho, totMEcg
      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg

      totME = 0.

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         totMEcg = 0.
         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie
                  if (cg%leafmap(i, j, k)) then
                     select case (dom%geometry_type)
                        case (GEO_XYZ)
                           select case (mode)
                              case (REL_CALC)
                                 totMEcg = totMEcg + rhoH(h(C, cg%sgp(i,j,k)))
                              case (REL_SET)
                                 rho = rhoH(h(C, cg%sgp(i,j,k)))
                                 cg%u(flind%ion%idn, i, j, k) = rho
                                 totMEcg = totMEcg + rho
                           end select
                        case (GEO_RPZ)
                           select case (mode)
                              case (REL_CALC)
                                 totMEcg = totMEcg + rhoH(h(C, cg%sgp(i,j,k))) *cg%x(i)
                              case (REL_SET)
                                 rho = rhoH(h(C, cg%sgp(i,j,k)))
                                 cg%u(flind%ion%idn, i, j, k) = rho
                                 totMEcg = totMEcg + rho *cg%x(i)
                           end select
                        case default
                           call die("[initproblem:totalMEnthalpic] Unsupported geometry")
                     end select
                  endif
               enddo
            enddo
         enddo

         totME = totME + totMEcg * cg%dvol

         cgl => cgl%nxt
      enddo

      call piernik_MPI_Allreduce(totME, pSUM)

   end subroutine totalMEnthalpic

!-------------------------------------------------------------------------------
! calculate enthalpy

   real function h(C, Phi)

      implicit none

      real, intent(in) :: C, Phi

      h = C - Phi ! rotation will be included here

   end function h

!-------------------------------------------------------------------------------
! find pressure corresponding to given density (EOS dependent)

   real function presrho(rho)

      use fluidindex, only: flind

      implicit none

      real, intent(in) :: rho

      presrho = clump_K *  rho ** flind%ion%gam

   end function presrho

!-------------------------------------------------------------------------------
! find density corresponding to given enthalpy (EOS dependent)

   real function rhoH(H)

      use fluidindex, only: flind
      use global,     only: smalld

      implicit none

      real, intent(in) :: H

      rhoH = smalld

      if (H > 0.) rhoH = ( (1. - 1./flind%ion%gam) * H / clump_K)**(1./flind%ion%gam_1)

      rhoH = max(smalld, rhoH)

   end function rhoH

!-------------------------------------------------------------------------------
! find enthalpy corresponding to given density (EOS dependent)

   real function hRho(rho)

      use fluidindex, only: flind

      implicit none

      real, intent(in) :: rho

      hRho = clump_K * flind%ion%gam / flind%ion%gam_1 * rho ** flind%ion%gam_1

   end function hRho

!> \brief Performa late initialization of the cg added after domain expansion

   subroutine sg_late_init

      use cg_list,        only: cg_list_element
      use cg_list_dataop, only: expanded_domain
      use constants,      only: xdim, ydim, zdim
      use dataio_pub,     only: die
      use fluidindex,     only: flind
      use func,           only: ekin, emag
      use global,         only: smalld, smallei

      implicit none

      type(cg_list_element),  pointer :: cgl

      integer :: p, i, j, k

      cgl => expanded_domain%first
      do while (associated(cgl))
         if (cgl%cg%is_old) call die("[initproblem:sg_late_init] Old piece on a new list")
         do p = 1, flind%energ
            associate ( fl => flind%ion )
            cgl%cg%u(fl%idn, :, :, :) = smalld
            cgl%cg%u(fl%imx, :, :, :) = clump_vel(xdim) * cgl%cg%u(fl%idn,:,:,:)
            cgl%cg%u(fl%imy, :, :, :) = clump_vel(ydim) * cgl%cg%u(fl%idn,:,:,:)
            cgl%cg%u(fl%imz, :, :, :) = clump_vel(zdim) * cgl%cg%u(fl%idn,:,:,:)
            cgl%cg%b(:,    :, :, :) = 0.
            ! cgl%cg%sgp ?
            do k = cgl%cg%ks, cgl%cg%ke
               do j = cgl%cg%js, cgl%cg%je
                  do i = cgl%cg%is, cgl%cg%ie
                     cgl%cg%u(fl%ien, i, j, k) = max(smallei, presrho(cgl%cg%u(fl%idn, i, j, k)) / fl%gam_1                              + &
                          &                     ekin(cgl%cg%u(fl%imx,i,j,k), cgl%cg%u(fl%imy,i,j,k), cgl%cg%u(fl%imz,i,j,k), cgl%cg%u(fl%idn,i,j,k)) + &
                          &                     emag(cgl%cg%b(xdim,i,j,k), cgl%cg%b(ydim,i,j,k), cgl%cg%b(zdim,i,j,k)))
                  enddo
               enddo
            enddo
            end associate
         enddo
         cgl => cgl%nxt
      enddo

   end subroutine sg_late_init

!> \brief Find hov close it the clump to the external edges and call expansion routine if necessary

   subroutine sg_dist_to_edge

      use cg_leaves,     only: leaves
      use cg_level_base, only: base
      use cg_list,       only: cg_list_element
      use constants,     only: xdim, ydim, zdim, LO, HI
      use domain,        only: dom
      use fluidindex,    only: iarr_all_dn

      implicit none

      type(cg_list_element),  pointer :: cgl
      real, dimension(xdim:zdim, LO:HI) :: ddist
      integer :: i
      integer, parameter :: iprox = 2

      if (dtrig < 0.) return

      ddist = huge(1.)
      cgl => leaves%first
      do while (associated(cgl))
         if (any(cgl%cg%ext_bnd)) then
            !> \todo roll it to a nested loop
            if (dom%has_dir(xdim)) then
               if (cgl%cg%ext_bnd(xdim, LO)) then
                  do i = cgl%cg%is, cgl%cg%ie
                     if (any(cgl%cg%u(iarr_all_dn, i, :, :) > dmax*dtrig)) then
                        ddist(xdim, LO) = min(ddist(xdim, LO), (cgl%cg%x(i) - cgl%cg%fbnd(xdim, LO))/cgl%cg%dx)
                        exit
                     endif
                  enddo
               endif
               if (cgl%cg%ext_bnd(xdim, HI)) then
                  do i = cgl%cg%ie, cgl%cg%is, -1
                     if (any(cgl%cg%u(iarr_all_dn, i, :, :) > dmax*dtrig)) then
                        ddist(xdim, HI) = min(ddist(xdim, HI), (cgl%cg%fbnd(xdim, HI) - cgl%cg%x(i))/cgl%cg%dx)
                        exit
                     endif
                  enddo
               endif
            endif

            if (dom%has_dir(ydim)) then
               if (cgl%cg%ext_bnd(ydim, LO)) then
                  do i = cgl%cg%js, cgl%cg%je
                     if (any(cgl%cg%u(iarr_all_dn, :, i, :) > dmax*dtrig)) then
                        ddist(ydim, LO) = min(ddist(ydim, LO), (cgl%cg%y(i) - cgl%cg%fbnd(ydim, LO))/cgl%cg%dy)
                        exit
                     endif
                  enddo
               endif
               if (cgl%cg%ext_bnd(ydim, HI)) then
                  do i = cgl%cg%je, cgl%cg%js, -1
                     if (any(cgl%cg%u(iarr_all_dn, :, i, :) > dmax*dtrig)) then
                        ddist(ydim, HI) = min(ddist(ydim, HI), (cgl%cg%fbnd(ydim, HI) - cgl%cg%y(i))/cgl%cg%dy)
                        exit
                     endif
                  enddo
               endif
            endif

            if (dom%has_dir(zdim)) then
               if (cgl%cg%ext_bnd(zdim, LO)) then
                  do i = cgl%cg%ks, cgl%cg%ke
                     if (any(cgl%cg%u(iarr_all_dn, :, :, i) > dmax*dtrig)) then
                        ddist(zdim, LO) = min(ddist(zdim, LO), (cgl%cg%z(i) - cgl%cg%fbnd(zdim, LO))/cgl%cg%dz)
                        exit
                     endif
                  enddo
               endif
               if (cgl%cg%ext_bnd(zdim, HI)) then
                  do i = cgl%cg%ke, cgl%cg%ks, -1
                     if (any(cgl%cg%u(iarr_all_dn, :, :, i) > dmax*dtrig)) then
                        ddist(zdim, HI) = min(ddist(zdim, HI), (cgl%cg%fbnd(zdim, HI) - cgl%cg%z(i))/cgl%cg%dz)
                        exit
                     endif
                  enddo
               endif
            endif
         endif
         cgl => cgl%nxt
      enddo

      !> \todo shrink the domain in the direction opposite to expansion (shift the domain in given direction)
      call base%expand(ddist(:,:) < iprox)

   end subroutine sg_dist_to_edge

!> \brief Write dmax to the restart file

   subroutine sg_attrs_wr(file_id)

      use dataio_pub, only: die
      use hdf5,       only: HID_T, SIZE_T
      use h5lt,       only: h5ltset_attribute_double_f

      implicit none

      integer(HID_T), intent(in) :: file_id

      integer(SIZE_T), parameter :: bufsize = 1
      integer(kind=4)            :: error

      call h5ltset_attribute_double_f(file_id, "/", "dmax", [dmax], bufsize, error)
      if (error<0) call die("[initproblem:sg_attrs_rd] Error writing dmax attribute")

   end subroutine sg_attrs_wr

!> \brief Read dmax from the restart file

   subroutine sg_attrs_rd(file_id)

      use constants,  only: I_ONE
      use dataio_pub, only: die
      use global,     only: smalld
      use hdf5,       only: HID_T
      use h5lt,       only: h5ltget_attribute_double_f

      implicit none

      integer(HID_T), intent(in) :: file_id

      integer(kind=4)        :: error
      real, dimension(I_ONE) :: rbuff

      call h5ltget_attribute_double_f(file_id, "/", "dmax", rbuff, error)
      dmax = rbuff(1)

      if (error<0) call die("[initproblem:sg_attrs_rd] Error reading dmax attribute")
      if (dmax <= smalld) call die("[initproblem:sg_attrs_rd] dmax <= smalld")

   end subroutine sg_attrs_rd

!> \brief

   subroutine sg_vars(var, tab, ierrh, cg)

      use constants,  only: pi, xdim, zdim
      use domain,     only: dom
      use fluidindex, only: flind
      use func,       only: ekin
      use grid_cont,  only: grid_container
      use units,      only: newtong

      implicit none

      character(len=*),               intent(in)    :: var
      real(kind=4), dimension(:,:,:), intent(inout) :: tab
      integer,                        intent(inout) :: ierrh
      type(grid_container), pointer,  intent(in)    :: cg

      real :: delx
      integer :: d, i, j, k

      ierrh = 0
      select case (trim(var))
         case ("nJ")
            ! Jeans wavelength divided by longest cell dimension n_J = lambda_J/dz = c_sound * sqrt(pi/(G rho))/dx
            delx = 0.
            do d = xdim, zdim
               if (dom%has_dir(d)) delx = max(delx, cg%dl(d))
            enddo
            do k = cg%ks, cg%ke
               do j = cg%js, cg%je
                  do i = cg%is, cg%ie
                     tab(i-cg%is+1, j-cg%js+1, k-cg%ks+1) = real(sqrt(pi/newtong/cg%u(flind%ion%idn, i, j, k)) * flind%ion%get_cs(i, j, k, cg%u, cg%b, cg%cs_iso2)/delx, kind=4)
                  enddo
               enddo
            enddo
         case ("cs")
            do k = cg%ks, cg%ke
               do j = cg%js, cg%je
                  do i = cg%is, cg%ie
                     tab(i-cg%is+1, j-cg%js+1, k-cg%ks+1) = real(flind%ion%get_cs(i, j, k, cg%u, cg%b, cg%cs_iso2), kind=4)
                  enddo
               enddo
            enddo
         case ("mach")
            do k = cg%ks, cg%ke
               do j = cg%js, cg%je
                  do i = cg%is, cg%ie
                     tab(i-cg%is+1, j-cg%js+1, k-cg%ks+1) = real(sqrt(2*ekin(cg%u(flind%ion%imx, i, j, k), cg%u(flind%ion%imy, i, j, k), cg%u(flind%ion%imz, i, j, k), cg%u(flind%ion%idn, i, j, k)**2)) / flind%ion%get_cs(i, j, k, cg%u, cg%b, cg%cs_iso2), kind=4)
                  enddo
               enddo
            enddo
         case default
            ierrh = -1
      end select

   end subroutine sg_vars

end module initproblem
