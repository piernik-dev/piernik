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

module initproblem

  use problem_pub, only: problem_name, run_id

  real               :: clump_mass, clump_pos_x, clump_pos_y, clump_pos_z, clump_vel_x, clump_vel_y, clump_vel_z, clump_K, clump_r, epsC, epsM
  logical            :: crashNotConv
  integer            :: maxitC, maxitM

  integer, parameter :: REL_CALC = 1, REL_SET = REL_CALC + 1

  namelist /PROBLEM_CONTROL/  problem_name, run_id, clump_mass, clump_vel_x, clump_vel_y, clump_vel_z, clump_K, clump_r, epsC, epsM, maxitC, maxitM, crashNotConv

contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use grid,     only : xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz
      use errh,     only : namelist_errh, die
      use mpisetup, only : cwd, ierr, rbuff, cbuff, ibuff, lbuff, proc, buffer_dim, comm, &
           &               MPI_CHARACTER, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_LOGICAL

      implicit none

      integer            :: ierrh
      character(LEN=100) :: par_file, tmp_log_file

      ! namelist default parameter values
      problem_name = 'selfgrav_clump'      !< The default problem name
      run_id       = '_'                   !< Auxiliary run identifier
      clump_mass   = 1.0e10                !< Mass of the clump
      clump_vel_x  = 0.                    !< X-velocity in the domain
      clump_vel_y  = 0.                    !< Y-velocity in the domain
      clump_vel_z  = 0.                    !< Z-velocity in the domain
      clump_K      = 1.                    !< polytropic constant K for p = K rho**gamma formula
      clump_r      = 0.                    !< initial radius of the clump
      epsC         = 1.e-5                 !< tolerance limit for energy level change
      epsM         = 1.e-10                !< tolerance limit for clump mass change
      maxitC       = 100                   !< iteration limit for energy level
      maxitM       = 100                   !< iteration limit for clump mass
      crashNotConv = .true.                !< Crash if unable to converge initial conditions
      !\todo add rotation

      if(proc == 0) then
         par_file = trim(cwd)//'/problem.par'
         tmp_log_file = trim(cwd)//'/tmp.log'

         open(1,file=par_file)
            read(unit=1,nml=PROBLEM_CONTROL,iostat=ierrh)
            call namelist_errh(ierrh,'PROBLEM_CONTROL')
         close(1)
         open(3, file=tmp_log_file, position='append')
            write(3,nml=PROBLEM_CONTROL)
            write(3,*)
         close(3)
      endif

      if (proc == 0) then

         cbuff(1) =  problem_name
         cbuff(2) =  run_id

         rbuff(1) = clump_mass
         rbuff(2) = clump_vel_x
         rbuff(3) = clump_vel_y
         rbuff(4) = clump_vel_z
         rbuff(5) = clump_K
         rbuff(6) = clump_r
         rbuff(7) = epsC
         rbuff(8) = epsM

         ibuff(1) = maxitC
         ibuff(2) = maxitM

         lbuff(1) = crashNotConv

      end if

      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      call MPI_BCAST(lbuff,    buffer_dim, MPI_LOGICAL,          0, comm, ierr)

      if (proc /= 0) then

         problem_name = cbuff(1)
         run_id       = cbuff(2)(1:3)

         clump_mass   = rbuff(1)
         clump_vel_x  = rbuff(2)
         clump_vel_y  = rbuff(3)
         clump_vel_z  = rbuff(4)
         clump_K      = rbuff(5)
         clump_r      = rbuff(6)
         epsC         = rbuff(7)
         epsM         = rbuff(8)

         maxitC       = ibuff(1)
         maxitM       = ibuff(2)

         crashNotConv = lbuff(1)

      endif

      if (clump_mass <= 0.) call die("[initproblem:read_problem_par] Negative mass of the clump.")
      if (clump_K <= 0.) call die("[initproblem:read_problem_par] Negative polytropic constant.")

#ifdef ISO
      call die("[initproblem:read_problem_par] Isothermal EOS not supported.")
#endif

      clump_pos_x = (xmax+xmin)/2.
      clump_pos_y = (ymax+ymin)/2.
      clump_pos_z = (zmax+zmin)/2.
      clump_r = max(clump_r, dx, dy, dz)

   end subroutine read_problem_par

!-----------------------------------------------------------------------------
!  Iterate density field to get it self-consistent with potential, rotation etc.
!  Based on: Hachisu, ApJS, 61, 479
!
! BEWARE: hardcoded numbers
!
   subroutine init_prob

      use mpisetup,      only : proc, smalld, smallei, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_MIN, MPI_MAX, MPI_SUM, comm, ierr
      use arrays,        only : u, b, mgp, gpot
      use constants,     only : fpiG, pi, newtong
      use grid,          only : xmin, xmax, ymin, ymax, zmin, zmax, x, y, z, nx, ny, nz, dx, dy, dz, is, ie, js, je, ks, ke
      use initionized,   only : gamma_ion, idni, imxi, imyi, imzi, ieni
      use dataio_public, only : tend
      use multigrid,     only : multigrid_solve
      use errh,          only : die

      implicit none

      integer, parameter             :: LOW=1, HIGH=LOW+1, TRY=3
      integer                        :: i, j, k, t, tmax, iC, iM, il, ih, jl, jh, kl, kh
      logical                        :: doneC, doneM
      real, dimension(LOW:HIGH)      :: Cint, totME
      character, dimension(LOW:HIGH) :: ind
      real, dimension(TRY)           :: Cint_try, totME_try
      character, dimension(TRY)      :: i_try
      real                           :: Cint_old = HUGE(1.)!, Cint_min

      u(imxi, :, :, :) = clump_vel_x
      u(imyi, :, :, :) = clump_vel_y
      u(imzi, :, :, :) = clump_vel_z
      b(:,    :, :, :) = 0.

      u(idni, :, :, :) = smalld
      u(ieni, :, :, :) = smallei

      ! Initialize with point source
      ! /todo use polytrope here to improve convergence for gamma_ion <= 1.4
      ! awk 'BEGIN {n=0; f=0.; t=1.; dx=1e-3; print 0., t, f; for (x=0.; x<=10. && t>0.; x+=dx) {if (x==0) f-=t**n*dx; else f-=(t**n + 2.*f/x)*dx; t+=f*dx; if (t>0) print x, t, f}}'
      ! rho = rho_c * theta**n
      ! x = sqrt(4 * pi * newtong / (K * (n + 1)) * rho_c**(1-1./n)) *r
      ! gamma = 1 + 1./n
      il = ie+1
      ih = is-1
      do i = is, ie
         if (abs(x(i) - clump_pos_x) <= clump_r) then
            il = min(i, il)
            ih = max(i, ih)
         end if
      end do
      jl = je+1
      jh = js-1
      do j = js, je
         if (abs(y(j) - clump_pos_y) <= clump_r) then
            jl = min(j, jl)
            jh = max(j, jh)
         end if
      end do
      kl = ke+1
      kh = ks-1
      do k = ks, ke
         if (abs(z(k) - clump_pos_z) <= clump_r) then
            kl = min(k, kl)
            kh = max(k, kh)
         end if
      end do

      iC = 0
      totME(1) = clump_mass / (4./3. * pi * clump_r**3)
      do k = kl, kh
         do j = jl, jh
            do i = il, ih
               if ((x(i)-clump_pos_x)**2 + (y(j)-clump_pos_y)**2 + (z(k)-clump_pos_z)**2 < clump_r**2) then
                  u(idni, i, j, k) = totME(1)
                  iC =iC + 1
               end if
            end do
         end do
      end do

      call MPI_AllReduce (MPI_IN_PLACE, iC, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
      if (proc == 0) write(*,'(a,es13.7,a,i7,a)')"[initproblem:init_prob] Starting with uniform sphere with M = ", iC*totME(1) * dx * dy * dz, " (", iC, " cells)"

      iC = 1
      doneC = .false.

      do while (.not. doneC)

         call multigrid_solve(u(idni,:,:,:))

         Cint = [ minval(mgp(is:ie,js:je,ks:ke)), maxval(mgp(is:ie,js:je,ks:ke)) ] ! rotation will modify this

         call MPI_AllReduce (MPI_IN_PLACE, Cint(LOW),  1, MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
         call MPI_AllReduce (MPI_IN_PLACE, Cint(HIGH), 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)
         !Cint_min = Cint(LOW)

         call totalMEnthalpic(Cint(LOW),  totME(LOW),  REL_CALC)
         call totalMEnthalpic(Cint(HIGH), totME(HIGH), REL_CALC)
         ind = [ '-', '+' ]

         if (iC > 1) then ! try previous C
            tmax = HIGH
            i_try(1:tmax) = [ '1', '2' ]
            do t = 1, tmax ! replicated code
               call totalMEnthalpic(Cint_try(t), totME_try(t), REL_CALC)
               if (totME_try(t) > clump_mass .and. totME_try(t) < totME(HIGH)) then
                  Cint(HIGH)  = Cint_try(t)
                  totME(HIGH) = totME_try(t)
                  ind(HIGH)   = i_try(t)
               else if (totME_try(t) < clump_mass .and. totME_try(t) > totME(LOW)) then
                  Cint(LOW)   = Cint_try(t)
                  totME(LOW)  = totME_try(t)
                  ind(LOW)    = i_try(t)
               end if
            end do
            if (Cint(LOW) > Cint(HIGH)) then
               Cint_try(LOW:HIGH) = Cint(LOW:HIGH)
               Cint(LOW:HIGH) = Cint_try(HIGH:LOW:-1)
               totME_try(LOW:HIGH) = totME(LOW:HIGH)
               totME(LOW:HIGH) = totME_try(HIGH:LOW:-1)
            end if
         end if

         if (proc == 0) write(*,'(2(a,i4),2(a,2es15.7),3a)')"[initproblem:init_prob] iter = ",iC,"/",0," dM= ",totME-clump_mass, " C= ", Cint, " ind = ",ind

         iM = 1
         doneM = .false.

         do while (.not. doneM)

            if (clump_mass > totME(LOW) .and. clump_mass < totME(HIGH)) then
               ind = [ '<', '>' ]
               tmax = TRY
               Cint_try(LOW)  = Cint(LOW) + (Cint(HIGH) - Cint(LOW))*(clump_mass - totME(LOW))/(totME(HIGH) - totME(LOW) ) ! secant
               Cint_try(HIGH) = (Cint(LOW) + Cint(HIGH))/2.                                                                ! bisection
               Cint_try(TRY)  = 2*Cint_try(1) - Cint(LOW)                                                                  ! 2*overshoot secant
               i_try = [ 's', 'b', 'S' ]
               do t = 1, tmax
                  call totalMEnthalpic(Cint_try(t), totME_try(t), REL_CALC)
                  if (totME_try(t) > clump_mass .and. totME_try(t) < totME(HIGH)) then
                     Cint(HIGH)  = Cint_try(t)
                     totME(HIGH) = totME_try(t)
                     ind(HIGH)   = i_try(t)
                  else if (totME_try(t) < clump_mass .and. totME_try(t) > totME(LOW)) then
                     Cint(LOW)  = Cint_try(t)
                     totME(LOW) = totME_try(t)
                     ind(LOW)   = i_try(t)
                  end if
               end do
            else
               tmax = HIGH
               if (clump_mass > totME(HIGH)) then     ! amoeba crawls up
                  Cint_try(1)   = Cint(HIGH)
                  Cint_try(2)   = 3*Cint(HIGH) - 2*Cint(LOW)
                  i_try(1:tmax) = [ '+', 'u' ]
               else if (clump_mass < totME(LOW)) then ! amoeba crawls down
                  Cint_try(2)   = Cint(LOW)
                  Cint_try(1)   = 3*Cint(LOW) - 2*Cint(HIGH)
                  i_try(1:tmax) = [ 'd', '-' ]
               end if
               do t = LOW, HIGH
                  Cint(t) = Cint_try(t)
                  call totalMEnthalpic(Cint(t), totME(t), REL_CALC)
                  ind(t) = i_try(t)
               end do
            end if

            if (proc == 0) write(*,'(2(a,i4),2(a,2es15.7),3a)')"[initproblem:init_prob] iter = ",iC,"/",iM," dM= ",totME-clump_mass, " C= ", Cint, " ind = ",ind

            t = LOW
            if (abs(1. - totME(LOW)/clump_mass) > abs(1. - totME(HIGH)/clump_mass)) t = HIGH
            if (abs(1. - totME(t)/clump_mass) < epsM) doneM = .true.

            iM = iM + 1
            if (iM > maxitM .and. .not. doneM) then
               if (crashNotConv) then
                  call die("[initproblem:init_prob] M-iterations not converged")
               else
                  write(*,'(a)')"[initproblem:init_prob] M-iterations not converged. Continue anyway."
                  doneM = .true.
               end if
            end if

         end do
         call totalMEnthalpic(Cint(t), totME(t), REL_SET)
         call virialCheck(huge(1.0))

         if (proc == 0) write(*,'(a,i4,2(a,es15.7))')"[initproblem:init_prob] iter = ",iC,"     M=",totME(t), " C=", Cint(t)

         if (abs(1. - Cint(t)/Cint_old) < epsC) doneC = .true.
         Cint_old = Cint(t)
         Cint_try(1:2) = Cint ! try them as first guesses in next iteration

         iC = iC + 1
         if (iC > maxitC .and. .not. doneC) then
            if (crashNotConv) then
               call die("[initproblem:init_prob] C-iterations not converged")
            else
               write(*,'(a)')"[initproblem:init_prob] C-iterations not converged. Continue anyway."
               doneC = .true.
            end if
         end if

      end do

      call virialCheck(0.01)

      ! final touch
      call multigrid_solve(u(idni,:,:,:))
      gpot = mgp
      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               u(idni,i,j,k) = max(u(idni,i,j,k), smalld)
               u(ieni,i,j,k) = presrho(u(idni, i, j, k)) / (gamma_ion-1.0) + &
                    &          0.5 * sum(u(imxi:imzi,i,j,k)**2,1) / u(idni,i,j,k)      + &
                    &          0.5 * sum(b(:,i,j,k)**2,1)
               u(ieni,i,j,k) = max(u(ieni,i,j,k), smallei)
            enddo
         enddo
      enddo

      if (proc == 0) write(*,'(a,g13.7)')"[initproblem:init_prob] Relaxation finished. Largest orbital period: ",2.*pi*sqrt( (min(xmax-xmin, ymax-ymin, zmax-zmin)/2.)**3/(newtong * clump_mass) )

   end subroutine init_prob

!-------------------------------------------------------------------------------
! check the value of | 2T + W + 3P | / | W |. Should be small

   subroutine virialCheck(tol)

      use grid,        only : is, ie, js, je, ks, ke, dx, dy, dz
      use arrays,      only : u, mgp
      use mpisetup,    only : proc, comm, ierr, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
      use initionized, only : idni
      use errh,        only : die

      implicit none

      real, intent(in)      :: tol

      integer               :: i, j, k
      integer, parameter    :: nTWP = 3
      real, dimension(nTWP) :: TWP
      real                  :: vc

      TWP(:) = 0.

      do k = ks, ke
         do j = js, je
            do i = is, ie
!               TWP(1) = TWP(1) + u(idni, i, j, k) * 0.                !T, will be /= 0. for rotating clump
               TWP(2) = TWP(2) + u(idni, i, j, k) * mgp(i, j, k) * 0.5 !W
               TWP(3) = TWP(3) + presrho(u(idni, i, j, k))         !P
            end do
         end do
      end do

      call MPI_AllReduce (MPI_IN_PLACE, TWP, nTWP, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)

      TWP = TWP * dx * dy * dz
      vc = abs(2.*TWP(1) + TWP(2) + 3*TWP(3))/abs(TWP(2))
      if (proc == 0) write(*,'(a,es15.7,a,3es15.7,a)')"[initproblem:virialCheck] VC=",vc, " TWP=(",TWP(:),")"

      if (vc > tol) then
         if (3*abs(TWP(3)) < abs(TWP(2))) then
            write(*,'(a)')"[initproblem:virialCheck] Virial imbalance occured because the clump is not resolved"
         else
            write(*,'(a)')"[initproblem:virialCheck] Virial imbalance occured because the clump overfills the domain"
         end if
         if (crashNotConv) call die("[initproblem:virialCheck] Virial defect too high.")
      end if

   end subroutine virialCheck

!-------------------------------------------------------------------------------
! Try two values of integral constant C and return corresponding masses

   subroutine totalMEnthalpic(C, totME, mode)

      use mpisetup,    only : smalld, comm, ierr, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
      use arrays,      only : mgp, u
      use grid,        only : is, ie, js, je, ks, ke, dx, dy, dz
      use initionized, only : idni

      implicit none

      real,    intent(in)  :: C
      real,    intent(out) :: totME
      integer, intent(in)  :: mode

      integer              :: i, j, k
      real                 :: rho

      totME = 0.

      do k = ks, ke
         do j = js, je
            do i = is, ie
               select case (mode)
               case (REL_CALC)
                  totME     = totME     + rhoH(h(C,     mgp(i,j,k)))
               case (REL_SET)
                  rho = rhoH(h(C, mgp(i,j,k)))
                  u(idni, i, j, k) = rho
                  totME = totME + rho
               end select
            end do
         end do
      end do

      call MPI_AllReduce (MPI_IN_PLACE, totME, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)

      totME     = totME     * dx * dy * dz

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

      use initionized, only : gamma_ion

      implicit none

      real, intent(in) :: rho

      presrho = clump_K *  rho ** gamma_ion

   end function presrho

!-------------------------------------------------------------------------------
! find density corresponding to given enthalpy (EOS dependent)

   real function rhoH(H)

      use initionized, only : gamma_ion
      use mpisetup,    only : smalld

      implicit none

      real, intent(in) :: H

      rhoH = smalld

      if (H > 0.) rhoH = ( (1. - 1./gamma_ion) * H / clump_K)**(1./(gamma_ion - 1.))

      rhoH = max(smalld, rhoH)

   end function rhoH

!-------------------------------------------------------------------------------
! find enthalpy corresponding to given density (EOS dependent)

   real function hRho(rho)

      use initionized,   only : gamma_ion

      implicit none

      real, intent(in) :: rho

      hRho = clump_K * gamma_ion / (gamma_ion - 1.) * rho ** (gamma_ion - 1.)

   end function hRho

end module initproblem
