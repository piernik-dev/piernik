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

  real              :: clump_mass, clump_vel_x, clump_vel_y, clump_vel_z, clump_K

  integer, parameter :: REL_CALC = 1, REL_SET = REL_CALC + 1

  namelist /PROBLEM_CONTROL/  problem_name, run_id, clump_mass, clump_vel_x, clump_vel_y, clump_vel_z, clump_K

contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par

!      use grid,     only : xmin, xmax, ymin, ymax, zmin, zmax, nx, ny, nz
      use errh,     only : namelist_errh, die
      use mpisetup, only : cwd, ierr, rbuff, cbuff, ibuff, proc, buffer_dim, comm, &
           &               MPI_CHARACTER, MPI_DOUBLE_PRECISION, MPI_INTEGER
      use constants, only: pi

      implicit none

      integer :: ierrh
      character(LEN=100) :: par_file, tmp_log_file

      ! namelist default parameter values
      problem_name = 'selfgrav_clump'      !< The default problem name
      run_id       = '_'                   !< Auxiliary run identifier
      clump_mass   = 1.0e10                !< Mass of the clump
      clump_vel_x  = 0.                    !< X-velocity in the domain
      clump_vel_y  = 0.                    !< Y-velocity in the domain
      clump_vel_z  = 0.                    !< Z-velocity in the domain
      clump_K      = 1.                    !< polytropic constant K for p = K rho**gamma formula

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

      end if

      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
!      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (proc /= 0) then

         problem_name = cbuff(1)
         run_id       = cbuff(2)(1:3)

         clump_mass   = rbuff(1)
         clump_vel_x  = rbuff(2)
         clump_vel_y  = rbuff(3)
         clump_vel_z  = rbuff(4)
         clump_K      = rbuff(5)

      endif

      if (clump_mass <= 0.) call die("[initproblem:read_problem_par] Negative mass of the clump.")
      if (clump_K <= 0.) call die("[initproblem:read_problem_par] Negative polytropic constant.")

#ifdef ISO
      call die("[initproblem:read_problem_par] Isothermal EOS not supported.")
#endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------
!  Iterate density field to get it self-consistent with potential, rotation etc.
!  Based on: Hachisu, ApJS, 61, 479

   subroutine init_prob

      use mpisetup,      only : proc, smalld, smallei
      use arrays,        only : u, b, mgp
      use constants,     only : fpiG, pi, newtong
      use grid,          only : xmin, xmax, ymin, ymax, zmin, zmax, x, y, z, nx, ny, nz, dx, dy, dz
      use initionized,   only : gamma_ion, idni, imxi, imyi, imzi, ieni
      use dataio_public, only : tend
      use multigrid,     only : multigrid_solve

      implicit none

      logical, parameter :: crashNotConv = .false.
      real, parameter    :: epsC = 1.e-5, epsM = 1.e-10
      integer, parameter :: maxitC=100, maxitM = 100
      integer :: i, j, k, iC, iM
      logical :: doneC, doneM
      real    :: Cint, Cint_aux, Cint_old, totME, totME_aux, dC, dM, dMdC, cfac

      u(imxi, :, :, :) = clump_vel_x
      u(imyi, :, :, :) = clump_vel_y
      u(imzi, :, :, :) = clump_vel_z
      b(:,    :, :, :) = 0.

      u(idni, :, :, :) = smalld
      u(ieni, :, :, :) = smallei

      ! Start with a single-cell clump
      u(idni, nx/2, ny/2, nz/2) = (clump_mass - smalld*((xmax-xmin)*(ymax-ymin)*(zmax-zmin) - (dx * dy * dz)))/(dx * dy * dz)

      Cint = - newtong * clump_mass / (sqrt(dx**2 + dy**2 + dz**2)/2.)

      iC = 1
      doneC = .false.

      do while (.not. doneC)

         call multigrid_solve(u(idni,:,:,:))

         Cint_old = Cint
!AJG WARNING: hardcoded numbers go to cfac
         cfac = -1e-5
         iM = 1
         doneM = .false.

         do while (.not. doneM)

            Cint_aux = (1+cfac)*Cint
            call sim_totalMEnthalpic(Cint, Cint_aux, totME, totME_aux, REL_CALC)
            dMdC = (totME_aux - totME)/(Cint_aux - Cint)

!            write(*,'(a,5g20.10)')"ii stmE ",dMdC,totME_aux,totME,Cint_aux,Cint

            if (proc == 0) write(*,'(2(a,i3),3(a,es15.7))')"[initproblem:init_prob] iter = ",iC,"/",iM," M= ",totME, " C= ", Cint, " dM/dC= ", dMdC

            if (dMdC == 0.) then
               dC = -0.1 * Cint
            else
               if (totME/clump_mass > 2.) then
                  dM = -0.5*totME
               else if (totME/clump_mass > 0.1) then
                  dC = totME/clump_mass
                  dM = (clump_mass - totME)* ( 3. * dC / ( dC**3 + 2 ) ) ! crude safety factor
                  cfac = max(-0.01, -0.001 * abs(1. - totME/clump_mass))
               else
                  dM = 0.1*clump_mass
                  cfac = min(-1e-10, cfac / 2.)
               end if
               dC = dM / dMdC
            end if
            if (abs(dC/Cint) < epsM .and. abs(1. - totME/clump_mass) < epsM) doneM = .true.
            Cint = Cint + dC
!            write(*,'(a,6g20.10)')"ii d ",Cint,dC,dM,dMdC,totME/clump_mass,cfac

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
         call sim_totalMEnthalpic(Cint, Cint_aux, totME, totME_aux, REL_SET)
         call sim_virialCheck

         if (proc == 0) write(*,'(a,i3,2(a,es15.7))')"[initproblem:init_prob] iter = ",iC,"     M=",totME, " C=", Cint

         if (abs(1. - Cint/Cint_old) < epsC) doneC = .true.

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
      call multigrid_solve(u(idni,:,:,:))
      if (proc == 0) write(*,'(a)')"[initproblem:init_prob] Relaxation finished."

      ! final touch
      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               if (u(idni,i,j,k) < smalld) u(idni,i,j,k) = smalld
               u(ieni,i,j,k) = sim_presrho(u(idni, i, j, k)) / (gamma_ion-1.0) + &
                    &          0.5 * sum(u(imxi:imzi,i,j,k)**2,1) / u(idni,i,j,k)      + &
                    &          0.5 * sum(b(:,i,j,k)**2,1)
!               write(*,'(a,3i3,4g20.10)')"ii ",i,j,k,u(idni,i,j,k),sim_presrho(u(idni, i, j, k)),u(ieni,i,j,k),mgp(i,j,k)
            enddo
         enddo
      enddo


   end subroutine init_prob

!-------------------------------------------------------------------------------
! check the value of | 2T + W + 3P | / | W |. Should be small

   subroutine sim_virialCheck

      use grid,        only: is, ie, js, je, ks, ke, dx, dy, dz
      use arrays,      only: u, mgp
      use mpisetup,    only: proc, comm, ierr, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
      use initionized, only: idni

      implicit none

      integer                           :: i, j, k

      integer, parameter :: nTWP = 3
      real, dimension(nTWP) :: TWP

      TWP(:) = 0.

      do k = ks, ke
         do j = js, je
            do i = is, ie
!               TWP(1) = TWP(1) + u(idni, i, j, k) * 0.                !T
               TWP(2) = TWP(2) + u(idni, i, j, k) * mgp(i, j, k) * 0.5 !W
               TWP(3) = TWP(3) + sim_presrho(u(idni, i, j, k))         !P
            end do
         end do
      end do

      call MPI_AllReduce (MPI_IN_PLACE, TWP, nTWP, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)

      TWP = TWP * dx * dy * dz

      if (proc == 0) write(*,'(a,2(es15.7,a),3es15.7,a)')"[sim_virialCheck] VC=",abs(2.*TWP(1) + TWP(2) + 3*TWP(3))/abs(TWP(2)), " T/|W|=", TWP(1)/abs(TWP(2))," TWP=(",TWP(:),")"

   end subroutine sim_virialCheck

!-------------------------------------------------------------------------------
! Try two values of integral constant C and return corresponding masses

   subroutine sim_totalMEnthalpic(C, C_aux, totME, totME_aux, mode)

      use mpisetup,    only: smalld, comm, ierr, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM
      use arrays,      only: mgp, u
      use grid,        only: is, ie, js, je, ks, ke, dx, dy, dz
      use initionized, only: idni

      implicit none

      real,    intent(in)  :: C, C_aux
      real,    intent(out) :: totME, totME_aux
      integer, intent(in)  :: mode

      integer              :: i, j, k
      real                 :: rho

      totME = 0.
      totME_aux = 0.

      do k = ks, ke
         do j = js, je
            do i = is, ie
               select case (mode)
               case (REL_CALC)
                  totME     = totME     + sim_rhoH(sim_h(C,     mgp(i,j,k)))
                  totME_aux = totME_aux + sim_rhoH(sim_h(C_aux, mgp(i,j,k)))
               case (REL_SET)
                  rho = sim_rhoH(sim_h(C, mgp(i,j,k)))
                  u(idni, i, j, k) = max(smalld, rho)
                  totME = totME + rho
               end select
!               write(*,'(a,3i3,5g15.7)')"i:stME ",i,j,k,totME,C,mgp(i,j,k),sim_h(C,     mgp(i,j,k)),sim_rhoH(sim_h(C,     mgp(i,j,k)))
            end do
         end do
      end do

      call MPI_AllReduce (MPI_IN_PLACE, totME,     1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
      call MPI_AllReduce (MPI_IN_PLACE, totME_aux, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)

      totME     = totME     * dx * dy * dz
      totME_aux = totME_aux * dx * dy * dz

   end subroutine sim_totalMEnthalpic

!-------------------------------------------------------------------------------
! calculate enthalpy

   real function sim_h(C, Phi)

      implicit none

      real, intent(in) :: C, Phi

      sim_h = C - Phi ! rotation will be included here

   end function sim_h

!-------------------------------------------------------------------------------
! find pressure corresponding to given density (EOS dependent)

   real function sim_presrho(rho)

      use initionized,   only : gamma_ion

      implicit none

      real, intent(in) :: rho

      sim_presrho = clump_K *  rho ** gamma_ion

   end function sim_presrho

!-------------------------------------------------------------------------------
! find density corresponding to given enthalpy (EOS dependent)

   real function sim_rhoH(H)

      use initionized, only : gamma_ion
      use mpisetup,    only : smalld

      implicit none

      real, intent(in) :: H

      if (H > 0.) then
         sim_rhoH = ( (1. - 1./gamma_ion) * H / clump_K)**(1./(gamma_ion - 1.))
      else
         sim_rhoH = smalld
      end if

   end function sim_rhoH

!-------------------------------------------------------------------------------
! find enthalpy corresponding to given density (EOS dependent)

   real function sim_hRho(rho)

      use initionized,   only : gamma_ion

      implicit none

      real, intent(in) :: rho

      sim_hRho = clump_K * gamma_ion / (gamma_ion - 1.) * rho ** (gamma_ion - 1.)

   end function sim_hRho

end module initproblem
