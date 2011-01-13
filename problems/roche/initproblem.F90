! $Id: initproblem.F90 3071 2010-11-05 12:25:40Z xarth $
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

! Initial condition for Keplerian disk
! Written by: M. Hanasz, March 2006

   use mpisetup,    only: cbuff_len

   implicit none

   private
   public :: read_problem_par, init_prob

   real :: dblob, xblob, yblob, zblob, dnamb, rclear, pblob, dnblob, dnin, vxfac, p0ambfac, dnambfac, taucool, Tblob, pamb

   real    :: ptmass1               !< Roche lobe: mass1
   real    :: ptmass2               !< Roche lobe: mass2
   real    :: ptm1_x                !< Roche lobe: location of mass1 (y=z=0)
   real    :: ptm2_x                !< Roche lobe: -//- mass2
   real    :: cmass_x               !< Roche lobe: -//- center of mass
   real    :: Omega                 !< Roche lobe: Rotational angular velocity


   namelist /PROBLEM_CONTROL/  dnblob,xblob,yblob,zblob,rclear,Tblob,dblob, vxfac, p0ambfac, taucool, dnambfac, &
      & ptmass1, ptmass2, ptm1_x, ptm2_x

   contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par
      use dataio_pub,    only: ierrh, par_file, namelist_errh, compare_namelist      ! QA_WARN required for diff_nml
      use mpisetup,      only: rbuff, buffer_dim, proc, comm, ierr
      use mpi,           only: MPI_DOUBLE_PRECISION
      use grid,          only: cg
!      use constants,     only: mH, kboltz, gasRconst
      use types,       only: problem_customize_solution

      implicit none

      xblob            = -0.5
      yblob            = 0.0
      zblob            = 0.0
      dblob            = 0.0005
      Tblob            = 3500
      dnblob           = 1.0
      dnambfac            = 1.e-5
      p0ambfac         = 1.e-4
      rclear           = (cg%x(2)-cg%x(1))
      vxfac            = 1.0
      taucool          = 1.0e30

      ptmass1 = 1.0
      ptmass2 = 1.0
      ptm1_x = 0.0
      ptm2_x = -1.0

      if (proc .eq. 0) then

         diff_nml(PROBLEM_CONTROL)

         rbuff(1) = xblob
         rbuff(2) = yblob
         rbuff(3) = zblob
         rbuff(4) = dblob
         rbuff(5) = Tblob
         rbuff(6) = dnblob
         rbuff(7) = dnambfac
         rbuff(8) = p0ambfac
         rbuff(9) = rclear
         rbuff(10) = vxfac
         rbuff(11) = taucool
         rbuff(12)  = ptmass1
         rbuff(13)  = ptmass2
         rbuff(14)  = ptm1_x
         rbuff(15)  = ptm2_x


      endif

      call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (proc /= 0) then

         xblob            = rbuff(1)
         yblob            = rbuff(2)
         zblob            = rbuff(3)
         dblob            = rbuff(4)
         Tblob            = rbuff(5)
         dnblob           = rbuff(6)
         dnambfac         = rbuff(7)
         p0ambfac         = rbuff(8)
         rclear           = rbuff(9)
         vxfac            = rbuff(10)
         taucool          = rbuff(11)
         ptmass1          = rbuff(12)
         ptmass2          = rbuff(13)
         ptm1_x           = rbuff(14)
         ptm2_x           = rbuff(15)


      endif

      !pressure following the temperature


      pblob = dnblob * Tblob / 7.457e-9
!      pblob = kboltz * dnblob * Tblob / mH / gasRconst
      pamb = p0ambfac * pblob
      dnamb = dnambfac * dnblob

      problem_customize_solution => impose_inflow

   end subroutine read_problem_par

!-----------------------------------------------------------------------------
   subroutine grav_roche
      use arrays,       only: gp
      use grid,         only: cg
      use mpisetup,     only: smalld
      use constants,    only: newtong
      use gravity,      only: r_smooth
      use dataio_pub,   only: printinfo

      implicit none

      integer :: i, j, k
      real    :: r_smooth2, gm, gmr, z2, yz2
      real    :: r1, r2, rc, GM1, GM2
      logical, save :: frun = .true.

      if (.not.frun) return
      call printinfo("[initprob:grav_roche] setting gravity for problem Roche")

       r_smooth2 = r_smooth**2
       GM1 =  newtong * ptmass1
       GM2 =  newtong * ptmass2

!       write(*,*) GM1,GM2,ptmass1,ptmass2,Omega,cmass_x  !QA_WARN

       gm = 0.0 ! to get rid of warnings
       gmr = 0.0
       z2 = 0.0
       yz2 = 0.0

       do k = 1, cg%nz
          do j = 1, cg%ny
             do i = 1, cg%nx
                r1 = (cg%x(i) - ptm1_x)**2  + (cg%y(j) - 0.0)**2 + (cg%z(k) - 0.0)**2
                r2 = (cg%x(i) - ptm2_x)**2  + (cg%y(j) - 0.0)**2 + (cg%z(k) - 0.0)**2
                rc = (cg%x(i) - cmass_x)**2 + (cg%y(j) - 0.0)**2 + (cg%z(k) - 0.0)**2
                    gp(i,j,k) =  - GM1 / dsqrt(r1 + r_smooth2) &
                     - GM2 / dsqrt(r2 + r_smooth2) &
                     - 0.5 * Omega**2 * rc
             enddo
          enddo
       enddo
       frun = .false.
   end subroutine grav_roche
!-----------------------------------------------------------------------------
   subroutine init_prob
      use arrays,      only: u, b, dprof
      use constants,   only: newtong
      use fluidindex,  only: ibx, iby, ibz
      use gravity,     only: r_smooth, r_grav, n_gravr, ptmass, grav_pot_3d, user_grav
      use grid,        only: cg
      use hydrostatic, only: hydrostatic_zeq_densmid
      use initneutral, only: idnn, imxn, imyn, imzn
      use initneutral, only: ienn, gamma_neu
      use mpisetup,    only: smalld, smallei


      implicit none

      integer :: i, j, k, imx,imy,imz,idn
      real    :: xi, yj, vx, vy, vz, zk

      if (user_grav) then
         grav_pot_3d => grav_roche
         call grav_pot_3d
      endif

      cmass_x = (ptmass1*ptm1_x + ptmass2*ptm2_x)/(ptmass1+ptmass2)
      Omega = dsqrt(newtong*(ptmass1+ptmass2)/(abs(ptm1_x-ptm2_x))**3)

      imx=imxn
      imy=imyn
      imz=imzn
      idn=idnn

!      write(*,*) dnamb,dnblob,xblob,dblob  !QA_WARN
      do j = 1,cg%ny
         yj = cg%y(j)
         do i = 1,cg%nx
            xi = cg%x(i)
            do k = 1,cg%nz
               zk = cg%z(k)
!blob
               u(idn,i,j,k) = dnamb + &
                    dnblob*exp(-((xi-xblob)**2+(yj-yblob)**2+zk**2)/dblob)
               u(idn,i,j,k) = max(u(idn,i,j,k), smalld)
               u(ienn,i,j,k) = pamb/(gamma_neu-1.0) + &
                    pblob/(gamma_neu-1.0)*exp(-((xi-xblob)**2+(yj-yblob)**2+zk**2)/dblob)

               vx = vxfac*dsqrt(gamma_neu*pblob/dnblob)
               vy = 0.0
               vz = 0.0

               u(imx,i,j,k) = vx*(u(idn,i,j,k)-dnamb)
               u(imy,i,j,k) = vy*(u(idn,i,j,k)-dnamb)
               u(imz,i,j,k) = vz*(u(idn,i,j,k)-dnamb)


            enddo
         enddo
      enddo


      return
   end subroutine init_prob
!-----------------------------------------------------------------------------

   subroutine impose_inflow

      use dataio_pub,     only: code_progress, PIERNIK_FINISHED, halfstep, msg, die, printinfo
      use arrays,         only: b, u
      use grid,           only: cg
      use initneutral, only: idnn, imxn, imyn, imzn
      use initneutral, only: ienn, gamma_neu
      use mpisetup,    only: smalld, smallei, t, dt

      implicit none

      integer            :: i, j, k, imx,imy,imz,idn,ien
      real    :: xi, yj,zk,r1,r2,dntemp,vx,vy,vz, dnold, enold, entemp, &
           csaim, enaim, coolfac

      imx=imxn
      imy=imyn
      imz=imzn
      idn=idnn
      ien=ienn

!clearing the compact star
      do j = 1,cg%ny
         yj = cg%y(j)
         do i = 1,cg%nx
            xi = cg%x(i)
            do k = 1,cg%nz
               zk = cg%z(k)
               r1=dsqrt((xi-ptm1_x)**2 + yj**2 + zk**2)
               r2=dsqrt((xi-ptm2_x)**2 + yj**2 + zk**2)



               if (r1<rclear) then

                  u(idn,i,j,k) = smalld + u(idn,i,j,k)*exp((r1/rclear-1.0)/1.0)
!                  u(idn,i,j,k) = dnamb + u(idn,i,j,k)*exp((r1/rclear-1.0)/1.0)
                  u(imx,i,j,k) = u(imx,i,j,k)*exp((r1/rclear-1.0)/1.0)
                  u(imy,i,j,k) = u(imy,i,j,k)*exp((r1/rclear-1.0)/1.0)
                  u(imz,i,j,k) = u(imz,i,j,k)*exp((r1/rclear-1.0)/1.0)
                  u(ien,i,j,k) = smallei + u(ien,i,j,k)*exp((r1/rclear-1.0)/1.0)

               endif

               if (r2<rclear) then

                  u(idn,i,j,k) = smalld + u(idn,i,j,k)*exp((r2/rclear-1.0)/1.0)
                  u(imx,i,j,k) = u(imx,i,j,k)*exp((r2/rclear-1.0)/1.0)
                  u(imy,i,j,k) = u(imy,i,j,k)*exp((r2/rclear-1.0)/1.0)
                  u(imz,i,j,k) = u(imz,i,j,k)*exp((r2/rclear-1.0)/1.0)
                  u(ien,i,j,k) = smallei + u(ien,i,j,k)*exp((r2/rclear-1.0)/1.0)

               endif

             enddo
         enddo
      enddo

!imposing inflow
      do j = 1,cg%ny
         yj = cg%y(j)
         do i = 1,cg%nx
            xi = cg%x(i)
            do k = 1,cg%nz
               zk = cg%z(k)
!blob
               dnold = u(idn,i,j,k)
               dntemp = dnamb + dnblob*exp(-((xi-xblob)**2+(yj-yblob)**2+zk**2)/dblob)
               u(idn,i,j,k) = max(dntemp, dnold)

!              cooling
!              ~sound speed ~temp
               csaim = pamb/dnamb/(gamma_neu-1.0)
               enaim = csaim * u(idnn,i,j,k)
               enold = u(ienn,i,j,k)
               coolfac = max((taucool / dt), 1.0)
!               u(ienn,i,j,k) = enold - (enold - enaim)/coolfac

!               u(ienn,i,j,k) = enold - (enold - enaim)/2.0
!               u(ienn,i,j,k) = max(enaim,u(ienn,i,j,k))
!               if ( zk == 0.0 .and. abs(yj)<= 0.01 .and. xi<-.4) then
!               write(*,*) xi,yk,enaim,u(ienn,i,j,k),u(idnn,i,j,k) !QA_WARN
!               endif

!hot blob
               enold = u(ienn,i,j,k)
               entemp = pamb/(gamma_neu-1.0) + &
                     pblob/(gamma_neu-1.0)*exp(-((xi-xblob)**2+(yj-yblob)**2+zk**2)/dblob)
               u(ienn,i,j,k) = max(entemp, enold)

!no inflow
!               u(idnn,i,j,k) = dnold
!               u(ienn,i,j,k) = enold

               vx = vxfac*dsqrt(gamma_neu*pblob/dnblob)
               vy = 0.0
               vz = 0.0

               u(imx,i,j,k) = u(imx,i,j,k) + vx*(u(idn,i,j,k)-dnold)
               u(imy,i,j,k) = u(imy,i,j,k) + vy*(u(idn,i,j,k)-dnold)
               u(imz,i,j,k) = u(imz,i,j,k) + vz*(u(idn,i,j,k)-dnold)

            enddo
         enddo
      enddo


   end subroutine impose_inflow

end module initproblem
