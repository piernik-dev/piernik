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

   implicit none

   private
   public :: read_problem_par, init_prob, problem_pointers

   real :: dblob, xblob, yblob, zblob, dnamb, rclear, pblob, dnblob, dnin, vxfac, p0ambfac, dnambfac, taucool, Tblob, pamb

   namelist /PROBLEM_CONTROL/  dnblob,xblob,yblob,zblob,rclear,Tblob,dblob, &
                               vxfac, p0ambfac, taucool, dnambfac

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use user_hooks, only: problem_customize_solution

      implicit none

      problem_customize_solution => impose_inflow

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use constants,  only: xdim
      use dataio_pub, only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml, lun, getlun      ! QA_WARN required for diff_nml
      use domain,     only: dom
      use mpi,        only: MPI_DOUBLE_PRECISION
      use mpisetup,   only: rbuff, buffer_dim, proc, comm, ierr, FIRST

      implicit none

      xblob            = -0.5
      yblob            = 0.0
      zblob            = 0.0
      dblob            = 0.0005
      Tblob            = 3500
      dnblob           = 1.0
      dnambfac            = 1.e-5
      p0ambfac         = 1.e-4
      rclear           = dom%L_(xdim)/dom%n_d(xdim) ! BEWARE: why the x-direction is so special here?
      vxfac            = 1.0
      taucool          = 1.0e30

      if (proc == 0) then

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

      endif

      call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)

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

      endif

      !pressure following the temperature


      pblob = dnblob * Tblob / 7.457e-9
!      pblob = kboltz * dnblob * Tblob / mH / gasRconst
      pamb = p0ambfac * pblob
      dnamb = dnambfac * dnblob

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use constants,   only: xdim, ydim, zdim
      use global,      only: smalld
      use grid,        only: all_cg
      use gc_list,     only: cg_list_element
      use grid_cont,   only: grid_container
      use initneutral, only: idnn, imxn, imyn, imzn, ienn, gamma_neu

      implicit none

      integer :: i, j, k, imx,imy,imz,idn
      real    :: xi, yj, vx, vy, vz, zk
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      imx=imxn
      imy=imyn
      imz=imzn
      idn=idnn

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         do j = 1,cg%n_(ydim)
            yj = cg%y(j)
            do i = 1,cg%n_(xdim)
               xi = cg%x(i)
               do k = 1,cg%n_(zdim)
                  zk = cg%z(k)
!blob
                  cg%u(idn,i,j,k) = dnamb + dnblob*exp(-((xi-xblob)**2+(yj-yblob)**2+zk**2)/dblob)
                  cg%u(idn,i,j,k) = max(cg%u(idn,i,j,k), smalld)
                  cg%u(ienn,i,j,k) = pamb/(gamma_neu-1.0) + pblob/(gamma_neu-1.0)*exp(-((xi-xblob)**2+(yj-yblob)**2+zk**2)/dblob)

                  vx = vxfac*sqrt(gamma_neu*pblob/dnblob)
                  vy = 0.0
                  vz = 0.0

                  cg%u(imx,i,j,k) = vx*(cg%u(idn,i,j,k)-dnamb)
                  cg%u(imy,i,j,k) = vy*(cg%u(idn,i,j,k)-dnamb)
                  cg%u(imz,i,j,k) = vz*(cg%u(idn,i,j,k)-dnamb)

               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine init_prob
!-----------------------------------------------------------------------------

   subroutine impose_inflow

      use constants,   only: xdim, ydim, zdim
      use global,      only: smalld, smallei, t, dt
      use grid,        only: all_cg
      use gc_list,     only: cg_list_element
      use grid_cont,   only: grid_container
      use gravity,     only: ptm_x,ptm2_x
      use initneutral, only: idnn, imxn, imyn, imzn, ienn, gamma_neu

      implicit none

      integer :: i, j, k, imx,imy,imz,idn,ien
      real :: xi, yj,zk,r1,r2,dntemp,vx,vy,vz, dnold, enold, entemp, csaim, enaim, coolfac
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg


      imx=imxn
      imy=imyn
      imz=imzn
      idn=idnn
      ien=ienn

!clearing the compact star
      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         do k = 1,cg%n_(zdim)
            zk = cg%z(k)
            do j = 1,cg%n_(ydim)
               yj = cg%y(j)
               do i = 1,cg%n_(xdim)
                  xi = cg%x(i)


                  r1=sqrt((xi-ptm_x)**2  + yj**2 + zk**2)
                  r2=sqrt((xi-ptm2_x)**2 + yj**2 + zk**2)


                  if (r1<rclear) then

                     cg%u(idn,i,j,k) = smalld + cg%u(idn,i,j,k)*exp((r1/rclear-1.0)/1.0)
! cg%u(idn,i,j,k) = dnamb + cg%u(idn,i,j,k)*exp((r1/rclear-1.0)/1.0)
                     cg%u(imx,i,j,k) = cg%u(imx,i,j,k)*exp((r1/rclear-1.0)/1.0)
                     cg%u(imy,i,j,k) = cg%u(imy,i,j,k)*exp((r1/rclear-1.0)/1.0)
                     cg%u(imz,i,j,k) = cg%u(imz,i,j,k)*exp((r1/rclear-1.0)/1.0)
                     cg%u(ien,i,j,k) = smallei + cg%u(ien,i,j,k)*exp((r1/rclear-1.0)/1.0)


                  endif

                  if (r2<rclear) then

                     cg%u(idn,i,j,k) = smalld + cg%u(idn,i,j,k)*exp((r2/rclear-1.0)/1.0)
                     cg%u(imx,i,j,k) = cg%u(imx,i,j,k)*exp((r2/rclear-1.0)/1.0)
                     cg%u(imy,i,j,k) = cg%u(imy,i,j,k)*exp((r2/rclear-1.0)/1.0)
                     cg%u(imz,i,j,k) = cg%u(imz,i,j,k)*exp((r2/rclear-1.0)/1.0)
                     cg%u(ien,i,j,k) = smallei + cg%u(ien,i,j,k)*exp((r2/rclear-1.0)/1.0)

                  endif

               enddo
            enddo
         enddo

!imposing inflow
         do k = 1,cg%n_(zdim)
            zk = cg%z(k)
            do j = 1,cg%n_(ydim)
               yj = cg%y(j)
               do i = 1,cg%n_(xdim)
                  xi = cg%x(i)
!blob
                  dnold = cg%u(idn,i,j,k)
                  dntemp = dnamb + dnblob*exp(-((xi-xblob)**2+(yj-yblob)**2+zk**2)/dblob)
                  cg%u(idn,i,j,k) = max(dntemp, dnold)
!              cooling
!              ~sound speed ~temp
                  csaim = pamb/dnamb/(gamma_neu-1.0)
                  enaim = csaim * cg%u(idnn,i,j,k)
                  enold = cg%u(ienn,i,j,k)
                  coolfac = max((taucool / dt), 1.0)
!               cg%u(ienn,i,j,k) = enold - (enold - enaim)/coolfac
!               cg%u(ienn,i,j,k) = max(enaim,cg%u(ienn,i,j,k))

!hot blob
                  enold = cg%u(ienn,i,j,k)
                  entemp = pamb/(gamma_neu-1.0) + pblob/(gamma_neu-1.0)*exp(-((xi-xblob)**2+(yj-yblob)**2+zk**2)/dblob)
                  cg%u(ienn,i,j,k) = max(entemp, enold)
!no inflow
!               cg%u(idnn,i,j,k) = dnold
!               cg%u(ienn,i,j,k) = enold

                  vx = vxfac*sqrt(gamma_neu*pblob/dnblob)
                  vy = 0.0
                  vz = 0.0

                  cg%u(imx,i,j,k) = cg%u(imx,i,j,k) + vx*(cg%u(idn,i,j,k)-dnold)
                  cg%u(imy,i,j,k) = cg%u(imy,i,j,k) + vy*(cg%u(idn,i,j,k)-dnold)
                  cg%u(imz,i,j,k) = cg%u(imz,i,j,k) + vz*(cg%u(idn,i,j,k)-dnold)

               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine impose_inflow

end module initproblem
