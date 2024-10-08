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

module initproblem

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real :: dblob, xblob, yblob, zblob, dnamb, rclear, pblob, dnblob, vxfac, p0ambfac, dnambfac, taucool, Tblob, pamb

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

      use bcast,      only: piernik_MPI_Bcast
      use constants,  only: xdim
      use dataio_pub, only: nh
      use domain,     only: dom
      use mpisetup,   only: rbuff, master, slave

      implicit none

      xblob        = -0.5
      yblob        = 0.0
      zblob        = 0.0
      dblob        = 0.0005
      Tblob        = 3500
      dnblob       = 1.0
      dnambfac     = 1.e-5
      p0ambfac     = 1.e-4
      rclear       = dom%L_(xdim)/dom%n_d(xdim) ! BEWARE: why the x-direction is so special here?
      vxfac        = 1.0
      taucool      = 1.0e30

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

         rbuff(1)  = xblob
         rbuff(2)  = yblob
         rbuff(3)  = zblob
         rbuff(4)  = dblob
         rbuff(5)  = Tblob
         rbuff(6)  = dnblob
         rbuff(7)  = dnambfac
         rbuff(8)  = p0ambfac
         rbuff(9)  = rclear
         rbuff(10) = vxfac
         rbuff(11) = taucool

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         xblob     = rbuff(1)
         yblob     = rbuff(2)
         zblob     = rbuff(3)
         dblob     = rbuff(4)
         Tblob     = rbuff(5)
         dnblob    = rbuff(6)
         dnambfac  = rbuff(7)
         p0ambfac  = rbuff(8)
         rclear    = rbuff(9)
         vxfac     = rbuff(10)
         taucool   = rbuff(11)

      endif

      !pressure following the temperature


      pblob = dnblob * Tblob / 7.457e-9
!      pblob = kboltz * dnblob * Tblob / mH / gasRconst
      pamb = p0ambfac * pblob
      dnamb = dnambfac * dnblob

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use constants,   only: xdim, ydim, zdim, LO, HI
      use fluidindex,  only: flind
      use fluidtypes,  only: component_fluid
      use global,      only: smalld
      use grid_cont,   only: grid_container

      implicit none

      class(component_fluid), pointer :: fl
      integer                         :: i, j, k
      real                            :: xi, yj, vx, vy, vz, zk
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

      fl => flind%neu

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
            xi = cg%x(i)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               yj = cg%y(j)
               do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
                  zk = cg%z(k)
!blob
                  cg%u(fl%idn,i,j,k) = dnamb + dnblob*exp(-((xi-xblob)**2+(yj-yblob)**2+zk**2)/dblob)
                  cg%u(fl%idn,i,j,k) = max(cg%u(fl%idn,i,j,k), smalld)
                  cg%u(fl%ien,i,j,k) = pamb/fl%gam_1 + pblob/fl%gam_1*exp(-((xi-xblob)**2+(yj-yblob)**2+zk**2)/dblob)

                  vx = vxfac*sqrt(fl%gam*pblob/dnblob)
                  vy = 0.0
                  vz = 0.0

                  cg%u(fl%imx,i,j,k) = vx*(cg%u(fl%idn,i,j,k)-dnamb)
                  cg%u(fl%imy,i,j,k) = vy*(cg%u(fl%idn,i,j,k)-dnamb)
                  cg%u(fl%imz,i,j,k) = vz*(cg%u(fl%idn,i,j,k)-dnamb)

               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine problem_initial_conditions
!-----------------------------------------------------------------------------

   subroutine impose_inflow(forward)

      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use constants,   only: xdim, ydim, zdim, LO, HI
      use fluidindex,  only: flind
      use fluidtypes,  only: component_fluid
      use global,      only: smalld, smallei, dt
      use grid_cont,   only: grid_container
      use gravity,     only: ptm_x, ptm2_x

      implicit none

      logical, intent(in)             :: forward
      class(component_fluid), pointer :: fl
      integer                         :: i, j, k
      real                            :: xi, yj, zk, r1, r2, dntemp, vx, vy, vz, dnold, enold, entemp, csaim, enaim, coolfac
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg


      fl => flind%neu

!clearing the compact star
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            zk = cg%z(k)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               yj = cg%y(j)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  xi = cg%x(i)


                  r1=sqrt((xi-ptm_x)**2  + yj**2 + zk**2)
                  r2=sqrt((xi-ptm2_x)**2 + yj**2 + zk**2)


                  if (r1<rclear) then

                     cg%u(fl%idn,i,j,k) = smalld + cg%u(fl%idn,i,j,k)*exp((r1/rclear-1.0)/1.0)
!                     cg%u(fl%idn,i,j,k) = dnamb + cg%u(fl%idn,i,j,k)*exp((r1/rclear-1.0)/1.0)
                     cg%u(fl%imx,i,j,k) = cg%u(fl%imx,i,j,k)*exp((r1/rclear-1.0)/1.0)
                     cg%u(fl%imy,i,j,k) = cg%u(fl%imy,i,j,k)*exp((r1/rclear-1.0)/1.0)
                     cg%u(fl%imz,i,j,k) = cg%u(fl%imz,i,j,k)*exp((r1/rclear-1.0)/1.0)
                     cg%u(fl%ien,i,j,k) = smallei + cg%u(fl%ien,i,j,k)*exp((r1/rclear-1.0)/1.0)


                  endif

                  if (r2<rclear) then

                     cg%u(fl%idn,i,j,k) = smalld + cg%u(fl%idn,i,j,k)*exp((r2/rclear-1.0)/1.0)
                     cg%u(fl%imx,i,j,k) = cg%u(fl%imx,i,j,k)*exp((r2/rclear-1.0)/1.0)
                     cg%u(fl%imy,i,j,k) = cg%u(fl%imy,i,j,k)*exp((r2/rclear-1.0)/1.0)
                     cg%u(fl%imz,i,j,k) = cg%u(fl%imz,i,j,k)*exp((r2/rclear-1.0)/1.0)
                     cg%u(fl%ien,i,j,k) = smallei + cg%u(fl%ien,i,j,k)*exp((r2/rclear-1.0)/1.0)

                  endif

               enddo
            enddo
         enddo

!imposing inflow
         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            zk = cg%z(k)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               yj = cg%y(j)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  xi = cg%x(i)
!blob
                  dnold = cg%u(fl%idn,i,j,k)
                  dntemp = dnamb + dnblob*exp(-((xi-xblob)**2+(yj-yblob)**2+zk**2)/dblob)
                  cg%u(fl%idn,i,j,k) = max(dntemp, dnold)
!              cooling
!              ~sound speed ~temp
                  csaim = pamb/dnamb/fl%gam_1
                  enaim = csaim * cg%u(fl%idn,i,j,k)
                  enold = cg%u(fl%ien,i,j,k)
                  coolfac = max((taucool / dt), 1.0)
!               cg%u(fl%ien,i,j,k) = enold - (enold - enaim)/coolfac
!               cg%u(fl%ien,i,j,k) = max(enaim,cg%u(fl%ien,i,j,k))

!hot blob
                  enold = cg%u(fl%ien,i,j,k)
                  entemp = pamb/fl%gam_1 + pblob/fl%gam_1*exp(-((xi-xblob)**2+(yj-yblob)**2+zk**2)/dblob)
                  cg%u(fl%ien,i,j,k) = max(entemp, enold)
!no inflow
!               cg%u(fl%idn,i,j,k) = dnold
!               cg%u(fl%ien,i,j,k) = enold

                  vx = vxfac*sqrt(fl%gam*pblob/dnblob)
                  vy = 0.0
                  vz = 0.0

                  cg%u(fl%imx,i,j,k) = cg%u(fl%imx,i,j,k) + vx*(cg%u(fl%idn,i,j,k)-dnold)
                  cg%u(fl%imy,i,j,k) = cg%u(fl%imy,i,j,k) + vy*(cg%u(fl%idn,i,j,k)-dnold)
                  cg%u(fl%imz,i,j,k) = cg%u(fl%imz,i,j,k) + vz*(cg%u(fl%idn,i,j,k)-dnold)

               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

      return
      if (forward) i = j ! suppress compiler warnings on unused arguments

   end subroutine impose_inflow

end module initproblem
