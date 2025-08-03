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
! ----------------------------------------- !
! Initial condition for a 1D MHD shock tube !
! See: Ryu, Jones, ApJ 442:228-258, (1995)  !
!    and reference therein                  !
! ----------------------------------------- !
! Initial condition for a 1D Sod shcok tube !
! See: Gary A. Sod, J.Comp.Physics 27,1-31(1978) !
! ------------------------------------------!

   use constants, only: I_TWO

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real                   :: dl,vxl,vyl,vzl,bxl,byl,bzl,el
   real                   :: dr,vxr,vyr,vzr,bxr,byr,bzr,er
   character(len=I_TWO)   :: face
   real                   :: theta, c
   namelist /PROBLEM_CONTROL/  dl,vxl,vyl,vzl,bxl,byl,bzl,el,dr,vxr,vyr,vzr,bxr,byr,bzr,er,face, theta ,c

contains

!-----------------------------------------------------------------------------
   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use constants,  only: cbuff_len
      use dataio_pub, only: nh
      use mpisetup,   only: rbuff, cbuff, master, slave

      implicit none

      dl = 1.
      vxl = 0.
      vyl = 0.
      vzl = 0.
      bxl = 0.
      byl = 0.
      bzl = 0.
      el = 0.

      dr = 1.
      vxr = 0.
      vyr = 0.
      vzr = 0.
      bxr = 0.
      byr = 0.
      bzr = 0.
      er = 0.
      c = 0.0
      face = "xy"
      theta = 0.785398163  ! pi/4

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

         rbuff(1)  = dl
         rbuff(2)  = vxl
         rbuff(3)  = vyl
         rbuff(4)  = vzl
         rbuff(5)  = bxl
         rbuff(6)  = byl
         rbuff(7)  = bzl
         rbuff(8)  = el
         rbuff(9)  = dr
         rbuff(10) = vxr
         rbuff(11) = vyr
         rbuff(12) = vzr
         rbuff(13) = bxr
         rbuff(14) = byr
         rbuff(15) = bzr
         rbuff(16) = er
         rbuff(17) = theta
         rbuff(18) = c
         cbuff(1)  = face

      endif

      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(cbuff, cbuff_len)

      if (slave) then

         dl  = rbuff(1)
         vxl = rbuff(2)
         vyl = rbuff(3)
         vzl = rbuff(4)
         bxl = rbuff(5)
         byl = rbuff(6)
         bzl = rbuff(7)
         el  = rbuff(8)
         dr  = rbuff(9)
         vxr = rbuff(10)
         vyr = rbuff(11)
         vzr = rbuff(12)
         bxr = rbuff(13)
         byr = rbuff(14)
         bzr = rbuff(15)
         er  = rbuff(16)
         theta = rbuff(17)
         c = rbuff(18)
         face = trim(cbuff(1))

      endif

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine problem_initial_conditions

      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use constants,   only: xdim, ydim, zdim, LO, HI, pi
      use dataio_pub,  only: die
      use fluidindex,  only: flind
      use fluidtypes,  only: component_fluid
      use func,        only: ekin, emag
      use grid_cont,   only: grid_container
#ifndef ISO
!      use global,      only: smallei
#endif /* !ISO */
      implicit none

      class(component_fluid), pointer :: fl
      integer                         :: i, j, k
      real                            :: xi, yj, zk, vx, vy, vz, rho, pre, bx, by, bz
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      integer                         :: p
      real                            :: m
      real, parameter                 :: tol = 1e-10
      real                            :: theta_mod
      theta_mod = modulo(theta, 2.0*pi)
      p = nint(2.0 * theta_mod / pi)  ! nearest integer to 2*theta/pi
      if (mod(p,2) == 1 .and. abs(theta_mod - p*pi/2.0) < tol) then
         call die("[initproblem:problem_initial_conditions] theta cannot be odd multiple of pi/2")
      endif
      m = tan(theta)

      !   Secondary parameters
      do p = 1, flind%fluids

         fl => flind%all_fluids(p)%fl

         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               yj = cg%y(j)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  xi = cg%x(i)
                  do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
                     zk = cg%z(k)

                     select case (face)
                        case ("xy","yx")
                           if ((-m*xi+yj < c)) then
                              rho = dl
                              pre = el
                              vx  = vxl
                              vy  = vyl
                              vz  = vzl
                              bx  = bxl
                              by  = byl
                              bz  = bzl
                           else
                              rho = dr
                              pre = er
                              vx  = vxr
                              vy  = vyr
                              vz  = vzr
                              bx  = bxr
                              by  = byr
                              bz  = bzr
                           endif
                        case ("yz","zy")
                           if ((-m*yj+zk < c)) then
                              rho = dl
                              pre = el
                              vx  = vxl
                              vy  = vyl
                              vz  = vzl
                              bx  = bxl
                              by  = byl
                              bz  = bzl
                           else
                              rho = dr
                              pre = er
                              vx  = vxr
                              vy  = vyr
                              vz  = vzr
                              bx  = bxr
                              by  = byr
                              bz  = bzr
                           endif
                        case ("zx","xz")
                           if ((-m*xi+zk < c)) then
                              rho = dl
                              pre = el
                              vx  = vxl
                              vy  = vyl
                              vz  = vzl
                              bx  = bxl
                              by  = byl
                              bz  = bzl
                           else
                              rho = dr
                              pre = er
                              vx  = vxr
                              vy  = vyr
                              vz  = vzr
                              bx  = bxr
                              by  = byr
                              bz  = bzr
                           endif
                        case default
                           call die("[initproblem:problem_initial_conditions] wrong problem direction")
                           rho = 0.; pre = 0.; vx = 0.; vy = 0.; vz = 0.; bx = 0.; by = 0.; bz = 0. ! suppress compiler warnings on possibly uninitialized variables
                     end select
                     cg%u(fl%idn,i,j,k) = rho
                     cg%u(fl%imx,i,j,k) = vx*cg%u(fl%idn,i,j,k)
                     cg%u(fl%imy,i,j,k) = vy*cg%u(fl%idn,i,j,k)
                     cg%u(fl%imz,i,j,k) = vz*cg%u(fl%idn,i,j,k)
                     if (fl%has_energy) then

                        cg%u(fl%ien,i,j,k) = pre ! p/(gamma - 1)
                        !cg%u(fl%ien,i,j,k) = max(cg%u(fl%ien,i,j,k, smallei))
                        cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k))

                        if (fl%is_magnetized) then

                           cg%b(xdim,i,j,k)   =  bx
                           cg%b(ydim,i,j,k)   =  by
                           cg%b(zdim,i,j,k)   =  bz
                           cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))

                        endif

                     endif

!#ifndef ISO
!                     cg%u(fl%ien,i,j,k) = pre ! pre here means eint
!                     cg%u(fl%ien,i,j,k) = max(cg%u(fl%ien,i,j,k), smallei)
!                     cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) +ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k))
!#endif /* !ISO */
!                     cg%b(xdim,i,j,k)   =  bx
!                     cg%b(ydim,i,j,k)   =  by
!                     cg%b(zdim,i,j,k)   =  bz

!#ifndef ISO
!                     cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
!#endif /* !ISO */
                  enddo
               enddo
            enddo
            cgl => cgl%nxt
         enddo
      enddo

   end subroutine problem_initial_conditions
!-----------------------------------------------------------------------------
end module initproblem
