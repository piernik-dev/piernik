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
! ----------------------------------------- !
! Initial condition for a 1D MHD shock tube !
! See: Ryu, Jones, ApJ 442:228-258, (1995)  !
!    and reference therein                  !
! ----------------------------------------- !

   implicit none

   private
   public :: read_problem_par, init_prob, problem_pointers

   real               :: dl,vxl,vyl,vzl,bxl,byl,bzl,el
   real               :: dr,vxr,vyr,vzr,bxr,byr,bzr,er
   integer, parameter :: one = 1
   character(len=one) :: dir

   namelist /PROBLEM_CONTROL/  dl,vxl,vyl,vzl,bxl,byl,bzl,el,dr,vxr,vyr,vzr,bxr,byr,bzr,er

contains

!-----------------------------------------------------------------------------
   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_pub,    only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml, lun      ! QA_WARN required for diff_nml
      use mpi,           only: MPI_DOUBLE_PRECISION
      use mpisetup,      only: rbuff, buffer_dim, master, slave, comm, mpi_err, FIRST

      implicit none

      if (master) then

         diff_nml(PROBLEM_CONTROL)

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

      endif

      call MPI_Bcast(rbuff, buffer_dim, MPI_DOUBLE_PRECISION, FIRST, comm, mpi_err)

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

      endif

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine init_prob

      use constants,  only: xdim, ydim, zdim
      use fluidindex, only: flind
      use fluidtypes, only: component_fluid
      use func,       only: ekin, emag
      use gc_list,    only: cg_list_element
      use grid,       only: leaves
      use grid_cont,  only: grid_container
#ifndef ISO
      use global,     only: smallei
#endif /* !ISO */
      implicit none

      class(component_fluid), pointer :: fl
      integer                         :: i, j, k
      real                            :: xi, yj, zk, vx, vy, vz, rho, pre, bx, by, bz
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

      call read_problem_par

!   Secondary parameters
      fl => flind%ion

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do j = 1, cg%n_(ydim)
            yj = cg%y(j)
            do i = 1, cg%n_(xdim)
               xi = cg%x(i)
               do k = 1, cg%n_(zdim)
                  zk = cg%z(k)

                  if ((xi <= 0.5)) then
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

                  cg%u(fl%idn,i,j,k) = rho
                  cg%u(fl%imx,i,j,k) = vx*cg%u(fl%idn,i,j,k)
                  cg%u(fl%imy,i,j,k) = vy*cg%u(fl%idn,i,j,k)
                  cg%u(fl%imz,i,j,k) = vz*cg%u(fl%idn,i,j,k)
#ifndef ISO
                  cg%u(fl%ien,i,j,k) = pre ! pre here means eint
                  cg%u(fl%ien,i,j,k) = max(cg%u(fl%ien,i,j,k), smallei)
                  cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) +ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k))
#endif /* !ISO */
                  cg%b(xdim,i,j,k)   =  bx
                  cg%b(ydim,i,j,k)   =  by
                  cg%b(zdim,i,j,k)   =  bz

#ifndef ISO
                  cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
#endif /* !ISO */
               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine init_prob
!-----------------------------------------------------------------------------
end module initproblem
