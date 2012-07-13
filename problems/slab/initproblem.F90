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

   implicit none

   private
   public :: read_problem_par, init_prob, problem_pointers

   real               :: d0, r0, bx0, by0, bz0
   integer, parameter :: one = 1
   character(len=one) :: dir

   namelist /PROBLEM_CONTROL/  d0, r0, bx0, by0, bz0

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

      d0      = 1.0
      r0      = 0.25

      if (master) then

         diff_nml(PROBLEM_CONTROL)

         rbuff(1) = d0
         rbuff(2) = r0
         rbuff(3) = bx0
         rbuff(4) = by0
         rbuff(5) = bz0

      endif

      call MPI_Bcast(rbuff, buffer_dim, MPI_DOUBLE_PRECISION, FIRST, comm, mpi_err)

      if (slave) then

         d0       = rbuff(1)
         r0       = rbuff(2)
         bx0      = rbuff(3)
         by0      = rbuff(4)
         bz0      = rbuff(5)

      endif

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine init_prob

      use constants,  only: xdim, ydim, zdim, half
      use fluidindex, only: flind
      use fluidtypes, only: component_fluid
      use func,       only: ekin, emag
      use cg_list,    only: cg_list_element
      use grid,       only: leaves
      use grid_cont,  only: grid_container
#ifndef ISO
      use global,     only: smallei
#endif /* !ISO */
#ifndef FFTW
      use shear,      only: qshear, omega
#endif /* !FFTW */

      implicit none

      class(component_fluid), pointer :: fl
      integer                         :: i, j, k
      real                            :: xi, yj, zk, vx, vy, vz
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
                  vx = 0.0
#ifdef FFTW
                  vy = 0.0
#else /* !FFTW */
                  vy = -qshear*omega*xi
#endif /* !FFTW */
                  vz = 0.0
                  if (abs(yj) <= r0 ) then
                     cg%u(fl%idn,i,j,k) = d0
                  else
                     cg%u(fl%idn,i,j,k) = half*d0
                  endif

                  cg%u(fl%imx,i,j,k) = vx*cg%u(fl%idn,i,j,k)
                  cg%u(fl%imy,i,j,k) = vy*cg%u(fl%idn,i,j,k)
                  cg%u(fl%imz,i,j,k) = vz*cg%u(fl%idn,i,j,k)
#ifndef ISO
                  cg%u(fl%ien,i,j,k) = 1.0/fl%gam_1 !*cg%u(fl%idn,i,j,k)
                  cg%u(fl%ien,i,j,k) = max(cg%u(fl%ien,i,j,k), smallei)
                  cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k))
#endif /* !ISO */
                  cg%b(xdim,i,j,k)   =  bx0
                  cg%b(ydim,i,j,k)   =  by0
                  cg%b(zdim,i,j,k)   =  bz0

#ifndef ISO
                  cg%u(fl%ien,i,j,k)   = cg%u(fl%ien,i,j,k) + emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
#endif /* !ISO */
               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine init_prob
!-----------------------------------------------------------------------------
end module initproblem
