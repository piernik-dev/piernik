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
   public :: read_problem_par, init_prob

   real :: d0,r0,bx0,by0,bz0
   integer, parameter :: one = 1
   character(len=one) :: dir

   namelist /PROBLEM_CONTROL/  d0, r0,bx0,by0,bz0

contains

!-----------------------------------------------------------------------------
   subroutine read_problem_par

      use dataio_pub,    only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml      ! QA_WARN required for diff_nml
      use mpisetup,      only: rbuff, buffer_dim, master, slave, comm, ierr
      use mpi,           only: MPI_DOUBLE_PRECISION

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

      call MPI_Bcast(rbuff, buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (slave) then

         d0           = rbuff(1)
         r0           = rbuff(2)
         bx0          = rbuff(3)
         by0          = rbuff(4)
         bz0          = rbuff(5)

      endif

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine init_prob

      use arrays,       only: u,b
      use grid,         only: cg
      use initionized,  only: idni,imxi,imyi,imzi
#ifndef ISO
      use initionized,  only: ieni, gamma_ion
      use mpisetup,     only: smallei
#endif /* !ISO */
      use shear,        only: qshear, omega
      implicit none

      integer  :: i,j,k
      real     :: xi,yj,zk
      real     :: vx,vy,vz

      call read_problem_par

!   Secondary parameters

      do j = 1, cg%ny
         yj = cg%y(j)
         do i = 1, cg%nx
            xi = cg%x(i)
            do k = 1, cg%nz
               zk = cg%z(k)
               vx = 0.0
#ifdef FFTW
               vy = 0.0
#else /* !FFTW */
               vy = -qshear*omega*xi
#endif /* !FFTW */
               vz = 0.0
               if (abs(yj) <= r0 ) then
                  u(idni,i,j,k) = d0
               else
                  u(idni,i,j,k) = 0.5*d0
               endif

               u(imxi,i,j,k) = vx*u(idni,i,j,k)
               u(imyi,i,j,k) = vy*u(idni,i,j,k)
               u(imzi,i,j,k) = vz*u(idni,i,j,k)
#ifndef ISO
               u(ieni,i,j,k) = 1.0/(gamma_ion-1.0)!*u(idni,i,j,k)
               u(ieni,i,j,k) = max(u(ieni,i,j,k), smallei)
               u(ieni,i,j,k) = u(ieni,i,j,k) +0.5*(vx**2+vy**2+vz**2)*u(idni,i,j,k)
#endif /* !ISO */
               b(1,i,j,k)   =  bx0
               b(2,i,j,k)   =  by0
               b(3,i,j,k)   =  bz0

#ifndef ISO
               u(ieni,i,j,k)   = u(ieni,i,j,k) +0.5*sum(b(:,i,j,k)**2,1)
#endif /* !ISO */
            enddo
         enddo
      enddo

   end subroutine init_prob
!-----------------------------------------------------------------------------
end module initproblem
