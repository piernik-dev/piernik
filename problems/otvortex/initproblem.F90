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

   real               :: d0,r0,bx0,by0,bz0

   namelist /PROBLEM_CONTROL/  d0, r0,bx0,by0,bz0

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_pub,    only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml, lun, getlun      ! QA_WARN required for diff_nml
      use mpi,           only: MPI_DOUBLE_PRECISION
      use mpisetup,      only: rbuff, comm, ierr, buffer_dim, master, slave, FIRST

      implicit none

      d0      = 1.0
      r0      = 0.25

      if (master) then

         diff_nml(PROBLEM_CONTROL)

         rbuff(1) = d0
         rbuff(2) = r0

      endif

      call MPI_Bcast(rbuff, buffer_dim, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)

      if (slave) then

         d0           = rbuff(1)
         r0           = rbuff(2)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use constants,   only: pi, dpi, fpi, xdim, ydim, zdim
      use global,      only: smallei
      use grid,        only: all_cg
      use gc_list,     only: cg_list_element
      use grid_cont,   only: grid_container
      use initionized, only: idni, imxi, imyi, imzi
#ifndef ISO
      use initionized, only: ieni, gamma_ion
#endif /* !ISO */

      implicit none

      integer :: i, j, k
      real    :: xi, yj, zk
      real    :: vx, vy, vz, rho, pre, bx, by, bz, b0
      real, dimension(:,:,:),allocatable :: A
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

!   Secondary parameters

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         if (.not.allocated(A)) allocate(A(cg%n_(xdim), cg%n_(ydim),1))

         rho = 25.0/(36.0*pi)
         pre =  5.0/(12.0*pi)
         b0  = 1./sqrt(fpi)
         vz = 0.0
         bz0 = 0.0

         do j=1, cg%n_(ydim)
            do i = 1, cg%n_(xdim)
               A(i,j,1) = b0*(cos(fpi*cg%xl(i))/fpi + cos(dpi*cg%yl(j))/dpi)
            enddo
         enddo

         do j = 1, cg%n_(ydim)
            yj = cg%y(j)
            do i = 1, cg%n_(xdim)
               xi = cg%x(i)
               do k = 1, cg%n_(zdim)
                  zk = cg%z(k)

                  vx  = -sin(dpi*yj)
                  vy  = sin(dpi*xi)
                  bx  = b0*vx
                  by  = b0*sin(fpi*xi)
                  bz  = 0.0

                  cg%u(idni,i,j,k) = rho
                  cg%u(imxi,i,j,k) = vx*cg%u(idni,i,j,k)
                  cg%u(imyi,i,j,k) = vy*cg%u(idni,i,j,k)
                  cg%u(imzi,i,j,k) = vz*cg%u(idni,i,j,k)
#ifndef ISO
                  cg%u(ieni,i,j,k) = pre/(gamma_ion-1.0)
                  cg%u(ieni,i,j,k) = max(cg%u(ieni,i,j,k), smallei)
                  cg%u(ieni,i,j,k) = cg%u(ieni,i,j,k) +0.5*(vx**2+vy**2+vz**2)*cg%u(idni,i,j,k)
#endif /* !ISO */
                  cg%b(1,i,j,k)  = bx
                  cg%b(2,i,j,k)  = by
                  cg%b(3,i,j,k)  = bz

#ifndef ISO
                  cg%u(ieni,i,j,k)   = cg%u(ieni,i,j,k) +0.5*sum(cg%b(:,i,j,k)**2,1)
#endif /* !ISO */
               enddo
            enddo
         enddo
         if (allocated(A)) deallocate(A)

         cgl => cgl%nxt
      enddo

   end subroutine init_prob

end module initproblem
