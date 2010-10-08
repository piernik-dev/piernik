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
#include "macros.h"

module initproblem

   use problem_pub, only: problem_name, run_id

   real              :: d0,r0,bx0,by0,bz0
   character(len=1)  :: dir

   namelist /PROBLEM_CONTROL/  problem_name, run_id, d0, r0,bx0,by0,bz0

   contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_public, only: ierrh, msg, par_file
      use errh,          only: namelist_errh
      use func,          only: compare_namelist
      use mpisetup,      only: MPI_CHARACTER, MPI_INTEGER, MPI_DOUBLE_PRECISION, &
           &                   cbuff_len, cbuff, ibuff, rbuff, comm, ierr, buffer_dim, proc

      implicit none

      problem_name = 'shock'
      run_id  = 'tst'
      d0      = 1.0
      r0      = 0.25

      if(proc .eq. 0) then

         diff_nml(PROBLEM_CONTROL)

         cbuff(1) =  problem_name
         cbuff(2) =  run_id

         rbuff(1) = d0
         rbuff(2) = r0

      end if

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (proc /= 0) then

         problem_name = cbuff(1)
         run_id       = cbuff(2)(1:3)

         d0           = rbuff(1)
         r0           = rbuff(2)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use arrays,       only: u,b
      use constants,    only: pi, dpi, fpi
      use grid,         only: x, y, z, nx, ny, nz, xl, yl, dx, dy
      use initionized,  only: idni, imxi, imyi, imzi
#ifndef ISO
      use initionized,  only: ieni, gamma_ion
#endif /* !ISO */
      use mpisetup,     only: smallei

      implicit none

      integer :: i, j, k
      real    :: xi, yj, zk
      real    :: vx, vy, vz, rho, pre, bx, by, bz, b0
      real, dimension(:,:,:),allocatable :: A

!   Secondary parameters

      if (.not.allocated(A)) allocate(A(nx,ny,1))

      rho = 25.0/(36.0*pi)
      pre =  5.0/(12.0*pi)
      b0  = 1./sqrt(fpi)
      vz = 0.0
      bz0 = 0.0

      do j=1,ny
         do i = 1,nx
            A(i,j,1) = b0*(dcos(fpi*xl(i))/fpi + dcos(dpi*yl(j))/dpi)
         enddo
      enddo

      do j = 1,ny
         yj = y(j)
         do i = 1,nx
            xi = x(i)
            do k = 1,nz
               zk = z(k)

               vx  = -dsin(dpi*yj)
               vy  = dsin(dpi*xi)
               bx  = b0*vx
               by  = b0*dsin(fpi*xi)
               bz  = 0.0

               u(idni,i,j,k) = rho
               u(imxi,i,j,k) = vx*u(idni,i,j,k)
               u(imyi,i,j,k) = vy*u(idni,i,j,k)
               u(imzi,i,j,k) = vz*u(idni,i,j,k)
#ifndef ISO
               u(ieni,i,j,k) = pre/(gamma_ion-1.0)
               u(ieni,i,j,k) = max(u(ieni,i,j,k), smallei)
               u(ieni,i,j,k) = u(ieni,i,j,k) +0.5*(vx**2+vy**2+vz**2)*u(idni,i,j,k)
#endif /* !ISO */
               b(1,i,j,k)  = bx
               b(2,i,j,k)  = by
               b(3,i,j,k)  = bz

#ifndef ISO
               u(ieni,i,j,k)   = u(ieni,i,j,k) +0.5*sum(b(:,i,j,k)**2,1)
#endif /* !ISO */
            enddo
         enddo
      enddo
      if (allocated(A)) deallocate(A)
      return
   end subroutine init_prob

end module initproblem
