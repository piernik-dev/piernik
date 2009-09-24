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

   use arrays,       only : u,b
   use grid,         only : x,y,z,nx,ny,nz
   use initionized,  only : idni,imxi,imyi,imzi
#ifndef ISO
   use initionized,  only : ieni, gamma_ion
#endif /* ISO */
   use shear,        only : qshear, omega
   use mpisetup

   real :: d0,r0,bx0,by0,bz0
   character ::  problem_name*32,run_id*3,dir*1

   namelist /PROBLEM_CONTROL/  problem_name, run_id, d0, r0,bx0,by0,bz0


   contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      implicit none

      problem_name = 'slab'
      run_id  = 'tst'
      d0      = 1.0
      r0      = 0.25

      if(proc .eq. 0) then
         open(1,file='problem.par')
         read(unit=1,nml=PROBLEM_CONTROL)
         write(*,nml=PROBLEM_CONTROL)
         close(1)
         open(3, file='tmp.log', position='append')
         write(3,nml=PROBLEM_CONTROL)
         write(3,*)
         close(3)
      endif

      if(proc .eq. 0) then

         cbuff(1) =  problem_name
         cbuff(2) =  run_id

         rbuff(1) = d0
         rbuff(2) = r0
         rbuff(3) = bx0
         rbuff(4) = by0
         rbuff(5) = bz0

         call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
         call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      else

         call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
         call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

         problem_name = cbuff(1)
         run_id       = cbuff(2)

         d0           = rbuff(1)
         r0           = rbuff(2)
         bx0          = rbuff(3)
         by0          = rbuff(4)
         bz0          = rbuff(5)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      implicit none

      integer  :: i,j,k
      real     :: xi,yj,zk
      real     :: vx,vy,vz
      real     :: kn,Lx,kJ,Ly,Lz,Ln

      call read_problem_par

!   Secondary parameters

      do j = 1,ny
         yj = y(j)
         do i = 1,nx
            xi = x(i)
            do k = 1,nz
               zk = z(k)
               vx = 0.0
#ifndef FFTW
               vy = -qshear*omega*xi
#else  /* FFTW */
               vy = 0.0
#endif /* FFTW */
               vz = 0.0
               if(abs(yj) <= r0 ) then
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
#endif /* ISO */
               b(1,i,j,k)   =  bx0
               b(2,i,j,k)   =  by0
               b(3,i,j,k)   =  bz0

#ifndef ISO
               u(ieni,i,j,k)   = u(ieni,i,j,k) +0.5*sum(b(:,i,j,k)**2,1)
#endif /* ISO */
            enddo
         enddo
      enddo
      return
   end subroutine init_prob

end module initproblem
