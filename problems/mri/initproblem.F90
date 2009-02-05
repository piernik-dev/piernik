!$Id$
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

! Initial condition for Parker instability in realistic galactic gravity
! Written by: M. Hanasz, February 2006

   use mpisetup
   use arrays
   use fluidindex,  only : ibx,iby,ibz
   use initionized, only : idni,imxi,imyi,imzi
#ifndef ISO
   use initionized, only : ieni
#endif /* ISO */
   use shear,       only : qshear,omega

   real ::  d0,nbx0,nby0,nbz0,beta,dv
   character problem_name*32,run_id*3

   namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                               d0,nbx0,nby0,nbz0,beta,dv,dp

   contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      implicit none

      problem_name = 'xxx'
      run_id  = 'aaa'
      d0      = 1.0
      nbx0    = 0.0
      nby0    = 1.0
      nbz0    = 0.0
      beta    = 1.0

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
         rbuff(2) = nbx0
         rbuff(3) = nby0
         rbuff(4) = nbz0
         rbuff(5) = beta
         rbuff(6) = dv
         rbuff(7) = dp

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
         nbx0         = rbuff(2)
         nby0         = rbuff(3)
         nbz0         = rbuff(4)
         beta         = rbuff(5)
         dv           = rbuff(6)
         dp           = rbuff(7)
      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob
      use grid,         only : x,nx,ny,nz
      use constants,    only : small
      use initionized,  only : cs_ion, gamma_ion

      implicit none

      integer :: i,j,k
      real :: b0,p0,c_si
      real, dimension(4) :: rand

      c_si = cs_ion

      call read_problem_par

!   Secondary parameters
      p0 = c_si**2*d0
      b0 = sqrt(2.0*p0/beta)

      call random_seed()
      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               u(idni,i,j,k) = d0
               call random_number(rand)
               u(imxi:imzi,i,j,k) = 0.0
               u(imyi,i,j,k) = -qshear*omega*x(i)*u(idni,i,j,k)

               u(imxi,i,j,k) = u(imxi,i,j,k) + dv*c_si*(rand(1)-0.5)/sqrt(gamma_ion)
               u(imyi,i,j,k) = u(imyi,i,j,k) + dv*c_si*(rand(2)-0.5)/sqrt(gamma_ion)
               u(imzi,i,j,k) = u(imzi,i,j,k) + dv*c_si*(rand(3)-0.5)/sqrt(gamma_ion)
#ifndef ISO
               u(ieni,i,j,k)   = p0/(gamma_ion-1.0)*(1.0 + dp*(rand(4)-0.5)) &
                               + 0.5*sum(u(imxi:imzi,i,j,k)**2,1)
#endif /* ISO */
            enddo
         enddo
      enddo


      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               b(ibx,i,j,k)   = b0*sqrt(u(idni,i,j,k)/d0)* nbx0/sqrt(nbx0**2+nby0**2+nbz0**2+small)
               b(iby,i,j,k)   = b0*sqrt(u(idni,i,j,k)/d0)* nby0/sqrt(nbx0**2+nby0**2+nbz0**2+small)
               b(ibz,i,j,k)   = b0*sqrt(u(idni,i,j,k)/d0)* nbz0/sqrt(nbx0**2+nby0**2+nbz0**2+small)
#ifndef ISO
               u(ieni,i,j,k)   = u(ieni,i,j,k) +0.5*sum(b(:,i,j,k)**2,1)
#endif /* ISO */
            enddo
         enddo
      enddo

      return
   end subroutine init_prob

end module initproblem
