! $Id:
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

! Initial condition for Sedov-Taylor explosion
! Written by: M. Hanasz, March 2006

#ifdef IONIZED
   use initionized
#endif /* IONIZED */
#ifdef NEUTRAL
   use initneutral
#endif /* NEUTRAL */

   use arrays
   use grid
   use mpisetup

   real :: t_sn
   integer :: n_sn
   real :: d0,p0,bx0,by0,bz0,Eexpl, x0,y0,z0,r0, dt_sn,r
   character(len=32) :: problem_name
   character(len=3)  :: run_id

   namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                               d0,p0, bx0,by0,bz0, Eexpl,  x0,y0,z0, r0, &
                               n_sn, dt_sn

   contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par
      use errh, only : namelist_errh

      implicit none
      integer :: ierrh

      t_sn = 0.0

      problem_name = 'aaa'
      run_id  = 'aaa'
      d0      = 1.0
      p0      = 1.e-3
      bx0     =   0.
      by0     =   0.
      bz0     =   0.
      Eexpl   = 1.e5
      x0      = 0.0
      y0      = 0.0
      z0      = 0.0
      r0      = dxmn/2.
      n_sn    = 1
      dt_sn   = 0.0

      if(proc .eq. 0) then
         open(1,file='problem.par')
         read(unit=1,nml=PROBLEM_CONTROL,iostat=ierrh)
         call namelist_errh(ierrh,'PROBLEM_CONTROL')
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
         rbuff(2) = p0
         rbuff(3) = bx0
         rbuff(4) = by0
         rbuff(5) = bz0
         rbuff(6) = Eexpl
         rbuff(7) = x0
         rbuff(8) = y0
         rbuff(9) = z0
         rbuff(10)= r0
         rbuff(11)= dt_sn

         ibuff(1) = n_sn

         call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
         call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      else

         call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
         call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

         problem_name = cbuff(1)
         run_id       = cbuff(2)(1:3)

         d0           = rbuff(1)
         p0           = rbuff(2)
         bx0          = rbuff(3)
         by0          = rbuff(4)
         bz0          = rbuff(5)
         Eexpl        = rbuff(6)
         x0           = rbuff(7)
         y0           = rbuff(8)
         z0           = rbuff(9)
         r0           = rbuff(10)
         dt_sn        = rbuff(11)

         n_sn         = ibuff(1)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

#ifdef IONIZED
      use initionized, only : gamma_ion
#endif /* IONIZED */
#ifdef NEUTRAL
      use initneutral, only : gamma_neu
#endif /* NEUTRAL */

      implicit none

      integer :: i,j,k

! Uniform equilibrium state

#ifdef NEUTRAL
      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               u(idnn,i,j,k) = d0
               u(imxn,i,j,k) = 0.0
               u(imyn,i,j,k) = 0.0
               u(imzn,i,j,k) = 0.0
               u(ienn,i,j,k) = p0/(gamma_neu-1.0)
               u(ienn,i,j,k) = u(ienn,i,j,k) + 0.5*(u(imxn,i,j,k)**2 +u(imyn,i,j,k)**2 &
                                                   +u(imzn,i,j,k)**2)/u(idnn,i,j,k)
            enddo
         enddo
      enddo

! Explosions

      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               r = dsqrt( (x(i)-x0)**2 + (y(j)-y0)**2 + (z(k)-z0)**2 )
               if( r**2 < r0**2) then
                  u(ienn,i,j,k)   = u(ienn,i,j,k) + Eexpl
               endif
            enddo
         enddo
      enddo

#endif /* NEUTRAL */

#ifdef IONIZED
      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               u(idni,i,j,k) = d0
               u(imxi,i,j,k) = 0.0
               u(imyi,i,j,k) = 0.0
               u(imzi,i,j,k) = 0.0
               u(ieni,i,j,k) = p0/(gamma_ion-1.0)
               u(ieni,i,j,k) = u(ieni,i,j,k) + 0.5*(u(imxi,i,j,k)**2 +u(imyi,i,j,k)**2 &
                                                   +u(imzi,i,j,k)**2)/u(idni,i,j,k)
               b(1,i,j,k)    = bx0
               b(2,i,j,k)    = by0
               b(3,i,j,k)    = bz0
               u(ieni,i,j,k) = u(ieni,i,j,k) + 0.5*sum(b(:,i,j,k)**2,1)
            enddo
         enddo
      enddo

! Explosions

      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               r = dsqrt( (x(i)-x0)**2 + (y(j)-y0)**2 + (z(k)-z0)**2 )
               if( r**2 < r0**2) then
                  u(ieni,i,j,k)   = u(ieni,i,j,k) + Eexpl
               endif
            enddo
         enddo
      enddo
#endif /* IONIZED */
      return
   end subroutine init_prob

end module initproblem
