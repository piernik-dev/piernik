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

! Initial condition for dust fronts
! Written by: M. Hanasz, January 2009
   use problem_pub, only: problem_name, run_id

   real      :: d0, v0, v1
   integer   :: m_x, m_y, m_z

   namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                               d0, v0, v1, m_x, m_y, m_z

   contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_public, only: cwd, msg, par_file
      use func,          only: compare_namelist
      use mpisetup,      only: MPI_CHARACTER, MPI_INTEGER, MPI_DOUBLE_PRECISION, &
                &              cbuff_len, cbuff, ibuff, rbuff, comm, ierr, buffer_dim, proc

      implicit none

      problem_name = 'aaa'
      run_id       = 'aaa'
      d0           = 1.0
      v0           = 0.0
      v1           = 0.01
      m_x          = 1
      m_y          = 0
      m_z          = 0

      if (proc .eq. 0) then

         diff_nml(PROBLEM_CONTROL)


         cbuff(1) =  problem_name
         cbuff(2) =  run_id

         rbuff(1) = d0
         rbuff(2) = v0
         rbuff(3) = v1

         ibuff(1) = m_x
         ibuff(2) = m_y
         ibuff(3) = m_z

      endif

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_Bcast(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (proc /= 0) then

         problem_name = cbuff(1)
         run_id       = cbuff(2)(1:3)

         d0           = rbuff(1)
         v0           = rbuff(2)
         v1           = rbuff(3)

         m_x          = ibuff(1)
         m_y          = ibuff(2)
         m_z          = ibuff(3)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use arrays,    only: u
      use constants, only: pi
      use grid,      only: Lx, Ly, Lz, nx, ny, nz, x, y, z
      use initdust,  only: idnd, imxd, imyd, imzd


      implicit none

      integer :: i,j,k
      real    :: k_x,k_y,k_z,k_a

      write(*,*) 'initproblem:', pi

      k_x = 2.*pi/Lx*real(m_x)
      k_y = 2.*pi/Ly*real(m_y)
      k_z = 2.*pi/Lz*real(m_z)
      k_a = sqrt(k_x**2+k_y**2+k_z**2)

      do i = 1,nx
         do j = 1,ny
            do k = 1,nz

               u(idnd,i,j,k) = d0
               u(imxd,i,j,k) = d0*k_x/k_a*(v0 +v1*sin(k_x*x(i)+k_y*y(j)+k_z*z(k)))
               u(imyd,i,j,k) = d0*k_y/k_a*(v0 +v1*sin(k_x*x(i)+k_y*y(j)+k_z*z(k)))
               u(imzd,i,j,k) = d0*k_z/k_a*(v0 +v1*sin(k_x*x(i)+k_y*y(j)+k_z*z(k)))

            enddo
         enddo
      enddo

      return
   end subroutine init_prob


end module initproblem
