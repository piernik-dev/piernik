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
#include "piernik.def"
#include "macros.h"

module initproblem

! Initial condition for tearing instability problem
! Written by: R.K. Pawlaszek July 2007

!       dimdir - sets the direction of magnetic field change
!                choose between: 'x', 'y', 'z'
!
!       magdir - sets the magnetic field component to change
!                choose between: 'x', 'y', 'z'
!
!       dimdir can't be equal magdir!!

   use problem_pub, only: problem_name, run_id

   real              :: beta, v0

   namelist /PROBLEM_CONTROL/ problem_name, run_id, beta, v0

   contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use mpisetup,      only: MPI_CHARACTER, MPI_DOUBLE_PRECISION, &
           &                   cbuff_len, cbuff, rbuff, buffer_dim, comm, ierr, proc
      use dataio_public, only: ierrh, msg, par_file, namelist_errh, compare_namelist
      use types,         only: idlen

      implicit none


      problem_name = 'tearing'
      run_id       = 'tst'
      beta         =  1.0
      v0           =  0.1

      if (proc == 0) then

         diff_nml(PROBLEM_CONTROL)

         cbuff(1) =  problem_name
         cbuff(2) =  run_id

         rbuff(1) = beta
         rbuff(2) = v0

      endif

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (proc /= 0) then

         problem_name = cbuff(1)
         run_id       = cbuff(2)(1:idlen)

         beta         = rbuff(1)
         v0           = rbuff(2)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob
      use arrays,       only: u, b
      use constants,    only: pi
      use fluidindex,   only: ibx, iby, ibz
      use grid,         only: x, y, nx, ny, nz, xmax, ymax, zmax
      use initionized,  only: idni, imxi, imyi, imzi
#ifndef ISO
      use initionized,  only: ieni
#endif /* !ISO */


      implicit none

      integer  :: i, j, k
      real     :: xmid, ymid, zmid, vzab, b0

      xmid = 0.5*xmax
      ymid = 0.5*ymax
      zmid = 0.5*zmax

      u(idni,:,:,:) = 1.0
      u(imyi,:,:,:) = 0.0
      u(imzi,:,:,:) = 0.0

      b(ibx,:,:,:)  = 0.0
      b(ibz,:,:,:)  = 0.0

      call read_problem_par

      b0 = 1.0

      do k = 1,nz
         do j = 1,ny
            do i = 1,nx

               vzab = v0*dcos(2.*pi*y(j))
               u(imxi,i,j,k) = u(idni,i,j,k)*vzab

               if (abs(x(i)) .LE. xmid) then
                  b(iby,i,j,k) = -b0
               else
                  b(iby,i,j,k) =  b0
               endif
            enddo
         enddo
      enddo

#ifndef ISO
      u(ieni,:,:,:)   = 0.5*beta &
                      + 0.5*(u(imxi,:,:,:)**2  + u(imyi,:,:,:)**2 &
                           + u(imzi,:,:,:)**2) / u(idni,:,:,:)

      u(ieni,:,:,:)   = u(ieni,:,:,:) + 0.5*sum(b(:,:,:,:)**2,1)
#endif /* !ISO */

      return
   end subroutine init_prob

end module initproblem
