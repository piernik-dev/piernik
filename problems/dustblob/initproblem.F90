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

! Initial condition for blob test
   use problem_pub, only: problem_name, run_id
   implicit none

   private
   public :: read_problem_par, init_prob

   real              :: d_gas, p_gas, v_gas, d_dust, v_dust, x0, y0, z0, r0

   namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                               d_gas, p_gas, v_gas, d_dust, v_dust, x0, y0, z0, r0

   contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par
      use dataio_pub,    only: ierrh, par_file, namelist_errh, compare_namelist     ! QA_WARN
      use mpisetup,      only: cbuff_len, cbuff, rbuff, buffer_dim, proc, comm, ierr, &
                               MPI_CHARACTER, MPI_DOUBLE_PRECISION
      use types,         only: idlen

      implicit none

      problem_name = 'aaa'
      run_id  = 'aaa'
      d_gas   =  1.0
      p_gas   =  1.0
      v_gas   =  1.0
      d_dust  =  1.0
      v_dust  =  0.0
      x0      =  0.0
      y0      =  0.0
      z0      =  0.0
      r0      =  1.0

      if (proc .eq. 0) then

         diff_nml(PROBLEM_CONTROL)

         cbuff(1)     =  problem_name
         cbuff(2)     =  run_id

         rbuff(1)     = d_gas
         rbuff(2)     = p_gas
         rbuff(3)     = v_gas
         rbuff(4)     = d_dust
         rbuff(5)     = v_dust
         rbuff(6)     = x0
         rbuff(7)     = y0
         rbuff(8)     = z0
         rbuff(9)     = r0

      endif

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (proc /= 0) then

         problem_name = cbuff(1)
         run_id       = cbuff(2)(1:idlen)

         d_gas        = rbuff(1)
         p_gas        = rbuff(2)
         v_gas        = rbuff(3)
         d_dust       = rbuff(4)
         v_dust       = rbuff(5)
         x0           = rbuff(6)
         y0           = rbuff(7)
         z0           = rbuff(8)
         r0           = rbuff(9)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob
      use arrays,       only: u
      use grid,         only: x, y, z, nx, ny, nz, nzd
      use initdust,     only: idnd, imxd, imyd, imzd
      use initneutral,  only: idnn, imxn, imyn, imzn, gamma_neu
      use mpisetup,     only: smalld
#ifndef ISO
      use initneutral,  only: ienn
#endif /* !ISO */
      implicit none

      real    :: xi,yj,zk, rc
      integer :: i, j, k

      do i = 1,nx
         xi = x(i)
         do j = 1,ny
            yj = y(j)
            do k = 1,nz
               if (nzd /= 1) then
                  zk = z(k)
                  rc = sqrt((xi-x0)**2+(yj-y0)**2+(zk-z0)**2)
               else
                  rc = sqrt((xi-x0)**2+(yj-y0)**2)
               endif

               u(idnn,i,j,k) = d_gas
               u(imxn,i,j,k) = 0.0
               u(imyn,i,j,k) = d_gas*v_gas
               u(imzn,i,j,k) = 0.0
               if (rc <= r0) then
                  u(idnd,i,j,k) = d_dust
                  u(imyd,i,j,k) = d_dust*v_dust
               else
                  u(idnd,i,j,k) = smalld
                  u(imyd,i,j,k) = 0.0
               endif
               u(imxd,i,j,k) = 0.0
               u(imzd,i,j,k) = 0.0
#ifndef ISO
               u(ienn,i,j,:) = p_gas/(gamma_neu-1.0) + 0.5*(u(imxn,i,j,k)**2 + &
                  u(imyn,i,j,k)**2+u(imzn,i,j,k)**2)/u(idnn,i,j,k)
#endif /* !ISO */
            enddo
         enddo
      enddo

      return
   end subroutine init_prob

!-----------------------------------------------------------------------------------------------------------------------------------

end module initproblem
