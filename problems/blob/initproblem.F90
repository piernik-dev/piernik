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

! Initial condition for blob test
! Blob test by Agertz et al., 2007, MNRAS, 380, 963.
   use problem_pub, only: problem_name, run_id

! ToDo: write support for original, SPH-noisy, initial conditions

   real              :: chi, rblob, blobxc, blobyc, blobzc, Mext, denv, tkh, vgal

   namelist /PROBLEM_CONTROL/  problem_name, run_id, chi, rblob, blobxc, blobyc, blobzc, Mext, denv, tkh, vgal

   contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_public, only: ierrh, msg, par_file
      use errh,          only: namelist_errh
      use func,          only: compare_namelist
      use mpisetup,      only: MPI_CHARACTER, MPI_DOUBLE_PRECISION, cbuff_len, cbuff, rbuff, buffer_dim, comm, ierr, proc

      implicit none


      problem_name = 'aaa'
      run_id  = 'aa'
      chi     = 10.0
      rblob   =  1.0
      blobxc  =  5.0
      blobyc  =  5.0
      blobzc  =  5.0
      Mext    =  2.7
      denv    =  1.0
      tkh     =  1.7
      vgal    =  0.0

      if (proc == 0) then

         diff_nml(PROBLEM_CONTROL)

         cbuff(1) =  problem_name
         cbuff(2) =  run_id

         rbuff(1) = chi
         rbuff(2) = rblob
         rbuff(3) = blobxc
         rbuff(4) = blobyc
         rbuff(5) = blobzc
         rbuff(6) = Mext
         rbuff(7) = denv
         rbuff(8) = tkh
         rbuff(9) = vgal

      endif

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (proc /= 0) then

         problem_name = cbuff(1)
         run_id       = cbuff(2)(1:3)

         chi          = rbuff(1)
         rblob        = rbuff(2)
         blobxc       = rbuff(3)
         blobyc       = rbuff(4)
         blobzc       = rbuff(5)
         Mext         = rbuff(6)
         denv         = rbuff(7)
         tkh          = rbuff(8)
         vgal         = rbuff(9)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use arrays,       only: u
      use grid,         only: x, y, z, nx, ny, nz, nzd, ymin, ymax
      use initneutral,  only: gamma_neu, idnn, imxn, imyn, imzn, ienn

      implicit none

      real    :: penv, rcx, rcy, rrel
      integer :: i, j, k

      penv = 3.2*rblob*sqrt(chi)/tkh/(Mext*gamma_neu/denv)

      u(imzn, :, :, :) = 0.0
      u(ienn, :, :, :) = penv/(gamma_neu-1.0)

      do i = 1,nx
         rcx = (x(i)-blobxc)**2
         do j = 1,ny
            rcy = (y(j)-blobyc)**2
            do k = 1,nz
               if(nzd /= 1) then
                  rrel = sqrt(rcx + rcy + (z(k)-blobzc)**2)
               else
                  rrel = sqrt(rcx + rcy)
               endif

               if (rblob >= rrel) then
                  u(idnn,i,j,k) = chi*denv
                  u(imxn,i,j,k) = chi*denv*vgal
                  u(imyn,i,j,k) = 0.0
               else
                  u(idnn,i,j,k) = denv
                  u(imxn,i,j,k) = denv*vgal
                  u(imyn,i,j,k) = Mext*gamma_neu*penv
               endif
            enddo
         enddo
      enddo

   end subroutine init_prob

!------------------------------------------------------------------------------------------

end module initproblem

