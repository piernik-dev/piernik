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

! Initial condition for blob test
! Blob test by Agertz et al., 2007, MNRAS, 380, 963.
! ToDo: write support for original, SPH-noisy, initial conditions

   implicit none

   private
   public :: read_problem_par, init_prob, problem_pointers

   real   :: chi, rblob, blobxc, blobyc, blobzc, Mext, denv, tkh, vgal

   namelist /PROBLEM_CONTROL/  chi, rblob, blobxc, blobyc, blobzc, Mext, denv, tkh, vgal

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_pub,    only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml, lun   ! QA_WARN required for diff_nml
      use mpi,           only: MPI_DOUBLE_PRECISION
      use mpisetup,      only: rbuff, buffer_dim, comm, mpi_err, master, slave, FIRST

      implicit none

      chi     = 10.0
      rblob   =  1.0
      blobxc  =  5.0
      blobyc  =  5.0
      blobzc  =  5.0
      Mext    =  2.7
      denv    =  1.0
      tkh     =  1.7
      vgal    =  0.0

      if (master) then

         diff_nml(PROBLEM_CONTROL)

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

      call MPI_Bcast(rbuff, buffer_dim, MPI_DOUBLE_PRECISION, FIRST, comm, mpi_err)

      if (slave) then

         chi      = rbuff(1)
         rblob    = rbuff(2)
         blobxc   = rbuff(3)
         blobyc   = rbuff(4)
         blobzc   = rbuff(5)
         Mext     = rbuff(6)
         denv     = rbuff(7)
         tkh      = rbuff(8)
         vgal     = rbuff(9)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use constants,  only: xdim, ydim, zdim
      use domain,     only: dom
      use fluidindex, only: flind
      use fluidtypes, only: component_fluid
      use gc_list,    only: cg_list_element
      use grid,       only: leaves
      use grid_cont,  only: grid_container

      implicit none

      class(component_fluid), pointer :: fl
      real                            :: penv, rcx, rcy, rrel
      integer                         :: i, j, k
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

      fl => flind%neu

      penv = 3.2*rblob*sqrt(chi)/tkh/(Mext*fl%gam/denv)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         cg%u(fl%imz, :, :, :) = 0.0
         cg%u(fl%ien, :, :, :) = penv/fl%gam_1

         do i = 1, cg%n_(xdim)
            rcx = (cg%x(i)-blobxc)**2
            do j = 1, cg%n_(ydim)
               rcy = (cg%y(j)-blobyc)**2
               do k = 1, cg%n_(zdim)
                  if (dom%has_dir(zdim)) then
                     rrel = sqrt(rcx + rcy + (cg%z(k)-blobzc)**2)
                  else
                     rrel = sqrt(rcx + rcy)
                  endif

                  if (rblob >= rrel) then
                     cg%u(fl%idn,i,j,k) = chi*denv
                     cg%u(fl%imx,i,j,k) = chi*denv*vgal
                     cg%u(fl%imy,i,j,k) = 0.0
                  else
                     cg%u(fl%idn,i,j,k) = denv
                     cg%u(fl%imx,i,j,k) = denv*vgal
                     cg%u(fl%imy,i,j,k) = Mext*fl%gam*penv
                  endif
               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine init_prob

!------------------------------------------------------------------------------------------

end module initproblem
