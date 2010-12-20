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

! Initial condition for the cosmic ray driven dynamo
! Based on Parker instability setup
! Written by: M. Hanasz, February 2006
! Modified by M.Hanasz for CR-driven dynamo

   implicit none

   private
   public :: read_problem_par, init_prob

   real :: d0, bxn,byn,bzn, alpha

   namelist /PROBLEM_CONTROL/  &
                               d0, bxn,byn,bzn, alpha
   contains

!-----------------------------------------------------------------------------

   subroutine read_problem_par
      use dataio_pub,    only: ierrh, par_file, namelist_errh, compare_namelist   ! QA_WARN required for diff_nml
      use mpisetup,      only: cbuff_len, cbuff, rbuff, buffer_dim, master, slave, comm, ierr
      use mpi,           only: MPI_CHARACTER, MPI_DOUBLE_PRECISION

      implicit none
      d0     = 1.0
      bxn    = 0.0
      byn    = 1.0
      bzn    = 0.0

      if (master) then

         diff_nml(PROBLEM_CONTROL)

         rbuff(1)  = d0
         rbuff(2)  = bxn
         rbuff(3)  = byn
         rbuff(4)  = bzn
         rbuff(5)  = alpha

      endif

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_Bcast(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (slave) then

         d0           = rbuff(1)
         bxn          = rbuff(2)
         byn          = rbuff(3)
         bzn          = rbuff(4)
         alpha        = rbuff(5)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob
      use arrays,      only: u, b, dprof
      use fluidindex,  only: ibx, iby, ibz
      use grid,        only: nx, ny, nz
      use hydrostatic, only: hydrostatic_zeq_densmid
      use initfluids,  only: cs_iso2
      use initionized, only: idni, imxi, imyi, imzi
      use mpisetup,    only: smalld
#ifdef SHEAR
      use shear,       only: qshear, omega
#endif /* SHEAR */

      implicit none

      integer :: i,j,k
      real :: b0, csim2
!   Secondary parameters

      b0 = sqrt(2.*alpha*d0*cs_iso2)

      csim2 = cs_iso2*(1.0+alpha)

      call hydrostatic_zeq_densmid(1, 1, d0, csim2)

      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               u(idni,i,j,k)   = max(smalld,dprof(k))

               u(imxi,i,j,k) = 0.0
               u(imyi,i,j,k) = 0.0
               u(imzi,i,j,k) = 0.0
#ifdef SHEAR
               u(imyi,i,j,k) = -qshear*omega*x(i)*u(idni,i,j,k)
#endif /* SHEAR */

#ifndef ISO
               u(ieni,i,j,k)   = cs_iso2/(gamma_ion-1.0) * u(idni,i,j,k) &
                               + 0.5*(u(imxi,i,j,k)**2 + u(imyi,i,j,k)**2 + &
                                      u(imzi,i,j,k)**2 ) / u(idni,i,j,k)
#endif /* !ISO */

            enddo
         enddo
      enddo

      do k = 1,nz
         do j = 1,ny
            do i = 1,nx
               b(ibx,i,j,k)   = b0*sqrt(u(idni,i,j,k)/d0)* bxn/sqrt(bxn**2+byn**2+bzn**2)
               b(iby,i,j,k)   = b0*sqrt(u(idni,i,j,k)/d0)* byn/sqrt(bxn**2+byn**2+bzn**2)
               b(ibz,i,j,k)   = b0*sqrt(u(idni,i,j,k)/d0)* bzn/sqrt(bxn**2+byn**2+bzn**2)
#ifndef ISO
               u(ieni,i,j,k)   = u(ieni,i,j,k) +0.5*sum(b(:,i,j,k)**2,1)
#endif /* !ISO */
            enddo
         enddo
      enddo

      return
   end subroutine init_prob

end module initproblem
