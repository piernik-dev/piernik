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

! Initial condition for tearing instability problem
! Written by: R.K. Pawlaszek July 2007

!       dimdir - sets the direction of magnetic field change
!                choose between: 'x', 'y', 'z'
!
!       magdir - sets the magnetic field component to change
!                choose between: 'x', 'y', 'z'
!
!       dimdir can't be equal magdir!!

   implicit none

   private
   public :: read_problem_par, init_prob, problem_pointers

   real   :: beta, v0, d0, alpha

   namelist /PROBLEM_CONTROL/ beta, v0, d0, alpha

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_pub,    only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml      ! QA_WARN required for diff_nml
      use mpi,           only: MPI_DOUBLE_PRECISION
      use mpisetup,      only: rbuff, buffer_dim, comm, ierr, master, slave

      implicit none

      beta         =  1.0
      d0           =  1.0
      v0           =  0.1
      alpha        =  1.0

      if (master) then

         diff_nml(PROBLEM_CONTROL)

         rbuff(1) = beta
         rbuff(2) = v0
         rbuff(3) = d0
         rbuff(4) = alpha

      endif

      call MPI_Bcast(rbuff, buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      if (slave) then

         beta         = rbuff(1)
         v0           = rbuff(2)
         d0           = rbuff(3)
         alpha        = rbuff(4)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use constants,   only: pi
      use domain,      only: dom
      use fluidindex,  only: ibx, iby, ibz, flind
      use grid,        only: cga
      use grid_cont,   only: cg_list_element, grid_container
      use initionized, only: idni, imxi, imyi, imzi
#ifndef ISO
      use initionized, only: ieni
#endif /* !ISO */

      implicit none

      integer  :: i, j, k
      real     :: xmid, ymid, zmid, vzab, b0
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      xmid = dom%x0
      ymid = dom%y0
      zmid = dom%z0


      cgl => cga%cg_leafs%cg_l(1)
      do while (associated(cgl))
         cg => cgl%cg

         cg%u%arr(idni,:,:,:) = d0
         cg%u%arr(imyi,:,:,:) = 0.0
         cg%u%arr(imzi,:,:,:) = 0.0

         cg%b%arr(ibx,:,:,:)  = 0.0
         cg%b%arr(ibz,:,:,:)  = 0.0

         call read_problem_par

         b0 = sqrt(2.*alpha*d0*flind%ion%cs2)

         do k = 1, cg%nz
            do j = 1, cg%ny
               do i = 1, cg%nx

                  vzab = v0*cos(2.*pi*cg%y(j)/dom%Ly)
                  cg%u%arr(imxi,i,j,k) = cg%u%arr(idni,i,j,k)*vzab

                  if (abs(cg%x(i)-xmid) <= 0.25*dom%Lx) then
                     cg%b%arr(iby,i,j,k) = -b0
                  else
                     cg%b%arr(iby,i,j,k) =  b0
                  endif
               enddo
            enddo
         enddo

#ifndef ISO
         cg%u%arr(ieni,:,:,:) = 0.5*beta + &
              &               0.5*(cg%u%arr(imxi,:,:,:)**2  + cg%u%arr(imyi,:,:,:)**2 + &
              &               cg%u%arr(imzi,:,:,:)**2) / cg%u%arr(idni,:,:,:)

         cg%u%arr(ieni,:,:,:)   = cg%u%arr(ieni,:,:,:) + 0.5*sum(cg%b%arr(:,:,:,:)**2,1)
#endif /* !ISO */

         cgl => cgl%nxt
      enddo

   end subroutine init_prob

end module initproblem
