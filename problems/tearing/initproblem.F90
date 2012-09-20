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

      use dataio_pub, only: nh      ! QA_WARN required for diff_nml
      use mpisetup,   only: rbuff, master, slave, piernik_MPI_Bcast

      implicit none

      beta        =  1.0
      d0          =  1.0
      v0          =  0.1
      alpha       =  1.0

      if (master) then

         diff_nml(PROBLEM_CONTROL)

         rbuff(1) = beta
         rbuff(2) = v0
         rbuff(3) = d0
         rbuff(4) = alpha

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         beta     = rbuff(1)
         v0       = rbuff(2)
         d0       = rbuff(3)
         alpha    = rbuff(4)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use cg_list,     only: cg_list_element
      use cg_leaves,   only: leaves
      use constants,   only: pi, xdim, ydim, zdim, half
      use domain,      only: dom
      use fluidindex,  only: flind
      use fluidtypes,  only: component_fluid
      use func,        only: ekin, emag
      use grid_cont,   only: grid_container

      implicit none

      class(component_fluid), pointer :: fl
      integer                         :: i, j, k
      real                            :: vzab, b0
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

      fl => flind%ion
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         cg%u(fl%idn,:,:,:) = d0
         cg%u(fl%imy,:,:,:) = 0.0
         cg%u(fl%imz,:,:,:) = 0.0

         cg%b(xdim,:,:,:)   = 0.0
         cg%b(zdim,:,:,:)   = 0.0

         call read_problem_par

         b0 = sqrt(2.*alpha*d0*fl%cs2)

         do k = 1, cg%n_(zdim)
            do j = 1, cg%n_(ydim)
               do i = 1, cg%n_(xdim)

                  vzab = v0*cos(2.*pi*cg%y(j)/dom%L_(ydim))
                  cg%u(fl%imx,i,j,k) = cg%u(fl%idn,i,j,k)*vzab

                  if (abs(cg%x(i)-dom%C_(xdim)) <= 0.25*dom%L_(xdim)) then
                     cg%b(ydim,i,j,k) = -b0
                  else
                     cg%b(ydim,i,j,k) =  b0
                  endif
               enddo
            enddo
         enddo

#ifndef ISO
         cg%u(fl%ien,:,:,:) = half*beta + ekin(cg%u(fl%imx,:,:,:), cg%u(fl%imy,:,:,:), cg%u(fl%imz,:,:,:), cg%u(fl%idn,:,:,:))
         cg%u(fl%ien,:,:,:) = cg%u(fl%ien,:,:,:) + emag(cg%b(xdim,:,:,:), cg%b(ydim,:,:,:), cg%b(zdim,:,:,:))
#endif /* !ISO */

         cgl => cgl%nxt
      enddo

   end subroutine init_prob

end module initproblem
