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

! Initial condition for dust fronts
! Written by: M. Hanasz, January 2009
   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real            :: d0, v0, v1
   integer(kind=4) :: m_x, m_y, m_z

   namelist /PROBLEM_CONTROL/  d0, v0, v1, m_x, m_y, m_z

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_pub, only: nh     ! QA_WARN required for diff_nml
      use mpisetup,   only: ibuff, rbuff, master, slave, piernik_MPI_Bcast

      implicit none

      d0           = 1.0
      v0           = 0.0
      v1           = 0.01
      m_x          = 1
      m_y          = 0
      m_z          = 0

      if (master) then

         diff_nml(PROBLEM_CONTROL)

         rbuff(1) = d0
         rbuff(2) = v0
         rbuff(3) = v1

         ibuff(1) = m_x
         ibuff(2) = m_y
         ibuff(3) = m_z

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         d0           = rbuff(1)
         v0           = rbuff(2)
         v1           = rbuff(3)

         m_x          = ibuff(1)
         m_y          = ibuff(2)
         m_z          = ibuff(3)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use constants,   only: pi, xdim, ydim, zdim, LO, HI
      use domain,      only: dom
      use grid_cont,   only: grid_container
      use fluidindex,  only: flind

      implicit none

      integer :: i,j,k
      real    :: k_x,k_y,k_z,k_a
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      k_x = 2.*pi/dom%L_(xdim)*real(m_x)
      k_y = 2.*pi/dom%L_(ydim)*real(m_y)
      k_z = 2.*pi/dom%L_(zdim)*real(m_z)
      k_a = sqrt(k_x**2+k_y**2+k_z**2)


      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)

                  cg%u(flind%dst%idn,i,j,k) = d0
                  cg%u(flind%dst%imx,i,j,k) = d0*k_x/k_a*(v0 +v1*sin(k_x*cg%x(i)+k_y*cg%y(j)+k_z*cg%z(k)))
                  cg%u(flind%dst%imy,i,j,k) = d0*k_y/k_a*(v0 +v1*sin(k_x*cg%x(i)+k_y*cg%y(j)+k_z*cg%z(k)))
                  cg%u(flind%dst%imz,i,j,k) = d0*k_z/k_a*(v0 +v1*sin(k_x*cg%x(i)+k_y*cg%y(j)+k_z*cg%z(k)))

               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine problem_initial_conditions

end module initproblem
