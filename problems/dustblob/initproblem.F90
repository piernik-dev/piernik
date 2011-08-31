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
   implicit none

   private
   public :: read_problem_par, init_prob, problem_pointers

   real              :: d_gas, p_gas, v_gas, d_dust, v_dust, x0, y0, z0, r0

   namelist /PROBLEM_CONTROL/  &
                               d_gas, p_gas, v_gas, d_dust, v_dust, x0, y0, z0, r0

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par
      use dataio_pub,    only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml, lun, getlun     ! QA_WARN required for diff_nml
      use mpisetup,      only: rbuff, buffer_dim, master, slave, comm, ierr, FIRST
      use mpi,           only: MPI_DOUBLE_PRECISION

      implicit none

      d_gas   =  1.0
      p_gas   =  1.0
      v_gas   =  1.0
      d_dust  =  1.0
      v_dust  =  0.0
      x0      =  0.0
      y0      =  0.0
      z0      =  0.0
      r0      =  1.0

      if (master) then

         diff_nml(PROBLEM_CONTROL)

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

      call MPI_Bcast(rbuff, buffer_dim, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)

      if (slave) then

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

      use constants,   only: xdim, ydim, zdim
      use domain,      only: has_dir
      use fluidindex,  only: flind
      use global,      only: smalld
      use grid,        only: all_cg
      use gc_list,     only: cg_list_element
      use grid_cont,   only: grid_container
      use initdust,    only: idnd, imxd, imyd, imzd
      use initneutral, only: idnn, imxn, imyn, imzn, gamma_neu
#ifndef ISO
      use initneutral, only: ienn
#endif /* !ISO */

      implicit none

      real    :: xi,yj,zk, rc
      integer :: i, j, k
      type(cg_list_element), pointer :: cgl
      type(grid_container), pointer :: cg

      cgl => all_cg%first
      do while (associated(cgl))
         cg => cgl%cg

         if (associated(cg%cs_iso2%arr)) cg%cs_iso2%arr(:,:,:) = flind%neu%cs2

         do i = 1, cg%n_(xdim)
            xi = cg%x(i)
            do j = 1, cg%n_(ydim)
               yj = cg%y(j)
               do k = 1, cg%n_(zdim)
                  if (has_dir(zdim)) then
                     zk = cg%z(k)
                     rc = sqrt((xi-x0)**2+(yj-y0)**2+(zk-z0)**2)
                  else
                     rc = sqrt((xi-x0)**2+(yj-y0)**2)
                  endif

                  cg%u%arr(idnn,i,j,k) = d_gas
                  cg%u%arr(imxn,i,j,k) = 0.0
                  cg%u%arr(imyn,i,j,k) = d_gas*v_gas
                  cg%u%arr(imzn,i,j,k) = 0.0
                  if (rc <= r0) then
                     cg%u%arr(idnd,i,j,k) = d_dust
                     cg%u%arr(imyd,i,j,k) = d_dust*v_dust
                  else
                     cg%u%arr(idnd,i,j,k) = smalld
                     cg%u%arr(imyd,i,j,k) = 0.0
                  endif
                  cg%u%arr(imxd,i,j,k) = 0.0
                  cg%u%arr(imzd,i,j,k) = 0.0
#ifndef ISO
                  cg%u%arr(ienn,i,j,:) = p_gas/(gamma_neu-1.0) + 0.5*(cg%u%arr(imxn,i,j,k)**2 + &
                       &                 cg%u%arr(imyn,i,j,k)**2+cg%u%arr(imzn,i,j,k)**2)/cg%u%arr(idnn,i,j,k)
#endif /* !ISO */
               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine init_prob

!-----------------------------------------------------------------------------------------------------------------------------------

end module initproblem
