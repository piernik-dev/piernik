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

!>
!! \brief HLLD Riemann solver for ideal magnetohydrodynamics
!!
!! Varadarajan Parthasarathy, CAMK, Warszawa. 2015.
!! Dr. Artur Gawryszczak, CAMK, Warszawa.
!!
!! RK(N) with N .GE. 3 could be helpful for WENO3 ( this statement to be tested )
!!
!! Reference:Relativistic Hydrodynamics, L. Rezzolla, O. Zanotti
!! ---------------------------------------------------------------------------
!! L (or dtodx)--> discretization of spatial differential operator (Eq. 9.135)
!! ---------------------------------------------------------------------------
!! RK2 (Eq. 9.140)
!! u^(1)   = u^(n) + \Delta t L(u^(n))
!! u^(n+1) = 1/2 ( u^(n) + u^(1) + \Delta t L(u^(1)  )
!! ---------------------------------------------------------------------------
!! RK3 (Eq. 9.141)
!! u^(1)   = u(n) + \Delta t L(u^(n))
!! u^(2)   = 1/4 ( 3 u^(n) + u^(1) + \Delta t L(u^(1) ) )
!! u^(n+1) = 1/3 u^(n) + 2/3 u^(2) + 2/3 \Delta t (u^(2))
!! ---------------------------------------------------------------------------
!!
!! Energy fix up routines for CT and its related comments are not used in the current version.
!! The algorithm is simply present for experimental purposes.
!<

module solvecg_unsplit

! pulled by ANY

   implicit none

   private
   public  :: solve_cg_unsplit

contains

! This routine has to conform to the interface defined in sweeps::sweep

   subroutine solve_cg_unsplit(cg,istep)

      use constants,        only: mag_n, GEO_XYZ
      use dataio_pub,       only: die
      use domain,           only: dom
      use fluidindex,       only: flind
      use grid_cont,        only: grid_container
      use global,           only: use_fargo
      use named_array_list, only: wna
      use sources,          only: prepare_sources

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      integer,                       intent(in) :: istep     ! stage in the time integration scheme
      integer :: nmag, i

      if (dom%geometry_type /= GEO_XYZ) call die("[solve_cg_riemann:solve_cg_riemann] Non-cartesian geometry is not implemented yet in this Riemann solver.")

      call prepare_sources(cg)

      if (wna%exists(mag_n)) then
         nmag = 0
         do i = 1, flind%fluids
            if (flind%all_fluids(i)%fl%is_magnetized) nmag = nmag + 1
         enddo
         if (nmag > 1) call die("[solve_cg_riemann:solve_cg_riemann] At most one magnetized fluid is implemented")
         !call solve_cg_ub(cg, ddim, istep)
      else
         !call solve_cg_u(cg, ddim, istep)
      endif

      cg%processed = .true.

   end subroutine solve_cg_unsplit


end module solvecg_unsplit
