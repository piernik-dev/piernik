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

! Initial condition for fluid flows for Kelvin-Helmholtz Instability
! based on Agertz et al. 2008
! Written by: D. Woltanski, February 2008

   implicit none

   private
   public :: read_problem_par, init_prob, problem_pointers
   real   :: chi, dbot, lpert, Mtop, Mbot, dpert, tkh, vtransf

   namelist /PROBLEM_CONTROL/  chi, dbot, lpert, Mtop, Mbot, dpert, tkh, vtransf

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_pub, only: ierrh, par_file, namelist_errh, compare_namelist, cmdl_nml, lun     ! QA_WARN required for diff_nml
      use mpisetup,   only: rbuff, master, slave, piernik_MPI_Bcast

      implicit none

      chi     = 8.0
      dbot    = 1.0
      lpert   = 0.05
      Mtop    = 0.11
      Mbot    = 0.34
      dpert   = 80.0
      tkh     = 1.70
      vtransf = 0.0

      if (master) then

         diff_nml(PROBLEM_CONTROL)

         rbuff(1) = chi
         rbuff(2) = dbot
         rbuff(3) = lpert
         rbuff(4) = Mtop
         rbuff(5) = Mbot
         rbuff(6) = dpert
         rbuff(7) = tkh
         rbuff(8) = vtransf

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         chi      = rbuff(1)
         dbot     = rbuff(2)
         lpert    = rbuff(3)
         Mtop     = rbuff(4)
         Mbot     = rbuff(5)
         dpert    = rbuff(6)
         tkh      = rbuff(7)
         vtransf  = rbuff(8)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine init_prob

      use cg_list,     only: cg_list_element
      use cg_list_bnd, only: leaves
      use constants,   only: dpi, xdim, ydim, zdim
      use domain,      only: dom
      use fluidindex,  only: flind
      use fluidtypes,  only: component_fluid
      use grid_cont,   only: grid_container

      implicit none

      class(component_fluid), pointer :: fl
      real                            :: dtop, lambda, p0, vtop, vbot, k0, vp, rcx, rcy, rc
      integer                         :: i, j
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

      fl => flind%neu

      dtop = dbot/chi
      lambda = 1./6.

      p0 = lambda**2 * (1.+chi)**2 /(chi*tkh**2) *dbot / ( (Mtop*sqrt(chi)+Mbot)**2 * fl%gam)
      vtop  =  1.*Mtop*sqrt(fl%gam*p0/dtop)
      vbot  = -1.*Mbot*sqrt(fl%gam*p0/dbot)
      k0    = dpi/lambda
      vp    = (Mtop*sqrt(chi)+Mbot)*sqrt(fl%gam*p0/dbot)/dpert

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do i = 1, cg%n_(xdim)
            rcx = cg%x(i)
            do j = 1, cg%n_(ydim)
               rcy = cg%y(j)
               rc  = rcy-0.5*dom%L_(ydim)
               if (rc > 0.0) then
                  cg%u(fl%idn,i,j,:) = dtop
                  cg%u(fl%imx,i,j,:) = vtop*dtop
               endif
               if (rc <= 0.0) then
                  cg%u(fl%idn,i,j,:) = dbot
                  cg%u(fl%imx,i,j,:) = vbot*dbot
               endif
               if (abs(rc) < lpert) then
                  cg%u(fl%imy,i,j,:) = vp*sin(k0*rcx)*cg%u(fl%idn,i,j,:)
               else
                  cg%u(fl%imy,i,j,:) = 0.0
               endif
               if (dom%has_dir(zdim)) then
                  cg%u(fl%imz,i,j,:) = vtransf*cg%u(fl%idn,i,j,:)
               else
                  cg%u(fl%imz,i,j,:) = 0.0
               endif
               cg%u(fl%ien,i,j,:) = p0/fl%gam_1
            enddo
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine init_prob

end module initproblem
