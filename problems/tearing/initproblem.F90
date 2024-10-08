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
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real   :: beta, v0, d0, alpha

   namelist /PROBLEM_CONTROL/ beta, v0, d0, alpha

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use dataio_pub, only: nh
      use mpisetup,   only: rbuff, master, slave

      implicit none

      beta        =  1.0
      d0          =  1.0
      v0          =  0.1
      alpha       =  1.0

      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=PROBLEM_CONTROL, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL")
         read(nh%cmdl_nml,nml=PROBLEM_CONTROL, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PROBLEM_CONTROL", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=PROBLEM_CONTROL)
         close(nh%lun)
         call nh%compare_namelist()

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

   subroutine problem_initial_conditions

      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use constants,   only: pi, xdim, ydim, zdim, LO, HI
      use domain,      only: dom
      use fluidindex,  only: flind
      use fluidtypes,  only: component_fluid
      use func,        only: ekin, emag
      use grid_cont,   only: grid_container
#ifndef ISO
      use constants,   only: half
#endif /* !ISO */

      implicit none

      class(component_fluid), pointer :: fl
      integer                         :: i, j
      real                            :: vzab, b0
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

      fl => flind%ion
      b0 = sqrt(2.*alpha*d0*fl%cs2)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         cg%u(fl%idn,:,:,:) = d0
         cg%u(fl%imy,:,:,:) = 0.0
         cg%u(fl%imz,:,:,:) = 0.0

         call cg%set_constant_b_field([0., 0., 0.]) ! slight overkill at ydim for simplicity

         do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
            vzab = v0*cos(2.*pi*cg%y(j)/dom%L_(ydim))
            cg%u(fl%imx,:,j,:) = cg%u(fl%idn,:,j,:)*vzab
         enddo

         do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
            if (abs(cg%x(i)-dom%C_(xdim)) <= 0.25*dom%L_(xdim)) then
               cg%b(ydim,i,:,:) = -b0
            else
               cg%b(ydim,i,:,:) =  b0
            endif
         enddo

#ifndef ISO
         cg%u(fl%ien,:,:,:) = half*beta + ekin(cg%u(fl%imx,:,:,:), cg%u(fl%imy,:,:,:), cg%u(fl%imz,:,:,:), cg%u(fl%idn,:,:,:))
         cg%u(fl%ien,:,:,:) = cg%u(fl%ien,:,:,:) + emag(cg%b(xdim,:,:,:), cg%b(ydim,:,:,:), cg%b(zdim,:,:,:))
#endif /* !ISO */

         cgl => cgl%nxt
      enddo

   end subroutine problem_initial_conditions

end module initproblem
