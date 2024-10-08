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

! Initial condition for the cosmic ray driven dynamo
! Based on Parker instability setup
! Written by: M. Hanasz, February 2006
! Modified by M.Hanasz for CR-driven dynamo

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real :: d0, bxn, byn, bzn, alpha

   namelist /PROBLEM_CONTROL/ d0, bxn, byn, bzn, alpha

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
      d0     = 1.0
      bxn    = 0.0
      byn    = 1.0
      bzn    = 0.0
      alpha  = 0.0

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

         rbuff(1) = d0
         rbuff(2) = bxn
         rbuff(3) = byn
         rbuff(4) = bzn
         rbuff(5) = alpha

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         d0       = rbuff(1)
         bxn      = rbuff(2)
         byn      = rbuff(3)
         bzn      = rbuff(4)
         alpha    = rbuff(5)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use constants,   only: xdim, ydim, zdim, LO, HI
      use fluidindex,  only: flind
      use fluidtypes,  only: component_fluid
      use func,        only: ekin, emag
      use global,      only: smalld
      use grid_cont,   only: grid_container
      use hydrostatic, only: hydrostatic_zeq_densmid, set_default_hsparams, dprof
#ifdef SHEAR
      use shear,       only: qshear, omega
#endif /* SHEAR */

      implicit none

      class(component_fluid), pointer :: fl
      integer                         :: i,j,k
      real                            :: b0, csim2
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

!   Secondary parameters
      fl => flind%ion

      b0 = sqrt(2.*alpha*d0*fl%cs2)

      csim2 = fl%cs2*(1.0+alpha)

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         cg%cs_iso2(:,:,:) = fl%cs2

         call set_default_hsparams(cg)
         i = cg%lhn(xdim,LO) ; j = cg%lhn(ydim,LO)
         call hydrostatic_zeq_densmid(i, j, d0, csim2)

         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  cg%u(fl%idn,i,j,k)   = max(smalld, dprof(k))

                  cg%u(fl%imx,i,j,k) = 0.0
                  cg%u(fl%imy,i,j,k) = 0.0
                  cg%u(fl%imz,i,j,k) = 0.0
#ifdef SHEAR
                  cg%u(fl%imy,i,j,k) = -qshear*omega*x(i)*cg%u(fl%idn,i,j,k)
#endif /* SHEAR */

#ifndef ISO
                  cg%u(fl%ien,i,j,k) = fl%cs2/fl%gam_1 * cg%u(fl%idn,i,j,k) + ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k))
#endif /* !ISO */

               enddo
            enddo
         enddo

         do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
            do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
               do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                  cg%b(xdim,i,j,k)   = b0*sqrt(cg%u(fl%idn,i,j,k)/d0)* bxn/sqrt(bxn**2+byn**2+bzn**2)
                  cg%b(ydim,i,j,k)   = b0*sqrt(cg%u(fl%idn,i,j,k)/d0)* byn/sqrt(bxn**2+byn**2+bzn**2)
                  cg%b(zdim,i,j,k)   = b0*sqrt(cg%u(fl%idn,i,j,k)/d0)* bzn/sqrt(bxn**2+byn**2+bzn**2)
#ifndef ISO
                  cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + emag(cg%b(xdim,i,j,k), cg%b(xdim,i,j,k), cg%b(xdim,i,j,k))
#endif /* !ISO */
               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine problem_initial_conditions

end module initproblem
