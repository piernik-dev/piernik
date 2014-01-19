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

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real   :: d0,r0,bx0,by0,bz0

   namelist /PROBLEM_CONTROL/  d0, r0,bx0,by0,bz0

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

      d0      = 1.0
      r0      = 0.25

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
         rbuff(2) = r0

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         d0       = rbuff(1)
         r0       = rbuff(2)

      endif

   end subroutine read_problem_par

!-----------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use constants,   only: pi, dpi, fpi, xdim, ydim, zdim, LO, HI, LEFT
      use fluidindex,  only: flind
      use fluidtypes,  only: component_fluid
      use func,        only: ekin, emag
      use global,      only: smallei
      use grid_cont,   only: grid_container

      implicit none

      class(component_fluid), pointer    :: fl
      integer                            :: i, j, k
      real                               :: xi, yj, zk, vx, vy, vz, rho, pre, bx, by, bz, b0
      real, dimension(:,:,:),allocatable :: A
      type(cg_list_element),  pointer    :: cgl
      type(grid_container),   pointer    :: cg

!   Secondary parameters
      fl => flind%ion

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         if (.not.allocated(A)) allocate(A(cg%lhn(xdim,LO):cg%lhn(xdim,HI), cg%lhn(ydim,LO):cg%lhn(ydim,HI), 1))

         rho = 25.0/(36.0*pi)
         pre =  5.0/(12.0*pi)
         b0  = 1./sqrt(fpi)
         vz  = 0.0
         bz0 = 0.0

         do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
            do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
               A(i,j,1) = b0*(cos(fpi*cg%coord(LEFT, xdim)%r(i))/fpi + cos(dpi*cg%coord(LEFT, ydim)%r(j))/dpi)
            enddo
         enddo

         do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
            yj = cg%y(j)
            do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
               xi = cg%x(i)
               do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
                  zk = cg%z(k)

                  vx  = -sin(dpi*yj)
                  vy  = sin(dpi*xi)
                  bx  = b0*vx
                  by  = b0*sin(fpi*xi)
                  bz  = 0.0

                  cg%u(fl%idn,i,j,k) = rho
                  cg%u(fl%imx,i,j,k) = vx*cg%u(fl%idn,i,j,k)
                  cg%u(fl%imy,i,j,k) = vy*cg%u(fl%idn,i,j,k)
                  cg%u(fl%imz,i,j,k) = vz*cg%u(fl%idn,i,j,k)
#ifndef ISO
                  cg%u(fl%ien,i,j,k) = pre/fl%gam_1
                  cg%u(fl%ien,i,j,k) = max(cg%u(fl%ien,i,j,k), smallei)
                  cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) +ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k))
#endif /* !ISO */
                  cg%b(xdim,i,j,k)  = bx
                  cg%b(ydim,i,j,k)  = by
                  cg%b(zdim,i,j,k)  = bz

#ifndef ISO
                  cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))
#endif /* !ISO */
               enddo
            enddo
         enddo
         if (allocated(A)) deallocate(A)

         cgl => cgl%nxt
      enddo

   end subroutine problem_initial_conditions

end module initproblem
