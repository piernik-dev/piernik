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

! Initial condition for Sedov-Taylor explosion
! Written by: M. Hanasz, March 2006

   implicit none

   private
   public  :: read_problem_par, problem_initial_conditions, problem_pointers

   real :: d0, T0, bx0, by0, bz0, pertamp

   namelist /PROBLEM_CONTROL/ d0, T0, bx0, by0, bz0, pertamp

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use dataio_pub, only: nh
      use mpisetup,   only: ibuff, rbuff, master, slave

      implicit none

      T0      = 1.0e0
      d0      = 1.0
      bx0     = 0.
      by0     = 0.
      bz0     = 0.
      pertamp = 0.

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

! store global parameters to share with other processors with mpi

         rbuff(1) = d0
         rbuff(2) = T0
         rbuff(3) = bx0
         rbuff(4) = by0
         rbuff(5) = bz0
         rbuff(6) = pertamp

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         d0           = rbuff(1)
         T0           = rbuff(2)
         bx0          = rbuff(3)
         by0          = rbuff(4)
         bz0          = rbuff(5)
         pertamp      = rbuff(6)

      endif

   end subroutine read_problem_par
!-----------------------------------------------------------------------------
   subroutine problem_initial_conditions

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use constants,  only: ION, xdim, ydim, zdim, LO, HI, pi
      use domain,     only: dom
      use fluidindex, only: flind
      use grid_cont,  only: grid_container
      use thermal,    only: itemp, thermal_active
      use units,      only: kboltz, mH

      implicit none

      integer                         :: i, j, k, p
      real                            :: cs, p0, kx, ky, kz
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      complex                         :: im

      im = (0,1)
      kx = 2.*pi/dom%L_(xdim)
      ky = 2.*pi/dom%L_(ydim)
      kz = 2.*pi/dom%L_(zdim)
      do p = 1, flind%energ
         associate(fl => flind%all_fluids(p)%fl)

         p0 = d0 * T0 * kboltz / mH
         cs = sqrt(fl%gam * T0 * kboltz / mH)

! Uniform equilibrium state

         cgl => leaves%first
         do while (associated(cgl))
            cg => cgl%cg

! Unperturbed initial state

            do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
               do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
                  do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                     cg%u(fl%idn,i,j,k) = d0
                     cg%u(fl%imx,i,j,k) = 0.0
                     cg%u(fl%imy,i,j,k) = 0.0
                     cg%u(fl%imz,i,j,k) = 0.0
                     cg%u(fl%ien,i,j,k) = p0/(fl%gam_1)
! Perturbation
                     cg%u(fl%imx,i,j,k) = cg%u(fl%imx,i,j,k) + pertamp*cg%u(fl%idn,i,j,k)*cs*sin(2.*pi*cg%x(i)/dom%L_(xdim))*cos(2.*pi*cg%y(j)/dom%L_(ydim))*cos(2.*pi*cg%z(k)/dom%L_(zdim))
                     cg%u(fl%imy,i,j,k) = cg%u(fl%imy,i,j,k) + pertamp*cg%u(fl%idn,i,j,k)*cs*cos(2.*pi*cg%x(i)/dom%L_(xdim))*sin(2.*pi*cg%y(j)/dom%L_(ydim))*cos(2.*pi*cg%z(k)/dom%L_(zdim))
                     cg%u(fl%imz,i,j,k) = cg%u(fl%imz,i,j,k) + pertamp*cg%u(fl%idn,i,j,k)*cs*cos(2.*pi*cg%x(i)/dom%L_(xdim))*cos(2.*pi*cg%y(j)/dom%L_(ydim))*sin(2.*pi*cg%z(k)/dom%L_(zdim))
                     !cg%u(fl%idn,i,j,k) = cg%u(fl%idn,i,j,k) + pertamp*cg%u(fl%idn,i,j,k)*cos(2.*pi*cg%x(i)/dom%L_(xdim))*cos(2.*pi*cg%y(j)/dom%L_(ydim))*sin(2.*pi*cg%z(k)/dom%L_(zdim))

                     !cg%u(fl%idn,i,j,k) = cg%u(fl%idn,i,j,k) + pertamp*cg%u(fl%idn,i,j,k)      * exp(im*(kx*cg%x(i)+ky*cg%y(j)+kz*cg%z(k)))
                     !cg%u(fl%imx,i,j,k) = cg%u(fl%imx,i,j,k) + pertamp*cg%u(fl%idn,i,j,k) * cs * exp(im*(kx*cg%x(i)+ky*cg%y(j)+kz*cg%z(k)))
                     !cg%u(fl%imy,i,j,k) = cg%u(fl%imy,i,j,k) + pertamp*cg%u(fl%idn,i,j,k) * cs * exp(im*(kx*cg%x(i)+ky*cg%y(j)+kz*cg%z(k)))
                     !cg%u(fl%imz,i,j,k) = cg%u(fl%imz,i,j,k) + pertamp*cg%u(fl%idn,i,j,k) * cs * exp(im*(kx*cg%x(i)+ky*cg%y(j)+kz*cg%z(k)))

                     !cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + pertamp*cg%u(fl%ien,i,j,k)      * exp(im*(kx*cg%x(i)+ky*cg%y(j)+kz*cg%z(k)))
                     cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + 0.5*(cg%u(fl%imx,i,j,k)**2 +cg%u(fl%imy,i,j,k)**2 + cg%u(fl%imz,i,j,k)**2)/cg%u(fl%idn,i,j,k)
                  enddo
               enddo
            enddo

            if (thermal_active) cg%q(itemp)%arr(:,:,:) = T0

            if (fl%tag == ION) then
               call cg%set_constant_b_field([bx0, by0, bz0])
               do k = cg%lhn(zdim,LO), cg%lhn(zdim,HI)
                  do j = cg%lhn(ydim,LO), cg%lhn(ydim,HI)
                     do i = cg%lhn(xdim,LO), cg%lhn(xdim,HI)
                        cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + 0.5*(cg%b(xdim,i,j,k)**2 + cg%b(ydim,i,j,k)**2 + cg%b(zdim,i,j,k)**2)
                     enddo
                  enddo
               enddo
            endif

            cgl => cgl%nxt
         enddo

         end associate
      enddo

   end subroutine problem_initial_conditions

end module initproblem
