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

   real    :: d0, T0, bx0, by0, bz0, pertamp

   namelist /PROBLEM_CONTROL/ d0, T0, bx0, by0, bz0, pertamp

contains

!-----------------------------------------------------------------------------

   subroutine problem_pointers

      use dataio_user, only: user_vars_hdf5
      use dataio_user, only: user_tsl

      implicit none

      user_tsl       => thermal_tsl
      user_vars_hdf5 => crtest_analytic_ecr1

   end subroutine problem_pointers

!-----------------------------------------------------------------------------

   subroutine read_problem_par

      use dataio_pub, only: nh
      use mpisetup,   only: ibuff, rbuff, master, slave, piernik_MPI_Bcast

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
         pertamp = rbuff(6)

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
      use units,      only: kboltz, mH

      implicit none

      integer                         :: i, j, k, p
      real                            :: cs, p0
      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg

      do p = 1, flind%energ
         associate(fl => flind%all_fluids(p)%fl)

         p0 = d0/mH*kboltz*T0
         cs = sqrt(fl%gam*T0*kboltz/mH)

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

                     cg%u(fl%ien,i,j,k) = cg%u(fl%ien,i,j,k) + 0.5*(cg%u(fl%imx,i,j,k)**2 +cg%u(fl%imy,i,j,k)**2 + cg%u(fl%imz,i,j,k)**2)/cg%u(fl%idn,i,j,k)
                  enddo
               enddo
            enddo

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

! Add perturbation

            cgl => cgl%nxt
         enddo

         end associate
      enddo

   end subroutine problem_initial_conditions
!-----------------------------------------------------------------------------
   subroutine thermal_tsl(user_vars, tsl_names)

      use constants,   only: pSUM
      use diagnostics, only: pop_vector
      use mpisetup,    only: proc, master, piernik_MPI_Allreduce

      implicit none

      real, dimension(:), intent(inout), allocatable                       :: user_vars
      character(len=*), dimension(:), intent(inout), allocatable, optional :: tsl_names
      real :: output

      if (present(tsl_names)) then
         call pop_vector(tsl_names, len(tsl_names(1)), ["foobar_thermal"])    !   add to header
      else
         ! do mpi stuff here...
         output = real(proc,8)
         call piernik_MPI_Allreduce(output, pSUM)
         if (master) call pop_vector(user_vars,[output])                 !   pop value
      endif

   end subroutine thermal_tsl

!-----------------------------------------------------------------------------

   subroutine crtest_analytic_ecr1(var, tab, ierrh, cg)

      use constants,        only: xdim, ydim, zdim
      use fluidindex,       only: flind
      use fluidtypes,       only: component_fluid
      use func,             only: emag, ekin
      use grid_cont,        only: grid_container
      use named_array_list, only: wna
      use units,            only: kboltz, mH

      implicit none

      character(len=*),               intent(in)             :: var
      real, dimension(:,:,:),         intent(inout)          :: tab
      integer,                        intent(inout)          :: ierrh
      type(grid_container), pointer,  intent(in)             :: cg
      real, dimension(cg%is:cg%ie, cg%js:cg%je, cg%ks:cg%ke) :: eint, kin_ener, mag_ener, temp
      class(component_fluid), pointer                        :: pfl

      !> \warning ONLY ONE FLUID IS USED!!!
      pfl => flind%all_fluids(1)%fl
      if (pfl%has_energy) then
            kin_ener = ekin(cg%w(wna%fi)%span(pfl%imx,cg%ijkse), cg%w(wna%fi)%span(pfl%imy,cg%ijkse), cg%w(wna%fi)%span(pfl%imz,cg%ijkse), cg%w(wna%fi)%span(pfl%idn,cg%ijkse))
            if (pfl%is_magnetized) then
               mag_ener = emag(cg%w(wna%bi)%span(xdim,cg%ijkse), cg%w(wna%bi)%span(ydim,cg%ijkse), cg%w(wna%bi)%span(zdim,cg%ijkse))
               eint = cg%w(wna%fi)%span(pfl%ien,cg%ijkse) - kin_ener - mag_ener
            else
               eint = cg%w(wna%fi)%span(pfl%ien,cg%ijkse) - kin_ener
            endif
      endif

      temp = (pfl%gam-1)*mH/kboltz*eint/cg%w(wna%fi)%span(pfl%idn,cg%ijkse)

      ierrh = 0
      select case (trim(var))
         case ("temp")
            tab(:,:,:) = real(temp, 4)
         case default
            ierrh = -1
      end select

   end subroutine crtest_analytic_ecr1

end module initproblem
