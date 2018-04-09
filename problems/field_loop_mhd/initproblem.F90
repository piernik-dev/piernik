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
!    Two-dimensional magentic field loop advection test problem
!    Implementation: Dr. Varadarajan Parthasarathy, CAMK (Warszawa)
!    Read info for details of the test problem.
!
#include "piernik.h"

module initproblem

  implicit none

  private
  public :: read_problem_par, problem_initial_conditions, problem_pointers

  real   :: uni_dens, uni_pres, v0, sinalpha, cosalpha, A0, R

  namelist /PROBLEM_CONTROL/ uni_dens, uni_pres, v0, sinalpha, cosalpha, A0, R

contains

!----------------------------------------------------------------------------------------------

  subroutine problem_pointers

    implicit none

  end subroutine problem_pointers

!----------------------------------------------------------------------------------------------

  subroutine read_problem_par

    use dataio_pub, only: nh
    use mpisetup, only: rbuff, master, slave, PIERNIK_MPI_Bcast

    implicit none

    uni_dens = 1.0
    uni_pres = 1.0
    v0       = sqrt(5.0)
    sinalpha = 1.0/sqrt(5.0)
    cosalpha = 2.0/sqrt(5.0)
    A0       = 1.e-3
    R        = 0.3

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

         rbuff(1) = uni_dens
         rbuff(2) = uni_pres
         rbuff(3) = v0
         rbuff(4) = sinalpha
         rbuff(5) = cosalpha
         rbuff(6) = A0
         rbuff(7) = R

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         uni_dens = rbuff(1)
         uni_pres = rbuff(2)
         v0       = rbuff(3)
         sinalpha = rbuff(4)
         cosalpha = rbuff(5)
         A0       = rbuff(6)
         R        = rbuff(7)

      endif

  end subroutine read_problem_par

!-------------------------------------------------------------------------------------------------------

  subroutine problem_initial_conditions

    use cg_leaves,   only: leaves
    use cg_list,     only: cg_list_element
    use constants,   only: xdim, ydim, zdim
    use grid_cont,   only: grid_container
    use fluidindex,  only: flind
    use fluidtypes,  only: component_fluid
    use func,        only: ekin, emag
    use constants,   only: zero

    implicit none

    type(cg_list_element),  pointer :: cgl
    type(grid_container),   pointer :: cg
    class(component_fluid), pointer :: fl

    integer :: i, j, k

    fl => flind%ion
    cgl => leaves%first
    do while (associated(cgl))
       cg => cgl%cg

       do k = cg%ks, cg%ke
          do j = cg%js, cg%je
             do i = cg%is, cg%ie

                ! Density
                cg%u(fl%idn,i,j,k) = uni_dens
                ! Velocity
                cg%u(fl%imx,i,j,k) = v0*cosalpha*cg%u(fl%idn,i,j,k)
                cg%u(fl%imy,i,j,k) = v0*sinalpha*cg%u(fl%idn,i,j,k)
                cg%u(fl%imz,i,j,k) = zero
                ! Mangetic field
                if ( sqrt(cg%x(i)*cg%x(i) + cg%y(j)*cg%y(j) ) .le. R ) then

                   cg%b(xdim,i,j,k) = -A0*cg%y(j)/(sqrt(cg%x(i)*cg%x(i) + cg%y(j)*cg%y(j) )) !  dA_z/dy
                   cg%b(ydim,i,j,k) =  A0*cg%x(i)/(sqrt(cg%x(i)*cg%x(i) + cg%y(j)*cg%y(j) )) ! -dA_z/dx
                else
                   cg%b(xdim,i,j,k) = zero
                   cg%b(ydim,i,j,k) = zero
                endif
                cg%b(zdim,i,j,k) = zero
                ! Pressure/Energy
                cg%u(fl%ien,i,j,k) = uni_pres/fl%gam_1 + ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k)) + &
                                             emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))

             enddo
          enddo
       enddo

       cgl => cgl%nxt
    enddo

  end subroutine problem_initial_conditions

!-------------------------------------------------------------------------------------------------------------------------------------
end module initproblem
