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
!! \brief Two-dimensional magentic field loop advection test problem
!! Implementation: Dr. Varadarajan Parthasarathy, CAMK (Warszawa)
!!
!! \details The two-dimensional magnetic field loop advection test problem is implemented to PIERNIK following
!! Journal of Computational Physics 229 (2010) 2117–2138, Mignone. A, Tzeferacos. P, Sec. 4.4. Such
!! an implementation can also be found in the PLUTO code.
!!
!! A weak magnetic field loop is advected in a uniform velocity field. As the total pressure is
!! dominated by the thermal component, the magnetic field is transported as a passive scalar.
!! The robustness of the scheme is tested by its ability to preserve the initial circular shape
!! of the loop.
!!
!! The box size is given by (x, y):([-1,1],[-0.5,0.5]) discretized on (N_x,N_x/2) grid cells,
!! where N_x=128. At the initial instant density and pressure are set to 1. The flow velocity
!! is defined as
!!
!!    \vec{V} = (v0*cosalpha,v0*sinalpha)
!!
!! where v0 = \sqrt(5), cosalpha = 2/\sqrt(5) and sinalpha=1/sqrt(5). The magnetic field is defined
!! through the vector potential as
!!
!!    A_z = A0*(R - r) if r <= R
!!        = 0          if r > R
!!
!! where A0 = 10^-3, R = 0.3, and r = \sqrt{x^2 + y^2}.
!!
!! Periodic boundary conditions are used in both the directions.
!<

module initproblem

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real   :: uni_dens, uni_pres, vx, vy, A0, R

   namelist /PROBLEM_CONTROL/ uni_dens, uni_pres, vx, vy, A0, R

contains

!----------------------------------------------------------------------------------------------

   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!----------------------------------------------------------------------------------------------

   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use dataio_pub, only: nh
      use mpisetup,   only: rbuff, master, slave

      implicit none

      uni_dens = 1.0
      uni_pres = 1.0
      vx       = 2.
      vy       = 1.
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
         rbuff(3) = vx
         rbuff(4) = vy
         rbuff(5) = A0
         rbuff(6) = R

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         uni_dens = rbuff(1)
         uni_pres = rbuff(2)
         vx       = rbuff(3)
         vy       = rbuff(4)
         A0       = rbuff(5)
         R        = rbuff(6)

      endif

   end subroutine read_problem_par

!-------------------------------------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use constants,   only: xdim, ydim, zdim, zero, LEFT
      use grid_cont,   only: grid_container
      use fluidindex,  only: flind
      use fluidtypes,  only: component_fluid
      use func,        only: ekin, emag
      use global,      only: cc_mag

      implicit none

      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      class(component_fluid), pointer :: fl

      integer :: i, j, k
      real :: r2

      fl => flind%ion
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         call cg%set_constant_b_field([0., 0., 0.])

         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie

                  ! Density
                  cg%u(fl%idn,i,j,k) = uni_dens
                  ! Velocity
                  cg%u(fl%imx,i,j,k) = vx*cg%u(fl%idn,i,j,k)
                  cg%u(fl%imy,i,j,k) = vy*cg%u(fl%idn,i,j,k)
                  cg%u(fl%imz,i,j,k) = zero
                  ! Mangetic field
                  if (cc_mag) then
                     if ( sqrt(cg%x(i)*cg%x(i) + cg%y(j)*cg%y(j) ) .le. R ) then
                        cg%b(xdim,i,j,k) = -A0*cg%y(j)/(sqrt(cg%x(i)*cg%x(i) + cg%y(j)*cg%y(j) )) !  dA_z/dy
                        cg%b(ydim,i,j,k) =  A0*cg%x(i)/(sqrt(cg%x(i)*cg%x(i) + cg%y(j)*cg%y(j) )) ! -dA_z/dx
                     endif
                  else  ! face-centered components
                     r2 = sum([cg%coord(LEFT, xdim)%r(i), cg%y(j)]**2)
                     if (r2 <= R**2) cg%b(xdim, i, j, k) = -A0 * cg%y(j) / sqrt(r2) !  dA_z/dy
                     r2 = sum([cg%x(i), cg%coord(LEFT, ydim)%r(j)]**2)
                     if (r2 <= R**2) cg%b(ydim, i, j, k) =  A0 * cg%x(i) / sqrt(r2) ! -dA_z/dx
                  endif

                  ! Pressure/Energy
                  cg%u(fl%ien,i,j,k) = uni_pres/fl%gam_1 + ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k)) + &
                       &               emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))

               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine problem_initial_conditions

!-------------------------------------------------------------------------------------------------------------------------------------
end module initproblem
