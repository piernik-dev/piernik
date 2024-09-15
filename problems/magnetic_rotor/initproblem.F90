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
!! \brief Two-dimensional magentic rotor problem
!! Author: Dr Varadarajan Parthasarathy (CAMK, Warszawa)
!!
!! \details The two dimensional rotor problem is implemented
!! to PIERNIK following Journal of Computational Physics 229 (2010) 2117–2138, Mignone. A, Tzeferacos. P, Sec. 4.6.
!! Such an implementation can be found in the PLUTO code.
!!
!! This test problem consists of a dense disk rotating in an ambient medium, which is initially static, and threaded
!! by an initially uniform magnetic field. Due to the spin of the rotor, the magnetic field is wrapped around the
!! disk creating torsional Alfven waves. The interaction of the torsional waves with the disk extracts the angular
!! momentum. Also, the buid-up of the magnetic pressure leads to the compression of the rotor.
!!
!! The box size is given by (x,y):([-0.5,0.5],[-0.5,0.5]) discretized on on N_x = N_y = 400 grid cells.
!! At t=0,
!! (\rho,v_{x},v_{y}) = \begin{cases}
!!                        (10,-\omega y, \omega x) & \text{if} r \leq r_{0} \\
!!                        (1 + 9f,-f \omega y \frac{r_{0}}{r}, f \omega x \frac{r_{0}}{r}) & \text{if} r_{0} \le r \le r_{1} \\
!!                        (1,0,0) & \text{if} r \geq r_{1},
!!                      \end{cases}
!! where \omega = 20, r_{0} = 1, r_{1} = 0.115, r = \sqrt{x^{2} + y^{2}}. The taper function is f = (r_{1} - r)/(r_{1} - r_{0}).
!! Thermal pressure is unifrom at t=0, and set to 1 (\gamma = 1.4 is used). The only non-vanishing component of the magnetic
!! field, B_{x} = 5/\sqrt{4\pi}.
!! Outflow boundary conditions are used in both the directions.
!<

module initproblem

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real   :: omega, r0, r1, uni_pres, Bx

   namelist /PROBLEM_CONTROL/ omega, r0, r1, uni_pres, Bx

contains

!-------------------------------------------------------------------------------------------

   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!-------------------------------------------------------------------------------------------

   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use constants,  only: pi
      use dataio_pub, only: nh
      use mpisetup,   only: rbuff, master, slave

      implicit none

      omega    = 20.0
      r0       = 0.1
      r1       = 0.115
      uni_pres = 1.0
      Bx       = 5.0/sqrt(4.0*pi)

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

         rbuff(1) = omega
         rbuff(2) = r0
         rbuff(3) = r1
         rbuff(4) = uni_pres
         rbuff(5) = Bx

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         omega    = rbuff(1)
         r0       = rbuff(2)
         r1       = rbuff(3)
         uni_pres = rbuff(4)
         Bx       = rbuff(5)

      endif

   end subroutine read_problem_par

!-------------------------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use cg_leaves,   only: leaves
      use cg_list,     only: cg_list_element
      use constants,   only: xdim, ydim, zdim, one, zero
      use grid_cont,   only: grid_container
      use fluidindex,  only: flind
      use fluidtypes,  only: component_fluid
      use func,        only: ekin, emag
!    use global,      only: cc_mag

      implicit none

      type(cg_list_element),  pointer :: cgl
      type(grid_container),   pointer :: cg
      class(component_fluid), pointer :: fl

      integer                         :: i,j,k
      real                            :: f_taper, r0_r

      fl  => flind%ion
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         call cg%set_constant_b_field([0., 0., 0.])

         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie

                  ! r = sqrt( cg%x(i)*cg%x(i)+cg%y(j)*cg%y(j) )

                  f_taper = ( r1 - sqrt( cg%x(i)*cg%x(i)+cg%y(j)*cg%y(j) ) )/( r1 - r0 )
                  r0_r    = r0/sqrt( cg%x(i)*cg%x(i)+cg%y(j)*cg%y(j) )

                  if ( sqrt( cg%x(i)*cg%x(i)+cg%y(j)*cg%y(j) ) <= r0 ) then ! r <= r0

                     ! Density
                     cg%u(fl%idn,i,j,k) = 10.0
                     ! Velocity
                     cg%u(fl%imx,i,j,k) = -omega*cg%y(j) * cg%u(fl%idn,i,j,k)
                     cg%u(fl%imy,i,j,k) =  omega*cg%x(i) * cg%u(fl%idn,i,j,k)
                     cg%u(fl%imz,i,j,k) = zero

                  else if ( sqrt( cg%x(i)*cg%x(i)+cg%y(j)*cg%y(j) ) < r1 ) then ! r< r1

                     ! Density
                     cg%u(fl%idn,i,j,k) = one + 9.0*f_taper
                     ! Velocity
                     cg%u(fl%imx,i,j,k) = -f_taper*omega*cg%y(j)*r0_r * cg%u(fl%idn,i,j,k)
                     cg%u(fl%imy,i,j,k) =  f_taper*omega*cg%x(i)*r0_r * cg%u(fl%idn,i,j,k)
                     cg%u(fl%imz,i,j,k) = zero

                  else ! r > r1

                     ! Density
                     cg%u(fl%idn,i,j,k) = 1.0
                     ! Velocity
                     cg%u(fl%imx,i,j,k) = zero
                     cg%u(fl%imy,i,j,k) = zero
                     cg%u(fl%imz,i,j,k) = zero

                  endif

                  ! Pressure/Energy
                  cg%u(fl%ien,i,j,k) = uni_pres/fl%gam_1 + ekin(cg%u(fl%imx,i,j,k), cg%u(fl%imy,i,j,k), cg%u(fl%imz,i,j,k), cg%u(fl%idn,i,j,k)) + &
                       &               emag(cg%b(xdim,i,j,k), cg%b(ydim,i,j,k), cg%b(zdim,i,j,k))

                  !if (cc_mag) then

                  cg%b(xdim,i,j,k) = Bx
                  cg%b(ydim,i,j,k) = zero
                  cg%b(zdim,i,j,k) = zero

                  !endif
                  ! For CT, with face-centered A_{z} = Bx*y

               enddo
            enddo
         enddo

         cgl => cgl%nxt
      enddo

   end subroutine problem_initial_conditions

end module initproblem
