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
!
! Initial condition for 2D Riemann problem with 4 constant states.
! Author: Varadarajan Parthasarathy, CAMK, Warszawa. 15 October 2014.
!
! Problem description
!
! ........................................
! .                 .                    .
! .                 .                    .
! .     *_mp        .         *_pp       .
! .                 .                    .
! .                 .                    .
! ........................................
! .                 .                    .
! .                 .                    .
! .     *_mm        .         *_pm       .
! .                 .                    .
! ........................................
!
! ---------------------------------------------------------------------------------------------------------------------------
! References: 1.  Eleuterio F. Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics-A Practical Introduction
!                 Thrid Edition, Springer.
!
!             2.  Randall J. Leveque, Finite Volume Methods for Hyperbolic Problems
!                 CAMBRIDGE TEXTS IN APPLIED MATHEMATICS, CAMBRIDGE UNIVERSITY PRESS.
! ----------------------------------------------------------------------------------------------------------------------------

   implicit none

   private
   public :: read_problem_par, problem_initial_conditions, problem_pointers

   real :: den_pp, p_pp, velx_pp, vely_pp, den_mp, p_mp, velx_mp, vely_mp, den_mm, p_mm, velx_mm, vely_mm, den_pm, p_pm, velx_pm, vely_pm

   namelist /PROBLEM_CONTROL/ den_pp, p_pp, velx_pp, vely_pp, den_mp, p_mp, velx_mp, vely_mp, den_mm, p_mm, velx_mm, vely_mm, den_pm, p_pm, velx_pm, vely_pm

contains

!----------------------------------------------------------------------------------------------------

   subroutine problem_pointers

      implicit none

   end subroutine problem_pointers

!-----------------------------------------------------------------------------------------------------

   subroutine read_problem_par

      use bcast,      only: piernik_MPI_Bcast
      use dataio_pub, only: nh
      use mpisetup,   only: rbuff, master, slave

      implicit none

      !top right state
      den_pp  = 1.5
      p_pp    = 1.5
      velx_pp = 0.0
      vely_pp = 0.0

      !top left state
      den_mp  = 0.5323
      p_mp    = 0.3
      velx_mp = 1.206
      vely_mp = 0.0

      !bottom left state
      den_mm  = 0.138
      p_mm    = 0.029
      velx_mm = 1.206
      vely_mm = 1.206

      !bottom right state
      den_pm  = 0.5323
      p_pm    = 0.3
      velx_pm = 0.0
      vely_pm = 1.206

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

         rbuff(1) = den_pp
         rbuff(2) = p_pp
         rbuff(3) = velx_pp
         rbuff(4) = vely_pp
         rbuff(5) = den_mp
         rbuff(6) = p_mp
         rbuff(7) = velx_mp
         rbuff(8) = vely_mp
         rbuff(9) = den_mm
         rbuff(10) = p_mm
         rbuff(11) = velx_mm
         rbuff(12) = vely_mm
         rbuff(13) = den_pm
         rbuff(14) = p_pm
         rbuff(15) = velx_pm
         rbuff(16) = vely_pm

      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         den_pp  = rbuff(1)
         p_pp    = rbuff(2)
         velx_pp = rbuff(3)
         vely_pp = rbuff(4)
         den_mp  = rbuff(5)
         p_mp    = rbuff(6)
         velx_mp = rbuff(7)
         vely_mp = rbuff(8)
         den_mm  = rbuff(9)
         p_mm    = rbuff(10)
         velx_mm = rbuff(11)
         vely_mm = rbuff(12)
         den_pm  = rbuff(13)
         p_pm    = rbuff(14)
         velx_pm = rbuff(15)
         vely_pm = rbuff(16)

      endif

   end subroutine read_problem_par

!------------------------------------------------------------------------------------------------------

   subroutine problem_initial_conditions

      use cg_list,    only: cg_list_element
      use cg_leaves,  only: leaves
      use fluidindex, only: flind, iarr_all_mz
      use fluidtypes, only: component_fluid
      use grid_cont,  only: grid_container

      implicit none

      type(cg_list_element), pointer :: cgl
      type(grid_container),  pointer :: cg
      class(component_fluid), pointer :: fl
      integer :: i, j, k, f

      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg

         do k = cg%ks, cg%ke
            do j = cg%js, cg%je
               do i = cg%is, cg%ie

                  do f = 1, flind%fluids
                     fl  => flind%all_fluids(f)%fl
                     if ((cg%x(i) .gt. 0.0) .and. (cg%y(j) .gt. 0.0)) then
                        cg%u(fl%idn, i, j, k) = den_pp
                        cg%u(fl%imx, i, j, k) = den_pp*velx_pp
                        cg%u(fl%imy, i, j, k) = den_pp*vely_pp
                        cg%u(fl%ien, i, j, k) = p_pp / fl%gam_1 + 0.5 / cg%u(fl%idn, i, j, k) * sum(cg%u([fl%imx,fl%imy], i, j, k)**2)
                     else if ((cg%x(i) .lt. 0.0) .and. (cg%y(j) .gt. 0.0)) then
                        cg%u(fl%idn, i, j, k) = den_mp
                        cg%u(fl%imx, i, j, k) = den_mp*velx_mp
                        cg%u(fl%imy, i, j, k) = den_mp*vely_mp
                        cg%u(fl%ien, i, j, k) = p_mp / fl%gam_1 + 0.5 / cg%u(fl%idn, i, j, k) * sum(cg%u([fl%imx,fl%imy], i, j, k)**2)
                     else if ((cg%x(i) .lt. 0.0) .and. (cg%y(j) .lt. 0.0)) then
                        cg%u(fl%idn, i, j, k) = den_mm
                        cg%u(fl%imx, i, j, k) = den_mm*velx_mm
                        cg%u(fl%imy, i, j, k) = den_mm*vely_mm
                        cg%u(fl%ien, i, j, k) = p_mm / fl%gam_1 + 0.5 / cg%u(fl%idn, i, j, k) * sum(cg%u([fl%imx,fl%imy], i, j, k)**2)
                     else if ((cg%x(i) .gt. 0.0) .and. (cg%y(j) .lt. 0.0)) then
                        cg%u(fl%idn, i, j, k) = den_pm
                        cg%u(fl%imx, i, j, k) = den_pm*velx_pm
                        cg%u(fl%imy, i, j, k) = den_pm*vely_pm
                        cg%u(fl%ien, i, j, k) = p_pm / fl%gam_1 + 0.5 / cg%u(fl%idn, i, j, k) * sum(cg%u([fl%imx,fl%imy], i, j, k)**2)
                     endif
                  enddo

               enddo
            enddo
         enddo

         call cg%set_constant_b_field([0., 0., 0.])
         cg%u(iarr_all_mz, :, :, :) = 0.

         cgl => cgl%nxt
      enddo

   end subroutine problem_initial_conditions

end module initproblem
