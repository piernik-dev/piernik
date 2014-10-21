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

! Initial condition for 2D Riemann problem with 4 constant states.
! P. Varadarajan, CAMK, Warszawa. 15 October 2014.

  implicit none

  private
  public :: read_problem_par, problem_initial_conditions, problem_pointers

  real :: den_pp, en_pp, velx_pp, vely_pp, den_mp, en_mp, velx_mp, vely_mp, den_mm, en_mm, velx_mm, vely_mm, den_pm, en_pm, velx_pm, vely_pm

  namelist /PROBLEM_CONTROL/ den_pp, en_pp, velx_pp, vely_pp, den_mp, en_mp, velx_mp, vely_mp, den_mm, en_mm, velx_mm, vely_mm, den_pm, en_pm, velx_pm, vely_pm

contains

!----------------------------------------------------------------------------------------------------

  subroutine problem_pointers

    implicit none

  end subroutine problem_pointers

!-----------------------------------------------------------------------------------------------------

  subroutine read_problem_par

    use dataio_pub, only: nh
    use mpisetup, only: rbuff, master, slave, piernik_MPI_Bcast

    implicit none

    den_pp  = 1.5
    en_pp   = 3.75
    velx_pp = 0.0
    vely_pp = 0.0
    den_mp  = 0.5323
    en_mp   = 1.1371
    velx_mp = 1.206
    vely_mp = 0.0
    den_mm  = 0.138
    en_mm   = 0.2732
    velx_mm = 1.206
    vely_mm = 1.206
    den_pm  = 0.5323
    en_pm   = 1.1371
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
         rbuff(2) = en_pp
         rbuff(3) = velx_pp
         rbuff(4) = vely_pp
         rbuff(5) = den_mp
         rbuff(6) = en_mp
         rbuff(7) = velx_mp
         rbuff(8) = vely_mp
         rbuff(9) = den_mm
         rbuff(10) = en_mm
         rbuff(11) = velx_mm
         rbuff(12) = vely_mm
         rbuff(13) = den_pm
         rbuff(14) = en_pm
         rbuff(15) = velx_pm
         rbuff(16) = vely_pm
         
      endif

      call piernik_MPI_Bcast(rbuff)

      if (slave) then 

         den_pp  = rbuff(1)
         en_pp   = rbuff(2)
         velx_pp = rbuff(3)
         vely_pp = rbuff(4)
         den_mp  = rbuff(5)
         en_mp   = rbuff(6)
         velx_mp = rbuff(7)
         vely_mp = rbuff(8)
         den_mm  = rbuff(9)
         en_mm   = rbuff(10)
         velx_mm = rbuff(11)
         vely_mm = rbuff(12)
         den_pm  = rbuff(13)
         en_pm   = rbuff(14)
         velx_pm = rbuff(15)
         vely_pm = rbuff(16)

      endif
      

  
  end subroutine read_problem_par

!------------------------------------------------------------------------------------------------------

  subroutine problem_initial_conditions

    use cg_list,    only: cg_list_element
    use cg_leaves,  only: leaves
    use fluidindex, only: flind
    use domain,     only: dom
    use grid_cont,  only: grid_container
    use constants,  only: xdim, ydim

    implicit none
    
    type(cg_list_element), pointer :: cgl
    type(grid_container),  pointer :: cg

    integer :: i, j, k

    cgl => leaves%first
    do while (associated(cgl))
       cg => cgl%cg

       do k = cg%ks, cg%ke
          do j = cg%js, cg%je
             do i = cg%is, cg%ie
                
                if ((cg%x(i) .gt. 0.0) .and. (cg%y(j) .gt. 0.0)) then
                   cg%u(flind%neu%idn, i, j, k) = den_pp
                   cg%u(flind%neu%ien, i, j, k) = en_pp
                   cg%u(flind%neu%imx, i, j, k) = den_pp*velx_pp
                   cg%u(flind%neu%imy, i, j, k) = den_pp*vely_pp
                else if ((cg%x(i) .lt. 0.0) .and. (cg%y(j) .gt. 0.0)) then
                   cg%u(flind%neu%idn, i, j, k) = den_mp
                   cg%u(flind%neu%ien, i, j, k) = en_mp
                   cg%u(flind%neu%imx, i, j, k) = den_mp*velx_mp
                   cg%u(flind%neu%imy, i, j, k) = den_mp*vely_mp
                else if ((cg%x(i) .lt. 0.0) .and. (cg%y(j) .lt. 0.0)) then
                   cg%u(flind%neu%idn, i, j, k) = den_mm
                   cg%u(flind%neu%ien, i, j, k) = en_mm
                   cg%u(flind%neu%imx, i, j, k) = den_mm*velx_mm
                   cg%u(flind%neu%imy, i, j, k) = den_mm*vely_mm
                else if ((cg%x(i) .gt. 0.0) .and. (cg%y(j) .lt. 0.0)) then
                   cg%u(flind%neu%idn, i, j, k) = den_pm
                   cg%u(flind%neu%ien, i, j, k) = en_pm
                   cg%u(flind%neu%imx, i, j, k) = den_pm*velx_pm
                   cg%u(flind%neu%imy, i, j, k) = den_pm*vely_pm
                
                endif
             
             enddo
          enddo
       enddo

       cgl => cgl%nxt
    enddo

  end subroutine problem_initial_conditions


end module initproblem
