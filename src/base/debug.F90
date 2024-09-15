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
!! \brief Module providing few parameters for debugging and developing the code.
!!
!! \details Keep this module as free of dependences as possible to allow its use early after start and eliminate risk of cyclic dependences.
!<
module piernikdebug
! pulled by DEBUG

   implicit none

   private
   public :: init_piernikdebug, has_const_dt, constant_dt

   real,    protected :: constant_dt               !< value of timestep regardless of fluid state
   logical, protected :: has_const_dt              !< true if piernikdebug::constant_dt > 0

   namelist /PIERNIK_DEBUG/ constant_dt

contains

!>
!! \brief Routine to set debug auxiliary parameters
!!
!! \n \n
!! @b PIERNIK_DEBUG
!! \n \n
!! <table border="+1">
!!   <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!!   <tr><td>constant_dt      </td><td>0.0    </td><td>real value             </td><td>\copydoc piernikdebug::constant_dt      </td></tr>
!! </table>
!! The list is active while \b "DEBUG" is defined.
!! \n \n
!<

   subroutine init_piernikdebug

      use bcast,      only: piernik_MPI_Bcast
      use constants,  only: PIERNIK_INIT_MPI  !, cbuff_len
      use dataio_pub, only: code_progress, die, nh
      use mpisetup,   only: master, slave, rbuff  ! lbuff, cbuff, ibuff, buffer_dim

      implicit none

      if (code_progress < PIERNIK_INIT_MPI) call die("[debug:init_piernikdebug] MPI not initialized.")

      constant_dt = 0.0

      if (master) then
         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=PIERNIK_DEBUG)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=PIERNIK_DEBUG, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PIERNIK_DEBUG")
         read(nh%cmdl_nml,nml=PIERNIK_DEBUG, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PIERNIK_DEBUG", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=PIERNIK_DEBUG)
         close(nh%lun)
         call nh%compare_namelist()

         rbuff(1) = constant_dt

      endif

!      call piernik_MPI_Bcast(cbuff, cbuff_len)
!      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)
!      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         constant_dt = rbuff(1)

      endif

      has_const_dt = (constant_dt > 0.0)

   end subroutine init_piernikdebug

end module piernikdebug
