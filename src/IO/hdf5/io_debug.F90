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

module piernikiodebug
! pulled by DEBUG

   implicit none

   private
   public :: init_piernikiodebug, force_dumps

   logical, protected :: force_hdf5_dump           !< dump hdf5 every sweep regardless of dataio_pub::dt_hdf
   logical, protected :: force_res_dump            !< dump restart every sweep regardless of dataio_pub::dt_res
   logical, protected :: force_allbnd_dump         !< dump restart with all boundaries every sweep regardless of dataio_pub::dt_res
   logical, protected :: force_log_dump            !< dump log every sweep regardless of dataio_pub:dt_log

   namelist /PIERNIK_IO_DEBUG/ force_hdf5_dump, force_res_dump, force_allbnd_dump, force_log_dump

contains

!>
!! \brief Routine to set debug auxiliary parameters
!!
!! \n \n
!! @b PIERNIK_IO_DEBUG
!! \n \n
!! <table border="+1">
!!   <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!!   <tr><td>force_hdf5_dump  </td><td>       </td><td>logical value          </td><td>\copydoc piernikiodebug::force_hdf5_dump  </td></tr>
!!   <tr><td>force_res_dump   </td><td>       </td><td>logical value          </td><td>\copydoc piernikiodebug::force_res_dump   </td></tr>
!!   <tr><td>force_allbnd_dump</td><td>       </td><td>logical value          </td><td>\copydoc piernikiodebug::force_allbnd_dump</td></tr>
!!   <tr><td>force_log_dump   </td><td>       </td><td>logical value          </td><td>\copydoc piernikiodebug::force_log_dump   </td></tr>
!! </table>
!! The list is active while \b "DEBUG" is defined.
!! \n \n
!<

   subroutine init_piernikiodebug

      use constants,  only: PIERNIK_INIT_MPI
      use dataio_pub, only: nh  ! QA_WARN required for diff_nml
      use dataio_pub, only: code_progress, die
      use mpisetup,   only: master, slave, lbuff, piernik_MPI_Bcast

      implicit none

      if (code_progress < PIERNIK_INIT_MPI) call die("[io_debug:init_piernikdebug] MPI not initialized.")

      if (master) then
         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=PIERNIK_IO_DEBUG)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=PIERNIK_IO_DEBUG, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "PIERNIK_IO_DEBUG")
         read(nh%cmdl_nml,nml=PIERNIK_IO_DEBUG, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "PIERNIK_IO_DEBUG", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=PIERNIK_IO_DEBUG)
         close(nh%lun)
         call nh%compare_namelist()

         lbuff(1) = force_hdf5_dump
         lbuff(2) = force_log_dump
         lbuff(3) = force_res_dump
         lbuff(4) = force_allbnd_dump

      endif

      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         force_hdf5_dump   = lbuff(1)
         force_log_dump    = lbuff(2)
         force_res_dump    = lbuff(3)
         force_allbnd_dump = lbuff(4)

      endif

   end subroutine init_piernikiodebug

   subroutine force_dumps

      use constants,    only: LOGF
      use dataio,       only: write_data
      use dataio_pub,   only: warn
      use data_hdf5,    only: write_hdf5
      use restart_hdf5, only: write_restart_hdf5

      implicit none

      if (force_hdf5_dump)   call write_hdf5
      if (force_res_dump)    call write_restart_hdf5
      if (force_allbnd_dump) call warn("[io_debug:make_sweep] force_allbnd_dump has no effect for single-file HDF5 restart files")
      if (force_log_dump)    call write_data(output=LOGF)

   end subroutine force_dumps

end module piernikiodebug
