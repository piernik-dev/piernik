! $Id$
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
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen 
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD 
!             for original source code "mhd.f90" 
!   
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.def"
module timer

   implicit none

   integer :: nzones, cpuhours, cpumins, cpusecs , wchours , wcmins  , wcsecs
   real    :: zcps,  cputot, cpuallp, wctot, cpu_start, cpu_stop
   integer :: iarray(3)
   real(kind=4), dimension(2) ::  tarray

contains

   subroutine timer_start
      implicit none
      real(kind=4) :: dtime
!
!  Initialise cpu and wall clocks.  "cputot" will be the total cpu
!  time (in seconds) consumed by this job.  "wctot" will be the total
!  elapsed wall clock time (in seconds) since the job began.
!

      call itime ( iarray )
      wctot  = iarray(1) * 3600. + iarray(2) * 60. + iarray(3)

      cputot  = dtime ( tarray )
      cpu_start = tarray(1) +tarray(2)

   end subroutine timer_start

!------------------------------------------------------------------------------------------

   subroutine timer_stop
      use start, only: t,dt, tend, nstep, nend
      use grid, only :nxd,nyd,nzd
      use mpisetup
      use dataio, only : log_file,log_lun

      implicit none
      real(kind=4) :: dtime

!      Final wall clock time, expressed in hours, minutes, and seconds.
!
      call itime ( iarray )
      wctot  = iarray(1) * 3600. + iarray(2) * 60. + iarray(3) - wctot
      wchours  =  int ( wctot / 3600.0 )
      wcmins   =  int ( wctot / 60.0   ) - 60   * wchours
      wcsecs   =  int ( wctot + 0.5    ) - 3600 * wchours &
                                         - 60   * wcmins
!
!      cpu usage, expressed in hours, minutes, and seconds.
!
      cputot  = dtime ( tarray )
      cpu_stop = tarray(1) +tarray(2)
      cputot = cpu_stop-cpu_start

      cpuhours =  int ( cputot / 3600.0 )
      cpumins  =  int ( cputot / 60.0   ) - 60   * cpuhours
      cpusecs  =  int ( cputot + 0.5    ) - 3600 * cpuhours &
                                        - 60   * cpumins
      nzones = nxd * nyd * nzd


      call MPI_REDUCE(cputot, cpuallp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)

      if(proc == 0) then

         zcps  = real(nstep) * real(nzones) / cpuallp

         write(*,*)
         write (*, 10) cpuallp
         write (*, 20) wctot
         write (*, 30) zcps
         write(*,*)

         open(log_lun, file=log_file, position='append')
            write(log_lun,*)
            write (log_lun, 10) cpuallp
            write (log_lun, 20) wctot
            write (log_lun, 30) zcps
            write(log_lun,*)
         close(log_lun)
      endif

10    format('CPU time        = ', f12.2,' s')
20    format('Wall clock time = ', f12.2,' s')
30    format('Zone-cycles / s = ',1pe12.5)

   end subroutine timer_stop

end module timer


