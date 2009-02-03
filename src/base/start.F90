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

module start

   use mpisetup

   implicit none

   real    :: t, dt, tend, cfl, smalld, smallei, rorder, dt_cr
   real    :: csi2, csim2, amp_ecr_sn, ethu, f_sn
   integer :: nstep, istep, integration_order, nend, nstep_start
   real, dimension(1) :: gamma   !!! do poprawy
   real, dimension(3,2)  :: cni

!-------------------------------------------------------------------------------
   contains

      subroutine read_params  
         use errh, only : namelist_errh
         implicit none
         integer :: ierrh
         character(len=100) :: par_file, tmp_log_file
         namelist /END_CONTROL/ nend, tend

         namelist /NUMERICAL_SETUP/  cfl, smalld, smallei, integration_order

         par_file = trim(cwd)//'/problem.par'
         tmp_log_file = trim(cwd)//'/tmp.log'

         nstep  = 0
         t      = 0.0
         dt     = 0.0

         tend   = 1.0
         nend   = 10

         cfl     = 0.7
         smalld  = 1.e-10
         smallei = 1.e-10
         integration_order  = 2

         if(proc == 0) then

            open(1,file=par_file)
               read(unit=1,nml=END_CONTROL,iostat=ierrh)
               call namelist_errh(ierrh,'END_CONTROL')
               read(unit=1,nml=NUMERICAL_SETUP,iostat=ierrh)
               call namelist_errh(ierrh,'NUMERICAL_SETUP')
            close(1)

            open(3, file=tmp_log_file, position='append')
               write(unit=3,nml=END_CONTROL)
               write(unit=3,nml=NUMERICAL_SETUP)
            close(3)

         endif


         if(proc == 0) then

!  namelist /END_CONTROL/ nend, tend

            ibuff(30) = nend

            rbuff(30) = tend

!  namelist /NUMERICAL_SETUP/  cfl, smalld, smallei,
!                              integration_order
            rbuff(80) = cfl
            rbuff(83) = smalld
            rbuff(84) = smallei

            ibuff(80) = integration_order


! Broadcasting parameters

            call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
            call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
            call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

         else

            call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
            call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
            call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

!  namelist /END_CONTROL/ nend, tend

            nend                = ibuff(30)

            tend                = rbuff(30)

!  namelist /NUMERICAL_SETUP/  cfl, smalld, smallei,
!                              integration_order,

            cfl                 = rbuff(80)
            smalld              = rbuff(83)
            smallei             = rbuff(84)

            integration_order   = ibuff(80)

         endif  ! (proc .eq. 0)

!         cn(1:3,1) = (/ 1. , 0.5 , 0.0 /)
!         cn(1:3,2) = (/ 1. , 1.  , 0.0 /)
!         if(integration_order == 1) cn(2,1) = 1.

!-------------------------
         if(integration_order > 2) then
            stop 'For "ORIG" scheme integration_order must be 1 or 2'
         endif

      end subroutine read_params

end module start








