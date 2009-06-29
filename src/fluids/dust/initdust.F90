! $Id: initdsttral.F90 594 2009-01-21 12:33:44Z xarth $
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

module initdust

  implicit none

    integer               :: idnd, imxd, imyd, imzd
    real                  :: dragc_gas_dust, taus, dalpha

    integer, allocatable, dimension(:)  :: iarr_dst 
    integer, allocatable, dimension(:)  :: iarr_dst_swpx, iarr_dst_swpy, iarr_dst_swpz

  contains


  subroutine init_dust

    use mpisetup
    use errh, only : namelist_errh
    implicit none
    integer :: ierrh      
    character(LEN=100) :: par_file, tmp_log_file

    namelist /FLUID_DUST/ dragc_gas_dust, dalpha
    
      dragc_gas_dust  = 1.0
      dalpha = 1.0

      if(proc .eq. 0) then
         par_file = trim(cwd)//'/problem.par'
         tmp_log_file = trim(cwd)//'/tmp.log'
         open(1,file=par_file)
            read(unit=1,nml=FLUID_DUST,iostat=ierrh)
            call namelist_errh(ierrh,'FLUID_DUST')
         close(1)
         open(3, file='tmp.log', position='append')
           write(3,nml=FLUID_DUST)
           write(3,*)
         close(3)
      endif

 
    if(proc .eq. 0) then

      rbuff(1)   = dragc_gas_dust
      rbuff(2)   = dalpha
    
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else
    
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      
      dragc_gas_dust  = rbuff(1)  
      dalpha          = rbuff(2)

    endif
    taus = 1. / dragc_gas_dust

  end subroutine init_dust



  subroutine dust_index(nvar,nvar_dst)
  
    implicit none
    integer :: nvar, nvar_dst

   
      idnd = nvar + 1
      imxd = nvar + 2
      imyd = nvar + 3
      imzd = nvar + 4
           
      nvar_dst      = 4
      nvar          = imzd

      allocate(iarr_dst(nvar_dst),iarr_dst_swpx(nvar_dst), iarr_dst_swpy(nvar_dst), iarr_dst_swpz(nvar_dst))

      iarr_dst      = [idnd,imxd,imyd,imzd] 
      iarr_dst_swpx = [idnd,imxd,imyd,imzd]
      iarr_dst_swpy = [idnd,imyd,imxd,imzd]
      iarr_dst_swpz = [idnd,imzd,imyd,imxd]

      
   end subroutine dust_index

end module initdust
