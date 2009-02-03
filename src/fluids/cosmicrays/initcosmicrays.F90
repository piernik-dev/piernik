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

module initcosmicrays

  implicit none

    real                  :: cfl_cr,smallecr,cr_active
    real                  :: gamma_cr,cr_eff,K_cr_paral,K_cr_perp 

    integer               :: iecr

    integer, allocatable, dimension(:)  :: iarr_crs 
    integer, allocatable, dimension(:)  :: iarr_crs_swpx, iarr_crs_swpy, iarr_crs_swpz

 contains


  subroutine init_cosmicrays
    use errh, only : namelist_errh
    use mpisetup

    implicit none
    integer :: ierrh
    character(LEN=100) :: par_file, tmp_log_file

    namelist /COSMIC_RAYS/ cfl_cr,smallecr,cr_active,gamma_cr,cr_eff,K_cr_paral,K_cr_perp
    
    cfl_cr     = 0.9
    smallecr   = 0.0
    cr_active  = 1.0

    gamma_cr   = 4./3.
    cr_eff     = 0.1       !  canonical conversion rate of SN en.-> CR
                           !  we fix E_SN=10**51 erg
    K_cr_paral = 0.0
    K_cr_perp  = 0.0
    
      if(proc .eq. 0) then
         par_file = trim(cwd)//'/problem.par'
         tmp_log_file = trim(cwd)//'/tmp.log'
         open(1,file=par_file)
            read(unit=1,nml=COSMIC_RAYS,iostat=ierrh)
            call namelist_errh(ierrh,'COSMIC_RAYS')
         close(1)
         open(3, file='tmp.log', position='append')
           write(3,nml=COSMIC_RAYS)
           write(3,*)
         close(3)
      endif 
 
 
    if(proc .eq. 0) then
      
      rbuff(1)   = cfl_cr     
      rbuff(2)   = smallecr  
      rbuff(3)   = cr_active  
      rbuff(4)   = gamma_cr  
      rbuff(5)   = cr_eff           
      rbuff(6)   = K_cr_paral 
      rbuff(7)   = K_cr_perp  
    
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else
    
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      
      cfl_cr     = rbuff(1)  
      smallecr   = rbuff(2)  
      cr_active  = rbuff(3)  
      gamma_cr   = rbuff(4)  
      cr_eff     = rbuff(5)        
      K_cr_paral = rbuff(6)  
      K_cr_perp  = rbuff(7)  

    endif

  end subroutine init_cosmicrays



  subroutine cosmicray_index(nvar,nvar_crs)
  
    implicit none
    integer :: nvar, nvar_crs

      iecr          = nvar + 1           
      nvar_crs      = 1
      nvar          = iecr
      
      allocate(iarr_crs(nvar_crs),iarr_crs_swpx(nvar_crs), iarr_crs_swpy(nvar_crs), iarr_crs_swpz(nvar_crs))

      iarr_crs      = [iecr] 
      iarr_crs_swpx = [iecr]
      iarr_crs_swpy = [iecr]
      iarr_crs_swpz = [iecr]
      
   end subroutine cosmicray_index

end module initcosmicrays

