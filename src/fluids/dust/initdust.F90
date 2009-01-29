! $Id: initdsttral.F90 594 2009-01-21 12:33:44Z xarth $
#include "piernik.def"

module initdust

  implicit none

    integer               :: idnd, imxd, imyd, imzd
    real                  :: dragc_gas_dust

    integer, allocatable, dimension(:)  :: iarr_dst 
    integer, allocatable, dimension(:)  :: iarr_dst_swpx, iarr_dst_swpy, iarr_dst_swpz

  contains


  subroutine init_dust

    use mpisetup

    implicit none
          
    character par_file*(100), tmp_log_file*(100)

    namelist /FLUID_DUST/ dragc_gas_dust
    
      dragc_gas_dust  = 0.0

      if(proc .eq. 0) then
         par_file = trim(cwd)//'/problem.par'
         tmp_log_file = trim(cwd)//'/tmp.log'
         open(1,file=par_file)
            read(unit=1,nml=FLUID_DUST)
         close(1)
         open(3, file='tmp.log', position='append')
           write(3,nml=FLUID_DUST)
           write(3,*)
         close(3)
      endif

 
    if(proc .eq. 0) then

      rbuff(1)   = dragc_gas_dust
    
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else
    
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      
      dragc_gas_dust  = rbuff(1)  

    endif

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
