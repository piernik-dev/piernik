! $Id$
#include "piernik.def"

module initionized

  implicit none

    real                  :: gamma_ion, cs_iso_ion,cs_iso_ion2, cs_ion

    integer               :: idni, imxi, imyi, imzi
#ifndef ISO
    integer               :: ieni
#endif /* ISO */    

    integer, allocatable, dimension(:)  :: iarr_ion 
    integer, allocatable, dimension(:)  :: iarr_ion_swpx, iarr_ion_swpy, iarr_ion_swpz

 contains


  subroutine init_ionized

    use mpi_setup

    implicit none
          
    character par_file*(100), tmp_log_file*(100)

    namelist /FLUID_IONIZED/ gamma_ion, cs_iso_ion, cs_ion
    
      gamma_ion  = 1.66666666
      cs_iso_ion = 1.0

      if(proc .eq. 0) then
         par_file = trim(cwd)//'/problem.par'
         tmp_log_file = trim(cwd)//'/tmp.log'
         open(1,file=par_file)
            read(unit=1,nml=FLUID_IONIZED,end=1234,err=1234)
         close(1)

	 goto 4567	 
 1234    write(*,*) 'Check namelist in "problem.par" :'
  	 write(*,*) ''
         write(*,nml=FLUID_IONIZED) 
	 write(*,*) ''
	 stop	 
 4567    continue
	 	 
         open(3, file='tmp.log', position='append')
           write(3,nml=FLUID_IONIZED)
           write(3,*)
         close(3)
	 
	 
      endif 
 
 
    if(proc .eq. 0) then

      rbuff(1)   = gamma_ion
      rbuff(2)   = cs_iso_ion
      rbuff(3)   = cs_ion
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      
      gamma_ion  = rbuff(1)  
      cs_iso_ion = rbuff(2)  
      cs_ion     = rbuff(3)  

    endif
    
    cs_iso_ion2  = cs_iso_ion**2

  end subroutine init_ionized



  subroutine ionized_index(nvar,nvar_ion)
  
    implicit none
    integer :: nvar, nvar_ion

   
      idni = nvar + 1
      imxi = nvar + 2
      imyi = nvar + 3
      imzi = nvar + 4
           
#ifdef ISO
      nvar_ion      = 4
      nvar          = imzi

      allocate(iarr_ion(nvar_ion),iarr_ion_swpx(nvar_ion), iarr_ion_swpy(nvar_ion), iarr_ion_swpz(nvar_ion))

      iarr_ion      = [idni,imxi,imyi,imzi] 
      iarr_ion_swpx = [idni,imxi,imyi,imzi]
      iarr_ion_swpy = [idni,imyi,imxi,imzi]
      iarr_ion_swpz = [idni,imzi,imyi,imxi]
#else
      ieni          = nvar + 5
      nvar_ion      = 5
      nvar          = ieni 

      allocate(iarr_ion(nvar_ion),iarr_ion_swpx(nvar_ion), iarr_ion_swpy(nvar_ion), iarr_ion_swpz(nvar_ion))     

      iarr_ion      = [idni,imxi,imyi,imzi,ieni] 
      iarr_ion_swpx = [idni,imxi,imyi,imzi,ieni]
      iarr_ion_swpy = [idni,imyi,imxi,imzi,ieni]
      iarr_ion_swpz = [idni,imzi,imyi,imxi,ieni]
#endif /* ISO */   

!      write(*,*) 'ionized_index', iarr_ion_swpx
      
   end subroutine ionized_index

end module initionized
