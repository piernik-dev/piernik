! $Id$
#include "piernik.def"

module initneutral

  implicit none

    real                  :: gamma_neu, cs_iso_neu,cs_iso_neu2

    integer               :: idnn, imxn, imyn, imzn
#ifndef ISO
    integer               :: ienn
#endif /* ISO */    

    integer, allocatable, dimension(:)  :: iarr_neu 
    integer, allocatable, dimension(:)  :: iarr_neu_swpx, iarr_neu_swpy, iarr_neu_swpz

  contains


  subroutine init_neutral

    use mpi_setup

    implicit none
          
    character par_file*(100), tmp_log_file*(100)

    namelist /FLUID_neutral/ gamma_neu, cs_iso_neu
    
      gamma_neu  = 1.66666666
      cs_iso_neu = 1.0

      if(proc .eq. 0) then
         par_file = trim(cwd)//'/problem.par'
         tmp_log_file = trim(cwd)//'/tmp.log'
         open(1,file=par_file)
            read(unit=1,nml=FLUID_neutral)
         close(1)
         open(3, file='tmp.log', position='append')
           write(3,nml=FLUID_neutral)
           write(3,*)
         close(3)
      endif	 
	 
	 
    if(proc .eq. 0) then

      rbuff(1)   = gamma_neu
      rbuff(2)   = cs_iso_neu
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      
      gamma_neu  = rbuff(1)  
      cs_iso_neu = rbuff(2)  

    endif

  end subroutine init_neutral



  subroutine neutral_index(nvar,nvar_neu)
  
    implicit none
    integer, nvar, nvar_neu

   
      idnn = nvar + 1
      imxn = nvar + 2
      imyn = nvar + 3
      imzn = nvar + 4
           
#ifdef ISO
      nvar_neu      = 4
      nvar          = imzn

      allocate(iarr_neu(nvar_neu),iarr_neu_swpx(nvar_neu), iarr_neu_swpy(nvar_neu), iarr_neu_swpz(nvar_neu))

      iarr_neu      = [idnn,imxn,imyn,imzn] 
      iarr_neu_swpx = [idnn,imxn,imyn,imzn]
      iarr_neu_swpy = [idnn,imyn,imxn,imzn]
      iarr_neu_swpz = [idnn,imzn,imyn,imxn]
#else
      ienn          = nvar + 5
      nvar_neu      = 5
      nvar          = ienn 

      allocate(iarr_neu(nvar_neu),iarr_neu_swpx(nvar_neu), iarr_neu_swpy(nvar_neu), iarr_neu_swpz(nvar_neu))     

      iarr_neu      = [idnn,imxn,imyn,imzn,ienn] 
      iarr_neu_swpx = [idnn,imxn,imyn,imzn,ienn]
      iarr_neu_swpy = [idnn,imyn,imxn,imzn,ienn]
      iarr_neu_swpz = [idnn,imzn,imyn,imxn,ienn]
#endif /* ISO */   

!      write(*,*) 'neutral_index', iarr_neu_swpx
      
   end subroutine neutral_index

end module initneutral
