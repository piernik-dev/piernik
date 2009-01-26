#include "piernik.def"

module init_problem
  
! Initial condition for tearing instability problem
! Written by: R.K. Pawlaszek July 2007

!       dimdir - sets the direction of magnetic field change
!                choose between: 'x', 'y', 'z'
!
!       magdir - sets the magnetic field component to change
!                choose between: 'x', 'y', 'z'
!
!       dimdir can't be equal magdir!!

!  use start
  use initionized, only : idni,imxi,imyi,imzi
#ifndef ISO
  use initionized, only : ieni
#endif 
  use fluidindex, only : ibx,iby,ibz
  use arrays, only : u,b
  use grid
  use constants
  use mpi_setup
!  use resistivity

  integer   ibi
  real      beta,v0
  character problem_name*32,run_id*3

  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                              beta,v0

contains

!-----------------------------------------------------------------------------

  subroutine read_problem_par

    implicit none
  
      problem_name =  'tearing'
      run_id       =  'tst'
      beta         =   1.0
      v0           =   0.1
      
    if(proc .eq. 0) then
    
      open(1,file='problem.par')
        read(unit=1,nml=PROBLEM_CONTROL)
        write(*,nml=PROBLEM_CONTROL)
      close(1)
      open(3, file='tmp.log', position='append')
        write(3,nml=PROBLEM_CONTROL)
        write(3,*)
      close(3)
    endif


    if(proc .eq. 0) then

      cbuff(1) =  problem_name
      cbuff(2) =  run_id
      
      rbuff(1) = beta
      rbuff(2) = v0
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      

      problem_name = cbuff(1)   
      run_id       = cbuff(2)   

      beta         = rbuff(1)  
      v0           = rbuff(2)
        
    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob

    implicit none

    integer i,j,k
    real xmid, ymid, zmid, vzab, b0
  
    xmid = 0.5*xmax
    ymid = 0.5*ymax
    zmid = 0.5*zmax
       

    u(idni,:,:,:) = 1.0
    u(imyi,:,:,:) = 0.0
    u(imzi,:,:,:) = 0.0

    b(ibx,:,:,:)  = 0.0
    b(ibz,:,:,:)  = 0.0

    call read_problem_par

    b0 = 1.0

    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          
           vzab = v0*dcos(2.*pi*y(j))
           u(imxi,i,j,k) = u(idni,i,j,k)*vzab
    
           if (abs(x(i)) .LE. xmid) then
               b(iby,i,j,k) = -b0
           else
               b(iby,i,j,k) =  b0
           endif
        enddo
      enddo
    enddo

#ifndef ISO
    u(ieni,:,:,:)   = 0.5*beta &
                      + 0.5*(u(imxi,:,:,:)**2 + u(imyi,:,:,:)**2 &
                             + u(imzi,:,:,:)**2) / u(idni,:,:,:)

    u(ieni,:,:,:)   = u(ieni,:,:,:) + 0.5*sum(b(:,:,:,:)**2,1)
#endif /* ISO */

    return
  end subroutine init_prob  

end module init_problem
