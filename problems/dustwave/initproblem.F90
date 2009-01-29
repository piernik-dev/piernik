#include "piernik.def"

module initproblem
  
! Initial condition for dust fronts
! Written by: D. Wolt, June 2007

  use arrays
  use start
  use grid

  real      :: d0, v0, v1
  integer   :: m_x, m_y, m_z
  character :: problem_name*32,run_id*3

  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                              d0, v0, v1, m_x, m_y, m_z


 contains

!-----------------------------------------------------------------------------

  subroutine read_problem_par

    implicit none
  

    problem_name = 'aaa'
    run_id   = 'aaa'
    d0       = 1.0
    v0       = 0.0
    v1       = 0.01
    m_x	     = 1
    m_y	     = 0
    m_z	     = 0
    
    if(proc .eq. 0) then
      open(1,file='problem.par')
        read(unit=1,nml=PROBLEM_CONTROL)
      close(1)
      open(3, file='tmp.log', position='append')
        write(3,nml=PROBLEM_CONTROL)
        write(3,*)
      close(3)
    endif

    if(proc .eq. 0) then


      cbuff(1) =  problem_name
      cbuff(2) =  run_id

      rbuff(1) = d0
      rbuff(2) = v0
      rbuff(3) = v1
      
      ibuff(1) = m_x
      ibuff(2) = m_y
      ibuff(3) = m_z
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      
      problem_name = cbuff(1)   
      run_id       = cbuff(2)   

      d0           = rbuff(1)  
      v0	   = rbuff(2)
      v1	   = rbuff(3)
      
      m_x	   = ibuff(1)
      m_y	   = ibuff(2)
      m_z	   = ibuff(3)
    
    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob

    use constants,   only : pi
    use initdust,    only : idnd,imxd,imyd,imzd

    implicit none
 
    integer :: i,j,k
    real    :: xi,yj,zk 
    real    :: k_x,k_y,k_z,k_a
    
    write(*,*) 'initproblem:', pi
    
    k_x = 2.*pi/Lx*real(m_x)
    k_y = 2.*pi/Ly*real(m_y)
    k_z = 2.*pi/Lz*real(m_z)  
    k_a = sqrt(k_x**2+k_y**2+k_z**2)
        
    do i = 1,nx
      do j = 1,ny
        do k = 1,nz

	  u(idnd,i,j,k) = d0
	  u(imxd,i,j,k) = d0*k_x/k_a*(v0 +v1*sin(k_x*x(i)+k_y*y(j)+k_z*z(k)))
          u(imyd,i,j,k) = d0*k_y/k_a*(v0 +v1*sin(k_x*x(i)+k_y*y(j)+k_z*z(k)))
          u(imzd,i,j,k) = d0*k_z/k_a*(v0 +v1*sin(k_x*x(i)+k_y*y(j)+k_z*z(k)))

        enddo
      enddo
    enddo
    
    return
  end subroutine init_prob  
  

end module initproblem

