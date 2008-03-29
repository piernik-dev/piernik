#include "mhd.def"

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

  use start
  use arrays
  use grid
  use constants
  use mpi_setup
  use resistivity

  integer   ibi
  real      d0,v0
  character problem_name*32,run_id*3,dimdir*1,magdir*1

  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                              d0,dimdir,magdir,v0

contains

!-----------------------------------------------------------------------------

  subroutine read_problem_par

    implicit none
  
      problem_name =  'tearing'
      run_id       =  'tst'
      d0           =   1.0
      dimdir       =  'x'
      magdir       =  'y'
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
      cbuff(3) =  dimdir
      cbuff(4) =  magdir
      
      rbuff(1) = d0
      rbuff(2) = v0
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      
!  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
!                              d0,p0, bx0,by0,bz0,dimdir,magdir,v0

      problem_name = cbuff(1)   
      run_id       = cbuff(2)   
      dimdir       = cbuff(3)
      magdir       = cbuff(4)

      d0           = rbuff(1)  
      v0           = rbuff(2)
        
!      write(*,*) 'proc=',proc  
!      write(*,nml=PROBLEM_CONTROL)  
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
       
    u(2:4,:,:,:) = 0.
    call read_problem_par

    b0 = sqrt(2.*alpha*d0*c_si**2)

    do k = 1,nz
      do j = 1,ny
        do i = 1,nx

           u(1,i,j,k)   = d0
        select case(dimdir)
          case('x')
          
          select case(magdir)
            case('y')
            
            b(1,i,j,k) = 0.0
            b(3,i,j,k) = 0.0
            vzab = v0*dcos(2.*pi*y(j))
            u(2,i,j,k) = u(1,i,j,k)*vzab
    
          if (abs(x(i)) .LE. xmid) then
            b(2,i,j,k) = -b0
          else
            b(2,i,j,k) =  b0
          endif
            
            case('z')
            
            b(1,i,j,k) = 0.0
            b(2,i,j,k) = 0.0
            vzab = v0*dcos(2.*pi*z(k))
            u(2,i,j,k) = u(1,i,j,k)*vzab
    
          if (abs(x(i)) .LE. xmid) then
            b(3,i,j,k) = -b0
          else
            b(3,i,j,k) =  b0
          endif
            
            case default
             write(*,*) 'Change in x-direction is improper'
            end select
            
          case('y')

            select case(magdir)
             case('x')
             
            b(2,i,j,k) = 0.0
            b(3,i,j,k) = 0.0
            vzab = v0*dcos(2.*pi*x(i))
            u(3,i,j,k) = u(1,i,j,k)*vzab
       
          if (abs(y(j)) .LE. ymid) then
            b(1,i,j,k) = -b0
          else
            b(1,i,j,k) =  b0
          endif

              case('z')
            
            b(2,i,j,k) = 0.0
            b(1,i,j,k) = 0.0
            vzab = v0*dcos(2.*pi*z(k))
            u(3,i,j,k) = u(1,i,j,k)*vzab
       
          if (abs(y(j)) .LE. ymid) then
            b(3,i,j,k) = -b0
          else
            b(3,i,j,k) =  b0
          endif

            case default
              write(*,*) 'Change in y-direction is improper'
            end select

          case('z')

            select case(magdir)
              case('x')
              
            b(3,i,j,k) = 0.0
            b(2,i,j,k) = 0.0
            vzab = v0*dcos(2.*pi*x(i))
            u(4,i,j,k) = u(1,i,j,k)*vzab
         
          if (abs(z(k)) .LE. zmid) then
            b(1,i,j,k) = -b0
          else
            b(1,i,j,k) =  b0
          endif

              case('y')
  
            b(3,i,j,k) = 0.
            b(1,i,j,k) = 0.
            vzab = v0*dcos(2.*pi*z(k))
            u(4,i,j,k) = u(1,i,j,k)*vzab
         
          if (abs(z(k)) .LE. zmid) then
            b(2,i,j,k) = -b0
          else
            b(2,i,j,k) =  b0
          endif

            case default
              write(*,*) 'Change in z-direction is improper'
            end select
          case default
           write(*,*) 'Something is wrong in init_problem.f90'
        end select
        
          u(5,i,j,k)   = c_si**2/(gamma-1.0)*u(1,i,j,k) &
                         + 0.5*sum(u(2:4,i,j,k)**2,1)/u(1,i,j,k)
          u(5,i,j,k)   = u(5,i,j,k) + 0.5*sum(b(:,i,j,k)**2,1)
        enddo
      enddo
    enddo


    return
  end subroutine init_prob  

end module init_problem
