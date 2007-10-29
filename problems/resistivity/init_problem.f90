module init_problem
  
! Initial condition for resistiviy test problem
! Written by: R. Pawlaszek November 2006

! Nie uzywajmy polskich znakow, bo jest wiele sposobow kodowania
!
! Nie zmieniać nic tutaj!! To jest testowy init_problem do sprawdzania symetrii
! kodu. Skompilować kodem Magneto. W pliku problem par są dwie wielkości: 
!       bla - określa w jakim kierunku chcemy sprawdzić zmianę
!       bli - określa, która składowa pola magnetycznego ma się zmieniać.
! Sprawdzanie symetrii polega teraz po prostu na określaniu, czy równa ilość
! komórek po obu stronach płaszczyzny prądowej ma przeciwne znaki. 

  use start
  use arrays
  use grid
  use constants
  use mpi_setup
  use resistivity

  integer ibi
  real d0,p0,bxy,bxz,byx,byz,bzx,bzy,v0
  character problem_name*32,run_id*3,bla*1,bli*1

  real radius
  
  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                              d0,p0, bx0,by0,bz0,bla,bli,v0

contains

!-----------------------------------------------------------------------------

  subroutine read_problem_par

    implicit none
  
      problem_name = 'aaa'
      run_id       = 'aaa'
      d0           = 1.0
      p0           = 1.e-3
      bx0          =   0.0
      by0          =   0.0
      bz0          =   0.0
      bla          =   'x'
      bli          =   'y'
      v0           = 1.0
      
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

      rbuff(1) = d0
      rbuff(2) = p0
      rbuff(3) = bx0
      rbuff(4) = by0
      rbuff(5) = bz0
      rbuff(6) = v0
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      
!  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
!                              d0,p0, bx0,by0,bz0,bla,bli,v0

      problem_name = cbuff(1)   
      run_id       = cbuff(2)   

      d0           = rbuff(1)  
      p0           = rbuff(2)  
      bx0          = rbuff(3)  
      by0          = rbuff(4)  
      bz0          = rbuff(5)  
      v0           = rbuff(6)
        
!      write(*,*) 'proc=',proc  
!      write(*,nml=PROBLEM_CONTROL)  
    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob

    implicit none

    integer i,j,k
    real xmid, ymid, zmid
  
    xmid = 0.5*xmax
    ymid = 0.5*ymax
    zmid = 0.5*zmax
       
    u(2:4,:,:,:) = 0.
    
    call read_problem_par
    call random_seed()

    call random_number(u(2:4,:,:,:))

    u(2:4,:,:,:)=v0*(u(2:4,:,:,:)-0.5)

    do k = 1,nz
      do j = 1,ny
        do i = 1,nx

        select case(bla)
          case('x')
          
          select case(bli)
	    case('y')
	    
            b(1,i,j,k) = bx0
            b(3,i,j,k) = bz0
    
          if (abs(x(i)) .LE. xmid) then
            b(2,i,j,k) = by0
          else
            b(2,i,j,k) = -by0
          endif
            
	    case('z')
	    
            b(1,i,j,k) = bx0
            b(2,i,j,k) = by0
    
          if (abs(x(i)) .LE. xmid) then
            b(3,i,j,k) = bz0
          else
            b(3,i,j,k) = -bz0
          endif
	    
	    case default
	     write(*,*) 'Zmiana w x jest zle zadana'
	    end select
	    
          case('y')

            select case(bli)
	     case('x')
	     
            b(2,i,j,k) = by0
	    b(3,i,j,k) = bz0
       
          if (abs(y(j)) .LE. ymid) then
            b(1,i,j,k) = bx0
          else
            b(1,i,j,k) = -bx0
          endif

	      case('z')
	    
            b(2,i,j,k) = by0
	    b(1,i,j,k) = bx0
       
          if (abs(y(j)) .LE. ymid) then
            b(3,i,j,k) = bz0
          else
            b(3,i,j,k) = -bz0
          endif

	    case default
	      write(*,*) 'Zmiana w y jest zle zadana'
	    end select

          case('z')

            select case(bli)
	      case('x')
	      
            b(3,i,j,k) = bz0
            b(2,i,j,k) = by0
         
          if (abs(z(k)) .LE. zmid) then
            b(1,i,j,k) = bx0
          else
            b(1,i,j,k) = -bx0
          endif

	      case('y')
  
            b(3,i,j,k) = bz0
            b(1,i,j,k) = bx0
         
          if (abs(z(k)) .LE. zmid) then
            b(2,i,j,k) = by0
          else
            b(2,i,j,k) = -by0
          endif

	    case default
	      write(*,*) 'Zmiana w z jest zle zadana'
	    end select
          case default
           write(*,*) 'Cos nie tak w init.problem'
        end select
	
	  u(1,i,j,k)   = d0
          u(5,i,j,k)   = p0/(gamma-1.0)
	  u(5,i,j,k)   = u(5,i,j,k) + 0.5*sum(u(2:4,i,j,k)**2,1)/u(1,i,j,k)
	  u(5,i,j,k)   = u(5,i,j,k) + 0.5*sum(b(:,i,j,k)**2,1)
        enddo
      enddo
    enddo


    return
  end subroutine init_prob  

end module init_problem
