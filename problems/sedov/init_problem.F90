module init_problem
  
! Initial condition for Sedov-Taylor explosion
! Written by: M. Hanasz, March 2006

  use start
  use arrays
  use grid
  use mpi_setup

  real t_sn
  integer n_sn
  real d0,p0,bx0,by0,bz0,Eexpl, x0,y0,z0,r0, dt_sn
  character problem_name*32,run_id*3

  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                              d0,p0, bx0,by0,bz0, Eexpl,  x0,y0,z0, r0, &
			      n_sn, dt_sn 

contains

!-----------------------------------------------------------------------------

  subroutine read_problem_par
    use comp_log
    implicit none
    integer :: i
    
      t_sn = 0.0
  
      problem_name = 'aaa'
      run_id  = 'aaa'
      d0      = 1.0
      p0      = 1.e-3
      bx0     =   0.
      by0     =   0.
      bz0     =   0.
      Eexpl   = 1.e5
      x0      = 0.0 
      y0      = 0.0 
      z0      = 0.0 
      r0      = dxmn/2.
      n_sn   = 1
      dt_sn = 0.0
         
    if(proc .eq. 0) then
    
      open(1,file='problem.par')
        read(unit=1,nml=PROBLEM_CONTROL)
        write(*,nml=PROBLEM_CONTROL)
      close(1)
      open(3, file='tmp.log', position='append')
        write(3,nml=PROBLEM_CONTROL)
        write(3,*)
        do i=1,nenv
           write(3,*) trim(env(i))
        enddo
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
      rbuff(6) = Eexpl
      rbuff(7) = x0
      rbuff(8) = y0
      rbuff(9) = z0
      rbuff(10)= r0
      rbuff(11)= dt_sn
      
      ibuff(1) = n_sn
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      
!  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
!                              d0,p0, bx0,by0,bz0, Eexpl,  x0,y0,z0, r0 &
!			       n_sn, dt_sn 

      problem_name = cbuff(1)   
      run_id       = cbuff(2)   

      d0           = rbuff(1)  
      p0           = rbuff(2)  
      bx0          = rbuff(3)  
      by0          = rbuff(4)  
      bz0          = rbuff(5)  
      Eexpl        = rbuff(6)  
      x0           = rbuff(7) 
      y0           = rbuff(8)
      z0           = rbuff(9)
      r0           = rbuff(10)
      dt_sn        = rbuff(11)
      
      n_sn         = ibuff(1)  

    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob

    implicit none

    integer i,j,k, n
    
    
    call read_problem_par
 
! Uniform equilibrium state

    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          u(1,i,j,k)   = d0 
          u(2:4,i,j,k) = 0.0
          u(5,i,j,k)   = p0/(gamma-1.0)
	  u(5,i,j,k)   = u(5,i,j,k) + 0.5*sum(u(2:4,i,j,k)**2,1)/u(1,i,j,k)
          b(1,i,j,k)   = bx0
          b(2,i,j,k)   = by0
          b(3,i,j,k)   = bz0
          u(5,i,j,k)   = u(5,i,j,k) + 0.5*sum(b(:,i,j,k)**2,1)
        enddo
      enddo
    enddo
    
! Explosions


  if(n_sn .eq. 1) then
    do k = ks,ke
      do j = nb+1,ny-nb
        do i = nb+1,nx-nb
          if(((x(i)-x0)**2+(y(j)-y0)**2+(z(k)-z0)**2) .lt. r0**2) then
            u(5,i,j,k)   = u(5,i,j,k) + Eexpl
          endif
        enddo
      enddo
    enddo
  else if (n_sn .gt. 1) then
  
!    call random_seed()
  
!    do n=2,n_sn
!      call random_explosion
!    enddo
!  else
!    write(*,*) 'n_sn =', n_sn
!    stop
  endif
    
    return
  end subroutine init_prob  


  subroutine random_explosion
  
  implicit none
    integer i,j,k, n, nexpl
    real rand(3)
    
    call random_number(rand)

    x0 = xmin + (xmax-xmin)*rand(1)
    y0 = ymin + (ymax-ymin)*rand(2)
    z0 = zmin + (zmax-zmin)*rand(3)

    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if(((x(i)-x0)**2+(y(j)-y0)**2+(z(k)-z0)**2) .lt. r0**2) then
            u(5,i,j,k)   = u(5,i,j,k) + Eexpl
          endif
        enddo
      enddo
    enddo


  end subroutine random_explosion

end module init_problem

