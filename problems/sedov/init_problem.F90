#define RNG nb+1:nx-nb, nb+1:ny-nb, nb+1:nz-nb
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
          u(idna,i,j,k)   = d0 
          u(imxa,i,j,k) = 0.0
          u(imya,i,j,k) = 0.0
          u(imza,i,j,k) = 0.0
          u(iena,i,j,k)   = p0/(gamma-1.0)
	  u(iena,i,j,k)   = u(iena,i,j,k) + 0.5*(u(imxa,i,j,k)**2 +u(imya,i,j,k)**2 &
	                                        +u(imza,i,j,k)**2)/u(idna,i,j,k)
          b(1,i,j,k)   = bx0
          b(2,i,j,k)   = by0
          b(3,i,j,k)   = bz0
          u(iena(fmagn),i,j,k)   = u(iena(fmagn),i,j,k) + spread(0.5*sum(b(:,i,j,k)**2,1),1,nfmagn)
        enddo
      enddo
    enddo
    
! Explosions


  if(n_sn .eq. 1) then
    do k = ks,ke
      do j = nb+1,ny-nb
        do i = nb+1,nx-nb
          if(((x(i)-x0)**2+(y(j)-y0)**2+(z(k)-z0)**2) .lt. r0**2) then
            u(iena,i,j,k)   = u(iena,i,j,k) + Eexpl
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
            u(iena,i,j,k)   = u(iena,i,j,k) + Eexpl
          endif
        enddo
      enddo
    enddo


  end subroutine random_explosion

  subroutine user_plt(var,ij,xn,tab)
    use start, only : nb
    use arrays, only : u,b,nyb,nzb,nxb,idna,imxa,imya,imza,ibx,iby,ibz
#ifdef ISO
    use arrays, only : iena
#endif
    implicit none
    character(LEN=4)     :: var
    character(LEN=2)     :: ij
    integer              :: xn
    real, dimension(:,:) :: tab

    select case(var)
      case ("dens")
         if(ij=="yz") tab(:,:) = u(idna,xn,nb+1:nyb+nb,nb+1:nzb+nb)
         if(ij=="xz") tab(:,:) = u(idna,nb+1:nxb+nb,xn,nb+1:nzb+nb)
         if(ij=="xy") tab(:,:) = u(idna,nb+1:nxb+nb,nb+1:nyb+nb,xn)
      case ("velx")
         if(ij=="yz") tab(:,:) = u(imxa,xn,nb+1:nyb+nb,nb+1:nzb+nb) / &
                        u(idna,xn,nb+1:nyb+nb,nb+1:nzb+nb)
         if(ij=="xz") tab(:,:) = u(imxa,nb+1:nxb+nb,xn,nb+1:nzb+nb) / &
                        u(idna,nb+1:nxb+nb,xn,nb+1:nzb+nb)
         if(ij=="xy") tab(:,:) = u(imxa,nb+1:nxb+nb,nb+1:nyb+nb,xn) / &
                        u(idna,nb+1:nxb+nb,nb+1:nyb+nb,xn)
      case ("vely")
         if(ij=="yz") tab(:,:) = u(imya,xn,nb+1:nyb+nb,nb+1:nzb+nb) / &
                        u(idna,xn,nb+1:nyb+nb,nb+1:nzb+nb)
         if(ij=="xz") tab(:,:) = u(imya,nb+1:nxb+nb,xn,nb+1:nzb+nb) / &
                        u(idna,nb+1:nxb+nb,xn,nb+1:nzb+nb)
         if(ij=="xy") tab(:,:) = u(imya,nb+1:nxb+nb,nb+1:nyb+nb,xn) / &
                        u(idna,nb+1:nxb+nb,nb+1:nyb+nb,xn)
      case ("velz")
         if(ij=="yz") tab(:,:) = u(imza,xn,nb+1:nyb+nb,nb+1:nzb+nb) / &
                        u(idna,xn,nb+1:nyb+nb,nb+1:nzb+nb)
         if(ij=="xz") tab(:,:) = u(imza,nb+1:nxb+nb,xn,nb+1:nzb+nb) / &
                        u(idna,nb+1:nxb+nb,xn,nb+1:nzb+nb)
         if(ij=="xy") tab(:,:) = u(imza,nb+1:nxb+nb,nb+1:nyb+nb,xn) / &
                        u(idna,nb+1:nxb+nb,nb+1:nyb+nb,xn)
      case ("ener")
#ifndef ISO
         if(ij=="yz") tab(:,:) = u(iena,xn,nb+1:nyb+nb,nb+1:nzb+nb)
         if(ij=="xz") tab(:,:) = u(iena,nb+1:nxb+nb,xn,nb+1:nzb+nb)
         if(ij=="xy") tab(:,:) = u(iena,nb+1:nxb+nb,nb+1:nyb+nb,xn)
#else
         if(ij=="yz") tab(:,:) = 0.5 * ( 
                          u(imxa,xn,nb+1:nyb+nb,nb+1:nzb+nb)**2 &
                        + u(imya,xn,nb+1:nyb+nb,nb+1:nzb+nb)**2 &
                        + u(imza,xn,nb+1:nyb+nb,nb+1:nzb+nb)**2) / &
                             u(idna,xn,nb+1:nyb+nb,nb+1:nzb+nb)
         if(ij=="xz") tab(:,:) = 0.5 * (
                          u(imxa,nb+1:nxb+nb,xn,nb+1:nzb+nb)**2 &
                         +u(imya,nb+1:nxb+nb,xn,nb+1:nzb+nb)**2 &
                         +u(imza,nb+1:nxb+nb,xn,nb+1:nzb+nb)**2) / &
                             u(idna,nb+1:nxb+nb,xn,nb+1:nzb+nb)
         if(ij=="xy") tab(:,:) = 0.5 * (
                          u(imxa,nb+1:nxb+nb,nb+1:nyb+nb,xn)**2 &
                          u(imya,nb+1:nxb+nb,nb+1:nyb+nb,xn)**2 &
                          u(imza,nb+1:nxb+nb,nb+1:nyb+nb,xn)**2) / &
                             u(idna,nb+1:nxb+nb,nb+1:nyb+nb,xn)
#endif /* ISO */
      case ("magx")
         if(ij=="yz") tab(:,:) = b(ibx,xn,nb+1:nyb+nb,nb+1:nzb+nb)
         if(ij=="xz") tab(:,:) = b(ibx,nb+1:nxb+nb,xn,nb+1:nzb+nb)
         if(ij=="xy") tab(:,:) = b(ibx,nb+1:nxb+nb,nb+1:nyb+nb,xn)
      case ("magy")
         if(ij=="yz") tab(:,:) = b(iby,xn,nb+1:nyb+nb,nb+1:nzb+nb)
         if(ij=="xz") tab(:,:) = b(iby,nb+1:nxb+nb,xn,nb+1:nzb+nb)
         if(ij=="xy") tab(:,:) = b(iby,nb+1:nxb+nb,nb+1:nyb+nb,xn)
      case ("magz")
         if(ij=="yz") tab(:,:) = b(ibz,xn,nb+1:nyb+nb,nb+1:nzb+nb)
         if(ij=="xz") tab(:,:) = b(ibz,nb+1:nxb+nb,xn,nb+1:nzb+nb)
         if(ij=="xy") tab(:,:) = b(ibz,nb+1:nxb+nb,nb+1:nyb+nb,xn)

      case default
         write(*,*) "Unknonw var = ",var," in user_plt!"
   end select
 
  end subroutine user_plt

  subroutine user_hdf5(var,tab)
    use arrays, only : u,b,idna,imxa,imya,imza,ibx,iby,ibz
#ifndef ISO
    use arrays, only : iena
#endif
    implicit none
    character(LEN=4)     :: var
    real*4, dimension(:,:,:) :: tab

    select case(var)
      case("dens")
         tab(:,:,:) = real(u(idna,RNG),4)
      case("velx")
         tab(:,:,:) = real(u(imxa,RNG) / u(idna,RNG),4)
      case("vely")
         tab(:,:,:) = real(u(imya,RNG) / u(idna,RNG),4)
      case("velz")
         tab(:,:,:) = real(u(imza,RNG) / u(idna,RNG),4)
      case("ener")
#ifdef ISO
         tab(:,:,:) = real(0.5 *( u(imxa,RNG)**2 + u(imya,RNG)**2 + u(imza,RNG)**2 ) &
            / u(idna,RNG),4)
#else
         tab(:,:,:) = real(u(iena,RNG),4)
#endif
      case("magx")
         tab(:,:,:) = real(b(ibx,RNG),4)
      case("magy")
         tab(:,:,:) = real(b(iby,RNG),4)
      case("magz")
         tab(:,:,:) = real(b(ibz,RNG),4)
      case default
         write(*,*) "Error in user_hdf5"
    end select

  end subroutine user_hdf5

end module init_problem

