#define RNG nb+1:nx-nb, nb+1:ny-nb, nb+1:nz-nb
module initproblem

! Initial condition for Sedov-Taylor explosion

  use start, only : c_si, gamma
  use arrays
  use grid
  use mpisetup

  real t_sn
  integer n_sn
  real d0, Mrms
  character problem_name*32,run_id*3

  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                              d0, c_si, Mrms

contains

!-----------------------------------------------------------------------------

  subroutine read_problem_par
    use version
    implicit none
    integer :: i

      t_sn = 0.0

      problem_name = 'aaa'
      run_id  = 'aaa'
      d0      = 1.0
      c_si    = 0.1
      Mrms    = 5.0

    if(proc .eq. 0) then

      open(1,file='problem.par')
        read(unit=1,nml=PROBLEM_CONTROL)
        write(*,nml=PROBLEM_CONTROL)
      close(1)
      open(3, file='tmp.log', position='append')
        write(3,nml=PROBLEM_CONTROL)
      close(3)
    endif


    if(proc .eq. 0) then


      cbuff(1) =  problem_name
      cbuff(2) =  run_id

      rbuff(1) = d0
      rbuff(2) = Mrms
      rbuff(3) = c_si

      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else

      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

!  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
!                              d0,p0, bx0,by0,bz0, Eexpl,  x0,y0,z0, r0 &
!                n_sn, dt_sn

      problem_name = cbuff(1)
      run_id       = cbuff(2)

      d0           = rbuff(1)
      Mrms         = rbuff(2)
      c_si         = rbuff(3)

    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob
    use initneutral, only : idnn,imxn,imyn,imzn,ienn, gamma_neu
    use constants, only : dpi
    implicit none

    integer i,j,k,m, n,l
    integer, parameter :: kp = 8
    real, dimension(6) :: mn
    real, dimension(3) :: deltav
    real, dimension(3,nx,ny,nz) :: dv
    real :: somx,somy,somz,rms,cma, vol

! Uniform equilibrium state
    gamma = gamma_neu

    call random_seed()
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          deltav(:) = 0.0
          do m=-kp,kp
            do n = -kp,kp
               call random_number(mn)
!               somx = dpi*(float(n)*y(j) + float(m)*z(k)) / Lx
!               somy = dpi*(float(n)*x(i) + float(m)*z(k)) / Ly
!               somz = dpi*(float(n)*x(i) + float(m)*y(j)) / Lz
!               deltav(1) = deltav(1) + mn(1)*dsin(somx) + mn(2)*dcos(somx)
!               deltav(2) = deltav(2) + mn(3)*dsin(somy) + mn(4)*dcos(somy)
!               deltav(3) = deltav(3) + mn(5)*dsin(somz) + mn(6)*dcos(somz)
               deltav(1) = deltav(1) + mn(1)*dsin(float(m)*z(k)) &
                                     + mn(2)*dcos(float(n)*y(j))
               deltav(2) = deltav(2) + mn(3)*dsin(float(m)*x(i)) &
                                     + mn(4)*dcos(float(n)*z(k))
               deltav(3) = deltav(3) + mn(5)*dsin(float(m)*y(j)) &
                                     + mn(6)*dcos(float(n)*x(i))
            enddo
          enddo
          dv(:,i,j,k) = deltav(:)
       enddo
     enddo
   enddo
   vol = nx*ny*nz
   rms = dsqrt( sum(dv**2) / vol )

   cma = 1.0
   if( rms/c_si < 0.1) then
      cma = rms/c_si
   else
      l=1
      do while (cma >= 0.1 .and. l <=10)
        cma = rms / c_si * (0.1)**l
        l=l+1
      enddo
   endif

    write(*,*) 'cma = ',cma, ' rms = ',rms, ' on ', proc
    write(*,*) 'c_si = ',c_si, ' l = ',l, ' on ', proc
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          u(idnn,i,j,k) = d0
          u(imxn,i,j,k) = u(idnn,i,j,k) * dv(1,i,j,k) * cma
          u(imyn,i,j,k) = u(idnn,i,j,k) * dv(2,i,j,k) * cma
          u(imzn,i,j,k) = u(idnn,i,j,k) * dv(3,i,j,k) * cma
          u(ienn,i,j,k) = c_si**2*d0/(gamma_neu*(gamma_neu-1.0))
          u(ienn,i,j,k) = u(ienn,i,j,k) + 0.5*(u(imxn,i,j,k)**2 &
                         + u(imyn,i,j,k)**2 + u(imzn,i,j,k) )/u(idnn,i,j,k)
          b(:,i,j,k)   = 0.0
          u(ienn,i,j,k)   = u(ienn,i,j,k) + 0.5*sum(b(:,i,j,k)**2,1)
        enddo
      enddo
    enddo

! Explosions

    return
  end subroutine init_prob

end module initproblem

