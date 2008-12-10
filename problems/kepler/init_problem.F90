#include "piernik.def"


module init_problem

! Initial condition for Keplerian disk
! Written by: M. Hanasz, March 2006

  use mpi_setup

  real d0, collf, r_max, alfasupp, rsup1, rsup2, rsup3, rsup4, dout
  character problem_name*32,run_id*3, mag_field_orient*32, suppkind*6

  namelist /PROBLEM_CONTROL/  problem_name, run_id, &
                              d0,dout,r_max,collf,alfasupp, &
			      rsup1,rsup2,rsup3,rsup4,suppkind,mag_field_orient    


contains

!-----------------------------------------------------------------------------

  subroutine read_problem_par

    implicit none


    problem_name = 'aaa'
    run_id  = 'aa'
    d0      = 1.0
    dout    = 1.0e-4
    r_max   = 1.0
    collf   = 0.0
    alfasupp= 0.0
    rsup1   = 0.0
    rsup2   = 1.0
    rsup3   = 2.0
    rsup4   = 3.0
    suppkind= 'full' !'smooth' !'sharp' !'smout' !'smpeak'
    mag_field_orient = 'none'
    
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
      cbuff(3) =  suppkind
      cbuff(4) =  mag_field_orient

      rbuff(1) = d0
      rbuff(2) = dout
      rbuff(3) = r_max
      rbuff(4) = collf
      rbuff(5) = alfasupp
      rbuff(6) = rsup1
      rbuff(7) = rsup2
      rbuff(8) = rsup3
      rbuff(9) = rsup4
    
      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    else

      call MPI_BCAST(cbuff, 32*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,    buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,    buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)

      problem_name = cbuff(1)
      run_id       = cbuff(2)
      suppkind     = cbuff(3)
      mag_field_orient = cbuff(4)

      d0           = rbuff(1)
      dout         = rbuff(2)
      r_max	   = rbuff(3)
      collf	   = rbuff(4)
      alfasupp     = rbuff(5)
      rsup1        = rbuff(6)
      rsup2        = rbuff(7)
      rsup3        = rbuff(8)
      rsup4        = rbuff(9)
    
    endif

  end subroutine read_problem_par

!-----------------------------------------------------------------------------

  subroutine init_prob
    use arrays,      only : u,b,x,y,z,nx,ny,nz,nfluid,idna,imxa,imya,imza
    use constants,   only : newtong
    use hydrostatic, only : hydrostatic_zeq
    use start,       only : alpha,c_si,dimensions,r_smooth,r_grav,n_gravr,ptmass,smalld,collfaq
#ifndef ISO
    use arrays,      only : iena,fadiab,fmagn
    use start,       only : smallei,gamma
#endif /* !ISO */
#ifdef KEPLER_SUPPRESSION
    use arrays,      only : alfsup,omx0,omy0
#endif /* KEPLER_SUPPRESSION */
    implicit none

    integer i,j,k,kmid
    real xi, yj, rc, vx, vy, vz, b0, sqr_gm
    real, allocatable ::dprof(:), gpot(:)

!    call read_problem_par

!   Secondary parameters

    sqr_gm = sqrt(newtong*ptmass)

    b0 = sqrt(2.*alpha*d0*c_si**2) 
    write(*,*) 'b0 = ',b0 
    allocate(dprof(nz),gpot(nz))

    do k=1, nz
      if(z(k) .lt. 0.0) kmid = k       ! the midplane is in between
    enddo                                  ! ksmid and ksmid+1

    collfaq = collf*3.4e15*(8000.0/8000.0)**0.375 !*cm**3/gram/sek

#ifdef KEPLER_SUPPRESSION
    allocate(alfsup(nx,ny))
    allocate(omx0(nfluid,nx,ny,nz))
    allocate(omy0(nfluid,nx,ny,nz))
#endif /* KEPLER_SUPPRESSION */
      
    do j = 1,ny
      yj = y(j)
      do i = 1,nx
        xi = x(i)
        rc = sqrt(xi**2+yj**2)

        if(dimensions .eq.'3d') then
          call hydrostatic_zeq(i, j, d0, dprof)
        endif

        do k = 1,nz

          vx = sqr_gm * (-yj)/(rc**2+r_smooth**2)**0.75
          vy = sqr_gm * ( xi)/(rc**2+r_smooth**2)**0.75
          vz = 0.0

          u(idna,i,j,k) = min((rc/r_grav)**n_gravr,100.0)
          if(dimensions .eq.'3d') then
            u(idna,i,j,k) = dout + (dprof(k)-dout)/cosh(u(idna,i,j,k))
          else
            u(idna,i,j,k) = dout + (d0 - dout)/cosh(u(idna,i,j,k))
          endif
          u(idna,i,j,k) = max(u(idna,i,j,k), smalld)
#ifdef KEPLER_SUPPRESSION
	  if(suppkind .eq. 'sharp') then
	    if(rc .le. rsup2) then
	      alfsup(i,j) = alfasupp*min((rc-rsup2)/(rsup1-rsup2),1.0)
	    endif
	    if(rc .ge. rsup3) then
	      alfsup(i,j) = alfasupp*min((rc-rsup3)/(rsup4-rsup3),1.0)
	    endif
	  endif
	  if(suppkind .eq. 'smooth') then
	    alfsup(i,j) = alfasupp*(tanh((rc-0.6)*12) + tanh((-rc+0.25)*20)+2.)/2.
	  endif
	  if(suppkind .eq. 'smout') then
	    alfsup(i,j) = alfasupp*(tanh((rc-0.6)*12) + 1.)/2.
	  endif
	  if(suppkind .eq. 'smpeak') then
	    alfsup(i,j) = alfasupp*((tanh((rc-0.4)*12) + 1.)/2.)*((tanh((-rc+0.9)*12)+1.)/2.)
	  endif
	  if(suppkind .eq. 'full') then
	    alfsup = alfasupp
	  endif
#endif /* KEPLER_SUPPRESSION */
                          
          u(imxa,i,j,k) = vx*u(idna,i,j,k)
          u(imya,i,j,k) = vy*u(idna,i,j,k)
          u(imza,i,j,k) = vz*u(idna,i,j,k)      
#ifdef KEPLER_SUPPRESSION
	    omx0(:,i,j,k) = vx
	    omy0(:,i,j,k) = vy
#endif /* KEPLER_SUPPRESSION */
#ifndef ISO	  
          u(iena(fadiab),i,j,k) = c_si**2/(gamma(fadiab)-1.0)*u(idna(fadiab),i,j,k)
          u(iena,i,j,k) = max(u(iena,i,j,k), smallei)
          u(iena(fadiab),i,j,k) = u(iena(fadiab),i,j,k) +0.5*(vx**2+vy**2+vz**2)*u(idna(fadiab),i,j,k)
#endif /* ISO */
          if(trim(mag_field_orient) .eq. 'toroidal') then
            b(1,i,j,k)   = -b0*sqrt(u(idna(1),i,j,k)/d0)*yj/rc
            b(2,i,j,k)   =  b0*sqrt(u(idna(1),i,j,k)/d0)*xi/rc
            b(3,i,j,k)   =  0.0
          else if(trim(mag_field_orient) .eq. 'vertical') then
            b(1,i,j,k)   =  0.0
            b(2,i,j,k)   =  0.0
            b(3,i,j,k)   =  b0
          else if(trim(mag_field_orient) .eq. 'none') then
            b(1,i,j,k)   =  0.0
            b(2,i,j,k)   =  0.0
            b(3,i,j,k)   =  0.0
          endif

#ifndef ISO
          u(iena(fmagn),i,j,k)   = u(iena(fmagn),i,j,k) +0.5*sum(b(:,i,j,k)**2,1)
#endif /* ISO */
        enddo
      enddo
    enddo

    deallocate(dprof,gpot)

    return
  end subroutine init_prob


end module init_problem

