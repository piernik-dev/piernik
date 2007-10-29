module mhdblock

! Based on: Pen Arras & Wong (2003)
! Modified extensively by M. Hanasz November 2005 - May 2006
! Optimalization of cpu-time efficiency of mhdflux, tvd1 by K. Kowalik
! Modification history:  see "changelog"  -- warning: out of date

  use mpi_setup
  use start
  use arrays
  use grid
  use init_problem 
  use fluid_boundaries
  use resistivity
  use mag_boundaries
  use dataio  
  use diagnostics
! use poisson_solver
  use thermal
  use resistivity
  
  implicit none
  real c

contains

!==========================================================================================

  subroutine mhdstep 

  implicit none
  

! wybuchy najlepiej dodac przed obliczeniem kroku czasowego
!DW+
   !if(sn_prob .ne. 'null') call supernovae_distribution(dt)
!DW-
!    if((dt_sn .ne. 0.0) .and. (t .gt. t_sn) ) then
!       call random_explosion
!       t_sn = t_sn + dt_sn
!    endif



    call timestep

    if (dt_log .gt. 0.0) then
      if (nlog .lt. (int(t / dt_log) + 1)) then
        call write_log
        nlog = nlog + 1
      endif
    endif
    
    if(proc.eq.0) write(*,900) nstep,dt,t
900 format('   nstep = ',i7,'   dt = ',f22.16,'   t = ',f22.16)   


    t=t+ 2*dt
    
!    write(*,*) 'magnetic=',magnetic

! x sweep
      call fluidx
      if(magfield) call advectbyzx
    
! y sweep
      call fluidy
      if(magfield) call advectbzxy

! z sweep
    if(dimensions .eq. '3d') then
      call fluidz
      if(magfield) call advectbxyz
    endif
    
!      call write_data(output='all')
!         call write_hdf        
!          nhdf = nhdf + 1
    
!if(selfgravity .ne. 'null') call poisson 

! back z
    if(dimensions .eq. '3d') then
      if(magfield) call advectbxyz   
      call fluidz
    endif

! back y
      if(magfield) call advectbzxy
      call fluidy

! back x
      if(magfield) call advectbyzx
      call fluidx


  end subroutine mhdstep

!==========================================================================================

  subroutine timestep
  
  
    implicit none
  

! Timestep computation

    call timestep_mhd
    dt=min(dt_mhd,(tend-t)/2.)

    if (coolheat .and. coolheat_active .eq. 'yes') then
      call timestep_coolheat
      dt = min(dt,dt_coolheat)
    endif

    if(heatcond) then
      dt_heatcond = cfl_heatcond * dxmn**2/(K_heatcond+small)
      dt = min(dt,dt_heatcond)
    endif

    if(bulk_viscosity) then
      dt_visc = cfl_visc * dxmn**2/(nu_bulk+small)
      dt = min(dt,dt_visc)
!      write(*,*) dt_visc
    endif
  
  end subroutine timestep

!------------------------------------------------------------------------------------------



  subroutine timestep_mhd

    real dt_mhd_proc, dt_mhd_all, c_max_all 
!DW+
    real dt_mhd_proc_x, dt_mhd_proc_y, dt_mhd_proc_z
    real cx, cy, cz, vx, vy, vz, cf
!DW-

!pawlaszek
    real dtdiff
!pawlaszek
    
! locals
    real pmag
    real v,ps,p,bx,by,bz
    integer i,j,k,ip,jp,kp

!DW+    
    cx    = 0.0
    cy    = 0.0
    cz    = 0.0
!DW-
    c     = 0.0


    do k=kb,ke
      do j=nb+1,nb+nyb
        do i=nb+1,nb+nxb

          vx=abs(u(2,i,j,k)/u(1,i,j,k))
          vy=abs(u(3,i,j,k)/u(1,i,j,k))
          vz=abs(u(4,i,j,k)/u(1,i,j,k))
          v = max(vx,vy,vz)
          
          pmag = sum(b(:,i,j,k)**2,1)/2.

          ps=(u(5,i,j,k)-sum(u(2:4,i,j,k)**2,1)/u(1,i,j,k)/2.)*(gamma-1.)+(2.-gamma)*pmag
          p=ps-pmag
          
          cf = sqrt(abs(  (2.*pmag+gamma*p)/u(1,i,j,k)) )
!DW+
          cx=max(cx,vx+cf)
          cy=max(cy,vy+cf)
          cz=max(cz,vz+cf)
          c =max(c,cx,cy,cz)
!DW-

        end do
      end do
    end do

!DW+
    dt_mhd_proc_x = cfl*dx0/cx
    dt_mhd_proc_y = cfl*dy0/cy
    dt_mhd_proc_z = cfl*dz0/cz
    dt_mhd_proc   = min(dt_mhd_proc_x, dt_mhd_proc_y, dt_mhd_proc_z)
!DW-

    call MPI_REDUCE(c, c_max_all, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)
    call MPI_BCAST(c_max_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    
    c = c_max_all
    

!pawlaszek
    if (resist) then
      call timestep_resist
      dt_mhd_proc=min(dt_mhd_proc,dt_resist)
    endif
!pawlaszek
        
    call MPI_REDUCE(dt_mhd_proc, dt_mhd_all, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr)
    call MPI_BCAST(dt_mhd_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    dt_mhd = dt_mhd_all


  end subroutine timestep_mhd

!==========================================================================================

  subroutine fluidx
    implicit none
    real, dimension(3,nx) :: b_x
    real, dimension(5,nx) :: u_x, d_x
    integer j,k,jp,kp
    real yj,zk
  
    d_x = spread(dx,1,5)
    b_x = 0.0

    do k=kb,ke
      kp=k+1
      do j=1,ny-1
        if(magfield)then
          jp=j+1
          b_x=b(:,:,j,k)/2
          b_x(1,:)=b_x(1,:)+eoshift(b_x(1,:),shift=1,boundary=big)
          b_x(2,:)=b_x(2,:)+b(2,:,jp,k)/2
          if(dimensions .eq. '3d')  b_x(3,:)=b_x(3,:)+b(3,:,j,kp)/2
        endif
        u_x=u(:,:,j,k)
        call relaxing_tvd(u_x,b_x,'xsweep',j,k,x,xl,xr,d_x,nx,nb,dt)
        u(:,:,j,k)=u_x
      end do
    end do

!    call maxu('fluidx')

    call bnd_u('xdim')
    call bnd_u('ydim')
    if(dimensions .eq. '3d') call bnd_u('zdim')

end subroutine fluidx

!------------------------------------------------------------------------------------------

  subroutine fluidy
    implicit none
    real, dimension(3,ny) :: b_y
    real, dimension(5,ny) :: u_y,d_y
    integer i,j,k,ip,jp,kp
    real xi,zk

    d_y = spread(dy,1,5)
    b_y = 0.0

    do k=kb,ke
      kp=k+1
      do i=1,nx-1
        ip=i+1
        if(magfield)then
          b_y(:,:)=b(:,i,:,k)/2
          b_y(1,:)=b_y(1,:)+b(1,ip,:,k)/2
          b_y(2,:)=b_y(2,:)+eoshift(b_y(2,:),shift=1,boundary=big)
          if (dimensions .eq. '3d') b_y(3,:)=b_y(3,:)+b(3,i,:,kp)/2
          b_y((/2,1,3/),:)=b_y(:,:)
        endif

        u_y((/1,3,2,4,5/),:)=u(:,i,:,k)
        call relaxing_tvd(u_y,b_y,'ysweep',k,i,y,yl,yr,d_y,ny,nb,dt)
        u(:,i,:,k)=u_y((/1,3,2,4,5/),:)
      end do
    end do

!    call maxu('fluidy')


    call bnd_u('xdim')
    call bnd_u('ydim')
    if(dimensions .eq. '3d') call bnd_u('zdim')
    

end subroutine fluidy

!------------------------------------------------------------------------------------------

  subroutine fluidz
    implicit none
    real, dimension(3,nz) :: b_z
    real, dimension(5,nz) :: u_z,d_z
    integer i,j,k,ip,jp,kp
    real xi,yj

    d_z = spread(dz,1,5)
    b_z = 0.0

    do j=1,ny-1
      jp=j+1
      do i=1,nx-1
        if(magfield)then
          ip=i+1
          b_z(:,:)=b(:,i,j,:)/2
          b_z(1,:)=b_z(1,:)+b(1,ip,j,:)/2
          b_z(2,:)=b_z(2,:)+b(2,i,jp,:)/2
          b_z(3,:)=b_z(3,:)+eoshift(b_z(3,:),shift=1,boundary=big)
          b_z((/3,2,1/),:)=b_z(:,:)
        endif
        u_z((/1,4,3,2,5/),:)=u(:,i,j,:)
        call relaxing_tvd(u_z,b_z,'zsweep',i,j,z,zl,zr,d_z,nz,nb,dt)
        u(:,i,j,:)=u_z((/1,4,3,2,5/),:)
      end do
    end do
    
!    call maxu('fluidz')
 

    call bnd_u('xdim')
    call bnd_u('ydim')
    if(dimensions .eq. '3d') call bnd_u('zdim')

  end subroutine fluidz

!==========================================================================================


  subroutine advectbyzx

IF (trim(magnetic_scheme) .eq. 'orig' .and. magnetic_order .eq. 2) THEN

    call advectby_x

    b(2,:,:,:) = b(2,:,:,:) - wa/dx0
    wa = cshift(wa,shift=-1,dim=1)
    b(2,:,:,:) = b(2,:,:,:) + wa/dx0
    b(1,:,:,:) = b(1,:,:,:) - wa/dy0
    wa = cshift(wa,shift=1,dim=2)
    b(1,:,:,:) = b(1,:,:,:) + wa/dy0

    call bnd_b('xdim')
    call bnd_b('ydim')
    if(dimensions .eq. '3d') call bnd_b('zdim')

  if (dimensions .eq. '3d') then
    call advectbz_x
    
    b(3,:,:,:) = b(3,:,:,:) - wa/dx0
    wa = cshift(wa,shift=-1,dim=1)
    b(3,:,:,:) = b(3,:,:,:) + wa/dx0
    b(1,:,:,:) = b(1,:,:,:) - wa/dz0
    wa = cshift(wa,shift=1,dim=3)
    b(1,:,:,:) = b(1,:,:,:) + wa/dz0
  
    call bnd_b('xdim')
    call bnd_b('ydim')
    call bnd_b('zdim')
  endif

ELSE IF(trim(magnetic_scheme) .eq. 'ssp' .and. magnetic_order .eq. 3) THEN

! b1 is now the basic magnetic field!!
    b1(:,:,:,:) = b(:,:,:,:)
    
! DIFFUSION FIRST STEP for by_x

    call advectby_x
    call diffuseby_x

    b(2,:,:,:) = b1(2,:,:,:) - wcu/dx0
    wcu = cshift(wcu,shift=1,dim=1)
    b(2,:,:,:) = b(2,:,:,:) + wcu/dx0
    wcu = cshift(wcu,shift=-1,dim=1)
    b(1,:,:,:) = b1(1,:,:,:) + wcu/dy0
    wcu = cshift(wcu,shift=1,dim=2)
    b(1,:,:,:) = b(1,:,:,:) - wcu/dy0

! FIRST STEP for by_x

    b(2,:,:,:) = b(2,:,:,:) - wa/dx0
    wa = cshift(wa,shift=-1,dim=1)
    b(2,:,:,:) = b(2,:,:,:) + wa/dx0
    b(1,:,:,:) = b(1,:,:,:) - wa/dy0
    wa = cshift(wa,shift=1,dim=2)
    b(1,:,:,:) = b(1,:,:,:) + wa/dy0

    call bnd_b('xdim')
    call bnd_b('ydim')
    if(dimensions .eq. '3d') call bnd_b('zdim')

! DIFFUSION SECOND STEP for by_x
    
    call advectby_x
    call diffuseby_x

    b(2,:,:,:) = 0.75*b1(2,:,:,:) + 0.25*(b(2,:,:,:) - wcu/dx0)
    wcu = cshift(wcu,shift=1,dim=1)
    b(2,:,:,:) = b(2,:,:,:) + 0.25*wcu/dx0
    wcu = cshift(wcu,shift=-1,dim=1)
    b(1,:,:,:) = 0.75*b1(1,:,:,:) + 0.25*(b(1,:,:,:) + wcu/dy0)
    wcu = cshift(wcu,shift=1,dim=2)
    b(1,:,:,:) = b(1,:,:,:) - 0.25*wcu/dy0

! SECOND STEP for by_x

    b(2,:,:,:) = b(2,:,:,:) - 0.25*wa/dx0
    wa = cshift(wa,shift=-1,dim=1)
    b(2,:,:,:) = b(2,:,:,:) + 0.25*wa/dx0
    b(1,:,:,:) = b(1,:,:,:) - 0.25*wa/dy0
    wa = cshift(wa,shift=1,dim=2)
    b(1,:,:,:) = b(1,:,:,:) + 0.25*wa/dy0

    call bnd_b('xdim')
    call bnd_b('ydim')
    if(dimensions .eq. '3d') call bnd_b('zdim') 

! DIFFUSION THIRD STEP for by_x

    call advectby_x
    call diffuseby_x

    b(2,:,:,:) = 1./3.*b1(2,:,:,:) + 2./3.*(b(2,:,:,:) - wcu/dx0)
    wcu = cshift(wcu,shift=1,dim=1)
    b(2,:,:,:) = b(2,:,:,:) + 2./3.*wcu/dx0
    wcu = cshift(wcu,shift=-1,dim=1)
    b(1,:,:,:) = 1./3.*b1(1,:,:,:) + 2./3.*(b(1,:,:,:) + wcu/dy0)
    wcu = cshift(wcu,shift=1,dim=2)
    b(1,:,:,:) = b(1,:,:,:) - 2./3.*wcu/dy0

! THIRD STEP for by_x

    b(2,:,:,:) = b(2,:,:,:) - 2./3.*wa/dx0
    wa = cshift(wa,shift=-1,dim=1)
    b(2,:,:,:) = b(2,:,:,:) + 2./3.*wa/dx0
    b(1,:,:,:) = b(1,:,:,:) - 2./3.*wa/dy0
    wa = cshift(wa,shift=1,dim=2)
    b(1,:,:,:) = b(1,:,:,:) + 2./3.*wa/dy0

    call bnd_b('xdim')
    call bnd_b('ydim')
    if(dimensions .eq. '3d') call bnd_b('zdim')

! b is the magnetic field back again!!
!-------------------------------------------------------------------------------

  if (dimensions .eq. '3d') then

! b1 is the initial magnetic field
    b1(:,:,:,:) = b(:,:,:,:)

! DIFFUSION FIRST STEP for bz_x

    call advectbz_x
    call diffusebz_x

    b(3,:,:,:) = b1(3,:,:,:) - wcu/dx0
    wcu = cshift(wcu,shift=1,dim=1)
    b(3,:,:,:) = b(3,:,:,:) + wcu/dx0
    wcu = cshift(wcu,shift=-1,dim=1)
    b(1,:,:,:) = b1(1,:,:,:) + wcu/dz0
    wcu = cshift(wcu,shift=1,dim=3)
    b(1,:,:,:) = b(1,:,:,:) - wcu/dz0

! FIRST STEP for bz_x
  
    b(3,:,:,:) = b(3,:,:,:) - wa/dx0
    wa = cshift(wa,shift=-1,dim=1)
    b(3,:,:,:) = b(3,:,:,:) + wa/dx0
    b(1,:,:,:) = b(1,:,:,:) - wa/dz0
    wa = cshift(wa,shift=1,dim=3)
    b(1,:,:,:) = b(1,:,:,:) + wa/dz0
  
    call bnd_b('xdim')
    call bnd_b('ydim')
    call bnd_b('zdim')

! DIFFUSE SECOND STEP for bz_x

    call advectbz_x
    call diffusebz_x

    b(3,:,:,:) = 0.75*b1(3,:,:,:) + 0.25*(b(3,:,:,:) - wcu/dx0)
    wcu = cshift(wcu,shift=1,dim=1)
    b(3,:,:,:) = b(3,:,:,:) + 0.25*wcu/dx0
    wcu = cshift(wcu,shift=-1,dim=1)
    b(1,:,:,:) = 0.75*b1(1,:,:,:) + 0.25*(b(1,:,:,:) + wcu/dz0)
    wcu = cshift(wcu,shift=1,dim=3)
    b(1,:,:,:) = b(1,:,:,:) - 0.25*wcu/dz0

! SECOND STEP for bz_x

    b(3,:,:,:) = b(3,:,:,:) - 0.25*wa/dx0
    wa = cshift(wa,shift=-1,dim=1)
    b(3,:,:,:) = b(3,:,:,:) + 0.25*wa/dx0
    b(1,:,:,:) = b(1,:,:,:) - 0.25*wa/dz0
    wa = cshift(wa,shift=1,dim=3)
    b(1,:,:,:) = b(1,:,:,:) + 0.25*wa/dz0
  
    call bnd_b('xdim')
    call bnd_b('ydim')
    call bnd_b('zdim')

! DIFFUSION THIRD STEP for bz_x

    call advectbz_x
    call diffusebz_x

    b(3,:,:,:) = 1./3.*b1(3,:,:,:) + 2./3.*(b(3,:,:,:) - wcu/dx0)
    wcu = cshift(wcu,shift=1,dim=1)
    b(3,:,:,:) = b(3,:,:,:) + 2./3.*wcu/dx0
    wcu = cshift(wcu,shift=-1,dim=1)
    b(1,:,:,:) = 1./3.*b1(1,:,:,:) + 2./3.*(b(1,:,:,:) + wcu/dz0)
    wcu = cshift(wcu,shift=1,dim=3)
    b(1,:,:,:) = b(1,:,:,:) - 2./3.*wcu/dz0

! THIRD STEP for bz_x

    b(3,:,:,:) = b(3,:,:,:) - 2./3.*wa/dx0
    wa = cshift(wa,shift=-1,dim=1)
    b(3,:,:,:) = b(3,:,:,:) + 2./3.*wa/dx0
    b(1,:,:,:) = b(1,:,:,:) - 2./3.*wa/dz0
    wa = cshift(wa,shift=1,dim=3)
    b(1,:,:,:) = b(1,:,:,:) + 2./3.*wa/dz0
  
    call bnd_b('xdim')
    call bnd_b('ydim')
    call bnd_b('zdim')

!b is main magnetic field back again

  endif

ELSE

900 format(' Scheme ',a4,' of order ',i2,' not implemented for magnetic field advection')   
    write(*,*) ' '
    write(*,900),trim(magnetic_scheme),magnetic_order
    write(*,*) ' '
    STOP
ENDIF

  end subroutine advectbyzx

!------------------------------------------------------------------------------------------

  subroutine advectbxyz
    implicit none

IF(trim(magnetic_scheme) .eq. 'orig' .and. magnetic_order .eq. 2) THEN

    call advectbx_z

    b(1,:,:,:) = b(1,:,:,:) - wa/dz0
    wa = cshift(wa,shift=-1,dim=3)
    b(1,:,:,:) = b(1,:,:,:) + wa/dz0
    b(3,:,:,:) = b(3,:,:,:) - wa/dx0
    wa = cshift(wa,shift=1,dim=1)
    b(3,:,:,:) = b(3,:,:,:) + wa/dx0
  
    call bnd_b('xdim')
    call bnd_b('ydim')
    call bnd_b('zdim')

    call advectby_z

    b(2,:,:,:) = b(2,:,:,:) - wa/dz0
    wa = cshift(wa,shift=-1,dim=3)
    b(2,:,:,:) = b(2,:,:,:) + wa/dz0
    b(3,:,:,:) = b(3,:,:,:) - wa/dy0
    wa = cshift(wa,shift=1,dim=2)
    b(3,:,:,:) = b(3,:,:,:) + wa/dy0
  
    call bnd_b('xdim')
    call bnd_b('ydim')
    call bnd_b('zdim')
    
ELSE IF (trim(magnetic_scheme) .eq. 'ssp' .and. magnetic_order .eq. 3) THEN

! Now b1 is the initial magnetic field
    b1(:,:,:,:) = b(:,:,:,:)

! DIFFUSION FIRST STEP for bx_z

    call advectbx_z
    call diffusebx_z
 
    b(1,:,:,:) = b1(1,:,:,:) - wcu/dz0
    wcu = cshift(wcu,shift=1,dim=3)
    b(1,:,:,:) = b(1,:,:,:) + wcu/dz0
    wcu = cshift(wcu,shift=-1,dim=3)
    b(3,:,:,:) = b1(3,:,:,:) + wcu/dx0
    wcu = cshift(wcu,shift=1,dim=1)
    b(3,:,:,:) = b(3,:,:,:) - wcu/dx0

! FIRST STEP for bx_z
 
    b(1,:,:,:) = b(1,:,:,:) - wa/dz0
    wa = cshift(wa,shift=-1,dim=3)
    b(1,:,:,:) = b(1,:,:,:) + wa/dz0
    b(3,:,:,:) = b(3,:,:,:) - wa/dx0
    wa = cshift(wa,shift=1,dim=1)
    b(3,:,:,:) = b(3,:,:,:) + wa/dx0
  
    call bnd_b('xdim')
    call bnd_b('ydim')
    call bnd_b('zdim')

! DIFFUSION SECOND STEP for bx_z

    call advectbx_z
    call diffusebx_z

    b(1,:,:,:) = 0.75*b1(1,:,:,:) + 0.25*(b(1,:,:,:) - wcu/dz0)
    wcu = cshift(wcu,shift=1,dim=3)
    b(1,:,:,:) = b(1,:,:,:) + 0.25*wcu/dz0
    wcu = cshift(wcu,shift=-1,dim=3)
    b(3,:,:,:) = 0.75*b1(3,:,:,:) + 0.25*(b(3,:,:,:) + wcu/dx0)
    wcu = cshift(wcu,shift=1,dim=1)
    b(3,:,:,:) = b(3,:,:,:) - 0.25*wcu/dx0

! SECOND STEP for bx_z

    b(1,:,:,:) = b(1,:,:,:) - 0.25*wa/dz0
    wa = cshift(wa,shift=-1,dim=3)
    b(1,:,:,:) = b(1,:,:,:) + 0.25*wa/dz0
    b(3,:,:,:) = b(3,:,:,:) - 0.25*wa/dx0
    wa = cshift(wa,shift=1,dim=1)
    b(3,:,:,:) = b(3,:,:,:) + 0.25*wa/dx0
 
    call bnd_b('xdim')
    call bnd_b('ydim')
    call bnd_b('zdim')

! DIFFUSION THIRD STEP for bx_z

    call advectbx_z
    call diffusebx_z

    b(1,:,:,:) = 1./3.*b1(1,:,:,:) + 2./3.*(b(1,:,:,:) - wcu/dz0)
    wcu = cshift(wcu,shift=1,dim=3)
    b(1,:,:,:) = b(1,:,:,:) + 2./3.*wcu/dz0
    wcu = cshift(wcu,shift=-1,dim=3)
    b(3,:,:,:) = 1./3.*b1(3,:,:,:) + 2./3.*(b(3,:,:,:) + wcu/dx0)
    wcu = cshift(wcu,shift=1,dim=1)
    b(3,:,:,:) = b(3,:,:,:) - 2./3.*wcu/dx0

! THIRD STEP for bx_z

    b(1,:,:,:) = b(1,:,:,:) - 2./3.*wa/dz0
    wa = cshift(wa,shift=-1,dim=3)
    b(1,:,:,:) = b(1,:,:,:) + 2./3.*wa/dz0
    b(3,:,:,:) = b(3,:,:,:) - 2./3.*wa/dx0
    wa = cshift(wa,shift=1,dim=1)
    b(3,:,:,:) = b(3,:,:,:) + 2./3.*wa/dx0
  
    call bnd_b('xdim')
    call bnd_b('ydim')
    call bnd_b('zdim')

!--------------------------------------------------------------------------

! Now b1 is the initial magnetic field
    b1(:,:,:,:) = b(:,:,:,:)

! DIFFUSION FIRST STEP for by_z

    call advectby_z
    call diffuseby_z

    b(2,:,:,:) = b1(2,:,:,:) - wcu/dz0
    wcu = cshift(wcu,shift=1,dim=3)
    b(2,:,:,:) = b(2,:,:,:) + wcu/dz0
    wcu = cshift(wcu,shift=-1,dim=3)
    b(3,:,:,:) = b1(3,:,:,:) + wcu/dy0
    wcu = cshift(wcu,shift=1,dim=2)
    b(3,:,:,:) = b(3,:,:,:) - wcu/dy0

! FIRST STEP for by_z

    b(2,:,:,:) = b(2,:,:,:) - wa/dz0
    wa = cshift(wa,shift=-1,dim=3)
    b(2,:,:,:) = b(2,:,:,:) + wa/dz0
    b(3,:,:,:) = b(3,:,:,:) - wa/dy0
    wa = cshift(wa,shift=1,dim=2)
    b(3,:,:,:) = b(3,:,:,:) + wa/dy0
 
    call bnd_b('xdim')
    call bnd_b('ydim')
    call bnd_b('zdim')

! DIFFUSION SECOND STEP for by_z

    call advectby_z
    call diffuseby_z

    b(2,:,:,:) = 0.75*b1(2,:,:,:) + 0.25*(b(2,:,:,:) - wcu/dz0)
    wcu = cshift(wcu,shift=1,dim=3)
    b(2,:,:,:) = b(2,:,:,:) + 0.25*wcu/dz0
    wcu = cshift(wcu,shift=-1,dim=3)
    b(3,:,:,:) = 0.75*b1(3,:,:,:) + 0.25*(b(3,:,:,:) + wcu/dy0)
    wcu = cshift(wcu,shift=1,dim=2)
    b(3,:,:,:) = b(3,:,:,:) - 0.25*wcu/dy0

! SECOND STEP for by_z

    b(2,:,:,:) = b(2,:,:,:) - 0.25*wa/dz0
    wa = cshift(wa,shift=-1,dim=3)
    b(2,:,:,:) = b(2,:,:,:) + 0.25*wa/dz0
    b(3,:,:,:) = b(3,:,:,:) - 0.25*wa/dy0
    wa = cshift(wa,shift=1,dim=2)
    b(3,:,:,:) = b(3,:,:,:) + 0.25*wa/dy0
 
    call bnd_b('xdim')
    call bnd_b('ydim')
    call bnd_b('zdim')

! DIFFUSION THIRD STEP for by_z

    call advectby_z
    call diffuseby_z

    b(2,:,:,:) = 1./3.*b1(2,:,:,:) + 2./3.*(b(2,:,:,:) - wcu/dz0)
    wcu = cshift(wcu,shift=1,dim=3)
    b(2,:,:,:) = b(2,:,:,:) + 2./3.*wcu/dz0
    wcu = cshift(wcu,shift=-1,dim=3)
    b(3,:,:,:) = 1./3.*b1(3,:,:,:) + 2./3.*(b(3,:,:,:) + wcu/dy0)
    wcu = cshift(wcu,shift=1,dim=2)
    b(3,:,:,:) = b(3,:,:,:) - 2./3.*wcu/dy0

! THIRD STEP for by_z

    b(2,:,:,:) = b(2,:,:,:) - 2./3.*wa/dz0
    wa = cshift(wa,shift=-1,dim=3)
    b(2,:,:,:) = b(2,:,:,:) + 2./3.*wa/dz0
    b(3,:,:,:) = b(3,:,:,:) - 2./3.*wa/dy0
    wa = cshift(wa,shift=1,dim=2)
    b(3,:,:,:) = b(3,:,:,:) + 2./3.*wa/dy0
  
    call bnd_b('xdim')
    call bnd_b('ydim')
    call bnd_b('zdim')

ELSE
900 format(' Scheme ',a4,' of order ',i2,' not implemented for magnetic field advection')   
    write(*,*) ' '
    write(*,900),trim(magnetic_scheme),magnetic_order
    write(*,*) ' '
    STOP
ENDIF
  
  end subroutine advectbxyz

!------------------------------------------------------------------------------------------

  subroutine advectbzxy
    implicit none

IF(trim(magnetic_scheme) .eq. 'orig' .and. magnetic_order .eq. 2) THEN

  if (dimensions .eq. '3d') then
    call advectbz_y
  
    b(3,:,:,:) = b(3,:,:,:) - wa/dy0
    wa = cshift(wa,shift=-1,dim=2)
    b(3,:,:,:) = b(3,:,:,:) + wa/dy0
    b(2,:,:,:) = b(2,:,:,:) - wa/dz0
    wa = cshift(wa,shift=1,dim=3)
    b(2,:,:,:) = b(2,:,:,:) + wa/dz0
   
    call bnd_b('xdim')
    call bnd_b('ydim')
    call bnd_b('zdim')

  endif

    call advectbx_y
   
    b(1,:,:,:) = b(1,:,:,:) - wa/dy0
    wa = cshift(wa,shift=-1,dim=2)
    b(1,:,:,:) = b(1,:,:,:) + wa/dy0
    b(2,:,:,:) = b(2,:,:,:) - wa/dx0
    wa = cshift(wa,shift=1,dim=1)
    b(2,:,:,:) = b(2,:,:,:) + wa/dx0

    call bnd_b('xdim')
    call bnd_b('ydim')
    if(dimensions .eq. '3d') call bnd_b('zdim')

ELSE IF(trim(magnetic_scheme) .eq. 'ssp' .and. magnetic_order .eq. 3) THEN

  if (dimensions .eq. '3d') then

! Now b1 is the initial magnetic field
    b1(:,:,:,:) = b(:,:,:,:)

! DIFFUSION FIRST STEP for bz_y

    call advectbz_y
    call diffusebz_y

    b(3,:,:,:) = b1(3,:,:,:) - wcu/dy0
    wcu = cshift(wcu,shift=1,dim=2)
    b(3,:,:,:) = b(3,:,:,:) + wcu/dy0
    wcu = cshift(wcu,shift=-1,dim=2)
    b(2,:,:,:) = b1(2,:,:,:) + wcu/dz0
    wcu = cshift(wcu,shift=1,dim=3)
    b(2,:,:,:) = b(2,:,:,:) - wcu/dz0

! FIRST STEP for bz_y

    b(3,:,:,:) = b(3,:,:,:) - wa/dy0
    wa = cshift(wa,shift=-1,dim=2)
    b(3,:,:,:) = b(3,:,:,:) + wa/dy0
    b(2,:,:,:) = b(2,:,:,:) - wa/dz0
    wa = cshift(wa,shift=1,dim=3)
    b(2,:,:,:) = b(2,:,:,:) + wa/dz0
   
    call bnd_b('xdim')
    call bnd_b('ydim')
    call bnd_b('zdim')

! DIFFUSION SECOND STEP for bz_y

    call advectbz_y
    call diffusebz_y

    b(3,:,:,:) = 0.75*b1(3,:,:,:) + 0.25*(b(3,:,:,:) - wcu/dy0)
    wcu = cshift(wcu,shift=1,dim=2)
    b(3,:,:,:) = b(3,:,:,:) + 0.25*wcu/dy0
    wcu = cshift(wcu,shift=-1,dim=2)
    b(2,:,:,:) = 0.75*b1(2,:,:,:) + 0.25*(b(2,:,:,:) + wcu/dz0)
    wcu = cshift(wcu,shift=1,dim=3)
    b(2,:,:,:) = b(2,:,:,:) - 0.25*wcu/dz0

! SECOND STEP for bz_y

    b(3,:,:,:) = b(3,:,:,:) - 0.25*wa/dy0
    wa = cshift(wa,shift=-1,dim=2)
    b(3,:,:,:) = b(3,:,:,:) + 0.25*wa/dy0
    b(2,:,:,:) = b(2,:,:,:) - 0.25*wa/dz0
    wa = cshift(wa,shift=1,dim=3)
    b(2,:,:,:) = b(2,:,:,:) + 0.25*wa/dz0
   
    call bnd_b('xdim')
    call bnd_b('ydim')
    call bnd_b('zdim')

! DIFFUSION THIRD STEP for bz_y

    call advectbz_y
    call diffusebz_y

    b(3,:,:,:) = 1./3.*b1(3,:,:,:) + 2./3.*(b(3,:,:,:) - wcu/dy0)
    wcu = cshift(wcu,shift=1,dim=2)
    b(3,:,:,:) = b(3,:,:,:) + 2./3.*wcu/dy0
    wcu = cshift(wcu,shift=-1,dim=2)
    b(2,:,:,:) = 1./3.*b1(2,:,:,:) + 2./3.*(b(2,:,:,:) + wcu/dz0)
    wcu = cshift(wcu,shift=1,dim=3)
    b(2,:,:,:) = b(2,:,:,:) - 2./3.*wcu/dz0

! THIRD STEP for bz_y

    b(3,:,:,:) = b(3,:,:,:) - 2./3.*wa/dy0
    wa = cshift(wa,shift=-1,dim=2)
    b(3,:,:,:) = b(3,:,:,:) + 2./3.*wa/dy0
    b(2,:,:,:) = b(2,:,:,:) - 2./3.*wa/dz0
    wa = cshift(wa,shift=1,dim=3)
    b(2,:,:,:) = b(2,:,:,:) + 2./3.*wa/dz0
   
    call bnd_b('xdim')
    call bnd_b('ydim')
    call bnd_b('zdim')

  endif

! Now b1 is the initial magnetic field 
    b1(:,:,:,:) = b(:,:,:,:)

! DIFFUSION FIRST STEP for bx_y

    call advectbx_y
    call diffusebx_y
    
    b(1,:,:,:) = b1(1,:,:,:) - wcu/dy0
    wcu = cshift(wcu,shift=1,dim=2)
    b(1,:,:,:) = b(1,:,:,:) + wcu/dy0
    wcu = cshift(wcu,shift=-1,dim=2)
    b(2,:,:,:) = b1(2,:,:,:) + wcu/dx0
    wcu = cshift(wcu,shift=1,dim=1)
    b(2,:,:,:) = b(2,:,:,:) - wcu/dx0

! FIRST STEP for bx_y

    b(1,:,:,:) = b(1,:,:,:) - wa/dy0
    wa = cshift(wa,shift=-1,dim=2)
    b(1,:,:,:) = b(1,:,:,:) + wa/dy0
    b(2,:,:,:) = b(2,:,:,:) - wa/dx0
    wa = cshift(wa,shift=1,dim=1)
    b(2,:,:,:) = b(2,:,:,:) + wa/dx0

    call bnd_b('xdim')
    call bnd_b('ydim')
    if(dimensions .eq. '3d') call bnd_b('zdim')

! DIFFUSION SECOND STEP for bx_y

    call advectbx_y
    call diffusebx_y

    b(1,:,:,:) = 0.75*b1(1,:,:,:) + 0.25*(b(1,:,:,:) - wcu/dy0)
    wcu = cshift(wcu,shift=1,dim=2)
    b(1,:,:,:) = b(1,:,:,:) + 0.25*wcu/dy0
    wcu = cshift(wcu,shift=-1,dim=2)
    b(2,:,:,:) = 0.75*b1(2,:,:,:) + 0.25*(b(2,:,:,:) + wcu/dx0)
    wcu = cshift(wcu,shift=1,dim=1)
    b(2,:,:,:) = b(2,:,:,:) - 0.25*wcu/dx0

! SECOND STEP for bx_y

    b(1,:,:,:) = b(1,:,:,:) - 0.25*wa/dy0
    wa = cshift(wa,shift=-1,dim=2)
    b(1,:,:,:) = b(1,:,:,:) + 0.25*wa/dy0
    b(2,:,:,:) = b(2,:,:,:) - 0.25*wa/dx0
    wa = cshift(wa,shift=1,dim=1)
    b(2,:,:,:) = b(2,:,:,:) + 0.25*wa/dx0

    call bnd_b('xdim')
    call bnd_b('ydim')
    if(dimensions .eq. '3d') call bnd_b('zdim')

! DIFFUSION THIRD STEP for bx_y

    call advectbx_y
    call diffusebx_y

    b(1,:,:,:) = 1./3.*b1(1,:,:,:) + 2./3.*(b(1,:,:,:) - wcu/dy0)
    wcu = cshift(wcu,shift=1,dim=2)
    b(1,:,:,:) = b(1,:,:,:) + 2./3.*wcu/dy0
    wcu = cshift(wcu,shift=-1,dim=2)
    b(2,:,:,:) = 1./3.*b1(2,:,:,:) + 2./3.*(b(2,:,:,:) + wcu/dx0)
    wcu = cshift(wcu,shift=1,dim=1)
    b(2,:,:,:) = b(2,:,:,:) - 2./3.*wcu/dx0

! THIRD STEP for bx_y

    b(1,:,:,:) = b(1,:,:,:) - 2./3.*wa/dy0
    wa = cshift(wa,shift=-1,dim=2)
    b(1,:,:,:) = b(1,:,:,:) + 2./3.*wa/dy0
    b(2,:,:,:) = b(2,:,:,:) - 2./3.*wa/dx0
    wa = cshift(wa,shift=1,dim=1)
    b(2,:,:,:) = b(2,:,:,:) + 2./3.*wa/dx0

    call bnd_b('xdim')
    call bnd_b('ydim')
    if(dimensions .eq. '3d') call bnd_b('zdim')

ELSE
900 format(' Scheme ',a4,' of order ',i2,' not implemented for magnetic field advection')   
    write(*,*) ' '
    write(*,900),trim(magnetic_scheme),magnetic_order
    write(*,*) ' '
    STOP
ENDIF

  end subroutine advectbzxy

!==========================================================================================

  subroutine tvdb(vibj,b,vg,n,dt,di)
    integer i,n,ip,ipp,im
    real dt
    real, dimension(n) :: vibj,b,vg,di
! locals
    real, dimension(n) :: b1,vibj1,vh
    real w,wp,wm,dw,v

  ! unlike the B field, the vibj lives on the right cell boundary

    vh=(vg+eoshift(vg,1,boundary=big))/2.

IF (trim(magnetic_scheme) .ne. 'ssp') then

    where(vh .gt. 0.)
      vibj1=b*vg
    elsewhere
      vibj1=eoshift(b*vg,1,boundary=big)
    end where
    b1=b-(vibj1-eoshift(vibj1,-1,boundary=big))*dt/di/2.

ELSE
    b1=b
ENDIF

    do i=3,n-3
      ip=i+1
      ipp=ip+1
      im=i-1
      v=vh(i)
      if (v .gt. 0.) then
        w=vg(i)*b1(i)
        wp=(vg(ip)*b1(ip)-w)/2.
        wm=(w-vg(im)*b1(im))/2.
      else
        w=vg(ip)*b1(ip)
        wp=(w-vg(ipp)*b1(ipp))/2.
        wm=(vg(i)*b1(i)-w)/2.
      end if
      dw=0.
      if(wm*wp .gt. 0.) dw=2.*wm*wp/(wm+wp)
      vibj(i)=(w+dw)*dt
    end do

  end subroutine tvdb

!pawlaszek

  subroutine advectby_x

    implicit none
    real, dimension(nx) :: vxby,by_x,vx
    integer                j,k,jm
    
    do k=1,nz
      do j=2,ny
        jm=j-1
        vx=(u(2,:,jm,k)+u(2,:,j,k))/(u(1,:,jm,k)+u(1,:,j,k))
        vx=(eoshift(vx,shift=-1,boundary=big) &
           +eoshift(vx,shift=1,boundary=big)+2.*vx)/4.
        by_x=b(2,:,j,k)
        call tvdb(vxby,by_x,vx,nx,dt,dx)
        wa(:,j,k) = vxby
      end do
    end do

    call bnd_emf(wa, 'vxby', 'xdim')
    call bnd_emf(wa, 'vxby', 'ydim')
    if(dimensions .eq. '3d') call bnd_emf(wa, 'vxby', 'zdim')


  end subroutine advectby_x

  subroutine advectbz_x

    implicit none
    real, dimension(nx) :: vxbz,bz_x,vx
    integer                j,k,km

    do k=2,nz
      km=k-1
      do j=1,ny
        vx=(u(2,:,j,km)+u(2,:,j,k))/(u(1,:,j,km)+u(1,:,j,k))
        vx=(eoshift(vx,-1,boundary=big) &
           +eoshift(vx,1,boundary=big)+2*vx)/4
        bz_x=b(3,:,j,k)
        call tvdb(vxbz,bz_x,vx,nx,dt,dx)
        wa(:,j,k) = vxbz
      end do
    end do

    call bnd_emf(wa, 'vxbz', 'xdim')
    call bnd_emf(wa, 'vxbz', 'ydim')
    if(dimensions .eq. '3d') call bnd_emf(wa, 'vxbz', 'zdim')

  end subroutine advectbz_x

  subroutine advectbx_z

    implicit none
    real, dimension(nz)  :: vzbx,bx_z,vz
    integer                 j,i,im

    do j=1,ny
      do i=2,nx
        im=i-1
        vz=(u(4,im,j,:)+u(4,i,j,:))/(u(1,im,j,:)+u(1,i,j,:))
        vz=(eoshift(vz,-1,boundary=big) &
           +eoshift(vz,1,boundary=big)+2*vz)/4
        bx_z=b(1,i,j,:)
        call tvdb(vzbx,bx_z,vz,nz,dt,dz)
        wa(i,j,:) = vzbx
      end do
    end do
    
    if(dimensions .eq. '3d') call bnd_emf(wa, 'vzbx', 'zdim')
    call bnd_emf(wa, 'vzbx', 'xdim')
    call bnd_emf(wa, 'vzbx', 'ydim')

  end subroutine advectbx_z

  subroutine advectby_z
    
    implicit none
    real, dimension(nz)  :: vzby,by_z,vz
    integer                 i,j,jm

    do j=2,ny
      jm=j-1
      do i=1,nx
        vz=(u(4,i,jm,:)+u(4,i,j,:))/(u(1,i,jm,:)+u(1,i,j,:))
        vz=(eoshift(vz,-1,boundary=big) &
           +eoshift(vz,1,boundary=big)+2*vz)/4
        by_z=b(2,i,j,:)
        call tvdb(vzby,by_z,vz,nz,dt,dz)
        wa(i,j,:) = vzby
      end do
    end do

     if(dimensions .eq. '3d') call bnd_emf(wa, 'vzby', 'zdim')
     call bnd_emf(wa, 'vzby', 'xdim')
     call bnd_emf(wa, 'vzby', 'ydim')

  end subroutine advectby_z

  subroutine advectbz_y

    implicit none
    real, dimension(ny)   :: vybz,bz_y,vy
    integer                  i,k,km

    do k=2,nz
      km=k-1
      do i=1,nx  
        vy=(u(3,i,:,km)+u(3,i,:,k))/(u(1,i,:,km)+u(1,i,:,k))
        vy=(eoshift(vy,-1,boundary=big) &
           +eoshift(vy,1,boundary=big)+2*vy)/4
        bz_y=b(3,i,:,k)
        call tvdb(vybz,bz_y,vy,ny,dt,dy)
        wa(i,:,k) = vybz
      end do
    end do

    call bnd_emf(wa, 'vybz', 'ydim')
    if(dimensions .eq. '3d') call bnd_emf(wa, 'vybz', 'zdim')
    call bnd_emf(wa, 'vybz', 'xdim')
  
  end subroutine advectbz_y

  subroutine advectbx_y

    implicit none
    real, dimension(ny) :: vybx,bx_y,vy
    integer                k,i,im

    do k=1,nz
      do i=2,nx
        im=i-1
        vy=(u(3,im,:,k)+u(3,i,:,k))/(u(1,im,:,k)+u(1,i,:,k))
        vy=(eoshift(vy,-1,boundary=big) &
           +eoshift(vy,1,boundary=big)+2*vy)/4
        bx_y=b(1,i,:,k)
        call tvdb(vybx,bx_y,vy,ny,dt,dy)
        wa(i,:,k) = vybx
      end do
    end do

    call bnd_emf(wa, 'vybx', 'ydim')
    if(dimensions .eq. '3d') call bnd_emf(wa, 'vybx', 'zdim')
    call bnd_emf(wa, 'vybx', 'xdim')

  end subroutine advectbx_y

!pawlaszek

!==========================================================================================

   subroutine relaxing_tvd(u,b,sweep,i1,i2,x,xl,xr,dx,n,nb,dt)
 ! Cooling and heating implemented following Rafal Kosinski

   use gravity

    implicit none
    integer i1,i2, i, n,nb
    real dt
    real, dimension(5,n) :: u,dx,cfr,ul,ur
    real, dimension(3,n) :: b
    real, dimension(n)   :: x,xl,xr,grav,gravl, gravr, drdx, dgdx 
    real, dimension(n)   :: velx, dvelx, flux_momx, dmomx
           
!    real coolheat_profile(n)
    character sweep*6
    integer nsubstep,substep
    real dt_substep, dt_coolheat_cell
          
!locals    
    real, dimension(5,n) :: w, u1,u2,u3,ul0,ur0,ul1,ur1,ul2,ur2,ul3,ur3, &
                            fr,fl,flux,dfrp,dfrm,dflm,dflp,dfl,dfr, &
                            dulf,durf,duls,durs
                            
                            
    real, dimension(n)   :: dens, dens1, ekin, eint, eint1, einth, emag, eint_src
    real, dimension(n)   :: temp, dtemp, deint, flux_heat
    real K_heat 
    logical bulk_viscosity1     
    
    bulk_viscosity1=.false.   

    w(:,:) = 0.0
    durs(:,:) = 0.0
    duls(:,:) = 0.0

!-------------------------------------------------------------------------------
!   RUNGE-KUTTA 1-ST STEP

    call mhdflux(w,cfr,u,b,n,nb)

    fr = (u*cfr+w)/2.
    fl=(u*cfr-w)/2.
    ur0 = fr/cfr
    ul0 = fl/cfr
    
    fl=cshift(fl,shift=1,dim=2)                 


    if(trim(integration_scheme) .eq. 'ssp') then

    if(bulk_viscosity) then
      velx = ur0(2,:)/ur0(1,:)
      dvelx = (cshift(velx,1) - velx)/dx(1,:)
      flux_momx = - nu_bulk * ur0(1,:) * dvelx
      fr(2,:)=fr(2,:) + flux_momx
    endif


      dfrp=(cshift(fr,shift=1,dim=2)-fr)*0.5    
      dfrm=cshift(dfrp,shift=-1,dim=2)          
      call flimiter(fr,dfrm,dfrp,n)
    
    if(bulk_viscosity) then
      velx = ul0(2,:)/ul0(1,:)
      dvelx = (cshift(velx,1) - velx)/dx(1,:)
      flux_momx = - nu_bulk * ul0(1,:) * dvelx
      fl(2,:)=fl(2,:) - flux_momx
    endif
    
      dflp=(fl-cshift(fl,shift=1,dim=2))*0.5    
      dflm=cshift(dflp,shift=-1,dim=2)          
      call flimiter(fl,dflm,dflp,n)

    endif 
    
    durf = (fr-cshift(fr,shift=-1,dim=2))/dx*dt 
    dulf = (fl-cshift(fl,shift=-1,dim=2))/dx*dt
    
    if(integration_order .eq. 2 ) then   
      ur1= ur0 - 0.5*durf               
      ul1= ul0 + 0.5*dulf       
    else            
      ur1= ur0 - durf           
      ul1= ul0 + dulf       
    endif          

!
    ur1(1,:) = max(ur1(1,:), smalld/2.)
    ul1(1,:) = max(ul1(1,:), smalld/2.)
    
! Gravity source terms -------------------------------------
    if(gravaccel) then

      call grav_pot2accel(sweep,i1,i2, n, gravr)
        
      gravl      = cshift(gravr,-1) 
      gravl(1)   = gravl(2)
      gravr(n)   = gravr(n-1)
 
      if(trim(integration_scheme) .eq. 'orig' .and. integration_order .eq. 2 ) then   
        duls(5,:) = gravr*(ul1(2,:))*dt 
        durs(5,:) = gravl*(ur1(2,:))*dt 
        duls(2,:) = gravr*(ul1(1,:))*dt 
        durs(2,:) = gravl*(ur1(1,:))*dt 
        ur1= ur1 + 0.5*durs             
        ul1= ul1 + 0.5*duls     
      else            
        duls(5,:) = gravr*(ul0(2,:))*dt 
        durs(5,:) = gravl*(ur0(2,:))*dt 
        duls(2,:) = gravr*(ul0(1,:))*dt 
        durs(2,:) = gravl*(ur0(1,:))*dt 
        ur1= ur1 + durs                 
        ul1= ul1 + duls             
      endif        
         
    endif
! ----------------------------------------------------------

    u1 = ul1 + ur1
    u1(1,:) = max(u1(1,:), smalld)
          
    dens = u1(1,:)
    ekin = 0.5*(u1(2,:)*u1(2,:) + u1(3,:)*u1(3,:) + u1(4,:)*u1(4,:))/u1(1,:)
    emag = 0.5*(b(1,:)*b(1,:) + b(2,:)*b(2,:) + b(3,:)*b(3,:)) 
    eint = u1(5,:)-ekin-emag
    eint = max(eint,smallei)
    u1(5,:) = eint+ekin+emag
    
    if(integration_order .eq. 1) then
      u(:,:) = u1(:,:)
      goto 999
    endif  
    
!-------------------------------------------------------------------------------
!   RUNGE-KUTTA 2-ND STEP

    call mhdflux(w,cfr,u1,b,n,nb)

    fr = (u1*cfr+w)/2.
    ur = fr/cfr

    if(bulk_viscosity) then
      velx = ur(2,:)/ur(1,:)
      dvelx = (cshift(velx,1) - velx)/dx(1,:)
      flux_momx = - nu_bulk * ur(1,:) * dvelx
      fr(2,:)=fr(2,:) + flux_momx
    endif

    dfrp=(cshift(fr,shift=1,dim=2)-fr)*0.5      
    dfrm=cshift(dfrp,shift=-1,dim=2)            
    call flimiter(fr,dfrm,dfrp,n)
    durf = (fr-cshift(fr,shift=-1,dim=2))/dx*dt 
    
    fl=(u1*cfr-w)/2.
    ul = fl/cfr
    fl=cshift(fl,shift=1,dim=2)                 

    if(bulk_viscosity) then
      velx = ul(2,:)/ul(1,:)
      dvelx = (cshift(velx,1) - velx)/dx(1,:)
      flux_momx = - nu_bulk * ul(1,:) * dvelx
      fl(2,:)=fl(2,:) - flux_momx
    endif
    
    dflp=(fl-cshift(fl,shift=1,dim=2))*0.5      
    dflm=cshift(dflp,shift=-1,dim=2)            
    call flimiter(fl,dflm,dflp,n)
    dulf = (fl-cshift(fl,shift=-1,dim=2))/dx*dt

    if(trim(integration_scheme) .eq. 'orig' .and. integration_order .eq. 2 ) then   
      ur2= ur0 - durf           
      ul2= ul0 + dulf   
    else            
      ur2 = 0.75*ur0 + 0.25*(ur1 - durf)        !
      ul2 = 0.75*ul0 + 0.25*(ul1 + dulf)        !
    endif          

    ur2(1,:) = max(ur2(1,:), smalld/2.)
    ul2(1,:) = max(ul2(1,:), smalld/2.)

! Gravity source terms (grav same as in the RK half-timestep)
    if(gravaccel) then

      if(trim(integration_scheme) .eq. 'orig' .and. integration_order .eq. 2 ) then   
        duls(5,:) = gravr*(ul0(2,:)+ul2(2,:))*dt/2. 
        durs(5,:) = gravl*(ur0(2,:)+ur2(2,:))*dt/2. 
        duls(2,:) = gravr*(ul0(1,:)+ul2(1,:))*dt/2. 
        durs(2,:) = gravl*(ur0(1,:)+ur2(1,:))*dt/2. 
        ur2= ur2 + durs         
        ul2= ul2 + duls 
      else            
        duls(5,:) = gravr*(ul1(2,:))*dt 
        durs(5,:) = gravl*(ur1(2,:))*dt 
        duls(2,:) = gravr*(ul1(1,:))*dt 
        durs(2,:) = gravl*(ur1(1,:))*dt 
        ur2= ur2 + 0.25 * durs          
        ul2= ul2 + 0.25 * duls      
      endif        

    endif

    u2 = ul2 + ur2
    dens = u2(1,:)
    ekin = 0.5*(u2(2,:)*u2(2,:)+u2(3,:)*u2(3,:)+u2(4,:)*u2(4,:))/u2(1,:)
    eint = u(5,:)-ekin-emag
    eint = max(eint,smallei)
    u(5,:) = eint+ekin+emag

    if(integration_order .eq. 2) then
      u(:,:) = u2(:,:)
      goto 999
    endif  

!   RUNGE-KUTTA 3-RD STEP

    call mhdflux(w,cfr,u2,b,n,nb)

    fr = (u2*cfr+w)/2.
    ur = fr/cfr

    if(bulk_viscosity) then
      velx = ur(2,:)/ur(1,:)
      dvelx = (cshift(velx,1) - velx)/dx(1,:)
      flux_momx = - nu_bulk * ur(1,:) * dvelx
      fr(2,:)=fr(2,:) + flux_momx
    endif

    dfrp=(cshift(fr,shift=1,dim=2)-fr)*0.5      
    dfrm=cshift(dfrp,shift=-1,dim=2)            
    call flimiter(fr,dfrm,dfrp,n)
    durf = (fr-cshift(fr,shift=-1,dim=2))/dx*dt 


    fl = (u2*cfr-w)/2.
    ul = fl/cfr
    
    fl=cshift(fl,shift=1,dim=2)                 

    if(bulk_viscosity) then
      velx = ul(2,:)/ul(1,:)
      dvelx = (cshift(velx,1) - velx)/dx(1,:)
      flux_momx = - nu_bulk * ul(1,:) * dvelx
      fl(2,:)=fl(2,:) - flux_momx
    endif
    
    
    
    dflp=(fl-cshift(fl,shift=1,dim=2))*0.5      
    dflm=cshift(dflp,shift=-1,dim=2)            
    call flimiter(fl,dflm,dflp,n)
    dulf = (fl-cshift(fl,shift=-1,dim=2))/dx*dt


!-----------------------------------------------------------------------!
    ur3 = 1./3.*ur0 + 2./3.*(ur2 - durf)        !
    ul3 = 1./3.*ul0 + 2./3.*(ul2 + dulf)        !
!-----------------------------------------------------------------------!


    ur3(1,:) = max(ur3(1,:), smalld/2.)
    ul3(1,:) = max(ul3(1,:), smalld/2.)

! Gravity source terms 
    if(gravaccel) then
    
      duls(5,:) = gravr*(ul2(2,:))*dt 
      durs(5,:) = gravl*(ur2(2,:))*dt 
      duls(2,:) = gravr*(ul2(1,:))*dt 
      durs(2,:) = gravl*(ur2(1,:))*dt 

      ul3 = ul3 + 2./3. * duls
      ur3 = ur3 + 2./3. * durs

    endif

    u3 = ul3 + ur3

    dens = u3(1,:)
    ekin = 0.5*(u3(2,:)*u3(2,:) + u3(3,:)*u3(3,:) + u3(4,:)*u3(4,:))/u3(1,:)
    eint = u3(5,:)-ekin-emag
    eint = max(eint,smallei)
    u3(5,:) = eint+ekin+emag

    u = u3


999 continue


    return
  end subroutine relaxing_tvd


!==========================================================================================
 
  subroutine mhdflux(v,cfr,u,b,n,nb)
    implicit none
    integer n,nb
    real, dimension(5,n)::v,u,cfr
    real, dimension(3,n):: b(3,n)
! locals
    real, dimension(n) :: vx,ps,p,pmag

    v   = 0.0
    cfr = 0.0

    pmag(2:n-1)=(b(1,2:n-1)*b(1,2:n-1)+b(2,2:n-1)*b(2,2:n-1)+b(3,2:n-1)*b(3,2:n-1))*0.5
    vx(2:n-1)=u(2,2:n-1)/u(1,2:n-1)
    ps(2:n-1)=(u(5,2:n-1)-(u(2,2:n-1)*u(2,2:n-1)+u(3,2:n-1)*u(3,2:n-1) & 
              +u(4,2:n-1)*u(4,2:n-1))/u(1,2:n-1)*0.5)*(gamma-1)+(2-gamma)*pmag(2:n-1)
    v(1,2:n-1)=u(2,2:n-1)
    v(2,2:n-1)=u(2,2:n-1)*vx(2:n-1)+ps(2:n-1)-b(1,2:n-1)*b(1,2:n-1)
    v(3,2:n-1)=u(3,2:n-1)*vx(2:n-1)-b(2,2:n-1)*b(1,2:n-1)
    v(4,2:n-1)=u(4,2:n-1)*vx(2:n-1)-b(3,2:n-1)*b(1,2:n-1)
    v(5,2:n-1)=(u(5,2:n-1)+ps(2:n-1))*vx(2:n-1)-b(1,2:n-1)*(b(1,2:n-1)*u(2,2:n-1) &
                +b(2,2:n-1)*u(3,2:n-1)+b(3,2:n-1)*u(4,2:n-1))/u(1,2:n-1)
    p(2:n-1)=ps(2:n-1)-pmag(2:n-1)

    select case (freezing_speed) 
      case ('local')
!       The freezing speed is now computed locally (in each cell) 
!       as in Trac & Pen (2003). This ensures much sharper shocks, 
!       but sometimes may lead to numerical instabilities      
        cfr(1,2:n-1) = abs(vx(2:n-1)) &
                      +max(sqrt( abs(2*pmag(2:n-1) + gamma*p(2:n-1))/u(1,2:n-1)),small)
        cfr(1,1) = cfr(1,2)
        cfr(1,n) = cfr(1,n-1)   
        cfr = spread(cfr(1,:),1,5)
      case ('global')
!       The freezing speed is now computed globally 
!       (c=const for the whole domain) in sobroutine 'timestep' 
!       Original computation of the freezing speed was done
!       for each sweep separately:
!       c=maxval(abs(vx(nb-2:n-1))+sqrt(abs(2*pmag(nb-2:n-1) &
!                            +gamma*p(nb-2:n-1))/u(1,nb-2:n-1)))  
        cfr(:,:) = c
      case default
        write(*,*) 'Freezing speed: ', freezing_speed, ' undefined'
        stop  
    end select


  end subroutine mhdflux

!==========================================================================================

  subroutine flimiter(f,a,b,n)
    implicit none
    integer n
    real, dimension(5,n) :: f,a,b,c
    
    select case(trim(flux_limiter))
      case('vanleer')
        c = a*b
        where (c .gt. 0)                                        
          f=f+2*c/(a+b)
        endwhere      
      case('minmod')
        f = f+(sign(1.,a)+sign(1.,b))*min(abs(a),abs(b))/2.
      case('superbee')
        where (abs(a) .gt. abs(b))
          f = f+(sign(1.,a)+sign(1.,b))*min(abs(a), abs(2*b))/2.
        elsewhere
          f = f+(sign(1.,a)+sign(1.,b))*min(abs(2*a), abs(b))/2.
        endwhere
      case default
        write(*,*) 'Flux limiter: ', flux_limiter,' not defined.'
        stop
    end select
    return
  end subroutine flimiter


!==========================================================================================

end module mhdblock









