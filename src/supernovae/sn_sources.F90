#include "mhd.def"
! Written by M. Hanasz, 2003, 2007 & K. Kowalik 2007 
! Adapted for this code by M. Hanasz, November 2007


module sn_sources


  use arrays
  use constants
  use grid
  use start
#ifdef SHEAR  
  use shear
#endif SHEAR  
  
  implicit none
  real xsn,ysn,zsn
  real epsi,epso
  real ysna, ysni, ysno
  real             :: phi,theta
  integer, save :: nsn, nsn_last
  real, save    :: dt_sn_prev, ecr_supl, decr_supl   
  
  real gasdev
  real,    save :: gset             
  integer, save :: irand, iset
  

 contains


 
  subroutine random_sn
! Written by: M. Hanasz  

    implicit none 
    real dt_sn, t_dw1
    integer isn, nsn_per_timestep
      
    dt_sn = 1./(f_sn+small)  
      
    nsn = t/dt_sn
    nsn_per_timestep = nsn - nsn_last
    nsn_last = nsn
    
    do isn = 0, nsn_per_timestep

      call rand_coords

#ifdef COSM_RAYS
      call cr_sn
#endif /* COSM_RAYS */

#ifdef DIPOLS
      call rand_angles !(phi, theta)
      call dipol_sn
#endif DIPOLS

    enddo ! isn

    return

  end subroutine random_sn

!--------------------------------------------------------------------------
 
  subroutine cr_sn
! Written by: M. Hanasz  
    implicit none
    integer i,j,k, ipm, jpm   
    real decr   

      do k=1,nz
        do j=1,ny
          do i=1,nx
 
            do ipm=-1,1

              if(ipm .eq. -1) ysna = ysno
              if(ipm .eq.  0) ysna = ysn
              if(ipm .eq.  1) ysna = ysni

              do jpm=-1,1


                decr = amp_ecr_sn * ethu  &
                       * EXP(-((x(i)-xsn+real(ipm)*Lx)**2  &
                       + (y(j)-ysna+real(jpm)*Ly)**2  &
                       + (z(k)-zsn)**2)/r_sn**2)  

                u(iecr,i,j,k) = u(iecr,i,j,k) + decr		
  
              enddo ! jpm
            enddo ! ipm
   
          enddo ! i
        enddo ! j
      enddo ! k

         
    return

  end subroutine cr_sn

!--------------------------------------------------------------------------

  subroutine dipol_sn !(amp,rmk)
! Writen by: K.Kowalik
  
    use constants
    use arrays, only : u,b,xl,yl,zl,ibx,iby,ibz,&
                        iena,imxa,imya,imza,idna
    use grid, only   : dx,dy,dz
    use start, only  : smallei
    implicit none
!    real, intent(in) :: phi,theta
!    real             :: phi,theta
!    real, intent(in) :: amp_dip_sn,rmk
    real             :: rmk
    real             :: temp1,temp2
    integer          :: ipm, jpm
    real, dimension(:,:,:,:), allocatable :: A
#ifndef ISO
    real, dimension(:,:,:), allocatable :: ekin,eint
#endif ISO


    real :: xx, yy, zz, x, y, z, r, rc, sint
    real :: Aphi
    integer :: nx,ny,nz,i,j,k

    rmk = r_sn
    
    nx = size(b,2)
    ny = size(b,3)
    nz = size(b,4)

    if(.not.allocated(A)) then
        allocate(A(3,nx,ny,nz))
    else
      write(*,*) 'erroneous array allocation in dipol'
      stop
    endif

#ifndef ISO
    if(.not.allocated(ekin).and. .not.allocated(eint)) then
        allocate(eint(nx,ny,nz))
        allocate(ekin(nx,ny,nz))
    else
      write(*,*) 'erroneous array allocation in dipol'
      stop
    endif

    ekin(:,:,:) = 0.5*( u(imxa,:,:,:)**2 + u(imya,:,:,:)**2 + u(imza,:,:,:)**2 ) / u(idna,:,:,:)
    eint(:,:,:) = max( u(iena,:,:,:) - ekin -  0.5*( b(ibx,:,:,:)**2 + b(iby,:,:,:)**2 + &
            b(ibz,:,:,:)**2) , smallei)
#endif ISO

!    call rand_angles(phi, theta)
     
    do i = 1,nx
      xx = xl(i)
      do j = 1,ny
        yy = yl(j)
        do k = 1,nz 
          zz = zl(k)
	  
	    A(:,i,j,k) = 0.0

            do ipm=-1,1

              if(ipm .eq. -1) ysna = ysno
              if(ipm .eq.  0) ysna = ysn
              if(ipm .eq.  1) ysna = ysni

              do jpm=-1,1

                x    = ((xx-xsn+real(ipm)*Lx)*cos(phi)+(yy-ysna+real(jpm)*Ly)*sin(phi)) &
	                       *cos(theta)-(zz-zsn)*sin(theta) 
                y    = ((yy-ysna+real(jpm)*Ly)*cos(phi)-(xx-xsn+real(ipm)*Lx)*sin(phi))
                z    = ((xx-xsn+real(ipm)*Lx)*cos(phi)+(yy-ysna+real(jpm)*Ly)*sin(phi)) &
	                       *sin(theta)+(zz-zsn)*cos(theta) 
                  
                r    = sqrt(x**2 + y**2 + z**2)
                rc   = sqrt(x**2 + y**2)
                sint = rc/(r+small)

                Aphi =  amp_dip_sn * r*sint / (rmk**2 + r**2 + 2.*rmk*r*sint)**1.5
                temp1 = -1.0 *Aphi* y / (rc+small)
                temp2 = Aphi * x / (rc+small)
	   
                A(1,i,j,k) = A(1,i,j,k) + temp1*cos(theta)*cos(phi) - temp2*sin(phi)
                A(2,i,j,k) = A(2,i,j,k) + temp1*cos(theta)*sin(phi) + temp2*cos(phi)
                A(3,i,j,k) = A(3,i,j,k) - temp1*sin(theta)

              enddo ! jpm
            enddo ! ipm


        enddo
      enddo
    enddo
     

    b(ibx,1:nx,  1:ny-1,1:nz-1) = b(ibx,1:nx,  1:ny-1,1:nz-1) + &
           (A(3,1:nx,2:ny,1:nz-1) - A(3,1:nx,1:ny-1,1:nz-1))/dy - &
           (A(2,1:nx,1:ny,2:nz  ) - A(2,1:nx,1:ny,  1:nz-1))/dz

    b(iby,1:nx-1,1:ny,  1:nz-1) = b(iby,1:nx-1,1:ny,  1:nz-1) + &
           (A(1,1:nx-1,1:ny,2:nz) - A(1,1:nx-1,1:ny,1:nz-1))/dz - &
           (A(3,2:nx,1:ny,1:nz-1) - A(3,1:nx-1,1:ny,1:nz-1))/dx

    b(ibz,1:nx-1,1:ny-1,1:nz  ) = b(ibz,1:nx-1,1:ny-1,1:nz  ) + & 
           (A(2,2:nx,1:ny-1,1:nz) - A(2,1:nx-1,1:ny-1,1:nz))/dx - &
           (A(1,1:nx-1,2:ny,1:nz) - A(1,1:nx-1,1:ny-1,1:nz))/dy
       
!     write(*,*) maxval(b(1,:,:,:))
!     write(*,*) maxval(b(2,:,:,:))
!     write(*,*) maxval(b(3,:,:,:))

!     stop 

#ifndef ISO
    u(iena,:,:,:) = eint + ekin + 0.5*(b(1,:,:,:)**2 + &
             b(2,:,:,:)**2 + b(3,:,:,:)**2)
    if (allocated(eint)) deallocate(eint)
    if (allocated(ekin)) deallocate(ekin)
#endif
    if (allocated(A)) deallocate(A)

  end subroutine dipol_sn


!--------------------------------------------------------------------------
 
  subroutine rand_coords
! Written by M. Hanasz

    real rand(4), znorm
    integer jsn,jremap 
    real dysn 

    call random_number(rand)
    xsn = xmin+ Lx*rand(1)
    ysn = ymin+ Ly*rand(2)

    if(dimensions == '3d') then
      irand = irand+4  
      znorm = gasdev(rand(3), rand(4)) 
      zsn = h_sn*znorm
    else
      zsn = 0.0
    endif

#ifdef SHEAR
      jsn  = js+int((ysn-ymin)/dy)
      dysn  = dmod(ysn,dy)
       
      epsi   = eps*dy 
      epso   = -epsi

!  outer boundary
      jremap = jsn - delj 
      jremap = mod(mod(jremap, nyd)+nyd,nyd)
      if (jremap .le. (js-1)) jremap = jremap + nyd

      ysno = y(jremap) + epso + dysn
 
!  inner boundary
      jremap = jsn + delj 
      jremap = mod(jremap, nyd)+nyd
      if (jremap .ge. (je+1)) jremap = jremap - nyd

      ysni = y(jremap) + epsi + dysn
#else SHEAR
      ysno = ysn
      ysni = ysn
#endif SHEAR
  
  
  end subroutine rand_coords

!-----------------------------------------------------------------------


  subroutine rand_angles !(phi,theta)
! Written by K. Kowalik
  
    use constants, only : pi
    implicit none
!    real, intent(out) :: phi, theta
!    real :: phi, theta
    real :: rand(2),rnz,rny,rnx

    call random_number(rand)
    rnz   = (1.0-2.0*rand(1))
    rnx = sqrt(1.0-rnz**2)*cos(2.*pi*rand(2))
    rny = sqrt(1.0-rnz**2)*sin(2.*pi*rand(2))
    phi = 2.0*pi*rand(2)
    theta = acos(rnz/sqrt(rnz**2+rny**2+rnx**2) )
    return
  end subroutine rand_angles

!=======================================================================
!
!      \\\\\\\         B E G I N   S U B R O U T I N E S        ///////
!      ///////    F R O M   N U M E R I C A L   R E C I P E S   \\\\\\\
!
!=======================================================================


      function gasdev(x,y)

      integer idum
      real x, y, x1, y1,  r
      real gasdev, rand(2)
      real fac,rsq
 
      if (iset.eq.0) then
1       x1=2.*x-1.
        y1=2.*y-1.
        r=x1**2+y1**2
        if(r.ge.1.) then
	  call random_number(rand)
          x = rand(1)
          y = rand(2)
          irand = irand+2
          go to 1
        endif
        fac=sqrt(-2.*log(r)/r)
        gset=x1*fac
        gasdev=y1*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      end function gasdev

!=======================================================================
!
!      \\\\\\\          E N D   S U B R O U T I N E S           ///////
!      ///////    F R O M   N U M E R I C A L   R E C I P E S   \\\\\\\
!
!=======================================================================
  

end module sn_sources
