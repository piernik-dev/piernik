#include "mhd.def"
module func
contains
  function tvdb_emf(vh,vg,b,dt)
    implicit none
    real, dimension(:), intent(in) :: vh,vg,b
    real, intent(in)               :: dt
    real, dimension(size(vh)) :: tvdb_emf
    integer :: i,ip,ipp,im
    real    :: w,wp,wm,dw,v

    tvdb_emf = 0.0

    do i = lbound(vh,1)+2, ubound(vh,1)-3
       ip  = i  + 1
       ipp = ip + 1
       im  = i  - 1
       v   = vh(i)
       if (v .gt. 0.) then
         w=vg(i)*b(i)
         wp=(vg(ip)*b(ip)-w)*0.5
         wm=(w-vg(im)*b(im))*0.5
       else
         w=vg(ip)*b(ip)
         wp=(w-vg(ipp)*b(ipp))*0.5
         wm=(vg(i)*b(i)-w)*0.5
       end if
       dw=0.0
       if(wm*wp > 0.0) dw=2.0*wm*wp/(wm+wp)
       tvdb_emf(i)=(w+dw)*dt
    enddo
    return
  end function tvdb_emf

   subroutine dipol(xi,yj,zk,phi,theta,amp,rmk)
     use arrays, only : u,b,xl,yl,zl,ibx,iby,ibz,&
                        iena,imxa,imya,imza,idna
     use grid, only   : dx,dy,dz
     use start, only  : smallei
     implicit none
     real, intent(in) :: phi,theta,amp,rmk
     real, intent(in) :: xi,yj,zk
     real  :: temp1,temp2
     real, dimension(:,:,:,:), allocatable :: A
#ifndef ISO
     real, dimension(:,:,:), allocatable :: ekin,eint
#endif ISO


     real :: xx, yy, zz, x, y, z, r, rc, sint
     real :: Aphi
     integer :: nx,ny,nz,i,j,k

     nx = size(u,1)
     ny = size(u,2)
     nz = size(u,3)

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
            
#endif
     do i = 1,nx
       xx = xl(i)
       do j = 1,ny
         yy = yl(j)
         do k = 1,nz 
           zz = zl(k)

           x    = ((xx-xi)*cos(phi)+(yy-yj)*sin(phi))*cos(theta)-(zz-zk)*sin(theta) 
           y    = ((yy-yj)*cos(phi)-(xx-xi)*sin(phi))
           z    = ((xx-xi)*cos(phi)+(yy-yj)*sin(phi))*sin(theta)+(zz-zk)*cos(theta) 
                  
           r    = sqrt(x**2 + y**2 + z**2)
           rc   = sqrt(x**2 + y**2)
           sint = rc/r

           Aphi =  amp * r*sint / (rmk**2 + r**2 + 2.*rmk*r*sint)**1.5
           temp1 = -1.0 *Aphi* y / rc
           temp2 = Aphi * x / rc

           A(1,i,j,k) =  temp1*cos(theta)*cos(phi) - temp2*sin(phi)
           A(2,i,j,k) =  temp1*cos(theta)*sin(phi) + temp2*cos(phi)
           A(3,i,j,k) = -temp1*sin(theta)
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
       
#ifndef ISO
     u(iena,:,:,:) = eint + ekin + 0.5*(b(1,:,:,:)**2 + &
             b(2,:,:,:)**2 + b(3,:,:,:)**2)
     if (allocated(eint)) deallocate(eint)
     if (allocated(ekin)) deallocate(ekin)
#endif
     if (allocated(A)) deallocate(A)

   end subroutine dipol

   subroutine rn_angles(phi,theta)
     use constants, only : pi
     implicit none
     real, intent(out) :: phi, theta
     real :: rn,rnz,rny,rnx

     call random_number(rn)
     rnz   = (1.0-2.0*rn)
     call random_number(rn)
     rnx = sqrt(1.0-rnz**2)*cos(2.*pi*rn)
     rny = sqrt(1.0-rnz**2)*sin(2.*pi*rn)
     phi = 2.0*pi*rn
     theta = acos(rnz/sqrt(rnz**2+rny**2+rnx**2) )
     return
   end subroutine rn_angles

end module func
