! $Id$
#include "piernik.def"
! Written by M. Hanasz, 2003, 2007 & K. Kowalik 2007 
! Adapted for this code by M. Hanasz, November 2007

module sn_sources


#ifdef SHEAR  
   use shear
#endif /* SHEAR */
  
   implicit none
   real epsi, epso
   real ysna, ysni, ysno
   integer, save :: nsn, nsn_last
   real, save    :: dt_sn_prev, ecr_supl, decr_supl
   real,    save :: gset
   integer, save :: irand, iset
  

   contains

   subroutine random_sn
! Written by: M. Hanasz  
      use start, only : f_sn,t
      use constants, only : small

      implicit none 
      real :: dt_sn
      real, dimension(2) :: orient
      real, dimension(3) :: snpos
      integer :: isn, nsn_per_timestep
      
      dt_sn = 1./(f_sn+small)  
        
      nsn = t/dt_sn
      nsn_per_timestep = nsn - nsn_last
      nsn_last = nsn
    
      do isn = 0, nsn_per_timestep

         call rand_coords(snpos)

#ifdef COSM_RAYS
         call cr_sn(snpos)
#endif /* COSM_RAYS */

#ifdef DIPOLS
!        call magn_multipole_sn(rand_angles(),snpos,)
#endif /* DIPOLS */
      enddo ! isn
      return
   end subroutine random_sn

!--------------------------------------------------------------------------
 
   subroutine cr_sn(pos)
! Written by: M. Hanasz 
      use arrays, only : nx,ny,nz,iecr,u,x,y,z
      use start,  only : amp_ecr_sn,r_sn,ethu
      use grid,   only : Lx,Ly
      implicit none
      real, dimension(3), intent(in) :: pos
      integer i,j,k, ipm, jpm   
      real decr, xsn,ysn,zsn
      xsn = pos(1)
      ysn = pos(2)
      zsn = pos(3)

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

   subroutine magn_multipole_sn(orient,pos,A)
      use start, only  : amp_dip_sn,howmulti,r_sn
      use grid, only   : Lx,Ly
      use constants
      use arrays, only : u,b,xl,yl,zl,ibx,iby,ibz,&
                          iena,imxa,imya,imza,idna
      use grid, only   : dx,dy,dz
      use start, only  : smallei

      implicit none
      real, intent(in), dimension(2) :: orient
      real, intent(in), dimension(3) :: pos
      real             :: phi,theta,xsn,ysn,zsn
      real             :: rmk
      real             :: temp1,temp2
      integer          :: ipm, jpm, kpm
      real, dimension(:,:,:,:) :: A
!     real, dimension(:,:,:,:), allocatable :: A
#ifndef ISO
!     real, dimension(:,:,:), allocatable :: ekin,eint
#endif /* ISO */
      real :: xx, yy, zz, x, y, z, r, rc, sint
      real :: Aphi, r_aux
      real :: sin_theta, cos_theta, sin_phi, cos_phi
      integer :: nx,ny,nz,i,j,k
      integer, dimension(2), parameter :: amp_fac = (/1,-1/)

      rmk = r_sn

      xsn = pos(1)
      ysn = pos(2)
      zsn = pos(3)

      phi   = orient(1)
      theta = orient(2)

      sin_theta = sin(theta)
      cos_theta = cos(theta)
      sin_phi = sin(phi)
      cos_phi = cos(phi)

      nx = size(b,2)
      ny = size(b,3)
      nz = size(b,4)

#ifndef ISO
!      if(fl == -1) then
!         if(.not.allocated(ekin).and. .not.allocated(eint)) then
!            allocate(eint(nx,ny,nz))
!            allocate(ekin(nx,ny,nz))
!         else
!            write(*,*) 'erroneous array allocation in dipol'
!            stop
!         endif

!
!         ekin(:,:,:) = 0.5*( u(imxa,:,:,:)**2 + u(imya,:,:,:)**2 + u(imza,:,:,:)**2 ) / u(idna,:,:,:)
!         eint(:,:,:) = max( u(iena,:,:,:) - ekin -  0.5*( b(ibx,:,:,:)**2 + b(iby,:,:,:)**2 + &
!               b(ibz,:,:,:)**2) , smallei)
!      endif
#endif /* ISO */

      do i = 1,nx
         xx = xl(i)
         do j = 1,ny
            yy = yl(j)
            do k = 1,nz 
               zz = zl(k)
 
!              A(:,i,j,k) = 0.0
             !>
             !! \todo Following lines are also needed for strict periodic domain,
             !! and corner-periodic boundaries. It's only temporary solution!!!
             !<
#ifdef SHEAR
               do ipm=-1,1

                  if(ipm .eq. -1) ysna = ysno
                  if(ipm .eq.  0) ysna = ysn
                  if(ipm .eq.  1) ysna = ysni

                  do jpm=-1,1
                     x    = ((xx-xsn+real(ipm)*Lx)*cos_phi+(yy-ysna+real(jpm)*Ly)*sin_phi) &
                             *cos_theta-(zz-zsn)*sin_theta 
                     y    = ((yy-ysna+real(jpm)*Ly)*cos_phi-(xx-xsn+real(ipm)*Lx)*sin_phi)
                     z    = ((xx-xsn+real(ipm)*Lx)*cos_phi+(yy-ysna+real(jpm)*Ly)*sin_phi) &
                             *sin_theta+(zz-zsn)*cos_theta 
                
#else
!                    ipm = 0; jpm=0; ysna = ysn
                     x    = ((xx-xsn)*cos_phi+(yy-ysn)*sin_phi) * cos_theta-(zz-zsn)*sin_theta 
                     y    = ((yy-ysn)*cos_phi-(xx-xsn)*sin_phi)
                     z    = ((xx-xsn)*cos_phi+(yy-ysn)*sin_phi) * sin_theta+(zz-zsn)*cos_theta          
#endif /* SHEAR */
                        do kpm=1,howmulti
                
                           r    = sqrt(x**2 + y**2 + (z+real(2*kpm-howmulti-1)*rmk)**2)
                           rc   = sqrt(x**2 + y**2)
                           sint = rc/(r+small)

                           r_aux = 1.0 / (rmk**2 + r**2 + 2.*rmk*r*sint)
                           r_aux = r_aux*r_aux*r_aux
                           r_aux = sqrt(r_aux)

!                          Aphi =  real(amp_fac(kpm)) * amp_dip_sn * r*sint / (rmk**2 + r**2 + 2.*rmk*r*sint)**1.5
                           Aphi =  real(amp_fac(kpm)) * amp_dip_sn * r*sint * r_aux
                           temp1 = -1.0 *Aphi* y / (rc+small)
                           temp2 = Aphi * x / (rc+small)
  
                           A(1,i,j,k) = A(1,i,j,k) + temp1*cos_theta*cos_phi - temp2*sin_phi
                           A(2,i,j,k) = A(2,i,j,k) + temp1*cos_theta*sin_phi + temp2*cos_phi
                           A(3,i,j,k) = A(3,i,j,k) - temp1*sin_theta

                        enddo ! kpm
#ifdef SHEAR
                  enddo ! jpm
               enddo ! ipm
#endif /* SHEAR */
            enddo
         enddo
      enddo
     
!      if(fl == 1) then
!         b(ibx,1:nx,  1:ny-1,1:nz-1) = b(ibx,1:nx,  1:ny-1,1:nz-1) + &
!            (A(3,1:nx,2:ny,1:nz-1) - A(3,1:nx,1:ny-1,1:nz-1))/dy - &
!            (A(2,1:nx,1:ny-1,2:nz) - A(2,1:nx,1:ny-1,1:nz-1))/dz
!
!         b(iby,1:nx-1,1:ny,  1:nz-1) = b(iby,1:nx-1,1:ny,  1:nz-1) + &
!            (A(1,1:nx-1,1:ny,2:nz) - A(1,1:nx-1,1:ny,1:nz-1))/dz - &
!            (A(3,2:nx,1:ny,1:nz-1) - A(3,1:nx-1,1:ny,1:nz-1))/dx
!
!         b(ibz,1:nx-1,1:ny-1,1:nz  ) = b(ibz,1:nx-1,1:ny-1,1:nz  ) + & 
!            (A(2,2:nx,1:ny-1,1:nz) - A(2,1:nx-1,1:ny-1,1:nz))/dx - &
!            (A(1,1:nx-1,2:ny,1:nz) - A(1,1:nx-1,1:ny-1,1:nz))/dy
       
#ifndef ISO
!         u(iena,:,:,:) = eint + ekin + 0.5*(b(1,:,:,:)**2 + &
!            b(2,:,:,:)**2 + b(3,:,:,:)**2)
!         if (allocated(eint)) deallocate(eint)
!         if (allocated(ekin)) deallocate(ekin)
#endif /* ISO */
!         if (allocated(A)) deallocate(A)
!     endif

   end subroutine magn_multipole_sn



!--------------------------------------------------------------------------
 
   subroutine rand_coords(pos)
! Written by M. Hanasz
      use start,  only : h_sn,dimensions,xmin,ymin
      use grid,   only : Lx,Ly
#ifdef SHEAR
      use arrays, only : y,js,je
      use grid,   only : dy
      use start,  only : nyd

      integer :: jsn,jremap 
      real :: dysn
#endif /* SHEAR */
    
      real, dimension(3), intent(out) :: pos
      real, dimension(4) :: rand
      real :: xsn,ysn,zsn,znorm

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
#else /* SHEAR */
      ysno = ysn
      ysni = ysn
#endif /* SHEAR */

      pos(1) = xsn
      pos(2) = ysn
      pos(3) = zsn
      return
  
   end subroutine rand_coords

!-----------------------------------------------------------------------

!>
!! \brief Function generate point on the surface of unit
!! sphere with uniform distribution and returns its latidude and longitude
!<
   function rand_angles()
  
      use constants, only : pi
      implicit none
      real :: rand(2)                  
      real :: rnx,rny,rnz               !> Point's position in cartesian coordinates
      real, dimension(2) :: rand_angles !> Latidue and longitude

      call random_number(rand)
      rnz = (1.0-2.0*rand(1))
      rnx = sqrt(1.0-rnz**2)*cos(2.*pi*rand(2))
      rny = sqrt(1.0-rnz**2)*sin(2.*pi*rand(2))
      rand_angles(1) = 2.0*pi*rand(2)
      rand_angles(2) = acos(rnz/sqrt(rnz**2+rny**2+rnx**2) )
 
   end function rand_angles

!=======================================================================
!
!      \\\\\\\         B E G I N   S U B R O U T I N E S        ///////
!      ///////    F R O M   N U M E R I C A L   R E C I P E S   \\\\\\\
!
!=======================================================================
! idum and rsq variables were moved out due to being unused

   function gasdev(x,y)

      implicit none
      real :: x, y, x1, y1,  r
      real :: gasdev
      real :: fac
      real, dimension(2) :: rand
      real, save :: gset
      integer, save :: iset, irand
 
      if (iset.eq.0) then
1        x1 = 2.0*x - 1.0
         y1 = 2.0*y - 1.0
         r  = x1**2 + y1**2
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
