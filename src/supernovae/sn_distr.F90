! $Id$
#include "piernik.def"

module sn_distr

! Supernovae distribution and drawing
! Written by: D. Woltanski, December 2007

   real, dimension(2)        :: SNtrest, SNheight, SNfreq
   real, dimension(2,1000,2) :: danta
   real    :: MexplSN,EexplSN,emagadd,tot_emagadd
   integer :: r0snx,r0sny,r0snz
   integer, dimension(2)   :: SNnohistory
   integer, dimension(3)   :: rintsn
   integer :: itype

   contains

!-----------------------------------------------------------------------------

   subroutine prepare_SNdistr
      use start, only: snemass,snenerg,sn1time,sn2time,r0sn
      use grid, only: dx,dy,dz
      use constants
      implicit none
      real :: RmaxI, RmaxII, rcl, rc
      integer :: i,imax

      emagadd     = 0.0
      tot_emagadd = 0.0
      imax        = 1e3
      RmaxI       = 50.0*kpc
      RmaxII      = 15.0*kpc
      SNheight(1) = 100.0*pc !325.0*pc 	!exp(-|z|/325*pc)
      SNheight(2) = 100.0*pc !266.0*pc	!exp(-|z|/266*pc)
      SNtrest(:)  = 0.0
      SNfreq(1)   = 1./sn1time/year
      SNfreq(2)   = 1./sn2time/year

      r0snx = int(r0sn/dx)+1
      r0sny = int(r0sn/dy)+1
      r0snz = int(r0sn/dz)+1

      MexplSN  = snemass*Msun
#ifndef ISO
      EexplSN  = snenerg*erg
#endif /* ISO */
      call random_seed()

      rcl=0.0
      rc=1./real(imax)*RmaxI
      danta(1,1,1)=(2.*pi*2.6/1.e6*exp(-((rc -8.5*kpc)/(4.9*kpc)))*rc &
         +(2.*pi*2.6/1.e6*exp(-((rcl-8.5*kpc)/(4.9*kpc)))*rcl))*(rc-rcl)*0.5
      danta(1,1,2)=rc
      rcl=rc
      do i=2,imax
         rc=real(i)/real(imax)*RmaxI
         danta(1,i,1)=danta(1,i-1,1)+(2.*pi*2.6/1.e6*exp(-((rc -8.5*kpc)/(4.9*kpc)))*rc &
            +(2.*pi*2.6/1.e6*exp(-((rcl-8.5*kpc)/(4.9*kpc)))*rcl))*(rc-rcl)*0.5
         danta(1,i,2)=rc
         rcl=rc
      enddo
      danta(1,:,1)=danta(1,:,1)/danta(1,imax,1)

      rcl=0.0
      rc=1./real(imax)*RmaxII
      rc=real(i)/real(imax)*RmaxII
      danta(2,1,1)=((2.*pi*19./1.e6*exp(-((rc -4.5*kpc)**2-(4.0*kpc)**2)/(2.9*kpc)**2)*rc )&
         +(2.*pi*19./1.e6*exp(-((rcl-4.5*kpc)**2-(4.0*kpc)**2)/(2.9*kpc)**2)*rcl))*(rc-rcl)*0.5
      danta(2,1,2)=rc
      rcl=rc
      do i=2,imax
         rc=real(i)/real(imax)*RmaxII
         danta(2,i,1)=danta(2,i-1,1)+((2.*pi*19./1.e6*exp(-((rc -4500.)**2-(4000.0)**2)/(2900.)**2)*rc )&
            +(2.*pi*19./1.e6*exp(-((rcl-4500.)**2-(4000.0)**2)/(2900.)**2)*rcl))*(rc-rcl)*0.5
         danta(2,i,2)=rc
         rcl=rc
      enddo
      danta(2,:,1)=danta(2,:,1)/danta(2,imax,1)
      SNnohistory(:)=0

   end subroutine prepare_SNdistr

!===============================================================================================
   subroutine supernovae_distribution
      use mpi_setup
      use start, only : dt,smallei,nb
      use arrays, only : u,b,nx,ny,nz,ibx,iby,ibz,idna,imxa,imya,imza,iena
      use grid, only : dx,dy,dz,dvol
      use sn_sources, only : rand_angles
      use mag_boundaries, only : compute_b_bnd, bnd_a
      implicit none

      integer :: isn
      real, dimension(3) :: snpos
      real, dimension(2) :: dtime
      integer, dimension(2) :: SNno, pot
      real, allocatable, dimension(:,:) :: snposarray,snangarray
      real, dimension(:,:,:,:), allocatable :: A
#ifndef ISO
      real, dimension(:,:,:), allocatable :: ekin,eint
#endif /* ISO */

      SNno(:)=0
      if(proc .eq. 0) then
         dtime=SNtrest+2*dt
         pot  = relato0((/0.0,0.0/),dtime*SNfreq)     ! zabezpiecza przed ujemna iloscia wybuchow
         SNno=pot*int(dtime*SNfreq)
         SNtrest=dtime-(real(SNno))/SNfreq
#ifdef VERBOSE
         call write_sninfo(SNno)
#endif /* VERBOSE */
      endif
      call MPI_BCAST(SNno, 2, MPI_INTEGER, 0, comm, ierr)

      allocate(snposarray(sum(SNno,1),3), snangarray(sum(SNno,1),2))
      if(proc .eq. 0) then
         do itype = 1,2
            if(SNno(itype) .gt. 0) then
               do isn=1,SNno(itype)
                  call rand_galcoord(snpos)
                  snposarray(isn+SNno(1)*(itype-1),:)=snpos
                  snangarray(isn+SNno(1)*(itype-1),:)=rand_angles()
               enddo
            endif
         enddo
      endif
      call MPI_BCAST(snposarray, 3*sum(SNno,1), MPI_DOUBLE_PRECISION, 0, comm, ierr)
      call MPI_BCAST(snangarray, 2*sum(SNno,1), MPI_DOUBLE_PRECISION, 0, comm, ierr)
#ifdef VERBOSE
      itype = 1
#endif /* VERBOSE */

      allocate(A(3,nx,ny,nz))
      A(:,:,:,:) = 0.0
#ifndef ISO
      allocate(eint(nx,ny,nz))
      allocate(ekin(nx,ny,nz))
      ekin(:,:,:) = 0.5*( u(imxa,:,:,:)**2 + u(imya,:,:,:)**2 + u(imza,:,:,:) )**2 / u(idna,:,:,:)
      eint(:,:,:) = max( u(iena,:,:,:) - ekin -  0.5*(b(ibx,:,:,:)**2 + b(iby,:,:,:)**2 + &
                         b(ibz,:,:,:)**2) , smallei)
      emagadd = sum(u(iena,nb+1:nx-nb,nb+1:ny-nb,nb+1:nz-nb))*dvol
#endif /* ISO */

      do isn=1,sum(SNno,1)
#ifdef VERBOSE
         if(isn .gt. SNno(1)) itype = 2
#endif /* VERBOSE */
         call add_explosion(snposarray(isn,:),snangarray(isn,:),A)
#ifdef VERBOSE
         if(proc .eq. 0) write(*,*) 'added ',isn,'. SN of ',sum(SNno,1)
#endif /* VERBOSE */
      enddo
      call MPI_BARRIER(comm,ierr)

      call bnd_a(A)

      b(ibx,1:nx,  1:ny-1,1:nz-1) = b(ibx,1:nx,  1:ny-1,1:nz-1) + &
         (A(3,1:nx,2:ny,1:nz-1) - A(3,1:nx,1:ny-1,1:nz-1))/dy - &
         (A(2,1:nx,1:ny-1,2:nz) - A(2,1:nx,1:ny-1,1:nz-1))/dz

      b(iby,1:nx-1,1:ny,  1:nz-1) = b(iby,1:nx-1,1:ny,  1:nz-1) + &
         (A(1,1:nx-1,1:ny,2:nz) - A(1,1:nx-1,1:ny,1:nz-1))/dz - &
         (A(3,2:nx,1:ny,1:nz-1) - A(3,1:nx-1,1:ny,1:nz-1))/dx

      b(ibz,1:nx-1,1:ny-1,1:nz  ) = b(ibz,1:nx-1,1:ny-1,1:nz  ) + &
         (A(2,2:nx,1:ny-1,1:nz) - A(2,1:nx-1,1:ny-1,1:nz))/dx - &
         (A(1,1:nx-1,2:ny,1:nz) - A(1,1:nx-1,1:ny-1,1:nz))/dy

      deallocate(A)
#ifndef ISO
      u(iena,:,:,:) = eint + ekin + 0.5*(b(1,:,:,:)**2 + &
           b(2,:,:,:)**2 + b(3,:,:,:)**2)
      emagadd = sum(u(iena,nb+1:nx-nb,nb+1:ny-nb,nb+1:nz-nb))*dvol - emagadd
      deallocate(ekin,eint)
#endif
      call compute_b_bnd
      deallocate(snposarray,snangarray)
      return

   end subroutine supernovae_distribution

!------------------------------------------------

   subroutine add_explosion(snpos,angles,A)
      use arrays, only: nx,ny,nz,x,y,z,nxb,nyb,nzb,xdim,ydim,zdim,dl,u
      use grid, only: xminb,xmaxb,yminb,ymaxb,zminb,zmaxb,dx,dy,dz,Lx,Ly,Lz
      use start, only: r0sn,nb,add_mass,add_ener,add_encr,add_magn
      use constants
#ifdef COSM_RAYS
      use arrays, only : iecr
      use start, only  : r_sn, cr_eff
#endif /* COSM_RAYS */
#ifdef DIPOLS
      use sn_sources, only : magn_multipole_sn
#endif /* DIPOLS */
      implicit none
      real, dimension(3) :: snpos
      real, dimension(2) :: angles
      real    :: r1sn
      real    :: massadd,eneradd,encradd,normscal,normvalu
      real, dimension(:,:,:,:) :: A
      integer :: ic,jc,kc,i,j,k
      real    :: e_sn, amp_sn, amp_cr


      if((snpos(1)+r0sn .ge. xminb-nb*dl(xdim)) .and. &
         (snpos(1)-r0sn .le. xmaxb+nb*dl(xdim))) then
         if((snpos(2)+r0sn .ge. yminb-nb*dl(ydim)) .and. &
            (snpos(2)-r0sn .le. ymaxb+nb*dl(ydim))) then
            if((snpos(3)+r0sn .ge. zminb-nb*dl(zdim)) .and. &
               (snpos(3)-r0sn .le. zmaxb+nb*dl(zdim))) then

               ic =nb+int((snpos(1)-xminb)/(xmaxb-xminb)*nxb)
               jc =nb+int((snpos(2)-yminb)/(ymaxb-yminb)*nyb)
               kc =nb+int((snpos(3)-zminb)/(zmaxb-zminb)*nzb)
               massadd=0.0
               eneradd=0.0
               encradd=0.0
               normscal=1./0.427796
               normvalu=normscal/(sqrt(pi)*r0sn)**3

               do i = ic-r0snx,ic+r0snx
                  do j = jc-r0sny,jc+r0sny
                     do k = kc-r0snz,kc+r0snz
                        if((i .ge. 1) .and. (i .le. nx) .and.&
                           (j .ge. 1) .and. (j .le. ny) .and.&
                           (k .ge. 1) .and. (k .le. nz)) then

                           r1sn = sqrt((snpos(1)-x(i))**2+(snpos(2)-y(j))**2+(snpos(3)-z(k))**2)
                           if(r1sn .lt. r0sn) then
                              if(add_mass .eq. 'yes') then
                                 u(2,i,j,k)=u(2,i,j,k)/u(1,i,j,k)*(u(1,i,j,k)+MexplSN*exp(-r1sn**2/r0sn**2))*normvalu
                                 u(3,i,j,k)=u(3,i,j,k)/u(1,i,j,k)*(u(1,i,j,k)+MexplSN*exp(-r1sn**2/r0sn**2))*normvalu
                                 u(1,i,j,k)=u(1,i,j,k)+MexplSN*exp(-r1sn**2/r0sn**2)*normvalu
                                 massadd=massadd+MexplSN*exp(-r1sn**2/r0sn**2)*normvalu*dx*dy*dz
                              endif
#ifndef ISO
                              if(add_ener .eq. 'yes') then
                                 u(5,i,j,k)=u(5,i,j,k)+EexplSN*exp(-r1sn**2/r0sn**2)*normvalu
                                 eneradd=eneradd+EexplSN*exp(-r1sn**2/r0sn**2)*normvalu*dx*dy*dz
                              endif
#endif /* ISO */
#ifdef COSM_RAYS
                              if(add_encr .eq. 'yes') then
                                 e_sn = 1.0e+51*erg
                                 amp_sn = e_sn * normvalu
                                 amp_cr = cr_eff * amp_sn

                                 u(iecr,i,j,k) = u(iecr,i,j,k) + amp_cr*exp(-r1sn**2/r0sn**2)

                                 encradd=encradd+amp_cr*exp(-r1sn**2/r0sn**2)*dx*dy*dz
                              endif
#endif /* COSM_RAYS */
                           endif
                        endif
                     enddo
                  enddo
               enddo
#ifdef VERBOSE
               write(*,'(a3,i1.1,a11,i4.4,i4.4,i4.4)') 'SN ',itype,' position: ',ic,jc,kc
               write(*,'(a15,f8.3,1x,f8.3,1x,f8.3)')   '     coords:   ',snpos(1),snpos(2),snpos(3)

               if((ic .ge. 1) .and. (ic .le. nx) .and. (jc .ge. 1) .and. &
                  (jc .le. ny) .and. (kc .ge. 1) .and. (kc .le. nz)) &
                     write(*,'(a15,f8.3,1x,f8.3,1x,f8.3)')   '     that is:  ',x(ic),y(jc),z(kc)
               if(add_mass .eq. 'yes') write(*,'(a19,e15.8,a5)') '   mass injection: ',massadd/Msun,' Msun'
               if(add_ener .eq. 'yes') write(*,'(a19,e15.8,a4)') ' energy injection: ',eneradd/erg,' erg'
               if(add_encr .eq. 'yes') write(*,'(a19,e15.8,a4)') 'CR energy inject.: ',encradd/erg,' erg'
#endif /* VERBOSE */
#ifdef DIPOLS
#ifndef DIP9
               if(add_magn == 'yes') call magn_multipole_sn(angles,snpos,A)
#endif /* DIP9 */
#endif /* DIPOLS */
            endif
         endif
      endif
#ifdef DIP9
        
      if(add_magn == 'yes') call magn_multipole_sn(angles,snpos,A)

#endif /* DIP9 */
      return
   end subroutine add_explosion

!--------------------------------------------------------------------------

   subroutine rand_galcoord(snpos)
      use constants, only: pi
      implicit none
      real, dimension(3) :: snpos
      real, dimension(4) :: los4
      real :: radius, azym
      integer ii
!      real rand
!      external rand

      call random_number(los4)

      ii=1
      do while(danta(itype,ii,1) .lt. los4(1))
         ii=ii+1
      enddo
      radius = danta(itype,ii+1,2)  -   (danta(itype,ii+1,2)-danta(itype,ii,2)) &
         *(danta(itype,ii+1,1)-los4(1))/(danta(itype,ii+1,1)-danta(itype,ii,1))
      azym = 2.*pi*los4(2)
      snpos(1) = radius*cos(azym)
      snpos(2) = radius*sin(azym)
      snpos(3) = gasdev(los4(3),los4(4))*SNheight(itype)

      return
   end subroutine rand_galcoord

!--------------------------------------------------

   subroutine write_sninfo(SNno)
      use start, only : t,dt
      implicit none
      integer, dimension(2) :: SNno

      write(*,'(a12,i8,a7,i8,a6)') 'explosions: ',SNno(1),' SN I, ',SNno(2),' SN II'
      SNnohistory = SNnohistory + SNno
      write(*,'(a22,f8.4,a9,f10.4)') ' SNE frequency: SN I: ',SNno(1)/2./dt,', SN II: ',SNno(2)/2./dt
      write(*,'(a22,f8.4,a9,f10.4)') 'mean frequency: SN I: ',SNnohistory(1)/t,', SN II: ',SNnohistory(2)/t
      return
   end subroutine write_sninfo

!--------------------------------------------------

   function relato0(a,b)
      implicit none
      real,dimension(2) :: a,b,relato0

      where(b > a)
         relato0 = 1
      elsewhere
         relato0 = 0
      endwhere
      return
   end function relato0

!=======================================================================
!
!      \\\\\\\         B E G I N   S U B R O U T I N E S        ///////
!      ///////    F R O M   N U M E R I C A L   R E C I P E S   \\\\\\\
!
!=======================================================================


   function gasdev(x,y)

      implicit none
      real x, y, x1, y1,  r
      real gasdev, rand(2)
      real fac
      real, save :: gset
      integer, save :: iset, irand

      if (iset.eq.0) then
1        x1=2.*x-1.
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

end module sn_distr

