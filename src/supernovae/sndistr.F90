! $Id$
!
! PIERNIK Code Copyright (C) 2006 Michal Hanasz
!
!    This file is part of PIERNIK code.
!
!    PIERNIK is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    PIERNIK is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with PIERNIK.  If not, see <http://www.gnu.org/licenses/>.
!
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen 
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD 
!             for original source code "mhd.f90" 
!   
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
#include "piernik.def"

module sndistr

! Supernovae distribution and drawing

   real, dimension(2)        :: SNtrest, SNheight, SNfreq
   real, dimension(2,1000,2) :: danta
   real    :: MexplSN,EexplSN,emagadd,tot_emagadd
   integer :: r0snx,r0sny,r0snz
   integer, dimension(2)   :: SNnohistory
   integer, dimension(3)   :: rintsn
   integer :: itype

   contains

!-----------------------------------------------------------------------------

   subroutine init_sndistr
      use start, only: snemass,snenerg,sn1time,sn2time,r0sn
      use grid, only: dx,dy,dz
      use constants
      implicit none
      integer :: seed
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

      call system_clock(seed)
      call rng_initialise(seed)

   end subroutine init_sndistr

!===============================================================================================
   subroutine supernovae_distribution
      use mpisetup
      use start, only : dt,smallei
      use arrays, only : u,b,nx,ny,nz,ibx,iby,ibz,idna,imxa,imya,imza,iena
      use grid, only : dx,dy,dz,dvol,nb
      use snsources, only : rand_angles
      use magboundaries, only : all_mag_boundaries, bnd_a
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
      call all_mag_boundaries
      deallocate(snposarray,snangarray)
      return

   end subroutine supernovae_distribution

!------------------------------------------------

   subroutine add_explosion(snpos,angles,A)
      use arrays, only: nx,ny,nz,x,y,z,nxb,nyb,nzb,xdim,ydim,zdim,dl,u
      use grid, only: xminb,xmaxb,yminb,ymaxb,zminb,zmaxb,dx,dy,dz,Lx,Ly,Lz,nb
      use start, only: r0sn,add_mass,add_ener,add_encr
      use constants
#ifdef COSM_RAYS
      use arrays, only : iecr
      use start, only  : r_sn, cr_eff
#endif /* COSM_RAYS */
#ifdef DIPOLS
      use snsources, only : magn_multipole_sn
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
                  call magn_multipole_sn(angles,snpos,A)
#endif /* DIP9 */
#endif /* DIPOLS */
            endif
         endif
      endif
#ifdef DIP9
      if((snpos(1) .ge. xminb-0.25*Lx) .and. &
         (snpos(1) .le. xmaxb+0.25*Lx)) then
         if((snpos(2) .ge. yminb-0.25*Ly) .and. &
            (snpos(2) .le. ymaxb+0.25*Ly)) then
            if((snpos(3) .ge. zminb-0.25*Lz) .and. &
               (snpos(3) .le. zmaxb+0.25*Lz))then

                  call magn_multipole_sn(angles,snpos,A)

            endif
         endif
      endif
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
      integer, parameter :: n=1
      real, dimension(n) :: x

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
      call rng_sample_gaussian(x,n,1d0)
      snpos(3) = x(1) * SNheight(itype)

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

end module sndistr

