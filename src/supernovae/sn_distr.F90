#include "mhd.def"

module sn_distr
  
! Supernovae distribution and drawing
! Written by: D. Woltanski, December 2007

  real,dimension(2)       :: SNtrest, SNheight, SNfreq
  real, dimension(2,1000,2) :: danta
  real EexplSN,r0snx,r0sny,r0snz
  integer, dimension(2)   :: SNnohistory
  integer, dimension(3)   :: rintsn
  integer itype


contains

!!-----------------------------------------------------------------------------


  subroutine prepare_SNdistr
  use start, only: snenerg,sn1time,sn2time,r0sn
  use grid, only: dx,dy,dz
  use constants, only: pi,kpc,pc,year
#ifndef ISO
  use constants, only: erg
#endif /* ISO */
  implicit none
  real RmaxI, RmaxII, rcl, rc
  integer i,imax

  imax       = 1e3
  RmaxI      = 50.0*kpc
  RmaxII     = 15.0*kpc
  SNheight(1) = 100.0*pc !325.0*pc 	!exp(-|z|/325*pc)
  SNheight(2) = 100.0*pc !266.0*pc	!exp(-|z|/266*pc)
  SNtrest(:) = 0.0
  SNfreq(1) = 1./sn1time/year
  SNfreq(2) = 1./sn2time/year
  r0snx = int(r0sn/dx)+1
  r0sny = int(r0sn/dy)+1
  r0snz = int(r0sn/dz)+1
#ifdef ISO
  EexplSN  = snenerg
#else /* ISO */
  EexplSN  = snenerg*erg/(4./3.*pi*r0sn**3)
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
    use start, only : dt
    implicit none

    integer i,j,k,ic,jc,kc, isn
    real, dimension(3) :: snpos
    real, dimension(2) :: dtime
    integer, dimension(2) :: SNno, pot
    real, allocatable, dimension(:,:) :: snposarray
    

    SNno(:)=0
    if(proc .eq. 0) then
      dtime=SNtrest+2*dt
      pot  = relato0((/0.0,0.0/),dtime*SNfreq)	   ! zabezpiecza przed ujemna iloscia wybuchow
      SNno=pot*int(dtime*SNfreq)
      SNtrest=dtime-(real(SNno))/SNfreq
      
!      call write_sninfo(SNno)
    endif
    call MPI_BCAST(SNno, 2, MPI_INTEGER, 0, comm, ierr)
    
    allocate(snposarray(sum(SNno,1),3))
    if(proc .eq. 0) then
      do itype = 1,2
        if(SNno(itype) .gt. 0) then
          do isn=1,SNno(itype)
	    call rand_galcoord(snpos)
	    snposarray(isn+SNno(1)*(itype-1),:)=snpos
          enddo
        endif
      enddo
    endif
    call MPI_BCAST(snposarray, 3*sum(SNno,1), MPI_DOUBLE_PRECISION, 0, comm, ierr)
    do isn=1,sum(SNno,1)
      call add_explosion(snpos)
    enddo

    return
  
  end subroutine supernovae_distribution  

!------------------------------------------------

  subroutine add_explosion(snpos)
  use arrays, only: nx,ny,nz,x,y,z,nxb,nyb,nzb,xdim,ydim,zdim,dl,u
  use grid, only: xminb,xmaxb,yminb,ymaxb,zminb,zmaxb
  use start, only: r0sn,nb
  use constants
#ifdef COSM_RAYS
  use arrays, only : iecr
  use start, only  : r_sn, cr_eff           
#endif COSM_RAYS
  implicit none
  real, dimension(3) :: snpos
  real r1sn
  integer ic,jc,kc,i,j,k
  real e_sn, amp_sn, amp_cr 
  
  
  if((snpos(1)+r0sn .ge. xminb-nb*dl(xdim)) .and. (snpos(1)-r0sn .le. xmaxb+nb*dl(xdim))) then
  if((snpos(2)+r0sn .ge. yminb-nb*dl(ydim)) .and. (snpos(2)-r0sn .le. ymaxb+nb*dl(ydim))) then
  if((snpos(3)+r0sn .ge. zminb-nb*dl(zdim)) .and. (snpos(3)-r0sn .le. zmaxb+nb*dl(zdim))) then
    ic =nb+int((snpos(1)-xminb)/(xmaxb-xminb)*nxb)
    jc =nb+int((snpos(2)-yminb)/(ymaxb-yminb)*nyb)
    kc =nb+int((snpos(3)-zminb)/(zmaxb-zminb)*nzb)
    do i = ic-r0snx,ic+r0snx
    do j = jc-r0sny,jc+r0sny
    do k = kc-r0snz,kc+r0snz
      if((i .ge. 1) .and. (i .le. nx) .and. (j .ge. 1) .and. (j .le. ny) .and. (k .ge. 1) .and. (k .le. nz)) then
!	      write(*,'(a3,i1.1,a11,i4.4,i4.4,i4.4)') 'SN ',itype,' position: ',i,j,k
!	      write(*,'(a15,f8.3,1x,f8.3,1x,f8.3)') '     coords:   ',snpos(1),snpos(2),snpos(3)
!	      write(*,'(a15,f8.3,1x,f8.3,1x,f8.3)') '     that is:  ',x(i),y(j),z(k)
        r1sn = sqrt((snpos(1)-x(i))**2+(snpos(2)-y(j))**2+(snpos(3)-z(k))**2)
	if(r1sn .lt. r0sn) then
#ifdef ISO
!	  u(2,i,j,k)=u(2,i,j,k)/u(1,i,j,k)*(u(1,i,j,k)+EexplSN*exp(-r1sn**2/r0sn**2))
!	  u(3,i,j,k)=u(3,i,j,k)/u(1,i,j,k)*(u(1,i,j,k)+EexplSN*exp(-r1sn**2/r0sn**2))
!          u(1,i,j,k)=u(1,i,j,k)+EexplSN*exp(-r1sn**2/r0sn**2)
#else /* ISO */
!          u(5,i,j,k)=u(5,i,j,k)+EexplSN*exp(-r1sn**2/r0sn**2)
#endif /* ISO */

#ifdef COSM_RAYS
         
	 e_sn = 10e+51*erg
	 
         amp_sn = e_sn/(sqrt(pi)*r0sn*pc)**3 
	 
	 amp_cr = cr_eff * amp_sn

         u(iecr,i,j,k) = u(iecr,i,j,k) + amp_cr*exp(-r1sn**2/r0sn**2)
	 
!	 write(*,*) e_sn, amp_cr, r0sn, pc, erg
!	 stop
	 
#endif COSM_RAYS

	endif
      endif
    enddo
    enddo
    enddo
  endif
  endif
  endif
  
  return
  end subroutine add_explosion

!--------------------------------------------------------------------------

 subroutine rand_galcoord(snpos)
 use constants, only: pi
 implicit none
 real, dimension(3) :: snpos
 real, dimension(4) :: los4
 real radius, azym
 integer ii, i
 real rand
 external rand

 call random_number(los4)

! do i=1,4
!   los4(i) =rand()
! enddo

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
 use start, only: t,dt
 implicit none
 integer, dimension(2) :: SNno
      write(*,'(a12,i8,a7,i8,a6)') 'explosions: ',SNno(1),' SN I, ',SNno(2),' SN II'
      SNnohistory = SNnohistory + SNno
      write(*,'(a22,f8.4,a9,f8.4)') ' SNE frequency: SN I: ',SNno(1)/2./dt,', SN II: ',SNno(2)/2./dt
      write(*,'(a22,f8.4,a9,f8.4)') 'mean frequency: SN I: ',SNnohistory(1)/t,', SN II: ',SNnohistory(2)/t
 return
 end subroutine write_sninfo
  
!--------------------------------------------------

  function relato0(a,b)
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
      integer idum
      real x, y, x1, y1,  r
      real gasdev, rand(2)
      real fac,rsq
      real, save :: gset
      integer, save :: iset, irand
 
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

end module sn_distr

