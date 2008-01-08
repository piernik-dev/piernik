#include "mhd.def"

module sn_distr
  
! Supernovae distribution and drawing
! Written by: D. Woltanski, December 2007

  use arrays
  use grid
  use constants

  real SNIrest, SNIIrest
  real, dimension(1000,2) :: danta1, danta2
  integer, dimension(2)   :: SNnohistory


contains

!-----------------------------------------------------------------------------


  subroutine compute_SNdistr
  implicit none
  real RmaxI, RmaxII, rcl, rc
  integer i,imax

  imax     = 1e3
  RmaxI    = 50000.0
  RmaxII   = 15000.0
  SNIrest  = 0.0
  SNIIrest = 0.0
  call random_seed()

  rcl=0.0
  rc=1./real(imax)*RmaxI
  danta1(1,1)=(2.*pi*2.6/1.e6*exp(-((rc-8500.)/(4900.)))*rc &
	   +(2.*pi*2.6/1.e6*exp(-((rcl-8500.)/(4900.)))*rcl))*(rc-rcl)*0.5
  danta1(1,2)=rc
  rcl=rc
  do i=2,imax
    rc=real(i)/real(imax)*RmaxI
    danta1(i,1)=danta1(i-1,1)+(2.*pi*2.6/1.e6*exp(-((rc-8500.)/(4900.)))*rc &
	   +(2.*pi*2.6/1.e6*exp(-((rcl-8500.)/(4900.)))*rcl))*(rc-rcl)*0.5
    danta1(i,2)=rc
    rcl=rc
  enddo
  danta1(:,1)=danta1(:,1)/danta1(imax,1)
    
  rcl=0.0
  rc=1./real(imax)*RmaxII
  rc=real(i)/real(imax)*RmaxII
  danta2(1,1)=((2.*pi*19./1.e6*exp(-((rc-4500.)**2-(4000.0)**2)/(2900.)**2)*rc)&
	   +(2.*pi*19./1.e6*exp(-((rcl-4500.)**2-(4000.0)**2)/(2900.)**2)*rcl))*(rc-rcl)*0.5
  danta2(1,2)=rc
  rcl=rc
  do i=2,imax
    rc=real(i)/real(imax)*RmaxII
    danta2(i,1)=danta2(i-1,1)+((2.*pi*19./1.e6*exp(-((rc-4500.)**2-(4000.0)**2)/(2900.)**2)*rc)&
	   +(2.*pi*19./1.e6*exp(-((rcl-4500.)**2-(4000.0)**2)/(2900.)**2)*rcl))*(rc-rcl)*0.5
    danta2(i,2)=rc
    rcl=rc
  enddo
  danta2(:,1)=danta2(:,1)/danta2(imax,1)
  SNnohistory(:)=0

  end subroutine compute_SNdistr


!===============================================================================================
  subroutine supernovae_distribution
    use mpi_setup
    use start, only : nb, t, dt, snenerg, sn1time, sn2time
    implicit none

    integer i,j,k, SNInum, SNIInum, isn, ii
    real iburst, burstdistr, rdist, zz, dtime, randomstart, SNIsum, SNIIsum, SNIsumall, SNIIsumall, los
    real radius, azym, SNIfreq, SNIIfreq, EexplSN
    real, dimension(3) :: snpos
    integer, dimension(2) :: SNno
    
    SNIfreq  = 1./sn1time/year
    SNIIfreq = 1./sn2time/year
#ifdef ISO
    EexplSN  = snenerg
#else /* ISO */
    EexplSN  = snenerg*erg
#endif /* ISO */

!===============obliczanie, jaka czesc SN z poprzedniego kroku trzeba doliczyc=======
    SNno(1)=0 					! zerowanie ilosci wylosowanych SNI w aktualnym kroku
    SNno(2)=0					! zerowanie ilosci wylosowanych SNII w aktualnym kroku
    if(proc .eq. 0) then			! losowanie przeprowadzone tylko w zerowym procesie
      dtime=SNIrest+2*dt			! czas dla losowania ilosci SN powiekszony o niewylosowane SN z poprzedniego kroku
      call random_number(los)
      do while(los .lt. dtime*SNIfreq)	! SN wylosowana jesli los jest mniejszy niz ilosc przewidywana srednia
        SNno(1)=SNno(1)+1			! wtedy wylosowana SN dodawana jest do ilosci wylosowanych w danym kroku
        SNIrest=dtime-1./SNIfreq		! obliczanie pozostalego czasu
!        SNIrest=dtime-2.*los/SNIfreq		! obliczanie pozostalego czasu - druga wersja
        dtime=SNIrest				! ustawianie nowego czasu losowania
        call random_number(los)			! powtorne losowanie i sprawdzenie warunku
      enddo
      dtime=SNIIrest+2*dt			!analogicznie jak dla SNI
      call random_number(los)
      do while(los .lt. dtime*SNIIfreq)
        SNno(2)=SNno(2)+1
	SNIIrest=dtime-1./SNIIfreq
!	SNIIrest=dtime-2.*los/SNIIfreq
	dtime=SNIIrest
	call random_number(los)
      enddo
      write(*,'(a12,i4.4,a7,i4.4,a6)') 'explosions: ',SNno(1),' SN I, ',SNno(2),' SN II'
      SNnohistory = SNnohistory + SNno
      write(*,'(a22,f8.4,a9,f8.4)') ' SNE frequency: SN I: ',SNno(1)/2./dt,', SN II: ',SNno(2)/2./dt
      write(*,'(a22,f8.4,a9,f8.4)') 'mean frequency: SN I: ',SNnohistory(1)/t,', SN II: ',SNnohistory(2)/t
    endif
    call MPI_BCAST(SNno, 2, MPI_INTEGER, 0, comm, ierr)	! rozeslanie informacji o ilosci wylosowanych SN
    
    if(SNno(1) .gt. 0) then				! jesli zostaly wylosowane jakies SNI, to losowane sa ich pozycje
      do isn=1,SNno(1)
        if(proc .eq. 0) then
	  call random_number(los)
!	  sigmaI=2.6/kpc**2/Myr*exp(-(radius-Rgc_sun)/4.9/kpc)
!	  radius = 2.4*5314.4*tan(los*1.5027)*exp(-los**2*1.4)			! ustalany jest promien cylindryczny
	  ii=1
	  do while(danta1(ii,1) .lt. los)
	  ii=ii+1
	  enddo
	  radius =danta1(ii+1,2)-(danta1(ii+1,2)-danta1(ii,2))*(danta1(ii+1,1)-los)/(danta1(ii+1,1)-danta1(ii,1))
	  call random_number(los)
	  call random_number(azym)				!exp(-|z|/325*pc)
	  snpos(3) = sqrt(-2.*log(azym))*cos(2.*pi*los)*100.0*pc !325.0*pc	! ustalana jest wysokosc nad plaszczyzna Galaktyki
	  call random_number(los)
	  azym = 2.*pi*los					! ustalany jest azymut
	  snpos(1) = radius*cos(azym)				! dokladne liczenie wspolrzednych x,y,z na podstawie promienia, z i azymutu
	  snpos(2) = radius*sin(azym)
!	  snpos(3) = z
	endif
	call MPI_BCAST(snpos, 3, MPI_DOUBLE_PRECISION, 0, comm, ierr)		! przekazywanie pozycji wybuchu SN do odpowiedniego procesu
	if((snpos(1) .ge. xminb-nb*dl(xdim)) .and. (snpos(1) .le. xmaxb+nb*dl(xdim))) then
	if((snpos(2) .ge. yminb-nb*dl(ydim)) .and. (snpos(2) .le. ymaxb+nb*dl(ydim))) then
	if((snpos(3) .ge. zminb-nb*dl(zdim)) .and. (snpos(3) .le. zmaxb+nb*dl(zdim))) then
	  i =nb+int((snpos(1)-xminb)/(xmaxb-xminb)*nxb)
	  j =nb+int((snpos(2)-yminb)/(ymaxb-yminb)*nyb)
	  k =nb+int((snpos(3)-zminb)/(zmaxb-zminb)*nzb)
	  write(*,'(a15,i4.4,i4.4,i4.4)') 'SN I position: ',i,j,k
	  write(*,'(a15,f8.3,1x,f8.3,1x,f8.3)') '     coords:   ',snpos(1),snpos(2),snpos(3)
	  write(*,'(a15,f8.3,1x,f8.3,1x,f8.3)') '     that is:  ',x(i),y(j),z(k)
#ifdef ISO
	  u(2,i,j,k)=u(2,i,j,k)/u(1,i,j,k)*(u(1,i,j,k)+EexplSN)
	  u(3,i,j,k)=u(3,i,j,k)/u(1,i,j,k)*(u(1,i,j,k)+EexplSN)
          u(1,i,j,k)=u(1,i,j,k)+EexplSN
#else /* ISO */
          u(5,i,j,k)=u(5,i,j,k)+EexplSN
#endif /* ISO */
	endif
	endif
	endif
      enddo
    endif
    
    if(SNno(2) .gt. 0) then			! analogicznie jak dla SNI
      do isn=1,SNno(2)
        if(proc .eq. 0) then
	  call random_number(los)
!	  sigmaII=19.0/kpc**2/Myr*exp(-((radius-4.5*kpc)**2-(r_gc_sun-4.5*kpc)**2)/(2.9*kpc)**2)
!	  radius = 1.2*5178.1*tan(sqrt(los**2+2*los)/1.2027)*exp(-los**2*1.5)
	  ii=1
	  do while(danta2(ii,1) .lt. los)
	  ii=ii+1
	  enddo
	  radius =danta2(ii+1,2)-(danta2(ii+1,2)-danta2(ii,2))*(danta2(ii+1,1)-los)/(danta2(ii+1,1)-danta2(ii,1))
	  call random_number(los)
	  call random_number(azym)
	  snpos(3) = sqrt(-2.*log(azym))*cos(2.*pi*los)*100.0*pc !266.0*pc		!exp(-|z|/266*pc)
	  call random_number(los)
	  azym = 2.*pi*los
	  snpos(1) = radius*cos(azym)
	  snpos(2) = radius*sin(azym)
!	  snpos(3) = z
	endif
	call MPI_BCAST(snpos, 3, MPI_DOUBLE_PRECISION, 0, comm, ierr)		! przekazywanie pozycji wybuchu SN do odpowiedniego procesu
	if((snpos(1) .ge. xminb-nb*dl(xdim)) .and. (snpos(1) .le. xmaxb+nb*dl(xdim))) then
	if((snpos(2) .ge. yminb-nb*dl(ydim)) .and. (snpos(2) .le. ymaxb+nb*dl(ydim))) then
	if((snpos(3) .ge. zminb-nb*dl(zdim)) .and. (snpos(3) .le. zmaxb+nb*dl(zdim))) then
	  i =nb+int((snpos(1)-xminb)/(xmaxb-xminb)*nxb)
	  j =nb+int((snpos(2)-yminb)/(ymaxb-yminb)*nyb)
	  k =nb+int((snpos(3)-zminb)/(zmaxb-zminb)*nzb)
	  write(*,'(a15,i4.4,i4.4,i4.4)') 'SN II position: ',i,j,k
	  write(*,'(a15,f8.3,1x,f8.3,1x,f8.3)') '     coords:   ',snpos(1),snpos(2),snpos(3)
	  write(*,'(a15,f8.3,1x,f8.3,1x,f8.3)') '     that is:  ',x(i),y(j),z(k)
#ifdef ISO
	  u(2,i,j,k)=u(2,i,j,k)/u(1,i,j,k)*(u(1,i,j,k)+EexplSN)
	  u(3,i,j,k)=u(3,i,j,k)/u(1,i,j,k)*(u(1,i,j,k)+EexplSN)
          u(1,i,j,k)=u(1,i,j,k)+EexplSN
#else /* ISO */
          u(5,i,j,k)=u(5,i,j,k)+EexplSN
#endif /* ISO */
	endif
	endif
	endif
      enddo
    endif

    
    return
  
  end subroutine supernovae_distribution  


end module sn_distr

