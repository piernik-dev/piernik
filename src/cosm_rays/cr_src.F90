#include "mhd.def"

! Written by M. Hanasz, 2003
! Adapted for this code by M. Hanasz, November 2007


module cr_src


  use arrays
  use constants
  use grid
  use start
#ifdef SHEAR  
  use shear
#endif SHEAR  
  use fluid_boundaries, only : compute_u_bnd 
  
  implicit none
  real xsn,ysn,zsn
  integer, save :: nsn
  real, save    :: dt_sn_prev, ecr_supl, decr_supl   
  
  real gasdev
  real,    save :: gset             
  integer, save :: irand, iset
  

 contains
  
 
  subroutine ran_sncr
    integer i,j,k, ipm, jpm
    integer jsn,jremap 
    real epsi,epso
    real ysna, ysni, ysno,dysn 
    real decr, part_cr_sn      
    real rand(4), znorm
    real dt_sn, t_dw1            
      
    dt_sn = 1./(f_sn+small) * t_arm/(t_dw+small)   
      
    t_dw1  = mod(t,t_dw)
      
    if(t_dw1 .lt. t_arm) then
      part_cr_sn = amp_ecr_sn * 2.*dt/(dt_sn+small)
    else
      part_cr_sn = 0.0
    endif

    if(dt_sn_prev .gt. dt_sn) then
      dt_sn_prev = 0.0
      decr_supl = 0.0
      nsn = nsn+1
      
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
    endif
 
!    if(nsn .GT. 0) then  

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

     
      do k=1,nz
        do j=1,ny
          do i=1,nx
 
            do ipm=-1,1

              if(ipm .eq. -1) ysna = ysno
              if(ipm .eq.  0) ysna = ysn
              if(ipm .eq.  1) ysna = ysni

              do jpm=-1,1
                decr = part_cr_sn * ethu  &
                       * EXP(-((x(i)-xsn+real(ipm)*Lx)**2  &
                       + (y(j)-ysna+real(jpm)*Ly)**2  &
                       + (z(k)-zsn)**2)/r_sn**2)  

                u(iecr,i,j,k) = u(iecr,i,j,k) + decr		
 
                if(((i .ge. is) .and. (i .le. ie)) .and. &
     	           ((j .ge. js) .and. (j .le. je)) .and. &
     	           ((k .ge. ks) .and. (k .le. ke))) then
                  decr_supl = decr_supl + decr*dvol
                  ecr_supl = ecr_supl + decr*dvol
                endif		
 
              enddo ! jpm
            enddo ! ipm
   
          enddo ! i
        enddo ! j
      enddo ! k
!    endif ! nsn
     
    dt_sn_prev = dt_sn_prev + 2.*dt 

    call compute_u_bnd
      
    return

  end subroutine ran_sncr


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
  

end module cr_src
