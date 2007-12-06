#include "mhd.def"


module diagnostics

! Written by: M. Hanasz, December 2005 - January 2006
! Subroutines: "transpose12" and "transpose12" from the original 
!              Pen, Arras & Wong (2003) code



contains


subroutine maxu(comment)
  use start, only  : nstep,proc
  use arrays, only : u,is,ie,js,je,ks,ke
  implicit none
  character comment*(*)

  write(*,*)    
  write(*,*) comment,':    nstep=', nstep, '       proc =',proc
  write(*,*) 'max(u1)=',maxval(u(1,is:ie,js:je,ks:ke)), &
                        maxloc(u(1,is:ie,js:je,ks:ke))
  write(*,*) 'max(u2)=',maxval(u(2,is:ie,js:je,ks:ke)), &
                        maxloc(u(2,is:ie,js:je,ks:ke))
  write(*,*) 'max(u3)=',maxval(u(3,is:ie,js:je,ks:ke)), &
                        maxloc(u(3,is:ie,js:je,ks:ke))
  write(*,*) 'max(u4)=',maxval(u(4,is:ie,js:je,ks:ke)), &
                        maxloc(u(4,is:ie,js:je,ks:ke))
#ifndef ISO
  write(*,*) 'max(u5)=',maxval(u(5,is:ie,js:je,ks:ke)), &
                        maxloc(u(5,is:ie,js:je,ks:ke))
#endif
  write(*,*)    
  
end subroutine maxu  


!------------------------------------------------------------------------------------------
subroutine maxb(comment)
  use start, only  : nstep,proc
  use arrays, only : b,is,ie,js,je,ks,ke
  implicit none
  character comment*(*)

  write(*,*)    
  write(*,*) comment,':    nstep=', nstep
  write(*,*) 'max(b1)=',maxval(b(1,is:ie+1,js:je+1,ks:ke+1)), &
                        maxloc(b(1,is:ie+1,js:je+1,ks:ke+1))
  write(*,*) 'max(b2)=',maxval(b(2,is:ie+1,js:je+1,ks:ke+1)), &
                        maxloc(b(2,is:ie+1,js:je+1,ks:ke+1))
  write(*,*) 'max(b3)=',maxval(b(3,is:ie+1,js:je+1,ks:ke+1)), &
                        maxloc(b(3,is:ie+1,js:je+1,ks:ke+1))
  write(*,*)    
  
end subroutine maxb  

!------------------------------------------------------------------------------------------


subroutine test_divb
  use grid, only   : dx,dy,dz
  use arrays, only : nx,ny,nz,b
  use constants, only : big
  implicit none
  real, allocatable   :: divb(:,:,:)
  allocate (divb(nx,ny,nz))
    
  divb = (eoshift(b(1,:,:,:),shift=1,dim=1,boundary=big) -b(1,:,:,:))*dy*dz &
        +(eoshift(b(2,:,:,:),shift=1,dim=2,boundary=big) -b(2,:,:,:))*dx*dz &
        +(eoshift(b(3,:,:,:),shift=1,dim=3,boundary=big) -b(3,:,:,:))*dx*dy 

  write(*,*)  
  write(*,*) 'MAX(|div B|) = ', maxval(abs(divb(2:nx-1,2:ny-1,2:nz-1)))
  write(*,*)

  deallocate (divb)

end subroutine test_divb  

!==========================================================================================
!      subroutine fpesetup
!      use iflport
!      implicit none
!      interface
!         subroutine fpe_handler(signo,siginfo)
!            integer(4), intent(in) :: signo, siginfo
!         end subroutine
!      end interface
!      integer ir
!      ir = ieee_handler('set', 'overflow',fpe_handler)
!      end subroutine fpesetup

!      subroutine fpe_handler(sig, code)
!      use iflport
!      implicit none
!      integer sig, code
!      print *, 'fpesetup: code =', code, '  sig=', sig
!      if(code .eq.fpe$zerodivide) print *, 'occured divide by zero'
!      if(code .eq.fpe$overflow)   print *, 'occured overflow'
!      if(code .eq.fpe$underflow)  print *, 'occured underflow'
!      if(code .eq.fpe$invalid)    print *, 'occured invalid'
!      end subroutine fpe_handler

end module diagnostics
