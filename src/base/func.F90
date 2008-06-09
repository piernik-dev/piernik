#include "mhd.def"
module func
contains

! Te procedury powinny sie znalezc docelowo w jakims innym module. 

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

!-----------------------------------------------------------------------------

  function pshift(tab,d)
    implicit none
    real, dimension(:,:,:) :: tab
    integer :: d
    integer :: lx,ly,lz
    real, dimension(SIZE(tab,1),SIZE(tab,2),SIZE(tab,3)) :: pshift

    lx = SIZE(tab,1)
    ly = SIZE(tab,2)
    lz = SIZE(tab,3)

    if(d==1) then
      pshift(1:lx-1,:,:) = tab(2:lx,:,:); tab(lx,:,:) = tab(1,:,:)
    else if(d==2) then
      pshift(:,1:ly-1,:) = tab(:,2:ly,:); tab(:,ly,:) = tab(:,1,:)
    else if(d==3) then
      pshift(:,:,1:lz-1) = tab(:,:,2:lz); tab(:,:,lz) = tab(:,:,1)
    else
       write(*,*) 'Dim ill defined in pshift!'
    endif

    return
  end function pshift

  function mshift(tab,d)
    implicit none
    real, dimension(:,:,:) :: tab
    integer :: d
    integer :: lx,ly,lz
    real, dimension(SIZE(tab,1),SIZE(tab,2),SIZE(tab,3)) :: mshift

    lx = SIZE(tab,1)
    ly = SIZE(tab,2)
    lz = SIZE(tab,3)

    if(d==1) then
      mshift(2:lx,:,:) = tab(1:lx-1,:,:); tab(1,:,:) = tab(lx,:,:)
    else if(d==2) then
      mshift(:,2:ly,:) = tab(:,1:ly-1,:); tab(:,1,:) = tab(:,ly,:)
    else if(d==3) then
      mshift(:,:,2:lz) = tab(:,:,1:lz-1); tab(:,:,1) = tab(:,:,lz)
    else
       write(*,*) 'Dim ill defined in mshift!'
    endif

    return
  end function mshift

end module func

!!! MH: procedury "dipol" i "rn_angles"
!!! zostaly przeniesione (z nieiwlkimi przrzerobkami) do modulu sn_sources.F90
!!! w katalogu src/supernovae
