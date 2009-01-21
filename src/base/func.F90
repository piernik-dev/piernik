! $Id$
#include "piernik.def"
!>
!! \brief Module that contains unclassified functions
!!
!! This module should be empty. Every function or subroutine placed here belong
!! elsewhere. We are yet unsure where to put them.
!! \todo Move all structures elsewhere
!! \warning Procedures \a dipol and \a rn_angles were moved to sn_sources.F90
!<
module func

implicit none

  contains

! Te procedury powinny sie znalezc docelowo w jakims innym module.

!>
!! \brief Function that evolves EMFs in time
!! \param vh velocity perpendicular to #b
!! \param vg velocity perpendicular to #b
!! \param b one component of magnetic field
!! \param dt timestep
!<
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

!>
!! \brief Function makes one-cell,foreward circular shift of 3D array in any direction
!! \param tab input array
!! \param d shift's direction, where 1,2,3 corresponds to \a x,\a y,\a z respectively
!! \return real, dimension(SIZE(tab,1),SIZE(tab,2),SIZE(tab,3))
!!
!! Functions was written in order to significantly improve
!! the performance at the cost of the flexibility of original \p CSHIFT.
!<
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
      pshift(1:lx-1,:,:) = tab(2:lx,:,:); pshift(lx,:,:) = tab(1,:,:)
    else if(d==2) then
      pshift(:,1:ly-1,:) = tab(:,2:ly,:); pshift(:,ly,:) = tab(:,1,:)
    else if(d==3) then
      pshift(:,:,1:lz-1) = tab(:,:,2:lz); pshift(:,:,lz) = tab(:,:,1)
    else
       write(*,*) 'Dim ill defined in pshift!'
    endif

    return
  end function pshift

!>
!! \brief Function makes one-cell, backward circular shift of 3D array in any direction
!! \param tab input array
!! \param d shift's direction, where 1,2,3 corresponds to \a x,\a y,\a z respectively
!! \return real, dimension(SIZE(tab,1),SIZE(tab,2),SIZE(tab,3))
!!
!! Functions was written in order to significantly improve
!! the performance at the cost of the flexibility of original \p CSHIFT.
!<
  function mshift(tab,d)
    implicit none
    real, dimension(:,:,:) :: tab
    integer :: d
    integer :: lx,ly,lz
    real, dimension(SIZE(tab,1) , SIZE(tab,2) , SIZE(tab,3)) :: mshift

    lx = SIZE(tab,1)
    ly = SIZE(tab,2)
    lz = SIZE(tab,3)

    if(d==1) then
      mshift(2:lx,:,:) = tab(1:lx-1,:,:); mshift(1,:,:) = tab(lx,:,:)
    else if(d==2) then
      mshift(:,2:ly,:) = tab(:,1:ly-1,:); mshift(:,1,:) = tab(:,ly,:)
    else if(d==3) then
      mshift(:,:,2:lz) = tab(:,:,1:lz-1); mshift(:,:,1) = tab(:,:,lz)
    else
       write(*,*) 'Dim ill defined in mshift!'
    endif

    return
  end function mshift

!-----------------------------------------------------------------------------

#if defined COSM_RAYS || defined PRESS_GRAD_EXCH
  subroutine div_v
    use fluidindex, only : idna,imxa,imya,imza
    use arrays, only : nx,ny,nz,u,nfluid,divvel
    use grid,   only : dx,dy,dz,nzd
    implicit none
!    real, dimension(nx)  :: dvx,dvy,dvz
    real, dimension(nfluid,nx) :: tmp
    integer j,k

    divvel(:,:,:,:) = 0.0
    if(nzd == 1) then
      k=1
      do j=2,ny-1
        tmp=u(imxa,:,j,k)/u(idna,:,j,k)
        divvel(:,:,j,k) = (eoshift(tmp,1,DIM=2) - eoshift(tmp,-1,DIM=2))/(2.*dx)
        tmp=(u(imya,:,j+1,k)/u(idna,:,j+1,k)-u(imya,:,j-1,k)/u(idna,:,j-1,k))/(2.*dy)
        divvel(:,:,j,k) = divvel(:,:,j,k) + tmp
      enddo
    else if (nzd /= 1) then
      do k=2,nz-1
        do j=2,ny-1
          tmp=u(imxa,:,j,k)/u(idna,:,j,k)
          divvel(:,:,j,k) = (eoshift(tmp,1,DIM=2) - eoshift(tmp,-1,DIM=2))/(2.*dx)
          tmp=(u(imya,:,j+1,k)/u(idna,:,j+1,k)-u(imya,:,j-1,k)/u(idna,:,j-1,k))/(2.*dy)
          divvel(:,:,j,k) = divvel(:,:,j,k) + tmp
          tmp=(u(imza,:,j,k+1)/u(idna,:,j,k+1)-u(imza,:,j,k-1)/u(idna,:,j,k-1))/(2.*dz)
          divvel(:,:,j,k) = divvel(:,:,j,k) + tmp
        enddo
      enddo
    endif
  end subroutine div_v


  subroutine div_vx(k,j)

    use arrays, only : divvel,nx,nfluid
    real,dimension(nfluid,nx) :: divv
    integer j,k
    divv = divvel(:,:,j,k)

  end subroutine div_vx

  subroutine div_vy(k,i)

    use arrays, only : divvel,ny,nfluid
    real,dimension(nfluid,ny) :: divv
    integer i,k
    divv = divvel(:,i,:,k)

  end subroutine div_vy

  subroutine  div_vz(j,i)

    use arrays, only : divvel,nz,nfluid
    real,dimension(nfluid,nz) :: divv
    integer i,j
    divv = divvel(:,i,j,:)
  end subroutine div_vz

  subroutine whichfaq(faq,i,j,n)
  implicit none
  real faq
  integer i,j,n

  faq = 0.5
  if(i .eq. 0) then
   i=1
   faq=1.0
  endif
  if(j-1 .eq. n) then
   j=n
   faq=1.0
  endif

  end subroutine whichfaq

#endif /* COSM_RAYS || PRESS_GRAD_EXCH */


end module func

!!! MH: procedury "dipol" i "rn_angles"
!!! zostaly przeniesione (z nieiwlkimi przrzerobkami) do modulu sn_sources.F90
!!! w katalogu src/supernovae
