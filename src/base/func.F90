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

#if defined COSM_RAYS || defined PRESS_GRAD_EXCH
#ifndef SPLIT
  real, allocatable, dimension(:,:) :: divv
#endif /* !SPLIT */
#endif /* COSM_RAYS || PRESS_GRAD_EXCH */

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
    use arrays, only : nx,ny,nz,u,idna,imxa,imya,imza,nfluid,divvel
    use grid,   only : dx,dy,dz
    use start,  only : dimensions
    implicit none
!    real, dimension(nx)  :: dvx,dvy,dvz
    real, dimension(nfluid,nx) :: tmp
    integer j,k

    divvel(:,:,:,:) = 0.0
    if(dimensions .eq. '2dxy') then
      k=1
      do j=2,ny-1
        tmp=u(imxa,:,j,k)/u(idna,:,j,k)
        divvel(:,:,j,k) = (eoshift(tmp,1,DIM=2) - eoshift(tmp,-1,DIM=2))/(2.*dx)
        tmp=(u(imya,:,j+1,k)/u(idna,:,j+1,k)-u(imya,:,j-1,k)/u(idna,:,j-1,k))/(2.*dy)
        divvel(:,:,j,k) = divvel(:,:,j,k) + tmp
      enddo
    else if (dimensions .eq. '3d') then
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

#ifdef SPLIT
    use arrays, only : divvel,nx,nfluid
    real,dimension(nfluid,nx) :: divv
    integer j,k
    divv = divvel(:,:,j,k)
#else /* SPLIT */
    use arrays, only : u,idna,imxa,imya,imza,nx,ny,nz,nfluid
    use grid,   only : dx,dy,dz
    use start,  only : dimensions

    implicit none
    real, dimension(nfluid,nx) :: aux, divv
    integer :: j,jm,jp,k,km,kp
    real rfaq

    jp = j+1; jm = j-1
    kp = k+1; km = k-1
    call whichfaq(rfaq,jm,jp,ny)

    aux = u(imxa,:,j,k)/u(idna,:,j,k)
    divv(:,2:nx-1) = 0.5*(aux(:,3:nx)-aux(:,1:nx-2))/dx
    aux = rfaq*(u(imya,:,jp,k)/u(idna,:,jp,k)-u(imya,:,jm,k)/u(idna,:,jm,k))/dy
    divv(:,2:nx-1) = divv(:,2:nx-1) + aux

    if(dimensions == '3d') then
      call whichfaq(rfaq,km,kp,nz)
      aux = rfaq*(u(imza,:,j,kp)/u(idna,:,j,kp)-u(imza,:,j,km)/u(idna,:,j,km))/dz
      divv(:,2:nx-1) = divv(:,2:nx-1) + aux
    endif

   divv(:,1) = divv(:,2) ; divv(:,nx) = divv(:,nx-1)
#endif /* SPLIT */

  end subroutine div_vx

  subroutine div_vy(k,i)

#ifdef SPLIT
    use arrays, only : divvel,ny,nfluid
    real,dimension(nfluid,ny) :: divv
    integer i,k
    divv = divvel(:,i,:,k)
#else /* SPLIT */
    use arrays, only : u,idna,imxa,imya,imza,nx,ny,nz,nfluid
    use grid,   only : dx,dy,dz
    use start,  only : dimensions

    implicit none
    real, dimension(nfluid,ny) :: aux, divv
    integer :: i,im,ip,k,km,kp
    real rfaq

    ip = i+1; im = i-1
    kp = k+1; km = k-1
    call whichfaq(rfaq,im,ip,nx)

    aux = u(imya,i,:,k)/u(idna,i,:,k)
    divv(:,2:ny-1) = 0.5*(aux(:,3:ny)-aux(:,1:ny-2))/dy
    aux = rfaq*(u(imxa,ip,:,k)/u(idna,ip,:,k)-u(imxa,im,:,k)/u(idna,im,:,k))/dx
    divv(:,2:ny-1) = divv(:,2:ny-1) + aux

    if(dimensions == '3d') then
      call whichfaq(rfaq,km,kp,nz)
      aux = rfaq*(u(imza,i,:,kp)/u(idna,i,:,kp)-u(imza,i,:,km)/u(idna,i,:,km))/dz
      divv(:,2:ny-1) = divv(:,2:ny-1) + aux
    endif

   divv(:,1) = divv(:,2) ; divv(:,ny) = divv(:,ny-1)
#endif /* SPLIT */

  end subroutine div_vy

  subroutine  div_vz(j,i)

#ifdef SPLIT
    use arrays, only : divvel,nz,nfluid
    real,dimension(nfluid,nz) :: divv
    integer i,j
    divv = divvel(:,i,j,:)
#else /* SPLIT */
    use arrays, only : u,idna,imxa,imya,imza,nx,ny,nz,nfluid
    use grid,   only : dx,dy,dz
    use start,  only : dimensions

    implicit none
    real, dimension(nfluid,nz) :: aux, divv
    integer :: j,jm,jp,i,im,ip
    real rfaq

    jp = j+1; jm = j-1
    ip = i+1; im = i-1
    call whichfaq(rfaq,im,ip,nx)

    aux = u(imza,i,j,:)/u(idna,i,j,:)
    divv(:,2:nz-1) = 0.5*(aux(:,3:nz)-aux(:,1:nz-2))/dz
    aux = rfaq*(u(imxa,ip,j,:)/u(idna,ip,j,:)-u(imxa,im,j,:)/u(idna,im,j,:))/dx
    divv(:,2:nz-1) = divv(:,2:nz-1) + aux
    call whichfaq(rfaq,jm,jp,ny)
    aux = rfaq*(u(imya,i,jp,:)/u(idna,i,jp,:)-u(imya,i,jm,:)/u(idna,i,jm,:))/dy
    divv(:,2:nz-1) = divv(:,2:nz-1) + aux

    divv(:,1) = divv(:,2) ; divv(:,nz) = divv(:,nz-1)
#endif /* SPLIT */
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
