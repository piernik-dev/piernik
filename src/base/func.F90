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
   subroutine namelist_errh(errh,nm)
      implicit none
      integer, intent(in) :: errh
      character(len=*), intent(in) :: nm
      
      select case (errh)
         case (19)
            write(*,*) "severe (19): Invalid reference to variable in ",trim(nm), " namelist"
            write(*,*) "One of the following conditions occurred: "
            write(*,*) "    * The variable was not a member of the namelist group."
            write(*,*) "    * An attempt was made to subscript a scalar variable."
            write(*,*) "    * A subscript of the array variable was out-of-bounds."
            write(*,*) "    * An array variable was specified with too many or too few subscripts for the variable."
            write(*,*) "    * An attempt was made to specify a substring of a noncharacter variable or array name."
            write(*,*) "    * A substring specifier of the character variable was out-of-bounds."
            write(*,*) "    * A subscript or substring specifier of the variable was not an integer constant."
            write(*,*) "    * An attempt was made to specify a substring by using an unsubscripted array variable."
            stop
         case (-1)
            write(*,*) "Namelist: ",trim(nm)," not found in problem.par"
            stop
         case (5010)
            write(*,*) "One of the variables found in problem.par doesn't belong to ",trim(nm), " namelist"
            stop
         case (239)
            write(*,*) "One of the variables found in problem.par doesn't belong to ",trim(nm), " namelist"
            stop
         case (0)
         case default
            write(*,*) 'Unknown error (', errh,') in namelist ',trim(nm)
            stop
      endselect

   end subroutine namelist_errh

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
    real    :: w,dwp,dwm,dw,v

    tvdb_emf = 0.0

    do i = lbound(vh,1)+2, ubound(vh,1)-3
       ip  = i  + 1
       ipp = ip + 1
       im  = i  - 1
       v   = vh(i)
       if (v .gt. 0.) then
         w=vg(i)*b(i)
         dwp=(vg(ip)*b(ip)-w)*0.5
         dwm=(w-vg(im)*b(im))*0.5
       else
         w=vg(ip)*b(ip)
         dwp=(w-vg(ipp)*b(ipp))*0.5
         dwm=(vg(i)*b(i)-w)*0.5
       end if
       dw=0.0
       if(dwm*dwp > 0.0) dw=2.0*dwm*dwp/(dwm+dwp)
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

#ifdef COSM_RAYS 
  subroutine div_v(ifluid)
    use mpi_setup
    use initionized, only : idni,imxi,imyi,imzi
    use fluidindex,  only : nfluid,iarr_all_dn,iarr_all_mx,iarr_all_my,iarr_all_mz
    use grid,        only : nx,ny,nz
    use grid,        only : dx,dy,dz,nxd,nyd,nzd
    use arrays,      only : u,divvel
    implicit none
!    real, dimension(nx)  :: dvx,dvy,dvz
!    real, dimension(nfluid,nx) :: tmp
    real, dimension(nx) :: vx
    real, dimension(ny) :: vy
    real, dimension(nz) :: vz
    integer             :: i,j,k,ifluid
    integer             :: idnf,imxf,imyf,imzf 
    
    
    idnf = iarr_all_dn(ifluid) 
    imxf = iarr_all_mx(ifluid) 
    imyf = iarr_all_my(ifluid) 
    imzf = iarr_all_mz(ifluid) 

    divvel(:,:,:) = 0.0

!    if(nyd == 1 .or. nxd == 1) thenddd
!         write(*,*) "div_v @ func.F90 does not support 'nyd' and/or 'nxd' == 1"
!        call mpistop
!    endif

    if(nxd /= 1) then
       do k = 1, nz
         do j = 1, ny
            vx = u(imxf,:,j,k) / u(idnf,:,j,k)
            divvel(2:nx-1,j,k) = ( vx(3:nx) - vx(1:nx-2) )  / (2.*dx)
         enddo
       enddo
       divvel(1,:,:) = divvel(2,:,:); divvel(nx,:,:) = divvel(nx-1,:,:) ! for sanity
    endif
           
    if(nyd /= 1) then
       do k = 1, nz
         do i = 1, nx
            vy = u(imxf,i,:,k) / u(idnf,i,:,k)
            divvel(i,2:ny-1,k) = ( vy(3:ny) - vy(1:ny-2) )  / (2.*dy)
         enddo
       enddo
       divvel(:,1,:) = divvel(:,2,:); divvel(:,ny,:) = divvel(:,ny-1,:) ! for sanity
    endif
           
    if(nzd /= 1) then
       do j = 1, ny
         do i = 1, nx
            vz = u(imxf,:,j,k) / u(idnf,:,j,k)
            divvel(i,j,2:nz-1) = ( vz(3:nz) - vz(1:nz-2) )  / (2.*dz)
         enddo
       enddo
       divvel(:,:,1) = divvel(:,:,2); divvel(:,:,nz) = divvel(:,:,nz-1) ! for sanity
    endif
           

!    if(nzd == 1) then
!      k=1
!      do j=2,ny-1
!        tmp = u(imxf,:,j,k)/u(idnf,:,j,k)
!        divvel(:,j,k) = (eoshift(tmp,1,DIM=1) - eoshift(tmp,-1,DIM=1))/(2.*dx)
!        tmp = (u(imyf,:,j+1,k)/u(idnf,:,j+1,k)-u(imyf,:,j-1,k)/u(idnf,:,j-1,k))/(2.*dy)
!        divvel(:,j,k) = divvel(:,j,k) + tmp
!      enddo
!    else if (nzd /= 1) then
!      do k=2,nz-1
!        do j=2,ny-1
!          tmp =   u(imxf,:,j,k)/u(idnf,:,j,k)
!          divvel(:,j,k) = (eoshift(tmp,1,DIM=1) - eoshift(tmp,-1,DIM=1))/(2.*dx)
!          tmp = ( u(imyf,:,j+1,k)/u(idnf,:,j+1,k)-u(imyf,:,j-1,k)/u(idnf,:,j-1,k) )/(2.*dy)
!          divvel(:,j,k) = divvel(:,j,k) + tmp
!          tmp = ( u(imzf,:,j,k+1)/u(idnf,:,j,k+1)-u(imzf,:,j,k-1)/u(idnf,:,j,k-1) )/(2.*dz)
!          divvel(:,j,k) = divvel(:,j,k) + tmp
!        enddo
!      enddo
!    endif
  end subroutine div_v


  subroutine div_vx(k,j)

    use grid,        only : nx
    use fluidindex,  only : nfluid
    use arrays,      only : divvel
    
    real,dimension(nx) :: divv
    integer j,k
    
    divv = divvel(:,j,k)

  end subroutine div_vx

  subroutine div_vy(k,i)

    use grid,        only : ny
    use fluidindex,  only : nfluid
    use arrays, only      : divvel
    
    real,dimension(ny) :: divv
    integer i,k
    
    divv = divvel(i,:,k)

  end subroutine div_vy

  subroutine  div_vz(j,i)

    use grid,        only : nz
    use fluidindex,  only : nfluid
    use arrays, only      : divvel
    
    real,dimension(nz) :: divv
    integer i,j
    
    divv = divvel(i,j,:)
    
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

#endif /* COSM_RAYS  */


end module func

!!! MH: procedury "dipol" i "rn_angles"
!!! zostaly przeniesione (z nieiwlkimi przrzerobkami) do modulu sn_sources.F90
!!! w katalogu src/supernovae
