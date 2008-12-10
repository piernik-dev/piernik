! $Id$
#include "piernik.def"

!>
!! \brief Module that collects source components specified for each fluid
!!
!! This module should not be changed by user.
!<
module allsources
#ifdef IONIZED
  use ionizeds, only : src_ionized
#endif /* IONIZED */
#ifdef NEUTRAL
  use neutrals, only : src_neutral
#endif /* NEUTRAL */
#ifdef MOLECULAR
  use moleculars, only : src_molecular
#endif /* MOLECULAR */
#ifdef DUST
  use dusts, only : src_dust
#endif /* DUST */
#ifdef COSM_RAYS
  use cosmic_rays, only : src_cr
#endif /* COSM_RAYS */

contains

!================GENERAL SUBROUTINE OF ALLSOURCES MODULE=======================
!>
!! \brief Subroutine which calls source subroutines for characteristic for each fluid or each couple of fluids
!! \param Duus timestep correction contribution from fluid source terms
!! \param uu currently used fluid table
!! \param bb magnetic field x,y,z-components table
!! \param sweep current sweep
!! \param i1 number of the row in the first direction after current sweep direction
!! \param i2 number of the row in the first direction before current sweep direction
!! \param n number of cells in the current sweep
!! \param dt current time step
!<
subroutine all_src(Duus,uu,bb,sweep,i1,i2,n,dt)
    use arrays, only : nu
    implicit none
    integer               :: i1,i2,n
    real                  :: dt
    real, dimension(nu,n) :: uu,Duus,Duuss
    real, dimension(3,n)  :: bb
    character sweep*6

    Duus=0.0
#ifdef COLLISIONS
    call collisions(Duuss,uu,sweep,i1,i2,n,dt)
    Duus=Duus+Duuss
#endif /* COLLISIONS */

#ifdef PRESS_GRAD_EXCH
    call press_grad_exchange(Duuss,uu,bb,sweep,i1,i2,n,dt)
    Duus=Duus+Duuss
#endif /* PRESS_GRAD_EXCH */

#ifdef IONIZED
   call src_ionized !(Duuss,uu,sweep,i1,i2,n,dt)
!   Duus=Duus+Duuss
#endif /* IONIZED */

#ifdef NEUTRAL
   call src_neutral !(Duuss,uu,sweep,i1,i2,n,dt)
!   Duus=Duus+Duuss
#endif /* NEUTRAL */

#ifdef MOLECULAR
   call src_molecular !(Duuss,uu,sweep,i1,i2,n,dt)
!   Duus=Duus+Duuss
#endif /* MOLECULAR */

#ifdef DUST
   call src_dust !(Duuss,uu,sweep,i1,i2,n,dt)
!   Duus=Duus+Duuss
#endif /* DUST */

#ifdef COSM_RAYS
   call src_cr !(Duuss,uu,sweep,i1,i2,n,dt)
!   Duus=Duus+Duuss
#endif /* COSM_RAYS */

#ifdef ANY_COOLING
   call all_cooling(eint,n)
#endif /* ANY_COOLING */

#ifdef KEPLER_SUPPRESSION
    call kepler_suppression(Duuss,uu,sweep,i1,i2,n,dt)
    Duus=Duus+Duuss
#endif /* KEPLER_SUPPRESSION */

end subroutine all_src

!===============END OF GENERAL SUBROUTINE OF ALLSOURCES MODULE=====================

!==========================ELEMENTARY SUBROUTINES==================================

#ifdef ANY_COOLING
!>
!! \brief Subroutine which calls cooling subroutines if specified.
!! \param eint internal energy density
!! \param n number of cells in the current sweep
!<
subroutine all_cooling(eint,n)
  use arrays, only : nadiab
  implicit none
  integer n
  real, dimension(nadiab,n)    :: eint

#if !defined ISO && defined SIMPLE_COOL
    call simple_cool(eint,n)
#endif /* !ISO && SIMPLE_COOL */
end subroutine all_cooling
#endif /* ANY_COOLING */

#ifdef COLLISIONS
!>
!! \brief Subroutine which computes collisions/friction between each pair of fluids.
!! \param Duus timestep correction contribution from fluid source terms
!! \param uu currently used fluid table
!! \param sweep current sweep
!! \param i1 number of the row in the first direction after current sweep direction
!! \param i2 number of the row in the first direction before current sweep direction
!! \param n number of cells in the current sweep
!! \param dt current time step
!<
subroutine collisions(Duus,uu,sweep,i1,i2,n,dt)
  use arrays, only : nu,x,y,z,idna,imxa,imya,nfluid
  use start,  only : maxxyz,collfaq
#ifndef ISO
  use arrays, only : iena,fadiab
#endif /* !ISO */
  implicit none
  character sweep*6
  integer n,i1,i2,ifl,jfl
  real, dimension(nfluid,nfluid,n) :: flch
  real, dimension(nfluid,n)        :: colls,velcoor
  real, dimension(nu,n)            :: Duus,uu
  real, dimension(maxxyz)          :: r1,r2
  real    :: a1,dt
  integer :: rend

    Duus=0.0
!-----collisions between fluids------
    select case(sweep)
       case('xsweep')
         a1    = y(i1)
         rend  = size(x)
         r1(1:rend) = x(:)                            ! r1   max(size(x),size(y),size(z))
         r2(1:rend) = a1 / (r1(1:rend)*r1(1:rend) + a1 * a1)
       case('ysweep')
         a1    = x(i2)
         rend  = size(y)
         r1(1:rend) = y(1:rend)
         r2(1:rend) = a1 / (r1(1:rend)*r1(1:rend) + a1 * a1)
       case('zsweep')
         a1    = 1.0
         rend  = size(z)
         r1(1:rend) = 0.0
         r2(1:rend) = 1.0
    end select

#ifdef COLLS_MODIF
    velcoor(:,:) = (uu(imxa,:)*a1 - uu(imya,:)*spread(r1(:),1,nfluid))/uu(idna,:)*spread(r2(:),1,nfluid)
#endif /* COLLS_MODIF */
    do ifl=1,nfluid
      do jfl=1,nfluid
        if(ifl .ne. jfl) then
#ifdef COLLS_PART
          flch(ifl,jfl,:)= uu(idna(2),:) &
#else /* COLLS_PART */
          flch(ifl,jfl,:)= uu(idna(ifl),:) * uu(idna(jfl),:) &
#endif /* COLLS_PART */
#ifdef COLLS_MODIF
              *( velcoor(jfl,:) - velcoor(ifl,:) )
#else /* COLLS_MODIF */
              *( uu(imxa(jfl),:)/uu(idna(jfl),:) - uu(imxa(ifl),:)/uu(idna(ifl),:))
#endif /* COLLS_MODIF */
        else
          flch(ifl,jfl,:)=0.0
        endif
      enddo
    enddo
    do ifl=1,nfluid
      colls(ifl,:)=collfaq*sum(flch(ifl,:,:),1)
    enddo
#ifndef ISO
!    Duus(iena(fadiab),:)=(uu(imxa(fadiab),:)*colls*dt+0.5*(colls*dt)**2)/uu(idna(fadiab),:)
!    Duus(iena(fadiab),:)=uu(imxa(fadiab),:)**2/uu(idna(fadiab),:)*colls*dt
    Duus(iena(fadiab),:)=uu(imxa(fadiab),:)*colls(fadiab,:)*dt
#endif /* ISO */
    Duus(imxa,:)=colls*dt

!------------------------------------
end subroutine collisions
#endif /* COLLISIONS */

#ifdef KEPLER_SUPPRESSION
!>
!! \brief Subroutine which helps to maintain rotation of fluids close to boundary regions.
!! \param Duus timestep correction contribution from fluid source terms
!! \param uu currently used fluid table
!! \param sweep current sweep
!! \param i1 number of the row in the first direction after current sweep direction
!! \param i2 number of the row in the first direction before current sweep direction
!! \param n number of cells in the current sweep
!! \param dt current time step
!<
subroutine kepler_suppression(Duus,uu,sweep,i1,i2,n,dt)
    use arrays, only : alfsup,omx0,omy0,nu,nfluid,x,y,z,nx,ny,idna,imxa,imya
    use start,  only : maxxyz
#ifndef ISO
  use arrays, only : nadiab,fadiab,iena
#endif /* !ISO */
  implicit none
  integer n,i,j,k,i1,i2,ifl,ni1,ni2
  character sweep*6
  real,dimension(nu,n)      :: Duus,uu
  real, dimension(nfluid,n) :: kplrsup,velcoor,vel0
  real, dimension(maxxyz)   :: r1,r2
  real :: a1,dt
  integer :: rend,ii

    Duus=0.0
    select case(sweep)
      case('xsweep')
        a1    = y(i1)
        rend  = size(x)
        r1(1:rend) = x(:)                            ! r1   max(size(x),size(y),size(z))
        r2(1:rend) = a1 / (r1(1:rend)*r1(1:rend) + a1 * a1)
        vel0(:,:)  = omx0(:,:,i1,i2)
      case('ysweep')
        a1    = x(i2)
        rend  = size(y)
        r1(1:rend) = y(1:rend)
        r2(1:rend) = a1 / (r1(1:rend)*r1(1:rend) + a1 * a1)
        vel0(:,:)  = omy0(:,i2,:,i1)
      case('zsweep')
        a1    = 1.0
        rend  = size(z)
        r1(1:rend) = 0.0
        r2(1:rend) = 1.0
    end select

    if(sweep .ne. 'zsweep') then
    do ifl=1,nfluid
#ifdef KEPL_SUPP_SIMX
        velcoor(ifl,:)=uu(imxa(ifl),:)/uu(idna(ifl),:)
        kplrsup(ifl,:)=-alfsup(:,i1)*(velcoor(ifl,:)-vel0(ifl,:)) !*uu(idna(ifl),:)
#else /* KEPL_SUPP_SIMX */
        velcoor(ifl,:)=(uu(imxa(ifl),:)*a1-uu(imya(ifl),:)*r1)/uu(idna(ifl),:)*r2
        kplrsup(ifl,:)=-alfsup(:,i1)*(velcoor(ifl,:)-vel0(ifl,:))!*uu(idna(ifl),:)
#endif /* KEPL_SUPP_SIMX */
    enddo
    else
      kplrsup(:,:)=0.0
    endif
#ifndef ISO
!    Duus(iena(fadiab),:)=(uu(imxa(fadiab),:)*kplrsup*uu(idna(fadiab),:)*dt &
!            +0.5*(kplrsup*uu(idna(fadiab),:)*dt)**2)/uu(idna(fadiab),:)
    Duus(iena(fadiab),:)=kplrsup(fadiab,:)*uu(imxa(fadiab),:)*dt
#endif /* ISO */
    Duus(imxa,:)=kplrsup(:,:)*uu(idna,:)*dt
!    if((i1 .eq. 50) .and. (sweep .eq. 'xsweep')) then
!    write(*,*) Duus(imxa(1),10),kplrsup(1,10),uu(idna(1),10),omx0(1,10,50,1),alfsup(10,50),uu(imxa(1),10)
!    do ii=50,62
!    if((sweep .eq. 'xsweep') .and. (abs(uu(imxa(1),ii)) .ge. abs(omx0(1,ii,i1,i2)))) then
!    if(i1 .eq. 40) then
!    write(*,*) 'omx0 smaller!: ',sweep,ii,i1,i2
!    write(*,*) Duus(imxa(1),ii),kplrsup(1,i1),uu(idna(1),ii),omx0(1,ii,i1,i2),alfsup(ii,i1),uu(imxa(1),ii)
!    endif
!    endif
!    if((sweep .eq. 'ysweep') .and. (abs(uu(imxa(1),ii)) .ge. abs(omy0(1,i2,ii,i1)))) then
!    write(*,*) 'omy0 smaller!: ',sweep,i2,ii,i1
!    write(*,*) Duus(imxa(1),ii),kplrsup(1,ii),uu(idna(1),ii),omx0(1,i2,ii,i1),alfsup(i2,ii),uu(imxa(1),ii)
!    endif
!    enddo

end subroutine kepler_suppression
#endif /* KEPLER_SUPPRESSION */

#ifdef PRESS_GRAD_EXCH
!>
!! \brief This subroutine guarantees an exchange of PdivV and gradP/rho between fluids.
!! \param Duus timestep correction contribution from fluid source terms
!! \param uu currently used fluid table
!! \param sweep current sweep
!! \param i1 number of the row in the first direction after current sweep direction
!! \param i2 number of the row in the first direction before current sweep direction
!! \param n number of cells in the current sweep
!! \param dt current time step
!! \todo To be extensively tested.
!! \warning NOT READY YET !!!
!<
subroutine press_grad_exchange(Duus,uu,bb,sweep,i1,i2,n,dt)
!    use arrays, only : alfsup,omx0,omy0,nu,nfluid,x,y,z,nx,ny,idna,imxa,imya
    use arrays, only : nfluid,nu,divvel,idna,imxa,imya,imza,ibx,iby,ibz
    use start,  only : gamma,c_si
    use grid,   only : dx
#ifndef ISO
  use arrays, only : nadiab,fadiab,iena
#endif /* !ISO */
  implicit none
  integer n,i,j,k,i1,i2,ifl,jfl,ni1,ni2
  character sweep*6
  real, dimension(nu,n)     :: Duus,uu
  real, dimension(3,n)      :: bb
  real, dimension(n)        :: emag
  real, dimension(nfluid,n) :: divv,pressgrad,press
#ifndef ISO
  real, dimension(nadiab,n) :: ekin,eint
#endif /* ISO */
  real dt


    select case (sweep)
      case('xsweep')
        divv = divvel(:,:,i1,i2)
      case('ysweep')
        divv = divvel(:,i2,:,i1)
      case('zsweep')
        divv = divvel(:,i1,i2,:)
    end select
    press(:,:) = uu(idna,:)*c_si**2
#ifdef DUST
    press(fdust,:) = 0.0
#endif /* DUST */
#ifndef ISO
    ekin(:,:) = 0.5*(uu(imxa(fadiab),:)*uu(imxa(fadiab),:) &
                        +uu(imya(fadiab),:)*uu(imya(fadiab),:) &
                        +uu(imza(fadiab),:)*uu(imza(fadiab),:))/uu(idna(fadiab),:)
    emag = 0.5*(bb(ibx,:)*bb(ibx,:) + bb(iby,:)*bb(iby,:) + bb(ibz,:)*bb(ibz,:))
    eint = uu(iena(fadiab),:)-ekin-spread(emag,1,nadiab)
    press(fadiab,:) = eint*spread(gamma(fadiab)-1.0,2,n)
#endif /* !ISO */
    pressgrad(:,2:n-1) = (press(:,3:n)-press(:,1:n-2))/2.0/dx
    pressgrad(:,1) = pressgrad(:,2)
    pressgrad(:,n) = pressgrad(:,n-1)
    do ifl=1,nfluid
      do jfl=1,nfluid
        if(ifl .ne. jfl) then
          Duus(imxa(ifl),:)   =-pressgrad(jfl,:)/uu(idna(ifl),:)*dt
#ifndef ISO
          if((ifl .le. nadiab) .and. (jfl .le. nadiab)) then
            Duus(iena(ifl,:)=-press(jfl,:)/uu(idna(ifl,:)*divv(ifl,:)*dt
          endif
#endif /* ISO */
        endif
      enddo
    enddo

end subroutine press_grad_exchange
#endif /* PRESS_GRAD_EXCH */

#if !defined ISO && defined SIMPLE_COOL
subroutine simple_cool(eint,n)
    use arrays, only : nadiab
    use start,  only : dt,tauc
    implicit none
    integer n
    real, dimension(nadiab,n)    :: eint

    eint = eint*(1.0-0.5*dt/tauc)
end subroutine simple_cool
#endif /* !ISO && SIMPLE_COOL */


!==========================END OF ELEMENTARY SUBROUTINES====================================

end module allsources
