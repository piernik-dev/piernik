! $Id$
#include "piernik.def"
#ifdef SSP
#define UL_TAB ul
#define UR_TAB ur
#endif /* SSP */

#ifdef ORIG
#define UL_TAB ul0
#define UR_TAB ur0
#endif /* ORIG */

#ifdef SPLIT
#define VIBJ_TAB vibj1
#define U_TABL uu(:,:)
#define U_IDNA uu(idna,:)
#define U_IMXA uu(imxa,:)
#define U_IMYA uu(imya,:)
#define U_IMZA uu(imza,:)
#define U_IENA uu(iena,:)
#define U_IECR uu(iecr,:)
#define B_IBX bb(ibx,:)
#define B_IBY bb(iby,:)
#define B_IBZ bb(ibz,:)
#define VZ_LIM_DN uu(idna,:)
#define VZ_LIM_VZ uu(imxa,:)
#else /* SPLIT */
#define VIBJ_TAB vibj
#define U_TABL u(:,:,:,:)
#define U_IDNA u(idna,:,:,:)
#define U_IMXA u(imxa,:,:,:)
#define U_IMYA u(imya,:,:,:)
#define U_IMZA u(imza,:,:,:)
#define U_IENA u(iena,:,:,:)
#define U_IECR u(iecr,:,:,:)
#define B_IBX b(ibx,:,:,:)
#define B_IBY b(iby,:,:,:)
#define B_IBZ b(ibz,:,:,:)
#define VZ_LIM_DN u(idna,:)
#define VZ_LIM_VZ u(imza,:)
#endif /* SPLIT */

#if defined ORIG && defined SPLIT
#define DIORIG ,di
#define B1ORIG b1
#else /* ORIG && SPLIT */
#define DIORIG
#define B1ORIG b
#endif /* ORIG && SPLIT */

module tv   ! both both
  contains

  subroutine tvdb(vibj,b,vg,n,dt DIORIG)
#ifdef ORIG
    use constants, only : big
#ifndef SPLIT
    use start, only     : istep
#endif /* SPLIT */
#endif /* ORIG */
    use func, only      : tvdb_emf
    implicit none
    integer, intent(in) :: n
    real, intent(in)    :: dt DIORIG
    real, dimension(n)  :: vibj,b,vg
! locals
    real, dimension(n)  :: vh
#if defined SPLIT && defined ORIG
    real, dimension(n)  :: b1,vibj1
    real :: dti
#endif /* SPLIT && ORIG */

  ! unlike the B field, the vibj lives on the right cell boundary
    vh = 0.0
    vh(1:n-1) =(vg(1:n-1)+ vg(2:n))*0.5;     vh(n) = vh(n-1)

#ifdef ORIG
#ifdef SPLIT
    dti = dt/di
#endif /* SPLIT */
#ifndef SPLIT
    if(istep .eq. 1) then
#endif /* SPLIT */
      where(vh > 0.)
        VIBJ_TAB=b*vg
      elsewhere
        VIBJ_TAB=eoshift(b*vg,1,boundary=big)
      end where
#ifdef SPLIT
    b1(2:n) = b(2:n) -(vibj1(2:n)-vibj1(1:n-1))*dti*0.5;    b1(1) = b(2)
#endif /* SPLIT */
#ifndef SPLIT
      vibj = vibj*dt

    else if(istep .eq. 2) then
#endif /* SPLIT */
#endif /* ORIG */
    vibj = tvdb_emf(vh,vg,B1ORIG,dt)
#ifndef SPLIT
#ifdef ORIG
    else
      stop
    endif
#endif /* ORIG */
#endif /* SPLIT */

  end subroutine tvdb

  subroutine relaxing_tvd(Duu,uu,bb,sweep,i1,i2,dx,n,dt)
    use arrays,  only : idna,imxa,imya,imza,ibx,iby,ibz,nu
    use fluxes,  only : flux_limit,grav_limit,mhdflux
#if defined SPLIT || defined ORIG
    use start,   only : istep
    use arrays,  only : ul0,ur0
#endif /* SPLIT || ORIG */

#ifndef ISO
    use arrays,  only : iena
#if defined GRAV || defined TIDAL_FORCES
    use arrays,  only : nadiab,fadiab
#endif /* GRAV || TIDAL_FORCES */
#endif /* ISO */
#ifdef GRAV
    use gravity, only : grav_pot2accel
    use arrays,  only : nfluid
#endif /* GRAV */
#ifdef GLOBAL_FR_SPEED
    use time,    only : c
#endif /* GLOBAL_FR_SPEED */
#if defined SHEAR && defined TIDAL_FORCES
    use start,  only : qshear, omega
#endif /* SHEAR && TIDAL_FORCES */
#ifdef ANY_SOURCES
    use allsources, only : all_src
#endif /* ANY_SOURCES */
#ifdef ANY_LIMITS
    use floor_ceil, only : floor_ceiling
#endif /* ANY_LIMITS */

    implicit none
    integer               :: i1,i2,n
    real                  :: dt,dx,dtx
    real, dimension(nu,n) :: Duu,uu,cfr,ul,ur,Duus
    real, dimension(3,n)  :: bb

    character sweep*6

!locals
    real, dimension(nu,n) :: w,fr,fl,dulf,durf
    real, dimension(nu,n) :: Dul,Dur
#ifdef GRAV
    real, dimension(n)    :: gravl,gravr
    real, dimension(nu,n) :: duls,durs
    integer               :: ifluid
#endif /* GRAV */
#if defined SHEAR && defined TIDAL_FORCES
    real, dimension(nfluid,n)    :: vxr,rotfr,rotfl
#endif /* SHEAR && TIDAL_FORCES */
#ifdef COSM_RAYS
    real, dimension(n)    :: divv,decr,gpcr,ecr,vx
#endif /* COSM_RAYS */

    w         = 0.0
    cfr       = 0.0
    dtx       = dt / dx

#ifdef GRAV
    duls = 0.0
    durs = 0.0
#endif /* GRAV */

    Dul  = 0.0
    Dur  = 0.0

    call mhdflux(w,cfr,uu,bb,n)

    fr = (uu*cfr+w)*0.5
    fl = (uu*cfr-w)*0.5
    ur = fr/cfr
    ul = fl/cfr

    fl(:,1:n-1) = fl(:,2:n) ; fl(:,n) = fl(:,n-1)

#if defined SPLIT || defined ORIG
    if(istep == 1) then
      ur0 = ur
      ul0 = ul
    endif
#endif /* SPLIT || ORIG */

    call flux_limit(fr,fl,nu,n)

    durf(:,2:n) = dtx*(fr(:,2:n) - fr(:,1:n-1));     durf(:,1) = durf(:,2)
    dulf(:,2:n) = dtx*(fl(:,2:n) - fl(:,1:n-1));     dulf(:,1) = dulf(:,2)

!Do przemyslenia!!
!#ifdef SPLIT
!#ifdef SSP
!    Dur = ur-durf
!    Dul = ul+dulf
!#else /* SSP */
!    Dur = -durf
!    Dul =  durf
!#endif /* SSP */
!#else /* SPLIT */
    Dur = -durf
    Dul =  dulf
!#endif /* SPLIT */

#if defined SHEAR && defined TIDAL_FORCES
    vxr(:,1:n-1) = 0.5*(uu(imya,2:n)/uu(idna,2:n) + uu(imya,1:n-1)/uu(idna,1:n-1))

#ifdef SHEAR_MY
    if(sweep .eq. 'xsweep') then
      rotfr(:,:) =   2.0*omega*vxr(:,:)
    else if(sweep .eq. 'ysweep')  then
      rotfr(:,:) = (- 2.0 + qshear)*omega*vxr(:,:)
    else
      rotfr(:,:) = 0.0
    endif
#endif /* SHEAR_MY */

#ifdef SHEAR_MPI
    if(sweep .eq. 'xsweep') then
      rotfr(:,:) =   2.0*omega*vxr(:,:)  + qshear*omega**2*xr(:)
    else if(sweep .eq. 'ysweep')  then
      rotfr(:,:) = - 2.0*omega*vxr(:,:)
    else
      rotfr(:,:) = 0.0
    endif

#endif /* SHEAR_MPI */
    rotfr(:,n)   = rotfr(:,n-1)
    rotfl(:,2:n) = rotfr(:,1:n-1)
    rotfl(:,1)   = rotfl(:,2)

    call grav_limit(rotfr,rotfl,n)

#ifndef ISO
    duls(iena(fadiab),:)  = rotfr*UL_TAB(imxa(fadiab),:)*dt
    durs(iena(fadiab),:)  = rotfl*UR_TAB(imxa(fadiab),:)*dt
#endif /* ISO */
    duls(imxa,:)  = rotfr*UL_TAB(idna,:)*dt
    durs(imxa,:)  = rotfl*UR_TAB(idna,:)*dt
    Dur = Dur + durs
    Dul = Dul + duls
#endif /* SHEAR && TIDAL_FORCES */

#ifdef GRAV
! Gravity source terms -------------------------------------

    call grav_pot2accel(sweep,i1,i2, n, gravr)
    gravr(n)   = gravr(n-1)
    gravl(2:n) = gravr(1:n-1)
    gravl(1)   = gravl(2)

    call grav_limit(gravr,gravl,n)

    do ifluid=1,nfluid
#ifndef ISO
    if (ifluid .le. nadiab) then
    duls(iena(fadiab(ifluid)),:)  = gravr*UL_TAB(imxa(fadiab(ifluid)),:)*dt
    durs(iena(fadiab(ifluid)),:)  = gravl*UR_TAB(imxa(fadiab(ifluid)),:)*dt
    endif
#endif /* ISO */
    duls(imxa(ifluid),:)  = gravr*UL_TAB(idna(ifluid),:)*dt
    durs(imxa(ifluid),:)  = gravl*UR_TAB(idna(ifluid),:)*dt
    enddo
    Dur = Dur + durs
    Dul = Dul + duls
#endif /* GRAV */

!--- user-defined left- and rightsided source terms
!-
!-   Dul = Dul + LEFT_SOUCRE_TERM
!-   Dur = Dur + RIGHT_SOUCRE_TERM
!-
!-   NOTE: there may be also needed an addition in integrate procedure (e.g. COSM_RAYS)
!---

    Duu = Dul + Dur

#ifdef ANY_LIMITS
! wolt: mozliwe, ze wywolanie procedury floor_ceiling powinno byc tylko w integrate
    call floor_ceiling(1,n)
#endif /* ANY_LIMITS */
#ifdef ANY_SOURCES
    call all_src(Duus,uu,bb,sweep,i1,i2,n,dt)
    Duu = Duu + Duus
#endif /* ANY_SOURCES */

!--- user-defined centered source terms
!-
!-   Duu = Duu + SOUCRE_TERM
!-
!-   NOTE: there may be also needed an addition in integrate procedure (e.g. COSM_RAYS)
!---

    return
  end subroutine relaxing_tvd

!========================================================================
#ifndef SPLIT
  subroutine initials
    use arrays, only : ui,bi,u,b
    implicit none
    ui(:,:,:,:) = u(:,:,:,:)
    bi(:,:,:,:) = b(:,:,:,:)

  end subroutine initials
#endif /* SPLIT */
!========================================================================
#ifdef SPLIT
  subroutine integrate(uu,Duu,bb,sweep,i1,i2,n)
#else /* SPLIT */
  subroutine integrate(sweep,n)
#endif /* SPLIT */

    use start,  only : istep,smalld,cn
    use func, only : pshift, mshift
#ifndef SPLIT
#ifdef SSP
    use start,  only : integration_order
#endif /* SSP */
    use fluid_boundaries, only : compute_u_bnd
    use mag_boundaries, only :   compute_b_bnd
#endif /* !SPLIT */
    use arrays, only : idna,imxa,imya,imza,ibx,iby,ibz, &
#ifdef SPLIT
                       ul0,ur0,nu
#else /* SPLIT */
                       ui,u,bi,b,Lu,Lb,nx,ny,nz
#endif /* SPLIT */
#ifndef ISO
    use arrays, only : iena,nadiab
    use start,  only : smallei
#endif /* ISO */
#ifdef ANY_LIMITS
    use floor_ceil, only : floor_ceiling
#endif /* ANY_LIMITS */
#ifdef SHEAR_MY
    use start, only : qshear, omega, dt
#endif /* SHEAR_MY */

  implicit none
  integer               :: n,i1,i2
  character sweep*6
#ifdef SPLIT
  real, dimension(nu,n) :: uu,Duu
  real, dimension(3,n)  :: bb
#endif /* SPLIT */
#ifndef ISO
#ifdef SPLIT
    real, dimension(nadiab,n)    :: ekin,eint
    real, dimension(n)    :: emag
#else /* SPLIT */
  real, allocatable :: ekin(:,:,:,:), emag(:,:,:), eint(:,:,:,:)

    allocate(ekin(nadiab,nx,ny,nz),emag(nx,ny,nz),eint(nadiab,nx,ny,nz))
#endif /* SPLIT */
#endif /* ISO */

#ifdef SPLIT
#ifdef ORIG
    uu = cn(1,istep)*(ul0+ur0) + cn(2,istep)*Duu
#endif /* ORIG */
#ifdef SSP
    uu = cn(1,istep)*(ul0+ur0) + cn(2,istep)*(uu+Duu)
#endif /* SSP */
#else /* SPLIT */
        u(:,:,:,:) = cn(1,istep)*ui(:,:,:,:) + cn(2,istep)*(cn(3,istep)*u(:,:,:,:)+Lu(:,:,:,:))

#endif /* SPLIT */
    U_IDNA = max(U_IDNA,smalld)

#ifndef ISO
    ekin = 0.5*(U_IMXA*U_IMXA+U_IMYA*U_IMYA &
               +U_IMZA*U_IMZA)/U_IDNA
#ifdef SPLIT
    emag = 0.5*(B_IBX*B_IBX+B_IBY*B_IBY &
               +B_IBZ*B_IBZ)
#else /* SPLIT */
    emag = 0.125*( (B_IBX+pshift(B_IBX,dim=1))**2 &
                  + (B_IBY+pshift(B_IBY,dim=2))**2 &
                  + (B_IBZ+pshift(B_IBZ,dim=3))**2)

#endif /* SPLIT */

    eint = U_IENA-ekin-spread(emag,1,nadiab)

#ifdef SHEAR_MY
    if (sweep == 'xsweep') then
      eint = eint + cn(2,istep)* qshear*Omega*U_IMYA*(U_IMXA/U_IDNA)*dt
    endif
#endif /* SHEAR_MY */

#ifdef ANY_COOLING
    call all_cooling(eint,n)
#endif /* ANY_COOLING */
    eint = max(eint,smallei)
    U_IENA = eint+ekin+spread(emag,1,nadiab)
#ifndef SPLIT
    deallocate(ekin,emag,eint)
#endif /* SPLIT */
#endif /* ISO */
#ifdef ANY_LIMITS
! Dwukrotne ograniczanie vz moze byc zbedne (patrz wyzej)
! Ten fragment bedzie mozna ewentualnie wykasowac po testach
! wolt: pytanie, czy ten nie powinien zostac, a wykasowac ten pierwszy w relaxing_tvd
    call floor_ceiling(2,n)
#endif /* ANY_LIMITS */

#ifndef SPLIT
    call compute_u_bnd
    b(:,:,:,:) = cn(1,istep)*bi(:,:,:,:) + cn(2,istep)*(cn(3,istep)*b(:,:,:,:)+Lb(:,:,:,:))

    call compute_b_bnd
#endif /* SPLIT */

  end subroutine integrate

end module tv
