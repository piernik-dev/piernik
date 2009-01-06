! $Id$
#include "piernik.def"
#ifdef IONIZED
#define NUMBION IONIZED
#else /* IONIZED */
#define NUMBION 0
#endif /* IONIZED */
#ifdef NEUTRAL
#define NUMBNEUT NEUTRAL
#else /* NEUTRAL */
#define NUMBNEUT 0
#endif /* NEUTRAL */
#ifdef DUST
#define NUMBDUST DUST
#else /* DUST */
#define NUMBDUST 0
#endif /* DUST */
#ifdef COSM_RAYS
#define NUMBCR COSM_RAYS
#else /* COSM_RAYS */
#define NUMBCR 0
#endif /* COSM_RAYS */
#define NUMBFLUID NUMBION+NUMBNEUT+NUMBDUST
#ifdef ISO
#define NUPF 4
#else /* ISO */
#define NUPF 5
#endif /* ISO */

module arrays

   implicit none
   integer :: nx, ny, nz
   integer,parameter  :: xdim=1, ydim=2, zdim=3
   integer,parameter  :: nrlscal=100, nintscal=100
   integer,allocatable,dimension(:) :: idna, imxa, imya, imza, iena, magn
   integer,allocatable,dimension(:) :: fmagn, fdust, findex, fadiab

   integer, parameter :: nfi=NUMBION
   integer, parameter :: nui=nfi*NUPF
   integer, parameter :: nfn=NUMBNEUT
   integer, parameter :: nun=nfn*NUPF
   integer, parameter :: nfd=NUMBDUST
   integer, parameter :: nud=nfd*4
   integer, parameter :: nfluid=NUMBFLUID
   integer, parameter :: nadiab=(NUMBFLUID-NUMBDUST)*(NUPF-4)
   integer, parameter :: nfion =NUMBION
   integer, parameter :: nfmagn=NUMBION
   integer, parameter :: ndust =NUMBDUST
   integer, parameter :: nuc=NUMBCR

   integer, parameter :: nu=nui+nun+nud+nuc
   integer, parameter :: nm=3

   integer, dimension(nu) :: iuswpx, iuswpy, iuswpz
   integer, dimension(nm) :: ibswpx, ibswpy, ibswpz

   integer, parameter  :: ibx=1, iby=2, ibz=3
   integer, parameter  :: icx=1, icy=2, icz=3

   integer nxb, nyb, nzb
   integer nxt, nyt, nzt
   integer is, ie, js, je, ks, ke

   real, allocatable :: dl(:)

   real, allocatable, dimension(:)  :: x, xl, xr
   real, allocatable, dimension(:)  :: y, yl, yr
   real, allocatable, dimension(:)  :: z, zl, zr

   real, allocatable, dimension(:,:,:,:) :: u, b
#ifdef GRAV
   real, allocatable, dimension(:,:,:) :: gp
   real, allocatable, dimension(:)     :: dprof, eprof
#endif /* GRAV */
#ifdef SPLIT
   real, allocatable, dimension(:,:,:) :: wa, wcu
#ifdef SSP
   real, allocatable, dimension(:,:,:,:) :: bi
#endif /* SSP */

#else /* ~SPLIT */
   real, allocatable, dimension(:,:,:,:) :: bi
   real, allocatable, dimension(:,:,:,:) :: ui
   real, allocatable, dimension(:,:,:,:) :: Lu, Lb
   real, allocatable, dimension(:,:,:)   :: wa
#ifdef FLX_BND
   real, allocatable, dimension(:,:,:,:) :: flx
   real, allocatable, dimension(:,:,:)   :: cfr
#endif /* FLX_BND */
#endif /* ~SPLIT */
#ifdef MASS_COMPENS
    real, allocatable, dimension(:,:,:,:)   :: dinit
#endif /* MASS_COMPENS */

!#ifdef COOL_HEAT
   real, allocatable, dimension(:)       :: coolheat_profile
!#endif /* COOL_HEAT */

   real(kind=4), allocatable, dimension(:,:,:)  :: outwa, outwb, outwc

#if defined COSM_RAYS || defined PRESS_GRAD_EXCH
   real, allocatable, dimension(:,:,:,:) :: divvel
#endif /* COSM_RAYS || PRESS_GRAD_EXCH */

   real,    allocatable, dimension(:)       :: rlscal
   integer, allocatable, dimension(:)       :: intscal
   real,    allocatable, dimension(:,:)     :: ul0, ur0
   real,    allocatable, dimension(:,:,:,:) :: bndxrar, bndyrar

#ifdef KEPLER_SUPPRESSION
   real, allocatable, dimension(:,:,:,:) :: omx0, omy0
   real, allocatable, dimension(:,:)     :: alfsup
#endif /* KEPLER_SUPPRESSION */

   contains

   subroutine arrays_allocate
      use start, only : nxd, nyd, nzd, nb, dimensions, gamma
      use mpi_setup
#ifdef IONIZED
      use ionizeds, only : add_ion_index, specify_ionized
#endif /* IONIZED */
#ifdef NEUTRAL
      use neutrals, only : add_neut_index, specify_neutral
#endif /* NEUTRAL */
#ifdef DUST
      use dusts, only : add_dust_index, specify_dust
#endif /* DUST */
#ifdef COSM_RAYS
      use cosmic_rays, only : add_cr_index, specify_cosmrays
#endif /* COSM_RAYS */

      integer ifluid,imagn,iadiab,idust,iter
#ifdef IONIZED
      integer nuipf
#endif /* IONIZED */
#ifdef NEUTRAL
      integer nunpf
#endif /* NEUTRAL */
#ifdef DUST
      integer nudpf
#endif /* DUST */
      character fluidtype*10

      ifluid = 1
      imagn  = 1
      idust  = 1
      iadiab = 1
#ifdef IONIZED
      call add_ion_index
#endif /* IONIZED */
#ifdef NEUTRAL
      call add_neut_index(nui)
#endif /* NEUTRAL */
#ifdef DUST
      call add_dust_index(nui+nun+num)
#endif /* DUST */
#ifdef COSM_RAYS
      call add_cr_index(nui+nun+num+nud)
#endif /* COSM_RAYS */


      allocate(idna(nfluid),imxa(nfluid),imya(nfluid),imza(nfluid))
      allocate(magn(nfluid),fmagn(nfmagn),fdust(ndust),findex(nfluid),fadiab(nadiab))
      magn=0
#ifndef ISO
      allocate(iena(nadiab))
#endif /* ISO */
#ifdef IONIZED
      call specify_ionized(nfluid,nadiab,ifluid,iadiab,imagn,idna,imxa,imya,imza,iena,findex,magn,fmagn,gamma,fadiab)
#endif /* IONIZED */
#ifdef NEUTRAL
      call specify_neutral(nfluid,nadiab,ifluid,iadiab,idna,imxa,imya,imza,iena,findex,gamma,fadiab)
#endif /* NEUTRAL */
#ifdef DUST
      call specify_dust(nfluid,ifluid,idna,imxa,imya,imza,findex,gamma,fdust,idust)
#endif /* DUST */
      iuswpx(idna)=idna ; iuswpx(imxa)=imxa ; iuswpx(imya)=imya ; iuswpx(imza)=imza
      iuswpy(idna)=idna ; iuswpy(imxa)=imya ; iuswpy(imya)=imxa ; iuswpy(imza)=imza
      iuswpz(idna)=idna ; iuswpz(imxa)=imza ; iuswpz(imya)=imya ; iuswpz(imza)=imxa
#ifndef ISO
      iuswpx(iena) = iena
      iuswpy(iena) = iena
      iuswpz(iena) = iena
#endif /* ISO */
#ifdef COSM_RAYS
      iuswpx(iecr) = iecr
      iuswpy(iecr) = iecr
      iuswpz(iecr) = iecr
      call specify_cosmrays
#endif /* COSM_RAYS */

      ibswpx(ibx:ibz)  =(/ibx,iby,ibz/)
      ibswpy(ibx:ibz)  =(/iby,ibx,ibz/)
      ibswpz(ibx:ibz)  =(/ibz,iby,ibx/)

      if((mod(nxd, pxsize) .ne. 0) .or. &
         (mod(nyd, pysize) .ne. 0) .or. &
         (mod(nzd, pzsize) .ne. 0)) then
         call mpistop
         if (proc .eq. 0) then
            write(*,*) 'One of: (mod(nxd,pxsize) .or. mod(nyd,pysize) .or. mod(nzd,pzsize)) .ne. 0'
         endif
         stop
      endif

      nxb = nxd/pxsize     !
      nyb = nyd/pysize     ! Block 'physical' grid sizes
      nzb = nzd/pzsize     !

      nx=nxb+2*nb          !
      ny=nyb+2*nb          ! Block total grid sizes
      nz=nzb+2*nb          !

      nxt=nxd+2*nb         !
      nyt=nyd+2*nb         ! Domain total grid sizes
      nzt=nzd+2*nb         !


      is = nb+1
      ie = nb+nxb
      js = nb+1
      je = nb+nyb

      if(dimensions .eq. '3d') then
         ks = nb+1
         ke = nb+nzb
      else if(dimensions .eq. '2dxy') then
         nz     = 1
         nzd    = 1
         nzb    = 1
         pzsize = 1
         ks     = 1
         ke     = 1
      endif

      allocate(dl(3))

      allocate(rlscal(nrlscal),intscal(nintscal)) ! arrays to store additional
                                                ! scalar quantities in restart files


      allocate(x(nx), xl(nx), xr(nx))
      allocate(y(ny), yl(ny), yr(ny))
      allocate(z(nz), zl(nz), zr(nz))

#ifdef SPLIT
#ifdef ORIG
      allocate(u(nu,nx,ny,nz),b(3,nx,ny,nz))
#else /* ~ORIG */
      allocate(u(nu,nx,ny,nz),b(3,nx,ny,nz),bi(3,nx,ny,nz))
#endif /* ~ORIG */
#else /* ~SPLIT */
      allocate(u(nu,nx,ny,nz),b(3,nx,ny,nz))
      allocate(ui(nu,nx,ny,nz),bi(3,nx,ny,nz))
      allocate(Lu(nu,nx,ny,nz),Lb(3,nx,ny,nz))
#ifdef FLX_BND
      allocate(flx(nu,nx,ny,nz))
      allocate(cfr(nx,ny,nz))
#endif /* FLX_BND */
#endif /* ~SPLIT */

#ifdef GRAV
      allocate(gp(nx,ny,nz))
      allocate(dprof(nz),eprof(nz))
#endif /* GRAV */
#ifdef SPLIT
      allocate(wa(nx,ny,nz),wcu(nx,ny,nz))
#else /* SPLIT */
      allocate(wa(nx,ny,nz))
#endif /* SPLIT */
      allocate(outwa(nx,ny,nz),outwb(nx,ny,nz),outwc(nx,ny,nz))
#ifdef MASS_COMPENS
      allocate(dinit(nfluid,nx,ny,nz))
#endif /* MASS_COMPENS */
!#ifdef COOL_HEAT
      allocate( coolheat_profile(nz))
!#endif /* COOL_HEAT */
#if defined COSM_RAYS || defined PRESS_GRAD_EXCH
      allocate(divvel(nfluid,nx,ny,nz))
#endif /* COSM_RAYS || PRESS_GRAD_EXCH */

   end subroutine arrays_allocate

   subroutine arrays_deallocate

      deallocate(dl)
      deallocate(rlscal,intscal)
      deallocate(x, xl, xr)
      deallocate(y, yl, yr)
      deallocate(z, zl, zr)

#ifdef SPLIT
#ifdef ORIG
      deallocate(u,b)
#else /* ~ORIG */
      deallocate(u,b,bi)
#endif /* ~ORIG */
#else /* ~SPLIT */
      deallocate(u,b)
      deallocate(ui,bi)
      deallocate(Lu,Lb)
#ifdef FLX_BND
      deallocate(cfr,flx)
#endif /* FLX_BND */
#endif /* ~SPLIT */
#ifdef GRAV
      deallocate(gp)
      deallocate(dprof,eprof)
#endif /* GRAV */
#ifdef SPLIT
      deallocate(wa,wcu)
#else /* SPLIT */
      deallocate(wa)
#endif /* SPLIT */
#ifdef MASS_COMPENS
      deallocate(dinit)
#endif /* MASS_COMPENS */
!#ifdef COOL_HEAT
      deallocate(coolheat_profile)
!#endif /* COOL_HEAT */
#if defined COSM_RAYS || defined PRESS_GRAD_EXCH
      deallocate(divvel)
#endif /* COSM_RAYS || PRESS_GRAD_EXCH */
#ifdef KEPLER_SUPPRESSION
      if(allocated(alfsup)) deallocate(alfsup)
      if(allocated(omx0)) deallocate(omx0)
      if(allocated(omy0)) deallocate(omy0)
#endif /* KEPLER_SUPPRESSION */
      if(allocated(bndxrar)) deallocate(bndxrar)
      if(allocated(bndyrar)) deallocate(bndyrar)
      deallocate(outwa,outwb,outwc)

   end subroutine arrays_deallocate

end module arrays

