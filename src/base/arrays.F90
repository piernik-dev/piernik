#include "mhd.def"

module arrays
  use start 
  use mpi_setup

  implicit none
  integer :: nx, ny, nz
  integer,parameter  :: xdim=1,ydim=2,zdim=3
  
#ifdef ISO
  integer,parameter  :: nua=4
#else ISO
  integer,parameter  :: nua=5
#endif ISO  

#ifdef COSM_RAYS
  integer,parameter  :: nuc=1
#else COSM_RAYS
  integer,parameter  :: nuc=0  
#endif COSM_RAYS

#ifdef DUST
  integer,parameter  :: nud=1
#else DUST
  integer,parameter  :: nud=0  
#endif DUST

  integer,parameter  :: nu=nua+nuc+nud
  
  integer,parameter  :: idna=1,imxa=2,imya=3,imza=4
#ifdef ISO
  integer,parameter  :: iena=0 
#else ISO
  integer,parameter  :: iena=5   
#endif ISO 

#ifdef COSM_RAYS
  integer,parameter  :: iecr = nua+1
#else COSM_RAYS
  integer,parameter  :: iecr = 0  
#endif COSM_RAYS

#ifdef DUST
  integer,parameter  :: idnd=nua+nuc+1,ivxd=nua+nuc+2,ivyd=nua+nuc+3,ivzd=nua+nuc+4
#else DUST
  integer,parameter  :: idnd=0,ivxd=0,ivyd=0,ivzd=0  
#endif DUST
   
  integer, dimension(nu) :: iuswpx,iuswpy,iuswpz    
  integer, dimension(3) :: ibswpx,ibswpy,ibswpz    
   
  integer,parameter  :: ibx=1,iby=2,ibz=3
  integer,parameter  :: icx=1,icy=2,icz=3

  integer nxb, nyb, nzb
  integer nxt, nyt, nzt
  integer is, ie, js, je, ks, ke 

  real, allocatable :: dl(:)

  real, allocatable  :: x(:), xl(:), xr(:)
  real, allocatable  :: y(:), yl(:), yr(:)
  real, allocatable  :: z(:), zl(:), zr(:)

  real, allocatable :: u(:,:,:,:)
  real, allocatable :: b(:,:,:,:)
#ifdef GRAV
  real, allocatable :: gp(:,:,:)     
  real, allocatable :: dprof(:),eprof(:)
#endif GRAV
#ifdef SPLIT
  real, allocatable :: wa(:,:,:),wcu(:,:,:)
#ifdef SSP
  real, allocatable :: bi(:,:,:,:)
#endif SSP

#else SPLIT
  real, allocatable :: bi(:,:,:,:)
  real, allocatable :: ui(:,:,:,:)
  real, allocatable :: Lu(:,:,:,:),Lb(:,:,:,:)
  real, allocatable :: wa(:,:,:)
#endif SPLIT


  real(kind=4), allocatable  :: outwa(:,:,:),outwb(:,:,:),outwc(:,:,:)

  real, allocatable  :: coolheat_profile(:)

contains
  
  subroutine arrays_allocate

    iuswpx(idna:imza)=(/idna,imxa,imya,imza/) 
    iuswpy(idna:imza)=(/idna,imya,imxa,imza/) 
    iuswpz(idna:imza)=(/idna,imza,imya,imxa/) 
    ibswpx(ibx:ibz)  =(/ibx,iby,ibz/)    
    ibswpy(ibx:ibz)  =(/iby,ibx,ibz/)    
    ibswpz(ibx:ibz)  =(/ibz,iby,ibx/)    
    
#ifndef ISO
    iuswpx(iena) = iena
    iuswpy(iena) = iena
    iuswpz(iena) = iena
#endif ISO
#ifdef COSM_RAYS
    iuswpx(iecr) = iecr
    iuswpy(iecr) = iecr
    iuswpz(iecr) = iecr
#endif COSM_RAYS

  
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

    allocate(x(nx), xl(nx), xr(nx))
    allocate(y(ny), yl(ny), yr(ny))
    allocate(z(nz), zl(nz), zr(nz))    

#ifdef SPLIT
#ifdef ORIG
    allocate(u(nu,nx,ny,nz),b(3,nx,ny,nz))
#else ORIG
    allocate(u(nu,nx,ny,nz),b(3,nx,ny,nz),bi(3,nx,ny,nz))
#endif ORIG
#else SPLIT
    allocate(u(nu,nx,ny,nz),b(3,nx,ny,nz))
    allocate(ui(nu,nx,ny,nz),bi(3,nx,ny,nz))
    allocate(Lu(nu,nx,ny,nz),Lb(3,nx,ny,nz))
#endif SPLIT

#ifdef GRAV
    allocate(gp(nx,ny,nz))
    allocate(dprof(nz),eprof(nz))
#endif GRAV
#ifdef SPLIT
    allocate(wa(nx,ny,nz),wcu(nx,ny,nz))
#else SPLIT
    allocate(wa(nx,ny,nz))
#endif SPLIT
    allocate(outwa(nx,ny,nz),outwb(nx,ny,nz),outwc(nx,ny,nz))

    if(coolheat) allocate( coolheat_profile(nz))

  end subroutine arrays_allocate
      
  subroutine arrays_deallocate

    deallocate(dl)
    deallocate(x, xl, xr)
    deallocate(y, yl, yr)
    deallocate(z, zl, zr)    

#ifdef SPLIT
#ifdef ORIG
    deallocate(u,b)
#else ORIG
    deallocate(u,b,bi)
#endif ORIG
#else SPLIT
    deallocate(u,b)
    deallocate(ui,bi)
    deallocate(Lu,Lb)
#endif SPLIT
#ifdef GRAV
    deallocate(gp)
    deallocate(dprof,eprof)
#endif GRAV
#ifdef SPLIT
    deallocate(wa,wcu)
#else SPLIT
    deallocate(wa)
#endif SPLIT
    deallocate(outwa,outwb,outwc)
#ifdef COOL_HEAT
    deallocate(coolheat_profile)
#endif COOL_HEAT

  end subroutine arrays_deallocate
      
end module arrays

