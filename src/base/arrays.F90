! $Id$
#include "piernik.def"

module arrays

   implicit none
   integer :: nx, ny, nz
   integer,parameter  :: xdim=1, ydim=2, zdim=3
   integer,parameter  :: nrlscal=100, nintscal=100

   integer nxb, nyb, nzb
   integer nxt, nyt, nzt
   integer is, ie, js, je, ks, ke

   real, allocatable :: dl(:)

   real, allocatable, dimension(:)  :: x, xl, xr
   real, allocatable, dimension(:)  :: y, yl, yr
   real, allocatable, dimension(:)  :: z, zl, zr

   real, allocatable, dimension(:,:,:,:) :: u, b
   real, allocatable, dimension(:,:,:) :: wa, wcu
#ifdef GRAV
   real, allocatable, dimension(:,:,:) :: gp
   real, allocatable, dimension(:)     :: dprof, eprof
#endif /* GRAV */

#ifdef COOL_HEAT
   real, allocatable, dimension(:)       :: coolheat_profile
#endif /* COOL_HEAT */

   real(kind=4), allocatable, dimension(:,:,:)  :: outwa, outwb, outwc

#if defined COSM_RAYS || defined PRESS_GRAD_EXCH
   real, allocatable, dimension(:,:,:,:) :: divvel
#endif /* COSM_RAYS || PRESS_GRAD_EXCH */

   real,    allocatable, dimension(:)       :: rlscal
   integer, allocatable, dimension(:)       :: intscal
   real,    allocatable, dimension(:,:)     :: ul0, ur0
   real,    allocatable, dimension(:,:,:,:) :: bndxrar, bndyrar

   contains
   
   subroutine arrays_allocate

      use mpi_setup
      use start, only : nxd, nyd, nzd, nb, dimensions, gamma
      use fluidindex, only : nvar
      

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
      allocate(u(nvar,nx,ny,nz),b(3,nx,ny,nz))
#ifdef GRAV
      allocate(gp(nx,ny,nz))
      allocate(dprof(nz),eprof(nz))
#endif /* GRAV */
      allocate(wa(nx,ny,nz),wcu(nx,ny,nz))
      allocate(outwa(nx,ny,nz),outwb(nx,ny,nz),outwc(nx,ny,nz))
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
      deallocate(u,b)
      deallocate(wa,wcu)
#ifdef GRAV
      deallocate(gp)
      deallocate(dprof,eprof)
#endif /* GRAV */
#ifdef COOL_HEAT
      deallocate(coolheat_profile)
#endif /* COOL_HEAT */
#if defined COSM_RAYS || defined PRESS_GRAD_EXCH
      deallocate(divvel)
#endif /* COSM_RAYS || PRESS_GRAD_EXCH */
      if(allocated(bndxrar)) deallocate(bndxrar)
      if(allocated(bndyrar)) deallocate(bndyrar)
      deallocate(outwa,outwb,outwc)

   end subroutine arrays_deallocate

end module arrays

