! $Id$
#include "piernik.def"

module arrays

   implicit none
   integer,parameter  :: nrlscal=100, nintscal=100

   real, allocatable, dimension(:,:,:,:) :: u, b
   real, allocatable, dimension(:,:,:) :: wa, wcu
#ifdef GRAV
   real, allocatable, dimension(:,:,:) :: gp
   real, allocatable, dimension(:)     :: dprof, eprof
#endif /* GRAV */

   real(kind=4), allocatable, dimension(:,:,:)  :: outwa, outwb, outwc

#if defined COSM_RAYS || defined PRESS_GRAD_EXCH
   real, allocatable, dimension(:,:,:,:) :: divvel
#endif /* COSM_RAYS || PRESS_GRAD_EXCH */

   real,    allocatable, dimension(:)       :: rlscal
   integer, allocatable, dimension(:)       :: intscal
   real,    allocatable, dimension(:,:)     :: ul0, ur0
   real,    allocatable, dimension(:,:,:,:) :: bndxrar, bndyrar

   contains
   
   subroutine arrays_allocate(nx,ny,nz,nvar)
      implicit none
      integer, intent(in) :: nx,ny,nz,nvar

      allocate(rlscal(nrlscal),intscal(nintscal)) ! arrays to store additional
                                                  ! scalar quantities in restart files
      allocate(u(nvar,nx,ny,nz),b(3,nx,ny,nz))

#ifdef GRAV
      allocate(gp(nx,ny,nz))
      allocate(dprof(nz),eprof(nz))
#endif /* GRAV */

      allocate(wa(nx,ny,nz),wcu(nx,ny,nz))
      allocate(outwa(nx,ny,nz),outwb(nx,ny,nz),outwc(nx,ny,nz))

!#if defined COSM_RAYS || defined PRESS_GRAD_EXCH
!      allocate(divvel(nfluid,nx,ny,nz))
!#endif /* COSM_RAYS || PRESS_GRAD_EXCH */

   end subroutine arrays_allocate

   subroutine arrays_deallocate

      deallocate(rlscal,intscal)
      deallocate(u,b)
      deallocate(wa,wcu)

#ifdef GRAV
      deallocate(gp)
      deallocate(dprof,eprof)
#endif /* GRAV */

#if defined COSM_RAYS || defined PRESS_GRAD_EXCH
      deallocate(divvel)
#endif /* COSM_RAYS || PRESS_GRAD_EXCH */

      if(allocated(bndxrar)) deallocate(bndxrar)
      if(allocated(bndyrar)) deallocate(bndyrar)
      deallocate(outwa,outwb,outwc)

   end subroutine arrays_deallocate

end module arrays

