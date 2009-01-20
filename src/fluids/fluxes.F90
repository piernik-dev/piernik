! $Id$
#include "piernik.def"
!>
!! \brief Module that collects all flux components from each fluid
!!
!! This module should not be changed by user.
!<
module fluxes

#ifdef IONIZED
  use initionized, only : iarr_ion
  use fluxesionized, only : flux_ion 
#endif /* IONIZED */
#ifdef NEUTRAL
  use initneutral, only : iarr_neu
  use fluxesneutral, only : flux_neu
#endif /* NEUTRAL */
#ifdef DUST
  use initdust, only : iarr_dst
  use fluxesdust, only : flux_dst
#endif /* DUST */
#ifdef COSM_RAYS
  use initcrs, only : iarr_crs
  use fluxescr, only : flux_crs
#endif /* COSM_RAYS */

contains

!>
!! \brief Subroutine which changes flux and cfr from mhdflux regarding specified fluids.
!! \param flux flux
!! \param cfr freezing speed
!! \param uu currently used fluid table
!! \param bb magnetic field x,y,z-components table
!! \param n number of cells in the current sweep
!<
subroutine all_fluxes(flux,cfr,uu,bb,n)

    use fluidindex, only : nvar

    use fluidindex, only : nmag

#ifdef IONIZED
    use fluidindex, only : nvar_ion
#endif /* IONIZED */

#ifdef NEUTRAL
    use fluidindex, only : nvar_neu
#endif /* NEUTRAL */

#ifdef DUST
    use fluidindex, only : nvar_dst
#endif /* DUST */

#ifdef COSM_RAYS
    use cosmic_rays, only : nvar_crs
#endif /* COSM_RAYS */

    implicit none
    integer n
    real, dimension(nvar,n)::flux,uu,cfr
    real, dimension(nmag,n)::bb

#ifdef IONIZED
    real, dimension(nvar_ion,n) :: fluxion,cfrion,uuion
#endif /* IONIZED */

#ifdef NEUTRAL
    real, dimension(nvar_neu,n) :: fluxneu,cfrneu,uuneu
#endif /* NEUTRAL */

#ifdef DUST
    real, dimension(nvar_dst,n) :: fluxdst,cfrdst,uudst
#endif /* DUST */

#ifdef COSM_RAYS
    real, dimension(nvar_crs,n) :: fluxcrs,uucrs
#endif /* COSM_RAYS */

#ifdef IONIZED
   uuion(:,:)=uu(iarr_ion,:)
   
!   write(*,*) 'ion:',nvar_ion
!   write(*,*) uuion
   
   call flux_ion(fluxion,cfrion,uuion,bb,n)
   flux(iarr_ion,:)=fluxion
   cfr(iarr_ion,:) =cfrion
   uu(iarr_ion,:)  =uuion
#endif /* IONIZED */

!   write(*,*) 'allfluxes:'
!   write(*,*) fluxion(1,:)
!   write(*,*) fluxion(2,:)
!   write(*,*) fluxion(5,:)
!   stop

#ifdef NEUTRAL
   uuneu(:,:)=uu(iarr_neu,:)

   write(*,*) 'neu:',nvar_neu
   write(*,*) uuneu
   stop

   call flux_neu(fluxneu,cfrneu,uuneu,n)
   flux(iarr_neu,:)=fluxneu
   cfr(iarr_neu,:)=cfrneu
   uu(iarr_neu,:)=uuneu
#endif /* NEUTRAL */

#ifdef DUST
   uudst=uu(iarr_dst,:)
   call flux_dst(fluxdst,cfrdst,uudst,n)
   flux(iarr_dst,:)=fluxdst
   cfr(iarr_dst,:)=cfrdst
   uu(iarr_dst,:)=uudst
#endif /* DUST */

#ifdef COSM_RAYS
   uucrs=uu(iarr_crs,:)
   call flux_crs(fluxcrs,uucrs,n)
   flux(iarr_crs,:)=fluxcrs
   uu(iarr_crs,:)=uucrs
#endif /* COSM_RAYS */

end subroutine all_fluxes

!==========================================================================================


  subroutine flimiter(f,a,b,m,n)
    implicit none
    integer m,n
    real, dimension(m,n) :: f,a,b,c
#ifdef VANLEER
      c = a*b
      where (c .gt. 0.0)
        f = f+2.0*c/(a+b)
      endwhere
#endif /* VANLEER */
#ifdef MONCEN
        f = f+(sign(1.0,a)+sign(1.0,b))*min(2.*abs(a),2.*abs(b),0.5*abs(a+b))*0.5
#endif /* MONCEN */
#ifdef MINMOD
        f = f+(sign(1.0,a)+sign(1.0,b))*min(abs(a),abs(b))*0.5
#endif /* MINMOD */
#ifdef SUPERBEE
      where (abs(a) .gt. abs(b))
        f = f+(sign(1.0,a)+sign(1.0,b))*min(abs(a), abs(2.0*b))*0.5
      elsewhere
        f = f+(sign(1.0,a)+sign(1.0,b))*min(abs(2.0*a), abs(b))*0.5
      endwhere
#endif /* SUPERBEE */

    return
  end subroutine flimiter
  
  subroutine flux_limit(fr,fl,m,n)
    use start, only : istep

    implicit none
    integer             :: m,n
    real,dimension(m,n) :: fr,fl,dfrp,dfrm,dflp,dflm

#ifdef ORIG
    if(istep == 2) then
      dfrp(:,1:n-1) = 0.5*(fr(:,2:n) - fr(:,1:n-1)); dfrp(:,n) = dfrp(:,n-1)
      dfrm(:,2:n)   = dfrp(:,1:n-1);                 dfrm(:,1) = dfrm(:,2)
      call flimiter(fr,dfrm,dfrp,m,n)

      dflp(:,1:n-1) = 0.5*(fl(:,1:n-1) - fl(:,2:n)); dflp(:,n) = dflp(:,n-1)
      dflm(:,2:n)   = dflp(:,1:n-1);                 dflm(:,1) = dflm(:,2)
      call flimiter(fl,dflm,dflp,m,n)
    endif
#endif /* ORIG */
#ifdef SSP
      dfrp(:,1:n-1) = 0.5*(fr(:,2:n) - fr(:,1:n-1)); dfrp(:,n) = dfrp(:,n-1)
      dfrm(:,2:n)   = dfrp(:,1:n-1);                 dfrm(:,1) = dfrm(:,2)
      call flimiter(fr,dfrm,dfrp,m,n)

      dflp(:,1:n-1) = 0.5*(fl(:,1:n-1) - fl(:,2:n)); dflp(:,n) = dflp(:,n-1)
      dflm(:,2:n)   = dflp(:,1:n-1);                 dflm(:,1) = dflm(:,2)
      call flimiter(fl,dflm,dflp,m,n)
#endif /* SSP */

  end subroutine flux_limit

  subroutine grav_limit(gravr,gravl,n)
   use start, only : istep

   implicit none
   integer             :: n
   real,dimension(n)   :: gravr,gravl,dgrp,dgrm,dglp,dglm

#ifdef ORIG
    if(istep == 2) then
      dgrp(1:n-1) = 0.5*(gravr(1:n-1) - gravr(2:n));  dgrp(n) = dgrp(n-1)
      dgrm(2:n) = dgrp(1:n-1)                      ;  dgrm(1) = dgrm(2)
      call flimiter(gravr,dgrm,dgrp,1,n)

      dglp(1:n-1) = 0.5*(gravl(2:n) - gravl(1:n-1));  dglp(n) = dglp(n-1)
      dglm(2:n)   = dglp(1:n-1)                    ;  dglm(1) = dglm(2)
      call flimiter(gravl,dglm,dglp,1,n)
    endif
#endif /* ORIG */
#ifdef SSP
      dgrp(1:n-1) = 0.5*(gravr(1:n-1) - gravr(2:n));  dgrp(n) = dgrp(n-1)
      dgrm(2:n) = dgrp(1:n-1)                      ;  dgrm(1) = dgrm(2)
      call flimiter(gravr,dgrm,dgrp,1,n)

      dglp(1:n-1) = 0.5*(gravl(2:n) - gravl(1:n-1));  dglp(n) = dglp(n-1)
      dglm(2:n)   = dglp(1:n-1)                    ;  dglm(1) = dglm(2)
      call flimiter(gravl,dglm,dglp,1,n)
#endif /* SSP */

  end subroutine grav_limit
  


end module fluxes
