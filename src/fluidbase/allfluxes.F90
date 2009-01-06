! $Id$
#include "piernik.def"
!>
!! \brief Module that collects all flux components from each fluid
!!
!! This module should not be changed by user.
!<
module allfluxes
#ifdef IONIZED
  use ionizeds, only : flux_ionized, aions
#endif /* IONIZED */
#ifdef NEUTRAL
  use neutrals, only : flux_neutral, aneut
#endif /* NEUTRAL */
#ifdef DUST
  use dusts, only : flux_dust, adust
#endif /* DUST */
#ifdef COSM_RAYS
  use cosmic_rays, only : flux_cr
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
subroutine all_flux(flux,cfr,uu,bb,n)
    use arrays, only : nu,nui,nun,nud
#ifdef COSM_RAYS
    use cosmic_rays, only : iecr
#endif /* COSM_RAYS */
    implicit none
    integer n
    real, dimension(nu,n)::flux,uu,cfr
    real, dimension(3,n):: bb
#ifdef IONIZED
    real, dimension(nui,n) :: iflux,icfr,iuu
#endif /* IONIZED */
#ifdef NEUTRAL
    real, dimension(nun,n) :: nflux,ncfr,nuu
#endif /* NEUTRAL */
#ifdef DUST
    real, dimension(nud,n) :: dflux,dcfr,duu
#endif /* DUST */

#ifdef IONIZED
   iflux=flux(aions,:)
   icfr=cfr(aions,:)
   iuu=uu(aions,:)
   call flux_ionized(iflux,icfr,iuu,nui,bb,n)
   flux(aions,:)=iflux
   cfr(aions,:) =icfr
   uu(aions,:)  =iuu
#endif /* IONIZED */
#ifdef NEUTRAL
   nflux=flux(aneut,:)
   ncfr=cfr(aneut,:)
   nuu=uu(aneut,:)
   call flux_neutral(nflux,ncfr,nuu,nun,n)
   flux(aneut,:)=nflux
   cfr(aneut,:)=ncfr
   uu(aneut,:)=nuu
#endif /* NEUTRAL */
#ifdef DUST
   dflux=flux(adust,:)
   dcfr=cfr(adust,:)
   duu=uu(adust,:)
   call flux_dust(dflux,dcfr,duu,nud,n)
   flux(adust,:)=dflux
   cfr(adust,:)=dcfr
   uu(adust,:)=duu
#endif /* DUST */
#ifdef COSM_RAYS
   call flux_cr(flux,cfr,uu,bb,nu,n)
#endif /* COSM_RAYS */

end subroutine all_flux

end module allfluxes
