! $Id$
#include "piernik.def"

module timestep

! Based on: Pen Arras & Wong (2003)
! Modified extensively by M. Hanasz November 2005 - May 2006
! Optimalization of cpu-time efficiency of mhdflux, tvd1 by K. Kowalik
! Modification history:  see "changelog"  -- warning: out of date

  use mpi_setup

#ifdef COSM_RAYS
  use cr_diffusion
#endif /* COSM_RAYS */

  implicit none
  real c_all

contains


  subroutine time_step
    use start, only : dt, tend, t
    use constants, only : small,big
    
#ifdef IONIZED
    use timestepionized, only : timestep_ion
    use timestepionized, only : dt_ion,c_ion
#endif /* IONIZED */

#ifdef NEUTRAL
    use timestepneutral, only : timestep_neut
    use timestepneutral, only : dt_neut,c_neut
#endif /* NEUTRAL */

#ifdef DUST
    use timestepdust, only : timestep_dust
    use timestepdust, only : dt_dust,c_dust
#endif /* DUST */

#ifdef COSM_RAYS
    use timestepcr, only : timestep_cr
    use timestepcr, only : dt_cr
#endif COSM_RAYS
    
    
#ifdef SIMPLE_COOL
    use start, only : tauc
#endif /* SIMPLE_COOL */
#ifdef RESIST
    use resistivity, only : dt_resist, timestep_resist
#endif /* RESIST */
#ifdef COLLISIONS
    use start, only : collfaq, dt_colls
#endif /* COLLISIONS */
#if defined COLLISIONS || defined KEPLER_SUPPRESSION
    use start, only : cfl_colls
    use fluidindex, only : idna
    use arrays, only : u
#endif /* COLLISIONS || KEPLER_SUPPRESSION */
#ifdef KEPLER_SUPPRESSION
    use fluidindex, only : idna,imxa,imya,imza
    use arrays, only : alfsup,nfluid,nz,omx0,omy0,x,y,nx,ny
    use start,  only : dt_supp
#endif /* KEPLER_SUPPRESSION */
    implicit none
#ifdef KEPLER_SUPPRESSION
    real,allocatable,dimension(:,:) :: velx,vely,dvx,dvy
    integer j,k,ifl
#endif /* KEPLER_SUPPRESSION */
! Timestep computation


    c_all = 0.0

#ifdef IONIZED
    call timestep_ion
    dt=min(dt_ion,(tend-t)/2.)
    c_all = max(c_all,c_ion)
#endif /* IONIZED */

#ifdef NEUTRAL
    call timestep_neut
    dt=min(dt_neut,(tend-t)/2.)
    c_all = max(c_all,c_neut)
#endif /* NEUTRAL */

#ifdef DUST
    call timestep_dust
    dt=min(dt_dust,(tend-t)/2.)
    c_all = max(c_all,c_dust)
#endif /* DUST */

#ifdef COSM_RAYS
    call timestep_cr
    dt=min(dt_cr,(tend-t)/2.)
#endif /* COSM_RAYS */

#ifdef RESIST
    call timestep_resist
    dt = min(dt,dt_resist)
#endif /* RESIST */

#ifdef COOL_HEAT
    call timestep_coolheat
    dt = min(dt,dt_coolheat)
#endif /* COOL_HEAT */

!#ifdef HEAT_COND
!    dt_heatcond = cfl_heatcond * 0.5*dxmn**2/(K_heatcond+small)
!    dt = min(dt,dt_heatcond)
!#endif /* HEAT_COND */

!#ifdef VISC
!    dt_visc = cfl_visc * 0.5*dxmn**2/(nu_bulk+small)
!    dt = min(dt,dt_visc)
!#endif /* VISC */

!#ifdef COSM_RAYS  !!! przeniesc do "timestep_cr"
!    dt_cr = cfl_cr * 0.5*dxmn**2/(K_cr_paral+K_cr_perp+small)
!    dt = min(dt,dt_cr)
!endif /* COSM_RAYS */

#ifdef SIMPLE_COOL
    dt = min(dt,0.01 * tauc)
#endif /* SIMPLE_COOL */
#ifdef COLLISIONS
    dt_colls = cfl_colls / (collfaq+small) / maxval(u(idna,:,:,:))
    dt = min(dt,dt_colls)
#endif /* COLLISIONS */
#ifdef KEPLER_SUPPRESSION

    dt_supp = big

    do ifl=1,nfluid
      do k=1,nz
        dt_supp = min(dt_supp, minval(abs(1./(alfsup+small)/u(idna(ifl),:,:,k))))
      enddo
    enddo
    dt = min(dt,dt_supp)
#endif /* KEPLER_SUPPRESSION */
  end subroutine time_step

!------------------------------------------------------------------------------------------

end module timestep
