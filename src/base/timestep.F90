! $Id$
#include "piernik.def"

module timestep

! Based on: Pen Arras & Wong (2003)
! Modified extensively by M. Hanasz November 2005 - May 2006
! Optimalization of cpu-time efficiency of mhdflux, tvd1 by K. Kowalik
! Modification history:  see "changelog"  -- warning: out of date

  use mpi_setup

!#ifdef COSM_RAYS
!  use cr_diffusion
!#endif /* COSM_RAYS */

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
    use timestepneutral, only : timestep_neu
    use timestepneutral, only : dt_neu,c_neu
#endif /* NEUTRAL */

#ifdef DUST
    use timestepdust, only : timestep_dst
    use timestepdust, only : dt_dst,c_dst
#endif /* DUST */

#ifdef COSM_RAYS
    use timestepcosmicrays, only : timestep_crs
    use timestepcosmicrays, only : dt_crs
#endif /* COSM_RAYS */
    
    
#ifdef SIMPLE_COOL
    use start, only : tauc
#endif /* SIMPLE_COOL */
#ifdef RESISTIVE
    use resistivity, only : dt_resist, timestep_resist
#endif /* RESISTIVE */
    implicit none
! Timestep computation


    c_all = 0.0

#ifdef IONIZED
    call timestep_ion
    dt=min(dt_ion,(tend-t)/2.)
    c_all = max(c_all,c_ion)
#endif /* IONIZED */

#ifdef NEUTRAL
    call timestep_neu
    dt=min(dt_neu,(tend-t)/2.)
    c_all = max(c_all,c_neu)
#endif /* NEUTRAL */

#ifdef DUST
    call timestep_dst
    dt=min(dt_dst,(tend-t)/2.)
    c_all = max(c_all,c_dst)
#endif /* DUST */

#ifdef COSM_RAYS
    call timestep_crs
    dt=min(dt_crs,(tend-t)/2.)
#endif /* COSM_RAYS */

#ifdef RESISTIVE
    call timestep_resist
    dt = min(dt,dt_resist)
#endif /* RESISTIVE */

!#ifdef COSM_RAYS  !!! przeniesc do "timestep_crs"
!    dt_crs = cfl_crs * 0.5*dxmn**2/(K_crs_paral+K_crs_perp+small)
!    dt = min(dt,dt_crs)
!endif /* COSM_RAYS */

#ifdef SIMPLE_COOL
    dt = min(dt,0.01 * tauc)
#endif /* SIMPLE_COOL */
  end subroutine time_step

!------------------------------------------------------------------------------------------

end module timestep
