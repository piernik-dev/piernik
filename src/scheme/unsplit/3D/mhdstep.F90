! $Id$
#include "piernik.def"
module mod_mhdstep ! UNSPLIT 3D

  implicit none

  contains

subroutine mhdstep
  use start,  only : dimensions,dt,dt_log,dt_tsl,nstep,t,magfield
  use dataio, only : nlog,ntsl,write_log,write_timeslice
  use time,   only : timestep
  use mpi_setup, only : proc
#ifdef RESIST
  use resistivity
#endif /* RESIST */
#ifdef DEBUG
  use dataio, only : nhdf,write_hdf
#endif /* DEBUG */
#ifdef COSM_RAYS
  use cr_diffusion
#endif /* COSM_RAYS */
#ifdef SHEAR
  use shear, only : yshift
  use fluid_boundaries, only : bnd_u
#endif /* SHEAR */
#ifdef SN_SRC
  use sn_sources
#endif /* SN_SRC */
#ifdef SNE_DISTR
  use sn_distr
#endif /* SNE_DISTR */
#ifdef SELF_GRAV
  use poisson_solver, only : poisson
#endif /* SELF_GRAV */

  implicit none
#ifdef DEBUG
  integer system, syslog
#endif /* DEBUG */

  call timestep

  if(dt_log .gt. 0.0) then
    if(nlog .lt. (int(t / dt_log) + 1)) then
      call write_log
      nlog = nlog + 1
    endif
  endif

  if(dt_tsl .gt. 0.0) then
    if(ntsl .lt. (int(t / dt_tsl) + 1)) then
      call write_timeslice
      ntsl = ntsl + 1
    endif
  endif

  if(proc.eq.0) write(*,900) nstep,dt,t
900      format('   nstep = ',i7,'   dt = ',e22.16,'   t = ',e22.16)

      t=t+2.*dt

#ifdef SHEAR
      call yshift(t,2.*dt)
      call bnd_u('xdim')
      call bnd_u('ydim')
#endif /* SHEAR */
#ifdef SELF_GRAV
      call poisson
#endif /* SELF_GRAV */

!------------------- X->Y->Z ---------------------

  call initials


  do istep=1,integration_order
#ifndef ONLYZSWEEP
#ifdef COSM_RAYS
      call cr_diff_x
#endif /* COSM_RAYS */
      if(magfield) call magfieldbyzx
      call fluidx
#ifdef DEBUG
      syslog = system('echo -n sweep x')
      call write_hdf
      nhdf = nhdf + 1
#endif /* DEBUG */
#endif /* ONLYZSWEEP */

#ifndef ONLYZSWEEP
#ifdef COSM_RAYS
      call cr_diff_y
#endif /* COSM_RAYS */
      if(magfield) call magfieldbzxy
      call fluidy
#ifdef DEBUG
      syslog = system('echo -n sweep y')
      call write_hdf
      nhdf = nhdf + 1
#endif /* DEBUG */
#endif /* ONLYZSWEEP */

    if(dimensions .eq. '3d') then
#ifdef COSM_RAYS
      call cr_diff_z
#endif  /* COSM_RAYS */
      if(magfield) call magfieldbxyz
      call fluidz
#ifdef DEBUG
      syslog = system('echo -n sweep z')
      call write_hdf
      nhdf = nhdf + 1
#endif /* DEBUG */
    endif
    call integrate
  enddo

! Sources ----------------------------------------

#ifdef SN_SRC
      call random_sn
      call dipol_sn
#endif /* SN_SRC */

      t=t+dt
#ifdef SHEAR
      call yshift(t,dt)
      call bnd_u('xdim')
      call bnd_u('ydim')
#endif /* SHEAR */

#ifdef SELF_GRAV
      call poisson
#endif /* SELF_GRAV */
!-------------------------------------------------



!------------------- Z->Y->X ---------------------

  call initials

  do istep=1,integration_order
    if(dimensions .eq. '3d') then
#ifdef COSM_RAYS
      call cr_diff_z
#endif /* COSM_RAYS */
      if(magfield) call magfieldbxyz
      call fluidz
#ifdef DEBUG
      syslog = system('echo -n sweep z')
      call write_hdf
      nhdf = nhdf + 1
#endif /* DEBUG */
    endif

#ifndef ONLYZSWEEP
#ifdef COSM_RAYS
      call cr_diff_y
#endif /* COSM_RAYS */
      if(magfield) call magfieldbzxy
      call fluidy
#ifdef DEBUG
      syslog = system('echo -n sweep y')
      call write_hdf
      nhdf = nhdf + 1
#endif /* DEBUG */
#endif /* ONLYZSWEEP */

#ifndef ONLYZSWEEP
#ifdef COSM_RAYS
      call cr_diff_x
#endif /* COSM_RAYS */
      if(magfield) call magfieldbyzx
      call fluidx
#ifdef DEBUG
      syslog = system('echo -n sweep x')
      call write_hdf
      nhdf = nhdf + 1
#endif /* DEBUG */
#endif /* ONLYZSWEEP */
    call integrate
  enddo

#ifdef SNE_DISTR
      call supernovae_distribution
#ifdef DEBUG
      syslog = system('echo -n sne_dis')
      call write_hdf
      nhdf = nhdf + 1
#endif /* DEBUG */
#endif /* SNE_DISTR */

end subroutine mhdstep


  subroutine magfieldbyzx
    use start,   only : dimensions
    use advects, only : advectby_x,advectbz_x
#ifdef RESIST
    use resistivity, only : diffuseby_x,diffusebz_x
#endif /* RESIST */

      call advectby_x
#ifdef RESIST
      call diffuseby_x
#endif /* RESIST */

    if(dimensions .eq. '3d') then

      call advectbz_x

#ifdef RESIST
      call diffusebz_x
#endif /* RESIST */

    endif

  end subroutine magfieldbyzx

!------------------------------------------------------------------------------------------

  subroutine magfieldbzxy
    use start,   only : dimensions
    use advects, only : advectbx_y,advectbz_y
#ifdef RESIST
    use resistivity, only : diffusebz_y,diffusebx_y
#endif /* RESIST */

    if(dimensions .eq. '3d') then

      call advectbz_y
#ifdef RESIST
      call diffusebz_y
#endif /* RESIST */

    endif

      call advectbx_y
#ifdef RESIST
      call diffusebx_y
#endif /* RESIST */

  end subroutine magfieldbzxy

!------------------------------------------------------------------------------------------

  subroutine magfieldbxyz
    use advects, only : advectbx_z,advectby_z
#ifdef RESIST
    use resistivity, only : diffusebx_z,diffuseby_z
#endif /* RESIST */

      call advectbx_z
#ifdef RESIST
      call diffusebx_z
#endif /* RESIST */

      call advectby_z
#ifdef RESIST
      call diffuseby_z
#endif /* RESIST */

  end subroutine magfieldbxyz

end module mod_mhdstep
