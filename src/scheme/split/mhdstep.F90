! $Id$
#include "piernik.def"
module mod_mhdstep   ! SPLIT

  implicit none

  contains

subroutine mhdstep
  use start,  only : dimensions,dt,dt_log,dt_tsl,nstep,t,magfield
  use dataio, only : nlog,ntsl,write_log,write_timeslice
  use time,   only : timestep
  use fluids, only : fluidx,fluidy,fluidz
  use mpi_setup, only : proc
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
#ifdef MASS_COMPENS
  use init_problem, only : mass_loss_compensate
#endif /* MASS_COMPENS */
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

      t=t+dt

#ifdef SHEAR
      call yshift(t,dt)
      call bnd_u('xdim')
      call bnd_u('ydim')
#endif /* SHEAR */
#ifdef SELF_GRAV
      call poisson
#endif /* SELF_GRAV */

!------------------- X->Y->Z ---------------------
#ifndef ONLYZSWEEP
      call fluidx
      if(magfield) call magfieldbyzx
#ifdef COSM_RAYS
      call cr_diff_x
#endif /* COSM_RAYS */
#ifdef DEBUG
      syslog = system('echo -n sweep x')
      call write_hdf
      nhdf = nhdf + 1
#endif /* DEBUG */
#endif /* ONLYZSWEEP */

#ifndef ONLYZSWEEP
      call fluidy
      if(magfield) call magfieldbzxy
#ifdef COSM_RAYS
      call cr_diff_y
#endif /* COSM_RAYS */
#ifdef DEBUG
      syslog = system('echo -n sweep y')
      call write_hdf
      nhdf = nhdf + 1
#endif /* DEBUG */
#endif /* ONLYZSWEEP */

    if(dimensions .eq. '3d') then
      call fluidz
      if(magfield) call magfieldbxyz
#ifdef COSM_RAYS
      call cr_diff_z
#endif  /* COSM_RAYS */
#ifdef DEBUG
      syslog = system('echo -n sweep z')
      call write_hdf
      nhdf = nhdf + 1
#endif /* DEBUG */
    endif

! Sources ----------------------------------------

#ifdef SN_SRC
#ifndef SNE_DISTR
      call random_sn
#endif /* SNE_DISTR */
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
    use arrays,  only : b,ibx,iby,ibz,xdim,ydim,zdim
    use advects, only : advectby_x,advectbz_x
#ifdef SSP
    use start,   only : integration_order,istep
    use arrays,  only : bi
#endif /* SSP */
#ifdef RESIST
    use resistivity, only : diffuseby_x,diffusebz_x
#endif /* RESIST */

#ifdef SSP
! Now bi is the initial magnetic field
    bi(:,:,:,:) = b(:,:,:,:)

    do istep=1, integration_order
#endif /* SSP */

      call advectby_x
#ifdef RESIST
      call diffuseby_x
#endif /* RESIST */
      call mag_add(iby,xdim,ibx,ydim)
#ifdef SSP
    enddo
#endif /* SSP */

    if(dimensions .eq. '3d') then

#ifdef SSP
! Now bi is the initial magnetic field
    bi(:,:,:,:) = b(:,:,:,:)

    do istep=1, integration_order
#endif /* SSP */
      call advectbz_x

#ifdef RESIST
      call diffusebz_x
#endif /* RESIST */

      call mag_add(ibz,xdim,ibx,zdim)
#ifdef SSP
    enddo
#endif /* SSP */

    endif

  end subroutine magfieldbyzx

!------------------------------------------------------------------------------------------

  subroutine magfieldbzxy
    use start,   only : dimensions
    use arrays,  only : b,ibx,iby,ibz,xdim,ydim,zdim
    use advects, only : advectbx_y,advectbz_y
#ifdef SSP
    use start,   only : integration_order,istep
    use arrays,  only : bi
#endif /* SSP */
#ifdef RESIST
    use resistivity, only : diffusebx_y,diffusebz_y
#endif /* RESIST */

    if(dimensions .eq. '3d') then

#ifdef SSP
! Now bi is the initial magnetic field
    bi(:,:,:,:) = b(:,:,:,:)

    do istep=1, integration_order
#endif /* SSP */
      call advectbz_y
#ifdef RESIST
      call diffusebz_y
#endif /* RESIST */
      call mag_add(ibz,ydim,iby,zdim)
#ifdef SSP
    enddo
#endif /* SSP */

    endif

#ifdef SSP
! Now bi is the initial magnetic field
    bi(:,:,:,:) = b(:,:,:,:)

    do istep=1,integration_order
#endif /* SSP */
      call advectbx_y
#ifdef RESIST
      call diffusebx_y
#endif /* RESIST */
      call mag_add(ibx,ydim,iby,xdim)
#ifdef SSP
    enddo
#endif /* SSP */

  end subroutine magfieldbzxy

!------------------------------------------------------------------------------------------

  subroutine magfieldbxyz
    use start,   only : dimensions
    use arrays,  only : b,ibx,iby,ibz,xdim,ydim,zdim
    use advects, only : advectbx_z,advectby_z
#ifdef SSP
    use start,   only : integration_order,istep
    use arrays,  only : bi
#endif /* SSP */
#ifdef RESIST
    use resistivity, only : diffusebx_z,diffuseby_z
#endif /* RESIST */

#ifdef SSP
! Now bi is the initial magnetic field
    bi(:,:,:,:) = b(:,:,:,:)

    do istep=1,integration_order
#endif /* SSP */
      call advectbx_z
#ifdef RESIST
      call diffusebx_z
#endif /* RESIST */
      call mag_add(ibx,zdim,ibz,xdim)
#ifdef SSP
    enddo
#endif /* SSP */

#ifdef SSP
! Now bi is the initial magnetic field
    bi(:,:,:,:) = b(:,:,:,:)

    do istep=1,integration_order
#endif /* SSP */
      call advectby_z
#ifdef RESIST
      call diffuseby_z
#endif /* RESIST */
      call mag_add(iby,zdim,ibz,ydim)
#ifdef SSP
    enddo
#endif /* SSP */

  end subroutine magfieldbxyz

  subroutine mag_add(ib1,dim1,ib2,dim2)
    use func,   only : pshift, mshift
    use arrays, only : b,dl,wa,wcu
    use mag_boundaries, only : compute_b_bnd
#ifdef SSP
    use start,   only : cn,istep
    use arrays,  only : bi
#endif /* SSP */

    implicit none
    integer             :: ib1,ib2,dim1,dim2

#ifdef ORIG
#ifdef RESIST
! DIFFUSION FULL STEP

    b(ib1,:,:,:) = b(ib1,:,:,:) - wcu/dl(dim1)
!   wcu = cshift(wcu,shift= 1,dim=dim1)
    wcu = pshift(wcu,dim1)
    b(ib1,:,:,:) = b(ib1,:,:,:) + wcu/dl(dim1)
!   wcu = cshift(wcu,shift=-1,dim=dim1)
    wcu = mshift(wcu,dim1)
    b(ib2,:,:,:) = b(ib2,:,:,:) + wcu/dl(dim2)
!   wcu = cshift(wcu,shift= 1,dim=dim2)
    wcu = pshift(wcu,dim2)
    b(ib2,:,:,:) = b(ib2,:,:,:) - wcu/dl(dim2)
#endif /* RESIST */
! ADVECTION FULL STEP

    b(ib1,:,:,:) = b(ib1,:,:,:) - wa/dl(dim1)
!   wa = cshift(wa,shift=-1,dim=dim1)
    wa = mshift(wa,dim1)
    b(ib1,:,:,:) = b(ib1,:,:,:) + wa/dl(dim1)
    b(ib2,:,:,:) = b(ib2,:,:,:) - wa/dl(dim2)
!   wa = cshift(wa,shift=1,dim=dim2)
    wa = pshift(wa,dim2)
    b(ib2,:,:,:) = b(ib2,:,:,:) + wa/dl(dim2)
#endif /* ORIG */

#ifdef SSP

    if (istep .ne. 1) then
      b(ib1,:,:,:) = cn(1,istep)*bi(ib1,:,:,:)+cn(2,istep)*b(ib1,:,:,:)
      b(ib2,:,:,:) = cn(1,istep)*bi(ib2,:,:,:)+cn(2,istep)*b(ib2,:,:,:)
    endif

#ifdef RESIST
! DIFFUSION STEP

    b(ib1,:,:,:) = b(ib1,:,:,:) - cn(2,istep)*wcu/dl(dim1)
!   wcu = cshift(wcu,shift=1,dim=dim1)
    wcu = pshift(wcu,dim1)
    b(ib1,:,:,:) = b(ib1,:,:,:) + cn(2,istep)*wcu/dl(dim1)
!   wcu = cshift(wcu,shift=-1,dim=dim1)
    wcu = mshift(wcu,dim1)
    b(ib2,:,:,:) = b(ib2,:,:,:) + cn(2,istep)*wcu/dl(dim2)
!   wcu = cshift(wcu,shift=1,dim=dim2)
    wcu = pshift(wcu,dim2)
    b(ib2,:,:,:) = b(ib2,:,:,:) - cn(2,istep)*wcu/dl(dim2)
#endif /* RESIST */
! ADVECTION STEP

    b(ib1,:,:,:) = b(ib1,:,:,:) - cn(2,istep)*wa/dl(dim1)
!   wa = cshift(wa,shift=-1,dim=dim1)
    wa = mshift(wa,dim1)
    b(ib1,:,:,:) = b(ib1,:,:,:) + cn(2,istep)*wa/dl(dim1)
    b(ib2,:,:,:) = b(ib2,:,:,:) - cn(2,istep)*wa/dl(dim2)
!   wa = cshift(wa,shift=1,dim=dim2)
    wa = pshift(wa,dim2)
    b(ib2,:,:,:) = b(ib2,:,:,:) + cn(2,istep)*wa/dl(dim2)
#endif /* SSP */

    call compute_b_bnd

  end subroutine mag_add

end module mod_mhdstep
