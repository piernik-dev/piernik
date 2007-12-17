#include "mhd.def"
module mod_mhdstep   ! SPLIT
  use start, only  : dt_log, magfield, dimensions, t, dt, &
      dt_tsl, nstep, proc, istep, integration_order, cn 
  use dataio, only : write_log, write_timeslice, nlog, ntsl
  use time, only   : timestep
  use arrays, only : dl,ibx,iby,ibz,xdim,ydim,zdim,b,wa
  use fluids, only : fluidx, fluidy, fluidz
  use advects,only : advectby_x, advectbz_x, advectbx_y, &
      advectbz_y, advectbx_z, advectby_z
  use mag_boundaries, only : compute_b_bnd
#ifdef SSP
  use arrays, only : bi
#endif /* SSP */
#ifdef RESIST  
  use resistivity
#endif /* RESIST */
#ifdef SHEAR  
  use shear
#endif /* SHEAR */
#ifdef SN_SRC
  use sn_sources
#endif /* SN_SRC */
#ifdef SNE_DISTR
  use sn_distr
#endif /* SNE_DISTR */
#ifdef SELF_GRAV
  use poisson_solver
#endif /* SELF_GRAV */
#ifdef GALAXY
  use init_problem, only : mass_loss_compensate
#endif /* GALAXY */

  implicit none

  contains
subroutine mhdstep
    
  implicit none
  real tmp

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
      call yshift(t)
#endif /* SHEAR */
#ifdef SELF_GRAV
      call poisson
#endif /* SELF_GRAV */

        
!------------------- X->Y->Z ---------------------
      call fluidx                         ! x sweep                      
      if(magfield) call magfieldbyzx 
#ifdef COSM_RAYS          
      call cr_diff_x 
#endif /* COSM_RAYS */
      call fluidy                         ! y sweep                      
      if(magfield) call magfieldbzxy        
#ifdef COSM_RAYS          
      call cr_diff_y   
#endif /* COSM_RAYS */
    if(dimensions .eq. '3d') then
      call fluidz                         ! z sweep                      
      if(magfield) call magfieldbxyz        
#ifdef COSM_RAYS          
      call cr_diff_z   
#endif  /* COSM_RAYS */
    endif
! Sources ----------------------------------------

#ifdef SN_SRC          
      call random_sn 
      call dipol_sn  
#endif /* SN_SRC */
#ifdef SNE_DISTR
      call supernovae_distribution(dt)
#endif /* SNE_DISTR */
      
      t=t+dt
#ifdef SHEAR
      call yshift(t)
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
      if(magfield) call magfieldbxyz      ! z sweep                       
      call fluidz                         
    endif
#ifdef COSM_RAYS          
      call cr_diff_y   
#endif /* COSM_RAYS */
      if(magfield) call magfieldbzxy      ! y sweep                       
      call fluidy                         
#ifdef COSM_RAYS          
      call cr_diff_x   
#endif /* COSM_RAYS */               
      if(magfield) call magfieldbyzx      ! x sweep                      
      call fluidx                         

end subroutine mhdstep


  subroutine magfieldbyzx

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

    implicit none
    integer             :: ib1,ib2,dim1,dim2

#ifdef ORIG
#ifdef RESIST
! DIFFUSION FULL STEP

    b(ib1,:,:,:) = b(ib1,:,:,:) - wcu/dl(dim1)
    wcu = cshift(wcu,shift=1,dim=dim1)
    b(ib1,:,:,:) = b(ib1,:,:,:) + wcu/dl(dim1)
    wcu = cshift(wcu,shift=-1,dim=dim1)
    b(ib2,:,:,:) = b(ib2,:,:,:) + wcu/dl(dim2)
    wcu = cshift(wcu,shift=1,dim=dim2)
    b(ib2,:,:,:) = b(ib2,:,:,:) - wcu/dl(dim2)
#endif /* RESIST */
! ADVECTION FULL STEP

    b(ib1,:,:,:) = b(ib1,:,:,:) - wa/dl(dim1)
    wa = cshift(wa,shift=-1,dim=dim1)
    b(ib1,:,:,:) = b(ib1,:,:,:) + wa/dl(dim1)
    b(ib2,:,:,:) = b(ib2,:,:,:) - wa/dl(dim2)
    wa = cshift(wa,shift=1,dim=dim2)
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
    wcu = cshift(wcu,shift=1,dim=dim1)
    b(ib1,:,:,:) = b(ib1,:,:,:) + cn(2,istep)*wcu/dl(dim1)
    wcu = cshift(wcu,shift=-1,dim=dim1)
    b(ib2,:,:,:) = b(ib2,:,:,:) + cn(2,istep)*wcu/dl(dim2)
    wcu = cshift(wcu,shift=1,dim=dim2)
    b(ib2,:,:,:) = b(ib2,:,:,:) - cn(2,istep)*wcu/dl(dim2)
#endif /* RESIST */
! ADVECTION STEP

    b(ib1,:,:,:) = b(ib1,:,:,:) - cn(2,istep)*wa/dl(dim1)
    wa = cshift(wa,shift=-1,dim=dim1)
    b(ib1,:,:,:) = b(ib1,:,:,:) + cn(2,istep)*wa/dl(dim1)
    b(ib2,:,:,:) = b(ib2,:,:,:) - cn(2,istep)*wa/dl(dim2)
    wa = cshift(wa,shift=1,dim=dim2)
    b(ib2,:,:,:) = b(ib2,:,:,:) + cn(2,istep)*wa/dl(dim2)
#endif /* SSP */ 

    call compute_b_bnd

  end subroutine mag_add

end module mod_mhdstep
