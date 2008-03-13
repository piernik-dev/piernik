#include "mhd.def" 
module mod_mhdstep ! UNSPLIT 3D
   use start
   use dataio
   use time, only : timestep
   use fluids
   use advects
#ifdef RESIST
   use resistivity
#endif /* RESIST */
#ifdef SNE_DISTR
  use sn_distr
#endif /* SNE_DISTR */
   contains

subroutine mhdstep
    
  implicit none
  real tmp

  call timestep

  if (dt_log .gt. 0.0) then
     if (nlog .lt. (int(t / dt_log) + 1)) then
       call write_log
       nlog = nlog + 1
     endif
  endif

  if(proc.eq.0) write(*,900) nstep,dt,t
900      format('   nstep = ',i7,'   dt = ',f22.16,'   t = ',f22.16)
  t=t+2.*dt

         
!------------------- X->Y->Z ---------------------

  call initials


  do istep=1,integration_order
    if(magfield) call magfieldbyzx
    call fluidx                           ! x sweep
    if(magfield) call magfieldbzxy
    call fluidy                           ! y sweep
    if(dimensions .eq. '3d') then
      if(magfield) call magfieldbxyz
      call fluidz                         ! z sweep
    endif
    call integrate
  enddo

#ifdef SNE_DISTR
      call supernovae_distribution
#endif /* SNE_DISTR */

!------------------- Z->Y->X ---------------------

  call initials

  do istep=1,integration_order
    if(dimensions .eq. '3d') then
      if(magfield) call magfieldbxyz
      call fluidz                         ! z sweep
    endif
    if(magfield) call magfieldbzxy
    call fluidy                           ! y sweep
    if(magfield) call magfieldbyzx
    call fluidx                           ! x sweep
    call integrate
  enddo

end subroutine mhdstep


  subroutine magfieldbyzx

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
