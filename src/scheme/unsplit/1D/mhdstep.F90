! $Id$
#include "mhd.def"
module mod_mhdstep ! UNSPLIT 1D
#ifdef RESIST
   use resistivity
#endif /* RESIST */
#ifdef SHEAR
   use shear, only : yshift
   use fluid_boundaries, only: bnd_u
#endif /* SHEAR */
#ifdef SNE_DISTR
  use sn_distr
#endif /* SNE_DISTR */
  implicit none

   contains
      subroutine mhdstep
        use start, only : dimensions,dt,dt_log,nstep,t
        use dataio, only : nlog, write_log
        use time, only : timestep
        use mpi_setup, only : proc
#ifdef SELF_GRAV
        use poisson_solver, only : poisson
#endif /* SELF_GRAV */

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
900      format('   nstep = ',i7,'   dt = ',f22.16,'    t = ',f22.16)
         t=t+dt
!------------------- X->Y->Z ---------------------
#ifdef SHEAR
         call yshift(t)
#endif /* SHEAR */
#ifdef SELF_GRAV
         call poisson
#endif /* SELF_GRAV */

         call sweepx
         call sweepy
         if(dimensions .eq. '3d') then
            call sweepz
         endif
!         stop
#ifdef SNE_DISTR
      call supernovae_distribution
#endif /* SNE_DISTR */
 
        t = t+dt
!------------------- Z->Y->X ---------------------
#ifdef SHEAR
         call yshift(t)
#endif /* SHEAR */
#ifdef SELF_GRAV
         call poisson
#endif /* SELF_GRAV */


         if(dimensions .eq. '3d') then
            call sweepz
         endif
         call sweepy
         call sweepx
      end subroutine mhdstep

!-------------------------------------------------

      subroutine sweepx
        use start, only: magfield,istep,integration_order
        use tv, only : initials,integrate
        use fluids, only : fluidx

         call initials
         do istep=1,integration_order
            if(magfield) call magfieldbyzx
            call fluidx
            call integrate
         enddo

      end subroutine sweepx

!-----------------------------------------------------------------

      subroutine sweepy
        use start, only: magfield,istep,integration_order
        use tv, only : initials,integrate
        use fluids, only : fluidy

         call initials
         do istep=1,integration_order
            if(magfield) call magfieldbzxy
            call fluidy
            call integrate
         enddo
#ifdef FLX_BND
         call bnd_u('xdim')
         call bnd_u('ydim')
#endif /* FLX_BND */

      end subroutine sweepy

!------------------------------------------------------------------

      subroutine sweepz
        use start, only: magfield,istep,integration_order
        use tv, only : initials,integrate
        use fluids, only : fluidz

         call initials
         do istep=1,integration_order
            if(magfield) call magfieldbxyz
            call fluidz
            call integrate
         enddo

      end subroutine sweepz


  subroutine magfieldbyzx
    use start, only : dimensions
    use advects, only : advectbz_x, advectby_x

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
    use start, only : dimensions
    use advects, only : advectbz_y, advectbx_y

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
    use advects, only : advectbx_z, advectby_z

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
