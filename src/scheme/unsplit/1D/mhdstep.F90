#include "mhd.def"
module mod_mhdstep ! UNSPLIT 1D
   use start
   use dataio
   use time, only : timestep
   use fluids
   use advects
#ifdef RESIST
   use resistivity
#endif RESIST
#ifdef SHEAR
   use shear, only : yshift
#endif

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
900      format('   nstep = ',i7,'   dt = ',f22.16,'    t = ',f22.16)
         t=t+2.*dt
!------------------- X->Y->Z ---------------------
#ifdef SHEAR
         call yshift(t-dt)
#endif
         call sweepx
         call sweepy
         if(dimensions .eq. '3d') then
            call sweepz
         endif

!------------------- Z->Y->X ---------------------
#ifdef SHEAR
         call yshift(t)
#endif

         if(dimensions .eq. '3d') then
            call sweepz
         endif
         call sweepy
         call sweepx
      end subroutine mhdstep

!-------------------------------------------------

      subroutine sweepx

         call initials
         do istep=1,integration_order
            if(magfield) call magfieldbyzx
            call fluidx
            call integrate
         enddo

      end subroutine sweepx

!-----------------------------------------------------------------

      subroutine sweepy

         call initials
         do istep=1,integration_order
            if(magfield) call magfieldbzxy
            call fluidy
            call integrate
         enddo

      end subroutine sweepy

!------------------------------------------------------------------

      subroutine sweepz

         call initials
         do istep=1,integration_order
            if(magfield) call magfieldbxyz
            call fluidz
            call integrate
         enddo

      end subroutine sweepz


  subroutine magfieldbyzx

      call advectby_x
#ifdef RESIST
      call diffuseby_x
#endif RESIST

    if(dimensions .eq. '3d') then

      call advectbz_x
#ifdef RESIST
      call diffusebz_x
#endif RESIST

    endif

  end subroutine magfieldbyzx

!------------------------------------------------------------------------------------------

  subroutine magfieldbzxy

    if(dimensions .eq. '3d') then

      call advectbz_y
#ifdef RESIST
      call diffusebz_y
#endif RESIST

    endif

      call advectbx_y
#ifdef RESIST
      call diffusebx_y
#endif RESIST

  end subroutine magfieldbzxy

!------------------------------------------------------------------------------------------

  subroutine magfieldbxyz

      call advectbx_z
#ifdef RESIST
      call diffusebx_z
#endif RESIST

      call advectby_z
#ifdef RESIST
      call diffuseby_z
#endif RESIST

  end subroutine magfieldbxyz

end module mod_mhdstep
