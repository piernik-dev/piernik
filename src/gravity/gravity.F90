#include "mhd.def"

module gravity

  use mpi_setup
  use start
  use arrays
  use grid
  use constants

  character gp_status*9

contains

!--------------------------------------------------------------------------
  subroutine grav_pot_3d
    use grid, only : z 
    implicit none
    integer i, j
    real, allocatable :: gpot(:)
    allocate(gpot(nz))

    do j = 1,ny
      do i = 1,nx
        call grav_pot('zsweep', i,j, z(:), nz, gpot, gp_status)
        if(gp_status .eq. 'undefined') exit
        gp(i,j,:) = gpot
      enddo
    enddo
    
    if(gp_status .eq. 'undefined') call grav_accel2pot

    deallocate(gpot)
    
  end subroutine grav_pot_3d


!--------------------------------------------------------------------------
  subroutine grav_pot(sweep, i1,i2, xsw, n, gpot,status)
   
    implicit none
    character, intent(in) :: sweep*6
    integer, intent(in)   :: i1, i2                   ! 
    integer,intent(in)    :: n
    real, dimension(n)    :: xsw, x3, rc,fr,gr
    real, dimension(n),intent(out) :: gpot
    character, intent(out) :: status*9
    real                  :: x1, x2
    
#if GRAV == 'null'
        gpot = 0.0
      ! do nothing
#elif GRAV == 'ptmass'
        select case (sweep)
          case('xsweep') 
            x1  = y(i1)-ptm_y
            x2  = z(i2)-ptm_z
            x3  = xsw - ptm_x
          case('ysweep') 
            x1  = z(i1)-ptm_z
            x2  = x(i2)-ptm_x
            x3  = xsw - ptm_y
          case('zsweep') 
            x1  = x(i1)-ptm_x
            x2  = y(i2)-ptm_y
            x3  = xsw - ptm_z
        end select
        rc = dsqrt(x1**2+x2**2)
        fr = min( (rc/r_grav)**n_gravr , 100.0)
        fr = max(1./cosh(fr),smalld/100.)
        gpot = -newtong*ptmass/dsqrt(x1**2+x2**2+x3**2+r_smooth**2)
        gpot = gpot - csim2*dlog(fr) ! *d0
#elif GRAV == 'ptflat'
        select case (sweep)
          case('xsweep')
            x1  = y(i1)-ptm_y
            x2  = 0.0
            x3  = xsw - ptm_x
          case('ysweep')
            x1  = 0.0
            x2  = x(i2)-ptm_x
            x3  = xsw - ptm_y
          case('zsweep')
            x1  = x(i1)-ptm_x
            x2  = y(i2)-ptm_y
            x3  = 0.0
        end select
        rc = dsqrt(x1**2+x2**2)
        fr(:) = min((rc(:)/r_grav)**n_gravr,100.0)
        fr = 1./cosh(fr)+smalld
        gpot = -newtong*ptmass/dsqrt(x1**2+x2**2+x3**2+r_smooth**2)
        gpot = gpot - csim2*dlog(fr) ! *d0
#else
        status = 'undefined'
#warning 'GRAV declared, but gravity model undefined in grav_pot'
#endif

  end subroutine grav_pot

!--------------------------------------------------------------------------
   
  subroutine grav_accel(sweep, i1,i2, xsw, n, grav)
   
    implicit none
    character, intent(in) :: sweep*6
    integer, intent(in)   :: i1, i2                   ! 
    integer, intent(in)   :: n
    real, dimension(:)    :: xsw
!DW+
    real, dimension(n)    :: rgc_vect
!DW-
    real, dimension(n),intent(out) :: grav
    
    real                  :: x1, x2
    real, dimension(n)    :: r, r32
!DW+
    real rstar
    real aconst, bconst, cconst, flat, omega, g_shear_rate
    real rgc_scal
!DW-
    
    integer i

    select case (sweep)
      case('xsweep') 
        x1  = y(i1)
        x2  = z(i2)
      case('ysweep') 
        x1  = z(i1)
        x2  = x(i2)
      case('zsweep') 
        x1  = x(i1)
        x2  = y(i2)
    end select
    
!    select case (grav_model)


#if GRAV == 'uniform'
!      case ('uniform')
       if (sweep == 'zsweep') then
          grav = g_z
       else
          grav = 0.0
       endif
#elif GRAV == 'linear'
       if (sweep == 'zsweep') then
          grav = dg_dz*xsw
       else
          grav = 0.0
       endif
#elif GRAV ==  'ptmass'
       if (sweep == 'xsweep') then
         r  = sqrt((xsw-ptm_x)*(xsw-ptm_x) + (x1-ptm_y)*(x1-ptm_y)  &
               + (x2-ptm_z)*(x2-ptm_z) + r_smooth*r_smooth)
       elseif(sweep == 'ysweep') then
         r  = sqrt((xsw-ptm_y)*(xsw-ptm_y) + (x1-ptm_z)*(x1-ptm_z)  &
               + (x2-ptm_x)*(x2-ptm_x) + r_smooth*r_smooth)
       else
         r  = sqrt((xsw-ptm_z)*(xsw-ptm_z) + (x1-ptm_x)*(x1-ptm_x)  &
               + (x2-ptm_y)*(x2-ptm_y) + r_smooth*r_smooth)
       endif
       r32  = r*r*r
       grav = -newtong*ptmass*xsw/r32
#elif GRAV ==  'ptflat'	
       if (sweep == 'xsweep') then
         r  = sqrt((xsw-ptm_x)*(xsw-ptm_x) + (x1-ptm_y)*(x1-ptm_y)  &
               + r_smooth*r_smooth)
         r32  = r*r*r
         grav = -newtong*ptmass*xsw/r32
       elseif(sweep == 'ysweep') then
         r  = sqrt((xsw-ptm_y)*(xsw-ptm_y) + (x2-ptm_x)*(x2-ptm_x)  &
               + r_smooth*r_smooth)
         r32  = r*r*r
         grav = -newtong*ptmass*xsw/r32
       else
         grav = 0.0
       endif
#elif GRAV ==  'gal_mali'
       if (sweep == 'zsweep') then
         grav = -0.87 * (6.8*tanh(3.2 * xsw/(pc*kpc)) + 1.7*xsw/(pc*kpc)) * 1.e-9
       else
         grav = 0.0
       endif
       
#elif GRAV ==  'galactic'
       if (sweep == 'zsweep') then
         grav = cmps2 * (  &
           (-4.4e-9 * exp(-(r_gc-r_gc_sun)/(4.9*kpc)) * xsw/sqrt(xsw**2+(0.2*kpc)**2)) &    
           -( 1.7e-9 * (r_gc_sun**2 + (2.2*kpc)**2)/(r_gc**2 + (2.2*kpc)**2)*xsw/kpc) )     
!          -Om*(Om+G) * Z * (kpc ?) ! in the transition region between rigid 
!                                   ! and flat rotation F'98: eq.(36)       
       else
         grav=0.0
       endif

!dw+
! galactic case as in ferriere'98
#elif GRAV ==  'galactic_dw'
        aconst = -vsun / (5000.0 * pc *(5000.0-3000.0) * pc)            ![pc/s/pc2]
        bconst = -3000.0 * pc * aconst                                  ![1/s]
        cconst =  vsun/5000.0*pc-aconst*5000.0*pc-bconst*log(5000.0*pc) ![1/s]
        flat   = aconst * 3000.0 * pc + bconst * log(3000.0*pc)+ cconst ![1/s]

        if (sweep == 'zsweep') then
          rgc_scal = sqrt(x(i1)**2+y(i2)**2)
          if(rgc_scal.le.3.0*kpc) then
            omega=flat
            g_shear_rate=0.0
          elseif(rgc_scal.ge.5.0*kpc) then
            omega=vsun/rgc_scal
            g_shear_rate=(-1.0)*vsun/rgc_scal
          else     ! 3.0 < x < 5.0
            omega=aconst*rgc_scal+bconst*log(rgc_scal)+cconst
            g_shear_rate=aconst*rgc_scal+bconst
          endif

          grav=-4.4e-9*cm/sek**2*exp((-1.0)*(rgc_scal-r_gc_sun)/(4.9*kpc))*xsw/sqrt(xsw**2+(0.2*kpc)**2)
          grav=grav-1.7e-9*cm/sek**2*(r_gc_sun**2+(2.2*kpc)**2)/(rgc_scal**2+(2.2*kpc)**2)*xsw/kpc
          grav=grav+2*omega*(omega+g_shear_rate)*xsw
          
        else  ! y or x
          if(sweep == 'xsweep') then 
            rgc_vect(:) = sqrt(xsw(:)**2 +   y(i1)**2)
          elseif(sweep == 'ysweep') then 
            rgc_vect(:) = sqrt(x(i2)**2  + xsw(i1)**2)
          endif

          grav = -1.0 * xsw
          where (rgc_vect <= 3.0*kpc)
             grav = grav * flat**2 
          elsewhere (rgc_vect >= 5.0*kpc)
             grav = grav * (vsun / rgc_vect)**2 
          elsewhere
             grav = grav * (aconst*rgc_vect + bconst*log(rgc_vect) + cconst)**2
          endwhere
        endif

!        select case (sweep)
!          case('xsweep')
!            do i=1,n 
!              rgc_vect(i) = sqrt(xsw(i)**2+y(i1)**2)
!              if(rgc_vect(i).le.3.0*kpc) then
!                omega=flat
!              else
!                if(rgc_vect(i).ge.5.0*kpc) then
!                  omega=vsun/rgc_vect(i)
!                else
!                  omega=aconst*rgc_vect(i)+bconst*log(rgc_vect(i))+cconst
!                endif
!              endif
!              grav(i)=(-1.0)*omega**2*rgc_vect(i)*xsw(i)/rgc_vect(i)
!            enddo
!          case('ysweep') 
!            do i=1,n
!              rgc_vect(i) = sqrt(x(i2)**2+xsw(i)**2)
!              if(rgc_vect(i).le.3.0*kpc) then
!                omega=flat
!              else
!                if(rgc_vect(i).ge.5.0*kpc) then
!                  omega=vsun/rgc_vect(i)
!                else
!                  omega=aconst*rgc_vect(i)+bconst*log(rgc_vect(i))+cconst
!                endif
!              endif
!              grav(i)=(-1.0)*omega**2*rgc_vect(i)*xsw(i)/rgc_vect(i)
!            enddo
!          case('zsweep') 
!            rgc_scal = sqrt(x(i1)**2+y(i2)**2)
!            if(rgc_scal.le.3.0*kpc) then
!              omega=flat
!              g_shear_rate=0.0
!            else
!              if(rgc_scal.ge.5.0*kpc) then
!                omega=vsun/rgc_scal
!                g_shear_rate=(-1.0)*vsun/rgc_scal
!              else
!                omega=aconst*rgc_scal+bconst*log(rgc_scal)+cconst
!                g_shear_rate=aconst*rgc_scal+bconst
!              endif
!            endif
!            grav=-4.4e-9*cm/sek**2*exp((-1.0)*(rgc_scal-r_gc_sun)/(4.9*kpc))*xsw/sqrt(xsw**2+(0.2*kpc)**2)
!            grav=grav-1.7e-9*cm/sek**2*(r_gc_sun**2+(2.2*kpc)**2)/(rgc_scal**2+(2.2*kpc)**2)*xsw/kpc
!            grav=grav+2*omega*(omega+g_shear_rate)*xsw
!        end select

! galactic case as in ferriere'98, but without "omega" in z component
#elif GRAV ==  'galwomeg'
!      case ('galwomeg')
        aconst=-vsun/(5000.0*pc*(5000.0-3000.0)*pc) ![pc/s/pc2]
        bconst=-3000.0*pc*aconst ![1/s]
        cconst=vsun/5000.0*pc-aconst*5000.0*pc-bconst*log(5000.0*pc) ![1/s]
        flat=aconst*3000.0*pc+bconst*log(3000.0*pc)+cconst ![1/s]
        select case (sweep)
          case('xsweep')
            do i=1,n 
              rgc_vect(i) = sqrt(xsw(i)**2+y(i1)**2)
              if(rgc_vect(i).le.3.0*kpc) then
                omega=flat
              else
                if(rgc_vect(i).ge.5.0*kpc) then
                  omega=vsun/rgc_vect(i)
                else
                  omega=aconst*rgc_vect(i)+bconst*log(rgc_vect(i))+cconst
                endif
              endif
              grav(i)=(-1.0)*omega**2*rgc_vect(i)*xsw(i)/rgc_vect(i)
            enddo
          case('ysweep') 
            do i=1,n
              rgc_vect(i) = sqrt(x(i2)**2+xsw(i)**2)
              if(rgc_vect(i).le.3.0*kpc) then
                omega=flat
              else
                if(rgc_vect(i).ge.5.0*kpc) then
                  omega=vsun/rgc_vect(i)
                else
                  omega=aconst*rgc_vect(i)+bconst*log(rgc_vect(i))+cconst
                endif
              endif
              grav(i)=(-1.0)*omega**2*rgc_vect(i)*xsw(i)/rgc_vect(i)
            enddo
          case('zsweep') 
            rgc_scal = sqrt(x(i1)**2+y(i2)**2)
            if(rgc_scal.le.3.0*kpc) then
              omega=flat
              g_shear_rate=0.0
            else
              if(rgc_scal.ge.5.0*kpc) then
                omega=vsun/rgc_scal
                g_shear_rate=(-1.0)*vsun/rgc_scal
              else
                omega=aconst*rgc_scal+bconst*log(rgc_scal)+cconst
                g_shear_rate=aconst*rgc_scal+bconst
              endif
            endif
            grav=-4.4e-9*cm/sek**2*exp((-1.0)*(rgc_scal-r_gc_sun)/(4.9*kpc))*xsw/sqrt(xsw**2+(0.2*kpc)**2)
            grav=grav-1.7e-9*cm/sek**2*(r_gc_sun**2+(2.2*kpc)**2)/(rgc_scal**2+(2.2*kpc)**2)*xsw/kpc
!           grav=grav+2*omega*(omega+g_shear_rate)*xsw
          end select

! galactic case as in ferriere'98 but constant within 5kpc radius 
#elif GRAV ==  'galcentr'
!      case ('galcentr')
        aconst=-vsun/(5000.0*pc*(5000.0-3000.0)*pc)
        bconst=-3000.0*pc*aconst
        cconst=vsun/5000.0*pc-aconst*5000.0*pc-bconst*log(5000.0*pc)
        flat=aconst*3000.0*pc+bconst*log(3000.0*pc)+cconst
        select case (sweep)
          case('xsweep')
            do i=1,n 
              rgc_vect(i) = sqrt(xsw(i)**2+y(i1)**2)
              if(rgc_vect(i).le.3.0*kpc) then
                omega=flat
              else
                if(rgc_vect(i).ge.5.0*kpc) then
                  omega=vsun/rgc_vect(i)
                else
                  omega=aconst*rgc_vect(i)+bconst*log(rgc_vect(i))+cconst
                endif
              endif
              grav(i)=(-1.0)*omega**2*rgc_vect(i)*xsw(i)/rgc_vect(i)
              enddo
            case('ysweep') 
              do i=1,n
                rgc_vect(i) = sqrt(x(i2)**2+xsw(i)**2)
                if(rgc_vect(i).le.3.0*kpc) then
                  omega=flat
                else
                  if(rgc_vect(i).ge.5.0*kpc) then
                    omega=vsun/rgc_vect(i)
                  else
                    omega=aconst*rgc_vect(i)+bconst*log(rgc_vect(i))+cconst
                  endif
                endif
                grav(i)=(-1.0)*omega**2*rgc_vect(i)*xsw(i)/rgc_vect(i)
              enddo
            case('zsweep') 
              rgc_scal = sqrt(x(i1)**2+y(i2)**2)
              if(rgc_scal.ge.5.0*kpc) then
                omega=vsun/rgc_scal
                g_shear_rate=(-1.0)*vsun/rgc_scal
                grav=-4.4e-9*cm/sek**2*exp((-1.0)*(rgc_scal-r_gc_sun)/(4.9*kpc))*xsw/sqrt(xsw**2+(0.2*kpc)**2)
                grav=grav-1.7e-9*cm/sek**2*(r_gc_sun**2+(2.2*kpc)**2)/(rgc_scal**2+(2.2*kpc)**2)*xsw/kpc
              else
                rgc_scal=5.0*kpc
                omega=aconst*rgc_scal+bconst*log(rgc_scal)+cconst
                g_shear_rate=aconst*rgc_scal+bconst
                grav=-4.4e-9*cm/sek**2*exp((-1.0)*(rgc_scal-r_gc_sun)/(4.9*kpc))*xsw/sqrt(xsw**2+(0.2*kpc)**2)
                grav=grav-1.7e-9*cm/sek**2*(r_gc_sun**2+(2.2*kpc)**2)/(rgc_scal**2+(2.2*kpc)**2)*xsw/kpc
              endif
          end select

! galactic case with z component constant as for sun galactic radius
#elif GRAV ==  'gallksun'
!      case ('gallksun')
        aconst=-vsun/(5000.0*pc*(5000.0-3000.0)*pc)
        bconst=-3000.0*pc*aconst
        cconst=vsun/5000.0*pc-aconst*5000.0*pc-bconst*log(5000.0*pc)
        flat=aconst*3000.0*pc+bconst*log(3000.0*pc)+cconst
        select case (sweep)
          case('xsweep')
            do i=1,n 
              rgc_vect(i) = sqrt(xsw(i)**2+y(i1)**2)
              if(rgc_vect(i).le.3.0*kpc) then
                omega=flat
              else
                if(rgc_vect(i).ge.5.0*kpc) then
                  omega=vsun/rgc_vect(i)
                else
                  omega=aconst*rgc_vect(i)+bconst*log(rgc_vect(i))+cconst
                endif
              endif
              grav(i)=(-1.0)*omega**2*rgc_vect(i)*xsw(i)/rgc_vect(i)
            enddo
          case('ysweep') 
            do i=1,n
              rgc_vect(i) = sqrt(x(i2)**2+xsw(i)**2)
              if(rgc_vect(i).le.3.0*kpc) then
                omega=flat
              else
                if(rgc_vect(i).ge.5.0*kpc) then
                  omega=vsun/rgc_vect(i)
                else
                  omega=aconst*rgc_vect(i)+bconst*log(rgc_vect(i))+cconst
                endif
              endif
              grav(i)=(-1.0)*omega**2*rgc_vect(i)*xsw(i)/rgc_vect(i)
            enddo
          case('zsweep') 
            rgc_scal = 8.5*kpc
            omega=vsun/rgc_scal
            g_shear_rate=(-1.0)*vsun/rgc_scal
            grav=-4.4e-9*cm/sek**2*exp((-1.0)*(rgc_scal-r_gc_sun) &
                    /(4.9*kpc))*xsw/sqrt(xsw**2+(0.2*kpc)**2)
            grav=grav-1.7e-9*cm/sek**2*(r_gc_sun**2+(2.2*kpc)**2) &
                    /(rgc_scal**2+(2.2*kpc)**2)*xsw/kpc
            grav=grav+2*omega*(omega+g_shear_rate)*xsw*3.0856e18
        end select


! galactic case for pseudo-2d simulations
#elif GRAV ==  'gal_flat'
!      case ('gal_flat')
        aconst=-vsun/(5000.0*pc*(5000.0-3000.0)*pc)
        bconst=-3000.0*pc*aconst
        cconst=vsun/5000.0*pc-aconst*5000.0*pc-bconst*log(5000.0*pc)
        flat=aconst*3000.0*pc+bconst*log(3000.0*pc)+cconst
          select case (sweep)
            case('xsweep')
              do i=1,n 
                rgc_vect(i) = sqrt(xsw(i)**2+y(i1)**2)
                if(rgc_vect(i).le.3.0*kpc) then
                  omega=flat
                else
                  if(rgc_vect(i).ge.5.0*kpc) then
                    omega=vsun/rgc_vect(i)
                  else
                    omega=aconst*rgc_vect(i)+bconst*log(rgc_vect(i))+cconst
                  endif
                endif
                grav(i)=(-1.0)*omega**2*rgc_vect(i)*xsw(i)/rgc_vect(i)
              enddo
            case('ysweep') 
              do i=1,n
                rgc_vect(i) = sqrt(x(i2)**2+xsw(i)**2)
                if(rgc_vect(i).le.3.0*kpc) then
                  omega=flat
                else
                  if(rgc_vect(i).ge.5.0*kpc) then
                    omega=vsun/rgc_vect(i)
                  else
                    omega=aconst*rgc_vect(i)+bconst*log(rgc_vect(i))+cconst
                  endif
                endif
                grav(i)=(-1.0)*omega**2*rgc_vect(i)*xsw(i)/rgc_vect(i)
              enddo
             case('zsweep') 
               grav=0.0
           end select
!dw-

#elif GRAV ==  'null'
       grav=0.0
#else
#error 'GRAV declared, but gravity model undefined'
#endif

    if (n_gravh .ne. 0) then
      grav(:)=grav(:)/cosh((xsw(:)/h_grav)**n_gravh )
    endif      
     
  end subroutine grav_accel

!--------------------------------------------------------------------------
   
!  subroutine grav_profile
!
!    implicit none 
!
!        if(n_gravh .ne. 0) then
!          gravity_profile = 1./cosh((z/(h_grav))**n_gravh)
!        else
!         gravity_profile(:) = 1.
!       endif
!       
!       write(*,*) gravity_profile
! 
!  end subroutine  grav_profile
!
!*****************************************************************************  
 
 
  subroutine grav_pot2accel(sweep, i1,i2, n, grav)
  
    implicit none
    character, intent(in)          :: sweep*6
    integer, intent(in)            :: i1, i2                   
    integer, intent(in)            :: n
    real, dimension(n),intent(out) :: grav
  
! Gravitational accelleration is computed on right cell boundaries 
  
    select case(sweep)
      case('xsweep')
        grav(1:n-1) = (gp(1:n-1,i1,i2) - gp(2:n,i1,i2))/dl(xdim)
      case('ysweep')
        grav(1:n-1) = (gp(i2,1:n-1,i1) - gp(i2,2:n,i1))/dl(ydim)
      case('zsweep')
        grav(1:n-1) = (gp(i1,i2,1:n-1) - gp(i1,i2,2:n))/dl(zdim)
    end select
  
  end subroutine grav_pot2accel

!--------------------------------------------------------------------------

  subroutine grav_accel2pot
  
    implicit none
    integer               :: i, j, k, ip, pgpmax
    real gravrx(nx), gravry(ny), gravrz(nz)
    real                  :: gp_max, dgp(3), dgp0
    integer, dimension(3) :: loc_gp_max, gploc
    integer               :: proc_gp_max
    integer               :: px, py, pz, pc(3)
    real dgpx_proc, dgpx_all(0:nproc-1), &
         dgpy_proc, dgpy_all(0:nproc-1), &
         dgpz_proc, dgpz_all(0:nproc-1), &
         dgpx(0:pxsize-1,0:pysize-1,0:pzsize-1), &
         dgpy(0:pxsize-1,0:pysize-1,0:pzsize-1), &
         dgpz(0:pxsize-1,0:pysize-1,0:pzsize-1), &
         ddgp(0:pxsize-1,0:pysize-1,0:pzsize-1)
    
    
! Gravitational potential gp(i,j,k) is defined in cell centers  
! First we find gp(:,:,:) in each block separately assuming:  
    
    gp(1,1,1) = 0.0

    call grav_accel("xsweep", 1, 1, xr(:), nx, gravrx)    
    do i = 1, nx-1      
      gp(i+1,1,1) = gp(i,1,1) - gravrx(i)*dl(xdim)
    enddo
    
    do i=1,nx
      call grav_accel('ysweep',1, i, yr(:), ny, gravry)
      do j = 1, ny-1      
        gp(i,j+1,1) = gp(i,j,1) - gravry(j)*dl(ydim)
      enddo
    enddo

    do i=1,nx
      do j=1,ny
        call grav_accel('zsweep', i, j, zr(:), nz, gravrz)
        do k = 1, nz-1      
          gp(i,j,k+1) = gp(i,j,k) - gravrz(k)*dl(zdim)
        enddo
      enddo
    enddo
            
    dgpx_proc = gp(is,1,1)-gp(1,1,1)
    dgpy_proc = gp(1,js,1)-gp(1,1,1)
    dgpz_proc = gp(1,1,ks)-gp(1,1,1)

    call MPI_GATHER ( dgpx_proc, 1, MPI_DOUBLE_PRECISION, &
                      dgpx_all,  1, MPI_DOUBLE_PRECISION, &
                      0, comm3d,err )
    call MPI_GATHER ( dgpy_proc, 1, MPI_DOUBLE_PRECISION, &
                      dgpy_all,  1, MPI_DOUBLE_PRECISION, &
                      0, comm3d,err )
    call MPI_GATHER ( dgpz_proc, 1, MPI_DOUBLE_PRECISION, &
                      dgpz_all,  1, MPI_DOUBLE_PRECISION, &
                      0, comm3d,err )
                      

    if(proc .eq. 0) then
    
      do ip = 0, nproc-1
        call MPI_CART_COORDS(comm3d, ip, ndims, pc, ierr)
       
        px = pc(1)
        py = pc(2)
        pz = pc(3)
      
        dgpx(px,py,pz) = dgpx_all(ip)
        dgpy(px,py,pz) = dgpy_all(ip)
        dgpz(px,py,pz) = dgpz_all(ip)

      enddo
    
      ddgp(0,0,0) = 0.0
      do i = 1, pxsize-1      
        ddgp(i,0,0) = ddgp(i-1,0,0) + dgpx(i-1,0,0)
      enddo
   
      do i=0, pxsize-1
        do j = 1, pysize-1      
          ddgp(i,j,0) = ddgp(i,j-1,0) + dgpy(i,j-1,0)
        enddo
      enddo

      do i=0, pxsize-1
        do j=0, pysize-1
          do k = 1, pzsize-1      
             ddgp(i,j,k)= ddgp(i,j,k-1) + dgpz(i,j,k-1)
          enddo
        enddo
      enddo
    
    endif

    call MPI_BCAST(ddgp, nproc, MPI_DOUBLE_PRECISION, 0, comm, ierr)       

    px = pcoords(1)
    py = pcoords(2)
    pz = pcoords(3)
    
    gp = gp + ddgp(px,py,pz)
      
    gp_max      = maxval(gp(is:ie,js:je,ks:ke)) 
    loc_gp_max  = maxloc(gp(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
                  
    call mpifind(gp_max, 'max', loc_gp_max, proc_gp_max)
    pgpmax = proc_gp_max
    
    call MPI_BCAST(gp_max, 1, MPI_DOUBLE_PRECISION, pgpmax, comm, ierr)       
    gp = gp - gp_max
   
  end subroutine grav_accel2pot


end module gravity
