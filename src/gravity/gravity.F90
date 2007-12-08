#include "mhd.def"

module gravity

  use constants


  character gp_status*9
#ifdef GALACTIC_DISK
    real, allocatable :: gpotdisk(:,:,:),gpothalo(:,:,:),gpotbulge(:,:,:)
#endif GALACTIC_DISK

contains

!--------------------------------------------------------------------------
  subroutine grav_pot_3d
    use arrays, only : nx,ny,nz,gp,z
    implicit none
    integer i, j
    real, allocatable :: gpot(:)
    allocate(gpot(nz))

#ifdef GALACTIC_DISK
allocate(gpotdisk(nx,ny,nz),gpothalo(nx,ny,nz),gpotbulge(nx,ny,nz))
#endif GALACTIC_DISK

    do j = 1,ny
      do i = 1,nx
        call grav_pot('zsweep', i,j, z(:), nz, gpot, gp_status)
        if(gp_status .eq. 'undefined') exit
        gp(i,j,:) = gpot
      enddo
    enddo
    
    if(gp_status .eq. 'undefined') then
      call grav_accel2pot('default')
#ifdef GALACTIC_DISK
        call grav_accel2pot('diskprt')
	call grav_accel2pot('haloprt')
	call grav_accel2pot('blgpart')
#endif GALACTIC_DISK
    endif

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
    
    real, dimension(n)	  :: rgc_vect
    real aconst, bconst, cconst, flat, omega, g_shear_rate
    real rgc_scal
    integer i
    real, allocatable :: gpdisk(:), gphalo(:), gpblg(:)
    real Mhalo,Mbulge,Mdisk,ahalo,bbulge,adisk,bdisk

    
#ifdef GRAV_NULL
        gpot = 0.0
      ! do nothing
#elif defined (GRAV_UNIFORM)
        select case (sweep)
	  case('xsweep')
	    gpot = -g_z*z(i2)
	  case('ysweep')
	    gpot = -g_z*z(i1)
	  case('zsweep')
	    gpot = -g_z*xsw
	end select
#elif defined (GRAV_LINEAR)
        select case (sweep)
	  case('xsweep') 
	    gpot = -0.5*dg_dz*z(i2)**2
	  case('ysweep') 
	    gpot = -0.5*dg_dz*z(i1)**2
	  case('zsweep') 
	    gpot = -0.5*dg_dz*xsw**2
        end select
	    
#elif defined (GRAV_PTMASS)
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
#elif defined (GRAV_PTFLAT)
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

! galactic case as in vollmer'01
#elif defined (GRAV_GAL_VOLLMER)
        allocate(gpdisk(n),gphalo(n),gpblg(n))
	Mhalo = 8.6e10*Msun 
	Mbulge= 5.6e9*Msun
	Mdisk = 2.6e10*Msun
	ahalo = 12.*kpc
	bbulge= 387.*pc
	adisk = 2.7*kpc
	bdisk = 250.*pc
        select case (sweep)
          case('xsweep')
            do i=1,n 
              rgc_vect(i) = sqrt(xsw(i)**2+y(i1)**2)
	    gpot(i)=-Mdisk/sqrt(rgc_vect(i)**2+(adisk+sqrt(z(i2)**2+bdisk**2))**2)
	    gpot(i)=gpot(i)-1./(rgc_vect(i)**2+z(i2)**2)*Mhalo*(sqrt(rgc_vect(i)**2+z(i2)**2)/ahalo)**2.02/(1.+(sqrt(rgc_vect(i)**2+z(i2)**2)/ahalo))
	    gpot(i)=gpot(i)-(Mhalo/1.02/ahalo)*(-1.02/(1.+(100.0*kpc/ahalo)**1.02)+log(1.0+(100.0*kpc/ahalo)**1.02))
	    gpot(i)=gpot(i)+(Mhalo/1.02/ahalo)*(-1.02/(1.+(sqrt(rgc_vect(i)**2+z(i2)**2)/ahalo)**1.02)+log(1.0+(sqrt(rgc_vect(i)**2+z(i2)**2)/ahalo)**1.02))
	    gpot(i)=gpot(i)-Mbulge/sqrt(rgc_vect(i)**2+z(i2)**2+bbulge**2)
            enddo
          case('ysweep') 
            do i=1,n
              rgc_vect(i) = sqrt(x(i2)**2+xsw(i)**2)
	    gpot(i)=-Mdisk/sqrt(rgc_vect(i)**2+(adisk+sqrt(z(i1)**2+bdisk**2))**2)
	    gpot(i)=gpot(i)-1./(rgc_vect(i)**2+z(i1)**2)*Mhalo*(sqrt(rgc_vect(i)**2+z(i1)**2)/ahalo)**2.02/(1.+(sqrt(rgc_vect(i)**2+z(i1)**2)/ahalo))
	    gpot(i)=gpot(i)-(Mhalo/1.02/ahalo)*(-1.02/(1.+(100.0*kpc/ahalo)**1.02)+log(1.0+(100.0*kpc/ahalo)**1.02))
	    gpot(i)=gpot(i)+(Mhalo/1.02/ahalo)*(-1.02/(1.+(sqrt(rgc_vect(i)**2+z(i1)**2)/ahalo)**1.02)+log(1.0+(sqrt(rgc_vect(i)**2+z(i1)**2)/ahalo)**1.02))
	    gpot(i)=gpot(i)-Mbulge/sqrt(rgc_vect(i)**2+z(i1)**2+bbulge**2)
            enddo
          case('zsweep') 
            rgc_scal = sqrt(x(i1)**2+y(i2)**2)
	    gpdisk=-Mdisk/sqrt(rgc_scal**2+(adisk+sqrt(xsw**2+bdisk**2))**2)
	    gphalo=-1./(rgc_scal**2+xsw**2)*Mhalo*(sqrt(rgc_scal**2+xsw**2)/ahalo)**2.02/(1.+(sqrt(rgc_scal**2+xsw**2)/ahalo))
	    gphalo=gphalo-(Mhalo/1.02/ahalo)*(-1.02/(1.+(100.0*kpc/ahalo)**1.02)+log(1.0+(100.0*kpc/ahalo)**1.02))
	    gphalo=gphalo+(Mhalo/1.02/ahalo)*(-1.02/(1.+(sqrt(rgc_scal**2+xsw**2)/ahalo)**1.02)+log(1.0+(sqrt(rgc_scal**2+xsw**2)/ahalo)**1.02))
	    gpblg=-Mbulge/sqrt(rgc_scal**2+xsw**2+bbulge**2)
            gpot=gpdisk+gpblg+gphalo
	    gpotdisk(i1,i2,:)=gpdisk
	    gpothalo(i1,i2,:)=gphalo
	    gpotbulge(i1,i2,:)=gpblg
        end select
	    deallocate(gpdisk,gphalo,gpblg)

! galactic case as in ferriere'98
#elif defined (GRAV_GAL_FERRIERE)
        allocate(gpdisk(n),gphalo(n),gpblg(n))
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
                g_shear_rate=0.0
              else
                if(rgc_vect(i).ge.5.0*kpc) then
                  omega=vsun/rgc_vect(i)
                  g_shear_rate=(-1.0)*vsun/rgc_scal
                else
                  omega=aconst*rgc_vect(i)+bconst*log(rgc_vect(i))+cconst
                  g_shear_rate=aconst*rgc_scal+bconst
                endif
              endif
	      gpot(i)=4.4e-9*cm/sek**2*exp((-1.0)*(rgc_vect(i)-r_gc_sun)/(4.9*kpc))*(sqrt(z(i2)**2+(0.2*kpc)**2)-sqrt((0.2*kpc)**2))
	      gpot(i)=gpot(i)+1.7e-9*cm/sek**2*(r_gc_sun**2+(2.2*kpc)**2)/(rgc_vect(i)**2+(2.2*kpc)**2)*z(i2)**2/2.0/kpc
	      gpot(i)=gpot(i)+omega*(omega+g_shear_rate)*z(i2)**2
	      if(rgc_vect(i).le.3.0*kpc) then
	        gpot(i)=gpot(i)+0.5*flat**2*rgc_vect(i)**2
	      else
	        if(rgc_vect(i).ge.5.0*kpc) then
	          gpot(i)=gpot(i)+0.5*flat**2*(3.0*kpc)**2 &
		      +8.0*aconst*kpc**2+bconst*(5.0*kpc*log(5.0*kpc/e)-3.0*kpc*log(3.0*kpc/e))+cconst*2.0*kpc &
		      +vsun*log(rgc_vect(i)/5.0/kpc)
	        else
	          gpot(i)=gpot(i)+0.5*flat**2*(3.0*kpc)**2 +0.5*aconst*(rgc_vect(i)**2-(3.0*kpc)**2) &
		      +bconst*(rgc_vect(i)*log(rgc_vect(i)/e)-3.0*kpc*log(3.0*kpc/e))+cconst*(rgc_vect(i)-3.0*kpc)
	        endif
	      endif
            enddo
          case('ysweep') 
            do i=1,n
              rgc_vect(i) = sqrt(x(i2)**2+xsw(i)**2)
              if(rgc_vect(i).le.3.0*kpc) then
                omega=flat
                g_shear_rate=0.0
              else
                if(rgc_vect(i).ge.5.0*kpc) then
                  omega=vsun/rgc_vect(i)
                  g_shear_rate=(-1.0)*vsun/rgc_scal
                else
                  omega=aconst*rgc_vect(i)+bconst*log(rgc_vect(i))+cconst
                g_shear_rate=aconst*rgc_scal+bconst
                endif
              endif
	      gpot(i)=4.4e-9*cm/sek**2*exp((-1.0)*(rgc_vect(i)-r_gc_sun)/(4.9*kpc))*(sqrt(z(i1)**2+(0.2*kpc)**2)-sqrt((0.2*kpc)**2))
	      gpot(i)=gpot(i)+1.7e-9*cm/sek**2*(r_gc_sun**2+(2.2*kpc)**2)/(rgc_vect(i)**2+(2.2*kpc)**2)*z(i1)**2/2.0/kpc
	      gpot(i)=gpot(i)+omega*(omega+g_shear_rate)*z(i1)**2
	      if(rgc_vect(i).le.3.0*kpc) then
	        gpot(i)=gpot(i)+0.5*flat**2*rgc_vect(i)**2
	      else
	        if(rgc_vect(i).ge.5.0*kpc) then
	          gpot(i)=gpot(i)+0.5*flat**2*(3.0*kpc)**2 &
		      +8.0*aconst*kpc**2+bconst*(5.0*kpc*log(5.0*kpc/e)-3.0*kpc*log(3.0*kpc/e))+cconst*2.0*kpc &
		      +vsun*log(rgc_vect(i)/5.0/kpc)
	        else
	          gpot(i)=gpot(i)+0.5*flat**2*(3.0*kpc)**2 +0.5*aconst*(rgc_vect(i)**2-(3.0*kpc)**2) &
		      +bconst*(rgc_vect(i)*log(rgc_vect(i)/e)-3.0*kpc*log(3.0*kpc/e))+cconst*(rgc_vect(i)-3.0*kpc)
	        endif
	      endif
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
	    gpdisk=4.4e-9*cm/sek**2*exp((-1.0)*(rgc_scal-r_gc_sun)/(4.9*kpc))*(sqrt(xsw**2+(0.2*kpc)**2)-sqrt((0.2*kpc)**2))
	    gphalo=1.7e-9*cm/sek**2*(r_gc_sun**2+(2.2*kpc)**2)/(rgc_scal**2+(2.2*kpc)**2)*xsw**2/2.0/kpc
	    gpblg=omega*(omega+g_shear_rate)*xsw**2
	    if(rgc_scal.le.3.0*kpc) then
	      gpblg=gpblg+0.5*flat**2*rgc_scal**2
	    else
	      if(rgc_scal.ge.5.0*kpc) then
	        gpblg=gpblg+0.5*flat**2*(3.0*kpc)**2 &
		    +8.0*aconst*kpc**2+bconst*(5.0*kpc*log(5.0*kpc/e)-3.0*kpc*log(3.0*kpc/e))+cconst*2.0*kpc &
		    +vsun*log(rgc_scal/5.0/kpc)
	      else
	        gpblg=gpblg+0.5*flat**2*(3.0*kpc)**2 +0.5*aconst*(rgc_scal**2-(3.0*kpc)**2) &
		    +bconst*(rgc_scal*log(rgc_scal/e)-3.0*kpc*log(3.0*kpc/e))+cconst*(rgc_scal-3.0*kpc)
	      endif
	    endif
            gpot=gpdisk+gpblg+gphalo
	    gpotdisk(i1,i2,:)=gpdisk
	    gpothalo(i1,i2,:)=gphalo
	    gpotbulge(i1,i2,:)=gpblg
        end select
	    deallocate(gpdisk,gphalo,gpblg)


#else
        status = 'undefined'
!#warning 'GRAV declared, but gravity model undefined in grav_pot'
! niektore modele grawitacji realizowane sa za pomoca przyspieszenia 
! (np. 'galactic') z ktorego liczony jest potencjal
#endif

  end subroutine grav_pot

!--------------------------------------------------------------------------
   
  subroutine grav_accel(sweep, i1,i2, xsw, n, grav, gravpart)
    use start, only  : h_grav, n_gravh,nb
    use arrays, only : gp,x,y,z,dl,xdim,ydim,zdim,nx,ny,nz, &
      is,ie,js,je,ks,ke, xr,yr,zr
   
    implicit none
    character, intent(in) :: sweep*6, gravpart*7
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


#ifdef GRAV_UNIFORM
!      case ('uniform')
       if (sweep == 'zsweep') then
          grav = g_z
       else
          grav = 0.0
       endif
#elif defined (GRAV_LINEAR)
       if (sweep == 'zsweep') then
          grav = dg_dz*xsw
       else
          grav = 0.0
       endif
#elif defined (GRAV_PTMASS)
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
#elif defined (GRAV_PTFLAT)
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
#elif defined (GRAV_GAL_MALI)
       if (sweep == 'zsweep') then
         grav = -0.87 * (6.8*tanh(3.2 * xsw/(pc*kpc)) + 1.7*xsw/(pc*kpc)) * 1.e-9
       else
         grav = 0.0
       endif
       
#elif defined (GRAV_GALACTIC)
       if (sweep == 'zsweep') then
	    select case(gravpart)
	      case('diskprt')
	        grav = cmps2*(-4.4e-9 * exp(-(r_gc-r_gc_sun)/(4.9*kpc)) * xsw/sqrt(xsw**2+(0.2*kpc)**2))
	      case('haloprt')
	        grav = cmps2*(-(1.7e-9 * (r_gc_sun**2 + (2.2*kpc)**2)/(r_gc**2 + (2.2*kpc)**2)*xsw/kpc) )
	      case('blgpart')
	        grav = 0.0
	      case('default')
         grav = cmps2 * (  &
           (-4.4e-9 * exp(-(r_gc-r_gc_sun)/(4.9*kpc)) * xsw/sqrt(xsw**2+(0.2*kpc)**2)) &    
           -( 1.7e-9 * (r_gc_sun**2 + (2.2*kpc)**2)/(r_gc**2 + (2.2*kpc)**2)*xsw/kpc) )     
!          -Om*(Om+G) * Z * (kpc ?) ! in the transition region between rigid 
!                                   ! and flat rotation F'98: eq.(36)       
        end select
       else
         grav=0.0
       endif

!dw+
! galactic case as in ferriere'98
#elif defined (GRAV_GALACTIC_DW)
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
          if((gravpart .eq. 'default') .or. (gravpart .eq. 'diskprt')) then
            grav=-4.4e-9*cm/sek**2*exp((-1.0)*(rgc_scal-r_gc_sun)/(4.9*kpc))*xsw/sqrt(xsw**2+(0.2*kpc)**2)
          endif
	  if((gravpart .eq. 'default') .or. (gravpart .eq. 'haloprt')) then
            grav=grav-1.7e-9*cm/sek**2*(r_gc_sun**2+(2.2*kpc)**2)/(rgc_scal**2+(2.2*kpc)**2)*xsw/kpc
	    endif
	  if((gravpart .eq. 'default') .or. (gravpart .eq. 'blgpart')) then
            grav=grav+2*omega*(omega+g_shear_rate)*xsw
	  endif

        else  ! y or x
          if(sweep == 'xsweep') then 
            rgc_vect(:) = sqrt(xsw(:)**2 +   y(i1)**2)
          elseif(sweep == 'ysweep') then 
            rgc_vect(:) = sqrt(x(i2)**2  +  xsw(:)**2)
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

! galactic case for pseudo-2d simulations
#elif defined (GRAV_GAL_FLAT)
        aconst = -vsun / (5000.0 * pc *(5000.0-3000.0) * pc)            ![pc/s/pc2]
        bconst = -3000.0 * pc * aconst                                  ![1/s]
        cconst =  vsun/5000.0*pc-aconst*5000.0*pc-bconst*log(5000.0*pc) ![1/s]
        flat   = aconst * 3000.0 * pc + bconst * log(3000.0*pc)+ cconst ![1/s]

        if (sweep == 'zsweep') then
            grav = 0.0
          
        else  ! y or x
          if(sweep == 'xsweep') then 
            rgc_vect(:) = sqrt(xsw(:)**2 +   y(i1)**2)
          elseif(sweep == 'ysweep') then 
            rgc_vect(:) = sqrt(x(i2)**2  +  xsw(:)**2)
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

! galactic case as in ferriere'98 but with smoothed omega
#elif defined (GRAV_GALSMOOTH)
        if (sweep == 'zsweep') then
	    rgc_scal = sqrt(x(i1)**2+y(i2)**2)
	    omega = vsun * tanh(rgc_scal/3.0/kpc)/rgc_scal
	    g_shear_rate = vsun/rgc_scal * ((1.0-(tanh(rgc_scal/3.0/kpc))**2)/(3.0*kpc)-tanh(rgc_scal/3.0/kpc))
	    if((gravpart .eq. 'default') .or. (gravpart .eq. 'diskprt')) then
              grav=-4.4e-9*cm/sek**2*exp((-1.0)*(rgc_scal-r_gc_sun)/(4.9*kpc))*xsw/sqrt(xsw**2+(0.2*kpc)**2)
	    endif
	    if((gravpart .eq. 'default') .or. (gravpart .eq. 'haloprt')) then
              grav=grav-1.7e-9*cm/sek**2*(r_gc_sun**2+(2.2*kpc)**2)/(rgc_scal**2+(2.2*kpc)**2)*xsw/kpc
	    endif
	    if((gravpart .eq. 'default') .or. (gravpart .eq. 'blgpart')) then
              grav=grav+2*omega*(omega+g_shear_rate)*xsw
	    endif
        else  ! y or x
          if(sweep == 'xsweep') then 
            rgc_vect(:) = sqrt(xsw(:)**2 +   y(i1)**2)
          elseif(sweep == 'ysweep') then 
            rgc_vect(:) = sqrt(x(i2)**2  +  xsw(:)**2)
          endif
          grav = -1.0 * xsw(:) * vsun*tanh(rgc_vect(:)/3.0/kpc)/rgc_vect(:)
	endif
! galactic case as in ferriere'98 but with smoothed omega with decreasing rotation with z
#elif defined (GRAV_GALSMTDEC)
        if (sweep == 'zsweep') then
	    rgc_scal = sqrt(x(i1)**2+y(i2)**2)
	    omega=vsun*tanh(rgc_scal/3.0/kpc)/rgc_scal !*exp(-xsw**2/2.0/h_grav**n_gravh)
	    g_shear_rate=(vsun/rgc_scal*((1.0-(tanh(rgc_scal/3.0/kpc))**2)/(3.0*kpc)-tanh(rgc_scal/3.0/kpc)))!*exp(-xsw**2/2.0/h_grav**n_gravh)
	    if((gravpart .eq. 'default') .or. (gravpart .eq. 'diskprt')) then
              grav=-4.4e-9*cm/sek**2*exp((-1.0)*(rgc_scal-r_gc_sun)/(4.9*kpc))*xsw/sqrt(xsw**2+(0.2*kpc)**2)
	    endif
	    if((gravpart .eq. 'default') .or. (gravpart .eq. 'haloprt')) then
              grav=grav-1.7e-9*cm/sek**2*(r_gc_sun**2+(2.2*kpc)**2)/(rgc_scal**2+(2.2*kpc)**2)*xsw/kpc
	    endif
	    if((gravpart .eq. 'default') .or. (gravpart .eq. 'blgpart')) then
              grav=grav+2*omega*(omega+g_shear_rate)*xsw*(exp(-xsw**2/2.0/h_grav**n_gravh))**2
	    endif
        else  ! y or x
          if(sweep == 'xsweep') then 
            rgc_vect(:) = sqrt(xsw(:)**2 +   y(i1)**2)
            grav(:)=vsun*tanh(rgc_vect(:)/3.0/kpc)/rgc_vect(:)*exp(-z(i2)**2/2.0/h_grav**n_gravh)
            grav(:)=grav(:)**2*(-1.0)*rgc_vect(:)*xsw(:)/rgc_vect(:)
          elseif(sweep == 'ysweep') then 
            rgc_vect(:) = sqrt(x(i2)**2  +  xsw(:)**2)
	    grav(:)=vsun*tanh(rgc_vect(i)/3000.0/pc)/rgc_vect(i)*exp(-z(i1)**2/2.0/h_grav**n_gravh)
            grav(:)=grav(:)**2*(-1.0)*rgc_vect(:)*xsw(:)/rgc_vect(:)
          endif
	endif
! galactic case as in ferriere'98 but with smoothed omega with decreasing rotation with x (x axis symmetry)
#elif defined (GRAV_GALDCZTOX)
        if (sweep == 'xsweep') then
          case('xsweep') 
            grav(:)=0.0
	    rgc_scal = sqrt(z(i2)**2+y(i1)**2)
	    omega=vsun*tanh(rgc_scal/3.0/kpc)/rgc_scal !*exp(-xsw**2/2.0/h_grav**n_gravh)
	    g_shear_rate=(vsun/rgc_scal*((1.0-(tanh(rgc_scal/3.0/kpc))**2)/(3.0*kpc)-tanh(rgc_scal/3.0/kpc)))!*exp(-xsw**2/2.0/h_grav**n_gravh)
	    if((gravpart .eq. 'default') .or. (gravpart .eq. 'diskprt')) then
              grav=-4.4e-9*cm/sek**2*exp((-1.0)*(rgc_scal-r_gc_sun)/(4.9*kpc))*xsw/sqrt(xsw**2+(0.2*kpc)**2)
	    endif
	    if((gravpart .eq. 'default') .or. (gravpart .eq. 'haloprt')) then
              grav=grav-1.7e-9*cm/sek**2*(r_gc_sun**2+(2.2*kpc)**2)/(rgc_scal**2+(2.2*kpc)**2)*xsw/kpc
	    endif
	    if((gravpart .eq. 'default') .or. (gravpart .eq. 'blgpart')) then
              grav=grav+2*omega*(omega+g_shear_rate)*xsw*(exp(-xsw**2/2.0/h_grav**n_gravh))**2
	    endif
        else  ! y or z
          if(sweep == 'zsweep') then 
            rgc_vect(:) = sqrt(xsw(:)**2+y(i2)**2)
	    grav(:)=vsun*tanh(rgc_vect(:)/3.0/kpc)/rgc_vect(:)*exp(-x(i1)**2/2.0/h_grav**n_gravh)
            grav(:)=grav(:)**2*(-1.0)*rgc_vect(:)*xsw(:)/rgc_vect(:)
          elseif(sweep == 'ysweep') then 
            rgc_vect(:) = sqrt(z(i1)**2+xsw(:)**2)
	    grav(:)=vsun*tanh(rgc_vect(:)/3000.0/pc)/rgc_vect(:)*exp(-x(i2)**2/2.0/h_grav**n_gravh)
            grav(:)=grav(:)**2*(-1.0)*rgc_vect(:)*xsw(:)/rgc_vect(:)
            rgc_vect(:) = sqrt(x(i2)**2  +  xsw(:)**2)
	    grav(:)=vsun*tanh(rgc_vect(i)/3000.0/pc)/rgc_vect(i)*exp(-z(i1)**2/2.0/h_grav**n_gravh)
            grav(:)=grav(:)**2*(-1.0)*rgc_vect(:)*xsw(:)/rgc_vect(:)
          endif
	endif
! galactic case as in ferriere'98 but with smoothed omega and increased at large radii
#elif defined (GRAV_GALSMINCR)
        if (sweep == 'zsweep') then
	    rgc_scal = sqrt(x(i1)**2+y(i2)**2)
	    omega=vsun*tanh(rgc_scal/3.0/kpc)/rgc_scal
	    g_shear_rate=vsun/rgc_scal*((1.0-(tanh(rgc_scal/3.0/kpc))**2)/(3.0*kpc)-tanh(rgc_scal/3.0/kpc))
	    if((gravpart .eq. 'default') .or. (gravpart .eq. 'diskprt')) then
              grav=-4.4e-9*cm/sek**2*exp((-1.0)*(rgc_scal-r_gc_sun)/(4.9*kpc))*xsw/sqrt(xsw**2+(0.2*kpc)**2)
	    endif
	    if((gravpart .eq. 'default') .or. (gravpart .eq. 'haloprt')) then
              grav=grav-1.7e-9*cm/sek**2*(r_gc_sun**2+(2.2*kpc)**2)/(rgc_scal**2+(2.2*kpc)**2)*xsw/kpc
	    endif
	    if((gravpart .eq. 'default') .or. (gravpart .eq. 'blgpart')) then
              grav=grav+2*omega*(omega+g_shear_rate)*xsw
!	      grav=grav*(2.-1./cosh((rc/r_grav)**n_gravr))
            endif
        else  ! y or x
          if(sweep == 'xsweep') then 
            rgc_vect(:) = sqrt(xsw(:)**2+y(i1)**2)
	    grav(:)=vsun*tanh(rgc_vect(:)/3.0/kpc)/rgc_vect(:)
            grav(:)=grav(:)**2*(-1.0)*xsw(i)
	    grav(:)=grav(:)*(2.-1./cosh((rgc_vect(:)/r_grav)**n_gravr))
          elseif(sweep == 'ysweep') then 
            rgc_vect(:) = sqrt(x(i2)**2+xsw(:)**2)
	    grav(:)=vsun*tanh(rgc_vect(:)/3000.0/pc)/rgc_vect(:)
            grav(:)=grav(:)**2*(-1.0)*xsw(:)
	    grav(:)=grav(:)*(2.-1./cosh((rgc_vect(:)/r_grav)**n_gravr))
          endif
	endif

!dw-

#elif defined (GRAV_NULL)
       grav=0.0
#else
!#error 'GRAV declared, but gravity model undefined'
#endif

#ifndef GRAV == 'galsmtdec'
#ifndef GRAV == 'galsmincr'
    if (n_gravh .ne. 0) then
      grav(:)=grav(:)/cosh((xsw(:)/h_grav)**n_gravh )
    endif      
#endif GRAV == 'galsmincr'
#endif GRAV == 'galsmtdec'
     
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
    use arrays, only : gp, dl, xdim, ydim, zdim
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

  subroutine grav_accel2pot(gravpart)
    use mpi_setup
    use arrays, only : gp,dl,xdim,ydim,zdim,is,ie,js,je,ks,ke,nx,ny,nz, &
      zr,yr,xr
    use start, only  : nb

    implicit none
    character*7 gravpart
    integer               :: i, j, k, ip, pgpmax
    real, allocatable     :: gpwork(:,:,:)
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
! Instead of gp, gpwork is used to not change already computed gp when
! disk, bulge and halo parts of galactic potential are computed.
! gravpart is used to distinguish which part has to be computed.
    
    allocate(gpwork(nx,ny,nz))
    gpwork(1,1,1) = 0.0

    call grav_accel('xsweep',1,1, xr(:), nx, gravrx, gravpart)    
    do i = 1, nx-1      
      gpwork(i+1,1,1) = gpwork(i,1,1) - gravrx(i)*dl(xdim)
    enddo
    
    do i=1,nx
      call grav_accel('ysweep',1, i, yr(:), ny, gravry, gravpart)
      do j = 1, ny-1      
        gpwork(i,j+1,1) = gpwork(i,j,1) - gravry(j)*dl(ydim)
      enddo
    enddo

    do i=1,nx
      do j=1,ny
        call grav_accel('zsweep', i, j, zr(:), nz, gravrz, gravpart)
        do k = 1, nz-1      
          gpwork(i,j,k+1) = gpwork(i,j,k) - gravrz(k)*dl(zdim)
        enddo
      enddo
    enddo
            
    dgpx_proc = gpwork(is,1,1)-gpwork(1,1,1)
    dgpy_proc = gpwork(1,js,1)-gpwork(1,1,1)
    dgpz_proc = gpwork(1,1,ks)-gpwork(1,1,1)

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
    
    gpwork = gpwork + ddgp(px,py,pz)
      
    gp_max      = maxval(gpwork(is:ie,js:je,ks:ke)) 
    loc_gp_max  = maxloc(gpwork(is:ie,js:je,ks:ke)) &
                  + (/nb,nb,nb/)
                  
    call mpifind(gp_max, 'max', loc_gp_max, proc_gp_max)
    pgpmax = proc_gp_max
    
    call MPI_BCAST(gp_max, 1, MPI_DOUBLE_PRECISION, pgpmax, comm, ierr)       
    gpwork = gpwork - gp_max
   
#ifdef GALACTIC_DISK
    select case(gravpart)
      case('default')
        gp=gpwork
      case('diskprt')
        gpotdisk=gpwork
      case('haloprt')
        gpothalo=gpwork
      case('blgpart')
        gpotbulge=gpwork
    end select
#else GALACTIC_DISK
        gp=gpwork
#endif GALACTIC_DISK
    deallocate(gpwork)
   
  end subroutine grav_accel2pot


end module gravity
