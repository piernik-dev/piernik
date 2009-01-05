! $Id$
#include "piernik.def"

module time

! Based on: Pen Arras & Wong (2003)
! Modified extensively by M. Hanasz November 2005 - May 2006
! Optimalization of cpu-time efficiency of mhdflux, tvd1 by K. Kowalik
! Modification history:  see "changelog"  -- warning: out of date

  use mpi_setup

#ifdef COSM_RAYS
  use cr_diffusion
#endif /* COSM_RAYS */

  implicit none
  real c

contains


  subroutine timestep
    use start, only : dt, tend, t, dt_mhd, dt_coolheat, &
         dt_visc, dt_cr
    use constants, only : small,big
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
    use arrays, only : u, idna
#endif /* COLLISIONS || KEPLER_SUPPRESSION */
#ifdef KEPLER_SUPPRESSION
    use arrays, only : alfsup,nfluid,nz,omx0,omy0,x,y,nx,ny,imxa,imya
    use start,  only : dt_supp
#endif /* KEPLER_SUPPRESSION */
    implicit none
#ifdef KEPLER_SUPPRESSION
    real,allocatable,dimension(:,:) :: velx,vely,dvx,dvy
    integer j,k,ifl
#endif /* KEPLER_SUPPRESSION */
! Timestep computation

    call timestep_mhd
    dt=min(dt_mhd,(tend-t)/2.)

#ifdef RESIST
    call timestep_resist
    dt = min(dt,dt_resist)
#endif /* RESIST */
#ifdef COOL_HEAT
    call timestep_coolheat
    dt = min(dt,dt_coolheat)
#endif /* COOL_HEAT */
#ifdef HEAT_COND
    dt_heatcond = cfl_heatcond * 0.5*dxmn**2/(K_heatcond+small)
    dt = min(dt,dt_heatcond)
#endif /* HEAT_COND */
#ifdef VISC
    dt_visc = cfl_visc * 0.5*dxmn**2/(nu_bulk+small)
    dt = min(dt,dt_visc)
#endif /* VISC */
#ifdef COSM_RAYS
    dt_cr = cfl_cr * 0.5*dxmn**2/(K_cr_paral+K_cr_perp+small)
    dt = min(dt,dt_cr)
#endif /* COSM_RAYS */
#ifdef SIMPLE_COOL
    dt = min(dt,0.01 * tauc)
#endif /* SIMPLE_COOL */
#ifdef COLLISIONS
    dt_colls = cfl_colls / (collfaq+small) / maxval(u(idna,:,:,:))
    dt = min(dt,dt_colls)
#endif /* COLLISIONS */
#ifdef KEPLER_SUPPRESSION
!!    allocate(velx(nx,ny),vely(nx,ny),dvx(nx,ny),dvy(nx,ny))
    dt_supp = big
!!    do ifl = 1,nfluid
!!    do k = 1,nz
!!!    dt_supp = min(dt_supp,cfl_colls / maxval(alfsup(:,:)/u(idna(ifl),:,:,k)+small))
!!!    vel1=sqrt(u(imxa(ifl),:,:,k)**2+u(imya(ifl),:,:,k)**2)/u(idna(ifl),:,:,k)
!!    velx=u(imxa(ifl),:,:,k)/u(idna(ifl),:,:,k)
!!#ifdef KEPL_SUPP_SIMX
!!    dvx(:,:)=u(imxa(ifl),:,:,k)-omx0(ifl,:,:,k)
!!else /* KEPL_SUPP_SIMX */
!!    vely=u(imya(ifl),:,:,k)/u(idna(ifl),:,:,k)
!!    do j=1,ny
!!    dvx(:,j)=(velx(:,j)*y(j)-vely(:,j)*x(:))*y(j)/(x(:)**2+y(j)**2)-omx0(ifl,:,j,k)
!!    enddo
!!    do j=1,nx
!!    dvy(j,:)=(vely(j,:)*x(j)-velx(j,:)*y(:))*x(j)/(y(:)**2+x(j)**2)-omy0(ifl,j,:,k)
!!    enddo
!!    dt_supp = min(dt_supp,cfl_colls*minval((abs(vely)+small)/(abs(dvy)+small)/(abs(alfsup)+small)))
!!#endif /* KEPL_SUPP_SIMX */
!!    dt_supp = min(dt_supp,cfl_colls*minval((abs(velx)+small)/(abs(dvx)+small)/(abs(alfsup)+small)))
!!!    dt_supp = min(dt_supp,cfl_colls*sqrt(u(imxa(ifl),:,:,k)**2+u(imya(ifl),:,:,k)**2)/u(idna(ifl),:,:,k)/alfsup(:,:)/
!!    enddo
!!    enddo
!!    deallocate(velx,vely,dvx,dvy)
    do ifl=1,nfluid
      do k=1,nz
        dt_supp = min(dt_supp, minval(abs(1./(alfsup+small)/u(idna(ifl),:,:,k))))
      enddo
    enddo
    dt = min(dt,dt_supp)
#endif /* KEPLER_SUPPRESSION */
  end subroutine timestep

!------------------------------------------------------------------------------------------

  subroutine timestep_mhd
    use grid, only   : dx,dy,dz
    use start, only  : nb,cfl, gamma,dt_mhd, csi2
#ifdef ISO
    use start, only : csi2
#endif /* ISO */
    use arrays, only : ks,ke,nyb,nxb,idna,imxa,imya,imza,u,b
    use arrays, only : iena,nfluid,nadiab,magn
#ifdef DUST
    use dusts, only : idnd
#endif /* DUST */

    implicit none

    real dt_mhd_proc, dt_mhd_all, c_max_all
    real dt_mhd_proc_x, dt_mhd_proc_y, dt_mhd_proc_z
    real cx, cy, cz, vx, vy, vz, cf

! locals
    real pmag
    real v,ps,p
    integer i,j,k,ifluid

    cx    = 0.0
    cy    = 0.0
    cz    = 0.0
    c     = 0.0

    do ifluid=1,nfluid
    do k=ks,ke
      do j=nb+1,nb+nyb
        do i=nb+1,nb+nxb

          vx=abs(u(imxa(ifluid),i,j,k)/u(idna(ifluid),i,j,k))
          vy=abs(u(imya(ifluid),i,j,k)/u(idna(ifluid),i,j,k))
          vz=abs(u(imza(ifluid),i,j,k)/u(idna(ifluid),i,j,k))
          v = max(vx,vy,vz)

          pmag = real(magn(ifluid))*sum(b(:,i,j,k)**2,1)/2.

          if(ifluid .gt. nadiab) then
            p = csi2*u(idna(ifluid),i,j,k)
            ps =p+pmag
            cf = sqrt(abs(  (2.*pmag+p)/u(idna(ifluid),i,j,k)) )
          else
#ifdef DUST
            if(idna(ifluid) .lt. idnd(1)) then
#endif /* DUST */
	      ps=(u(iena(ifluid),i,j,k)-sum(u(imxa(ifluid):imza(ifluid),i,j,k)**2,1) &
	         /u(idna(ifluid),i,j,k)/2.)*(gamma(ifluid)-1.)+(2.-gamma(ifluid))*pmag
              p=ps-pmag
              cf = sqrt(abs(  (2.*pmag+gamma(ifluid)*p)/u(idna(ifluid),i,j,k)) )
#ifdef DUST
	    else
	      cf = 0.0
	    endif
#endif /* DUST */
          endif

          cx=max(cx,vx+cf)
          cy=max(cy,vy+cf)
          cz=max(cz,vz+cf)
          c =max(c,cx,cy,cz)

        end do
      end do
    end do
    enddo

    dt_mhd_proc_x = dx/cx
    dt_mhd_proc_y = dy/cy
    dt_mhd_proc_z = dz/cz
    dt_mhd_proc   = min(dt_mhd_proc_x, dt_mhd_proc_y, dt_mhd_proc_z)

    call MPI_REDUCE(c, c_max_all, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)
    call MPI_BCAST(c_max_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    c = c_max_all

    call MPI_REDUCE(dt_mhd_proc, dt_mhd_all, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr)
    call MPI_BCAST(dt_mhd_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    dt_mhd = cfl*dt_mhd_all

  end subroutine timestep_mhd

end module time
