#include "mhd.def"

module time

! Based on: Pen Arras & Wong (2003)
! Modified extensively by M. Hanasz November 2005 - May 2006
! Optimalization of cpu-time efficiency of mhdflux, tvd1 by K. Kowalik
! Modification history:  see "changelog"  -- warning: out of date

  use mpi_setup
  use start
  use arrays
  use grid
  use init_problem 
!  use diagnostics
  use thermal

#ifdef COSM_RAYS
  use cr_diffusion
#endif COSM_RAYS  
#ifdef RESIST
  use resistivity
#endif
  
  implicit none
  real c

contains


  subroutine timestep
  
! Timestep computation

    call timestep_mhd
    dt=min(dt_mhd,(tend-t)/2.)

#ifdef RESIST
    call timestep_resist
    dt = min(dt,dt_resist)
#endif RESIST
#ifdef COOL_HEAT
    call timestep_coolheat
    dt = min(dt,dt_coolheat)
#endif COOL_HEAT
#ifdef HEAT_COND
    dt_heatcond = cfl_heatcond * 0.5*dxmn**2/(K_heatcond+small)
    dt = min(dt,dt_heatcond)
#endif HEAT_COND
#ifdef VISC
    dt_visc = cfl_visc * 0.5*dxmn**2/(nu_bulk+small)
    dt = min(dt,dt_visc)
#endif VISC   
#ifdef COSM_RAYS
    dt_cr = cfl_cr * 0.5*dxmn**2/(K_cr_paral+K_cr_perp+small)
    dt = min(dt,dt_cr)
#endif COSM_RAYS 

  end subroutine timestep

!------------------------------------------------------------------------------------------

  subroutine timestep_mhd

    real dt_mhd_proc, dt_mhd_all, c_max_all 
    real dt_mhd_proc_x, dt_mhd_proc_y, dt_mhd_proc_z
    real cx, cy, cz, vx, vy, vz, cf

! locals
    real pmag
    real v,ps,p,bx,by,bz
    integer i,j,k,ip,jp,kp

    cx    = 0.0
    cy    = 0.0
    cz    = 0.0
    c     = 0.0

    do k=ks,ke
      do j=nb+1,nb+nyb
        do i=nb+1,nb+nxb

          vx=abs(u(imxa,i,j,k)/u(idna,i,j,k))
          vy=abs(u(imya,i,j,k)/u(idna,i,j,k))
          vz=abs(u(imza,i,j,k)/u(idna,i,j,k))
          v = max(vx,vy,vz)
          
          pmag = sum(b(:,i,j,k)**2,1)/2.

#ifdef ISO
          p = csi2*u(idna,i,j,k)
          ps =p+pmag
          cf = sqrt(abs(  (2.*pmag+p)/u(idna,i,j,k)) )
#else ISO
          ps=(u(iena,i,j,k)-sum(u(imxa:imza,i,j,k)**2,1)/u(idna,i,j,k)/2.)*(gamma-1.)+(2.-gamma)*pmag
          p=ps-pmag          
          cf = sqrt(abs(  (2.*pmag+gamma*p)/u(idna,i,j,k)) )
#endif ISO

          cx=max(cx,vx+cf)
          cy=max(cy,vy+cf)
          cz=max(cz,vz+cf)
          c =max(c,cx,cy,cz)

        end do
      end do
    end do

    dt_mhd_proc_x = cfl*dx/cx
    dt_mhd_proc_y = cfl*dy/cy
    dt_mhd_proc_z = cfl*dz/cz
    dt_mhd_proc   = min(dt_mhd_proc_x, dt_mhd_proc_y, dt_mhd_proc_z)

    call MPI_REDUCE(c, c_max_all, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)
    call MPI_BCAST(c_max_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    
    c = c_max_all
    
    call MPI_REDUCE(dt_mhd_proc, dt_mhd_all, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr)
    call MPI_BCAST(dt_mhd_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    dt_mhd = dt_mhd_all

  end subroutine timestep_mhd

end module time
