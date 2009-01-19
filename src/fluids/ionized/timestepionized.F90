! $Id$
#include "piernik.def"

module timestepionized

  real :: dt_ion,c_ion

contains

  subroutine timestep_ion
    use mpi_setup
    use grid, only   : dx,dy,dz
    use start, only  : nb,cfl
!    use timestep, only : dt_ion
    use initionized, only : gamma_ion, cs_iso_ion,cs_iso_ion2   
    use initionized, only : idni,imxi,imyi,imzi
#ifndef ISO
    use initionized, only : ieni
#endif /* ISO */

    use arrays, only : ks,ke,nyb,nxb,u,b

    implicit none

    real dt_ion_proc, dt_ion_all, c_max_all
    real dt_ion_proc_x, dt_ion_proc_y, dt_ion_proc_z
    real cx, cy, cz, vx, vy, vz, cf
    

! locals
    real pmag
    real v,ps,p
    integer i,j,k
    real csi_sq


    cx    = 0.0
    cy    = 0.0
    cz    = 0.0
    c_ion     = 0.0

    do k=ks,ke
      do j=nb+1,nb+nyb
        do i=nb+1,nb+nxb

          vx=abs(u(imxi,i,j,k)/u(idni,i,j,k))
          vy=abs(u(imyi,i,j,k)/u(idni,i,j,k))
          vz=abs(u(imzi,i,j,k)/u(idni,i,j,k))

          pmag = sum(b(:,i,j,k)**2,1)/2.

#ifdef ISO
            p = csi_sq*u(idni,i,j,k)
            ps =p+pmag
            cf = sqrt(abs(  (2.*pmag+p)/u(idni,i,j,k)) )
#else /* ISO */
	    ps=(u(ieni,i,j,k)-sum(u(imxi:imzi,i,j,k)**2,1) &
	        /u(idni,i,j,k)/2.)*(gamma_ion-1.)+(2.-gamma_ion)*pmag
            p=ps-pmag
            cf = sqrt(abs(  (2.*pmag+gamma_ion*p)/u(idni,i,j,k)) )
#endif /* ISO */

          cx=max(cx,vx+cf)
          cy=max(cy,vy+cf)
          cz=max(cz,vz+cf)
          c_ion =max(c_ion,cx,cy,cz)

        end do
      end do
    end do

    dt_ion_proc_x = dx/cx
    dt_ion_proc_y = dy/cy
    dt_ion_proc_z = dz/cz
    dt_ion_proc   = min(dt_ion_proc_x, dt_ion_proc_y, dt_ion_proc_z)

    call MPI_REDUCE(c_ion, c_max_all, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)
    call MPI_BCAST(c_max_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    c_ion = c_max_all

    call MPI_REDUCE(dt_ion_proc, dt_ion_all, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr)
    call MPI_BCAST(dt_ion_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    dt_ion = cfl*dt_ion_all
    
!    write(*,*) 'timestep_ion:', dt_ion

  end subroutine timestep_ion

!-------------------------------------------------------------------------------
end module timestepionized

