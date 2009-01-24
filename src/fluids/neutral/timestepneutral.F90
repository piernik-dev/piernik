! $Id$
#include "piernik.def"

module timestepneutral

  real :: dt_neu,c_neu

contains

  subroutine timestep_neu
    use mpi_setup
    use grid, only     : dx,dy,dz,nb,ks,ke,nyb,nxb
    use arrays, only   : u,b
    use start, only  : cfl
    use initneutral, only : gamma_neu, cs_iso_neu,cs_iso_neu2   
    use initneutral, only : idnn,imxn,imyn,imzn
#ifndef ISO
    use initneutral, only : ienn
#endif /* ISO */


    implicit none

    real dt_neu_proc, dt_neu_all, c_max_all
    real dt_neu_proc_x, dt_neu_proc_y, dt_neu_proc_z
    real cx, cy, cz, vx, vy, vz, cs
    

! locals

    real p
    integer i,j,k


    cx    = 0.0
    cy    = 0.0
    cz    = 0.0
    c_neu = 0.0

    do k=ks,ke
      do j=nb+1,nb+nyb
        do i=nb+1,nb+nxb

          vx=abs(u(imxn,i,j,k)/u(idnn,i,j,k))
          vy=abs(u(imyn,i,j,k)/u(idnn,i,j,k))
          vz=abs(u(imzn,i,j,k)/u(idnn,i,j,k))

#ifdef ISO
            p = cs_iso_neu2*u(idnn,i,j,k)
            cs = sqrt(abs( p/u(idnn,i,j,k)) )
#else /* ISO */
            p=(u(ienn,i,j,k)-sum(u(imxn:imzn,i,j,k)**2,1) &
              /u(idnn,i,j,k)/2.)*(gamma_neu-1.)

            cs = sqrt(abs(  (gamma_neu*p)/u(idnn,i,j,k)) )
#endif /* ISO */

          cx=max(cx,vx+cs)
          cy=max(cy,vy+cs)
          cz=max(cz,vz+cs)
          c_neu =max(c_neu,cx,cy,cz)

        end do
      end do
    end do

    dt_neu_proc_x = dx/cx
    dt_neu_proc_y = dy/cy
    dt_neu_proc_z = dz/cz
    dt_neu_proc   = min(dt_neu_proc_x, dt_neu_proc_y, dt_neu_proc_z)

    call MPI_REDUCE(c_neu, c_max_all, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)
    call MPI_BCAST(c_max_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    c_neu = c_max_all

    call MPI_REDUCE(dt_neu_proc, dt_neu_all, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr)
    call MPI_BCAST(dt_neu_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    dt_neu = cfl*dt_neu_all
    
!    write(*,*) 'timestep_neu:', dt_neu

  end subroutine timestep_neu

!-------------------------------------------------------------------------------
end module timestepneutral

