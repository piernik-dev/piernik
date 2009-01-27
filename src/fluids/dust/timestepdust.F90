! $Id$
#include "piernik.def"

module timestepdust

  real :: dt_dst,c_dst

contains

  subroutine timestep_dst
    use mpi_setup
    use grid, only      : dx,dy,dz,nb,ks,ke,is,ie,js,je
    use arrays, only    : u,b
    use start, only     : cfl  
    use initdust, only  : idnd,imxd,imyd,imzd
    use constants, only : small

    implicit none

    real dt_dst_proc, dt_dst_all, c_max_all
    real dt_dst_proc_x, dt_dst_proc_y, dt_dst_proc_z
    real cx, cy, cz, vx, vy, vz
    

! locals

    real p
    integer i,j,k


    cx    = 0.0
    cy    = 0.0
    cz    = 0.0
    c_dst     = 0.0

    do k=ks,ke
      do j=js,je
        do i=is,ie

          vx=abs(u(imxd,i,j,k)/u(idnd,i,j,k))
          vy=abs(u(imyd,i,j,k)/u(idnd,i,j,k))
          vz=abs(u(imzd,i,j,k)/u(idnd,i,j,k))

          cx=max(cx,vx)
          cy=max(cy,vy)
          cz=max(cz,vz)
          c_dst =max(c_dst,cx,cy,cz)

        end do
      end do
    end do

    dt_dst_proc_x = dx/max(cx,small)
    dt_dst_proc_y = dy/max(cy,small)
    dt_dst_proc_z = dz/max(cz,small)
    dt_dst_proc   = min(dt_dst_proc_x, dt_dst_proc_y, dt_dst_proc_z)

    call MPI_REDUCE(c_dst, c_max_all, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)
    call MPI_BCAST(c_max_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    c_dst = c_max_all

    call MPI_REDUCE(dt_dst_proc, dt_dst_all, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr)
    call MPI_BCAST(dt_dst_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    dt_dst = cfl*dt_dst_all
    
!    write(*,*) 'timestep_dst:', dt_dst

  end subroutine timestep_dst

!-------------------------------------------------------------------------------
end module timestepdust

