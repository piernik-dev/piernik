
  subroutine timestep_dust
    use grid, only   : dx,dy,dz,nb
    use start, only  : cfl, dt_dust
    use fluidindex, only : idnd,imxd,imyd,imzd

    use arrays, only : ks,ke,nyb,nxb,u,b

    implicit none

    real dt_dust_proc, dt_dust_all, c_max_all
    real dt_dust_proc_x, dt_dust_proc_y, dt_dust_proc_z
    real cx, cy, cz, vx, vy, vz, cf

! locals
    real v,ps,p
    integer i,j,k

    cx    = 0.0
    cy    = 0.0
    cz    = 0.0
    c_dust     = 0.0

    do k=ks,ke
      do j=nb+1,nb+nyb
        do i=nb+1,nb+nxb

          vx=abs(u(imxd,i,j,k)/u(idnd,i,j,k))
          vy=abs(u(imyd,i,j,k)/u(idnd,i,j,k))
          vz=abs(u(imzd,i,j,k)/u(idnd,i,j,k))

          cx=max(cx,vx)
          cy=max(cy,vy)
          cz=max(cz,vz)
          c_dust =max(c_dust,cx,cy,cz)

        end do
      end do
    end do

    dt_dust_proc_x = dx/cx
    dt_dust_proc_y = dy/cy
    dt_dust_proc_z = dz/cz
    dt_dust_proc   = min(dt_dust_proc_x, dt_dust_proc_y, dt_dust_proc_z)

    call MPI_REDUCE(c_dust, c_max_all, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)
    call MPI_BCAST(c_max_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    c_dust = c_max_all

    call MPI_REDUCE(dt_dust_proc, dt_dust_all, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr)
    call MPI_BCAST(dt_dust_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    dt_dust = cfl*dt_dust_all

  end subroutine timestep_dust

!-------------------------------------------------------------------------------
