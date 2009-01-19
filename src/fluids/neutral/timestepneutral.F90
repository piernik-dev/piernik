
  subroutine timestep_neutral
    use grid, only   : dx,dy,dz
    use start, only  : nb,cfl, gamma,dt_neutal
    use fluidindex, only : nu, idnn,imxn,imyn,imzn
#ifndef ISO
    use fluidindex, only : ienn
#else /* ISO */
    use start, only : csi2
#endif /* ISO */

    use arrays, only : ks,ke,nyb,nxb,u,b

    implicit none

    real dt_neutral_proc, dt_neutral_all, c_max_all
    real dt_neutral_proc_x, dt_neutral_proc_y, dt_neutral_proc_z
    real cx, cy, cz, vx, vy, vz, cf

! locals
    real v,ps,p
    integer i,j,k

    cx    = 0.0
    cy    = 0.0
    cz    = 0.0
    c_neutral     = 0.0

    do k=ks,ke
      do j=nb+1,nb+nyb
        do i=nb+1,nb+nxb

          vx=abs(u(imxn,i,j,k)/u(idnn,i,j,k))
          vy=abs(u(imyn,i,j,k)/u(idnn,i,j,k))
          vz=abs(u(imzn,i,j,k)/u(idnn,i,j,k))

#ifdef ISO
            cf = csi
#else /* ISO */
	    p=(u(ienn,i,j,k)-sum(u(imxn:imzn,i,j,k)**2,1) &
	        /u(idnn,i,j,k)/2.)*(gamma(ineutral)-1.)
            cf = sqrt(abs( (gamma(ineutral)*p)/u(idnn,i,j,k)) )
#endif /* ISO */

          cx=max(cx,vx+cf)
          cy=max(cy,vy+cf)
          cz=max(cz,vz+cf)
          c_neutral =max(c_neutral,cx,cy,cz)

        end do
      end do
    end do

    dt_neutral_proc_x = dx/cx
    dt_neutral_proc_y = dy/cy
    dt_neutral_proc_z = dz/cz
    dt_neutral_proc   = min(dt_neutral_proc_x, dt_neutral_proc_y, dt_neutral_proc_z)

    call MPI_REDUCE(c_ionized, c_max_all, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)
    call MPI_BCAST(c_max_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    c_ionized = c_max_all

    call MPI_REDUCE(dt_neutral_proc, dt_neutral_all, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr)
    call MPI_BCAST(dt_neutral_all, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    dt_neutral = cfl*dt_neutral_all

  end subroutine timestep_neutral

!-------------------------------------------------------------------------------
