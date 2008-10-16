! $Id$
module mpi_bnd

contains
   subroutine mpi_bnd_prep
    use mpi_setup
    use arrays, only : nu,nx,ny,nz,u,nxb,nyb,nzb
    use start, only  : nb

    implicit none
    integer, dimension(4) :: sizes, subsizes, starts
    integer*4 :: ord
    integer*4 :: old

    ord = MPI_ORDER_FORTRAN
    old = MPI_DOUBLE_PRECISION

!------------------------!
!   X dimension - fluid  !
!------------------------!
    sizes    = (/nu,nx,ny,nz/)
    subsizes = (/nu,nb,ny,nz/)
    starts   = (/0,0,0,0/)
   
    call MPI_TYPE_CREATE_SUBARRAY(4,sizes,subsizes,starts,ord,&
           old,MPI_YZ_LEFT_BND,ierr)
    call MPI_TYPE_COMMIT(MPI_YZ_LEFT_BND,ierr)

    starts(2) = nb
    call MPI_TYPE_CREATE_SUBARRAY(4,sizes,subsizes,starts,ord,&
           old,MPI_YZ_LEFT_DOM,ierr)
    call MPI_TYPE_COMMIT(MPI_YZ_LEFT_DOM,ierr)

    starts(2) = nxb
    call MPI_TYPE_CREATE_SUBARRAY(4,sizes,subsizes,starts,ord,&
           old,MPI_YZ_RIGHT_DOM,ierr)
    call MPI_TYPE_COMMIT(MPI_YZ_RIGHT_DOM,ierr)

    starts(2) = nxb+nb
    call MPI_TYPE_CREATE_SUBARRAY(4,sizes,subsizes,starts,ord,&
           old,MPI_YZ_RIGHT_BND,ierr)
    call MPI_TYPE_COMMIT(MPI_YZ_RIGHT_BND,ierr)

!------------------------!
!   X dimension - Bfield !
!------------------------!
    sizes    = (/3,nx,ny,nz/)
    subsizes = (/3,nb,ny,nz/)
    starts   = (/0,0,0,0/)

    call MPI_TYPE_CREATE_SUBARRAY(4,sizes,subsizes,starts,ord,&
           old,MAG_YZ_LEFT_BND,ierr)
    call MPI_TYPE_COMMIT(MAG_YZ_LEFT_BND,ierr)

    starts(2) = nb
    call MPI_TYPE_CREATE_SUBARRAY(4,sizes,subsizes,starts,ord,&
           old,MAG_YZ_LEFT_DOM,ierr)
    call MPI_TYPE_COMMIT(MAG_YZ_LEFT_DOM,ierr)

    starts(2) = nxb
    call MPI_TYPE_CREATE_SUBARRAY(4,sizes,subsizes,starts,ord,&
           old,MAG_YZ_RIGHT_DOM,ierr)
    call MPI_TYPE_COMMIT(MAG_YZ_RIGHT_DOM,ierr)

    starts(2) = nxb+nb
    call MPI_TYPE_CREATE_SUBARRAY(4,sizes,subsizes,starts,ord,&
           old,MAG_YZ_RIGHT_BND,ierr)
    call MPI_TYPE_COMMIT(MAG_YZ_RIGHT_BND,ierr)

!------------------------!
!   Y dimension - fluid  !
!------------------------!
    sizes    = (/nu,nx,ny,nz/)
    subsizes = (/nu,nx,nb,nz/)
    starts   = (/0,0,0,0/)
   
    call MPI_TYPE_CREATE_SUBARRAY(4,sizes,subsizes,starts,ord,&
           old,MPI_XZ_LEFT_BND,ierr)
    call MPI_TYPE_COMMIT(MPI_XZ_LEFT_BND,ierr)

    starts(3) = nb
    call MPI_TYPE_CREATE_SUBARRAY(4,sizes,subsizes,starts,ord,&
           old,MPI_XZ_LEFT_DOM,ierr)
    call MPI_TYPE_COMMIT(MPI_XZ_LEFT_DOM,ierr)

    starts(3) = nyb
    call MPI_TYPE_CREATE_SUBARRAY(4,sizes,subsizes,starts,ord,&
           old,MPI_XZ_RIGHT_DOM,ierr)
    call MPI_TYPE_COMMIT(MPI_XZ_RIGHT_DOM,ierr)

    starts(3) = nyb+nb
    call MPI_TYPE_CREATE_SUBARRAY(4,sizes,subsizes,starts,ord,&
           old,MPI_XZ_RIGHT_BND,ierr)
    call MPI_TYPE_COMMIT(MPI_XZ_RIGHT_BND,ierr)

!------------------------!
!   Y dimension - Bfield !
!------------------------!
    sizes    = (/3,nx,ny,nz/)
    subsizes = (/3,nx,nb,nz/)
    starts   = (/0,0,0,0/)
   
    call MPI_TYPE_CREATE_SUBARRAY(4,sizes,subsizes,starts,ord,&
           old,MAG_XZ_LEFT_BND,ierr)
    call MPI_TYPE_COMMIT(MAG_XZ_LEFT_BND,ierr)

    starts(3) = nb
    call MPI_TYPE_CREATE_SUBARRAY(4,sizes,subsizes,starts,ord,&
           old,MAG_XZ_LEFT_DOM,ierr)
    call MPI_TYPE_COMMIT(MAG_XZ_LEFT_DOM,ierr)

    starts(3) = nyb
    call MPI_TYPE_CREATE_SUBARRAY(4,sizes,subsizes,starts,ord,&
           old,MAG_XZ_RIGHT_DOM,ierr)
    call MPI_TYPE_COMMIT(MAG_XZ_RIGHT_DOM,ierr)

    starts(3) = nyb+nb
    call MPI_TYPE_CREATE_SUBARRAY(4,sizes,subsizes,starts,ord,&
           old,MAG_XZ_RIGHT_BND,ierr)
    call MPI_TYPE_COMMIT(MAG_XZ_RIGHT_BND,ierr)

!------------------------!
!   Z dimension - fluid  !
!------------------------!
    sizes    = (/nu,nx,ny,nz/)
    subsizes = (/nu,nx,ny,nb/)
    starts   = (/0,0,0,0/)
   
    call MPI_TYPE_CREATE_SUBARRAY(4,sizes,subsizes,starts,ord,&
           old,MPI_XY_LEFT_BND,ierr)
    call MPI_TYPE_COMMIT(MPI_XY_LEFT_BND,ierr)

    starts(4) = nb
    call MPI_TYPE_CREATE_SUBARRAY(4,sizes,subsizes,starts,ord,&
           old,MPI_XY_LEFT_DOM,ierr)
    call MPI_TYPE_COMMIT(MPI_XY_LEFT_DOM,ierr)

    starts(4) = nzb
    call MPI_TYPE_CREATE_SUBARRAY(4,sizes,subsizes,starts,ord,&
           old,MPI_XY_RIGHT_DOM,ierr)
    call MPI_TYPE_COMMIT(MPI_XY_RIGHT_DOM,ierr)

    starts(4) = nzb+nb
    call MPI_TYPE_CREATE_SUBARRAY(4,sizes,subsizes,starts,ord,&
           old,MPI_XY_RIGHT_BND,ierr)
    call MPI_TYPE_COMMIT(MPI_XY_RIGHT_BND,ierr)

!------------------------!
!   Z dimension - Bfield !
!------------------------!
    sizes    = (/3,nx,ny,nz/)
    subsizes = (/3,nx,ny,nb/)
    starts   = (/0,0,0,0/)
   
    call MPI_TYPE_CREATE_SUBARRAY(4,sizes,subsizes,starts,ord,&
           old,MAG_XY_LEFT_BND,ierr)
    call MPI_TYPE_COMMIT(MAG_XY_LEFT_BND,ierr)

    starts(4) = nb
    call MPI_TYPE_CREATE_SUBARRAY(4,sizes,subsizes,starts,ord,&
           old,MAG_XY_LEFT_DOM,ierr)
    call MPI_TYPE_COMMIT(MAG_XY_LEFT_DOM,ierr)

    starts(4) = nzb
    call MPI_TYPE_CREATE_SUBARRAY(4,sizes,subsizes,starts,ord,&
           old,MAG_XY_RIGHT_DOM,ierr)
    call MPI_TYPE_COMMIT(MAG_XY_RIGHT_DOM,ierr)

    starts(4) = nzb+nb
    call MPI_TYPE_CREATE_SUBARRAY(4,sizes,subsizes,starts,ord,&
           old,MAG_XY_RIGHT_BND,ierr)
    call MPI_TYPE_COMMIT(MAG_XY_RIGHT_BND,ierr)
    

  end subroutine mpi_bnd_prep
end module mpi_bnd
