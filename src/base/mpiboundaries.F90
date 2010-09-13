! $Id$
!
! PIERNIK Code Copyright (C) 2006 Michal Hanasz
!
!    This file is part of PIERNIK code.
!
!    PIERNIK is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    PIERNIK is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with PIERNIK.  If not, see <http://www.gnu.org/licenses/>.
!
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
!    Ue-Li Pen
!        see: Pen, Arras & Wong (2003) for algorithm and
!             http://www.cita.utoronto.ca/~pen/MHD
!             for original source code "mhd.f90"
!
!    For full list of developers see $PIERNIK_HOME/license/pdt.txt
!
module mpiboundaries
!>
!! \brief (KK)
!<

contains
   subroutine mpi_boundaries_prep
      use mpisetup,   only : ierr, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
         MPI_YZ_LEFT_BND, MPI_YZ_RIGHT_BND, MPI_YZ_LEFT_DOM, MPI_YZ_RIGHT_DOM, &
         MPI_XZ_LEFT_BND, MPI_XZ_RIGHT_BND, MPI_XZ_LEFT_DOM, MPI_XZ_RIGHT_DOM, &
         MPI_XY_LEFT_BND, MPI_XY_RIGHT_BND, MPI_XY_LEFT_DOM, MPI_XY_RIGHT_DOM, &
         MAG_YZ_LEFT_BND, MAG_YZ_RIGHT_BND, MAG_YZ_LEFT_DOM, MAG_YZ_RIGHT_DOM, &
         MAG_XZ_LEFT_BND, MAG_XZ_RIGHT_BND, MAG_XZ_LEFT_DOM, MAG_XZ_RIGHT_DOM, &
         MAG_XY_LEFT_BND, MAG_XY_RIGHT_BND, MAG_XY_LEFT_DOM, MAG_XY_RIGHT_DOM, &
         ARR_YZ_LEFT_BND, ARR_YZ_RIGHT_BND, ARR_YZ_LEFT_DOM, ARR_YZ_RIGHT_DOM, &
         ARR_XZ_LEFT_BND, ARR_XZ_RIGHT_BND, ARR_XZ_LEFT_DOM, ARR_XZ_RIGHT_DOM, &
         ARR_XY_LEFT_BND, ARR_XY_RIGHT_BND, ARR_XY_LEFT_DOM, ARR_XY_RIGHT_DOM
      use fluidindex, only : nvar
      use arrays,     only : u
      use grid,       only : nb, nx, ny, nz, nxb, nyb, nzb, nxd, nyd, nzd

      implicit none
      integer, dimension(:), allocatable :: sizes, subsizes, starts
      integer(kind=4) :: ord
      integer(kind=4) :: old

      ord = MPI_ORDER_FORTRAN
      old = MPI_DOUBLE_PRECISION

!------------------------!
!   X dimension - fluid  !
!------------------------!
      if(nxd /= 1) then
         allocate(sizes(4), subsizes(4), starts(4))
         sizes    = (/nvar%all,nx,ny,nz/)
         subsizes = (/nvar%all,nb,ny,nz/)
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
         deallocate(sizes,subsizes,starts)

!--------------------------------!
!   X dimension - nx*ny*nz array !
!--------------------------------!
         allocate(sizes(3), subsizes(3), starts(3))
         sizes    = (/nx,ny,nz/)
         subsizes = (/nb,ny,nz/)
         starts   = (/0,0,0/)

         call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,ord,&
            old,ARR_YZ_LEFT_BND,ierr)
         call MPI_TYPE_COMMIT(ARR_YZ_LEFT_BND,ierr)

         starts(1) = nb
         call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,ord,&
            old,ARR_YZ_LEFT_DOM,ierr)
         call MPI_TYPE_COMMIT(ARR_YZ_LEFT_DOM,ierr)

         starts(1) = nxb
         call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,ord,&
            old,ARR_YZ_RIGHT_DOM,ierr)
         call MPI_TYPE_COMMIT(ARR_YZ_RIGHT_DOM,ierr)

         starts(1) = nxb+nb
         call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,ord,&
            old,ARR_YZ_RIGHT_BND,ierr)
         call MPI_TYPE_COMMIT(ARR_YZ_RIGHT_BND,ierr)
         deallocate(sizes,subsizes,starts)
      endif

!------------------------!
!   Y dimension - fluid  !
!------------------------!
      if(nyd /= 1) then
         allocate(sizes(4), subsizes(4), starts(4))
         sizes    = (/nvar%all,nx,ny,nz/)
         subsizes = (/nvar%all,nx,nb,nz/)
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
         deallocate(sizes,subsizes,starts)

!--------------------------------!
!   Y dimension - nx*ny*nz array !
!--------------------------------!
         allocate(sizes(3), subsizes(3), starts(3))
         sizes    = (/nx,ny,nz/)
         subsizes = (/nx,nb,nz/)
         starts   = (/0,0,0/)

         call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,ord,&
            old,ARR_XZ_LEFT_BND,ierr)
         call MPI_TYPE_COMMIT(ARR_XZ_LEFT_BND,ierr)

         starts(2) = nb
         call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,ord,&
            old,ARR_XZ_LEFT_DOM,ierr)
         call MPI_TYPE_COMMIT(ARR_XZ_LEFT_DOM,ierr)

         starts(2) = nyb
         call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,ord,&
            old,ARR_XZ_RIGHT_DOM,ierr)
         call MPI_TYPE_COMMIT(ARR_XZ_RIGHT_DOM,ierr)

         starts(2) = nyb+nb
         call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,ord,&
            old,ARR_XZ_RIGHT_BND,ierr)
         call MPI_TYPE_COMMIT(ARR_XZ_RIGHT_BND,ierr)
         deallocate(sizes,subsizes,starts)
      endif

!------------------------!
!   Z dimension - fluid  !
!------------------------!
      if(nzd /= 1) then
         allocate(sizes(4), subsizes(4), starts(4))
         sizes    = (/nvar%all,nx,ny,nz/)
         subsizes = (/nvar%all,nx,ny,nb/)
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
         deallocate(sizes,subsizes,starts)

!--------------------------------!
!   Z dimension - nx*ny*nz array !
!--------------------------------!
         allocate(sizes(3), subsizes(3), starts(3))
         sizes    = (/nx,ny,nz/)
         subsizes = (/nx,ny,nb/)
         starts   = (/0,0,0/)

         call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,ord,&
            old,ARR_XY_LEFT_BND,ierr)
         call MPI_TYPE_COMMIT(ARR_XY_LEFT_BND,ierr)

         starts(3) = nb
         call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,ord,&
            old,ARR_XY_LEFT_DOM,ierr)
         call MPI_TYPE_COMMIT(ARR_XY_LEFT_DOM,ierr)

         starts(3) = nzb
         call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,ord,&
            old,ARR_XY_RIGHT_DOM,ierr)
         call MPI_TYPE_COMMIT(ARR_XY_RIGHT_DOM,ierr)

         starts(3) = nzb+nb
         call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,ord,&
            old,ARR_XY_RIGHT_BND,ierr)
         call MPI_TYPE_COMMIT(ARR_XY_RIGHT_BND,ierr)
         deallocate(sizes,subsizes,starts)
      endif
   end subroutine mpi_boundaries_prep

end module mpiboundaries
