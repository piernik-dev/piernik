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
!    Initial implementation of PIERNIK code was based on TVD split MHD code by
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
   implicit none
   private
   public :: mpi_boundaries_prep
contains
   subroutine mpi_boundaries_prep
      use arrays,     only: u
      use fluidindex, only: nvar
      use grid,       only: cg
      use mpi,        only: MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION
      use mpisetup,   only: ierr, xdim, ydim, zdim, has_dir, &
                            MPI_YZ_LEFT_BND, MPI_YZ_RIGHT_BND, MPI_YZ_LEFT_DOM, MPI_YZ_RIGHT_DOM, &
                            MPI_XZ_LEFT_BND, MPI_XZ_RIGHT_BND, MPI_XZ_LEFT_DOM, MPI_XZ_RIGHT_DOM, &
                            MPI_XY_LEFT_BND, MPI_XY_RIGHT_BND, MPI_XY_LEFT_DOM, MPI_XY_RIGHT_DOM, &
                            MAG_YZ_LEFT_BND, MAG_YZ_RIGHT_BND, MAG_YZ_LEFT_DOM, MAG_YZ_RIGHT_DOM, &
                            MAG_XZ_LEFT_BND, MAG_XZ_RIGHT_BND, MAG_XZ_LEFT_DOM, MAG_XZ_RIGHT_DOM, &
                            MAG_XY_LEFT_BND, MAG_XY_RIGHT_BND, MAG_XY_LEFT_DOM, MAG_XY_RIGHT_DOM, &
                            ARR_YZ_LEFT_BND, ARR_YZ_RIGHT_BND, ARR_YZ_LEFT_DOM, ARR_YZ_RIGHT_DOM, &
                            ARR_XZ_LEFT_BND, ARR_XZ_RIGHT_BND, ARR_XZ_LEFT_DOM, ARR_XZ_RIGHT_DOM, &
                            ARR_XY_LEFT_BND, ARR_XY_RIGHT_BND, ARR_XY_LEFT_DOM, ARR_XY_RIGHT_DOM

      implicit none
      integer, dimension(:), allocatable :: sizes, subsizes, starts
      integer(kind=4) :: ord
      integer(kind=4) :: old

      ord = MPI_ORDER_FORTRAN
      old = MPI_DOUBLE_PRECISION

!------------------------!
!   X dimension - fluid  !
!------------------------!
      if (has_dir(xdim)) then

         allocate(sizes(4), subsizes(4), starts(4))

         sizes    = [ nvar%all, cg%nx, cg%ny, cg%nz ]
         subsizes = [ nvar%all, cg%nb, cg%ny, cg%nz ]
         starts   = [ 0, 0, 0, 0 ]

         call MPI_Type_create_subarray(4,sizes,subsizes,starts,ord,&
            old,MPI_YZ_LEFT_BND,ierr)
         call MPI_Type_commit(MPI_YZ_LEFT_BND,ierr)

         starts(2) = cg%nb
         call MPI_Type_create_subarray(4,sizes,subsizes,starts,ord,&
            old,MPI_YZ_LEFT_DOM,ierr)
         call MPI_Type_commit(MPI_YZ_LEFT_DOM,ierr)

         starts(2) = cg%nxb
         call MPI_Type_create_subarray(4,sizes,subsizes,starts,ord,&
            old,MPI_YZ_RIGHT_DOM,ierr)
         call MPI_Type_commit(MPI_YZ_RIGHT_DOM,ierr)

         starts(2) = cg%nxb+cg%nb
         call MPI_Type_create_subarray(4,sizes,subsizes,starts,ord,&
            old,MPI_YZ_RIGHT_BND,ierr)
         call MPI_Type_commit(MPI_YZ_RIGHT_BND,ierr)

!------------------------!
!   X dimension - Bfield !
!------------------------!
         sizes    = [ 3, cg%nx, cg%ny, cg%nz ]
         subsizes = [ 3, cg%nb, cg%ny, cg%nz ]
         starts   = [ 0, 0, 0, 0 ]

         call MPI_Type_create_subarray(4,sizes,subsizes,starts,ord,&
            old,MAG_YZ_LEFT_BND,ierr)
         call MPI_Type_commit(MAG_YZ_LEFT_BND,ierr)

         starts(2) = cg%nb
         call MPI_Type_create_subarray(4,sizes,subsizes,starts,ord,&
            old,MAG_YZ_LEFT_DOM,ierr)
         call MPI_Type_commit(MAG_YZ_LEFT_DOM,ierr)

         starts(2) = cg%nxb
         call MPI_Type_create_subarray(4,sizes,subsizes,starts,ord,&
            old,MAG_YZ_RIGHT_DOM,ierr)
         call MPI_Type_commit(MAG_YZ_RIGHT_DOM,ierr)

         starts(2) = cg%nxb+cg%nb
         call MPI_Type_create_subarray(4,sizes,subsizes,starts,ord,&
            old,MAG_YZ_RIGHT_BND,ierr)
         call MPI_Type_commit(MAG_YZ_RIGHT_BND,ierr)

         deallocate(sizes,subsizes,starts)

!---------------------------------------!
!   X dimension - cg%nx*cg%ny*cg%nz array (grav) !
!---------------------------------------!
         allocate(sizes(3), subsizes(3), starts(3))

         sizes    = [ cg%nx, cg%ny, cg%nz ]
         subsizes = [ cg%nb, cg%ny, cg%nz ]
         starts   = [ 0, 0, 0 ]

         call MPI_Type_create_subarray(3, sizes, subsizes, starts, ord, old, ARR_YZ_LEFT_BND,  ierr)
         call MPI_Type_commit(ARR_YZ_LEFT_BND,  ierr)

         starts(1) = cg%nb
         call MPI_Type_create_subarray(3, sizes, subsizes, starts, ord, old, ARR_YZ_LEFT_DOM,  ierr)
         call MPI_Type_commit(ARR_YZ_LEFT_DOM,  ierr)

         starts(1) = cg%nxb
         call MPI_Type_create_subarray(3, sizes, subsizes, starts, ord, old, ARR_YZ_RIGHT_DOM, ierr)
         call MPI_Type_commit(ARR_YZ_RIGHT_DOM, ierr)

         starts(1) = cg%nxb+cg%nb
         call MPI_Type_create_subarray(3, sizes, subsizes, starts, ord, old, ARR_YZ_RIGHT_BND, ierr)
         call MPI_Type_commit(ARR_YZ_RIGHT_BND, ierr)

         deallocate(sizes,subsizes,starts)

      endif

!------------------------!
!   Y dimension - fluid  !
!------------------------!
      if (has_dir(ydim)) then

         allocate(sizes(4), subsizes(4), starts(4))

         sizes    = [ nvar%all, cg%nx, cg%ny, cg%nz ]
         subsizes = [ nvar%all, cg%nx, cg%nb, cg%nz ]
         starts   = [ 0, 0, 0, 0 ]

         call MPI_Type_create_subarray(4,sizes,subsizes,starts,ord,&
            old,MPI_XZ_LEFT_BND,ierr)
         call MPI_Type_commit(MPI_XZ_LEFT_BND,ierr)

         starts(3) = cg%nb
         call MPI_Type_create_subarray(4,sizes,subsizes,starts,ord,&
            old,MPI_XZ_LEFT_DOM,ierr)
         call MPI_Type_commit(MPI_XZ_LEFT_DOM,ierr)

         starts(3) = cg%nyb
         call MPI_Type_create_subarray(4,sizes,subsizes,starts,ord,&
            old,MPI_XZ_RIGHT_DOM,ierr)
         call MPI_Type_commit(MPI_XZ_RIGHT_DOM,ierr)

         starts(3) = cg%nyb+cg%nb
         call MPI_Type_create_subarray(4,sizes,subsizes,starts,ord,&
            old,MPI_XZ_RIGHT_BND,ierr)
         call MPI_Type_commit(MPI_XZ_RIGHT_BND,ierr)

!------------------------!
!   Y dimension - Bfield !
!------------------------!
         sizes    = [ 3, cg%nx, cg%ny, cg%nz ]
         subsizes = [ 3, cg%nx, cg%nb, cg%nz ]
         starts   = [ 0, 0, 0, 0 ]

         call MPI_Type_create_subarray(4,sizes,subsizes,starts,ord,&
            old,MAG_XZ_LEFT_BND,ierr)
         call MPI_Type_commit(MAG_XZ_LEFT_BND,ierr)

         starts(3) = cg%nb
         call MPI_Type_create_subarray(4,sizes,subsizes,starts,ord,&
            old,MAG_XZ_LEFT_DOM,ierr)
         call MPI_Type_commit(MAG_XZ_LEFT_DOM,ierr)

         starts(3) = cg%nyb
         call MPI_Type_create_subarray(4,sizes,subsizes,starts,ord,&
            old,MAG_XZ_RIGHT_DOM,ierr)
         call MPI_Type_commit(MAG_XZ_RIGHT_DOM,ierr)

         starts(3) = cg%nyb+cg%nb
         call MPI_Type_create_subarray(4,sizes,subsizes,starts,ord,&
            old,MAG_XZ_RIGHT_BND,ierr)
         call MPI_Type_commit(MAG_XZ_RIGHT_BND,ierr)

         deallocate(sizes,subsizes,starts)

!---------------------------------------!
!   Y dimension - cg%nx*cg%ny*cg%nz array (grav) !
!---------------------------------------!
         allocate(sizes(3), subsizes(3), starts(3))

         sizes    = [ cg%nx, cg%ny, cg%nz ]
         subsizes = [ cg%nx, cg%nb, cg%nz ]
         starts   = [ 0, 0, 0 ]

         call MPI_Type_create_subarray(3, sizes, subsizes, starts, ord, old, ARR_XZ_LEFT_BND,  ierr)
         call MPI_Type_commit(ARR_XZ_LEFT_BND,  ierr)

         starts(2) = cg%nb
         call MPI_Type_create_subarray(3, sizes, subsizes, starts, ord, old, ARR_XZ_LEFT_DOM,  ierr)
         call MPI_Type_commit(ARR_XZ_LEFT_DOM,  ierr)

         starts(2) = cg%nyb
         call MPI_Type_create_subarray(3, sizes, subsizes, starts, ord, old, ARR_XZ_RIGHT_DOM, ierr)
         call MPI_Type_commit(ARR_XZ_RIGHT_DOM, ierr)

         starts(2) = cg%nyb+cg%nb
         call MPI_Type_create_subarray(3, sizes, subsizes, starts, ord, old, ARR_XZ_RIGHT_BND, ierr)
         call MPI_Type_commit(ARR_XZ_RIGHT_BND, ierr)

         deallocate(sizes,subsizes,starts)

      endif

!------------------------!
!   Z dimension - fluid  !
!------------------------!
      if (has_dir(zdim)) then

         allocate(sizes(4), subsizes(4), starts(4))

         sizes    = [ nvar%all, cg%nx, cg%ny, cg%nz ]
         subsizes = [ nvar%all, cg%nx, cg%ny, cg%nb ]
         starts   = [ 0, 0, 0, 0 ]

         call MPI_Type_create_subarray(4,sizes,subsizes,starts,ord,&
            old,MPI_XY_LEFT_BND,ierr)
         call MPI_Type_commit(MPI_XY_LEFT_BND,ierr)

         starts(4) = cg%nb
         call MPI_Type_create_subarray(4,sizes,subsizes,starts,ord,&
            old,MPI_XY_LEFT_DOM,ierr)
         call MPI_Type_commit(MPI_XY_LEFT_DOM,ierr)

         starts(4) = cg%nzb
         call MPI_Type_create_subarray(4,sizes,subsizes,starts,ord,&
            old,MPI_XY_RIGHT_DOM,ierr)
         call MPI_Type_commit(MPI_XY_RIGHT_DOM,ierr)

         starts(4) = cg%nzb+cg%nb
         call MPI_Type_create_subarray(4,sizes,subsizes,starts,ord,&
            old,MPI_XY_RIGHT_BND,ierr)
         call MPI_Type_commit(MPI_XY_RIGHT_BND,ierr)

!------------------------!
!   Z dimension - Bfield !
!------------------------!
         sizes    = [ 3, cg%nx, cg%ny, cg%nz ]
         subsizes = [ 3, cg%nx, cg%ny, cg%nb ]
         starts   = [ 0, 0, 0, 0 ]

         call MPI_Type_create_subarray(4,sizes,subsizes,starts,ord,&
            old,MAG_XY_LEFT_BND,ierr)
         call MPI_Type_commit(MAG_XY_LEFT_BND,ierr)

         starts(4) = cg%nb
         call MPI_Type_create_subarray(4,sizes,subsizes,starts,ord,&
            old,MAG_XY_LEFT_DOM,ierr)
         call MPI_Type_commit(MAG_XY_LEFT_DOM,ierr)

         starts(4) = cg%nzb
         call MPI_Type_create_subarray(4,sizes,subsizes,starts,ord,&
            old,MAG_XY_RIGHT_DOM,ierr)
         call MPI_Type_commit(MAG_XY_RIGHT_DOM,ierr)

         starts(4) = cg%nzb+cg%nb
         call MPI_Type_create_subarray(4,sizes,subsizes,starts,ord,&
            old,MAG_XY_RIGHT_BND,ierr)
         call MPI_Type_commit(MAG_XY_RIGHT_BND,ierr)

         deallocate(sizes,subsizes,starts)

!---------------------------------------!
!   Z dimension - cg%nx*cg%ny*cg%nz array (grav) !
!---------------------------------------!
         allocate(sizes(3), subsizes(3), starts(3))

         sizes    = [ cg%nx, cg%ny, cg%nz ]
         subsizes = [ cg%nx, cg%ny, cg%nb ]
         starts   = [ 0, 0, 0 ]

         call MPI_Type_create_subarray(3, sizes, subsizes, starts, ord, old, ARR_XY_LEFT_BND,  ierr)
         call MPI_Type_commit(ARR_XY_LEFT_BND,  ierr)

         starts(3) = cg%nb
         call MPI_Type_create_subarray(3, sizes, subsizes, starts, ord, old, ARR_XY_LEFT_DOM,  ierr)
         call MPI_Type_commit(ARR_XY_LEFT_DOM,  ierr)

         starts(3) = cg%nzb
         call MPI_Type_create_subarray(3, sizes, subsizes, starts, ord, old, ARR_XY_RIGHT_DOM, ierr)
         call MPI_Type_commit(ARR_XY_RIGHT_DOM, ierr)

         starts(3) = cg%nzb+cg%nb
         call MPI_Type_create_subarray(3, sizes, subsizes, starts, ord, old, ARR_XY_RIGHT_BND, ierr)
         call MPI_Type_commit(ARR_XY_RIGHT_BND, ierr)

         deallocate(sizes,subsizes,starts)

      endif

   end subroutine mpi_boundaries_prep

end module mpiboundaries
