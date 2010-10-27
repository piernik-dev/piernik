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
#include "piernik.def"

module magboundaries

contains

   subroutine bnd_a(A)
      use mpisetup, only: MAG_YZ_LEFT_DOM, MAG_YZ_RIGHT_DOM, MAG_YZ_LEFT_BND, MAG_YZ_RIGHT_BND, &
         MAG_XY_LEFT_DOM, MAG_XY_RIGHT_DOM, MAG_XY_LEFT_BND, MAG_XY_RIGHT_BND, &
         MAG_XZ_LEFT_DOM, MAG_XZ_RIGHT_DOM, MAG_XZ_LEFT_BND, MAG_XZ_RIGHT_BND, &
         ierr, req, comm3d, procxl, procxr, procyl, procyr, proczl, proczr, status, &
         pxsize, pysize, pzsize
      implicit none
      real, dimension(:,:,:,:) :: A

      if (pxsize .gt. 1) then

         CALL MPI_Isend  (A(1,1,1,1), 1, MAG_YZ_LEFT_DOM,  procxl, 10, comm3d, req(1), ierr)
         CALL MPI_Isend  (A(1,1,1,1), 1, MAG_YZ_RIGHT_DOM, procxr, 20, comm3d, req(3), ierr)
         CALL MPI_Irecv  (A(1,1,1,1), 1, MAG_YZ_LEFT_BND,  procxl, 20, comm3d, req(2), ierr)
         CALL MPI_Irecv  (A(1,1,1,1), 1, MAG_YZ_RIGHT_BND, procxr, 10, comm3d, req(4), ierr)

         call MPI_Waitall(4,req(:),status(:,:),ierr)
      endif

      if (pysize .gt. 1) then
         CALL MPI_Isend  (A(1,1,1,1), 1, MAG_XZ_LEFT_DOM,  procyl, 30, comm3d, req(1), ierr)
         CALL MPI_Isend  (A(1,1,1,1), 1, MAG_XZ_RIGHT_DOM, procyr, 40, comm3d, req(3), ierr)
         CALL MPI_Irecv  (A(1,1,1,1), 1, MAG_XZ_LEFT_BND,  procyl, 40, comm3d, req(2), ierr)
         CALL MPI_Irecv  (A(1,1,1,1), 1, MAG_XZ_RIGHT_BND, procyr, 30, comm3d, req(4), ierr)

         call MPI_Waitall(4,req(:),status(:,:),ierr)
      endif

      if (pzsize .gt. 1) then
         CALL MPI_Isend  (A(1,1,1,1), 1, MAG_XY_LEFT_DOM,  proczl, 50, comm3d, req(1), ierr)
         CALL MPI_Isend  (A(1,1,1,1), 1, MAG_XY_RIGHT_DOM, proczr, 60, comm3d, req(3), ierr)
         CALL MPI_Irecv  (A(1,1,1,1), 1, MAG_XY_LEFT_BND,  proczl, 60, comm3d, req(2), ierr)
         CALL MPI_Irecv  (A(1,1,1,1), 1, MAG_XY_RIGHT_BND, proczr, 50, comm3d, req(4), ierr)

         call MPI_Waitall(4,req(:),status(:,:),ierr)
      endif
   end subroutine bnd_a

   subroutine bnd_b(dim)
      use arrays,        only: b
      use dataio_pub, only: msg, warn
      use fluidindex,    only: ibx, iby, ibz
      use grid,          only: nb, nx, ny, nz, nxb, nyb, nzb
      use mpisetup,      only: MPI_DOUBLE_PRECISION, bnd_xl, bnd_xr, bnd_yl, bnd_yr, bnd_zl, bnd_zr, &
                               ierr, req, comm3d, procxl, procxr, procyl, procyr, proczl, proczr, status, &
                               pxsize, pysize, pzsize, procxyl, procyxl, pcoords, comm, &
                               MAG_YZ_LEFT_DOM, MAG_YZ_RIGHT_DOM, MAG_YZ_LEFT_BND, MAG_YZ_RIGHT_BND, &
                               MAG_XY_LEFT_DOM, MAG_XY_RIGHT_DOM, MAG_XY_LEFT_BND, MAG_XY_RIGHT_BND, &
                               MAG_XZ_LEFT_DOM, MAG_XZ_RIGHT_DOM, MAG_XZ_LEFT_BND, MAG_XZ_RIGHT_BND
#ifdef SHEAR
      use shear,         only: eps,delj
#endif /* SHEAR */

      implicit none
      character(len=*)  :: dim
      integer           :: i, j
      real, allocatable :: send_left(:,:,:,:),recv_left(:,:,:,:)
#ifdef SHEAR
      real, allocatable :: send_right(:,:,:,:),recv_right(:,:,:,:)
#endif /* SHEAR */
! MPI block comunication

      select case (dim)
         case ("xdim")
#ifdef SHEAR
            allocate(send_left(3,nb,ny,nz),send_right(3,nb,ny,nz), &
                     recv_left(3,nb,ny,nz),recv_right(3,nb,ny,nz))

            send_left (:,:,:,:)  = b(:,nb+1:2*nb,:,:)
            send_right(:,:,:,:)  = b(:,nxb+1:nxb+nb,:,:)

            if (bnd_xl == "she") then
!
! przesuwamy o calkowita liczbe komorek + periodyczny wb w kierunku y
!
               send_left (:,:,nb+1:nb+nyb,:)         = cshift(send_left (:,:,nb+1:nb+nyb,:),dim=3,shift= delj)
               send_left (:,:,1:nb,:)                = send_left  (:,:,nyb+1:nyb+nb,:)
               send_left (:,:,nb+nyb+1:nyb+2*nb,:)   = send_left  (:,:,nb+1:2*nb,:)
!
! remapujemy  - interpolacja kwadratowa
!
               send_left (:,:,:,:)  = (1.+eps)*(1.-eps) * send_left (:,:,:,:) &
                                      -0.5*eps*(1.-eps) * cshift(send_left (:,:,:,:),shift=-1,dim=3) &
                                      +0.5*eps*(1.+eps) * cshift(send_left (:,:,:,:),shift=1,dim=3)
            endif ! (bnd_xl == "she")


            if (bnd_xr == "she") then
!
! przesuwamy o calkowita liczbe komorek + periodyczny wb w kierunku y
!
               send_right (:,:,nb+1:nb+nyb,:)        = cshift(send_right(:,:,nb+1:nb+nyb,:),dim=3,shift=-delj)
               send_right (:,:,1:nb,:)               = send_right (:,:,nyb+1:nyb+nb,:)
               send_right (:,:,nb+nyb+1:nyb+2*nb,:)  = send_right (:,:,nb+1:2*nb,:)
!
! remapujemy - interpolacja kwadratowa
!
               send_right (:,:,:,:) = (1.+eps)*(1.-eps) * send_right (:,:,:,:) &
                                      -0.5*eps*(1.-eps) * cshift(send_right (:,:,:,:),shift=1,dim=3) &
                                      +0.5*eps*(1.+eps) * cshift(send_right (:,:,:,:),shift=-1,dim=3)
            endif ! (bnd_xr == "she")
!
! wysylamy na drugi brzeg
!
            CALL MPI_Isend   (send_left , 3*ny*nz*nb, MPI_DOUBLE_PRECISION, procxl, 10, comm, req(1), ierr)
            CALL MPI_Isend   (send_right, 3*ny*nz*nb, MPI_DOUBLE_PRECISION, procxr, 20, comm, req(3), ierr)
            CALL MPI_Irecv   (recv_left , 3*ny*nz*nb, MPI_DOUBLE_PRECISION, procxl, 20, comm, req(2), ierr)
            CALL MPI_Irecv   (recv_right, 3*ny*nz*nb, MPI_DOUBLE_PRECISION, procxr, 10, comm, req(4), ierr)

            call MPI_Waitall(4,req(:),status(:,:),ierr)

            b(:,1:nb-1,:,:)               = recv_left(:,1:nb-1,:,:)
            b(:,nxb+nb+1+1:nxb+2*nb,:,:)  = recv_right(:,1+1:nb,:,:)

            if (allocated(send_left))  deallocate(send_left)
            if (allocated(send_right)) deallocate(send_right)
            if (allocated(recv_left))  deallocate(recv_left)
            if (allocated(recv_right)) deallocate(recv_right)

!===============================================================================
#else /* !SHEAR */

            if (pxsize .gt. 1) then

               CALL MPI_Isend  (b(1,1,1,1), 1, MAG_YZ_LEFT_DOM,  procxl, 10, comm3d, req(1), ierr)
               CALL MPI_Isend  (b(1,1,1,1), 1, MAG_YZ_RIGHT_DOM, procxr, 20, comm3d, req(3), ierr)
               CALL MPI_Irecv  (b(1,1,1,1), 1, MAG_YZ_LEFT_BND,  procxl, 20, comm3d, req(2), ierr)
               CALL MPI_Irecv  (b(1,1,1,1), 1, MAG_YZ_RIGHT_BND, procxr, 10, comm3d, req(4), ierr)

               call MPI_Waitall(4,req(:),status(:,:),ierr)

            endif
#endif /* !SHEAR */


         case ("ydim")
            if (pysize .gt. 1) then

               CALL MPI_Isend  (b(1,1,1,1), 1, MAG_XZ_LEFT_DOM,  procyl, 30, comm3d, req(1), ierr)
               CALL MPI_Isend  (b(1,1,1,1), 1, MAG_XZ_RIGHT_DOM, procyr, 40, comm3d, req(3), ierr)
               CALL MPI_Irecv  (b(1,1,1,1), 1, MAG_XZ_LEFT_BND,  procyl, 40, comm3d, req(2), ierr)
               CALL MPI_Irecv  (b(1,1,1,1), 1, MAG_XZ_RIGHT_BND, procyr, 30, comm3d, req(4), ierr)

               call MPI_Waitall(4,req(:),status(:,:),ierr)
            endif

         case ("zdim")
            if (pzsize .gt. 1) then
               CALL MPI_Isend  (b(1,1,1,1), 1, MAG_XY_LEFT_DOM,  proczl, 50, comm3d, req(1), ierr)
               CALL MPI_Isend  (b(1,1,1,1), 1, MAG_XY_RIGHT_DOM, proczr, 60, comm3d, req(3), ierr)
               CALL MPI_Irecv  (b(1,1,1,1), 1, MAG_XY_LEFT_BND,  proczl, 60, comm3d, req(2), ierr)
               CALL MPI_Irecv  (b(1,1,1,1), 1, MAG_XY_RIGHT_BND, proczr, 50, comm3d, req(4), ierr)

               call MPI_Waitall(4,req(:),status(:,:),ierr)
            endif
      end select ! (dim)

! MPI + non-MPI corner-periodic boundary condition

      if (bnd_xl .eq. "cor") then
!   - lower to left
         if (pcoords(1) .eq. 0 .and. pcoords(2) .eq. 0) then
            do i=1,nb
               do j=nb+1,ny
                  b(ibx,i,j,:) = -b(iby,j,2*nb+1-i,:)
                  b(iby,i,j,:) =  b(ibx,j,2*nb+1-i,:)
                  b(ibz,i,j,:) =  b(ibz,j,2*nb+1-i,:)
               enddo
            enddo
         endif

         if (procxyl .gt. 0) then
            allocate(send_left(3,nb,ny,nz), recv_left(3,nx,nb,nz))

            send_left(:,:,:,:) = b(:,nb+1:2*nb,:,:)

            CALL MPI_Isend   (send_left , 3*nb*ny*nz, MPI_DOUBLE_PRECISION, procxyl, 70, comm, req(1), ierr)
            CALL MPI_Irecv   (recv_left , 3*nx*nb*nz, MPI_DOUBLE_PRECISION, procxyl, 80, comm, req(2), ierr)

            call MPI_Waitall(2,req(:),status(:,:),ierr)

            do i=1,nb
               do j=1,ny
                  b(ibx,i,j,:) = -recv_left(iby,j,nb+1-i,:)
                  b(iby,i,j,:) =  recv_left(ibx,j,nb+1-i,:)
                  b(ibz,i,j,:) =  recv_left(ibz,j,nb+1-i,:)
               enddo
            enddo

            if (allocated(send_left))  deallocate(send_left)
            if (allocated(recv_left))  deallocate(recv_left)
         endif
      endif

      if (bnd_yl .eq. "cor") then
!   - left to lower
         if (pcoords(2) .eq. 0 .and. pcoords(1) .eq. 0 ) then
            do j=1,nb
               do i=nb+1,nx
                  b(ibx,i,j,:) =  b(iby,2*nb+1-j,i,:)
                  b(iby,i,j,:) = -b(ibx,2*nb+1-j,i,:)
                  b(ibz,i,j,:) =  b(ibz,2*nb+1-j,i,:)
               enddo
            enddo
!   - interior to corner
            do j=1,nb
               do i=1,nb
                  b(ibx,i,j,:) =  -b(ibx,2*nb+1-i,2*nb+1-j,:)
                  b(iby,i,j,:) =  -b(iby,2*nb+1-i,2*nb+1-j,:)
                  b(ibz,i,j,:) =   b(ibz,2*nb+1-i,2*nb+1-j,:)
               enddo
            enddo
         endif

         if (procyxl .gt. 0) then
            allocate(send_left(3,nx,nb,nz), recv_left(3,nb,ny,nz))

            send_left(:,:,:,:) = b(:,:,nb+1:2*nb,:)

            CALL MPI_Isend   (send_left , 3*nx*nb*nz, MPI_DOUBLE_PRECISION, procyxl, 80, comm, req(1), ierr)
            CALL MPI_Irecv   (recv_left , 3*nb*ny*nz, MPI_DOUBLE_PRECISION, procyxl, 70, comm, req(2), ierr)

            call MPI_Waitall(2,req(:),status(:,:),ierr)

            do j=1,nb
               do i=1,nx
                  b(ibx,i,j,:) =  recv_left(iby,nb+1-j,i,:)
                  b(iby,i,j,:) = -recv_left(ibx,nb+1-j,i,:)
                  b(ibz,i,j,:) =  recv_left(ibz,nb+1-j,i,:)
               enddo
            enddo

            if (allocated(send_left))  deallocate(send_left)
            if (allocated(recv_left))  deallocate(recv_left)
         endif
      endif


! Non-MPI boundary conditions

      select case (dim)
         case ("xdim")

            select case (bnd_xl(1:3))
               case ("she")
!              Do nothing if "she"
               case ("mpi")
!              Do nothing if "mpi"
               case ("per")
                  b(:,1:nb,:,:)              = b(:,nxb+1:nxb+nb,:,:)
               case ("cor")
!              Do nothing if "cor"
               case ("inf")
!              Do nothing if "inf"
               case ("ref")
!              Do nothing if "ref"
               case ("out")
                  b(:,1,:,:) = b(:,2,:,:)
               case default
                  write(msg,'(4a)') "[magboundaries:bnd_b]: Boundary condition ",bnd_xl," not implemented in ",dim
                  call warn(msg)
            end select  ! (bnd_xl)

            select case (bnd_xr(1:3))
               case ("she")
!              Do nothing if "she"
               case ("mpi")
!              Do nothing if "mpi"
               case ("per")
                  b(:,nxb+nb+1:nxb+2*nb,:,:) = b(:,nb+1:2*nb,:,:)
               case ("cor")
!              Do nothing if "cor"
               case ("inf")
!              Do nothing if "inf"
               case ("ref")
!              Do nothing if "ref"
               case ("out")
                  b(:,nx,:,:) = b(:,nx-1,:,:)
               case default
                  write(msg,'(4a)') "[magboundaries:bnd_b]: Boundary condition ",bnd_xr," not implemented in ",dim
                  call warn(msg)
            end select  ! (bnd_xr)

         case ("ydim")

            select case (bnd_yl(1:3))
               case ("mpi")
!              Do nothing if "mpi"
               case ("per")
                  b(:,:,1:nb,:)              = b(:,:,nyb+1:nyb+nb,:)
               case ("cor")
!              Do nothing if "cor"
               case ("inf")
!              Do nothing if "inf"
               case ("ref")
!              Do nothing if "ref"
               case ("out")
                  b(:,:,1,:) = b(:,:,2,:)
               case default
                  write(msg,'(4a)') "[magboundaries:bnd_b]: Boundary condition ",bnd_yl," not implemented in ",dim
                  call warn(msg)
            end select  ! (bnd_yl)

            select case (bnd_yr(1:3))
               case ("mpi")
!              Do nothing if "mpi"
               case ("per")
                  b(:,:,nyb+nb+1:nyb+2*nb,:) = b(:,:,nb+1:2*nb,:)
               case ("cor")
!              Do nothing if "cor"
               case ("inf")
!              Do nothing if "inf"
               case ("ref")
!              Do nothing if "ref"
               case ("out")
                  b(:,:,ny,:) = b(:,:,ny-1,:)
               case default
                  write(msg,'(4a)') "[magboundaries:bnd_b]: Boundary condition ",bnd_yr," not implemented in ",dim
                  call warn(msg)

            end select  ! (bnd_yr)


         case ("zdim")

            select case (bnd_zl(1:3))
               case ("mpi")
!              Do nothing if "mpi"
               case ("per")
                  b(:,:,:,1:nb)              = b(:,:,:,nzb+1:nzb+nb)
               case ("ref")
!              Do nothing if "ref"
               case ("out")
                  b(:,:,:,1) = b(:,:,:,2)
               case default
                  write(msg,'(4a)') "[magboundaries:bnd_b]: Boundary condition ",bnd_zl," not implemented in ",dim
                  call warn(msg)
            end select  ! (bnd_zl)

            select case (bnd_zr(1:3))
               case ("mpi")
!              Do nothing if "mpi"
               case ("per")
                  b(:,:,:,nzb+nb+1:nzb+2*nb) = b(:,:,:,nb+1:2*nb)
               case ("ref")
!              Do nothing if "ref"
               case ("out")
                  b(:,:,:,nz) = b(:,:,:,nz-1)
               case default
                  write(msg,'(4a)') "[magboundaries:bnd_b]: Boundary condition ",bnd_zr," not implemented in ",dim
                  call warn(msg)
            end select  ! (bnd_zr)

      end select  ! (dim)

   end subroutine bnd_b

!=====================================================================================================

   subroutine bnd_emf(var, name, dim)
      use dataio_pub, only: msg, warn
      use grid,          only: nx, ny, nz, nb, nxb, nyb, nzb
      use mpisetup,      only: bnd_xl, bnd_xr, bnd_yl, bnd_yr, bnd_zl, bnd_zr

      implicit none
      real, dimension(nx,ny,nz) :: var
      real, dimension(ny,nz)    :: dvarx
      real, dimension(nx,nz)    :: dvary
      real, dimension(nx,ny)    :: dvarz
      character(len=4)          :: name, dim
      integer                   :: ib


      select case (dim)
         case ("xdim")

            select case (name)
               case ("vxby","vxbz")

                  select case (bnd_xl(1:3))
                     case ("she")
!                    Do nothing if "she"
                     case ("mpi")
!                    Do nothing if "mpi"
                     case ("per")
!                    Do nothing if "per"
                     case ("cor")
!                    Do nothing if "cor"
                     case ("inf")
!                    Do nothing if "inf"
                     case ("ref")
                        var(nb,:,:)             = 0.0
                        do ib=1,nb-1
                           var(nb-ib,:,:)       = -var(nb+ib,:,:)
                        enddo
                     case ("out")
                        dvarx = var(nb+1,:,:)-var(nb,:,:)
                        do ib=1,nb-1
                           var(ib,:,:)          = var(nb+1,:,:)  - real(nb+1-ib)*dvarx
                        enddo
                     case default
                        write(msg,'(6a)') "[magboundaries:bnd_emf]: Boundary condition ",bnd_xl," not implemented for ",name, " in ", dim
                        call warn(msg)
                  end select  ! (bnd_xl)

                  select case (bnd_xr(1:3))
                     case ("she")
!                    Do nothing if "she"
                     case ("mpi")
!                    Do nothing if "mpi"
                     case ("per")
!                    Do nothing if "per"
                     case ("cor")
!                    Do nothing if "cor"
                     case ("inf")
!                    Do nothing if "inf"
                     case ("ref")
                        var(nb+nxb,:,:)         = 0.0
                        do ib=1,nb-1
                           var(nb+nxb+ib,:,:)   = -var(nb+nxb-ib,:,:)
                        enddo
                     case ("out")
                        dvarx = var(nb+nxb,:,:)-var(nb+nxb-1,:,:)
                        do ib=1,nb-1
                           var(nb+nxb+ib,:,:)   = var(nb+nxb,:,:) + real(ib)*dvarx
                        enddo
                     case default
                        write(msg,'(6a)') "[magboundaries:bnd_emf]: Boundary condition ",bnd_xr," not implemented for ",name, " in ", dim
                        call warn(msg)
                  end select ! (bnd_xr)

               case ("vybx","vzbx","emfy","emfz")

                  select case (bnd_xl(1:3))
                     case ("she")
!                    Do nothing if "she"
                     case ("mpi")
!                    Do nothing if "mpi"
                     case ("per")
!                    Do nothing if "per"
                     case ("cor")
!                    Do nothing if "cor"
                     case ("inf")
!                    Do nothing if "inf"
                     case ("ref")
                        var(nb+1,:,:)           = 0.0
                        do ib=1,nb
                           var(nb+1-ib,:,:)     = -var(nb+1+ib,:,:)
                        enddo
                     case ("out")
                        dvarx = var(nb+2,:,:)-var(nb+1,:,:)
                        do ib=1,nb
                           var(ib,:,:)          = var(nb+1,:,:)  - real(nb+1-ib)*dvarx
                        enddo
                     case default
                        write(msg,'(6a)') "[magboundaries:bnd_emf]: Boundary condition ",bnd_xl," not implemented for ",name, " in ", dim
                        call warn(msg)
                  end select   ! (bnd_xl)

                  select case (bnd_xr(1:3))
                     case ("she")
!                    Do nothing if "she"
                     case ("mpi")
!                    Do nothing if "mpi"
                     case ("per")
!                    Do nothing if "per"
                     case ("cor")
!                    Do nothing if "cor"
                     case ("inf")
!                    Do nothing if "inf"
                     case ("ref")
                        var(nb+nxb+1,:,:)       = 0.0
                        do ib=1,nb-1
                           var(nb+nxb+1+ib,:,:) = -var(nb+nxb+1-ib,:,:)
                        enddo
                     case ("out")
                        dvarx = var(nb+nxb+1,:,:)-var(nb+nxb,:,:)
                        do ib=1,nb-1
                           var(nb+nxb+1+ib,:,:) =  var(nb+nxb+1,:,:) + real(ib)*dvarx
                        enddo
                     case default
                        write(msg,'(6a)') "[magboundaries:bnd_emf]: Boundary condition ",bnd_xr," not implemented for ",name, " in ", dim
                        call warn(msg)
                  end select  ! (bnd_xr)

               case ("vybz","vzby","emfx")

                  select case (bnd_xl(1:3))
                     case ("she")
!                    Do nothing if "she"
                     case ("mpi")
!                    Do nothing if "mpi"
                     case ("per")
!                    Do nothing if "per"
                     case ("cor")
!                    Do nothing if "cor"
                     case ("inf")
!                    Do nothing if "inf"
                     case ("ref")
                        do ib=1,nb
                           var(nb+1-ib,:,:)     = var(nb+ib,:,:)
                        enddo
                     case ("out")
                        dvarx = var(nb+1,:,:)-var(nb,:,:)
                        do ib=1,nb
                           var(ib,:,:)          = var(nb+1,:,:)  - real(nb+1-ib)*dvarx
                        enddo
                     case default
                        write(msg,'(6a)') "[magboundaries:bnd_emf]: Boundary condition ",bnd_xl," not implemented for ",name, " in ", dim
                        call warn(msg)
                  end select   ! (bnd_xl)

                  select case (bnd_xr(1:3))
                     case ("she")
!                    Do nothing if "she"
                     case ("mpi")
!                    Do nothing if "mpi"
                     case ("per")
!                    Do nothing if "per"
                     case ("cor")
!                    Do nothing if "cor"
                     case ("inf")
!                    Do nothing if "inf"
                     case ("ref")
                        do ib=1,nb
                           var(nb+nxb+ib,:,:)   = var(nb+nxb+1-ib,:,:)
                        enddo
                     case ("out")
                        dvarx = var(nb+nxb+1,:,:)-var(nb+nxb,:,:)
                        do ib=1,nb
                           var(nb+nxb+ib,:,:)   = var(nb+nxb,:,:) + real(ib)*dvarx
                        enddo
                     case default
                        write(msg,'(6a)') "[magboundaries:bnd_emf]: Boundary condition ",bnd_xr," not implemented for ",name, " in ", dim
                        call warn(msg)
                  end select   ! (bnd_xr)

            end select  ! (name)

         case ("ydim")

            select case (name)
               case ("vybz","vybx")

                  select case (bnd_yl(1:3))
                     case ("mpi")
!                    Do nothing if "mpi"
                     case ("per")
!                    Do nothing if "per"
                     case ("cor")
!                    Do nothing if "cor"
                     case ("inf")
!                    Do nothing if "inf"
                     case ("ref")
                        var(:,nb,:)             = 0.0
                        do ib=1,nb-1
                           var(:,nb-ib,:)       = -var(:,nb+ib,:)
                        enddo
                     case ("out")
                        dvary = var(:,nb+1,:)-var(:,nb,:)
                        do ib=1,nb-1
                           var(:,ib,:)          = var(:,nb+1,:) - real(nb+1-ib)*dvary
                        enddo
                     case default
                        write(msg,'(6a)') "[magboundaries:bnd_emf]: Boundary condition ",bnd_yl," not implemented for ",name, " in ", dim
                        call warn(msg)
                  end select  ! (bnd_yl)

                  select case (bnd_yr(1:3))
                     case ("mpi")
!                    Do nothing if "mpi"
                     case ("per")
!                    Do nothing if "per"
                     case ("cor")
!                    Do nothing if "cor"
                     case ("inf")
!                    Do nothing if "inf"
                     case ("ref")
                        var(:,nb+nyb,:)         = 0.0
                        do ib=1,nb-1
                           var(:,nb+nyb+ib,:)   = -var(:,nb+nyb-ib,:)
                        enddo
                     case ("out")
                        dvary = var(:,nb+nyb,:)-var(:,nb+nyb-1,:)
                        do ib=1,nb-1
                           var(:,nb+nyb+ib,:)   = var(:,nb+nyb,:) + real(ib)*dvary
                        enddo
                     case default
                        write(msg,'(6a)') "[magboundaries:bnd_emf]: Boundary condition ",bnd_yr," not implemented for ",name, " in ", dim
                        call warn(msg)
                  end select ! (bnd_yr)

               case ("vzby","vxby","emfz","emfx")

                  select case (bnd_yl(1:3))
                     case ("mpi")
!                    Do nothing if "mpi"
                     case ("per")
!                    Do nothing if "per"
                     case ("cor")
!                    Do nothing if "cor"
                     case ("inf")
!                    Do nothing if "inf"
                     case ("ref")
                        var(:,nb+1,:)           = 0.0
                        do ib=1,nb
                           var(:,nb+1-ib,:)     = -var(:,nb+1+ib,:)
                        enddo
                     case ("out")
                        dvary = var(:,nb+2,:)-var(:,nb+1,:)
                        do ib=1,nb
                           var(:,ib,:)          = var(:,nb+1,:) - real(nb+1-ib)*dvary
                        enddo
                     case default
                        write(msg,'(6a)') "[magboundaries:bnd_emf]: Boundary condition ",bnd_yl," not implemented for ",name, " in ", dim
                        call warn(msg)
                  end select   ! (bnd_yl)

                  select case (bnd_yr(1:3))
                     case ("mpi")
!                    Do nothing if "mpi"
                     case ("per")
!                    Do nothing if "per"
                     case ("cor")
!                    Do nothing if "cor"
                     case ("inf")
!                    Do nothing if "inf"
                     case ("ref")
                        var(:,nb+nyb+1,:)       = 0.0
                        do ib=1,nb-1
                           var(:,nb+nyb+1+ib,:) = -var(:,nb+nyb+1-ib,:)
                        enddo
                     case ("out")
                        dvary = var(:,nb+nyb+1,:)-var(:,nb+nyb,:)
                        do ib=1,nb-1
                           var(:,nb+nyb+1+ib,:) =  var(:,nb+nyb+1,:) + real(ib)*dvary
                        enddo
                     case default
                        write(msg,'(6a)') "[magboundaries:bnd_emf]: Boundary condition ",bnd_yr," not implemented for ",name, " in ", dim
                        call warn(msg)
                  end select  ! (bnd_yr)

               case ("vzbx","vxbz","emfy")

                  select case (bnd_yl(1:3))
                     case ("mpi")
!                    Do nothing if "mpi"
                     case ("per")
!                    Do nothing if "per"
                     case ("cor")
!                    Do nothing if "cor"
                     case ("inf")
!                    Do nothing if "inf"
                     case ("ref")
                        do ib=1,nb
                           var(:,nb+1-ib,:)     = var(:,nb+ib,:)
                        enddo
                     case ("out")
                        dvary = var(:,nb+1,:)-var(:,nb,:)
                        do ib=1,nb
                           var(:,ib,:)          = var(:,nb+1,:) - real(nb+1-ib)*dvary
                        enddo
                     case default
                        write(msg,'(6a)') "[magboundaries:bnd_emf]: Boundary condition ",bnd_yl," not implemented for ",name, " in ", dim
                        call warn(msg)
                  end select   ! (bnd_yl)

                  select case (bnd_yr(1:3))
                     case ("mpi")
!                    Do nothing if "mpi"
                     case ("per")
!                    Do nothing if "per"
                     case ("cor")
!                    Do nothing if "cor"
                     case ("inf")
!                    Do nothing if "inf"
                     case ("ref")
                        do ib=1,nb
                           var(:,nb+nyb+ib,:)   = var(:,nb+nyb+1-ib,:)
                        enddo
                     case ("out")
                        dvary = var(:,nb+nyb+1,:)-var(:,nb+nyb,:)
                        do ib=1,nb
                           var(:,nb+nyb+ib,:)   = var(:,nb+nyb,:) + real(ib)*dvary
                        enddo
                     case default
                        write(msg,'(6a)') "[magboundaries:bnd_emf]: Boundary condition ",bnd_yr," not implemented for ",name, " in ", dim
                        call warn(msg)
                  end select   ! (bnd_yr)

            end select  ! (name)


         case ("zdim")

            select case (name)

               case ("vzbx","vzby")

                  select case (bnd_zl(1:3))
                     case ("mpi")
!                    Do nothing if "mpi"
                     case ("per")
!                    Do nothing if "per"
                     case ("cor")
!                    Do nothing if "cor"
                     case ("inf")
!                    Do nothing if "inf"
                     case ("ref")
                        var(:,:,nb)             = 0.0
                        do ib=1,nb-1
                           var(:,:,nb-ib)       = -var(:,:,nb+ib)
                        enddo
                     case ("out")
                        dvarz = var(:,:,nb+1)-var(:,:,nb)
                        do ib=1,nb-1
                           var(:,:,ib)          = var(:,:,nb+1) - real(nb+1-ib)*dvarz
                        enddo
                     case default
                        write(msg,'(6a)') "[magboundaries:bnd_emf]: Boundary condition ",bnd_zl," not implemented for ",name, " in ", dim
                        call warn(msg)
                  end select  ! (bnd_zl)

                  select case (bnd_zr(1:3))
                     case ("mpi")
!                    Do nothing if "mpi"
                     case ("per")
!                    Do nothing if "per"
                     case ("cor")
!                    Do nothing if "cor"
                     case ("inf")
!                    Do nothing if "inf"
                     case ("ref")
                        var(:,:,nb+nzb)         = 0.0
                        do ib=1,nb-1
                           var(:,:,nb+nzb+ib)   = -var(:,:,nb+nzb-ib)
                        enddo
                     case ("out")
                        dvarz = var(:,:,nb+nzb)-var(:,:,nb+nzb-1)
                        do ib=1,nb-1
                           var(:,:,nb+nzb+ib)   = var(:,:,nb+nzb) + real(ib)*dvarz
                        enddo
                     case default
                        write(msg,'(6a)') "[magboundaries:bnd_emf]: Boundary condition ",bnd_zr," not implemented for ",name, " in ", dim
                        call warn(msg)
                  end select ! (bnd_zr)


               case ("vxbz","vybz","emfy","emfx")

                  select case (bnd_zl(1:3))
                     case ("mpi")
!                    Do nothing if "mpi"
                     case ("per")
!                    Do nothing if "per"
                     case ("cor")
!                    Do nothing if "cor"
                     case ("inf")
!                    Do nothing if "inf"
                     case ("ref")
                        var(:,:,nb+1)           = 0.0
                        do ib=1,nb
                           var(:,:,nb+1-ib)     = -var(:,:,nb+1+ib)
                        enddo
                     case ("out")
                        dvarz = var(:,:,nb+2)-var(:,:,nb+1)
                        do ib=1,nb
                           var(:,:,ib)          = var(:,:,nb+1) - real(nb+1-ib)*dvarz
                        enddo
                     case default
                        write(msg,'(6a)') "[magboundaries:bnd_emf]: Boundary condition ",bnd_zl," not implemented for ",name, " in ", dim
                        call warn(msg)
                  end select   ! (bnd_zl)

                  select case (bnd_zr(1:3))
                     case ("mpi")
!                    Do nothing if "mpi"
                     case ("per")
!                    Do nothing if "per"
                     case ("cor")
!                    Do nothing if "cor"
                     case ("inf")
!                    Do nothing if "inf"
                     case ("ref")
                        var(:,:,nb+nzb+1)       = 0.0
                        do ib=1,nb-1
                           var(:,:,nb+nzb+1+ib) = -var(:,:,nb+nzb+1-ib)
                        enddo
                     case ("out")
                        dvarz = var(:,:,nb+nzb+1)-var(:,:,nb+nzb)
                        do ib=1,nb-1
                           var(:,:,nb+nzb+1+ib) =  var(:,:,nb+nzb+1) + real(ib)*dvarz
                        enddo
                     case default
                        write(msg,'(6a)') "[magboundaries:bnd_emf]: Boundary condition ",bnd_zr," not implemented for ",name, " in ", dim
                        call warn(msg)
                  end select   ! (bnd_zr)

               case ("vxby","vybx","emfz")

                  select case (bnd_zl(1:3))
                     case ("mpi")
!                    Do nothing if "mpi"
                     case ("per")
!                    Do nothing if "per"
                     case ("cor")
!                    Do nothing if "cor"
                     case ("inf")
!                    Do nothing if "inf"
                     case ("ref")
                        do ib=1,nb
                           var(:,:,nb+1-ib)          = var(:,:,nb+ib)
                        enddo
                     case ("out")
                        dvarz = var(:,:,nb+1)-var(:,:,nb)
                        do ib=1,nb
                           var(:,:,ib)               = var(:,:,nb+1) - real(nb+1-ib)*dvarz
                        enddo
                     case default
                        write(msg,'(6a)') "[magboundaries:bnd_emf]: Boundary condition ",bnd_zl," not implemented for ",name, " in ", dim
                        call warn(msg)
                  end select   ! (bnd_zl)

                  select case (bnd_zr(1:3))
                     case ("mpi")
!                    Do nothing if "mpi"
                     case ("per")
!                    Do nothing if "per"
                     case ("cor")
!                    Do nothing if "cor"
                     case ("inf")
!                    Do nothing if "inf"
                     case ("ref")
                        do ib=1,nb
                           var(:,:,nb+nzb+ib)        = var(:,:,nb+nzb+1-ib)
                        enddo
                     case ("out")
                        dvarz = var(:,:,nb+nzb+1)-var(:,:,nb+nzb)
                        do ib=1,nb
                           var(:,:,nb+nzb+ib)        = var(:,:,nb+nzb) + real(ib)*dvarz
                        enddo
                     case default
                        write(msg,'(6a)') "[magboundaries:bnd_emf]: Boundary condition ",bnd_zr," not implemented for ",name, " in ", dim
                        call warn(msg)
                  end select ! (bnd_zr)

            end select  ! (name)

      end select ! (dim)
   end subroutine bnd_emf

   subroutine all_mag_boundaries
      use grid, only: nxd, nyd, nzd
      implicit none

      if (nxd /= 1) call bnd_b("xdim")
      if (nyd /= 1) call bnd_b("ydim")
      if (nzd /= 1) call bnd_b("zdim")

   end subroutine all_mag_boundaries

end module magboundaries
