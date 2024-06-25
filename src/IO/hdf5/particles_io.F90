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
#include "piernik.h"

module particles_io
! pulled by NBODY && HDF5

   use constants, only: dsetnamelen

   implicit none

   private
   public :: init_nbody_hdf5, pvarn, serial_nbody_datafields, parallel_nbody_datafields

   character(len=dsetnamelen), dimension(*), parameter  :: pvarn = ['ppid', 'mass', 'ener', 'posx', 'posy', 'posz', 'velx', 'vely', 'velz', 'accx', 'accy', 'accz', 'tfor', 'tdyn']
   logical,                    dimension(size(pvarn))   :: pvarl = .false.
   integer(kind=4), dimension(:), allocatable           :: tabi
   real,            dimension(:), allocatable           :: tabr

contains

   subroutine init_nbody_hdf5(pvars)

      use common_hdf5, only: pdsets
      use dataio_pub,  only: msg, warn
      use mpisetup,    only: master

      implicit none

      character(len=dsetnamelen), dimension(:), intent(inout) :: pvars  !< quantities to be plotted, see dataio::vars
      integer                                                 :: ie, il, k
      logical                                                 :: var_found

      var_found = .false.
      do il = lbound(pvars, 1), ubound(pvars, 1)
         if (len(trim(pvars(il))) == 0) cycle
         var_found = .true.
      enddo
      if (.not. var_found) pvars(1:size(pvarn)) = pvarn

      k = 0
      do il = lbound(pvars, 1), ubound(pvars, 1)
         if (len(trim(pvars(il))) == 0) cycle
         var_found = .false.
         do ie = lbound(pvarn, 1), ubound(pvarn, 1)
            if (trim(pvars(il)) == trim(pvarn(ie))) then
               pvarl(ie) = .true.
               var_found = .true.
               k = k + 1
            endif
         enddo
         if (.not. var_found) then
            write(msg,'(2a)')'[particles_io::init_nbody_hdf5]: unknown particle var: ', pvars(il)
            if (master) call warn(msg)
         endif
      enddo

      allocate(pdsets(k))
      pdsets = pack(pvarn, pvarl)

   end subroutine init_nbody_hdf5

   subroutine parallel_nbody_datafields(group_id, pvars, ncg, cg)

      use grid_cont, only: grid_container
      use hdf5,      only: HID_T

      implicit none

      integer(HID_T), dimension(:,:), intent(in)    :: group_id       !< File identifier
      character(len=*), dimension(:), intent(in)    :: pvars
      integer(kind=4),                intent(in)    :: ncg
      type(grid_container), pointer,  intent(inout) :: cg
      integer(kind=4)                               :: ivar, n_part

      n_part = cg%count_particles()
      if (n_part == 0) return

      allocate(tabi(n_part), tabr(n_part))
      do ivar = lbound(pvars, dim=1, kind=4), ubound(pvars, dim=1, kind=4)
         select case (pvars(ivar))
            case ('id')
               call collect_intr(pvars(ivar), cg)
               call write_nbody_h5_intr(group_id(ncg, ivar))
            case default
               call collect_real(pvars(ivar), cg)
               call write_nbody_h5_real(group_id(ncg, ivar))
         end select
      enddo
      deallocate(tabi, tabr)

   end subroutine parallel_nbody_datafields

   subroutine serial_nbody_datafields(group_id, pvars, ncg, cg_src_ncg, proc_ncg, tot_cg_n)

      use common_hdf5,    only: get_nth_cg, hdf_vars
      use constants,      only: I_ONE
      use grid_cont,      only: grid_container
      use hdf5,           only: HID_T
      use MPIF,           only: MPI_INTEGER, MPI_STATUS_IGNORE, MPI_COMM_WORLD
      use MPIFUN,         only: MPI_Recv, MPI_Send
      use mpisetup,       only: master, FIRST, proc, err_mpi

      implicit none

      integer(HID_T), dimension(:,:), intent(in) :: group_id       !< File identifier
      character(len=*), dimension(:), intent(in) :: pvars
      integer(kind=4),                intent(in) :: ncg, cg_src_ncg, proc_ncg, tot_cg_n
      integer(kind=4)                            :: ptag, ivar, n_part
      type(grid_container), pointer              :: cg

      if (.not. master .and. proc_ncg /= proc) return

      if (proc_ncg == proc) then
         cg => get_nth_cg(cg_src_ncg)
         n_part = cg%count_particles()
      endif

      ptag = ncg + tot_cg_n * (ubound(hdf_vars, 1, kind=4) + I_ONE)
      if (master) then
         if (proc_ncg /= proc) call MPI_Recv(n_part, I_ONE, MPI_INTEGER, proc_ncg, ptag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
      else
         if (proc_ncg == proc) call MPI_Send(n_part, I_ONE, MPI_INTEGER, FIRST, ptag, MPI_COMM_WORLD, err_mpi)
      endif

      if (n_part == 0) return

      allocate(tabi(n_part), tabr(n_part))
      do ivar = lbound(pvars, 1, kind=4), ubound(pvars, 1, kind=4)
         select case (pvars(ivar))
            case ('id')
               call serial_write_intr(group_id, ncg, ivar, pvars(ivar), n_part, cg_src_ncg, proc_ncg, ptag+ivar)
            case default
               call serial_write_real(group_id, ncg, ivar, pvars(ivar), n_part, cg_src_ncg, proc_ncg, ptag+ivar)
         end select
      enddo
      deallocate(tabi, tabr)

   end subroutine serial_nbody_datafields

   subroutine collect_intr(pvar, cg)

      use grid_cont,      only: grid_container
      use particle_types, only: particle

      implicit none

      character(len=*),               intent(in) :: pvar
      type(grid_container), pointer,  intent(in) :: cg

      integer                                    :: cgnp
      type(particle), pointer                    :: pset

      cgnp = 0
      select case (pvar)
         case ('id')
            pset => cg%pset%first
            do while (associated(pset))
               if (pset%pdata%phy) then
                  cgnp = cgnp + 1
                  tabi(cgnp) = pset%pdata%pid
               endif
               pset => pset%nxt
            enddo
         case default
      end select

   end subroutine collect_intr

   subroutine collect_real(pvar, cg)

      use constants,      only: xdim, ydim, zdim
      use grid_cont,      only: grid_container
      use particle_types, only: particle

      implicit none

      character(len=*),              intent(in) :: pvar
      type(grid_container), pointer, intent(in) :: cg

      integer                                   :: cgnp
      type(particle), pointer                   :: pset

      cgnp = 0
      pset => cg%pset%first
      do while (associated(pset))
         if (pset%pdata%phy) then
            cgnp = cgnp + 1
            select case (pvar)
               case ('mass')
                  tabr(cgnp) = pset%pdata%mass
               case ('energy')
                  tabr(cgnp) = pset%pdata%energy
               case ('position_x')
                  tabr(cgnp) = pset%pdata%pos(xdim)
               case ('position_y')
                  tabr(cgnp) = pset%pdata%pos(ydim)
               case ('position_z')
                  tabr(cgnp) = pset%pdata%pos(zdim)
               case ('velocity_x')
                  tabr(cgnp) = pset%pdata%vel(xdim)
               case ('velocity_y')
                  tabr(cgnp) = pset%pdata%vel(ydim)
               case ('velocity_z')
                  tabr(cgnp) = pset%pdata%vel(zdim)
               case ('acceleration_x')
                  tabr(cgnp) = pset%pdata%acc(xdim)
               case ('acceleration_y')
                  tabr(cgnp) = pset%pdata%acc(ydim)
               case ('acceleration_z')
                  tabr(cgnp) = pset%pdata%acc(zdim)
               case ('formation_time')
                  tabr(cgnp) = pset%pdata%tform
               case ('dynamical_time')
                  tabr(cgnp) = pset%pdata%tdyn
               case default
            end select
         endif
         pset => pset%nxt
      enddo

   end subroutine collect_real

   subroutine serial_write_intr(group_id, ncg, ivar, pvar, n_part, cg_src_ncg, proc_ncg, ptag)

      use common_hdf5, only: get_nth_cg
      use hdf5,        only: HID_T
      use MPIF,        only: MPI_INTEGER, MPI_STATUS_IGNORE, MPI_COMM_WORLD
      use MPIFUN,      only: MPI_Recv, MPI_Send
      use mpisetup,    only: master, FIRST, proc, err_mpi

      implicit none

      integer(HID_T), dimension(:,:), intent(in) :: group_id       !< group identifiers
      character(len=*),               intent(in) :: pvar
      integer(kind=4),                intent(in) :: ncg, ivar, n_part, cg_src_ncg, proc_ncg, ptag

      if (proc_ncg == proc) then
         call collect_intr(pvar, get_nth_cg(cg_src_ncg))
         if (.not. master) call MPI_Send(tabi, n_part, MPI_INTEGER, FIRST, ptag, MPI_COMM_WORLD, err_mpi)
      endif

      if (master) then
         if (proc_ncg /= proc) call MPI_Recv(tabi, n_part, MPI_INTEGER, proc_ncg, ptag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
         call write_nbody_h5_intr(group_id(ncg, ivar))
      endif

   end subroutine serial_write_intr

   subroutine serial_write_real(group_id, ncg, ivar, pvar, n_part, cg_src_ncg, proc_ncg, ptag)

      use common_hdf5, only: get_nth_cg
      use hdf5,        only: HID_T
      use MPIF,        only: MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, MPI_COMM_WORLD
      use MPIFUN,      only: MPI_Recv, MPI_Send
      use mpisetup,    only: master, FIRST, proc, err_mpi

      implicit none

      integer(HID_T), dimension(:,:), intent(in) :: group_id       !< group identifiers
      character(len=*),               intent(in) :: pvar
      integer(kind=4),                intent(in) :: ncg, ivar, n_part, cg_src_ncg, proc_ncg, ptag

      if (proc_ncg == proc) then
         call collect_real(pvar, get_nth_cg(cg_src_ncg))
         if (.not. master) call MPI_Send(tabr, n_part, MPI_DOUBLE_PRECISION, FIRST, ptag, MPI_COMM_WORLD, err_mpi)
      endif

      if (master) then
         if (proc_ncg /= proc) call MPI_Recv(tabr, n_part, MPI_DOUBLE_PRECISION, proc_ncg, ptag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
         call write_nbody_h5_real(group_id(ncg, ivar))
      endif

   end subroutine serial_write_real

   subroutine write_nbody_h5_intr(dataset_id)

      use dataio_pub, only: die
      use hdf5,       only: h5dwrite_f, HID_T, HSIZE_T, H5T_NATIVE_INTEGER

      implicit none

      integer(HID_T),     intent(in) :: dataset_id
      integer(HSIZE_T), dimension(1) :: dimm
      integer(kind=4)                :: error

      if (all(kind(dataset_id) /= [4, 8])) call die("[particles_io:write_nbody_h5_intr] HID_T doesn't fit to MPI_INTEGER8")

      dimm = shape(tabi)
      call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, tabi, dimm, error)  ! beware: 64-bit tabi(:) produces "no specific subroutine for the generic ‘h5dwrite_f'" error

   end subroutine write_nbody_h5_intr

   subroutine write_nbody_h5_real(dataset_id)

      use dataio_pub, only: die
      use hdf5,       only: h5dwrite_f, HID_T, HSIZE_T, H5T_NATIVE_DOUBLE

      implicit none

      integer(HID_T),     intent(in) :: dataset_id
      integer(HSIZE_T), dimension(1) :: dimm
      integer(kind=4)                :: error

      if (all(kind(dataset_id) /= [4, 8])) call die("[particles_io:write_nbody_h5_real] HID_T doesn't fit to MPI_INTEGER8")

      dimm = shape(tabr)
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, tabr, dimm, error)

   end subroutine write_nbody_h5_real

end module particles_io
