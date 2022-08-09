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

module cg_particles_io
! pulled by NBODY && HDF5

   use constants, only: dsetnamelen

   implicit none

   private
   public :: init_nbody_hdf5, pvarn, serial_nbody_datafields, parallel_nbody_datafields

   character(len=dsetnamelen), dimension(*), parameter  :: pvarn = ['ppid', 'mass', 'ener', 'posx', 'posy', 'posz', 'velx', 'vely', 'velz', 'accx', 'accy', 'accz', 'tfor', 'tdyn']
   logical,                    dimension(size(pvarn))   :: pvarl = .false.

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
            write(msg,'(2a)')'[particles_io_hdf5::init_nbody_hdf5]: unknown particle var: ', pvars(il)
            if (master) call warn(msg)
         endif
      enddo

#ifdef NBODY_1FILE
      allocate(pdsets(k))
      pdsets = pack(pvarn, pvarl)
#endif /* NBODY_1FILE */

   end subroutine init_nbody_hdf5

   subroutine nbody_datafields(group_id, pvar, n_part, cg)

      use grid_cont, only: grid_container
      use hdf5,      only: HID_T

      implicit none

      integer(HID_T),                intent(in)    :: group_id       !< File identifier
      character(len=*),              intent(in)    :: pvar
      integer(kind=4),               intent(in)    :: n_part
      type(grid_container), pointer, intent(inout) :: cg

      select case (pvar)
         case ('id')
            call parallel_write_intr1(group_id, pvar, n_part, cg)
         case default
            call parallel_write_rank1(group_id, pvar, n_part, cg)
      end select

   end subroutine nbody_datafields

   subroutine parallel_nbody_datafields(group_id, pvars, ncg, cg)

      use grid_cont,      only: grid_container
      use hdf5,           only: HID_T
      use particle_utils, only: count_cg_particles

      implicit none

      integer(HID_T), dimension(:,:), intent(in)    :: group_id       !< File identifier
      character(len=*), dimension(:), intent(in)    :: pvars
      integer(kind=4),                intent(in)    :: ncg
      type(grid_container), pointer,  intent(inout) :: cg
      integer(kind=4)                               :: i, n_part

      n_part = count_cg_particles(cg)
      if (n_part == 0) return

      do i = lbound(pvars, dim=1, kind=4), ubound(pvars, dim=1, kind=4)
         call nbody_datafields(group_id(ncg, i), pvars(i), n_part, cg)
      enddo

   end subroutine parallel_nbody_datafields

   subroutine serial_nbody_datafields(group_id, pvars, ncg, cg_src_ncg, proc_ncg, tot_cg_n)

      use common_hdf5, only: hdf_vars
      use hdf5,        only: HID_T

      implicit none

      integer(HID_T), dimension(:,:), intent(in) :: group_id       !< File identifier
      character(len=*), dimension(:), intent(in) :: pvars
      integer(kind=4),                intent(in) :: ncg, cg_src_ncg, proc_ncg, tot_cg_n
      integer(kind=4)                            :: ptag, ivar

      do ivar = lbound(pvars, 1, kind=4), ubound(pvars, 1, kind=4)
         ptag = ncg + tot_cg_n * (ubound(hdf_vars, 1, kind=4) + 2*ivar)
         select case (pvars(ivar))
            case ('id')
               call serial_write_intr1(group_id, ncg, ivar, pvars(ivar), cg_src_ncg, proc_ncg, ptag)
            case default
               call serial_write_rank1(group_id, ncg, ivar, pvars(ivar), cg_src_ncg, proc_ncg, ptag)
         end select
      enddo

   end subroutine serial_nbody_datafields

   subroutine collect_intr1(pvar, cg, tabi1)

      use grid_cont,      only: grid_container
      use particle_types, only: particle

      implicit none

      character(len=*),              intent(in)    :: pvar
      type(grid_container), pointer, intent(inout) :: cg

      integer                                      :: cgnp
      integer(kind=4), dimension(:), allocatable   :: tabi1
      type(particle), pointer                      :: pset

      cgnp = 0
      select case (pvar)
         case ('id')
            pset => cg%pset%first
            do while (associated(pset))
               if (pset%pdata%phy) then
                  cgnp = cgnp + 1
                  tabi1(cgnp) = pset%pdata%pid
               endif
               pset => pset%nxt
            enddo
         case default
      end select

   end subroutine collect_intr1

   subroutine collect_rank1(pvar, cg, tabr1)

      use constants,      only: xdim, ydim, zdim
      use grid_cont,      only: grid_container
      use particle_types, only: particle

      implicit none

      character(len=*),              intent(in)    :: pvar
      type(grid_container), pointer, intent(inout) :: cg

      integer                                      :: cgnp
      real, dimension(:), allocatable              :: tabr1
      type(particle), pointer                      :: pset

      cgnp = 0
      pset => cg%pset%first
      do while (associated(pset))
         if (pset%pdata%phy) then
            cgnp = cgnp + 1
            select case (pvar)
               case ('mass')
                  tabr1(cgnp) = pset%pdata%mass
               case ('energy')
                  tabr1(cgnp) = pset%pdata%energy
               case ('position_x')
                  tabr1(cgnp) = pset%pdata%pos(xdim)
               case ('position_y')
                  tabr1(cgnp) = pset%pdata%pos(ydim)
               case ('position_z')
                  tabr1(cgnp) = pset%pdata%pos(zdim)
               case ('velocity_x')
                  tabr1(cgnp) = pset%pdata%vel(xdim)
               case ('velocity_y')
                  tabr1(cgnp) = pset%pdata%vel(ydim)
               case ('velocity_z')
                  tabr1(cgnp) = pset%pdata%vel(zdim)
               case ('acceleration_x')
                  tabr1(cgnp) = pset%pdata%acc(xdim)
               case ('acceleration_y')
                  tabr1(cgnp) = pset%pdata%acc(ydim)
               case ('acceleration_z')
                  tabr1(cgnp) = pset%pdata%acc(zdim)
               case ('formation_time')
                  tabr1(cgnp) = pset%pdata%tform
               case ('dynamical_time')
                  tabr1(cgnp) = pset%pdata%tdyn
               case default
            end select
         endif
         pset => pset%nxt
      enddo

   end subroutine collect_rank1

   subroutine serial_write_intr1(group_id, ncg, ivar, pvar, cg_src_ncg, proc_ncg, ptag)

      use common_hdf5,    only: get_nth_cg
      use constants,      only: I_ONE
      use dataio_pub,     only: can_i_write, die
      use grid_cont,      only: grid_container
      use hdf5,           only: HID_T
      use MPIF,           only: MPI_INTEGER, MPI_STATUS_IGNORE, MPI_COMM_WORLD
      use MPIFUN,         only: MPI_Recv, MPI_Send
      use mpisetup,       only: master, FIRST, proc, err_mpi
      use particle_utils, only: count_cg_particles

      implicit none

      integer(HID_T), dimension(:,:), intent(in) :: group_id       !< group identifiers
      character(len=*),               intent(in) :: pvar
      integer(kind=4),                intent(in) :: ncg, ivar, cg_src_ncg, proc_ncg, ptag

      type(grid_container), pointer              :: cg
      integer(kind=4)                            :: n_part
      integer(kind=4), dimension(:), allocatable :: tabi

      if (proc_ncg == proc) then
         cg => get_nth_cg(cg_src_ncg)
         n_part = count_cg_particles(cg)
         allocate(tabi(n_part))
         if (n_part > 0) call collect_intr1(pvar, cg, tabi)
      endif

      if (master) then
         if (.not. can_i_write) call die("[cg_particles_io] Master can't write")
         if (proc_ncg /= proc) then
            call MPI_Recv(n_part, I_ONE, MPI_INTEGER, proc_ncg, ptag-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
            allocate(tabi(n_part))
            call MPI_Recv(tabi, n_part, MPI_INTEGER, proc_ncg, ptag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
         endif
         if (n_part > 0) call write_nbody_h5_int_rank1(group_id(ncg, ivar), tabi)
         deallocate(tabi)
      else
         if (can_i_write) call die("[cg_particles_io] Slave can write")
         if (proc_ncg == proc) then
            call MPI_Send(n_part, I_ONE, MPI_INTEGER, FIRST, ptag-1, MPI_COMM_WORLD, err_mpi)
            call MPI_Send(tabi, n_part, MPI_INTEGER, FIRST, ptag, MPI_COMM_WORLD, err_mpi)
            deallocate(tabi)
         endif
      endif

   end subroutine serial_write_intr1

   subroutine serial_write_rank1(group_id, ncg, ivar, pvar, cg_src_ncg, proc_ncg, ptag)

      use common_hdf5,    only: get_nth_cg
      use constants,      only: I_ONE
      use dataio_pub,     only: can_i_write, die
      use grid_cont,      only: grid_container
      use hdf5,           only: HID_T
      use MPIF,           only: MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_STATUS_IGNORE, MPI_COMM_WORLD
      use MPIFUN,         only: MPI_Recv, MPI_Send
      use mpisetup,       only: master, FIRST, proc, err_mpi
      use particle_utils, only: count_cg_particles

      implicit none

      integer(HID_T), dimension(:,:), intent(in) :: group_id       !< group identifiers
      character(len=*),               intent(in) :: pvar
      integer(kind=4),                intent(in) :: ncg, ivar, cg_src_ncg, proc_ncg, ptag

      type(grid_container), pointer              :: cg
      integer(kind=4)                            :: n_part
      real, dimension(:), allocatable            :: tabr

      if (proc_ncg == proc) then
         cg => get_nth_cg(cg_src_ncg)
         n_part = count_cg_particles(cg)
         allocate(tabr(n_part))
         if (n_part > 0) call collect_rank1(pvar, cg, tabr)
      endif

      if (master) then
         if (.not. can_i_write) call die("[cg_particles_io] Master can't write")
         if (proc_ncg /= proc) then
            call MPI_Recv(n_part, I_ONE, MPI_INTEGER, proc_ncg, ptag-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
            allocate(tabr(n_part))
            call MPI_Recv(tabr, n_part, MPI_DOUBLE_PRECISION, proc_ncg, ptag, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
         endif
         if (n_part > 0) call write_nbody_h5_rank1(group_id(ncg, ivar), tabr)
         deallocate(tabr)
      else
         if (can_i_write) call die("[cg_particles_io] Slave can write")
         if (proc_ncg == proc) then
            call MPI_Send(n_part, I_ONE, MPI_INTEGER, FIRST, ptag-1, MPI_COMM_WORLD, err_mpi)
            call MPI_Send(tabr, n_part, MPI_DOUBLE_PRECISION, FIRST, ptag, MPI_COMM_WORLD, err_mpi)
            deallocate(tabr)
         endif
      endif

   end subroutine serial_write_rank1

   subroutine parallel_write_intr1(group_id, pvar, n_part, cg)

      use grid_cont, only: grid_container
      use hdf5,      only: HID_T

      implicit none

      integer(HID_T),                intent(in)    :: group_id       !< File identifier
      character(len=*),              intent(in)    :: pvar
      integer(kind=4),               intent(in)    :: n_part
      type(grid_container), pointer, intent(inout) :: cg

      integer(kind=4), dimension(:), allocatable   :: tabi

      allocate(tabi(n_part))
      call collect_intr1(pvar, cg, tabi)
      call write_nbody_h5_int_rank1(group_id, tabi)
      deallocate(tabi)

   end subroutine parallel_write_intr1

   subroutine parallel_write_rank1(group_id, pvar, n_part, cg)

      use grid_cont, only: grid_container
      use hdf5,      only: HID_T

      implicit none

      integer(HID_T),                intent(in)    :: group_id       !< File identifier
      character(len=*),              intent(in)    :: pvar
      integer(kind=4),               intent(in)    :: n_part
      type(grid_container), pointer, intent(inout) :: cg

      real, dimension(:), allocatable              :: tabr

      allocate(tabr(n_part))
      call collect_rank1(pvar, cg, tabr)
      call write_nbody_h5_rank1(group_id, tabr)
      deallocate(tabr)

   end subroutine parallel_write_rank1

   subroutine write_nbody_h5_int_rank1(dataset_id, tab)

      use dataio_pub, only: die
      use hdf5,       only: h5dwrite_f, HID_T, HSIZE_T, H5T_NATIVE_INTEGER

      implicit none

      integer(HID_T),                intent(in) :: dataset_id
      integer(kind=4), dimension(:), intent(in) :: tab
      integer(HSIZE_T), dimension(1)            :: dimm
      integer(kind=4)                           :: error

      if (all(kind(dataset_id) /= [4, 8])) call die("[cg_particles_io:write_nbody_h5_int_rank1] HID_T doesn't fit to MPI_INTEGER8")

      dimm = shape(tab)
      call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, tab, dimm, error)  ! beware: 64-bit tab(:) produces "no specific subroutine for the generic ‘h5dwrite_f'" error

   end subroutine write_nbody_h5_int_rank1

   subroutine write_nbody_h5_rank1(dataset_id, tab)

      use dataio_pub, only: die
      use hdf5,       only: h5dwrite_f, HID_T, HSIZE_T, H5T_NATIVE_DOUBLE

      implicit none

      integer(HID_T),     intent(in) :: dataset_id
      real, dimension(:), intent(in) :: tab
      integer(HSIZE_T), dimension(1) :: dimm
      integer(kind=4)                :: error

      if (all(kind(dataset_id) /= [4, 8])) call die("[cg_particles_io:write_nbody_h5_rank1] HID_T doesn't fit to MPI_INTEGER8")

      dimm = shape(tab)
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, tab, dimm, error)

   end subroutine write_nbody_h5_rank1

end module cg_particles_io

