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
  public :: dump_cg_particles, init_nbody_hdf5, pdsets, nbody_datafields

   character(len=dsetnamelen), dimension(*), parameter  :: pvarn = ['ppid', 'mass', 'ener', 'posx', 'posy', 'posz', 'velx', 'vely', 'velz', 'accx', 'accy', 'accz']
   logical,                    dimension(size(pvarn))   :: pvarl = .false.
   character(len=dsetnamelen), allocatable, dimension(:) ::pdsets

   contains

   subroutine init_nbody_hdf5(pvars)

      use dataio_pub, only: msg, warn
      use mpisetup,   only: master

      implicit none

      character(len=dsetnamelen), dimension(:), intent(in) :: pvars  !< quantities to be plotted, see dataio::vars
      integer                                              :: ie, il, k, l
      logical                                              :: var_found

      k=0
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
         if (.not.var_found) then
            write(msg,'(2a)')'[particles_io_hdf5::init_nbody_hdf5]: unknown particle var: ', pvars(il)
            if (master) call warn(msg)
         endif
      enddo

#ifdef NBODY_1FILE
      allocate(pdsets(k))
      l=1
      do il = lbound(pvars, 1), ubound(pvars, 1)
         if (len(trim(pvars(il))) == 0) cycle
         if (pvarl(il) .eqv. .true.) then
            pdsets(l)=pvars(il)
            l=l+1
         endif
      enddo
#endif /* NBODY_1FILE */

    end subroutine init_nbody_hdf5

  subroutine dump_cg_particles(group_id)

    use particle_utils, only: count_all_particles
    use hdf5, only: HID_T

    implicit none
    integer(HID_T), intent(in)    :: group_id
    integer(kind=4)               :: n_part

    n_part = count_all_particles()
    call nbody_datasets(n_part, group_id)

  end subroutine dump_cg_particles

  subroutine nbody_datasets(n_part, group_id)
     use hdf5, only: HID_T

     implicit none

     integer(HID_T)              :: group_id       !< File identifier
     integer(kind=4), intent(in) :: n_part
     integer                     :: i

     do i = lbound(pvarl, 1), ubound(pvarl, 1)
        if (pvarl(i)) then
           call nbody_datafields(group_id, trim(pvarn(i)), n_part)
        endif
     enddo

   end subroutine nbody_datasets

   subroutine nbody_datafields(group_id, pvar, n_part)

     use hdf5, only: HID_T

      implicit none

      integer(HID_T),    intent(in) :: group_id       !< File identifier
      character(len=*), intent(in) :: pvar
      integer(kind=4),  intent(in) :: n_part

      select case (pvar)
         case ('id')
            call collect_and_write_intr1(group_id, pvar, n_part)
         case default
            call collect_and_write_rank1(group_id, pvar, n_part)
      end select

   end subroutine nbody_datafields

   subroutine collect_and_write_intr1(group_id, pvar, n_part)

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: I_ONE
      use dataio_pub,     only: nproc_io, can_i_write, die
      use domain,         only: is_multicg
      use hdf5,           only: HID_T
      use MPIF,           only: MPI_INTEGER, MPI_STATUS_IGNORE, MPI_DOUBLE_INT, MPI_Recv, MPI_Send
      use mpisetup,       only: master, FIRST, LAST, proc, comm, mpi_err
      use particle_types, only: particle

      implicit none

      integer(HID_T),   intent(in)       :: group_id       !< File identifier
      character(len=*), intent(in)       :: pvar
      integer(kind=4),  intent(in)       :: n_part
      integer                            :: cgnp, recnp
      integer(kind=4)                    :: ncg
      integer(kind=4), dimension(:), allocatable :: tabi1, tabi2
      type(cg_list_element), pointer     :: cgl
      type(particle), pointer        :: pset

      if (is_multicg) call die("[cg_particles_io:collect_and_write_intr1] several cg per processor not implemented yet")

      allocate(tabi1(n_part))
      recnp = 0

      cgl => leaves%first
      do while (associated(cgl))
         cgnp = 0
         select case (pvar)
            case ('id')
               pset => cgl%cg%pset%first
               do while (associated(pset))
                  if (pset%pdata%phy) then
                     cgnp = cgnp + 1
                     tabi1(recnp+cgnp) = pset%pdata%pid
                  endif
                  pset => pset%nxt
               enddo
            case default
         end select
         recnp = recnp+cgnp
         cgl => cgl%nxt
      enddo

      ! Not compatible with AMR or several cg per processor
      if (nproc_io == 1) then !perform serial write
         ! write all cg, one by one
         do ncg = FIRST, LAST
            if (master) then
               if (.not. can_i_write) call die("[cg_particles_io] Master can't write")
               if (ncg == proc) then
                  call write_nbody_h5_int_rank1(group_id, pvar, tabi1)
               else
                  call MPI_Recv(n_part, I_ONE, MPI_INTEGER, ncg, ncg, comm, MPI_STATUS_IGNORE, mpi_err)
                  allocate(tabi2(n_part))
                  call MPI_Recv(tabi2, n_part, MPI_INTEGER, ncg, ncg, comm, MPI_STATUS_IGNORE, mpi_err)
                  call MPI_Recv(group_id, I_ONE, MPI_DOUBLE_INT, ncg, ncg, comm, MPI_STATUS_IGNORE, mpi_err)
                  call write_nbody_h5_int_rank1(group_id, pvar, tabi2)
                  deallocate(tabi2)
               endif
            else
               if (can_i_write) call die("[cg_particles_io] Slave can write")
               if (ncg == proc) then
                  call MPI_Send(n_part, I_ONE, MPI_INTEGER, FIRST, ncg, comm, mpi_err)
                  call MPI_Send(tabi1, n_part, MPI_INTEGER, FIRST, ncg, comm, mpi_err)
                  call MPI_Send(group_id, I_ONE, MPI_DOUBLE_INT, FIRST, ncg, comm, mpi_err)
               endif
            endif
         enddo
      else
         call write_nbody_h5_int_rank1(group_id, pvar, tabi1)
      endif

      deallocate(tabi1)

   end subroutine collect_and_write_intr1

   subroutine collect_and_write_rank1(group_id, pvar, n_part)

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: xdim, ydim, zdim, I_ONE
      use dataio_pub,     only: nproc_io, can_i_write, die
      use domain,         only: is_multicg
      use hdf5,           only: HID_T
      use MPIF,           only: MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_STATUS_IGNORE, MPI_DOUBLE_INT, MPI_Recv, MPI_Send
      use mpisetup,       only: master, FIRST, LAST, proc, comm, mpi_err
      use particle_types, only: particle

      implicit none

      integer(HID_T),   intent(in)    :: group_id       !< File identifier
      character(len=*), intent(in)    :: pvar
      integer(kind=4),  intent(in)    :: n_part

      integer                         :: cgnp, recnp, i
      integer(kind=4)                 :: ncg
      real, dimension(:), allocatable :: tabr1, tabr2
      type(cg_list_element), pointer  :: cgl
      type(particle), pointer         :: pset

      if (is_multicg) call die("[cg_particles_io:collect_and_write_rank1] several cg per processor not implemented yet")

      allocate(tabr1(n_part))
      recnp = 0

      cgl => leaves%first
      do while (associated(cgl))
         cgnp = 0
         pset => cgl%cg%pset%first
         do while (associated(pset))
            if (pset%pdata%phy) then
               cgnp = cgnp + 1
               select case (pvar)
                  case ('mass')
                     tabr1(recnp+cgnp) = pset%pdata%mass
                     i=2
                  case ('energy')
                     tabr1(recnp+cgnp) = pset%pdata%energy
                     i=3
                  case ('position_x')
                     tabr1(recnp+cgnp) = pset%pdata%pos(xdim)
                     i=4
                  case ('position_y')
                     tabr1(recnp+cgnp) = pset%pdata%pos(ydim)
                     i=5
                  case ('position_z')
                     tabr1(recnp+cgnp) = pset%pdata%pos(zdim)
                     i=6
                  case ('velocity_x')
                     tabr1(recnp+cgnp) = pset%pdata%vel(xdim)
                     i=7
                  case ('velocity_y')
                     tabr1(recnp+cgnp) = pset%pdata%vel(ydim)
                     i=8
                  case ('velocity_z')
                     tabr1(recnp+cgnp) = pset%pdata%vel(zdim)
                     i=9
                  case ('acceleration_x')
                     tabr1(recnp+cgnp) = pset%pdata%acc(xdim)
                     i=10
                  case ('acceleration_y')
                     tabr1(recnp+cgnp) = pset%pdata%acc(ydim)
                     i=11
                  case ('acceleration_z')
                     tabr1(recnp+cgnp) = pset%pdata%acc(zdim)
                     i=12
                  case default
               end select
            endif
            pset => pset%nxt
         enddo
         recnp = recnp+cgnp
         cgl => cgl%nxt
      enddo

      ! Not compatible with AMR or several cg per processor
      if (nproc_io == 1) then !perform serial write
         ! write all cg, one by one
         do ncg = FIRST, LAST
            if (master) then
               if (.not. can_i_write) call die("[cg_particles_io] Master can't write")
               if (ncg == proc) then
                  call write_nbody_h5_rank1(group_id, pvar, tabr1)
               else
                  call MPI_Recv(n_part, I_ONE, MPI_INTEGER, ncg, ncg, comm, MPI_STATUS_IGNORE, mpi_err)
                  allocate(tabr2(n_part))
                  call MPI_Recv(tabr2, n_part, MPI_DOUBLE_PRECISION, ncg, ncg, comm, MPI_STATUS_IGNORE, mpi_err)
                  call MPI_Recv(group_id, I_ONE, MPI_DOUBLE_INT, ncg, ncg, comm, MPI_STATUS_IGNORE, mpi_err)
                  call write_nbody_h5_rank1(group_id, pvar, tabr2)
                  deallocate(tabr2)
               endif
            else
               if (can_i_write) call die("[cg_particles_io] Slave can write")
               if (ncg == proc) then
                  call MPI_Send(n_part, I_ONE, MPI_INTEGER, FIRST, ncg, comm, mpi_err)
                  call MPI_Send(tabr1, n_part, MPI_DOUBLE_PRECISION, FIRST, ncg, comm, mpi_err)
                  call MPI_Send(group_id, I_ONE, MPI_DOUBLE_INT, FIRST, ncg, comm, mpi_err)
               endif
            endif
         enddo
      else ! perform parallel write
         call write_nbody_h5_rank1(group_id, pvar, tabr1)
      endif

      deallocate(tabr1)

   end subroutine collect_and_write_rank1

   subroutine write_nbody_h5_int_rank1(group_id, vvar, tab)

      use hdf5, only: h5dcreate_f, h5dclose_f, h5dwrite_f, h5screate_simple_f, h5sclose_f, HID_T, HSIZE_T, H5T_NATIVE_INTEGER

      implicit none

      character(len=*),      intent(in) :: vvar
      integer(HID_T),        intent(in) :: group_id
      integer(kind=4), dimension(:), intent(in) :: tab
      integer(HSIZE_T), dimension(1)    :: dimm
      integer(kind=4)                   :: error
      integer(HID_T)                    :: dataset_id
#ifndef NBODY_1FILE
      integer(HID_T)                    :: dataspace_id
      integer(kind=4)                   :: rank1 = 1
#endif /* !NBODY_1FILE */

      dimm = shape(tab)
      dataset_id = group_id !For 1 file writing group_id is the cg_g_id
#ifndef NBODY_1FILE
      call h5screate_simple_f(rank1, dimm, dataspace_id, error)
      call h5dcreate_f(group_id, vvar, H5T_NATIVE_INTEGER, dataspace_id, dataset_id, error)
#endif /* !NBODY_1FILE */
      call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, tab, dimm, error)  ! beware: 64-bit tab(:) produces "no specific subroutine for the generic ‘h5dwrite_f'" error
#ifndef NBODY_1FILE
      call h5dclose_f(dataset_id, error)
      call h5sclose_f(dataspace_id, error)
#endif /* !NBODY_1FILE */

#ifdef NBODY_1FILE
      if (.false.) error = len(vvar, kind=4)  ! suppress -Wunused-dummy-argument
#endif /* NBODY_1FILE */

   end subroutine write_nbody_h5_int_rank1

   subroutine write_nbody_h5_rank1(group_id, vvar, tab)

     use hdf5, only: h5dcreate_f, h5dclose_f, h5dwrite_f, h5screate_simple_f, h5sclose_f, HID_T, HSIZE_T, H5T_NATIVE_DOUBLE

      implicit none

      character(len=*),   intent(in) :: vvar
      integer(HID_T),     intent(in) :: group_id
      real, dimension(:), intent(in) :: tab
      integer(HSIZE_T), dimension(1) :: dimm
      integer(HID_T)                 :: dataset_id
      integer(kind=4)                :: error
#ifndef NBODY_1FILE
      integer(HID_T)                 :: dataspace_id
      integer(kind=4)                :: rank1 = 1
#endif /* !NBODY_1FILE */

      dimm = shape(tab)
      dataset_id = group_id
#ifndef NBODY_1FILE
      call h5screate_simple_f(rank1, dimm, dataspace_id, error)
      call h5dcreate_f(group_id, vvar, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error)
#endif /* !NBODY_1FILE */
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, tab, dimm, error)
#ifndef NBODY_1FILE
      call h5dclose_f(dataset_id, error)
      call h5sclose_f(dataspace_id, error)
#endif /* !NBODY_1FILE */

#ifdef NBODY_1FILE
      if (.false.) error = len(vvar, kind=4)  ! suppress -Wunused-dummy-argument
#endif /* NBODY_1FILE */

   end subroutine write_nbody_h5_rank1


end module cg_particles_io

