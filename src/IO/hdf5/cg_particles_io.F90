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

   character(len=dsetnamelen), dimension(*), parameter  :: pvarn = ['ppid', 'mass', 'ener', 'ppos', 'pvel', 'pacc']
   logical,                    dimension(size(pvarn))   :: pvarl = .false.
   character(len=dsetnamelen), allocatable, dimension(:) ::pdsets

   contains

   subroutine init_nbody_hdf5(pvars)

      use dataio_pub, only: msg, warn

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
            write(msg,'(2a)')'[particles_io_hdf5::init_nbody_hdf5]: unknown particle var: ', pvars(il) ; call warn(msg)
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

      integer(HID_T),   intent(in) :: group_id       !< File identifier
      character(len=*), intent(in) :: pvar
      integer(kind=4),  intent(in) :: n_part

      select case (pvar)
         case ('ppid')
            call collect_and_write_intr1(group_id, pvar, n_part)
         case ('mass', 'ener')
            call collect_and_write_rank1(group_id, pvar, n_part)
         case ('ppos', 'pvel', 'pacc')
            call collect_and_write_rank2(group_id, pvar, n_part)
         case default
      end select

   end subroutine nbody_datafields

   subroutine collect_and_write_intr1(group_id, pvar, n_part)

      use cg_leaves, only: leaves
      use cg_list,   only: cg_list_element
      use hdf5,      only: HID_T

      implicit none

      integer(HID_T),   intent(in)       :: group_id       !< File identifier
      character(len=*), intent(in)       :: pvar
      integer(kind=4),  intent(in)       :: n_part
      integer                            :: cgnp, recnp, i
      integer, dimension(:), allocatable :: tabi1
      type(cg_list_element), pointer     :: cgl

      allocate(tabi1(n_part))
      recnp = 0

      cgl => leaves%first
      do while (associated(cgl))
         cgnp = 0
         select case (pvar)
            case ('ppid')
               do i = 1, size(cgl%cg%pset%p, dim=1)
                  if (cgl%cg%pset%p(i)%phy) then
                     cgnp = cgnp + 1
                     tabi1(recnp+cgnp) = cgl%cg%pset%p(i)%pid
                  endif
               enddo
            case default
         end select
         recnp = recnp+cgnp
         cgl => cgl%nxt
      enddo

      call write_nbody_h5_int_rank1(group_id, pvar, tabi1)
      deallocate(tabi1)

   end subroutine collect_and_write_intr1

   subroutine collect_and_write_rank1(group_id, pvar, n_part)

      use cg_leaves, only: leaves
      use cg_list,   only: cg_list_element
      use hdf5,      only: HID_T

      implicit none

      integer(HID_T),   intent(in)    :: group_id       !< File identifier
      character(len=*), intent(in)    :: pvar
      integer(kind=4),  intent(in)    :: n_part
      integer                         :: cgnp, recnp, i
      real, dimension(:), allocatable :: tabr1
      type(cg_list_element), pointer  :: cgl

      allocate(tabr1(n_part))
      recnp = 0

      cgl => leaves%first
      do while (associated(cgl))
         !cgnp = size(cgl%cg%pset%p, dim=1)
         cgnp = 0
         do i = 1, size(cgl%cg%pset%p, dim=1)
            if (cgl%cg%pset%p(i)%phy) then
               cgnp = cgnp + 1
               select case (pvar)
                  case ('mass')
                     tabr1(recnp+cgnp) = cgl%cg%pset%p(i)%mass
                  case ('ener')
                     tabr1(recnp+cgnp) = cgl%cg%pset%p(i)%energy
                  case default
               end select
            endif
         enddo
         recnp = recnp+cgnp
         cgl => cgl%nxt
      enddo

      call write_nbody_h5_rank1(group_id, pvar, tabr1)
      deallocate(tabr1)

   end subroutine collect_and_write_rank1

   subroutine collect_and_write_rank2(group_id, pvar, n_part)

      use cg_leaves, only: leaves
      use cg_list,   only: cg_list_element
      use constants, only: ndims
      use hdf5,      only: HID_T

      implicit none

      integer(HID_T),   intent(in)      :: group_id       !< File identifier
      character(len=*), intent(in)      :: pvar
      integer(kind=4),  intent(in)      :: n_part
      integer                           :: cgnp, recnp, i, j
      real, dimension(:,:), allocatable :: tabr2
      type(cg_list_element), pointer    :: cgl

      allocate(tabr2(n_part, ndims))

      recnp = 0

      cgl => leaves%first
      do while (associated(cgl))
         cgnp = size(cgl%cg%pset%p, dim=1)
         j=1
         select case (pvar)
            case ('ppos')
               do i = 1, cgnp
                  if (cgl%cg%pset%p(i)%phy) then
                     tabr2(recnp+j,:) = cgl%cg%pset%p(i)%pos(:)
                     j = j + 1
                  endif
               enddo
            case ('pvel')
               do i = 1, cgnp
                  if (cgl%cg%pset%p(i)%phy) then
                     tabr2(recnp+j,:) = cgl%cg%pset%p(i)%vel(:)
                     j = j + 1
                  endif
               enddo
            case ('pacc')
               do i = 1, cgnp
                  if (cgl%cg%pset%p(i)%phy) then
                     tabr2(recnp+j,:) = cgl%cg%pset%p(i)%acc(:)
                     j = j + 1
                  endif
               enddo
            case default
         end select
         cgl => cgl%nxt
      enddo

      call write_nbody_h5_rank2(group_id, pvar, tabr2)

      deallocate(tabr2)

   end subroutine collect_and_write_rank2

   subroutine write_nbody_h5_int_rank1(group_id, vvar, tab)

      use hdf5, only: h5dcreate_f, h5dclose_f, h5dwrite_f, h5screate_simple_f, h5sclose_f, HID_T, HSIZE_T, H5T_NATIVE_INTEGER

      implicit none

      character(len=*),      intent(in) :: vvar
      integer(HID_T),        intent(in) :: group_id
      integer, dimension(:), intent(in) :: tab
      integer(HSIZE_T), dimension(1)    :: dimm
      integer(HID_T)                    :: dataspace_id, dataset_id
      integer(kind=4)                   :: error, rank1 = 1

      dimm = shape(tab)
      dataset_id = group_id !For 1 file writing group_id is the cg_g_id
#ifndef NBODY_1FILE
      call h5screate_simple_f(rank1, dimm, dataspace_id, error)
      call h5dcreate_f(group_id, vvar, H5T_NATIVE_INTEGER, dataspace_id, dataset_id, error)
#endif /* NBODY_1FILE */
      call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, tab, dimm, error)
#ifndef NBODY_1FILE
      call h5dclose_f(dataset_id, error)
      call h5sclose_f(dataspace_id, error)
#endif /* NBODY_1FILE */

   end subroutine write_nbody_h5_int_rank1

   subroutine write_nbody_h5_rank1(group_id, vvar, tab)

      use hdf5, only: h5dcreate_f, h5dclose_f, h5dwrite_f, h5screate_simple_f, h5sclose_f, HID_T, HSIZE_T, H5T_NATIVE_DOUBLE

      implicit none

      character(len=*),   intent(in) :: vvar
      integer(HID_T),     intent(in) :: group_id
      real, dimension(:), intent(in) :: tab
      integer(HSIZE_T), dimension(1) :: dimm
      integer(HID_T)                 :: dataspace_id, dataset_id
      integer(kind=4)                :: error, rank1 = 1

      dimm = shape(tab)
      dataset_id = group_id
#ifndef NBODY_1FILE
      call h5screate_simple_f(rank1, dimm, dataspace_id, error)
      call h5dcreate_f(group_id, vvar, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error)
#endif /* NBODY_1FILE */
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, tab, dimm, error)
#ifndef NBODY_1FILE
      call h5dclose_f(dataset_id, error)
      call h5sclose_f(dataspace_id, error)
#endif /* NBODY_1FILE */

   end subroutine write_nbody_h5_rank1

   subroutine write_nbody_h5_rank2(group_id, vvar, tab)

      use hdf5, only: h5dcreate_f, h5dclose_f, h5dwrite_f, h5screate_simple_f, h5sclose_f, HID_T, HSIZE_T, H5T_NATIVE_DOUBLE

      implicit none

      character(len=*),     intent(in) :: vvar
      integer(HID_T),       intent(in) :: group_id
      real, dimension(:,:), intent(in) :: tab
      integer(HSIZE_T), dimension(2)   :: dimv
      integer(HID_T)                   :: dataspace_id, dataset_id
      integer(kind=4)                  :: error, rank2 = 2

      dimv = shape(tab)
      dataset_id = group_id
#ifndef NBODY_1FILE
      call h5screate_simple_f(rank2, dimv, dataspace_id, error)
      call h5dcreate_f(group_id, vvar, H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error)
#endif /* NBODY_1FILE */
      call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, tab, dimv, error)
#ifndef NBODY_1FILE
      call h5dclose_f(dataset_id, error)
      call h5sclose_f(dataspace_id, error)
#endif /* NBODY_1FILE */

    end subroutine write_nbody_h5_rank2

end module cg_particles_io

