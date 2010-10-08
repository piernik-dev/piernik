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

#include "piernik.def"

module multigridio

   implicit none
   private
#ifdef NEW_HDF5
   public :: multigrid_add_hdf5
#endif /* NEW_HDF5 */
   public :: multigrid_write_hdf5

contains

!!$ ============================================================================
!!
!! I/O
!!

#ifdef NEW_HDF5
   subroutine multigrid_add_hdf5

      use multigridvars, only: level_min, level_max, lvl, source, solution, hdf5levels, &
           &                   defect, correction, XDIR, YDIR, ZDIR, NDIM, XLO, XHI, YLO, YHI, ZLO, ZHI
      use list_hdf5,     only: add_lhdf5, lhdf5_info

      implicit none

      type(lhdf5_info) :: item
      integer :: i

      if (.not. hdf5levels) return

      item%p    => get_lvl

      do i = level_min, level_max
         item%sz   = [ lvl(i)%nxb, lvl(i)%nyb, lvl(i)%nzb ]
         item%ind  = [ lvl(i)%is,  lvl(i)%ie,  lvl(i)%js, lvl(i)%je, lvl(i)%ks, lvl(i)%ke ]

         write(item%key, '(A,I1)')  "mg_src", i
         item%opti = [ i, source, 0, 0, 0 ]       !< optional integer passed to func
         call add_lhdf5(item)

         write(item%key, '(A,I1)')  "mg_sln", i
         item%opti = [ i, solution, 0, 0, 0 ]     !< optional integer passed to func
         call add_lhdf5(item)

         write(item%key, '(A,I1)')  "mg_def", i
         item%opti = [ i, defect, 0, 0, 0 ]       !< optional integer passed to func
         call add_lhdf5(item)

         write(item%key, '(A,I1)')  "mg_cor", i
         item%opti = [ i, correction, 0, 0, 0 ]   !< optional integer passed to func
         call add_lhdf5(item)
      enddo

      contains
         function get_lvl(i, sz, opti) result (outtab)
            implicit none
            integer, dimension(XLO:ZHI), intent(in)     :: i
            integer, dimension(NDIM), intent(in)        :: sz
            integer, dimension(5), intent(in)           :: opti
            real, dimension(sz(XDIR), sz(YDIR), sz(ZDIR)) :: outtab
            outtab(:,:,:) =  lvl(opti(1))%mgvar(i(XLO):i(XHI), i(YLO):i(YHI), i(ZLO):i(ZHI), opti(2))
         end function get_lvl

   end subroutine multigrid_add_hdf5
#endif /* NEW_HDF5 */

!!$ ============================================================================

   subroutine multigrid_write_hdf5(file_id)

      use multigridvars,   only: level_min, level_max, lvl
      use hdf5,            only: HID_T

      implicit none

      integer(HID_T), intent(in) :: file_id

      character(len=8) :: dname
      integer          :: i

      do i = level_min, level_max
         write(dname, '(A3,I1)') 'lvl', i
         call multigrid_write_arr(lvl(i), trim(dname), file_id)
      enddo

   end subroutine multigrid_write_hdf5

!!$ ============================================================================

   subroutine multigrid_write_arr(llvl, dsetname, file_id)

      use mpisetup, only: pcoords, psize
      use hdf5,     only: HID_T, HSIZE_T, HSSIZE_T, &
           &              H5FD_MPIO_INDEPENDENT_F, H5P_DATASET_CREATE_F, H5P_DATASET_XFER_F, H5S_SELECT_SET_F, H5T_NATIVE_DOUBLE, &
           &              h5screate_simple_f, h5sclose_f, h5sselect_hyperslab_f, &
           &              h5pcreate_f, h5pset_chunk_f, h5pset_dxpl_mpio_f, h5pset_dxpl_mpio_f, &
           &              h5dcreate_f, h5dget_space_f, h5dwrite_f, h5dclose_f
      use multigridvars, only : plvl, ngridvars, source, solution, defect, &
           &                    correction, hdf5levels, XDIR, YDIR, ZDIR, NDIM

      implicit none

      type(plvl), intent(in) :: llvl

      integer :: rank = NDIM ! Dataset rank

      CHARACTER(LEN=*) :: dsetname    !< Dataset name

      integer(HID_T) :: file_id       !< File identifier
      integer(HID_T) :: dset_id       !< Dataset identifier
      integer(HID_T) :: plist_id      !< Property list identifier
      integer(HID_T) :: filespace     !< Dataspace identifier in file
      integer(HID_T) :: memspace      !< Dataspace identifier in memory

      integer(HSIZE_T),  DIMENSION(NDIM) :: count
      integer(HSSIZE_T), DIMENSION(NDIM) :: offset
      integer(HSIZE_T),  DIMENSION(NDIM) :: stride
      integer(HSIZE_T),  DIMENSION(NDIM) :: block
      integer(HSIZE_T),  DIMENSION(NDIM) :: dimsf, dimsfi, chunk_dims
      integer :: error

      integer :: nxd, nyd, nzd
      integer :: v
      integer, parameter :: dnamelen = 32
      character(len=dnamelen) :: dname

      if (.not. hdf5levels) return

      nxd = llvl%nxb*psize(XDIR)
      nyd = llvl%nyb*psize(YDIR)
      nzd = llvl%nzb*psize(ZDIR)

      dimsf = [ nxd, nyd, nzd ] ! Dataset dimensions
      dimsfi = dimsf
      chunk_dims = [ llvl%nxb, llvl%nyb, llvl%nzb ] ! Chunks dimensions

      do v = 1, ngridvars
         select case(v)
            case(source)
               write(dname,'(2a)') trim(dsetname),"_src"
            case(solution)
               write(dname,'(2a)') trim(dsetname),"_soln"
            case(defect)
               write(dname,'(2a)') trim(dsetname),"_def"
            case(correction)
               write(dname,'(2a)') trim(dsetname),"_corr"
            case default
               write(dname,'(2a)') trim(dsetname),"_???"
         end select
         !
         ! Create the data space for the  dataset.
         !
         CALL h5screate_simple_f(rank, dimsf, filespace, error)
         CALL h5screate_simple_f(rank, chunk_dims, memspace, error)

         !
         ! Create chunked dataset.
         !
         CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
         CALL h5pset_chunk_f(plist_id, rank, chunk_dims, error)
         CALL h5dcreate_f(file_id, dname, H5T_NATIVE_DOUBLE, filespace, dset_id, error, plist_id)
         CALL h5sclose_f(filespace, error)

         !
         ! Each process defines dataset in memory and writes it to the hyperslab
         ! in the file.
         !
         stride(:) = 1
         count(:) =  1
         block(:) = chunk_dims(:)

         offset(:) = pcoords(:)*chunk_dims(:)
         !
         ! Select hyperslab in the file.
         !
         CALL h5dget_space_f(dset_id, filespace, error)
         CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error, stride, block)
         !
         ! Create property list for collective dataset write
         !
         CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
         CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)

         !
         ! Write the dataset collectively.
         !

         CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, llvl%mgvar(llvl%is:llvl%ie, llvl%js:llvl%je, llvl%ks:llvl%ke, v), dimsfi, error, &
              &          file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

         !
         ! Close dataspaces.
         !
         CALL h5sclose_f(filespace, error)
         CALL h5sclose_f(memspace, error)
         !
         ! Close the dataset.
         !
         CALL h5dclose_f(dset_id, error)

      end do

   end subroutine multigrid_write_arr

end module multigridio
