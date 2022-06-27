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

!>
!! \brief This module contains I/O instructions and handlers for CRESP (NR) data.
!<
module cresp_io
! pulled by CRESP

   implicit none

   private
   public   :: save_smap_to_open, save_cresp_smap_h5, write_cresp_to_restart, cresp_read_smaps_from_hdf, &
            &  read_cresp_smap_fields, read_real_arr2d_dset, read_smap_header_h5, &
            &  save_NR_smap, check_NR_smaps_headers, read_NR_smap, read_NR_smap_header, got_smaps_from_restart

   logical :: got_smaps_from_restart = .false.

contains

   subroutine write_cresp_to_restart(file_id)

      use hdf5,           only: HID_T, SIZE_T
      use h5lt,           only: h5ltset_attribute_int_f, h5ltset_attribute_double_f
      use initcosmicrays, only: ncrb
      use initcrspectrum, only: e_small, p_min_fix, p_max_fix, q_big, use_cresp

      implicit none

      integer(HID_T), intent(in) :: file_id
      integer(SIZE_T)            :: bufsize
      integer(kind=4)            :: error
      integer(kind=4), dimension(1) :: lnsnbuf_i
      real,    dimension(1)      :: lnsnbuf_r

      bufsize = 1
      lnsnbuf_i(bufsize) = ncrb
      call h5ltset_attribute_int_f(file_id,    "/", "ncrb",      lnsnbuf_i, bufsize, error)

      lnsnbuf_r(bufsize) = p_min_fix
      call h5ltset_attribute_double_f(file_id, "/", "p_min_fix", lnsnbuf_r, bufsize, error)

      lnsnbuf_r(bufsize) = p_max_fix
      call h5ltset_attribute_double_f(file_id, "/", "p_max_fix", lnsnbuf_r, bufsize, error)

      lnsnbuf_r(bufsize) = q_big
      call h5ltset_attribute_double_f(file_id, "/", "q_big",     lnsnbuf_r, bufsize, error)

      lnsnbuf_r(bufsize) = e_small
      call h5ltset_attribute_double_f(file_id, "/", "e_small",   lnsnbuf_r, bufsize, error)

      if (use_cresp) then
         call create_cresp_smap_fields(file_id) ! create "/cresp/smaps_{LO,UP}/..."
         call cresp_write_smaps_to_hdf(file_id) ! create "/cresp/smaps_{LO,UP}/{p_f}_ratio"
      endif

   end subroutine write_cresp_to_restart

!----------------------------------------------------------------------------------------------------
   subroutine cresp_write_smaps_to_hdf(file_id)

      use constants,     only: LO, HI
      use cresp_helpers, only: n_g_smaps, dset_attrs, p_ratios_lo, f_ratios_lo, p_ratios_up, f_ratios_up
      use hdf5,          only: HID_T

      implicit none

      integer(HID_T), intent(in) :: file_id

      call save_smap_to_open(file_id, n_g_smaps(LO), dset_attrs(1), p_ratios_lo)
      call save_smap_to_open(file_id, n_g_smaps(LO), dset_attrs(2), f_ratios_lo)
      call save_smap_to_open(file_id, n_g_smaps(HI), dset_attrs(1), p_ratios_up)
      call save_smap_to_open(file_id, n_g_smaps(HI), dset_attrs(2), f_ratios_up)

   end subroutine cresp_write_smaps_to_hdf
!----------------------------------------------------------------------------------------------------
!> \brief Create fields containing parameters of cresp solution maps, which are saved later in these groups.
!
   subroutine create_cresp_smap_fields(file_id)

      use constants,       only: LO, HI
      use cresp_helpers,   only: n_g_cresp, n_g_smaps, n_a_dims, n_a_esmall, n_a_max_p_r, &
      &  n_a_clight, n_a_qbig, n_a_qbig, n_a_amin, n_a_amax, n_a_nmin, n_a_nmax, hdr_io
      use hdf5,            only: HID_T, h5gcreate_f, h5gclose_f
      use helpers_hdf5,    only: create_attribute

      implicit none

      integer(HID_T), intent(in) :: file_id
      integer(HID_T)             :: g_id, g_id_cresp
      integer(kind=4)            :: error, i

      call h5gcreate_f(file_id, n_g_cresp, g_id_cresp, error)
         do i = LO, HI
            call h5gcreate_f(file_id, n_g_smaps(i), g_id, error)
               call create_attribute(g_id, n_a_dims,    [ hdr_io(i)%s_dim1, hdr_io(i)%s_dim2 ] )
               call create_attribute(g_id, n_a_esmall,  [ hdr_io(i)%s_es   ])
               call create_attribute(g_id, n_a_max_p_r, [ hdr_io(i)%s_pr   ])
               call create_attribute(g_id, n_a_qbig,    [ hdr_io(i)%s_qbig ])
               call create_attribute(g_id, n_a_clight,  [ hdr_io(i)%s_c    ])
               call create_attribute(g_id, n_a_amin,    [ hdr_io(i)%s_amin ])
               call create_attribute(g_id, n_a_amax,    [ hdr_io(i)%s_amax ])
               call create_attribute(g_id, n_a_nmin,    [ hdr_io(i)%s_nmin ])
               call create_attribute(g_id, n_a_nmax,    [ hdr_io(i)%s_nmax ])
            call h5gclose_f(g_id, error)
         enddo
      call h5gclose_f(g_id_cresp, error)

   end subroutine create_cresp_smap_fields
!---------------------------------------------------------------------------------------------------
   subroutine cresp_read_smaps_from_hdf(file_id)

      use constants,     only: LO, HI, I_ZERO
      use cresp_helpers, only: map_header, dset_attrs, n_g_smaps, p_ratios_lo, f_ratios_lo, p_ratios_up, f_ratios_up, hdr_res, allocate_smaps, deallocate_smaps
      use dataio_pub,    only: warn
      use hdf5,          only: HID_T

      implicit none

      integer(HID_T), intent(in) :: file_id

      call read_smap_header_h5(file_id, hdr_res)

      if (hdr_res(1)%s_dim1 .eq. I_ZERO .or. hdr_res(1)%s_dim2 .eq. I_ZERO) then
         call warn("[cresp_NR_method:cresp_read_smaps_from_hdf] Got solution map dimension = 0. Will solve for new maps.")
      else
         call deallocate_smaps ! TODO just in case. Reading should be called before "fill_guess_grids"

         call allocate_smaps(hdr_res(1)%s_dim1, hdr_res(1)%s_dim2) ! TODO decide whether the same dim is forced onto all maps (rather so)

         call read_real_arr2d_dset(file_id, n_g_smaps(LO)//"/"//dset_attrs(1), p_ratios_lo)
         call read_real_arr2d_dset(file_id, n_g_smaps(LO)//"/"//dset_attrs(2), f_ratios_lo)
         call read_real_arr2d_dset(file_id, n_g_smaps(HI)//"/"//dset_attrs(1), p_ratios_up)
         call read_real_arr2d_dset(file_id, n_g_smaps(HI)//"/"//dset_attrs(2), f_ratios_up)
         got_smaps_from_restart = .true.
      endif

   end subroutine cresp_read_smaps_from_hdf
!----------------------------------------------------------------------------------------------------
!> \brief Create an external file to save ONLY cresp solution maps, optional, WIP TODO, FIXME
!
   subroutine save_cresp_smap_h5(smap_data, bound, dsname, filename)

      use constants,     only: I_ZERO, I_ONE
      use cresp_helpers, only: n_g_smaps, cresp_gname
      use dataio_pub,    only: msg, printinfo
      use hdf5,          only: h5open_f, h5close_f, h5fcreate_f, h5fopen_f, h5fclose_f, h5dopen_f, h5dclose_f, h5eset_auto_f, h5gopen_f, h5gclose_f, &
           &                   HID_T, H5F_ACC_RDWR_F, H5F_ACC_TRUNC_F, H5F_ACC_RDONLY_F

      implicit none

      real, dimension(:,:), target, intent(inout) :: smap_data
      integer,                      intent(in)    :: bound
      character(len=*),             intent(in)    :: dsname
      character(len=*),             intent(in)    :: filename
      logical                                     :: file_exist, has_group, has_dset
      integer(HID_T)                              :: file_id, dset_id, group_id
      integer(kind=4)                             :: error

      file_exist = .false.
      has_group  = .false.

      inquire(file = trim(filename), exist = file_exist)    ! check if file "filename" exists

      if (file_exist) then
         has_group = .false.

         call h5open_f(error)
         call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)
         call h5eset_auto_f(I_ZERO, error)   ! Turn off printing messages
         call h5gopen_f(file_id, cresp_gname, group_id, error)
         if (error .eq. 0) has_group = .true.
         call h5eset_auto_f(I_ONE, error)    ! Turn on  printing messages
         call h5gclose_f(group_id, error)
         call h5fclose_f(file_id, error)
         call h5close_f(error)
      else
         write(msg,"(A,A,A)") "[cresp_io:save_cresp_smap_h5] File '", trim(filename),"' does not exist; must be created and initialized."
         call printinfo(msg)

         call h5open_f(error)
         call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, error)
         call h5fclose_f(file_id, error)
         call h5close_f(error)
      endif

      if (.not. has_group) then
         call h5open_f(error)
         call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, error)
         call create_cresp_smap_fields(file_id)   ! Creates "/cresp" group and associated attrs
         call h5fclose_f(file_id, error)
         call h5close_f(error)
      endif

      has_dset = .false.

      call h5open_f(error)
      call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)
      call h5eset_auto_f(I_ZERO, error)   ! Turn off printing messages
      call h5dopen_f(file_id, "/"//n_g_smaps(bound)//"/"//dsname, dset_id, error)
      if (error .eq. 0) has_dset = .true.
      call h5dclose_f(dset_id, error)
      call h5eset_auto_f(I_ONE, error)    ! Turn on  printing messages
      call h5fclose_f(file_id, error)
      call h5close_f(error)

      if (.not. has_dset) then
         call h5open_f(error)
         call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, error)
         call save_smap_to_open(file_id, n_g_smaps(bound), dsname, smap_data)
      else
         write(msg,"(A,A,A)") "[cresp_io:save_cresp_smap_h5] Dataset ", "/"//n_g_smaps(bound)//"/"//dsname, " already present in file; omitting."
         call printinfo(msg)
      endif

   end subroutine save_cresp_smap_h5
!---------------------------------------------------------------------------------------------------
!> \brief Prepare and save rank-2 real array dataset in given file_id
!
   subroutine save_smap_to_open(file_id, path, dsname, dsdata)

      use helpers_hdf5,    only: create_dataset
      use hdf5,            only: HID_T, h5gopen_f, h5gclose_f

      implicit none

      real, dimension(:,:), target, intent(in) :: dsdata
      character(len=*),             intent(in) :: dsname, path
      integer(HID_T),               intent(in) :: file_id
      integer(HID_T)                           :: group_id
      integer(kind=4)                          :: error
      real, dimension(:,:), pointer            :: pdsdata

      call h5gopen_f(file_id, path, group_id, error)
      pdsdata => dsdata
      call create_dataset(group_id, dsname, pdsdata)
      call h5gclose_f(group_id, error)
      pdsdata => null()


   end subroutine save_smap_to_open
!----------------------------------------------------------------------------------------------------
!> \brief Save NR solution map to an ASCII file DEPRECATED
!
   subroutine save_NR_smap(NR_smap, hdr, vname, bc)

      use constants,        only: I_ONE
      use cresp_helpers,    only: bound_name, extension, flen, map_header

      implicit none

      integer(kind=4),      intent(in) :: bc
      integer(kind=4), parameter       :: flun = 31
      character(len=flen)              :: fname
      type(map_header),     intent(in) :: hdr
      integer(kind=4)                  :: j
      real, dimension(:,:), intent(in) :: NR_smap
      character(len=*),     intent(in) :: vname

      fname = vname // bound_name(bc) // extension
      open(flun, file=fname, status="unknown", position="rewind")
         write(flun,"(A56,A2,A26)") "This is a storage file for NR init grid, boundary case: ", bound_name(bc), &
            & ". Do not append this file."

         write(flun, "(1E15.8, 2I10,10E22.15)") hdr%s_es, hdr%s_dim1, hdr%s_dim2, hdr%s_pr, hdr%s_qbig, hdr%s_c, hdr%s_amin, hdr%s_amax, hdr%s_nmin, hdr%s_nmax
         write(flun, "(A1)") " "                             ! Blank line
         do j = I_ONE, size(NR_smap, dim=2, kind=4)
            write(flun, "(*(E24.15E3))") NR_smap(:,j)  ! WARNING - MIGHT NEED EXPLICIT ELEMENT COUNT IN LINE IN OLDER COMPILERS
         enddo
         close(flun)

   end subroutine save_NR_smap
!----------------------------------------------------------------------------------------------------
!> \brief Read solution map from an ASCII file ! DEPRECATED
   subroutine read_NR_smap(NR_smap, vname, bc, exit_code)

      use constants,       only: I_ONE
      use cresp_helpers,   only: bound_name, extension, flen
      use dataio_pub,      only: msg, warn

      implicit none

      real, dimension(:,:), intent(inout) :: NR_smap
      character(len=*),     intent(in)    :: vname
      integer(kind=4),      intent(in)    :: bc
      logical,              intent(out)   :: exit_code
      integer(kind=4)                     :: j, rstat = 0, flun = 31
      character(len=flen)                 :: fname

      fname = vname // bound_name(bc) // extension
      open(flun, file=fname, status="old", position="rewind", IOSTAT=rstat)
      if (rstat > 0) then
         write(msg, "(A8,I4,A8,2A20)") "IOSTAT:", rstat, ": file ", vname//bound_name(bc)//extension," does not exist!"
         call warn(msg)
         exit_code = .true.
         return
      else
         read(flun, *) ! Skipping comment line
         read(flun, *) ! Skipping header
         read(flun, *) ! Skipping blank line
         do j = I_ONE, size(NR_smap, dim=2, kind=4)
            read(flun, "(*(E24.15E3))", IOSTAT=rstat) NR_smap(:,j)  ! WARNING - THIS MIGHT NEED EXPLICIT INDICATION OF ELEMENTS COUNT IN LINE IN OLDER COMPILERS
         enddo
         exit_code = .false.
      endif
      if (rstat > 0) exit_code = .true.
      close(flun)

   end subroutine read_NR_smap
!---------------------------------------------------------------------------------------------------
!> \brief Read parameters used to construct solution map from an ASCII file ! DEPRECATED
!
   subroutine read_NR_smap_header(var_name, hdr, exit_code)

      use dataio_pub,      only: msg, warn
      use constants,       only: fmt_len
      use cresp_helpers,   only: extension, flen, map_header

      implicit none

      logical                          :: exit_code
      integer(kind=4), parameter       :: flun = 31
      character(len=flen)              :: f_name
      type(map_header), intent(inout)  :: hdr
      character(len=fmt_len)           :: fmt
      character(len=*), intent(in)     :: var_name
      integer(kind=4)                  :: fstat, rstat

      fstat = 0
      rstat = 0
      f_name = var_name // extension
      fmt = "(1E15.8,2I10,10E22.15)"

      open(flun, file=f_name, status="old", position="rewind", IOSTAT=fstat)

      if (fstat > 0) then
         write(msg,"(A8,I4,A8,2A20)") "IOSTAT:", fstat, ": file ", f_name, " does not exist!"
         call warn(msg)
         exit_code = .true.
         return
      endif

      read(flun, fmt, IOSTAT=rstat) hdr%s_es, hdr%s_dim1, hdr%s_dim2, hdr%s_pr, hdr%s_qbig, hdr%s_c, hdr%s_amin, hdr%s_amax, hdr%s_nmin, hdr%s_nmax
      if (rstat > 0 ) then  ! should work for older files using the same format
         read(flun, fmt, IOSTAT=rstat) hdr%s_es, hdr%s_dim1, hdr%s_dim2, hdr%s_pr, hdr%s_qbig, hdr%s_c
         hdr%s_amin = 0.
         hdr%s_amax = 0.
         hdr%s_nmin = 0.
         hdr%s_nmax = 0.
      endif

      exit_code = .false.
      close(flun)

   end subroutine read_NR_smap_header
!---------------------------------------------------------------------------------------------------
!> \brief Check parameters used to construct solution map against the loadable file.
!
   subroutine check_NR_smap_header(hdr, hdr_std, hdr_equal)

      use constants,     only: zero
      use cresp_helpers, only: map_header
      use dataio_pub,    only: msg, printinfo, warn
      use func,          only: operator(.equals.)

      implicit none

      type(map_header), intent(in)  :: hdr, hdr_std
      logical                       :: hdr_equal

      hdr_equal = .true.

      hdr_equal = hdr_equal .and. (hdr%s_es   .equals.   hdr_std%s_es)
      hdr_equal = hdr_equal .and. (hdr%s_dim1 .eq.       hdr_std%s_dim1)
      hdr_equal = hdr_equal .and. (hdr%s_dim2 .eq.       hdr_std%s_dim2)
      hdr_equal = hdr_equal .and. (hdr%s_qbig .equals.   hdr_std%s_qbig)
      hdr_equal = hdr_equal .and. (hdr%s_pr   .equals.   hdr_std%s_pr)
      hdr_equal = hdr_equal .and. (hdr%s_c    .equals.   hdr_std%s_c)

!  WARNING allowing to read old solution maps; without saved a_tab and n_tab limits
      hdr_equal = hdr_equal .and. ((hdr%s_amin .equals. hdr_std%s_amin) .or. (hdr%s_amin .equals. zero))
      hdr_equal = hdr_equal .and. ((hdr%s_amax .equals. hdr_std%s_amax) .or. (hdr%s_amax .equals. zero))
      hdr_equal = hdr_equal .and. ((hdr%s_nmin .equals. hdr_std%s_nmin) .or. (hdr%s_nmin .equals. zero))
      hdr_equal = hdr_equal .and. ((hdr%s_nmax .equals. hdr_std%s_nmax) .or. (hdr%s_nmax .equals. zero))

      if (hdr_equal) then
         write(msg,"(A109)") "[cresp_io:check_NR_smap_header] Headers match (provided in ratios files vs. values resulting from parameters)"
         call printinfo(msg)
      else
         write(msg,"(A110)") "[cresp_io:check_NR_smap_header] Headers differ (provided in ratios files vs. values resulting from parameters)"
         call warn(msg)
      endif

   end subroutine check_NR_smap_header

   subroutine check_NR_smaps_headers(hdr_res, hdr, hdr_match_res)

      use constants,     only: LO, HI
      use cresp_helpers, only: map_header

      implicit none

      type(map_header), dimension(2), intent(in)  :: hdr_res, hdr
      logical,          dimension(2), intent(out) :: hdr_match_res

      call check_NR_smap_header(hdr_res(LO), hdr(LO), hdr_match_res(LO))
      call check_NR_smap_header(hdr_res(HI), hdr(HI), hdr_match_res(HI))

   end subroutine check_NR_smaps_headers

!---------------------------------------------------------------------------------------------------
!> \brief Prepare and read cresp solution map fields from hdf5 file. !DEPRECATED
!
   subroutine read_cresp_smap_fields(read_error, filename_opt)

      use constants,  only: cwdlen
      use dataio_pub, only: msg, printinfo
      use hdf5,       only: HID_T, H5F_ACC_RDONLY_F, h5close_f, h5fclose_f, h5open_f, h5fopen_f

      implicit none

      integer(kind=4)                            :: error
      integer(HID_T)                             :: file_id
      character(len=cwdlen)                      :: filename
      character(len=*), optional,     intent(in) :: filename_opt
      logical                                    :: file_exist
      logical,                       intent(out) :: read_error

      if (present(filename_opt)) then
         filename = filename_opt
      else
         filename = "CRESP_smaps.h5"   ! FIXME - I should not be here! TODO arrange the code to read smaps and header in common_hdf5
         write (msg,"(A,A)") "[cresp_io:read_cresp_smap_fields] Name of file to open: ", trim(filename)
         call printinfo(msg)
      endif

      inquire(file = trim(filename), exist = file_exist)

      if (.not. file_exist) then
         write(msg, "(A,A,A)") "[cresp_io:read_cresp_smap_fields] file '", trim(filename), "' does not exist."
         call printinfo(msg)
         read_error = .true.
         return
      else
         write(msg, "(A,A,A)") "[cresp_io:read_cresp_smap_fields] Proceeding with '", trim(filename), "' file."
         call printinfo(msg)
         read_error = .false.
      endif

      call h5open_f(error)
      call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)

! WARNING format version omitted - CRESP is not included in older versions of piernik
      call read_smap_header_h5(file_id)

      call h5fclose_f(file_id, error)
      call h5close_f(error)

   end subroutine read_cresp_smap_fields
!---------------------------------------------------------------------------------------------------
!> \brief Read parameters used to construct solution map from an h5 file
!
   subroutine read_smap_header_h5(file_id, hdr_out)

      use constants,          only: LO, HI
      use cresp_helpers,      only: hdr_io, map_header, n_g_smaps, n_a_clight, n_a_dims, n_a_esmall, &
            &  n_a_max_p_r, n_a_qbig, n_a_amax, n_a_amin, n_a_nmax, n_a_nmin, real_attrs, int_attrs
      use dataio_pub,         only: die
      use hdf5,               only: HID_T
      use h5lt,               only: h5ltget_attribute_double_f, h5ltget_attribute_int_f

      implicit none

      integer                                               :: i, ia
      integer(kind=4)                                       :: error
      integer(HID_T),                            intent(in) :: file_id
      type(map_header), dimension(2), optional, intent(out) :: hdr_out
      integer(kind=4),  dimension(2)                        :: ibuf
      real,             dimension(1)                        :: rbuf

      do i = LO, HI              ! use statement for LO, HI in the upper level
         do ia = lbound(int_attrs, dim=1), ubound(int_attrs, dim=1)
            call h5ltget_attribute_int_f(file_id, n_g_smaps(i), trim(int_attrs(ia)), ibuf, error)
            if ((int_attrs(ia)) .eq. n_a_dims ) then
               hdr_io(i)%s_dim1    = ibuf(1)
               hdr_io(i)%s_dim2    = ibuf(2)
            else
               call die("[cresp_io:read_cresp_smap_fields] Non-recognized area_type.")
            endif
         enddo
         do ia = lbound(real_attrs, dim=1), ubound(real_attrs, dim=1)
            call h5ltget_attribute_double_f(file_id, n_g_smaps(i), trim(real_attrs(ia)), rbuf, error)
            select case (real_attrs(ia))
               case (n_a_esmall)
                  hdr_io(i)%s_es      = rbuf(1)
               case (n_a_max_p_r)
                  hdr_io(i)%s_pr      = rbuf(1)
               case (n_a_qbig)
                  hdr_io(i)%s_qbig    = rbuf(1)
               case (n_a_clight)
                  hdr_io(i)%s_c       = rbuf(1)
               case (n_a_amin)
                  hdr_io(i)%s_amin    = rbuf(1)
               case (n_a_amax)
                  hdr_io(i)%s_amax    = rbuf(1)
               case (n_a_nmin)
                  hdr_io(i)%s_nmin    = rbuf(1)
               case (n_a_nmax)
                  hdr_io(i)%s_nmax    = rbuf(1)
               case default
                  call die("[cresp_io:read_cresp_smap_fields] Non-recognized area_type.")
            end select
         enddo
      enddo

      if (present(hdr_out)) hdr_out = hdr_io

   end subroutine read_smap_header_h5
!---------------------------------------------------------------------------------------------------
!
!> \brief  Read 2D double precision dataset 'dsdata' of provided name 'dsetname' with 'file_id'
!  Requires file to be open. WARNING Does not contain fail-safe instructions!
!
   subroutine read_real_arr2d_dset(file_id, dsetname, dsdata)

      use hdf5,      only: HID_T, h5dopen_f, h5dread_f, h5dopen_f, &
                     &  h5dclose_f, H5T_NATIVE_DOUBLE, HSIZE_T
      implicit none

      real, dimension(:,:), target, intent(in) :: dsdata
      integer(HID_T)                           :: dset_id
      character(len=*),             intent(in) :: dsetname
      integer(HID_T),               intent(in) :: file_id
      integer(kind=4)                          :: error
      real, dimension(:,:),            pointer :: pa2d
      integer(HSIZE_T), dimension(2)           :: pa2ddims

      pa2d => dsdata(:,:)
      if (associated(pa2d)) then
         pa2ddims = shape(pa2d)

         call h5dopen_f(file_id, dsetname, dset_id, error)
         call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, pa2d, pa2ddims, error)
         call h5dclose_f(dset_id, error)
      endif
      pa2d => null()

   end subroutine read_real_arr2d_dset

!---------------------------------------------------------------------------------------------------


end module cresp_io
