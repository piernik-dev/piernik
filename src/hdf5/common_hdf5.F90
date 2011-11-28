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

#include "piernik.h"

!>
!! \brief Module that contains common I/O routines for HDF5 outputs (restart, data dumps and slices)
!<
module common_hdf5

! pulled by ANY

   use constants, only: singlechar, ndims

   implicit none

   private
   public :: init_hdf5, cleanup_hdf5, set_common_attributes, common_shortcuts, write_to_hdf5_v2
   public :: nhdf_vars, hdf_vars, d_gname, base_d_gname, d_fc_aname, d_size_aname, d_edge_apname, d_bnd_apname, cg_gname, &
         & cg_cnt_aname, cg_lev_aname, cg_size_aname, cg_offset_aname, n_cg_name, dir_pref

   integer, parameter :: S_LEN = 30

   character(len=S_LEN), allocatable, dimension(:), protected :: hdf_vars  !< dataset names for hdf files
   integer, protected :: nhdf_vars !< number of quantities plotted to hdf files
   character(len=*), parameter :: d_gname = "domains", base_d_gname = "base", d_fc_aname = "fine_count", &
        &                         d_size_aname = "n_d", d_edge_apname = "-edge_position", d_bnd_apname = "-boundary_type", &
        &                         cg_gname = "cg", cg_cnt_aname = "cg_count", cg_lev_aname = "level", cg_size_aname = "n_b", &
        &                         cg_offset_aname = "off"
   character(len=singlechar), dimension(ndims), parameter :: dir_pref = [ "x", "y", "z" ]

!> \brief Add an attribute (1D array) to the given group and initialize its value
   interface create_attribute
      module procedure create_int_attribute
      module procedure create_real_attribute
   end interface

contains

!>
!! \brief Procedure initializing HDF5 module
!<

   subroutine init_hdf5(vars)

      use constants,   only: varlen
      use fluids_pub,  only: has_ion, has_dst, has_neu
      use fluidindex,  only: iarr_all_dn, iarr_all_mx, iarr_all_my, iarr_all_mz
#ifdef COSM_RAYS
      use dataio_pub,  only: warn, msg
      use fluidindex,  only: iarr_all_crs
#endif /* COSM_RAYS */

      implicit none

      character(len=varlen), dimension(:), intent(in) :: vars  !< quantities to be plotted, see dataio::vars

      integer               :: nvars, i, j
#if defined COSM_RAYS
      integer               :: k
      character(len=varlen) :: aux
#endif /* COSM_RAYS */

      nvars = 1
      do while ( len(trim(vars(nvars))) > 1)
         nvars = nvars + 1
      enddo
      nvars = nvars - 1

      nhdf_vars = 0
      do i = 1, nvars
         select case (vars(i))
            case ('dens')
               nhdf_vars = nhdf_vars + size(iarr_all_dn,1)
            case ('velx')
               nhdf_vars = nhdf_vars + size(iarr_all_mx,1)
            case ('vely')
               nhdf_vars = nhdf_vars + size(iarr_all_my,1)
            case ('velz')
               nhdf_vars = nhdf_vars + size(iarr_all_mz,1)
            case ('ener')
               nhdf_vars = nhdf_vars + size(iarr_all_mz,1)
               if (has_dst) nhdf_vars = nhdf_vars - 1
#ifdef GRAV
            case ('gpot')
               nhdf_vars = nhdf_vars + 1
#ifdef MULTIGRID
            case ('mgso') ! multigrid solution
               nhdf_vars = nhdf_vars + 1
#endif /* MULTIGRID */
#endif /* GRAV */
            case ('magx', 'magy', 'magz', 'pres')
               nhdf_vars = nhdf_vars + 1
#ifdef COSM_RAYS
            case ('encr')
               nhdf_vars = nhdf_vars + size(iarr_all_crs,1)
#endif /* COSM_RAYS */
            case default
               nhdf_vars = nhdf_vars + 1
         end select
      enddo
      allocate(hdf_vars(nhdf_vars)); j = 1
      do i = 1, nvars
         select case (vars(i))
            case ('dens')
               if (has_dst) then ; hdf_vars(j) = 'dend' ; j = j + 1 ; endif
               if (has_neu) then ; hdf_vars(j) = 'denn' ; j = j + 1 ; endif
               if (has_ion) then ; hdf_vars(j) = 'deni' ; j = j + 1 ; endif
            case ('velx')
               if (has_dst) then ; hdf_vars(j) = 'vlxd' ; j = j + 1 ; endif
               if (has_neu) then ; hdf_vars(j) = 'vlxn' ; j = j + 1 ; endif
               if (has_ion) then ; hdf_vars(j) = 'vlxi' ; j = j + 1 ; endif
            case ('vely')
               if (has_dst) then ; hdf_vars(j) = 'vlyd' ; j = j + 1 ; endif
               if (has_neu) then ; hdf_vars(j) = 'vlyn' ; j = j + 1 ; endif
               if (has_ion) then ; hdf_vars(j) = 'vlyi' ; j = j + 1 ; endif
            case ('velz')
               if (has_dst) then ; hdf_vars(j) = 'vlzd' ; j = j + 1 ; endif
               if (has_neu) then ; hdf_vars(j) = 'vlzn' ; j = j + 1 ; endif
               if (has_ion) then ; hdf_vars(j) = 'vlzi' ; j = j + 1 ; endif
            case ('ener')
               if (has_neu) then ; hdf_vars(j) = 'enen' ; j = j + 1 ; endif
               if (has_ion) then ; hdf_vars(j) = 'enei' ; j = j + 1 ; endif
            case ("magx", "magy", "magz")
               hdf_vars(j) = vars(i) ; j = j + 1
#ifdef COSM_RAYS
            case ('encr')
               do k = 1, size(iarr_all_crs,1)
                  if (k<=9) then
                     write(aux,'(A2,I1)') 'cr', k
                     hdf_vars(j) = aux ; j = j + 1
                  else
                     write(msg, '(a,i3)')"[common_hdf5:init_hdf5] Cannot create name for CR energy component #", k
                     call warn(msg)
                  endif
               enddo
#endif /* COSM_RAYS */
#ifdef GRAV
            case ('gpot')
               hdf_vars(j) = 'gpot' ; j = j + 1
#ifdef MULTIGRID
            case ('mgso')
               hdf_vars(j) = 'mgso' ; j = j + 1
#endif /* MULTIGRID */
#endif /* GRAV */
            case ('pres')
               if (has_neu) then ; hdf_vars(j) = 'pren' ; j = j + 1 ; endif
               if (has_ion) then ; hdf_vars(j) = 'prei' ; j = j + 1 ; endif
            case default
               hdf_vars(j) = trim(vars(i)) ; j = j + 1
         end select
      enddo

   end subroutine init_hdf5

!>
!! \brief Procedure finalizing HDF5 module
!<
   subroutine cleanup_hdf5
      implicit none

      if (allocated(hdf_vars)) deallocate(hdf_vars)

   end subroutine cleanup_hdf5

!-----------------------------------------------------------------------------
!>
!! \brief decode some useful indices from variable name, if possible
!<
   subroutine common_shortcuts(var, fl_dni, i_xyz)

      use constants,  only: varlen, singlechar
      use fluidindex, only: flind
      use fluidtypes, only: component_fluid

      implicit none

      character(len=varlen),          intent(in)    :: var
      type(component_fluid), pointer, intent(inout) :: fl_dni
      integer,                        intent(out)   :: i_xyz

      character(len=singlechar) :: dc

      nullify(fl_dni)
      if (any([ "den", "vlx", "vly", "vlz", "ene" ] == var(1:3))) then
         select case (var(4:4))
            case ("d")
               fl_dni => flind%dst
            case ("n")
               fl_dni => flind%neu
            case ("i")
               fl_dni => flind%ion
         end select
      endif

      i_xyz = huge(1)
      if (var(1:2) == "vl") then
         dc = var(3:3)
      else if (var(1:3) == "mag") then
         dc = var(4:4)
      else
         dc = '_'
      endif
      if (any([ "x", "y", "z" ] == dc)) i_xyz = ichar(dc) - ichar("x")

   end subroutine common_shortcuts

!>
!! \brief This routine writes all attributes that are common to restart and output files.
!!
!! \details Write real, integer and character attributes. Store contents of problem.par and env files.
!! Other common elements may also be moved here.
!<
   subroutine set_common_attributes(filename)

      use constants,   only: I_ONE, I_NINE
      use dataio_pub,  only: use_v2_io, parfile, parfilelines
      use dataio_user, only: additional_attrs
      use global,      only: magic_mass, local_magic_mass
      use hdf5,        only: HID_T, SIZE_T, HSIZE_T, H5F_ACC_TRUNC_F, H5T_NATIVE_CHARACTER, H5Z_FILTER_DEFLATE_F, H5P_DATASET_CREATE_F, &
           &                 h5open_f, h5fcreate_f, h5fclose_f, H5Zfilter_avail_f, H5Pcreate_f, H5Pset_deflate_f, H5Pset_chunk_f, &
           &                 h5tcopy_f, h5tset_size_f, h5screate_simple_f, H5Dcreate_f, H5Dwrite_f, H5Dclose_f, H5Sclose_f, H5Tclose_f, H5Pclose_f, h5close_f
      use mpi,         only: MPI_DOUBLE_PRECISION, MPI_SUM
      use mpisetup,    only: slave, comm, ierr, FIRST
      use version,     only: env, nenv

      implicit none

      character(len=*), intent(in)   :: filename        !> HDF File name

      integer(HID_T)                 :: file_id         !> File identifier
      integer(HID_T)                 :: type_id, dspace_id, dset_id, prp_id
      integer(HSIZE_T), dimension(1) :: dimstr
      logical(kind=4)                :: Z_avail
      integer(SIZE_T)                :: maxlen
      integer(kind=4)                :: error
      real                           :: magic_mass0

      call MPI_Reduce(local_magic_mass, magic_mass0, I_ONE, MPI_DOUBLE_PRECISION, MPI_SUM, FIRST, comm, ierr)
      magic_mass       = magic_mass + magic_mass0
      local_magic_mass = 0.0

      if (slave) return ! This data need not be written in parallel.

      call h5open_f(error)
      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

      if (use_v2_io) then
         call set_common_attributes_v2(file_id)
      else
         call set_common_attributes_v1(file_id)
      endif

      ! Store a compressed copy of the problem.par file.
      maxlen = int(maxval(len_trim(parfile(:parfilelines))), kind=4)
      dimstr = [parfilelines]
      call H5Zfilter_avail_f(H5Z_FILTER_DEFLATE_F, Z_avail, error)
      ! call H5Zget_filter_info_f ! everything should be always fine for gzip
      call H5Pcreate_f(H5P_DATASET_CREATE_F, prp_id, error)
      if (Z_avail) then
         call H5Pset_deflate_f(prp_id, I_NINE, error)
         call H5Pset_chunk_f(prp_id, I_ONE, dimstr, error)
      endif
      call H5Tcopy_f(H5T_NATIVE_CHARACTER, type_id, error)
      call H5Tset_size_f(type_id, maxlen, error)
      call H5Screate_simple_f(I_ONE, dimstr, dspace_id, error)
      call H5Dcreate_f(file_id, "problem.par", type_id,  dspace_id, dset_id, error, dcpl_id = prp_id)
      call H5Dwrite_f(dset_id, type_id, parfile(:)(:maxlen), dimstr, error)
      call H5Dclose_f(dset_id, error)
      call H5Sclose_f(dspace_id, error)

      ! Store a compressed copy of the piernik.def file and Id lines from source files.
      ! We recycle type_id and prp_id, so we don't close them yet.
      maxlen = int(maxval(len_trim(env(:nenv))), kind=4)
      dimstr = [nenv]
      if (Z_avail) call H5Pset_chunk_f(prp_id, I_ONE, dimstr, error)
      call H5Tset_size_f(type_id, maxlen, error)
      call H5Screate_simple_f(I_ONE, dimstr, dspace_id, error)
      call H5Dcreate_f(file_id, "env", type_id,  dspace_id, dset_id, error, dcpl_id = prp_id)
      call H5Dwrite_f(dset_id, type_id, env(:)(:maxlen), dimstr, error)
      call H5Dclose_f(dset_id, error)
      call H5Sclose_f(dspace_id, error)
      call H5Tclose_f(type_id, error)
      call H5Pclose_f(prp_id, error)

      if (associated(additional_attrs)) call additional_attrs(file_id)

      call h5fclose_f(file_id, error)
      call h5close_f(error)

      ! only master process exits here

   end subroutine set_common_attributes

!>
!! \brief Common attributes for v2 files
!!
!! The rr1 marks critical attributes that are read by read_restart_hdf5 and compared against value read from the problem.par file.
!! The rr2 marks runtime values that are read by read_restart_hdf5 and assigned to something in the code.
!> \todo Set up an universal table(s) of attribute names for use by both set_common_attributes and read_restart_hdf5.
!! Provide indices for critical attributes (rr1) and for runtime attributes (rr2).
!!
!<

   subroutine set_common_attributes_v1(file_id)

      use constants,   only: cbuff_len, xdim, ydim, zdim, I_ONE
      use dataio_pub,  only: require_init_prob, piernik_hdf5_version, problem_name, run_id, last_hdf_time, next_t_tsl, next_t_log, nres, nhdf, domain_dump, step_hdf
      use domain,      only: dom
      use global,      only: magic_mass, t, dt, nstep
      use hdf5,        only: HID_T, SIZE_T
      use h5lt,        only: h5ltset_attribute_double_f, h5ltset_attribute_int_f, h5ltset_attribute_string_f

      implicit none

      integer(HID_T), intent(in)                   :: file_id       !> File identifier

      integer(kind=4)                              :: fe
      integer(SIZE_T)                              :: i
      integer(SIZE_T), parameter                   :: bufsize = I_ONE
      integer(kind=4)                              :: error
      integer, parameter                           :: buf_len = 50
      integer(kind=4),          dimension(buf_len) :: ibuffer
      real,                     dimension(buf_len) :: rbuffer
      character(len=cbuff_len), dimension(buf_len) :: ibuffer_name = ''
      character(len=cbuff_len), dimension(buf_len) :: rbuffer_name = ''

      rbuffer(1)   = t                       ; rbuffer_name(1)   = "time" !rr2
      rbuffer(2)   = dt                      ; rbuffer_name(2)   = "timestep" !rr2
      rbuffer(3)   = last_hdf_time           ; rbuffer_name(3)   = "last_hdf_time" !rr2
      rbuffer(4:5) = dom%edge(xdim, :)       ; rbuffer_name(4:5) = [ "xmin", "xmax" ] !rr1
      rbuffer(6:7) = dom%edge(ydim, :)       ; rbuffer_name(6:7) = [ "ymin", "ymax" ] !rr1
      rbuffer(8:9) = dom%edge(zdim, :)       ; rbuffer_name(8:9) = [ "zmin", "zmax" ] !rr1
      rbuffer(10)  = piernik_hdf5_version    ; rbuffer_name(10)  = "piernik" !rr1, rr2
      rbuffer(11)  = magic_mass              ; rbuffer_name(11)  = "magic_mass" !rr2
      rbuffer(12)  = next_t_tsl              ; rbuffer_name(12)  = "next_t_tsl" !rr2
      rbuffer(13)  = next_t_log              ; rbuffer_name(13)  = "next_t_log" !rr2

      ibuffer(1)   = nstep                   ; ibuffer_name(1)   = "nstep" !rr2
      ibuffer(2)   = nres + I_ONE            ; ibuffer_name(2)   = "nres" !rr2
      ibuffer(3)   = nhdf                    ; ibuffer_name(3)   = "nhdf" !rr2
      ibuffer(4)   = nstep                   ; ibuffer_name(4)   = "step_res" !rr2
      ibuffer(5)   = step_hdf                ; ibuffer_name(5)   = "step_hdf" !rr2
      ibuffer(6:8) = dom%n_d(:)              ; ibuffer_name(6:8) = [ "nxd", "nyd", "nzd" ] !rr1
      ibuffer(9)   = dom%nb                  ; ibuffer_name(9)   = "nb" ! BEWARE: assuming cga%cg_all(:)%nb equal everywhere
      ibuffer(10)  = require_init_prob       ; ibuffer_name(10)  = "require_init_prob" !rr2

      i = 1
      do while (rbuffer_name(i) /= "")
         call h5ltset_attribute_double_f(file_id, "/", rbuffer_name(i), rbuffer(i), bufsize, error)
         i = i + bufsize
      enddo

      i = 1
      do while (ibuffer_name(i) /= "")
         call h5ltset_attribute_int_f(file_id, "/", ibuffer_name(i), ibuffer(i), bufsize, error)
         i = i + bufsize
      enddo

      fe = len_trim(problem_name, kind=4)
      call h5ltset_attribute_string_f(file_id, "/", "problem_name", problem_name(1:fe), error) !rr2
      fe = len_trim(domain_dump, kind=4)
      call h5ltset_attribute_string_f(file_id, "/", "domain", domain_dump(1:fe), error) !rr2
      fe = len_trim(run_id, kind=4)
      call h5ltset_attribute_string_f(file_id, "/", "run_id", run_id(1:fe), error) !rr2

   end subroutine set_common_attributes_v1

!>
!! \brief Common attributes for v2 files
!!
!! \warning Do not remove redundancy between set_common_attributes_v1 and set_common_attributes_v2 until we get mature state of the v2 outputs
!<

   subroutine set_common_attributes_v2(file_id)

      use constants,   only: cbuff_len, I_ONE
      use dataio_pub,  only: require_init_prob, piernik_hdf5_version2, problem_name, run_id, last_hdf_time, next_t_tsl, next_t_log, nres, nhdf, domain_dump, step_hdf
      use global,      only: magic_mass, t, dt, nstep
      use hdf5,        only: HID_T, SIZE_T
      use h5lt,        only: h5ltset_attribute_double_f, h5ltset_attribute_int_f, h5ltset_attribute_string_f

      implicit none

      integer(HID_T), intent(in)                   :: file_id       !> File identifier

      integer(kind=4)                              :: fe
      integer(SIZE_T)                              :: i
      integer(SIZE_T), parameter                   :: bufsize = I_ONE
      integer(kind=4)                              :: error
      integer, parameter                           :: buf_len = 50
      integer(kind=4),          dimension(buf_len) :: ibuffer
      real,                     dimension(buf_len) :: rbuffer
      character(len=cbuff_len), dimension(buf_len) :: ibuffer_name = ''
      character(len=cbuff_len), dimension(buf_len) :: rbuffer_name = ''

      rbuffer(1) = t                     ; rbuffer_name(1) = "time" !rr2
      rbuffer(2) = dt                    ; rbuffer_name(2) = "timestep" !rr2
      rbuffer(3) = last_hdf_time         ; rbuffer_name(3) = "last_hdf_time" !rr2
      rbuffer(4) = piernik_hdf5_version2 ; rbuffer_name(4) = "piernik" !rr1, rr2
      rbuffer(5) = magic_mass            ; rbuffer_name(5) = "magic_mass" !rr2
      rbuffer(6) = next_t_tsl            ; rbuffer_name(6) = "next_t_tsl" !rr2
      rbuffer(7) = next_t_log            ; rbuffer_name(7) = "next_t_log" !rr2

      ibuffer(1) = nstep                 ; ibuffer_name(1) = "nstep" !rr2
      ibuffer(2) = nres + I_ONE          ; ibuffer_name(2) = "nres" !rr2
      ibuffer(3) = nhdf                  ; ibuffer_name(3) = "nhdf" !rr2
      ibuffer(4) = nstep                 ; ibuffer_name(4) = "step_res" !rr2
      ibuffer(5) = step_hdf              ; ibuffer_name(5) = "step_hdf" !rr2
      ibuffer(6) = require_init_prob     ; ibuffer_name(6) = "require_init_prob" !rr2

      !> \todo  add number of pieces in the restart point/data dump

      i = 1
      do while (rbuffer_name(i) /= "")
         call h5ltset_attribute_double_f(file_id, "/", rbuffer_name(i), rbuffer(i), bufsize, error)
         i = i + bufsize
      enddo

      i = 1
      do while (ibuffer_name(i) /= "")
         call h5ltset_attribute_int_f(file_id, "/", ibuffer_name(i), ibuffer(i), bufsize, error)
         i = i + bufsize
      enddo

      fe = len_trim(problem_name, kind=4)
      call h5ltset_attribute_string_f(file_id, "/", "problem_name", problem_name(1:fe), error) !rr2
      fe = len_trim(domain_dump, kind=4)
      call h5ltset_attribute_string_f(file_id, "/", "domain", domain_dump(1:fe), error) !rr2
      fe = len_trim(run_id, kind=4)
      call h5ltset_attribute_string_f(file_id, "/", "run_id", run_id(1:fe), error) !rr2

      ! these values will go do  base domain description
!!$      rbuffer(4:5) = dom%edge(xdim, :)       ; rbuffer_name(4:5) = [ "xmin", "xmax" ] !rr1
!!$      rbuffer(6:7) = dom%edge(ydim, :)       ; rbuffer_name(6:7) = [ "ymin", "ymax" ] !rr1
!!$      rbuffer(8:9) = dom%edge(zdim, :)       ; rbuffer_name(8:9) = [ "zmin", "zmax" ] !rr1
!!$      ibuffer(6:8) = dom%n_d(:)              ; ibuffer_name(6:8) = [ "nxd", "nyd", "nzd" ] !rr1
!!$      ibuffer(9)   = dom%nb                  ; ibuffer_name(9)   = "nb" ! BEWARE: assuming cga%cg_all(:)%nb equal everywhere
!!$      external boundary types

   end subroutine set_common_attributes_v2

!> \brief Attach an 32-bit integer attribute (scalar or rank-1 small array) to the given group.

   subroutine create_int_attribute(g_id, name, int_array)

     use constants, only: I_ONE
     use hdf5,      only: H5T_NATIVE_INTEGER, HID_T, HSIZE_T, &
          &               h5acreate_f, h5aclose_f, h5awrite_f, h5screate_simple_f, h5sclose_f

     implicit none

     integer(HID_T), intent(in)                :: g_id      !< group id where to create the attribute
     character(len=*), intent(in)              :: name      !< name
     integer(kind=4), dimension(:), intent(in) :: int_array !< the data

     integer(HID_T)  :: aspace_id, attr_id
     integer(kind=4) :: error

     call h5screate_simple_f(I_ONE, [ size(int_array, kind=HSIZE_T) ], aspace_id, error)
     call h5acreate_f(g_id, name, H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
     call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, int_array, [ size(int_array, kind=HSIZE_T) ], error)
     call h5aclose_f(attr_id, error)
     call h5sclose_f(aspace_id, error)

   end subroutine create_int_attribute

!> \brief Attach an 64-bit real attribute (scalar or rank-1 small array) to the given group.

   subroutine create_real_attribute(g_id, name, real_array)

     use constants, only: I_ONE
     use hdf5,      only: H5T_NATIVE_DOUBLE, HID_T, HSIZE_T, &
          &               h5acreate_f, h5aclose_f, h5awrite_f, h5screate_simple_f, h5sclose_f

     implicit none

     integer(HID_T), intent(in)     :: g_id       !< group id where to create the attribute
     character(len=*), intent(in)   :: name       !< name
     real, dimension(:), intent(in) :: real_array !< the data

     integer(HID_T)  :: aspace_id, attr_id
     integer(kind=4) :: error

     call h5screate_simple_f(I_ONE, [ size(real_array, kind=HSIZE_T) ], aspace_id, error)
     call h5acreate_f(g_id, name, H5T_NATIVE_DOUBLE, aspace_id, attr_id, error)
     call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, real_array, [ size(real_array, kind=HSIZE_T) ], error)
     call h5aclose_f(attr_id, error)
     call h5sclose_f(aspace_id, error)

   end subroutine create_real_attribute

!> \brief Generate numbered cg group name
   function n_cg_name(g)
      use constants, only: dsetnamelen
      implicit none
      integer, intent(in)        :: g !< group number
      character(len=dsetnamelen) :: n_cg_name
      write(n_cg_name,'(2a,i8.8)')trim(cg_gname), "_", g
   end function n_cg_name
!>

!! \brief Write a multi-file, multi-domain HDF5 file
!!
!! \details There are three approaches to be implemented:
!! - Single-file, serial I/O. The easiest way. Master writes everything, slaves send their data to the master. Does not take advantage of parallel filesystems. Best choice for non-parallel filesystems.
!! - Multi-file, serial I/O. An extension of the above approach. Selected processes (can be all) write to their files, other processes send them their data.
!!   Can take advantage of parallel filesystems. Can use local scratch filesystems. Requires additional work on reading.
!! - Single-file, parallel I/O. The most ambitious approach. Selected processes (can be all) write to the files, other processes send them their data.
!!   Can take advantage of parallel filesystems. Currently does not allow for compression during write.
!!   Requires a lot of pseudo-collective operations. The "flexible PHDF5" would simplify the code, but it needs to be implemented.first.
!!
!! \warning Partial implementation: Single-file, serial I/O works for non-AMR setups.
!<

   subroutine write_to_hdf5_v2(filename, create_empty_cg_datasets, write_cg_to_hdf5)

      use constants,   only: cwdlen, dsetnamelen, xdim, zdim, ndims, I_ONE, I_TWO, INT4
      use dataio_pub,  only: die, nproc_io, can_i_write
      use dataio_user, only: problem_write_restart
      use domain,      only: dom
      use gc_list,     only: cg_list_element
      use grid,        only: all_cg
      use hdf5,        only: HID_T, H5F_ACC_RDWR_F, H5P_FILE_ACCESS_F, H5Z_FILTER_DEFLATE_F, &
           &                 h5open_f, h5close_f, h5fopen_f, h5fclose_f, h5gcreate_f, h5gopen_f, h5gclose_f, &
           &                 h5pcreate_f, h5pclose_f, h5pset_fapl_mpio_f, h5zfilter_avail_f
      use mpi,         only: MPI_INFO_NULL, MPI_INTEGER, MPI_INTEGER8, MPI_STATUS_IGNORE
      use mpisetup,    only: comm, proc, FIRST, LAST, master

      implicit none
      character(len=cwdlen), intent(in)             :: filename
      interface
         subroutine create_empty_cg_datasets(cgl_g_id, cg_n_b, Z_avail, g)
            use hdf5,     only: HID_T
            implicit none
            integer(HID_T), intent(in)                           :: cgl_g_id
            integer(kind=4), dimension(:,:), pointer, intent(in) :: cg_n_b
            logical(kind=4), intent(in)                          :: Z_avail
            integer, intent(in)                                  :: g
         end subroutine create_empty_cg_datasets

         subroutine write_cg_to_hdf5(cgl_g_id, cg_n, cg_all_n_b)
            use hdf5,     only: HID_T
            implicit none
            integer(HID_T), intent(in)                           :: cgl_g_id
            integer(kind=4), dimension(:),   pointer, intent(in) :: cg_n
            integer(kind=4), dimension(:,:), pointer, intent(in) :: cg_all_n_b
         end subroutine write_cg_to_hdf5
      end interface

      integer(HID_T)                                :: file_id                                  !> File identifier
      integer(HID_T)                                :: plist_id                                 !> Property list identifier
      integer(HID_T)                                :: cgl_g_id,  cg_g_id                       !> cg list and cg group identifiers
      integer(HID_T)                                :: doml_g_id, dom_g_id                      !> domain list and domain group identifiers
      integer(kind=4)                               :: error, cg_cnt
      integer                                       :: g, p, i
      integer, parameter                            :: tag = I_ONE
      integer(kind=4),  dimension(:),   pointer     :: cg_n                                     !> offset for cg group numbering
      integer(kind=4),  dimension(:,:), pointer     :: cg_all_n_b                               !> sizes of all cg
      integer(kind=4),  dimension(:),   allocatable :: cg_rl                                    !> list of refinement levels from all cgs/procs
      integer(kind=4),  dimension(:,:), pointer     :: cg_n_b                                   !> list of n_b from all cgs/procs
      integer(kind=8),  dimension(:,:), allocatable :: cg_off                                   !> list of offsets from all cgs/procs
      type(cg_list_element), pointer                :: cgl
      logical(kind=4)                               :: Z_avail                                  !> .true. if HDF5 was compiled with zlib support
      character(len=dsetnamelen)                    :: d_label

      ! Create a new file and initialize it

      ! Prepare groups and datasets for grid containers on the master
      allocate(cg_n(FIRST:LAST))
      call MPI_Allgather(all_cg%cnt, I_ONE, MPI_INTEGER, cg_n, I_ONE, MPI_INTEGER, comm, error)
      cg_cnt = sum(cg_n(:))
      allocate(cg_all_n_b(cg_cnt, ndims))

      if (master) then

         ! Open the HDF5 file only in master process and create all groups required for cg storage.
         ! Create also all related datasets and attributes. Do not write big datasets yet.

         call h5open_f(error)
         call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

         call h5gcreate_f(file_id, cg_gname, cgl_g_id, error)                                 ! create "/cg"

         call create_attribute(cgl_g_id, cg_cnt_aname, [ cg_cnt ])                            ! create "/cg/cg_count"

         Z_avail = .false.
         if (nproc_io == 1) call h5zfilter_avail_f(H5Z_FILTER_DEFLATE_F, Z_avail, error)
         !> \todo test it thoroughly before enabling for > 1

         ! Do not assume that the master knows all the lists
         do p = FIRST, LAST
            allocate(cg_rl(cg_n(p)), cg_n_b(cg_n(p), ndims), cg_off(cg_n(p), ndims))
            if (p == FIRST) then
               g = 1
               cgl => all_cg%first
               do while (associated(cgl))
                  cg_rl(g)     = 1  !> \deprecated change this once we introduce many levels
                  cg_n_b(g, :) = cgl%cg%n_b(:)
                  cg_off(g, :) = cgl%cg%off(:)
                  g = g + 1
                  cgl => cgl%nxt
               enddo
            else
               call MPI_Recv(cg_rl,  size(cg_rl),  MPI_INTEGER,  p, tag,       comm, MPI_STATUS_IGNORE, error)
               call MPI_Recv(cg_n_b, size(cg_n_b), MPI_INTEGER,  p, tag+I_ONE, comm, MPI_STATUS_IGNORE, error)
               call MPI_Recv(cg_off, size(cg_off), MPI_INTEGER8, p, tag+I_TWO, comm, MPI_STATUS_IGNORE, error)
            endif

            do g = 1, cg_n(p)
               call h5gcreate_f(cgl_g_id, n_cg_name(sum(cg_n(:p))-cg_n(p)+g), cg_g_id, error) ! create "/cg/cg_%08d"

               call create_attribute(cg_g_id, cg_lev_aname, [ cg_rl(g) ] )                    ! create "/cg/cg_%08d/level"
               call create_attribute(cg_g_id, cg_size_aname, cg_n_b(g, :))                    ! create "/cg/cg_%08d/n_b"
               call create_attribute(cg_g_id, cg_offset_aname, int(cg_off(g, :), kind=4))     ! create "/cg/cg_%08d/off"

               cg_all_n_b(sum(cg_n(:p))-cg_n(p)+g, :) = cg_n_b(g, :)

               if (any(cg_off(g, :) > 2.**31)) call die("[restart_hdf5:write_restart_hdf5_v2] large offsets require better treatment")

               call create_empty_cg_datasets(cg_g_id, cg_n_b, Z_avail, g) !!!!!

               call h5gclose_f(cg_g_id, error)
            enddo

            deallocate(cg_rl, cg_n_b, cg_off)
         enddo

         call h5gclose_f(cgl_g_id, error)

         ! describe_domains
         call h5gcreate_f(file_id, d_gname, doml_g_id, error)                    ! create "/domains"

         call h5gcreate_f(doml_g_id, base_d_gname, dom_g_id, error)              ! create "/domains/base"
         call create_attribute(dom_g_id, d_size_aname, dom%n_d(:))               ! create "/domains/base/n_d"
         do i = xdim, zdim
            write(d_label, '(2a)') dir_pref(i), d_edge_apname
            call create_attribute(dom_g_id, d_label, dom%edge(i, :))             ! create "/domains/base/[xyz]-edge_position"
            write(d_label, '(2a)') dir_pref(i), d_bnd_apname
            call create_attribute(dom_g_id, d_label, int(dom%bnd(i, :), kind=4)) ! create "/domains/base/[xyz]-boundary_type"
         enddo

         call h5gclose_f(dom_g_id, error)

         call create_attribute(doml_g_id, d_fc_aname, [ 0_INT4 ] )               ! create "/domains/fine_count"  ! we have only base domain at the moment

         !> \todo add here all fine domains
         ! name "fine_00000001"
         ! attributes: n_d(:), off(:), refinement
         call h5gclose_f(doml_g_id, error)

         call h5fclose_f(file_id, error)
         call h5close_f(error)

      else ! send all the necessary information to the master
         allocate(cg_rl(all_cg%cnt), cg_n_b(all_cg%cnt, ndims), cg_off(all_cg%cnt, ndims))
         g = 1
         cgl => all_cg%first
         do while (associated(cgl))
            cg_rl(g)     = 1  !> \deprecated change this once we introduce many levels
            cg_n_b(g, :) = cgl%cg%n_b(:)
            cg_off(g, :) = cgl%cg%off(:)
            g = g + 1
            cgl => cgl%nxt
         enddo
         call MPI_Send(cg_rl,  size(cg_rl),  MPI_INTEGER,  FIRST, tag,       comm, error)
         call MPI_Send(cg_n_b, size(cg_n_b), MPI_INTEGER,  FIRST, tag+I_ONE, comm, error)
         call MPI_Send(cg_off, size(cg_off), MPI_INTEGER8, FIRST, tag+I_TWO, comm, error)
         deallocate(cg_rl, cg_n_b, cg_off)
      endif

      call MPI_Bcast(cg_all_n_b, size(cg_all_n_b), MPI_INTEGER, FIRST, comm, error)

      call MPI_Barrier(comm, error)
      ! Reopen the HDF5 file for parallel write
      call h5open_f(error)
      if (can_i_write) then
         if (nproc_io > 1) then
            call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
            ! when nproc_io < nproc we'll probably need another communicator for subset of processes that have can_i_write flag set
            call h5pset_fapl_mpio_f(plist_id, comm, MPI_INFO_NULL, error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
            call h5pclose_f(plist_id, error)
         else
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
         endif
         call h5gopen_f(file_id, cg_gname, cgl_g_id, error)
      endif

      call write_cg_to_hdf5(cgl_g_id, cg_n, cg_all_n_b) !!!!!

      if (can_i_write) then
         call h5gclose_f(cgl_g_id, error)
         if (associated(problem_write_restart)) call problem_write_restart(file_id) ! may be called by not all processes
         call h5fclose_f(file_id, error)  ! Close the file
      endif

      call h5close_f(error)            ! Close HDF5 stuff

      deallocate(cg_n, cg_all_n_b)

   end subroutine write_to_hdf5_v2

! This routine will become useful when we begin to use multiple domain containers (AMR or non-rectangular compound domains)
! This routine will become obsolete or will need serious rework with HDF5-1.8.8 ~xarth
!!$   subroutine write_grid_containter(cg, file_id, plist_id)
!!$
!!$      use constants, only: xdim, ydim, zdim, ndims, dsetnamelen, LO, HI, I_ONE, I_TWO, I_THREE, I_FOUR, INT4
!!$      use grid_cont, only: grid_container
!!$      use hdf5,      only: HID_T, SIZE_T, HSIZE_T, H5T_NATIVE_INTEGER, H5T_STD_I8LE, H5T_NATIVE_DOUBLE, H5T_COMPOUND_F, &
!!$           &               h5screate_simple_f, h5tarray_create_f, h5tget_size_f, h5tcreate_f, h5tinsert_f, h5dwrite_f, h5sclose_f, h5tclose_f, h5dclose_f, h5dcreate_f
!!$
!!$      implicit none
!!$
!!$      type(grid_container), pointer, intent(in) :: cg
!!$      integer(HID_T), intent(in)                :: file_id, plist_id
!!$
!!$      integer(SIZE_T), parameter :: n_int4 = 19, n_r8 = 14, n_nxarr_r8 = 4, n_nyarr_r8 = 4, n_nzarr_r8 = 4, &
!!$         & n_ndims_r8 = 2, n_ndims_i4 =1, n_ndims_i8 = 1, n_ndims_lohi_i4 = 2
!!$
!!$      integer(SIZE_T) :: n_arr3d_r8, n_ndims_arr4d_r8, n_u_arr4d_r8, n_stub
!!$      integer :: total_no
!!$
!!$      integer(SIZE_T) :: int4_ts, r8_ts, nxarr_r8_ts, nyarr_r8_ts, nzarr_r8_ts, arr3d_r8_ts, ndims_r8_ts, &
!!$         & ndims_i4_ts, ndims_i8_ts, ndims_lohi_i4_ts, ndims_arr4d_r8_ts, u_arr4d_r8_ts, type_size, offset
!!$      integer(HID_T)  :: ndims_r8_t, ndims_i4_t, ndims_i8_t, nxarr_r8_t, nyarr_r8_t, nzarr_r8_t, arr3d_r8_t, ndims_lohi_i4_t, &
!!$         & ndims_arr4d_r8_t, u_arr4d_r8_t, dtype_id, dspace_id, dset_id
!!$      integer(HSIZE_T),  dimension(1) :: dims
!!$      integer(HID_T),    dimension(:), allocatable :: types, dmem_id
!!$      integer(SIZE_T),   dimension(:), allocatable :: types_sizes
!!$      character(len=dsetnamelen), dimension(:), allocatable :: types_names
!!$      character(len=dsetnamelen) :: dset_name
!!$
!!$      integer(kind=4) :: error
!!$      integer :: i
!!$
!!$      dims = 1
!!$      call h5screate_simple_f(I_ONE, dims, dspace_id, error)
!!$
!!$      call h5tarray_create_f(H5T_NATIVE_DOUBLE,  I_ONE, [integer(HSIZE_T):: ndims],              ndims_r8_t, error)
!!$      call h5tarray_create_f(H5T_NATIVE_INTEGER, I_ONE, [integer(HSIZE_T):: ndims],              ndims_i4_t, error)
!!$      call h5tarray_create_f(H5T_NATIVE_INTEGER, I_TWO, [integer(HSIZE_T):: ndims, HI-LO+1],     ndims_lohi_i4_t, error)
!!$      call h5tarray_create_f(H5T_STD_I8LE,       I_ONE, [integer(HSIZE_T):: ndims],              ndims_i8_t, error)
!!$      call h5tarray_create_f(H5T_NATIVE_DOUBLE,  I_ONE, [integer(HSIZE_T):: cg%n_(xdim)],        nxarr_r8_t, error)
!!$      call h5tarray_create_f(H5T_NATIVE_DOUBLE,  I_ONE, [integer(HSIZE_T):: cg%n_(ydim)],        nyarr_r8_t, error)
!!$      call h5tarray_create_f(H5T_NATIVE_DOUBLE,  I_ONE, [integer(HSIZE_T):: cg%n_(zdim)],        nzarr_r8_t, error)
!!$      call h5tarray_create_f(H5T_NATIVE_DOUBLE,  I_THREE, [integer(HSIZE_T):: cg%n_(:)   ],        arr3d_r8_t, error)
!!$      call h5tarray_create_f(H5T_NATIVE_DOUBLE,  I_FOUR, [integer(HSIZE_T):: ndims, cg%n_(:)],    ndims_arr4d_r8_t, error)
!!$      call h5tarray_create_f(H5T_NATIVE_DOUBLE,  I_FOUR, [integer(HSIZE_T):: size(cg%u,1), cg%n_(:)], u_arr4d_r8_t, error)
!!$
!!$      n_arr3d_r8 = 10  ! gc_{x,y,z}dim
!!$      n_stub     = 0
!!$      if (associated(cg%cs_iso2))  n_stub = n_stub + I_ONE
!!$      if (associated(cg%wa))       n_stub = n_stub + I_ONE
!!$      if (associated(cg%gpot))     n_stub = n_stub + I_ONE
!!$      if (associated(cg%hgpot))    n_stub = n_stub + I_ONE
!!$      if (associated(cg%gp))       n_stub = n_stub + I_ONE
!!$      if (associated(cg%sgp))      n_stub = n_stub + I_ONE
!!$      if (associated(cg%sgpm))     n_stub = n_stub + I_ONE
!!$      n_arr3d_r8 = n_arr3d_r8 - n_stub
!!$
!!$      n_ndims_arr4d_r8 = 1  ! b
!!$      if (associated(cg%b0%arr)) then
!!$         n_ndims_arr4d_r8 = n_ndims_arr4d_r8 + I_ONE
!!$      else
!!$         n_stub = n_stub + I_ONE
!!$      endif
!!$
!!$      n_u_arr4d_r8 = 2 ! u,uh
!!$      if (associated(cg%u0%arr)) then
!!$         n_u_arr4d_r8 = n_u_arr4d_r8 + I_ONE
!!$      else
!!$         n_stub = n_stub + I_ONE
!!$      endif
!!$
!!$      call h5tget_size_f(H5T_NATIVE_INTEGER, int4_ts,error)
!!$      call h5tget_size_f(H5T_NATIVE_DOUBLE,  r8_ts, error)
!!$      call h5tget_size_f(ndims_r8_t,         ndims_r8_ts, error)
!!$      call h5tget_size_f(ndims_i4_t,         ndims_i4_ts, error)
!!$      call h5tget_size_f(ndims_i8_t,         ndims_i8_ts, error)
!!$      call h5tget_size_f(ndims_lohi_i4_t,    ndims_lohi_i4_ts, error)
!!$      call h5tget_size_f(nxarr_r8_t,         nxarr_r8_ts, error)
!!$      call h5tget_size_f(nyarr_r8_t,         nyarr_r8_ts, error)
!!$      call h5tget_size_f(nzarr_r8_t,         nzarr_r8_ts, error)
!!$      call h5tget_size_f(arr3d_r8_t,         arr3d_r8_ts, error)
!!$      call h5tget_size_f(ndims_arr4d_r8_t,   ndims_arr4d_r8_ts, error)
!!$      call h5tget_size_f(u_arr4d_r8_t,       u_arr4d_r8_ts, error)
!!$
!!$      type_size = (n_int4+n_stub)*int4_ts + n_r8*r8_ts + n_ndims_r8*ndims_r8_ts + n_ndims_i4*ndims_i4_ts + n_ndims_i8*ndims_i8_ts + n_ndims_lohi_i4*ndims_lohi_i4_ts &
!!$            & + n_nxarr_r8*nxarr_r8_ts + n_nyarr_r8*nyarr_r8_ts + n_nzarr_r8*nzarr_r8_ts + n_arr3d_r8*arr3d_r8_ts + n_ndims_arr4d_r8*ndims_arr4d_r8_ts + n_u_arr4d_r8*u_arr4d_r8_ts
!!$
!!$      call h5tcreate_f(H5T_COMPOUND_F, type_size, dtype_id, error)
!!$
!!$      total_no =  int(n_int4 + n_r8 + n_ndims_r8 + n_ndims_i4  + n_ndims_i8 + n_ndims_lohi_i4 + n_nxarr_r8 + n_nyarr_r8 + &
!!$         n_nzarr_r8 + n_arr3d_r8 + n_ndims_arr4d_r8 + n_u_arr4d_r8 + n_stub, kind(total_no))
!!$
!!$      allocate(types(total_no), types_sizes(total_no), types_names(total_no), dmem_id(total_no))
!!$
!!$      types(1)  = H5T_NATIVE_INTEGER;  types_sizes(1)  = int4_ts;        types_names(1)  = "nx"
!!$      types(2)  = H5T_NATIVE_INTEGER;  types_sizes(2)  = int4_ts;        types_names(2)  = "ny"
!!$      types(3)  = H5T_NATIVE_INTEGER;  types_sizes(3)  = int4_ts;        types_names(3)  = "nz"
!!$      types(4)  = H5T_NATIVE_INTEGER;  types_sizes(4)  = int4_ts;        types_names(4)  = "nxb"
!!$      types(5)  = H5T_NATIVE_INTEGER;  types_sizes(5)  = int4_ts;        types_names(5)  = "nyb"
!!$      types(6)  = H5T_NATIVE_INTEGER;  types_sizes(6)  = int4_ts;        types_names(6)  = "nzb"
!!$      types(7)  = H5T_NATIVE_INTEGER;  types_sizes(7)  = int4_ts;        types_names(7)  = "is"
!!$      types(8)  = H5T_NATIVE_INTEGER;  types_sizes(8)  = int4_ts;        types_names(8)  = "ie"
!!$      types(9)  = H5T_NATIVE_INTEGER;  types_sizes(9)  = int4_ts;        types_names(9)  = "js"
!!$      types(10) = H5T_NATIVE_INTEGER;  types_sizes(10) = int4_ts;        types_names(10) = "je"
!!$      types(11) = H5T_NATIVE_INTEGER;  types_sizes(11) = int4_ts;        types_names(11) = "ks"
!!$      types(12) = H5T_NATIVE_INTEGER;  types_sizes(12) = int4_ts;        types_names(12) = "ke"
!!$      types(13) = H5T_NATIVE_INTEGER;  types_sizes(13) = int4_ts;        types_names(13) = "maxxyz"
!!$      types(14) = H5T_NATIVE_INTEGER;  types_sizes(14) = int4_ts;        types_names(14) = "isb"
!!$      types(15) = H5T_NATIVE_INTEGER;  types_sizes(15) = int4_ts;        types_names(15) = "ieb"
!!$      types(16) = H5T_NATIVE_INTEGER;  types_sizes(16) = int4_ts;        types_names(16) = "jsb"
!!$      types(17) = H5T_NATIVE_INTEGER;  types_sizes(17) = int4_ts;        types_names(17) = "jeb"
!!$      types(18) = H5T_NATIVE_INTEGER;  types_sizes(18) = int4_ts;        types_names(18) = "ksb"
!!$      types(19) = H5T_NATIVE_INTEGER;  types_sizes(19) = int4_ts;        types_names(19) = "keb"
!!$
!!$      types(20) = H5T_NATIVE_DOUBLE;   types_sizes(20) = r8_ts;          types_names(20) = "dx"
!!$      types(21) = H5T_NATIVE_DOUBLE;   types_sizes(21) = r8_ts;          types_names(21) = "dy"
!!$      types(22) = H5T_NATIVE_DOUBLE;   types_sizes(22) = r8_ts;          types_names(22) = "dz"
!!$      types(23) = H5T_NATIVE_DOUBLE;   types_sizes(23) = r8_ts;          types_names(23) = "idx"
!!$      types(24) = H5T_NATIVE_DOUBLE;   types_sizes(24) = r8_ts;          types_names(24) = "idy"
!!$      types(25) = H5T_NATIVE_DOUBLE;   types_sizes(25) = r8_ts;          types_names(25) = "idz"
!!$      types(26) = H5T_NATIVE_DOUBLE;   types_sizes(26) = r8_ts;          types_names(26) = "dxmn"
!!$      types(27) = H5T_NATIVE_DOUBLE;   types_sizes(27) = r8_ts;          types_names(27) = "dvol"
!!$      types(28) = H5T_NATIVE_DOUBLE;   types_sizes(28) = r8_ts;          types_names(28) = "xminb"
!!$      types(29) = H5T_NATIVE_DOUBLE;   types_sizes(29) = r8_ts;          types_names(29) = "xmaxb"
!!$      types(30) = H5T_NATIVE_DOUBLE;   types_sizes(30) = r8_ts;          types_names(30) = "yminb"
!!$      types(31) = H5T_NATIVE_DOUBLE;   types_sizes(31) = r8_ts;          types_names(31) = "ymaxb"
!!$      types(32) = H5T_NATIVE_DOUBLE;   types_sizes(32) = r8_ts;          types_names(32) = "zminb"
!!$      types(33) = H5T_NATIVE_DOUBLE;   types_sizes(33) = r8_ts;          types_names(33) = "zmaxb"
!!$
!!$      types(34) = ndims_r8_t;          types_sizes(34) = ndims_r8_ts;    types_names(34) = "dl"
!!$      types(35) = ndims_r8_t;          types_sizes(35) = ndims_r8_ts;    types_names(35) = "idl"
!!$
!!$      types(36) = nxarr_r8_t;          types_sizes(36) = nxarr_r8_ts;    types_names(36) = "x"
!!$      types(37) = nxarr_r8_t;          types_sizes(37) = nxarr_r8_ts;    types_names(37) = "xl"
!!$      types(38) = nxarr_r8_t;          types_sizes(38) = nxarr_r8_ts;    types_names(38) = "xr"
!!$      types(39) = nxarr_r8_t;          types_sizes(39) = nxarr_r8_ts;    types_names(39) = "inv_x"
!!$
!!$      types(40) = nyarr_r8_t;          types_sizes(40) = nyarr_r8_ts;    types_names(40) = "y"
!!$      types(41) = nyarr_r8_t;          types_sizes(41) = nyarr_r8_ts;    types_names(41) = "yl"
!!$      types(42) = nyarr_r8_t;          types_sizes(42) = nyarr_r8_ts;    types_names(42) = "yr"
!!$      types(43) = nyarr_r8_t;          types_sizes(43) = nyarr_r8_ts;    types_names(43) = "inv_y"
!!$
!!$      types(44) = nzarr_r8_t;          types_sizes(44) = nzarr_r8_ts;    types_names(44) = "z"
!!$      types(45) = nzarr_r8_t;          types_sizes(45) = nzarr_r8_ts;    types_names(45) = "zl"
!!$      types(46) = nzarr_r8_t;          types_sizes(46) = nzarr_r8_ts;    types_names(46) = "zr"
!!$      types(47) = nzarr_r8_t;          types_sizes(47) = nzarr_r8_ts;    types_names(47) = "inv_z"
!!$
!!$      types(48) = ndims_i4_t;          types_sizes(48) = ndims_i4_ts;    types_names(48) = "n_b"
!!$
!!$      types(49) = ndims_i8_t;          types_sizes(49) = ndims_i8_ts;    types_names(49) = "off"
!!$
!!$      types(50) = ndims_lohi_i4_t;     types_sizes(50) = ndims_lohi_i4_ts; types_names(50) = "ijkse"
!!$      types(51) = ndims_lohi_i4_t;     types_sizes(51) = ndims_lohi_i4_ts; types_names(51) = "bnd"
!!$
!!$      types(52) = arr3d_r8_t;          types_sizes(52) = arr3d_r8_ts;    types_names(52) = "gc_xdim"
!!$      types(53) = arr3d_r8_t;          types_sizes(53) = arr3d_r8_ts;    types_names(53) = "gc_ydim"
!!$      types(54) = arr3d_r8_t;          types_sizes(54) = arr3d_r8_ts;    types_names(54) = "gc_zdim"
!!$
!!$      types_names(55) = dname(WA)
!!$      if (associated(cg%wa)) then
!!$         types(55) = arr3d_r8_t;          types_sizes(55) = arr3d_r8_ts
!!$      else
!!$         types(55) = H5T_NATIVE_INTEGER;  types_sizes(55) = int4_ts
!!$      endif
!!$      types_names(56) = dname(GPOT)
!!$      if (associated(cg%gpot)) then
!!$         types(56) = arr3d_r8_t;          types_sizes(56) = arr3d_r8_ts
!!$      else
!!$         types(56) = H5T_NATIVE_INTEGER;  types_sizes(56) = int4_ts
!!$      endif
!!$      types_names(57) = dname(HGPOT)
!!$      if (associated(cg%hgpot)) then
!!$         types(57) = arr3d_r8_t;          types_sizes(57) = arr3d_r8_ts
!!$      else
!!$         types(57) = H5T_NATIVE_INTEGER;  types_sizes(57) = int4_ts
!!$      endif
!!$      types_names(58) = dname(GP)
!!$      if (associated(cg%gp)) then
!!$         types(58) = arr3d_r8_t;          types_sizes(58) = arr3d_r8_ts
!!$      else
!!$         types(58) = H5T_NATIVE_INTEGER;  types_sizes(58) = int4_ts
!!$      endif
!!$      types_names(59) = dname(SGP)
!!$      if (associated(cg%sgp)) then
!!$         types(59) = arr3d_r8_t;          types_sizes(59) = arr3d_r8_ts
!!$      else
!!$         types(59) = H5T_NATIVE_INTEGER;  types_sizes(59) = int4_ts
!!$      endif
!!$      types_names(60) = dname(SGPM)
!!$      if (associated(cg%sgpm)) then
!!$         types(60) = arr3d_r8_t;          types_sizes(60) = arr3d_r8_ts
!!$      else
!!$         types(60) = H5T_NATIVE_INTEGER;  types_sizes(60) = int4_ts
!!$      endif
!!$      types_names(61) = dname(CS_ISO2)
!!$      if (associated(cg%cs_iso2)) then
!!$         types(61) = arr3d_r8_t;          types_sizes(61) = arr3d_r8_ts
!!$      else
!!$         types(61) = H5T_NATIVE_INTEGER;  types_sizes(61) = int4_ts
!!$      endif
!!$
!!$      types(62) = ndims_arr4d_r8_t;    types_sizes(62) = ndims_arr4d_r8_ts;    types_names(62) = "b"
!!$      types_names(63) = dname(B0)
!!$      if (associated(cg%b0%arr)) then
!!$         types(63) = ndims_arr4d_r8_t;    types_sizes(63) = ndims_arr4d_r8_ts
!!$      else
!!$         types(63) = H5T_NATIVE_INTEGER;  types_sizes(63) = int4_ts
!!$      endif
!!$
!!$      types(64) = u_arr4d_r8_t;        types_sizes(64) = u_arr4d_r8_ts;        types_names(64) = "u"
!!$      types(65) = u_arr4d_r8_t;        types_sizes(65) = u_arr4d_r8_ts;        types_names(65) = "uh"
!!$      types_names(66) = dname(U0)
!!$      if (associated(cg%u0%arr)) then
!!$         types(66) = u_arr4d_r8_t;        types_sizes(66) = u_arr4d_r8_ts
!!$      else
!!$         types(66) = H5T_NATIVE_INTEGER;  types_sizes(66) = int4_ts
!!$      endif
!!$
!!$      offset = 0
!!$      do i = 1, total_no
!!$         call h5tinsert_f(dtype_id, types_names(i),  offset, types(i), error)
!!$         offset = offset + types_sizes(i)
!!$      enddo
!!$
!!$      write(dset_name,'("cg",i4.4)') 1
!!$      call h5dcreate_f(file_id, dset_name, dtype_id, dspace_id, dset_id, error)
!!$
!!$      do i = 1, total_no
!!$         call h5tcreate_f(H5T_COMPOUND_F, types_sizes(i), dmem_id(i), error)
!!$         offset = 0
!!$         call h5tinsert_f(dmem_id(i), types_names(i), offset, types(i), error)
!!$      enddo
!!$
!!$      dims = 1
!!$      call h5dwrite_f(dset_id, dmem_id(1),  int(cg%n_(xdim), kind=4), dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(2),  int(cg%n_(ydim), kind=4), dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(3),  int(cg%n_(zdim), kind=4), dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(4),  int(cg%nxb, kind=4),    dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(5),  int(cg%nyb, kind=4),    dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(6),  int(cg%nzb, kind=4),    dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(7),  int(cg%is, kind=4),     dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(8),  int(cg%ie, kind=4),     dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(9),  int(cg%js, kind=4),     dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(10), int(cg%je, kind=4),     dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(11), int(cg%ks, kind=4),     dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(12), int(cg%ke, kind=4),     dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(13), int(cg%maxxyz, kind=4), dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(14), int(cg%isb, kind=4),    dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(15), int(cg%ieb, kind=4),    dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(16), int(cg%jsb, kind=4),    dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(17), int(cg%jeb, kind=4),    dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(18), int(cg%ksb, kind=4),    dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(19), int(cg%keb, kind=4),    dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(20), cg%dx,     dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(21), cg%dy,     dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(22), cg%dz,     dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(23), cg%idx,    dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(24), cg%idy,    dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(25), cg%idz,    dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(26), cg%dxmn,   dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(27), cg%dvol,   dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(28), cg%fbnd(xdim, LO), dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(29), cg%fbnd(xdim, HI), dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(30), cg%fbnd(ydim, LO), dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(31), cg%fbnd(ydim, HI), dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(32), cg%fbnd(zdim, LO), dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(33), cg%fbnd(zdim, HI), dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(34), cg%dl,     dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(35), cg%idl,    dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(36), cg%x,      dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(37), cg%xl,     dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(38), cg%xr,     dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(39), cg%inv_x,  dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(40), cg%y,      dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(41), cg%yl,     dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(42), cg%yr,     dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(43), cg%inv_y,  dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(44), cg%z,      dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(45), cg%zl,     dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(46), cg%zr,     dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(47), cg%inv_z,  dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(48), int(cg%n_b, kind=4),    dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(49), int(cg%n_b, kind=4),    dims, error, xfer_prp=plist_id) !!!
!!$      call h5dwrite_f(dset_id, dmem_id(50), int(cg%ijkse, kind=4),  dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(51), int(cg%bnd, kind=4),    dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(52), int(cg%gc_xdim, kind=4),dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(53), int(cg%gc_ydim, kind=4),dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(54), int(cg%gc_zdim, kind=4),dims, error, xfer_prp=plist_id)
!!$      if (associated(cg%wa)) then
!!$         call h5dwrite_f(dset_id, dmem_id(55), cg%wa     ,dims, error, xfer_prp=plist_id)
!!$      else
!!$         call h5dwrite_f(dset_id, dmem_id(55), -999_INT4      ,dims, error, xfer_prp=plist_id)
!!$      endif
!!$      if (associated(cg%gpot)) then
!!$         call h5dwrite_f(dset_id, dmem_id(56), cg%gpot   ,dims, error, xfer_prp=plist_id)
!!$      else
!!$         call h5dwrite_f(dset_id, dmem_id(56), -999_INT4      ,dims, error, xfer_prp=plist_id)
!!$      endif
!!$      if (associated(cg%hgpot)) then
!!$         call h5dwrite_f(dset_id, dmem_id(57), cg%hgpot  ,dims, error, xfer_prp=plist_id)
!!$      else
!!$         call h5dwrite_f(dset_id, dmem_id(57), -999_INT4      ,dims, error, xfer_prp=plist_id)
!!$      endif
!!$      if (associated(cg%gp)) then
!!$         call h5dwrite_f(dset_id, dmem_id(58), cg%gp     ,dims, error, xfer_prp=plist_id)
!!$      else
!!$         call h5dwrite_f(dset_id, dmem_id(58), -999_INT4      ,dims, error, xfer_prp=plist_id)
!!$      endif
!!$      if (associated(cg%sgp)) then
!!$         call h5dwrite_f(dset_id, dmem_id(59), cg%sgp    ,dims, error, xfer_prp=plist_id)
!!$      else
!!$         call h5dwrite_f(dset_id, dmem_id(59), -999_INT4      ,dims, error, xfer_prp=plist_id)
!!$      endif
!!$      if (associated(cg%sgpm)) then
!!$         call h5dwrite_f(dset_id, dmem_id(60), cg%sgpm   ,dims, error, xfer_prp=plist_id)
!!$      else
!!$         call h5dwrite_f(dset_id, dmem_id(60), -999_INT4      ,dims, error, xfer_prp=plist_id)
!!$      endif
!!$      if (associated(cg%cs_iso2)) then
!!$         call h5dwrite_f(dset_id, dmem_id(61), cg%cs_iso2,dims, error, xfer_prp=plist_id)
!!$      else
!!$         call h5dwrite_f(dset_id, dmem_id(61), -999_INT4      ,dims, error, xfer_prp=plist_id)
!!$      endif
!!$      call h5dwrite_f(dset_id, dmem_id(62), cg%b,dims, error, xfer_prp=plist_id)
!!$      if (associated(cg%b0%arr)) then
!!$         call h5dwrite_f(dset_id, dmem_id(63), cg%b0%arr     ,dims, error, xfer_prp=plist_id)
!!$      else
!!$         call h5dwrite_f(dset_id, dmem_id(63), -999_INT4      ,dims, error, xfer_prp=plist_id)
!!$      endif
!!$      call h5dwrite_f(dset_id, dmem_id(64), cg%u,  dims, error, xfer_prp=plist_id)
!!$      call h5dwrite_f(dset_id, dmem_id(65), cg%uh%arr, dims, error, xfer_prp=plist_id)
!!$      if (associated(cg%u0%arr)) then
!!$         call h5dwrite_f(dset_id, dmem_id(66), cg%u0%arr     ,dims, error, xfer_prp=plist_id)
!!$      else
!!$         call h5dwrite_f(dset_id, dmem_id(66), -999_INT4      ,dims, error, xfer_prp=plist_id)
!!$      endif
!!$
!!$      call h5dclose_f(dset_id, error)
!!$      do i = 1, total_no
!!$         call h5tclose_f(dmem_id(i),error)
!!$      enddo
!!$      CALL h5sclose_f(dspace_id, error)
!!$      CALL h5tclose_f(dtype_id, error)
!!$      call h5tclose_f(ndims_r8_t, error)
!!$      call h5tclose_f(ndims_i4_t, error)
!!$      call h5tclose_f(ndims_i8_t, error)
!!$      call h5tclose_f(ndims_lohi_i4_t, error)
!!$      call h5tclose_f(nxarr_r8_t, error)
!!$      call h5tclose_f(nyarr_r8_t, error)
!!$      call h5tclose_f(nzarr_r8_t, error)
!!$      call h5tclose_f(arr3d_r8_t, error)
!!$      call h5tclose_f(ndims_arr4d_r8_t, error)
!!$      call h5tclose_f(u_arr4d_r8_t, error)
!!$
!!$      deallocate(types,types_sizes,types_names,dmem_id)
!!$
!!$   end subroutine write_grid_containter

end module common_hdf5
