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

! pulled by HDF5

   use constants, only: singlechar, ndims, dsetnamelen
   use hdf5,      only: HID_T

   implicit none

   private
   public :: init_hdf5, cleanup_hdf5, set_common_attributes, common_shortcuts, write_to_hdf5_v2, set_h5_properties
   public :: hdf_vars, hdf_vars_avail, cancel_hdf_var, d_gname, base_d_gname, d_fc_aname, d_size_aname, &
        d_edge_apname, d_bnd_apname, cg_gname, cg_cnt_aname, cg_lev_aname, cg_size_aname, cg_offset_aname, &
        n_cg_name, dir_pref, cg_ledge_aname, cg_redge_aname, cg_dl_aname, O_OUT, O_RES, STAT_OK, STAT_INV, &
        create_empty_cg_dataset, get_nth_cg, data_gname, output_fname, cg_output, enable_all_hdf_var, &
        dump_announcement, dump_announce_time
#ifdef NBODY
   public :: part_gname, st_gname, pdsets

   character(len=*), parameter :: part_types_gname = "particle_types", part_gname = "particles", st_gname = "stars"
   character(len=dsetnamelen), allocatable, dimension(:) :: pdsets
#endif /* NBODY */
   character(len=dsetnamelen), allocatable, dimension(:), protected :: hdf_vars  !< dataset names for hdf files
   logical,                    allocatable, dimension(:), protected :: hdf_vars_avail
   character(len=*), parameter :: d_gname = "domains", base_d_gname = "base", d_fc_aname = "fine_count", &
        & d_size_aname = "n_d", d_edge_apname = "-edge_position", d_bnd_apname = "-boundary_type", &
        & cg_gname = "grid", cg_cnt_aname = "cg_count", cg_lev_aname = "level", cg_size_aname = "n_b", &
        & cg_offset_aname = "off", cg_ledge_aname = "left_edge", cg_redge_aname = "right_edge", &
        & cg_dl_aname = "dl", data_gname = "data"
   character(len=singlechar), dimension(ndims), parameter :: dir_pref = [ "x", "y", "z" ]

   ! enumerator for 'otype' used in various functions to distinguish different
   ! types of output
   enum, bind(c)
      enumerator :: O_RES
      enumerator :: O_OUT
   end enum

   ! indicate success or failure
   enum, bind(C)
      enumerator :: STAT_OK = 0
      enumerator :: STAT_INV = -1
   end enum

   enum, bind(c)
      enumerator :: cg_le !< index list of left edges from all cgs/procs
      enumerator :: cg_re !< index list of right edges from all cgs/procs
      enumerator :: cg_dl !< index list of cell sizes from all cgs/procs
   end enum

   type :: cg_output
      integer(HID_T), dimension(:), allocatable   :: cg_g_id
      integer(HID_T), dimension(:), allocatable   :: part_g_id
      integer(HID_T), dimension(:), allocatable   :: st_g_id
      integer(HID_T), dimension(:,:), allocatable :: dset_id
      integer(HID_T), dimension(:,:), allocatable :: pdset_id
      integer(HID_T)                              :: xfer_prp
      integer(kind=4), allocatable, dimension(:)  :: offsets
      integer(kind=4), allocatable, dimension(:)  :: cg_src_p
      integer(kind=4), allocatable, dimension(:)  :: cg_src_n
      integer(kind=4)                             :: tot_cg_n
   contains
      procedure :: init => initialize_write_cg
      procedure :: clean => finalize_write_cg
   end type cg_output

contains

!>
!! \brief Procedure initializing HDF5 module
!<

   subroutine init_hdf5(vars)

      use constants,      only: dsetnamelen, singlechar
      use dataio_pub,     only: warn
      use fluids_pub,     only: has_ion, has_dst, has_neu
      use global,         only: cc_mag
      use mpisetup,       only: master
#ifdef COSM_RAYS
      use cr_data,        only: cr_names, cr_spectral
      use dataio_pub,     only: msg
#endif /* COSM_RAYS */
#ifdef CRESP
      use cr_data,        only: icr_spc
      use initcosmicrays, only: ncrb, nspc
#endif /* CRESP */

      implicit none

      character(len=dsetnamelen), dimension(:), intent(in) :: vars  !< quantities to be plotted, see dataio::vars

      integer                                              :: i
      character(len=singlechar)                            :: fc, ord
      character(len=dsetnamelen)                           :: aux
#ifdef COSM_RAYS
      integer                                              :: k, ke
#endif /* COSM_RAYS */
#ifdef CRESP
      integer                                              :: ks
#endif /* CRESP */

      if (.not. allocated(hdf_vars)) allocate(hdf_vars(0))

      do i = lbound(vars, 1), ubound(vars, 1)
         select case (trim(vars(i)))
            case ('')
            case ('dens')
               if (has_dst) call append_var('dend')
               if (has_neu) call append_var('denn')
               if (has_ion) call append_var('deni')
            case ('velx')
               if (has_dst) call append_var('vlxd')
               if (has_neu) call append_var('vlxn')
               if (has_ion) call append_var('vlxi')
            case ('vely')
               if (has_dst) call append_var('vlyd')
               if (has_neu) call append_var('vlyn')
               if (has_ion) call append_var('vlyi')
            case ('velz')
               if (has_dst) call append_var('vlzd')
               if (has_neu) call append_var('vlzn')
               if (has_ion) call append_var('vlzi')
            case ('momx')
               if (has_dst) call append_var('momxd')
               if (has_neu) call append_var('momxn')
               if (has_ion) call append_var('momxi')
            case ('momy')
               if (has_dst) call append_var('momyd')
               if (has_neu) call append_var('momyn')
               if (has_ion) call append_var('momyi')
            case ('momz')
               if (has_dst) call append_var('momzd')
               if (has_neu) call append_var('momzn')
               if (has_ion) call append_var('momzi')
            case ('ener')
               if (has_neu) call append_var('enen')
               if (has_ion) call append_var('enei')
            case ('ethr')
               if (has_neu) call append_var('ethn')
               if (has_ion) call append_var('ethi')
            case ('pres')
               if (has_neu) call append_var('pren')
               if (has_ion) call append_var('prei')
            case ('temp')
               if (has_neu) call append_var('temn')
               if (has_ion) call append_var('temi')
            case ("divb", "divB")
               if (cc_mag) then
                  call append_var("divbc")
               else
                  call append_var("divbf")
               endif
            case ("divb4", "divb6", "divb8")
               if (cc_mag) then
                  fc = "c"
               else
                  fc = "f"
               endif
               read(vars(i), '(a4,a1)') aux, ord
               write(aux, '(3a)') "divb", fc, ord
               call append_var(aux)
            case ("divb_diml")
               call append_var("divb_norm")
#ifdef COSM_RAYS
            case ('encr')
               do k = 1, size(cr_names)
                  if (cr_spectral(k)) cycle
                  if (len(trim(cr_names(k))) > 0) then
                     write(aux,'(a3,a)') 'cr_', trim(cr_names(k))
                     call append_var(aux)
                  else
                     ke = k - count(cr_spectral)
                     if (ke <= 99) then
                        write(aux,'(A2,I2.2)') 'cr', ke
                        call append_var(aux)
                     else
                        write(msg, '(a,i3)')"[common_hdf5:init_hdf5] Cannot create name for CR energy component #", ke
                        call warn(msg)
                     endif
                  endif
               enddo
#endif /* COSM_RAYS */
#ifdef CRESP
            case ('cren') !< CRESP number density fields ! TODO change labels, allow saving just the picked species
               do ks = 1, nspc
                  do k = 1, ncrb
                     if (k<=99) then
                        write(aux,'(A3,a,A1,I2.2)') 'cr_', trim(cr_names(icr_spc(ks))), 'n', k
                        call append_var(aux)
                     else
                        write(msg, '(a,i3)')"[common_hdf5:init_hdf5] Cannot create name for CRESP number density component #", k
                        call warn(msg)
                     endif
                  enddo
                  do k = lbound(vars, 1), ubound(vars, 1)
                     if (vars(k) .eq. 'cree') exit
                     if (k .eq. ubound(vars, 1)) then
                        write(msg, '(a)')"[common_hdf5:init_hdf5] CRESP 'cren' field created, but 'cree' not defined: reconstruction of spectrum from hdf files requires both."
                        call warn(msg)
                     endif
                  enddo
               enddo
            case ('cree') !< CRESP energy density fields ! TODO change labels, allow saving just the picked species
               do ks = 1, nspc
                  do k = 1, ncrb
                     if (k<=99) then
                        write(aux,'(A3,a,A1,I2.2)') 'cr_', trim(cr_names(icr_spc(ks))), 'e', k
                        call append_var(aux)
                     else
                        write(msg, '(a,i3)')"[common_hdf5:init_hdf5] Cannot create name for CRESP energy density component #", k
                        call warn(msg)
                     endif
                  enddo
                  do k = lbound(vars, 1), ubound(vars, 1)
                     if (vars(k) .eq. 'cren') exit
                     if (k .eq. ubound(vars, 1)) then
                        write(msg, '(a)')"[common_hdf5:init_hdf5] CRESP 'cree' field created, but 'cren' not defined: reconstruction of spectrum from hdf files requires both."
                        call warn(msg)
                     endif
                  enddo
               enddo
            case ('cref') !< CRESP distribution function
               do k = 1, ncrb+1
                  if (k<=99) then
                     write(aux,'(A4,I2.2)') 'cref', k
                     call append_var(aux)
                  else
                     write(msg, '(a,i3)')"[common_hdf5:init_hdf5] Cannot create name for CRESP distribution function component #", k
                     call warn(msg)
                  endif
               enddo
            case ('crep') !< CRESP cutoff momenta
               do k = 1, 2
                  if (k<=99) then
                     write(aux,'(A4,I2.2)') 'crep', k
                     call append_var(aux)
                  else
                     write(msg, '(a,i3)')"[common_hdf5:init_hdf5] Cannot create name for CRESP cutoff momentum component #", k
                     call warn(msg)
                  endif
               enddo
            case ('creq') !< CRESP spectrum index
               do k = 1, ncrb
                  if (k<=99) then
                     write(aux,'(A4,I2.2)') 'creq', k
                     call append_var(aux)
                  else
                     write(msg, '(a,i3)')"[common_hdf5:init_hdf5] Cannot create name for CRESP spectrum index component #", k
                     call warn(msg)
                  endif
               enddo
#endif /* CRESP */
            case default
               if (.not. has_ion .and. (any(trim(vars(i)) == ["deni", "vlxi", "vlyi", "vlzi", "enei", "ethi", "prei"]) .or. any(trim(vars(i)) == ["momxi", "momyi", "momzi"]))) then
                  if (master) call warn("[common_hdf5:init_hdf5] Cannot safely use plot variable '" // trim(vars(i)) // "' without ionized fluid")
               else if (.not. has_neu .and. (any(trim(vars(i)) == ["denn", "vlxn", "vlyn", "vlzn", "enen", "ethn", "pren"]) .or. any(trim(vars(i)) == ["momxn", "momyn", "momzn"]))) then
                  if (master) call warn("[common_hdf5:init_hdf5] Cannot safely use plot variable '" // trim(vars(i)) // "' without neutral fluid")
               else if (.not. has_dst .and. (any(trim(vars(i)) == ["dend", "vlxd", "vlyd", "vlzd"]) .or. any(trim(vars(i)) == ["momxd", "momyd", "momzd"]))) then
                  if (master) call warn("[common_hdf5:init_hdf5] Cannot safely use plot variable '" // trim(vars(i)) // "' without dust fluid")
               else
                  call append_var(vars(i)) ! all other known and unknown field descriptions
               endif
         end select
      enddo

      allocate(hdf_vars_avail(size(hdf_vars)))
      call enable_all_hdf_var

   contains

      subroutine append_var(n)

         use dataio_pub, only: warn, die
         use mpisetup,   only: master

         implicit none

         character(len=*), intent(in) :: n

         character(len=dsetnamelen), allocatable, dimension(:) :: tmp

         if (len_trim(n) < 1) then
            if (master) call warn("[common_hdf5:init_hdf5:append_var] empty name")
            return
         endif

         if (.not. allocated(hdf_vars)) then
            call die("[common_hdf5:init_hdf5:append_var] hdf_vars not allocated")
         else
            if (.not. any(trim(n) == hdf_vars)) then
               allocate(tmp(lbound(hdf_vars, dim=1):ubound(hdf_vars, dim=1) + 1))
               tmp(:ubound(hdf_vars, dim=1)) = hdf_vars
               call move_alloc(from=tmp, to=hdf_vars)
               hdf_vars(ubound(hdf_vars, dim=1)) = trim(n)
            else
               if (master) call warn("[common_hdf5:init_hdf5:append_var] duplicated name: '" // trim(n) // "'")
            endif
         endif

      end subroutine append_var

   end subroutine init_hdf5

!> \brief Procedure finalizing HDF5 module

   subroutine cleanup_hdf5

      implicit none

      if (allocated(hdf_vars))       deallocate(hdf_vars)
      if (allocated(hdf_vars_avail)) deallocate(hdf_vars_avail)
#ifdef NBODY
      if (allocated(pdsets))         deallocate(pdsets)
#endif /* NBODY */

   end subroutine cleanup_hdf5

!> \brief Mark a plot variable as faulty (don't spam warnings unnecessarily)

   subroutine cancel_hdf_var(var)

      use constants, only: I_ONE

      implicit none

      character(len=dsetnamelen), intent(in) :: var
      integer(kind=4)                        :: i

      do i = I_ONE, size(hdf_vars, kind=4) !> \todo: introduce integer(kind=4), parameter :: first_hdf_var = 1 in constants.h and use it everywhere instead
         if (hdf_vars(i) == var) hdf_vars_avail(i) = .false.
      enddo

   end subroutine cancel_hdf_var

!> \brief Mark all plot variables as good

   subroutine enable_all_hdf_var

      implicit none

      if (size(hdf_vars_avail) > 0) hdf_vars_avail = .true.

   end subroutine enable_all_hdf_var

!-----------------------------------------------------------------------------
!>
!! \brief decode some useful indices from variable name, if possible
!<
   subroutine common_shortcuts(var, fl_dni, i_xyz)

      use constants,  only: dsetnamelen, singlechar, INT4
      use dataio_pub, only: warn
      use fluidindex, only: flind
      use fluids_pub, only: has_ion, has_dst, has_neu
      use fluidtypes, only: component_fluid

      implicit none

      character(len=dsetnamelen),      intent(in)    :: var
      class(component_fluid), pointer, intent(inout) :: fl_dni
      integer(kind=4),                 intent(out)   :: i_xyz

      character(len=singlechar)                      :: dc

      nullify(fl_dni)
      if (any([ "den", "vlx", "vly", "vlz", "ene" ] == var(1:3))) then
         select case (var(4:4))
            case ("d")
               if (has_dst) then
                  fl_dni => flind%dst
               else
                  call warn("[common_hdf5:common_shortcuts] cannot assign fluid to " // trim(var) // "' because we have no dust fluid")
               endif
            case ("n")
               if (has_neu) then
                  fl_dni => flind%neu
               else
                  call warn("[common_hdf5:common_shortcuts] cannot assign fluid to " // trim(var) // "' because we have no neutral fluid")
               endif
            case ("i")
               if (has_ion) then
                  fl_dni => flind%ion
               else
                  call warn("[common_hdf5:common_shortcuts] cannot assign fluid to " // trim(var) // "' because we have no ionized fluid")
               endif
            case default
               call warn("[common_hdf5:common_shortcuts] cannot assign fluid to '" // trim(var) // "'")
         end select
      else if (any([ "momx", "momy", "momz" ] == var(1:4))) then
         select case (var(5:5))
            case ("d")
               if (has_dst) then
                  fl_dni => flind%dst
               else
                  call warn("[common_hdf5:common_shortcuts] cannot assign fluid to " // trim(var) // "' because we have no dust fluid")
               endif
            case ("n")
               if (has_neu) then
                  fl_dni => flind%neu
               else
                  call warn("[common_hdf5:common_shortcuts] cannot assign fluid to " // trim(var) // "' because we have no neutral fluid")
               endif
            case ("i")
               if (has_ion) then
                  fl_dni => flind%ion
               else
                  call warn("[common_hdf5:common_shortcuts] cannot assign fluid to " // trim(var) // "' because we have no ionized fluid")
               endif
            case default
               call warn("[common_hdf5:common_shortcuts] cannot assign fluid to '" // trim(var) // "'")
         end select
      endif

      i_xyz = huge(1_INT4)
      if (var(1:2) == "vl") then
         dc = var(3:3)
      else if (var(1:3) == "mag" .or. var(1:3) == "mom") then
         dc = var(4:4)
      else
         dc = '_'
      endif
      if (any([ "x", "y", "z" ] == dc)) i_xyz = ichar(dc, kind=4) - ichar("x", kind=4)

   end subroutine common_shortcuts

   function tight_chararray(arrin, ndim, maxn) result(arrout)

      implicit none

      integer,                           intent(in) :: ndim
      integer(kind=4),                   intent(in) :: maxn
      character(len=*), dimension(ndim), intent(in) :: arrin
      character(len=maxn), dimension(ndim)          :: arrout

      arrout = arrin(:)(:maxn)

   end function tight_chararray
!>
!! \brief This routine writes all attributes that are common to restart and output files.
!!
!! \details Write real, integer and character attributes. Store contents of problem.par and env files.
!! Other common elements may also be moved here.
!!
!! \ToDo figure out if it is of use for us:
!! http://computation.llnl.gov/projects/floating-point-compression/zfp-and-derivatives
!! https://github.com/LLNL/H5Z-ZFP
!<
   subroutine set_common_attributes(filename)

      use constants,         only: I_ONE
      use dataio_pub,        only: maxenvlen, maxparlen, use_v2_io, parfile, parfilelines, gzip_level
      use dataio_user,       only: user_attrs_wr, user_attrs_pre
      use hdf5,              only: HID_T, SIZE_T, HSIZE_T, H5F_ACC_TRUNC_F, H5T_NATIVE_CHARACTER, H5Z_FILTER_DEFLATE_F, H5P_DATASET_CREATE_F, &
         & h5open_f, h5fcreate_f, h5fclose_f, H5Zfilter_avail_f, H5Pcreate_f, H5Pset_deflate_f, H5Pset_chunk_f, h5tcopy_f, h5tset_size_f, &
         & h5screate_simple_f, H5Dcreate_f, H5Dwrite_f, H5Dclose_f, H5Sclose_f, H5Tclose_f, H5Pclose_f, h5close_f
      use mpisetup,          only: master, slave
      use version,           only: env, nenv
#ifdef CRESP
      use cresp_io,          only: write_cresp_to_restart
#endif /* CRESP */
#ifdef RANDOMIZE
      use randomization,     only: write_current_seed_to_restart
#endif /* RANDOMIZE */
#ifdef SN_SRC
      use snsources,         only: write_snsources_to_restart
#endif /* SN_SRC */
#if defined(MULTIGRID) && defined(SELF_GRAV)
      use multigrid_gravity, only: write_oldsoln_to_restart
#endif /* MULTIGRID && SELF_GRAV */

      implicit none

      character(len=*), intent(in)   :: filename        !< HDF File name

      integer(HID_T)                 :: file_id         !< File identifier
      integer(HID_T)                 :: type_id, dspace_id, dset_id, prp_id
      integer(HSIZE_T), dimension(1) :: dimstr
      logical(kind=4)                :: Z_avail         !< Z_avail perhaps should be of type integer(HID_T)
      integer(kind=4)                :: error           !< error perhaps should be of type integer(HID_T)

      ! \ToDo Set up a stack of routines registered by appropriate modules
      if (associated(user_attrs_pre)) call user_attrs_pre

      if (master) then
         call h5open_f(error)
         call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
      endif

#if defined(MULTIGRID) && defined(SELF_GRAV)
      call write_oldsoln_to_restart(file_id)  ! Old solution fields have to be selectively marked for restart also on slaves
#endif /* MULTIGRID && SELF_GRAV */

      if (slave) return ! This data need not be written in parallel.

      if (use_v2_io) then
         call set_common_attributes_v2(file_id)
      else
         call set_common_attributes_v1(file_id)
      endif

      ! Store a compressed copy of the problem.par file.
      dimstr = [parfilelines]
      call H5Zfilter_avail_f(H5Z_FILTER_DEFLATE_F, Z_avail, error)
      ! call H5Zget_filter_info_f ! everything should be always fine for gzip
      call H5Pcreate_f(H5P_DATASET_CREATE_F, prp_id, error)
      if (Z_avail) then
         call H5Pset_deflate_f(prp_id, gzip_level, error)
         call H5Pset_chunk_f(prp_id, I_ONE, dimstr, error)
      endif
      call H5Tcopy_f(H5T_NATIVE_CHARACTER, type_id, error)
      call H5Tset_size_f(type_id, int(maxparlen, SIZE_T), error)
      call H5Screate_simple_f(I_ONE, dimstr, dspace_id, error)
      call H5Dcreate_f(file_id, "problem.par", type_id,  dspace_id, dset_id, error, dcpl_id = prp_id)
      call H5Dwrite_f(dset_id, type_id, tight_chararray(parfile, parfilelines, maxparlen), dimstr, error)
      call H5Dclose_f(dset_id, error)
      call H5Sclose_f(dspace_id, error)

      ! Store a compressed copy of the piernik.def file and Id lines from source files.
      ! We recycle type_id and prp_id, so we don't close them yet.
      dimstr = [nenv]
      if (Z_avail) call H5Pset_chunk_f(prp_id, I_ONE, dimstr, error)
      call H5Tset_size_f(type_id, int(maxenvlen, SIZE_T), error)
      call H5Screate_simple_f(I_ONE, dimstr, dspace_id, error)
      call H5Dcreate_f(file_id, "env", type_id,  dspace_id, dset_id, error, dcpl_id = prp_id)
      call H5Dwrite_f(dset_id, type_id, tight_chararray(env, nenv, maxenvlen), dimstr, error)
      call H5Dclose_f(dset_id, error)
      call H5Sclose_f(dspace_id, error)
      call H5Tclose_f(type_id, error)
      call H5Pclose_f(prp_id, error)

      ! \ToDo Set up a stack of routines registered by appropriate modules
#ifdef CRESP
      call write_cresp_to_restart(file_id)
#endif /* CRESP */
#ifdef RANDOMIZE
      call write_current_seed_to_restart(file_id)
#endif /* RANDOMIZE */
#ifdef SN_SRC
      call write_snsources_to_restart(file_id)
#endif /* SN_SRC */
      if (associated(user_attrs_wr)) call user_attrs_wr(file_id)

      call h5fclose_f(file_id, error)
      call h5close_f(error)

      ! only master process exits here

   end subroutine set_common_attributes

!>
!! \brief Common attributes for v2 files
!!
!! The rr1 marks critical attributes that are read by read_restart_hdf5 and compared against value read from the
!! problem.par file.
!! The rr2 marks runtime values that are read by read_restart_hdf5 and assigned to something in the code.
!> \todo Set up an universal table(s) of attribute names for use by both set_common_attributes and read_restart_hdf5.
!! Provide indices for critical attributes (rr1) and for runtime attributes (rr2).
!!
!! \ToDo convert to use "call set_attr", like in set_common_attributes_v2, then despaghettify.
!<

   subroutine set_common_attributes_v1(file_id)

      use cg_level_finest, only: finest
      use constants,       only: cbuff_len, xdim, ydim, zdim, I_ONE
      use dataio_pub,      only: require_problem_IC, piernik_hdf5_version, problem_name, run_id, last_hdf_time, &
         &                       last_res_time, last_tsl_time, last_log_time, nres, nhdf, domain_dump
      use domain,          only: dom
      use fluidindex,      only: flind
      use global,          only: t, dt, nstep
      use hdf5,            only: HID_T, SIZE_T
      use h5lt,            only: h5ltset_attribute_double_f, h5ltset_attribute_int_f, h5ltset_attribute_string_f
      use mass_defect,     only: magic_mass
      use timestep_pub,    only: c_all_old, cfl_c, stepcfl
      use units,           only: cm, gram, sek, kelvin, miu0

      implicit none

      integer(HID_T), intent(in)                   :: file_id       !< File identifier

      integer(kind=4)                              :: fe
      integer(SIZE_T)                              :: i
      integer(kind=4)                              :: error         !< error perhaps should be of type integer(HID_T)
      integer, parameter                           :: buf_len = 50
      integer(SIZE_T), parameter                   :: bufsize = I_ONE
      integer(SIZE_T),          dimension(buf_len) :: rbuffer_size
      integer(kind=4),          dimension(buf_len) :: ibuffer
      real,                     dimension(buf_len) :: rbuffer
      character(len=cbuff_len), dimension(buf_len) :: ibuffer_name = ''
      character(len=cbuff_len), dimension(buf_len) :: rbuffer_name = ''

      rbuffer_size = bufsize
      rbuffer(1)   = t                       ; rbuffer_name(1)   = "time" !rr2
      rbuffer(2)   = dt                      ; rbuffer_name(2)   = "timestep" !rr2
      rbuffer(3:4) = dom%edge(xdim, :)       ; rbuffer_name(3:4) = [ "xmin", "xmax" ] !rr1
      rbuffer(5:6) = dom%edge(ydim, :)       ; rbuffer_name(5:6) = [ "ymin", "ymax" ] !rr1
      rbuffer(7:8) = dom%edge(zdim, :)       ; rbuffer_name(7:8) = [ "zmin", "zmax" ] !rr1
      rbuffer(9)   = piernik_hdf5_version    ; rbuffer_name(9)   = "piernik" !rr1, rr2
      rbuffer(10)  = last_log_time           ; rbuffer_name(10)  = "last_log_time" !rr2
      rbuffer(11)  = last_tsl_time           ; rbuffer_name(11)  = "last_tsl_time" !rr2
      rbuffer(12)  = last_hdf_time           ; rbuffer_name(12)  = "last_hdf_time" !rr2
      rbuffer(13)  = last_res_time           ; rbuffer_name(13)  = "last_res_time" !rr2
      rbuffer(14)  = -99999.9                ; rbuffer_name(14)  = "last_plt_time" !rr2 ! FIXME
      rbuffer(15)  = c_all_old               ; rbuffer_name(15)  = "c_all_old" !rr2
      rbuffer(16)  = stepcfl                 ; rbuffer_name(16)  = "stepcfl" !rr2
      rbuffer(17)  = cfl_c                   ; rbuffer_name(17)  = "cfl_c" !rr2
      rbuffer(18)  = cm                      ; rbuffer_name(18)  = "cm" !rr2
      rbuffer(19)  = gram                    ; rbuffer_name(19)  = "gram" !rr2
      rbuffer(20)  = sek                     ; rbuffer_name(20)  = "sek" !rr2
      rbuffer(21)  = miu0                    ; rbuffer_name(21)  = "miu0" !rr2
      rbuffer(22)  = kelvin                  ; rbuffer_name(22)  = "kelvin" !rr2
      rbuffer_size(23) = flind%fluids
      rbuffer(23:22+rbuffer_size(23)) = magic_mass ; rbuffer_name(23:22+rbuffer_size(23)) = "magic_mass" !rr2

      ibuffer(1)   = nstep                   ; ibuffer_name(1)   = "nstep" !rr2
      ibuffer(2)   = nres                    ; ibuffer_name(2)   = "nres" !rr2
      ibuffer(3)   = nhdf                    ; ibuffer_name(3)   = "nhdf" !rr2
      ibuffer(4)   = -1                      ; ibuffer_name(4)   = "nimg" !rr2 !FIXME
      !>
      !! \todo check if finest is complete, if not then find finest complete level
      !! (see data_hdf5::h5_write_to_single_file_v1)
      !<
      ibuffer(5:7) = int(finest%level%l%n_d(:), kind=4) ; ibuffer_name(5:7) = [ "nxd", "nyd", "nzd" ] !rr1
      ibuffer(8)   = dom%nb                  ; ibuffer_name(8)   = "nb"
      ibuffer(9)   = require_problem_IC      ; ibuffer_name(9)   = "require_problem_IC" !rr2

      i = 1
      do while (rbuffer_name(i) /= "")
         call h5ltset_attribute_double_f(file_id, "/", rbuffer_name(i), rbuffer(i:i-I_ONE+rbuffer_size(i)), rbuffer_size(i), error)
         i = i + rbuffer_size(i)
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
!<

   subroutine set_common_attributes_v2(file_id)

      use constants,          only: I_ONE
      use dataio_pub,         only: require_problem_IC, piernik_hdf5_version2, problem_name, run_id, last_hdf_time, &
         &                          last_res_time, last_log_time, last_tsl_time, nres, nhdf, domain_dump
      use global,             only: t, dt, nstep
      use hdf5,               only: HID_T
      use mass_defect,        only: magic_mass
      use set_get_attributes, only: set_attr
      use timestep_pub,       only: c_all_old, cfl_c, stepcfl

      implicit none

      integer(HID_T), intent(in)                   :: file_id       !< File identifier

      ! real attributes
      call set_attr(file_id, "time",          [t                     ]) !rr2
      call set_attr(file_id, "timestep",      [dt                    ]) !rr2
      call set_attr(file_id, "piernik",       [piernik_hdf5_version2 ]) !rr1, rr2
      call set_attr(file_id, "last_log_time", [last_log_time         ]) !rr2
      call set_attr(file_id, "last_tsl_time", [last_tsl_time         ]) !rr2
      call set_attr(file_id, "last_hdf_time", [last_hdf_time         ]) !rr2
      call set_attr(file_id, "last_res_time", [last_res_time         ]) !rr2
      call set_attr(file_id, "last_plt_time", [-99999.99999          ]) !rr2 !FIXME
      call set_attr(file_id, "c_all_old",     [c_all_old             ]) !rr2
      call set_attr(file_id, "stepcfl",       [stepcfl               ]) !rr2
      call set_attr(file_id, "cfl_c",         [cfl_c                 ]) !rr2
      call set_attr(file_id, "magic_mass",     magic_mass)

      ! integer attributes
      call set_attr(file_id, "nstep",              [nstep                 ]) !rr2
      call set_attr(file_id, "nres",               [nres                  ]) !rr2
      call set_attr(file_id, "nhdf",               [nhdf                  ]) !rr2
      call set_attr(file_id, "nimg",               [-I_ONE                ]) !rr2 !FIXME
      call set_attr(file_id, "require_problem_IC", [require_problem_IC    ]) !rr2
      !> \todo  add number of pieces in the restart point/data dump

      ! string attributes
      call set_attr(file_id, "problem_name", [trim(problem_name)]) !rr2
      call set_attr(file_id, "domain",       [trim(domain_dump) ]) !rr2
      call set_attr(file_id, "run_id",       [trim(run_id)      ]) !rr2

   end subroutine set_common_attributes_v2

!> \brief Generate numbered cg group name
   function n_cg_name(g)

      use constants, only: dsetnamelen

      implicit none

      integer(kind=4), intent(in) :: g !< group number
      character(len=dsetnamelen) :: n_cg_name

      write(n_cg_name,'(2a,i10.10)')trim(cg_gname), "_", g-1

   end function n_cg_name

!>
!! \brief Find a n-th grid container mentioned in cg_desc structure.
!!
!! The loop here has to match the loop in collect_cg_data!
!<

   function get_nth_cg(n) result(cg)

      use cg_leaves,  only: leaves
      use cg_list,    only: cg_list_element
      use dataio_pub, only: die
      use grid_cont,  only: grid_container

      implicit none

      integer(kind=4), intent(in) :: n

      type(grid_container),  pointer :: cg
      type(cg_list_element), pointer :: cgl

      integer :: i

      nullify(cg)
      i = 1
      cgl => leaves%first
      do while (associated(cgl))
         if (i == n) then
            cg => cgl%cg
            exit
         endif
         i = i + 1
         cgl => cgl%nxt
      enddo

      if (.not. associated(cg)) call die("[common_hdf5:get_nth_cg] cannot find n-th cg")

   end function get_nth_cg


!> \brief Create an empty double precision dataset of given dimensions. Use compression if available.
   subroutine create_empty_cg_dataset(cg_g_id, name, ddims, Z_avail, otype)

      use dataio_pub, only: enable_compression, gzip_level, die, h5_64bit
      use hdf5,       only: HID_T, HSIZE_T, H5P_DATASET_CREATE_F, H5T_NATIVE_REAL, H5T_NATIVE_DOUBLE, &
           &                h5dcreate_f, h5dclose_f, h5screate_simple_f, h5sclose_f, h5pcreate_f, h5pclose_f, h5pset_deflate_f, &
           &                h5pset_shuffle_f, h5pset_chunk_f

      implicit none

      integer(HID_T),                 intent(in) :: cg_g_id !< group id where to create the dataset
      character(len=*),               intent(in) :: name    !< name
      integer(HSIZE_T), dimension(:), intent(in) :: ddims   !< dimensionality
      logical(kind=4),                intent(in) :: Z_avail !< can use compression?
      integer(kind=4),                intent(in) :: otype   !< output type

      integer(HID_T)                             :: prp_id, filespace, dset_id, dtype
      integer(kind=4)                            :: error   !< error perhaps should be of type integer(HID_T)

      call h5pcreate_f(H5P_DATASET_CREATE_F, prp_id, error)
      if (enable_compression .and. Z_avail) then
         call h5pset_shuffle_f(prp_id, error)
         call h5pset_deflate_f(prp_id, gzip_level, error)
         call h5pset_chunk_f(prp_id, size(ddims, kind=4), ddims, error)
      endif

      if (otype == O_RES) then
         dtype = H5T_NATIVE_DOUBLE
      else if (otype == O_OUT) then
         if (h5_64bit) then
            dtype = H5T_NATIVE_DOUBLE
         else
            dtype = H5T_NATIVE_REAL
         endif
      else
         call die("[common_hdf5:create_empty_cg_dataset] Unknown output time")
      endif

      call h5screate_simple_f(size(ddims, kind=4), ddims, filespace, error)
      call h5dcreate_f(cg_g_id, name, dtype, filespace, dset_id, error, dcpl_id = prp_id)
      call h5dclose_f(dset_id, error)
      call h5sclose_f(filespace, error)
      call h5pclose_f(prp_id, error)

   end subroutine create_empty_cg_dataset

!! \brief Write a multi-file, multi-domain HDF5 file
!!
!! \details There are three approaches to be implemented:
!! - Single-file, serial I/O. The easiest way. Master writes everything, slaves send their data to the master. Does not
!!   take advantage of parallel filesystems. Best choice for non-parallel filesystems.
!! - Multi-file, serial I/O. An extension of the above approach. Selected processes (can be all) write to their files,
!!   other processes send them their data.
!!   Can take advantage of parallel filesystems. Can use local scratch filesystems. Requires additional work on reading.
!! - Single-file, parallel I/O. The most ambitious approach. Selected processes (can be all) write to the files, other
!!   processes send them their data.
!!   Can take advantage of parallel filesystems. Currently does not allow for compression during write.
!!   Requires a lot of pseudo-collective operations. The "flexible PHDF5" would simplify the code, but it needs to be
!!   implemented.first.
!!
!! \warning Partial implementation: Single-file, serial I/O works for non-AMR setups.
!!
!! \param create_empty_cg_datasets
!!    Function responsible for creating empty datasets, called by master
!! \param write_cg_to_hdf5
!!    Function that performs actual I/O, called by all
!<

   subroutine write_to_hdf5_v2(filename, otype, create_empty_cg_datasets, write_cg_to_hdf5)

      use cg_leaves,    only: leaves
      use constants,    only: cwdlen, dsetnamelen, xdim, zdim, ndims, I_ONE, I_TWO, I_THREE, I_FOUR, INT4, LO, HI, &
         &                    GEO_XYZ, GEO_RPZ
      use dataio_pub,   only: die, nproc_io, can_i_write, domain_dump, msg
      use domain,       only: dom
      use gdf,          only: gdf_create_format_stamp, gdf_create_simulation_parameters, gdf_create_root_datasets, &
         &                    gdf_root_datasets_t, gdf_parameters_t, GDF_CARTESIAN, GDF_POLAR
      use global,       only: t
      use hdf5,         only: HID_T, H5F_ACC_RDWR_F, H5P_FILE_ACCESS_F, H5P_GROUP_ACCESS_F, H5Z_FILTER_DEFLATE_F, &
         &                    h5open_f, h5close_f, h5fopen_f, h5fclose_f, h5gcreate_f, h5gopen_f, h5gclose_f, h5pclose_f, &
         &                    h5zfilter_avail_f
      use helpers_hdf5, only: create_attribute!, create_corefile
      use MPIF,         only: MPI_INTEGER, MPI_INTEGER8, MPI_STATUS_IGNORE, MPI_REAL8, MPI_COMM_WORLD
      use MPIFUN,       only: MPI_Allgather, MPI_Recv, MPI_Send
      use mpisetup,     only: FIRST, LAST, master, err_mpi, piernik_MPI_Bcast
#ifdef NBODY
      use constants,    only: I_FIVE, I_SIX
      !use particle_utils, only: count_all_particles
#else /* !NBODY */
      use constants,    only: I_ZERO
#endif /* NBODY */

      implicit none

      character(len=cwdlen), intent(in) :: filename         !< Name of the HDF5 file
      integer(kind=4),       intent(in) :: otype            !< Output type (restart, data)
      interface
         !>
         !! Function responsible for creating empty datasets, called by master
         !<
         subroutine create_empty_cg_datasets(cgl_g_id, cg_n_b, cg_n_o, Z_avail, n_part, st_g_id)
            use hdf5, only: HID_T
            implicit none
            integer(HID_T),                intent(in) :: cgl_g_id
            integer(kind=4), dimension(:), intent(in) :: cg_n_b
            integer(kind=4), dimension(:), intent(in) :: cg_n_o
            logical(kind=4),               intent(in) :: Z_avail
            integer(kind=8),               intent(in) :: n_part
            integer(HID_T),                intent(in) :: st_g_id
         end subroutine create_empty_cg_datasets

         !>
         !! Function that performs actual I/O, called by all
         !<
         subroutine write_cg_to_hdf5(cgl_g_id, cg_n, cg_all_n_b, cg_all_n_o)
            use hdf5, only: HID_T
            implicit none
            integer(HID_T),                           intent(in) :: cgl_g_id
            integer(kind=4), dimension(:),   pointer, intent(in) :: cg_n
            integer(kind=4), dimension(:,:), pointer, intent(in) :: cg_all_n_b
            integer(kind=4), dimension(:,:), pointer, intent(in) :: cg_all_n_o
         end subroutine write_cg_to_hdf5
      end interface

      integer(HID_T)                                :: file_id          !< File identifier
      integer(HID_T)                                :: plist_id         !< Property list identifier
      integer(HID_T)                                :: cgl_g_id         !< cg list identifiers
      integer(HID_T)                                :: cg_g_id          !< cg group identifiers
#ifdef NBODY
      integer(HID_T)                                :: pt_g_id          !< particle_types identifiers
      integer(HID_T)                                :: ptst_g_id        !< particle_types stars identifiers
      integer(HID_T)                                :: part_g_id        !< particles identifiers
#endif /* NBODY */
      integer(HID_T)                                :: st_g_id          !< stars identifiers
      integer(HID_T)                                :: doml_g_id        !< domain list identifier
      integer(HID_T)                                :: dom_g_id         !< domain group identifier
      integer(kind=4)                               :: error, cg_cnt, p !< error perhaps should be of type integer(HID_T)
      integer                                       :: g, i
      integer(kind=4), parameter                    :: tag = I_ONE
      integer(kind=4),  dimension(:),   pointer     :: cg_n             !< offset for cg group numbering
      integer(kind=4),  dimension(:,:), pointer     :: cg_all_n_b       !< sizes of all cg
      integer(kind=4),  dimension(:,:), pointer     :: cg_all_n_o       !< sizes of all cg, expanded by external boundaries
      integer(kind=4),  dimension(:),   pointer     :: cg_rl            !< list of refinement levels from all cgs/procs
      integer(kind=4),  dimension(:,:), pointer     :: cg_n_b           !< list of n_b from all cgs/procs
      integer(kind=4),  dimension(:,:), pointer     :: cg_n_o           !< list of grid dimensions with external guardcells from all cgs/procs
      integer(kind=8),  dimension(:,:), pointer     :: cg_off           !< list of offsets from all cgs/procs
      integer(kind=8),  dimension(:),   pointer     :: cg_npart         !< list of particle count from all cgs/procs
      integer(kind=8),  dimension(:),   pointer     :: cg_pid_max       !< list of maximum particle id from all cgs/procs

      !>
      !! auxiliary array for communication of {cg_le, cg_re, cg_dl} lists
      !<
      real(kind=8), dimension(:,:,:), pointer       :: dbuf
      logical(kind=4)                               :: Z_avail !< .true. if HDF5 was compiled with zlib support
      character(len=dsetnamelen)                    :: d_label
      integer(kind=4)                               :: indx
      real, dimension(LO:HI)                        :: edge
      real, dimension(ndims)                        :: temp

      type(gdf_root_datasets_t)                     :: rd
      type(gdf_parameters_t)                        :: gdf_sp

      ! Create a new file and initialize it

      ! Prepare groups and datasets for grid containers on the master
      allocate(cg_n(FIRST:LAST))
      call MPI_Allgather(leaves%cnt, I_ONE, MPI_INTEGER, cg_n, I_ONE, MPI_INTEGER, MPI_COMM_WORLD, err_mpi)
      cg_cnt = sum(cg_n(:))
      allocate(cg_all_n_b(ndims, cg_cnt), cg_all_n_o(ndims, cg_cnt))

      if (master) then
         call rd%init(cg_cnt)

         ! Open the HDF5 file only in master process and create all groups required for cg storage.
         ! Create also all related datasets and attributes. Do not write big datasets yet.

         call h5open_f(error)
         call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
!         call create_corefile(filename, file_id, bstore=.True.)

         if (otype == O_OUT) then
            call gdf_create_format_stamp(file_id)
            call gdf_sp%init()
            gdf_sp%current_time = t
            select case (dom%geometry_type)
               case (GEO_XYZ)
                  gdf_sp%domain_left_edge = dom%edge(:, LO)
                  gdf_sp%domain_right_edge = dom%edge(:, HI)
                  gdf_sp%geometry = GDF_CARTESIAN
               case (GEO_RPZ)
                  gdf_sp%domain_left_edge = dom%edge(:, LO)
                  gdf_sp%domain_right_edge = dom%edge(:, HI)
                  gdf_sp%geometry = GDF_POLAR
               case default
                  write(msg,'(a,i3)') "[common_hdf5:write_to_hdf5_v2] Unknown system of coordinates ", dom%geometry_type
                  call die(msg)
            end select
            gdf_sp%field_ordering = 1
            gdf_sp%boundary_conditions = int([0,0,0,0,0,0], kind=4)  !! \todo fix hardcoded integers
            gdf_sp%refine_by = int([2], kind=8) !! \todo fix hardcoded integers
            gdf_sp%cosmological_simulation = int([0], kind=8) !! \todo fix hardcoded integers
            gdf_sp%dimensionality = int([dom%eff_dim], kind=8)
            gdf_sp%unique_identifier = "Piernik"
            select case (trim(domain_dump))
               case ('phys_domain')
                  gdf_sp%num_ghost_zones = int(0, kind=8)
                  gdf_sp%domain_dimensions = dom%n_d
               case ('full_domain')
                  gdf_sp%num_ghost_zones = int(dom%nb, kind=8)
                  gdf_sp%domain_dimensions = dom%n_d + dom%nb*2
            end select

            call gdf_create_simulation_parameters(file_id, gdf_sp)
            call gdf_sp%cleanup()
         endif
         call h5gcreate_f(file_id, data_gname, cgl_g_id, error)     ! create "/data"
#ifdef NBODY
         call h5gcreate_f(file_id, part_types_gname, pt_g_id, error)     ! create "/particle_types"
         call h5gcreate_f(pt_g_id, st_gname, ptst_g_id, error)     ! create "/particle_types/stars"
#endif /* NBODY */

         call create_attribute(cgl_g_id, cg_cnt_aname, [ cg_cnt ])  ! create "/data/cg_count"

         Z_avail = .false.
         if (nproc_io == 1) call h5zfilter_avail_f(H5Z_FILTER_DEFLATE_F, Z_avail, error)
         call h5zfilter_avail_f(H5Z_FILTER_DEFLATE_F, Z_avail, error)
         !> \todo test it thoroughly before enabling for > 1

         ! Do not assume that the master knows all the lists
         do p = FIRST, LAST
            allocate(cg_rl(cg_n(p)), cg_n_b(cg_n(p), ndims), cg_off(cg_n(p), ndims), cg_npart(cg_n(p)))
#ifdef NBODY
            allocate(cg_pid_max(cg_n(p)))
#else /* !NBODY */
            cg_npart = I_ZERO
#endif /* !NBODY */
            nullify(dbuf)
            if (otype == O_OUT) allocate(dbuf(cg_le:cg_dl, cg_n(p), ndims))
            nullify(cg_n_o)
            if (otype == O_RES) allocate(cg_n_o(cg_n(p), ndims))
            if (p == FIRST) then
               call collect_cg_data(cg_rl, cg_n_b, cg_n_o, cg_off, cg_npart, cg_pid_max, dbuf, otype)
#ifdef NBODY
!               n_part = count_all_particles()
!               pid_max = pid_gen
#endif /* NBODY */
            else
               call MPI_Recv(cg_rl,  size(cg_rl, kind=4),  MPI_INTEGER,  p, tag,         MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
               call MPI_Recv(cg_n_b, size(cg_n_b, kind=4), MPI_INTEGER,  p, tag+I_ONE,   MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
               call MPI_Recv(cg_off, size(cg_off, kind=4), MPI_INTEGER8, p, tag+I_TWO,   MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
               if (otype == O_OUT) &
                  & call MPI_Recv(dbuf,   size(dbuf, kind=4),   MPI_REAL8,    p, tag+I_THREE, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
               if (otype == O_RES) &
                  & call MPI_Recv(cg_n_o, size(cg_n_o, kind=4), MPI_INTEGER,  p, tag+I_FOUR,  MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
#ifdef NBODY
               call MPI_Recv(cg_npart,   size(cg_npart,   kind=4), MPI_INTEGER8, p, tag+I_FIVE, MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
               call MPI_Recv(cg_pid_max, size(cg_pid_max, kind=4), MPI_INTEGER8, p, tag+I_SIX,  MPI_COMM_WORLD, MPI_STATUS_IGNORE, err_mpi)
#endif /* NBODY */
            endif

            do g = 1, cg_n(p)
               call h5gcreate_f(cgl_g_id, n_cg_name(int(sum(cg_n(:p))-cg_n(p)+g, kind=4)), cg_g_id, error) ! create "/data/grid_%08d"

#ifdef NBODY
               call h5gcreate_f(cg_g_id, part_gname, part_g_id, error) ! create "/data/grid_%08d/particles
               call h5gcreate_f(part_g_id, st_gname, st_g_id, error) ! create "/data/grid_%08d/particles/stars"

               if (int(cg_npart(g),   kind=4) /= cg_npart(g))   call die("[common_hdf5:write_to_hdf5_v2] cg_npart needs to be 64-bit")
               if (int(cg_pid_max(g), kind=4) /= cg_pid_max(g)) call die("[common_hdf5:write_to_hdf5_v2] pid_max needs to be 64-bit")

               call create_attribute(st_g_id, "n_part",  [ int(cg_npart(g),   kind=4) ])  ! create "/data/grid_%08d/particles/stars/n_part"
               call create_attribute(st_g_id, "pid_max", [ int(cg_pid_max(g), kind=4) ])  ! create "/data/grid_%08d/particles/stars/pid_max"
#endif /* NBODY */

               call create_attribute(cg_g_id, cg_lev_aname, [ cg_rl(g) ] )                ! create "/data/grid_%08d/level"
               temp = cg_n_b(g, :)
               call create_attribute(cg_g_id, cg_size_aname, temp)                        ! create "/data/grid_%08d/n_b"
               call create_attribute(cg_g_id, cg_offset_aname, int(cg_off(g, :), kind=4)) ! create "/data/grid_%08d/off"

               if (otype == O_OUT) then
                  temp(:) = dbuf(cg_le, g, :)
                  call create_attribute(cg_g_id, cg_ledge_aname, temp)  ! create "/data/grid_%08d/left_edge"
                  temp(:) = dbuf(cg_re, g, :)
                  call create_attribute(cg_g_id, cg_redge_aname, temp)  ! create "/data/grid_%08d/right_edge"
                  temp(:) = dbuf(cg_dl, g, :)
                  call create_attribute(cg_g_id, cg_dl_aname, temp)     ! create "/data/grid_%08d/dl"
               endif

               cg_all_n_b(:, sum(cg_n(:p))-cg_n(p)+g) = cg_n_b(g, :)
               if (associated(cg_n_o)) cg_all_n_o(:, sum(cg_n(:p))-cg_n(p)+g) = cg_n_o(g, :)
               if (otype == O_OUT) then
                  indx = int(sum(cg_n(:p))-cg_n(p)+g, kind=4)
                  rd%grid_level(indx) = cg_rl(g)
                  rd%grid_left_index(:,indx) = cg_off(g,:)
                  rd%grid_parent_id(indx)     = -1
                  rd%grid_particle_count(1,indx) = int(cg_npart(g), kind=4)
               endif

               if (any(cg_off(g, :) > 2.**31)) &
                  & call die("[common_hdf5:write_to_hdf5_v2] large offsets require better treatment")

               if (otype == O_OUT .and. .not. associated(cg_n_o)) allocate(cg_n_o(cg_n(p), ndims))  ! Prevent crash due to passing unallocated array.
               ! Ugliness of the above fix suggests overcomplication somewhere.
               call create_empty_cg_datasets(cg_g_id, cg_n_b(g, :), cg_n_o(g, :), Z_avail, cg_npart(g), st_g_id) !!!!!
#ifdef NBODY
               call h5gclose_f(st_g_id, error)
               call h5gclose_f(part_g_id, error)
#endif /* NBODY */

               call h5gclose_f(cg_g_id, error)
            enddo

            deallocate(cg_rl, cg_n_b, cg_off, cg_npart)
#ifdef NBODY
            deallocate(cg_pid_max)
#endif /* NBODY */
            if (associated(dbuf)) deallocate(dbuf)
            if (associated(cg_n_o)) deallocate(cg_n_o)
         enddo
         rd%grid_dimensions = cg_all_n_b

#ifdef NBODY
         call h5gclose_f(pt_g_id, error)
         call h5gclose_f(ptst_g_id, error)
#endif /* NBODY */

         call h5gclose_f(cgl_g_id, error)

         if (otype == O_OUT) &
            & call gdf_create_root_datasets(file_id, rd)

         ! describe_domains
         call h5gcreate_f(file_id, d_gname, doml_g_id, error) ! create "/domains"

         call h5gcreate_f(doml_g_id, base_d_gname, dom_g_id, error) ! create "/domains/base"
         call create_attribute(dom_g_id, d_size_aname, dom%n_d(:)) ! create "/domains/base/n_d"
         do i = xdim, zdim
            write(d_label, '(2a)') dir_pref(i), d_edge_apname
            edge(:) = dom%edge(i, :)
            call create_attribute(dom_g_id, d_label, edge) ! create "/domains/base/[xyz]-edge_position"
            write(d_label, '(2a)') dir_pref(i), d_bnd_apname
            ! create "/domains/base/[xyz]-boundary_type"
            call create_attribute(dom_g_id, d_label, int(dom%bnd(i, :), kind=4))
         enddo

         call h5gclose_f(dom_g_id, error)

         ! create "/domains/fine_count"
         ! we have only base domain at the moment
         call create_attribute(doml_g_id, d_fc_aname, [ 0_INT4 ] )

         !> \todo add here all fine domains
         ! name "fine_00000001"
         ! attributes: n_d(:), off(:), refinement
         call h5gclose_f(doml_g_id, error)

         call h5fclose_f(file_id, error)
         call h5close_f(error)

         call rd%cleanup()
      else ! send all the necessary information to the master
         allocate(cg_rl(leaves%cnt), cg_n_b(leaves%cnt, ndims), cg_off(leaves%cnt, ndims))
#ifdef NBODY
         allocate(cg_npart(leaves%cnt), cg_pid_max(leaves%cnt))
#endif /* NBODY */
         nullify(dbuf)
         if (otype == O_OUT) allocate(dbuf(cg_le:cg_dl, leaves%cnt, ndims))
         nullify(cg_n_o)
         if (otype == O_RES) allocate(cg_n_o(leaves%cnt, ndims))
         call collect_cg_data(cg_rl, cg_n_b, cg_n_o, cg_off, cg_npart, cg_pid_max, dbuf, otype)
         call MPI_Send(cg_rl,  size(cg_rl, kind=4),  MPI_INTEGER,  FIRST, tag,         MPI_COMM_WORLD, err_mpi)
         call MPI_Send(cg_n_b, size(cg_n_b, kind=4), MPI_INTEGER,  FIRST, tag+I_ONE,   MPI_COMM_WORLD, err_mpi)
         call MPI_Send(cg_off, size(cg_off, kind=4), MPI_INTEGER8, FIRST, tag+I_TWO,   MPI_COMM_WORLD, err_mpi)
         if (otype == O_OUT) call MPI_Send(dbuf,   size(dbuf, kind=4),   MPI_REAL8,    FIRST, tag+I_THREE, MPI_COMM_WORLD, err_mpi)
         if (otype == O_RES) call MPI_Send(cg_n_o, size(cg_n_o, kind=4), MPI_INTEGER,  FIRST, tag+I_FOUR,  MPI_COMM_WORLD, err_mpi)
#ifdef NBODY
         call MPI_Send(cg_npart,   size(cg_npart,   kind=4), MPI_INTEGER8, FIRST, tag+I_FIVE, MPI_COMM_WORLD, err_mpi)
         call MPI_Send(cg_pid_max, size(cg_pid_max, kind=4), MPI_INTEGER8, FIRST, tag+I_SIX,  MPI_COMM_WORLD, err_mpi)
         deallocate(cg_npart, cg_pid_max)
#endif /* NBODY */
         deallocate(cg_rl, cg_n_b, cg_off)
         if (associated(dbuf)) deallocate(dbuf)
         if (associated(cg_n_o)) deallocate(cg_n_o)
      endif

      call piernik_MPI_Bcast(cg_all_n_b)
      call piernik_MPI_Bcast(cg_all_n_o)
      ! Reopen the HDF5 file for parallel write
      call h5open_f(error)
      if (can_i_write) then
         plist_id = set_h5_properties(H5P_FILE_ACCESS_F, nproc_io)
         call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
         call h5pclose_f(plist_id, error)
         plist_id = set_h5_properties(H5P_GROUP_ACCESS_F, nproc_io)
         call h5gopen_f(file_id, data_gname, cgl_g_id, error, gapl_id = plist_id)
         call h5pclose_f(plist_id, error)
      endif

      call write_cg_to_hdf5(cgl_g_id, cg_n, cg_all_n_b, cg_all_n_o) !!!!!

      if (can_i_write) then
         call h5gclose_f(cgl_g_id, error)
         call h5fclose_f(file_id, error)  ! Close the file
      endif

      call h5close_f(error)            ! Close HDF5 stuff

      deallocate(cg_n, cg_all_n_b, cg_all_n_o)

      return
      if (.false.) error = int(st_g_id, kind=4)

   end subroutine write_to_hdf5_v2

!>
!! \brief Collect data to fill cg_desc structure
!!
!! The loop here has to match the loop in get_nth_cg!
!>
   subroutine collect_cg_data(cg_rl, cg_n_b, cg_n_o, cg_off, cg_npart, cg_pid_max, dbuf, otype)

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: LO, HI, I_ONE
#ifdef NBODY
      use star_formation, only: pid_gen
#endif /* NBODY */

      implicit none

      integer(kind=4), dimension(:),     pointer, intent(inout) :: cg_rl            !< list of refinement levels from all cgs/procs
      integer(kind=4), dimension(:,:),   pointer, intent(inout) :: cg_n_b           !< list of n_b from all cgs/procs
      integer(kind=4), dimension(:,:),   pointer, intent(inout) :: cg_n_o           !< list of grid dimensions with external guardcells from all cgs/procs
      integer(kind=8), dimension(:,:),   pointer, intent(inout) :: cg_off           !< list of offsets from all cgs/procs with respect to level offset (lose level offset in the restart)
      integer(kind=8), dimension(:),     pointer, intent(inout) :: cg_npart         !< list of particle count from all cgs/procs
      integer(kind=8), dimension(:),     pointer, intent(inout) :: cg_pid_max       !< list of maximum particle id from all cgs/procs
      real(kind=8),    dimension(:,:,:), pointer, intent(inout) :: dbuf
      integer(kind=4),                            intent(in)    :: otype            !< Output type (restart, data)

      type(cg_list_element), pointer :: cgl
      integer :: g

      g = 1
      cgl => leaves%first
      do while (associated(cgl))
         cg_rl (g   ) = int(cgl%cg%l%id, kind=4)
         cg_n_b(g, :) = cgl%cg%n_b(:)
         cg_off(g, :) = cgl%cg%my_se(:, LO) - cgl%cg%l%off(:)
#ifdef NBODY
         cg_npart(g)   = cgl%cg%count_particles()
         cg_pid_max(g) = pid_gen
#endif /* NBODY */
         if (otype == O_OUT) then
            dbuf(cg_le, g, :)  = cgl%cg%fbnd(:, LO)
            dbuf(cg_re, g, :)  = cgl%cg%fbnd(:, HI)
            dbuf(cg_dl, g, :)  = cgl%cg%dl
         endif
         if (otype == O_RES) cg_n_o(g, :) = cgl%cg%lh_out(:, HI) - cgl%cg%lh_out(:, LO) + I_ONE
         g = g + 1
         cgl => cgl%nxt
      enddo

      return
      if (.false.) g = kind(cg_npart) + kind(cg_pid_max)

   end subroutine collect_cg_data

   function set_h5_properties(h5p, nproc_io) result (plist_id)
      use hdf5,     only: HID_T, H5P_FILE_ACCESS_F, h5pcreate_f, h5pset_fapl_mpio_f, &
         &                H5FD_MPIO_COLLECTIVE_F, h5pset_dxpl_mpio_f, H5P_DATASET_XFER_F
      use MPIF,     only: MPI_INFO_NULL, MPI_COMM_WORLD

      implicit none
      integer(HID_T),  intent(in) :: h5p
      integer(kind=4), intent(in) :: nproc_io
      integer(HID_T)              :: plist_id
      integer(kind=4)             :: error    !< error perhaps should be of type integer(HID_T)

      call h5pcreate_f(h5p, plist_id, error)
      if (nproc_io > 1) then
         ! when nproc_io < nproc we'll probably need another communicator for subset of processes that
         ! have can_i_write flag set
         if (h5p == H5P_FILE_ACCESS_F) then
#ifdef MPIF08
            call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD%mpi_val, MPI_INFO_NULL%mpi_val, error)  ! really?
#else /* !MPIF08 */
            call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
#endif /* !MPIF08 */
         else if (h5p == H5P_DATASET_XFER_F) then
            call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
         endif
      endif
   end function set_h5_properties

   function output_fname(wr_rd, ext, no, allproc, bcast, prefix) result(filename)

      use constants,  only: cwdlen, idlen, RD, WR, I_FOUR, domlen, fnamelen
      use dataio_pub, only: problem_name, run_id, res_id, wd_wr, wd_rd, warn, die, msg
      use mpisetup,   only: master, piernik_MPI_Bcast, proc

      implicit none

      integer(kind=4),       intent(in)           :: wr_rd, no
      character(len=I_FOUR), intent(in)           :: ext
      logical,               intent(in), optional :: allproc
      logical,               intent(in), optional :: bcast
      character(len=*),      intent(in), optional :: prefix
      character(len=cwdlen)                       :: filename, temp  ! File name
      character(len=domlen)                       :: fullext
      character(len=idlen)                        :: file_id
      logical                                     :: exec_allproc


      ! Sanity checks go here
      if (present(prefix)) then
         if (len_trim(prefix) >= fnamelen) then
            write(msg,*) "[common_hdf5:output_fname]:", trim(prefix), "is longer than the allowed filename"
            call die(msg)
         endif
         if (len_trim(prefix) > fnamelen/2) then
            write(msg,*) "[common_hdf5:output_fname]: There is high chance that ", trim(prefix), &
               & " will overflow the filename"
            call warn(msg)
         endif
      endif

      if ((wr_rd == RD) .and. (res_id /= '')) then
         file_id = res_id
      else
         file_id = run_id
      endif

      exec_allproc = .false.
      fullext = ext
      if (present(allproc)) then
         if (allproc) then
            exec_allproc = .true.
            write(fullext,'(".cpu",i5.5,a)') proc, ext
         endif
      endif

      if (master .or. exec_allproc) then
         if (present(prefix)) then
            write(temp,'(2(a,"_"),a3,"_",i4.4,a)') trim(prefix), trim(problem_name), file_id, no, fullext
         else
            write(temp,'(a,"_",a3,"_",i4.4,a)') trim(problem_name), file_id, no, fullext
         endif
         select case (wr_rd)
            case (RD)
               write(filename,'(2a)') trim(wd_rd),trim(temp)
            case (WR)
               write(filename,'(2a)') trim(wd_wr),trim(temp)
            case default
               write(filename,'(2a)') './',trim(temp)
         end select
      endif

      if (present(bcast)) then
         if (bcast) call piernik_MPI_Bcast(filename, cwdlen)
      endif

   end function output_fname

   subroutine dump_announcement(dumptype, nio, fname, last_dump_time, sequential)

      use constants,  only: cwdlen, tmr_hdf, I_ONE, I_FOUR, I_SEVEN, RES, WR
      use dataio_pub, only: msg, printio, thdf, multiple_h5files, piernik_hdf5_version, piernik_hdf5_version2, use_v2_io
      use global,     only: t
      use mpisetup,   only: master
      use timer,      only: set_timer

      implicit none

      integer(kind=4),                     intent(in)    :: dumptype
      integer(kind=4),                     intent(inout) :: nio
      character(len=cwdlen),               intent(out)   :: fname
      real,                                intent(in)    :: last_dump_time
      logical,                             intent(in)    :: sequential
      character(len=I_SEVEN), dimension(2), parameter    :: dumpname = ['restart', 'dataset']
      character(len=I_FOUR),  dimension(2), parameter    :: extname  = ['.res', '.h5 ']
      real                                               :: phv
      logical                                            :: exec_mh5f

      thdf = set_timer(tmr_hdf,.true.)

      nio = nio + I_ONE
! restart
      exec_mh5f = (multiple_h5files .and. dumptype /= RES)
      fname = output_fname(WR, extname(dumptype), nio, allproc=exec_mh5f, bcast=(.not.exec_mh5f))

      if (.not. master) return

      phv = piernik_hdf5_version
      if (use_v2_io .and. .not. exec_mh5f) phv = piernik_hdf5_version2

      if (sequential) then
         write(msg,'(a,es23.16,a,a,a,f5.2,1x,2a)') 'ordered t ', last_dump_time, ': Writing ', dumpname(dumptype), ' v', phv, trim(fname), " ... "
      else
         write(msg,'(a,es23.16,a,a,a,f5.2,1x,2a)') 'requested at t ', t, ': Writing ', dumpname(dumptype), ' v', phv, trim(fname), " ... "
      endif
      call printio(msg, .true.)

   end subroutine dump_announcement

   subroutine dump_announce_time

      use constants,  only: tmr_hdf
      use dataio_pub, only: msg, printinfo, thdf
      use mpisetup,   only: master
      use timer,      only: set_timer

      implicit none

      thdf = set_timer(tmr_hdf)
      if (master) then
         write(msg,'(a6,f10.2,a2)') ' done ', thdf, ' s'
         call printinfo(msg, .true.)
      endif

   end subroutine dump_announce_time

#ifdef NBODY
   subroutine initialize_write_cg(this, cgl_g_id, cg_n, nproc_io, ntags, dsets, pdsets)
#else /* !NBODY */
   subroutine initialize_write_cg(this, cgl_g_id, cg_n, nproc_io, ntags, dsets)
#endif /* !NBODY */

      use constants,  only: dsetnamelen, I_ONE
      use dataio_pub, only: can_i_write
      use hdf5,       only: HID_T, H5P_GROUP_ACCESS_F, H5P_DATASET_ACCESS_F, H5P_DATASET_XFER_F, &
          &                 h5gopen_f, h5pclose_f, h5dopen_f
      use mpisetup,   only: FIRST, LAST

      implicit none

      class(cg_output),                         intent(inout) :: this
      integer(HID_T),                           intent(in)    :: cgl_g_id
      integer(kind=4),   pointer, dimension(:), intent(in)    :: cg_n
      integer(kind=4),                          intent(in)    :: nproc_io
      integer(kind=4),                          intent(out)   :: ntags
      character(len=dsetnamelen), dimension(:), intent(in)    :: dsets
#ifdef NBODY
      character(len=dsetnamelen), dimension(:), intent(in)    :: pdsets
#endif /* NBODY */

      integer(kind=4)                                         :: i, ncg
      integer(HID_T)                                          :: plist_id
      integer(kind=4)                                         :: error    !< error perhaps should be of type integer(HID_T)

      ntags = ubound(dsets, 1, kind=4) + I_ONE
      this%tot_cg_n = sum(cg_n)
      allocate(this%cg_src_p(1:this%tot_cg_n))
      allocate(this%cg_src_n(1:this%tot_cg_n))
      allocate(this%cg_g_id(1:this%tot_cg_n))
#ifdef NBODY
      allocate(this%part_g_id(1:this%tot_cg_n))
      allocate(this%st_g_id(1:this%tot_cg_n))
      ntags = ntags + ubound(pdsets, 1, kind=4) + I_ONE
#endif /* NBODY */
      allocate(this%offsets(0:nproc_io-1))

      ! construct source addresses of the cg to be written
      do i = FIRST, LAST
         this%cg_src_p(sum(cg_n(:i))-cg_n(i)+1:sum(cg_n(:i))) = i
         do ncg = 1, cg_n(i)
            this%cg_src_n(sum(cg_n(:i))-cg_n(i)+ncg) = ncg
         enddo
      enddo

      !> \todo silent assumption that nproc_io == nproc FIXME
      this%offsets(:) = 0
      if (nproc_io > 0) then
         do i = 1, nproc_io - I_ONE
            this%offsets(i) = sum(cg_n(:i-1))
         enddo
      endif

      !> \todo Do a consistency check
      if (can_i_write) then

         plist_id = set_h5_properties(H5P_GROUP_ACCESS_F, nproc_io)
         do ncg = I_ONE, this%tot_cg_n
            call h5gopen_f(cgl_g_id, n_cg_name(ncg), this%cg_g_id(ncg), error, gapl_id = plist_id)
#ifdef NBODY
            call h5gopen_f(this%cg_g_id(ncg), part_gname, this%part_g_id(ncg), error, gapl_id = plist_id)
            call h5gopen_f(this%part_g_id(ncg), st_gname, this%st_g_id(ncg),   error, gapl_id = plist_id)
#endif /* NBODY */
         enddo
         call h5pclose_f(plist_id, error)

         plist_id = set_h5_properties(H5P_DATASET_ACCESS_F, nproc_io)
         allocate(this%dset_id(1:this%tot_cg_n, lbound(dsets, dim=1):ubound(dsets, dim=1)))
#ifdef NBODY
         allocate(this%pdset_id(1:this%tot_cg_n, lbound(pdsets, dim=1):ubound(pdsets, dim=1)))
#endif /* NBODY */
         do ncg = I_ONE, this%tot_cg_n
            do i = lbound(dsets, dim=1, kind=4), ubound(dsets, dim=1, kind=4)
               call h5dopen_f(this%cg_g_id(ncg), dsets(i), this%dset_id(ncg,i), error, dapl_id = plist_id)
            enddo

#ifdef NBODY
            do i = lbound(pdsets, dim=1, kind=4), ubound(pdsets, dim=1, kind=4)
               call h5dopen_f(this%st_g_id(ncg), pdsets(i), this%pdset_id(ncg,i), error, dapl_id = plist_id)
            enddo
#endif /* NBODY */
         enddo
         call h5pclose_f(plist_id, error)
      endif

      this%xfer_prp = set_h5_properties(H5P_DATASET_XFER_F, nproc_io)

   end subroutine initialize_write_cg

   subroutine finalize_write_cg(this)

      use hdf5,       only: h5dclose_f, h5gclose_f, h5pclose_f
      use dataio_pub, only: can_i_write

      implicit none

      class(cg_output), intent(inout) :: this

      integer                         :: ncg, i
      integer(kind=4)                 :: error    !< error perhaps should be of type integer(HID_T)

      if (can_i_write) then
         do ncg = lbound(this%cg_g_id, 1), ubound(this%cg_g_id, 1)
            do i = lbound(this%dset_id, 2), ubound(this%dset_id, 2)
               call h5dclose_f(this%dset_id(ncg, i), error)
            enddo
#ifdef NBODY
            do i = lbound(this%pdset_id, 2), ubound(this%pdset_id, 2)
               call h5dclose_f(this%pdset_id(ncg, i), error)
            enddo
            call h5gclose_f(this%st_g_id(ncg), error)
            call h5gclose_f(this%part_g_id(ncg), error)
#endif /* NBODY */
            call h5gclose_f(this%cg_g_id(ncg), error)
         enddo
      endif
      call h5pclose_f(this%xfer_prp, error)
      if (allocated(this%dset_id))   deallocate(this%dset_id)
#ifdef NBODY
      if (allocated(this%pdset_id))  deallocate(this%pdset_id)
      if (allocated(this%st_g_id))   deallocate(this%st_g_id)
      if (allocated(this%part_g_id)) deallocate(this%part_g_id)
#endif /* NBODY */
      if (allocated(this%cg_g_id))   deallocate(this%cg_g_id)
      if (allocated(this%cg_src_p))  deallocate(this%cg_src_p)
      if (allocated(this%cg_src_n))  deallocate(this%cg_src_n)

   end subroutine finalize_write_cg

end module common_hdf5
