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
#include "macros.h"

!!$ ============================================================================
!>
!! \brief Multigrid Poisson solver
!<

module multigrid

#ifdef MULTIGRID

   implicit none

   private
   public :: init_multigrid, cleanup_multigrid

contains

!!$ ============================================================================
!!
!! Initializations and cleanup
!!

   subroutine init_multigrid(cgrid)

      use multigridvars,      only: lvl, level_max, level_min, level_gb, roof, base, gb, gb_cartmap, mg_nb, ngridvars, correction, &
           &                        is_external, eff_dim, NDIM, has_dir, XDIR, YDIR, ZDIR, XLO, XHI, YLO, YHI, ZLO, ZHI, LOW, HIGH, D_x, D_y, D_z, &
           &                        ord_prolong, ord_prolong_face, stdout, verbose_vcycle, hdf5levels, tot_ts
      use types,              only: grid_container
      use mpisetup,           only: buffer_dim, comm, comm3d, ierr, proc, nproc, ndims, pxsize, pysize, pzsize, &
           &                        ibuff, rbuff, lbuff, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_LOGICAL
      use multigridhelpers,   only: mg_write_log, dirtyH, do_ascii_dump, dirty_debug, multidim_code_3D, &
           &                        aux_par_I0, aux_par_I1, aux_par_I2, aux_par_R0, aux_par_R1, aux_par_R2
      use multigridmpifuncs,  only: mpi_multigrid_prep
      use dataio_pub,         only: msg, par_file, die, warn, namelist_errh, compare_namelist
#ifdef GRAV
      use multigrid_gravity,  only: init_multigrid_grav, init_multigrid_grav_post
#endif /* GRAV */
#ifdef COSM_RAYS
      use multigrid_diffusion, only: init_multigrid_diff, init_multigrid_diff_post
#endif /* COSM_RAYS */
#ifdef NEW_HDF5
      use multigridio,        only: multigrid_add_hdf5
#endif /* NEW_HDF5 */

      implicit none

      type(grid_container), intent(in) :: cgrid                  !< copy of grid variables

      integer                          :: ierrh, div, idx, i, j, nxc, nx
      logical, save                    :: frun = .true.          !< First run flag
      real                             :: mb_alloc               !< Allocation counter
      integer, dimension(6)            :: aerr                   !BEWARE: hardcoded magic integer. Update when you change number of simultaneous error checks

      namelist /MULTIGRID_SOLVER/ level_max, ord_prolong, ord_prolong_face, &
           &                      stdout, verbose_vcycle, do_ascii_dump, dirty_debug, hdf5levels, multidim_code_3D, &
           &                      aux_par_I0, aux_par_I1, aux_par_I2, aux_par_R0, aux_par_R1, aux_par_R2

      if (.not.frun) call die("[multigrid:init_multigrid] Called more than once.")
      frun = .false.

      if (ndims /= NDIM) call die("[multigrid:init_multigrid] broken dimensional constants")

      ! Default values for namelist variables
      level_max         = 1
      ord_prolong       = 0
      ord_prolong_face  = 0
      ! May all the logical parameters be .false. by default
      stdout                 = .false.
      verbose_vcycle         = .false.
      do_ascii_dump          = .false.
      dirty_debug            = .false.
      hdf5levels             = .false.
      multidim_code_3D       = .false.

      aux_par_I0 = 0 ; aux_par_I1 = 0 ; aux_par_I2 = 0
      aux_par_R0 = 0.; aux_par_R1 = 0.; aux_par_R2 = 0.

      if (proc == 0) then

         diff_nml(MULTIGRID_SOLVER)

         ibuff(1) = level_max
         ibuff(2) = ord_prolong
         ibuff(3) = ord_prolong_face

         lbuff(1) = stdout
         lbuff(2) = verbose_vcycle
         lbuff(3) = do_ascii_dump
         lbuff(4) = dirty_debug
         lbuff(5) = hdf5levels
         lbuff(6) = multidim_code_3D

         rbuff(buffer_dim  ) = aux_par_R0
         rbuff(buffer_dim-1) = aux_par_R1
         rbuff(buffer_dim-2) = aux_par_R2

         ibuff(buffer_dim  ) = aux_par_I0
         ibuff(buffer_dim-1) = aux_par_I1
         ibuff(buffer_dim-2) = aux_par_I2

      endif

      call MPI_Bcast(ibuff, buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_Bcast(rbuff, buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      call MPI_Bcast(lbuff, buffer_dim, MPI_LOGICAL,          0, comm, ierr)

      if (proc /= 0) then

         level_max        = ibuff(1)
         ord_prolong      = ibuff(2)
         ord_prolong_face = ibuff(3)

         stdout           = lbuff(1)
         verbose_vcycle   = lbuff(2)
         do_ascii_dump    = lbuff(3)
         dirty_debug      = lbuff(4)
         hdf5levels       = lbuff(5)
         multidim_code_3D = lbuff(6)

         aux_par_R0       = rbuff(buffer_dim)
         aux_par_R1       = rbuff(buffer_dim-1)
         aux_par_R2       = rbuff(buffer_dim-2)

         aux_par_I0       = ibuff(buffer_dim)
         aux_par_I1       = ibuff(buffer_dim-1)
         aux_par_I2       = ibuff(buffer_dim-2)

      endif

      ngridvars = correction  !< 4 variables are required for basic use of the multigrid solver

! ToDo: Make array of subroutine pointers
#ifdef GRAV
      call init_multigrid_grav
#endif /* GRAV */
#ifdef COSM_RAYS
      call init_multigrid_diff
#endif /* COSM_RAYS */

      !! Sanity checks
      if (abs(ord_prolong) > 2*mg_nb) call die("[multigrid:init_multigrid] not enough guardcells for given prolongation operator order")
      if (allocated(lvl)) call die("[multigrid:init_multigrid] lvl already allocated")
      allocate(lvl(level_gb:level_max), stat=aerr(1))                                 ! level_gb = level_min-1 contains some global base level data
      if (aerr(1) /= 0) call die("[multigrid:init_multigrid] Allocation error: lvl")
      mb_alloc = size(lvl)

      has_dir(XDIR) = (cgrid%nxb > 1)
      has_dir(YDIR) = (cgrid%nyb > 1)
      has_dir(ZDIR) = (cgrid%nzb > 1)
      eff_dim = count(has_dir(:))
      if (eff_dim < 1 .or. eff_dim > 3) call die("[multigrid:init_multigrid] Unsupported number of dimensions.")
      if (has_dir(XDIR)) D_x = 1
      if (has_dir(YDIR)) D_y = 1
      if (has_dir(ZDIR)) D_z = 1

      !! Initialization of all regular levels (all but global base)
      !! Following loop gives us:
      !!    * SHAPE (lvl(level_max  )) = (nxb  , nyb  , nzd  ) + (2*nb, 2*nb, 2*nb)
      !!    * SHAPE (lvl(level_max-1)) = (nxb/2, nyb/2, nzd/2) + (2*nb, 2*nb, 2*nb)
      !!    * SHAPE (lvl(level_max-2)) = (nxb/4, nyb/4, nzd/4) + (2*nb, 2*nb, 2*nb)
      !!    * ...
      !!    * SHAPE (lvl(1)) = (nxb/2**(level_max-1)0, nyb/2**(level_max-1), nzd/2**(level_max-1)) + (2*nb, 2*nb, 2*nb)
      do idx = level_max, level_min, -1

         lvl(idx)%level = idx                                      ! level number

         div = 2**(level_max -idx)                                 ! derefinement factor with respect to the top level
         lvl(idx)%nb    = mg_nb                                    ! number of guardcells

         nxc = 1  ! suppres warning on possibly uninitialized variables
         do i = XDIR, ZDIR ! this can be rewritten as a three subroutine/function calls
            select case (i)
               case (XDIR)
                  nxc = cgrid%nxb
               case (YDIR)
                  nxc = cgrid%nyb
               case (ZDIR)
                  nxc = cgrid%nzb
            end select

            nx = 1
            if (has_dir(i)) then
               nx = nxc / div ! number of interior cells in direction i
               if (nx < lvl(idx)%nb) then
                  write(msg, '(2(a,i1),a,i4,2(a,i2))')"[multigrid:init_multigrid] Number of guardcells exceeds number of interior cells in the ",i," direction, ", &
                       lvl(idx)%nb, " > ", nx, " at level ", idx, ". You may try to set level_max <=", level_max-idx
                  call die(msg)
               endif
               if (nx * div /= nxc) then
                  write(msg, '(a,i1,a,3f6.1,2(a,i2))')"[multigrid:init_multigrid] Fractional number of cells in ",i," direction ", &
                       nxc/real(div), " at level ", idx, ". You may try to set level_max <=", level_max-idx
                  call die(msg)
               endif
            endif

            select case (i)
               case (XDIR)
                  lvl(idx)%nxb = nx
               case (YDIR)
                  lvl(idx)%nyb = nx
               case (ZDIR)
                  lvl(idx)%nzb = nx
            end select
         enddo

         lvl(idx)%dvol = 1.
         lvl(idx)%vol  = 1.

         ! /todo: check if these are correctly defined for multipole solver
         lvl(idx)%dxy = 1.
         lvl(idx)%dxz = 1.
         lvl(idx)%dyz = 1.

         if (has_dir(XDIR)) then
            lvl(idx)%nx    = lvl(idx)%nxb + 2*lvl(idx)%nb             ! total number of cells in x, y and z directions
            lvl(idx)%dx    = (cgrid%xmaxb-cgrid%xminb) / lvl(idx)%nxb ! cell size in x, y and z directions
            lvl(idx)%is    = lvl(idx)%nb + 1                          ! lowest and highest indices for interior cells
            lvl(idx)%ie    = lvl(idx)%nb + lvl(idx)%nxb
            lvl(idx)%idx2  = 1. / lvl(idx)%dx**2                      ! auxiliary invariants
            lvl(idx)%dvol  = lvl(idx)%dvol * lvl(idx)%dx              ! cell volume
            lvl(idx)%vol   = lvl(idx)%vol * (cgrid%xmaxb-cgrid%xminb)
            lvl(idx)%dxy   = lvl(idx)%dxy * lvl(idx)%dx
            lvl(idx)%dxz   = lvl(idx)%dxz * lvl(idx)%dx
         else
            lvl(idx)%nx    = 1
            lvl(idx)%dx    = huge(1.0)
            lvl(idx)%is    = 1
            lvl(idx)%ie    = 1
            lvl(idx)%idx2  = 0.
         endif

         if (has_dir(YDIR)) then
            lvl(idx)%ny    = lvl(idx)%nyb + 2*lvl(idx)%nb
            lvl(idx)%dy    = (cgrid%ymaxb-cgrid%yminb) / lvl(idx)%nyb
            lvl(idx)%js    = lvl(idx)%nb + 1
            lvl(idx)%je    = lvl(idx)%nb + lvl(idx)%nyb
            lvl(idx)%idy2  = 1. / lvl(idx)%dy**2
            lvl(idx)%dvol  = lvl(idx)%dvol * lvl(idx)%dy
            lvl(idx)%vol   = lvl(idx)%vol * (cgrid%ymaxb-cgrid%yminb)
            lvl(idx)%dxy   = lvl(idx)%dxy * lvl(idx)%dy
            lvl(idx)%dyz   = lvl(idx)%dyz * lvl(idx)%dy
         else
            lvl(idx)%ny    = 1
            lvl(idx)%dy    = huge(1.0)
            lvl(idx)%js    = 1
            lvl(idx)%je    = 1
            lvl(idx)%idy2  = 0.
         endif

         if (has_dir(ZDIR)) then
            lvl(idx)%nz    = lvl(idx)%nzb + 2*lvl(idx)%nb
            lvl(idx)%dz    = (cgrid%zmaxb-cgrid%zminb) / lvl(idx)%nzb
            lvl(idx)%ks    = lvl(idx)%nb + 1
            lvl(idx)%ke    = lvl(idx)%nb + lvl(idx)%nzb
            lvl(idx)%idz2  = 1. / lvl(idx)%dz**2
            lvl(idx)%dvol  = lvl(idx)%dvol * lvl(idx)%dz
            lvl(idx)%vol   = lvl(idx)%vol * (cgrid%zmaxb-cgrid%zminb)
            lvl(idx)%dxz   = lvl(idx)%dxz * lvl(idx)%dz
            lvl(idx)%dyz   = lvl(idx)%dyz * lvl(idx)%dz
         else
            lvl(idx)%nz    = 1
            lvl(idx)%dz    = huge(1.0)
            lvl(idx)%ks    = 1
            lvl(idx)%ke    = 1
            lvl(idx)%idz2  = 0.
         endif

         lvl(idx)%dvol2 = lvl(idx)%dvol**2

         ! data storage
         ! BEWARE prolong_x and %prolong_xy are used only with RBGS relaxation when ord_prolong /= 0
         if ( allocated(lvl(idx)%prolong_x) .or. allocated(lvl(idx)%prolong_xy) .or. allocated(lvl(idx)%mgvar) .or. &
              allocated(lvl(idx)%x) .or. allocated(lvl(idx)%y) .or. allocated(lvl(idx)%z) ) call die("[multigrid:init_multigrid] multigrid arrays already allocated")
         allocate( lvl(idx)%mgvar     (lvl(idx)%nx, lvl(idx)%ny,                  lvl(idx)%nz,                  ngridvars), stat=aerr(1) )
         allocate( lvl(idx)%prolong_x (lvl(idx)%nx, lvl(idx)%nyb/2+2*lvl(idx)%nb, lvl(idx)%nzb/2+2*lvl(idx)%nb),            stat=aerr(2) )
         allocate( lvl(idx)%prolong_xy(lvl(idx)%nx, lvl(idx)%ny,                  lvl(idx)%nzb/2+2*lvl(idx)%nb),            stat=aerr(3) )
         allocate( lvl(idx)%x         (lvl(idx)%nx),                                                                        stat=aerr(4) )
         allocate( lvl(idx)%y         (lvl(idx)%ny),                                                                        stat=aerr(5) )
         allocate( lvl(idx)%z         (lvl(idx)%nz),                                                                        stat=aerr(6) )
         if (any(aerr(1:6) /= 0)) call die("[multigrid:init_multigrid] Allocation error: lvl(idx)%*")
         if ( .not. allocated(lvl(idx)%prolong_x) .or. .not. allocated(lvl(idx)%prolong_xy) .or. .not. allocated(lvl(idx)%mgvar) .or. &
              .not. allocated(lvl(idx)%x) .or. .not. allocated(lvl(idx)%y) .or. .not. allocated(lvl(idx)%z) ) &
              call die("[multigrid:init_multigrid] some multigrid arrays not allocated")
         mb_alloc  = mb_alloc + size(lvl(idx)%prolong_x) + size(lvl(idx)%prolong_xy) + size(lvl(idx)%mgvar) + size(lvl(idx)%x)  + size(lvl(idx)%y) + size(lvl(idx)%z)

         if ( allocated(lvl(idx)%bnd_x) .or. allocated(lvl(idx)%bnd_y) .or. allocated(lvl(idx)%bnd_z)) call die("[multigrid:init_multigrid] multigrid boundary arrays already allocated")
         allocate( lvl(idx)%bnd_x(lvl(idx)%js:lvl(idx)%je, lvl(idx)%ks:lvl(idx)%ke, LOW:HIGH), stat=aerr(1) )
         allocate( lvl(idx)%bnd_y(lvl(idx)%is:lvl(idx)%ie, lvl(idx)%ks:lvl(idx)%ke, LOW:HIGH), stat=aerr(2) )
         allocate( lvl(idx)%bnd_z(lvl(idx)%is:lvl(idx)%ie, lvl(idx)%js:lvl(idx)%je, LOW:HIGH), stat=aerr(3) )
         if (any(aerr(1:3) /= 0)) call die("[multigrid:init_multigrid] Allocation error: lvl(idx)%bnd_?")
         mb_alloc  = mb_alloc + size(lvl(idx)%bnd_x) + size(lvl(idx)%bnd_y) + size(lvl(idx)%bnd_z)

         ! array initialization
         if (dirty_debug) then
            lvl(idx)%mgvar     (:, :, :, :) = dirtyH
            lvl(idx)%prolong_x (:, :, :)    = dirtyH
            lvl(idx)%prolong_xy(:, :, :)    = dirtyH
            lvl(idx)%bnd_x     (:, :, :)    = dirtyH
            lvl(idx)%bnd_y     (:, :, :)    = dirtyH
            lvl(idx)%bnd_z     (:, :, :)    = dirtyH
         else
            lvl(idx)%mgvar     (:, :, :, :) = 0.0 ! should not be necessary if dirty_debug shows nothing suspicious
         endif

         if (has_dir(XDIR)) then
            do j = 1, lvl(idx)%nx
               lvl(idx)%x(j)  = cgrid%xminb + 0.5*lvl(idx)%dx + (j-lvl(idx)%nb-1)*lvl(idx)%dx
            enddo
         else
            lvl(idx)%x(:) = (cgrid%xminb + cgrid%xmaxb) / 2.
         endif

         if (has_dir(YDIR)) then
            do j = 1, lvl(idx)%ny
               lvl(idx)%y(j)  = cgrid%yminb + 0.5*lvl(idx)%dy + (j-lvl(idx)%nb-1)*lvl(idx)%dy
            enddo
         else
            lvl(idx)%y(:) = (cgrid%yminb + cgrid%ymaxb) / 2.
         endif

         if (has_dir(ZDIR)) then
            do j = 1, lvl(idx)%nz
               lvl(idx)%z(j)  = cgrid%zminb + 0.5*lvl(idx)%dz + (j-lvl(idx)%nb-1)*lvl(idx)%dz
            enddo
         else
            lvl(idx)%z(:) = (cgrid%zminb + cgrid%zmaxb) / 2.
         endif

      enddo

      ! handy shortcuts
      base => lvl(level_min)
      roof => lvl(level_max)
      gb   => lvl(level_gb)

      call mpi_multigrid_prep

#ifdef NEW_HDF5
      call multigrid_add_hdf5
#endif /* NEW_HDF5 */

      tot_ts = 0.

      ! construct global PE mapping
      if (allocated(gb_cartmap)) call die("[multigrid:init_multigrid] gb_cartmap array already allocated")
      allocate(gb_cartmap(0:nproc-1), stat=aerr(1))
      if (aerr(1) /= 0) call die("[multigrid:init_multigrid] Allocation error: gb_cartmap")
      mb_alloc = mb_alloc + size(gb_cartmap) !may be inaccurate
      do j=0, nproc-1
         call MPI_Cart_coords(comm3d, j, NDIM, gb_cartmap(j)%proc, ierr)
         gb_cartmap(j)%lo(XDIR) = gb_cartmap(j)%proc(XDIR) * base%nxb + 1 ! starting x, y and z indices of interior cells from
         gb_cartmap(j)%lo(YDIR) = gb_cartmap(j)%proc(YDIR) * base%nyb + 1 ! coarsest level on the gb_src array
         gb_cartmap(j)%lo(ZDIR) = gb_cartmap(j)%proc(ZDIR) * base%nzb + 1
         gb_cartmap(j)%up(XDIR) = gb_cartmap(j)%lo(XDIR)   + base%nxb - 1 ! ending indices
         gb_cartmap(j)%up(YDIR) = gb_cartmap(j)%lo(YDIR)   + base%nyb - 1
         gb_cartmap(j)%up(ZDIR) = gb_cartmap(j)%lo(ZDIR)   + base%nzb - 1
      enddo

      ! mark external faces
      is_external(:) = .false.

      ! BEWARE: ignores periodicity of the grid
      if (gb_cartmap(proc)%proc(XDIR) == 0)        is_external(XLO) = .true.
      if (gb_cartmap(proc)%proc(XDIR) == pxsize-1) is_external(XHI) = .true.
      if (gb_cartmap(proc)%proc(YDIR) == 0)        is_external(YLO) = .true.
      if (gb_cartmap(proc)%proc(YDIR) == pysize-1) is_external(YHI) = .true.
      if (gb_cartmap(proc)%proc(ZDIR) == 0)        is_external(ZLO) = .true.
      if (gb_cartmap(proc)%proc(ZDIR) == pzsize-1) is_external(ZHI) = .true.

      if (.not. has_dir(XDIR)) is_external(XLO:XHI) = .false.
      if (.not. has_dir(YDIR)) is_external(YLO:YHI) = .false.
      if (.not. has_dir(ZDIR)) is_external(ZLO:ZHI) = .false.

#ifdef GRAV
      call init_multigrid_grav_post(cgrid, mb_alloc)
#endif /* !GRAV */
#ifdef COSM_RAYS
      call init_multigrid_diff_post(cgrid, mb_alloc)
#endif /* COSM_RAYS */

      ! summary
      if (proc == 0) then
         write(msg, '(a,i2,a,3(i4,a),f6.1,a)')"[multigrid:init_multigrid] Initialized ", level_max, " levels, coarsest resolution [ ", &
            lvl(1)%nxb, ",", lvl(1)%nyb, ",", lvl(1)%nzb, " ] per processor, allocated", mb_alloc*8./1048576., "MiB" ! sizeof(double)/2.**20
         call mg_write_log(msg)
      endif

   end subroutine init_multigrid

!!$ ============================================================================
!!
!! Deallocate, destroy, demolish ...
!!

   subroutine cleanup_multigrid

      use multigridvars,      only: lvl, level_gb, level_min, level_max, has_dir, XDIR, YDIR, ZDIR, tot_ts, gb_cartmap
      use mpisetup,           only: proc, nproc, MPI_DOUBLE_PRECISION, comm3d, ierr
      use multigridhelpers,   only: mg_write_log
      use dataio_pub,         only: msg
#ifdef GRAV
      use multigrid_gravity,  only: cleanup_multigrid_grav
#endif /* GRAV */
#ifdef COSM_RAYS
      use multigrid_diffusion, only: cleanup_multigrid_diff
#endif /* COSM_RAYS */

      implicit none

      integer :: i, ib
      real, allocatable, dimension(:) :: all_ts

#ifdef GRAV
      call cleanup_multigrid_grav
#endif /* GRAV */
#ifdef COSM_RAYS
      call cleanup_multigrid_diff
#endif /* COSM_RAYS */

      if (allocated(lvl)) then
         do i=level_gb, level_max
            if (allocated(lvl(i)%prolong_xy)) deallocate(lvl(i)%prolong_xy)
            if (allocated(lvl(i)%prolong_x))  deallocate(lvl(i)%prolong_x)
            if (allocated(lvl(i)%mgvar))      deallocate(lvl(i)%mgvar)
            if (allocated(lvl(i)%x))          deallocate(lvl(i)%x)
            if (allocated(lvl(i)%y))          deallocate(lvl(i)%y)
            if (allocated(lvl(i)%z))          deallocate(lvl(i)%z)
            if (allocated(lvl(i)%bnd_x))      deallocate(lvl(i)%bnd_x)
            if (allocated(lvl(i)%bnd_y))      deallocate(lvl(i)%bnd_y)
            if (allocated(lvl(i)%bnd_z))      deallocate(lvl(i)%bnd_z)

            if (i >= level_min) then
               do ib = 1, lvl(i)%nb
                  if (has_dir(XDIR)) then
                     call MPI_Type_free(lvl(i)%MPI_YZ_LEFT_BND(ib), ierr)
                     call MPI_Type_free(lvl(i)%MPI_YZ_LEFT_DOM(ib), ierr)
                     call MPI_Type_free(lvl(i)%MPI_YZ_RIGHT_DOM(ib), ierr)
                     call MPI_Type_free(lvl(i)%MPI_YZ_RIGHT_BND(ib), ierr)
                  endif

                  if (has_dir(YDIR)) then
                     call MPI_Type_free(lvl(i)%MPI_XZ_LEFT_BND(ib), ierr)
                     call MPI_Type_free(lvl(i)%MPI_XZ_LEFT_DOM(ib), ierr)
                     call MPI_Type_free(lvl(i)%MPI_XZ_RIGHT_DOM(ib), ierr)
                     call MPI_Type_free(lvl(i)%MPI_XZ_RIGHT_BND(ib), ierr)
                  endif

                  if (has_dir(ZDIR)) then
                     call MPI_Type_free(lvl(i)%MPI_XY_LEFT_BND(ib), ierr)
                     call MPI_Type_free(lvl(i)%MPI_XY_LEFT_DOM(ib), ierr)
                     call MPI_Type_free(lvl(i)%MPI_XY_RIGHT_DOM(ib), ierr)
                     call MPI_Type_free(lvl(i)%MPI_XY_RIGHT_BND(ib), ierr)
                  endif

               enddo
            endif

         enddo
         deallocate(lvl)
      endif

      if (allocated(gb_cartmap))  deallocate(gb_cartmap)

      if (allocated(all_ts)) deallocate(all_ts)
      allocate(all_ts(0:nproc-1))

      call MPI_Gather(tot_ts, 1, MPI_DOUBLE_PRECISION, all_ts, 1, MPI_DOUBLE_PRECISION, 0, comm3d, ierr)

      if (proc == 0) then
         write(msg, '(a,3(g11.4,a))')"[multigrid] Spent ", sum(all_ts)/nproc, " seconds in multigrid_solve_* (min= ",minval(all_ts)," max= ",maxval(all_ts),")."
         call mg_write_log(msg, .false.)
      endif

      if (allocated(all_ts)) deallocate(all_ts)

   end subroutine cleanup_multigrid

#else /* MULTIGRID */
#warning This should not happen. Probably the multigrid.F90 file is included in object directory by mistake.
#endif /* MULTIGRID */

end module multigrid

!!$ ============================================================================
