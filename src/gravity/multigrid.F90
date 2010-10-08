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
#include "macros.h"

!!$ ============================================================================
!>
!! \brief Multigrid poisson solver
!<

module multigrid

#ifdef MULTIGRID

   use multigridvars    ! QA: ignoring only check

   implicit none

   private
   public :: init_multigrid, cleanup_multigrid, multigrid_solve

contains

!!$ ============================================================================
!!
!! Initializations and cleanup
!!

   subroutine init_multigrid(cgrid)

      use types,              only: grid_container
      use errh,               only: namelist_errh, die, warn
      use arrays,             only: sgp
      use constants,          only: pi, dpi
      use mpisetup,           only: buffer_dim, comm, comm3d, ierr, proc, nproc, ndims, cbuff_len, &
           &                        bnd_xl_dom, bnd_xr_dom, bnd_yl_dom, bnd_yr_dom, bnd_zl_dom, bnd_zr_dom, &
           &                        bnd_xl, bnd_xr, bnd_yl, bnd_yr, bnd_zl, bnd_zr, &
           &                        ibuff, cbuff, rbuff, lbuff, &
           &                        pxsize, pysize, pzsize, &
           &                        MPI_CHARACTER, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_LOGICAL
      use grid,               only: xmin, xmax, ymin, ymax, zmin, zmax
      use func,               only: compare_namelist
      use multigridhelpers,   only: mg_write_log, dirtyH, do_ascii_dump, dirty_debug, multidim_code_3D, &
           &                        aux_par_I0, aux_par_I1, aux_par_I2, aux_par_R0, aux_par_R1, aux_par_R2
#ifdef NEW_HDF5
      use multigridio,        only: multigrid_add_hdf5
#endif /* NEW_HDF5 */
      use multipole,          only: init_multipole, use_point_monopole, lmax, mmax, ord_prolong_mpole, coarsen_multipole, interp_pt2mom, interp_mom2pot
      use multigridmpifuncs,  only: mpi_multigrid_prep
      use dataio_public,      only: msg, par_file, cwd

      implicit none

      type(grid_container), intent(in) :: cgrid                  !< copy of grid variables

      integer                          :: ierrh, div, idx, i, j, nxc=1, nx
      logical, save                    :: frun = .true.          !< First run flag
      real                             :: mb_alloc               !< Allocation counter
      integer, dimension(6)            :: aerr                   !BEWARE: hardcoded magic integer. Update when you change number of simultaneous error checks
      real, allocatable, dimension(:)  :: kx, ky, kz             !< FFT kernel directional components for convolution
      type(soln_history), pointer      :: os

      namelist /MULTIGRID_SOLVER/ norm_tol, overrelax, overrelax_x, overrelax_y, overrelax_z, Jacobi_damp, vcycle_abort, L4_strength, &
           &                      level_max, coarsen_multipole, lmax, mmax, max_cycles, nsmool, nsmoob, nsmoof, &
           &                      ord_laplacian, ord_prolong, ord_prolong_face, ord_prolong_mpole, ord_time_extrap, &
           &                      use_point_monopole, trust_fft_solution, stdout, verbose_vcycle, gb_no_fft, prefer_rbgs_relaxation, &
           &                      fft_full_relax, prefer_modified_norm, gb_solve_gather, fft_patient, do_ascii_dump, dirty_debug, hdf5levels, multidim_code_3D, &
           &                      interp_pt2mom, interp_mom2pot, &
           &                      grav_bnd_str, &
           &                      aux_par_I0, aux_par_I1, aux_par_I2, aux_par_R0, aux_par_R1, aux_par_R2

      if (.not.frun) call die("[multigrid:init_multigrid] Called more than once.")
      frun = .false.

      if (ndims /= NDIM) call die("[multigrid:init_multigrid] broken dimensional constants")

      ! Default values for namelist variables
      norm_tol      = 1.e-6
      overrelax     = 1.
      overrelax_x   = 1.
      overrelax_y   = 1.
      overrelax_z   = 1.
      Jacobi_damp   = 1.
      vcycle_abort  = 2.
      L4_strength   = 1.0

      level_max         = 1
      coarsen_multipole = 1
      lmax              = 16
      mmax              = -1 ! will be automatically set to lmax unless explicitly limited in problem.par
      max_cycles        = 20
      nsmool            = 4
      nsmoob            = 100
      nsmoof            = 1
      ord_laplacian     = 2
      ord_prolong       = 0
      ord_prolong_face  = 0
      ord_prolong_mpole = -2
      ord_time_extrap   = 1

      ! May all the logical parameters be .false. by default
      use_point_monopole     = .false.
      trust_fft_solution     = .false.
      stdout                 = .false.
      verbose_vcycle         = .false.
      gb_no_fft              = .false.
      prefer_rbgs_relaxation = .false.
      fft_full_relax         = .false.
      prefer_modified_norm   = .false.
      gb_solve_gather        = .false.
      fft_patient            = .false.
      do_ascii_dump          = .false.
      dirty_debug            = .false.
      hdf5levels             = .false.
      multidim_code_3D       = .false.
      interp_pt2mom          = .false.
      interp_mom2pot         = .false.

      if (bnd_xl_dom /= 'per' .or. bnd_xr_dom /= 'per' .or. bnd_yl_dom /= 'per' .or. bnd_yr_dom /= 'per' .or. bnd_zl_dom /= 'per' .or. bnd_zr_dom /= 'per') then
         grav_bnd_str = "isolated" ! /todo make this a default, correct default problem.par where necessary
      else
         grav_bnd_str = "periodic"
      end if

      aux_par_I0 = 0 ; aux_par_I1 = 0 ; aux_par_I2 = 0
      aux_par_R0 = 0.; aux_par_R1 = 0.; aux_par_R2 = 0.

      if (proc == 0) then

         diff_nml(MULTIGRID_SOLVER)

         rbuff(1) = norm_tol
         rbuff(2) = overrelax
         rbuff(3) = overrelax_x
         rbuff(4) = overrelax_y
         rbuff(5) = overrelax_z
         rbuff(6) = Jacobi_damp
         rbuff(7) = vcycle_abort
         rbuff(8) = L4_strength

         ibuff( 1) = level_max
         ibuff( 2) = coarsen_multipole
         ibuff( 3) = lmax
         ibuff( 4) = mmax
         ibuff( 5) = max_cycles
         ibuff( 6) = nsmool
         ibuff( 7) = nsmoob
         ibuff( 8) = nsmoof
         ibuff( 9) = ord_laplacian
         ibuff(10) = ord_prolong
         ibuff(11) = ord_prolong_face
         ibuff(12) = ord_prolong_mpole
         ibuff(13) = ord_time_extrap

         lbuff( 1) = use_point_monopole
         lbuff( 2) = trust_fft_solution
         lbuff( 3) = stdout
         lbuff( 4) = verbose_vcycle
         lbuff( 5) = gb_no_fft
         lbuff( 6) = prefer_rbgs_relaxation
         lbuff( 7) = fft_full_relax
         lbuff( 8) = prefer_modified_norm
         lbuff( 9) = gb_solve_gather
         lbuff(10) = fft_patient
         lbuff(11) = do_ascii_dump
         lbuff(12) = dirty_debug
         lbuff(13) = hdf5levels
         lbuff(14) = multidim_code_3D
         lbuff(15) = interp_pt2mom
         lbuff(16) = interp_mom2pot

         cbuff(1) = grav_bnd_str

         rbuff(buffer_dim  ) = aux_par_R0
         rbuff(buffer_dim-1) = aux_par_R1
         rbuff(buffer_dim-2) = aux_par_R2

         ibuff(buffer_dim  ) = aux_par_I0
         ibuff(buffer_dim-1) = aux_par_I1
         ibuff(buffer_dim-2) = aux_par_I2

      end if

      call MPI_BCAST(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_BCAST(ibuff,           buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_BCAST(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      call MPI_BCAST(lbuff,           buffer_dim, MPI_LOGICAL,          0, comm, ierr)

      if (proc /= 0) then

         norm_tol       = rbuff(1)
         overrelax      = rbuff(2)
         overrelax_x    = rbuff(3)
         overrelax_y    = rbuff(4)
         overrelax_z    = rbuff(5)
         Jacobi_damp    = rbuff(6)
         vcycle_abort   = rbuff(7)
         L4_strength    = rbuff(8)

         level_max         = ibuff( 1)
         coarsen_multipole = ibuff( 2)
         lmax              = ibuff( 3)
         mmax              = ibuff( 4)
         max_cycles        = ibuff( 5)
         nsmool            = ibuff( 6)
         nsmoob            = ibuff( 7)
         nsmoof            = ibuff( 8)
         ord_laplacian     = ibuff( 9)
         ord_prolong       = ibuff(10)
         ord_prolong_face  = ibuff(11)
         ord_prolong_mpole = ibuff(12)
         ord_time_extrap   = ibuff(13)

         use_point_monopole      = lbuff( 1)
         trust_fft_solution      = lbuff( 2)
         stdout                  = lbuff( 3)
         verbose_vcycle          = lbuff( 4)
         gb_no_fft               = lbuff( 5)
         prefer_rbgs_relaxation  = lbuff( 6)
         fft_full_relax          = lbuff( 7)
         prefer_modified_norm    = lbuff( 8)
         gb_solve_gather         = lbuff( 9)
         fft_patient             = lbuff(10)
         do_ascii_dump           = lbuff(11)
         dirty_debug             = lbuff(12)
         hdf5levels              = lbuff(13)
         multidim_code_3D        = lbuff(14)
         interp_pt2mom           = lbuff(15)
         interp_mom2pot          = lbuff(16)

         grav_bnd_str   = cbuff(1)(1:len(grav_bnd_str))

         aux_par_R0     = rbuff(buffer_dim)
         aux_par_R1     = rbuff(buffer_dim-1)
         aux_par_R2     = rbuff(buffer_dim-2)

         aux_par_I0     = ibuff(buffer_dim)
         aux_par_I1     = ibuff(buffer_dim-1)
         aux_par_I2     = ibuff(buffer_dim-2)

      endif

      ! boundaries
      grav_bnd = bnd_invalid
      select case(grav_bnd_str)
         case("isolated", "iso")
            grav_bnd = bnd_isolated
         case("periodic", "per")
            if (bnd_xl_dom /= 'per' .or. bnd_xr_dom /= 'per' .or. bnd_yl_dom /= 'per' .or. bnd_yr_dom /= 'per' .or. bnd_zl_dom /= 'per' .or. bnd_zr_dom /= 'per') &
                 call die("[multigrid:init_multigrid] cannot use periodic boundaries for gravity on nonperiodic domain")
            grav_bnd = bnd_periodic
         case("dirichlet", "dir")
            grav_bnd = bnd_dirichlet
         case default
            call die("[multigrid:init_multigrid] Non-recognized boundary description.")
      end select

      if (.not. (grav_bnd == bnd_periodic .or. grav_bnd == bnd_dirichlet .or. grav_bnd == bnd_isolated) .and. .not. gb_no_fft) then
         gb_no_fft = .true.
         if (proc == 0) call warn("[multigrid:init_multigrid] Use of FFT not allowed by current boundary type/combination.")
      end if

      if (.not. prefer_rbgs_relaxation .and. any([ overrelax, overrelax_x, overrelax_y, overrelax_z ] /= 1.)) then
         if (proc == 0) call warn("[multigrid:init_multigrid] Overrelaxation is disabled for FFT local solver.")
         overrelax = 1.
         overrelax_x   = 1.
         overrelax_y   = 1.
         overrelax_z   = 1.
      end if

      if ((Jacobi_damp <= 0. .or. Jacobi_damp>1.) .and. proc == 0) then
         write(msg, '(a,g12.5,a)')"[multigrid:init_multigrid] Jacobi_damp = ",Jacobi_damp," is outside (0, 1] interval."
         call warn(msg)
      end if

      if (fft_patient) fftw_flags = FFTW_PATIENT

      !! Sanity checks
      if (max(abs(ord_laplacian), abs(ord_prolong)) > 2*mg_nb) call die("[multigrid:init_multigrid] not enough guardcells for given operator order")
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

         do i = XDIR, ZDIR ! this can be rewritten as a three subroutine/function calls
            select case(i)
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
               end if
               if (nx * div /= nxc) then
                  write(msg, '(a,i1,a,3f6.1,2(a,i2))')"[multigrid:init_multigrid] Fractional number of cells in ",i," direction ", &
                       nxc/real(div), " at level ", idx, ". You may try to set level_max <=", level_max-idx
                  call die(msg)
               end if
            end if

            select case(i)
               case (XDIR)
                  lvl(idx)%nxb = nx
               case (YDIR)
                  lvl(idx)%nyb = nx
               case (ZDIR)
                  lvl(idx)%nzb = nx
            end select
         end do

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
         end if

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
         end if

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
         end if

         lvl(idx)%dvol2 = lvl(idx)%dvol**2

         ! this should work correctly also when eff_dim < 3
         lvl(idx)%r  = overrelax   / 2.
         lvl(idx)%rx = lvl(idx)%dvol2 * lvl(idx)%idx2
         lvl(idx)%ry = lvl(idx)%dvol2 * lvl(idx)%idy2
         lvl(idx)%rz = lvl(idx)%dvol2 * lvl(idx)%idz2
         lvl(idx)%r  = lvl(idx)%r  / (lvl(idx)%rx + lvl(idx)%ry + lvl(idx)%rz)
         lvl(idx)%rx = overrelax_x * lvl(idx)%rx * lvl(idx)%r
         lvl(idx)%ry = overrelax_y * lvl(idx)%ry * lvl(idx)%r
         lvl(idx)%rz = overrelax_z * lvl(idx)%rz * lvl(idx)%r
         lvl(idx)%r  = lvl(idx)%r  * lvl(idx)%dvol2

         ! BEWARE: some of the above invariants may be not optimally defined - the convergence ratio drops when dx /= dy or dy /= dz or dx /= dz
         ! and overrelaxation factors are required to get any convergence (often poor)

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
         end if

         if (has_dir(XDIR)) then
            do j = 1, lvl(idx)%nx
               lvl(idx)%x(j)  = cgrid%xminb + 0.5*lvl(idx)%dx + (j-lvl(idx)%nb-1)*lvl(idx)%dx
            enddo
         else
            lvl(idx)%x(:) = (cgrid%xminb + cgrid%xmaxb) / 2.
         end if

         if (has_dir(YDIR)) then
            do j = 1, lvl(idx)%ny
               lvl(idx)%y(j)  = cgrid%yminb + 0.5*lvl(idx)%dy + (j-lvl(idx)%nb-1)*lvl(idx)%dy
            enddo
         else
            lvl(idx)%y(:) = (cgrid%yminb + cgrid%ymaxb) / 2.
         end if

         if (has_dir(ZDIR)) then
            do j = 1, lvl(idx)%nz
               lvl(idx)%z(j)  = cgrid%zminb + 0.5*lvl(idx)%dz + (j-lvl(idx)%nb-1)*lvl(idx)%dz
            enddo
         else
            lvl(idx)%z(:) = (cgrid%zminb + cgrid%zmaxb) / 2.
         end if

         if (prefer_rbgs_relaxation) then
            lvl(idx)%fft_type = fft_none
         else if (grav_bnd == bnd_periodic .and. nproc == 1) then
            lvl(idx)%fft_type = fft_rcr
         else if (grav_bnd == bnd_periodic .or. grav_bnd == bnd_dirichlet .or. grav_bnd == bnd_isolated) then
            lvl(idx)%fft_type = fft_dst
         else
            lvl(idx)%fft_type = fft_none
         end if

      end do

      ! handy shortcuts
      base => lvl(level_min)
      roof => lvl(level_max)
      gb   => lvl(level_gb)

      ! solution recycling

      ord_time_extrap = min(nold_max-1, max(-1, ord_time_extrap))
      nold = ord_time_extrap + 1
      if (nold > 0) then
         do j = 1, 2 ! inner and outer arrays
            if (associated(os)) nullify(os)
            select case (j)
               case(1)
                  os => inner
               case(2)
                  if (grav_bnd == bnd_isolated) os => outer
            end select
            if (associated(os)) then
               do i = 1, nold
                  if ( allocated(os%old(i)%soln) ) call die("[multigrid:init_multigrid] os%old(:)%soln arrays already allocated")
                  allocate( os%old(i)%soln( roof%nx, roof%ny, roof%nz), stat=aerr(i) )
                  mb_alloc = mb_alloc + size(os%old(i)%soln)
                  os%old(i)%time= -HUGE(1.0)
                  if (dirty_debug) os%old(i)%soln(:, :, :) = dirtyH
               end do
               if (any(aerr(1:nold) /= 0)) call die("[multigrid:init_multigrid] Allocation error: os%old(:)%soln")
               os%valid = .false.
               os%last  = 1
            end if
         end do
      end if

      sgp(:,:,:) = 0. !Initialize all the guardcells, even those which does not impact the solution

      call mpi_multigrid_prep

#ifdef NEW_HDF5
      call multigrid_add_hdf5
#endif /* NEW_HDF5 */

      ! data related to local and global base-level FFT solver
      if (gb_no_fft) then
         gb%fft_type = fft_none
      else
         select case (grav_bnd)
            case (bnd_periodic)
               gb%fft_type = fft_rcr
            case (bnd_dirichlet, bnd_isolated)
               gb%fft_type = fft_dst
            case default
               gb%fft_type = fft_none
               if (proc == 0) call warn("[multigrid:init_multigrid] gb_no_fft set but no suitable boundary conditions found. Reverting to RBGS relaxation.")
         end select
      end if

      !special initialization of global base-level FFT-related data
      if (gb%fft_type /= fft_none) then

         !BEWARE only small subset of gb% members is ever initialized

         gb%nxb = base%nxb * pxsize
         gb%nyb = base%nyb * pysize
         gb%nzb = base%nzb * pzsize

         gb%level = level_gb

         gb%idx2  = base%idx2
         gb%idy2  = base%idy2
         gb%idz2  = base%idz2

         if (gb_solve_gather) then
            if (allocated(gb_src_temp)) call die("[multigrid:init_multigrid] gb_src_temp already allocated")
            allocate(gb_src_temp(base%nxb, base%nyb, base%nzb, 0:nproc-1), stat=aerr(1))
            if(aerr(1) /= 0) call die("[multigrid:init_multigrid] Allocation error: gb_src_temp")
            mb_alloc = mb_alloc + size(gb_src_temp)
         endif

      end if

      ! FFT solver storage and data
      do idx = level_max, level_gb, -1

         lvl(idx)%planf = 0
         lvl(idx)%plani = 0

         if (lvl(idx)%fft_type /= fft_none) then

            select case (lvl(idx)%fft_type)
               case (fft_rcr)
                  lvl(idx)%nxc = lvl(idx)%nxb / 2 + 1
               case (fft_dst)
                  lvl(idx)%nxc = lvl(idx)%nxb
               case default
                  call die("[multigrid:init_multigrid] Unknown FFT type.")
            end select

            if (allocated(lvl(idx)%Green3D) .or. allocated(lvl(idx)%src)) call die("[multigrid:init_multigrid] Green3D or src arrays already allocated")
            allocate(lvl(idx)%Green3D(lvl(idx)%nxc, lvl(idx)%nyb, lvl(idx)%nzb), stat=aerr(1))
            allocate(lvl(idx)%src    (lvl(idx)%nxb, lvl(idx)%nyb, lvl(idx)%nzb), stat=aerr(2))
            if (any(aerr(1:2) /= 0)) call die("[multigrid:init_multigrid] Allocation error: lvl(idx)%Green3D or lvl(idx)%src.")
            mb_alloc = mb_alloc + size(lvl(idx)%Green3D) + size(lvl(idx)%src)

            if (allocated(kx)) deallocate(kx)
            if (allocated(ky)) deallocate(ky)
            if (allocated(kz)) deallocate(kz)
            allocate(kx(lvl(idx)%nxc), stat=aerr(1))
            allocate(ky(lvl(idx)%nyb), stat=aerr(2))
            allocate(kz(lvl(idx)%nzb), stat=aerr(3))
            if (any(aerr(1:3) /= 0)) call die("[multigrid:init_multigrid] Allocation error: k[xyz]")

            select case (lvl(idx)%fft_type)

              ! lvl(idx)%fft_norm is set such that the following sequence gives identity:
              ! call dfftw_execute(lvl(idx)%planf); lvl(idx)%fftr(:, :, :) = lvl(idx)%fftr(:, :, :) * lvl(idx)%fft_norm ; call dfftw_execute(lvl(idx)%plani)

               case (fft_rcr)
                  if (allocated(lvl(idx)%fft)) call die("[multigrid:init_multigrid] fft or Green3D array already allocated")
                  allocate(lvl(idx)%fft(lvl(idx)%nxc, lvl(idx)%nyb, lvl(idx)%nzb), stat=aerr(1))
                  if (aerr(1) /= 0) call die("[multigrid:init_multigrid] Allocation error: fft.")
                  mb_alloc  = mb_alloc + 2*size(lvl(idx)%fft)

                  lvl(idx)%fft_norm = 1. / real( lvl(idx)%nxb * lvl(idx)%nyb * lvl(idx)%nzb ) ! No 4 pi G factor here because the source was already multiplied by it

                  ! FFT local solver initialization for 2nd order (3-point) Laplacian
                  ! sin(k*x-d) - 2.*sin(k*x) + sin(k*x+d) = 2 * (cos(d)-1) * sin(k*x) = -4 * sin(d/2)**2 * sin(k*x)
                  ! For 4th order: a*sin(k*x) + b*(sin(k*x-d) + sin(k*x+d)) + c*(sin(k*x-2*d) + sin(k*x+2*d)), a+2*b+2*c == 0 it would be:
                  ! 4*(a+b+(a+2*b)*cos(d)) * sin(d/2)**2 * sin(k*x)
                  ! For 6th order: a*sin(k*x) + b*(sin(k*x-d) + sin(k*x+d)) + c*(sin(k*x-2*d) + sin(k*x+2*d)) + e*(sin(k*x-3*d) + sin(k*x+3*d)), a+2*b+2*c+2*e == 0 it would be:
                  ! 2*(3*a+4*b+2*c+4*(a+2*b+c)*cos(d)+2*(a+2*(b+c))*cos(2*d)) * sin(d/2)**2 * sin(k*x)
                  ! asymptotically: -d**2/2 for d<pi

                  kx(:) = lvl(idx)%idx2 * (cos(dpi/lvl(idx)%nxb*(/(j, j=0, lvl(idx)%nxc-1)/)) - 1.)
                  ky(:) = lvl(idx)%idy2 * (cos(dpi/lvl(idx)%nyb*(/(j, j=0, lvl(idx)%nyb-1)/)) - 1.)
                  kz(:) = lvl(idx)%idz2 * (cos(dpi/lvl(idx)%nzb*(/(j, j=0, lvl(idx)%nzb-1)/)) - 1.)
                  call dfftw_plan_dft_r2c_3d(lvl(idx)%planf, lvl(idx)%nxb, lvl(idx)%nyb, lvl(idx)%nzb, lvl(idx)%src, lvl(idx)%fft, fftw_flags)
                  call dfftw_plan_dft_c2r_3d(lvl(idx)%plani, lvl(idx)%nxb, lvl(idx)%nyb, lvl(idx)%nzb, lvl(idx)%fft, lvl(idx)%src, fftw_flags)

               case (fft_dst)

                  if (allocated(lvl(idx)%fftr)) call die("[multigrid:init_multigrid] fftr array already allocated")
                  allocate(lvl(idx)%fftr(lvl(idx)%nxc, lvl(idx)%nyb, lvl(idx)%nzb), stat=aerr(1))
                  if (aerr(1) /= 0) call die("[multigrid:init_multigrid] Allocation error: fftr.")
                  mb_alloc  = mb_alloc + size(lvl(idx)%fftr)

                  lvl(idx)%fft_norm = 1. / (8. * real( lvl(idx)%nxb * lvl(idx)%nyb * lvl(idx)%nzb ))
                  kx(:) = lvl(idx)%idx2 * (cos(pi/lvl(idx)%nxb*(/(j, j=1, lvl(idx)%nxc)/)) - 1.)
                  ky(:) = lvl(idx)%idy2 * (cos(pi/lvl(idx)%nyb*(/(j, j=1, lvl(idx)%nyb)/)) - 1.)
                  kz(:) = lvl(idx)%idz2 * (cos(pi/lvl(idx)%nzb*(/(j, j=1, lvl(idx)%nzb)/)) - 1.)
                  call dfftw_plan_r2r_3d(lvl(idx)%planf, lvl(idx)%nxb, lvl(idx)%nyb, lvl(idx)%nzb, lvl(idx)%src,  lvl(idx)%fftr, &
                       &                 FFTW_RODFT10, FFTW_RODFT10, FFTW_RODFT10, fftw_flags)
                  call dfftw_plan_r2r_3d(lvl(idx)%plani, lvl(idx)%nxb, lvl(idx)%nyb, lvl(idx)%nzb, lvl(idx)%fftr, lvl(idx)%src,  &
                       &                 FFTW_RODFT01, FFTW_RODFT01, FFTW_RODFT01, fftw_flags)

               case default
                  call die("[multigrid:init_multigrid] Unknown FFT type.")
            end select

            !// compute Green's function for 7-point 3D discrete laplacian
            do i = 1, lvl(idx)%nxc
               do j = 1, lvl(idx)%nyb
                  where( (kx(i) + ky(j) + kz(:)) /= 0 )
                     lvl(idx)%Green3D(i,j,:) = 0.5 * lvl(idx)%fft_norm / (kx(i) + ky(j) + kz(:))
                  elsewhere
                     lvl(idx)%Green3D(i,j,:) = 0.0
                  endwhere
               enddo
            enddo

         end if
      end do

      if (roof%fft_type == fft_none .and. trust_fft_solution) then
         if (proc == 0) call warn("[multigrid:init_multigrid] cannot trust FFT solution on the roof.")
         trust_fft_solution = .false.
      end if

      if (allocated(kx)) deallocate(kx)
      if (allocated(ky)) deallocate(ky)
      if (allocated(kz)) deallocate(kz)

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
      end do

      ! mark external faces
      is_external(:) = .false.
      if (grav_bnd /= bnd_periodic) then
         ! BEWARE I am not sure if the checks are complete. It probably will crash on shear boundaries
         if (bnd_xl(1:3) /= "mpi") is_external(XLO) = .true.
         if (bnd_xr(1:3) /= "mpi") is_external(XHI) = .true.
         if (bnd_yl(1:3) /= "mpi") is_external(YLO) = .true.
         if (bnd_yr(1:3) /= "mpi") is_external(YHI) = .true.
         if (bnd_zl(1:3) /= "mpi") is_external(ZLO) = .true.
         if (bnd_zr(1:3) /= "mpi") is_external(ZHI) = .true.

         ! force domain boundaries for gravity if requested even when the domain is periodic
         if (gb_cartmap(proc)%proc(XDIR) == 0)        is_external(XLO) = .true.
         if (gb_cartmap(proc)%proc(XDIR) == pxsize-1) is_external(XHI) = .true.
         if (gb_cartmap(proc)%proc(YDIR) == 0)        is_external(YLO) = .true.
         if (gb_cartmap(proc)%proc(YDIR) == pysize-1) is_external(YHI) = .true.
         if (gb_cartmap(proc)%proc(ZDIR) == 0)        is_external(ZLO) = .true.
         if (gb_cartmap(proc)%proc(ZDIR) == pzsize-1) is_external(ZHI) = .true.
      end if

      if (.not. has_dir(XDIR)) is_external(XLO:XHI) = .false.
      if (.not. has_dir(YDIR)) is_external(YLO:YHI) = .false.
      if (.not. has_dir(ZDIR)) is_external(ZLO:ZHI) = .false.

      if (grav_bnd == bnd_isolated) call init_multipole(mb_alloc,cgrid)

      tot_ts = 0.

      if (allocated(vcycle_factors)) deallocate(vcycle_factors)
      allocate(vcycle_factors(0:max_cycles, 2), stat=aerr(1))
      if (aerr(1) /= 0) call die("[multigrid:init_multigrid] Allocation error: vcycle_factors")
      mb_alloc = mb_alloc + size(vcycle_factors)

      cprefix=""

      ! summary
      if (proc == 0) then
         write(msg, '(a,i2,a,3(i4,a),f6.1,a)')"[multigrid:init_multigrid] Initialized ", level_max, " levels, coarsest resolution [ ", &
            lvl(1)%nxb, ",", lvl(1)%nyb, ",", lvl(1)%nzb, " ] per processor, allocated", mb_alloc*8./1048576., "MiB" ! sizeof(double)/2.**20
         call mg_write_log(msg)
         if (overrelax /= 1. .or. overrelax_x /= 1. .or. overrelax_y /= 1. .or. overrelax_z /= 1.) then
            write(msg, '(a,f8.5,a,3f8.5,a)')"[multigrid:init_multigrid] Overrelaxation factors: global = ", overrelax, ", directional = [", overrelax_x, overrelax_y, overrelax_z, "]"
            call mg_write_log(msg)
         end if
      endif

   end subroutine init_multigrid

!!$ ============================================================================
!!
!! Deallocate, destroy, demolish ...
!!

   subroutine cleanup_multigrid

      use mpisetup,           only: proc, nproc, MPI_DOUBLE_PRECISION, comm3d, ierr
      use multigridhelpers,   only: mg_write_log
      use multipole,          only: cleanup_multipole
      use dataio_public,      only: msg

      implicit none

      integer :: i, ib
      real, allocatable, dimension(:) :: all_ts

      if(allocated(lvl)) then
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
            if (allocated(lvl(i)%fft))        deallocate(lvl(i)%fft)
            if (allocated(lvl(i)%fftr))       deallocate(lvl(i)%fftr)
            if (allocated(lvl(i)%src))        deallocate(lvl(i)%src)
            if (allocated(lvl(i)%Green3D))    deallocate(lvl(i)%Green3D)

            if (lvl(i)%planf /= 0) call dfftw_destroy_plan(lvl(i)%planf)
            if (lvl(i)%plani /= 0) call dfftw_destroy_plan(lvl(i)%plani)

            if (i >= level_min) then
               do ib = 1, lvl(i)%nb
                  if (has_dir(XDIR)) then
                     call MPI_Type_free(lvl(i)%MPI_YZ_LEFT_BND(ib), ierr)
                     call MPI_Type_free(lvl(i)%MPI_YZ_LEFT_DOM(ib), ierr)
                     call MPI_Type_free(lvl(i)%MPI_YZ_RIGHT_DOM(ib), ierr)
                     call MPI_Type_free(lvl(i)%MPI_YZ_RIGHT_BND(ib), ierr)
                  end if

                  if (has_dir(YDIR)) then
                     call MPI_Type_free(lvl(i)%MPI_XZ_LEFT_BND(ib), ierr)
                     call MPI_Type_free(lvl(i)%MPI_XZ_LEFT_DOM(ib), ierr)
                     call MPI_Type_free(lvl(i)%MPI_XZ_RIGHT_DOM(ib), ierr)
                     call MPI_Type_free(lvl(i)%MPI_XZ_RIGHT_BND(ib), ierr)
                  end if

                  if (has_dir(ZDIR)) then
                     call MPI_Type_free(lvl(i)%MPI_XY_LEFT_BND(ib), ierr)
                     call MPI_Type_free(lvl(i)%MPI_XY_LEFT_DOM(ib), ierr)
                     call MPI_Type_free(lvl(i)%MPI_XY_RIGHT_DOM(ib), ierr)
                     call MPI_Type_free(lvl(i)%MPI_XY_RIGHT_BND(ib), ierr)
                  end if

               end do
            end if

         enddo
         deallocate(lvl)
      endif

      call dfftw_cleanup

      do i = 1, nold
         if (allocated(inner%old(i)%soln)) deallocate(inner%old(i)%soln)
         if (allocated(outer%old(i)%soln)) deallocate(outer%old(i)%soln)
      end do

      if (allocated(gb_cartmap))  deallocate(gb_cartmap)
      if (allocated(gb_src_temp)) deallocate(gb_src_temp)

      call cleanup_multipole

      if (allocated(all_ts)) deallocate(all_ts)
      allocate(all_ts(0:nproc-1))

      call MPI_Gather(tot_ts, 1, MPI_DOUBLE_PRECISION, all_ts, 1, MPI_DOUBLE_PRECISION, 0, comm3d, ierr)

      if (proc == 0) then
         write(msg, '(a,3(g11.4,a))')"[multigrid] Spent ", sum(all_ts)/nproc, " seconds in multigrid_solve (min= ",minval(all_ts)," max= ",maxval(all_ts),")."
         call mg_write_log(msg, .false.)
      endif

      if (allocated(vcycle_factors)) deallocate(vcycle_factors)
      if (allocated(all_ts)) deallocate(all_ts)

   end subroutine cleanup_multigrid

!!$ ============================================================================
!!
!! This routine tries to construct first guess of potential based on previously obtained solution, if any.
!! If for some reason it is not desired to do the extrapolation (i.e. timestep varies too much) it would be good to have more control on this behavior.
!!
!! Quadratic extrapolation in time often gives better guess for smooth potentials, but is more risky for sharp-peaked potential fields (like moving self-bound clumps)
!! Rational extrapolation (1 + a t)(b + c t) can give 2 times better and 10 times worse guess depending on timestep.
!!
!! Set history%valid to .false. to force start from scratch.
!!

   subroutine init_solution(history)

      use errh,             only: die
      use mpisetup,         only: proc, t
      use multigridhelpers, only: set_dirty, check_dirty, mg_write_log
      use dataio_public,    only: msg

      implicit none

      type(soln_history), intent(inout) :: history !< inner or outer potential history for recycling

      integer :: l, p0, p1, p2, ordt
      real, dimension(3)  :: dt_fac

      call set_dirty(solution)

      p0 = history%last
      if (nold > 0) then
         p1 = 1 + mod(p0 + nold - 2, nold) ! index of previous save
         p2 = 1 + mod(p1 + nold - 2, nold)
      else
         p1 = p0
         p2 = p0
      end if

      ordt = ord_time_extrap
      if (history%valid) then
         if ( history%old(p2)%time /= history%old(p1)%time .and. &       ! quadratic interpolation
              history%old(p2)%time /= history%old(p0)%time .and. &
              history%old(p1)%time /= history%old(p0)%time) then
            ordt = min(2, ord_time_extrap)
         else if (history%old(p0)%time /= history%old(p1)%time) then     ! linear extrapolation
            ordt = min(1, ord_time_extrap)
         else                                                            ! simple recycling
            ordt = min(0, ord_time_extrap)
         end if
      else                                                               ! coldstart
         ordt = min(-1, ord_time_extrap)
      end if

      select case(ordt)
         case (:-1)
            if (proc == 0 .and. ord_time_extrap > -1) then
               write(msg, '(3a)')"[multigrid:init_solution] Clearing ",trim(cprefix),"solution."
               call mg_write_log(msg, stdout)
            endif
            do l = level_min, level_max
               lvl(l)%mgvar(:, :, :, solution) = 0.
            end do
            history%old(:)%time = -HUGE(1.0)
         case (0)
            roof%mgvar(:, :, :, solution) = history%old(p0)%soln(:, :, :)
            if (proc == 0 .and. ord_time_extrap > 0) then
               write(msg, '(3a)')"[multigrid:init_solution] No extrapolation of ",trim(cprefix),"solution."
               call mg_write_log(msg, stdout)
            endif
         case (1)
            dt_fac(1) = (t - history%old(p0)%time) / (history%old(p0)%time - history%old(p1)%time)
            roof%mgvar(:, :, :, solution) = (1. + dt_fac(1)) * history%old(p0)%soln(:, :, :) - dt_fac(1) *  history%old(p1)%soln(:, :, :)
            if (proc == 0 .and. ord_time_extrap > 1) then
               write(msg, '(3a)')"[multigrid:init_solution] Linear extrapolation of ",trim(cprefix),"solution."
               call mg_write_log(msg, stdout)
            end if
         case (2)
            dt_fac(1) = (t - history%old(p0)%time) / (history%old(p1)%time - history%old(p2)%time)
            dt_fac(2) = (t - history%old(p1)%time) / (history%old(p2)%time - history%old(p0)%time)
            dt_fac(3) = (t - history%old(p2)%time) / (history%old(p0)%time - history%old(p1)%time)
            roof%mgvar(:, :, :, solution) = - dt_fac(2) * dt_fac(3) * history%old(p0)%soln(:, :, :) &
                 &                          - dt_fac(1) * dt_fac(3) * history%old(p1)%soln(:, :, :) &
                 &                          - dt_fac(1) * dt_fac(2) * history%old(p2)%soln(:, :, :)
         case default
            call die("[multigrid:init_solution] Extrapolation order not implemented")
      end select

      call check_dirty(level_max, solution, "init_soln")

   end subroutine init_solution

!!$ ============================================================================
!!
!! Solve finest level if allowed (typically on single CPU)
!!

   subroutine fft_solve_roof

      implicit none

      if (roof%fft_type == fft_none) return

      roof%src(:, :, :) = roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, source)
      call fft_convolve(roof%level)
      roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, solution) = roof%src(:, :, :)

   end subroutine fft_solve_roof

!!$ ============================================================================
!!
!! Make a local copy of source (density) and multiply by 4 pi G
!!
   subroutine init_source(dens)

#ifdef JEANS_PROBLEM
      use initproblem,        only: d0, mode ! hack for tests
#endif /* JEANS_PROBLEM */
      use constants,          only: fpiG
      use grid,               only: is, ie, js, je, ks, ke
      use errh,               only: die
      use multigridhelpers,   only: set_dirty, check_dirty
      use multigridbasefuncs, only: norm_sq, substract_average

      implicit none

      real, optional, dimension(:,:,:), intent(in)  :: dens !< input source field or nothing for empty space

      call set_dirty(source)

      if (present(dens)) then
         roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, source) = fpiG * dens(is:ie, js:je, ks:ke)
         call norm_sq(source, norm_rhs_orig) ! The norm of corrected source will be calculated in multigrid_solve
      else
         if (grav_bnd /= bnd_givenval) call die("[multigrid:init_source] empty space allowed only for given value boundaries.")
         roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, source) = 0.
         norm_rhs_orig = 0.
      end if

      select case (grav_bnd)
         case (bnd_periodic) ! probably also bnd_neumann
            call substract_average(level_max, source)
         case (bnd_dirichlet)
#ifdef JEANS_PROBLEM
            if (mode == 1) roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, source) = &
                 &         roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, source) - fpiG * d0 ! remove density bias
#endif /* JEANS_PROBLEM */
         case (bnd_givenval) ! convert potential into a layer of imaginary mass (substract second derivative normal to computational domain boundary)
            if (is_external(XLO)) roof%mgvar(roof%is,         roof%js:roof%je, roof%ks:roof%ke, source) = &
                 &                roof%mgvar(roof%is,         roof%js:roof%je, roof%ks:roof%ke, source) - &
                 &                roof%bnd_x(                 roof%js:roof%je, roof%ks:roof%ke, LOW)  * 2. * roof%idx2 / fpiG
            if (is_external(XHI)) roof%mgvar(        roof%ie, roof%js:roof%je, roof%ks:roof%ke, source) = &
                 &                roof%mgvar(        roof%ie, roof%js:roof%je, roof%ks:roof%ke, source) - &
                 &                roof%bnd_x(                 roof%js:roof%je, roof%ks:roof%ke, HIGH) * 2. * roof%idx2 / fpiG
            if (is_external(YLO)) roof%mgvar(roof%is:roof%ie, roof%js,         roof%ks:roof%ke, source) = &
                 &                roof%mgvar(roof%is:roof%ie, roof%js,         roof%ks:roof%ke, source) - &
                 &                roof%bnd_y(roof%is:roof%ie,                  roof%ks:roof%ke, LOW)  * 2. * roof%idy2 / fpiG
            if (is_external(YHI)) roof%mgvar(roof%is:roof%ie,         roof%je, roof%ks:roof%ke, source) = &
                 &                roof%mgvar(roof%is:roof%ie,         roof%je, roof%ks:roof%ke, source) - &
                 &                roof%bnd_y(roof%is:roof%ie,                  roof%ks:roof%ke, HIGH) * 2. * roof%idy2 / fpiG
            if (is_external(ZLO)) roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks,         source) = &
                 &                roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks,         source) - &
                 &                roof%bnd_z(roof%is:roof%ie, roof%js:roof%je,                  LOW)  * 2. * roof%idz2 / fpiG
            if (is_external(ZHI)) roof%mgvar(roof%is:roof%ie, roof%js:roof%je,         roof%ke, source) = &
                 &                roof%mgvar(roof%is:roof%ie, roof%js:roof%je,         roof%ke, source) - &
                 &                roof%bnd_z(roof%is:roof%ie, roof%js:roof%je,                  HIGH) * 2. * roof%idz2 / fpiG
         case default
            call die("[multigrid:init_source] Unknown boundary type")
      end select

      call check_dirty(level_max, source, "init_src")

   end subroutine init_source

!!$ ============================================================================
!!
!! This routine manages old copies of potential for recycling.
!!

   subroutine store_solution(history)

      use mpisetup,          only: proc, t
      use multigridmpifuncs, only: mpi_multigrid_bnd

      implicit none

      type(soln_history), intent(inout) :: history !< inner or outer potential history to store recent solution
      logical :: extrapolate_bnd

      if (nold <= 0) return

      if (history%valid) then
         history%last = 1 + mod(history%last, nold)
      else
         history%old(:)%time = t ! prevents extrapolation too early
      end if

      history%old(history%last)%soln(:, :, :) = roof%mgvar(:, :, :, solution)
      history%old(history%last)%time = t
      history%valid = .true.

      ! Update guardcells of the solution before leaving. This can be done in higher-level routines that collect all the gravity contributions, but would be less safe.
      ! Extrapolate isolated boundaries, remember that grav_bnd is messed up by multigrid_solve
      extrapolate_bnd = (grav_bnd == bnd_isolated .or. grav_bnd == bnd_givenval)
      call mpi_multigrid_bnd(level_max, solution, mg_nb, extrapolate_bnd)

   end subroutine store_solution

!!$ ============================================================================
!!
!! Multigrid driver. This is the only multigrid routine intended to be called from the gravity module.
!! This routine is also responsible for communicating the solution to the rest of world via sgp array.
!!

   subroutine multigrid_solve(dens)

      use timer,     only: timer_
      use arrays,    only: sgp
      use grid,      only: is, ie, js, je, ks, ke
      use errh,      only: die
      use multipole, only: multipole_solver

      implicit none

      real, dimension(:,:,:), intent(in) :: dens !< input source field

      logical :: isolated
      integer :: isb, ieb, jsb, jeb, ksb, keb

      ts =  timer_("multigrid", .true.)
      if ( (has_dir(XDIR) .and. is-mg_nb <= 0) .or. &
           (has_dir(YDIR) .and. js-mg_nb <= 0) .or. &
           (has_dir(ZDIR) .and. ks-mg_nb <= 0) )    &
           call die("[multigrid:multigrid_solve] Current implementation requires at least 2 guardcells in the hydro part")

      isolated = (grav_bnd == bnd_isolated) ! BEWARE: not elegant; probably there should be two global grav_bnd variables

      if (isolated) then
         grav_bnd = bnd_dirichlet
         cprefix = "i-"
      end if

      call init_source(dens)

      call vcycle(inner)

      ! /todo: move to multigridvars and init_multigrid
      if (has_dir(XDIR)) then
         isb = is-mg_nb
         ieb = ie+mg_nb
      else
         isb = 1
         ieb = 1
      end if

      if (has_dir(YDIR)) then
         jsb = js-mg_nb
         jeb = je+mg_nb
      else
         jsb = 1
         jeb = 1
      end if

      if (has_dir(ZDIR)) then
         ksb = ks-mg_nb
         keb = ke+mg_nb
      else
         ksb = 1
         keb = 1
      end if

      sgp(isb:ieb, jsb:jeb, ksb:keb) = roof%mgvar(:, :, :, solution)

      if (isolated) then
         grav_bnd = bnd_givenval

         cprefix = "o-"
         call multipole_solver
         call init_source

         call vcycle(outer)
         sgp(isb:ieb, jsb:jeb, ksb:keb) = sgp(isb:ieb, jsb:jeb, ksb:keb) + roof%mgvar(:, :, :, solution)

         grav_bnd = bnd_isolated ! restore
      end if

      ts = timer_("multigrid")
      tot_ts = tot_ts + ts

   end subroutine multigrid_solve

!!$ ============================================================================
!!
!! The solver. Here we choose an adaptation of the Huang-Greengard V-cycle.
!! For more difficult problems, like variable coefficient diffusion equation a more sophisticated V-cycle may be more effective.
!!

   subroutine vcycle(history)

      use errh,               only: die, warn
      use mpisetup,           only: proc, nproc, cbuff_len
      use timer,              only: timer_
      use multigridhelpers,   only: set_dirty, check_dirty, mg_write_log, brief_v_log, do_ascii_dump, numbered_ascii_dump
      use multigridbasefuncs, only: norm_sq, residual, restrict_all, substract_average
      use dataio_public,      only: msg

      implicit none

      type(soln_history), intent(inout) :: history !< inner or outer potential history used for initializing first guess

      real,    parameter :: suspicious_factor = 1.05 ! If the norm decreases too slowly then dump diagnostic output (BEWARE: this option is for tests only)
      integer            :: l, v
      real               :: norm_rhs, norm_lhs, norm_old, norm_lowest
      logical            :: dump_every_step, dump_result
      integer, parameter       :: fmtlen = 32
      character(len=fmtlen)    :: fmt
      character(len=cbuff_len) :: dname

      inquire(file = "_dump_every_step_", EXIST=dump_every_step) ! use for debug only
      inquire(file = "_dump_result_", EXIST=dump_result)

      do_ascii_dump = do_ascii_dump .or. dump_every_step .or. dump_result

      ! On single CPU use FFT if possible because it is faster. Can be disabled by prefer_rbgs_relaxation = .true.
      if (nproc == 1 .and. roof%fft_type /= fft_none) then
         call set_dirty(solution)
         call fft_solve_roof
         if (trust_fft_solution) then
            write(msg, '(3a)')"[multigrid:vcycle] FFT solution trusted, skipping ", trim(cprefix), "cycle."
            call mg_write_log(msg, stdout)
            return
         end if
      else
         call init_solution(history)
      end if

      if ((grav_bnd /= bnd_givenval) .and. .not. prefer_modified_norm) then
         norm_rhs = norm_rhs_orig
      else
         call norm_sq(source, norm_rhs)
         if (proc == 0 .and. norm_rhs<(1.-1e-6)*norm_rhs_orig) then
            write(msg, '(a,f8.5)')"[multigrid:vcycle] norm_rhs/norm_rhs_orig = ", norm_rhs/norm_rhs_orig
            call mg_write_log(msg, stdout)
         endif
      end if
      norm_old = norm_rhs
      norm_lowest = norm_rhs

      if (norm_rhs == 0.) then ! empty domain => potential == 0.
         if (proc == 0) then
            write(msg, '(a)')"[multigrid:vcycle] source == 0"
            call mg_write_log(msg)
         endif
         return
      end if

!      if (.not. history%valid .and. prefer_rbgs_relaxation) call approximate_solution(level_max, source, solution) ! not necessary when init_solution called FFT
! difficult statement: for approximate_solution_fft it requires to pass a flag to use guardcells instead of prolonging faces.
! how much does it improve? (make a benchmark at some point)

      ! iterations
      do v = 0, max_cycles

         call set_dirty(defect)
         call residual(level_max, source, solution, defect)
         if (grav_bnd == bnd_periodic) call substract_average(level_max, defect)
         call check_dirty(level_max, defect, "residual")

         call norm_sq(defect, norm_lhs)
         ts = timer_("multigrid")
         tot_ts = tot_ts + ts
         if (proc == 0 .and. verbose_vcycle) then
            if (norm_old/norm_lhs < 1.e5) then
               fmt='(3a,i3,a,f12.9,a,f8.2,a,f7.3)'
            else
               fmt='(3a,i3,a,f12.9,a,es8.2,a,f7.3)'
            end if
            write(msg, fmt)"[multigrid] ", trim(cprefix), "Cycle:", v, " norm/rhs= ", norm_lhs/norm_rhs, " reduction factor= ", norm_old/norm_lhs, "   dt_wall= ", ts
            call mg_write_log(msg, stdout)
         end if
         vcycle_factors(v,:) = [ norm_old/norm_lhs, ts ]

         if (norm_old/norm_lhs <= suspicious_factor .or. dump_every_step .or. (norm_lhs/norm_rhs <= norm_tol .and. dump_result)) then
            write(dname,'(2a)')trim(cprefix),"mdump"
            if (dump_result .and. norm_lhs/norm_rhs <= norm_tol) then
               call numbered_ascii_dump(dname)
            else
               call numbered_ascii_dump(dname, v)
            end if
         end if

         if (norm_lhs/norm_rhs <= norm_tol) exit

         if (v<1) then ! forgive poor convergence in some first V-cycles
            norm_lowest = norm_lhs
         else
            if (norm_lhs < norm_lowest) then
               norm_lowest = norm_lhs
            else
               if (norm_lhs/norm_lowest > vcycle_abort) then
                  if (.not. verbose_vcycle) call brief_v_log(v, norm_lhs/norm_rhs)
                  call die("[multigrid:vcycle] Serious nonconvergence detected.")
                  !In such case one may increase nsmool, decrease refinement depth or use FFT
               end if
            end if
         end if

         norm_old = norm_lhs

         ! the Huang-Greengard V-cycle
         call restrict_all(defect)

         call set_dirty(correction)
         base%mgvar(:, :, :, correction) = 0.

         do l = level_min, level_max
            call approximate_solution(l, defect, correction)
            call check_dirty(l, correction, "Vup relax+")
         end do
         roof%mgvar     (roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, solution) = &
              roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, solution) + &
              roof%mgvar(roof%is:roof%ie, roof%js:roof%je, roof%ks:roof%ke, correction)

      end do

      if (v == max_cycles + 1) then
         if (proc == 0 .and. norm_lhs/norm_rhs > norm_tol) call warn("[multigrid:vcycle] Not enough V-cycles to achieve convergence.")
         v = max_cycles
      end if

      if (.not. verbose_vcycle) call brief_v_log(v, norm_lhs/norm_rhs)

      call check_dirty(level_max, solution, "final_solution")

      call store_solution(history)

   end subroutine vcycle

!!$ ============================================================================
!!
!! This routine has to find an approximate solution for given source field and implemented differential operator
!!

   subroutine approximate_solution(lev, src, soln)

      use errh,               only: die
      use multigridhelpers,   only: check_dirty
      use multigridbasefuncs, only: prolong_level

      implicit none

      integer, intent(in) :: lev  !< level for which approximate the solution
      integer, intent(in) :: src  !< index of source in lvl()%mgvar
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar

      if (any( [ src, soln ] <= 0) .or. any( [ src, soln ] > ngridvars)) call die("[multigrid:approximate_solution] Invalid variable index.")
      if (lev < level_min .or. lev > level_max) call die("[multigrid:approximate_solution] Invalid level number.")

      call check_dirty(lev, src, "approx_soln src-")

      if (lev == level_min .and. .not. gb_no_fft) then
         call gb_fft_solve(src, soln)
      else
         if (prefer_rbgs_relaxation) then
            call check_dirty(lev, soln, "approx_soln soln-")
            call approximate_solution_rbgs(lev, src, soln)
         else
            call approximate_solution_fft(lev, src, soln)
         end if
      end if

      if (prefer_rbgs_relaxation .and. soln == correction .and. lev <  level_max) call prolong_level(lev, correction)
      ! BEWARE other implementations of the multigrid algorithm may be incompatible with prolongation called from here

      call check_dirty(lev, soln, "approx_soln soln+")

   end subroutine approximate_solution

!!$ ============================================================================
!!
!! FFT given-boundary Poisson solver applied to local domain. Should require less communication than RBGS implementation.
!!
!! \todo test a configuration with wider area being subjected to FFT (sizes would no longer be 2**n) to avoid the need of relaxation
!!

   subroutine approximate_solution_fft(lev, src, soln)

      use errh,               only: die, warn
      use mpisetup,           only: nproc
      use multigridhelpers,   only: dirty_debug, check_dirty, dirtyL, multidim_code_3D
      use multigridmpifuncs,  only: mpi_multigrid_bnd

      implicit none

      integer, intent(in) :: lev  !< level for which approximate the solution
      integer, intent(in) :: src  !< index of source in lvl()%mgvar
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar

      integer :: nf, n

      do nf = 1, nsmoof
         lvl(lev)%src(:, :, :) = lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie, lvl(lev)%js:lvl(lev)%je, lvl(lev)%ks:lvl(lev)%ke, src)

         if (lvl(lev)%fft_type == fft_dst) then !correct boundaries on nonperiodic local domain
            if (nf == 1) then
               call make_face_boundaries(lev, soln)
            else
               call mpi_multigrid_bnd(lev, soln, 1, .false.)
               if (has_dir(XDIR)) then
                  lvl(lev)%bnd_x(:, :, LOW)  = 0.5* sum (lvl(lev)%mgvar(lvl(lev)%is-1:lvl(lev)%is, lvl(lev)%js:lvl(lev)%je, lvl(lev)%ks:lvl(lev)%ke, soln), 1)
                  lvl(lev)%bnd_x(:, :, HIGH) = 0.5* sum (lvl(lev)%mgvar(lvl(lev)%ie:lvl(lev)%ie+1, lvl(lev)%js:lvl(lev)%je, lvl(lev)%ks:lvl(lev)%ke, soln), 1)
               end if
               if (has_dir(YDIR)) then
                  lvl(lev)%bnd_y(:, :, LOW)  = 0.5* sum (lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie, lvl(lev)%js-1:lvl(lev)%js, lvl(lev)%ks:lvl(lev)%ke, soln), 2)
                  lvl(lev)%bnd_y(:, :, HIGH) = 0.5* sum (lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie, lvl(lev)%je:lvl(lev)%je+1, lvl(lev)%ks:lvl(lev)%ke, soln), 2)
               end if
               if (has_dir(ZDIR)) then
                  lvl(lev)%bnd_z(:, :, LOW)  = 0.5* sum (lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie, lvl(lev)%js:lvl(lev)%je, lvl(lev)%ks-1:lvl(lev)%ks, soln), 3)
                  lvl(lev)%bnd_z(:, :, HIGH) = 0.5* sum (lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie, lvl(lev)%js:lvl(lev)%je, lvl(lev)%ke:lvl(lev)%ke+1, soln), 3)
               end if
            end if

            if (dirty_debug) then
               if (has_dir(XDIR) .and. any(abs(lvl(lev)%bnd_x(:, :, :)) > dirtyL)) call warn("approximate_solution_fft dirty bnd_x")
               if (has_dir(YDIR) .and. any(abs(lvl(lev)%bnd_y(:, :, :)) > dirtyL)) call warn("approximate_solution_fft dirty bnd_y")
               if (has_dir(ZDIR) .and. any(abs(lvl(lev)%bnd_z(:, :, :)) > dirtyL)) call warn("approximate_solution_fft dirty bnd_z")
            end if

            if (has_dir(XDIR)) then
               lvl(lev)%src(1,            :, :) = lvl(lev)%src(1,            :, :) - lvl(lev)%bnd_x(:, :, LOW)  * 2. * lvl(lev)%idx2
               lvl(lev)%src(lvl(lev)%nxb, :, :) = lvl(lev)%src(lvl(lev)%nxb, :, :) - lvl(lev)%bnd_x(:, :, HIGH) * 2. * lvl(lev)%idx2
            end if
            if (has_dir(YDIR)) then
               lvl(lev)%src(:, 1,            :) = lvl(lev)%src(:, 1,            :) - lvl(lev)%bnd_y(:, :, LOW)  * 2. * lvl(lev)%idy2
               lvl(lev)%src(:, lvl(lev)%nyb, :) = lvl(lev)%src(:, lvl(lev)%nyb, :) - lvl(lev)%bnd_y(:, :, HIGH) * 2. * lvl(lev)%idy2
            end if
            if (has_dir(ZDIR)) then
               lvl(lev)%src(:, :, 1           ) = lvl(lev)%src(:, :, 1           ) - lvl(lev)%bnd_z(:, :, LOW)  * 2. * lvl(lev)%idz2
               lvl(lev)%src(:, :, lvl(lev)%nzb) = lvl(lev)%src(:, :, lvl(lev)%nzb) - lvl(lev)%bnd_z(:, :, HIGH) * 2. * lvl(lev)%idz2
            end if
         end if

         call fft_convolve(lev)

         lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie, lvl(lev)%js:lvl(lev)%je, lvl(lev)%ks:lvl(lev)%ke, soln) = lvl(lev)%src(:, :, :)

         call check_dirty(lev, soln, "approx_soln fft+")

         !BEWARE use has_dir() here in a way that does not degrade performance

         !relax the boundaries
         do n = 1, nsmool
            call mpi_multigrid_bnd(lev, soln, 1, .false.)
            ! Possible optimization: This is a quite costly part of the local FFT solver
            if (fft_full_relax) then
               if (eff_dim == NDIM .and. .not. multidim_code_3D) then
                  lvl                    (lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     soln)  = &
                       lvl(lev)%rx * (lvl(lev)%mgvar(lvl(lev)%is-1:lvl(lev)%ie-1, lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is+1:lvl(lev)%ie+1, lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     soln)) + &
                       lvl(lev)%ry * (lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js-1:lvl(lev)%je-1, lvl(lev)%ks:lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js+1:lvl(lev)%je+1, lvl(lev)%ks:lvl(lev)%ke,     soln)) + &
                       lvl(lev)%rz * (lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks-1:lvl(lev)%ke-1, soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks+1:lvl(lev)%ke+1, soln)) - &
                       lvl(lev)%r  *  lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     src)
               else

                  call die("[multigrid:approximate_solution_fft] fft_full_relax is allowed only for 3D at the moment")

                  ! An additional array (lvl(lev)%src would be good enough) is required here to assemble partial results or use red-black passes
                  lvl                    (lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     soln)  = &
                       - lvl(lev)%r * lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     src)
                  if (has_dir(XDIR)) &
                       lvl               (lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     soln)  = &
                       lvl               (lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     soln)  + &
                       lvl(lev)%rx * (lvl(lev)%mgvar(lvl(lev)%is-1:lvl(lev)%ie-1, lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is+1:lvl(lev)%ie+1, lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     soln))
                  if (has_dir(YDIR)) &
                       lvl               (lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     soln)  = &
                       lvl               (lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     soln)  + &
                       lvl(lev)%ry * (lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js-1:lvl(lev)%je-1, lvl(lev)%ks:lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js+1:lvl(lev)%je+1, lvl(lev)%ks:lvl(lev)%ke,     soln))
                  if (has_dir(ZDIR)) &
                       lvl               (lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     soln)  = &
                       lvl               (lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks:lvl(lev)%ke,     soln)  + &
                       lvl(lev)%rz * (lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks-1:lvl(lev)%ke-1, soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is:lvl(lev)%ie,     lvl(lev)%js:lvl(lev)%je,     lvl(lev)%ks+1:lvl(lev)%ke+1, soln))
               end if
            else
               ! relax only two layers of cells (1 is  significantly worse, 3 does not improve much)
               ! edges are relaxed twice, corners are relaxed three times which seems to be good

               if (has_dir(XDIR)) then
                  lvl                    (lev)%mgvar(lvl(lev)%is    :lvl(lev)%is+D_x,   lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ke,     soln)  = & ! -X
                       lvl(lev)%rx * (lvl(lev)%mgvar(lvl(lev)%is-D_x:lvl(lev)%is,       lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is+D_x:lvl(lev)%is+2*D_x, lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ke,     soln)) + &
                       lvl(lev)%ry * (lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%is+D_x,   lvl(lev)%js-D_y:lvl(lev)%je-D_y, lvl(lev)%ks    :lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%is+D_x,   lvl(lev)%js+D_y:lvl(lev)%je+D_y, lvl(lev)%ks    :lvl(lev)%ke,     soln)) + &
                       lvl(lev)%rz * (lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%is+D_x,   lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks-D_z:lvl(lev)%ke-D_z, soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%is+D_x,   lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks+D_z:lvl(lev)%ke+D_z, soln)) - &
                       lvl(lev)%r  *  lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%is+D_x,   lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ke,     src)

                  lvl                    (lev)%mgvar(lvl(lev)%ie-D_x  :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ke,     soln)  = & ! +X
                       lvl(lev)%rx * (lvl(lev)%mgvar(lvl(lev)%ie-2*D_x:lvl(lev)%ie-D_x, lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%ie      :lvl(lev)%ie+D_x, lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ke,     soln)) + &
                       lvl(lev)%ry * (lvl(lev)%mgvar(lvl(lev)%ie-D_x  :lvl(lev)%ie,     lvl(lev)%js-D_y:lvl(lev)%je-D_y, lvl(lev)%ks    :lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%ie-D_x  :lvl(lev)%ie,     lvl(lev)%js+D_y:lvl(lev)%je+D_y, lvl(lev)%ks    :lvl(lev)%ke,     soln)) + &
                       lvl(lev)%rz * (lvl(lev)%mgvar(lvl(lev)%ie-D_x  :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks-D_z:lvl(lev)%ke-D_z, soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%ie-D_x  :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks+D_z:lvl(lev)%ke+D_z, soln)) - &
                       lvl(lev)%r  *  lvl(lev)%mgvar(lvl(lev)%ie-D_x  :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ke,     src)
               end if

               if (has_dir(YDIR)) then
                  lvl                    (lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%js+D_y,   lvl(lev)%ks    :lvl(lev)%ke,     soln)  = & ! -Y
                       lvl(lev)%rx * (lvl(lev)%mgvar(lvl(lev)%is-D_x:lvl(lev)%ie-D_x, lvl(lev)%js    :lvl(lev)%js+D_y,   lvl(lev)%ks    :lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is+D_x:lvl(lev)%ie+D_x, lvl(lev)%js    :lvl(lev)%js+D_y,   lvl(lev)%ks    :lvl(lev)%ke,     soln)) + &
                       lvl(lev)%ry * (lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js-D_y:lvl(lev)%js,       lvl(lev)%ks    :lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js+D_y:lvl(lev)%js+2*D_y, lvl(lev)%ks    :lvl(lev)%ke,     soln)) + &
                       lvl(lev)%rz * (lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%js+D_y,   lvl(lev)%ks-D_z:lvl(lev)%ke-D_z, soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%js+D_y,   lvl(lev)%ks+D_z:lvl(lev)%ke+D_z, soln)) - &
                       lvl(lev)%r  *  lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%js+D_y,   lvl(lev)%ks    :lvl(lev)%ke,     src)

                  lvl                    (lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%je-D_y  :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ke,     soln)  = & ! +Y
                       lvl(lev)%rx * (lvl(lev)%mgvar(lvl(lev)%is-D_x:lvl(lev)%ie-D_x, lvl(lev)%je-D_y  :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is+D_x:lvl(lev)%ie+D_x, lvl(lev)%je-D_y  :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ke,     soln)) + &
                       lvl(lev)%ry * (lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%je-2*D_y:lvl(lev)%je-D_y, lvl(lev)%ks    :lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%je      :lvl(lev)%je+D_y, lvl(lev)%ks    :lvl(lev)%ke,     soln)) + &
                       lvl(lev)%rz * (lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%je-D_y  :lvl(lev)%je,     lvl(lev)%ks-D_z:lvl(lev)%ke-D_z, soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%je-D_y  :lvl(lev)%je,     lvl(lev)%ks+D_z:lvl(lev)%ke+D_z, soln)) - &
                       lvl(lev)%r  *  lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%je-D_y  :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ke,     src)
               end if

               if (has_dir(ZDIR)) then
                  lvl                    (lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ks+D_z,   soln)  = & ! -Z
                       lvl(lev)%rx * (lvl(lev)%mgvar(lvl(lev)%is-D_x:lvl(lev)%ie-D_x, lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ks+D_z,   soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is+D_x:lvl(lev)%ie+D_x, lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ks+D_z,   soln)) + &
                       lvl(lev)%ry * (lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js-D_y:lvl(lev)%je-D_y, lvl(lev)%ks    :lvl(lev)%ks+D_z,   soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js+D_y:lvl(lev)%je+D_y, lvl(lev)%ks    :lvl(lev)%ks+D_z,   soln)) + &
                       lvl(lev)%rz * (lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks-D_z:lvl(lev)%ks,       soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks+D_z:lvl(lev)%ks+2*D_z, soln)) - &
                       lvl(lev)%r  *  lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ks    :lvl(lev)%ks+D_z,   src)

                  lvl                    (lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ke-D_z  :lvl(lev)%ke ,    soln)  = & ! +Z
                       lvl(lev)%rx * (lvl(lev)%mgvar(lvl(lev)%is-D_x:lvl(lev)%ie-D_x, lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ke-D_z  :lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is+D_x:lvl(lev)%ie+D_x, lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ke-D_z  :lvl(lev)%ke,     soln)) + &
                       lvl(lev)%ry * (lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js-D_y:lvl(lev)%je-D_y, lvl(lev)%ke-D_z  :lvl(lev)%ke,     soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js+D_y:lvl(lev)%je+D_y, lvl(lev)%ke-D_z  :lvl(lev)%ke,     soln)) + &
                       lvl(lev)%rz * (lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ke-2*D_z:lvl(lev)%ke-D_z, soln)  + &
                       &              lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ke      :lvl(lev)%ke+D_z, soln)) - &
                       lvl(lev)%r  *  lvl(lev)%mgvar(lvl(lev)%is    :lvl(lev)%ie,     lvl(lev)%js    :lvl(lev)%je,     lvl(lev)%ke-D_z  :lvl(lev)%ke,     src)
               end if

            end if
         end do

         call check_dirty(lev, soln, "approx_soln relax+")

      end do

   end subroutine approximate_solution_fft

!!$ ============================================================================
!!
!! This routine prepares boundary values for local-FFT solver
!!

   subroutine make_face_boundaries(lev, soln)

      use mpisetup,           only: nproc
      use multigridbasefuncs, only: zero_boundaries
      use errh,               only: warn

      implicit none

      integer, intent(in) :: lev  !< level for which approximate the solution
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar

      if (grav_bnd == bnd_periodic .and. nproc == 1) then
         call zero_boundaries(lev)
      else
         if (lev > level_min) then
            call prolong_faces(lev, soln)
         else
            if (grav_bnd /= bnd_givenval) call zero_boundaries(lev)
            call warn("m:mfb WTF?")
         end if
      end if

   end subroutine make_face_boundaries

!!$ ============================================================================
!!
!! Prolong solution data at level (lev-1) to faces at level lev
!!

   subroutine prolong_faces(lev, soln)

      use mpisetup,           only: proc
      use errh,               only: die, warn
      use multigridhelpers,   only: check_dirty
      use multigridmpifuncs,  only: mpi_multigrid_bnd

      implicit none

      integer, intent(in) :: lev  !< level for which approximate the solution
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar

      integer                       :: i, j, k
      type(plvl), pointer           :: coarse, fine
      real, parameter, dimension(3) :: p0  = [ 0.,       1.,     0.     ] ! injection
      real, parameter, dimension(3) :: p1  = [ 0.,       3./4.,  1./4.  ] ! 1D linear prolongation stencil
      real, parameter, dimension(3) :: p2i = [ -1./8.,   1.,     1./8.  ] ! 1D integral cubic prolongation stencil
      real, parameter, dimension(3) :: p2d = [ -3./32., 30./32., 5./32. ] ! 1D direct cubic prolongation stencil
      real, dimension(-1:1)         :: p
      real, dimension(-1:1,-1:1,2,2):: pp   ! 2D prolongation stencil
      real                          :: pp_norm

      if (lev < level_min .or. lev > level_max) call die("[multigrid:prolong_faces] Invalid level")

      if (lev == level_min) then
         call warn("[multigrid:prolong_faces] Cannot prolong anything to base level")
         return
      end if

      select case(ord_prolong_face)
         case(0)
            p(:) = p0(:)
         case(1,-1)
            p(:) = p1(:)
         case(2)
            p(:) = p2i(:)
         case(-2)
            p(:) = p2d(:)
         case default
            p(:) = p0(:)
      end select

      do i = -1, 1
         pp(i,:,1,1) = 0.5*p( i)*p(:)       ! 0.5 because of face averaging
         pp(i,:,1,2) = 0.5*p( i)*p(1:-1:-1) ! or use matmul()
         pp(i,:,2,1) = 0.5*p(-i)*p(:)
         pp(i,:,2,2) = 0.5*p(-i)*p(1:-1:-1)
      end do

      call mpi_multigrid_bnd(lev-1, soln, 1, .false.) !BEWARE for higher prolongation order more guardcell will be required
      call check_dirty(lev-1, soln, "prolong_faces", 1)

      coarse => lvl(lev - 1)
      fine   => lvl(lev)

      if (has_dir(XDIR)) then
         pp_norm = 2.*sum(pp(-D_y:D_y, -D_z:D_z, 1, 1)) ! normalization is required for ord_prolong_face == 1 and -2
         do j = coarse%js, coarse%je
            do k = coarse%ks, coarse%ke
               fine%bnd_x(-fine%js+2*j,    -fine%ks+2*k,    LOW) =sum(pp(-D_y:D_y, -D_z:D_z, 1, 1) * (coarse%mgvar(coarse%is,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_x(-fine%js+2*j+D_y,-fine%ks+2*k,    LOW) =sum(pp(-D_y:D_y, -D_z:D_z, 2, 1) * (coarse%mgvar(coarse%is,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_x(-fine%js+2*j,    -fine%ks+2*k+D_z,LOW) =sum(pp(-D_y:D_y, -D_z:D_z, 1, 2) * (coarse%mgvar(coarse%is,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_x(-fine%js+2*j+D_y,-fine%ks+2*k+D_z,LOW) =sum(pp(-D_y:D_y, -D_z:D_z, 2, 2) * (coarse%mgvar(coarse%is,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%is-1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_x(-fine%js+2*j,    -fine%ks+2*k,    HIGH)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 1) * (coarse%mgvar(coarse%ie,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_x(-fine%js+2*j+D_y,-fine%ks+2*k,    HIGH)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 1) * (coarse%mgvar(coarse%ie,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_x(-fine%js+2*j,    -fine%ks+2*k+D_z,HIGH)=sum(pp(-D_y:D_y, -D_z:D_z, 1, 2) * (coarse%mgvar(coarse%ie,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_x(-fine%js+2*j+D_y,-fine%ks+2*k+D_z,HIGH)=sum(pp(-D_y:D_y, -D_z:D_z, 2, 2) * (coarse%mgvar(coarse%ie,j-D_y:j+D_y,k-D_z:k+D_z,soln) + coarse%mgvar(coarse%ie+1,j-D_y:j+D_y,k-D_z:k+D_z,soln))) / pp_norm
            end do
         end do
      end if

      if (has_dir(YDIR)) then
         pp_norm = 2.*sum(pp(-D_x:D_x, -D_z:D_z, 1, 1))
         do i = coarse%is, coarse%ie
            do k = coarse%ks, coarse%ke
               fine%bnd_y(-fine%is+2*i,    -fine%ks+2*k,    LOW) =sum(pp(-D_x:D_x, -D_z:D_z, 1, 1) * (coarse%mgvar(i-D_x:i+D_x,coarse%js,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-1,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_y(-fine%is+2*i+D_x,-fine%ks+2*k,    LOW) =sum(pp(-D_x:D_x, -D_z:D_z, 2, 1) * (coarse%mgvar(i-D_x:i+D_x,coarse%js,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-1,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_y(-fine%is+2*i,    -fine%ks+2*k+D_z,LOW) =sum(pp(-D_x:D_x, -D_z:D_z, 1, 2) * (coarse%mgvar(i-D_x:i+D_x,coarse%js,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-1,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_y(-fine%is+2*i+D_x,-fine%ks+2*k+D_z,LOW) =sum(pp(-D_x:D_x, -D_z:D_z, 2, 2) * (coarse%mgvar(i-D_x:i+D_x,coarse%js,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%js-1,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_y(-fine%is+2*i,    -fine%ks+2*k,    HIGH)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 1) * (coarse%mgvar(i-D_x:i+D_x,coarse%je,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+1,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_y(-fine%is+2*i+D_x,-fine%ks+2*k,    HIGH)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 1) * (coarse%mgvar(i-D_x:i+D_x,coarse%je,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+1,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_y(-fine%is+2*i,    -fine%ks+2*k+D_z,HIGH)=sum(pp(-D_x:D_x, -D_z:D_z, 1, 2) * (coarse%mgvar(i-D_x:i+D_x,coarse%je,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+1,k-D_z:k+D_z,soln))) / pp_norm
               fine%bnd_y(-fine%is+2*i+D_x,-fine%ks+2*k+D_z,HIGH)=sum(pp(-D_x:D_x, -D_z:D_z, 2, 2) * (coarse%mgvar(i-D_x:i+D_x,coarse%je,k-D_z:k+D_z,soln) + coarse%mgvar(i-D_x:i+D_x,coarse%je+1,k-D_z:k+D_z,soln))) / pp_norm
            end do
         end do
      end if

      if (has_dir(ZDIR)) then
         pp_norm = 2.*sum(pp(-D_x:D_x, -D_y:D_y, 1, 1))
         do i = coarse%is, coarse%ie
            do j = coarse%js, coarse%je
               fine%bnd_z(-fine%is+2*i,    -fine%js+2*j,    LOW) =sum(pp(-D_x:D_x, -D_y:D_y, 1, 1) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-1,soln))) / pp_norm
               fine%bnd_z(-fine%is+2*i+D_x,-fine%js+2*j,    LOW) =sum(pp(-D_x:D_x, -D_y:D_y, 2, 1) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-1,soln))) / pp_norm
               fine%bnd_z(-fine%is+2*i,    -fine%js+2*j+D_y,LOW) =sum(pp(-D_x:D_x, -D_y:D_y, 1, 2) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-1,soln))) / pp_norm
               fine%bnd_z(-fine%is+2*i+D_x,-fine%js+2*j+D_y,LOW) =sum(pp(-D_x:D_x, -D_y:D_y, 2, 2) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ks-1,soln))) / pp_norm
               fine%bnd_z(-fine%is+2*i,    -fine%js+2*j,    HIGH)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 1) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+1,soln))) / pp_norm
               fine%bnd_z(-fine%is+2*i+D_x,-fine%js+2*j,    HIGH)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 1) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+1,soln))) / pp_norm
               fine%bnd_z(-fine%is+2*i,    -fine%js+2*j+D_y,HIGH)=sum(pp(-D_x:D_x, -D_y:D_y, 1, 2) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+1,soln))) / pp_norm
               fine%bnd_z(-fine%is+2*i+D_x,-fine%js+2*j+D_y,HIGH)=sum(pp(-D_x:D_x, -D_y:D_y, 2, 2) * (coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke,soln) + coarse%mgvar(i-D_x:i+D_x,j-D_y:j+D_y,coarse%ke+1,soln))) / pp_norm
            end do
         end do
      end if

   end subroutine prolong_faces

!!$ ============================================================================
!!
!! Red-Black Gauss-Seidel relaxation.
!!
!! This is the most costly routine in a serial run. Try to find optimal values for nsmool and nsmoob.
!! This routine also depends a lot on communication so it  may limit scalability of the multigrid.
!! \todo Implement convergence check on base level (not very important since we have a FFT solver for base level)
!!

   subroutine approximate_solution_rbgs(lev, src, soln)

      use multigridhelpers,   only: dirty_debug, check_dirty, multidim_code_3D, dirty_label
      use multigridmpifuncs,  only: mpi_multigrid_bnd
      use errh,               only: die

      implicit none

      integer, intent(in) :: lev  !< level for which approximate the solution
      integer, intent(in) :: src  !< index of source in lvl()%mgvar
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar

      integer, parameter :: RED_BLACK = 2 !< the checkerboard requires two sweeps

      integer :: n, j, k, i1, j1, k1, id, jd, kd
      integer :: nsmoo

      if (lev == level_min) then
         nsmoo = nsmoob
      else
         nsmoo = nsmool
      end if

      do n = 1, RED_BLACK*nsmoo
         call mpi_multigrid_bnd(lev, soln, 1, .false.) ! no corners are required here

         if (dirty_debug) then
            write(dirty_label, '(a,i5)')"relax soln- smoo=", n
            call check_dirty(lev, soln, dirty_label)
         end if

         ! Possible optimization: this is the most costly part of the RBGS relaxation (instruction count, read and write data, L1 and L2 read cache miss)
         ! do n = 1, nsmoo
         !    call mpi_multigrid_bnd(lev, soln, 1, .false.)
         !    relax single layer of red cells at all faces
         !    call mpi_multigrid_bnd(lev, soln, 1, .false.)
         !    relax interior cells (except for single layer of cells at all faces), first red, then 1-cell behind black one.
         !    relax single layer of black cells at all faces
         ! end do

         ! with explicit outer loops it is easier to describe a 3-D checkerboard :-)

         if (eff_dim==NDIM .and. .not. multidim_code_3D) then
            do k = lvl(lev)%ks, lvl(lev)%ke
               do j = lvl(lev)%js, lvl(lev)%je
                  i1 = lvl(lev)%is + mod(n+j+k, RED_BLACK)
                  lvl(          lev          )%mgvar(i1  :lvl(lev)%ie  :2, j,   k,   soln) = &
                       lvl(lev)%rx * (lvl(lev)%mgvar(i1-1:lvl(lev)%ie-1:2, j,   k,   soln) + lvl(lev)%mgvar(i1+1:lvl(lev)%ie+1:2, j,   k,   soln)) + &
                       lvl(lev)%ry * (lvl(lev)%mgvar(i1  :lvl(lev)%ie  :2, j-1, k,   soln) + lvl(lev)%mgvar(i1:  lvl(lev)%ie:  2, j+1, k,   soln)) + &
                       lvl(lev)%rz * (lvl(lev)%mgvar(i1  :lvl(lev)%ie  :2, j,   k-1, soln) + lvl(lev)%mgvar(i1:  lvl(lev)%ie:  2, j,   k+1, soln)) - &
                       lvl(lev)%r  *  lvl(lev)%mgvar(i1  :lvl(lev)%ie  :2, j,   k,   src)
               end do
            end do
         else
            ! In 3D this variant significantly increases instruction count and also some data read
            i1 = lvl(lev)%is; id = 1 ! mv to multigridvars, init_multigrid
            j1 = lvl(lev)%js; jd = 1
            k1 = lvl(lev)%ks; kd = 1
            if (has_dir(XDIR)) then
               id = RED_BLACK
            else if (has_dir(YDIR)) then
               jd = RED_BLACK
            else if (has_dir(ZDIR)) then
               kd = RED_BLACK
            end if

            if (kd == RED_BLACK) k1 = lvl(lev)%ks + mod(n, RED_BLACK)
            do k = k1, lvl(lev)%ke, kd
               if (jd == RED_BLACK) j1 = lvl(lev)%js + mod(n+k, RED_BLACK)
               do j = j1, lvl(lev)%je, jd
                  if (id == RED_BLACK) i1 = lvl(lev)%is + mod(n+j+k, RED_BLACK)
                  lvl(      lev)%mgvar(i1  :lvl(lev)%ie  :id, j,   k,   soln) = &
                       & (1. - Jacobi_damp)* lvl(lev)%mgvar(i1  :lvl(lev)%ie  :id, j,   k,   soln) &
                       &     - Jacobi_damp * lvl(lev)%mgvar(i1  :lvl(lev)%ie  :id, j,   k,   src)  * lvl(lev)%r
                  if (has_dir(XDIR)) &
                       lvl (lev)%mgvar(i1  :lvl(lev)%ie  :id, j,   k,   soln) = lvl(lev)%mgvar(i1:  lvl(lev)%ie:  id, j,   k,   soln)  + &
                       &       Jacobi_damp *(lvl(lev)%mgvar(i1-1:lvl(lev)%ie-1:id, j,   k,   soln) + lvl(lev)%mgvar(i1+1:lvl(lev)%ie+1:id, j,   k,   soln)) * lvl(lev)%rx
                  if (has_dir(YDIR)) &
                       lvl (lev)%mgvar(i1  :lvl(lev)%ie  :id, j,   k,   soln) = lvl(lev)%mgvar(i1:  lvl(lev)%ie:  id, j,   k,   soln)  + &
                       &       Jacobi_damp *(lvl(lev)%mgvar(i1  :lvl(lev)%ie  :id, j-1, k,   soln) + lvl(lev)%mgvar(i1:  lvl(lev)%ie:  id, j+1, k,   soln)) * lvl(lev)%ry
                  if (has_dir(ZDIR)) &
                       lvl (lev)%mgvar(i1  :lvl(lev)%ie  :id, j,   k,   soln) = lvl(lev)%mgvar(i1:  lvl(lev)%ie:  id, j,   k,   soln)  + &
                       &       Jacobi_damp *(lvl(lev)%mgvar(i1  :lvl(lev)%ie  :id, j,   k-1, soln) + lvl(lev)%mgvar(i1:  lvl(lev)%ie:  id, j,   k+1, soln)) * lvl(lev)%rz
               end do
            end do
         end if

         if (dirty_debug) then
            write(dirty_label, '(a,i5)')"relax soln+ smoo=", n
            call check_dirty(lev, soln, dirty_label)
         end if

      end do

   end subroutine approximate_solution_rbgs

!!$ ============================================================================
!!
!! Obtain an exact solution on base level. (wrapper)
!!
   subroutine gb_fft_solve(src, soln)
      implicit none

      integer, intent(in) :: src  !< index of source in lvl()%mgvar
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar

      if (gb_solve_gather) then
         call gb_fft_solve_gather(src, soln)
      else
         call gb_fft_solve_sendrecv(src, soln)
      end if

   end subroutine gb_fft_solve

!!$ ============================================================================
!!
!! Obtain an exact solution on base level.
!! The source is gathered on the master PE (possible bottleneck), solution is obtained through FFT, then its parts are sent to all PEs.
!! Unfortunately non-blocking communication probably will not change much here.
!! Alternative implementation: set gb_src to 0, copy local part to it, then call MPI_Allreduce and solve with FFT on each PE (no need to communicate the solution)
!!

   subroutine gb_fft_solve_sendrecv(src, soln)

      use mpisetup, only: nproc, proc, ierr, comm3d, status, MPI_DOUBLE_PRECISION

      implicit none

      integer, intent(in) :: src  !< index of source in lvl()%mgvar
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar

      integer :: p

      if (proc == 0) then

         ! collect the source on the base level
         gb%src(gb_cartmap(0)%lo(XDIR):gb_cartmap(0)%up(XDIR), &
              & gb_cartmap(0)%lo(YDIR):gb_cartmap(0)%up(YDIR), &
              & gb_cartmap(0)%lo(ZDIR):gb_cartmap(0)%up(ZDIR)) = &
              base%mgvar(base%is:base%ie, base%js:base%je, base%ks:base%ke, src)
         do p = 1, nproc -1
            call MPI_Recv(gb%src(gb_cartmap(p)%lo(XDIR):gb_cartmap(p)%up(XDIR), &
                 &               gb_cartmap(p)%lo(YDIR):gb_cartmap(p)%up(YDIR), &
                 &               gb_cartmap(p)%lo(ZDIR):gb_cartmap(p)%up(ZDIR)), &
                 &        base%nxb*base%nyb*base%nzb, MPI_DOUBLE_PRECISION, p, p, comm3d, status, ierr)
         end do

         call fft_convolve(gb%level)

         ! send all the pieces away
         base%mgvar(base%is:base%ie, base%js:base%je, base%ks:base%ke, soln) = &
              gb%src(gb_cartmap(0)%lo(XDIR):gb_cartmap(0)%up(XDIR), &
              &      gb_cartmap(0)%lo(YDIR):gb_cartmap(0)%up(YDIR), &
              &      gb_cartmap(0)%lo(ZDIR):gb_cartmap(0)%up(ZDIR))
         do p = 1, nproc -1
            call MPI_Send(gb%src(gb_cartmap(p)%lo(XDIR):gb_cartmap(p)%up(XDIR), &
                 &               gb_cartmap(p)%lo(YDIR):gb_cartmap(p)%up(YDIR), &
                 &               gb_cartmap(p)%lo(ZDIR):gb_cartmap(p)%up(ZDIR)), &
                 &        base%nxb*base%nyb*base%nzb, MPI_DOUBLE_PRECISION, p, p, comm3d, ierr)
         end do

      else
         call MPI_Send(base%mgvar(base%is:base%ie, base%js:base%je, base%ks:base%ke, src),  base%nxb*base%nyb*base%nzb, MPI_DOUBLE_PRECISION, 0, proc, comm3d, ierr)
         call MPI_Recv(base%mgvar(base%is:base%ie, base%js:base%je, base%ks:base%ke, soln), base%nxb*base%nyb*base%nzb, MPI_DOUBLE_PRECISION, 0, proc, comm3d, status, ierr)
      end if

   end subroutine gb_fft_solve_sendrecv

!!$ ============================================================================
!!
!! Alternative version using MPI_Gather
!!
!! \todo using 1d ffts and gather to reduce dimensions one by one could prevent possible bottleneck and could be (?) faster:
!!       1. N_x * N_y independent gathers reduce Z dim
!!       2. fft in Z dim on N_x*N_y proc
!!       3. N_x independent gathers of (2.) reduce X dim
!!       4. fft in Y dim on N_x
!!       5. gather (4.) on master proc
!!       6. fft in X dim on master
!!       7. Scatter

   subroutine gb_fft_solve_gather(src, soln)

      use mpisetup, only: nproc, proc, ierr, comm3d, status, MPI_DOUBLE_PRECISION

      implicit none

      integer, intent(in) :: src  !< index of source in lvl()%mgvar
      integer, intent(in) :: soln !< index of solution in lvl()%mgvar

      integer :: p

      call MPI_Gather(base%mgvar(base%is:base%ie, base%js:base%je, base%ks:base%ke, src),  base%nxb*base%nyb*base%nzb, MPI_DOUBLE_PRECISION, &
           &          gb_src_temp, base%nxb*base%nyb*base%nzb, MPI_DOUBLE_PRECISION, 0, comm3d, ierr)

      if (proc == 0) then
         do p = 0, nproc-1
            gb%src(gb_cartmap(p)%lo(XDIR):gb_cartmap(p)%up(XDIR), &
                   gb_cartmap(p)%lo(YDIR):gb_cartmap(p)%up(YDIR), &
                   gb_cartmap(p)%lo(ZDIR):gb_cartmap(p)%up(ZDIR)) = gb_src_temp(:,:,:,p)
         enddo

         call fft_convolve(gb%level)

         do p = 0, nproc-1
            gb_src_temp(:,:,:,p) = &
               gb%src(gb_cartmap(p)%lo(XDIR):gb_cartmap(p)%up(XDIR), &
                      gb_cartmap(p)%lo(YDIR):gb_cartmap(p)%up(YDIR), &
                      gb_cartmap(p)%lo(ZDIR):gb_cartmap(p)%up(ZDIR))
         enddo
      endif

      call MPI_Scatter(gb_src_temp,  base%nxb*base%nyb*base%nzb, MPI_DOUBLE_PRECISION, &
                       base%mgvar(base%is:base%ie, base%js:base%je, base%ks:base%ke, soln), base%nxb*base%nyb*base%nzb, MPI_DOUBLE_PRECISION, &
                       0, comm3d, ierr)

   end subroutine gb_fft_solve_gather

!!$ ============================================================================
!!
!! Do the FFT convolution
!!

   subroutine fft_convolve(level)

      use errh, only: die

      implicit none

      integer, intent(in) :: level !< level at which make the convolution

      ! do the convolution in Fourier space; lvl(level)%src(:,:,:) -> lvl(level)%fft{r}(:,:,:)
      call dfftw_execute(lvl(level)%planf)

      select case (lvl(level)%fft_type)
         case (fft_rcr)
            lvl(level)%fft  = lvl(level)%fft  * lvl(level)%Green3D
         case (fft_dst)
            lvl(level)%fftr = lvl(level)%fftr * lvl(level)%Green3D
         case default
            call die("[multigrid:fft_convolve] Unknown FFT type.")
      end select

      call dfftw_execute(lvl(level)%plani) ! lvl(level)%fft{r}(:,:,:) -> lvl(level)%src(:,:,:)

    end subroutine fft_convolve

#else /* MULTIGRID */
#warning This should not happen. Probably the multigrid.F90 file is included in object directory by mistake.
#endif /* MULTIGRID */

end module multigrid

!!$ ============================================================================
