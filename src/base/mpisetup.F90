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
#include "macros.h"
!>
!! \brief (KK)
!!
!! In this module following namelists of parameters are specified:
!! \copydetails mpisetup::init_mpi
!<
module mpisetup

   use types,      only: domain_container
   use constants,  only: ndims, cbuff_len, LO, HI, BLK, BND

   implicit none

   private
   public :: cleanup_mpi, init_mpi, mpifind, translate_bnds_to_ints_dom, is_neigh, &
        &    buffer_dim, cbuff, ibuff, lbuff, rbuff, comm, comm3d, req, status, ierr, info, &
        &    master, slave, have_mpi, has_dir, eff_dim, is_uneven, is_mpi_noncart, &
        &    nproc, pcoords, proc, procn, psize, procxyl, procyxl, &
        &    dom, geometry_type, procmask, &
        &    cfl, cfl_max, cflcontrol, cfl_violated, &
        &    dt, dt_initial, dt_max_grow, dt_min, dt_old, dtm, t, nstep, &
        &    integration_order, limiter, smalld, smallei, smallp, use_smalld, magic_mass, local_magic_mass, &
        &    relax_time, grace_period_passed, cfr_smooth, repeat_step

   integer, protected :: nproc, proc, ierr, info

   integer, parameter :: nreq = size([LO, HI]) * size([BLK, BND]) ! just another way of defining '4' ;-)
   integer, allocatable, dimension(:)   :: req
   integer, allocatable, dimension(:,:) :: status

   logical, protected    :: master, slave

   real, parameter       :: dt_default_grow = 2.
   logical               :: cfl_violated             !< True when cfl condition is violated
   real                  :: t, dt, dt_old, dtm
   real, save            :: magic_mass = 0.0
   real, save            :: local_magic_mass = 0.0
   integer               :: nstep

   integer, protected        :: comm, comm3d
   integer, dimension(ndims), protected :: pcoords
   integer, protected :: procxyl, procyxl
   integer, protected, dimension(ndims, LO:HI) :: procn   !< array of neighbours proc numbers
   logical, protected, dimension(ndims) :: has_dir   !< .true. for existing directions
   integer, protected    :: eff_dim                  !< effective dimensionality of the simulation

   type(domain_container), protected :: dom

   integer, parameter :: buffer_dim = 200                !< size of [cilr]buff arrays used to exchange namelist parameters
   character(len=cbuff_len), dimension(buffer_dim) :: cbuff
   integer,                  dimension(buffer_dim) :: ibuff
   real,                     dimension(buffer_dim) :: rbuff
   logical,                  dimension(buffer_dim) :: lbuff

   logical, protected :: have_mpi           !< .true. when run on more than one processor
   logical, protected :: is_uneven          !< .true. when n[xyz]b depend on process number
   logical, protected :: is_mpi_noncart     !< .true. when there exist a process that has more than one neighbour in any direction

   integer, allocatable, dimension(:) :: primes
   real :: ideal_bsize
   integer, protected :: geometry_type  !< the type of geometry: cartesian: GEO_XYZ, cylindrical: GEO_RPZ, other: GEO_INVALID
   integer, dimension(:), allocatable :: procmask !< (0:nproc-1)-sized auxiliary array for searching overlaps, neighbours etc. BEWARE: antiparallel

   ! Namelist variables

   integer, dimension(ndims), protected :: psize !< desired number of MPI blocks in x, y and z-dimension
   logical :: reorder                 !< allows processes reordered for efficiency (a parameter of MPI_Cart_create and MPI_graph_create)
   logical :: allow_uneven    !< allows different values of n[xyz]b on divverent processes
   logical :: allow_noncart   !< allows more than one neighbour on a boundary
   real    :: dd_unif_quality !< uniform domain decomposition may be rejected it its quality is below this threshold (e.g. very elongated local domains are found)
   real    :: dd_rect_quality !< rectilinear domain decomposition may be rejected it its quality is below this threshold (not used yet)
   logical :: use_comm3d      !< If .false. then do not call any MPI_Cart_* functions

   namelist /MPI_BLOCKS/ psize, reorder, allow_uneven, allow_noncart, dd_unif_quality, dd_rect_quality, use_comm3d

   integer, protected :: nxd  !< number of %grid cells in physical domain (without boundary cells) in x-direction (if == 1 then x-dimension is reduced to a point with no boundary cells)
   integer, protected :: nyd  !< number of %grid cells in physical domain (without boundary cells) in y-direction (-- || --)
   integer, protected :: nzd  !< number of %grid cells in physical domain (without boundary cells) in z-direction (-- || --)
   integer, protected :: nb   !< number of boundary cells surrounding the physical domain, same for all directions
   character(len=cbuff_len) :: bnd_xl     !< type of boundary conditions for the left  x-boundary
   character(len=cbuff_len) :: bnd_xr     !< type of boundary conditions for the right x-boundary
   character(len=cbuff_len) :: bnd_yl     !< type of boundary conditions for the left  y-boundary
   character(len=cbuff_len) :: bnd_yr     !< type of boundary conditions for the right y-boundary
   character(len=cbuff_len) :: bnd_zl     !< type of boundary conditions for the left  z-boundary
   character(len=cbuff_len) :: bnd_zr     !< type of boundary conditions for the right z-boundary
   real    :: xmin                           !< physical domain left x-boundary position
   real    :: xmax                           !< physical domain right x-boundary position
   real    :: ymin                           !< physical domain left y-boundary position
   real    :: ymax                           !< physical domain right y-boundary position
   real    :: zmin                           !< physical domain left z-boundary position
   real    :: zmax                           !< physical domain right z-boundary position
   character(len=cbuff_len) :: geometry      !< define system of coordinates: "cartesian" or "cylindrical"

   namelist /DOMAIN/ nxd, nyd, nzd, nb, bnd_xl, bnd_xr, bnd_yl, bnd_yr, bnd_zl, bnd_zr, xmin, xmax, ymin, ymax, zmin, zmax, geometry
   !namelist /DOMAIN/ dom, geometry, nb

   real    :: dt_initial               !< initial timestep
   real    :: dt_max_grow              !< maximum timestep growth rate
   real    :: dt_min                   !< minimum allowed timestep
   real    :: cfl                      !< desired Courant–Friedrichs–Lewy number
   real    :: cfl_max                  !< warning threshold for the effective CFL number achieved
   logical :: use_smalld               !< correct denisty when it gets lower than smalld
   real    :: smallp                   !< artificial infimum for pressure
   real    :: smalld                   !< artificial infimum for density
   real    :: smallc                   !< artificial infimum for freezing speed
   real    :: smallei                  !< artificial infimum for internal energy density
   !>
   !! small number used to smooth freezing speed, especially handy in dust with random noise in velocity field.
   !! \f$c_{\textrm{fr}} = \sqrt{v^2 + \frac{1}{2}(\max{v} - \min{v})c_{\textrm{fr}}^{\textrm{smooth}}} + \ldots\f$
   !<
   real    :: cfr_smooth
   integer, protected :: integration_order           !< Runge-Kutta time integration order (1 - 1st order, 2 - 2nd order)
   character(len=cbuff_len) :: limiter     !< type of flux limiter
   character(len=cbuff_len) :: cflcontrol  !< type of cfl control just before each sweep (possibilities: 'none', 'main', 'user')
   logical                  :: repeat_step !< repeat fluid step if cfl condition is violated (significantly increases mem usage)
   real    :: relax_time                   !< relaxation/grace time, additional physics will be turned off until mpisetup::t >= mpisetup::relax_time

   namelist /NUMERICAL_SETUP/  cfl, smalld, smallei, integration_order, cfr_smooth, dt_initial, dt_max_grow, dt_min, smallc, smallp, limiter, cflcontrol, use_smalld, cfl_max, relax_time, repeat_step

   interface is_neigh
      module procedure is_neigh_simple, is_neigh_per
   end interface

contains

!-----------------------------------------------------------------------------
!>
!! \brief Routine to start MPI
!!
!! \n \n
!! @b MPI_BLOCKS
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>psize(3)       </td><td>1      </td><td>integer</td><td>\copydoc mpisetup::psize          </td></tr>
!! <tr><td>reorder        </td><td>.false.</td><td>logical</td><td>\copydoc mpisetup::reorder        </td></tr>
!! <tr><td>allow_uneven   </td><td>.false.</td><td>logical</td><td>\copydoc mpisetup::allow_uneven   </td></tr>
!! <tr><td>allow_noncart  </td><td>.false.</td><td>logical</td><td>\copydoc mpisetup::allow_noncart  </td></tr>
!! <tr><td>dd_unif_quality</td><td>0.9    </td><td>real   </td><td>\copydoc mpisetup::dd_unif_quality</td></tr>
!! <tr><td>dd_rect_quality</td><td>0.9    </td><td>real   </td><td>\copydoc mpisetup::dd_rect_quality</td></tr>
!! </table>
!! \n \n
!! @b DOMAIN
!! \n \n
!! <table border="+1">
!!   <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!!   <tr><td>nxd</td><td>1</td><td>positive integer    </td><td>\copydoc mpisetup::nxd</td></tr>
!!   <tr><td>nyd</td><td>1</td><td>positive integer    </td><td>\copydoc mpisetup::nyd</td></tr>
!!   <tr><td>nzd</td><td>1</td><td>positive integer    </td><td>\copydoc mpisetup::nzd</td></tr>
!!   <tr><td>nb </td><td>4</td><td>non-negative integer</td><td>\copydoc mpisetup::nb </td></tr>
!!   <tr><td>bnd_xl</td><td>'per'</td><td>'per', 'ref', 'out', 'outd', 'outh', 'cor'</td><td>\copydoc mpisetup::bnd_xl</td></tr>
!!   <tr><td>bnd_xr</td><td>'per'</td><td>'per', 'ref', 'out', 'outd', 'outh'       </td><td>\copydoc mpisetup::bnd_xr</td></tr>
!!   <tr><td>bnd_yl</td><td>'per'</td><td>'per', 'ref', 'out', 'outd', 'outh', 'cor'</td><td>\copydoc mpisetup::bnd_yl</td></tr>
!!   <tr><td>bnd_yr</td><td>'per'</td><td>'per', 'ref', 'out', 'outd', 'outh'       </td><td>\copydoc mpisetup::bnd_yr</td></tr>
!!   <tr><td>bnd_zl</td><td>'per'</td><td>'per', 'ref', 'out', 'outd', 'outh'       </td><td>\copydoc mpisetup::bnd_zl</td></tr>
!!   <tr><td>bnd_zr</td><td>'per'</td><td>'per', 'ref', 'out', 'outd', 'outh'       </td><td>\copydoc mpisetup::bnd_zr</td></tr>
!!   <tr><td> xmin     </td><td> 0.          </td><td> real                     </td><td> physical domain left x-boundary position  </td></tr>
!!   <tr><td> xmax     </td><td> 1.          </td><td> real                     </td><td> physical domain right x-boundary position </td></tr>
!!   <tr><td> ymin     </td><td> 0.          </td><td> real                     </td><td> physical domain left y-boundary position  </td></tr>
!!   <tr><td> ymax     </td><td> 1.          </td><td> real                     </td><td> physical domain right y-boundary position </td></tr>
!!   <tr><td> zmin     </td><td> 0.          </td><td> real                     </td><td> physical domain left z-boundary position  </td></tr>
!!   <tr><td> zmax     </td><td> 1.          </td><td> real                     </td><td> physical domain right z-boundary position </td></tr>
!!   <tr><td> geometry </td><td> "cartesian" </td><td> character(len=cbuff_len) </td><td> \copydoc mpisetup::geometry                   </td></tr>
!! </table>
!! \n \n
!! @b NUMERICAL_SETUP
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>cfl              </td><td>0.7   </td><td>real value between 0.0 and 1.0       </td><td>\copydoc mpisetup::cfl              </td></tr>
!! <tr><td>cfl_max          </td><td>0.9   </td><td>real value between cfl and 1.0       </td><td>\copydoc mpisetup::cfl_max          </td></tr>
!! <tr><td>cflcontrol       </td><td>       </td><td>string                              </td><td>\copydoc mpisetup::cflcontrol       </td></tr>
!! <tr><td>repeat_step      </td><td>.true.</td><td>logical value                        </td><td>\copydoc mpisetup::use_smalld       </td></tr>
!! <tr><td>smallp           </td><td>1.e-10</td><td>real value                           </td><td>\copydoc mpisetup::smallp           </td></tr>
!! <tr><td>smalld           </td><td>1.e-10</td><td>real value                           </td><td>\copydoc mpisetup::smalld           </td></tr>
!! <tr><td>use_smalld       </td><td>.true.</td><td>logical value                        </td><td>\copydoc mpisetup::use_smalld       </td></tr>
!! <tr><td>smallei          </td><td>1.e-10</td><td>real value                           </td><td>\copydoc mpisetup::smallei          </td></tr>
!! <tr><td>smallc           </td><td>1.e-10</td><td>real value                           </td><td>\copydoc mpisetup::smallc           </td></tr>
!! <tr><td>integration_order</td><td>2     </td><td>1 or 2 (or 3 - currently unavailable)</td><td>\copydoc mpisetup::integration_order</td></tr>
!! <tr><td>cfr_smooth       </td><td>0.0   </td><td>real value                           </td><td>\copydoc mpisetup::cfr_smooth       </td></tr>
!! <tr><td>dt_initial       </td><td>-1.   </td><td>positive real value or -1.           </td><td>\copydoc mpisetup::dt_initial       </td></tr>
!! <tr><td>dt_max_grow      </td><td>2.    </td><td>real value > 1.1                     </td><td>\copydoc mpisetup::dt_max_grow      </td></tr>
!! <tr><td>dt_min           </td><td>0.    </td><td>positive real value                  </td><td>\copydoc mpisetup::dt_min           </td></tr>
!! <tr><td>limiter          </td><td>vanleer</td><td>string                              </td><td>\copydoc mpisetup::limiter          </td></tr>
!! <tr><td>relax_time       </td><td>0.0   </td><td>real value                           </td><td>\copydoc mpisetup::relax_time       </td></tr>
!! </table>
!! \n \n
!<
   subroutine init_mpi

      use constants,     only: cwdlen, xdim, ydim, zdim, LO, HI, big_float, dpi, GEO_XYZ, GEO_RPZ, GEO_INVALID, BND_PER, BND_COR, BND_REF, DD_CART, DD_UE
      use mpi,           only: MPI_COMM_WORLD, MPI_COMM_NULL, MPI_INFO_NULL, MPI_PROC_NULL, MPI_STATUS_SIZE, MPI_CHARACTER, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_LOGICAL
      use dataio_pub,    only: die, printinfo, msg, cwd, ansi_white, ansi_black, warn, tmp_log_file
      use dataio_pub,    only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml  ! QA_WARN required for diff_nml

      implicit none

      real :: xmno, ymno, ymxo
      integer :: iproc

      integer, parameter :: hnlen = 32             !< hostname length limit
      character(len=cwdlen) :: cwd_proc
      character(len=hnlen)  :: host_proc
      integer               :: pid_proc

      character(len=cwdlen), allocatable, dimension(:) :: cwd_all
      character(len=hnlen) , allocatable, dimension(:) :: host_all
      integer              , allocatable, dimension(:) :: pid_all

      integer(kind=1)       :: getcwd, hostnm
      integer(kind=4)       :: getpid
      integer :: cwd_status, host_status
      logical :: par_file_exist
      logical :: tmp_log_exist
      integer :: p, px, py
      real :: maxcnt
      integer, dimension(:), allocatable :: pz_slab, py_slab
      integer, dimension(ndims) :: pc

      call MPI_Init( ierr )
      call MPI_Comm_rank(MPI_COMM_WORLD, proc, ierr)
      master = (proc == 0)
      slave  = (proc /= 0)
      comm = MPI_COMM_WORLD
      info = MPI_INFO_NULL
      call MPI_Comm_size(comm, nproc, ierr)
      have_mpi = (nproc > 1)

      !Assume that any tmp_log_file existed before Piernik was started and contains invalid/outydated/... data.
      !Delete it now and keep in mind that any warn, die, printinfo or printio messages issued before this point will be lost as well.
      if (master) then
         inquire(file = tmp_log_file, exist = tmp_log_exist)
         if (tmp_log_exist) then
            open(3, file=tmp_log_file)
            close(3, status="delete")
         endif
      endif
#ifdef VERBOSE
      if (master) call printinfo("[mpisetup:init_mpi]: commencing...")
#endif /* VERBOSE */

      if (allocated(cwd_all) .or. allocated(host_all) .or. allocated(pid_all)) call die("[mpisetup:init_mpi] cwd_all, host_all or pid_all already allocated")

      !> \deprecated BEWARE on slave it is probably enough to allocate only one element or none at all (may depend on MPI implementation)
      allocate(cwd_all(0:nproc))
      allocate(host_all(0:nproc))
      allocate(pid_all(0:nproc))

      pid_proc    = getpid()
      host_status = hostnm(host_proc)
      cwd_status  = getcwd(cwd_proc)

      if (cwd_status /= 0) call die("[mpisetup:init_mpi] problems accessing current working directory.")
#ifdef DEBUG
      write(msg,'(3a,i6,3a)') 'mpisetup: host="',trim(host_proc),'", PID=',pid_proc,' CWD="',trim(cwd_proc),'"'
      call printinfo(msg)
#endif /* DEBUG */

      call MPI_Gather(cwd_proc,  cwdlen, MPI_CHARACTER, cwd_all,  cwdlen, MPI_CHARACTER, 0, comm, ierr)
      call MPI_Gather(host_proc, hnlen,  MPI_CHARACTER, host_all, hnlen,  MPI_CHARACTER, 0, comm, ierr)
      call MPI_Gather(pid_proc,  1,      MPI_INTEGER,   pid_all,  1,      MPI_INTEGER,   0, comm, ierr)

      ! cwd = trim(cwd_proc)  !> \deprecated BEWARE: It's redundant, we get cwd for command line in init_piernik subroutine

      if (master) then
         inquire(file=par_file, exist=par_file_exist)
         if (.not. par_file_exist) call die('[mpisetup:init_mpi] Cannot find "problem.par" in the working directory',0)

         call printinfo("------------------------------------------------------------------------------------------------------", .false.)
         call printinfo("###############     Environment     ###############", .false.)
         call printinfo("", .false.)
         call printinfo("PROCESSES:", .false.)
         do iproc = 0, nproc-1
            write(msg,"(a,i4,a,i6,4a)") " proc=", iproc, ", pid= ",pid_all(iproc), " @",trim(host_all(iproc)), " cwd=",trim(cwd_all(iproc))
            call printinfo(msg, .false.)
         enddo
         call printinfo("", .true.)
         write(msg,"(5a,i5)") 'Start of the',ansi_white,' PIERNIK ',ansi_black,'code. No. of procs = ', nproc
         call printinfo(msg, .true.)
         call printinfo("", .true.)
         call printinfo("###############     Namelist parameters     ###############", .false.)
      endif

      deallocate(host_all)
      deallocate(pid_all)
      deallocate(cwd_all)

      dt_old = -1.

      ! Begin processing of namelist parameters

      psize(:) = [1, 1, 1]

      nxd    = 1
      nyd    = 1
      nzd    = 1
      nb     = 4
      xmin = 0.; xmax = 1.
      ymin = -big_float; ymax = big_float
      zmin = 0.; zmax = 1.
      geometry = "cartesian"

      reorder   = .false.      !< \todo test it!
      allow_uneven = .true.
      allow_noncart = .false.  !< experimental implementation
      use_comm3d = .true.      !< \todo make a big benchmark with and without comm3d

      bnd_xl = 'per'
      bnd_xr = 'per'
      bnd_yl = 'per'
      bnd_yr = 'per'
      bnd_zl = 'per'
      bnd_zr = 'per'

      limiter     = 'vanleer'
      cflcontrol  = 'warn'
      repeat_step = .true.

      cfl         = 0.7
      cfl_max     = 0.9
      cfr_smooth  = 0.0
      smallp      = big_float
      smalld      = big_float
      use_smalld  = .true.
      smallc      = 1.e-10
      smallei     = 1.e-10
      dt_initial  = -1.              !< negative value indicates automatic choice of initial timestep
      dt_max_grow = dt_default_grow  !< for sensitive setups consider setting this as low as 1.1
      dt_min      = tiny(1.)
      dd_unif_quality = 0.9
      dd_rect_quality = 0.9

      integration_order  = 2

      if (master) then
         diff_nml(MPI_BLOCKS)
         diff_nml(NUMERICAL_SETUP)
         diff_nml(DOMAIN)

         ! Sanitize input parameters, if possible
         dom%n_d(:) = max(1, [nxd, nyd, nzd])
         cfl_max = min(max(cfl_max, min(cfl*1.1, cfl+0.05, (1.+cfl)/2.) ), 1.0) ! automatically sanitize cfl_max
         if (integration_order > 2) call die ('[mpisetup:init_mpi]: "ORIG" scheme integration_order must be 1 or 2')

         if (dt_max_grow < 1.01) then
            if (master) then
               write(msg,'(2(a,g10.3))')"[mpisetup:init_mpi] dt_max_grow = ",dt_max_grow," is way too low. Resetting to ",dt_default_grow
               call warn(msg)
            endif
            dt_max_grow = dt_default_grow
         endif

      endif

      if (master) then

         cbuff(1) = bnd_xl
         cbuff(2) = bnd_xr
         cbuff(3) = bnd_yl
         cbuff(4) = bnd_yr
         cbuff(5) = bnd_zl
         cbuff(6) = bnd_zr
         cbuff(7) = limiter
         cbuff(8) = cflcontrol
         cbuff(9) = geometry

         ibuff(xdim:zdim) = psize(:)
         ibuff(4) = integration_order
         ibuff(5:7) = dom%n_d(:)
         ibuff(8) = nb

         rbuff( 1) = smalld
         rbuff( 2) = smallc
         rbuff( 3) = smallp
         rbuff( 4) = smallei
         rbuff( 5) = cfl
         rbuff( 6) = cfr_smooth
         rbuff( 7) = dt_initial
         rbuff( 8) = dt_max_grow
         rbuff( 9) = dt_min
         rbuff(10) = cfl_max
         rbuff(11) = relax_time
         rbuff(12) = xmin
         rbuff(13) = xmax
         rbuff(14) = ymin
         rbuff(15) = ymax
         rbuff(16) = zmin
         rbuff(17) = zmax
         rbuff(18) = dd_unif_quality
         rbuff(19) = dd_rect_quality

         lbuff(1) = use_smalld
         lbuff(2) = reorder
         lbuff(3) = repeat_step
         lbuff(4) = allow_uneven
         lbuff(5) = allow_noncart
         lbuff(6) = use_comm3d

      endif

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
      call MPI_Bcast(ibuff,           buffer_dim, MPI_INTEGER,          0, comm, ierr)
      call MPI_Bcast(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      call MPI_Bcast(lbuff,           buffer_dim, MPI_LOGICAL,          0, comm, ierr)

      if (slave) then

         use_smalld    = lbuff(1)
         reorder       = lbuff(2)
         repeat_step   = lbuff(3)
         allow_uneven  = lbuff(4)
         allow_noncart = lbuff(5)
         use_comm3d    = lbuff(6)

         smalld      = rbuff( 1)
         smallc      = rbuff( 2)
         smallp      = rbuff( 3)
         smallei     = rbuff( 4)
         cfl         = rbuff( 5)
         cfr_smooth  = rbuff( 6)
         dt_initial  = rbuff( 7)
         dt_max_grow = rbuff( 8)
         dt_min      = rbuff( 9)
         cfl_max     = rbuff(10)
         relax_time  = rbuff(11)
         xmin        = rbuff(12)
         xmax        = rbuff(13)
         ymin        = rbuff(14)
         ymax        = rbuff(15)
         zmin        = rbuff(16)
         zmax        = rbuff(17)
         dd_unif_quality = rbuff(18)
         dd_rect_quality = rbuff(19)

         bnd_xl     = cbuff(1)
         bnd_xr     = cbuff(2)
         bnd_yl     = cbuff(3)
         bnd_yr     = cbuff(4)
         bnd_zl     = cbuff(5)
         bnd_zr     = cbuff(6)
         limiter    = cbuff(7)
         cflcontrol = cbuff(8)
         geometry   = cbuff(9)

         psize(:)   = ibuff(xdim:zdim)
         integration_order = ibuff(4)
         dom%n_d(:) = ibuff(5:7)
         nb         = ibuff(8)

      endif

      has_dir(:) = dom%n_d(:) > 1
      eff_dim = count(has_dir(:))

      select case (geometry)
         case ("cartesian", "cart", "xyz", "XYZ")
            geometry_type = GEO_XYZ
         case ("cylindrical", "cyl", "rpz", "RPZ")
            geometry_type = GEO_RPZ
         case default
            geometry_type = GEO_INVALID
      end select

      dom%bnd(:,:) = translate_bnds_to_ints_dom()

      ! sanitize domain
      xmno = xmin
      ymno = ymin
      ymxo = ymax
      select case (geometry_type)
         case (GEO_XYZ)
            if (ymin <= -big_float .or. ymax >= big_float) call warn("[mpisetup:init_mpi] y range not specified. Defaulting to [0..1]")
            if (ymin <= -big_float) ymin = 0.
            if (ymax >= big_float) ymax = 1.
         case (GEO_RPZ)
            if (ymin <= -big_float) ymin = 0.
            if (ymax >= big_float) ymax = dpi
            if (ymax-ymin > dpi) then
               call warn("[mpisetup:init_mpi] Hyperbolic spaces are not implemented. Setting azimuthal span to 2pi.")
               if (abs(ymin) < 1./epsilon(1.)) then
                  ymax = ymin + dpi
               else
                  ymin = ymax - dpi
               endif
               if (abs(ymax-ymin - dpi) > 100*epsilon(1.)) call die("[mpisetup:init_mpi] absolute values for both ymax and ymin too high.") ! magic number
            endif
            if (xmin <= 0.) then
               xmin = 0.
               if (dom%bnd(xdim, LO) /= BND_REF) call warn("[mpisetup:init_mpi] Enforcing dom%bnd(xdim, LO) = 'ref'.")
               dom%bnd(xdim, LO) = BND_REF
            endif
            if (dom%bnd(xdim, HI) == BND_PER) call die("[mpisetup:init_mpi] Periodicity in radial direction is not allowed in cylindrical coordinates")
         case default
            call die("[mpisetup:init_mpi] Invalid geometry type.")
      end select
      if (xmno /= xmin) then
         write(msg,'(2(a,g20.12))')"[mpisetup:init_mpi] Sanitized xmin: ",xmno," -> ",xmin
         call warn(msg)
      endif
      if (ymno /= ymin .and. ymno /= -big_float) then
         write(msg,'(2(a,g20.12))')"[mpisetup:init_mpi] Sanitized ymin: ",ymno," -> ",ymin
         call warn(msg)
      endif
      if (ymxo /= ymax .and. ymxo /= big_float) then
         write(msg,'(2(a,g20.12))')"[mpisetup:init_mpi] Sanitized ymax: ",ymno," -> ",ymax
         call warn(msg)
      endif
      if (xmin > xmax) call die("[[mpisetup:init_mpi] Negative span in X-direction")
      if (ymin > ymax) call die("[[mpisetup:init_mpi] Negative span in Y-direction")
      if (zmin > zmax) call die("[[mpisetup:init_mpi] Negative span in Z-direction")

      ! set up the global domain
      dom%nb = nb

      dom%xmin = xmin
      dom%ymin = ymin
      dom%zmin = zmin
      dom%xmax = xmax
      dom%ymax = ymax
      dom%zmax = zmax

      call dom%set_derived ! finish up with the rest of domain_container members

      if (master .and. have_mpi) then
         if (allow_uneven) call warn("[mpisetup:init_mpi] Uneven domain decomposition is experimental.")
         if (allow_noncart) call warn("[mpisetup:init_mpi] Non-cartesian domain decomposition is highly experimental.")
      endif
      is_uneven = .false.
      is_mpi_noncart = .false.

      if (allocated(procmask)) call die("[mpisetup:init_mpi] procmask already allocated")
      allocate(procmask(0:nproc-1))

      where (.not. has_dir(:)) psize(:) = 1
      call divide_domain

      if (is_mpi_noncart) is_uneven = .true.

      ! most of the code below is incompatible with noncartesian domain decomposition and shoould be called from divide_domain_uniform and divide_domain_rectlinear

      if (any(dom%bnd(xdim:ydim, LO:HI) == BND_COR)) then
         if (is_mpi_noncart) call die("[mpisetup:init_mpi] Corner BC with noncartesian domain division not implemented")
         if (psize(xdim) /= psize(ydim) .or. dom%n_d(xdim) /= dom%n_d(ydim)) then
            write(msg, '(a,4(i4,a))')"[mpisetup:init_mpi] Corner BC require psize(xdim) equal to psize(ydim) and nxd equal to nyd. Detected: [",psize(xdim),",",psize(ydim),&
                 &                   "] and [",dom%n_d(xdim),",",dom%n_d(ydim),"]"
            call die(msg)
         endif
      endif
      if (any(dom%bnd(zdim, LO:HI) == BND_COR)) call die("[mpisetup:init_mpi] Corner BC not allowed for z-direction")

      if (allocated(dom%se)) call die("[mpisetup:init_mpi] dom%se already allocated")
      allocate(dom%se(0:nproc-1, xdim:zdim, LO:HI))
      dom%se(:, :, :) = 0

      procn(:,:) = MPI_PROC_NULL
      comm3d = MPI_COMM_NULL

      select case (dom%pdiv_type)
         case (DD_CART)

            if (is_mpi_noncart) call die("[mpisetup:init_mpi] inconsistent decomposition specs")

            if (use_comm3d) then
               if (is_mpi_noncart) call die("[mpisetup:init_mpi] MPI_Cart_create cannot be used for non-rectilinear or AMR domains")
               if (master) call printinfo("[mpisetup:init_mpi] Cartesian decomposition with comm3d")

               call MPI_Cart_create(comm, ndims, psize, dom%periodic, reorder, comm3d, ierr)
               call MPI_Cart_coords(comm3d, proc, ndims, pcoords, ierr)

               do p = 0, nproc-1
                  call MPI_Cart_coords(comm3d, p, ndims, pc, ierr)
                  where (has_dir(:))
                     dom%se(p, :, LO) = (dom%n_d(:) *  pc(:) ) / psize(:)     ! offset of low boundaries of the local domain (0 at low external boundaries)
                     dom%se(p, :, HI) = (dom%n_d(:) * (pc(:)+1))/psize(:) - 1 ! offset of high boundaties of the local domain (n_d(:) - 1 at right external boundaries)
                  endwhere
               enddo

               ! Compute neighbors
               do p = xdim, zdim
                  call MPI_Cart_shift(comm3d, p-xdim, 1, procn(p, LO), procn(p, HI), ierr)
               enddo

               if (any(dom%bnd(xdim:ydim, LO) == BND_COR)) then
                  if (pcoords(xdim) == 0 .and. pcoords(ydim) > 0) then
                     pc = (/pcoords(ydim),pcoords(xdim),pcoords(zdim)/)
                     call MPI_Cart_rank(comm3d,pc,procxyl,ierr)
                  else
                     procxyl = MPI_PROC_NULL
                  endif
                  if (pcoords(ydim) == 0 .and. pcoords(xdim) > 0 ) then
                     pc = (/pcoords(ydim),pcoords(xdim),pcoords(zdim)/)
                     call MPI_Cart_rank(comm3d,pc,procyxl,ierr)
                  else
                     procyxl = MPI_PROC_NULL
                  endif
               endif

               if (any(dom%bnd(xdim:ydim, HI) == BND_COR)) call die("[mpisetup:init_mpi] Corner boundary on the right side not implemented anywhere")

            else
               if (master) call printinfo("[mpisetup:init_mpi] Cartesian decomposition without comm3d")
               do p = 0, nproc-1
                  if (product(psize(:)) /= nproc) call die("[mpisetup:init_mpi] product(psize(:)) /= nproc")
                  pc(xdim) = mod(p, psize(xdim))
                  pc(ydim) = mod(p/psize(xdim), psize(ydim))
                  pc(zdim) = p/product(psize(xdim:ydim))
                  where (has_dir(:))
                     dom%se(p, :, LO) = (dom%n_d(:) *  pc(:) ) / psize(:)     ! offset of low boundaries of the local domain (0 at low external boundaries)
                     dom%se(p, :, HI) = (dom%n_d(:) * (pc(:)+1))/psize(:) - 1 ! offset of high boundaties of the local domain (n_d(:) - 1 at right external boundaries)
                  endwhere
               enddo

               if (any(dom%bnd(xdim:ydim, LO:HI) == BND_COR)) call die("[mpisetup:init_mpi] BND_COR not implemented for use_comm3d == .false.")
            endif
         case (DD_UE)
            if (master) call printinfo("[mpisetup:init_mpi] Non-cartesian decomposition (no comm3d possible)")
            allocate(pz_slab(psize(zdim) + 1))
            do p = 0, psize(zdim)
               pz_slab(p+1) = (p * nproc) / psize(zdim)
            enddo
            do p = 1, psize(zdim)
               dom%se(pz_slab(p):pz_slab(p+1)-1, zdim, LO) = nint((dom%n_d(zdim) *  pz_slab(p)   ) / real(nproc))
               dom%se(pz_slab(p):pz_slab(p+1)-1, zdim, HI) = nint((dom%n_d(zdim) *  pz_slab(p+1) ) / real(nproc)) - 1

               allocate(py_slab(psize(ydim) + 1))
               do py = 0, psize(ydim)
                  py_slab(py+1) = (py * (pz_slab(p+1)-pz_slab(p))) / psize(ydim) !> \todo try to sort lengths
               enddo
               do py = 1, psize(ydim)
                  dom%se(pz_slab(p)+py_slab(py):pz_slab(p)+py_slab(py+1)-1, ydim, LO) = nint((dom%n_d(ydim) *  py_slab(py)   ) / real(pz_slab(p+1)-pz_slab(p)))
                  dom%se(pz_slab(p)+py_slab(py):pz_slab(p)+py_slab(py+1)-1, ydim, HI) = nint((dom%n_d(ydim) *  py_slab(py+1) ) / real(pz_slab(p+1)-pz_slab(p))) - 1
                  do px = 0, py_slab(py+1)-py_slab(py) - 1
                     dom%se(pz_slab(p)+py_slab(py)+px, xdim, LO) = (dom%n_d(xdim) *  px    ) / (py_slab(py+1)-py_slab(py))
                     dom%se(pz_slab(p)+py_slab(py)+px, xdim, HI) = (dom%n_d(xdim) * (px+1) ) / (py_slab(py+1)-py_slab(py)) - 1
                  enddo
               enddo
               if (allocated(py_slab)) deallocate(py_slab)
            enddo
            if (allocated(pz_slab)) deallocate(pz_slab)
         case default
            call die("[mpisetup:init_mpi] unknown strategy for generating domain division")
      end select

      if (allocated(req) .or. allocated(status)) call die("[mpisetup:init_mpi] req or status already allocated")
      if (comm3d == MPI_COMM_NULL) then
         allocate(req(nreq))
         allocate(status(MPI_STATUS_SIZE, nreq))
      else
         allocate(req(4*nproc)) ! 4 = count([i_bnd, o_bnd]) * two sides
         allocate(status(MPI_STATUS_SIZE, size(req)))
      endif

      if (any(dom%bnd(:, :) == BND_COR) .and. comm3d == MPI_COMM_NULL) call die("[mpisetup:init_mpi] Corner BC not implemented without comm3d")

      if ((dom%periodic(xdim) .and. dom%se(proc, xdim, HI) /= dom%n_d(xdim) - 1) .or. dom%se(proc, xdim, LO) /= 0)                 bnd_xl = 'mpi'
      if ((dom%periodic(xdim) .and. dom%se(proc, xdim, LO) /= 0)                 .or. dom%se(proc, xdim, HI) /= dom%n_d(xdim) - 1) bnd_xr = 'mpi'
      if ((dom%periodic(ydim) .and. dom%se(proc, ydim, HI) /= dom%n_d(ydim) - 1) .or. dom%se(proc, ydim, LO) /= 0)                 bnd_yl = 'mpi'
      if ((dom%periodic(ydim) .and. dom%se(proc, ydim, LO) /= 0)                 .or. dom%se(proc, ydim, HI) /= dom%n_d(ydim) - 1) bnd_yr = 'mpi'
      if ((dom%periodic(zdim) .and. dom%se(proc, zdim, HI) /= dom%n_d(zdim) - 1) .or. dom%se(proc, zdim, LO) /= 0)                 bnd_zl = 'mpi'
      if ((dom%periodic(zdim) .and. dom%se(proc, zdim, LO) /= 0)                 .or. dom%se(proc, zdim, HI) /= dom%n_d(zdim) - 1) bnd_zr = 'mpi'

      ! For shear boundaries and some domain decompositions it is possible that a boundary can be mixed 'per' with 'mpi'
#ifdef SHEAR_BND
      if (comm3d == MPI_COMM_NULL) call die("[mpisetup:init_mpi] SHEAR_BND not implemented without comm3d")

      if (psize(ydim) > 1) call die("[mpisetup:initmpi] Shear-pediodic boundary conditions do not permit psize(ydim) > 1")
      ! This will be possible when is_mpi_noncart will be fully implemented

#ifndef FFTW
      if (pcoords(xdim) == 0) then
         bnd_xl = 'she'
      else
         bnd_xl = 'mpi'
      endif

      if (pcoords(xdim) == psize(xdim)-1) then
         bnd_xr = 'she'
      else
         bnd_xr = 'mpi'
      endif
#endif /* !FFTW */
#endif /* SHEAR_BND */

#ifdef DEBUG
      if (comm3d /= MPI_COMM_NULL) then
         do p = xdim, zdim
            write(msg,*) 'dir',p,': ',procn(p, LO), proc, procn(p, HI)
            call printinfo(msg)
         enddo
      endif
#endif /* DEBUG */

!#ifdef VERBOSE
      if (master) then
         maxcnt = 0
         do p = 0, nproc - 1
            write(msg,'(a,i4,a,2(3i6,a),i8,a)')"[mpisetup:init_mpi] segment @",p," : [",dom%se(p, :, LO),"] : [",dom%se(p, :, HI),"] #",product(dom%se(p, :, HI)-dom%se(p, :, LO)+1)," cells"
            call printinfo(msg)
            maxcnt = max(maxcnt, product(real(dom%se(p, :, HI)-dom%se(p, :, LO)+1)))
         enddo
         write(msg,'(a,f8.5)')"[mpisetup:init_mpi] Load balance: ",product(real(dom%n_d(:)))/(nproc*maxcnt) !> \todo add calculation of total internal boundary surface in cells
         call printinfo(msg)
      endif
!      call identify_neighbours(dom%se(proc, :, :), dom)
!#endif /* VERBOSE */

   end subroutine init_mpi

!-----------------------------------------------------------------------------
!
! \todo: prepare some useful lists here, so routines grid:: arr3d_boundaries, grid::grid_mpi_boundaries_prep, fluidboundaries, magboundaries, multigridmpifuncs
! can use them and no longer rely on pcoords(:), psize(:), proc[xyz][lr], proc{xy,yx}l, ...
! note that corner boundaries will require transposing on the acceptor side
!
   subroutine identify_neighbours(area, domain)

      use constants,  only: xdim, zdim, LO, HI
#ifdef VERBOSE
      use dataio_pub, only: printinfo, msg
#endif /* VERBOSE */

      implicit none

      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: area
      type(domain_container), intent(in) :: domain

      integer :: p, pp
      logical, dimension(xdim:zdim, LO:HI) :: neigh
      logical :: sharing, cor, fac
      integer(kind=8), dimension(xdim:zdim) :: per

      pp = 0 !nproc/3

      do p = 0, nproc-1
         per(:) = 0
         where (dom%periodic(:)) per(:) = dom%n_d(:)
         call is_neigh(area(:,:), domain%se(p, :, :), neigh(:,:), sharing, cor, fac, per)
#ifdef VERBOSE
         if (proc == pp) then
            write(msg,'(a,i4,a,6i4,a,l2,a,6l2,2(a,l2))')"m:in <@>",p," [",domain%se(p,:,:),"] * ",sharing, " > [",neigh(:,:),"] ?f",fac," ?c ",cor
            call printinfo(msg)
         endif
#endif /* VERBOSE */
      enddo
#ifdef VERBOSE
      if (proc == pp) then
         write(msg,'(a,i4,a,6i4)')"m:in <@>",pp," <",area(:,:)
         call printinfo(msg)
      endif
#endif /* VERBOSE */

   end subroutine identify_neighbours

!!-----------------------------------------------------------------------------
!!
!> \brief is_neigh_per checks if two given blocks placed within a periodic domain are overlapping and determines whether the second block can provide face or corner guardcells
!!
!! \details to handle shearing box which is divided in y-direction atthe edges, one has to provide another subroutine (is_neigh_per_shear) and add it to interface is_neigh
!!
!! \deprecated corner and face aren't used
!<
   subroutine is_neigh_per(this, other, neigh, share, corner, face, periods)

      use constants,  only: xdim, ydim, zdim, LO, HI

      implicit none

      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: this, other !< two boxes
      logical, dimension(xdim:zdim, LO:HI), intent(out)        :: neigh       !< on which sides the other can offer useful guardcells?
      logical, intent(out)                                     :: share       !< is there overlap between this and the other?
      logical, intent(out)                                     :: corner      !< is there overlap in edge or corner guardcells? (anything but face guardcells)
      logical, intent(out)                                     :: face        !< are there any directly (face) neighbouring grid cells?
      integer(kind=8), dimension(xdim:zdim), intent(in)                :: periods     !< where >0 then the direction is periodic with the given number of cells

      integer :: i, j, k
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: oth
      logical, dimension(xdim:zdim, LO:HI) :: nei
      logical :: sha, cor, fac

      neigh(:,:) = .false. ; share = .false. ; corner = .false. ; face = .false.
      do i = -1, 1
         if ((has_dir(xdim) .or. periods(xdim)>0) .or. i==0) then
            do j = -1, 1
               if ((has_dir(ydim) .or. periods(ydim)>0) .or. j==0) then
                  do k = -1, 1
                     if ((has_dir(zdim) .or. periods(zdim)>0) .or. k==0) then
                        oth(:,:) = other(:,:) + reshape([i*periods(xdim), j*periods(ydim), k*periods(zdim), i*periods(xdim), j*periods(ydim), k*periods(zdim)], [3,2])
                        call is_neigh_simple(this, oth, nei, sha, cor, fac)
                        neigh(:,:) = neigh(:,:) .or. nei(:,:)
                        share = share .or. sha
                        face = face .or. fac
                        corner = (corner .or. cor) .and. .not. face
                     endif
                  enddo
               endif
            enddo
         endif
      enddo

   end subroutine is_neigh_per

!!-----------------------------------------------------------------------------
!!
!> \brief is_neigh_simple checks if two given blocks placed within a nonperiodic domain are overlapping and determines whether the second block can provide face or corner guardcells
!!

   subroutine is_neigh_simple(this, other, neigh, share, corner, face)

      use constants,  only: xdim, zdim, LO, HI

      implicit none

      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: this, other !< two boxes
      logical, dimension(xdim:zdim, LO:HI), intent(out)        :: neigh       !< on which sides the other can offer useful guardcells?
      logical, intent(out)                                     :: share       !< is there overlap between this and the other?
      logical, intent(out)                                     :: corner      !< is there overlap in edge or corner guardcells? (anything but face guardcells)
      logical, intent(out)                                     :: face        !< are there any directly (face) neighbouring grid cells?

      logical, dimension(xdim:zdim, LO:HI) :: d_neigh
      logical, dimension(xdim:zdim)        :: d_share
      logical                              :: aux
      integer                              :: d

      ! periodicity ignored

      do d = xdim, zdim
         if (has_dir(d)) then
            d_neigh(d, LO) = (other(d, HI) + 1 >= this(d, LO)) .and. (other(d, LO) <  this(d, LO))
            d_neigh(d, HI) = (other(d, LO) - 1 <= this(d, HI)) .and. (other(d, HI) >  this(d, HI))
            d_share(d)     = (other(d, LO)     <= this(d, HI)) .and. (other(d, HI) >= this(d, LO))
         else
            d_neigh(d, :) = .false.
            d_share(d)    = .true.
         endif
      enddo

      aux = all(d_neigh(:, LO) .or. d_neigh(:, HI) .or. d_share(:))

      neigh(:,:) = d_neigh(:, :) .and. aux
      share = all(d_share(:))
      corner = count(d_share(:))<=1 .and. count((d_neigh(:, LO) .or. d_neigh(:, HI)) .and. aux) > 1
      face = count(neigh(:, :)) > 0 .and. .not. corner

   end subroutine is_neigh_simple

!-----------------------------------------------------------------------------

   subroutine cleanup_mpi

      use dataio_pub, only: printinfo
      use mpi,        only: MPI_COMM_NULL

      implicit none

      if (allocated(dom%se)) deallocate(dom%se)
      if (allocated(primes)) deallocate(primes)
      if (allocated(procmask)) deallocate(procmask)
      if (allocated(req)) deallocate(req)
      if (allocated(status)) deallocate(status)

      if (comm3d /= MPI_COMM_NULL) call MPI_Comm_free(comm3d, ierr)

      if (master) call printinfo("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++", .false.)
      call MPI_Barrier(comm,ierr)
      if (nproc > 1) call sleep(1) ! Prevent random SIGSEGVs in openmpi's MPI_Finalize
      call MPI_Finalize(ierr)

   end subroutine cleanup_mpi

!-----------------------------------------------------------------------------

   subroutine mpifind(var, what)

      use types,         only: value
      use constants,     only: MINL, MAXL
      use dataio_pub,    only: die
      use mpi,           only: MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_MAXLOC, MPI_IN_PLACE

      implicit none

      type(value), intent(inout) :: var
      integer, intent(in)       :: what

      real, dimension(2)  :: v_red
      integer, dimension(MINL:MAXL), parameter :: op = [ MPI_MINLOC, MPI_MAXLOC ]

      v_red(:) = [ var%val, real(proc) ]

      if (any([MINL, MAXL] == what)) then
         call MPI_Allreduce(MPI_IN_PLACE, v_red, 1, MPI_2DOUBLE_PRECISION, op(what), comm, ierr)
      else
         call die("[mpisetup:mpifind] invalid extremum type")
      endif

      var%val = v_red(1)
      var%proc = int(v_red(2))

   end subroutine mpifind

!-----------------------------------------------------------------------------
!
!> \brief This routine computes optimal allowed domain division
!
!> \todo make this a member of types::domain_container
!
   subroutine divide_domain

      use constants,  only: DD_CART, DD_UE
      use dataio_pub, only: die, warn, printinfo, msg

      implicit none

      real :: quality

      dom%pdiv_type = DD_CART
      dom%pdiv(:) = psize(:)

      if (product(psize) == nproc) then
         if (all(mod(dom%n_d(:), psize(:)) == 0)) then
            if (master .and. have_mpi) then
               write(msg,'(a,3i4,a,3i6,a)')"[mpisetup:divide_domain] Domain divided to [",psize(:)," ] pieces, each of [",dom%n_d(:)/psize(:)," ] cells."
               call printinfo(msg)
            endif
            return
         else
            write(msg,'(a,3i6,a,3i4,a)')"[mpisetup:divide_domain] Cannot divide domain with [",dom%n_d(:)," ] cells to [",psize(:)," ] piecess. "
            if (master) call warn(msg)
            psize(:) = 1
         endif
      endif

      call Eratosthenes_sieve(primes, nproc) ! it is possible to use primes only to sqrt(nproc), but it is easier to have the full table. Cheap for any reasonable nproc.

      ! this is the minimal total area of internal boundaries (periodic case), achievable for some perfect domain divisions
      ideal_bsize = eff_dim * (nproc * product(real(dom%n_d(:)))**(eff_dim-1))**(1./eff_dim)

      call divide_domain_uniform
      if (product(psize) == nproc) then
         quality = ideal_bsize / sum(psize(:)/real(dom%n_d(:)) * product(real(dom%n_d(:))), MASK = dom%n_d(:) > 1)
         if (quality >= dd_unif_quality .or. .not. (allow_uneven .or. allow_noncart)) then
            dom%pdiv(:) = psize(:)
            return
         endif
         write(msg,'(2(a,f6.3),a)')"[mpisetup:divide_domain] Quality of uniform division = ",quality," is below threshold ",dd_unif_quality, ", trying harder ..."
         if (master) call warn(msg)
      endif

      if (allow_uneven) then
         call divide_domain_rectlinear
         quality = 1 !< \todo make an estimate
         if (product(psize) == nproc) then
            if (quality > dd_rect_quality .or. .not. allow_noncart) then
               dom%pdiv(:) = psize(:)
               return
            endif
         endif
      else
         if (master) call warn("[mpisetup:divide_domain] Did not try uneven domain division")
      endif

      if (allow_noncart) then
         dom%pdiv_type = DD_UE
         psize(:) = dom%pdiv(:)
         call divide_domain_slices
         ! if good_enough then return
         dom%pdiv(:) = psize(:)
         return
      else
         if (master) call warn("[mpisetup:divide_domain] Did not try non-cartesian domain division")
      endif

      write(msg,'(a,3i6,a,i4,a)') "[mpisetup:divide_domain] Cannot divide domain with [",dom%n_d(:)," ] cells to ",nproc," piecess. "
      call die(msg)

   end subroutine divide_domain

!-----------------------------------------------------------------------------
!>
!! \brief This routine tries to divide the computational domain into local domains.
!! \details The goal is to minimize the ratio of longest to shortest edge to minimize the amount of inter-process communication.
!! If the benchmarks show that some direction should be partitioned in more pieces than other directions, implement appropriate weighting in j1, j2 and j3 calculations.
!!
!! For some weird domains and PE counts this routine may find tiling that does not satisfy multigrid restrictions even if there exists an acceptable tiling.
!! In such case the user must divide domain manually by providing p[xyz]size parameters through problem.par.
!!
!! \attention Must be called by all procs to avoid communication and ensure that every proc has proper psize
!<
   subroutine divide_domain_uniform

      use constants,     only: xdim, zdim
      use dataio_pub,    only: warn, printinfo, msg

      implicit none

      integer :: j1, j2, j3, jj, n, p
      integer, dimension(ndims) :: ldom, tmp

      ldom(xdim:zdim) = dom%n_d(zdim:xdim:-1) ! Maxloc returns first occurrence of max, reversing direction order (to ZYX) gives better cache utilization.
      n = nproc
      psize(:) = 1
      if (n == 1) return

      do p = size(primes), 1, -1 ! start from largest defined primes, continue down to 2
         do while (mod(n, primes(p))==0)

            jj = 0
            j1 = sum(maxloc(ldom), 1) ! First try the longest edge; note the trick to make a scalar from 1-element vector without assignment to another variable
            if (mod(ldom(j1), primes(p))==0) then
               jj = j1
            else
               j2 = 1 + mod(j1 + 0, ndims)
               j3 = 1 + mod(j1 + ndims -2, ndims)
               if (ldom(j2) > ldom(j3)) then
                  j2 = 1 + mod(j1 + ndims -2, ndims)
                  j3 = 1 + mod(j1 + 0, ndims)
               endif
               if (mod(ldom(j2), primes(p))==0) jj = j2 ! middle edge ...
               if (jj == 0 .and. mod(ldom(j3), primes(p))==0) jj = j3 ! try the shortest edge on last resort
            endif

            if (jj == 0) then
               if (master) call warn("[mpisetup:divide_domain_uniform]: Can't find divisible edge")
               psize(:) = 1
               return
            else
               psize(jj) = psize(jj) * primes(p)
               n         = n         / primes(p)
               ldom(jj)  = ldom(jj)  / primes(p)
            endif

         enddo
      enddo

      if (n /= 1) then
         if (master) call warn("[mpisetup:divide_domain_uniform]: I am not that intelligent") ! nproc has too big prime factors
         psize(:) = 1
         return
      endif

      tmp(xdim:zdim) = psize(zdim:xdim:-1) ! directions were reverted at ldom assignment
      psize(:) = tmp(:)

      if (master) then
         write(msg,'(a,3i4,a,3i6,a)')"[mpisetup:divide_domain_uniform] Domain divided to [",psize(:)," ] pieces, each of [",ldom(zdim:xdim:-1)," ] cells."
         call printinfo(msg)
      endif

   end subroutine divide_domain_uniform

!-----------------------------------------------------------------------------
!>
!! \brief Divide the computational domain into local domains. Allow their size to change by +/- 1 depending on CPU rank (this will introduce some load imbalance)
!! if it is not possible to divide an edge evenly. Try to minimize the imbalance and total internal boundaries size.
!<
   subroutine divide_domain_rectlinear

      use constants,     only: xdim, ydim
      use dataio_pub,    only: printinfo, msg

      implicit none

      real, parameter :: b_load_fac = 0.25 ! estimated increase of execution time after doubling the total size of internal boundaries.
      ! \todo estimate this factor for massively parallel runs and for Intel processors

      integer, allocatable, dimension(:) :: ppow
      integer, allocatable, dimension(:,:) :: fac
      integer, dimension(ndims) :: ldom
      integer :: p, i, j, k, n, nf, ii, bsize
      real :: load_balance, best, quality

      psize(:) = 1
      if (nproc == 1) return

      allocate(ppow(size(primes)))

      p = nproc
      do i = 1, size(primes)
         ppow(i) = 0
         do while (mod(p, primes(i)) == 0)
            ppow(i) = ppow(i) + 1
            p = p / primes(i)
         enddo
      enddo

      nf = count(ppow(:) > 0)
      allocate(fac(nf,3))
      j = 1
      do i = 1, size(primes)
         if (ppow(i)>0) then
            fac(j,:) = [ primes(i), ppow(i), (ppow(i)+1)*(ppow(i)+2)/2 ] ! prime, its power and number of different decompositions in three dimensions for this prime
            j = j + 1
         endif
      enddo
      deallocate(ppow)

      best = 0.
      ii = 0
      do while (all(fac(:,3) > 0))
         ldom(:) = 1
         do n = 1, nf ! find an unique decomposition of fac(n,2) into [i,j,k], all([i,j,k] >= 0) && i+j+k = fac(n,2). The decompositions are enumerated with fac(n,3).
            i = int(sqrt(1./4.+2.*(fac(n,3)-1)) - 1./2.) ! i and k enumerate a point in a triangle: (i>=0 && k>=0 && i+k<=fac(n,2))
            k = fac(n,3) - 1 - i*(i+1)/2
            i = fac(n,2) - i
            j = fac(n,2) - (i + k)
            ldom(:) = ldom(:) * fac(n,1)**[i, j, k]
         enddo

         ii = ii + 1
         bsize = int(sum(ldom(:)/dble(dom%n_d(:)) * product(dom%n_d(:)), MASK = dom%n_d(:) > 1)) !ldom(1)*dom%n_d(2)*dom%n_d(3) + ldom(2)*dom%n_d(1)*dom%n_d(3) + ldom(3)*dom%n_d(1)*dom%n_d(2)
         load_balance = product(real(dom%n_d(:))) / ( real(nproc) * product( int((dom%n_d(:)-1)/ldom(:)) + 1 ) )

         quality = load_balance/ (1 + b_load_fac*(bsize/ideal_bsize - 1.))
         ! \todo add a factor that estimates lower cost when x-direction is not chopped too much
         quality = quality * (1. - (0.001 * ldom(xdim) + 0.0001 * ldom(ydim))/nproc) ! \deprecated estimate these magic numbers

#ifdef DEBUG
         if (master) then
            write(msg,'(a,i3,a,3i4,a,i10,2(a,f10.7))')"m:ddr ",ii," psize= [",ldom(:)," ], bcells= ", bsize, ", balance = ", load_balance, ", est_quality = ", quality
            call printinfo(msg)
         endif
#endif /* DEBUG */
         if (quality > best) then
            best = quality
            psize(:) = ldom(:)
         endif
         do j = 1, nf ! search for next unique combination
            if (fac(j,3) > 1) then
               fac(j,3) = fac(j,3) - 1
               exit
            else
               if (j<nf) then
                  fac(j,3) = (fac(j,2)+1)*(fac(j,2)+2)/2
               else
                  fac(:,3) = 0 ! no more combinations to try
               endif
            endif
         enddo
      enddo

      deallocate(fac)

      is_uneven = any(mod(dom%n_d(:), psize(:)) /= 0)

      if (master) then
#ifdef DEBUG
         write(msg,'(a,3f10.2,a,i10)')"m:ddr id psize = [",(nproc/product(dble(dom%n_d(:))))**(1./eff_dim)*dom%n_d(:),"], bcells= ", int(ideal_bsize)
         call printinfo(msg)
#endif /* DEBUG */
         write(msg,'(a,3i4,a)')      "[mpisetup:divide_domain_rectlinear] Domain divided to [",psize(:)," ] pieces"
         call printinfo(msg)
         if (is_uneven) then
            write(msg,'(2(a,3i5),a)')"                                    Sizes are from [", int(dom%n_d(:)/psize(:))," ] to [",int((dom%n_d(:)-1)/psize(:))+1," ] cells."
            call printinfo(msg)
            write(msg,'(a,f8.5)')    "                                    Load balance is ",product(dom%n_d(:)) / ( dble(nproc) * product( int((dom%n_d(:)-1)/psize(:)) + 1 ) )
         else
            write(msg,'(a,3i5,a)')   "                                    Size is [", int(dom%n_d(:)/psize(:))," ] cells."
         endif
         call printinfo(msg)
      endif

   end subroutine divide_domain_rectlinear

!>
!! \brief Divide the computational domain into local domains. Allow their size to depend significantly on CPU rank and allow for more than one neighbour on a single boundary.
!! Try to minimize the imbalance and total internal boundaries size.
!<
   subroutine divide_domain_slices

      use constants,  only: xdim, ydim, zdim
      use dataio_pub, only: msg, printinfo, warn

      implicit none

      real, parameter :: minfac = 1.3 ! prevent domain division to halves if cell count in a given direction is too low. (not verified for optimality)
      real :: optc

      !> \todo Try to make an intelligent guess for slicing, then go down to the local minimum and explore neighbourhood. Exploring all possibilities is an O(nproc)**2 task
      ! The best solution is probably near (nproc/product(dble(dom%n_d(:))))**(1./eff_dim)*dom%n_d(:)

      is_mpi_noncart = .true.
      is_uneven = .true.

      if (all(psize(ydim:zdim) == 1)) then
         if (has_dir(zdim)) then
            optc = (product(dom%n_d(:))/real(nproc)) ** (1./eff_dim) ! number of cells for ideal cubes
            if (dom%n_d(zdim) > minfac*optc) psize(zdim) = ceiling(dom%n_d(zdim)/optc)
         endif
         if (has_dir(ydim)) then
            optc = (product(dom%n_d(xdim:ydim))*psize(zdim)/real(nproc)) ** (1./count(has_dir(xdim:ydim)))
            if (dom%n_d(ydim) > minfac*optc) psize(ydim) = ceiling(dom%n_d(ydim)/optc)
         endif
      endif
      if (has_dir(xdim)) psize(xdim) = (nproc - 1)/(psize(ydim)*psize(zdim)) + 1 !sometimes it might be less by 1

      where (.not. has_dir(:)) psize(:) = 1 ! just in case
      do while (product(psize(:)) < nproc)
         write(msg,'(a,3i4,a)') "[mpisetup:divide_domain_slices] imperfect noncartesian division to [",psize(:)," ] pieces"
         if (master) call warn(msg)
         if (has_dir(xdim)) then
            psize(xdim) = psize(xdim) + 1
         else if (has_dir(ydim)) then
            psize(ydim) = psize(ydim) + 1
         else
            psize(zdim) = psize(zdim) + 1
         endif
      enddo
      write(msg,'(a,3i4,a)') "[mpisetup:divide_domain_slices] performed noncartesian division to [",psize(:)," ] pieces"
      if (master) call printinfo(msg)

   end subroutine divide_domain_slices

!-----------------------------------------------------------------------------

   logical function grace_period_passed()
      implicit none
      grace_period_passed = (t >= relax_time)
   end function grace_period_passed

!-----------------------------------------------------------------------------

   subroutine Eratosthenes_sieve(tab,n)

      use dataio_pub, only: die
#ifdef DEBUG
      use dataio_pub, only: msg, printinfo
#endif /* DEBUG */

      implicit none

      integer, intent(inout), allocatable, dimension(:) :: tab
      integer, intent(in)                               :: n
      integer, dimension(n)                             :: numb

      integer :: i, no_primes

      if (allocated(tab)) call die("[mpisetup:Eratosthenes_sieve] tab already allocated")

      numb = [0, (i, i = 2, n)]

      do i = 2, n
         if (numb(i) /= 0) numb( 2*i : n : i ) = 0
      enddo

      no_primes = count(numb /= 0)
      allocate(tab(no_primes))
      tab       = pack(numb, numb /= 0)
#ifdef DEBUG
      if (master) then
         write(msg,'(2(A,I5))') "There are ", no_primes, " prime numbers less than", n
         call printinfo(msg)
         print "(5i7)", tab
      endif
#endif /* DEBUG */
   end subroutine Eratosthenes_sieve

!-----------------------------------------------------------------------------
!
!  public wrapper
!
   function translate_bnds_to_ints_dom() result(tab)

      use constants, only: ndims, LO, HI

      implicit none

      integer, dimension(ndims, LO:HI) :: tab

      tab(:,:) = translate_bnds_to_ints([bnd_xl, bnd_xr, bnd_yl, bnd_yr, bnd_zl, bnd_zr])

   end function translate_bnds_to_ints_dom

!-----------------------------------------------------------------------------

   function translate_bnds_to_ints(bnds) result(tab)

      use constants, only: xdim, zdim, ndims, LO, HI, &
         &                 BND_MPI, BND_PER, BND_REF, BND_OUT, BND_OUTD, BND_OUTH, BND_COR, BND_SHE, BND_USER, BND_INF, BND_INVALID

      implicit none

      character(len=*), dimension(HI*ndims) :: bnds ! expect exactly 6 descriptions
      integer, dimension(ndims, LO:HI) :: tab

      integer :: d, lh

      do d = xdim, zdim
         do lh = LO, HI
            select case (bnds(HI*(d-xdim)+lh))
               case ('per', 'periodic')
                  tab(d, lh) = BND_PER
               case ('ref', 'refl', 'reflecting')
                  tab(d, lh) = BND_REF
               case ('out', 'free')
                  tab(d, lh) = BND_OUT
               case ('outd', 'diode')
                  tab(d, lh) = BND_OUTD
               case ('outh')
                  tab(d, lh) = BND_OUTH
               case ('she', 'shear', 'shearing')
                  tab(d, lh) = BND_SHE
               case ('cor', 'corner')
                  tab(d, lh) = BND_COR
               case ('mpi')
                  tab(d, lh) = BND_MPI
               case ('user')
                  tab(d, lh) = BND_USER
               case ('inf') ! what is this?
                  tab(d, lh) = BND_INF
               case default
                  tab(d, lh) = BND_INVALID
            end select
         enddo
      enddo
   end function translate_bnds_to_ints

end module mpisetup
