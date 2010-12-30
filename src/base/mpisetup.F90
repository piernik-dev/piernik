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
!   use mpi, only: MPI_STATUS_SIZE
   use dataio_pub, only: cbuff_len

   implicit none

   integer, parameter :: MPI_STATUS_SIZE = 5  ! taken from mpi to silence warnings

   private
   public :: &
        & ARR_XY_LEFT_BND, ARR_XY_LEFT_DOM, ARR_XY_RIGHT_BND, ARR_XY_RIGHT_DOM, ARR_XZ_LEFT_BND, ARR_XZ_LEFT_DOM, ARR_XZ_RIGHT_BND, ARR_XZ_RIGHT_DOM, &
        & ARR_YZ_LEFT_BND, ARR_YZ_LEFT_DOM, ARR_YZ_RIGHT_BND, ARR_YZ_RIGHT_DOM, MAG_XY_LEFT_BND, MAG_XY_LEFT_DOM, MAG_XY_RIGHT_BND, MAG_XY_RIGHT_DOM, &
        & MAG_XZ_LEFT_BND, MAG_XZ_LEFT_DOM, MAG_XZ_RIGHT_BND, MAG_XZ_RIGHT_DOM, MAG_YZ_LEFT_BND, MAG_YZ_LEFT_DOM, MAG_YZ_RIGHT_BND, MAG_YZ_RIGHT_DOM, &
        & MPI_XY_LEFT_BND, MPI_XY_LEFT_DOM, MPI_XY_RIGHT_BND, MPI_XY_RIGHT_DOM, MPI_XZ_LEFT_BND, MPI_XZ_LEFT_DOM, &
        & MPI_XZ_RIGHT_BND, MPI_XZ_RIGHT_DOM, MPI_YZ_LEFT_BND, MPI_YZ_LEFT_DOM, MPI_YZ_RIGHT_BND, MPI_YZ_RIGHT_DOM, bnd_xl, bnd_xl_dom, bnd_xr, &
        & bnd_xr_dom, bnd_yl, bnd_yl_dom, bnd_yr, bnd_yr_dom, bnd_zl, bnd_zl_dom, bnd_zr, bnd_zr_dom, buffer_dim, cbuff, cbuff_len, cfl, cfl_max, cflcontrol, &
        & cfr_smooth, cleanup_mpi, comm, comm3d, dt, dt_initial, dt_max_grow, dt_min, dt_old, dtm, err, ibuff, ierr, info, init_mpi, &
        & integration_order, lbuff, limiter, mpifind, ndims, nproc, nstep, pcoords, proc, procxl, procxr, procxyl, procyl, procyr, procyxl, proczl, &
        & proczr, psize, pxsize, pysize, pzsize, rbuff, req, smalld, smallei, smallp, status, t, use_smalld, magic_mass, local_magic_mass, master, slave, &
        & nxd, nyd, nzd, nb, big_float

   integer :: nproc, proc, ierr , rc, info
   integer :: status(MPI_STATUS_SIZE,4)
   integer, dimension(4) :: req, err

   logical, protected    :: master, slave

   real, parameter       :: dt_default_grow = 2.
   real, parameter       :: big_float =  huge(real(1.0,4))
   real                  :: t, dt, dt_old, dtm
   real, save            :: magic_mass = 0.0
   real, save            :: local_magic_mass = 0.0
   integer               :: nstep

   integer, parameter    :: ndims = 3       ! 3D grid
   integer               :: comm, comm3d
   logical               :: reorder
   integer, dimension(ndims) :: psize, pcoords, coords
   logical, dimension(ndims) :: periods
   integer               ::   procxl, procxr, procyl, procyr, proczl, proczr, procxyl, procyxl, procxyr, procyxr

   integer, parameter               :: buffer_dim=200
   character(len=cbuff_len), dimension(buffer_dim) :: cbuff
   integer,   dimension(buffer_dim) :: ibuff
   real,      dimension(buffer_dim) :: rbuff
   logical,   dimension(buffer_dim) :: lbuff

   logical :: mpi_magic
   integer :: pxsize    !< number of MPI blocks in x-dimension
   integer :: pysize    !< number of MPI blocks in y-dimension
   integer :: pzsize    !< number of MPI blocks in z-dimension

   namelist /MPI_BLOCKS/ pxsize, pysize, pzsize, mpi_magic

   integer, protected :: nxd  !< number of %grid cells in physical domain (without boundary cells) in x-direction (if equal to 1 then x-dimension is reduced to a point with no boundary cells)
   integer, protected :: nyd  !< number of %grid cells in physical domain (without boundary cells) in y-direction (if equal to 1 then y-dimension is reduced to a point with no boundary cells)
   integer, protected :: nzd  !< number of %grid cells in physical domain (without boundary cells) in z-direction (if equal to 1 then z-dimension is reduced to a point with no boundary cells)
   integer, protected :: nb   !< number of boundary cells surrounding the physical domain, same for all directions

   namelist /DOMAIN_SIZES/ nxd, nyd, nzd, nb

   integer, parameter    :: bndlen = 4 !< length of boundary names
   character(len=bndlen) :: bnd_xl     !< type of boundary conditions for the left  x-boundary
   character(len=bndlen) :: bnd_xr     !< type of boundary conditions for the right x-boundary
   character(len=bndlen) :: bnd_yl     !< type of boundary conditions for the left  y-boundary
   character(len=bndlen) :: bnd_yr     !< type of boundary conditions for the right y-boundary
   character(len=bndlen) :: bnd_zl     !< type of boundary conditions for the left  z-boundary
   character(len=bndlen) :: bnd_zr     !< type of boundary conditions for the right z-boundary
   character(len=bndlen) :: bnd_xl_dom, bnd_xr_dom, bnd_yl_dom, bnd_yr_dom, bnd_zl_dom, bnd_zr_dom !< computational domain boundaries

   namelist /BOUNDARIES/ bnd_xl, bnd_xr, bnd_yl, bnd_yr, bnd_zl, bnd_zr

   real    :: dt_initial
   real    :: dt_max_grow
   real    :: dt_min
   real    :: cfl
   real    :: cfl_max
   logical :: use_smalld               !< correct denisty when it gets lower than smalld
   real    :: smallp                   !< artificial infimum for pressure
   real    :: smalld                   !< artificial infimum for density
   real    :: smallc                   !< artificial infimum for freezing speed
   real    :: smallei                  !< artificial infimum for internal energy density
   real    :: cfr_smooth
   integer :: integration_order
   character(len=cbuff_len) :: limiter !< type of flux limiter
   character(len=cbuff_len) :: cflcontrol !< type of cfl control just before each sweep (possibilities: 'none', 'main', 'user')

   namelist /NUMERICAL_SETUP/  cfl, smalld, smallei, integration_order, cfr_smooth, dt_initial, dt_max_grow, dt_min, smallc, smallp, limiter, cflcontrol, use_smalld, cfl_max

   integer, dimension(3) :: domsize   !< local copy of nxd, nyd, nzd

   logical     :: have_mpi

   integer, save :: MPI_XZ_LEFT_BND=-1, MPI_XZ_RIGHT_BND=-1
   integer, save :: MPI_XZ_LEFT_DOM=-1, MPI_XZ_RIGHT_DOM=-1
   integer, save :: MPI_XY_LEFT_BND=-1, MPI_XY_RIGHT_BND=-1
   integer, save :: MPI_XY_LEFT_DOM=-1, MPI_XY_RIGHT_DOM=-1
   integer, save :: MPI_YZ_LEFT_BND=-1, MPI_YZ_RIGHT_BND=-1
   integer, save :: MPI_YZ_LEFT_DOM=-1, MPI_YZ_RIGHT_DOM=-1

   integer, save :: MAG_XZ_LEFT_BND=-1, MAG_XZ_RIGHT_BND=-1
   integer, save :: MAG_XZ_LEFT_DOM=-1, MAG_XZ_RIGHT_DOM=-1
   integer, save :: MAG_XY_LEFT_BND=-1, MAG_XY_RIGHT_BND=-1
   integer, save :: MAG_XY_LEFT_DOM=-1, MAG_XY_RIGHT_DOM=-1
   integer, save :: MAG_YZ_LEFT_BND=-1, MAG_YZ_RIGHT_BND=-1
   integer, save :: MAG_YZ_LEFT_DOM=-1, MAG_YZ_RIGHT_DOM=-1

   integer, save :: ARR_XZ_LEFT_BND=-1, ARR_XZ_RIGHT_BND=-1
   integer, save :: ARR_XZ_LEFT_DOM=-1, ARR_XZ_RIGHT_DOM=-1
   integer, save :: ARR_XY_LEFT_BND=-1, ARR_XY_RIGHT_BND=-1
   integer, save :: ARR_XY_LEFT_DOM=-1, ARR_XY_RIGHT_DOM=-1
   integer, save :: ARR_YZ_LEFT_BND=-1, ARR_YZ_RIGHT_BND=-1
   integer, save :: ARR_YZ_LEFT_DOM=-1, ARR_YZ_RIGHT_DOM=-1

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
!! <tr><td>pxsize   </td><td>1     </td><td>integer</td><td>\copydoc mpisetup::pxsize   </td></tr>
!! <tr><td>pysize   </td><td>1     </td><td>integer</td><td>\copydoc mpisetup::pysize   </td></tr>
!! <tr><td>pzsize   </td><td>1     </td><td>integer</td><td>\copydoc mpisetup::pzsize   </td></tr>
!! <tr><td>mpi_magic</td><td>.true.</td><td>logical</td><td>\copydoc mpisetup::mpi_magic</td></tr>
!! </table>
!! \n \n
!! @b DOMAIN_SIZES
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>nxd</td><td>1</td><td>positive integer    </td><td>\copydoc grid::nxd</td></tr>
!! <tr><td>nyd</td><td>1</td><td>positive integer    </td><td>\copydoc grid::nyd</td></tr>
!! <tr><td>nzd</td><td>1</td><td>positive integer    </td><td>\copydoc grid::nzd</td></tr>
!! <tr><td>nb </td><td>4</td><td>non-negative integer</td><td>\copydoc grid::nb </td></tr>
!! </table>
!! \n \n
!! @b BOUNDARIES
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>bnd_xl</td><td>'per'</td><td>'per', 'ref', 'out', 'outd', 'outh', 'cor'</td><td>\copydoc mpisetup::bnd_xl</td></tr>
!! <tr><td>bnd_xr</td><td>'per'</td><td>'per', 'ref', 'out', 'outd', 'outh'       </td><td>\copydoc mpisetup::bnd_xr</td></tr>
!! <tr><td>bnd_yl</td><td>'per'</td><td>'per', 'ref', 'out', 'outd', 'outh', 'cor'</td><td>\copydoc mpisetup::bnd_yl</td></tr>
!! <tr><td>bnd_yr</td><td>'per'</td><td>'per', 'ref', 'out', 'outd', 'outh'       </td><td>\copydoc mpisetup::bnd_yr</td></tr>
!! <tr><td>bnd_zl</td><td>'per'</td><td>'per', 'ref', 'out', 'outd', 'outh'       </td><td>\copydoc mpisetup::bnd_zl</td></tr>
!! <tr><td>bnd_zr</td><td>'per'</td><td>'per', 'ref', 'out', 'outd', 'outh'       </td><td>\copydoc mpisetup::bnd_zr</td></tr>
!! </table>
!! \n \n
!! @b NUMERICAL_SETUP
!! \n \n
!! <table border="+1">
!! <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!! <tr><td>cfl              </td><td>0.7   </td><td>real value between 0.0 and 1.0       </td><td>\copydoc mpisetup::cfl              </td></tr>
!! <tr><td>cfl_max          </td><td>0.9   </td><td>real value between cfl and 1.0       </td><td>\copydoc mpisetup::cfl_max          </td></tr>
!! <tr><td>cflcontrol       </td><td>       </td><td>string                              </td><td>\copydoc mpisetup::cflcontrol       </td></tr>
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
!! </table>
!! \n \n
!<
      subroutine init_mpi
         use mpi,           only: MPI_COMM_WORLD, MPI_INFO_NULL, MPI_INFO_NULL, MPI_CHARACTER, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_LOGICAL, MPI_PROC_NULL
         use dataio_pub,    only: die, printinfo, msg, cwdlen, hnlen, cwd, ansi_white, ansi_black, warn
         use dataio_pub,    only: par_file, ierrh, namelist_errh, compare_namelist  ! QA_WARN required for diff_nml

         implicit none

         integer :: iproc

         character(LEN=cwdlen) :: cwd_proc
         character(LEN=hnlen)  :: host_proc
         integer               :: pid_proc

         character(LEN=cwdlen), allocatable, dimension(:) :: cwd_all
         character(LEN=hnlen) , allocatable, dimension(:) :: host_all
         integer              , allocatable, dimension(:) :: pid_all

         integer(kind=1)       :: getcwd, hostnm
         integer(kind=4)       :: getpid
         integer :: cwd_status
         logical :: par_file_exist

         call MPI_Init( ierr )
#ifdef VERBOSE
         call printinfo("[mpisetup:init_mpi]: commencing...")
#endif /* VERBOSE */
         call MPI_Comm_rank(MPI_COMM_WORLD, proc, ierr)
         master = (proc == 0)
         slave  = (proc /= 0)
         comm = MPI_COMM_WORLD
         info = MPI_INFO_NULL
         call MPI_Comm_size(comm, nproc, ierr)

         if (allocated(cwd_all) .or. allocated(host_all) .or. allocated(pid_all)) call die("[mpisetup:init_mpi] cwd_all, host_all or pid_all already allocated")

         ! BEWARE if slave it is probably enough to allocate only one element or none at all (may depend on MPI implementation)
         allocate(cwd_all(0:nproc))
         allocate(host_all(0:nproc))
         allocate(pid_all(0:nproc))

         pid_proc   = getpid()
         status     = hostnm(host_proc)
         cwd_status = getcwd(cwd_proc)

         if (cwd_status /= 0) call die("[mpisetup:init_mpi] problems accessing current working directory.")
#ifdef DEBUG
         write(msg,'(3a,i6,3a)') 'mpisetup: host="',trim(host_proc),'", PID=',pid_proc,' CWD="',trim(cwd_proc),'"'
         call printinfo(msg)
#endif /* DEBUG */

         call MPI_Gather(cwd_proc,  cwdlen, MPI_CHARACTER, cwd_all,  cwdlen, MPI_CHARACTER, 0, comm, err)
         call MPI_Gather(host_proc, hnlen,  MPI_CHARACTER, host_all, hnlen,  MPI_CHARACTER, 0, comm, err)
         call MPI_Gather(pid_proc,  1,      MPI_INTEGER,   pid_all,  1,      MPI_INTEGER,   0, comm, err)

         ! cwd = trim(cwd_proc)  !BEWARE: It's redundant, we get cwd for command line in init_piernik subroutine

         if (master) then
            inquire(file=par_file, exist=par_file_exist)
            if (.not. par_file_exist) call die('[mpisetup:init_mpi] Cannot find "problem.par" in the working directory',0)

            call printinfo("------------------------------------------------------------------------------------------------------", .false.)
            call printinfo("###############     Environment     ###############", .false.)
            call printinfo("", .false.)
            call printinfo("PROCESSES:", .false.)
            do iproc = 0, nproc-1
               write(msg,"(a6,i2,a7,i6,a1,a,a7,a)") " proc=",iproc, &
                    ", pid= ",pid_all(iproc), "@",trim(host_all(iproc)), &
                    ", cwd=",trim(cwd)
               call printinfo(msg, .false.)
            enddo
            call printinfo("", .true.)
            write(msg,"(5a,i5)") 'Start of the',ansi_white,' PIERNIK ',ansi_black,'code. No. of procs = ', nproc
            call printinfo(msg, .true.)
            call printinfo("", .true.)
            call printinfo("###############     Namelist parameters     ###############", .false.)
         endif

         pxsize = 1
         pysize = 1
         pzsize = 1
         nxd    = 1
         nyd    = 1
         nzd    = 1
         nb     = 4

         mpi_magic = .true.

         bnd_xl = 'per'
         bnd_xr = 'per'
         bnd_yl = 'per'
         bnd_yr = 'per'
         bnd_zl = 'per'
         bnd_zr = 'per'

         ! Provide backward compatibility for choosing limiter via preprocessor flag
         ! ToDo: Remove it when all problems are fixed

         limiter     = 'vanleer'
         cflcontrol  = 'warn'

         cfl         = 0.7
         cfl_max     = 0.9
         cfr_smooth  = 0.0
         smallp      = big_float
         smalld      = big_float
         use_smalld  = .true.
         smallc      = 1.e-10
         smallei     = 1.e-10
         dt_initial  = -1.              !< negative value indicates automatic choice of initial timestep
         dt_max_grow = dt_default_grow !< for sensitive setups consider setting this as low as 1.1
         dt_min      = tiny(1.)

         integration_order  = 2

         if (master) then
            diff_nml(MPI_BLOCKS)
            diff_nml(BOUNDARIES)
            diff_nml(NUMERICAL_SETUP)
            diff_nml(DOMAIN_SIZES)
         endif

         nxd = max(1, nxd)
         nyd = max(1, nyd)
         nzd = max(1, nzd)

         cfl_max = min(max(cfl_max, min(cfl*1.1, cfl+0.05, (1.+cfl)/2.) ), 1.0) ! automatically sanitize cfl_max

         if (master) then

            cbuff(1) = bnd_xl
            cbuff(2) = bnd_xr
            cbuff(3) = bnd_yl
            cbuff(4) = bnd_yr
            cbuff(5) = bnd_zl
            cbuff(6) = bnd_zr
            cbuff(7) = limiter
            cbuff(8) = cflcontrol

            ibuff(1) = pxsize
            ibuff(2) = pysize
            ibuff(3) = pzsize
            ibuff(4) = integration_order
            ibuff(5) = nxd
            ibuff(6) = nyd
            ibuff(7) = nzd
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

            lbuff(1) = mpi_magic
            lbuff(2) = use_smalld

         endif

         call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        0, comm, ierr)
         call MPI_Bcast(ibuff,           buffer_dim, MPI_INTEGER,          0, comm, ierr)
         call MPI_Bcast(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, 0, comm, ierr)
         call MPI_Bcast(lbuff,           buffer_dim, MPI_LOGICAL,          0, comm, ierr)

         if (slave) then

            mpi_magic   = lbuff(1)
            use_smalld  = lbuff(2)

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

            bnd_xl     = cbuff(1)(1:4)
            bnd_xr     = cbuff(2)(1:4)
            bnd_yl     = cbuff(3)(1:4)
            bnd_yr     = cbuff(4)(1:4)
            bnd_zl     = cbuff(5)(1:4)
            bnd_zr     = cbuff(6)(1:4)
            limiter    = cbuff(7)
            cflcontrol = cbuff(8)

            pxsize = ibuff(1)
            pysize = ibuff(2)
            pzsize = ibuff(3)

            integration_order = ibuff(4)

            nxd = ibuff(5)
            nyd = ibuff(6)
            nzd = ibuff(7)
            nb  = ibuff(8)

         endif

         domsize(:) = [nxd, nyd, nzd]

         bnd_xl_dom = bnd_xl
         bnd_xr_dom = bnd_xr
         bnd_yl_dom = bnd_yl
         bnd_yr_dom = bnd_yr
         bnd_zl_dom = bnd_zl
         bnd_zr_dom = bnd_zr

         psize(:) = [ pxsize, pysize, pzsize ]

         if (pxsize*pysize*pzsize /= nproc) then
            if (mpi_magic) then
               call divide_domain_voodoo(nproc)
            else
               if (master) then
                  write(msg,'(A,I5,A,I10)') 'nproc =',nproc,' MUST BE EQUAL TO   pxsize*pysize*pzsize =',pxsize*pysize*pzsize
                  call die(msg,0)
               endif
               call MPI_Barrier(MPI_COMM_WORLD, ierr)
               call MPI_Finalize(ierr)
            endif
         endif

         if ( (bnd_xl(1:3) == 'cor' .or. bnd_yl(1:3) == 'cor' .or. bnd_xr(1:3) == 'cor' .or. bnd_yr(1:3) == 'cor') .and. (pxsize /= pysize) ) then
            write(msg, '(a,2(i3,a))')"[mpisetup:init_mpi] Corner BC require equal pxsize and pysize. Detected: [",pxsize,",",pysize,"]"
            call die(msg)
         endif

         if (pxsize*pysize*pzsize /= 1) have_mpi = .true.

         periods(:) = .false.

         if (bnd_xl(1:3) == 'per' .or. bnd_xl(1:3) == 'she') then
            periods(1) = .true.  ! x periodic
            if (bnd_xr(1:3) /= bnd_xl(1:3)) call die("[mpisetup:init_mpi] Periodic or shear BC do not match in X-direction")
         endif

         if (bnd_yl(1:3) == 'per') then
            periods(2) = .true.  ! y periodic
            if (bnd_yr(1:3) /= bnd_yl(1:3)) call die("[mpisetup:init_mpi] Periodic BC do not match in Y-direction")
         endif

         if (bnd_zl(1:3) == 'per') then
            periods(3) = .true.  ! z periodic
            if (bnd_zr(1:3) /= bnd_zl(1:3)) call die("[mpisetup:init_mpi] Periodic BC do not match in Z-direction")
         endif

         reorder = .false.     ! allows processes reordered for efficiency

         call MPI_Cart_create(comm, ndims, psize, periods, reorder, comm3d, ierr)
         call MPI_Cart_coords(comm3d, proc, ndims, pcoords, ierr)

         deallocate(host_all)
         deallocate(pid_all)
         deallocate(cwd_all)

! Compute neighbors

         call MPI_Cart_shift(comm3d,0,1,procxl,procxr,ierr)   ! x dim
         call MPI_Cart_shift(comm3d,1,1,procyl,procyr,ierr)   ! y dim
         call MPI_Cart_shift(comm3d,2,1,proczl,proczr,ierr)   ! z dim

         if (bnd_xl(1:3) == 'cor' .and. bnd_yl(1:3) == 'cor' ) then
            if (pcoords(1) == 0 .and. pcoords(2) > 0) then
               coords = (/pcoords(2),pcoords(1),pcoords(3)/)
               call MPI_Cart_rank(comm3d,coords,procxyl,ierr)
            else
               procxyl = MPI_PROC_NULL
            endif
            if (pcoords(2) == 0 .and. pcoords(1) > 0 ) then
               coords = (/pcoords(2),pcoords(1),pcoords(3)/)
               call MPI_Cart_rank(comm3d,coords,procyxl,ierr)
            else
               procyxl = MPI_PROC_NULL
            endif
         endif

         if (bnd_xr(1:3) == 'cor' .and. bnd_yr(1:3) == 'cor' ) then
            if (pcoords(1) == psize(1)-1 .and. pcoords(2) < psize(2)-1) then
               coords = (/pcoords(2),pcoords(1),pcoords(3)/)
               call MPI_Cart_rank(comm3d,coords,procxyr,ierr)
            else
               procxyr = MPI_PROC_NULL
            endif
            if (pcoords(2) == psize(2)-1 .and. pcoords(1) < psize(2)-1 ) then
               coords = (/pcoords(2),pcoords(1),pcoords(3)/)
               call MPI_Cart_rank(comm3d,coords,procyxr,ierr)
            else
               procyxr = MPI_PROC_NULL
            endif
         endif

#ifdef SHEAR_BND
         if (pysize > 1) stop 'Shear-pediodic boundary conditions do not permit pysize > 1'

#ifndef FFTW
         if (pcoords(1) == 0) then
            bnd_xl = 'she'
         else
            bnd_xl = 'mpi'
         endif

         if (pcoords(1) == pxsize-1) then
            bnd_xr = 'she'
         else
            bnd_xr = 'mpi'
         endif
#endif /* !FFTW */

#else /* !SHEAR_BND */
         if (procxl /= MPI_PROC_NULL .and. procxl /= proc) bnd_xl = 'mpi'
         if (procxr /= MPI_PROC_NULL .and. procxr /= proc) bnd_xr = 'mpi'
#endif /* !SHEAR_BND */

         if (procyl /= MPI_PROC_NULL .and. procyl /= proc) bnd_yl = 'mpi'
         if (procyr /= MPI_PROC_NULL .and. procyr /= proc) bnd_yr = 'mpi'

         if (proczl /= MPI_PROC_NULL .and. proczl /= proc) bnd_zl = 'mpi'
         if (proczr /= MPI_PROC_NULL .and. proczr /= proc) bnd_zr = 'mpi'

#ifdef DEBUG
         write(msg,*) 'xdir: ',procxl, proc, procxr
         call printinfo(msg)
         write(msg,*) 'ydir: ',procyl, proc, procyr
         call printinfo(msg)
         write(msg,*) 'zdir: ',proczl, proc, proczr
         call printinfo(msg)
#endif /* DEBUG */

         if (integration_order > 2) call die ('[mpisetup:init_mpi]: "ORIG" scheme integration_order must be 1 or 2')

         dt_old = -1.
         if (dt_max_grow < 1.01) then
            if (master) then
               write(msg,'(2(a,g10.3))')"[mpisetup:init_mpi] dt_max_grow = ",dt_max_grow," is way too low. Resetting to ",dt_default_grow
               call warn(msg)
            endif
            dt_max_grow = dt_default_grow
         endif
#ifdef VERBOSE
         call printinfo("[mpisetup:init_mpi]: finished. \o/")
#endif /* VERBOSE */

      end subroutine init_mpi

!-----------------------------------------------------------------------------

      subroutine cleanup_mpi

         use dataio_pub,    only: printinfo

         implicit none

         call MPI_Comm_free(comm3d, ierr)

         if (domsize(1) /= 1) then
            call MPI_Type_free(MPI_YZ_LEFT_BND, ierr)
            call MPI_Type_free(MPI_YZ_LEFT_DOM, ierr)
            call MPI_Type_free(MPI_YZ_RIGHT_DOM, ierr)
            call MPI_Type_free(MPI_YZ_RIGHT_BND, ierr)
            call MPI_Type_free(MAG_YZ_LEFT_BND, ierr)
            call MPI_Type_free(MAG_YZ_LEFT_DOM, ierr)
            call MPI_Type_free(MAG_YZ_RIGHT_DOM, ierr)
            call MPI_Type_free(MAG_YZ_RIGHT_BND, ierr)
            call MPI_Type_free(ARR_YZ_LEFT_BND, ierr)
            call MPI_Type_free(ARR_YZ_LEFT_DOM, ierr)
            call MPI_Type_free(ARR_YZ_RIGHT_DOM, ierr)
            call MPI_Type_free(ARR_YZ_RIGHT_BND, ierr)
         endif

         if (domsize(2) /= 1) then
            call MPI_Type_free(MPI_XZ_LEFT_BND, ierr)
            call MPI_Type_free(MPI_XZ_LEFT_DOM, ierr)
            call MPI_Type_free(MPI_XZ_RIGHT_DOM, ierr)
            call MPI_Type_free(MPI_XZ_RIGHT_BND, ierr)
            call MPI_Type_free(MAG_XZ_LEFT_BND, ierr)
            call MPI_Type_free(MAG_XZ_LEFT_DOM, ierr)
            call MPI_Type_free(MAG_XZ_RIGHT_DOM, ierr)
            call MPI_Type_free(MAG_XZ_RIGHT_BND, ierr)
            call MPI_Type_free(ARR_XZ_LEFT_BND, ierr)
            call MPI_Type_free(ARR_XZ_LEFT_DOM, ierr)
            call MPI_Type_free(ARR_XZ_RIGHT_DOM, ierr)
            call MPI_Type_free(ARR_XZ_RIGHT_BND, ierr)
         endif

         if (domsize(3) /= 1) then
            call MPI_Type_free(MPI_XY_LEFT_BND, ierr)
            call MPI_Type_free(MPI_XY_LEFT_DOM, ierr)
            call MPI_Type_free(MPI_XY_RIGHT_DOM, ierr)
            call MPI_Type_free(MPI_XY_RIGHT_BND, ierr)
            call MPI_Type_free(MAG_XY_LEFT_BND, ierr)
            call MPI_Type_free(MAG_XY_LEFT_DOM, ierr)
            call MPI_Type_free(MAG_XY_RIGHT_DOM, ierr)
            call MPI_Type_free(MAG_XY_RIGHT_BND, ierr)
            call MPI_Type_free(ARR_XY_LEFT_BND, ierr)
            call MPI_Type_free(ARR_XY_LEFT_DOM, ierr)
            call MPI_Type_free(ARR_XY_RIGHT_DOM, ierr)
            call MPI_Type_free(ARR_XY_RIGHT_BND, ierr)
         endif

         if (master) call printinfo("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++", .false.)
         call MPI_Barrier(comm,ierr)
         if (nproc > 1) call sleep(1) ! Prevent random SIGSEGVs in openmpi's MPI_Finalize
         call MPI_Finalize(ierr)

      end subroutine cleanup_mpi

!-----------------------------------------------------------------------------

      subroutine mpifind(var, what, loc_arr, loc_proc)
         use mpi,           only: MPI_2DOUBLE_PRECISION, MPI_INTEGER, MPI_MINLOC, MPI_MAXLOC
         use dataio_pub,    only: msg, warn

         implicit none

         character(len=*), intent(in) :: what
         real                         :: var
         real, dimension(2)           :: rsend, rrecv
         integer, dimension(3)        :: loc_arr
         integer                      :: loc_proc

         rsend(1) = var
         rsend(2) = proc

         select case (what(1:3))
            case ('min')
               CALL MPI_Reduce(rsend, rrecv, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, 0, comm, ierr)
            case ('max')
               CALL MPI_Reduce(rsend, rrecv, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, 0, comm, ierr)
            case default
               write(msg,*) '[mpisetup:mpifind] actual parameter "', what, '"is not allowed'
               call warn(msg)
         end select

         if (master) then
            var = rrecv(1)
            loc_proc = rrecv(2)
         endif

         call MPI_Bcast(loc_proc, 1, MPI_INTEGER, 0, comm, ierr)

         if (loc_proc /= 0) then
            if (proc == loc_proc) then
               CALL MPI_Send  (loc_arr, 3, MPI_INTEGER,     0, 11, comm, ierr)
            else if (master) then
               CALL MPI_Recv  (loc_arr, 3, MPI_INTEGER, loc_proc, 11, comm, status, ierr)
            endif
         endif

      end subroutine mpifind

!------------------------------------------------------------------------------------------
! Must be called by all procs to avoid communication and ensure that every proc has
! proper psize, pxsize, pysize, pzsize
!
! This routine tries to divide the computational domain into local domains.
! The goal is to minimize the ratio of longest to shortest edge to minimize the amount of inter-process communication.
! If the benchmarks show that some direction should be partitioned in more pieces than other directions, implement appropriate weighting in j1, j2 and j3 calculations.
!
! For some weird domains and PE counts this routine may find tiling that does not satisfy multigrid restrictions even if there is some. In such case divide domain manually.
!

      subroutine divide_domain_voodoo(np)

         use dataio_pub,    only: die, printinfo, msg

         implicit none
         integer, parameter, dimension(26) :: some_primes = [ 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101 ]

         integer, intent(in) :: np

         integer :: j1, j2, j3, jj, n, p
         integer, dimension(3) :: ldom

         ldom(1:3) = domsize(3:1:-1) ! Maxloc returns first occurrence of max, reversing direction order (to ZYX) gives better cache utilization.
         n = np
         psize(:) = 1

         do p = size(some_primes), 1, -1 ! start from largest defined primes, continue down to 2
            do while (mod(n, some_primes(p))==0)

               jj = 0
               j1 = sum(maxloc(ldom), 1) ! First try the longest edge; note the trick to make a scalar from 1-element vector without assignment to another variable
               if (mod(ldom(j1), some_primes(p))==0) then
                  jj = j1
               else
                  j2 = 1 + mod(j1 + 0, ndims)
                  j3 = 1 + mod(j1 + ndims -2, ndims)
                  if (ldom(j2) > ldom(j3) .and. mod(ldom(j2), some_primes(p))==0) jj = j2 ! middle edge ...
                  if (jj == 0 .and. mod(ldom(j3), some_primes(p))==0) jj = j3 ! try the shortest edge on last resort
               endif

               if (jj == 0) then
                  call die("[divide_domain_voodoo]: Can't find divisible edge")
               else
                  psize(jj) = psize(jj) * some_primes(p)
                  n         = n         / some_primes(p)
                  ldom(jj)  = ldom(jj)  / some_primes(p)
               endif

            enddo
         enddo

         if (n /= 1) call die("[divide_domain_voodoo]: I am not that intelligent") ! np has too big prime factors

         pxsize = psize(3) ! directions were reverted at ldom assignment
         pysize = psize(2)
         pzsize = psize(1)

         psize = [ pxsize, pysize, pzsize ]

         if (master .and. np > 1) then
            write(msg,'(a,3i4,a,3i6,a)')"[mpisetup:divide_domain_voodoo] Domain divided to [",psize(:)," ] pieces, each of [",ldom(3:1:-1)," ] cells."
            call printinfo(msg)
         endif

      end subroutine divide_domain_voodoo

end module mpisetup
