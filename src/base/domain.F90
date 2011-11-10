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
!! This module contains everything closely related to management of computational domain.
!!
!! In this module following namelists of parameters are specified:
!! \copydetails domain::init_domain
!<
module domain

   use constants, only: ndims, cbuff_len, LO, HI

   implicit none

   private
   public :: cleanup_domain, init_domain, is_overlap, domain_container, user_divide_domain, allocate_pse, deallocate_pse, set_pse_sel, &
        &    pdom, dom, is_uneven, is_mpi_noncart, is_refined, is_multicg, cdd

! AMR: There will be at least one domain container for the base grid.
!      It will be possible to host one or more refined domains on the base container and on the refined containers.
!      The refined domains may cover whole parent domain or only a tiny part of it.
!      The multigrid solver will operate also on a stack of coarser domains - parents of the base domain.
!      The coarser domains must be no smaller than the base domain.

   type :: cuboids
      integer(kind=8), dimension(:,:,:), allocatable :: sel !< list of grid chunks (:, xdim:zdim, LO:HI)
   end type cuboids

   type :: domain_container
      ! primary parameters, read from /DOMAIN_SIZES/, /BOUNDARIES/ and /DOMAIN_LIMITS/ namelists
      real, dimension(ndims, LO:HI) :: edge     !< physical domain boundary positions
      integer(kind=4), dimension(ndims) :: n_d  !< number of grid cells in physical domain in x-, y- and z-direction (where equal to 1, the dimension is reduced to a point with no boundary cells)
      integer(kind=4)                  :: nb    !< number of boundary cells surrounding the physical domain, same for all directions
      integer, dimension(ndims, LO:HI) :: bnd   !< type of boundary conditions coded in integers

      ! derived parameters
      real, dimension(ndims) :: L_              !< span of the physical domain [ xmax-xmin, ymax-ymin, zmax-zmin ]
      real, dimension(ndims) :: C_              !< center of the physical domain [ (xmax+xmin)/2., (ymax+ymin)/2., (zmax+zmin)/2. ]
      real    :: Vol                            !< total volume of the physical domain

      logical, dimension(ndims) :: periodic     !< .true. for periodic and shearing boundary pairs

      type(cuboids), dimension(:), allocatable :: pse  !< lists of grid chunks on each process (FIRST:LAST); Use with care, because this is an antiparallel thing

      ! Do not use n_t(:) in the Piernik source tree without a good reason
      ! Avoid as a plague allocating buffers of that size because it negates benefits of parallelization
      ! \todo move them to another type, that extends domain_container?
      integer(kind=4), dimension(ndims) :: n_t          !< total number of %grid cells in the whole domain in every direction (n_d(:) + 2* nb for existing directions)

      integer(kind=8) :: total_ncells !< total number of %grid cells
      integer :: geometry_type  !< the type of geometry: cartesian: GEO_XYZ, cylindrical: GEO_RPZ, other: GEO_INVALID
      integer :: D_x                  !< set to 1 when x-direction exists, 0 otherwise
      integer :: D_y                  !< set to 1 when y-direction exists, 0 otherwise.
      integer :: D_z                  !< set to 1 when z-direction exists, 0 otherwise.

      integer, dimension(ndims) :: D_ !< set to 1 for existing directions, 0 otherwise. Useful for dimensionally-safe indices for difference operators on arrays,

      logical, dimension(ndims) :: has_dir   !< .true. for existing directions
      integer :: eff_dim                     !< effective dimensionality of the simulation

    contains

      procedure :: set_derived
      procedure :: translate_bnds_to_ints
      procedure :: print_me
      procedure :: init => init_domain_container

   end type domain_container

   type cart_decomposition
      integer(kind=4)                          :: comm3d  !< cartesian communicator
      integer(kind=4), dimension(ndims)        :: psize   !< number of divisions in each direction
      integer(kind=4), dimension(ndims)        :: pcoords !< own process coordinates within psize(:)-shaped array of processes
      integer(kind=4), dimension(ndims, LO:HI) :: procn   !< array of neighbours proc numbers
      integer(kind=4)                          :: procxyl !< neighbour in corner boundaries
      integer(kind=4)                          :: procyxl !< neighbour in corner boundaries
   end type cart_decomposition

   type(cart_decomposition), protected :: cdd !< Cartesian Domain Decomposition stuff

   type(domain_container), protected, target :: dom !< complete description of base level domain
   type(domain_container), pointer :: pdom

   logical, protected :: is_uneven          !< .true. when n_b(:) depend on process rank
   logical, protected :: is_mpi_noncart     !< .true. when there exist a process that has more than one neighbour in any direction
   logical, protected :: is_refined         !< .true. when AMR or static refinement is employed
   logical, protected :: is_multicg         !< .true. when at leas one process has more than one grid container

   ! Private variables

   logical :: dom_divided  !< Flags that domain decomposition was succesfull and quality of solution meets the criteria
   integer(kind=4), allocatable, dimension(:) :: primes
   real :: ideal_bsize

   ! Namelist variables

   integer(kind=4), dimension(ndims) :: psize !< desired number of MPI blocks in x, y and z-dimension
   integer(kind=4), dimension(ndims) :: bsize !< the size of cg for multiblock decomposition
   logical :: reorder                 !< allows processes reordered for efficiency (a parameter of MPI_Cart_create and MPI_graph_create)
   logical :: allow_uneven            !< allows different values of n_b(:) on different processes
   logical :: allow_noncart           !< allows more than one neighbour on a boundary
   logical :: allow_AMR               !< allows AMR
   real    :: dd_unif_quality         !< uniform domain decomposition may be rejected it its quality is below this threshold (e.g. very elongated local domains are found)
   real    :: dd_rect_quality         !< rectilinear domain decomposition may be rejected it its quality is below this threshold (not used yet)
   logical :: use_comm3d              !< If .false. then do not call any MPI_Cart_* functions

   namelist /MPI_BLOCKS/ psize, bsize, reorder, allow_uneven, allow_noncart, allow_AMR, dd_unif_quality, dd_rect_quality, use_comm3d

   integer(kind=4), dimension(ndims) :: n_d !< number of %grid cells in physical domain without boundary cells (where  == 1 then that dimension is reduced to a point with no boundary cells)
   !< \deprecated nxd, nyd and nzd are obsoleted by the above definition. Will be removed soon
   integer(kind=4), protected :: nxd    !< number of %grid cells in physical domain (without boundary cells) in x-direction (if == 1 then x-dimension is reduced to a point with no boundary cells)
   integer(kind=4), protected :: nyd    !< number of %grid cells in physical domain (without boundary cells) in y-direction (-- || --)
   integer(kind=4), protected :: nzd    !< number of %grid cells in physical domain (without boundary cells) in z-direction (-- || --)
   integer(kind=4), protected :: nb     !< number of boundary cells surrounding the physical domain, same for all directions
   character(len=cbuff_len) :: bnd_xl   !< type of boundary conditions for the left  x-boundary
   character(len=cbuff_len) :: bnd_xr   !< type of boundary conditions for the right x-boundary
   character(len=cbuff_len) :: bnd_yl   !< type of boundary conditions for the left  y-boundary
   character(len=cbuff_len) :: bnd_yr   !< type of boundary conditions for the right y-boundary
   character(len=cbuff_len) :: bnd_zl   !< type of boundary conditions for the left  z-boundary
   character(len=cbuff_len) :: bnd_zr   !< type of boundary conditions for the right z-boundary
   character(len=cbuff_len), dimension(HI*ndims) :: bnds  !< Six strings, describing boundary conditions
   real, dimension(ndims, LO:HI) :: edges
   real :: xmin                         !< physical domain left x-boundary position
   real :: xmax                         !< physical domain right x-boundary position
   real :: ymin                         !< physical domain left y-boundary position
   real :: ymax                         !< physical domain right y-boundary position
   real :: zmin                         !< physical domain left z-boundary position
   real :: zmax                         !< physical domain right z-boundary position
   character(len=cbuff_len) :: geometry !< define system of coordinates: "cartesian" or "cylindrical"

   namelist /BASE_DOMAIN/ n_d, nxd, nyd, nzd, nb, bnd_xl, bnd_xr, bnd_yl, bnd_yr, bnd_zl, bnd_zr, xmin, xmax, ymin, ymax, zmin, zmax, geometry
   !namelist /BASE_DOMAIN/ dom, geometry, nb

   interface is_overlap
      module procedure is_overlap_simple, is_overlap_per
   end interface

   interface
      subroutine divide_domain_template(dom_divided)

         implicit none

         logical, intent(inout) :: dom_divided

      end subroutine divide_domain_template
   end interface

   procedure(divide_domain_template), pointer :: user_divide_domain => Null() ! It should at least allocate and set up dom%pse(:)%sel(:,:,:)

contains

!-----------------------------------------------------------------------------
!>
!! \brief Routine to set up computational domain and define its decomposition
!!
!! \n \n
!! @b MPI_BLOCKS
!! \n \n
!! <table border="+1">
!!   <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!!   <tr><td>psize(3)       </td><td>1      </td><td>integer</td><td>\copydoc domain::psize          </td></tr>
!!   <tr><td>bsize(3)       </td><td>0      </td><td>integer</td><td>\copydoc domain::bsize          </td></tr>
!!   <tr><td>reorder        </td><td>.false.</td><td>logical</td><td>\copydoc domain::reorder        </td></tr>
!!   <tr><td>allow_uneven   </td><td>.false.</td><td>logical</td><td>\copydoc domain::allow_uneven   </td></tr>
!!   <tr><td>allow_noncart  </td><td>.false.</td><td>logical</td><td>\copydoc domain::allow_noncart  </td></tr>
!!   <tr><td>allow_AMR      </td><td>.false.</td><td>logical</td><td>\copydoc domain::allow_AMR      </td></tr>
!!   <tr><td>dd_unif_quality</td><td>0.9    </td><td>real   </td><td>\copydoc domain::dd_unif_quality</td></tr>
!!   <tr><td>dd_rect_quality</td><td>0.9    </td><td>real   </td><td>\copydoc domain::dd_rect_quality</td></tr>
!!   <tr><td>use_comm3d     </td><td>.true. </td><td>logical</td><td>\copydoc domain::use_comm3d     </td></tr>
!! </table>
!! \n \n
!! @b BASE_DOMAIN
!! \n \n
!! <table border="+1">
!!   <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!!   <tr><td>nxd     </td><td>1          </td><td>positive integer                          </td><td>\copydoc domain::nxd     </td></tr>
!!   <tr><td>nyd     </td><td>1          </td><td>positive integer                          </td><td>\copydoc domain::nyd     </td></tr>
!!   <tr><td>nzd     </td><td>1          </td><td>positive integer                          </td><td>\copydoc domain::nzd     </td></tr>
!!   <tr><td>n_d(3)  </td><td>1          </td><td>positive integer                          </td><td>\copydoc domain::n_d     </td></tr>
!!   <tr><td>nb      </td><td>4          </td><td>positive integer                          </td><td>\copydoc domain::nb      </td></tr>
!!   <tr><td>bnd_xl  </td><td>'per'      </td><td>'per', 'ref', 'out', 'outd', 'outh', 'cor'</td><td>\copydoc domain::bnd_xl  </td></tr>
!!   <tr><td>bnd_xr  </td><td>'per'      </td><td>'per', 'ref', 'out', 'outd', 'outh'       </td><td>\copydoc domain::bnd_xr  </td></tr>
!!   <tr><td>bnd_yl  </td><td>'per'      </td><td>'per', 'ref', 'out', 'outd', 'outh', 'cor'</td><td>\copydoc domain::bnd_yl  </td></tr>
!!   <tr><td>bnd_yr  </td><td>'per'      </td><td>'per', 'ref', 'out', 'outd', 'outh'       </td><td>\copydoc domain::bnd_yr  </td></tr>
!!   <tr><td>bnd_zl  </td><td>'per'      </td><td>'per', 'ref', 'out', 'outd', 'outh'       </td><td>\copydoc domain::bnd_zl  </td></tr>
!!   <tr><td>bnd_zr  </td><td>'per'      </td><td>'per', 'ref', 'out', 'outd', 'outh'       </td><td>\copydoc domain::bnd_zr  </td></tr>
!!   <tr><td>xmin    </td><td>0.         </td><td>real                                      </td><td>\copydoc domain::xmin    </td></tr>
!!   <tr><td>xmax    </td><td>1.         </td><td>real                                      </td><td>\copydoc domain::xmax    </td></tr>
!!   <tr><td>ymin    </td><td>0.         </td><td>real                                      </td><td>\copydoc domain::ymin    </td></tr>
!!   <tr><td>ymax    </td><td>1.         </td><td>real                                      </td><td>\copydoc domain::ymax    </td></tr>
!!   <tr><td>zmin    </td><td>0.         </td><td>real                                      </td><td>\copydoc domain::zmin    </td></tr>
!!   <tr><td>zmax    </td><td>1.         </td><td>real                                      </td><td>\copydoc domain::zmax    </td></tr>
!!   <tr><td>geometry</td><td>"cartesian"</td><td>character(len=cbuff_len)                  </td><td>\copydoc domain::geometry</td></tr>
!! </table>
!! \n \n
!<
   subroutine init_domain

      use constants,  only: xdim, zdim, LO, HI, PIERNIK_INIT_MPI, I_ONE, I_ZERO, big_float
      use dataio_pub, only: die, printinfo, msg, warn, code_progress
      use dataio_pub, only: par_file, ierrh, namelist_errh, compare_namelist, cmdl_nml, lun  ! QA_WARN required for diff_nml
      use mpi,        only: MPI_COMM_NULL, MPI_PROC_NULL, MPI_CHARACTER, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_LOGICAL, MPI_IN_PLACE, MPI_LOR, MPI_LAND
      use mpisetup,   only: buffer_dim, cbuff, ibuff, lbuff, rbuff, master, slave, proc, FIRST, LAST, nproc, comm, ierr, have_mpi

      implicit none

      integer :: p, i
      real, allocatable, dimension(:) :: maxcnt

      if (code_progress < PIERNIK_INIT_MPI) call die("[domain:init_domain] MPI not initialized.")

      ! Begin processing of namelist parameters

      psize(:) = I_ONE
      bsize(:) = I_ZERO

      nxd    = 1
      nyd    = 1
      nzd    = 1
      n_d(:) = I_ONE
      nb     = 4
      xmin = 0.; xmax = 1.
      ymin = -big_float; ymax = big_float
      zmin = 0.; zmax = 1.
      geometry = "cartesian"

      reorder   = .false.      !< \todo test it!
      allow_uneven = .true.
      allow_noncart = .false.  !< experimental implementation
      allow_AMR = .false.      !< not implemented yet
      use_comm3d = .true.      !< \todo make a big benchmark with and without comm3d

      bnd_xl = 'per'
      bnd_xr = 'per'
      bnd_yl = 'per'
      bnd_yr = 'per'
      bnd_zl = 'per'
      bnd_zr = 'per'

      dd_unif_quality = 0.9
      dd_rect_quality = 0.9

      if (master) then
         diff_nml(MPI_BLOCKS)
         diff_nml(BASE_DOMAIN)

         if (any([nxd, nyd, nzd] > 1)) call warn("[domain:init_domain] Use n_d instead of nxd, nyd and nzd")
         n_d(:) = max(n_d(:), [nxd, nyd, nzd])

         if (any(bsize(:) > 0 .and. bsize(:) < nb .and. n_d(:) > 1)) call die("[domain:init_domain] bsize(:) is too small.")

         cbuff(1) = bnd_xl
         cbuff(2) = bnd_xr
         cbuff(3) = bnd_yl
         cbuff(4) = bnd_yr
         cbuff(5) = bnd_zl
         cbuff(6) = bnd_zr
         cbuff(7) = geometry

         ibuff(         xdim:zdim) = psize(:)
         ibuff(  zdim+xdim:2*zdim) = n_d(:)
         ibuff(2*zdim+xdim:3*zdim) = bsize(:)
         ibuff(3*zdim+1)           = nb

         rbuff(1) = xmin
         rbuff(2) = xmax
         rbuff(3) = ymin
         rbuff(4) = ymax
         rbuff(5) = zmin
         rbuff(6) = zmax
         rbuff(7) = dd_unif_quality
         rbuff(8) = dd_rect_quality

         lbuff(1) = reorder
         lbuff(2) = allow_uneven
         lbuff(3) = allow_noncart
         lbuff(4) = allow_AMR
         lbuff(5) = use_comm3d

      endif

      call MPI_Bcast(cbuff, cbuff_len*buffer_dim, MPI_CHARACTER,        FIRST, comm, ierr)
      call MPI_Bcast(ibuff,           buffer_dim, MPI_INTEGER,          FIRST, comm, ierr)
      call MPI_Bcast(rbuff,           buffer_dim, MPI_DOUBLE_PRECISION, FIRST, comm, ierr)
      call MPI_Bcast(lbuff,           buffer_dim, MPI_LOGICAL,          FIRST, comm, ierr)

      if (slave) then

         reorder       = lbuff(1)
         allow_uneven  = lbuff(2)
         allow_noncart = lbuff(3)
         allow_AMR     = lbuff(4)
         use_comm3d    = lbuff(5)

         xmin            = rbuff(1)
         xmax            = rbuff(2)
         ymin            = rbuff(3)
         ymax            = rbuff(4)
         zmin            = rbuff(5)
         zmax            = rbuff(6)
         dd_unif_quality = rbuff(7)
         dd_rect_quality = rbuff(8)

         bnd_xl     = cbuff(1)
         bnd_xr     = cbuff(2)
         bnd_yl     = cbuff(3)
         bnd_yr     = cbuff(4)
         bnd_zl     = cbuff(5)
         bnd_zr     = cbuff(6)
         geometry   = cbuff(7)

         psize(:)   = int(ibuff(         xdim:zdim), kind=4)
         n_d(:)     = int(ibuff(  zdim+xdim:2*zdim), kind=4)
         bsize(:)   = int(ibuff(2*zdim+xdim:3*zdim), kind=4)
         nb         = int(ibuff(3*zdim+1),           kind=4)

      endif

      bnds = [bnd_xl, bnd_xr, bnd_yl, bnd_yr, bnd_zl, bnd_zr]
      edges = reshape( [xmin, ymin, zmin, xmax, ymax, zmax], shape=[ndims,HI-LO+I_ONE] )

      pdom => dom

      call dom%init(nb, n_d, bnds, edges, geometry)

#ifdef MULTIGRID
      if (allow_AMR .and. master) call warn("[domain:init_domain] Multigrid solver is not yet capable of using AMR domains.")
      allow_AMR = .false.
#endif /* MULTIGRID */
      if (master) then
         if (have_mpi) then
            if (allow_uneven) call warn("[domain:init_domain] Uneven domain decomposition is experimental.")
            if (allow_noncart) call warn("[domain:init_domain] Non-cartesian domain decomposition is highly experimental.")
         endif
         if (allow_AMR) call warn("[domain:init_domain allow_AMR is not implemented")
      endif
      is_uneven = .false.
      is_mpi_noncart = .false.
      is_refined = .false.
      is_multicg = .false.

      where (.not. dom%has_dir(:)) psize(:) = I_ONE

      ! cdd% will contain valid values if and only if comm3d becomes valid communicator
      cdd%procn(:,:) = MPI_PROC_NULL
      cdd%psize(:) = -I_ONE
      cdd%pcoords(:) = -I_ONE
      cdd%procxyl = MPI_PROC_NULL
      cdd%procyxl = MPI_PROC_NULL
      cdd%comm3d = MPI_COMM_NULL

      dom_divided = .false.
      ! \todo domain division should be converted to type-bound procedure or
      ! at least operate on pdom
      if (associated(user_divide_domain)) call user_divide_domain(dom_divided)
      call is_too_small(dom_divided, "user_divide_domain")
      if (.not. dom_divided) call divide_domain(dom_divided)
      call is_too_small(dom_divided, "divide_domain")
      call MPI_Allreduce(MPI_IN_PLACE, dom_divided, I_ONE, MPI_LOGICAL, MPI_LAND, comm, ierr)
      if (.not. dom_divided .and. master) call die("[domain:init_domain] Domain decomposition failed")

      !\todo Analyze the decomposition and set up [ is_uneven, is_mpi_noncart, is_refined, ... ]
      is_multicg = (ubound(dom%pse(proc)%sel(:, :, :), dim=1) > 1)
      call MPI_Allreduce(MPI_IN_PLACE, is_multicg, I_ONE, MPI_LOGICAL, MPI_LOR, comm, ierr)
      if (is_multicg .and. cdd%comm3d /= MPI_COMM_NULL) call die("[domain:init_domain] is_multicg cannot be used with comm3d")
      if (is_refined) then
         is_mpi_noncart = .true.
         is_multicg = .true.
      endif
      if (is_mpi_noncart) is_uneven = .true.

      ! bnd_[xyz][lr] now become altered according to local topology of processes
      if (is_multicg .and. master) call warn("[domain:init_domain] Multiple blocks per process not fully implemented yet")

!#ifdef VERBOSE
      if (master) then
         !call dom%print_me
         allocate(maxcnt(FIRST:LAST))
         maxcnt(:) = 0
         do p = FIRST, LAST
            i = 1
            write(msg,'(a,i4,a,2(3i6,a),i8,a)') "[domain:init_domain] segment @",p," : [",dom%pse(p)%sel(i, :, LO),"] : [",dom%pse(p)%sel(i, :, HI),"] #", &
                 &                              product(dom%pse(p)%sel(i, :, HI)-dom%pse(p)%sel(i, :, LO)+1)," cells"
            call printinfo(msg)
            if (ubound(dom%pse(p)%sel(:, :, :), dim=1) > 1) then
               do i= 2, ubound(dom%pse(p)%sel(:, :, :), dim=1)
                  write(msg,'(a,2(3i6,a),i8,a)')"                                   : [",dom%pse(p)%sel(i, :, LO),"] : [",dom%pse(p)%sel(i, :, HI),"] #", &
                       &                        product(dom%pse(p)%sel(i, :, HI)-dom%pse(p)%sel(i, :, LO)+1)," cells"
                  call printinfo(msg)
               enddo
            endif
            do i= 1, ubound(dom%pse(p)%sel(:, :, :), dim=1)
               maxcnt(p) = maxcnt(p) + product(real(dom%pse(p)%sel(i, :, HI)-dom%pse(p)%sel(i, :, LO)+1))
            enddo
         enddo
         write(msg,'(a,f8.5)')"[domain:init_domain] Load balance: ",product(real(dom%n_d(:)))/(nproc*maxval(maxcnt(:))) !> \todo add calculation of total internal boundary surface in cells
         call printinfo(msg)
         deallocate(maxcnt)
      endif
!#endif /* VERBOSE */

      if (is_refined) call die("[domain:init_domain] Refinements are not implemented")

   end subroutine init_domain

!>
!! \brief Prevent domain decompositions into pieces that are narrower than number of guardcells
!!
!! \details When a piece of grid is narrower than number of guardcells we may expect the following problems:
!!  * Complicated boundary exchange routines because a single neighbour cannot provide valid boundary data in one step.
!!  * Huge memory overhead because number of guardcells is much larger than number of active cells.
!!  * Huge pergormance penalty because everything becomes dominated by guardcell operations.
!! If this routine prevents you runnina a simulation then probably you either try to use too many processors or you use wrong domain decomposition scheme.
!<

   subroutine is_too_small(dom_divided, label)

      use constants,  only: I_ONE
      use dataio_pub, only: warn, msg
      use mpi,        only: MPI_IN_PLACE, MPI_LOGICAL, MPI_LAND
      use mpisetup,   only: proc, comm, ierr

      implicit none

      logical, intent(inout) :: dom_divided
      character(len=*), intent(in) :: label

      integer :: p

      if (.not. dom_divided) return

      if (allocated(dom%pse)) then
         if (allocated(dom%pse(proc)%sel)) then
            do p = lbound(dom%pse(proc)%sel(:, :, :), dim=1), ubound(dom%pse(proc)%sel(:, :, :), dim=1)
               dom_divided = dom_divided .and. all(dom%pse(proc)%sel(p, :, HI) - dom%pse(proc)%sel(p, :, LO) >= dom%nb - 1 .or. .not. dom%has_dir(:))
               if (.not. all(dom%pse(proc)%sel(p, :, HI) - dom%pse(proc)%sel(p, :, LO) >= dom%nb - 1 .or. .not. dom%has_dir(:))) then
                  write(msg,'(3a,2(3i6,a),i4,a,3l2)')"[domain:is_too_small] ",label," [",dom%pse(proc)%sel(p, :, LO),"]:[",dom%pse(proc)%sel(p, :, HI), "] nb",dom%nb," h_d",dom%has_dir(:)
                  call warn(msg)
               endif
            enddo
         else
            dom_divided = .false.
            write(msg,'(3a)')"[domain:is_too_small] ",label," no dom%pse(proc)%sel"
            call warn(msg)
         endif
      else
         dom_divided = .false.
         write(msg,'(3a)')"[domain:is_too_small] ",label," no dom%pse"
         call warn(msg)
      endif
      call MPI_Allreduce(MPI_IN_PLACE, dom_divided, I_ONE, MPI_LOGICAL, MPI_LAND, comm, ierr)

      if (allocated(dom%pse) .and. .not. dom_divided) call deallocate_pse

   end subroutine is_too_small

!!-----------------------------------------------------------------------------

   subroutine allocate_pse(n_cg)

      use constants,  only: xdim, zdim, LO, HI, I_ONE
      use dataio_pub, only: die
      use mpisetup,   only: FIRST, LAST, nproc

      implicit none

      integer, dimension(nproc), intent(in), optional :: n_cg
      integer :: p

      if (allocated(dom%pse)) call die("[domain:allocate_pse] dom%pse already allocated")
      allocate(dom%pse(FIRST:LAST))
      do p = FIRST, LAST
         if (present(n_cg)) then
            allocate(dom%pse(p)%sel(n_cg(p+1), xdim:zdim, LO:HI))
         else
            allocate(dom%pse(p)%sel(I_ONE, xdim:zdim, LO:HI))
         endif
         dom%pse(p)%sel(:, :, :) = 0
      enddo

   end subroutine allocate_pse

!!-----------------------------------------------------------------------------

   subroutine deallocate_pse

      use mpisetup, only: FIRST, LAST

      implicit none

      integer :: p

      if (allocated(dom%pse)) then
         do p = FIRST, LAST
            deallocate(dom%pse(p)%sel)
         enddo
         deallocate(dom%pse)
      endif

   end subroutine deallocate_pse

!!-----------------------------------------------------------------------------

   subroutine cartesian_tiling(p_size)

      use constants,  only: xdim, ydim, zdim, ndims, LO, HI, BND_COR, I_ONE
      use dataio_pub, only: printinfo, die
      use mpi,        only: MPI_COMM_NULL
      use mpisetup,   only: master, FIRST, LAST, nproc, comm, proc, ierr

      implicit none

      integer(kind=4), dimension(ndims), intent(in) :: p_size

      integer(kind=4) :: p
      integer(kind=4), dimension(ndims) :: pc

      if (is_mpi_noncart) call die("[domain:cartesian_tiling] inconsistent decomposition specs")
      if (product(p_size(:)) /= nproc) call die("[domain:cartesian_tiling] product(p_size(:)) /= nproc")

      call allocate_pse

      if (use_comm3d) then
         cdd%psize(:) = p_size(:)
         if (is_mpi_noncart .or. is_refined) call die("[domain:cartesian_tiling] MPI_Cart_create cannot be used for non-rectilinear or AMR domains")

         call MPI_Cart_create(comm, ndims, cdd%psize, dom%periodic, reorder, cdd%comm3d, ierr)
         call MPI_Cart_coords(cdd%comm3d, proc, ndims, cdd%pcoords, ierr)

         ! Compute neighbors
         do p = xdim, zdim
            call MPI_Cart_shift(cdd%comm3d, p-xdim, I_ONE, cdd%procn(p, LO), cdd%procn(p, HI), ierr)
         enddo

         if (any(dom%bnd(xdim:ydim, LO) == BND_COR)) then
            if (cdd%pcoords(xdim) == 0 .and. cdd%pcoords(ydim) > 0) then
               pc(:) = [ cdd%pcoords(ydim), cdd%pcoords(xdim), cdd%pcoords(zdim) ]
               call MPI_Cart_rank(cdd%comm3d, pc, cdd%procxyl, ierr)
            endif
            if (cdd%pcoords(ydim) == 0 .and. cdd%pcoords(xdim) > 0 ) then
               pc(:) = [ cdd%pcoords(ydim), cdd%pcoords(xdim), cdd%pcoords(zdim) ]
               call MPI_Cart_rank(cdd%comm3d, pc, cdd%procyxl, ierr)
            endif
         endif

         if (any(dom%bnd(xdim:ydim, HI) == BND_COR)) call die("[domain:cartesian_tiling] Corner boundary on the right side not implemented anywhere")

         if (master) call printinfo("[domain:cartesian_tiling] Cartesian decomposition with comm3d")
      else
         if (master) call printinfo("[domain:cartesian_tiling] Cartesian decomposition without comm3d")
      endif

      do p = FIRST, LAST
         if (cdd%comm3d == MPI_COMM_NULL) then
            if (use_comm3d) call die("[domain:cartesian_tiling] MPI_Cart_create failed")
            pc(:) = [ mod(p, p_size(xdim)), mod(p/p_size(xdim), p_size(ydim)), p/product(p_size(xdim:ydim)) ]
         else
            call MPI_Cart_coords(cdd%comm3d, p, ndims, pc, ierr)
         endif
         where (dom%has_dir(:))
            dom%pse(p)%sel(I_ONE, :, LO) = (dom%n_d(:) *  pc(:) ) / p_size(:)     ! offset of low boundaries of the local domain (0 at low external boundaries)
            dom%pse(p)%sel(I_ONE, :, HI) = (dom%n_d(:) * (pc(:)+1))/p_size(:) - 1 ! offset of high boundaries of the local domain (n_d(:) - 1 at right external boundaries)
         endwhere
      enddo

      if (ubound(dom%pse(proc)%sel(:,:,:), dim=1) > 1) call die("[domain:cartesian_tiling] No multiblock support.")

   end subroutine cartesian_tiling

!!-----------------------------------------------------------------------------

   subroutine choppy_tiling(p_size)

      use constants,  only: xdim, ydim, zdim, LO, HI, I_ZERO, I_ONE
      use dataio_pub, only: printinfo, die
      use mpisetup,   only: master, nproc, proc

      implicit none

      integer(kind=4), dimension(ndims), intent(in) :: p_size

      integer(kind=4) :: p, px, py
      integer(kind=4), dimension(:), allocatable :: pz_slab, py_slab

      call allocate_pse

      if (master) call printinfo("[domain:choppy_tiling] Non-cartesian decomposition (no comm3d)")
      allocate(pz_slab(p_size(zdim) + 1))
      pz_slab(1) = I_ZERO
      do p = I_ONE, p_size(zdim)
         pz_slab(p+1) = pz_slab(p) + nproc / p_size(zdim)
         if (p <= mod(nproc, p_size(zdim))) pz_slab(p+1) = pz_slab(p+1) + I_ONE ! longer slabs go first
      enddo
      do p = I_ONE, p_size(zdim)
         do px = pz_slab(p), pz_slab(p+1)-I_ONE
            dom%pse(px)%sel(I_ONE, zdim, LO) = nint((dom%n_d(zdim) *  pz_slab(p)   ) / real(nproc))
            dom%pse(px)%sel(I_ONE, zdim, HI) = nint((dom%n_d(zdim) *  pz_slab(p+1) ) / real(nproc)) - 1
         enddo
         allocate(py_slab(p_size(ydim) + 1))
         py_slab(1) = I_ZERO
         do py = I_ONE, p_size(ydim)
            py_slab(py+1) = py_slab(py) + (pz_slab(p+1)-pz_slab(p)) / p_size(ydim)
            if (py <= mod((pz_slab(p+1)-pz_slab(p)), p_size(ydim))) py_slab(py+1) = py_slab(py+1) + I_ONE ! longer slabs go first
         enddo
         do py = I_ONE, p_size(ydim)
            do px = pz_slab(p)+py_slab(py), pz_slab(p)+py_slab(py+1) - I_ONE
               dom%pse(px)%sel(I_ONE, ydim, LO) = nint((dom%n_d(ydim) *  py_slab(py)   ) / real(pz_slab(p+1)-pz_slab(p)))
               dom%pse(px)%sel(I_ONE, ydim, HI) = nint((dom%n_d(ydim) *  py_slab(py+1) ) / real(pz_slab(p+1)-pz_slab(p))) - I_ONE
            enddo
            do px = I_ZERO, py_slab(py+1)-py_slab(py) - I_ONE
               dom%pse(pz_slab(p)+py_slab(py)+px)%sel(I_ONE, xdim, LO) = (dom%n_d(xdim) *  px    ) / (py_slab(py+1)-py_slab(py))
               dom%pse(pz_slab(p)+py_slab(py)+px)%sel(I_ONE, xdim, HI) = (dom%n_d(xdim) * (px+1) ) / (py_slab(py+1)-py_slab(py)) - 1 ! no need to sort lengths here
            enddo
         enddo
         if (allocated(py_slab)) deallocate(py_slab)
      enddo
      if (allocated(pz_slab)) deallocate(pz_slab)

      if (ubound(dom%pse(proc)%sel(:,:,:), dim=1) > 1) call die("[domain:choppy_tiling] No multiblock support.")

   end subroutine choppy_tiling

!!-----------------------------------------------------------------------------
!!
!> \brief public routine for setting user decompositions
!!
   subroutine set_pse_sel(p, n, se)

      use constants,  only: xdim, zdim, LO, HI

      implicit none

      integer(kind=4), intent(in) :: p !< process
      integer, intent(in) :: n !< block number
      integer, dimension(xdim:zdim, LO:HI), intent(in) :: se !< segment

      dom%pse(p)%sel(n,:,:) = se(:,:)

   end subroutine set_pse_sel

!!-----------------------------------------------------------------------------
!!
!> \brief is_overlap_per checks if two given blocks placed within a periodic domain are overlapping.
!! \details to handle shearing box which is divided in y-direction at the edges, one has to provide another subroutine (is_overlap_per_shear) and add it to interface is_overlap
!<
   logical function is_overlap_per(this, other, periods) result(share)

      use constants,  only: xdim, ydim, zdim, ndims, LO, HI

      implicit none

      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: this, other !< two boxes
      integer(kind=8), dimension(xdim:zdim), intent(in)        :: periods     !< where >0 then the direction is periodic with the given number of cells

      integer :: i, j, k
      integer(kind=8), dimension(xdim:zdim, LO:HI) :: oth

      share = .false.
      do i = -1, 1
         if ((dom%has_dir(xdim) .or. periods(xdim)>0) .or. i==0) then
            do j = -1, 1
               if ((dom%has_dir(ydim) .or. periods(ydim)>0) .or. j==0) then
                  do k = -1, 1
                     if ((dom%has_dir(zdim) .or. periods(zdim)>0) .or. k==0) then
                        oth(:,:) = other(:,:) + reshape([i*periods(xdim), j*periods(ydim), k*periods(zdim), i*periods(xdim), j*periods(ydim), k*periods(zdim)], [ndims, HI])
                        share = share .or. is_overlap_simple(this, oth)
                     endif
                  enddo
               endif
            enddo
         endif
      enddo

   end function is_overlap_per

!!-----------------------------------------------------------------------------
!!
!> \brief is_overlap_simple checks if two given blocks placed within a nonperiodic domain are overlapping.
!! This routine is not supposed to take care of periodic domain - use is_overlap_per when you check overlap for boxes that cross the periodic domain boundary
!<
   logical function is_overlap_simple(this, other) result(share)

      use constants,  only: xdim, zdim, LO, HI

      implicit none

      integer(kind=8), dimension(xdim:zdim, LO:HI), intent(in) :: this, other !< two boxes

      integer :: d

      share = .true.
      do d = xdim, zdim
         if (dom%has_dir(d)) share = share .and. (other(d, LO) <= this(d, HI)) .and. (other(d, HI) >= this(d, LO))
      enddo

   end function is_overlap_simple

!-----------------------------------------------------------------------------

   subroutine cleanup_domain

      use mpi,      only: MPI_COMM_NULL
      use mpisetup, only: ierr

      implicit none

      if (cdd%comm3d /= MPI_COMM_NULL) call MPI_Comm_free(cdd%comm3d, ierr)

      call deallocate_pse
      if (allocated(primes)) deallocate(primes)

   end subroutine cleanup_domain

!-----------------------------------------------------------------------------
!
!> \brief This routine computes optimal allowed domain division
!
!> \todo make this a member of types::domain_container
!
   subroutine divide_domain(dom_divided)

      use dataio_pub, only: die, warn, printinfo, msg
      use mpisetup,   only: nproc, master, have_mpi

      implicit none

      logical, intent(inout) :: dom_divided

      real :: quality
      integer(kind=4), dimension(ndims) :: p_size

      if (dom_divided) then
         call warn("[domain:divide_domain] Domain already decomposed")
         return
      endif

      if (all(bsize(:) > 0 .or. .not. dom%has_dir(:))) then
         call stamp_cg(dom_divided)
         call is_too_small(dom_divided, "stamp_cg")
      endif
      if (dom_divided) return

      if (product(psize(:)) == nproc) then
         if (all(mod(dom%n_d(:), psize(:)) == 0)) then
            dom_divided = .true.
            if (master .and. have_mpi) then
               write(msg,'(a,3i4,a,3i6,a)')"[domain:divide_domain] Domain divided to [",psize(:)," ] pieces, each of [",dom%n_d(:)/psize(:)," ] cells."
               call printinfo(msg)
            endif
            call cartesian_tiling(psize(:))
            call is_too_small(dom_divided, "cartesian_tiling")
            if (dom_divided) return
         else
            write(msg,'(a,3i6,a,3i4,a)')"[domain:divide_domain] Cannot divide domain with [",dom%n_d(:)," ] cells to [",psize(:)," ] piecess. "
            if (master) call warn(msg)
         endif
      endif

      call Eratosthenes_sieve(primes, nproc) ! it is possible to use primes only to sqrt(nproc), but it is easier to have the full table. Cheap for any reasonable nproc.

      ! this is the minimal total area of internal boundaries (periodic case), achievable for some perfect domain divisions
      ideal_bsize = dom%eff_dim * (nproc * product(real(dom%n_d(:)))**(dom%eff_dim-1))**(1./dom%eff_dim)

      call divide_domain_uniform(p_size(:))
      if (product(p_size(:)) == nproc) then
         quality = ideal_bsize / sum(p_size(:)/real(dom%n_d(:)) * product(real(dom%n_d(:))), MASK = dom%n_d(:) > 1)
         if (quality >= dd_unif_quality .or. .not. (allow_uneven .or. allow_noncart)) then
            call cartesian_tiling(p_size(:))
            dom_divided = .true.
            call is_too_small(dom_divided, "divide_domain_uniform")
            if (dom_divided) return
         else
            write(msg,'(2(a,f6.3),a)')"[domain:divide_domain] Quality of uniform division = ",quality," is below threshold ",dd_unif_quality, ", trying harder ..."
            if (master) call warn(msg)
         endif
      endif

      if (allow_uneven) then
         call divide_domain_rectlinear(p_size(:))
         quality = 1 !< \todo make an estimate
         if (product(p_size(:)) == nproc) then
            if (quality > dd_rect_quality .or. .not. allow_noncart) then
               call cartesian_tiling(p_size(:))
               dom_divided = .true.
               call is_too_small(dom_divided, "divide_domain_rectlinear")
               if (dom_divided) return
            endif
         endif
         if (master) call warn("[domain:divide_domain] divide_domain_rectlinear failed")
      else
         if (master) call warn("[domain:divide_domain] Did not try uneven domain division")
      endif

      if (allow_noncart) then
         p_size(:) = psize(:)
         call divide_domain_slices(p_size(:))
         ! if good_enough then return
         call choppy_tiling(p_size(:))
         dom_divided = .true.
         call is_too_small(dom_divided, "divide_domain_slices")
         if (dom_divided) return
         if (master) call warn("[domain:divide_domain] divide_domain_slices failed")
      else
         if (master) call warn("[domain:divide_domain] Did not try non-cartesian domain division")
      endif

      write(msg,'(a,3i6,a,i4,a)') "[domain:divide_domain] Cannot divide domain with [",dom%n_d(:)," ] cells to ",nproc," piecess. "
      if (master) call die(msg)

   end subroutine divide_domain

!-----------------------------------------------------------------------------
!>
!! \brief This routine tries to divide the computational domain into local domains.
!! \details The goal is to minimize the ratio of longest to shortest edge to minimize the amount of inter-process communication.
!! If the benchmarks show that some direction should be partitioned in more pieces than other directions, implement appropriate weighting in j1, j2 and j3 calculations.
!!
!! For some weird domains and PE counts this routine may find tiling that does not satisfy multigrid restrictions even if there exists an acceptable tiling.
!! In such case the user must divide domain manually by providing psize(:) parameters through problem.par.
!!
!! \attention Must be called by all procs to avoid communication and ensure that every proc has proper psize
!<
   subroutine divide_domain_uniform(p_size)

      use constants,  only: xdim, zdim
      use dataio_pub, only: warn, printinfo, msg
      use mpisetup,   only: master, nproc

      implicit none

      integer(kind=4), dimension(ndims), intent(out) :: p_size

      integer(kind=4) :: n
      integer :: j1, j2, j3, jj, p
      integer(kind=4), dimension(ndims) :: ldom, tmp

      ldom(xdim:zdim) = dom%n_d(zdim:xdim:-1) ! Maxloc returns first occurrence of max, reversing direction order (to ZYX) gives better cache utilization.
      n = nproc
      p_size(:) = 1
      if (n == 1) return

      do p = size(primes), 1, -1 ! start from largest defined primes, continue down to 2
         do while (mod(n, primes(p))==0)

            jj = 0
            j1 = sum(maxloc(ldom), 1) ! First try the longest edge; note the trick to make a scalar from 1-element vector without assignment to another variable
            if (mod(ldom(j1), primes(p))==0) then
               jj = j1
            else
               j2 = 1 + mod(j1 + 0, int(ndims))
               j3 = 1 + mod(j1 + ndims -2, int(ndims))
               if (ldom(j2) > ldom(j3)) then
                  j2 = 1 + mod(j1 + ndims -2, int(ndims))
                  j3 = 1 + mod(j1 + 0, int(ndims))
               endif
               if (mod(ldom(j2), primes(p))==0) jj = j2 ! middle edge ...
               if (jj == 0 .and. mod(ldom(j3), primes(p))==0) jj = j3 ! try the shortest edge on last resort
            endif

            if (jj == 0) then
               if (master) call warn("[domain:divide_domain_uniform]: Can't find divisible edge")
               p_size(:) = 1
               return
            else
               p_size(jj) = p_size(jj) * primes(p)
               n          = n          / primes(p)
               ldom(jj)   = ldom(jj)   / primes(p)
            endif

         enddo
      enddo

      if (any(ldom(:) < dom%nb .and. dom%has_dir(:)) .or. n /= 1) then
         if (master) call warn("[domain:divide_domain_uniform]: I am not that intelligent") ! nproc has too big prime factors
         p_size(:) = 1
         return
      endif

      tmp(xdim:zdim) = p_size(zdim:xdim:-1) ! directions were reverted at ldom assignment
      p_size(:) = tmp(:)

      if (master) then
         write(msg,'(a,3i4,a,3i6,a)')"[domain:divide_domain_uniform] Domain divided to [",p_size(:)," ] pieces, each of [",ldom(zdim:xdim:-1)," ] cells."
         call printinfo(msg)
      endif

   end subroutine divide_domain_uniform

!-----------------------------------------------------------------------------
!>
!! \brief Divide the computational domain into local domains. Allow their size to change by +/- 1 depending on CPU rank (this will introduce some load imbalance)
!! if it is not possible to divide an edge evenly. Try to minimize the imbalance and total internal boundaries size.
!<
   subroutine divide_domain_rectlinear(p_size)

      use constants,  only: xdim, ydim, I_ZERO, I_ONE
      use dataio_pub, only: printinfo, msg
      use mpisetup,   only: master, nproc

      implicit none

      integer(kind=4), dimension(ndims), intent(out) :: p_size

      real, parameter :: b_load_fac = 0.25 ! estimated increase of execution time after doubling the total size of internal boundaries.
      ! \todo estimate this factor for massively parallel runs and for Intel processors

      integer(kind=4), allocatable, dimension(:) :: ppow
      integer(kind=4), allocatable, dimension(:,:) :: fac
      integer(kind=4), dimension(ndims) :: ldom
      integer(kind=4) :: p, i, j, k, nf
      integer :: n, ii, bsize
      real :: load_balance, best, quality

      p_size(:) = 1
      if (nproc == 1) return

      allocate(ppow(size(primes)))

      p = nproc
      do i = I_ONE, int(size(primes), kind=4)
         ppow(i) = 0
         do while (mod(p, primes(i)) == 0)
            ppow(i) = ppow(i) + I_ONE
            p = p / primes(i)
         enddo
      enddo

      nf = int(count(ppow(:) > 0), kind=4)
      allocate(fac(nf,3))
      j = I_ONE
      do i = I_ONE, int(size(primes), kind=4)
         if (ppow(i)>0) then
            fac(j,:) = [ primes(i), ppow(i), int((ppow(i)+1)*(ppow(i)+2)/2, kind=4) ] ! prime, its power and number of different decompositions in three dimensions for this prime
            j = j + I_ONE
         endif
      enddo
      deallocate(ppow)

      best = 0.
      ii = 0
      do while (all(fac(:,3) > 0))
         ldom(:) = 1
         do n = 1, nf ! find an unique decomposition of fac(n,2) into [i,j,k], all([i,j,k] >= 0) && i+j+k = fac(n,2). The decompositions are enumerated with fac(n,3).
            i = int(sqrt(1./4.+2.*(fac(n,3)-1)) - 1./2., kind=4) ! i and k enumerate a point in a triangle: (i>=0 && k>=0 && i+k<=fac(n,2))
            k = fac(n,3) - int(1 + i*(i+1)/2, kind=4)
            i = fac(n,2) - i
            j = fac(n,2) - (i + k)
            ldom(:) = ldom(:) * fac(n,1)**[i, j, k]
         enddo

         bsize = int(sum(ldom(:)/dble(dom%n_d(:)) * product(int(dom%n_d(:), kind=8)), MASK = dom%n_d(:) > 1)) !ldom(1)*dom%n_d(2)*dom%n_d(3) + ldom(2)*dom%n_d(1)*dom%n_d(3) + ldom(3)*dom%n_d(1)*dom%n_d(2)
         load_balance = product(real(dom%n_d(:))) / ( real(nproc) * product( int((dom%n_d(:)-1)/ldom(:)) + 1 ) )

         quality = load_balance/ (1 + b_load_fac*(bsize/ideal_bsize - 1.))
         ! \todo add a factor that estimates lower cost when x-direction is not chopped too much
         quality = quality * (1. - (0.001 * ldom(xdim) + 0.0001 * ldom(ydim))/nproc) ! \deprecated estimate these magic numbers

         if (any(ldom(:) > dom%n_d(:))) quality = 0
         if (any(dom%n_d(:)/ldom(:) < dom%nb .and. dom%has_dir(:))) quality = 0

#ifdef DEBUG
         if (quality > 0 .and. master) then
            ii = ii + 1
            write(msg,'(a,i3,a,3i4,a,i10,2(a,f10.7))')"m:ddr ",ii," p_size= [",ldom(:)," ], bcells= ", bsize, ", balance = ", load_balance, ", est_quality = ", quality
            call printinfo(msg)
         endif
#endif /* DEBUG */
         if (quality > best) then
            best = quality
            p_size(:) = ldom(:)
         endif
         do j = I_ONE, nf ! search for next unique combination
            if (fac(j,3) > 1) then
               fac(j,3) = fac(j,3) - I_ONE
               exit
            else
               if (j<nf) then
                  fac(j,3) = int((fac(j,2)+1)*(fac(j,2)+2)/2, kind=4)
               else
                  fac(:,3) = I_ZERO ! no more combinations to try
               endif
            endif
         enddo
      enddo

      deallocate(fac)

      is_uneven = any(mod(dom%n_d(:), p_size(:)) /= 0)

      if (master) then
#ifdef DEBUG
         write(msg,'(a,3f10.2,a,i10)')"m:ddr id p_size = [",(nproc/product(dble(dom%n_d(:))))**(1./dom%eff_dim)*dom%n_d(:),"], bcells= ", int(ideal_bsize)
         call printinfo(msg)
#endif /* DEBUG */
         write(msg,'(a,3i4,a)')      "[domain:divide_domain_rectlinear] Domain divided to [",p_size(:)," ] pieces"
         call printinfo(msg)
         if (is_uneven) then
            write(msg,'(2(a,3i5),a)')"                                    Sizes are from [", int(dom%n_d(:)/p_size(:))," ] to [",int((dom%n_d(:)-1)/p_size(:))+1," ] cells."
            call printinfo(msg)
            write(msg,'(a,f8.5)')    "                                    Load balance is ",product(int(dom%n_d(:), kind=8)) / ( dble(nproc) * product( int((dom%n_d(:)-1)/p_size(:)) + 1 ) )
         else
            write(msg,'(a,3i5,a)')   "                                    Size is [", int(dom%n_d(:)/p_size(:))," ] cells."
         endif
         call printinfo(msg)
      endif

   end subroutine divide_domain_rectlinear

!-----------------------------------------------------------------------------
!>
!! \brief Divide the computational domain into local domains. Allow their size to depend significantly on CPU rank and allow for more than one neighbour on a single boundary.
!! Try to minimize the imbalance and total internal boundaries size.
!<
   subroutine divide_domain_slices(p_size)

      use constants,  only: xdim, ydim, zdim, I_ONE
      use dataio_pub, only: msg, printinfo, warn
      use mpisetup,   only: master, nproc

      implicit none

      integer(kind=4), dimension(ndims), intent(inout) :: p_size

      real, parameter :: minfac = 1.3 ! prevent domain division to halves if cell count in a given direction is too low. (not verified for optimality)
      real :: optc

      !> \todo Try to make an intelligent guess for slicing, then go down to the local minimum and explore neighbourhood. Exploring all possibilities is an O(nproc)**2 task
      ! The best solution is probably near (nproc/product(dble(dom%n_d(:))))**(1./dom%eff_dim)*dom%n_d(:)

      is_mpi_noncart = .true.
      is_uneven = .true.

      if (all(p_size(ydim:zdim) == 1)) then
         if (dom%has_dir(zdim)) then
            optc = max(real(dom%nb), (product(int(dom%n_d(:), kind=8))/real(nproc)) ** (1./dom%eff_dim)) ! number of cells for ideal cubes
            if (dom%n_d(zdim) > minfac*optc) p_size(zdim) = int(ceiling(dom%n_d(zdim)/optc), kind=4)
         endif
         if (dom%has_dir(ydim)) then
            optc = max(real(dom%nb), (product(int(dom%n_d(xdim:ydim), kind=8))*p_size(zdim)/real(nproc)) ** (1./count(dom%has_dir(xdim:ydim))))
            if (dom%n_d(ydim) > minfac*optc) p_size(ydim) = int(ceiling(dom%n_d(ydim)/optc), kind=4)
         endif
      endif
      if (dom%has_dir(xdim)) p_size(xdim) = (nproc - I_ONE)/(p_size(ydim)*p_size(zdim)) + I_ONE !sometimes it might be less by 1

      where (.not. dom%has_dir(:)) p_size(:) = 1 ! just in case
      do while (product(p_size(:)) < nproc)
         write(msg,'(a,3i4,a)') "[domain:divide_domain_slices] imperfect noncartesian division to [",p_size(:)," ] pieces"
         if (master) call warn(msg)
         if (dom%has_dir(xdim)) then
            p_size(xdim) = p_size(xdim) + I_ONE
         else if (dom%has_dir(ydim)) then
            p_size(ydim) = p_size(ydim) + I_ONE
         else
            p_size(zdim) = p_size(zdim) + I_ONE
         endif
      enddo
      write(msg,'(a,3i4,a)') "[domain:divide_domain_slices] performed noncartesian division to [",p_size(:)," ] pieces"
      if (master) call printinfo(msg)

   end subroutine divide_domain_slices

!-----------------------------------------------------------------------------
!>
!! \brief Divide the domain into lots of identical blocks
!<
   subroutine stamp_cg(dom_divided)

      use constants,    only: xdim, zdim, LO, HI
      use dataio_pub,   only: warn, msg, die
      use mpisetup,     only: master, proc, nproc

      implicit none

      logical, intent(inout) :: dom_divided
      integer, dimension(xdim:zdim) :: n_bl, b_loc
      integer, dimension(xdim:zdim, LO:HI) :: se
      integer :: tot_bl, loc_bl, bl_s, bl_e
      integer, dimension(nproc) :: n_cg ! BEWARE: antiparallel
      integer, allocatable, dimension(:) :: pb ! BEWARE: antiparallel
      integer :: p, b

      dom_divided = .false.

      if (any(bsize(xdim:zdim) <=0)) then
         if (master) call warn("[initproblem:stamp_cg] some(bsize(1:3)) <=0")
         return
      endif

      if (any(mod(dom%n_d(:), bsize(xdim:zdim)) /= 0 .and. dom%has_dir(:))) then
         write(msg,'(a,3f10.3,a)')"[initproblem:stamp_cg] Fractional number of blocks: dom%n_d(:)/bsize(1:3) = [",dom%n_d(:)/real(bsize(xdim:zdim)),"]"
         if (master) call warn(msg)
         return
      endif

      where (dom%has_dir(:))
         n_bl(:) = dom%n_d(:) / bsize(xdim:zdim)
      elsewhere
         n_bl(:) = 1
      endwhere
      tot_bl = product(n_bl(:), mask=dom%has_dir(:))
      if (allocated(pb)) call die("[initproblem:stamp_cg] pb already allocated")
      allocate(pb(tot_bl))

      if (tot_bl < nproc) then
         if (master) call warn("[initproblem:stamp_cg] tot_bl < nproc")
         return
      endif

      do p = 0, nproc-1
         n_cg(p+1) = (((p + 1) * tot_bl)/nproc) - ((p * tot_bl)/nproc) ! number of blocks on process p
      enddo
      call allocate_pse(n_cg(:))

      bl_s = (proc * tot_bl)/nproc + 1
      bl_e = ((proc + 1) * tot_bl)/nproc
      loc_bl = bl_e - bl_s + 1
      if (loc_bl /= n_cg(proc + 1)) call die("[initproblem:stamp_cg] loc_bl /= n_cg(proc + 1)")

      do p = 0, nproc-1
         n_cg(p+1) = (p * tot_bl)/nproc + 1 ! first block on process p
         if (p>0) pb(n_cg(p):n_cg(p+1)-1) = p-1
      enddo
      pb(n_cg(nproc):) = nproc-1

      do p = 0, nproc-1
         if (size(dom%pse(p)%sel(:,:,:), dim=1) /= count(pb(:) == p)) call die("[initproblem:stamp_cg] size(dom%pse(p)%sel(:,:,:), dim=1) /= count(pb(:) == p)")
      enddo

      se(:,:) = 0
      do b = 1, tot_bl
         ! pb(b) = min((b*nproc)/tot_bl, nproc-1)
         call simple_ordering(b, n_bl(:), b_loc(:)) !> \todo implement Morton and Hilbert ordering

         where (dom%has_dir(:))
            se(:, LO) = b_loc(:) * bsize(xdim:zdim)
            se(:, HI) = se(:, LO) + bsize(xdim:zdim) - 1
         endwhere
         dom%pse(pb(b))%sel(b-n_cg(pb(b)+1)+1,:,:) = se(:,:) ! equivalent to call set_pse_sel(pb(b), b-n_cg(pb(b)+1)+1, se)
      enddo
      dom_divided = .true.
      !> \todo implement merging to larger cuboids

      if (master) call warn("[initproblem:stamp_cg] Experimental implementation")
      if (.not. dom_divided) call deallocate_pse
      deallocate(pb)

   end subroutine stamp_cg

!-----------------------------------------------------------------------------
!>
!! \brief Simple block ordering
!<

   subroutine simple_ordering(b, n_bl, b_loc)

      use constants, only: xdim, ydim, zdim

      implicit none

      integer, intent(in) :: b                            !< block number
      integer, dimension(xdim:zdim), intent(in)  :: n_bl  !< extents of block array
      integer, dimension(xdim:zdim), intent(out) :: b_loc !< location of block b

      b_loc(:) = [ mod(b-1, n_bl(xdim)), mod((b-1)/n_bl(xdim), n_bl(ydim)), (b-1)/product(n_bl(xdim:ydim)) ]

   end subroutine simple_ordering

!-----------------------------------------------------------------------------
!>
!! \brief Create an array of prime numbers smaller than given integer
!<
   subroutine Eratosthenes_sieve(tab, n)

      use constants,  only: I_ZERO, I_TWO
      use dataio_pub, only: die
#ifdef DEBUG
      use dataio_pub, only: msg, printinfo
      use mpisetup,   only: master
#endif /* DEBUG */

      implicit none

      integer(kind=4), intent(inout), allocatable, dimension(:) :: tab !< array of found prime numbers
      integer(kind=4), intent(in) :: n                                 !< max value of prime number

      integer(kind=4), dimension(n) :: numb
      integer(kind=4) :: i
      integer :: no_primes

      if (allocated(tab)) call die("[domain:Eratosthenes_sieve] tab already allocated")

      numb = [I_ZERO, (i, i = I_TWO, n)]

      do i = 2, n
         if (numb(i) /= 0) numb( 2*i : n : i ) = 0
      enddo

      no_primes = count(numb /= 0)
      allocate(tab(no_primes))
      tab = pack(numb, numb /= 0)
#ifdef DEBUG
      if (master) then
         write(msg,'(2(A,I5))') "There are ", no_primes, " prime numbers less than", n
         call printinfo(msg)
         print "(15i7)", tab
      endif
#endif /* DEBUG */
   end subroutine Eratosthenes_sieve

!-----------------------------------------------------------------------------
!>
!! \brief An interpreter of string-defined boundary types
!!
!! \todo Write a routine that does the reverse transformation - can be useful for easy-to-interpret messages
!<

   subroutine translate_bnds_to_ints(this, bnds)

      use constants, only: xdim, zdim, ndims, LO, HI, &
         &                 BND_MPI, BND_PER, BND_REF, BND_OUT, BND_OUTD, BND_OUTH, BND_COR, BND_SHE, BND_USER, BND_INF, BND_INVALID

      implicit none

      class(domain_container), intent(inout) :: this
      character(len=*), dimension(HI*ndims), intent(in) :: bnds  !< Six strings, describing boundary conditions

      integer :: d, lh

      do d = xdim, zdim
         do lh = LO, HI
            select case (bnds(HI*(d-xdim)+lh))
               case ('per', 'periodic')
                  this%bnd(d, lh) = BND_PER
               case ('ref', 'refl', 'reflecting')
                  this%bnd(d, lh) = BND_REF
               case ('out', 'free')
                  this%bnd(d, lh) = BND_OUT
               case ('outd', 'diode')
                  this%bnd(d, lh) = BND_OUTD
               case ('outh')
                  this%bnd(d, lh) = BND_OUTH
               case ('she', 'shear', 'shearing')
                  this%bnd(d, lh) = BND_SHE
               case ('cor', 'corner')
                  this%bnd(d, lh) = BND_COR
               case ('mpi')
                  this%bnd(d, lh) = BND_MPI
               case ('user')
                  this%bnd(d, lh) = BND_USER
               case ('inf') ! what is this?
                  this%bnd(d, lh) = BND_INF
               case default
                  this%bnd(d, lh) = BND_INVALID
            end select
         enddo
      enddo
   end subroutine translate_bnds_to_ints

!-----------------------------------------------------------------------------

   subroutine set_derived(this)

      use constants,  only: xdim, ydim, zdim, LO, HI, BND_PER, BND_SHE, I_TWO
      use dataio_pub, only: die

      implicit none

      class(domain_container), intent(inout) :: this
      integer :: d

      this%has_dir(:) = this%n_d(:) > 1  ! redundant

      this%periodic(:) = .false.
      do d = xdim, zdim
         if ((any(this%bnd(d, :) == BND_PER) .or. (d==xdim .and. any(this%bnd(d, :) == BND_SHE))) .and. this%has_dir(d)) then
            this%periodic(d) = .true.
            if (this%bnd(d, LO) /= this%bnd(d, HI)) call die("[types:set_derived] Periodic BC do not match")
         endif
      enddo
      if (any(this%bnd(ydim:zdim, :) == BND_SHE)) call die("[types:set_derived] Shearing BC not allowed for y- and z-direction")

      ! auxiliary lengths
      this%L_(:) = this%edge(:, HI) - this%edge(:, LO)
      this%C_(:) = (this%edge(:, HI) + this%edge(:, LO))/2

      !volume and total grid sizes
      this%Vol = product(this%L_(:), mask=this%has_dir(:))
      !> \deprecated BEWARE: Vol computed above is not true for non-cartesian geometry

      where (this%has_dir(:))
         this%n_t(:) = this%n_d(:) + I_TWO * this%nb
      elsewhere
         this%n_t(:) = 1
      endwhere

   end subroutine set_derived

!-----------------------------------------------------------------------------

   subroutine print_me(this)
      use constants,  only: LO, HI
      use dataio_pub, only: printinfo, msg

      implicit none

      class(domain_container), intent(inout) :: this
      integer :: i, j

      write(msg,'(a,3(F5.1,1X))') "LO edge: ", this%edge(:,LO); call printinfo(msg)
      write(msg,'(a,3(F5.1,1X))') "HI edge: ", this%edge(:,HI); call printinfo(msg)
      write(msg,'(a,3(I4,1X))')   "n_d    : ", this%n_d(:); call printinfo(msg)
      write(msg,'(a,3(I4,1X))')   "n_t    : ", this%n_t(:); call printinfo(msg)
      write(msg,'(a,1(I4,1X))')   "nb     : ", this%nb; call printinfo(msg)
      write(msg,'(a,3(I4,1X))')   "LO bnd : ", this%bnd(:,LO); call printinfo(msg)
      write(msg,'(a,3(I4,1X))')   "HI bnd : ", this%bnd(:,HI); call printinfo(msg)
      write(msg,'(a,3(F5.1,1X))') "L_     : ", this%L_(:); call printinfo(msg)
      write(msg,'(a,3(F5.1,1X))') "C_     : ", this%C_(:); call printinfo(msg)
      write(msg,'(a,1(F5.1,1X))') "Vol    : ", this%Vol; call printinfo(msg)
      write(msg,'(a,3(L1,1X))')   "period : ", this%periodic(:); call printinfo(msg)
      write(msg,'(a,I4)')         "size(pse) = ", size(this%pse)

      if (allocated(this%pse)) then
         do i = lbound(this%pse,1), ubound(this%pse,1)
            if (allocated(this%pse(i)%sel)) then
               write(msg,'(a,I4)')         "size(pse%sel,1) = ", size(this%pse(i)%sel,dim=1)
               do j = lbound(this%pse(i)%sel,dim=1), ubound(this%pse(i)%sel,dim=1)
                  write(msg, '(a,3(I6,1X))')   "LO sel : ", this%pse(i)%sel(j,:,LO)
                  write(msg, '(a,3(I6,1X))')   "HI sel : ", this%pse(i)%sel(j,:,HI)
               enddo
            endif
         enddo
      endif
   end subroutine print_me

   subroutine init_domain_container(this, nb, n_d, bnds, edges, geometry)
      use constants,    only: ndims, LO, HI, big_float, dpi, xdim, ydim, zdim, &
           &                GEO_XYZ, GEO_RPZ, GEO_INVALID, BND_PER, BND_REF, BND, I_ONE
      use dataio_pub,   only: die, warn, msg

      implicit none

      class(domain_container), intent(inout) :: this
      integer(kind=4), intent(in) :: nb    !< number of boundary cells surrounding the physical domain, same for all directions
      integer(kind=4), dimension(ndims) :: n_d !< number of %grid cells in physical domain without boundary cells (where  == 1 then that dimension is reduced to a point with no boundary cells)
      character(len=*), dimension(HI*ndims), intent(in) :: bnds  !< Six strings, describing boundary conditions
      real, dimension(ndims, LO:HI), intent(inout)      :: edges  !< physical domain boundaries position
      character(len=*), intent(in) :: geometry !< define system of coordinates: "cartesian" or "cylindrical"

      real :: xmno, ymno, ymxo

      ! Sanitize input parameters, if possible
      this%n_d(:) = max(I_ONE, n_d(:))
      this%has_dir(:) = this%n_d(:) > 1
      this%eff_dim = count(this%has_dir(:))
      where (this%has_dir(:))
         this%D_(:) = 1
      elsewhere
         this%D_(:) = 0
      endwhere

      ! shortcuts
      this%D_x = this%D_(xdim)
      this%D_y = this%D_(ydim)
      this%D_z = this%D_(zdim)

      this%total_ncells = product(int(this%n_d(:), kind=8))
      if (any(this%total_ncells < this%n_d(:))) call die("[domain:init_domain] Integer overflow: too many cells")

      select case (geometry)
         case ("cartesian", "cart", "xyz", "XYZ")
            this%geometry_type = GEO_XYZ
         case ("cylindrical", "cyl", "rpz", "RPZ")
            this%geometry_type = GEO_RPZ
         case default
            this%geometry_type = GEO_INVALID
      end select

      call this%translate_bnds_to_ints(bnds)

      ! sanitize domain
      xmno = edges(xdim,LO)
      ymno = edges(ydim,LO)
      ymxo = edges(ydim,HI)
      select case (this%geometry_type)
         case (GEO_XYZ)
            if (edges(ydim,LO) <= -big_float .or. edges(ydim,HI) >= big_float) call warn("[domain:init_domain] y range not specified. Defaulting to [0..1]")
            if (edges(ydim,LO) <= -big_float) edges(ydim,LO) = 0.
            if (edges(ydim,HI) >= big_float) edges(ydim,HI) = 1.
         case (GEO_RPZ)
            if (edges(ydim,LO) <= -big_float) edges(ydim,LO) = 0.
            if (edges(ydim,HI) >= big_float) edges(ydim,HI) = dpi
            if (edges(ydim,HI)-edges(ydim,LO) > dpi) then
               call warn("[domain:init_domain] Hyperbolic spaces are not implemented. Setting azimuthal span to 2pi.")
               if (abs(edges(ydim,LO)) < 1./epsilon(1.)) then
                  edges(ydim,HI) = edges(ydim,LO) + dpi
               else
                  edges(ydim,LO) = edges(ydim,HI) - dpi
               endif
               if (abs(edges(ydim,HI)-edges(ydim,LO) - dpi) > 100*epsilon(1.)) call die("[domain:init_domain] absolute values for both edges(ydim,HI) and edges(ydim,LO) too high.") ! magic number
            endif
            if (edges(xdim,LO) <= 0.) then
               edges(xdim,LO) = 0.
               if (this%bnd(xdim, LO) /= BND_REF) call warn("[domain:init_domain] Enforcing this%bnd(xdim, LO) = 'ref'.")
               this%bnd(xdim, LO) = BND_REF
            endif
            if (this%bnd(xdim, HI) == BND_PER) call die("[domain:init_domain] Periodicity in radial direction is not allowed in cylindrical coordinates")
         case default
            call die("[domain:init_domain] Invalid geometry type.")
      end select
      if (xmno /= edges(xdim,LO)) then
         write(msg,'(2(a,g20.12))')"[domain:init_domain] Sanitized edges(xdim,LO): ",xmno," -> ",edges(xdim,LO)
         call warn(msg)
      endif
      if (ymno /= edges(ydim,LO) .and. ymno /= -big_float) then
         write(msg,'(2(a,g20.12))')"[domain:init_domain] Sanitized edges(ydim,LO): ",ymno," -> ",edges(ydim,LO)
         call warn(msg)
      endif
      if (ymxo /= edges(ydim,HI) .and. ymxo /= big_float) then
         write(msg,'(2(a,g20.12))')"[domain:init_domain] Sanitized edges(ydim,HI): ",ymno," -> ",edges(ydim,HI)
         call warn(msg)
      endif
      if (edges(xdim,LO) > edges(xdim,HI)) call die("[[domain:init_domain] Negative span in X-direction")
      if (edges(ydim,LO) > edges(ydim,HI)) call die("[[domain:init_domain] Negative span in Y-direction")
      if (edges(zdim,LO) > edges(zdim,HI)) call die("[[domain:init_domain] Negative span in Z-direction")
      if (nb < 1) call die("[[domain:init_domain] no guardcells")

      ! set up the global domain
      this%nb = nb

      this%edge(:,:) = edges(:,:)

      call this%set_derived ! finish up with the rest of domain_container members

   end subroutine init_domain_container

end module domain
