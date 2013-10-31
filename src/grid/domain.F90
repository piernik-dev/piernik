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
!! \brief This module contains everything closely related to management of computational domain.
!!
!! In this module following namelists of parameters are specified:
!! \copydetails domain::init_domain
!<
module domain

   use constants, only: ndims, cbuff_len, LO, HI

   implicit none

   private
   public :: cleanup_domain, init_domain, translate_ints_to_bnds, domain_container, dom, is_uneven, is_mpi_noncart, is_refined, is_multicg, &
        &    psize, AMR_bsize, minsize, allow_noncart, allow_uneven, dd_unif_quality, dd_rect_quality! temporary export

! AMR: There will be at least one domain container for the base grid.
!      It will be possible to host one or more refined domains on the base container and on the refined containers.
!      The refined domains may cover whole parent domain or only a tiny part of it.
!      The multigrid solver operates also on a stack of coarser domains - parents of the base domain.
!      The coarser domains must be no smaller than the base domain.

   type :: domain_container
      ! primary parameters, read from /DOMAIN_SIZES/, /BOUNDARIES/ and /DOMAIN_LIMITS/ namelists
      real,            dimension(ndims, LO:HI) :: edge !< physical domain boundary positions
      integer(kind=4), dimension(ndims)        :: n_d  !< number of grid cells in physical domain in x-, y- and z-direction (where equal to 1, the dimension is reduced to a point with no boundary cells)
      integer(kind=4)                          :: nb   !< number of boundary cells surrounding the physical domain, same for all directions
      integer(kind=4), dimension(ndims, LO:HI) :: bnd  !< type of boundary conditions coded in integers

      ! derived parameters
      real, dimension(ndims)    :: L_           !< span of the physical domain [ xmax-xmin, ymax-ymin, zmax-zmin ]
      real, dimension(ndims)    :: C_           !< center of the physical domain [ (xmax+xmin)/2., (ymax+ymin)/2., (zmax+zmin)/2. ]
      real                      :: Vol          !< total volume of the physical domain

      logical, dimension(ndims) :: periodic     !< .true. for periodic and shearing boundary pairs

      ! Do not use n_t(:) in the Piernik source tree without a good reason
      ! Avoid as a plague allocating buffers of that size because it negates benefits of parallelization
      ! \todo move them to another type, that extends domain_container?
      integer(kind=4), dimension(ndims) :: n_t  !< total number of %grid cells in the whole domain in every direction (n_d(:) + 2* nb for existing directions)

      integer(kind=8) :: total_ncells           !< total number of %grid cells
      integer         :: geometry_type          !< the type of geometry: cartesian: GEO_XYZ, cylindrical: GEO_RPZ, other: GEO_INVALID
      integer         :: D_x                    !< set to 1 when x-direction exists, 0 otherwise
      integer         :: D_y                    !< set to 1 when y-direction exists, 0 otherwise.
      integer         :: D_z                    !< set to 1 when z-direction exists, 0 otherwise.

      integer(kind=4), dimension(ndims) :: D_   !< set to 1 for existing directions, 0 otherwise. Useful for dimensionally-safe indices for difference operators on arrays,

      logical, dimension(ndims) :: has_dir      !< .true. for existing directions
      integer                   :: eff_dim      !< effective dimensionality of the simulation
      integer(kind=8), dimension(ndims) :: off  !< offset of the base level

    contains

      procedure :: translate_bnds_to_ints     !< Convert strings to integer-coded boundary types
      procedure :: print_me                   !< Print computational domain details
      procedure :: init                       !< Initialize all variables of domain_container type
      procedure :: modify_side                !< Required when base level changes

   end type domain_container

   type(domain_container), protected :: dom !< complete description of base level domain

   ! This set of flags can be useful to gently disable modules that does not support certain grid improvements yet
   !> \todo Get rid of them as soon as all important modules are updated
   logical :: is_uneven          !< .true. when n_b(:) depend on process rank \todo protect it
   logical :: is_mpi_noncart     !< .true. when there exist a process that has more than one neighbour in any direction \todo protect it
   logical :: is_refined         !< .true. when AMR or static refinement is employed \todo protect it
   logical :: is_multicg         !< .true. when at leas one process has more than one grid container \todo protect it

   ! Namelist variables

   integer(kind=4), dimension(ndims) :: psize     !< desired number of MPI blocks in x, y and z-dimension
   integer(kind=4), dimension(ndims) :: AMR_bsize !< the size of cg for multiblock decomposition
   integer(kind=4), dimension(ndims) :: minsize   !< minimum size of cg, default is dom$nb
   !! \todo Implement maximum size of a cg (in cells) for use with GPGPU kernels. The minimum size id nb**dom%eff_dim

   logical :: allow_uneven            !< allows different values of n_b(:) on different processes
   logical :: allow_noncart           !< allows more than one neighbour on a boundary
   real    :: dd_unif_quality         !< uniform domain decomposition may be rejected it its quality is below this threshold (e.g. very elongated local domains are found)
   real    :: dd_rect_quality         !< rectilinear domain decomposition may be rejected it its quality is below this threshold (not used yet)

   namelist /MPI_BLOCKS/ psize, AMR_bsize, minsize, allow_uneven, allow_noncart, dd_unif_quality, dd_rect_quality

   integer(kind=4), dimension(ndims) :: n_d               !< number of %grid cells in physical domain without boundary cells (where  == 1 then that dimension is reduced to a point with no boundary cells)
   integer(kind=4), protected        :: nb                !< number of boundary cells surrounding the physical domain, same for all directions
   character(len=cbuff_len)          :: bnd_xl            !< type of boundary conditions for the left  x-boundary
   character(len=cbuff_len)          :: bnd_xr            !< type of boundary conditions for the right x-boundary
   character(len=cbuff_len)          :: bnd_yl            !< type of boundary conditions for the left  y-boundary
   character(len=cbuff_len)          :: bnd_yr            !< type of boundary conditions for the right y-boundary
   character(len=cbuff_len)          :: bnd_zl            !< type of boundary conditions for the left  z-boundary
   character(len=cbuff_len)          :: bnd_zr            !< type of boundary conditions for the right z-boundary
   real                              :: xmin              !< physical domain left x-boundary position
   real                              :: xmax              !< physical domain right x-boundary position
   real                              :: ymin              !< physical domain left y-boundary position
   real                              :: ymax              !< physical domain right y-boundary position
   real                              :: zmin              !< physical domain left z-boundary position
   real                              :: zmax              !< physical domain right z-boundary position
   character(len=cbuff_len)          :: geometry          !< define system of coordinates: "cartesian" or "cylindrical"
   integer(kind=8), dimension(ndims) :: dom_off           !< offset of the base level

   namelist /BASE_DOMAIN/ n_d, nb, bnd_xl, bnd_xr, bnd_yl, bnd_yr, bnd_zl, bnd_zr, xmin, xmax, ymin, ymax, zmin, zmax, geometry, dom_off
   !namelist /BASE_DOMAIN/ dom, geometry, nb

contains

!-----------------------------------------------------------------------------
!>
!! \brief Routine to set up computational domain and define its decomposition
!!
!! \n \n
!! @b BASE_DOMAIN
!! \n \n
!! <table border="+1">
!!   <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
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
!!   <tr><td>dom_off </td><td>0, 0, 0    </td><td>integer                                   </td><td>\copydoc domain::dom_off </td></tr>
!! </table>
!! \n \n
!! @b MPI_BLOCKS
!! \n \n
!! <table border="+1">
!!   <tr><td width="150pt"><b>parameter</b></td><td width="135pt"><b>default value</b></td><td width="200pt"><b>possible values</b></td><td width="315pt"> <b>description</b></td></tr>
!!   <tr><td>psize(3)       </td><td>1      </td><td>integer</td><td>\copydoc domain::psize          </td></tr>
!!   <tr><td>AMR_bsize(3)   </td><td>0      </td><td>integer</td><td>\copydoc domain::AMR_bsize      </td></tr>
!!   <tr><td>minsize(3)     </td><td>domain::nb </td><td>integer</td><td>\copydoc domain::minsize        </td></tr>
!!   <tr><td>allow_uneven   </td><td>.false.</td><td>logical</td><td>\copydoc domain::allow_uneven   </td></tr>
!!   <tr><td>allow_noncart  </td><td>.false.</td><td>logical</td><td>\copydoc domain::allow_noncart  </td></tr>
!!   <tr><td>dd_unif_quality</td><td>0.9    </td><td>real   </td><td>\copydoc domain::dd_unif_quality</td></tr>
!!   <tr><td>dd_rect_quality</td><td>0.9    </td><td>real   </td><td>\copydoc domain::dd_rect_quality</td></tr>
!! </table>
!! \n \n
!<
   subroutine init_domain

      use constants,  only: xdim, zdim, ndims, LO, HI, PIERNIK_INIT_MPI, I_ONE, I_ZERO, INVALID
      use dataio_pub, only: die, warn, code_progress
      use dataio_pub, only: nh  ! QA_WARN required for diff_nml
      use mpisetup,   only: cbuff, ibuff, lbuff, rbuff, master, slave, piernik_MPI_Bcast, have_mpi

      implicit none

      character(len=cbuff_len), dimension(HI*ndims) :: bnds  !< Six strings, describing boundary conditions
      real, dimension(ndims, LO:HI)     :: edges

      if (code_progress < PIERNIK_INIT_MPI) call die("[domain:init_domain] MPI not initialized.")

      ! Begin processing of namelist parameters

      psize(:)     = I_ONE
      AMR_bsize(:) = I_ZERO
      minsize(:)   = INVALID
      dom_off(:)   = I_ZERO

      n_d(:)   = I_ONE
      nb       = 4
      xmin     = 0.; xmax = 1.
      ymin     = -huge(1.); ymax = huge(1.) ! required for detection of problems in cylindrical coordinates
      zmin     = 0.; zmax = 1.
      geometry = "cartesian"

      allow_uneven  = .true.
      allow_noncart = .false.     !< experimental implementation

      bnd_xl = 'per'
      bnd_xr = 'per'
      bnd_yl = 'per'
      bnd_yr = 'per'
      bnd_zl = 'per'
      bnd_zr = 'per'

      dd_unif_quality = 0.9
      dd_rect_quality = 0.9

      if (master) then
         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=MPI_BLOCKS)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=MPI_BLOCKS, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "MPI_BLOCKS")
         read(nh%cmdl_nml,nml=MPI_BLOCKS, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "MPI_BLOCKS", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=MPI_BLOCKS)
         close(nh%lun)
         call nh%compare_namelist()
         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=BASE_DOMAIN)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=BASE_DOMAIN, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "BASE_DOMAIN")
         read(nh%cmdl_nml,nml=BASE_DOMAIN, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "BASE_DOMAIN", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=BASE_DOMAIN)
         close(nh%lun)
         call nh%compare_namelist()

         if (any(AMR_bsize(:) > 0 .and. AMR_bsize(:) < nb .and. n_d(:) > 1)) call die("[domain:init_domain] AMR_bsize(:) is too small.")

         cbuff(1) = bnd_xl
         cbuff(2) = bnd_xr
         cbuff(3) = bnd_yl
         cbuff(4) = bnd_yr
         cbuff(5) = bnd_zl
         cbuff(6) = bnd_zr
         cbuff(7) = geometry

         ibuff(         xdim:zdim) = psize(:)
         ibuff(  zdim+xdim:2*zdim) = n_d(:)
         ibuff(2*zdim+xdim:3*zdim) = AMR_bsize(:)
         ibuff(3*zdim+xdim:4*zdim) = minsize(:)
         ibuff(4*zdim+xdim:5*zdim) = int(dom_off(:), kind=4)
         ibuff(5*zdim+1)           = nb

         rbuff(1) = xmin
         rbuff(2) = xmax
         rbuff(3) = ymin
         rbuff(4) = ymax
         rbuff(5) = zmin
         rbuff(6) = zmax
         rbuff(7) = dd_unif_quality
         rbuff(8) = dd_rect_quality

         lbuff(1) = allow_uneven
         lbuff(2) = allow_noncart

      endif

      call piernik_MPI_Bcast(cbuff, cbuff_len)
      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(rbuff)
      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         allow_uneven  = lbuff(1)
         allow_noncart = lbuff(2)

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

         psize(:)     = int(ibuff(         xdim:zdim), kind=4)
         n_d(:)       = int(ibuff(  zdim+xdim:2*zdim), kind=4)
         AMR_bsize(:) = int(ibuff(2*zdim+xdim:3*zdim), kind=4)
         minsize(:)   = int(ibuff(3*zdim+xdim:4*zdim), kind=4)
         dom_off(:)   = int(ibuff(4*zdim+xdim:5*zdim), kind=4)
         nb           = int(ibuff(5*zdim+1),           kind=4)

      endif

      bnds = [bnd_xl, bnd_xr, bnd_yl, bnd_yr, bnd_zl, bnd_zr]
      edges = reshape( [xmin, ymin, zmin, xmax, ymax, zmax], shape=[ndims,HI-LO+I_ONE] )

      call dom%init(nb, n_d, bnds, edges, geometry, dom_off)

      where (dom%has_dir(:))
         minsize(:) = max(minsize(:), dom%nb)
      elsewhere
         minsize(:) = 1
      endwhere

      if (have_mpi .and. master) then
         if (allow_uneven) call warn("[domain:init_domain] Uneven domain decomposition is experimental.")
         if (allow_noncart) call warn("[domain:init_domain] Non-cartesian domain decomposition is experimental.")
      endif
      is_uneven = .false.
      is_mpi_noncart = .false.
      is_refined = .false.
      is_multicg = .false.

      where (.not. dom%has_dir(:)) psize(:) = I_ONE

   end subroutine init_domain

!> \brief Cleanup. Nothing special to do here at the moment.

   subroutine cleanup_domain

      implicit none

   end subroutine cleanup_domain

!> \brief An interpreter of string-defined boundary types

   subroutine translate_bnds_to_ints(this, bnds)

      use constants, only: xdim, zdim, ndims, LO, HI, BND_PER, BND_REF, BND_OUT, BND_OUTD, BND_OUTH, BND_OUTHD, BND_COR, BND_SHE, BND_USER, BND_INVALID

      implicit none

      class(domain_container),               intent(inout) :: this  !< object invoking type-bound procedure
      character(len=*), dimension(HI*ndims), intent(in)    :: bnds  !< Six strings, describing boundary conditions

      integer                                              :: d, lh

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
               case ('outhd')
                  this%bnd(d, lh) = BND_OUTHD
               case ('she', 'shear', 'shearing')
                  this%bnd(d, lh) = BND_SHE
               case ('cor', 'corner')
                  this%bnd(d, lh) = BND_COR
               case ('user')
                  this%bnd(d, lh) = BND_USER
               case default
                  this%bnd(d, lh) = BND_INVALID
            end select
         enddo
      enddo
   end subroutine translate_bnds_to_ints

!> \brief Convert integer-defined boundary type to human-readable string

   elemental function translate_ints_to_bnds(ibnd) result(bstr)

      use constants, only: dsetnamelen, BND_MPI, BND_FC, BND_MPI_FC,BND_PER, BND_REF, BND_OUT, BND_OUTD, BND_OUTH, BND_OUTHD, BND_COR, BND_SHE, BND_USER, BND_INVALID

      implicit none

      integer, intent(in)        :: ibnd !< integer boundary

      character(len=dsetnamelen) :: bstr !< output string boundary

      select case (ibnd)
         case (BND_PER)
            bstr = 'periodic'
         case (BND_REF)
            bstr = 'reflecting'
         case (BND_OUT)
            bstr = 'out'
         case (BND_OUTD)
            bstr = 'outd'
         case (BND_OUTH)
            bstr = 'outh'
         case (BND_OUTHD)
            bstr = 'outhd'
         case (BND_SHE)
            bstr = 'shearing'
         case (BND_COR)
            bstr = 'corner'
         case (BND_MPI)
            bstr = 'mpi'
         case (BND_FC)
            bstr = "FC"
         case (BND_MPI_FC)
            bstr = "mpi/FC"
         case (BND_USER)
            bstr = 'user'
         case (BND_INVALID)
            bstr = 'invalid'
         case default
            bstr = 'HORRIBLE!'
      end select

   end function translate_ints_to_bnds

!> \brief Print computational domain details

   subroutine print_me(this)

      use constants,  only: LO, HI
      use dataio_pub, only: printinfo, msg

      implicit none

      class(domain_container), intent(inout) :: this  !< object invoking type-bound procedure

      write(msg,'(a,3(F5.1,1X))') "LO edge: ", this%edge(:,LO)  ; call printinfo(msg)
      write(msg,'(a,3(F5.1,1X))') "HI edge: ", this%edge(:,HI)  ; call printinfo(msg)
      write(msg,'(a,3(I4,1X))')   "n_d    : ", this%n_d(:)      ; call printinfo(msg)
      write(msg,'(a,3(I4,1X))')   "n_t    : ", this%n_t(:)      ; call printinfo(msg)
      write(msg,'(a,1(I4,1X))')   "nb     : ", this%nb          ; call printinfo(msg)
      write(msg,'(a,3(I4,1X))')   "LO bnd : ", this%bnd(:,LO)   ; call printinfo(msg)
      write(msg,'(a,3(I4,1X))')   "HI bnd : ", this%bnd(:,HI)   ; call printinfo(msg)
      write(msg,'(a,3(F5.1,1X))') "L_     : ", this%L_(:)       ; call printinfo(msg)
      write(msg,'(a,3(F5.1,1X))') "C_     : ", this%C_(:)       ; call printinfo(msg)
      write(msg,'(a,1(F5.1,1X))') "Vol    : ", this%Vol         ; call printinfo(msg)
      write(msg,'(a,3(L1,1X))')   "period : ", this%periodic(:) ; call printinfo(msg)

   end subroutine print_me

!> \brief Initialize all variables of domain_container type

   subroutine init(this, nb, n_d, bnds, edges, geometry, off)

      use constants,  only: ndims, LO, HI, big_float, dpi, xdim, ydim, zdim, GEO_XYZ, GEO_RPZ, GEO_INVALID, BND_PER, BND_REF, BND_SHE, I_ONE, I_TWO
      use dataio_pub, only: die, warn, msg
      use func,       only: operator(.notequals.)
      use mpisetup,   only: master

      implicit none

      class(domain_container),                   intent(inout) :: this     !< object invoking type-bound procedure
      integer(kind=4),                           intent(in)    :: nb       !< number of boundary cells surrounding the physical domain, same for all directions
      integer(kind=4),  dimension(ndims)                       :: n_d      !< number of %grid cells in physical domain without boundary cells (where  == 1 then that dimension is reduced to a point with no boundary cells)
      character(len=*), dimension(HI*ndims),     intent(in)    :: bnds     !< Six strings, describing boundary conditions
      real,             dimension(ndims, LO:HI), intent(inout) :: edges    !< physical domain boundaries position
      character(len=*),                          intent(in)    :: geometry !< define system of coordinates: "cartesian" or "cylindrical"
      integer(kind=8), dimension(ndims),         intent(in)    :: off      !< offset of the base level

      real    :: xmno, ymno, ymxo
      integer :: d

      ! Sanitize input parameters, if possible
      this%n_d(:) = max(I_ONE, n_d(:))
      this%has_dir(:) = this%n_d(:) > 1
      this%eff_dim = count(this%has_dir(:))
      where (this%has_dir(:))
         this%D_(:) = 1
         this%off = off
      elsewhere
         this%D_(:) = 0
         this%off = 0
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
            if (edges(ydim,LO) <= -0.99*huge(1.) .or. edges(ydim,HI) >= 0.99*huge(1.)) then
               if (master) call warn("[domain:init_domain] y range not specified. Defaulting to [0..1]")
               if (edges(ydim,LO) <= -0.99*huge(1.)) edges(ydim,LO) = 0.
               if (edges(ydim,HI) >=  0.99*huge(1.)) edges(ydim,HI) = 1.
            endif
         case (GEO_RPZ)
            if (edges(ydim,LO) <= -big_float) edges(ydim,LO) = 0.
            if (edges(ydim,HI) >= big_float) edges(ydim,HI) = dpi
            if (edges(ydim,HI)-edges(ydim,LO) > dpi) then
               if (master) call warn("[domain:init_domain] Hyperbolic spaces are not implemented. Setting azimuthal span to 2pi.")
               if (abs(edges(ydim,LO)) < 1./epsilon(1.)) then
                  edges(ydim,HI) = edges(ydim,LO) + dpi
               else
                  edges(ydim,LO) = edges(ydim,HI) - dpi
               endif
               if (abs(edges(ydim,HI)-edges(ydim,LO) - dpi) > 100*epsilon(1.)) call die("[domain:init_domain] absolute values for both edges(ydim,HI) and edges(ydim,LO) too high.") ! magic number
            endif
            if (edges(xdim,LO) <= 0.) then
               edges(xdim,LO) = 0.
               if (this%bnd(xdim, LO) /= BND_REF) then
                  if (master) call warn("[domain:init_domain] Enforcing this%bnd(xdim, LO) = 'ref'.")
               endif
               this%bnd(xdim, LO) = BND_REF
            endif
            if (this%bnd(xdim, HI) == BND_PER) call die("[domain:init_domain] Periodicity in radial direction is not allowed in cylindrical coordinates")
         case default
            call die("[domain:init_domain] Invalid geometry type.")
      end select
      if (xmno .notequals. edges(xdim,LO)) then
         write(msg,'(2(a,g20.12))')"[domain:init_domain] Sanitized edges(xdim,LO): ",xmno," -> ",edges(xdim,LO)
         if (master) call warn(msg)
      endif
      if ((ymno .notequals. edges(ydim,LO)) .and. (ymno .notequals. -big_float)) then
         write(msg,'(2(a,g20.12))')"[domain:init_domain] Sanitized edges(ydim,LO): ",ymno," -> ",edges(ydim,LO)
         if (master) call warn(msg)
      endif
      if ((ymxo .notequals. edges(ydim,HI)) .and. (ymno .notequals. -big_float)) then
         write(msg,'(2(a,g20.12))')"[domain:init_domain] Sanitized edges(ydim,HI): ",ymno," -> ",edges(ydim,HI)
         if (master) call warn(msg)
      endif
      if (edges(xdim,LO) > edges(xdim,HI)) call die("[[domain:init_domain] Negative span in X-direction")
      if (edges(ydim,LO) > edges(ydim,HI)) call die("[[domain:init_domain] Negative span in Y-direction")
      if (edges(zdim,LO) > edges(zdim,HI)) call die("[[domain:init_domain] Negative span in Z-direction")
      if (nb < 1) call die("[[domain:init_domain] no guardcells")

      ! set up the global domain
      this%nb = nb

      this%edge(:,:) = edges(:,:)

      ! finish up with the rest of domain_container members

      this%periodic(:) = .false.
      do d = xdim, zdim
         if ((any(this%bnd(d, :) == BND_PER) .or. (d==xdim .and. any(this%bnd(d, :) == BND_SHE))) .and. this%has_dir(d)) then
            this%periodic(d) = .true.
            if (this%bnd(d, LO) /= this%bnd(d, HI)) call die("[domain:set_derived] Periodic BC do not match")
         endif
      enddo
      if (any(this%bnd(ydim:zdim, :) == BND_SHE)) call die("[domain:set_derived] Shearing BC not allowed for y- and z-direction")

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

   end subroutine init

!> \brief Call this routine after you update base level shape, before you add new blocks

   subroutine  modify_side(this, d, lh, n)

      use constants,  only: cbuff_len, LO, HI, xdim, ydim, zdim, ndims, I_ONE
      use dataio_pub, only: warn, die, msg
      use mpisetup,   only: master

      implicit none

      class(domain_container), intent(in) :: this     !< object invoking type-bound procedure
      integer,                 intent(in) :: d        !< direction to be updated
      integer,                 intent(in) :: lh       !< side to be updated
      integer(kind=4),         intent(in) :: n        !< how many dells are added/removed

      character(len=cbuff_len), dimension(HI*ndims) :: bnds  !< Six strings, describing boundary conditions
      real, dimension(ndims, LO:HI)     :: edges
      real :: d_edge

      if (n==0) then
         call warn("[domain:modify_side] nothing to change")
         return
      endif

      if (.not. this%has_dir(d)) call die("[domain:modify_side] Cannot increase dimenionality of the simulation")
      if (this%periodic(d)) call die("[domain:modify_side] Cannot change periodic direction")
      if (this%nb > this%n_d(d)+n) call die("[domain:modify_side] Cannot shrink that much")

      bnds = [bnd_xl, bnd_xr, bnd_yl, bnd_yr, bnd_zl, bnd_zr]

      d_edge = n/real(this%n_d(d))*this%L_(d)
      select case (HI*(d-xdim)+lh)
         case (HI*(xdim-xdim)+LO)
            xmin = xmin - d_edge
            write(msg, '(a,g15.7)')"New xmin = ", xmin
         case (HI*(xdim-xdim)+HI)
            xmax = xmax + d_edge
            write(msg, '(a,g15.7)')"New xmax = ", xmax
         case (HI*(ydim-xdim)+LO)
            ymin = ymin - d_edge
            write(msg, '(a,g15.7)')"New ymin = ", ymin
         case (HI*(ydim-xdim)+HI)
            ymax = ymax + d_edge
            write(msg, '(a,g15.7)')"New ymax = ", ymax
         case (HI*(zdim-xdim)+LO)
            zmin = zmin - d_edge
            write(msg, '(a,g15.7)')"New zmin = ", zmin
         case (HI*(zdim-xdim)+HI)
            zmax = zmax + d_edge
            write(msg, '(a,g15.7)')"New zmax = ", zmax
      end select
      n_d(d) = n_d(d) + n
      write(msg(len_trim(msg)+1:),'(a,3i6,a)')". New n_d = [",n_d,"]"
      if (master) call warn(msg) ! As long as the restart file does not automagically recognize changed parameters, this message should be easily visible

      if (lh==LO) dom_off(d) = dom_off(d) - n

      edges = reshape( [xmin, ymin, zmin, xmax, ymax, zmax], shape=[ndims,HI-LO+I_ONE] )

      call dom%init(nb, n_d, bnds, edges, geometry, dom_off)

   end subroutine modify_side

end module domain
