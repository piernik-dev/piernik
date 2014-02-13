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

!> \brief This module contains variables and initialization routines related to refinement

module refinement

   use constants,       only: ndims, LO, HI
   use refinement_flag, only: level_min, level_max

   implicit none

   private
   public :: n_updAMR, allow_face_rstep, allow_corner_rstep, oop_thr, refine_points, &
        &    refine_boxes, init_refinement, emergency_fix, set_n_updAMR, strict_SFC_ordering, prefer_n_bruteforce

   integer(kind=4), protected :: n_updAMR            !< how often to update the refinement structure
   logical,         protected :: allow_face_rstep    !< Allows >1 refinement step across faces (do not use it for any physical problems)
   logical,         protected :: allow_corner_rstep  !< Allows >1 refinement step across edges and corners (do not use it for any physical problems)
   logical,         protected :: strict_SFC_ordering !< Enforce strict SFC ordering to allow optimized neighbour search
   real,            protected :: oop_thr             !< Maximum allowed ratio of Out-of-Place grid pieces (according to current ordering scheme)
   logical,         protected :: prefer_n_bruteforce !< if .false. then try DFC algorithms for neighbor searches

   ! some refinement primitives
   integer, parameter :: nshapes = 10

   !> \brief Refinement point
   type :: ref_point
      integer(kind=4)        :: level  !> desired level of refinement
      real, dimension(ndims) :: coords !> coordinates, where to refine
   end type ref_point
   type(ref_point), dimension(nshapes), protected :: refine_points

   !> \brief Refinement box
   type :: ref_box
      integer(kind=4)               :: level  !> desired level of refinement
      real, dimension(ndims, LO:HI) :: coords !> coordinates, where to refine
   end type ref_box
   type(ref_box), dimension(nshapes), protected :: refine_boxes

   logical :: emergency_fix !< set to .true. if you want to call update_refinement ASAP

   namelist /AMR/ level_min, level_max, n_updAMR, allow_face_rstep, allow_corner_rstep, strict_SFC_ordering, &
        &         prefer_n_bruteforce, oop_thr, refine_points, refine_boxes

contains

!> \brief Initialization of parameters of refinement mechanics

   subroutine init_refinement

      use constants,  only: base_level_id, PIERNIK_INIT_DOMAIN, xdim, ydim, zdim, I_ONE, LO, HI
      use dataio_pub, only: nh      ! QA_WARN required for diff_nml
      use dataio_pub, only: die, code_progress, warn
      use domain,     only: AMR_bsize, dom
      use mpisetup,   only: ibuff, lbuff, rbuff, master, slave, piernik_MPI_Bcast

      implicit none

      integer :: d
      logical :: allow_AMR

      if (code_progress < PIERNIK_INIT_DOMAIN) call die("[refinement:init_refinement] Domain not initialized.")

      level_min = base_level_id
      level_max = level_min
      n_updAMR  = huge(I_ONE)
      allow_face_rstep    = .false.
      allow_corner_rstep  = .false.
      strict_SFC_ordering = .false.
      allow_AMR = .true.
      prefer_n_bruteforce = .false.
      oop_thr = 0.1
      do d = xdim, zdim
         if (dom%has_dir(d))  then
            if (AMR_bsize(d) < dom%nb) then
               if (allow_AMR .and. master) call warn("[refinement:init_refinement] Refinements disabled (AMR_bsize too small)")
               allow_AMR = .false.
            else
               if (mod(dom%n_d(d), AMR_bsize(d)) /= 0) then
                  if (allow_AMR .and. master) call warn("[refinement:init_refinement] Refinements disabled (domain not divisible by AMR_bsize)")
                  allow_AMR = .false.
               endif
            endif
         endif
      enddo
      refine_points(:) = ref_point(base_level_id-1, [ 0., 0., 0.] )
      refine_boxes (:) = ref_box  (base_level_id-1, reshape([ 0., 0., 0., 0., 0., 0.], [ndims, HI-LO+I_ONE] ) )

      if (1 + 9*nshapes > ubound(rbuff, dim=1)) call die("[refinement:init_refinement] increase rbuff size") ! should be detected at compile time but it is only a warning
      if (master) then

         if (.not.nh%initialized) call nh%init()
         open(newunit=nh%lun, file=nh%tmp1, status="unknown")
         write(nh%lun,nml=AMR)
         close(nh%lun)
         open(newunit=nh%lun, file=nh%par_file)
         nh%errstr=""
         read(unit=nh%lun, nml=AMR, iostat=nh%ierrh, iomsg=nh%errstr)
         close(nh%lun)
         call nh%namelist_errh(nh%ierrh, "AMR")
         read(nh%cmdl_nml,nml=AMR, iostat=nh%ierrh)
         call nh%namelist_errh(nh%ierrh, "AMR", .true.)
         open(newunit=nh%lun, file=nh%tmp2, status="unknown")
         write(nh%lun,nml=AMR)
         close(nh%lun)
         call nh%compare_namelist()

         ! sanitizing
         if (allow_AMR) then
            level_min = max(level_min, base_level_id)
            level_max = max(level_max, level_min)
         else
            level_min = base_level_id
            level_max = base_level_id
            n_updAMR  = huge(I_ONE)
         endif
         where (.not. dom%has_dir(:)) AMR_bsize(:) = huge(I_ONE)

         ibuff(1) = level_min
         ibuff(2) = level_max
         ibuff(3) = n_updAMR
         ibuff(4        :3+  nshapes) = refine_points(:)%level
         ibuff(4+nshapes:3+2*nshapes) = refine_boxes (:)%level

         lbuff(1) = allow_face_rstep
         lbuff(2) = allow_corner_rstep
         lbuff(3) = allow_AMR
         lbuff(4) = strict_SFC_ordering
         lbuff(5) = prefer_n_bruteforce

         rbuff(1) = oop_thr
         rbuff(2          :1+  nshapes) = refine_points(:)%coords(xdim)
         rbuff(2+  nshapes:1+2*nshapes) = refine_points(:)%coords(ydim)
         rbuff(2+2*nshapes:1+3*nshapes) = refine_points(:)%coords(zdim)
         rbuff(2+3*nshapes:1+4*nshapes) = refine_boxes (:)%coords(xdim, LO)
         rbuff(2+4*nshapes:1+5*nshapes) = refine_boxes (:)%coords(xdim, HI)
         rbuff(2+5*nshapes:1+6*nshapes) = refine_boxes (:)%coords(ydim, LO)
         rbuff(2+6*nshapes:1+7*nshapes) = refine_boxes (:)%coords(ydim, HI)
         rbuff(2+7*nshapes:1+8*nshapes) = refine_boxes (:)%coords(zdim, LO)
         rbuff(2+8*nshapes:1+9*nshapes) = refine_boxes (:)%coords(zdim, HI)

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(lbuff)
      call piernik_MPI_Bcast(rbuff)

      if (slave) then

         level_min = ibuff(1)
         level_max = ibuff(2)
         n_updAMR  = ibuff(3)
         refine_points(:)%level = ibuff(4        :3+  nshapes)
         refine_boxes (:)%level = ibuff(4+nshapes:3+2*nshapes)

         allow_face_rstep    = lbuff(1)
         allow_corner_rstep  = lbuff(2)
         allow_AMR           = lbuff(3)
         strict_SFC_ordering = lbuff(4)
         prefer_n_bruteforce = lbuff(5)

         oop_thr = rbuff(1)
         refine_points(:)%coords(xdim)     = rbuff(2          :1+  nshapes)
         refine_points(:)%coords(ydim)     = rbuff(2+  nshapes:1+2*nshapes)
         refine_points(:)%coords(zdim)     = rbuff(2+2*nshapes:1+3*nshapes)
         refine_boxes (:)%coords(xdim, LO) = rbuff(2+3*nshapes:1+4*nshapes)
         refine_boxes (:)%coords(xdim, HI) = rbuff(2+4*nshapes:1+5*nshapes)
         refine_boxes (:)%coords(ydim, LO) = rbuff(2+5*nshapes:1+6*nshapes)
         refine_boxes (:)%coords(ydim, HI) = rbuff(2+6*nshapes:1+7*nshapes)
         refine_boxes (:)%coords(zdim, LO) = rbuff(2+7*nshapes:1+8*nshapes)
         refine_boxes (:)%coords(zdim, HI) = rbuff(2+8*nshapes:1+9*nshapes)

      endif

      if (.not. allow_AMR) AMR_bsize=0

      ! Such large refinements may require additional work in I/O routines, visualization, computing MPI tags and so on.
      if (level_max > 40) call warn("[refinement:init_refinement] BEWARE: At such large refinements, integer overflows may happen under certain conditions.")

      emergency_fix = .false.

   end subroutine init_refinement

!> \brief Change the protected parameter n_updAMR

   subroutine set_n_updAMR(n)

      use dataio_pub, only: die
      use mpisetup,   only: piernik_MPI_Bcast

      implicit none

      integer(kind=4), intent(in) :: n

      integer(kind=4) :: nn

      nn = n
      call piernik_MPI_Bcast(nn)
      if (nn /= n) call die("[refinement:set_n_updAMR] n_updAMR goes out of sync")
      n_updAMR = nn

   end subroutine set_n_updAMR

end module refinement
