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

!> \brief This module contains variables and initialization routines related to refinement

module refinement

   implicit none

   private
   public :: init_refinement, ref_flag, level_max, level_min, n_updAMR, allow_face_rstep, allow_corner_rstep, oop_thr

   type :: ref_flag
      logical :: refine   !> a request to refine
      logical :: derefine !> a request to derefine
   contains
      procedure :: sanitize
   end type ref_flag

   integer(kind=4), protected :: level_min          !< minimum allowed refinement
   integer(kind=4), protected :: level_max          !< maximum allowed refinement (don't need to be reached if not necessary)
   integer(kind=4), protected :: n_updAMR           !< how often to update the refinement structure
   logical,         protected :: allow_face_rstep   !< Allows >1 refinement step across faces (do not use it for any physical problems)
   logical,         protected :: allow_corner_rstep !< Allows >1 refinement step across edges and corners (do not use it for any physical problems)
   real,            protected :: oop_thr            !< Maximum allowed ratio of Out-of-Place grid pieces (according to current ordering scheme)

   namelist /AMR/ level_min, level_max, n_updAMR, allow_face_rstep, allow_corner_rstep, oop_thr

contains

!> \brief Initialization of parameters of refinement mechanics

   subroutine init_refinement

      use constants,  only: base_level_id, PIERNIK_INIT_DOMAIN, xdim, zdim, I_ONE
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
      allow_face_rstep   = .false.
      allow_corner_rstep = .false.
      allow_AMR = .true.
      oop_thr = 0.1
      do d = xdim, zdim
         if (dom%has_dir(d))  then
            if (AMR_bsize(d) < dom%nb) then
               if (allow_AMR) call warn("[refinement:init_refinement] Refinements disabled (AMR_bsize too small)")
               allow_AMR = .false.
            else
               if (mod(dom%n_d(d), AMR_bsize(d)) /= 0) then
                  if (allow_AMR) call warn("[refinement:init_refinement] Refinements disabled (domain not divisible by AMR_bsize)")
                  allow_AMR = .false.
               endif
            endif
         endif
      enddo

      if (master) then

         diff_nml(AMR)

         ! sanitizing
         if (allow_AMR) then
            level_min = max(level_min, base_level_id)
            level_max = max(level_max, level_min)
         else
            level_min = base_level_id
            level_max = base_level_id
            n_updAMR  = huge(I_ONE)
         endif

         ibuff(1) = level_min
         ibuff(2) = level_max
         ibuff(3) = n_updAMR

         lbuff(1) = allow_face_rstep
         lbuff(2) = allow_corner_rstep
         lbuff(3) = allow_AMR

         rbuff(1) = oop_thr

      endif

      call piernik_MPI_Bcast(ibuff)
      call piernik_MPI_Bcast(lbuff)

      if (slave) then

         level_min = ibuff(1)
         level_max = ibuff(2)
         n_updAMR  = ibuff(3)

         allow_face_rstep   = lbuff(1)
         allow_corner_rstep = lbuff(2)
         allow_AMR          = lbuff(3)

         oop_thr = rbuff(1)

      endif

      if (.not. allow_AMR) AMR_bsize=0

      ! Such large refinements may require additional work in I/O routines, visualization, computing MPI tags and so on.
      if (level_max > 40) call warn("[refinement:init_refinement] BEWARE: At such large refinements, integer overflows may happen under certain conditions.")

   end subroutine init_refinement

!> \brief sanitize the refinement flags

   subroutine sanitize(this, my_level)

      implicit none

      class(ref_flag), intent(inout) :: this     ! object invoking this procedure
      integer,         intent(in)    :: my_level ! refinement level at which the flag has to be sanitized

      if (my_level >= level_max) this%refine   = .false.
      if (my_level <  level_min) this%refine   = .true.

      if (this%refine) this%derefine = .false.

      if (my_level >  level_max) this%derefine = .true.
      if (my_level <= level_min) this%derefine = .false.

   end subroutine sanitize

end module refinement
