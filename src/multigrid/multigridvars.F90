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

!!$ ============================================================================
!>
!! \brief Variables and data structures required by multigrid routines.
!!
!! \details These variables are not meant to be accessed (or worse: altered) from the outside of multigrid routines (Fortran has no friends :-( )
!!
!! The multigrid solver relies typically at most on 2 guardcells, but the impact on performance is not big so we've decided to have it uniform with the rest of Piernik
!<

module multigridvars
! pulled by MULTIGRID

   use constants,   only: xdim, zdim, LO, HI, dsetnamelen
   use cg_list_lev, only: cg_list_level

   implicit none

   public ! QA_WARN no secrets are kept here
   private :: xdim, zdim, LO, HI, dsetnamelen, cg_list_level ! QA_WARN prevent re-exporting

   ! these 4 variables are required for basic use of the multigrid solver
   character(len=dsetnamelen), parameter :: source_n     = "source"     !< density field
   character(len=dsetnamelen), parameter :: solution_n   = "solution"   !< iterated solution (potential) fields
   character(len=dsetnamelen), parameter :: defect_n     = "defect"     !< defect field (effectively the density not accounted in current solution)
   character(len=dsetnamelen), parameter :: correction_n = "correction" !< correction to the potential to be applied at the end of V-cycle
   integer :: source, solution, defect, correction !< indices to the fields described above

   ! these constants should be moved to constants module

   ! namelist parameters
   integer            :: ord_prolong                                  !< Prolongation operator order; allowed values are -4 -2, 0 (default), 2 and 4; -2 is often fast
   integer            :: ord_prolong_face_norm                        !< Face prolongation operator order in the direction normal to the face; allowed values are 0, 1  and 2
   integer            :: ord_prolong_face_par                         !< Face prolongation operator order in the directions parallel to the face; allowed values are -2 .. 2
   logical            :: stdout                                       !< print verbose messages to stdout
   logical            :: verbose_vcycle                               !< Print one line of log per V-cycle, summary otherwise

   ! boundaries
   enum, bind(C)                                                      !< constants for enumerating multigrid boundary types
      enumerator :: bnd_periodic                                      !< periodic
      enumerator :: bnd_dirichlet                                     !< 0-value boundary type (uniform Dirichlet)
      enumerator :: bnd_isolated                                      !< isolated boundary type
      enumerator :: bnd_neumann                                       !< 0-gradient boundary type (uniform Neumann)
      enumerator :: bnd_givenval                                      !< given value boundary type (general Dirichlet)
      enumerator :: bnd_invalid = bnd_periodic - 1                    !< invalid
   end enum
   logical, dimension(xdim:zdim, LO:HI) :: is_external                !< .true. for non-"mpi" local domain boundaries
   enum, bind(C)
      enumerator :: extbnd_donothing                                  !< Do not touch external boundaries
      enumerator :: extbnd_zero                                       !< Fill external boundaries with zeroes
      enumerator :: extbnd_extrapolate                                !< Perform extrapolation in external boundaries
      enumerator :: extbnd_mirror                                     !< Zero-gradient, mirroring external boundaries
      enumerator :: extbnd_antimirror = - extbnd_mirror               !< mirroring external boundaries with opposite sign
   end enum

   ! miscellaneous
   real                    :: ts                                      !< time for runtime profiling
   real                    :: tot_ts                                  !< total multigrid time
   logical                 :: is_mg_uneven                            !< .true. when domain shapes differ across procesors, even on the coarsest grids
   logical                 :: single_base                             !< .true. when the whole base level is located on a single cpu
   logical                 :: need_general_pf                         !< .false. only for most regular domain decomposition

   integer, parameter :: prefix_len = 3                               !< length of prefix for distinguishing V-cycles in the log
   type :: vcycle_stats
      real, allocatable, dimension(:) :: factor                       !< norm reduction factor
      real, allocatable, dimension(:) :: time                         !< time spent
      integer                         :: count                        !< number of executed V-cycles
      real                            :: norm_rhs                     !< norm of the source
      real                            :: norm_final                   !< norm of the defect relative to the source
      character(len=prefix_len)       :: cprefix                      !< prefix for distinguishing V-cycles in the log (e.g inner or outer potential, CR component)
    contains
       procedure :: init
       procedure :: brief_v_log
   end type vcycle_stats

   type(cg_list_level), pointer                           :: base     !< pointer to coarsest level
   type(cg_list_level), pointer                           :: roof     !< pointer to finest level

contains

!> \brief Initialize vcycle_stats

   subroutine init(this, size)

      use dataio_pub,    only: die

      implicit none

      class(vcycle_stats), intent(out) :: this   !< V-cycle statistics variable to be created or reset
      integer,             intent(in)  :: size !< size of the vs structure (usually max_cycles); for nonpositive value perform reset only

      if (size > 0) then
         if (allocated(this%factor) .or. allocated(this%time)) call die("[multigridhelpers:vcycle_stats_init] vcycle_stats already allocated.")
         allocate(this%factor(0:size), this%time(0:size))
      endif

      this%factor(:)  = 0.
      this%time(:)    = 0.
      this%count      = 0
      this%norm_rhs   = 0.
      this%norm_final = 0.
      this%cprefix    = ""

   end subroutine init

!> \brief Assembles one-line log of V-cycle achievements

   subroutine brief_v_log(this)

      use constants,     only: fplen, fmt_len
      use mpisetup,      only: slave
      use dataio_pub,    only: msg, warn, printinfo

      implicit none

      class(vcycle_stats), intent(in) :: this

      real                   :: at
      integer                :: i, lm, ftype
      character(len=fplen)   :: normred
      character(len=fmt_len), parameter, dimension(2) :: fmt_norm = [ '(a,i3,1x,2a,f7.3,a,i3,a,f7.3,a,f13.10,a)', '(a,i3,1x,2a,f7.3,a,i3,a,f7.3,a,e13.6,a) ' ]

      if (slave) return

      if (this%count > ubound(this%factor, 1)) call warn("[multigridhelpers:brief_v_log] Trying to read beyond upper bound of vcycle_stats.")

      at = 0.
      if (this%count > 0) at = sum(this%time(1:this%count))/this%count ! average V-cycle time on PE# 0

      ftype = 1
      if (this%norm_final < 1e-8) ftype = 2
      write(msg, fmt_norm(ftype))"[multigrid] ", this%count, trim(this%cprefix), "cycles, dt_wall=", this%time(0), " +", this%count, "*", at, ", norm/rhs= ", this%norm_final, " : "

      do i = 0, min(this%count, ubound(this%factor, 1))
         if (this%factor(i) < 1.0e4) then
            write(normred, '(f8.2)') this%factor(i)
         else if (this%factor(i) < 1.0e7) then
            write(normred, '(f8.0)') this%factor(i)
         else
            write(normred, '(es9.2)') this%factor(i)
         endif
         lm = len_trim(msg)
         if (len(msg) >= lm + 9) msg(lm+2:lm+9) = normred(1:8)
      enddo

      call printinfo(msg, stdout)

   end subroutine brief_v_log

end module multigridvars
