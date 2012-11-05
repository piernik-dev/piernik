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

!> \brief Simple structure that records progress of the V-cycle

module multigrid_vstats
! pulled by MULTIGRID

   implicit none

   private
   public :: vcycle_stats

   integer, parameter :: prefix_len = 3              !< length of prefix for distinguishing V-cycles in the log

   type :: vcycle_stats

      real, allocatable, dimension(:) :: factor      !< norm reduction factor
      real, allocatable, dimension(:) :: time        !< time spent
      integer                         :: count       !< number of executed V-cycles
      real                            :: norm_rhs    !< norm of the source
      real                            :: norm_final  !< norm of the defect relative to the source
      character(len=prefix_len)       :: cprefix     !< prefix for distinguishing V-cycles in the log (e.g inner or outer potential, CR component)

    contains

       procedure :: init                             !< Initialize vcycle_stats
       procedure :: brief_v_log                      !< Assembles one-line log of V-cycle achievements
       procedure :: cleanup                          !< Free the memory

   end type vcycle_stats

contains

!> \brief Initialize vcycle_stats

   subroutine init(this, size)

      use dataio_pub,    only: die

      implicit none

      class(vcycle_stats), intent(out) :: this !< V-cycle statistics variable to be created or reset
      integer(kind=4),     intent(in)  :: size !< size of the vs structure (usually max_cycles); for nonpositive value perform reset only

      if (size >= 0) then
         if (allocated(this%factor) .or. allocated(this%time)) call die("[multigrid_vstats:vcycle_stats_init] vcycle_stats already allocated.")
         allocate(this%factor(0:size), this%time(0:size))
      else
         call die("[multigrid_vstats:vcycle_stats_init] size < 0")
      endif

      this%factor(:)  = 0.
      this%time(:)    = 0.
      this%count      = 0
      this%norm_rhs   = 0.
      this%norm_final = 0.
      this%cprefix    = ""

   end subroutine init

!> \brief Free the memory

   subroutine cleanup(this)

      implicit none

      class(vcycle_stats), intent(inout) :: this

      if (allocated(this%factor)) deallocate(this%factor)
      if (allocated(this%time)) deallocate(this%time)

   end subroutine cleanup

!> \brief Assembles one-line log of V-cycle achievements

   subroutine brief_v_log(this)

      use constants,     only: fplen, fmt_len
      use mpisetup,      only: slave
      use dataio_pub,    only: msg, warn, printinfo
      use multigridvars, only: stdout

      implicit none

      class(vcycle_stats), intent(in) :: this

      real                   :: at
      integer                :: i, lm, ftype
      character(len=fplen)   :: normred
      character(len=fmt_len), parameter, dimension(2) :: fmt_norm = [ '(a,i3,1x,2a,f7.3,a,i3,a,f7.3,a,f13.10,a)', '(a,i3,1x,2a,f7.3,a,i3,a,f7.3,a,e13.6,a) ' ]

      if (slave) return

      if (this%count > ubound(this%factor, 1)) call warn("[multigrid_vstats:brief_v_log] Trying to read beyond upper bound of vcycle_stats.")

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
            write(normred, '(es8.1)') this%factor(i)
         endif
         lm = len_trim(msg)
         if (len(msg) >= lm + 9) msg(lm+2:lm+9) = normred(1:8)
      enddo

      call printinfo(msg, stdout)

   end subroutine brief_v_log

end module multigrid_vstats
