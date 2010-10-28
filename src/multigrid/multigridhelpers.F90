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

!$ ============================================================================
!!
!! This module contains variables and routines useful mostly for developing and debugging
!!

module multigridhelpers

   implicit none

   private

   public :: set_dirty, check_dirty, vcycle_stats_init, brief_v_log, mg_write_log, numbered_ascii_dump
   public :: do_ascii_dump, multidim_code_3D, dirty_debug, dirty_label, dirtyH, dirtyL, aux_par_I0, aux_par_I1, aux_par_I2, aux_par_R0, aux_par_R1, aux_par_R2

   ! namelist parameters
   logical            :: do_ascii_dump                      !< to dump, or not to dump: that is a question (ascii)
   logical            :: multidim_code_3D                   !< prefer code written for any 1D and 2D configuration even in 3D for benchmarking and debugging
   logical            :: dirty_debug                        !< Initialize everything with some insane values (dirtyH, defined below) and check if they can propagate
   integer            :: aux_par_I0, aux_par_I1, aux_par_I2 !< auxiliary integer parameters
   real               :: aux_par_R0, aux_par_R1, aux_par_R2 !< auxiliary real parameters
   integer, parameter    :: dl_len = 64                     !<
   character(len=dl_len) :: dirty_label                     !< buffer for label for check_dirty subroutine

   real, parameter    :: dirtyH = 1e200, dirtyL = 1e50      !< If dirty_debug, initialize arrays with dirtyH and check if the solution contains anything above dirtyL

contains

!/todo write also general set_dirty_array and check_dirty_array routines so there will be no need to use dirtyH and dirtyL outside multigridhelpers module.

!$ ============================================================================
!!
!! This routine pollutes selected multigrid variable with an insane value dirtyH.
!! If anything in the multigrid works by accident, through compiler-dependent initialization or unintentional relying on outdated values,
!! the insane value should pollute the solution in an easily visible way.
!!

   subroutine set_dirty(iv)

      use dataio_pub,    only: die
      use multigridvars, only: ngridvars, lvl, level_min, level_max

      implicit none

      integer, intent(in) :: iv   !< index of variable in lvl()%mgvar which we want to pollute

      integer :: l

      if (.not. dirty_debug) return

      if (iv < 1 .or. iv > ngridvars) call die("[multigridhelpers:set_dirty] Invalid variable index.")

      do l = level_min, level_max
         lvl(l)%mgvar(:, :, :, iv) = dirtyH
      enddo

   end subroutine set_dirty

!!$ ============================================================================
!!
!! This routine checks for detectable traces of set_dirty calls.
!!

   subroutine check_dirty(lev, iv, label, expand)

      use dataio_pub,    only: die, warn, msg
      use mpisetup,      only: proc
      use multigridvars, only: ngridvars, lvl, level_min, level_max, mg_nb, eff_dim, NDIM

      implicit none

      integer,           intent(in) :: iv     !< index of variable in lvl()%mgvar which we want to pollute
      integer,           intent(in) :: lev    !< level which we are checking. Invalid means all levels
      character(len=*),  intent(in) :: label  !< label to indicate the origin of call
      integer, optional, intent(in) :: expand !< also check guardcells

      integer :: l, l1, l2, i, j, k, ng

      if (.not. dirty_debug) return

      if (iv < 1 .or. iv > ngridvars) call die("[multigridhelpers:check_dirty] Invalid variable index.")

      if (lev < level_min .or. lev > level_max) then
         l1 = level_min
         l2 = level_max
      else
         l1 = lev
         l2 = lev
      endif

      if (present(expand) .and. eff_dim==NDIM) then ! for 1D and 2D one should define ng_x,ng_y and ng_z
         if (expand > mg_nb) then
            ng = mg_nb
         else
            ng = expand
         endif
      else
         ng = 0
      endif

      do l = l1, l2
         do k = lvl(l)%ks-ng, lvl(l)%ke+ng
            do j = lvl(l)%js-ng, lvl(l)%je+ng
               do i = lvl(l)%is-ng, lvl(l)%ie+ng
                  if (abs(lvl(l)%mgvar(i, j, k, iv)) > dirtyL) then
                     write(msg, '(3a,i4,a,i2,a,3(i3,a),i2,a,g20.12)') &
                        "[multigridhelpers:check_dirty] ", label, "@", proc, " lvl(", l, ")%mgvar(", i, ",", j, ",", k, ",", iv, ") = ", lvl(l)%mgvar(i, j, k, iv)
                     call warn(msg)
                  endif
               enddo
            enddo
         enddo
      enddo

   end subroutine check_dirty

!!$ ============================================================================
!!
!! Initialize vcycle_stats
!!

   subroutine vcycle_stats_init(vs, size)

      use dataio_pub,    only: die, ierrh
      use multigridvars, only: vcycle_stats

      implicit none

      type(vcycle_stats), intent(out) :: vs   !< V-cycle statistics variable to be created or reset
      integer, intent(in)             :: size !< size of the vs structure (usually max_cycles); for nonpositive value perform reset only

      if (size > 0) then
         if (allocated(vs%factor) .or. allocated(vs%time)) call die("[multigridhelpers:vcycle_stats_init] vcycle_stats already allocated.")
         allocate(vs%factor(0:size), stat=ierrh)
         if (ierrh /= 0) call die("[multigridhelpers:vcycle_stats_init] Allocation error: vs%factor.")
         allocate(vs%time(0:size), stat=ierrh)
         if (ierrh /= 0) call die("[multigridhelpers:vcycle_stats_init] Allocation error: vs%time.")
      endif

      vs%factor(:)     = 0.
      vs%time(:)       = 0.
      vs%count         = 0
      vs%norm_rhs_orig = 0.
      vs%norm_final    = 0.
      vs%cprefix       = ""

   end subroutine vcycle_stats_init

!!$ ============================================================================
!!
!! Assembles one-line log of V-cycle achievements
!!

   subroutine brief_v_log(vs)

      use multigridvars, only: stdout, vcycle_stats
      use mpisetup,      only: proc
      use dataio_pub,    only: msg, fplen, warn

      implicit none

      type(vcycle_stats), intent(in) :: vs

      real                 :: at
      integer              :: i, lm
      character(len=fplen) :: normred
      character(len=1)     :: dash

      if (proc /= 0) return

      if (vs%count > ubound(vs%factor, 1)) call warn("[multigridhelpers:brief_v_log] Trying to read beyond upper bound of vcycle_stats.")

      at = 0.
      if (vs%count > 0) at = sum(vs%time(1:vs%count))/vs%count ! average V-cycle time on PE# 0

      dash = ""
      if (len_trim(vs%cprefix) > 0) dash="-"

      write(msg, '(a,i3,1x,3a,f7.3,a,i3,a,f7.3,a,f11.9,a)') &
           "[multigrid] ", vs%count, trim(vs%cprefix), dash, "cycles, dt_wall=", vs%time(0), " +", vs%count, "*", at, ", norm/rhs= ", vs%norm_final, " : "

      do i = 0, min(vs%count, ubound(vs%factor, 1))
         if (vs%factor(i) < 1.0e4) then
            write(normred, '(f8.2)') vs%factor(i)
         else if (vs%factor(i) < 1.0e7) then
            write(normred, '(f8.0)') vs%factor(i)
         else
            write(normred, '(es8.2)') vs%factor(i)
         endif
         lm = len_trim(msg)
         if (len(msg) >= lm + 9) msg(lm+2:lm+9) = normred(1:8)
      enddo

      call mg_write_log(msg, stdout)

   end subroutine brief_v_log

!!$ ============================================================================
!!
!! Logfile
!!

   subroutine mg_write_log(msg, stdout_cntrl)

      use dataio_pub,    only: printinfo
      use mpisetup,      only: proc

      implicit none

      character(len=*), intent(in)  :: msg          !< message to be logged
      logical, optional, intent(in) :: stdout_cntrl !< flag to suppress printing to stdout

      logical            :: to_stdout

      to_stdout = (proc == 0)
      if (present(stdout_cntrl)) to_stdout = to_stdout .and. stdout_cntrl
      call printinfo(msg, to_stdout)

   end subroutine mg_write_log

!!$ ============================================================================
!!
!! Emergency routine for quick ASCII dumps
!!
!! Absolute integer coordinates also allow seamless concatenation of dumps made by all PEs.
!!

   subroutine ascii_dump(filename)

      use dataio_pub,    only: msg
      use mpisetup,      only: proc
      use multigridvars, only: level_min, level_max, lvl, gb_cartmap, ngridvars, XDIR, YDIR, ZDIR

      implicit none

      character(len=*), intent(in) :: filename !< name to write the emergency dump

      integer, parameter :: fu=30
      integer            :: l, i, j, k

      if (.not. do_ascii_dump) return

      open(fu, file=filename, status="unknown")

      write(fu, '("#",a3,2a4,a6,10a20,/)')"i", "j", "k", "level", "x(i)", "y(j)", "z(k)", "source", "solution", "defect", "correction", "bx", "by", "bz"
      do l = level_min, level_max
         do i = lvl(l)%is, lvl(l)%ie
            do j = lvl(l)%js, lvl(l)%je
               do k = lvl(l)%ks, lvl(l)%ke
                  write(fu, '(3i4,i6,10es20.11e3)')i-lvl(l)%is+gb_cartmap(proc)%proc(XDIR)*lvl(l)%nxb, &
                       &                          j-lvl(l)%js+gb_cartmap(proc)%proc(YDIR)*lvl(l)%nyb, &
                       &                          k-lvl(l)%ks+gb_cartmap(proc)%proc(ZDIR)*lvl(l)%nzb, &
                       &                          l, lvl(l)%x(i), lvl(l)%y(j), lvl(l)%z(k), lvl(l)%mgvar(i, j, k, 1:ngridvars)
               enddo
               write(fu, '(/)')
            enddo
            write(fu, '(/,/)')
         enddo
         write(fu, '(/,/,/)')
      enddo

      close(fu)

      if (proc == 0) then
         write(msg,'(3a)') "[multigridhelpers:ascii_dump] Wrote dump '",filename,"'"
         call mg_write_log(msg)
      endif

   end subroutine ascii_dump

!!$ ============================================================================
!!
!! Construct name of emergency ASCII dump
!!

   subroutine numbered_ascii_dump(basename, a)

      use mpisetup,      only: proc, nstep
      use dataio_pub,    only: halfstep, msg

      implicit none

      character(len=*), intent(in)  :: basename !< first part of the filename
      integer, optional, intent(in) :: a        !< additional number

      integer             :: l, n

      if (.not. do_ascii_dump) return

      n = 2 * (nstep - 1)
      if (halfstep) n = n + 1

      if (present(a)) then
         write(msg, '(a,i4,i6,i3)') trim(basename), proc, n, a
      else
         write(msg, '(a,i4,i6)')    trim(basename), proc, n
      endif
      do l = 1, len_trim(msg)
         if (msg(l:l) == " ") msg(l:l) = "_"
      enddo
      call ascii_dump(trim(msg))

   end subroutine numbered_ascii_dump

end module multigridhelpers
