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
!! \brief This module contains variables and routines useful mostly for developing and debugging
!<

module multigridhelpers
! pulled by MULTIGRID
   implicit none

   private

   public :: set_dirty, check_dirty, vcycle_stats_init, brief_v_log, mg_write_log, numbered_ascii_dump
   public :: do_ascii_dump, multidim_code_3D, dirty_debug, dirty_label, dirtyH, dirtyL

   ! namelist parameters
   logical            :: do_ascii_dump                      !< to dump, or not to dump: that is a question (ascii)
   logical            :: multidim_code_3D                   !< prefer code written for any 1D and 2D configuration even in 3D for benchmarking and debugging
   logical            :: dirty_debug                        !< Initialize everything with some insane values (dirtyH, defined below) and check if they can propagate
   integer, parameter    :: dl_len = 64                     !< length of label buffer
   character(len=dl_len) :: dirty_label                     !< buffer for label for check_dirty subroutine
   real, parameter    :: big_s =  huge(real(1.0,4))         !< largest single-precision number
   real, parameter    :: dirtyH = big_s                     !< If dirty_debug then initialize arrays with this insane value
   real, parameter    :: dirtyL = sqrt(big_s)               !< If dirty_debug then check if the solution got contaminated by dirtyH by checking if it contains anything above dirtyL

contains

!> \todo write also general set_dirty_array and check_dirty_array routines so there will be no need to use dirtyH and dirtyL outside multigridhelpers module.

!!$ ============================================================================
!>
!! \brief This routine pollutes selected multigrid variable with an insane value dirtyH.
!! \details If anything in the multigrid works by accident, through compiler-dependent initialization or unintentional relying on outdated values,
!! the insane value should pollute the solution in an easily visible way.
!<

   subroutine set_dirty(iv)

      use dataio_pub,    only: die
      use multigridvars, only: ngridvars, plvl, base

      implicit none

      integer(kind=4), intent(in) :: iv   !< index of variable in lvl()%mgvar which we want to pollute

      type(plvl), pointer :: curl

      if (.not. dirty_debug) return

      if (iv < 1 .or. iv > ngridvars) call die("[multigridhelpers:set_dirty] Invalid variable index.")

      curl => base
      do while (associated(curl))
         curl%mgvar(:, :, :, iv) = dirtyH
         curl => curl%finer
      enddo

   end subroutine set_dirty

!!$ ============================================================================
!>
!! \brief This routine checks for detectable traces of set_dirty calls.
!<

   subroutine check_dirty(curl, iv, label, expand)

      use dataio_pub,    only: die, warn, msg
      use domain,        only: eff_dim
      use mpisetup,      only: proc
      use multigridvars, only: ngridvars, plvl, mg_nb
      use constants,     only: ndims

      implicit none

      integer(kind=4),   intent(in) :: iv     !< index of variable in lvl()%mgvar which we want to pollute
      type(plvl), pointer, intent(in) :: curl !< level which we are checking
      character(len=*),  intent(in) :: label  !< label to indicate the origin of call
      integer, optional, intent(in) :: expand !< also check guardcells

      integer :: i, j, k, ng

      if (.not. dirty_debug) return
      if (iv < 1 .or. iv > ngridvars) call die("[multigridhelpers:check_dirty] Invalid variable index.")
      if (curl%empty) return

      if (present(expand) .and. eff_dim==ndims) then ! for 1D and 2D one should define ng_x,ng_y and ng_z
         if (expand > mg_nb) then
            ng = mg_nb
         else
            ng = expand
         endif
      else
         ng = 0
      endif

      do k = curl%ks-ng, curl%ke+ng
         do j = curl%js-ng, curl%je+ng
            do i = curl%is-ng, curl%ie+ng
               if (abs(curl%mgvar(i, j, k, iv)) > dirtyL) then
!                        if (count([i<curl%is .or. i>curl%ie, j<curl%js .or. j>curl%je, k<curl%ks .or. k>curl%ke]) <=1) then ! excludes corners
                  write(msg, '(3a,i4,a,i2,a,3(i3,a),i2,a,g20.12)') &
                          "[multigridhelpers:check_dirty] ", trim(label), "@", proc, " lvl(", curl%level, ")%mgvar(", i, ",", j, ",", k, ",", iv, ") = ", curl%mgvar(i, j, k, iv)
                  call warn(msg)
!                        endif
               endif
            enddo
         enddo
      enddo

   end subroutine check_dirty

!!$ ============================================================================
!>
!! \brief Initialize vcycle_stats
!<

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

      vs%factor(:)  = 0.
      vs%time(:)    = 0.
      vs%count      = 0
      vs%norm_rhs   = 0.
      vs%norm_final = 0.
      vs%cprefix    = ""

   end subroutine vcycle_stats_init

!!$ ============================================================================
!>
!! \brief Assembles one-line log of V-cycle achievements
!<

   subroutine brief_v_log(vs)

      use constants,     only: fplen, fmt_len
      use multigridvars, only: stdout, vcycle_stats
      use mpisetup,      only: slave
      use dataio_pub,    only: msg, warn

      implicit none

      type(vcycle_stats), intent(in) :: vs

      real                   :: at
      integer                :: i, lm, ftype
      character(len=fplen)   :: normred
      character(len=fmt_len), parameter, dimension(2) :: fmt_norm = [ '(a,i3,1x,2a,f7.3,a,i3,a,f7.3,a,f13.10,a)', '(a,i3,1x,2a,f7.3,a,i3,a,f7.3,a,e13.6,a) ' ]

      if (slave) return

      if (vs%count > ubound(vs%factor, 1)) call warn("[multigridhelpers:brief_v_log] Trying to read beyond upper bound of vcycle_stats.")

      at = 0.
      if (vs%count > 0) at = sum(vs%time(1:vs%count))/vs%count ! average V-cycle time on PE# 0

      ftype = 1
      if (vs%norm_final < 1e-8) ftype = 2
      write(msg, fmt_norm(ftype))"[multigrid] ", vs%count, trim(vs%cprefix), "cycles, dt_wall=", vs%time(0), " +", vs%count, "*", at, ", norm/rhs= ", vs%norm_final, " : "

      do i = 0, min(vs%count, ubound(vs%factor, 1))
         if (vs%factor(i) < 1.0e4) then
            write(normred, '(f8.2)') vs%factor(i)
         else if (vs%factor(i) < 1.0e7) then
            write(normred, '(f8.0)') vs%factor(i)
         else
            write(normred, '(es9.2)') vs%factor(i)
         endif
         lm = len_trim(msg)
         if (len(msg) >= lm + 9) msg(lm+2:lm+9) = normred(1:8)
      enddo

      call mg_write_log(msg, stdout)

   end subroutine brief_v_log

!!$ ============================================================================
!>
!! \brief Logfile
!<

   subroutine mg_write_log(msg, stdout_cntrl)

      use dataio_pub,    only: printinfo
      use mpisetup,      only: master

      implicit none

      character(len=*), intent(in)  :: msg          !< message to be logged
      logical, optional, intent(in) :: stdout_cntrl !< flag to suppress printing to stdout

      logical            :: to_stdout

      to_stdout = master
      if (present(stdout_cntrl)) to_stdout = to_stdout .and. stdout_cntrl
      call printinfo(msg, to_stdout)

   end subroutine mg_write_log

!!$ ============================================================================
!>
!! \brief Emergency routine for quick ASCII dumps
!!
!! \details Absolute integer coordinates also allow seamless concatenation of dumps made by all PEs.
!<

   subroutine ascii_dump(filename)

      use constants,     only: xdim, ydim, zdim
      use dataio_pub,    only: msg
      use mpisetup,      only: master
      use multigridvars, only: base, plvl, ngridvars

      implicit none

      character(len=*), intent(in) :: filename !< name to write the emergency dump

      integer, parameter :: fu=30
      integer            :: i, j, k
      type(plvl), pointer :: curl

      if (.not. do_ascii_dump) return

      open(fu, file=filename, status="unknown")

      write(fu, '("#",a3,2a4,a6,10a20,/)')"i", "j", "k", "level", "x(i)", "y(j)", "z(k)", "source", "solution", "defect", "correction", "bx", "by", "bz"
      curl => base
      do while (associated(curl))
         do i = curl%is, curl%ie
            do j = curl%js, curl%je
               do k = curl%ks, curl%ke
                  write(fu, '(3i4,i6,10es20.11e3)')i-curl%is+curl%off(xdim), j-curl%js+curl%off(ydim), k-curl%ks+curl%off(zdim), &
                       &                           curl%level, curl%x(i), curl%y(j), curl%z(k), curl%mgvar(i, j, k, 1:ngridvars)
               enddo
               write(fu, '(/)')
            enddo
            write(fu, '(/,/)')
         enddo
         write(fu, '(/,/,/)')
         curl => curl%finer
      enddo

      close(fu)

      if (master) then
         write(msg,'(3a)') "[multigridhelpers:ascii_dump] Wrote dump '",filename,"'"
         call mg_write_log(msg)
      endif

   end subroutine ascii_dump

!!$ ============================================================================
!>
!! \brief Construct name of emergency ASCII dump
!<

   subroutine numbered_ascii_dump(basename, a)

      use dataio_pub, only: halfstep, msg
      use global,     only: nstep
      use mpisetup,   only: proc

      implicit none

      character(len=*), intent(in)  :: basename !< first part of the filename
      integer, optional, intent(in) :: a        !< additional number

      integer             :: l, n

      if (.not. do_ascii_dump) return

      n = 2 * nstep
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
