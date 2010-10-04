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
!    Initial implemetation of PIERNIK code was based on TVD split MHD code by
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

   ! namelist parameters
   logical            :: do_ascii_dump                      !< to dump, or not to dump: that is a question (ascii)
   logical            :: multidim_code_3D                   !< prefer code written for any 1D and 2D configuration even in 3D for benchmarking and debugging
   logical            :: dirty_debug                        !< Initialize everything with some insane values (dirtyH, defined below) and check if they can propagate
   integer            :: aux_par_I0, aux_par_I1, aux_par_I2 !< auxiliary integer parameters
   real               :: aux_par_R0, aux_par_R1, aux_par_R2 !< auxiliary real parameters

   real, parameter    :: dirtyH = 1e200, dirtyL = 1e50      !< If dirty_debug, initialize arrays with dirtyH and check if the solution contains anything above dirtyL

   character(len=256) :: str                                !< string for messages, filenames, etc.

contains

!/todo write also general set_dirty_array and check_dirty_array routines so there will be no need to use dirtyH and dirtyL outside multigridhelpers module.

!$ ============================================================================
!!
!! This routine pollutes selected multigrid variable with an insane value dirtyH.
!! If anything in the multigrid works by accident, through compiler-dependent initialization or unintentional relying on outdated values,
!! the insane value should pollute the solution in an easily visible way.
!!

   subroutine set_dirty(iv)

      use errh,          only: die
      use multigridvars, only: ngridvars, lvl, level_min, level_max

      implicit none

      integer, intent(in) :: iv   !< index of variable in lvl()%mgvar which we want to pollute

      integer :: l

      if (.not. dirty_debug) return

      if (iv < 1 .or. iv > ngridvars) call die("[multigridhelpers:set_dirty] Invalid variable index.")

      do l = level_min, level_max
         lvl(l)%mgvar(:, :, :, iv) = dirtyH
      end do

   end subroutine set_dirty

!!$ ============================================================================
!!
!! This routine checks for detectable traces of set_dirty calls.
!!

   subroutine check_dirty(lev, iv, label, expand)

      use errh,          only: die
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
      end if

      if (present(expand) .and. eff_dim==NDIM) then ! for 1D and 2D one should define ng_x,ng_y and ng_z
         if (expand > mg_nb) then
            ng = mg_nb
         else
            ng = expand
         end if
      else
         ng = 0
      end if

      do l = l1, l2
         do k = lvl(l)%ks-ng, lvl(l)%ke+ng
            do j = lvl(l)%js-ng, lvl(l)%je+ng
               do i = lvl(l)%is-ng, lvl(l)%ie+ng
                  if (abs(lvl(l)%mgvar(i, j, k, iv)) > dirtyL) &
                       write (*, '(3a,i4,a,i2,a,3(i3,a),i2,a,g20.12)') &
                       "[multigridhelpers:check_dirty] ", label, "@", proc, " lvl(", l, ")%mgvar(", i, ",", j, ",", k, ",", iv, ") = ", lvl(l)%mgvar(i, j, k, iv)
               end do
            end do
         end do
      end do

   end subroutine check_dirty

!!$ ============================================================================
!!
!! Assembles one-line log of V-cycle achevments
!!

   subroutine brief_v_log(v, norm)

      use multigridvars, only: vcycle_factors, cprefix, stdout
      use mpisetup,      only: proc

      implicit none

      integer, intent(in) :: v    !< number of V-cycles taken
      real, intent(in)    :: norm !< norm reduction factor achieved

      real              :: at
      integer           :: i, lm
      character(len=64) :: normred

      if (proc /= 0) return

      at = 0.
      if (v>0) at = sum(vcycle_factors(1:v,2))/v

      write(str, '(a,i3,1x,2a,f7.3,a,i3,a,f7.3,a,f11.9,a)') &
           "[multigrid] ", v, trim(cprefix), "Cycles, dt_wall=", vcycle_factors(0,2), " +", v, "*", at, ", norm/rhs= ", norm, " : "

      do i = 0, v
         if (vcycle_factors(i,1) < 1.0e4) then
            write(normred, '(f8.2)') vcycle_factors(i,1)
         else if (vcycle_factors(i,1) < 1.0e7) then
            write(normred, '(f8.0)') vcycle_factors(i,1)
         else
            write(normred, '(es8.2)') vcycle_factors(i,1)
         end if
         lm = len_trim(str)
         if (len(str) >= lm + 9) str(lm+2:lm+9) = normred(1:8)
      end do

      call mg_write_log(str, stdout)

   end subroutine brief_v_log

!!$ ============================================================================
!!
!! Logfile
!!

   subroutine mg_write_log(msg, stdout_cntrl)

      use errh,          only: printinfo
      use mpisetup,      only: proc

      implicit none

      character(len=*), intent(in)  :: msg          !< message to be logged
      logical, optional, intent(in) :: stdout_cntrl !< flag to supress printing to stdout

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

      use mpisetup,      only: proc
      use multigridvars, only: level_min, level_max, lvl, gb_cartmap, ngridvars, XDIR, YDIR, ZDIR

      implicit none

      character(len=*), intent(in) :: filename !< name to write the emergency dump

      integer, parameter :: fu=30
      integer            :: l, i, j, k

      if (.not. do_ascii_dump) return

      open(fu, file=filename, status="unknown")

      write(fu, '("#",a3,2a4,a6,7a20,/)')"i", "j", "k", "level", "x(i)", "y(j)", "z(k)", "source", "solution", "defect", "correction"
      do l = level_min, level_max
         do i = lvl(l)%is, lvl(l)%ie
            do j = lvl(l)%js, lvl(l)%je
               do k = lvl(l)%ks, lvl(l)%ke
                  write(fu, '(3i4,i6,7es20.11e3)')i-lvl(l)%is+gb_cartmap(proc)%proc(XDIR)*lvl(l)%nxb, &
                       &                          j-lvl(l)%js+gb_cartmap(proc)%proc(YDIR)*lvl(l)%nyb, &
                       &                          k-lvl(l)%ks+gb_cartmap(proc)%proc(ZDIR)*lvl(l)%nzb, &
                       &                          l, lvl(l)%x(i), lvl(l)%y(j), lvl(l)%z(k), lvl(l)%mgvar(i, j, k, 1:ngridvars)
               end do
               write(fu, '(/)')
            end do
            write(fu, '(/,/)')
         end do
         write(fu, '(/,/,/)')
      end do

      close(fu)

      if (proc == 0) then
         write(str, '(3a)')"[multigridhelpers:ascii_dump] Wrote dump '", filename, "'"
         call mg_write_log(str)
      endif

   end subroutine ascii_dump

!!$ ============================================================================
!!
!! Construct name of emergency ASCII dump
!!

   subroutine numbered_ascii_dump(basename, a)

      use mpisetup,      only: proc, nstep
      use dataio_public, only: halfstep

      implicit none

      character(len=*), intent(in)  :: basename !< first part of the filename
      integer, optional, intent(in) :: a        !< additional number

      integer             :: l, n

      if (.not. do_ascii_dump) return

      n = 2 * (nstep - 1)
      if (halfstep) n = n + 1

      if (present(a)) then
         write(str, '(a,i4,i6,i3)') trim(basename), proc, n, a
      else
         write(str, '(a,i4,i6)')    trim(basename), proc, n
      end if
      do l = 1, len_trim(str)
         if (str(l:l) == " ") str(l:l) = "_"
      end do
      call ascii_dump(trim(str))

   end subroutine numbered_ascii_dump

end module multigridhelpers
