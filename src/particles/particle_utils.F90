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

!>
!! \brief Diagnostics of particle module
!!
!! This module contains particle diagnostics related things that are needed by I/O and debug routines.
!<

module particle_utils
! pulled by NBODY

   implicit none

   private
   public :: max_pvel_1d, add_part_in_proper_cg, print_all_particles, is_part_in_cg, npf
   public :: max_pacc_3d, particle_diagnostics, twodtscheme, dump_diagnose, tot_energy, d_energy, tot_angmom, d_angmom, count_cg_particles, count_all_particles, global_count_all_particles, part_leave_cg

   real    :: tot_angmom           !< angular momentum of set of the particles
   real    :: tot_energy           !< total energy of set of the particles
   real    :: d_energy             !< error of energy of set of the particles in succeeding timesteps
   real    :: d_angmom             !< error of angular momentum in succeeding timesteps
   logical :: twodtscheme
   logical :: dump_diagnose        !< dump diagnose for each particle to a seperate log file
   integer(kind=4), parameter :: npf = 14  !< number of single particle fields
   integer(kind=4), parameter :: npb = 1   !< number of cells between in and phy or between phy and out boundaries

contains

   subroutine print_all_particles

      !use cg_leaves, only: leaves
      !use cg_list,   only: cg_list_element

      implicit none

      !type(cg_list_element), pointer :: cgl

      !!!!!!!!!Does nothing for now

      !cgl => leaves%first
      !do while (associated(cgl))
      !   call cgl%cg%pset%print
      !   cgl => cgl%nxt
      !enddo

   end subroutine print_all_particles

   subroutine max_pvel_1d(cg, max_v)

      use constants,      only: ndims, xdim, zdim, zero
      use grid_cont,      only: grid_container
      use particle_types, only: particle

      implicit none

      type(grid_container), pointer, intent(in)  :: cg
      type(particle), pointer                    :: pset
      real, dimension(ndims),        intent(out) :: max_v
      integer(kind=4)                            :: cdim
      real                                       :: v_tmp

      !Better way to do this? Was easier with arrays
      max_v = zero
      do cdim = xdim, zdim
         pset => cg%pset%first
         do while (associated(pset))
            v_tmp = abs(pset%pdata%vel(cdim))
            if (v_tmp > max_v(cdim)) max_v(cdim) = v_tmp
            pset => pset%nxt
         enddo
      enddo

   end subroutine max_pvel_1d

   subroutine max_pacc_3d(cg, max_pacc)

      use constants, only: big, CENTER, half, LO, xdim, zdim, zero
      use grid_cont, only: grid_container
      use mpisetup,  only: proc
      use types,     only: value
      use particle_types, only: particle

      implicit none

      type(grid_container), pointer, intent(in)  :: cg
      type(value),                   intent(out) :: max_pacc
      type(particle), pointer                    :: pset
      real                                       :: acc2, max_acc
      integer(kind=4)                            :: cdim

      max_pacc%assoc = big

      max_acc  = zero
      pset => cg%pset%first
      do while (associated(pset))
         acc2 = sum(pset%pdata%acc(:)**2)
         if (acc2 > max_acc) then
            max_acc = acc2
            max_pacc%coords(:) = pset%pdata%pos(:)
            !max_pacc%proc = i !> \todo it might be an information about extremum particle, but the scheme of log file is to print the process number
         endif
         pset => pset%nxt
      enddo
      max_pacc%val = sqrt(max_acc)
      do cdim = xdim, zdim
         max_pacc%loc(cdim) = int( half + (max_pacc%coords(cdim) - cg%coord(CENTER,cdim)%r(cg%ijkse(cdim, LO))) * cg%idl(cdim) )
      enddo
      max_pacc%proc = proc

   end subroutine max_pacc_3d

   subroutine particle_diagnostics(regular)

#ifdef VERBOSE
      use dataio_pub, only: msg, printinfo
      use mpisetup,   only: master
#endif /* VERBOSE */

      implicit none

      logical, intent(in) :: regular              !< govern shorter diagnostics
      real, save          :: init_energy          !< total initial energy of set of the particles
      real, save          :: init_angmom          !< initial angular momentum of set of the particles
      logical, save       :: first_run_lf = .true.

      if (dump_diagnose .and. regular) call dump_particles_to_textfile

#ifndef VERBOSE
      if (regular) return
#endif /* !VERBOSE */

      call get_angmom_totener(tot_angmom, tot_energy)

      if (first_run_lf) then
         init_energy = tot_energy
         init_angmom = tot_angmom
         first_run_lf = .false.
      endif
      d_energy = log_error(tot_energy, init_energy)
      d_angmom = log_error(tot_angmom, init_angmom)

#ifdef VERBOSE
      if (master) then
         write(msg,'(a,3(1x,e12.5))') '[particle_utils:particle_diagnostics] Total energy: initial, current, error ', init_energy, tot_energy, d_energy
         call printinfo(msg)
         write(msg,'(a,3(1x,e12.5))') '[particle_utils:particle_diagnostics] ang_momentum: initial, current, error ', init_angmom, tot_angmom, d_angmom
         call printinfo(msg)
      endif
#endif /* VERBOSE */

   end subroutine particle_diagnostics

   real function log_error(vcurr,vinit)

      use constants, only: zero
      use func,      only: operator(.equals.)

      implicit none

      real, intent(in) :: vcurr, vinit

      if (vinit .equals. zero) then
         log_error = vcurr
         return
      endif
      log_error = (vcurr - vinit) / vinit
      if (log_error .equals. zero) then
         log_error = zero
      else
         log_error = log(abs(log_error))
      endif

   end function log_error

   subroutine get_angmom_totener(ang_momentum, total_energy)

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: pSUM, xdim, ydim, zdim, zero
      use mpisetup,       only: piernik_MPI_Allreduce
      use particle_types, only: particle

      implicit none

      real, intent(out)              :: ang_momentum
      real, intent(out)              :: total_energy !< total energy of set of particles
      real                           :: L1, L2, L3
      type(cg_list_element), pointer :: cgl
      type(particle), pointer        :: pset

      ang_momentum = zero
      total_energy = zero

      cgl => leaves%first
      do while (associated(cgl))

         pset => cgl%cg%pset%first
         do while (associated(pset))
            associate( part => pset%pdata )
               L1 = part%pos(ydim) * part%vel(zdim) - part%pos(zdim) * part%vel(ydim)
               L2 = part%pos(zdim) * part%vel(xdim) - part%pos(xdim) * part%vel(zdim)
               L3 = part%pos(xdim) * part%vel(ydim) - part%pos(ydim) * part%vel(xdim)
               ang_momentum = ang_momentum + part%mass * sqrt(L1**2 + L2**2 + L3**2)
               total_energy = total_energy + part%energy
            end associate
            pset => pset%nxt
         enddo

         cgl => cgl%nxt
      enddo

      call piernik_MPI_Allreduce(total_energy, pSUM)
      call piernik_MPI_Allreduce(ang_momentum, pSUM)

   end subroutine get_angmom_totener

   subroutine is_part_in_cg(cg, pos, in, phy, out)

      use constants,     only: LO, HI, ndims, xdim, ydim, zdim, LEFT, RIGHT
      use domain,        only: dom
      use grid_cont,     only: grid_container
      use particle_func, only: particle_in_area

      implicit none

      type(grid_container), pointer, intent(in)  :: cg
      real, dimension(ndims),        intent(in)  :: pos
      logical,                       intent(out) :: in, phy, out
      real, dimension(ndims,2)                   :: bnd_in, bnd_out

      !There is probably a better way to write this
      bnd_out(:,LO) = [cg%coord(LEFT, xdim)%r(cg%ijkse(xdim,LO)-npb), cg%coord(LEFT, ydim)%r(cg%ijkse(ydim,LO)-npb), cg%coord(LEFT, zdim)%r(cg%ijkse(zdim,LO)-npb)]
      bnd_out(:,HI) = [cg%coord(RIGHT,xdim)%r(cg%ijkse(xdim,HI)+npb), cg%coord(RIGHT,ydim)%r(cg%ijkse(ydim,HI)+npb), cg%coord(RIGHT,zdim)%r(cg%ijkse(zdim,HI)+npb)]

      bnd_in(:,LO)  = [cg%coord(LEFT, xdim)%r(cg%ijkse(xdim,LO)+npb), cg%coord(LEFT, ydim)%r(cg%ijkse(ydim,LO)+npb), cg%coord(LEFT, zdim)%r(cg%ijkse(zdim,LO)+npb)]
      bnd_in(:,HI)  = [cg%coord(RIGHT,xdim)%r(cg%ijkse(xdim,HI)-npb), cg%coord(RIGHT,ydim)%r(cg%ijkse(ydim,HI)-npb), cg%coord(RIGHT,zdim)%r(cg%ijkse(zdim,HI)-npb)]

      in  = particle_in_area(pos, bnd_in)
      phy = particle_in_area(pos, cg%fbnd)
      out = particle_in_area(pos, bnd_out)   ! Ghost particle

      if (particle_in_area(pos, dom%edge)) return

      phy = cg_outside_dom(pos, cg%fbnd)
      if (phy) then
         in  = .true.
         out = .true.
      endif

   end subroutine is_part_in_cg

   logical function cg_outside_dom(pos, fbnd) result (phy)

      use constants, only: LO, HI, ndims, xdim, zdim
      use domain,    only: dom
      use func,      only: operator(.equals.), operator(.notequals.)

      implicit none

      real, dimension(ndims),       intent(in) :: pos
      real, dimension(ndims,LO:HI), intent(in) :: fbnd
      integer                                  :: cdim, k, count1, count2, count3
      integer, dimension(ndims)                :: tmp

      phy = .false.
      count2 = 0
      do cdim = xdim, zdim
         if (pos(cdim) < dom%edge(cdim, LO)) then
            count2 = count2 + 1
            tmp(cdim) = LO
         else if (pos(cdim) > dom%edge(cdim, HI)) then
            count2 = count2 + 1
            tmp(cdim) = HI
         else
            tmp(cdim) = 0
         endif
      enddo

      !CORNER
      if (count2 == 3) then
         count3 = 0
         do cdim = xdim, zdim
            if (fbnd(cdim,tmp(cdim)) .equals. dom%edge(cdim,tmp(cdim))) count3 = count3 + 1
         enddo
         if (count3 == 3) phy = .true.
         return

      !EDGE
      else if (count2 == 2) then
         do cdim = xdim, zdim
            if (tmp(cdim) /= 0) then
               if (fbnd(cdim,tmp(cdim)) .notequals. dom%edge(cdim,tmp(cdim))) return
            else
               if ( (pos(cdim) < fbnd(cdim,LO)) .or. (pos(cdim) > fbnd(cdim,HI)) ) return
            endif
         enddo
         phy = .true.
         return

      !FACE
      else
         do cdim = xdim, zdim
            if ( ((pos(cdim) < dom%edge(cdim,LO)) .and. (fbnd(cdim,LO) .equals. dom%edge(cdim,LO))) .or. (((pos(cdim) > dom%edge(cdim,HI)) .and. (fbnd(cdim,HI) .equals. dom%edge(cdim,HI)))) ) then
               count1 = 0
               do k = xdim, zdim
                  if (k /= cdim) then
                     if ( (pos(k) > fbnd(k,LO)) .and. (pos(k) < fbnd(k,HI)) ) then
                        count1 = count1 + 1
                     endif
                  endif
               enddo
               if (count1 == 2) then
                  phy = .true.
                  return
               endif
            endif

         enddo
      endif

   end function cg_outside_dom

   subroutine add_part_in_proper_cg(pid, mass, pos, vel, acc, ener, tform, tdyn, success)

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: ndims
      use mpisetup,       only: proc
#ifdef NBODY_CHECK_PID
      use particle_types, only: particle
#endif /* NBODY_CHECK_PID */

      implicit none

      integer(kind=4),        intent(in)  :: pid
      real, dimension(ndims), intent(in)  :: pos, vel
      real, dimension(ndims), intent(in)  :: acc
      real,                   intent(in)  :: ener
      real,                   intent(in)  :: mass
      real, optional,         intent(in)  :: tform
      real, optional,         intent(in)  :: tdyn
      logical, optional,      intent(out) :: success
      type(cg_list_element), pointer      :: cgl
      logical                             :: in, phy, out, cgfound
      real                                :: tform1, tdyn1
#ifdef NBODY_CHECK_PID
      type(particle), pointer             :: pset
#endif /* NBODY_CHECK_PID */

      tform1 = 0.0
      tdyn1  = 0.0
      if (present(tform)) tform1 = tform
      if (present(tdyn))  tdyn1  = tdyn
      cgfound = .false.

      cgl => leaves%first
      do while (associated(cgl))
         call is_part_in_cg(cgl%cg, pos, in, phy, out)
#ifdef NBODY_CHECK_PID
         if (phy .or. out) then
            pset => cgl%cg%pset%first
            do while (associated(pset))
               if (pset%pdata%pid == pid) then
                  phy = .false.
                  out = .false.
                  exit
               endif
               pset => pset%nxt
            enddo
         endif
#endif /* NBODY_CHECK_PID */
         if (phy .or. out) call cgl%cg%pset%add(pid, mass, pos, vel, acc, ener, in, phy, out, tform1, tdyn1)
         if (phy .or. out) write(*,*) proc, 'added: ', pid, int(pos), in, phy, out, int(cgl%cg%fbnd)
         cgfound = cgfound .or. (phy .or. out)
         cgl => cgl%nxt
      enddo
      if (present(success)) success = cgfound
      write(*,*) proc, 'adding particle: ', cgfound, pid, pos

   end subroutine add_part_in_proper_cg

   logical function attribute_to_proc(pset, j, cgdl, se) result(to_send)

      use cg_level_base,  only: base
      use constants,      only: I_ONE, LO, HI
      use domain,         only: dom
      use mpisetup,       only: proc
      use particle_func,  only: particle_in_area
      use particle_types, only: particle

      implicit none

      integer,                 intent(in) :: j
      real, dimension(:),      intent(in) :: cgdl
      type(particle), pointer, intent(in) :: pset
      integer, dimension(:,:), intent(in) :: se
      integer :: b

      to_send = .false.
      do b = lbound(base%level%dot%gse(j)%c(:), dim=1), ubound(base%level%dot%gse(j)%c(:), dim=1)
         if (particle_in_area(pset%pdata%pos, [dom%edge(:,LO) + (base%level%dot%gse(j)%c(b)%se(:,LO) - npb) * cgdl(:), dom%edge(:,LO) + (base%level%dot%gse(j)%c(b)%se(:,HI) + I_ONE + npb) * cgdl(:)])) then
            to_send = .true.
         else if (pset%pdata%outside) then
            if (cg_outside_dom(pset%pdata%pos, [dom%edge(:,LO) + base%level%dot%gse(j)%c(b)%se(:,LO) * cgdl(:), dom%edge(:,LO) + (base%level%dot%gse(j)%c(b)%se(:,HI) + I_ONE) * cgdl(:)])) then
               to_send = .true.
            endif
         endif
#ifdef NBODY_CHECK_PID
         if (to_send) then
            if (j == proc .and. all(base%level%dot%gse(j)%c(b)%se == se)) to_send = .false.
         endif
#endif /* NBODY_CHECK_PID */
         if (to_send) return
      enddo

      return
      if (.false. .and. proc == sum(se)) return

   end function attribute_to_proc

   ! Sends leaving particles between processors, and creates ghosts
   subroutine part_leave_cg()

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: I_ONE, PPP_PART
      use dataio_pub,     only: die
      use domain,         only: is_refined
      use grid_cont,      only: grid_container
      use MPIF,           only: MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_COMM_WORLD
      use MPIFUN,         only: MPI_Alltoall, MPI_Alltoallv
      use mpisetup,       only: proc, err_mpi, FIRST, LAST
      use ppp,            only: ppp_main
      use particle_types, only: particle

      implicit none

      integer(kind=4), dimension(FIRST:LAST) :: nsend, nrecv, counts, countr, disps, dispr
      integer                                :: i, j, ind, inc
      integer(kind=4)                        :: nchcg
      real, dimension(:), allocatable        :: part_send, part_recv, part_chcg
      type(cg_list_element), pointer         :: cgl
      type(grid_container),  pointer         :: cg
      type(particle), pointer                :: pset, pset2
      character(len=*), parameter            :: ts_label = "leave_cg"

      if (is_refined) call die("[particle_utils:part_leave_cg] AMR not implemented yet")

      call ppp_main%start(ts_label, PPP_PART)

      nsend = 0
      nrecv = 0

      !Count number of particles to be sent
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         do j = FIRST, LAST
            pset => cg%pset%first
            do while (associated(pset))
               if (.not. pset%pdata%in) then
                  ! TO CHECK: PARTICLES CHANGING CG OUTSIDE DOMAIN?
                  if (attribute_to_proc(pset, j, cg%dl, cg%ijkse)) nsend(j) = nsend(j) + I_ONE ! WON'T WORK in AMR!!!
               endif
               pset => pset%nxt
            enddo
         enddo
         cgl => cgl%nxt
      enddo
      nchcg = nsend(proc)
      nsend(proc) = 0

      !Exchange information about particles numbers to be sent / received
      call MPI_Alltoall(nsend, I_ONE, MPI_INTEGER, nrecv, I_ONE, MPI_INTEGER, MPI_COMM_WORLD, err_mpi)

      !Store data of particles to be sent
      allocate(part_send(sum(nsend(:))*npf), part_chcg(nchcg*npf))
      ind = 1
      inc = 1
      cgl => leaves%first
      do while (associated(cgl))
         cg => cgl%cg
         do j = FIRST, LAST
            pset => cg%pset%first
            do while (associated(pset))
               if (.not. pset%pdata%in) then
                  if (attribute_to_proc(pset, j, cg%dl, cg%ijkse)) then
                     if (j == proc) then
                        part_chcg(inc:inc+npf-1) = collect_single_part_fields(inc, pset%pdata)
                     else
                        part_send(ind:ind+npf-1) = collect_single_part_fields(ind, pset%pdata)
                     endif
                  endif
               endif
               pset => pset%nxt
            enddo
         enddo

         !Remove particles out of cg
         pset => cg%pset%first
         do while (associated(pset))
#ifdef NBODY_CHECK_PID
            if (.not. pset%pdata%out) then
#else /* !NBODY_CHECK_PID */
            if (.not. pset%pdata%in) then
#endif /* !NBODY_CHECK_PID */
               pset2 => pset%nxt
               call cg%pset%remove(pset)
               pset => pset2
               cycle
            endif
            pset => pset%nxt
         enddo

         cgl => cgl%nxt
      enddo

      !Send / receive particle data
      counts = npf * nsend
      allocate(part_recv(sum(nrecv) * npf))
      countr = npf * nrecv
      disps(FIRST) = 0
      dispr(FIRST) = 0
      do j = FIRST + I_ONE, LAST
         disps(j) = disps(j-1) + counts(j-1)
         dispr(j) = dispr(j-1) + countr(j-1)
      enddo

      call MPI_Alltoallv(part_send, counts, disps, MPI_DOUBLE_PRECISION, part_recv, countr, dispr, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, err_mpi)

      !Add particles in cgs
      ind = 1
      do j = FIRST, LAST
         if (nrecv(j) /= 0) then
            do i = 1, nrecv(j)
               call unpack_single_part_fields(ind, part_recv(ind:ind+npf-1))
            enddo
         endif
      enddo
      inc = 1
      if (nchcg /= 0) then
         do i = 1, nchcg
            call unpack_single_part_fields(inc, part_chcg(inc:inc+npf-1))
         enddo
      endif

      deallocate(part_send, part_recv, part_chcg)

      call ppp_main%stop(ts_label, PPP_PART)

   end subroutine part_leave_cg

   function collect_single_part_fields(ind, p) result(pinfo)

      use particle_types, only: particle_data

      implicit none

      real, dimension(npf)               :: pinfo
      integer,             intent(inout) :: ind
      type(particle_data), intent(in)    :: p

      pinfo(1)    = p%pid
      if (.not. (pinfo(1) >=1)) print *, 'error, id cannot be zero', ind, pinfo(1)
      pinfo(2)    = p%mass
      pinfo(3:5)  = p%pos
      pinfo(6:8)  = p%vel
      pinfo(9:11) = p%acc
      pinfo(12)   = p%energy
      pinfo(13)   = p%tform
      pinfo(14)   = p%tdyn

      ind = ind + npf

   end function collect_single_part_fields

   subroutine unpack_single_part_fields(ind, pinfo)

      use constants, only: ndims

      implicit none

      integer,              intent(inout) :: ind
      real, dimension(npf), intent(in)    :: pinfo

      integer(kind=4)                     :: pid
      real, dimension(ndims)              :: pos, vel, acc
      real                                :: mass, ener, tform, tdyn
      logical                             :: attributed

      pid   = nint(pinfo(1), kind=4)
      mass  = pinfo(2)
      pos   = pinfo(3:5)
      vel   = pinfo(6:8)
      acc   = pinfo(9:11)
      ener  = pinfo(12)
      tform = pinfo(13)
      tdyn  = pinfo(14)
      call add_part_in_proper_cg(pid, mass, pos, vel, acc, ener, tform, tdyn, attributed) ! TO DO IN AMR USE GRID_ID TO CUT THE SEARCH SHORT
      if (.not. attributed) print *, 'error, particle', pid, 'cannot be attributed!' ! NON-AMR ONLY
      ind = ind + npf

   end subroutine unpack_single_part_fields

   integer(kind=4) function count_cg_particles(cg) result(n_part)

      use grid_cont,      only: grid_container
      use particle_types, only: particle

      implicit none

      type(grid_container), pointer, intent(in) :: cg
      type(particle),       pointer             :: pset

      n_part = 0
      pset => cg%pset%first
      do while (associated(pset))
         if (pset%pdata%phy) n_part = n_part + 1
      pset => pset%nxt
      enddo

   end function count_cg_particles

   integer(kind=4) function count_all_particles() result(pcount)

      use cg_leaves,      only: leaves
      use cg_list,        only: cg_list_element
      use constants,      only: I_ONE
      use particle_types, only: particle

      implicit none

      type(cg_list_element), pointer :: cgl
      type(particle), pointer    :: pset

      pcount = 0
      cgl => leaves%first
      do while (associated(cgl))
         pset => cgl%cg%pset%first
         do while (associated(pset))
            if (pset%pdata%phy) pcount = pcount + I_ONE
            pset => pset%nxt
         enddo
         cgl => cgl%nxt
      enddo

   end function count_all_particles

   integer(kind=4) function global_count_all_particles() result(gpcount)

      use constants, only: pSUM
      use mpisetup,  only: piernik_MPI_Allreduce

      implicit none

      gpcount = count_all_particles()
      call piernik_MPI_Allreduce(gpcount, pSUM)

   end function global_count_all_particles

   subroutine dump_particles_to_textfile

      use cg_leaves,        only: leaves
      use cg_list,          only: cg_list_element
      use constants,        only: cwdlen, half, ndims
      use dataio_pub,       only: log_wr, nrestart, problem_name, run_id
      use global,           only: t, dt
      use particle_gravity, only: get_acc_model
      use particle_types,   only: particle

      implicit none

      real                           :: kdt
      real, dimension(ndims)         :: acc2
      integer                        :: i, lun_out
      character(len=cwdlen)          :: plog_file
      type(cg_list_element), pointer :: cgl
      type(particle), pointer        :: pset

      kdt = dt ; if (.not.twodtscheme) kdt = half * dt

      write(plog_file,'(6a,i3.3,a)') trim(log_wr),'/',trim(problem_name),'_',trim(run_id),'_',nrestart,'_out.log'
      open(newunit=lun_out, file=plog_file, status='unknown',  position='append')
      cgl => leaves%first
      do while (associated(cgl))

         pset => cgl%cg%pset%first
         do while (associated(pset))
            associate( part => pset%pdata )
               call get_acc_model(part%pos, part%mass, acc2)
               write(lun_out, '(a,I3.3,1X,19(E13.6,1X))') 'particle', i, t, kdt, part%mass, part%pos, part%vel, part%acc, acc2(:), part%energy
            end associate
            pset => pset%nxt
         enddo

         cgl => cgl%nxt
      enddo
      close(lun_out)

   end subroutine dump_particles_to_textfile

end module particle_utils
