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
!! \brief Utilities for particle module
!!
!! This module contains particle utilties that are needed by other particle routines.
!<

module particle_utils
! pulled by NBODY

   implicit none

   private
   public :: add_part_in_proper_cg, is_part_in_cg, part_refresh_ghosts, part_leave_cg, detach_particle, &
        &    count_all_particles, global_count_all_particles, global_balance_particles, part_refresh_flags

contains

   subroutine is_part_in_cg(cg, pos, indomain, in, phy, out, fin)

      use constants,     only: LO, HI, ndims, xdim, ydim, zdim
      use domain,        only: dom
      use grid_cont,     only: grid_container
      use particle_func, only: particle_in_area, ijk_of_particle

      implicit none

      type(grid_container), pointer, intent(in)  :: cg
      real, dimension(ndims),        intent(in)  :: pos
      logical,                       intent(in)  :: indomain
      logical,                       intent(out) :: in, phy, out, fin
      integer(kind=8), dimension(ndims)          :: ijk

      in  = particle_in_area(pos, cg%bnd_in)
      phy = particle_in_area(pos, cg%fbnd)
      out = particle_in_area(pos, cg%bnd_out)   ! Ghost particle

      ijk = ijk_of_particle(pos, dom%edge(:,LO), cg%idl)
      ijk = max(ijk, int(cg%ijkse(:,LO), kind=8))
      ijk = min(ijk, int(cg%ijkse(:,HI), kind=8))
      fin = cg%leafmap(ijk(xdim), ijk(ydim), ijk(zdim))

      if (indomain) return

      if (outdom_part_in_cg(pos, cg%fbnd, cg%ext_bnd)) then
         in  = .true.
         phy = .true.
         out = .true.
         ijk = ijk_of_particle(pos, dom%edge(:,LO), cg%idl)
         ijk = max(ijk, cg%l%off)
         ijk = min(ijk, cg%l%off + cg%l%n_d - 1)
         fin = cg%leafmap(ijk(xdim), ijk(ydim), ijk(zdim))
      endif

   end subroutine is_part_in_cg

   logical function outdom_part_in_cg(pos, fbnd, ext_bnd) result (phy)

      use constants, only: LO, HI, ndims, xdim, zdim
      use domain,    only: dom

      implicit none

      real,    dimension(ndims),        intent(in) :: pos
      real,    dimension(ndims, LO:HI), intent(in) :: fbnd
      logical, dimension(ndims, LO:HI), intent(in) :: ext_bnd
      logical, dimension(ndims)                    :: fulfilled
      integer                                      :: cdim

      do cdim = xdim, zdim
         if (pos(cdim) < dom%edge(cdim, LO)) then
            fulfilled(cdim) = ext_bnd(cdim, LO)
         else if (pos(cdim) >= dom%edge(cdim, HI)) then
            fulfilled(cdim) = ext_bnd(cdim, HI)
         else
            fulfilled(cdim) = pos(cdim) >= fbnd(cdim, LO) .and. pos(cdim) < fbnd(cdim, HI)
         endif
      enddo
      phy = all(fulfilled)

   end function outdom_part_in_cg

   subroutine add_part_in_proper_cg(pid, mass, pos, vel, acc, ener, tform, tdyn, success)

      ! use cg_cost_data,   only: I_PARTICLE
      use cg_list,        only: cg_list_element
      use cg_list_global, only: all_cg
      use constants,      only: ndims, base_level_id
      use domain,         only: dom
      use particle_func,  only: particle_in_area
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
      logical                             :: in, phy, out, fin, indomain, cgfound, toadd, met, cgmet
      real                                :: tform1, tdyn1
#ifdef NBODY_CHECK_PID
      type(particle), pointer             :: pset
#endif /* NBODY_CHECK_PID */

      tform1 = 0.0
      tdyn1  = 0.0
      if (present(tform)) tform1 = tform
      if (present(tdyn))  tdyn1  = tdyn
      cgfound = .false.
      cgmet = .false.

      ! ToDo: OPT: Precompute list of possible cg using SFC_id, use direct loop as a fallback, perhaps with warning.
      indomain = particle_in_area(pos, dom%edge)
      cgl => all_cg%first
      do while (associated(cgl))
         if (cgl%cg%l%id >= base_level_id) then
            ! call cgl%cg%costs%start  ! Too big overhead when we call it for every particle
            call is_part_in_cg(cgl%cg, pos, indomain, in, phy, out, fin)
            toadd = out .and. fin
            met = out
#ifdef NBODY_CHECK_PID
            if (toadd) then
               pset => cgl%cg%pset%first
               do while (associated(pset))
                  if (pset%pdata%pid == pid) then
                     toadd = .false.
                     met = .false.
                     exit
                  endif
                  pset => pset%nxt
               enddo
            endif
#endif /* NBODY_CHECK_PID */
            if (toadd) call cgl%cg%pset%add(pid, mass, pos, vel, acc, ener, in, phy, out, fin, tform1, tdyn1)
            cgfound = cgfound .or. toadd
            cgmet = cgmet .or. met
            ! call cgl%cg%costs%stop(I_PARTICLE, ppp_c = .false.)
            ! There are too many calls to include this contribution as cg_cost:particles in the PPP output.
            ! It will be covered by add_part cumulative counter in part_leave_cg() instead.
            ! Collecting of the cg costs for load balancing putposes will still work.
         endif
         cgl => cgl%nxt
      enddo
      !if (present(success)) success = cgfound
      if (present(success)) success = cgmet

   end subroutine add_part_in_proper_cg

   function fbnd_npb(se, off, dl, nexp) result(fbnd)

      use constants, only: I_ONE, LO, HI, ndims
      use domain,    only: dom

      implicit none

      integer(kind=8), dimension(ndims, LO:HI), intent(in) :: se
      integer(kind=8), dimension(ndims),        intent(in) :: off
      real,            dimension(ndims),        intent(in) :: dl
      integer(kind=4),                          intent(in) :: nexp
      real,            dimension(ndims, LO:HI)             :: fbnd

      fbnd(:, LO)  = dom%edge(:, LO) + dl(:) * (se(:, LO)         - off(:) - nexp)
      fbnd(:, HI)  = dom%edge(:, LO) + dl(:) * (se(:, HI) + I_ONE - off(:) + nexp)

   end function fbnd_npb

   logical function attribute_to_proc(pset, j, se) result(to_send)

      use cg_level_base,      only: base
      use cg_level_connected, only: cg_level_connected_t
      use constants,          only: LO, HI, ndims, xdim, zdim, I_ZERO
      use domain,             only: dom
      use mpisetup,           only: proc
      use particle_func,      only: particle_in_area
      use particle_types,     only: particle, npb

      implicit none

      type(particle), pointer,         intent(in) :: pset
      integer,                         intent(in) :: j
      integer(kind=4), dimension(:,:), intent(in) :: se

      integer                             :: b
      integer(kind=4)                     :: cdim
      real,    dimension(ndims)           :: ldl
      logical, dimension(ndims, LO:HI)    :: ext_bnd
      type(cg_level_connected_t), pointer :: ll

      to_send = .false.
      if (pset%pdata%in .and. pset%pdata%fin) return

      ll => base%level
      do while (associated(ll))
         ldl(:) = dom%L_(:) / ll%l%n_d(:)
         do b = lbound(ll%dot%gse(j)%c(:), dim=1), ubound(ll%dot%gse(j)%c(:), dim=1)
            if (particle_in_area(pset%pdata%pos, fbnd_npb(ll%dot%gse(j)%c(b)%se, ll%l%off, ldl, npb))) then
               to_send = .true.
            else if (pset%pdata%outside) then
               do cdim = xdim, zdim
                  ext_bnd(cdim, LO) = ll%l%is_ext_bnd(ll%dot%gse(j)%c(b)%se, cdim, LO)
                  ext_bnd(cdim, HI) = ll%l%is_ext_bnd(ll%dot%gse(j)%c(b)%se, cdim, HI)
               enddo
               if (outdom_part_in_cg(pset%pdata%pos, fbnd_npb(ll%dot%gse(j)%c(b)%se, ll%l%off, ldl, I_ZERO), ext_bnd)) to_send = .true.
            endif
#ifdef NBODY_CHECK_PID
            if (to_send) then
               if (j == proc .and. all(ll%dot%gse(j)%c(b)%se == se)) to_send = .false.
            endif
#endif /* NBODY_CHECK_PID */
            if (to_send) return
         enddo
         ll => ll%finer
      enddo


      return
      if (.false. .and. proc == sum(se)) return

   end function attribute_to_proc

!>
!! \brief Refresh particle flags
!! \details This may be necessary after refinement update or when the domain is expanded.
!<

   subroutine part_refresh_flags

      use cg_cost_data,   only: I_PARTICLE
      use cg_list,        only: cg_list_element
      use cg_list_global, only: all_cg
      use domain,         only: dom
      use particle_func,  only: particle_in_area
      use particle_types, only: particle

      implicit none

      type(cg_list_element), pointer :: cgl
      type(particle), pointer :: pset

      cgl => all_cg%first
      do while (associated(cgl))
         call cgl%cg%costs%start
         pset => cgl%cg%pset%first
         do while (associated(pset))
            call is_part_in_cg(cgl%cg, pset%pdata%pos, particle_in_area(pset%pdata%pos, dom%edge), pset%pdata%in, pset%pdata%phy, pset%pdata%out, pset%pdata%fin)
            pset => pset%nxt
         enddo
         call cgl%cg%costs%stop(I_PARTICLE)
         cgl => cgl%nxt
      enddo

   end subroutine part_refresh_flags

!> Remove ghosts and reassign all physical particles. Useful for sanitizing after a change of refinement.

   subroutine part_refresh_ghosts

      use cg_cost_data,   only: I_PARTICLE
      use cg_list,        only: cg_list_element
      use cg_list_global, only: all_cg
      use particle_types, only: particle

      implicit none

      type(cg_list_element), pointer :: cgl
      type(particle), pointer :: pset

      cgl => all_cg%first
      do while (associated(cgl))
         call cgl%cg%costs%start
         pset => cgl%cg%pset%first
         do while (associated(pset))
            if (pset%pdata%phy) then
               pset => pset%nxt
            else !Remove ghosts
               call detach_particle(cgl%cg, pset)
            endif
         enddo
         call cgl%cg%costs%stop(I_PARTICLE)
         cgl => cgl%nxt
      enddo
      call part_leave_cg()

   end subroutine part_refresh_ghosts

   ! Sends leaving particles between processors, and creates ghosts
   ! OPT: It looks like operations here will cost about O(proc * number_of_cg * number_of_particles).
   !      At least O(proc * number_of_cg) can be reduced by efficient use of SFC_id and oct-trees.
   subroutine part_leave_cg()

      use cg_cost_data,   only: I_PARTICLE
      use cg_list,        only: cg_list_element
      use cg_list_global, only: all_cg
      use constants,      only: I_ONE, PPP_PART, base_level_id
      use grid_cont,      only: grid_container
      use MPIF,           only: MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_COMM_WORLD
      use MPIFUN,         only: MPI_Alltoall, MPI_Alltoallv
      use mpisetup,       only: proc, err_mpi, FIRST, LAST
      use ppp,            only: ppp_main
      use particle_types, only: particle, npf

      implicit none

      integer(kind=4), dimension(FIRST:LAST) :: nsend, nrecv, counts, countr, disps, dispr
      integer                                :: i, j, ind, inc
      integer(kind=4)                        :: nchcg
      real, dimension(:), allocatable        :: part_send, part_recv, part_chcg
      type(cg_list_element), pointer         :: cgl
      type(grid_container),  pointer         :: cg
      type(particle), pointer                :: pset
      character(len=*), parameter            :: ts_label = "leave_cg", cnt_label = "cnt_part", snd_label = "send_part_prep", &
           &                                    del_label = "detach_part", add_label = "add_part", addl_label = "chcg_part"

      call ppp_main%start(ts_label, PPP_PART)

      nsend = 0
      nrecv = 0

      ! ToDo: OPT: Precompute list of possible cg and processes using SFC_id, use direct loop as a fallback, perhaps with warning.
      ! It would be useful to have browsable oct-tree of cg in dot or somewhere.
      !Count number of particles to be sent
      call ppp_main%start(cnt_label, PPP_PART)
      do j = FIRST, LAST
         cgl => all_cg%first
         do while (associated(cgl))
            if (cgl%cg%l%id >= base_level_id) then
               call cgl%cg%costs%start
               cg => cgl%cg
               pset => cg%pset%first
               do while (associated(pset))
                  ! TO CHECK: PARTICLES CHANGING CG OUTSIDE DOMAIN?
                  if (attribute_to_proc(pset, j, cg%ijkse)) nsend(j) = nsend(j) + I_ONE ! WON'T WORK in AMR!!!
                  pset => pset%nxt
               enddo
               call cgl%cg%costs%stop(I_PARTICLE)
            endif
            cgl => cgl%nxt
         enddo
      enddo
      nchcg = nsend(proc)
      nsend(proc) = 0
      call ppp_main%stop(cnt_label, PPP_PART)

      !Store data of particles to be sent
      call ppp_main%start(snd_label, PPP_PART)
      allocate(part_send(sum(nsend) * npf), part_chcg(nchcg * npf))
      ind = 1
      inc = 1
      do j = FIRST, LAST
         cgl => all_cg%first
         do while (associated(cgl))
            if (cgl%cg%l%id >= base_level_id) then
               call cgl%cg%costs%start
               cg => cgl%cg
               pset => cg%pset%first
               do while (associated(pset))
                  if (attribute_to_proc(pset, j, cg%ijkse)) then
                     if (j == proc) then
                        part_chcg(inc:inc+npf-1) = collect_single_part_fields(inc, pset%pdata)
                     else
                        part_send(ind:ind+npf-1) = collect_single_part_fields(ind, pset%pdata)
                     endif
                  endif
                  pset => pset%nxt
               enddo
               call cgl%cg%costs%stop(I_PARTICLE)
            endif
            cgl => cgl%nxt
         enddo
      enddo
      call ppp_main%stop(snd_label, PPP_PART)

      !Remove particles out of cg
      call ppp_main%start(del_label, PPP_PART)
      cgl => all_cg%first
      do while (associated(cgl))
         if (cgl%cg%l%id >= base_level_id) then
            call cgl%cg%costs%start
            cg => cgl%cg
            pset => cg%pset%first
            do while (associated(pset))
#ifdef NBODY_CHECK_PID
               if (.not. pset%pdata%out .or. (pset%pdata%in .and. .not. pset%pdata%fin)) then
#else /* !NBODY_CHECK_PID */
               if (.not. (pset%pdata%in .and. pset%pdata%fin)) then
#endif /* !NBODY_CHECK_PID */
                  call detach_particle(cg, pset)
                  cycle
               endif
               pset => pset%nxt
            enddo
            call cgl%cg%costs%stop(I_PARTICLE)
         endif
         cgl => cgl%nxt
      enddo
      call ppp_main%stop(del_label, PPP_PART)

      !Exchange information about particles numbers to be sent / received
      call MPI_Alltoall(nsend, I_ONE, MPI_INTEGER, nrecv, I_ONE, MPI_INTEGER, MPI_COMM_WORLD, err_mpi)

      !Send / receive particle data
      allocate(part_recv(sum(nrecv) * npf))
      counts = npf * nsend
      countr = npf * nrecv
      disps(FIRST) = 0
      dispr(FIRST) = 0
      do j = FIRST + I_ONE, LAST
         disps(j) = disps(j-1) + counts(j-1)
         dispr(j) = dispr(j-1) + countr(j-1)
      enddo

      call MPI_Alltoallv(part_send, counts, disps, MPI_DOUBLE_PRECISION, part_recv, countr, dispr, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, err_mpi)

      call ppp_main%start(add_label, PPP_PART)
      !Add particles in cgs
      ind = 1
      do j = FIRST, LAST
         if (nrecv(j) /= 0) then
            do i = 1, nrecv(j)
               call unpack_single_part_fields(ind, part_recv(ind:ind+npf-1))
            enddo
         endif
      enddo
      call ppp_main%stop(add_label, PPP_PART)
      inc = 1
      if (nchcg /= 0) then
         call ppp_main%start(addl_label, PPP_PART)
         do i = 1, nchcg
            call unpack_single_part_fields(inc, part_chcg(inc:inc+npf-1))
         enddo
         call ppp_main%stop(addl_label, PPP_PART)
      endif

      deallocate(part_send, part_recv, part_chcg)

      call ppp_main%stop(ts_label, PPP_PART)

   end subroutine part_leave_cg

   function collect_single_part_fields(ind, p) result(pinfo)

      use particle_types, only: particle_data, P_ID, P_MASS, P_POS_X, P_POS_Z, P_VEL_X, P_VEL_Z, P_ACC_X, P_ACC_Z, P_ENER, P_TFORM, P_TDYN, npf

      implicit none

      real, dimension(npf)               :: pinfo
      integer,             intent(inout) :: ind
      type(particle_data), intent(in)    :: p

      pinfo(P_ID)            = p%pid
      if (.not. (pinfo(P_ID) >=1)) print *, 'error, id cannot be zero', ind, pinfo(P_ID)
      pinfo(P_MASS)          = p%mass
      pinfo(P_POS_X:P_POS_Z) = p%pos
      pinfo(P_VEL_X:P_VEL_Z) = p%vel
      pinfo(P_ACC_X:P_ACC_Z) = p%acc
      pinfo(P_ENER)          = p%energy
      pinfo(P_TFORM)         = p%tform
      pinfo(P_TDYN)          = p%tdyn

      ind = ind + npf

   end function collect_single_part_fields

   subroutine unpack_single_part_fields(ind, pinfo)

      use constants,      only: ndims
      use particle_types, only: P_ID, P_MASS, P_POS_X, P_POS_Z, P_VEL_X, P_VEL_Z, P_ACC_X, P_ACC_Z, P_ENER, P_TFORM, P_TDYN, npf
      implicit none

      integer,              intent(inout) :: ind
      real, dimension(npf), intent(in)    :: pinfo

      integer(kind=4)                     :: pid
      real, dimension(ndims)              :: pos, vel, acc
      real                                :: mass, ener, tform, tdyn
      logical                             :: attributed

      pid   = nint(pinfo(P_ID), kind=4)
      mass  = pinfo(P_MASS)
      pos   = pinfo(P_POS_X:P_POS_Z)
      vel   = pinfo(P_VEL_X:P_VEL_Z)
      acc   = pinfo(P_ACC_X:P_ACC_Z)
      ener  = pinfo(P_ENER)
      tform = pinfo(P_TFORM)
      tdyn  = pinfo(P_TDYN)
      call add_part_in_proper_cg(pid, mass, pos, vel, acc, ener, tform, tdyn, attributed) ! TO DO IN AMR USE GRID_ID TO CUT THE SEARCH SHORT
      if (.not. attributed) print *, 'error, particle', pid, 'cannot be attributed! ', pos ! NON-AMR CHECK ONLY
      ind = ind + npf

   end subroutine unpack_single_part_fields

!> \brief Count all particles on current thread

   integer(kind=4) function count_all_particles() result(pcount)

      use cg_cost_data,   only: I_PARTICLE
      use cg_list,        only: cg_list_element
      use cg_list_global, only: all_cg
      use constants,      only: base_level_id
      use dataio_pub,     only: die

      implicit none

      type(cg_list_element), pointer :: cgl

      pcount = 0
      cgl => all_cg%first
      do while (associated(cgl))
         call cgl%cg%costs%start

         if (associated(cgl%cg%pset%first) .and. (cgl%cg%l%id < base_level_id)) &
              call die("[particle_utils:count_all_particles] particles found below base level")
         pcount = pcount + cgl%cg%count_particles()

         call cgl%cg%costs%stop(I_PARTICLE)
         cgl => cgl%nxt
      enddo

   end function count_all_particles

!> \brief Count all particles on all threads

   integer(kind=8) function global_count_all_particles(message) result(gpcount)

      use allreduce,  only: piernik_MPI_Allreduce
      use constants,  only: pSUM, fmt_len, V_INFO
      use dataio_pub, only: printinfo, msg
      use mpisetup,   only: master

      implicit none

      character(len=*), optional, intent(in) :: message

      character(len=fmt_len) :: fmt

      gpcount = count_all_particles()
      call piernik_MPI_Allreduce(gpcount, pSUM)

      if (present(message)) then
         if (master) then
            select case (gpcount)
               case (:999999999)
                  fmt = '(2a,i9)'
               case default
                  fmt = '(2a,i20)'
            end select
            write(msg, fmt) trim(message), " Total number of particles is ", gpcount
            call printinfo(msg, V_INFO)
         endif
      endif

   end function global_count_all_particles

!> \brief Estimate load balance for particles

   real function global_balance_particles() result(lb_particles)

      use allreduce, only: piernik_MPI_Allreduce
      use constants, only: pMAX, pSUM, one
      use mpisetup,  only: nproc

      implicit none

      integer(kind=8) :: l_part_cnt, max_part_cnt, g_part_cnt

      l_part_cnt = count_all_particles()
      max_part_cnt = l_part_cnt
      call piernik_MPI_Allreduce(max_part_cnt, pMAX)
      if (max_part_cnt > 0.) then
         g_part_cnt = l_part_cnt
         call piernik_MPI_Allreduce(g_part_cnt, pSUM)
         lb_particles = g_part_cnt/(nproc * real(max_part_cnt))
      else
         lb_particles = one  ! 0/0 is balanced
      endif

   end function global_balance_particles

!> \brief Remove particle from cg%pset

   subroutine detach_particle(cg, part)

      use grid_cont,      only: grid_container
      use particle_types, only: particle

      implicit none

      type(grid_container),  pointer, intent(inout) :: cg
      type(particle),        pointer, intent(inout) :: part
      type(particle),        pointer                :: part2

      part2 => part%nxt
      call cg%pset%remove(part)
      part => part2

   end subroutine detach_particle

end module particle_utils
